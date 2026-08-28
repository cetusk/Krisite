// Krisite — 中核（`SPEC-phase3.md` §5 / §14 の CP2）
//
//     boolean : PolySoup × PolySoup → PolySoup     ★ CSG について閉じる
//
// **CP2 の範囲**は「入口と出口の分離 + 分類の台の移動」です（§14）。
// 分類は WNV ではなく**真偽値版**で、断片は各 source への内外ビットで分類します。
// 局所 BSP はまだ入れず、過剰分割のままです（CP4）。
//
// ---
//
// ## なぜ分類の台を移すのか
//
// 現在の二項実装は「断片の代表点から**相手メッシュ**へレイキャストする」形で、
// **相手が整数座標の三角メッシュであること**に依存しています。連鎖の 2 段目では
// 相手が構成点のスープになり、レイの平面のオフセットが有理数になるので
// **固定幅で表せません**（EMBER §3.2 と同じ壁）。
//
// **台を生成 0 の整数メッシュに移すと、何段連鎖してもレイキャストが厳密に行えます。**
//
// > **代償: 分類が一時的に「より大域的」になります。** 相手 1 枚へのレイキャストから
// > sources 全部へのレイキャストに変わるので、**連鎖の段数に比例してコストが増えます。**
// > **これは通過点であって目的地ではありません。** CP3 のセグメントトレース（§5.5）が
// > 参照点の伝播で局所化します。
#ifndef KRISITE_CSG_SOUP_BOOLEAN_HPP
#define KRISITE_CSG_SOUP_BOOLEAN_HPP

#include <algorithm>
#include <array>
#include <map>
#include <vector>
#if defined(KRISITE_DEBUG_SOUP)
#include <cstdio>
#endif

#include "krisite/csg/boolean.hpp"  // BoolOp / BoolOptions / BoolStats
#include "krisite/csg/interior.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/trace.hpp"
#include "krisite/octree/adaptive.hpp"

namespace krisite::csg {

namespace detail {

/// 平面表を移し替える。**係数は変わりません**（`PlaneTable` は代表を正準化するだけ）。
inline std::vector<PlaneRef> remap_planes(const PlaneTable& from, PlaneTable& to) {
    std::vector<PlaneRef> m(from.size());
    for (std::size_t i = 0; i < from.size(); ++i) {
        m[i] = to.intern(from.at(static_cast<PlaneId>(i)));
    }
    return m;
}

}  // namespace detail

/// $n$ 項のブール演算（`SPEC-phase3.md` §1 の契約）。**出力もスープなので連鎖できます。**
///
/// 入力スープ $i$ の多角形は **WNV の成分 $i$** を動かします（§5.1）。
/// 出力は再び **1 成分**のスープです（自分の多角形が自分の曲面）。
inline PolySoup boolean_nary(const std::vector<const PolySoup*>& inputs, const Indicator& ind,
                             const BoolOptions& opt, BoolStats* stats = nullptr) {
    KRISITE_CHECK(!inputs.empty(), "boolean_nary: 入力が空");
    BoolStats st;
    PointCache point_cache;
    PointCache* const cache = opt.cache_points ? &point_cache : nullptr;

    PolySoup out;
    out.components = 1;
    out.indicator = indicator_source(0);
    const std::size_t n_comp = inputs.size();

    // ---- 1. 平面表をまとめ、多角形を移す ------------------------------------
    std::vector<Poly> polys;
    for (std::size_t i = 0; i < n_comp; ++i) {
        const PolySoup& s = *inputs[i];
        const std::vector<PlaneRef> m = detail::remap_planes(s.table, out.table);
        for (const Poly& q : s.polys) {
            Poly r = q;
            r.frag.support = m[q.frag.support].id;
            r.frag.flipped = (q.frag.flipped != m[q.frag.support].flipped);
            for (PlaneId& e : r.frag.edge) e = m[e].id;
            r.comp = static_cast<std::uint32_t>(i);
            r.frag.owner = static_cast<int>(i);
            polys.push_back(std::move(r));
        }
    }

    // 分割に使う平面（過剰分割。CP4 で局所 BSP に置き換わる）
    std::vector<PlaneId> all_split;
    for (const Poly& q : polys) all_split.push_back(q.frag.support);
    std::sort(all_split.begin(), all_split.end());
    all_split.erase(std::unique(all_split.begin(), all_split.end()), all_split.end());

    // ---- 2. トレースの対象（多角形そのものが曲面）----------------------------
    std::vector<TracePoly> tps;
    std::vector<octree::Aabb> tps_aabb;
    tps.reserve(polys.size());
    tps_aabb.reserve(polys.size());
    for (const Poly& q : polys) {
        TracePoly tp;
        tp.support = q.frag.support;
        tp.edge = q.frag.edge;
        tp.comp = q.comp;
        tp.orient = q.frag.flipped ? -1 : +1;
        // 辺ごとの内側の符号（載っていない頂点の符号）
        const std::size_t nv = vertex_count(q.frag);
        std::vector<geom::HPointD> vs(nv);
        for (std::size_t i = 0; i < nv; ++i) vs[i] = fragment_vertex(out.table, q.frag, i, cache);
        tp.inward.resize(q.frag.edge.size());
        for (std::size_t k = 0; k < q.frag.edge.size(); ++k) {
            std::int8_t sgn = 0;
            for (const geom::HPointD& v : vs) {
                const int sv = geom::side(out.table.at(q.frag.edge[k]), v);
                if (sv != 0) {
                    sgn = static_cast<std::int8_t>(sv);
                    break;
                }
            }
            tp.inward[k] = sgn;
        }
        tps.push_back(std::move(tp));
        tps_aabb.push_back(q.aabb);
    }

    // ---- 3. 根の参照点（**全多角形の外**。だから WNV は 0）------------------------
    //
    // 以降、分割のたびに子の箱へ**伝播**させます（§5.3.1）。
    // **伝播は最適化ではなく成功率の前提です。** 大域の参照点から辿ると経路が長くなり、
    // 軸平行な入力では途中で必ず辺や頂点に当たります。
    geom::IPoint root_ref{};
    {
        octree::Aabb all{};
        for (int k = 0; k < 3; ++k) {
            all.lo[k] = krisite::kCoordMax;
            all.hi[k] = krisite::kCoordMin;
        }
        for (const Poly& q : polys) {
            for (int k = 0; k < 3; ++k) {
                all.lo[k] = std::min(all.lo[k], q.aabb.lo[k]);
                all.hi[k] = std::max(all.hi[k], q.aabb.hi[k]);
            }
        }
        // 外接箱の外（奇数ずらしで整列を崩す）。**整数点**であることが主経路の条件
        std::int64_t c[3];
        bool ok = false;
        for (int axis = 0; axis < 3 && !ok; ++axis) {
            for (int sidek = 0; sidek < 2 && !ok; ++sidek) {
                for (int k = 0; k < 3; ++k) {
                    c[k] = (k == axis) ? ((sidek == 0) ? all.lo[k] - 1 : all.hi[k] + 1)
                                       : ((all.lo[k] + all.hi[k]) / 2 + 1 + 2 * k);
                }
                ok = true;
                for (int k = 0; k < 3; ++k) {
                    if (c[k] < krisite::kCoordMin || c[k] > krisite::kCoordMax) ok = false;
                }
            }
        }
        KRISITE_CHECK(ok, "boolean_nary: 参照点を外接箱の外に取れない");
        root_ref = geom::IPoint{static_cast<std::int32_t>(c[0]), static_cast<std::int32_t>(c[1]),
                                static_cast<std::int32_t>(c[2])};
    }

    // ---- 4. セルごとの arrangement（過剰分割。CP2 と同じ）--------------------
    const octree::SubdivisionPolicy policy{opt.depth, !opt.adaptive, opt.leaf_threshold};

    // ---- 分割の再帰 + 参照点の伝播（§5.3.1）------------------------------------
    //
    // `build_leaves` と同じ規則で降りながら、**整数の参照点とその WNV** を運びます。
    // 子の箱に親の参照点が入っていなければ、**箱の中の整数点へ軸平行にトレース**して
    // 更新します（§3.3.0 の主経路）。経路は親の箱の中に収まるので短く済みます。
    struct Sub {
        octree::Cell cell;
        geom::IPoint ref;
        std::vector<std::int32_t> w;
    };
    std::vector<Sub> stack{Sub{octree::Cell{0, 0, 0, 0}, root_ref,
                               std::vector<std::int32_t>(n_comp, 0)}};
    std::vector<Sub> leaves;

    auto in_box = [](const geom::IPoint& p, const octree::CellBox& b) {
        const std::int64_t c[3] = {p.x, p.y, p.z};
        for (int k = 0; k < 3; ++k) {
            if (c[k] < b.lo[k] || c[k] >= b.hi[k]) return false;
        }
        return true;
    };
    // その箱に重なる多角形だけを見る（経路は箱の中に収まる）
    auto polys_in = [&](const octree::CellBox& b) {
        std::vector<TracePoly> v;
        for (std::size_t i = 0; i < tps.size(); ++i) {
            if (octree::assign_to_cell(tps_aabb[i], b)) v.push_back(tps[i]);
        }
        return v;
    };

    while (!stack.empty()) {
        const Sub cur = std::move(stack.back());
        stack.pop_back();
        const octree::CellBox cb = octree::box_of(cur.cell);

        bool split = false;
        if (cur.cell.depth < policy.max_depth) {
            if (policy.uniform) {
                split = true;
            } else {
                std::size_t na = 0, nb = 0;
                for (const Poly& q : polys) {
                    if (!octree::assign_to_cell(q.aabb, cb)) continue;
                    if (q.comp == 0) {
                        ++na;
                    } else {
                        ++nb;
                    }
                }
                split = (na > 0 && nb > 0) && (na + nb > policy.leaf_threshold);
            }
        }
        if (!split) {
            leaves.push_back(cur);
            continue;
        }

        const std::vector<TracePoly> here_polys = polys_in(cb);
        for (int t2 = 0; t2 < 8; ++t2) {
            const octree::Cell child{cur.cell.depth + 1,
                                     cur.cell.i * 2 + static_cast<std::uint32_t>(t2 & 1),
                                     cur.cell.j * 2 + static_cast<std::uint32_t>((t2 >> 1) & 1),
                                     cur.cell.k * 2 + static_cast<std::uint32_t>((t2 >> 2) & 1)};
            const octree::CellBox chb = octree::box_of(child);
            if (in_box(cur.ref, chb)) {
                stack.push_back(Sub{child, cur.ref, cur.w});
                continue;
            }
            // 子の箱の中の整数点へ移す。**まず射影**（EMBER §4.2.2 の第一候補）
            bool done = false;
            for (int attempt = 0; attempt < 8 && !done; ++attempt) {
                std::int64_t c[3];
                for (int k = 0; k < 3; ++k) {
                    const std::int64_t lo = chb.lo[k], hi = chb.hi[k] - 1;
                    std::int64_t v = (attempt == 0) ? std::int64_t{(&cur.ref.x)[k]}
                                                    : (lo + hi) / 2 + attempt * (1 + 2 * k);
                    if (v < lo) v = lo;
                    if (v > hi) v = hi;
                    c[k] = v;
                }
                const geom::IPoint np{static_cast<std::int32_t>(c[0]),
                                      static_cast<std::int32_t>(c[1]),
                                      static_cast<std::int32_t>(c[2])};
                std::vector<std::int32_t> w = cur.w;
                std::size_t segs = 0;
                ++st.ref_moves;
                if (trace_axis_path(out.table, here_polys, cur.ref, np, &w, &segs)) {
                    st.ref_segments += segs;
                    stack.push_back(Sub{child, np, std::move(w)});
                    done = true;
                } else {
                    ++st.ref_rejects;
                }
            }
            KRISITE_CHECK(done, "boolean_nary: 子の参照点を作れない（経路がすべて退化）");
        }
    }
    std::sort(leaves.begin(), leaves.end(),
              [](const Sub& a, const Sub& b) { return a.cell < b.cell; });

    st.leaf_depth_min = opt.depth;
    for (const Sub& l : leaves) {
        st.leaf_depth_min = std::min(st.leaf_depth_min, l.cell.depth);
        st.leaf_depth_max = std::max(st.leaf_depth_max, l.cell.depth);
    }
    st.total_cells = leaves.size();

    std::vector<Fragment> frags;
    std::vector<octree::Cell> frag_cell;
    std::vector<std::uint32_t> frag_comp, frag_tag;

    std::vector<std::vector<std::int32_t>> frag_wref;
    std::vector<geom::IPoint> frag_ref;
    for (const Sub& leaf : leaves) {
        const octree::Cell& cell = leaf.cell;
        const octree::CellBox cbox = octree::box_of(cell);
        const std::size_t frags_before = frags.size();
        std::vector<Fragment> local;
        std::vector<std::uint32_t> local_comp, local_tag;

        std::vector<std::size_t> here;
        for (std::size_t i = 0; i < polys.size(); ++i) {
            if (octree::assign_to_cell(polys[i].aabb, cbox)) here.push_back(i);
        }
        if (here.empty()) {
            ++st.empty_cells;
            continue;
        }

        struct CellPlane {
            PlaneId id;
            int keep;
        };
        std::vector<CellPlane> cps;
        if (cell.depth > 0) {
            const auto ps = octree::cell_planes(cell);
            for (int k = 0; k < 6; ++k) {
                const PlaneRef r = out.table.intern(ps[k]);
                const int base = (k % 2 == 0) ? +1 : -1;
                cps.push_back({r.id, r.flipped ? -base : base});
            }
        }

        std::vector<PlaneId> culled;
        if (opt.cull_planes) {
            for (PlaneId q : all_split) {
                if (geom::plane_crosses_box(out.table.at(q), cbox.lo, cbox.hi)) culled.push_back(q);
            }
        }
        const std::vector<PlaneId>& split_planes = opt.cull_planes ? culled : all_split;
        st.split_plane_slots += all_split.size();
        st.split_planes_used += split_planes.size();
        st.max_planes_per_cell = std::max(st.max_planes_per_cell, split_planes.size());

        for (std::size_t idx : here) {
            Fragment frag = polys[idx].frag;
            bool alive = true;
            for (const CellPlane& cp : cps) {
                if (cp.id == frag.support) continue;
                if (!clip_fragment(out.table, frag, cp.id, cp.keep, cache)) {
                    alive = false;
                    break;
                }
            }
            if (!alive) continue;
            std::vector<Fragment> pieces{frag};
            for (PlaneId q : split_planes) {
                std::vector<Fragment> next;
                next.reserve(pieces.size());
                for (const Fragment& p : pieces) {
                    if (q == p.support) {
                        next.push_back(p);
                        continue;
                    }
                    const SplitResult r = split_fragment(out.table, p, q, cache);
                    if (r.has_pos) next.push_back(r.pos);
                    if (r.has_neg) next.push_back(r.neg);
                }
                pieces.swap(next);
            }
            for (Fragment& p : pieces) {
                local.push_back(std::move(p));
                local_comp.push_back(polys[idx].comp);
                local_tag.push_back(polys[idx].tag);
            }
        }

        // 共平面重複を互いの辺平面で揃える（EMBER §4.3 の C4。CP2 と同じ）
        {
            std::map<PlaneId, std::vector<std::size_t>> by_sup;
            for (std::size_t i = 0; i < local.size(); ++i) by_sup[local[i].support].push_back(i);
            std::vector<Fragment> nl;
            std::vector<std::uint32_t> nc, nt;
            for (const auto& g : by_sup) {
                bool multi = false;
                for (std::size_t i : g.second) {
                    if (local_comp[i] != local_comp[g.second.front()] ||
                        local_tag[i] != local_tag[g.second.front()]) {
                        multi = true;
                        break;
                    }
                }
                if (g.second.size() < 2 || !multi) {
                    for (std::size_t i : g.second) {
                        nl.push_back(local[i]);
                        nc.push_back(local_comp[i]);
                        nt.push_back(local_tag[i]);
                    }
                    continue;
                }
                std::vector<PlaneId> es;
                for (std::size_t i : g.second) {
                    es.insert(es.end(), local[i].edge.begin(), local[i].edge.end());
                }
                std::sort(es.begin(), es.end());
                es.erase(std::unique(es.begin(), es.end()), es.end());
                for (std::size_t i : g.second) {
                    std::vector<Fragment> pieces{local[i]};
                    for (PlaneId q : es) {
                        std::vector<Fragment> next;
                        for (const Fragment& p : pieces) {
                            if (q == p.support) {
                                next.push_back(p);
                                continue;
                            }
                            const SplitResult r = split_fragment(out.table, p, q, cache);
                            if (r.has_pos) next.push_back(r.pos);
                            if (r.has_neg) next.push_back(r.neg);
                        }
                        pieces.swap(next);
                    }
                    for (Fragment& p : pieces) {
                        nl.push_back(std::move(p));
                        nc.push_back(local_comp[i]);
                        nt.push_back(local_tag[i]);
                    }
                }
            }
            local.swap(nl);
            local_comp.swap(nc);
            local_tag.swap(nt);
        }

        for (std::size_t i = 0; i < local.size(); ++i) {
            frags.push_back(std::move(local[i]));
            frag_cell.push_back(cell);
            frag_comp.push_back(local_comp[i]);
            frag_tag.push_back(local_tag[i]);
            frag_wref.push_back(leaf.w);
            frag_ref.push_back(leaf.ref);
        }
        if (frags.size() != frags_before) ++st.active_cells;
    }
    st.raw_fragments = frags.size();

    // ---- 5. 縫合（重複の仕分けに要る。CP2 と同じ）----------------------------
    std::map<std::array<PlaneId, 3>, std::uint32_t> by_key;
    std::vector<geom::HPointD> points;
    std::vector<std::vector<std::uint32_t>> raw(frags.size());
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        const Fragment& f = frags[fi];
        const std::size_t nv = vertex_count(f);
        raw[fi].reserve(nv);
        for (std::size_t i = 0; i < nv; ++i) {
            const auto k = detail::vertex_key(f, i);
            auto it = by_key.find(k);
            if (it == by_key.end()) {
                const auto id = static_cast<std::uint32_t>(points.size());
                points.push_back(fragment_vertex(out.table, f, i, cache));
                it = by_key.emplace(k, id).first;
            }
            raw[fi].push_back(it->second);
        }
    }
    st.constructed_points = points.size();
    std::vector<std::uint32_t> order(points.size());
    for (std::uint32_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
        return geom::lex_less(points[a], points[b]);
    });
    std::vector<std::uint32_t> remap(points.size());
    std::size_t merged_count = 0;
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        const auto id = static_cast<std::uint32_t>(merged_count++);
        while (j < order.size() && geom::h_equal(points[order[i]], points[order[j]])) {
            remap[order[j]] = id;
            ++j;
        }
        if (j - i > 1) st.merged_by_value += (j - i - 1);
        i = j;
    }
    st.merged_points = merged_count;

    std::map<detail::RegionKey, std::vector<std::size_t>> regions;
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        if (raw[fi].size() < 3) continue;
        std::vector<std::uint32_t> ids;
        ids.reserve(raw[fi].size());
        for (std::uint32_t v : raw[fi]) ids.push_back(remap[v]);
        regions[detail::region_key(frags[fi].support, std::move(ids))].push_back(fi);
    }

    // ---- 6. 分類（セグメントトレース。§5.5）----------------------------------
    //
    // **参照点から断片の内部点まで経路を張り、跨いだ多角形の Δw を足します。**
    // レイキャストと違い、無限遠まで数える必要がありません。
    for (const auto& kv : regions) {
        const std::size_t pick = *std::min_element(
            kv.second.begin(), kv.second.end(), [&](std::size_t a, std::size_t b) {
                return std::make_tuple(frag_comp[a], frag_tag[a], a) <
                       std::make_tuple(frag_comp[b], frag_tag[b], b);
            });
        const Fragment& f = frags[pick];
        // **経路は葉の箱に収まる**ので、その箱に重なる多角形だけを見ます（局所化）
        const std::vector<TracePoly> leaf_polys = polys_in(octree::box_of(frag_cell[pick]));

        // **この領域に載っている面（シート）の寄与。**
        //
        // トレースは始点の面を跨がないので、得られるのは**表裏に共通の分**です。
        // シートの内部がどちら側にあるかで、表と裏に別々に足します。
        //
        //   外向き法線が $+N$（orient = +1）→ 内部は $-N$ 側 → **裏**に +1
        //   外向き法線が $-N$（orient = -1）→ **表**に +1
        std::vector<std::int32_t> d_front(n_comp, 0), d_back(n_comp, 0);
        {
            // **同じ入力多角形の断片を二重に数えないこと。** 同一領域が複数のセルに
            // 割り当てられると、同じ多角形の断片が複数回現れます（重複割り当て）。
            // 真偽値なら代入なので無害でしたが、**巻き数では足し算なので効きます。**
            std::vector<std::pair<std::uint32_t, std::uint32_t>> seen;
            for (std::size_t fi : kv.second) {
                const auto key = std::make_pair(frag_comp[fi], frag_tag[fi]);
                if (std::find(seen.begin(), seen.end(), key) != seen.end()) continue;
                seen.push_back(key);
                if (frags[fi].flipped) {
                    ++d_front[frag_comp[fi]];
                } else {
                    ++d_back[frag_comp[fi]];
                }
            }
        }

        // **経路が退化したら、別の内部点で試します**（§3.3 の順序 6 通り x 内部点）。
        //
        // 軸平行な入力では、内部点を通る軸平行線が他の面の辺をちょうど通ることが
        // あります。順序を変えても直らないので、**点そのものを変えます**。
        // ---- 分類（§5.5）--------------------------------------------------
        //
        // **主経路**（§3.3.0）: 内部点 x は整数アンカー c を通る軸平行線の上にあるので、
        //   x → c は軸平行 1 本。c → 葉の参照点も整数どうしの軸平行 3 本まで。
        //   **経路は葉の箱に収まります**（参照点の伝播があるから）。
        //
        // **予備**（§3.3.1）: 平面を 1 枚ずつ置き換える 3 セグメント経路。
        const std::vector<std::int32_t>& w_leaf = frag_wref[pick];
        const geom::IPoint& ref_leaf = frag_ref[pick];
        std::vector<std::int32_t> w_at_x;
        bool traced = false;

        const std::size_t nv_f = vertex_count(f);
        for (unsigned variant = 0; variant <= nv_f && !traced; ++variant) {
            std::array<PlaneId, 3> xp{};
            geom::IPoint anchor{};
            geom::Axis aaxis = geom::Axis::X;
            const geom::HPointD x =
                interior_point(out.table, f, cache, &st.interior, &xp, variant, &anchor, &aaxis);
            std::vector<std::int32_t> w = w_leaf;

            if (variant == 0) {
                // x → アンカー（軸平行 1 本。x は面の上にある）
                const geom::Axis u = (aaxis == geom::Axis::X)   ? geom::Axis::Y
                                     : (aaxis == geom::Axis::Y) ? geom::Axis::Z
                                                                : geom::Axis::X;
                const geom::Axis v = (aaxis == geom::Axis::X)   ? geom::Axis::Z
                                     : (aaxis == geom::Axis::Y) ? geom::Axis::X
                                                                : geom::Axis::Y;
                const std::int64_t ac[3] = {anchor.x, anchor.y, anchor.z};
                const PlaneId l0 =
                    out.table.intern(geom::plane_axis_aligned(u, ac[static_cast<int>(u)])).id;
                const PlaneId l1 =
                    out.table.intern(geom::plane_axis_aligned(v, ac[static_cast<int>(v)])).id;
                const PlaneId eb =
                    out.table.intern(
                            geom::plane_axis_aligned(aaxis, ac[static_cast<int>(aaxis)]))
                        .id;
                std::size_t segs = 0;
                // **アンカーが支持平面の上にあることがあります**（軸平行な面では
                // 重心の丸めがちょうど面に乗る）。その場合 x == アンカーで、
                // 最初のセグメントは長さ 0。**面の上から出発する扱いを引き継ぎます。**
                const bool anchor_on_plane =
                    geom::side(out.table.at(f.support), geom::to_homogeneous(anchor)) == 0;
                const bool seg0 =
                    anchor_on_plane ||
                    trace_segment(out.table, leaf_polys, l0, l1, f.support, eb, x,
                                  geom::to_homogeneous(anchor), &w, true, false);
                if (seg0 && trace_axis_path(out.table, leaf_polys, anchor, ref_leaf, &w, &segs,
                                            anchor_on_plane)) {
                    traced = true;
                    if (segs == 0) ++st.path_single;
                    w_at_x = std::move(w);
                    break;
                }
                ++st.midpoint_retries;
                continue;
            }
            // 予備経路
            const std::array<PlaneId, 3> rp{
                out.table.intern(geom::plane_with_normal(1, 1, 0, ref_leaf)).id,
                out.table.intern(geom::plane_with_normal(0, 1, 1, ref_leaf)).id,
                out.table.intern(geom::plane_with_normal(1, 0, 1, ref_leaf)).id};
            const TraceResult tr = trace_path(out.table, leaf_polys, xp, rp, w_leaf, cache);
            if (tr.ok) {
                traced = true;
                ++st.path_fallback;
                w_at_x = tr.dw;
                break;
            }
            ++st.midpoint_retries;
        }
#if defined(KRISITE_TRACE_COUNT_ONLY)
        if (!traced) {
            ++st.path_fallback;  // 集計用: 失敗した領域
            continue;
        }
#else
        KRISITE_CHECK(traced, "boolean_nary: 経路がすべて退化した（内部点も順序も尽きた）");
#endif
        ++st.raycasts;

        std::vector<std::int32_t> w_front = w_at_x, w_back = w_at_x;
        for (std::size_t i = 0; i < n_comp; ++i) {
            w_front[i] += d_front[i];
            w_back[i] += d_back[i];
        }

        const bool in_front = ind.eval(w_front);
        const bool in_back = ind.eval(w_back);
        if (in_front == in_back) continue;

        Poly q;
        q.frag = f;
        q.frag.flipped = false;
        q.comp = 0;
        q.tag = frag_tag[pick];
        const octree::CellBox cb2 = octree::box_of(frag_cell[pick]);
        for (int k = 0; k < 3; ++k) {
            q.aabb.lo[k] = cb2.lo[k];
            q.aabb.hi[k] = cb2.hi[k];
        }
        if (in_front) {
            std::vector<std::uint32_t> dummy(vertex_count(q.frag), 0);
            std::vector<PlaneId> e = q.frag.edge;
            detail::reverse_polygon(dummy, e);
            q.frag.edge = std::move(e);
            q.frag.flipped = true;
        }
        out.polys.push_back(std::move(q));
    }
    st.fragments = out.polys.size();
    st.cache_hits = point_cache.hits();
    st.cache_misses = point_cache.misses();
    st.cache_entries = point_cache.entries();
    st.cache_bytes = point_cache.bytes();
    if (stats != nullptr) *stats = st;
    return out;
}

/// 二項の呼び出し形（$n = 2$ の特殊ケース）。
inline PolySoup boolean(const PolySoup& X, const PolySoup& Y, BoolOp op, const BoolOptions& opt,
                        BoolStats* stats = nullptr) {
    const Compose how = (op == BoolOp::Union)          ? Compose::Union
                        : (op == BoolOp::Intersection) ? Compose::Intersection
                                                       : Compose::Difference;
    const Indicator ind = compose(indicator_source(0), indicator_source(0), 1, how);
    return boolean_nary({&X, &Y}, ind, opt, stats);
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_SOUP_BOOLEAN_HPP

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
    tps.reserve(polys.size());
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
    }

    // ---- 3. 参照点（**全多角形の外**。だから WNV は 0）--------------------------
    //
    // **点は平面 3 つ組で持ちます**（経路の構成に要る。§3.3）。
    //
    // 外接箱の外に取れば「すべての立体の外」が保証され、巻き数は 0 です。
    // **座標は奇数ずらしにします。** 入力が軸平行だと、整列した参照点への経路が
    // 他の面の辺をちょうど通り、経路が退化します（実際に踏みました）。
    std::vector<PPoint> ref_candidates;
    const std::vector<std::int32_t> w_ref(n_comp, 0);
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
        const std::int64_t lim = -static_cast<std::int64_t>(krisite::kCoordMin);
        const std::int64_t odd[3] = {1, 3, 5};
        for (int axis = 0; axis < 3; ++axis) {
            for (int sidek = 0; sidek < 2; ++sidek) {
                std::int64_t c[3];
                bool ok = true;
                for (int k = 0; k < 3; ++k) {
                    if (k == axis) {
                        c[k] = (sidek == 0) ? all.lo[k] - 1 : all.hi[k] + 1;
                        if (c[k] < krisite::kCoordMin || c[k] > lim) ok = false;
                    } else {
                        // 箱の中央から奇数ぶんずらす（整列を崩す）
                        c[k] = (all.lo[k] + all.hi[k]) / 2 + odd[k];
                        if (c[k] < krisite::kCoordMin) c[k] = krisite::kCoordMin;
                        if (c[k] > lim) c[k] = lim;
                    }
                }
                if (!ok) continue;
                // **斜めの平面で定義します。** 軸平行だと、軸平行な入力との組で
                // 平面が平行になり、経路の点が作れません（実際に踏みました）。
                // 法線 (1,1,0),(0,1,1),(1,0,1) は行列式 2 で独立です。
                const geom::IPoint rp{static_cast<std::int32_t>(c[0]),
                                      static_cast<std::int32_t>(c[1]),
                                      static_cast<std::int32_t>(c[2])};
                ref_candidates.push_back(
                    {out.table.intern(geom::plane_with_normal(1, 1, 0, rp)).id,
                     out.table.intern(geom::plane_with_normal(0, 1, 1, rp)).id,
                     out.table.intern(geom::plane_with_normal(1, 0, 1, rp)).id});
            }
        }
        KRISITE_CHECK(!ref_candidates.empty(),
                      "boolean_nary: 参照点を外接箱の外に取れない（領域全体を覆う入力）");
    }

    // ---- 4. セルごとの arrangement（過剰分割。CP2 と同じ）--------------------
    const octree::SubdivisionPolicy policy{opt.depth, !opt.adaptive, opt.leaf_threshold};
    const std::vector<octree::Cell> leaves =
        octree::build_leaves(policy, [&](const octree::Cell& c, std::size_t* na, std::size_t* nb) {
            const octree::CellBox cb = octree::box_of(c);
            *na = 0;
            *nb = 0;
            for (const Poly& q : polys) {
                if (!octree::assign_to_cell(q.aabb, cb)) continue;
                if (q.comp == 0) {
                    ++*na;
                } else {
                    ++*nb;
                }
            }
        });
    st.leaf_depth_min = opt.depth;
    for (const octree::Cell& c : leaves) {
        st.leaf_depth_min = std::min(st.leaf_depth_min, c.depth);
        st.leaf_depth_max = std::max(st.leaf_depth_max, c.depth);
    }
    st.total_cells = leaves.size();

    std::vector<Fragment> frags;
    std::vector<octree::Cell> frag_cell;
    std::vector<std::uint32_t> frag_comp, frag_tag;

    for (const octree::Cell& cell : leaves) {
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
        TraceResult tr;
        const std::size_t nv_f = vertex_count(f);
        for (unsigned variant = 0; variant <= nv_f && !tr.ok; ++variant) {
            PPoint xp{};
            interior_point(out.table, f, cache, &st.interior, &xp, variant);
            for (const PPoint& xr : ref_candidates) {
                tr = trace_path(out.table, tps, xp, xr, w_ref, cache);
                if (tr.ok) break;
                ++st.midpoint_retries;
            }
        }
        KRISITE_CHECK(tr.ok, "boolean_nary: 経路がすべて退化した（内部点も順序も尽きた）");
        ++st.raycasts;

        std::vector<std::int32_t> w_front = tr.dw, w_back = tr.dw;
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

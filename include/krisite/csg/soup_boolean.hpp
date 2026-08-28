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
#include "krisite/csg/raycast.hpp"
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

/// 2 つのスープのブール演算（CP2）。**出力もスープなので連鎖できます。**
inline PolySoup boolean(const PolySoup& X, const PolySoup& Y, BoolOp op, const BoolOptions& opt,
                        BoolStats* stats = nullptr) {
    BoolStats st;
    PointCache point_cache;
    PointCache* const cache = opt.cache_points ? &point_cache : nullptr;

    PolySoup out;
    // ---- 1. sources と指示関数の合成（§5.2）----------------------------------
    out.sources = X.sources;
    out.sources.insert(out.sources.end(), Y.sources.begin(), Y.sources.end());
    const auto off = static_cast<std::uint32_t>(X.sources.size());
    const Compose how = (op == BoolOp::Union)          ? Compose::Union
                        : (op == BoolOp::Intersection) ? Compose::Intersection
                                                       : Compose::Difference;
    out.indicator = compose(X.indicator, Y.indicator, off, how);

    // ---- 2. 平面表を 1 つにまとめ、多角形を移す ------------------------------
    std::vector<Poly> polys;
    polys.reserve(X.polys.size() + Y.polys.size());
    for (int which = 0; which < 2; ++which) {
        const PolySoup& s = (which == 0) ? X : Y;
        const std::vector<PlaneRef> m = detail::remap_planes(s.table, out.table);
        for (const Poly& q : s.polys) {
            Poly r = q;
            r.frag.support = m[q.frag.support].id;
            r.frag.flipped = (q.frag.flipped != m[q.frag.support].flipped);
            for (PlaneId& e : r.frag.edge) e = m[e].id;
            r.src = q.src + (which == 0 ? 0u : off);
            r.frag.owner = which;
            polys.push_back(std::move(r));
        }
    }
    const std::size_t n_src = out.sources.size();

    // source ごとの支持平面（分類のキーに使う。§4.3.2 の符号ベクトルの一般化）。
    //
    // **残った多角形からではなく、source のメッシュ全体から作ります。** 連鎖では
    // 前段で消えた面があるので、多角形から集めると平面が抜け、符号ベクトルが
    // 粗くなって別々の領域が同じキーに落ちます（**実際に踏みました**）。
    //
    // 分類が断片上で一定であるためには、断片が**各 source の平面配置の 1 セル**に
    // 収まっている必要があります。だから分割にも全 source の平面が要ります。
    std::vector<std::vector<PlaneId>> planes_of_src(n_src);
    for (std::size_t i = 0; i < n_src; ++i) {
        const mesh::TriMesh& m = out.sources[i];
        for (const mesh::Tri& t : m.triangles) {
            const geom::PlaneD pl =
                geom::plane_from_triangle(m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]]);
            if (geom::is_degenerate(pl)) continue;
            planes_of_src[i].push_back(out.table.intern(pl).id);
        }
    }
    for (std::vector<PlaneId>& v : planes_of_src) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    std::vector<PlaneId> all_split;
    for (const std::vector<PlaneId>& v : planes_of_src) {
        all_split.insert(all_split.end(), v.begin(), v.end());
    }
    std::sort(all_split.begin(), all_split.end());
    all_split.erase(std::unique(all_split.begin(), all_split.end()), all_split.end());

    // source ごとの三角形の AABB。**early-out の判定に使います。**
    //
    // **多角形の有無で判定してはいけません。** 連鎖では前段で消えた面があるので、
    // 「この source の多角形が無い = この source の曲面が無い」は成り立ちません。
    // 曲面が横切っているのに「内外が一定」と決めつけると分類が壊れます
    // （**実際に踏みました**）。
    std::vector<std::vector<octree::Aabb>> src_aabb(n_src);
    for (std::size_t i = 0; i < n_src; ++i) {
        const mesh::TriMesh& m = out.sources[i];
        src_aabb[i].reserve(m.triangles.size());
        for (const mesh::Tri& t : m.triangles) {
            octree::Aabb r{};
            for (int k = 0; k < 3; ++k) {
                r.lo[k] = krisite::kCoordMax;
                r.hi[k] = krisite::kCoordMin;
            }
            for (int v = 0; v < 3; ++v) {
                const geom::IPoint& p = m.vertices[t[v]];
                const std::int64_t cc[3] = {p.x, p.y, p.z};
                for (int k = 0; k < 3; ++k) {
                    r.lo[k] = std::min(r.lo[k], cc[k]);
                    r.hi[k] = std::max(r.hi[k], cc[k]);
                }
            }
            src_aabb[i].push_back(r);
        }
    }

    // ---- 3. 葉の列挙（§3.1。固定深度は「常に最大深度」の特別な場合）----------
    const octree::SubdivisionPolicy policy{opt.depth, !opt.adaptive, opt.leaf_threshold};
    const std::vector<octree::Cell> leaves =
        octree::build_leaves(policy, [&](const octree::Cell& c, std::size_t* na, std::size_t* nb) {
            const octree::CellBox cb = octree::box_of(c);
            *na = 0;
            *nb = 0;
            for (const Poly& q : polys) {
                if (!octree::assign_to_cell(q.aabb, cb)) continue;
                if (q.frag.owner == 0) {
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

    // ---- 4. セルごとの arrangement -------------------------------------------
    std::vector<Fragment> frags;
    std::vector<octree::Cell> frag_cell;
    std::vector<std::uint32_t> frag_src, frag_tag;
    /// セルで「多角形が 1 枚も無かった source」の内外（-1 = 未確定）。§3.2 の early-out
    std::vector<std::vector<std::int8_t>> frag_forced;

    for (const octree::Cell& cell : leaves) {
        const octree::CellBox cbox = octree::box_of(cell);
        const std::size_t frags_before = frags.size();
        std::vector<Fragment> local;
        std::vector<std::uint32_t> local_src, local_tag;

        std::vector<std::size_t> here;
        for (std::size_t i = 0; i < polys.size(); ++i) {
            if (octree::assign_to_cell(polys[i].aabb, cbox)) here.push_back(i);
        }
        if (here.empty()) {
            ++st.empty_cells;
            continue;
        }

        // **曲面が 1 枚も横切らない source は、セル全体で内外が一定です。**
        // 隅 1 点（整数点）で決まります（§3.2 の early-out の一般化）。
        std::vector<std::size_t> count(n_src, 0);
        for (std::size_t i = 0; i < n_src; ++i) {
            for (const octree::Aabb& r : src_aabb[i]) {
                if (octree::assign_to_cell(r, cbox)) ++count[i];
            }
        }
        std::vector<std::int8_t> forced(n_src, -1);
        if (opt.early_out) {
            const geom::IPoint corner{static_cast<std::int32_t>(cbox.lo[0]),
                                      static_cast<std::int32_t>(cbox.lo[1]),
                                      static_cast<std::int32_t>(cbox.lo[2])};
            bool any = false;
            for (std::size_t i = 0; i < n_src; ++i) {
                if (count[i] > 0) continue;
                int wo = 0, cf = 0, cb = 0;
                winding_split(out.sources[i], geom::to_homogeneous(corner),
                              out.table.at(all_split.empty() ? PlaneId{0} : all_split.front()), &wo,
                              &cf, &cb);
                KRISITE_CHECK(cf == 0 && cb == 0,
                              "early-out: セルの隅が source の曲面に載っている"
                              "（曲面が無いはずのセル）");
                forced[i] = static_cast<std::int8_t>(wo);
                ++st.early_out_raycasts;
                any = true;
            }
            if (any) ++st.early_out_cells;
        }

        // セル面の平面（保持側つき）: lo 面は +、hi 面は -
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

        // 分割平面の絞り込み（`SPEC-phase2.md` §2.3）
        std::vector<PlaneId> culled;
        if (opt.cull_planes) {
            const std::int64_t* clo = cbox.lo;
            const std::int64_t* chi = cbox.hi;
            for (PlaneId q : all_split) {
                if (geom::plane_crosses_box(out.table.at(q), clo, chi)) culled.push_back(q);
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

            // 【全支持平面】で分割（過剰分割。§4.3.1。CP4 で局所 BSP に置き換わる）
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
                local_src.push_back(polys[idx].src);
                local_tag.push_back(polys[idx].tag);
            }
        }

        // ---- 共平面重複を【互いの辺平面で切って揃える】（EMBER §4.3 の C4）------
        //
        // **面併合をやめた代償です**（§3.1.4）。同じ平面に載る面が別々の三角形分割を
        // 持つと、重なりの領域が一致せず、`region_key` で潰せません。
        //
        // 支持平面が同じ断片どうしを、**相手の辺平面**で切ると領域が揃います。
        // 自分の辺平面で切っても no-op なので、まとめて適用して構いません。
        {
            std::map<PlaneId, std::vector<std::size_t>> by_sup;
            for (std::size_t i = 0; i < local.size(); ++i) by_sup[local[i].support].push_back(i);

            // **新しい配列を組み立てます。** 元の配列を消しながら回すと、
            // 2 つ目のグループ以降で添字が無効になります（実際に踏みました）。
            std::vector<Fragment> nl;
            std::vector<std::uint32_t> ns, nt;
            for (const auto& g : by_sup) {
                // 由来が 1 つだけなら、同じ多角形の断片どうしなので揃っています
                bool multi = false;
                for (std::size_t i : g.second) {
                    if (local_src[i] != local_src[g.second.front()] ||
                        local_tag[i] != local_tag[g.second.front()]) {
                        multi = true;
                        break;
                    }
                }
                if (g.second.size() < 2 || !multi) {
                    for (std::size_t i : g.second) {
                        nl.push_back(local[i]);
                        ns.push_back(local_src[i]);
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
                        ns.push_back(local_src[i]);
                        nt.push_back(local_tag[i]);
                    }
                }
            }
            local.swap(nl);
            local_src.swap(ns);
            local_tag.swap(nt);
        }

        for (std::size_t i = 0; i < local.size(); ++i) {
            frags.push_back(std::move(local[i]));
            frag_cell.push_back(cell);
            frag_src.push_back(local_src[i]);
            frag_tag.push_back(local_tag[i]);
            frag_forced.push_back(forced);
        }
        if (frags.size() != frags_before) ++st.active_cells;
    }
    st.raw_fragments = frags.size();

    // ---- 5. 縫合（重複の仕分けに要る）----------------------------------------
    std::map<std::array<PlaneId, 3>, std::uint32_t> by_key;
    std::vector<geom::HPointD> points;
    std::vector<std::vector<std::uint32_t>> raw(frags.size());
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        const Fragment& f = frags[fi];
        const std::size_t n = vertex_count(f);
        raw[fi].reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
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

    // ---- 6. 重複の仕分け（§4.3.3 / §5.5）-------------------------------------
    //
    // 同じ領域の断片は、(支持平面, 頂点集合) が一致します。**共平面重複は
    // 全順序で 1 枚だけ残します**（`SPEC-phase3.md` §5.4.1 / EMBER §4.3 の C4）。
    std::map<detail::RegionKey, std::vector<std::size_t>> regions;
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        if (raw[fi].size() < 3) continue;
        std::vector<std::uint32_t> ids;
        ids.reserve(raw[fi].size());
        for (std::uint32_t v : raw[fi]) ids.push_back(remap[v]);
        regions[detail::region_key(frags[fi].support, std::move(ids))].push_back(fi);
    }

    // ---- 7. 分類（WNV。`SPEC-phase3.md` §5.1 / §14 の CP3 の変更 1）--------------
    //
    // **内外の 1 ビットではなく巻き数で持ちます。** 同じ曲面を 2 度跨ぐ配置
    // （同じメッシュを 2 度使う連鎖、入れ子、自己交差）は、真偽値では表せません。
    //
    // 断片の表裏の巻き数は `winding_split` が 3 つに分けて返します。
    //
    //   w_other  断片に載っていない面から決まる分（表裏で共通）
    //   c_front  法線の側（+N）で、載っている面が寄与する分
    //   c_back   反対側（-N）で寄与する分
    //
    // **載っている面を除いてからレイキャストする**のが要点です。境界上から飛ばすと
    // 契約（`point_on_boundary` でないこと）が破れます。
    struct Wnd {
        int w_other, c_front, c_back;
    };
    std::map<std::pair<std::uint32_t, std::vector<std::int8_t>>, Wnd> wcache;

    for (const auto& kv : regions) {
        // 同じ領域に複数の断片が載っていても、出力するのは 1 枚です（§5.4.1）。
        // **巻き数は source のメッシュから決まる**ので、どの断片を代表に取っても同じです。
        const std::size_t pick = *std::min_element(
            kv.second.begin(), kv.second.end(), [&](std::size_t a, std::size_t b) {
                return std::make_tuple(frag_src[a], frag_tag[a], a) <
                       std::make_tuple(frag_src[b], frag_tag[b], b);
            });
        const Fragment& f = frags[pick];
        const geom::PlaneD& refpl = out.table.at(f.support);

        bool need_point = false;
        for (std::size_t i2 = 0; i2 < n_src; ++i2) {
            if (frag_forced[pick][i2] < 0) need_point = true;
        }
        geom::HPointD rep{};
        if (need_point) rep = interior_point(out.table, f, cache, &st.interior);

        std::vector<std::int32_t> w_front(n_src, 0), w_back(n_src, 0);
        for (std::size_t i2 = 0; i2 < n_src; ++i2) {
            if (frag_forced[pick][i2] >= 0) {
                // セルに曲面が無い source。隅で決めた巻き数が全体で一定
                w_front[i2] = frag_forced[pick][i2];
                w_back[i2] = frag_forced[pick][i2];
                continue;
            }
            std::vector<std::int8_t> sig(planes_of_src[i2].size());
            for (std::size_t k = 0; k < planes_of_src[i2].size(); ++k) {
                sig[k] = static_cast<std::int8_t>(
                    fragment_sign(out.table, f, planes_of_src[i2][k], cache));
            }
            const auto key = std::make_pair(static_cast<std::uint32_t>(i2), std::move(sig));
            auto it = wcache.find(key);
            if (it == wcache.end()) {
                Wnd v{};
                winding_split(out.sources[i2], rep, refpl, &v.w_other, &v.c_front, &v.c_back);
                ++st.raycasts;
                ++st.regions;
                it = wcache.emplace(key, v).first;
            }
            w_front[i2] = it->second.w_other + it->second.c_front;
            w_back[i2] = it->second.w_other + it->second.c_back;
        }

        const bool in_front = out.indicator.eval(w_front);
        const bool in_back = out.indicator.eval(w_back);
#if defined(KRISITE_DEBUG_SOUP)
        {
            std::fprintf(stderr, "region src=%u tag=%u n=%zu wF=[", frag_src[pick], frag_tag[pick],
                         kv.second.size());
            for (std::int32_t v : w_front) std::fprintf(stderr, "%d ", v);
            std::fprintf(stderr, "] wB=[");
            for (std::int32_t v : w_back) std::fprintf(stderr, "%d ", v);
            std::fprintf(stderr, "] → %d/%d %s\n", (int)in_front, (int)in_back,
                         in_front == in_back ? "捨てる" : "出力");
        }
#endif
        if (in_front == in_back) continue;  // (in,in) / (out,out) は捨てる（§5.2）

        Poly q;
        q.frag = f;
        q.frag.flipped = false;  // 基準は支持平面の法線
        q.src = frag_src[pick];
        q.tag = frag_tag[pick];
        const octree::CellBox cb2 = octree::box_of(frag_cell[pick]);
        for (int k = 0; k < 3; ++k) {
            q.aabb.lo[k] = cb2.lo[k];
            q.aabb.hi[k] = cb2.hi[k];
        }
        // **立体は in の側にあります。** 外向き法線は立体から外を向くので、
        // +N 側が in なら法線は -N、つまり反転して出力します（§5.2）。
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

}  // namespace krisite::csg

#endif  // KRISITE_CSG_SOUP_BOOLEAN_HPP

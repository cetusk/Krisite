// Krisite — 中核（`SPEC-phase3.md` §5 / §14 の CP2）
//
//     boolean : PolySoup × PolySoup → PolySoup     ★ CSG について閉じる
//
// **CP2 で入口と出口を分離し、CP3 で分類を WNV（$\mathbb{Z}^n$）へ移しました。**
// **CP4 で過剰分割を局所 BSP に置き換えます**（§5.4）。
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
#include <chrono>
#include <map>
#include <mutex>
#include <vector>
#if defined(KRISITE_DEBUG_SOUP)
#include <cstdio>
#endif

#include "krisite/csg/boolean.hpp"  // BoolOp / BoolOptions / BoolStats
#include "krisite/csg/interior.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/raycast.hpp"
#include "krisite/octree/adaptive.hpp"
#include "krisite/par/thread_pool.hpp"

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

/// 三角形が平面に触れる（またはまたぐ）か。**厳密に片側なら false。**
///
/// 局所 BSP の切断条件です（`SPEC-phase3.md` §5.4）。平面 `pl` の上に載る断片 $t$ と
/// 三角形 $v$ が交わるのは、$v$ が `pl` に触れるときだけです。厳密に片側に居る $v$ は
/// $t$ と交わらないので、$t$ を $v$ の平面で切る必要がありません。
///
/// **新しい述語ではありません。** `geom::side` を整数点に対して呼ぶだけで、
/// 幅は $3b+6$（`widths.hpp` の `kSide` に収まります）。
inline bool tri_touches_plane(const geom::PlaneD& pl, const mesh::TriMesh& m, const mesh::Tri& t) {
    int acc = 0;
    for (int v = 0; v < 3; ++v) {
        const int sg = geom::side(pl, geom::to_homogeneous(m.vertices[t[v]]));
        if (sg == 0) return true;  // 平面上の頂点がある = 触れている
        if (acc == 0) {
            acc = sg;
        } else if (acc != sg) {
            return true;  // またいでいる
        }
    }
    return false;
}

#if defined(KRISITE_EXPERIMENT_BSP_SKIP_DISJOINT)
/// 断片が箱と**分離している**か（どれかの軸で厳密に外側にあるか）。
///
/// **最適化の候補**であって変異ではありません（`ROADMAP.md`「切断候補の絞り込み」）。
/// **軸平行平面に対する `side` だけ**で判定できます。
inline bool fragment_outside_box(const PlaneTable& t, const Fragment& f, const octree::Aabb& box,
                                 PointCache* cache) {
    const std::size_t n = vertex_count(f);
    for (int k = 0; k < 3; ++k) {
        const auto ax = static_cast<geom::Axis>(k);
        // 箱の lo 面より手前（side < 0）に全頂点があるか
        for (int side_sel = 0; side_sel < 2; ++side_sel) {
            const geom::PlaneD pl =
                geom::plane_axis_aligned(ax, side_sel == 0 ? box.lo[k] : box.hi[k]);
            const int want = (side_sel == 0) ? -1 : +1;
            bool all = true;
            for (std::size_t i = 0; i < n && all; ++i) {
                if (geom::side(pl, fragment_vertex(t, f, i, cache)) != want) all = false;
            }
            if (all) return true;
        }
    }
    return false;
}
#endif

/// スレッド局所に貯めた統計を集約する（`SPEC-phase4.md` §1.1）。
///
/// **和と最大を取り違えないこと。** 最大の項目を足すと、スレッド数に比例して
/// 増える値になり、**決定性の検査を素通りしたまま数字だけが壊れます。**
inline void merge_stats(BoolStats& a, const BoolStats& b) {
    // ---- 和 ----
    a.active_cells += b.active_cells;
    a.empty_cells += b.empty_cells;
    a.early_out_cells += b.early_out_cells;
    a.early_out_raycasts += b.early_out_raycasts;
    a.split_plane_slots += b.split_plane_slots;
    a.split_planes_used += b.split_planes_used;
    a.bsp_cut_slots += b.bsp_cut_slots;
    a.bsp_cuts_used += b.bsp_cuts_used;
    a.bsp_cuts_skipped += b.bsp_cuts_skipped;
    a.regions += b.regions;
    a.raycasts += b.raycasts;
    a.interior.axis_line += b.interior.axis_line;
    a.interior.corner_offset += b.interior.corner_offset;
    a.interior.axis_failed += b.interior.axis_failed;
    a.interior.corner_tries += b.interior.corner_tries;
    a.cache_hits += b.cache_hits;
    a.cache_misses += b.cache_misses;
    a.cache_entries += b.cache_entries;
    a.cache_bytes += b.cache_bytes;
    a.leaf_input_total += b.leaf_input_total;
    a.leaf_nonempty += b.leaf_nonempty;
    a.ray_tri_tests += b.ray_tri_tests;
    // ---- 最大 ----
    a.max_planes_per_cell = std::max(a.max_planes_per_cell, b.max_planes_per_cell);
    a.leaf_input_max = std::max(a.leaf_input_max, b.leaf_input_max);
}

}  // namespace detail

/// 2 つのスープのブール演算（CP2）。**出力もスープなので連鎖できます。**
inline PolySoup boolean(const PolySoup& X, const PolySoup& Y, BoolOp op, const BoolOptions& opt,
                        BoolStats* stats = nullptr) {
    BoolStats st;
    // **段ごとの時間を測ります**（`SPEC-phase4.md` §9）。
    using Clock = std::chrono::steady_clock;
    const auto t_begin = Clock::now();
    auto lap = [](Clock::time_point& t0) {
        const auto t1 = Clock::now();
        const double ms =
            std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(t1 - t0).count();
        t0 = t1;
        return ms;
    };
    auto t_stage = t_begin;

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
    // **符号ベクトルは代表点で取ります**（CP4 で変更）。過剰分割の下では断片全体が
    // 各 source の平面配置の 1 セルに収まっていたので、断片の頂点から取れました。
    // 局所 BSP は**シートが横切らない**ことしか保証しないので、断片が平面をまたぎ得ます。
    //
    // **代表点で取れば、どちらの分割方式でも意味が変わりません。**
    // 分類が断片上で一定であること（＝代表点 1 点で足りること）は、
    // 「断片の相対内部をシートが横切らない」ことから従います。平面をまたいでも、
    // その平面上にシートが無ければ巻き数は変わりません。
    //
    // キーとしての妥当性: 閉多様体なら三角形の辺は隣接三角形との共有辺なので、
    // **平面配置は三角形分割を細分します。** よって同じ符号ベクトルの点は
    // 「同じ三角形の内側／外側」まで一致し、$w$ も $c_{front}, c_{back}$ も一致します。
    std::vector<std::vector<PlaneId>> planes_of_src(n_src);
    // 三角形ごとの平面 ID（退化は `kNoPlane`）。**局所 BSP の切断候補**（§5.4）。
    // `planes_of_src` と違って**重複を潰しません。** どの三角形がどのセルに居るかで
    // 切断候補が決まるので、三角形の添字で引ける必要があります。
    std::vector<std::vector<PlaneId>> tri_plane(n_src);
    for (std::size_t i = 0; i < n_src; ++i) {
        const mesh::TriMesh& m = out.sources[i];
        tri_plane[i].reserve(m.triangles.size());
        for (const mesh::Tri& t : m.triangles) {
            const geom::PlaneD pl =
                geom::plane_from_triangle(m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]]);
            if (geom::is_degenerate(pl)) {
                tri_plane[i].push_back(kNoPlane);
                continue;
            }
            const PlaneId id = out.table.intern(pl).id;
            tri_plane[i].push_back(id);
            planes_of_src[i].push_back(id);
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

    st.ms_prepare = lap(t_stage);

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

    st.ms_leaves = lap(t_stage);

    // ---- 4. セルごとの arrangement（**並列**。`SPEC-phase4.md` §2）-------------
    //
    // **葉は先に列挙されているので、タスクは最初から全部そろっています。**
    // 中央同期キューはアトミックな添字 1 本で足り、**burn-in がありません**
    // （EMBER §5.3 の burn-in は再帰でタスクを生む構造に由来します）。
    //
    // **決定性のために、各葉は自分のスロットにだけ書きます**（§4.2）。
    // 結合は葉の順（`build_leaves` が返す正準な順序）で行うので、
    // **スレッド数に依らず出力はビット単位で同一**です。
    std::vector<Fragment> frags;
    std::vector<octree::Cell> frag_cell;
    std::vector<std::uint32_t> frag_src, frag_tag;
    /// セルで「多角形が 1 枚も無かった source」の内外（-1 = 未確定）。§3.2 の early-out
    std::vector<std::vector<std::int8_t>> frag_forced;

    // **セル面の平面は先に登録します。** 平面表は葉の処理中に伸ばせません
    // （共有される可変状態になります）。セルだけで決まるので前に出せます。
    struct CellPlane {
        PlaneId id;
        int keep;
    };
    std::vector<std::vector<CellPlane>> cell_planes_of(leaves.size());
    for (std::size_t li = 0; li < leaves.size(); ++li) {
        if (leaves[li].depth == 0) continue;
        const auto ps = octree::cell_planes(leaves[li]);
        for (int k = 0; k < 6; ++k) {
            const PlaneRef r = out.table.intern(ps[k]);
            const int base = (k % 2 == 0) ? +1 : -1;
            cell_planes_of[li].push_back({r.id, r.flipped ? -base : base});
        }
    }

    /// 1 つの葉が出すもの。**スロットに書いて、あとで葉の順に結合します。**
    struct LeafOut {
        std::vector<Fragment> frags;
        std::vector<std::uint32_t> src, tag;
        std::vector<std::int8_t> forced;
        bool empty_cell = false;
        bool active = false;
    };
    std::vector<LeafOut> leaf_out(leaves.size());

    // **プールは持ち回します。** 渡されなければこの呼び出しの間だけ作ります
    // （生成コストは 8 スレッドで 0.2 ms。呼び出しが多い場面では `opt.pool` を渡すこと）
    const unsigned nthreads =
        (opt.pool != nullptr) ? opt.pool->size() : ((opt.threads <= 1) ? 1u : opt.threads);
    par::ThreadPool local_pool(opt.pool != nullptr ? 1u : nthreads);
    par::ThreadPool& pool = (opt.pool != nullptr) ? *opt.pool : local_pool;
    // **可変な器はスレッド局所に持ちます**（§1.1）。共有した瞬間に競合が入ります
    std::vector<PointCache> tl_cache(nthreads);
    std::vector<BoolStats> tl_stats(nthreads);

    pool.run(leaves.size(), [&](std::size_t li, unsigned tid) {
        const octree::Cell& cell = leaves[li];
        LeafOut& outl = leaf_out[li];
#if defined(KRISITE_MUTATION_SHARE_STATS)
        // 変異 22: **統計を共有カウンタにする**（§7.5）。TSan が検出します
        BoolStats& st = tl_stats[0];
#else
        BoolStats& st = tl_stats[tid];
#endif
#if defined(KRISITE_MUTATION_SHARE_CACHE)
        // 変異 21: **構成点キャッシュをスレッド間で共有する**（§7.5）。
        // 命中率は上がりますが、`std::map` への同時挿入で壊れます。TSan が検出します
        PointCache* const cache = opt.cache_points ? &tl_cache[0] : nullptr;
#else
        PointCache* const cache = opt.cache_points ? &tl_cache[tid] : nullptr;
#endif
        const octree::CellBox cbox = octree::box_of(cell);
        std::vector<Fragment> local;
        std::vector<std::uint32_t> local_src, local_tag;

        std::vector<std::size_t> here;
        for (std::size_t i = 0; i < polys.size(); ++i) {
            if (octree::assign_to_cell(polys[i].aabb, cbox)) here.push_back(i);
        }
        if (here.empty()) {
            ++st.empty_cells;
            outl.empty_cell = true;
            return;
        }

        // **曲面が 1 枚も横切らない source は、セル全体で内外が一定です。**
        // 隅 1 点（整数点）で決まります（§3.2 の early-out の一般化）。
        // **このセルに居る source の三角形**を、平面ごとに集めます。
        // early-out の判定（`count`）と局所 BSP の切断候補（§5.4）の両方に使います。
        //
        // **閉領域で見ます**（`overlaps_cell`）。割り当てに使う `assign_to_cell` は
        // 重複を避けるため半開区間ですが、ここで問うているのは
        // 「この曲面はこのセルに存在するか」で、**閉じたセルの上での話**です。
        //
        // 半開区間で見ると、面がセルの上側境界にちょうど乗る配置で
        // 「存在しない」と答えてしまい、early-out が**面の上の点を「曲面なし」として
        // 分類します**（`IMPL-phase3.md` §7.1。実際に踏みました）。
        std::vector<std::size_t> count(n_src, 0);
        std::map<PlaneId, std::vector<std::pair<std::uint32_t, std::uint32_t>>> cell_tri_by_plane;
        for (std::size_t i = 0; i < n_src; ++i) {
            for (std::size_t j = 0; j < src_aabb[i].size(); ++j) {
#if defined(KRISITE_MUTATION_HALF_OPEN_PRESENCE)
                // 変異 17: 「存在するか」を**割り当てと同じ半開区間**で見る（実際の誤り）。
                // 面がセルの上側境界に乗る配置でしか出ません。
                if (!octree::assign_to_cell(src_aabb[i][j], cbox)) continue;
#else
                if (!octree::overlaps_cell(src_aabb[i][j], cbox)) continue;
#endif
                ++count[i];
                const PlaneId pid = tri_plane[i][j];
                if (pid != kNoPlane) {
                    cell_tri_by_plane[pid].push_back(
                        {static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j)});
                }
            }
        }
        // **葉に入った三角形の数を記録します**（`SPEC-phase5.md` の CP1.5）。
        // **EMBER §4.5.3 が最適化しているのはこの量**（部分問題の多角形数）で、
        // $P$ でも深度でもありません。ここで持たないと比較になりません。
        {
            std::size_t in_leaf = 0;
            for (std::size_t i = 0; i < n_src; ++i) in_leaf += count[i];
            if (in_leaf > 0) {
                st.leaf_input_total += in_leaf;
                ++st.leaf_nonempty;
                st.leaf_input_max = std::max(st.leaf_input_max, in_leaf);
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
                // **`winding_split` は source の全三角形を走査します**（枝刈りなし）
                st.ray_tri_tests += out.sources[i].triangles.size();
                any = true;
            }
            if (any) ++st.early_out_cells;
        }

        // セル面の平面（保持側つき）。**登録は並列区間の前で済ませています**
        const std::vector<CellPlane>& cps = cell_planes_of[li];

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

        // ---- 局所 BSP の切断集合（`SPEC-phase3.md` §5.4）★ ----------------------
        //
        // 断片 $t$ を切るのは、**このセルに実際に居る三角形**のうち
        // $t$ の支持平面に触れるものの平面だけです。過剰分割は「箱を横切る全平面」で
        // 切っていました（`SPEC-phase1.md` §4.3.1）。
        //
        // **多角形からではなく source の三角形から作ります。** 連鎖では前段で消えた面が
        // あるので、多角形から集めると**共平面シートの境界を切り落とせません**
        // （`planes_of_src` を source から作っているのと同じ理由）。共平面シートの境界は
        // 隣接三角形との共有辺なので、その隣接三角形の平面で切られます。
        //
        // **保守的に平面まるごとで切ります**（EMBER Fig. 7）。交差線分だけを入れると
        // 断片が凸でなくなります。
        //
        // 切断集合は**支持平面だけで決まります。** セル内でメモ化することで、
        // 同じ平面に載る別々の多角形が必ず同じ切り方をすることが保証されます
        // （`region_key` で潰すために要ります）。
        std::map<PlaneId, std::vector<PlaneId>> cuts_memo;
        auto cuts_for = [&](PlaneId sup) -> const std::vector<PlaneId>& {
            const auto it = cuts_memo.find(sup);
            if (it != cuts_memo.end()) return it->second;
            const geom::PlaneD& sp = out.table.at(sup);
            std::vector<PlaneId> cs;
            for (const auto& kv : cell_tri_by_plane) {
                if (kv.first == sup) continue;  // 共平面は §5.4.1（全順序）の担当
                ++st.bsp_cut_slots;
                bool touch = false;
#if defined(KRISITE_MUTATION_BSP_NEVER_CUT)
                // 変異 16: **局所 BSP の切断条件を常に偽にする**（1 枚も切らない）。
                // 断片の相対内部をシートが横切るので、代表点 1 点の分類が壊れます。
#else
                for (const auto& ref : kv.second) {
                    const mesh::TriMesh& m = out.sources[ref.first];
                    if (detail::tri_touches_plane(sp, m, m.triangles[ref.second])) {
                        touch = true;
                        break;
                    }
                }
#endif
                if (touch) {
                    cs.push_back(kv.first);
                    ++st.bsp_cuts_used;
                } else {
                    ++st.bsp_cuts_skipped;
                }
            }
            return cuts_memo.emplace(sup, std::move(cs)).first->second;
        };

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

            // **局所 BSP**（§5.4）。`local_bsp` を切ると過剰分割に戻ります
            // （CP3 までの挙動 = §10.1 の正解器）。
            const std::vector<PlaneId>& cut_planes =
                opt.local_bsp ? cuts_for(frag.support) : split_planes;
            // 変異「局所 BSP の交差線分を 1 本落とす」（`SPEC-phase3.md` §10.5）。
            //
            // **落とすと断片が相手の曲面をまたいだまま残ります。** 代表点 1 点で
            // 分類できる前提（定理 7.2）が崩れるので、位相か体積に出るはずです。
            // **過剰分割の側では落としません**（そちらは正解器なので）。
            std::size_t cut_begin = 0;
#if defined(KRISITE_MUTATION_BSP_DROP_ONE_CUT)
            if (opt.local_bsp && !cut_planes.empty()) cut_begin = 1;
#endif
            std::vector<Fragment> pieces{frag};
            for (std::size_t ci = cut_begin; ci < cut_planes.size(); ++ci) {
                const PlaneId q = cut_planes[ci];
#if defined(KRISITE_EXPERIMENT_BSP_SKIP_DISJOINT)
                // **切断候補の絞り込み（最適化の候補。Phase 5）。変異ではありません。**
                //
                // 三角形が断片と交わらないなら、その断片をその平面で切る必要はありません。
                // **分類については健全です**（定理 7.2 の前提は保たれます）。
                // 全コーパスで出力は 1 ビットも変わらず、8,169 回多く省きました。
                //
                // **採用には証明か専用のコーパスケースが要ります。** 残る経路は
                // 「共平面に載る別々の多角形が違う切り方をして `region_key` が潰せなく
                // なる」で、**「コーパスに配置が無い」だけでは否定できません**
                // （`SPEC-phase2.md` §2.6）。Phase 3 は性能のフェーズではないので保留です。
                {
                    bool touches = false;
                    const auto it = cell_tri_by_plane.find(q);
                    if (it != cell_tri_by_plane.end()) {
                        for (const auto& ref : it->second) {
                            if (!detail::fragment_outside_box(
                                    out.table, frag, src_aabb[ref.first][ref.second], cache)) {
                                touches = true;
                                break;
                            }
                        }
                    }
                    // **変異が発火したことを数えます。** 空回りの変異は変異ではありません
                    if (!touches) {
                        ++st.bsp_cuts_skipped;
                        continue;
                    }
                }
#endif
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

        outl.frags = std::move(local);
        outl.src = std::move(local_src);
        outl.tag = std::move(local_tag);
        outl.forced = forced;
        outl.active = !outl.frags.empty();
        if (outl.active) ++st.active_cells;
    });

    // **葉の順に結合します**（§4.2 の正準な順序）。
    //
    // > **中核の出力は、実はこの順序に依存しません。** 縫合（§5）が構成点を
    // > `lex_less` で整列してから番号を振り直すので、**下流で正準化されます。**
    // > それでも葉の順で結合するのは、**依存していないことを確かめずに済ませない**
    // > ためです（順序を変えても同じ、を偶然に頼らない）。
    // >
    // > **依存が残っているのは分類の結合のほうです**（出力の多角形の並びが変わります）。
    for (std::size_t li = 0; li < leaves.size(); ++li) {
        LeafOut& o = leaf_out[li];
        for (std::size_t i = 0; i < o.frags.size(); ++i) {
            frags.push_back(std::move(o.frags[i]));
            frag_cell.push_back(leaves[li]);
            frag_src.push_back(o.src[i]);
            frag_tag.push_back(o.tag[i]);
            frag_forced.push_back(o.forced);
        }
    }
    for (const BoolStats& t : tl_stats) detail::merge_stats(st, t);
    st.raw_fragments = frags.size();

    st.ms_arrange = lap(t_stage);

    // ---- 5. 縫合（重複の仕分けに要る）----------------------------------------
    //
    // **ここが中核の逐次部分です**（実測 14.7%）。全構成点を大域の表に入れて
    // 値で併合するので、そのままでは並列にできません。**中核だけを並列化したときの
    // 上限は $1/0.161 = 6.2$ 倍**で、これは出口の頂点併合（`SPEC-phase4.md` §3.1）と
    // 同じ形の問題です。**並列整列 + 区分ごとの併合**で崩せますが、CP1 では扱いません。
    PointCache stitch_cache;
    PointCache* const cache = opt.cache_points ? &stitch_cache : nullptr;
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

    st.ms_stitch = lap(t_stage);

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
    using WCache = std::map<std::pair<std::uint32_t, std::vector<std::int8_t>>, Wnd>;

    // **順序非依存の検査**（§14 の CP3 の判定）。分類が可変な共有状態に依存していれば、
    // 逆順にすると結果が変わります。依存していなければ幾何の多重集合は同じです。
    std::vector<const std::pair<const detail::RegionKey, std::vector<std::size_t>>*> region_order;
    region_order.reserve(regions.size());
    for (const auto& kv : regions) region_order.push_back(&kv);
    if (opt.reverse_regions) std::reverse(region_order.begin(), region_order.end());

    // **分類も並列です**（`SPEC-phase4.md` §2）。領域は互いに独立で、
    // source メッシュへの読み取り専用アクセスしかありません。
    //
    // **メモ化（`wcache` / `PointCache`）はスレッド局所に持ちます**（§1.1）。
    // 命中率のために共有した瞬間に競合が入ります。**出力は 1 ビットも変わりません**
    // （キャッシュの有無で結果が変わらないことは Phase 2 で確かめてあります）。
    std::vector<WCache> tl_wcache(nthreads);
    std::vector<PointCache> tl_cache2(nthreads);
    std::vector<BoolStats> tl_stats2(nthreads);
    // **領域ごとのスロット。** 結合は `region_order` の順で行うので、
    // スレッド数に依らず出力はビット単位で同一になります（§4.2）
    std::vector<Poly> region_poly(region_order.size());
    std::vector<char> region_emit(region_order.size(), 0);
#if defined(KRISITE_MUTATION_JOIN_ORDER)
    // 変異 20: **結合時の正準な整列を外す**（`SPEC-phase4.md` §7.5）。
    //
    // スロットに書かず、**完了した順に共有の器へ直接積みます。** 排他は取るので
    // 競合はありませんが、**並びがスレッドのスケジューリングに依存**します。
    // 幾何は同じなので位相も体積も区別しません。**決定性の検査でしか捕まりません。**
    std::mutex emit_mutex;
#endif

    pool.run(region_order.size(), [&](std::size_t ri, unsigned tid) {
        const auto* kvp = region_order[ri];
#if defined(KRISITE_MUTATION_SHARE_STATS)
        BoolStats& st = tl_stats2[0];
#else
        BoolStats& st = tl_stats2[tid];
#endif
#if defined(KRISITE_MUTATION_SHARE_CACHE)
        WCache& wcache = tl_wcache[0];
        PointCache* const cache = opt.cache_points ? &tl_cache2[0] : nullptr;
#else
        WCache& wcache = tl_wcache[tid];
        PointCache* const cache = opt.cache_points ? &tl_cache2[tid] : nullptr;
#endif
        const auto& kv = *kvp;
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
#if defined(KRISITE_EXPERIMENT_NO_WCACHE)
            // **実験: 巻き数のメモ化をやめる**（`SPEC-phase5.md` の CP1.5）。
            //
            // **キーは「代表点が source の全平面のどちら側か」の符号ベクトル**でした。
            // これは点をシーン全体で同定するのに等しいので、**異なる領域はほぼ必ず
            // 異なるキーになります。** 実測で命中率 0%（エントリ 112,084 / 領域 112,084）。
            //
            // それでいて費用は領域あたり $O(\text{シーン})$ です。
            //
            //   空間  8,500 バイトのキー × 70 万件 = 6.1 GB
            //   時間  領域ごとに全平面へ side() = 1.2 × 10^10 回
            //
            // **節約ゼロでシーン規模の費用を両方払う構造**でした。
            Wnd v{};
            winding_split(out.sources[i2], rep, refpl, &v.w_other, &v.c_front, &v.c_back);
            ++st.raycasts;
            ++st.regions;
            st.ray_tri_tests += out.sources[i2].triangles.size();
            w_front[i2] = v.w_other + v.c_front;
            w_back[i2] = v.w_other + v.c_back;
            (void)wcache;
#else
            std::vector<std::int8_t> sig(planes_of_src[i2].size());
            for (std::size_t k = 0; k < planes_of_src[i2].size(); ++k) {
                sig[k] =
                    static_cast<std::int8_t>(geom::side(out.table.at(planes_of_src[i2][k]), rep));
            }
            const auto key = std::make_pair(static_cast<std::uint32_t>(i2), std::move(sig));
            auto it = wcache.find(key);
            if (it == wcache.end()) {
                Wnd v{};
                winding_split(out.sources[i2], rep, refpl, &v.w_other, &v.c_front, &v.c_back);
                ++st.raycasts;
                ++st.regions;
                st.ray_tri_tests += out.sources[i2].triangles.size();
                it = wcache.emplace(key, v).first;
            }
            w_front[i2] = it->second.w_other + it->second.c_front;
            w_back[i2] = it->second.w_other + it->second.c_back;
#endif
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
        if (in_front == in_back) return;  // (in,in) / (out,out) は捨てる（§5.2）

        Poly q;
        q.frag = f;
        // **断片の巻き順は「元の三角形の法線」基準**で、`f.flipped` がそれが支持平面の
        // 法線と逆かを記録しています（`from_mesh` が設定し、切断は保存します）。
        // 出力の巻き順は**外向き法線**基準にしたいので、
        //
        //     反転が要るのは (目標 = in_front) と (現在 = f.flipped) が食い違うとき
        //
        // です。**`in_front` だけで判定すると `f.flipped` が真の面で向きが狂います。**
        // (C, χ) は向きの反転で変わらないので、**体積でしか出ません**（§10.1）。
#if defined(KRISITE_MUTATION_ORIENT_IGNORE_FLIPPED)
        // 変異 15: **`f.flipped` を見ずに `in_front` だけで反転する**（実際の誤り）。
        // 位相は (C, χ) が向きの反転で変わらないので素通りします。**体積でしか出ません。**
        const bool need_reverse = in_front;
#else
        const bool need_reverse = (in_front != f.flipped);
#endif
        q.frag.flipped = in_front;  // 外向き法線が支持平面の法線と逆か
        q.src = frag_src[pick];
        q.tag = frag_tag[pick];
        const octree::CellBox cb2 = octree::box_of(frag_cell[pick]);
        for (int k = 0; k < 3; ++k) {
            q.aabb.lo[k] = cb2.lo[k];
            q.aabb.hi[k] = cb2.hi[k];
        }
        // **立体は in の側にあります。** 外向き法線は立体から外を向くので、
        // +N 側が in なら外向き法線は -N です（§5.2）。
        if (need_reverse) {
            std::vector<std::uint32_t> dummy(vertex_count(q.frag), 0);
            std::vector<PlaneId> e = q.frag.edge;
            detail::reverse_polygon(dummy, e);
            q.frag.edge = std::move(e);
        }
#if defined(KRISITE_MUTATION_JOIN_ORDER)
        {
            const std::lock_guard<std::mutex> g(emit_mutex);
            out.polys.push_back(std::move(q));
        }
        (void)ri;
#else
        region_poly[ri] = std::move(q);
        region_emit[ri] = 1;
#endif
    });

    // **領域の順に結合します**（§4.2）。**ここが出力の多角形の並びを決めます。**
#if !defined(KRISITE_MUTATION_JOIN_ORDER)
    for (std::size_t ri = 0; ri < region_order.size(); ++ri) {
        if (region_emit[ri] != 0) out.polys.push_back(std::move(region_poly[ri]));
    }
#endif
    for (unsigned k = 0; k < nthreads; ++k) {
        tl_stats2[k].cache_hits = tl_cache2[k].hits();
        tl_stats2[k].cache_misses = tl_cache2[k].misses();
        tl_stats2[k].cache_entries = tl_cache2[k].entries();
        tl_stats2[k].cache_bytes = tl_cache2[k].bytes();
        detail::merge_stats(st, tl_stats2[k]);
    }
    st.ms_classify = lap(t_stage);
    st.fragments = out.polys.size();
    if (stats != nullptr) *stats = st;
    return out;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_SOUP_BOOLEAN_HPP

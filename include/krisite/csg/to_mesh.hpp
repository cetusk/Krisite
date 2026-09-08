// Krisite — 出口（`SPEC-phase3.md` §6）
//
//     to_mesh : PolySoup → TriMesh
//
// **Phase 1 / 2 で作った機構がここに移動します。新規実装ではありません。**
//
//   1. 頂点の併合（平面3つ組 + 値ベースの第2段）   SPEC-phase1 §5.3
//   2. T 頂点の解決（案 D）                        SPEC-phase2 §2.4
//   3. 接触の分裂（辺・頂点）                      SPEC-phase2 §5
//   4. 三角形化                                    SPEC-phase1 §2.4.4
//
// **共平面の再併合（§6.4）と snap rounding（§6.5）は CP2 の範囲外です。**
//
// ---
//
// **縫合はここで独立に実装しています。** `boolean.hpp` の二項実装（§10.1 の正解器）が
// 同じ規則を持っていますが、**同じコードを共有すると正解器になりません**
// （`CLAUDE.md`「正解器は被検体と別経路で書く」）。規則は仕様（§5.3）が正です。
#ifndef KRISITE_CSG_TO_MESH_HPP
#define KRISITE_CSG_TO_MESH_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <map>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/tjunction.hpp"
#include "krisite/geom/predicates.hpp"
#include "krisite/mesh/split.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace krisite::csg {

/// 出力メッシュ（構成点 + 三角形）。`boolean.hpp` の `BoolMesh` と同じ形です。
struct SoupMesh {
    std::vector<geom::HPointD> vertices;
    std::vector<mesh::Tri> triangles;
    /// **三角形ごとの由来 source**（`SPEC-phase3.md` §4.3 の由来タグ）。
    ///
    /// `triangles` と同じ長さ。**接触の分裂は三角形の順序と個数を変えない**ので
    /// （`split_contacts` は `out = tris` から始めて置き換えるだけ）、
    /// 分裂の後もそのまま対応します。
    ///
    /// **次数 4 の辺で「4 枚がどの source から来たか」を数えるのに要ります** —
    /// 自己接触（A/A・B/B）か 2 立体の接触（A–B）かは、これでしか分かりません。
    std::vector<int> tri_src;
    /// **三角形ごとの由来タグ**（`SPEC-phase3.md` §4.3。**元の多角形 ID**）。
    ///
    /// `tri_src` は「どちらの入力メッシュか」しか区別しません。
    /// **「元の何番目の三角形か」は、ここでしか分かりません。**
    ///
    /// **異常な辺に接する面が、同じ入力三角形から 3 枚以上出ていないか**を
    /// 調べるのに要ります（`IMPL-phase5.md` §50）。
    std::vector<std::uint32_t> tri_tag;
    /// **三角形ごとの、スープの多角形の添字**（`triangles` と同じ長さ）。
    ///
    /// **これがあると、出力の三角形から支持平面・向き・領域まで辿れます。**
    /// 由来タグ（元の多角形 ID）とは別で、**こちらは演算後の多角形**を指します。
    ///
    /// **異常な辺に接する面が、どの領域から出たか**を調べるのに要ります
    /// （`IMPL-phase5.md` §52）。
    std::vector<std::uint32_t> tri_poly;
    /// **頂点ごとの平面 3 つ組**（`vertices` と同じ長さ）。
    ///
    /// 頂点は「支持平面 1 枚 + 隣り合う辺平面 2 枚」の交点として作られます。
    /// **その 3 つ組が、頂点の同一性の第 1 段の鍵**です。
    ///
    /// **値で併合された頂点では、代表の 3 つ組だけが残ります。**
    /// 何個の 3 つ組がその頂点に落ちたかは `vertex_merged` が持ちます。
    std::vector<std::array<PlaneId, 3>> vertex_key;
    /// **その頂点に落ちた平面 3 つ組の数**（1 なら併合されていない）。
    ///
    /// **2 以上なら「4 枚以上の平面が 1 点で交わった」**ということです。
    /// **別位置の 2 本の辺が 1 本に束ねられていないか**を調べるのに要ります。
    std::vector<std::uint32_t> vertex_merged;
    bool empty() const noexcept { return triangles.empty(); }
};

/// §11 の記録。
struct ToMeshStats {
    /// **段ごとの時間**（`SPEC-phase4.md` §3.2 / §9）。ミリ秒。
    ///
    /// **バリアの数が並列効率の上限を決めます。** 各段の実行時間が偏ると、
    /// バリアで待つ時間が増えます。**まず内訳を測ってから並列化すること。**
    /// **`to_mesh` の入り口から出口まで**（下の 5 段の和と比べるため）。
    ///
    /// > **内部で整合するはずの量を並べて出すと、自分で自分を検査します**
    /// > （`CLAUDE.md`）。**5 段の和と一致しなければ、計時されていない段があります。**
    double ms_total = 0;
    double ms_construct = 0;  ///< 構成点を作る（平面3つ組でメモ化）
    double ms_merge = 0;      ///< 値で併合する（整列 + 区分）
    double ms_index = 0;      ///< 平面ごとの頂点索引（T 解決の下ごしらえ）
    double ms_tri = 0;        ///< T 頂点の解決 + 三角形化
    double ms_split = 0;      ///< 接触の分裂

    std::size_t constructed_points = 0;  ///< 第1段（平面3つ組）で作った点
    std::size_t merged_points = 0;       ///< 第2段（値）の併合後
    std::size_t merged_by_value = 0;
    /// §4.4 の正準化で、**最初に来た点以外が代表になった回数**。
    ///
    /// **0 なら機構が空回りしています**（`CLAUDE.md`「足した機構が実際に発火した
    /// ことを、テスト自身に確かめさせること」）。
    std::size_t canonical_swaps = 0;  ///< 第2段が併合した数
    /// **A-3（セルで区切った索引）の (葉, 支持平面) の組の数**（`ToMeshOptions::cell_index`）。
    ///
    /// **0 なら退避しています**（箱がセルの箱でなかった、または旗が偽）。
    /// **スープ経路（`boolean` の出力）で 0 なら、機構が空回りしています。**
    std::size_t cell_index_groups = 0;
    /// **A-3 の位置決め（頂点をセルに割り振る二分探索）が走った回数**。
    /// **Morton なら不要になります**（`DESIGN-phase5-hotspots.md` §9.3 の恩恵 3）。
    std::size_t cell_index_locate_tests = 0;
    /// **A-3 の (葉, 平面) ごとの照合が走った回数**（索引の本体）。
    std::size_t cell_index_group_tests = 0;
    TJunctionStats t{};
    mesh::SplitStats split{};
};

struct ToMeshOptions {
    /// **radial sort で次数 4 の辺を分ける**（`SPEC-phase2.md` §5.1.2.1）。
    ///
    /// **偽にすると完全に外れ、連結成分だけの従来の挙動に戻ります**（比較の基準側）。
    bool radial_sort = true;
    /// **スレッド数**（`SPEC-phase4.md` §3）。0 か 1 なら逐次。
    ///
    /// 出口は**段ごと + バリア**で並列化します。中核（再帰タスク木）とは
    /// **構造が違います** — work-stealing のプールはそのままでは当たりません。
#if defined(KRISITE_DEFAULT_THREADS)
    unsigned threads = KRISITE_DEFAULT_THREADS;
#else
    unsigned threads = 1;
#endif
    /// 持ち回すプール。`nullptr` なら呼び出しごとに作ります。
    par::ThreadPool* pool = nullptr;
    bool split_contacts = true;  ///< §6.3。既定 ON、フラグで無効化可
    /// **扇の計算を逆順で回す**（`SPEC-phase4.md` §7.5）。`mesh::SplitOptions` へ渡します。
    ///
    /// **出力はバイト単位で変わってはいけません。** 変わるなら、ID の割り当てが
    /// 扇の計算順に依存しています（＝変異 23 と同じ欠陥）。
    /// **スケジューラに依存しない番人**で、1 スレッドでも効きます。
    bool reverse_fan = false;
    bool resolve_t = true;  ///< §6.2 の T 頂点解決
    /// **三角形化の一般解**（`SPEC-phase2.md` §2.4.4 (2)）。**退化三角形を 1 枚も作りません。**
    ///
    /// **偽にすると完全に外れ、従来の扇分割（退化を残す）に戻ります**（比較の基準側）。
    bool general_triangulation = true;
    /// **分裂後の多様体性を検査する**（`SplitOptions::verify_manifold`）。
    ///
    /// **`unresolved` の事後の加算に要ります**（`SPEC-phase5.md` §3.-1 の失敗判定）。
    /// **既定は真です。** 中身は同値な増分計算なので、費用は小さい。
    bool verify_split_manifold = true;
    /// **§5.5 の予測との突き合わせ**（`SplitOptions::verify_delta`）。**純粋な診断**。
    ///
    /// > **既定は偽です**（`HANDOVER.md` §5.1）。**`SPEC-phase5.md` §3.2 は
    /// > CP1〜CP3 で ON を要求している**ので、**実データの駆動プログラムは
    /// > 明示的に真にしてください。**
    bool verify_split_delta = false;
    /// **検証に従来の経路（`check_topology` を 2 回）を使う**（`SplitOptions::verify_naive`）。
    ///
    /// **正解器です。** 既定の増分計算と答えが一致することを検査するために残しています。
    bool verify_split_naive = false;
    /// **T 字接合の索引をセルで区切る**（`DESIGN-phase5-hotspots.md` §6.3 の A-3）。
    ///
    /// 平面ごとに全頂点を走査する代わりに、**多角形が属する葉の【閉じた箱】に
    /// 入る頂点だけ**を候補にします。**厳密な絞り込みで、落ちる T 頂点はありません。**
    ///
    /// **偽にすると完全に外れます**（`CLAUDE.md`「機構を追加したら、それを外す経路も
    /// 用意してください」）。**正しさの検査は真偽の両方で同じ出力を要求します。**
    ///
    /// **効くのはスープ経路だけです。二項メッシュ経路は意図的に素朴なまま**で、
    /// この旗を見ません（正解器は被検体と別経路で書く）。
    bool cell_index = true;
    /// **T 解決の照合で、候補集合を【辺ごとに】走査する**（`SPEC-phase5.md` §5.11）。
    ///
    /// **既定は偽で、多角形あたり 1 回だけ走査します。** 候補集合は
    /// (葉, 支持平面) 群で共有されるので**多角形の中で同じ**だからです。
    ///
    /// **判定は 1 つも変わりません。走査の順序だけが違います。**
    /// **真にすると従来の形に戻ります**（**比較の基準側**）。
    /// **両方で出力がバイト一致することを検査してください。**
    bool scan_per_edge = false;
    /// **候補を X 軸で整列し、多角形の区間を二分探索で絞る**
    /// （`DESIGN-phase5-hotspots.md` §13 の案 (b2)）。
    ///
    /// **厳密な絞り込みです。** 辺は多角形の境界上にあるので、
    /// **多角形の X 区間の外にある候補は、どの辺の相対内部にも載りません。**
    ///
    /// **偽にすると完全に外れます**（整列もしません。**比較の基準側**）。
    /// **出力はバイト単位で変わってはいけません。**
    bool sort_candidates = true;
    /// **案 (b)（多角形の箱で候補を落とす）の効果を数える**（`SPEC-phase5.md` §5.11）。
    ///
    /// **数えるだけで、実際には落としません。** 出力は 1 ビットも変わりません。
    /// **既定は偽です** — 箱を作るのに `cmp_h` を多角形あたり $3(n-1)$ 回払うので、
    /// **計測の費用が小さくありません**（`CLAUDE.md`「計測の費用を本番に持ち込まない」）。
    bool count_box_reject = false;
    /// **段の内訳を計時する**（`TJunctionStats::ms_insert_t` / `ms_fan_tri`）。
    ///
    /// **既定は偽です。** 多角形ごとに時計を 2 回読むので、
    /// **計測の費用を本番の経路に持ち込みません**（`CLAUDE.md`「計測の機構にも
    /// 外す経路を用意してください」）。**測るときだけ真にしてください。**
    ///
    /// > **採れるのは CPU 時間です。** 壁時計に直すには並列効率で割ること。
    bool time_stages = false;
};

/// スープを三角メッシュにする（`SPEC-phase3.md` §6）。
inline SoupMesh to_mesh(const PolySoup& s, const ToMeshOptions& opt = {},
                        ToMeshStats* stats = nullptr) {
    ToMeshStats st;
    using Clock = std::chrono::steady_clock;
    const auto t_enter = Clock::now();
    auto t_stage = t_enter;
    const auto lap = [](Clock::time_point& t0) {
        const auto t1 = Clock::now();
        const double ms =
            std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(t1 - t0).count();
        t0 = t1;
        return ms;
    };
    // **プールは持ち回します**（生成コストは 8 スレッドで 0.2 ms）
    const unsigned nthreads =
        (opt.pool != nullptr) ? opt.pool->size() : ((opt.threads <= 1) ? 1u : opt.threads);
    par::ThreadPool local_pool(opt.pool != nullptr ? 1u : nthreads);
    par::ThreadPool& pool = (opt.pool != nullptr) ? *opt.pool : local_pool;
    SoupMesh out;
    if (s.polys.empty()) {
        st.ms_total = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
                          Clock::now() - t_enter)
                          .count();
        if (stats != nullptr) *stats = st;
        return out;
    }

    // ---- 1. 縫合（§5.3）------------------------------------------------------
    //
    // 第1段: 平面3つ組をキーに引く。第2段: 値が厳密に等しい点を併合する。
    // **4 平面以上が 1 点で交わると 3つ組が違っても同じ点になる**ので、第2段が要ります。
    std::map<std::array<PlaneId, 3>, std::uint32_t> by_key;
    std::vector<geom::HPointD> points;
    std::vector<std::array<PlaneId, 3>> point_key;
    std::vector<std::vector<std::uint32_t>> raw(s.polys.size());

    for (std::size_t pi = 0; pi < s.polys.size(); ++pi) {
        const Fragment& f = s.polys[pi].frag;
        const std::size_t n = vertex_count(f);
        raw[pi].reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            std::array<PlaneId, 3> k{f.support, f.edge[(i + n - 1) % n], f.edge[i]};
            std::sort(k.begin(), k.end());
            auto it = by_key.find(k);
            if (it == by_key.end()) {
                const auto id = static_cast<std::uint32_t>(points.size());
                points.push_back(fragment_vertex(s.table, f, i));
                // **3 つ組を控えます**（`IMPL-phase5.md` §50）。
                // 頂点の同一性の第 1 段の鍵で、**異常な辺の端点を調べるのに要ります**
                point_key.push_back(k);
                it = by_key.emplace(k, id).first;
            }
            raw[pi].push_back(it->second);
        }
    }
    st.constructed_points = points.size();
    st.ms_construct = lap(t_stage);

    std::vector<std::uint32_t> order(points.size());
    for (std::uint32_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
        return geom::lex_less(points[a], points[b]);
    });
    std::vector<std::uint32_t> remap(points.size());
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        const auto id = static_cast<std::uint32_t>(out.vertices.size());
        // **代表は組の中で表現の辞書順が最小のものを選びます**（`SPEC-phase4.md` §4.4）。
        //
        // `lex_less` は同値な点に順序を付けないので、**何もしないと「どれが残るか」が
        // 入力の並び次第**になり、並びに触れる変更のたびに出力のバイト列が漂います。
        // **構成点の集合が同じ限り、これで代表は一意です**（集合が変われば変わり得ます）。
        std::size_t best = i;
        while (j < order.size() && geom::h_equal(points[order[i]], points[order[j]])) {
#if !defined(KRISITE_MUTATION_NO_CANONICAL_REPR)
            // 変異 24: **正準化をやめる**（最初に来た点を代表にする）。
            // 幾何も位相も体積も変わらないので、**並べ替え不変性でしか捕まりません**
            if (geom::repr_less(points[order[j]], points[order[best]])) best = j;
#endif
            remap[order[j]] = id;
            ++j;
        }
        if (best != i) ++st.canonical_swaps;
        out.vertices.push_back(points[order[best]]);
        // **代表の 3 つ組と、この頂点に落ちた 3 つ組の数**（§50）。
        // **2 以上なら 4 枚以上の平面が 1 点で交わっています**
        out.vertex_key.push_back(point_key[order[best]]);
        out.vertex_merged.push_back(static_cast<std::uint32_t>(j - i));
        if (j - i > 1) st.merged_by_value += (j - i - 1);
        i = j;
    }
    st.merged_points = out.vertices.size();
    st.ms_merge = lap(t_stage);

    // ---- 2. T 頂点の解決（§6.2）+ 3. 三角形化 --------------------------------
    PlaneVertexIndex index;
    CellPlaneVertexIndex cell_index;
    bool used_cell_index = false;
    if (opt.cell_index) {
        // **A-3（`DESIGN-phase5-hotspots.md` §6.3）。セルで区切って候補を絞ります。**
        // 実測の削減は 3〜391 倍で、**出力が占めるセルの数**で決まります（`BENCH.md`）。
        std::vector<octree::Aabb> box;
        std::vector<PlaneId> sup;
        // **★ 葉の深度は多角形が持っています**（`Poly::cell_depth`。§5.10.6）。
        // **箱の辺の長さから逆算する形はやめました**（箱を狭めたので逆算できません）。
        std::vector<std::uint8_t> cdep;
        box.reserve(s.polys.size());
        sup.reserve(s.polys.size());
        cdep.reserve(s.polys.size());
        for (const Poly& q : s.polys) {
            box.push_back(q.aabb);
            sup.push_back(q.frag.support);
            cdep.push_back(q.cell_depth);
        }
        used_cell_index = cell_index.build(s.table, out.vertices, box, sup, cdep, &pool,
                                           &st.cell_index_locate_tests, &st.cell_index_group_tests,
                                           opt.sort_candidates);
    }
    if (!used_cell_index) {
        std::vector<PlaneId> sup;
        sup.reserve(s.polys.size());
        for (const Poly& q : s.polys) sup.push_back(q.frag.support);
        std::sort(sup.begin(), sup.end());
        sup.erase(std::unique(sup.begin(), sup.end()), sup.end());
        // **平面ごとに独立**（§3）。**A-3 を外したときの経路です。**
        index.build(s.table, out.vertices, sup, &pool);
    }
    st.cell_index_groups = used_cell_index ? cell_index.groups() : 0;
    st.ms_index = lap(t_stage);

    // **多角形ごとに独立**（§3）。**スロットに書いて、あとで多角形の順に結合します**
    // （§4.2。スレッド数に依らず同じ三角形の列になります）。
    std::vector<std::vector<mesh::Tri>> poly_tris(s.polys.size());
    std::vector<TJunctionStats> tl_t(nthreads);
    pool.run(s.polys.size(), [&](std::size_t pi, unsigned tid) {
        const Fragment& f = s.polys[pi].frag;
        std::vector<std::uint32_t> poly;
        poly.reserve(raw[pi].size());
        for (std::uint32_t v : raw[pi]) poly.push_back(remap[v]);
        if (poly.size() < 3) return;
        std::vector<PlaneId> edge = f.edge;
        TJunctionStats& t = tl_t[tid];
        if (opt.resolve_t) {
            static const std::vector<std::uint32_t> kEmptyCand;
            const std::vector<std::uint32_t>* cand = nullptr;
            if (used_cell_index) {
                cand = &cell_index.candidates(pi);
            } else {
                const std::vector<std::uint32_t>* c = index.find(f.support);
                cand = (c == nullptr) ? &kEmptyCand : c;
            }
            // **計時は旗で囲みます**（`ToMeshOptions::time_stages`）。
            // **本番の経路に時計を持ち込みません。**
            const auto t0 = opt.time_stages ? Clock::now() : Clock::time_point{};
            const TPolygon tp = insert_t_vertices_with(
                s.table, out.vertices, *cand, edge, poly, &t, nullptr, opt.scan_per_edge,
                opt.count_box_reject, used_cell_index && cell_index.sorted());
            const auto t1 = opt.time_stages ? Clock::now() : Clock::time_point{};
            fan_triangulate(tp, poly_tris[pi], &t, opt.general_triangulation);
            if (opt.time_stages) {
                const auto t2 = Clock::now();
                using ms = std::chrono::duration<double, std::milli>;
                t.ms_insert_t += std::chrono::duration_cast<ms>(t1 - t0).count();
                t.ms_fan_tri += std::chrono::duration_cast<ms>(t2 - t1).count();
            }
        } else {
            TPolygon tp;
            tp.corners = static_cast<std::uint32_t>(poly.size());
            tp.vertex = poly;
            tp.is_corner.assign(poly.size(), 1);
            for (std::uint32_t i = 0; i < poly.size(); ++i) tp.orig.push_back(i);
            fan_triangulate(tp, poly_tris[pi], &t, opt.general_triangulation);
        }
    });
    for (std::size_t pi = 0; pi < s.polys.size(); ++pi) {
        for (const mesh::Tri& t : poly_tris[pi]) out.triangles.push_back(t);
        out.tri_src.insert(out.tri_src.end(), poly_tris[pi].size(),
                           static_cast<int>(s.polys[pi].src));
        // **由来タグ（元の多角形 ID）も引き継ぎます**（§4.3）。
        // `tri_src` だけでは「同じ入力三角形から何枚出たか」が分かりません
        out.tri_tag.insert(out.tri_tag.end(), poly_tris[pi].size(), s.polys[pi].tag);
        out.tri_poly.insert(out.tri_poly.end(), poly_tris[pi].size(),
                            static_cast<std::uint32_t>(pi));
    }
    for (const TJunctionStats& t : tl_t) detail::merge_tjunction_stats(st.t, t);

    st.ms_tri = lap(t_stage);

    // ---- 4. 接触の分裂（§6.3）------------------------------------------------
    if (opt.split_contacts && !out.triangles.empty()) {
        std::vector<std::uint32_t> origin;
        // **頂点ごとに独立**（§3）。ID の割り当ては逐次なので決定的です
        mesh::SplitOptions sopt;
        sopt.reverse_fan = opt.reverse_fan;
        // **radial sort に要る幾何を渡します**（`SPEC-phase2.md` §5.1.2.1）。
        //
        // **外向き法線をここで揃えます** — `Fragment::flipped` は
        // 「外向き法線が support の法線と逆か」なので、真なら反転します。
        // **規約を 2 箇所に分けないため、`split_contacts` の側では触りません。**
        std::vector<geom::PlaneD> tri_normal;
        mesh::RadialGeom rg;
        if (opt.radial_sort) {
            tri_normal.resize(out.triangles.size());
            for (std::size_t t = 0; t < out.triangles.size(); ++t) {
                const Poly& q = s.polys[out.tri_poly[t]];
                geom::PlaneD pl = s.table.at(q.frag.support);
                if (q.frag.flipped) {
                    pl.a = arith::neg(pl.a);
                    pl.b = arith::neg(pl.b);
                    pl.c = arith::neg(pl.c);
                    pl.d = arith::neg(pl.d);
                }
                tri_normal[t] = pl;
            }
            rg.vertices = &out.vertices;
            rg.normal = &tri_normal;
        }
        sopt.radial_sort = opt.radial_sort;
        sopt.verify_manifold = opt.verify_split_manifold;
        sopt.verify_delta = opt.verify_split_delta;
        sopt.verify_naive = opt.verify_split_naive;
        out.triangles =
            mesh::split_contacts(out.triangles, out.vertices.size(), &origin, &st.split, nullptr,
                                 nullptr, &pool, sopt, opt.radial_sort ? &rg : nullptr);
        // **分裂で複製された頂点は、元の頂点と同じ位置・同じ 3 つ組**です。
        // **添えた情報も一緒に複製しないと、長さが合わなくなります**
        for (std::uint32_t o : origin) {
            out.vertices.push_back(out.vertices[o]);
            out.vertex_key.push_back(out.vertex_key[o]);
            out.vertex_merged.push_back(out.vertex_merged[o]);
        }
    }

    st.ms_split = lap(t_stage);
    st.ms_total = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
                      Clock::now() - t_enter)
                      .count();
    if (stats != nullptr) *stats = st;
    return out;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_TO_MESH_HPP

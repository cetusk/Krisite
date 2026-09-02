// Krisite — ブール演算のパイプライン
//
// SPEC-phase1.md §4（パイプライン）, §4.3（局所 arrangement）, §5（縫合）, §6（分類）
//
//   1. 平面抽出・ID 付与（§3.1）
//   2. 固定深度で八分木を構築。各セルに三角形を割り当てる（§4.2）
//   3. セルごとに局所 arrangement を計算（§4.3）
//   4. 頂点の同一性を解決して大域メッシュに縫合（§5）
//   5. 断片を正準化し、重複割り当てと共平面重複を仕分ける（§4.2, §4.3.2）
//   6. 領域に分割し、内外を分類（§6）
//   7. 演算に応じて選択・向き付けして出力
//
// **手順 4 は 5 より先です。** 断片の同一性は「辺の平面 ID の列」では決まりません。
// 相異なる 2 平面が支持平面と同一の交線を持ち得るためです（例: 支持平面 z=0 に対し
// x=0 と x+z=0 は同じ交線を与える）。したがって
//
//   - 断片の重複判定
//   - 共平面重複の検出（A の断片と B の断片が同一領域を占めるか）
//
// はどちらも**縫合後の大域頂点 ID の集合**で行います。第2段（値ベースの併合）を
// 通した ID は正準なので、平面 ID の別名に影響されません。
#ifndef KRISITE_CSG_BOOLEAN_HPP
#define KRISITE_CSG_BOOLEAN_HPP

#include <algorithm>
#include <cstdint>
#include <map>
#include <vector>

#include "krisite/csg/faces.hpp"
#include "krisite/csg/fragment.hpp"
#include "krisite/csg/interior.hpp"
#include "krisite/csg/plane_table.hpp"
#include "krisite/csg/raycast.hpp"
#include "krisite/csg/tjunction.hpp"
#include "krisite/mesh/split.hpp"
#include "krisite/mesh/topology.hpp"
#include "krisite/mesh/tri_mesh.hpp"
#include "krisite/octree/adaptive.hpp"
#include "krisite/octree/uniform_grid.hpp"
#include "krisite/par/thread_pool.hpp"

namespace krisite::csg {

enum class BoolOp { Union, Intersection, Difference };

/// 出力メッシュ。頂点は構成点（有理数）なので `mesh::TriMesh` にはできません。
/// 位相検査は三角形の添字だけで行えるので `mesh::check_topology(triangles)` が使えます。
struct BoolMesh {
    std::vector<geom::HPointD> vertices;
    std::vector<mesh::Tri> triangles;
    bool empty() const noexcept { return triangles.empty(); }
};

/// 断片の分類（SPEC-phase1 §4.3.2, §6.1）。
///
/// 共平面重複は「相手の境界に載っている」ことなので、内外のどちらでもありません。
/// 相手の面との**向きの一致・不一致**で扱いが変わるので 2 値に分けます。
enum class FragClass {
    Outside,           ///< 相手の外部
    Inside,            ///< 相手の内部
    CoplanarSame,      ///< 相手の面と同一領域を占め、外向き法線が同じ向き
    CoplanarOpposite,  ///< 同一領域だが外向き法線が逆向き
};

/// SPEC-phase1 §4.3.3 と §5.4 が要求する計数。
struct BoolStats {
    std::size_t fragments = 0;      ///< 正準化後の断片数（§4.3.3）
    std::size_t raw_fragments = 0;  ///< 正準化前（重複を含む）
    /// 幾何を含むセルの数（§9.0）。**分割が働いているかの直接的な指標**です。
    /// 断片数は間接的で、ケース 5 のように面がセル境界に乗ると動きません。
    std::size_t active_cells = 0;
    std::size_t total_cells = 0;  ///< 8^depth
    /// **併合グループの空間的な広がり**（§5.4、第8版）。
    ///
    /// 並列化で問題になるのは頻度ではなく局所性です。値が一致する構成点は幾何的に
    /// 同一の点なので、同じセルか、その点を共有するセルにしか存在し得ません。
    /// **すべての併合グループが 1 セルとその面・辺・頂点隣接に収まる**なら、
    /// Phase 3 は「セル並列 + 境界併合」で組めます。
    ///
    /// セル添字の差の最大（軸ごとの最大）。0 は同一セル内、1 は隣接まで。
    /// **2 以上が出たら併合の誤り**（本来別の点を同一視した）を疑ってください。
    std::size_t max_merge_span = 0;
    std::size_t merge_groups = 0;              ///< 2 点以上が併合されたグループの数
    std::size_t duplicate_fragments = 0;       ///< 重複割り当てが生んだ重複断片（§5.4）
    std::size_t coplanar_same = 0;             ///< 共平面重複のうち向きが同じ対の数
    std::size_t coplanar_opposite = 0;         ///< 向きが逆の対の数
    std::size_t constructed_points = 0;        ///< 第1段が作った構成点の総数（§5.4 の分母）
    std::size_t merged_points = 0;             ///< 第2段の併合後の点数
    std::size_t merged_by_value = 0;           ///< 第1段が取りこぼし第2段が併合した数（§5.4）
    std::size_t max_planes_at_point = 0;       ///< 1 点に集まる平面の最大枚数（§5.4、セル面込み）
    std::size_t max_mesh_planes_at_point = 0;  ///< 同上、メッシュ平面のみ（対照）
    std::size_t planes_total = 0;              ///< 総当たりの分母（表に載った平面の総数）
    std::size_t mesh_planes = 0;               ///< うちメッシュ由来
    std::size_t regions = 0;                   ///< 相異なる符号ベクトルの数（§6.1）
    std::size_t raycasts = 0;                  ///< レイキャスト回数
    /// 代表点の構成（SPEC-phase3 §2.1 の段 0）。**どちらの経路で決まったか。**
    InteriorStats interior{};
    /// `side` と `intersect3` の呼び出し数（SPEC-phase1 §12）。
    ///
    /// **`KRISITE_COUNT_PREDICATES` を定義したビルドでのみ埋まります。**
    /// 既定ビルドでは 0 のままです（計数のコストを本番に持ち込まないため）。
    std::uint64_t side_calls = 0;
    std::uint64_t intersect3_calls = 0;

    /// §2.3 の絞り込み（SPEC-phase2）。
    ///
    /// `split_plane_slots` は「セル × 分割平面」の総当たり数、
    /// `split_planes_used` は絞り込んだ後に実際に使った数です。
    /// **比が Phase 2 の断片数削減の直接的な指標になります**（§11）。
    std::size_t split_plane_slots = 0;
    std::size_t split_planes_used = 0;
    std::size_t max_planes_per_cell = 0;  ///< 1 セルで使った分割平面の最大数

    // ---- 葉の粒度（`SPEC-phase5.md` の CP1.5。**EMBER §4.5.3 と直接比較する量**）----
    //
    // > EMBER: 分割の終盤では、部分問題はより少なく、より大きな多角形を含む。
    // > これは分割の効率を不可避に下げる。**25 多角形に最適点がある。**
    //
    // **EMBER が最適化しているのは「部分問題に入る多角形の数」です。**
    // $P$ でも深度でもないので、この量で持たないと比較になりません。
    /// 葉に入った三角形の総数（両入力の合計。**空でない葉のみ**）
    std::size_t leaf_input_total = 0;
    /// 同上の最大
    std::size_t leaf_input_max = 0;
    /// 空でない葉の数（平均を出すための分母）
    std::size_t leaf_nonempty = 0;
    /// **単一 source の葉**（`SPEC-phase5.md` の (c)）。
    ///
    /// 分割の判定は `(na > 0 && nb > 0)` なので、**片方しか居ない領域は
    /// どれだけ三角形が入っていても割られません。** そこで局所 BSP が全走行すると
    /// $O(P_\text{葉}^2)$ を払います。**NSI フラグ（`SPEC-phase3.md` §5.6）が
    /// あれば省ける費用**なので、いくら払っているかを数えます。
    std::size_t leaf_single_src = 0;
    std::size_t leaf_single_src_input_max = 0;
    std::size_t leaf_both_input_max = 0;
    std::size_t bsp_cut_slots_single = 0;  ///< 単一 source の葉での切断候補
    std::size_t bsp_cuts_used_single = 0;  ///< 同上、実際に切った枚数
    /// **NSI の宣言で局所 BSP を省いたセル**（`SPEC-phase3.md` §5.6）。
    /// **0 なら機構が空回りしています**（`CLAUDE.md`）。
    std::size_t bsp_cells_skipped_nsi = 0;

    /// **レイキャストが検査した三角形の総数**（`SPEC-phase5.md` の CP1.5）。
    ///
    /// `winding_split` は **source の全三角形を走査します**（枝刈りなし）。
    /// したがって 1 レイあたりの検査数は source の三角形数そのものです。
    ///
    /// **これが $O(\text{シーン})$ の実体で、$P^{1.31}$ の出どころです。**
    /// 加速構造を入れれば指数が下がるかどうかの判断材料になります。
    std::size_t ray_tri_tests = 0;

    /// §5.4 の局所 BSP（CP4）。**切断候補のうち何枚を実際に切ったか。**
    ///
    /// `bsp_cut_slots` は「セル内の候補平面 × 支持平面」の総当たり数。
    /// **`bsp_cuts_used` と `bsp_cuts_skipped` の両方が非零であること**を
    /// テストで確かめます。**片方が 0 なら判定が空回りしています**（一方に倒れている）。
    std::size_t bsp_cut_slots = 0;
    std::size_t bsp_cuts_used = 0;     ///< 支持平面に触れるので切った枚数
    std::size_t bsp_cuts_skipped = 0;  ///< 厳密に片側なので切らなかった枚数

    /// §3.1 の適応分割。葉の深度の分布（**深さの差が §2.4 の前提**）。
    unsigned leaf_depth_min = 0;
    unsigned leaf_depth_max = 0;

    /// §3.2 の early-out。
    std::size_t early_out_cells = 0;      ///< 相手が居ないので arrangement を省いたセル
    std::size_t empty_cells = 0;          ///< どちらも居ないので丸ごと飛ばしたセル
    std::size_t early_out_fragments = 0;  ///< そのセルで作った断片（分類を省いた数）
    std::size_t early_out_raycasts = 0;   ///< セルの隅（**整数点**）1 つで済ませた判定

    /// §4 の構成点キャッシュ（§4.4 の記録）。
    std::size_t cache_hits = 0;
    std::size_t cache_misses = 0;
    std::size_t cache_entries = 0;
    std::size_t cache_bytes = 0;

    /// §5 の接触の分裂。**予測と実測の突き合わせもここに入ります**（§5.5）。
    mesh::SplitStats split{};

    /// **段ごとの時間**（`SPEC-phase4.md` §9）。ミリ秒。
    ///
    /// **並列化の対象を選ぶために測ります。** どこに時間が行くかを知らずに
    /// 並列化しても無駄になります（`CLAUDE.md`「最適化の前に比率を測る」）。
    double ms_prepare = 0;   ///< 平面表の統合、source ごとの平面・AABB
    double ms_leaves = 0;    ///< 葉の列挙（八分木の構築）
    double ms_arrange = 0;   ///< セルごとの arrangement（§4）
    double ms_stitch = 0;    ///< 縫合と重複の仕分け（§5 / §6）
    double ms_classify = 0;  ///< 分類（§7）

    /// §2.4.3 の T 頂点の解決（SPEC-phase2）。
    ///
    /// **固定深度では `t_inserted` が 0 でなければなりません**（§2.4.4 (3) の負の対照）。
    /// 大域平面集合で切るので、線分の内部に載る頂点はその平面での切断点として既に
    /// 多角形の角になっているためです。**0 でなければ前提が崩れています。**
    TJunctionStats t{};
};

namespace detail {

/// 多角形の頂点順を反転し、辺の対応も付け替える。
///
/// `Difference` で採用した B の断片に使います（SPEC-phase1 §3.4）。辺 j は頂点 j と
/// j+1 を結ぶので、反転すると `edge'[k] = edge[(2n-2-k) mod n]` になります。
/// **T 頂点を入れる前に反転してください。** 入れた後だと `TPolygon::orig` の
/// 意味（元の辺の添字）が保てません。
inline void reverse_polygon(std::vector<std::uint32_t>& poly, std::vector<PlaneId>& edge) {
    const std::size_t n = poly.size();
    KRISITE_CHECK(n == edge.size(), "reverse_polygon: 頂点数と辺数が違う");
    std::vector<std::uint32_t> p(n);
    std::vector<PlaneId> e(n);
    for (std::size_t k = 0; k < n; ++k) {
        p[k] = poly[n - 1 - k];
        e[k] = edge[(2 * n - 2 - k) % n];
    }
    poly.swap(p);
    edge.swap(e);
}

}  // namespace detail

namespace detail {

/// 頂点の平面 3 つ組（昇順に正規化）。§5.3 の第1段のキー。
inline std::array<PlaneId, 3> vertex_key(const Fragment& f, std::size_t i) {
    const std::size_t n = f.edge.size();
    std::array<PlaneId, 3> k{f.support, f.edge[(i + n - 1) % n], f.edge[i]};
    std::sort(k.begin(), k.end());
    return k;
}

/// 断片の正準キー。支持平面と、縫合後の大域頂点 ID を昇順に並べたもの。
///
/// 凸多角形は頂点集合で一意に定まるので、順序を捨てても領域は復元できます。
using RegionKey = std::pair<PlaneId, std::vector<std::uint32_t>>;

inline RegionKey region_key(PlaneId support, std::vector<std::uint32_t> ids) {
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return RegionKey{support, std::move(ids)};
}

/// 演算と分類から、その断片を出力に採るかを決める（SPEC-phase1 §2.3 の正則化）。
///
/// 共平面重複は「A 側だけを残す」ことで二重出力を避けます。どちらを残すかは
/// 恣意的ですが、向きが問題になる `Difference` では A の向きがそのまま答えなので
/// A に固定するのが自然です。
///
///  | | ∪ | ∩ | A\B |
///  |同方向| A のみ | A のみ | 両方落とす |
///  |逆方向| 両方落とす | 両方落とす | A のみ |
///
/// 同方向の面は「両者の外側が同じ側」なので、∪ と ∩ の境界には残り、A\B では消えます。
/// 逆方向の面は「一方の外側が他方の内側」なので、∪ と ∩ では内部面になって消え、
/// A\B では A の表面として残ります。
inline bool select_fragment(BoolOp op, int owner, FragClass c) noexcept {
    switch (c) {
        case FragClass::Outside:
            return (op == BoolOp::Union)          ? true
                   : (op == BoolOp::Intersection) ? false
                                                  : (owner == 0);
        case FragClass::Inside:
            return (op == BoolOp::Union)          ? false
                   : (op == BoolOp::Intersection) ? true
                                                  : (owner == 1);
        case FragClass::CoplanarSame:
            return owner == 0 && op != BoolOp::Difference;
        case FragClass::CoplanarOpposite:
            return owner == 0 && op == BoolOp::Difference;
    }
    return false;
}

}  // namespace detail

/// ブール演算の設定。**すべて実行時パラメータです**（SPEC-phase1 §2.2、SPEC-phase2 §3.1）。
///
/// **既定は Phase 1 の意味論そのもの**（固定深度・early-out なし）で、絞り込みだけが
/// 有効です。適応分割を有効にした側と無効にした側を同一プロセス内で比較すれば、
/// §9.1 の分割戦略不変性がそのまま検査になります（§0.1）。
struct BoolOptions {
    unsigned depth = 0;  ///< 最大深度
    /// §2.3 の分割平面の絞り込み。**無効側が Phase 1 の挙動**（§0.1 の正解器）
    bool cull_planes = true;
    /// **レイキャストの 2 次元索引**（`SPEC-phase5.md` §6.3、CP1.5 の 2）。
    ///
    /// レイは軸平行なので、投影面のセルで候補を絞れます。**厳密な絞り込み**で、
    /// 落ちる三角形はありません（根拠は `ray_index.hpp` 冒頭の単調性）。
    ///
    /// **偽にすると完全に外れます。** 正しさの検査は真偽の両方で同じ出力を
    /// 要求します（`CLAUDE.md`「機構を追加したら、それを外す経路も用意してください」）。
    ///
    /// **効くのはスープ経路（`soup_boolean.hpp`）だけです。**
    /// **二項メッシュ経路（`boolean_op`）は意図的に素朴なまま**で、この旗を見ません
    /// （正解器は被検体と別経路で書く。`IMPL-phase5.md` §12）。
    bool ray_index = true;
    /// **適応分割**（§3.1）。偽なら常に最大深度まで分割する固定深度モード。
    /// **固定深度モードを消さないこと。** §9.1 の正解器です
    ///
    /// `KRISITE_DEFAULT_ADAPTIVE` を定義すると既定が反転します（§9.4 の CI ジョブ用）。
    /// **明示的に指定すれば従来どおり**なので、比較のテストは影響を受けません。
#if defined(KRISITE_DEFAULT_ADAPTIVE)
    bool adaptive = true;
#else
    bool adaptive = false;
#endif
    /// 分割を打ち切る三角形数の閾値（§3.1）。0 なら閾値では打ち切らない
    std::size_t leaf_threshold = 0;
    /// **接触の分裂**（§5）。**既定で有効です**（§5.2）。
    ///
    /// 正則化ブールは一般に多様体出力を保証できません。辺だけ・頂点だけを共有する接触が
    /// 出力に残るので、分裂させて多様体化します。Phase 5 のメッシュ化と GWN が多様体を
    /// 前提にできると扱いが楽なため、既定を分裂側にしています。
    ///
    /// **無効にしたときは $g$ で比較してはいけません**（§5.4）。非多様体出力では
    /// $\chi$ が奇数になり得ます。
    bool split_contacts = true;
    /// **構成点の保持**（§4）。平面3つ組をキーにメモ化する。
    ///
    /// **無効側が Phase 1 の挙動 = §9.1 の正解器です。** キャッシュの有無で出力は
    /// 1 ビットも変わらないはずなので、比較は値の完全一致で行えます。
#if defined(KRISITE_DEFAULT_ADAPTIVE)
    bool cache_points = true;
#else
    bool cache_points = false;
#endif
    /// **early-out**（§3.2）。相手の三角形が 1 つも無いセルで arrangement と分類を省く。
    ///
    /// **無効側が正解器です。** 有効・無効を同一プロセス内で比較すれば、省いたことで
    /// 答えが変わっていないことが直接検査できます（§0.1 と同じ構図）。
#if defined(KRISITE_DEFAULT_ADAPTIVE)
    bool early_out = true;
#else
    bool early_out = false;
#endif
    /// **局所 BSP**（`SPEC-phase3.md` §5.4。CP4）。**既定で有効です。**
    ///
    /// 断片を切るのは「このセルに居る三角形のうち、支持平面に触れるもの」の平面だけ。
    /// **無効側が過剰分割 = CP3 までの挙動で、§10.1 の正解器です。**
    ///
    /// 一致するのは $(C, \chi)$ と体積で、**三角形の集合と断片数は一致しません**
    /// （それが目的です）。比較のテストでは**両側を明示的に指定してください。**
    ///
    /// **`soup_boolean` でのみ効きます**（二項の `boolean` は過剰分割のままです）。
    bool local_bsp = true;
    /// **スレッド数**（`SPEC-phase4.md` §2 / §7.1）。**0 か 1 なら逐次**。
    ///
    /// **逐次経路はそのまま残ります。** `parallel_for` は `threads <= 1` で
    /// スレッドを作らないので、§7.1 の「逐次実装との一致」を同一プロセスで
    /// 比較できます（`SPEC-phase2.md` §0.1 と同じ構図）。
    ///
    /// **出力はスレッド数に依らずビット単位で同一**でなければなりません（§4）。
    /// **`soup_boolean` でのみ効きます**（二項の `boolean_op` は逐次のままです）。
#if defined(KRISITE_DEFAULT_THREADS)
    unsigned threads = KRISITE_DEFAULT_THREADS;
#else
    unsigned threads = 1;
#endif
    /// **持ち回すスレッドプール**（`SPEC-phase4.md` §2）。`nullptr` なら呼び出しごとに作ります。
    ///
    /// **呼び出しが多い場面では必ず渡してください。** 生成・破棄のコストは
    /// 8 スレッドで 1 回あたり 0.2 ms あり、**タスクが小さいと支配します**
    /// （実測 132 回で 26.2 ms = 中核の 40%。`IMPL-phase4.md` §1.2）。
    par::ThreadPool* pool = nullptr;
    /// **領域を逆順に分類する**（`SPEC-phase3.md` §14 の CP3 の判定）。
    ///
    /// 葉の分類が可変な共有状態に依存していれば、順序を変えると結果が変わります。
    /// **並列化できるかを機構で確かめるためのフラグ**です。出力の順序は変わりますが、
    /// **幾何の多重集合は変わらないはず**です。
    bool reverse_regions = false;
};

/// ブール演算。`depth` は八分木の深度（実行時パラメータ。SPEC-phase1 §2.2）。
///
/// `cull_planes` は §2.3（SPEC-phase2）の分割平面の絞り込み。**既定で有効です。**
///
/// **無効にした側が Phase 1 の挙動そのもの**なので、両者の出力を同一プロセス内で
/// 比較すれば §9.1 の分割戦略不変性がそのまま検査になります（§0.1）。
/// 深度と同じく**実行時パラメータにしてあります。** コンパイル時定数にすると
/// 1 回の実行で比較できません。
inline BoolMesh boolean_op(const mesh::TriMesh& A, const mesh::TriMesh& B, BoolOp op,
                           const BoolOptions& opt, BoolStats* stats = nullptr) {
    const unsigned depth = opt.depth;
    const bool cull_planes = opt.cull_planes;
    BoolStats st;
    // §4.2: **キャッシュはグローバルに持ちません。** ここが「CSG の文脈オブジェクト」で、
    // 以下すべての呼び出しに明示的に引き回します（`STYLE.md` と Phase 3 の並列化のため）。
    PointCache point_cache;
    PointCache* const cache = opt.cache_points ? &point_cache : nullptr;
#if defined(KRISITE_COUNT_PREDICATES)
    geom::counters::reset();
#endif

    // ---- 1. 平面抽出・ID 付与 ----
    PlaneTable table;
    const std::vector<Face> faces_a = build_faces(A, 0, table);
    const std::vector<Face> faces_b = build_faces(B, 1, table);
    const std::size_t n_mesh_planes = table.size();
    std::vector<PlaneId> mesh_planes(n_mesh_planes);
    for (std::size_t i = 0; i < n_mesh_planes; ++i) mesh_planes[i] = static_cast<PlaneId>(i);

    // 各面の AABB（セル割り当て用）
    auto face_aabb = [&](const mesh::TriMesh& m, const Face& f) {
        octree::Aabb r{};
        for (int t = 0; t < 3; ++t) {
            r.lo[t] = krisite::kCoordMax;
            r.hi[t] = krisite::kCoordMin;
        }
        for (mesh::VertexId vid : f.loop) {
            const geom::IPoint& p = m.vertices[vid];
            const std::int64_t c[3] = {p.x, p.y, p.z};
            for (int t = 0; t < 3; ++t) {
                r.lo[t] = std::min(r.lo[t], c[t]);
                r.hi[t] = std::max(r.hi[t], c[t]);
            }
        }
        return r;
    };

    // ---- 2 + 3. セルごとに局所 arrangement ----
    //
    // **ここでは重複を落としません。** 同一領域の断片が A と B の両方から出てくるのが
    // 共平面重複そのものなので、辺の平面 ID で潰すと区別がつかなくなります（§4.3.2）。
    std::vector<Fragment> frags;
    // §5.4: 併合の局所性を測るため、断片がどのセルで生まれたかを覚えておく
    std::vector<octree::Cell> frag_cell;
    /// §3.2 の early-out で分類が確定している断片（-1 = 未確定）。
    std::vector<std::int8_t> frag_forced;

    // 面の AABB は分割判定でも割り当てでも使うので、先に作っておく
    std::vector<octree::Aabb> aabb_a(faces_a.size()), aabb_b(faces_b.size());
    for (std::size_t i = 0; i < faces_a.size(); ++i) aabb_a[i] = face_aabb(A, faces_a[i]);
    for (std::size_t i = 0; i < faces_b.size(); ++i) aabb_b[i] = face_aabb(B, faces_b[i]);

    // §3.1 の分割判定で葉を列挙する。**固定深度は「常に最大深度まで分割する」特別な場合**
    const octree::SubdivisionPolicy policy{depth, !opt.adaptive, opt.leaf_threshold};
    const std::vector<octree::Cell> leaves =
        octree::build_leaves(policy, [&](const octree::Cell& c, std::size_t* na, std::size_t* nb) {
            const octree::CellBox cb = octree::box_of(c);
            *na = 0;
            *nb = 0;
            for (const octree::Aabb& r : aabb_a) {
                if (octree::assign_to_cell(r, cb)) ++*na;
            }
            for (const octree::Aabb& r : aabb_b) {
                if (octree::assign_to_cell(r, cb)) ++*nb;
            }
        });
    st.leaf_depth_min = depth;
    for (const octree::Cell& c : leaves) {
        st.leaf_depth_min = std::min(st.leaf_depth_min, c.depth);
        st.leaf_depth_max = std::max(st.leaf_depth_max, c.depth);
    }

    {
        {
            for (const octree::Cell& cell : leaves) {
                const octree::CellBox cbox = octree::box_of(cell);
                const std::size_t frags_before = frags.size();

                // §3.2 の early-out。**相手が居ないセルは arrangement を計算しません。**
                //
                // $C$ が $B$ の三角形を 1 つも含まないなら、$C$ 全体が $B$ の内側か外側かの
                // **どちらか一方**です。セルの隅 1 つで判定して、そのセルの $A$ の断片
                // すべてに同じ分類を与えれば済みます。
                //
                // **判定点は `lo` 隅（整数点）です。** `side(plane, IPoint)` は 2.35 ns、
                // `side(plane, HPoint)` は 7.80 ns（`BENCH.md`）。構成点を一切作りません。
                // `hi` 隅は `+2^(b-1)` になり得て `IPoint` に入りませんが、`lo` 隅は
                // 必ず範囲内です（`lo <= kCoordMax + 1 - セル幅`）。
                //
                // **`lo` 隅が $\partial B$ の上に無いことは、割り当てが保証します。**
                // $\partial B$ が隅を含むなら、その面の AABB は隅を含むので割り当てられ、
                // $n_B > 0$ になります。したがって `point_inside` の契約は満たされます。
                //
                // **閉領域で数えます**（`overlaps_cell`）。割り当ての `assign_to_cell` は
                // 重複を避けるため半開区間ですが、ここで問うているのは
                // 「相手の曲面がこのセルに存在するか」で、**閉じたセルの上での話**です。
                // 半開区間で数えると、面がセルの**上側境界**にちょうど乗る配置で
                // 「存在しない」と答え、early-out が面の上の点を「曲面なし」として
                // 分類します（`IMPL-phase3.md` §7.1）。
                std::size_t na = 0, nb = 0;
                for (const octree::Aabb& r : aabb_a) {
                    if (octree::overlaps_cell(r, cbox)) ++na;
                }
                for (const octree::Aabb& r : aabb_b) {
                    if (octree::overlaps_cell(r, cbox)) ++nb;
                }
                // -1 = 分類を省かない、0 = 相手の外、1 = 相手の内
                int forced[2] = {-1, -1};
                if (opt.early_out) {
#if defined(KRISITE_MUTATION_LOOSE_EARLY_OUT)
                    // SPEC-phase2 §9.3 の変異 6: 判定を緩め、**相手が居るセルでも省く。**
                    // 相手の境界がセルを横切っているので、隅 1 点の判定は断片ごとの正しい
                    // 分類と一致しません。**分類の誤りとして位相・体積に出るはずです。**
                    const bool skip_a = true, skip_b = true;
#else
                    const bool skip_a = (nb == 0), skip_b = (na == 0);
#endif
                    if (na == 0 && nb == 0) {
                        ++st.empty_cells;
                        continue;  // どちらも居ないセルは出力に寄与しません
                    }
                    const geom::IPoint corner{static_cast<std::int32_t>(cbox.lo[0]),
                                              static_cast<std::int32_t>(cbox.lo[1]),
                                              static_cast<std::int32_t>(cbox.lo[2])};
                    if (na > 0 && skip_a) {
                        forced[0] = point_inside(B, corner) ? 1 : 0;
                        ++st.early_out_raycasts;
                    }
                    if (nb > 0 && skip_b) {
                        forced[1] = point_inside(A, corner) ? 1 : 0;
                        ++st.early_out_raycasts;
                    }
                    if (forced[0] >= 0 || forced[1] >= 0) ++st.early_out_cells;
                }
                // セル面の平面 ID（保持側つき）: lo 面は +、hi 面は -
                struct CellPlane {
                    PlaneId id;
                    int keep;
                };
                std::vector<CellPlane> cps;
                if (cell.depth > 0) {
                    const auto ps = octree::cell_planes(cell);
                    for (int k = 0; k < 6; ++k) {
                        const PlaneRef r = table.intern(ps[k]);
                        // 平面の代表が裏返っているなら保持側も反転する
                        const int base = (k % 2 == 0) ? +1 : -1;
                        cps.push_back({r.id, r.flipped ? -base : base});
                    }
                }

                // §10.5 の変異 3: 分割平面を【割り当て集合】に変える。
                // 大域集合を使う理由は §4.3.1（継ぎ目の T 字接合）。
#if defined(KRISITE_MUTATION_LOCAL_PLANES)
                std::vector<PlaneId> local_planes;
                if (cell.depth > 0) {
                    for (int which = 0; which < 2; ++which) {
                        const mesh::TriMesh& lm = (which == 0) ? A : B;
                        for (const Face& f : (which == 0) ? faces_a : faces_b) {
                            if (!octree::assign_to_cell(face_aabb(lm, f), cbox)) continue;
                            if (std::find(local_planes.begin(), local_planes.end(), f.support) ==
                                local_planes.end()) {
                                local_planes.push_back(f.support);
                            }
                        }
                    }
                    std::sort(local_planes.begin(), local_planes.end());
                }
                const std::vector<PlaneId>& all_split =
                    (cell.depth > 0) ? local_planes : mesh_planes;
#else
                const std::vector<PlaneId>& all_split = mesh_planes;
#endif
                // SPEC-phase2 §2.3: **平面がセルの閉包と交わるときだけ切る。**
                //
                // 「三角形が届くか」で絞ってはいけません。plane(T) は無限に延びるので
                // T が届かないセルにも切断点を生みます（Phase 1 の変異 3）。
                //
                // **閉包で判定するので継ぎ目は割れません。** 共有面 F 上に切断点を
                // 生む平面は F を横切り、F は両セルの閉包に含まれるので、両セルが
                // 同じ平面で切ります（§2.3 の証明）。
                //
                // **深度 0 でも意味があります。** セルは 1 個ですが、座標範囲の外に
                // ある平面は落ちます。
                const std::int64_t* clo = cbox.lo;
                const std::int64_t* chi = cbox.hi;
                std::vector<PlaneId> culled;
                if (cull_planes) {
                    culled.reserve(all_split.size());
#if defined(KRISITE_MUTATION_REACHING_TRIANGLE)
                    // SPEC-phase2 §9.3 の変異 4: 判定基準を
                    // 「**平面が**セルを横切るか」から「**三角形が**セルに届くか」に戻す。
                    //
                    // Phase 1 の変異 3 と同型ですが、こちらは絞り込みの**基準**を
                    // 差し替えます。plane(T) は無限に延びるので、T が届かないセルにも
                    // 切断点を生みます。落とすと継ぎ目に T 字接合が出るはずです。
                    (void)clo;
                    (void)chi;
                    for (int which = 0; which < 2; ++which) {
                        const mesh::TriMesh& rm = (which == 0) ? A : B;
                        for (const Face& f : (which == 0) ? faces_a : faces_b) {
                            if (!octree::assign_to_cell(face_aabb(rm, f), cbox)) continue;
                            if (std::find(culled.begin(), culled.end(), f.support) ==
                                culled.end()) {
                                culled.push_back(f.support);
                            }
                        }
                    }
                    std::sort(culled.begin(), culled.end());
#else
                    for (PlaneId q : all_split) {
                        if (geom::plane_crosses_box(table.at(q), clo, chi)) culled.push_back(q);
                    }
#endif
                }
                const std::vector<PlaneId>& split_planes = cull_planes ? culled : all_split;
                st.split_plane_slots += all_split.size();
                st.split_planes_used += split_planes.size();
                st.max_planes_per_cell = std::max(st.max_planes_per_cell, split_planes.size());

                for (int which = 0; which < 2; ++which) {
                    const mesh::TriMesh& m = (which == 0) ? A : B;
                    const std::vector<Face>& fs = (which == 0) ? faces_a : faces_b;
                    for (const Face& f : fs) {
                        // §10.5 の変異 2: 割り当てを開領域に変える。
                        // 閉領域で行う理由は §4.2（共有面上の三角形が両側から落ちる）。
#if defined(KRISITE_MUTATION_OPEN_CELLS)
                        if (cell.depth > 0 && !octree::assign_to_cell_open(face_aabb(m, f), cbox)) {
                            continue;
                        }
#else
                        if (cell.depth > 0 && !octree::assign_to_cell(face_aabb(m, f), cbox)) {
                            continue;
                        }
#endif
                        Fragment frag = face_to_fragment(f);
                        // セルの 6 面でクリップ
                        bool alive = true;
                        for (const CellPlane& cp : cps) {
                            if (cp.id == frag.support) continue;
                            if (!clip_fragment(table, frag, cp.id, cp.keep, cache)) {
                                alive = false;
                                break;
                            }
                        }
                        if (!alive) continue;

                        // 【両メッシュの全平面】で分割（§4.3.1）。
                        //
                        // **early-out したセルでは分割しません**（§3.2「arrangement を
                        // 計算しない」）。継ぎ目は §2.4.3 の T 解決が埋めます。
                        std::vector<Fragment> pieces{frag};
                        const std::vector<PlaneId> no_split;
                        for (PlaneId q : (forced[which] >= 0) ? no_split : split_planes) {
                            std::vector<Fragment> next;
                            next.reserve(pieces.size());
                            for (const Fragment& p : pieces) {
                                if (q == p.support) {
                                    next.push_back(p);
                                    continue;
                                }
                                const SplitResult r = split_fragment(table, p, q, cache);
                                if (r.has_pos) next.push_back(r.pos);
                                if (r.has_neg) next.push_back(r.neg);
                            }
                            pieces.swap(next);
                        }
                        for (Fragment& p : pieces) {
                            frags.push_back(std::move(p));
                            frag_cell.push_back(cell);
                            frag_forced.push_back(static_cast<std::int8_t>(forced[which]));
                            if (forced[which] >= 0) ++st.early_out_fragments;
                        }
                    }
                }
                if (frags.size() != frags_before) ++st.active_cells;
            }
        }
    }
    st.raw_fragments = frags.size();
    st.total_cells = leaves.size();

    // ---- 4. 縫合（§5）----
    //
    // 第1段: 平面3つ組をキーに引く
    // 第2段: 全構成点を lex_less で整列し、値が厳密に等しいものを併合して再写像する
    //
    // **選択より先に、全断片について行います。** 断片の正準化（手順 5）が縫合後の
    // 頂点 ID を使うためです。
    std::map<std::array<PlaneId, 3>, std::uint32_t> by_key;
    std::vector<geom::HPointD> points;
    /// **その構成点が「覚えていた側」から返されたか**（§13 の CP5 の相互作用）。
    std::vector<char> point_from_cache;
    // 構成点ごとの、生成に関わったセル添字の範囲（§5.4 の局所性）
    struct CellRange {
        std::uint32_t lo[3], hi[3];
    };
    std::vector<CellRange> point_cells;

    auto vertex_id = [&](const Fragment& f, std::size_t i, const octree::Cell& c) {
        const auto k = detail::vertex_key(f, i);
        // **深さが混ざると添字はそのままでは比べられません。** 最大深度の格子に写します。
        // 固定深度では恒等写像なので、Phase 1 の数値がそのまま再現されます
        std::uint32_t ci[3];
        octree::normalized_index(c, depth, ci);
        auto it = by_key.find(k);
        if (it != by_key.end()) {
            CellRange& r = point_cells[it->second];
            for (int t = 0; t < 3; ++t) {
                r.lo[t] = std::min(r.lo[t], ci[t]);
                r.hi[t] = std::max(r.hi[t], ci[t]);
            }
            return it->second;
        }
        const auto id = static_cast<std::uint32_t>(points.size());
        const geom::HPointD v = fragment_vertex(table, f, i, cache);
        KRISITE_CHECK(arith::sign(v.w) != 0,
                      "boolean_op: 構成点の w が 0（3 平面が一点で交わっていない）");
        points.push_back(v);
        // **`fragment_vertex` を呼んだ「あと」に見ます。** 分割の途中で作られた点が
        // ここで再利用されたときだけ旗が立ちます（それが CP5 の見たい経路です）。
        point_from_cache.push_back(cache != nullptr && cache->served_from_cache(k) ? 1 : 0);
        point_cells.push_back(CellRange{{ci[0], ci[1], ci[2]}, {ci[0], ci[1], ci[2]}});
        by_key.emplace(k, id);
        return id;
    };

    std::vector<std::vector<std::uint32_t>> raw_polys(frags.size());
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        const Fragment& f = frags[fi];
        raw_polys[fi].reserve(vertex_count(f));
        for (std::size_t i = 0; i < vertex_count(f); ++i) {
            raw_polys[fi].push_back(vertex_id(f, i, frag_cell[fi]));
        }
    }
    st.constructed_points = points.size();

    // 第2段: 値ベースの併合
    std::vector<std::uint32_t> order(points.size());
    for (std::uint32_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
        return geom::lex_less(points[a], points[b]);
    });
    std::vector<std::uint32_t> remap(points.size());
    std::vector<geom::HPointD> merged;
    std::vector<char> merged_from_cache;
#if !defined(KRISITE_MUTATION_NO_STAGE2)
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        const auto id = static_cast<std::uint32_t>(merged.size());
        merged.push_back(points[order[i]]);
        merged_from_cache.push_back(0);
        std::uint32_t lo[3] = {point_cells[order[i]].lo[0], point_cells[order[i]].lo[1],
                               point_cells[order[i]].lo[2]};
        std::uint32_t hi[3] = {point_cells[order[i]].hi[0], point_cells[order[i]].hi[1],
                               point_cells[order[i]].hi[2]};
        while (j < order.size() && geom::h_equal(points[order[i]], points[order[j]])) {
            remap[order[j]] = id;
            if (point_from_cache[order[j]] != 0) merged_from_cache[id] = 1;
            for (int t = 0; t < 3; ++t) {
                lo[t] = std::min(lo[t], point_cells[order[j]].lo[t]);
                hi[t] = std::max(hi[t], point_cells[order[j]].hi[t]);
            }
            ++j;
        }
        if (j - i > 1) {
            st.merged_by_value += (j - i - 1);
            ++st.merge_groups;
        }
        // §5.4: 併合グループの空間的な広がり。同一セルなら 0、面・辺・頂点隣接なら 1
        for (int t = 0; t < 3; ++t) {
            st.max_merge_span = std::max(st.max_merge_span, std::size_t{hi[t] - lo[t]});
        }
        i = j;
    }
#else
    // §10.5 の変異 1: 第2段を無効化する
    merged = points;
    merged_from_cache = point_from_cache;
    for (std::uint32_t i = 0; i < remap.size(); ++i) remap[i] = i;
#endif
    st.merged_points = merged.size();

    // §2.4.3 の手順 1: **辺の平面対 → その交線上の大域頂点** の索引。
    //
    // 3つ組をそのまま登録すれば 3 つの対が張れます。線上に載ることは 3つ組から構造的に
    // 決まるので、幾何的な判定は「区間の内部か」だけで済みます。
    //
    // **平面3つ組をキーにしてはいけません。** 3つ組は正準ではなく、相異なる 2 平面が
    // 支持平面と同一の交線を与え得ます（`IMPL-phase1.md` §2.9）。別名で記録された頂点を
    // 取りこぼし、T 字接合がそのまま残ります。**実際に踏みました。**
    PlaneVertexIndex vertex_index;

    // 縫合後の多角形（大域 ID）
    std::vector<std::vector<std::uint32_t>> polys(frags.size());
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        polys[fi].reserve(raw_polys[fi].size());
        for (std::uint32_t v : raw_polys[fi]) polys[fi].push_back(remap[v]);
    }

    // ---- 5. 断片の正準化: 重複割り当てと共平面重複の仕分け ----
    //
    // 同一領域を占める断片をまとめます。
    //   - 同じ owner が複数 → §4.2 の重複割り当て。1 つ残して数える
    //   - 両方の owner が居る → 共平面重複。§4.3.2 の符号 0 の分岐
    std::map<detail::RegionKey, std::array<std::vector<std::size_t>, 2>> groups;
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        groups[detail::region_key(frags[fi].support, polys[fi])][frags[fi].owner].push_back(fi);
    }

    std::vector<std::size_t> reps;        // 代表となる断片の添字
    std::vector<FragClass> forced_class;  // 共平面重複で確定した分類
    std::vector<char> is_coplanar;        // forced が有効か（偽なら手順 6 で決める）
    reps.reserve(groups.size() * 2);
    for (const auto& kv : groups) {
        const auto& by_owner = kv.second;
        const bool both = !by_owner[0].empty() && !by_owner[1].empty();
        FragClass cls = FragClass::Outside;
        if (both) {
            // 外向き法線が同じ向きか。`flipped` は代表平面の法線に対する反転なので、
            // 同じ支持平面どうしなら `flipped` の一致がそのまま向きの一致になります。
            const bool same =
                frags[by_owner[0].front()].flipped == frags[by_owner[1].front()].flipped;
            cls = same ? FragClass::CoplanarSame : FragClass::CoplanarOpposite;
            if (same) {
                ++st.coplanar_same;
            } else {
                ++st.coplanar_opposite;
            }
        }
        for (int o = 0; o < 2; ++o) {
            if (by_owner[o].empty()) continue;
            st.duplicate_fragments += by_owner[o].size() - 1;
            reps.push_back(by_owner[o].front());
            forced_class.push_back(cls);
            is_coplanar.push_back(both ? 1 : 0);
        }
    }
    st.fragments = reps.size();

    // ---- 6. 符号ベクトルによる分類（§4.3.2, §6.1）----
    //
    // 相手メッシュの平面に対する符号ベクトルが等しい断片は、同じ凸領域に属するので
    // 内外が一致します。領域ごとに 1 回だけレイキャストします。
    auto other_planes = [&](int owner) {
        std::vector<PlaneId> r;
        for (const Face& f : (owner == 0) ? faces_b : faces_a) {
            if (std::find(r.begin(), r.end(), f.support) == r.end()) r.push_back(f.support);
        }
        std::sort(r.begin(), r.end());
        return r;
    };
    const std::vector<PlaneId> planes_of[2] = {other_planes(0), other_planes(1)};

    std::map<std::pair<int, std::vector<std::int8_t>>, bool> region_inside;
    std::vector<char> keep(reps.size(), 0);

    for (std::size_t ri = 0; ri < reps.size(); ++ri) {
        const Fragment& f = frags[reps[ri]];
        FragClass cls = forced_class[ri];
        // §3.2: early-out したセルの断片は、セルの隅 1 点の判定をそのまま使います。
        // **符号ベクトルもレイキャストも作りません。**
        if (!is_coplanar[ri] && frag_forced[reps[ri]] >= 0) {
            cls = (frag_forced[reps[ri]] == 1) ? FragClass::Inside : FragClass::Outside;
        } else if (!is_coplanar[ri]) {
            const std::vector<PlaneId>& qs = planes_of[f.owner];
            std::vector<std::int8_t> sig(qs.size());
            for (std::size_t k = 0; k < qs.size(); ++k) {
                sig[k] = static_cast<std::int8_t>(fragment_sign(table, f, qs[k], cache));
            }
            const auto key = std::make_pair(f.owner, sig);
            auto it = region_inside.find(key);
            if (it == region_inside.end()) {
                // 領域の代表点: **断片の相対内部の点**（SPEC-phase3 §2.1 の段 0）。
                //
                // Phase 1 / 2 は「頂点 → 対角線の中点 → 3 頂点の重心」の 3 段でした。
                // 2 段目以降は点を組み合わせる操作で、被符号値が 487 ビットに達します。
                // 段 0 で構成を平面ベースに戻し、**相対内部の点を直接作ります。**
                //
                // 相対内部の点は相手の平面配置のセルの内部にあるので、**$\partial B$ に
                // 載りません**（$\partial B$ は $B$ の平面の和集合に含まれ、セルの内部は
                // どの平面とも交わらない）。したがって境界判定が要りません。
                const mesh::TriMesh& other = (f.owner == 0) ? B : A;
                const geom::HPointD rep = interior_point(table, f, cache, &st.interior);
                KRISITE_CHECK(!point_on_boundary(other, rep),
                              "boolean_op: 相対内部の代表点が相手の境界上にある"
                              "（SPEC-phase3 §2.1 の前提が破れている）");
                // **二項経路は意図的に索引を使いません**（`IMPL-phase5.md` §12）。
                // **正解器なので素朴に保ちます。** 配線漏れではありません。
                const bool inside = point_inside(other, rep);
                ++st.raycasts;
                it = region_inside.emplace(key, inside).first;
                ++st.regions;
            }
            cls = it->second ? FragClass::Inside : FragClass::Outside;
        }
        keep[ri] = detail::select_fragment(op, f.owner, cls) ? 1 : 0;
    }

    // ---- 7. 出力（扇状三角形化）----
    //
    // 頂点順は「外から見て CCW」= 所有メッシュの外向き（§3.4）。`Difference` で
    // 採用する B の断片は、A\B の境界としては法線が逆になるので**順序を反転**します。
    //
    // **T 頂点の解決は選択のあと、三角形化の直前**（§2.4.4 (1)、`IMPL-phase2.md` §2.6.1）。
    // 正準化を先にするのは、同一領域が複数セルに重複割り当てされたとき、隣接の細かさが
    // 違えば入る T 頂点も違い得るためです。先に入れると `region_key` が食い違い、
    // 重複が併合されずに同じ面を二重に出力します。
    //
    // **一律に適用します。** 「隣が持っているから入れる」ではなく「線分の内部に載る
    // 大域頂点はすべて入れる」。片側だけに入れると T 字接合を作ってしまいます。
    BoolMesh out;
    std::vector<int> tri_owner;
    /// **early-out で arrangement を省いたセル由来か**（§13 の CP5 の相互作用）。
    std::vector<char> tri_from_early;
    out.vertices = std::move(merged);
    {
        // 出力に残る断片の支持平面だけ索引を張る（平面数 x 頂点数 回の `side`）
        std::vector<PlaneId> sup;
        for (std::size_t ri = 0; ri < reps.size(); ++ri) {
            if (keep[ri]) sup.push_back(frags[reps[ri]].support);
        }
        std::sort(sup.begin(), sup.end());
        sup.erase(std::unique(sup.begin(), sup.end()), sup.end());
        vertex_index.build(table, out.vertices, sup);
    }
    for (std::size_t ri = 0; ri < reps.size(); ++ri) {
        if (!keep[ri]) continue;
        const Fragment& f = frags[reps[ri]];
        std::vector<std::uint32_t> poly = polys[reps[ri]];
        std::vector<PlaneId> edge = f.edge;
        if (poly.size() < 3) continue;
        if (op == BoolOp::Difference && f.owner == 1) detail::reverse_polygon(poly, edge);
        const TPolygon tp = insert_t_vertices(table, out.vertices, vertex_index, f.support, edge,
                                              poly, &st.t, &merged_from_cache);
        const std::size_t before_n = out.triangles.size();
        fan_triangulate(tp, out.triangles, &st.t);
        // **owner は出力時に直接記録します。** 多角形の頂点数から再構成すると、
        // T 頂点が入ったときに数が合いません（§9.3 の変異 9 で使う）
        tri_owner.insert(tri_owner.end(), out.triangles.size() - before_n, f.owner);
        tri_from_early.insert(tri_from_early.end(), out.triangles.size() - before_n,
                              frag_forced[reps[ri]] >= 0 ? 1 : 0);
    }

    // ---- §5.4: 1 点に集まる平面の最大枚数（総当たり）----
    //
    // **セル面も数えます。** §5.2 が壊れる形は {P, F, Q} と {P, F, R}（F はセル面）
    // なので、セル面を外すと肝心の同時交差を数え落とします。
    st.planes_total = table.size();
    st.mesh_planes = n_mesh_planes;
    for (const geom::HPointD& v : out.vertices) {
        std::size_t cnt = 0, mcnt = 0;
        for (PlaneId q = 0; q < static_cast<PlaneId>(table.size()); ++q) {
            if (geom::side(table.at(q), v) == 0) {
                ++cnt;
                if (q < n_mesh_planes) ++mcnt;
            }
        }
        st.max_planes_at_point = std::max(st.max_planes_at_point, cnt);
        st.max_mesh_planes_at_point = std::max(st.max_mesh_planes_at_point, mcnt);
    }

#if defined(KRISITE_COUNT_PREDICATES)
    st.side_calls = geom::counters::side_calls;
    st.intersect3_calls = geom::counters::intersect3_calls;
#endif
    // ---- 8. 接触の分裂（§5）----
    //
    // **辺（次数 4）と頂点（k 個の扇）の両方を、頂点まわりの扇で一度に扱います。**
    // 次数 4 の辺は扇を繋がないので、辺の両端で扇が分かれ、辺そのものが 2 本になります。
    //
    // **面は増えません。** 境界の点集合は変わらず、組合せ的な表現だけが変わります（§5.1.3）。
#if defined(KRISITE_MUTATION_NO_SPLIT)
    // SPEC-phase2 §9.3 の変異 8: **接触の分裂を無効化**（辺・頂点の両方）。
    // §9.3 の除外 3 件が再び落ちるはずです。
    const bool do_split = false;
    (void)opt;
#else
    const bool do_split = opt.split_contacts;
#endif
    if (do_split && !out.triangles.empty()) {
        std::vector<std::uint32_t> origin;
#if defined(KRISITE_MUTATION_SPLIT_BY_OWNER)
        // SPEC-phase2 §9.3 の変異 9: §5.1.2 の対応付けを owner に戻す。
        // **$\chi$ と体積では捕まりません。** §5.5.1 の $C$ 不変性でのみ落ちます。
        // **ケース 16（自己接触）が無いと素通りします**（§8.1）。
        const std::vector<int>* owner_ptr = &tri_owner;
#else
        const std::vector<int>* owner_ptr = nullptr;
        (void)tri_owner;
#endif
        out.triangles = mesh::split_contacts(out.triangles, out.vertices.size(), &origin, &st.split,
                                             owner_ptr, &tri_from_early);
        for (std::uint32_t o : origin) out.vertices.push_back(out.vertices[o]);
    }

    st.cache_hits = point_cache.hits();
    st.cache_misses = point_cache.misses();
    st.cache_entries = point_cache.entries();
    st.cache_bytes = point_cache.bytes();
    if (stats) *stats = st;
    return out;
}

/// 互換の呼び出し形（Phase 1 からの呼び出し側を変えないため）。
inline BoolMesh boolean_op(const mesh::TriMesh& A, const mesh::TriMesh& B, BoolOp op,
                           unsigned depth, BoolStats* stats = nullptr, bool cull_planes = true) {
    BoolOptions opt;
    opt.depth = depth;
    opt.cull_planes = cull_planes;
    return boolean_op(A, B, op, opt, stats);
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_BOOLEAN_HPP

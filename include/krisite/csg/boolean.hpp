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
#include <tuple>
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
    /// **2 点以上が併合されたグループの数**（＝広がりが 1 以上の点の数）。
    ///
    /// > **★ この計数と `max_merge_span` は、2026-09-05 まで二項メッシュ経路にしか
    /// > 実装されていませんでした。** スープ経路（実際に使う経路）では常に 0 で、
    /// > **番人が空回りしていました**（`DESIGN-phase5-hotspots.md` §9.4）。
    std::size_t merge_groups = 0;
    std::size_t duplicate_fragments = 0;  ///< 重複割り当てが生んだ重複断片（§5.4）
    /// **★ 出力の外接箱が、セル箱より狭まった多角形の数**（`SPEC-phase5.md` §5.10.6）。
    ///
    /// **番人です。** 箱を「元の多角形の箱 $\cap$ セル箱」にする機構が、
    /// **1 個も狭められていないなら空回りしています。**
    /// **「機構を足したら、それが空回りしていないことを別に検査する」**（`CLAUDE.md`）。
    std::size_t out_aabb_narrowed = 0;
    // ---- ★ early-out の到達可能性解析（`SPEC-phase5.md` §5.10.8）の【測定】--------
    //
    // **判定するだけで、まだ捨てていません。** 出力は 1 ビットも変わりません。
    // **「捨てられる葉が 0 なら、案 A を実装する意味がありません。」**
    /// **確定した source が 1 つ以上ある葉の数**（＝ 機会の分母）。
    std::size_t eo_forced_leaves = 0;
    /// **★ 抽象評価が定値になった葉の数**（＝ 捨てられる葉。効果の上限）。
    std::size_t eo_const_leaves = 0;
    /// 同、その葉に届いた **source の三角形数**の和（$\sum P_\ell$）。
    ///
    /// **葉の数だけでは効きが分かりません。** 小さい葉ばかりなら効きません。
    std::size_t eo_const_input = 0;
    /// **★ 同、$\sum P_\ell^2$。** 局所 BSP が $O(P_\ell^2)$ なので**これが本当の効き**です。
    std::size_t eo_const_input_sq = 0;
    /// 同、その葉に**割り当てられた多角形数**の和と、その 2 乗和。
    ///
    /// **`leaf_input_*` は source の三角形、こちらは `polys`。両方が費用に効きます**
    /// （切断集合は三角形から作り、切る対象は多角形）。
    std::size_t eo_const_polys = 0;
    std::size_t eo_const_poly_sq = 0;
    /// **★ 定値の葉で実際に使われた仕事**（代理ではなく、省ける量そのもの）。
    ///
    /// > **`CLAUDE.md`「数えている量が、費用の代理になっているかを確かめてください」。**
    /// > **$\sum P_\ell^2$ は三角形で数えるか多角形で数えるかで大きく食い違いました**
    /// > **（実測 43.4% 対 0.9%）。そこで実際の演算回数を直接数えます。**
    ///
    /// `bsp_cut_slots`（局所 BSP の切断候補の走査）と、
    /// `frag_edges_count`（作った断片の数）の、定値の葉での増分です。
    /// **分母は同名の全体の計数です。**
    std::size_t eo_const_bsp_slots = 0;
    std::size_t eo_const_frags = 0;
    /// **★ 実際に捨てた葉の数**（`early_out_reachability` が真のとき）。
    ///
    /// **番人**: **0 ならバイト一致の検査が何も言っていません**（`SPEC-phase5.md` §5.10.8.5）。
    std::size_t eo_dropped_leaves = 0;
    std::size_t coplanar_same = 0;       ///< 共平面重複のうち向きが同じ対の数
    std::size_t coplanar_opposite = 0;   ///< 向きが逆の対の数
    std::size_t constructed_points = 0;  ///< 第1段が作った構成点の総数（§5.4 の分母）
    std::size_t merged_points = 0;       ///< 第2段の併合後の点数
    std::size_t merged_by_value = 0;     ///< 第1段が取りこぼし第2段が併合した数（§5.4）
    // ---- ★ 実験: 共平面重複の仕分けを、切断の符号列で行えるか --------------------
    //
    // **`RESEARCH-perf.md` §S3.5。`KRISITE_EXPERIMENT_REGION_HIST` で有効になります。**
    // **既定では 0 のままです**（実験の経路が存在しません）。
    /// **突き合わせたグループの数**（頂点 ID による仕分けのグループ数）
    std::size_t region_cmp_groups = 0;
    /// **★ そのうち断片が 2 個以上のグループ**（= 共平面重複が実際にある証拠）。
    ///
    /// **0 なら、両方の鍵が自明に一致しています。番人としてこれを数えます。**
    std::size_t region_cmp_multi = 0;
    /// **★ 2 つの鍵で仕分けが食い違ったグループの数。0 でなければ実験は失敗です。**
    ///
    /// **グループの【数】ではなく【中身】（断片の添字集合）を比べます。**
    /// 数が同じでも、違う断片が同じグループに入っている可能性があります。
    std::size_t region_cmp_mismatch = 0;
    /// **切断の履歴の長さ**（符号列の費用。断片あたりの平均を出すための分子）
    std::size_t region_hist_total = 0;
    std::size_t region_hist_max = 0;
    /// **符号列を数えた断片の数**（`region_hist_total` の分母）。
    /// **群の数で割ってはいけません** — 群と断片は 1 対 1 ではありません
    std::size_t region_hist_count = 0;
    // ---- ビット列版（`DESIGN-phase5-hotspots.md` §14.2）------------------------
    /// **ビット列（符号だけ。平面 ID を落とす）による仕分けが食い違った群の数**。
    ///
    /// **成立するのは「同じセル・同じ支持平面のグループの中」だけです**（§14.2.3）。
    /// **セルをまたぐ比較が起きるなら、ビット列は使えません。**
    std::size_t region_bits_mismatch = 0;
    /// **★ 番人: 同じ `region_key` の断片が、複数のセルに分かれている群の数。**
    ///
    /// **0 なら §14.2.3 の論証は検証されていません**（セルをまたぐ機会が無かっただけ）。
    /// **非零で食い違い 0 なら、論証が実測で裏づけられます。**
    std::size_t region_cross_cell = 0;
    /// **★ セルまたぎの【機構】を切り分けるための計数**（`SPEC-phase5.md` §5.10.5.7）。
    ///
    /// **`region_cross_cell` は「またいだ」ことしか言いません。**
    /// **`.claude/rules/deduction.md` §2.1 が要求する「なぜまたぐか」に答えるには、
    /// またいだ群の【中身】を見る必要があります。**
    ///
    /// 群の中の断片の辺平面の集合が**完全に一致する**群の数。
    /// **一致するなら「同じ多角形が 2 つのセルに割り当てられた」**ということで、
    /// 切断の結果が違ったのではありません。
    std::size_t region_cross_cell_same_edges = 0;
    /// 同上のうち、**支持平面が軸平行**（= セル境界面と同じ向き）の群の数。
    ///
    /// **セル境界平面に乗る多角形は、半開区間の割り当てでは片側だけに入るはず**です。
    /// **にもかかわらず両側に入っているなら、割り当てに使った外接箱が
    /// 多角形のものではありません**（`polysoup.hpp` の `Poly::aabb` は「保守的」）。
    std::size_t region_cross_cell_axis_support = 0;
    /// **符号列が 256 ビットを超えた断片の数**（§14.2.4 の「あふれ」）。
    std::size_t region_bits_overflow = 0;
    // ---- ビット列版が失敗した原因の切り分け（2026-09-08）------------------------
    /// **群の中で `hist` の【長さ】が揃っていない群の数。**
    ///
    /// **揃っていなければ、「$i$ 番目の切断」が群の中で同じ平面を指しません。**
    std::size_t region_bits_len_differ = 0;
    /// **群の中で `hist` の【平面の列】が揃っていない群の数**（長さは同じでも中身が違う）。
    std::size_t region_bits_planes_differ = 0;
    /// **同じ頂点 ID の群が、違うビット群に散った件数**（ビット列が細かすぎる側）。
    std::size_t region_bits_split = 0;
    /// **★ 違う頂点 ID の群が、同じビット群に混ざった件数**（ビット列が粗すぎる側）。
    std::size_t region_bits_merge = 0;
    /// **セルを鍵に加えたビット列**（セル, 支持平面, 符号列）での食い違い。
    std::size_t region_cellbits_mismatch = 0;
    std::size_t max_planes_at_point = 0;       ///< 1 点に集まる平面の最大枚数（§5.4、セル面込み）
    std::size_t max_mesh_planes_at_point = 0;  ///< 同上、メッシュ平面のみ（対照）
    std::size_t planes_total = 0;              ///< 総当たりの分母（表に載った平面の総数）
    std::size_t mesh_planes = 0;               ///< うちメッシュ由来
    std::size_t regions = 0;                   ///< 相異なる符号ベクトルの数（§6.1）
    /// **巻き数が負だった領域の数**（`SPEC-phase5.md` §2.9 の第四段階）。
    ///
    /// **指示関数を $w \ne 0$ から $w > 0$ に訂正したのが効いたかの、直接の証拠です。**
    /// **0 なら、その入力では旧定義と新定義が一致します。**
    ///
    /// 表と裏のどちらかで、いずれかの source の巻き数が負だった領域を数えます。
    std::size_t regions_negative_w = 0;
    /// **巻き数が 2 以上だった領域の数**（同上）。
    ///
    /// **入力の符号付き体積 $V_\text{in} = \int w\,dV$ と出力の体積が食い違うのは、
    /// 巻き数が 0 と 1 以外の値を取る場所があるときだけ**です。その内訳の片方。
    std::size_t regions_w_ge2 = 0;
    std::size_t raycasts = 0;  ///< レイキャスト回数
    /// 代表点の構成（SPEC-phase3 §2.1 の段 0）。**どちらの経路で決まったか。**
    InteriorStats interior{};
    /// `side` と `intersect3` の呼び出し数（SPEC-phase1 §12）。
    ///
    /// **`KRISITE_COUNT_PREDICATES` を定義したビルドでのみ埋まります。**
    /// 既定ビルドでは 0 のままです（計数のコストを本番に持ち込まないため）。
    std::uint64_t side_calls = 0;
    std::uint64_t intersect3_calls = 0;
    /// **★ 段ごとの `side` の内訳**（`KRISITE_COUNT_PREDICATES` のときだけ非零）。
    ///
    /// > **スープ経路には述語の計数がありませんでした**（2026-09-08 に発覚）。
    /// > **二項メッシュ経路（`boolean.hpp` の 1282 行）にしか無く、
    /// > 実際に使う経路では常に 0 でした。**
    /// > **`CLAUDE.md`「計装がどちらにあるかも確かめてください」の 2 度目です。**
    ///
    /// **`side_calls_arrange` は arrange（葉ごと）、
    /// `side_calls_classify` は分類（領域ごと）の合計**です。
    /// **代表点の構成（`interior.side_tests`）は後者の内側にあります。**
    std::uint64_t side_calls_arrange = 0;
    std::uint64_t side_calls_classify = 0;
    std::uint64_t intersect3_arrange = 0;
    std::uint64_t intersect3_classify = 0;
    // ---- ★ `side` の被符号値の【実際の】幅（`SPEC-phase5.md` §5.10.10 の案 E）------
    //
    // **区分は 64 / 128 / 192 / それ以上の 4 つ**（リム数の段が 64 ビットごとなので、
    // **192 ビットで済めば 4 リムから 1 リム減ります**）。
    // **`KRISITE_COUNT_PREDICATES` のときだけ非零です。**
    std::uint64_t side_w64_arrange = 0, side_w128_arrange = 0;
    std::uint64_t side_w192_arrange = 0, side_wmore_arrange = 0;
    std::uint64_t side_w64_classify = 0, side_w128_classify = 0;
    std::uint64_t side_w192_classify = 0, side_wmore_classify = 0;
    /// **観測した最大幅**（上界 `bits::kSide` と並べるために要ります）。
    std::uint64_t side_wmax = 0;
    /// **★ E2 の判定が選ぶリム数**（1 / 2 / 3 / 4）。段は分けず、全体で数えます。
    /// **被符号値の実際の幅（上の 4 区分）と並べると、判定の保守性が分かります。**
    std::uint64_t side_disp1 = 0, side_disp2 = 0, side_disp3 = 0, side_disp4 = 0;
    /// **★ 見積もり用のリム乗算の回数**（いま / E2）。**比が乗算の削減の見積もりです。**
    std::uint64_t side_mul_now = 0, side_mul_e2 = 0;
    /// **判定が読んだリムの数**（前判定の費用）。
    std::uint64_t side_disp_limbreads = 0;

    /// §2.3 の絞り込み（SPEC-phase2）。
    ///
    /// `split_plane_slots` は「セル × 分割平面」の総当たり数、
    /// `split_planes_used` は絞り込んだ後に実際に使った数です。
    /// **比が Phase 2 の断片数削減の直接的な指標になります**（§11）。
    /// **10-a が落とした割り当ての数**（`BoolOptions::exact_assign`）。
    /// **0 なら機構が空回りしています**（`CLAUDE.md`）。
    std::size_t assign_rejected_plane = 0;

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
    /// **$\sum_\ell P_\ell^2$**（`PERF.md` #1′）。局所 BSP の費用が $O(P_\ell^2)$ なら、
    /// **これが arrange の説明変数**になります。**最大や平均では代用できません。**
    ///
    /// 2 通り持ちます — **入力三角形**（葉に届いた source の三角形）と
    /// **多角形**（葉に割り当てられた `polys`）。**どちらが効くかは測って決めます。**
    std::size_t leaf_input_sq = 0;
    std::size_t leaf_poly_sq = 0;
    /// $\sum_\ell$（葉に割り当てられた多角形数）。**`leaf_poly_sq` の分母側**。
    std::size_t leaf_poly_total = 0;
    // ---- 無次元群（`PERF.md` §1.8。**借りてよいのは無次元群が一致するときだけ**）----
    //
    // **$\sum_\ell P_\ell^2$ が同じでも「単位の中身」が違えば、
    // 単位あたりの費用は違って当然です。** 合成と実データで揃っているかを見ます。
    /// 出来た断片の辺の総数（**頂点数 = 辺数**）と、その断片の個数。
    ///
    /// > **分母に `fragments` を使ってはいけません。** あれは
    /// > **`out.polys.size()`（最終出力の多角形数）**で、arrange が作った断片の数とは
    /// > 別の量です。**平均が最大を超えて気づきました**（`IMPL-phase5.md` §38）。
    std::size_t frag_edges_total = 0;
    std::size_t frag_edges_count = 0;
    std::size_t frag_edges_max = 0;
    /// 葉ごとの**相異なる支持平面の数**の総和。÷ `leaf_input_total` で
    /// 「葉あたりの平面数 / 多角形数」
    std::size_t leaf_planes_total = 0;

    // ---- arrange を段に刻む（`PERF.md` §1.1）------------------------------------
    //
    // **合成と実データで単位あたりの費用が 2 倍違う理由**を、当てにいかずに絞ります。
    // **無次元群では説明がつきませんでした**（差の向きが逆）。
    //
    // > **原因を推測するより、測定の粒度を上げるほうが速い**（`CLAUDE.md`）。
    //
    // **葉の粒度で計ります。** 断片ごとに時計を読むと、断片が数百万個ある実データで
    // 計時自体が支配します。**粗いところから刻むこと。**
    //
    // **スレッドをまたいで合算するので、これは CPU 時間**です（壁時計ではありません）。
    // 見るのは**割合**で、絶対値ではありません。
    double ms_arr_gather = 0.0;    ///< 葉に届く多角形の収集（$O(L \cdot P)$ の走査）
    double ms_arr_present = 0.0;   ///< 存在判定と平面ごとの三角形表（$O(L \cdot n)$ + `std::map`）
    double ms_arr_prep = 0.0;      ///< early-out / 平面の絞り込み / 切断集合の準備
    double ms_arr_frag = 0.0;      ///< **断片の生成**（セル平面のクリップ + 局所 BSP + 確定）
    double ms_arr_coplanar = 0.0;  ///< 共平面重複の突き合わせ（EMBER §4.3 の C4）
    double ms_arr_stitch =
        0.0;  ///< **縫合の葉ごとの部分**（鍵の重複除去 + 値の整列 + 境界の印。§5.10.13.3）
    // **断片の生成の内訳**（`measure_frag` のときだけ。全スレッドの CPU の和）
    double ms_fr_clip = 0.0;    ///< セル平面でのクリップ
    double ms_fr_prep = 0.0;    ///< 切断集合の取得（`cuts_for`）と器の準備
    double ms_fr_cut = 0.0;     ///< 切断ループ（O3 の判定 + `split_fragment_into`）
    double ms_fr_commit = 0.0;  ///< 確定（葉の配列へ移す）
    // **分類の内訳**（`measure_classify` のときだけ。CPU の和 = 全スレッドぶん。壁時計は
    // `ms_cl_par_wall`）
    double ms_cl_prep = 0.0;      ///< 領域の代表断片の選択、強制値の確認
    double ms_cl_rep = 0.0;       ///< 代表点の構成（`interior_point`）
    double ms_cl_ray = 0.0;       ///< レイキャスト（`winding_split` × source）
    double ms_cl_out = 0.0;       ///< 指示関数の評価と出力多角形の生成
    double ms_cl_par_wall = 0.0;  ///< 領域ループ（並列）の壁時計。`ms_classify` との差が逐次部分
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
    /// **単一 source の葉を割る規則（§3.1 の訂正）が発火した回数**。
    /// **0 なら `single_src_sq` は空回りしています**（`IMPL-phase5.md` §81）
    std::size_t single_src_splits = 0;

    /// **レイキャストが検査した三角形の総数**（`SPEC-phase5.md` の CP1.5）。
    ///
    /// **★ このコメントは CP1.5 で古くなっていました（2026-09-05 訂正）。**
    ///
    /// 旧: 「`winding_split` は source の全三角形を走査します（枝刈りなし）。
    /// したがって 1 レイあたりの検査数は source の三角形数そのものです」
    ///
    /// **`ray_index`（2 次元階層格子、既定 ON）が入ったので、いまは枝刈り後の数です。**
    /// **旗を OFF にすると全数に戻ります**（外せる形にしてあるので直接測れます）。
    ///
    /// | 模型 | ON | OFF |
    /// |---|---:|---:|
    /// | `425318`（8,352 三角形） | **210.3 / レイ** | 8,352.0 |
    /// | `78535`（10,218 三角形） | **733.7 / レイ** | 10,218.0 |
    ///
    /// **分類の時間は 27 倍 / 13 倍 速くなっています**（`BENCH.md`）。
    ///
    /// **古いコメントを信じて「枝刈りなし」と報告した事例があります。**
    /// **機構を足したら、それを説明しているコメントも直してください。**
    std::size_t ray_tri_tests = 0;
    /// **そのうち実際に寄与した三角形の数**（`DESIGN-phase5-hotspots.md` §11）。
    ///
    /// **`ray_tri_tests` との比が、索引の絞り込みの効きです。**
    /// **A-3 と同じ形の問い**で、A-3 では「索引が返す候補の 99.9% が捨てられている」でした。
    /// **D の前判定を通った候補の数**（`ray_tri_tests` との比が絞り込みの効き）。
    std::size_t ray_tri_kept = 0;
    std::size_t ray_tri_hits = 0;
    /// **投影 AABB の前判定を通った候補の数**（同上）。
    /// **`ray_tri_tests` との比が、安い前判定で落とせる分です。**
    std::size_t ray_tri_aabb = 0;
    /// **投影 AABB を通り、かつレイの前方にある候補の数**（D-1 の後の D-2）。
    std::size_t ray_tri_fwd = 0;
    /// **D-2 単独**（D-1 を掛けない）。
    std::size_t ray_tri_fwd_only = 0;
    /// **安い前判定が走った回数**（早期打ち切りあり）。
    std::size_t ray_cheap_tests = 0;
    /// **索引の段ごとの候補数と項目数**（`record_ray_levels` の下。段 0 が最も細かい）。
    /// 項目数は source × 3 軸の和。
    std::size_t ray_cand_level[12] = {};
    std::size_t ray_items_level[12] = {};
    std::size_t ray_levels_max = 0;
    std::size_t ray_fine_cap_min = 0,
                ray_fine_cap_max = 0;  ///< 使った K（source × 軸の最小 / 最大）
    /// **粒度**（同上。三角形の数 / その段のセル 1 つに収まる数 / 段 0 のセルで数えた項目数）
    std::size_t ray_tri_level[12] = {};
    std::size_t ray_fit1_level[12] = {};
    std::size_t ray_cells0_level[12] = {};
    /// 覆うセル数の区分（≤4 / ≤16 / ≤64 / ≤256 / ≤1024 / >1024）ごとの三角形数と $c_u c_v$ の和
    std::size_t ray_hist_tri[6] = {};
    std::size_t ray_hist_cells[6] = {};

    /// §5.4 の局所 BSP（CP4）。**切断候補のうち何枚を実際に切ったか。**
    ///
    /// `bsp_cut_slots` は「セル内の候補平面 × 支持平面」の総当たり数。
    /// **`bsp_cuts_used` と `bsp_cuts_skipped` の両方が非零であること**を
    /// テストで確かめます。**片方が 0 なら判定が空回りしています**（一方に倒れている）。
    std::size_t bsp_cut_slots = 0;
    /// **局所 BSP の切断のうち、実際に断片を 2 つに分けた数**（`SPEC-phase3.md` §5.4
    /// の未測定の問い。 「返した候補のうち何割使うか」の 6 例目）。`bsp_split_attempts` は（断片 ×
    /// 切断平面）の試行数
    std::size_t bsp_split_attempts = 0;
    std::size_t bsp_split_actual = 0;
    /// **切断平面に一度も分けられなかった断片**（セル平面のクリップだけで分かれた断片。O2
    /// の前提）と、 そのような断片を持つ元の三角形の数（`frags_uncut - frags_uncut_tris`
    /// がレイの削減の上界）
    std::size_t frags_uncut = 0;
    std::size_t frags_uncut_tris = 0;
    std::size_t bsp_cuts_used = 0;     ///< 支持平面に触れるので切った枚数
    std::size_t bsp_cuts_skipped = 0;  ///< 厳密に片側なので切らなかった枚数
    std::size_t bsp_skip_box = 0;      ///< O3: 安い前判定（整数の箱）で飛ばした（断片 × 平面）
    std::size_t bsp_skip_exact = 0;    ///< O3: 厳密な比較で飛ばした（同）
    std::size_t bsp_skip_boxside = 0;  ///< 多角形の箱が平面の片側なので飛ばした（断片 × 平面）

    /// §3.1 の適応分割。葉の深度の分布（**深さの差が §2.4 の前提**）。
    unsigned leaf_depth_min = 0;
    unsigned leaf_depth_max = 0;

    /// §3.2 の early-out。
    std::size_t early_out_cells = 0;      ///< 相手が居ないので arrangement を省いたセル
    std::size_t empty_cells = 0;          ///< どちらも居ないので丸ごと飛ばしたセル
    std::size_t early_out_fragments = 0;  ///< そのセルで作った断片（分類を省いた数）
    std::size_t early_out_raycasts = 0;   ///< セルの隅（**整数点**）1 つで済ませた判定
    /// **early-out が求めた巻き数が負だった回数**（`IMPL-phase5.md` §55.2）。
    ///
    /// 以前は `-1` を「未確定」の印にしていたので、**ここが非零なら誤認が起きていました。**
    /// **番人としても使えます** — 修正後にこれが非零なら、修正が効いた証拠です。
    std::size_t early_out_negative = 0;

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
    double ms_prepare = 0;  ///< 平面表の統合、source ごとの平面・AABB
    double ms_leaves = 0;   ///< 葉の列挙（八分木の構築）
    // **前処理と葉の列挙の内訳**（§5.10.14.33。逐次なので時計を置くだけ）
    double ms_pre_copy = 0;       ///< 多角形の複製と平面表の統合
    double ms_pre_intern = 0;     ///< 三角形の平面の intern（`plane_from_triangle` + 表）
    double ms_pre_rayplanes = 0;  ///< レイ用の平面の複製
    double ms_pre_rayindex = 0;   ///< レイ索引の構築（source × 3 軸）
    double ms_pre_split = 0;      ///< 切断平面の集合（整列・一意化）
    double ms_pre_aabb = 0;       ///< 三角形の箱
    double ms_leaves_count = 0;   ///< 葉の列挙のうち「数える」呼び出し（セルごとに全多角形を走査）
    std::size_t leaves_count_calls = 0;  ///< 数える呼び出しの回数（訪れたセル）
    std::size_t leaves_count_tests = 0;  ///< 同、多角形の箱とセルの判定の回数
    double ms_arrange = 0;               ///< セルごとの arrangement（§4）
    double ms_stitch = 0;                ///< 縫合と重複の仕分け（§5 / §6）
    double ms_classify = 0;              ///< 分類（§7）
    // **縫合の内訳**（壁時計。§5.10.13。逐次なので CPU と同じ）
    double ms_st_points = 0;   ///< 第 1 段: 平面 3 つ組の表 + 構成点の計算
    double ms_st_sort = 0;     ///< 第 2 段: `lex_less` の整列
    double ms_st_remap = 0;    ///< 第 2 段: 値の等しい区分に番号を振る
    double ms_st_regions = 0;  ///< 仕分け: `regions` の map
    // **並列化の前提の確認**（`measure_stitch` のときだけ）
    std::size_t pt_multi_leaf = 0;   ///< 2 つ以上の葉から参照された構成点（第 1 段の鍵で）
    std::size_t pt_mixed_depth = 0;  ///< うち、参照した葉の深度が 2 種類以上
    std::size_t pt_on_boundary =
        0;  ///< 最初に参照した葉の閉じた箱の境界に載る構成点（大域併合の候補の上界）
    std::size_t pt_multi_not_boundary = 0;  ///< 葉をまたぐのに境界に無い点（補題の対偶。0 のはず）
    std::size_t stitch_boundary_classes =
        0;                           ///< 縫合（葉ごと）: 境界に載る類の数（大域の併合の入力）
    std::size_t stitch_classes = 0;  ///< 縫合（葉ごと）: 葉の中の類の総数（境界 + 内部）
    std::size_t leaf_adj_pairs = 0;  ///< 閉じた箱が接する葉の対
    std::size_t leaf_adj_mixed = 0;  ///< うち深度の違う対
    std::size_t leaf_adj_max = 0;    ///< 1 つの葉が接する葉の最大数

    /// §2.4.3 の T 頂点の解決（SPEC-phase2）。
    ///
    /// **★★ この記述は誤りでした（2026-09-05 実測で否定）。**
    ///
    /// 旧: 「**固定深度では `t_inserted` が 0 でなければなりません**（負の対照）。
    /// 大域平面集合で切るので、線分の内部に載る頂点はその平面での切断点として既に
    /// 多角形の角になっているためです」
    ///
    /// **実測（`425318`、self-union、固定深度、絞り込み OFF）: 走査した辺 399,826 に対し
    /// 挿入 31,292（7.83%）。0 ではありません。**
    ///
    /// **「大域平面集合で切る」が成り立っていません。** 切断に使うのは
    /// **そのセルに割り当たった三角形の支持平面**であって、大域の平面集合ではありません。
    ///
    /// **そして深度を上げると T 頂点は【減ります】**（56,442 → 20,484、深度 2 → 8）。
    /// **セル境界が原因ではありません。機構は未特定です**（`DESIGN-phase5-hotspots.md` §6.5）。
    ///
    /// **`ray_index` のコメントに続いて 2 件目の陳腐化です。**
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

/// **断片の正準キー（v2）** — `DESIGN-phase5-hotspots.md` §14。
///
/// **(セル, 支持平面, 切断の符号列, ビット数)。**
/// **頂点 ID を使いません。したがって中核に縫合が要りません。**
///
/// **セルが要る理由**（§14.7.2。実測で確定）:
/// **第 2 段の切断（共平面揃え）はセルごと・支持平面ごとに決まる**ので、
/// **符号だけでは違うセルの領域が混ざります**（「粗すぎ」が全件）。
///
/// **ビット数を持つ理由**: `cutbits` の長さだけでは区別できません
/// （65 ビットと 128 ビットは、どちらも 2 ワード）。
using RegionKey2 = std::tuple<std::uint64_t, PlaneId, std::vector<std::uint64_t>, std::uint32_t>;

/// セルを 64 ビットに詰める（深度 20 まで）。
inline std::uint64_t cell_key(const octree::Cell& c) noexcept {
    return (static_cast<std::uint64_t>(c.depth) << 60) | (static_cast<std::uint64_t>(c.i) << 40) |
           (static_cast<std::uint64_t>(c.j) << 20) | static_cast<std::uint64_t>(c.k);
}

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
    /// **割り当てを厳しくする**（`DESIGN-phase5-hotspots.md` §10 の 10-a）。
    ///
    /// **いまの割り当ては三角形の AABB とセルの重なりだけを見ています。**
    /// 斜めの三角形の AABB は、セルが細かくなるほど過剰になります
    /// （実測: 深度 8 で割り当ての 95.4% が実際には交わらない）。
    ///
    /// **真なら、支持平面がセルの閉領域を横切るかも見ます**（既存の `plane_crosses_box`）。
    /// **平面は無限なので、これは保守的な絞り込みです。**
    /// 平面が横切らなければ三角形も交わらないので、**落ちる断片はありません。**
    ///
    /// > **★ 触るのは「この多角形をどのセルで処理するか」だけです。**
    /// > **「このセルを何で切るか」（切断平面の候補）も
    /// > 「この曲面はこのセルに在るか」（early-out）も変えません。**
    /// > `plane(T)` は無限に延びるので、T が届かないセルにも切断点を生みます。
    /// > **それが Phase 1 の変異 3 です**（`predicates.hpp` の警告）。
    /// > **同じ述語を、役割の違う問いに使っています。混同しないこと。**
    ///
    /// **偽にすると完全に外れます。** 正しさの検査は真偽の両方で同じ出力を要求します。
    bool exact_assign = true;
    /// **葉ごとの切断平面の数を書き出す**（検査用。`nullptr` なら何もしません）。
    ///
    /// **10-a が「切断平面の側」を変えていないことの直接の検査に使います**
    /// （`DESIGN-phase5-hotspots.md` §10.10）。
    ///
    /// **葉の列挙（`build_leaves`）は 10-a で変わらない**ので、
    /// **`leaves` の並びは旗の ON / OFF で同一**です。したがって**添字で比べられます。**
    ///
    /// 呼び出し側が `leaves.size()` を知らないので、中で `resize` します。
    /// **到達しなかった葉は `kNotReached`**（`here` が空で早期に戻った葉）。
    std::vector<std::size_t>* leaf_cull_out = nullptr;
    static constexpr std::size_t kNotReached = static_cast<std::size_t>(-1);
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
    /// **★ 分類のレイキャストに、安い前判定を掛ける**（`DESIGN-phase5-hotspots.md` §11 の D）。
    ///
    /// **索引が返した候補のうち、寄与するのはレイあたり 2.5〜2.8 だけ**（候補は 210〜734）。
    /// **D-2（前方か）→ D-1（投影 AABB）の順で、厳密に絞ります。**
    ///
    /// **偽にすると完全に外れます。** 正しさの検査は真偽の両方で同じ出力を要求します。
    bool ray_prefilter = true;
    /// **段の境界で進捗を出す**（既定 偽。診断用）。
    ///
    /// **対の中で「いまどの段にいるか」が見えないと、長時間かかったときに
    /// 「進んでいるか」を判断できません**（`IMPL-phase5.md` §79）。
    /// **段の時間は完了後にしか読めない**ので、走っている最中には使えません。
    bool verbose_stages = false;
    /// **前判定の効きを計測する**（既定 偽。真にすると計数のぶん遅くなります）。
    bool record_ray_filter = false;
    /// **索引の粒度の計測**（`SPEC-phase5.md` §5.10.14.11。既定 偽）。
    /// 候補が索引のどの段から来たかと、段ごとの項目数を数えます。
    bool record_ray_levels = false;
    /// **★ レイ索引の細かい割り当て（A。§5.10.14.15）**: 段 0 で覆うセル数がこれ以下の三角形は
    /// 覆うセル全部に入れる。**0 で従来どおり**（正解器）。値は記憶の上限から決める（`BENCH.md`）。
    std::size_t ray_index_fine_cells = 0;
    /// **同、記憶の上限で決める形**: 三角形 1 枚あたりの項目数の上限（例 16）。0 で使わない。
    /// 非零なら模型ごと・軸ごとに K を導き、`ray_index_fine_cells` は使わない。
    /// 既定 0（2026-09-10 判断: 上限は絶対量で置く。下の `ray_index_fine_bytes`）。
    std::size_t ray_index_fine_budget = 0;
    /// **★ レイ索引の記憶の上限（絶対量。source × 軸ごとのバイト数）**（§5.10.14.17。既定 16 MB）。
    /// この上限を満たす最大の K を模型・軸ごとに導く。**記憶の制約は絶対量なので、こちらが既定。**
    /// 0 で使わない（`ray_index_fine_budget` / `ray_index_fine_cells` へ）。
    std::size_t ray_index_fine_bytes = std::size_t{16} << 20;
    /// **★ O3: 切る三角形が触れない断片は切らない**（`SPEC-phase5.md` §5.10.14.24。**既定 0 =
    /// 切らない判定をしない**）。
    ///
    /// | 値 | 判定 |
    /// |---|---|
    /// | 0 | 従来（セル内の候補平面すべてで切る） |
    /// | 1 | **安い前判定だけ**: 多角形の箱 ∩
    /// セルの箱（整数）が、切る三角形の箱と交わらなければ飛ばす | | 2 | 1
    /// に加えて、交わるときは断片の頂点と三角形の箱の厳密な比較（Phase 3 の実験と同じ）で飛ばす |
    ///
    /// **既定 2**（2026-09-11。条件 1 は `DESIGN-phase5-hotspots.md` §36.8: 共平面の揃え（C4）が
    /// `region_key` の前に 走るので、切り方の違いは重なりの分割に残らない。突く構成
    /// `tests/csg/test_o3_coplanar.cpp` で発火と一致を確認）。
    int bsp_skip_disjoint = 2;
    /// **★ 多角形の箱が切断平面の片側に完全にあれば、その多角形の全断片でその平面を飛ばす**
    /// （§5.10.14.30 の 7 例目。既定 真）。箱は「多角形の箱 ∩ セルの箱」（整数）で、
    /// 断片はその中にあるので
    /// **早期 return と同じ判定を頂点を評価せずに行うだけ**（出力は不変）。偽で従来（正解器）。
    bool bsp_skip_boxside = true;
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
    /// **単一 source のセルを割る閾値**（$P_\ell^2$ で比べる。`IMPL-phase5.md` §80）。
    ///
    /// **既定 $128^2 = 16384$ は実測から決めました**（`docs/BENCH.md`）。
    /// `934258x111599` を深度 6 / 7 / 8 で、$P \in \{64,\dots,2048\}$ を振った結果、
    /// **3 つの深度すべてで 128 が最小**でした（64 とは 1% 差で並び、256 以降は悪化）。
    /// **体積は 18 設定すべてで厳密に一致**しています。
    ///
    /// `octree::kNoSingleSplit` にすると従来どおり割りません（§9 の比較用）。
    std::size_t single_src_sq = 16384;
    /// **接触の分裂**（§5）。**既定で有効です**（§5.2）。
    ///
    /// 正則化ブールは一般に多様体出力を保証できません。辺だけ・頂点だけを共有する接触が
    /// 出力に残るので、分裂させて多様体化します。Phase 5 のメッシュ化と GWN が多様体を
    /// 前提にできると扱いが楽なため、既定を分裂側にしています。
    ///
    /// **無効にしたときは $g$ で比較してはいけません**（§5.4）。非多様体出力では
    /// $\chi$ が奇数になり得ます。
    /// **共平面重複の仕分けに、切断の符号列を使う**（`DESIGN-phase5-hotspots.md` §14）。
    ///
    /// **鍵は `(セル, 支持平面, 切断の符号列)`。** 真にすると**中核の縫合が要らなくなります**
    /// （実測で中核の 29.2%）。
    ///
    /// > **★ 既定は偽です。連鎖で成立しないことが分かりました**（§14.8）。
    /// > **単発の演算では 3 入力すべてで食い違い 0 でしたが、
    /// > 連鎖（$(A\cup B)\setminus D$）では群がセルをまたぎ、
    /// > 同じ群が 2 つの鍵に割れます**（実測 371 件）。
    ///
    /// **真にすると、従来の鍵との突き合わせも同時に走ります**（検査のため）。
    /// **★ 出力の外接箱を「元の多角形の箱 $\cap$ セル箱」に狭めるか**
    /// （`SPEC-phase5.md` §5.10.6）。**既定は真。**
    ///
    /// > **偽にすると従来どおりセル箱を入れます。**
    /// > **`CLAUDE.md`「正しさの検査では、性能のための機構を無効化できること。
    /// > 『実質的に無効』ではなく『完全に外れる』形にすること」。**
    ///
    /// **実際に要りました。** 箱を狭めると格子が変わり、
    /// **変異 17（存在判定を半開区間で見る）が観測可能になる配置が
    /// コーパスから消えました**（`DESIGN-phase5-hotspots.md` §17.5）。
    /// **唯一の検出器だったので、外す経路が無ければ網が縮みます。**
    /// **★ 到達可能性による葉の除去**（`SPEC-phase5.md` §5.10.8。EMBER §4.5.2）。**既定は真。**
    ///
    /// **指示関数を 3 値論理で抽象評価し、
    /// 「このセルに居る source」に依存しないなら葉ごと捨てます。**
    ///
    /// **その葉のどの断片も `in_front == in_back` になるので、出力は変わりません。**
    /// **`early_out` が偽なら `forced` が無いので、この機構も自動的に無効になります。**
    ///
    /// > **偽にすると完全に外れます**（`CLAUDE.md`「正しさの検査では、
    /// > 性能のための機構を無効化できること」）。**比較の正解器側です。**
    bool early_out_reachability = true;
    /// **★ 構成点キャッシュを `std::map` に戻す**（`SPEC-phase5.md` §5.10.12）。**既定は偽。**
    ///
    /// **既定は開番地法のハッシュ表です。** 真にすると従来の `std::map` に戻ります。
    ///
    /// > **正解器として残しています**（`CLAUDE.md`「従来の経路を旗で残し、
    /// > 両者が一致することを検査してください」）。
    /// > **同一実行の中で A/B を取れるので、時間の比較が時間帯の交絡を受けません。**
    bool point_cache_map = false;
    /// **★ 断片の切断を従来の形（`SplitResult` を値で返す）に戻す**（§5.10.12.4）。**既定は偽。**
    ///
    /// **既定は `split_fragment_into`（切らないなら移動、確保 0 回）です。**
    /// **従来の経路を正解器として残し、旗の ON / OFF でバイト一致を検査します。**
    bool split_legacy = false;
    /// **★ 器の使い回し**（`SPEC-phase5.md` §5.10.12.4 の G1〜G3）。**ビットの集合。既定は
    /// 7（全部）。**
    ///
    /// | ビット | 何を | 何が消えるか |
    /// |---|---|---|
    /// | 1（G1） | 切断ループの `pieces` / `next` を葉の外で使い回す | 多角形 × 切断ごとの確保 |
    /// | 2（G2） | 分類の `w_front` / `w_back` と指示関数の作業配列をその場の器に | 領域ごとの確保
    /// | | 4（G3） | 縫合の `raw[fi]` を平坦な配列に | 断片ごとの確保 |
    ///
    /// **0 で従来どおり。** 確保の回数は項目ごとに、時間は合計で 1 回測ります（仕様側の条件）。
    unsigned alloc_reuse = 7;
    /// **★ 縫合の並列化に入る前の確認のための計測**（`SPEC-phase5.md` §5.10.13.3。既定 偽）。
    /// 真にすると、構成点ごとに「参照した葉が 2 つ以上か」「その葉の深度が混ざるか」を数え、
    /// 葉の閉じた箱が接する対を全部数えます（$O(L^2)$ なので既定では切ります）。
    bool measure_stitch = false;
    /// **★ 縫合を葉ごとに始める**（`SPEC-phase5.md` §5.10.13.3。既定 真）。
    ///
    /// 葉の中で平面 3 つ組の重複を除き、値で整列・併合し（arrange の並列部分）、
    /// **葉の閉じた箱の境界に載る類だけ**を大域で値で併合します。番号（札）は
    /// `(葉の番号 << 32) | 葉の中の類の番号` で、境界で併合された類は最小の札を取ります。
    /// **札は葉の列挙順と葉の中の整列だけで決まるので、スレッド数に依りません。**
    /// **偽で従来（大域の `std::map` と整列。逐次）。** 従来とは領域の順序が変わるので、
    /// 出力の一致は「順序を除いた鍵」で検査します（`chain_frag` の縫合の A/B）。
    bool stitch_parallel = true;
    /// **分類の中を刻む計測**（§5.10.14.7。既定 偽）。領域ごとに時計を 4 回読みます。
    bool measure_classify = false;
    /// **断片の生成の中を刻む計測**（§5.10.14.30。既定 偽）。多角形ごとに時計を 4 回読みます。
    bool measure_frag = false;
    bool tight_out_aabb = true;
    bool region_key_cuts = false;
    /// **仕分けは従来の鍵で行いつつ、切断の符号列とも突き合わせる**（**検査だけ**）。
    ///
    /// **出力は変わりません。** `region_key_cuts` と違い、仕分けには使いません。
    /// **コーパスでは常時真にしてください**（§14.3.2。置き換えの根拠を守る唯一の検査）。
    bool verify_region_key = false;
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
    /// **領域の巻き数を出力に残す**（`IMPL-phase5.md` §52）。**既定は偽。**
    ///
    /// **診断のためだけの旗**です。多角形の数だけ確保するので、
    /// **切っておけば費用はゼロ**です。
    bool record_winding = false;
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
    PointCache point_cache(opt.point_cache_map);
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
    const std::vector<octree::Cell> leaves = octree::build_leaves(
        policy, [&](const octree::Cell& c, std::size_t* na, std::size_t* nb, bool* bsp_skipped) {
            // 二項経路には NSI による局所 BSP の省略がありません（soup 経路のみ）。
            *bsp_skipped = false;
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
                            for (Fragment& p : pieces) {
                                if (opt.split_legacy) {
                                    if (q == p.support) {
                                        next.push_back(p);
                                        continue;
                                    }
                                    const SplitResult r = split_fragment(table, p, q, cache);
                                    if (r.has_pos) next.push_back(r.pos);
                                    if (r.has_neg) next.push_back(r.neg);
                                } else {
                                    split_fragment_into(table, std::move(p), q, cache, next);
                                }
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
        // **二項経路は正解器なので、検算を常に有効にします**（§5.5）。
        //
        // **スープ経路（`to_mesh`）は既定で切ります** — 純粋な診断を本番の経路に
        // 置かないため（`HANDOVER.md` §5.1）。**こちらはコーパスでしか使わないので、
        // 費用より検出力を取ります。** `CLAUDE.md`「正解器は被検体と別経路で書く」。
        mesh::SplitOptions sopt;
        sopt.verify_delta = true;
        out.triangles = mesh::split_contacts(out.triangles, out.vertices.size(), &origin, &st.split,
                                             owner_ptr, &tri_from_early, nullptr, sopt);
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

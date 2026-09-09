// Krisite — CP5 全域（$n$ 項 / スープ経路）。`SPEC-phase3.md` §14 の CP5
//
// **CP1〜CP4 で入った機構をすべて有効にして全域を回します。**
//
//     全ケース × 全演算 × 固定深度 0〜3 + 適応分割 × 分裂 on/off
//
// 個々の機構は CP1〜CP4 で確かめました。**CP5 が新しく見るのは機構どうしの相互作用**で、
// スープ経路ではどのチェックポイントでも単独では通らなかった経路が 4 つあります。
//
//   局所 BSP × 適応分割      葉ごとに切断集合が変わる（`test_soup` は固定深度のみ）
//   局所 BSP × early-out     曲面の無い source を隅で決めたセルでも切る経路
//   WNV × 接触の分裂         巻き数で選んだ断片が分裂に回る経路
//   代表点の段 0 / 段 1      速い経路と一般的な経路の両方
//
// **「全部有効にしたら通った」だけでは足りません。** 経路を**実際に踏んだこと**を数え、
// 踏んでいなければ落とします（`CLAUDE.md`「機構を足したら、その機構が空回りして
// いないことを別に検査してください」）。
//
// ---
//
// ## `test_soup.cpp` との分担
//
// | | `test_soup` | ここ |
// |---|---|---|
// | 構成 | `phase1_options`（**すべての機構を切った基準**） | **すべて有効** |
// | 分割 | 固定深度 0〜3 | 固定深度 + **適応分割** |
// | 分裂 | OFF 固定 | **on / off の両方** |
// | 目的 | 機構ごとの正解器比較 | **相互作用** |
//
// **基準側を明示するのは向こうの役目です。** ここは出荷時の構成そのものを回します。
//
// ## 分裂 on / off をまたぐ検査は「分裂で変わってはいけない量」に限る
//
// 分裂は頂点 ID を付け替えるだけで幾何を動かしません（`SPEC-phase2.md` §5.1.3 の定理）。
// したがって**頂点の値で見た三角形の多重集合**は一致しなければなりません。
// 一方 $\chi$ や頂点数は変わるので、またいで比べてはいけません。
#include <algorithm>
#include <array>
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"

#include "corpus.hpp"
#include "corpus_expect.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kMaxDepth = 3;
constexpr std::size_t kPasses = kMaxDepth + 2;  ///< 固定深度 0〜3 + 適応分割

/// **CP1〜CP4 の機構をすべて有効にした構成。**
///
/// **既定値に依存しません**（§9.4 の CI ジョブで既定が反転するため）。
BoolOptions all_on(unsigned depth, bool adaptive) {
    BoolOptions o;
    o.depth = depth;
    o.cull_planes = true;   // Phase 2 §2.3
    o.adaptive = adaptive;  // Phase 2 §3.1
    o.leaf_threshold = 0;
    o.early_out = true;     // Phase 2 §3.2
    o.cache_points = true;  // Phase 2 §4
    o.local_bsp = true;     // Phase 3 §5.4（CP4）
    o.reverse_regions = false;
    o.split_contacts = true;  // スープ経路では `to_mesh` 側が持つ
    return o;
}

const char* op_name(BoolOp op) {
    switch (op) {
        case BoolOp::Union:
            return "∪";
        case BoolOp::Intersection:
            return "∩";
        default:
            return "\\";
    }
}

/// 頂点の**値**で見た三角形の多重集合。
///
/// 値が等しい頂点に同じ番号を振り直し、三角形は**巡回順を保ったまま**最小の番号が
/// 先頭に来るよう回します（向きを落とさないため、並べ替えはしません）。
std::vector<std::array<std::uint32_t, 3>> geometric_key(const SoupMesh& m) {
    const std::size_t n = m.vertices.size();
    std::vector<std::uint32_t> order(n);
    for (std::uint32_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
        return krisite::geom::lex_less(m.vertices[a], m.vertices[b]);
    });
    std::vector<std::uint32_t> canon(n, 0);
    std::uint32_t next = 0;
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        while (j < order.size() &&
               krisite::geom::h_equal(m.vertices[order[i]], m.vertices[order[j]])) {
            canon[order[j]] = next;
            ++j;
        }
        ++next;
        i = j;
    }
    std::vector<std::array<std::uint32_t, 3>> out;
    out.reserve(m.triangles.size());
    for (const krisite::mesh::Tri& t : m.triangles) {
        const std::uint32_t c[3] = {canon[t[0]], canon[t[1]], canon[t[2]]};
        int s = 0;
        if (c[1] < c[s]) s = 1;
        if (c[2] < c[s]) s = 2;
        out.push_back({c[s], c[(s + 1) % 3], c[(s + 2) % 3]});
    }
    std::sort(out.begin(), out.end());
    return out;
}

struct Totals {
    std::size_t configs = 0;
    std::size_t order_checks = 0;
    std::size_t early_out_checks = 0;
    std::size_t nary_configs = 0;
    // ---- 相互作用の計数（**空回りの番人**）----
    std::size_t bsp_x_adaptive = 0;   ///< 適応分割で葉に深さの差があり、かつ BSP が切った
    std::size_t bsp_x_early_out = 0;  ///< early-out が発火し、かつ BSP が切った
    std::size_t wnv_x_split = 0;      ///< 分裂が実際に起きた構成
    std::size_t uneven_leaves = 0;    ///< 葉の深さに差が出た構成
    // ---- 機構が発火した総数 ----
    std::size_t bsp_cuts_used = 0;
    std::size_t bsp_cuts_skipped = 0;
    std::size_t early_out_cells = 0;
    std::size_t cache_hits = 0;
    std::size_t interior_axis = 0;    ///< 代表点の段 0
    std::size_t interior_corner = 0;  ///< 代表点の段 1
    std::size_t split_vertices = 0;
    std::size_t unresolved = 0;
    // ---- radial sort（`SPEC-phase2.md` §5.1.2.1）----
    std::size_t radial_attempted = 0;  ///< 連結成分では分けられなかった辺
    std::size_t radial_resolved = 0;   ///< そのうち角度順で分けられた辺
    std::size_t unresolved_off = 0;    ///< **radial sort を外したときの `unresolved`**
    // ---- 検証の増分計算（`IMPL-v2.md` §2）----
    /// **増分計算と従来経路（正解器）を突き合わせた構成の数。** 0 なら空回りです
    std::size_t verify_agree = 0;
    /// **事後の検査（分裂の【後】に非多様体）が発火した回数。**
    /// **`unresolved` のうち、対応付け不能な辺ではないほう**（`SPEC-phase5.md` §1.5.0.1）
    std::size_t unresolved_post = 0;
    std::size_t unsplit_edges = 0;  ///< 対応付けできず、分裂させずに残した辺
    // ---- T 解決の走査順（`SPEC-phase5.md` §5.11）----
    std::size_t scan_per_edge_total = 0;  ///< 辺ごとに走査したときの走査回数
    std::size_t scan_per_poly_total = 0;  ///< 多角形あたり 1 回にしたときの走査回数
    std::size_t side_per_edge_total = 0;  ///< 同、`side` の評価回数
    std::size_t side_per_poly_total = 0;
    // ---- 候補の区間の絞り込み（案 (b2)。`DESIGN-phase5-hotspots.md` §13）----
    std::size_t range_skipped = 0;        ///< 二分探索で飛ばした候補（**0 なら空回り**）
    std::size_t side_unsorted_total = 0;  ///< 絞り込みなしの `side`
    std::size_t side_sorted_total = 0;    ///< 絞り込みありの `side`
    // ---- 共平面重複の仕分けの鍵（`DESIGN-phase5-hotspots.md` §14）--------------
    std::size_t rk_groups = 0;    ///< 突き合わせた群
    std::size_t rk_multi = 0;     ///< **断片が 2 個以上の群**（重なりがある証拠）
    std::size_t rk_mismatch = 0;  ///< **食い違い。0 でなければ置き換えは誤り**
    std::size_t rk_cross = 0;     ///< **群がセルをまたいだ件数。0 でなければ鍵が割れる**
    std::size_t rk_overflow = 0;  ///< 符号列が 256 ビットを超えた断片
    std::size_t rk_bits_total = 0, rk_bits_max = 0, rk_bits_count = 0;
    // **連鎖の側は別に数えます**（§14.8。**連鎖では成立しません**）
    std::size_t rkn_groups = 0, rkn_mismatch = 0, rkn_cross = 0;
    /// **★ 出力の外接箱が狭まった多角形の数**（`SPEC-phase5.md` §5.10.6 の番人）。
    /// **0 なら機構が空回りしています。**
    std::size_t aabb_narrowed = 0;
    /// **★ セルまたぎの【機構】の内訳**（`DESIGN-phase5-hotspots.md` §15.3）。
    /// **「またいだ」は 1 段目です。中身を見ないと機構は決まりません。**
    std::size_t rkn_cross_same = 0, rkn_cross_axis = 0;
    /// **箱を狭める機構を【外した】側のセルまたぎ**（§17.3 の論証の対偶）。
    std::size_t rkn_cross_loose = 0;
};

Totals g;

/// 1 つの構成を回す。
void run_config(const kritest::Case& c, const TriMesh& a, const TriMesh& b, BoolOp op,
                const BoolOptions& opt, const std::string& tag) {
    const PolySoup sa = from_mesh(a), sb = from_mesh(b);

    BoolStats st{};
    // **★ コーパスでは、従来の鍵（頂点 ID）との突き合わせを常時走らせます**
    // （`DESIGN-phase5-hotspots.md` §14.3.2。**置き換えが正しいことを示す唯一の検査**）
    BoolOptions opt_v = opt;
    opt_v.verify_region_key = true;  // **仕分けは従来どおり。突き合わせだけ**（§14）
    const PolySoup r = boolean(sa, sb, op, opt_v, &st);
    g.aabb_narrowed += st.out_aabb_narrowed;
    g.rk_groups += st.region_cmp_groups;
    g.rk_multi += st.region_cmp_multi;
    g.rk_mismatch += st.region_cmp_mismatch;
    g.rk_cross += st.region_cross_cell;
    g.rk_overflow += st.region_bits_overflow;
    g.rk_bits_total += st.region_hist_total;
    g.rk_bits_max = std::max(g.rk_bits_max, st.region_hist_max);
    g.rk_bits_count += st.region_hist_count;

    // ---- 分裂 off / on -------------------------------------------------------
    ToMeshOptions off, on;
    off.split_contacts = false;
    on.split_contacts = true;
    // **§5.5 の検算は既定で切ってあります**（純粋な診断。`HANDOVER.md` §5.1）。
    // **検査するテストは明示的に立ててください。**
    on.verify_split_delta = true;
    ToMeshStats t_off{}, t_on{};
    const SoupMesh m_off = to_mesh(r, off, &t_off);
    const SoupMesh m_on = to_mesh(r, on, &t_on);

    const TopologyReport r_off = check_topology(m_off.triangles);
    const TopologyReport r_on = check_topology(m_on.triangles);

    // 分裂 ON は多様体でなければなりません（§5.1.2 / Phase 2 の判断 3）。
    //
    // **例外は §5.1.2.2 の「対応付けできなかった辺」だけです**（案 A）。
    // 判定は識別子ではなく `unresolved > 0` という機械的な事実で行います。
    const kritest::Exclusion ex_on = kritest::exclusion_when_split(t_on.split.unresolved, r_on);
    if (ex_on != kritest::Exclusion::None) {
        std::string why;
        KRI_CHECK_MSG(kritest::exclusion_conditions_ok(ex_on, r_on, &why),
                      tag + ": §9.3.1 の適用条件を満たさない（" + why + "）");
        KRI_CHECK_MSG(r_on.edges_deficient == 0, tag + ": **辺が裂けました**（次数 1 が " +
                                                     std::to_string(r_on.edges_deficient) +
                                                     " 本）");
    } else {
        KRI_CHECK_MSG(r_on.empty || r_on.ok(), tag + ": 分裂 ON で多様体になっていない");
    }
    // 分裂 OFF は接触辺が残るので `ok()` は使えません。向きは常に整合すべきです
    KRI_CHECK_MSG(r_off.empty || r_off.oriented, tag + ": 分裂 OFF で向きが整合していない");

    // **分裂で幾何は動きません**（§5.1.3 の定理）。値で見た三角形の多重集合が一致すること
    KRI_CHECK_MSG(geometric_key(m_off) == geometric_key(m_on),
                  tag + ": **分裂で幾何が動いた**（頂点 ID の付け替えのはず）");

    // 分裂の予測と実測（§9.4.2）
    KRI_CHECK_MSG(t_on.split.predicted_delta_v == t_on.split.actual_delta_v,
                  tag + ": ΔV の予測と実測が違う" +
                      kritest::pair_msg(t_on.split.predicted_delta_v, t_on.split.actual_delta_v));
    KRI_CHECK_MSG(t_on.split.predicted_delta_e == t_on.split.actual_delta_e,
                  tag + ": ΔE の予測と実測が違う" +
                      kritest::pair_msg(t_on.split.predicted_delta_e, t_on.split.actual_delta_e));

    // ---- ★ 検証の増分計算が、従来経路（正解器）と一致すること ----------------
    //
    // **`detail::split_topology` は `check_topology` を 2 回呼ぶ従来の経路を
    // 置き換えたものです**（`IMPL-v2.md` §2）。**検査を弱めていないことを、
    // ここで直接確かめます。**
    //
    // **正解器は被検体と別経路です** — 従来経路は大域の `std::map` で辺と
    // 頂点リンクを作る、まったく別の実装です（`CLAUDE.md`）。
    //
    // **`unresolved` は CP の失敗判定そのものなので**（`SPEC-phase5.md` §3.-1）、
    // **一致しなければ、その差は失敗件数の差として実データに出ます。**
    {
        ToMeshOptions on_naive = on;
        on_naive.verify_split_naive = true;
        ToMeshStats t_naive{};
        const SoupMesh m_naive = to_mesh(r, on_naive, &t_naive);
        KRI_CHECK_MSG(t_naive.split.unresolved == t_on.split.unresolved,
                      tag + ": unresolved が増分計算と従来経路で違う" +
                          kritest::pair_msg(t_naive.split.unresolved, t_on.split.unresolved));
        KRI_CHECK_MSG(
            t_naive.split.actual_delta_v == t_on.split.actual_delta_v,
            tag + ": ΔV が増分計算と従来経路で違う" +
                kritest::pair_msg(t_naive.split.actual_delta_v, t_on.split.actual_delta_v));
        KRI_CHECK_MSG(
            t_naive.split.actual_delta_e == t_on.split.actual_delta_e,
            tag + ": ΔE が増分計算と従来経路で違う" +
                kritest::pair_msg(t_naive.split.actual_delta_e, t_on.split.actual_delta_e));
        // **出力は 1 ビットも変わってはいけません**（検証は出力に触りません）
        KRI_CHECK_MSG(m_naive.triangles == m_on.triangles,
                      tag + ": 検証の経路を変えたら出力が変わった");
        // **前提（退化三角形が無い）が破れていないこと。** 一般解が入っているので 0
        KRI_CHECK_MSG(t_on.split.verify_fallback == 0,
                      tag + ": 増分計算が前提の破れで正解器に落ちた（退化三角形がある）");
        g.unresolved_post += t_on.split.unresolved_post;
        g.unsplit_edges += t_on.split.unsplit_edges;
        ++g.verify_agree;
    }

    // ---- ★ T 解決の走査順を変えても、出力がバイト一致すること ------------------
    //
    // **候補集合は (葉, 支持平面) 群で共有されるので、多角形の中で同じです。**
    // **辺ごとに走査し直す必要がありません**（`SPEC-phase5.md` §5.11）。
    //
    // **判定は 1 つも変えていません。走査の順序だけです。**
    // **だからバイト一致で守れます** — 絞り込みではないので、
    // 「落ちる T 頂点がないこと」を別に示す必要がありません。
    {
        ToMeshOptions per_edge = on;
        per_edge.scan_per_edge = true;
        ToMeshStats t_pe{};
        const SoupMesh m_pe = to_mesh(r, per_edge, &t_pe);
        KRI_CHECK_MSG(m_pe.triangles == m_on.triangles,
                      tag + ": **走査順を変えたら三角形が変わった**（判定は同じはず）");
        KRI_CHECK_MSG(m_pe.vertices.size() == m_on.vertices.size(),
                      tag + ": 走査順を変えたら頂点数が変わった");
        KRI_CHECK_MSG(t_pe.t.inserted == t_on.t.inserted,
                      tag + ": 走査順で挿入した T 頂点の数が変わった" +
                          kritest::pair_msg(t_pe.t.inserted, t_on.t.inserted));
        // **走査の回数は減っていなければなりません**（多角形あたり 1 回）。
        // **等しいなら機構が空回りしています**（`CLAUDE.md`）
        g.scan_per_edge_total += t_pe.t.cand_scans;
        g.scan_per_poly_total += t_on.t.cand_scans;
        g.side_per_edge_total += t_pe.t.side_tests;
        g.side_per_poly_total += t_on.t.side_tests;
    }

    // ---- ★ 候補の区間の絞り込み（案 (b2)）が出力を変えないこと ----------------
    //
    // **厳密な絞り込みです**（`DESIGN-phase5-hotspots.md` §13）。
    // **辺は多角形の境界上にあるので、辺の区間は多角形の区間に含まれます。**
    // **だから多角形の X 区間の外にある候補は、どの辺の相対内部にも載りません。**
    //
    // **`CLAUDE.md`「前判定の効果は『◯◯の数が 1 個も変わらない』で守れます」。**
    // **挿入した T 頂点の数が 1 個も変わらないことが、厳密であることの直接の検査です。**
    {
        ToMeshOptions unsorted = on;
        unsorted.sort_candidates = false;
        ToMeshStats t_us{};
        const SoupMesh m_us = to_mesh(r, unsorted, &t_us);
        KRI_CHECK_MSG(m_us.triangles == m_on.triangles,
                      tag + ": **候補の絞り込みで三角形が変わった**（厳密なはず）");
        KRI_CHECK_MSG(t_us.t.inserted == t_on.t.inserted,
                      tag + ": **絞り込みで T 頂点の数が変わった**（1 個も変わってはいけない）" +
                          kritest::pair_msg(t_us.t.inserted, t_on.t.inserted));
        KRI_CHECK_MSG(t_us.t.degenerate_kept == t_on.t.degenerate_kept,
                      tag + ": 絞り込みで残した退化三角形の枚数が変わった");
        g.range_skipped += t_on.t.cand_skipped_by_range;
        g.side_unsorted_total += t_us.t.side_tests;
        g.side_sorted_total += t_on.t.side_tests;
    }

    // ---- 二項正解器との一致（分裂 OFF どうし。§10.1）--------------------------
    BoolOptions bo = opt;
    bo.split_contacts = false;
    const BoolMesh ref = boolean_op(a, b, op, bo);
    const TopologyReport r_ref = check_topology(ref.triangles);
    KRI_CHECK_MSG(
        r_ref.components == r_off.components,
        tag + ": C が二項正解器と違う" + kritest::pair_msg(r_ref.components, r_off.components));
    KRI_CHECK_MSG(r_ref.chi == r_off.chi,
                  tag + ": χ が二項正解器と違う" + kritest::pair_msg(r_ref.chi, r_off.chi));

    // ---- early-out 有無で結果が変わらないこと ★ -----------------------------
    //
    // **`test_early_out.cpp` は二項経路しか見ていません。** スープ経路、とくに
    // source が 3 つ以上ある構成では、「この source の曲面がこのセルに存在するか」の
    // 判定が半開区間だと**面がセルの上側境界に乗る配置で壊れます**
    // （`IMPL-phase3.md` §7.1。実際に踏みました）。
    //
    // **これが変異 17 の主検出器です。**
    BoolOptions noeo = opt;
    noeo.early_out = false;
    const SoupMesh m_noeo = to_mesh(boolean(sa, sb, op, noeo), off);
    KRI_CHECK_MSG(geometric_key(m_off) == geometric_key(m_noeo),
                  tag + ": **early-out の有無で結果が変わった**（省いたことで分類が壊れている）");
    ++g.early_out_checks;

    // ---- 順序非依存（全機構を有効にした状態で。§14 の CP3 の判定）-------------
    //
    // **`test_soup` は機構を全部切った構成でしか見ていません。** メモ化と early-out を
    // 有効にした状態でも順序に依らないことは、ここでしか確かめられません。
    BoolOptions rev = opt;
    rev.reverse_regions = true;
    const SoupMesh m_rev = to_mesh(boolean(sa, sb, op, rev), off);
    KRI_CHECK_MSG(geometric_key(m_off) == geometric_key(m_rev),
                  tag + ": **領域を回す順序で結果が変わった**（並列化の前提）");
    ++g.order_checks;

    // ---- 計数（空回りの番人）--------------------------------------------------
    const bool uneven = st.leaf_depth_max > st.leaf_depth_min;
    if (uneven) ++g.uneven_leaves;
    if (uneven && st.bsp_cuts_used > 0) ++g.bsp_x_adaptive;
    if (st.early_out_cells > 0 && st.bsp_cuts_used > 0) ++g.bsp_x_early_out;
    if (t_on.split.split_vertices > 0) ++g.wnv_x_split;
    g.bsp_cuts_used += st.bsp_cuts_used;
    g.bsp_cuts_skipped += st.bsp_cuts_skipped;
    g.early_out_cells += st.early_out_cells;
    g.cache_hits += st.cache_hits;
    g.interior_axis += st.interior.axis_line;
    g.interior_corner += st.interior.corner_offset;
    g.split_vertices += t_on.split.split_vertices;
    g.unresolved += t_on.split.unresolved;
    g.radial_attempted += t_on.split.radial_attempted;
    g.radial_resolved += t_on.split.radial_resolved;
    // **radial sort を外した基準側**（§5.1.2.1 の配置に到達していることの番人）。
    // **機構を足したら、それを外す経路も用意する** — 外した側で到達を数えます。
    {
        ToMeshOptions no_radial = on;
        no_radial.radial_sort = false;
        ToMeshStats t_nr{};
        (void)to_mesh(r, no_radial, &t_nr);
        g.unresolved_off += t_nr.split.unresolved;
    }
    ++g.configs;
    (void)c;
}

/// $n$ 項（3 source）を全機構有効で回す（§10.2）。
///
/// **中間結果をメッシュに戻さないこと**が要点なので、二項の連鎖とは比べられません
/// （二項は `BoolMesh` を返し、それを入力に戻せない）。ここでは
/// **分割戦略を変えても $(C, \chi)$ が変わらないこと**を見ます。
void run_nary(const kritest::Case& c, const TriMesh& a, const TriMesh& b, const TriMesh& d) {
    const PolySoup sa = from_mesh(a), sb = from_mesh(b), sd = from_mesh(d);
    ToMeshOptions tm;
    tm.split_contacts = true;

    long long chi0 = 0;
    std::size_t comp0 = 0;
    // **★ 外接箱を狭める機構を、外した側でも回します**（`SPEC-phase5.md` §5.10.6）。
    //
    // > **`CLAUDE.md`「機構を追加したら、それを外す経路も用意してください。
    // > 外せないと、検査がその機構の影に入ります」。**
    //
    // **実際に影に入りました。** 箱を狭めると格子が変わり、
    // **変異 17（存在判定を半開区間で見る）が観測可能になる配置が消えました。**
    // **cp5n が唯一の検出器**なので、外した側を回さないと網が縮みます。
    for (std::size_t pass = 0; pass < kPasses * 2; ++pass) {
        const bool tight = (pass < kPasses);
        const std::size_t p2 = pass % kPasses;
        const bool adaptive = (p2 == kPasses - 1);
        const unsigned depth = adaptive ? kMaxDepth : static_cast<unsigned>(p2);
        BoolOptions o = all_on(depth, adaptive);
        o.tight_out_aabb = tight;
        // **★ 連鎖でも鍵の突き合わせを走らせます**（`DESIGN` §14）
        // **★ 連鎖でも突き合わせます**（§14.8 の食い違いは、ここで測ったものです）
        o.verify_region_key = true;
        BoolStats st1{}, st2{};
        const PolySoup s1 = boolean(sa, sb, BoolOp::Union, o, &st1);
        const PolySoup s2 = boolean(s1, sd, BoolOp::Difference, o, &st2);
        // **★ 連鎖の側は別に数えます**（§14.8）。**単発の結果と混ぜないこと**
        // **★ 鍵の突き合わせは、狭めた側でだけ数えます。**
        // **外した側は「昔の格子」なので、混ぜると §17.3 の論証が検査できません。**
        if (tight) {
            g.rkn_groups += st1.region_cmp_groups + st2.region_cmp_groups;
            g.rkn_mismatch += st1.region_cmp_mismatch + st2.region_cmp_mismatch;
            g.rkn_cross += st1.region_cross_cell + st2.region_cross_cell;
        } else {
            // **外した側では、またぎが復活するはずです**（§17.3 の論証の対偶）。
            g.rkn_cross_loose += st1.region_cross_cell + st2.region_cross_cell;
        }
        g.rkn_cross_same += st1.region_cross_cell_same_edges + st2.region_cross_cell_same_edges;
        g.rkn_cross_axis += st1.region_cross_cell_axis_support + st2.region_cross_cell_axis_support;
        const std::string tag = std::string("ケース ") + c.id + " (A∪B)\\D（" +
                                (adaptive ? "適応" : "深度 " + std::to_string(depth)) +
                                (tight ? "" : "・箱を狭めない") + "）";
        KRI_CHECK_MSG(s2.source_count() == 3, tag + ": source 数が 3 でない");
        const TopologyReport rep = check_topology(to_mesh(s2, tm).triangles);
        KRI_CHECK_MSG(rep.empty || rep.ok(), tag + ": n 項の出力が多様体でない");
        if (pass == 0) {
            chi0 = rep.chi;
            comp0 = rep.components;
        } else {
            KRI_CHECK_MSG(rep.chi == chi0,
                          tag + ": χ が分割戦略で変わった" + kritest::pair_msg(chi0, rep.chi));
            KRI_CHECK_MSG(rep.components == comp0, tag + ": C が分割戦略で変わった" +
                                                       kritest::pair_msg(comp0, rep.components));
        }
        ++g.nary_configs;
    }
}

void run_case(const kritest::Case& c) {
    const TriMesh a = c.make_a(), b = c.make_b();
    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        for (std::size_t pass = 0; pass < kPasses; ++pass) {
            const bool adaptive = (pass == kPasses - 1);
            const unsigned depth = adaptive ? kMaxDepth : static_cast<unsigned>(pass);
            const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) + "（" +
                                    (adaptive ? "適応" : "深度 " + std::to_string(depth)) + "）";
            run_config(c, a, b, op, all_on(depth, adaptive), tag);
        }
    }
    run_nary(c, a, b, kritest::corpus()[1].make_b());
}

/// **空回りの番人。** 期待値は「実測した数」ではなく**式**で持たせます。
void check_not_vacuous() {
    const std::size_t cases = kritest::corpus().size();
    const std::size_t want_configs = cases * 3 * kPasses;
    KRI_CHECK_MSG(g.configs == want_configs,
                  "構成の数が式と合わない" + kritest::pair_msg(want_configs, g.configs));
    KRI_CHECK_MSG(g.order_checks == want_configs, "順序非依存の検査数が構成数と合わない");
    KRI_CHECK_MSG(g.early_out_checks == want_configs, "early-out 比較の数が構成数と合わない");
    KRI_CHECK_MSG(g.nary_configs == cases * kPasses * 2, "n 項の構成数が式と合わない");
    // **★ 外した側では、セルまたぎが復活しなければなりません**（§17.3 の論証の対偶）。
    // **復活しないなら、狭めたことが効いた証拠になっていません。**
    KRI_CHECK_MSG(g.rkn_cross_loose > 0,
                  "**箱を狭めない側でもセルまたぎが 0 でした。**"
                  "§17.3 の論証は「箱が緩いから両側に入る」なので、"
                  "外した側では復活するはずです");

    // 機構が実際に発火したこと
    // **★ 外接箱を狭める機構の番人**（`SPEC-phase5.md` §5.10.6）。
    // **狭まった多角形が 1 個も無いなら、交差を取る意味がありません。**
    KRI_CHECK_MSG(g.aabb_narrowed > 0,
                  "出力の外接箱が 1 個も狭まっていない（`Poly::aabb` の交差が空回り）");
    KRI_CHECK_MSG(g.bsp_cuts_used > 0, "局所 BSP が 1 枚も切っていない");
    KRI_CHECK_MSG(g.bsp_cuts_skipped > 0, "局所 BSP が 1 枚も省いていない");
    KRI_CHECK_MSG(g.early_out_cells > 0, "early-out が 1 度も発火していない");
    KRI_CHECK_MSG(g.cache_hits > 0, "構成点の保持が 1 度も当たっていない");
    KRI_CHECK_MSG(g.split_vertices > 0, "接触の分裂が 1 度も起きていない");
    // **代表点は両経路とも踏むこと**（フォールバック連鎖の空回りの番人。CP1 と同じ）
    KRI_CHECK_MSG(g.interior_axis > 0, "代表点の段 0（軸平行）が 1 度も使われていない");
    KRI_CHECK_MSG(g.interior_corner > 0, "代表点の段 1（角のオフセット）が 1 度も使われていない");

    // 相互作用を**実際に踏んだこと**
    KRI_CHECK_MSG(g.uneven_leaves > 0, "適応分割で葉の深さに差が出た構成が無い");
    KRI_CHECK_MSG(g.bsp_x_adaptive > 0, "局所 BSP × 適応分割 を踏んだ構成が無い");
    KRI_CHECK_MSG(g.bsp_x_early_out > 0, "局所 BSP × early-out を踏んだ構成が無い");
    KRI_CHECK_MSG(g.wnv_x_split > 0, "WNV × 接触の分裂 を踏んだ構成が無い");

    // **0 でなければ報告すること**（§5.1.2.1。radial sort の必要性の判断材料）
    // **§5.1.2.2: 到達すること自体は欠陥ではありません。**
    // 分裂させずに次数 4 のまま残すので、閉じたまま非多様体になります。
    // **件数は記録**で、判定は上の `exclusion_when_split` + 適用条件が行います。
    //
    // **空回り防止**: ケース 24 を入れたので、到達 0 なら検査が効いていません
    // **★ radial sort を入れたので、この配置は解けるようになりました**（§5.1.2.1）。
    // **番人は「到達したか」を、外した側（`unresolved_off`）で数えます。**
    // **入れた側で 0 になることが、機構が効いていることの検査です。**
    KRI_CHECK_MSG(g.unresolved_off > 0,
                  "**§5.1.2.1 の配置に一度も到達していません。** ケース 24 が"
                  "コーパスから消えたか、対応付けの判定が変わっています");
    KRI_CHECK_MSG(g.radial_attempted > 0, "**radial sort が 1 度も呼ばれていません。空回りです**");
    KRI_CHECK_MSG(g.radial_resolved > 0, "**radial sort が 1 本も解けていません**");
    KRI_CHECK_MSG(g.unresolved == 0,
                  "**radial sort を入れたのに解けていない辺が残っています。**"
                  "コーパスに新しい形が入ったか、実装に穴があります");
    // **検証の増分計算の番人**（`IMPL-v2.md` §2）。
    // **突き合わせを 1 度も回していなければ、一致は何も言っていません。**
    KRI_CHECK_MSG(g.verify_agree > 0,
                  "**検証の増分計算と従来経路の突き合わせを 1 度も回していません。空回りです**");
    // ---- 共平面重複の仕分けの鍵（`DESIGN-phase5-hotspots.md` §14）--------------
    //
    // **中核から頂点 ID による仕分け（縫合）を外しました。**
    // **鍵は (セル, 支持平面, 切断の符号列) です。**
    // **置き換えが正しいことを示すのは、この突き合わせだけです。**
    // **単発の演算では成立します**（3 入力で実測。§14.7）
#if defined(KRISITE_FRAGMENT_CUTBITS)
    // **★ 符号列の検査は `KRISITE_FRAGMENT_CUTBITS` のビルドでだけ回ります**（§5.10.5）。
    //
    // > **仕様側の判断で「採りません」と決まった実験なので、
    // > 既定の経路から外しました**（2026-09-09。F1）。
    // > **「連鎖では成立しない」ことの固定は、この別ビルドが受け持ちます。**
    // > **CI に専用のジョブがあります**（`.github/workflows/ci.yml`）。
    KRI_CHECK_MSG(g.rk_mismatch == 0,
                  "**単発の演算で、新しい鍵（セル, 支持平面, 符号列）が"
                  "頂点 ID による仕分けと食い違いました**" +
                      kritest::pair_msg(g.rk_mismatch, 0));
    // **★ 群がセルをまたぐと、新しい鍵では 2 つに割れます。**
    KRI_CHECK_MSG(g.rk_cross == 0,
                  "**単発の演算で群がセルをまたぎました。**"
                  "新しい鍵では同じ群が 2 つに割れます" +
                      kritest::pair_msg(g.rk_cross, 0));
    // ---- ★ 連鎖では成立しないことを、検査として固定します（§14.8）--------------
    //
    // **`CLAUDE.md`「素通りする組合せも固定する」。**
    // **「連鎖でも通った」ことにして進むと、実データで静かに壊れます。**
    // **★ セルまたぎの番人は、証明とセットで置き換えました**（2026-09-08）。
    //
    // **`CLAUDE.md`「番人を retire するときは、証明とセットにしてください」。**
    // **以前は「連鎖で 371 件出る」ことを固定していました。**
    // **`Poly::aabb` を「元の多角形の箱 ∩ セル箱」に狭めたら 0 件になりました**
    // （`DESIGN-phase5-hotspots.md` §15.5 / §17）。
    //
    // **構造的に到達不能であることの論証**:
    //
    //   1. 371 件すべてが「支持平面が軸平行」かつ「辺平面の集合が同一」だった（実測）
    //   2. 多角形が軸平行平面 $x = c$ に乗るなら、その真の外接箱は $x$ で潰れている
    //   3. 箱は「元の箱 ∩ セル箱」で、`from_mesh` の箱は tight なので、
    //      **帰納的に、乗っている多角形の箱は必ず潰れている**
    //   4. 潰れた箱には半開区間 $[\mathrm{lo}, \mathrm{hi})$ の割り当てが片側だけを選ぶ
    //   5. 支持平面が軸平行でない多角形は、セル面と線か点でしか接しないので、
    //      隣の葉でクリップすると頂点数が 3 未満になり、仕分けから外れる
    //
    // **したがって「またぎが消えたこと」を固定します。**
    // **箱を再び緩めると、この検査が落ちます。**
    KRI_CHECK_MSG(g.rkn_cross == 0,
                  "**連鎖で群がセルをまたぎました。** `Poly::aabb` が緩んでいませんか"
                  "（§15.5 の論証。箱は「元の箱 ∩ セル箱」であるべきです）" +
                      kritest::pair_msg(g.rkn_cross, 0));
    // **連鎖で成立しないことは変わりません。** ただし**理由が変わりました**。
    //
    //   以前の理解  セルまたぎ（→ 箱を狭めて消えました）
    //   **いまの理解  適応分割の格子が段ごとに違う**（固定深度なら食い違い 0。§15.4）
    KRI_CHECK_MSG(g.rkn_mismatch > 0,
                  "**連鎖で食い違いが消えました。** §15.4 の制約が変わったなら、"
                  "新しい鍵を連鎖でも使えるか再検討してください");
    // **番人**: **共平面重複が 1 件も無ければ、両方の鍵は自明に一致します。**
    KRI_CHECK_MSG(g.rk_multi > 0,
                  "**断片が 2 個以上の群が 1 つもありません。**"
                  "共平面重複の無い入力だけで比較しており、一致は何も言っていません");
#endif
#if defined(KRISITE_FRAGMENT_CUTBITS)
    KRI_CHECK_MSG(g.rk_groups > 0, "**突き合わせを 1 度も回していません。空回りです**");
#endif
    std::printf(
        "    ★ 仕分けの鍵: 群 %zu / **重なりのある群 %zu** / **食い違い %zu** / "
        "セルまたぎ %zu\n",
        g.rk_groups, g.rk_multi, g.rk_mismatch, g.rk_cross);
    std::printf(
        "       **連鎖では成立しません**（適応の格子が段で違うため。§15.4）: "
        "群 %zu / 食い違い %zu / セルまたぎ %zu\n",
        g.rkn_groups, g.rkn_mismatch, g.rkn_cross);
    std::printf(
        "    外接箱を狭めた多角形 %zu（`Poly::aabb` = 元の箱 ∩ セル箱）"
        "／**狭めない側のセルまたぎ %zu**\n",
        g.aabb_narrowed, g.rkn_cross_loose);
    std::printf("       セルまたぎの内訳: 辺平面が同一 %zu / 支持平面が軸平行 %zu\n",
                g.rkn_cross_same, g.rkn_cross_axis);
    std::printf("       符号列: 断片あたり 平均 %.1f / 最大 %zu / 256 ビット超 %zu\n",
                g.rk_bits_count > 0 ? static_cast<double>(g.rk_bits_total) / g.rk_bits_count : 0.0,
                g.rk_bits_max, g.rk_overflow);
    // **走査順の入れ替えが実際に走査を減らしていること**（空回りの番人）。
    // **等しければ、多角形がすべて 1 辺しかないか、機構が効いていません。**
    // **絞り込みが実際に候補を飛ばしていること**（空回りの番人）。
    // **0 なら、整列が効いていないか旗が通っていません。**
    KRI_CHECK_MSG(g.range_skipped > 0,
                  "**候補の区間の絞り込みが 1 個も飛ばしていません。空回りです**");
    KRI_CHECK_MSG(g.side_sorted_total < g.side_unsorted_total,
                  "**絞り込みで `side` の評価が減っていません**" +
                      kritest::pair_msg(g.side_sorted_total, g.side_unsorted_total));
    KRI_CHECK_MSG(g.scan_per_poly_total < g.scan_per_edge_total,
                  "**T 解決の走査回数が減っていません。** 多角形あたり 1 回の走査が"
                  "効いていないか、旗が通っていません" +
                      kritest::pair_msg(g.scan_per_poly_total, g.scan_per_edge_total));

    std::printf("    構成 %zu（順序非依存 %zu / early-out 比較 %zu / n 項 %zu）\n", g.configs,
                g.order_checks, g.early_out_checks, g.nary_configs);
    std::printf("    局所 BSP: 切った %zu / 省いた %zu、early-out セル %zu、キャッシュ命中 %zu\n",
                g.bsp_cuts_used, g.bsp_cuts_skipped, g.early_out_cells, g.cache_hits);
    std::printf("    代表点: 段 0 %zu / 段 1 %zu、分裂した頂点 %zu\n", g.interior_axis,
                g.interior_corner, g.split_vertices);
    std::printf("    radial sort: 到達 %zu（外した側 %zu）/ 試行 %zu / 解決 %zu\n",
                g.unresolved_off, g.unresolved_off, g.radial_attempted, g.radial_resolved);
    std::printf("    相互作用: BSP×適応 %zu、BSP×early-out %zu、WNV×分裂 %zu（葉に差 %zu）\n",
                g.bsp_x_adaptive, g.bsp_x_early_out, g.wnv_x_split, g.uneven_leaves);
    // **`unresolved` の内訳**（`SPEC-phase5.md` §1.5.0.1。「解けなかった数」だけでは
    // 機構が働いたか分かりません）。**事後の検査が 0 なら、その経路は未検査です**
    std::printf(
        "    検証の増分計算: 突き合わせ %zu 構成、事後の非多様体 %zu、"
        "分裂させず残した辺 %zu\n",
        g.verify_agree, g.unresolved_post, g.unsplit_edges);
    // **走査順の効果は演算回数で出します**（`CLAUDE.md`「効果は演算回数で測る」）。
    // **`side` はほとんど減りません** — 減るのは候補集合の走査だけです（§5.11）
    std::printf("    候補の絞り込み: 飛ばした %zu、`side` %zu → %zu（**%.2f 倍**）\n",
                g.range_skipped, g.side_unsorted_total, g.side_sorted_total,
                g.side_sorted_total > 0
                    ? static_cast<double>(g.side_unsorted_total) / g.side_sorted_total
                    : 0.0);
    std::printf(
        "    T 解決の走査: 辺ごと %zu → 多角形あたり %zu（%.2f 倍）、"
        "`side` %zu → %zu（%.3f 倍）\n",
        g.scan_per_edge_total, g.scan_per_poly_total,
        g.scan_per_poly_total > 0
            ? static_cast<double>(g.scan_per_edge_total) / g.scan_per_poly_total
            : 0.0,
        g.side_per_edge_total, g.side_per_poly_total,
        g.side_per_poly_total > 0
            ? static_cast<double>(g.side_per_edge_total) / g.side_per_poly_total
            : 0.0);
}

}  // namespace

int main() {
    std::printf("\n  CP5 全域（n 項 / スープ経路）— 機構どうしの相互作用（SPEC-phase3 §14）\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) run_case(c);
    check_not_vacuous();
    std::printf("\n");
    return kritest::finish("csg/cp5n");
}

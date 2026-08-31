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
};

Totals g;

/// 1 つの構成を回す。
void run_config(const kritest::Case& c, const TriMesh& a, const TriMesh& b, BoolOp op,
                const BoolOptions& opt, const std::string& tag) {
    const PolySoup sa = from_mesh(a), sb = from_mesh(b);

    BoolStats st{};
    const PolySoup r = boolean(sa, sb, op, opt, &st);

    // ---- 分裂 off / on -------------------------------------------------------
    ToMeshOptions off, on;
    off.split_contacts = false;
    on.split_contacts = true;
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
    for (std::size_t pass = 0; pass < kPasses; ++pass) {
        const bool adaptive = (pass == kPasses - 1);
        const unsigned depth = adaptive ? kMaxDepth : static_cast<unsigned>(pass);
        const BoolOptions o = all_on(depth, adaptive);
        const PolySoup s1 = boolean(sa, sb, BoolOp::Union, o);
        const PolySoup s2 = boolean(s1, sd, BoolOp::Difference, o);
        const std::string tag = std::string("ケース ") + c.id + " (A∪B)\\D（" +
                                (adaptive ? "適応" : "深度 " + std::to_string(depth)) + "）";
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
    KRI_CHECK_MSG(g.nary_configs == cases * kPasses, "n 項の構成数が式と合わない");

    // 機構が実際に発火したこと
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
    KRI_CHECK_MSG(g.unresolved > 0,
                  "**§5.1.2.1 の配置に一度も到達していません。** ケース 24 が"
                  "コーパスから消えたか、対応付けの判定が変わっています");

    std::printf("    構成 %zu（順序非依存 %zu / early-out 比較 %zu / n 項 %zu）\n", g.configs,
                g.order_checks, g.early_out_checks, g.nary_configs);
    std::printf("    局所 BSP: 切った %zu / 省いた %zu、early-out セル %zu、キャッシュ命中 %zu\n",
                g.bsp_cuts_used, g.bsp_cuts_skipped, g.early_out_cells, g.cache_hits);
    std::printf("    代表点: 段 0 %zu / 段 1 %zu、分裂した頂点 %zu\n", g.interior_axis,
                g.interior_corner, g.split_vertices);
    std::printf("    相互作用: BSP×適応 %zu、BSP×early-out %zu、WNV×分裂 %zu（葉に差 %zu）\n",
                g.bsp_x_adaptive, g.bsp_x_early_out, g.wnv_x_split, g.uneven_leaves);
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

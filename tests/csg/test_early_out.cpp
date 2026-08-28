// Krisite — early-out（SPEC-phase2.md §3.2、§9.1）
//
// > $C$ が $B$ の三角形を 1 つも含まないなら、$C$ 全体が $B$ の内側か外側かの
// > **どちらか一方**です。セルの隅 1 つで判定して、arrangement を計算しない。
//
// **early-out は「仕事を省く」機構なので、省いたことで答えが変わっていないことを
// 直接検査します。** 無効側が正解器です（§0.1 と同じ構図）。
//
// ---
//
// **失う検出器を先に書き出しておきます**（`CLAUDE.md` の恒久ルール）。
//
// | 失うもの | 内容 |
// |---|---|
// | 分類の検査 | early-out
// したセルの断片は符号ベクトルもレイキャストも通りません。**分類経路のバグがそのセルでは露出しません**
// | | arrangement の検査 | 相手の平面で切らないので、切断の誤りがそのセルでは出ません |
//
// **残る検出器**: このテスト（early-out 有無の比較）、§10.3 の体積、§9.3 の変異 6。
// **このテストが主検出器です。**
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kMaxDepth = 3;

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

struct Totals {
    std::size_t leaves = 0;
    std::size_t early_cells = 0, empty_cells = 0;
    std::size_t frags_off = 0, frags_on = 0;
    std::size_t early_frags = 0;
    std::size_t raycasts_off = 0, raycasts_on = 0;
    std::size_t early_raycasts = 0;
    std::size_t regions_off = 0, regions_on = 0;
};

Totals g;

void run_case(const kritest::Case& c) {
    const TriMesh a = c.make_a(), b = c.make_b();
    std::printf("\n  ケース %-4s %s\n", c.id, c.what);

    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        // **既定値に依存しないこと**（§9.4 の CI ジョブで既定が反転する）
        BoolOptions off = kritest::phase1_options(kMaxDepth);
        off.adaptive = true;
        BoolOptions on = off;
        on.early_out = true;

        BoolStats s_off, s_on;
        const BoolMesh r_off = boolean_op(a, b, op, off, &s_off);
        const BoolMesh r_on = boolean_op(a, b, op, on, &s_on);
        const TopologyReport t_off = check_topology(r_off.triangles);
        const TopologyReport t_on = check_topology(r_on.triangles);

        const std::string tag =
            std::string("ケース ") + c.id + " " + op_name(op) + "（early-out 有無）";

        // ---- §9.1: 省いても答えは変わらない ----
        KRI_CHECK_MSG(t_on.empty == t_off.empty, tag + ": 空かどうかが変わった");
        KRI_CHECK_MSG(t_on.components == t_off.components,
                      tag + ": C が変わった（" + std::to_string(t_off.components) + " → " +
                          std::to_string(t_on.components) + "）");
        KRI_CHECK_MSG(t_on.chi == t_off.chi, tag + ": χ が変わった（" + std::to_string(t_off.chi) +
                                                 " → " + std::to_string(t_on.chi) + "）");
        KRI_CHECK_MSG(t_on.chi_even == t_off.chi_even, tag + ": χ の偶奇が変わった");
        KRI_CHECK_MSG(t_on.edge_manifold == t_off.edge_manifold,
                      tag + ": 辺多様体の判定が変わった（**T 字接合の疑い**）");
        KRI_CHECK_MSG(t_on.vertex_manifold == t_off.vertex_manifold,
                      tag + ": 頂点多様体が変わった");
        KRI_CHECK_MSG(t_on.oriented == t_off.oriented, tag + ": 向きの整合が変わった");
        KRI_CHECK_MSG(t_on.edges_odd_degree == t_off.edges_odd_degree,
                      tag + ": 次数が奇数の辺の数が変わった（" +
                          std::to_string(t_off.edges_odd_degree) + " → " +
                          std::to_string(t_on.edges_odd_degree) + "）");
        KRI_CHECK_MSG(t_on.no_degenerate == t_off.no_degenerate,
                      tag + ": 退化三角形の有無が変わった");

        g.leaves += s_on.total_cells;
        g.early_cells += s_on.early_out_cells;
        g.empty_cells += s_on.empty_cells;
        g.frags_off += s_off.fragments;
        g.frags_on += s_on.fragments;
        g.early_frags += s_on.early_out_fragments;
        g.raycasts_off += s_off.raycasts;
        g.raycasts_on += s_on.raycasts;
        g.early_raycasts += s_on.early_out_raycasts;
        g.regions_off += s_off.regions;
        g.regions_on += s_on.regions;

        if (op == BoolOp::Union) {
            std::printf(
                "    %-2s 葉 %4zu（省 %3zu 空 %3zu） 断片 %4zu→%-4zu "
                "レイ %3zu→%-3zu(隅%2zu) | C=%zu χ=%lld\n",
                op_name(op), s_on.total_cells, s_on.early_out_cells, s_on.empty_cells,
                s_off.fragments, s_on.fragments, s_off.raycasts, s_on.raycasts,
                s_on.early_out_raycasts, t_on.components, t_on.chi);
        }
    }
}

/// **空回りの番人。** early-out が一度も発火しないなら、この検査は何も検証していません。
void check_not_vacuous() {
    const double cell_pct = g.leaves ? 100.0 * static_cast<double>(g.early_cells + g.empty_cells) /
                                           static_cast<double>(g.leaves)
                                     : 0.0;
    const double frag_pct =
        g.frags_on ? 100.0 * static_cast<double>(g.early_frags) / static_cast<double>(g.frags_on)
                   : 0.0;
    std::printf("\n  §3.3 early-out の記録（全ケース × 3 演算 × 最大深度 %u の合計）\n", kMaxDepth);
    std::printf("    葉 %zu のうち 省略 %zu / 空 %zu = %.1f%%（セル数ベース）\n", g.leaves,
                g.early_cells, g.empty_cells, cell_pct);
    std::printf("    断片 %zu → %zu（うち分類を省いた断片 %zu = %.1f%%）\n", g.frags_off,
                g.frags_on, g.early_frags, frag_pct);
    std::printf("    レイキャスト %zu → %zu（うち隅の整数点 %zu）／領域 %zu → %zu\n",
                g.raycasts_off, g.raycasts_on, g.early_raycasts, g.regions_off, g.regions_on);

    KRI_CHECK_MSG(g.early_cells > 0, "early-out が一度も発火していない。**空回りです**");
    KRI_CHECK_MSG(g.early_frags > 0, "分類を省いた断片が 0。**空回りです**");
    // **省いたぶん仕事が減っていること。** 減らないなら §3.2 の主張が成り立ちません
    KRI_CHECK_MSG(g.frags_on < g.frags_off,
                  "断片が減っていない。early-out は arrangement を省くはずです");
    KRI_CHECK_MSG(g.regions_on < g.regions_off,
                  "領域（相異なる符号ベクトル）が減っていない。分類を省けていません");
}

}  // namespace

int main() {
    std::printf("\n  early-out — SPEC-phase2 §3.2 / §3.3 / CP2\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) run_case(c);
    check_not_vacuous();
    std::printf("\n");
    return kritest::finish("csg/early_out");
}

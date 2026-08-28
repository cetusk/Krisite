// Krisite — 適応分割（SPEC-phase2.md §3、§9.1）
//
// > **同一ケース・同一演算について、固定深度と適応分割の出力で、体積と位相 $(C, \chi)$ が
// > 厳密に一致すること。**
//
// **Phase 1 の実装（固定深度）がそのまま正解器です**（§0.1）。三角形分割は変わりますが
// 立体は同じなので、$(V,E,F)$ ではなく **$(C, \chi)$ で比較**します。
// 体積は GMP が要るので `csg/test_volume_gmp.cpp` が受け持ちます。
//
// **T 字接合はここに出ます。** 継ぎ目が割れれば位相が変わるためです。
//
// ---
//
// **空回りしないことを併せて検査します**（§9.0 と同じ規律）。適応分割が何も変えなければ、
// この検査は「常に同じ結果」を確かめているだけになります。
//
//   - 葉の数が固定深度より減っていること（**節約が実際に出ているか**。§13 の CP2 判定）
//   - **葉の深さに差が出ていること**（§2.4 の前提。差が無いなら T 字接合が生じない）
//   - **T 頂点が実際に入っていること**（案 D が効くべき場所で効いているか）
//   - **ケース 13 の 1 本の辺に T 頂点が 2 個以上**（§8 の要件。**変異 11 の前提**）
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "corpus_expect.hpp"
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
    std::size_t leaves_fixed = 0, leaves_adaptive = 0;
    std::size_t frags_fixed = 0, frags_adaptive = 0;
    std::size_t t_inserted = 0, t_dropped = 0;
    std::size_t max_depth_gap = 0;
    std::size_t max_t_per_edge = 0;
    std::size_t cases_with_gap = 0;
};

Totals g;

/// ケース 13 の要件（§8）を測るための保持。
std::size_t g_case13_max_t_per_edge = 0;
std::size_t g_case13_depth_gap = 0;

void run_case(const kritest::Case& c) {
    const TriMesh a = c.make_a(), b = c.make_b();
    std::printf("\n  ケース %-4s %s\n", c.id, c.what);
    bool gap_seen = false;

    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        // **既定値に依存しないこと**（§9.4 の CI ジョブで既定が反転する）
        BoolOptions fixed_opt = kritest::phase1_options(kMaxDepth);
        BoolOptions adap_opt = fixed_opt;
        adap_opt.adaptive = true;

        BoolStats sf, sa;
        const BoolMesh rf = boolean_op(a, b, op, fixed_opt, &sf);
        const BoolMesh ra = boolean_op(a, b, op, adap_opt, &sa);
        const TopologyReport tf = check_topology(rf.triangles);
        const TopologyReport ta = check_topology(ra.triangles);

        const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) +
                                "（適応分割 vs 固定深度 " + std::to_string(kMaxDepth) + "）";

        // ---- §9.1 が要求する形: (C, χ) の一致 ----
        KRI_CHECK_MSG(ta.empty == tf.empty, tag + ": 空かどうかが変わった");
        KRI_CHECK_MSG(ta.components == tf.components, tag + ": C が変わった（" +
                                                          std::to_string(tf.components) + " → " +
                                                          std::to_string(ta.components) + "）");
        KRI_CHECK_MSG(ta.chi == tf.chi, tag + ": χ が変わった（" + std::to_string(tf.chi) + " → " +
                                            std::to_string(ta.chi) + "）");
        KRI_CHECK_MSG(ta.chi_even == tf.chi_even, tag + ": χ の偶奇が変わった");

        // ---- §9.2: Phase 1 の検査体系が適応分割モードでも通ること ----
        //
        // **絶対値ではなく固定深度と比べます。** §9.3 の除外に該当する構成
        // （4T の ∪ など）は固定深度でも落ちるので、そこは「同じように落ちる」が正解です。
        KRI_CHECK_MSG(ta.edge_manifold == tf.edge_manifold,
                      tag + ": 辺多様体の判定が変わった（**T 字接合の疑い**）");
        KRI_CHECK_MSG(ta.vertex_manifold == tf.vertex_manifold, tag + ": 頂点多様体が変わった");
        KRI_CHECK_MSG(ta.oriented == tf.oriented, tag + ": 向きの整合が変わった");
        // **空メッシュでは `no_degenerate` が立ちません**（`check_topology` は空を
        // V=E=F=0 で見る）。他の旗と同じく固定深度と比べます
        KRI_CHECK_MSG(ta.no_degenerate == tf.no_degenerate,
                      tag + ": 退化三角形（同一頂点を 2 度使う面）の有無が変わった");
        KRI_CHECK_MSG(ta.edges_odd_degree == tf.edges_odd_degree,
                      tag + ": 次数が奇数の辺の数が変わった（" +
                          std::to_string(tf.edges_odd_degree) + " → " +
                          std::to_string(ta.edges_odd_degree) + "）。**面の過不足の直接的な証拠**");

        const unsigned gap = sa.leaf_depth_max - sa.leaf_depth_min;
        if (gap >= 2) gap_seen = true;
        g.leaves_fixed += sf.total_cells;
        g.leaves_adaptive += sa.total_cells;
        g.frags_fixed += sf.fragments;
        g.frags_adaptive += sa.fragments;
        g.t_inserted += sa.t.inserted;
        g.t_dropped += sa.t.degenerate_kept;
        g.max_depth_gap = std::max<std::size_t>(g.max_depth_gap, gap);
        g.max_t_per_edge = std::max(g.max_t_per_edge, sa.t.max_per_edge);
        if (std::string(c.id) == "13") {
            g_case13_max_t_per_edge = std::max(g_case13_max_t_per_edge, sa.t.max_per_edge);
            g_case13_depth_gap = std::max<std::size_t>(g_case13_depth_gap, gap);
        }

        std::printf(
            "    %-2s 葉 %4zu→%-4zu 断片 %4zu→%-4zu 深さ %u..%u T頂点 %zu(最大%zu/辺) "
            "退化%zu | C=%zu χ=%lld\n",
            op_name(op), sf.total_cells, sa.total_cells, sf.fragments, sa.fragments,
            sa.leaf_depth_min, sa.leaf_depth_max, sa.t.inserted, sa.t.max_per_edge,
            sa.t.degenerate_kept, ta.components, ta.chi);
    }
    if (gap_seen) ++g.cases_with_gap;
}

/// **空回りの番人。** 適応分割が何も変えていないなら §9.1 は何も検証していません。
void check_not_vacuous() {
    std::printf("\n  §13 CP2 の判定（全ケース × 3 演算 × 深度 %u の合計）\n", kMaxDepth);
    std::printf("    葉        %zu → %zu\n", g.leaves_fixed, g.leaves_adaptive);
    std::printf("    断片      %zu → %zu\n", g.frags_fixed, g.frags_adaptive);
    std::printf("    T 頂点    %zu 個（1 辺の最大 %zu 個）／残した退化三角形 %zu 枚\n",
                g.t_inserted, g.max_t_per_edge, g.t_dropped);
    std::printf("    深さの差  最大 %zu 段（2 段以上が出たケース %zu 件）\n", g.max_depth_gap,
                g.cases_with_gap);

    KRI_CHECK_MSG(g.leaves_adaptive < g.leaves_fixed,
                  "葉が 1 つも減っていない。適応分割が空回りしています");
    // §13: **減っていなければ D の複雑さに見合いません。**「決定したから使う」ではなく
    // 「効いたから使う」
    KRI_CHECK_MSG(g.frags_adaptive < g.frags_fixed,
                  "断片が減っていない。§13 のとおり、これでは案 D の複雑さに見合いません");
    KRI_CHECK_MSG(g.max_depth_gap >= 2,
                  "葉の深さの差が 2 段に届かない。**§2.4 の前提が作れていません**");
    KRI_CHECK_MSG(g.t_inserted > 0,
                  "T 頂点が 1 個も入っていない。**案 D が効くべき場所で効いていません**");
}

/// §8 のケース 13 の要件。**変異 11 の前提です。**
void check_case13_requirement() {
    std::printf("\n  §8 ケース 13 の要件\n");
    std::printf("    深さの差 %zu 段 / 1 本の辺の T 頂点の最大 %zu 個\n", g_case13_depth_gap,
                g_case13_max_t_per_edge);
    KRI_CHECK_MSG(g_case13_depth_gap >= 2,
                  "ケース 13 で隣接セルの深さが 2 段以上違っていない（§8 の本丸）");
    KRI_CHECK_MSG(g_case13_max_t_per_edge >= 2,
                  "ケース 13 の 1 本の辺に T 頂点が 2 個以上載っていない。"
                  "**§9.3 の変異 11 が検出できず、変異が無意味になります**（§8 の要件）");
}

}  // namespace

int main() {
    std::printf("\n  適応分割 — SPEC-phase2 §3 / §9.1 / CP2\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) run_case(c);
    check_not_vacuous();
    check_case13_requirement();
    std::printf("\n");
    return kritest::finish("csg/adaptive");
}

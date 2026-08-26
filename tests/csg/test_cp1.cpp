// Krisite — CP1: ケース 1（一般位置）、union のみ、深度 0 と 2
//
// SPEC-phase1.md §11 CP1
//
//   通れば: パイプラインの骨格は成立
//   落ちれば: 縫合以前の問題。設計を見直す
#include <cstdio>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;
using kritest::cube;

namespace {

/// 位相検査（§10.1）を通し、結果を返す。
TopologyReport run(const TriMesh& a, const TriMesh& b, BoolOp op, unsigned depth, BoolStats& st) {
    const BoolMesh r = boolean_op(a, b, op, depth, &st);
    return check_topology(r.triangles);
}

void report(const char* what, const TopologyReport& t, const BoolStats& st, unsigned depth) {
    std::printf("  %-22s 深度%u  V=%-5zu E=%-5zu F=%-5zu C=%zu g=%lld  断片=%zu 構成点=%zu\n", what,
                depth, t.v, t.e, t.f, t.components, t.genus_total, st.fragments,
                st.constructed_points);
}

// ---- CP1 -------------------------------------------------------------------

void test_cp1_union_general_position() {
    // ケース 1: 立方体 2 個、一般位置（重なるが面・辺・頂点は一致しない）
    //
    // **座標範囲に見合った大きさにすること。** 八分木は座標範囲全体を分割するので、
    // 隅に小さく置くと深度を上げてもセル境界が幾何をまたがない。
    // b=21 の深度 2 の内部境界は -2^19, 0, +2^19 で、下の 2 つはそれを実際にまたぐ。
    const TriMesh a = kritest::grid_scale_box(-600000, 200000);
    const TriMesh b = kritest::grid_scale_box(-100000, 700000);

    KRI_CHECK(check_topology(a).ok() && check_topology(b).ok());
    KRI_CHECK(krisite::mesh::is_outward_oriented(a));
    KRI_CHECK(krisite::mesh::is_outward_oriented(b));

    std::printf("\n  CP1: ケース 1（一般位置）× union\n");
    TopologyReport ref;
    bool first = true;
    for (unsigned d : {0u, 2u}) {
        BoolStats st;
        const TopologyReport t = run(a, b, BoolOp::Union, d, st);
        report("立方体2個 ∪", t, st, d);

        KRI_CHECK_MSG(t.f > 0, "union の出力が空");
        KRI_CHECK_MSG(t.edge_manifold, "辺多様体でない（深度 " + std::to_string(d) + "）");
        KRI_CHECK_MSG(t.vertex_manifold, "頂点多様体でない（深度 " + std::to_string(d) + "）");
        KRI_CHECK_MSG(t.oriented, "向きが整合しない（深度 " + std::to_string(d) + "）");
        KRI_CHECK_MSG(t.no_degenerate, "退化三角形がある（深度 " + std::to_string(d) + "）");
        KRI_CHECK_MSG(t.ok(), "位相検査に落ちた（深度 " + std::to_string(d) + "）");

        // 重なる 2 立方体の和は単一シェル・種数 0
        KRI_CHECK_MSG(t.components == 1, "C = " + std::to_string(t.components) + "（期待 1、深度 " +
                                             std::to_string(d) + "）");
        KRI_CHECK_MSG(t.genus_total == 0, "g = " + std::to_string(t.genus_total) +
                                              "（期待 0、深度 " + std::to_string(d) + "）");

        // §10.2.1 深度不変性: (C, g_total) が深度によらないこと
        if (first) {
            ref = t;
            first = false;
        } else {
            KRI_CHECK_MSG(t.components == ref.components, "深度不変性: C が深度で変わった");
            KRI_CHECK_MSG(t.genus_total == ref.genus_total, "深度不変性: g_total が深度で変わった");
        }
    }
}

/// 分割が実際に効いていること（深度 2 では断片が増える）。
void test_depth_actually_splits() {
    const TriMesh a = kritest::grid_scale_box(-600000, 200000);
    const TriMesh b = kritest::grid_scale_box(-100000, 700000);
    BoolStats s0, s2;
    boolean_op(a, b, BoolOp::Union, 0, &s0);
    boolean_op(a, b, BoolOp::Union, 2, &s2);
    KRI_CHECK_MSG(s2.fragments > s0.fragments,
                  "深度を上げても断片が増えない（分割が効いていない）");
    std::printf("  断片数: 深度0 = %zu → 深度2 = %zu\n", s0.fragments, s2.fragments);
    std::printf("  1 点に集まる平面の最大枚数: 深度0 = %zu / 深度2 = %zu\n", s0.max_planes_at_point,
                s2.max_planes_at_point);
    std::printf("  第2段が併合した点: 深度0 = %zu / 深度2 = %zu\n", s0.merged_by_value,
                s2.merged_by_value);
    std::printf("  領域数（符号ベクトル）: 深度0 = %zu / 深度2 = %zu\n", s0.regions, s2.regions);
}

/// 離れた 2 立方体の union は 2 シェルになる（分類が効いていることの確認）。
void test_disjoint_union() {
    const TriMesh a = kritest::grid_scale_box(-900000, -600000);
    const TriMesh b = kritest::grid_scale_box(600000, 900000);
    for (unsigned d : {0u, 2u}) {
        BoolStats st;
        const TopologyReport t = run(a, b, BoolOp::Union, d, st);
        report("離れた2立方体 ∪", t, st, d);
        KRI_CHECK(t.ok());
        KRI_CHECK_MSG(t.components == 2, "離れた 2 立方体の union は 2 シェル");
        KRI_CHECK(t.genus_total == 0);
    }
}

}  // namespace

int main() {
    test_cp1_union_general_position();
    test_depth_actually_splits();
    test_disjoint_union();
    std::printf("\n");
    return kritest::finish("csg/cp1");
}

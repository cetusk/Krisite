// Krisite — CP2: ケース 5（面がセル境界と一致）、union、深度 0〜3
//
// SPEC-phase1.md §11 CP2
//
//   通れば: EMBER 続行の判断が立つ
//   落ちれば: §5.3 の第2段を実装して再挑戦
//
// **あわせて §5.4 の数値を報告します。** とくに「1 点に集まる平面の最大枚数」が
// 5 枚以上になるかが、Phase 3（並列化）の設計を左右します。
//
// 順序について。ケース 5 は共平面重複と 4 平面同時交差を**同時に**突くので、
// 先に共平面重複だけを突くケース 2 / 8 を通してから当たります。そうしないと
// 落ちたときに原因が切り分けられません。
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

TopologyReport run(const TriMesh& a, const TriMesh& b, BoolOp op, unsigned depth, BoolStats& st) {
    const BoolMesh r = boolean_op(a, b, op, depth, &st);
    return check_topology(r.triangles);
}

/// 位相が閉じた向き付き多様体で、(C, g) が期待どおりか。
void expect(const std::string& what, const TopologyReport& t, std::size_t comps, long long genus,
            unsigned d) {
    const std::string tag = what + "（深度 " + std::to_string(d) + "）";
    KRI_CHECK_MSG(t.edge_manifold, tag + ": 辺多様体でない");
    KRI_CHECK_MSG(t.vertex_manifold, tag + ": 頂点多様体でない");
    KRI_CHECK_MSG(t.oriented, tag + ": 向きが整合しない");
    KRI_CHECK_MSG(t.no_degenerate, tag + ": 退化三角形がある");
    KRI_CHECK_MSG(t.ok(), tag + ": 位相検査に落ちた");
    KRI_CHECK_MSG(t.components == comps, tag + ": C = " + std::to_string(t.components) + "（期待 " +
                                             std::to_string(comps) + "）");
    KRI_CHECK_MSG(t.genus_total == genus, tag + ": g = " + std::to_string(t.genus_total) +
                                              "（期待 " + std::to_string(genus) + "）");
}

// ---- §9.0 (1) サイズ規律 ----------------------------------------------------

void test_size_discipline() {
    std::printf("\n  §9.0 (1) サイズ規律（両入力の AABB の和 >= 座標範囲の半分）\n");
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        const bool ok = kritest::size_discipline_ok(a, b);
        std::printf("    ケース %-2s %-24s %s\n", c.id, c.what, ok ? "OK" : "**違反**");
        KRI_CHECK_MSG(ok, std::string("ケース ") + c.id + ": §9.0 のサイズ規律に違反");
        // 入力そのものの健全性
        KRI_CHECK_MSG(check_topology(a).ok() && check_topology(b).ok(),
                      std::string("ケース ") + c.id + ": 入力が閉じた多様体でない");
        KRI_CHECK_MSG(
            krisite::mesh::is_outward_oriented(a) && krisite::mesh::is_outward_oriented(b),
            std::string("ケース ") + c.id + ": 入力の向きが外向きでない");
        KRI_CHECK_MSG(krisite::mesh::coords_in_range(a) && krisite::mesh::coords_in_range(b),
                      std::string("ケース ") + c.id + ": 入力が座標範囲外");
    }
}

// ---- 共平面重複（§4.3.2 の符号 0）------------------------------------------
//
// ケース 5 に当たる前に、共平面重複**だけ**を突くケースで処理が効くことを確かめます。

/// ケース 8: 同一の立方体。全 6 平面を共有し、すべて同方向。
///
/// $A \cup A = A$、$A \cap A = A$、$A \setminus A = \emptyset$（§9.1）。
void test_case8_identical() {
    const TriMesh a = kritest::cases::case8();
    std::printf("\n  ケース 8（同一の立方体）— 共平面重複が全面\n");
    for (unsigned d : {0u, 1u, 2u, 3u}) {
        BoolStats su, si, sd;
        const TopologyReport u = run(a, a, BoolOp::Union, d, su);
        const TopologyReport i = run(a, a, BoolOp::Intersection, d, si);
        const BoolMesh diff = boolean_op(a, a, BoolOp::Difference, d, &sd);

        std::printf(
            "    深度%u  ∪:C=%zu g=%lld  ∩:C=%zu g=%lld  \\:F=%zu  "
            "共平面(同/逆)=%zu/%zu  断片=%zu\n",
            d, u.components, u.genus_total, i.components, i.genus_total, diff.triangles.size(),
            su.coplanar_same, su.coplanar_opposite, su.fragments);

        expect("ケース8 ∪", u, 1, 0, d);
        expect("ケース8 ∩", i, 1, 0, d);
        KRI_CHECK_MSG(diff.triangles.empty(),
                      "ケース8: A\\A が空でない（深度 " + std::to_string(d) + "）");
        // すべての断片が共平面の対をなす（A と B が同一なので当然）
        KRI_CHECK_MSG(su.coplanar_opposite == 0, "ケース8: 逆方向の共平面が出た");
        KRI_CHECK_MSG(
            su.coplanar_same * 2 == su.fragments,
            "ケース8: 共平面の対が全断片を覆っていない（深度 " + std::to_string(d) + "）");
        // ∪ と ∩ は A そのもの。三角形数が一致すること
        KRI_CHECK_MSG(u.f == i.f, "ケース8: ∪ と ∩ の面数が違う");
    }
}

/// ケース 2: 面が完全共平面（$z = -2^{b-1}$ を共有、同方向）。
void test_case2_coplanar() {
    const TriMesh a = kritest::cases::case2_a(), b = kritest::cases::case2_b();
    std::printf("\n  ケース 2（面が完全共平面）\n");
    for (unsigned d : {0u, 1u, 2u, 3u}) {
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            BoolStats st;
            const TopologyReport t = run(a, b, op, d, st);
            std::printf("    深度%u %s  C=%zu g=%-2lld F=%-5zu  共平面(同/逆)=%zu/%zu 断片=%zu\n",
                        d, op_name(op), t.components, t.genus_total, t.f, st.coplanar_same,
                        st.coplanar_opposite, st.fragments);
            expect(op == BoolOp::Union          ? "ケース2 ∪"
                   : op == BoolOp::Intersection ? "ケース2 ∩"
                                                : "ケース2 \\",
                   t, 1, 0, d);
            KRI_CHECK_MSG(st.coplanar_same > 0,
                          "ケース2: 同方向の共平面重複が 1 つも検出されていない（深度 " +
                              std::to_string(d) + "）");
        }
    }
}

// ---- CP2 本体 --------------------------------------------------------------

void test_cp2_case5() {
    const TriMesh a = kritest::cases::case5_a(), b = kritest::cases::case5_b();
    // **3 演算すべてを回します**（第7版 §11）。第6版までは union のみでしたが、
    // それが選択規則と向き付けのバグ 2 件を隠していました。
    std::printf("\n  CP2: ケース 5（面がセル境界と完全一致）× 3 演算\n");

    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        TopologyReport ref;
        bool first = true;
        for (unsigned d : {0u, 1u, 2u, 3u}) {
            BoolStats st;
            const TopologyReport t = run(a, b, op, d, st);
            std::printf(
                "    深度%u %s V=%-5zu E=%-5zu F=%-5zu C=%zu g=%lld  断片=%-5zu "
                "重複=%-4zu 有効セル=%zu/%zu\n",
                d, op_name(op), t.v, t.e, t.f, t.components, t.genus_total, st.fragments,
                st.duplicate_fragments, st.active_cells, st.total_cells);
            KRI_CHECK_MSG(t.f > 0, std::string("ケース5 ") + op_name(op) + " の出力が空");
            expect(std::string("ケース5 ") + op_name(op), t, 1, 0, d);

            // §10.2.1 深度不変性
            if (first) {
                ref = t;
                first = false;
            } else {
                KRI_CHECK_MSG(t.components == ref.components, "深度不変性: C が深度で変わった");
                KRI_CHECK_MSG(t.genus_total == ref.genus_total, "深度不変性: g が深度で変わった");
            }
        }
    }
}

/// ケース 6（1 格子ずらした対照）。ケース 5 と同じ位相になること。
void test_case6_contrast() {
    const TriMesh a = kritest::cases::case6_a(), b = kritest::cases::case6_b();
    std::printf("\n  ケース 6（セル境界から 1 格子ずれ）— ケース 5 との対照\n");
    for (unsigned d : {0u, 1u, 2u, 3u}) {
        BoolStats st;
        const TopologyReport t = run(a, b, BoolOp::Union, d, st);
        std::printf("    深度%u  C=%zu g=%lld F=%-5zu 断片=%-5zu 重複=%zu 共平面=%zu\n", d,
                    t.components, t.genus_total, t.f, st.fragments, st.duplicate_fragments,
                    st.coplanar_same + st.coplanar_opposite);
        expect("ケース6 ∪", t, 1, 0, d);
    }
}

// ---- §5.4 の数値 -----------------------------------------------------------

void report_54() {
    std::printf("\n  §5.4 の実測（ケース別・深度別）\n");
    std::printf("    %-3s %-3s %-8s %-8s %-8s %-8s %-6s %-6s %-8s %s\n", "#", "深度", "断片",
                "重複断片", "構成点", "併合後", "値併合", "最大枚数", "うちmesh", "有効セル");
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (unsigned d : {0u, 1u, 2u, 3u}) {
            BoolStats st;
            boolean_op(a, b, BoolOp::Union, d, &st);
            std::printf("    %-3s %-4u %-8zu %-8zu %-8zu %-8zu %-8zu %-8zu %-8zu %zu/%zu\n", c.id,
                        d, st.fragments, st.duplicate_fragments, st.constructed_points,
                        st.merged_points, st.merged_by_value, st.max_planes_at_point,
                        st.max_mesh_planes_at_point, st.active_cells, st.total_cells);
            // 理論上界の見張り: 1 点に集まる平面は総平面数を超えない
            KRI_CHECK(st.max_planes_at_point <= st.planes_total);
            // §5.4（第8版）: 併合グループは 1 セルとその面・辺・頂点隣接に収まること。
            // **2 以上が出たら併合の誤り**（本来別の点を同一視した）を疑う。
            // これが成り立つなら Phase 3 は「セル並列 + 境界併合」で組める。
            KRI_CHECK_MSG(st.max_merge_span <= 1,
                          std::string("ケース ") + c.id + "（深度 " + std::to_string(d) +
                              "）: 併合グループがセル隣接を越えて広がっている（広がり " +
                              std::to_string(st.max_merge_span) +
                              "）。Phase 3 の「セル並列 + 境界併合」の前提が崩れます");
        }
    }
}

}  // namespace

int main() {
    test_size_discipline();
    test_case8_identical();
    test_case2_coplanar();
    test_cp2_case5();
    test_case6_contrast();
    report_54();
    std::printf("\n");
    return kritest::finish("csg/cp2");
}

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
void expect(const char* what, const TopologyReport& t, std::size_t comps, long long genus,
            unsigned d) {
    const std::string tag = std::string(what) + "（深度 " + std::to_string(d) + "）";
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
    std::printf("\n  CP2: ケース 5（面がセル境界と完全一致）× union\n");

    TopologyReport ref;
    bool first = true;
    for (unsigned d : {0u, 1u, 2u, 3u}) {
        BoolStats st;
        const TopologyReport t = run(a, b, BoolOp::Union, d, st);
        std::printf("    深度%u  V=%-5zu E=%-5zu F=%-5zu C=%zu g=%lld  断片=%-5zu 重複=%zu\n", d,
                    t.v, t.e, t.f, t.components, t.genus_total, st.fragments,
                    st.duplicate_fragments);
        KRI_CHECK_MSG(t.f > 0, "ケース5 ∪ の出力が空");
        expect("ケース5 ∪", t, 1, 0, d);

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

// ---- ケース 5T: §5.4 の計測が 4 枚以上を検出できることの証拠 ★ ---------------
//
// **軸平行な立方体だけでは 1 点に集まる平面が 3 枚を超えられません。**
// 1 軸につき点を通る平面は高々 1 枚（平行な 2 平面は交わらない）で、軸は 3 本だけ
// だからです。セル面が立方体の面と一致しても、同一の幾何平面なので `PlaneTable` で
// 同じ ID になり、枚数は増えません。
//
// したがって「ケース 5 で最大 3 枚だった」という報告は、それだけでは
// **計測器が壊れていても同じ結果になります。** 4 枚を実際に作って検出を確かめます。
//
// **位相は合否に使いません**（§9.3 のケース 11b と同じ扱い）。斜面を持つ入力では
// 現状の実装が非多様体を出しますが、これは CP2 の問いではなく CP3（ケース 7 / 9）の
// 問題であり、**この改訂の前から同じです**（旧実装は同じ入力で停止します）。
void test_case5t_detects_four_planes() {
    const TriMesh a = kritest::cases::case5t_a(), b = kritest::cases::case5t_b();
    std::printf("\n  ケース 5T（セル角を斜面が通る）— 計測の検出力の確認\n");
    std::size_t best_mesh = 0, best_all = 0, merged = 0;
    for (unsigned d : {0u, 1u, 2u, 3u}) {
        BoolStats st;
        const TopologyReport t = run(a, b, BoolOp::Union, d, st);
        std::printf(
            "    深度%u  最大枚数=%zu（mesh のみ %zu / 平面総数 %zu）値併合=%zu "
            "中点=%zu  ★位相 C=%zu g=%lld ok=%d（合否に使わない）\n",
            d, st.max_planes_at_point, st.max_mesh_planes_at_point, st.planes_total,
            st.merged_by_value, st.midpoint_raycasts, t.components, t.genus_total,
            static_cast<int>(t.ok()));
        best_mesh = std::max(best_mesh, st.max_mesh_planes_at_point);
        best_all = std::max(best_all, st.max_planes_at_point);
        merged = std::max(merged, st.merged_by_value);
    }
    KRI_CHECK_MSG(best_mesh >= 4,
                  "§5.4 の計測が 4 平面同時交差を検出できていない（意図的に作ったのに 3 枚以下）");
    KRI_CHECK_MSG(best_all >= best_mesh, "セル面込みの枚数が mesh のみを下回っている");
    KRI_CHECK_MSG(merged > 0,
                  "第2段（値ベースの併合）が発火していない。4 平面同時交差があれば "
                  "平面3つ組が相異なるのに同じ点になる頂点が出るはず（§5.3）");
}

// ---- §5.4 の数値 -----------------------------------------------------------

void report_54() {
    std::printf("\n  §5.4 の実測（ケース別・深度別）\n");
    std::printf("    %-3s %-3s %-8s %-8s %-8s %-8s %-6s %-6s %-6s\n", "#", "深度", "断片",
                "重複断片", "構成点", "併合後", "値併合", "最大枚数", "うちmesh");
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (unsigned d : {0u, 1u, 2u, 3u}) {
            BoolStats st;
            boolean_op(a, b, BoolOp::Union, d, &st);
            std::printf("    %-3s %-4u %-8zu %-8zu %-8zu %-8zu %-8zu %-8zu %-6zu\n", c.id, d,
                        st.fragments, st.duplicate_fragments, st.constructed_points,
                        st.merged_points, st.merged_by_value, st.max_planes_at_point,
                        st.max_mesh_planes_at_point);
            // 理論上界の見張り: 1 点に集まる平面は総平面数を超えない
            KRI_CHECK(st.max_planes_at_point <= st.planes_total);
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
    test_case5t_detects_four_planes();
    report_54();
    std::printf("\n");
    return kritest::finish("csg/cp2");
}

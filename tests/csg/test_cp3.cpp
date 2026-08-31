// Krisite — CP3: 全ケース × 3 演算 × 深度 0〜3
//
// SPEC-phase1.md §11 CP3
//
// **§9.3 の除外は §9.3.1 の条件を満たすときだけ適用します。** 条件を満たさない
// なら、それは除外すべき退化ではなくバグです。除外の件数も数えて報告します
// （多数に及ぶなら §2.1「出力は多様体」という前提のほうが誤り）。
#include <algorithm>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "corpus_expect.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;
using kritest::Exclusion;

namespace {

constexpr unsigned kMaxDepth = 3;

/// §5.3 の第2段が効いていること: **出力に値の重複する頂点が無い。**
///
/// 第2段は「全構成点を lex_less で整列し、値が厳密に等しいものを併合する」ので、
/// 併合後に等しい値の頂点が 2 つ以上残ることはありません。
///
/// **これは冪等性の比較では検出できません。** 比較側も値で同一視するため、
/// 第2段を止めても一致してしまいます（`test_idempotence.cpp` の注記）。
/// ここは値の重複そのものを見るので、第2段を止めれば必ず落ちます。
std::size_t duplicate_vertex_values(const BoolMesh& m) {
    std::vector<std::uint32_t> ord(m.vertices.size());
    for (std::uint32_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](std::uint32_t a, std::uint32_t b) {
        return krisite::geom::lex_less(m.vertices[a], m.vertices[b]);
    });
    std::size_t dup = 0;
    for (std::size_t i = 1; i < ord.size(); ++i) {
        if (krisite::geom::h_equal(m.vertices[ord[i - 1]], m.vertices[ord[i]])) ++dup;
    }
    return dup;
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

/// §9.3.1 の番人: 除外した **(ケース, 演算, 接触の次元)** の組。
///
/// 深度ごとに数えると 4 倍に見えてしまうので、構成の数で数えます。
/// **§13 は接触の次元別に記録することを求めています**（第11版）。
/// **2〜3 件なら特殊ケース。多数に及ぶなら §2.1 の前提を見直します。**
std::vector<std::pair<std::string, Exclusion>> excluded_configs;

std::size_t checked = 0;

void run_case(const kritest::Case& c) {
    const TriMesh a = c.make_a(), b = c.make_b();

    KRI_CHECK_MSG(kritest::size_discipline_ok(a, b),
                  std::string("ケース ") + c.id + ": §9.0 (1) のサイズ規律に違反");
    KRI_CHECK_MSG(check_topology(a).ok() && check_topology(b).ok(),
                  std::string("ケース ") + c.id + ": 入力が閉じた多様体でない");
    KRI_CHECK_MSG(krisite::mesh::is_outward_oriented(a) && krisite::mesh::is_outward_oriented(b),
                  std::string("ケース ") + c.id + ": 入力の向きが外向きでない");
    KRI_CHECK_MSG(krisite::mesh::coords_in_range(a) && krisite::mesh::coords_in_range(b),
                  std::string("ケース ") + c.id + ": 入力が座標範囲外");

    std::printf("\n  ケース %-4s %s\n", c.id, c.what);
    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        const kritest::ExpectedTopo want = kritest::expected_topo(c.id, op);
        // **除外は出力の性質から決めます**（識別子を使わない。`corpus_expect.hpp`）
        Exclusion ex = Exclusion::None;
        bool ex_recorded = false;

        TopologyReport ref;
        bool first = true;
        for (unsigned d = 0; d <= kMaxDepth; ++d) {
            BoolStats st;
            // **Phase 1 の意味論のまま（分裂 OFF）で見ます**（SPEC-phase2 §5.2）。
            // 分裂は設計どおり頂点の値を複製するので、§5.3 の第2段を突く
            // 「値の重複が無いこと」の検査と衝突します。**分裂側は
            // `tests/csg/test_split_semantics.cpp` が受け持ちます。**
            const BoolOptions opt = kritest::corpus_options(d);
            const BoolMesh r = boolean_op(a, b, op, opt, &st);
            const TopologyReport t = check_topology(r.triangles);
            const Exclusion ex_d = kritest::exclusion_from_topology(t);
            if (ex_d != Exclusion::None) {
                ex = ex_d;
                if (!ex_recorded) {
                    excluded_configs.emplace_back(std::string(c.id) + " " + op_name(op), ex);
                    ex_recorded = true;
                }
            }
            const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) + "（深度 " +
                                    std::to_string(d) + "）";

            std::printf(
                "    d%u %s C=%zu g=%-2lld χ=%-3lld F=%-5zu %-16s 断片%-5zu "
                "枚数%zu(m%zu) 値併合%zu 広%zu\n",
                d, op_name(op), t.components, t.genus_total, t.chi, t.f,
                ex == Exclusion::VertexContact ? "頂点接触（§9.3）"
                : ex == Exclusion::EdgeContact ? "辺接触（§9.3）"
                : t.ok()                       ? "ok"
                                               : "**NG**",
                st.fragments, st.max_planes_at_point, st.max_mesh_planes_at_point,
                st.merged_by_value, st.max_merge_span);

            // §5.3: 第2段が効いていること（値の重複した頂点が残っていない）
            KRI_CHECK_MSG(duplicate_vertex_values(r) == 0,
                          tag + ": 出力に値の重複する頂点が " +
                              std::to_string(duplicate_vertex_values(r)) +
                              " 組ある。§5.3 の第2段が効いていません");

            if (ex_d != Exclusion::None) {
                std::string why;
                KRI_CHECK_MSG(kritest::exclusion_conditions_ok(ex, t, &why),
                              tag + ": §9.3.1 の適用条件を満たさない（" + why +
                                  "）。除外すべき退化ではなくバグです");
            } else {
                KRI_CHECK_MSG(t.ok(), tag + ": 位相検査に落ちた");
                if (!t.empty) {
                    KRI_CHECK_MSG(t.edge_manifold, tag + ": 辺多様体でない");
                    KRI_CHECK_MSG(t.vertex_manifold, tag + ": 頂点多様体でない");
                    KRI_CHECK_MSG(t.oriented, tag + ": 向きが整合しない");
                    KRI_CHECK_MSG(t.no_degenerate, tag + ": 退化三角形がある");
                }
                ++checked;
            }

            // §9.2 / §9.1 が明示している期待値
            if (want.known) {
                KRI_CHECK_MSG(t.empty == want.empty, tag + ": 空かどうかが期待と違う（得 " +
                                                         std::to_string(t.empty) + " 期待 " +
                                                         std::to_string(want.empty) + "）");
                if (!want.empty) {
                    KRI_CHECK_MSG(t.components == want.components,
                                  tag + ": C = " + std::to_string(t.components) + "（期待 " +
                                      std::to_string(want.components) + "）");
                    KRI_CHECK_MSG(t.genus_total == want.genus_total,
                                  tag + ": g = " + std::to_string(t.genus_total) + "（期待 " +
                                      std::to_string(want.genus_total) + "）");
                }
            }

            // §10.2.1 深度不変性
            if (first) {
                ref = t;
                first = false;
            } else {
                KRI_CHECK_MSG(t.components == ref.components,
                              tag + ": 深度不変性 — C が深度で変わった（" +
                                  std::to_string(ref.components) + " → " +
                                  std::to_string(t.components) + "）");
                KRI_CHECK_MSG(t.genus_total == ref.genus_total,
                              tag + ": 深度不変性 — g が深度で変わった（" +
                                  std::to_string(ref.genus_total) + " → " +
                                  std::to_string(t.genus_total) + "）");
            }
        }
    }
}

}  // namespace

/// **§9.4 のジョブが空回りしていないことの番人。**
///
/// 既定を反転させても、テストが構成を明示していれば何も変わりません。**実際に一度、
/// 全テストを固定深度に固定したせいでジョブが無意味になりかけました。**
/// 適応分割モードでは**葉の深さに差が出ている**ことを確かめます。
void check_adaptive_mode_effective() {
    if (!kritest::adaptive_mode()) return;
    unsigned gap = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        BoolStats st;
        boolean_op(a, b, BoolOp::Union, kritest::corpus_options(3), &st);
        gap = std::max(gap, st.leaf_depth_max - st.leaf_depth_min);
    }
    std::printf("\n  §9.4 適応分割モード: 葉の深さの差 最大 %u 段\n", gap);
    KRI_CHECK_MSG(gap >= 2, "適応分割モードなのに葉の深さに差が出ていない。**空回りです**");
}

int main() {
    std::printf("\n  CP3 — 全ケース × 3 演算 × 深度 0〜3（SPEC-phase1 §11）\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) run_case(c);

    // §13: **接触の次元別に**記録する（第11版）
    std::size_t n_dim0 = 0, n_dim1 = 0;
    for (const auto& kv : excluded_configs) {
        if (kv.second == Exclusion::VertexContact) ++n_dim0;
        if (kv.second == Exclusion::EdgeContact) ++n_dim1;
    }
    std::printf("\n  §9.3.1 の除外（%zu 構成 = 次元 0 が %zu、次元 1 が %zu）\n",
                excluded_configs.size(), n_dim0, n_dim1);
    for (const auto& kv : excluded_configs) {
        std::printf("    %-12s 接触の次元 %d\n", kv.first.c_str(),
                    kv.second == Exclusion::VertexContact ? 0 : 1);
    }
    std::printf("  合否に使った組数: %zu\n", checked);

    // ---- §9.3.1 の番人: **件数の上限ではなく、増分の説明を要求します** ★ ----
    //
    // **件数の閾値を番人にすると、上げ続けることになります**（`CLAUDE.md`）。
    // コーパスは退化を狙って設計しているので、接触を持つケースを足せば除外は増えます。
    // 「上限を超えたら落ちる」ではなく「**増えたら理由を書かせる**」形にします。
    //
    // **由来のケースが「接触を意図して作ったもの」なら正当です。**
    // 意図していないケースから増えたなら、§2.1「出力は多様体」の前提を見直します。
    //
    // > **前提そのものの妥当性は、コーパスではなく実データで測ります。**
    // > CP1 の実測: 336 対中 3 件が接触由来 = **0.9%**。1% 未満なら §2.1 は正しい既定です。
    //
    // **接触を意図して作ったケース**（`SPEC-phase1.md` §9.3.1 の由来表）。
    // ここに無い ID から除外が出たら、それは説明の付かない増分です。
    static const std::set<std::string> kContactByDesign = {
        "4T",   // 頂点だけを共有（次元 0）
        "4T'",  // 頂点共有 + 体積も重なる
        "11b",  // 辺だけを共有（次元 1）
        "16",   // 自己接触（平板 − 菱形プリズム）
        "24",   // 連結成分で分けられない自己接触（Phase 5 の CP1 が実データで到達）
        "24'",  // 24 の対照（橋なし）
    };
    for (const auto& kv : excluded_configs) {
        const std::string id = kv.first.substr(0, kv.first.find(' '));
        KRI_CHECK_MSG(kContactByDesign.count(id) != 0,
                      "§9.3 の除外が **接触を意図していないケース** " + id + " から出ました（" +
                          kv.first +
                          "）。ケースを除外するのではなく、"
                          "§2.1「出力は多様体」という前提のほうを見直してください");
    }
    // **件数は記録として残します。** 判定材料は由来のほうです
    std::printf("  除外の由来: すべて接触を意図したケース（%zu 種）\n", kContactByDesign.size());
    // 空回り防止: 全ケース × 3 演算 × 深度 4 − 除外
    const std::size_t expect =
        kritest::corpus().size() * 3 * (kMaxDepth + 1) - excluded_configs.size() * (kMaxDepth + 1);
    KRI_CHECK_MSG(checked == expect, "合否に使った組数が期待と違う（" + std::to_string(checked) +
                                         " 対 " + std::to_string(expect) + "）");
    std::printf("\n");
    check_adaptive_mode_effective();
    return kritest::finish("csg/cp3");
}

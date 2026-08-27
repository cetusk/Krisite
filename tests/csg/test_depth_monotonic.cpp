// Krisite — §9.0 の番人: 断片数が深度 0→1→2→3 で狭義単調増加すること
//
// SPEC-phase1.md §9.0 (2), §13
//
// **これは「空回りしている緑」を検出するための検査です。**
//
// 八分木は座標範囲全体を分割するので、入力が範囲に比して小さいとすべてが 1 セルに
// 収まり、深度を上げても断片が増えません。このとき §10.2.1 の深度不変性は**通ります。
// しかし何も検証していません。** CP1 で実際に踏みました（IMPL-phase1 §6.1）。
//
// サイズ規律を手で守るのは必ずどこかで間違えるので、機械化します。
// **配置をどう間違えても、断片が増えなければここで落ちます。**
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kMaxDepth = 3;

/// **判定基準を「各段で狭義増加」から「非減少 かつ 深度 0 → 3 で狭義増加」に緩めています。**
///
/// §9.0 は各段での狭義単調増加を要求していますが、正準化後の断片数では成立しません。
/// 面がセル境界にちょうど乗るケース（5, 8, および 2 の一部）では
///
///   - セル境界平面が入力メッシュの平面と**同一**なので `PlaneTable` で同じ ID になり、
///     §4.3.1 の大域分割が既にその平面で切っています。深度を上げても
///     **新しい切断面が 1 枚も増えません**
///   - 増えるのは、面が境界に乗ることで両側のセルに割り当てられる**重複断片**だけ
///
/// ケース 5 は深度 1 の境界（$0$）も深度 2 の境界（$\pm 2^{b-2}$）も入力の面と一致するので、
/// 実際に新しい切断が入るのは深度 3 からです: 正準 48 → 48 → 48 → 192。
/// **これはケース 5 の定義そのものから来る現象で、配置の誤りではありません。**
///
/// **生の断片数（重複を落とす前）で判定するのは誤りです。** 規律違反の配置でも、
/// 面がたまたま境界に乗っていれば重複だけで数が増え、空回りを見逃します
/// （立方体を 0〜633 に置いた配置は $x=0$ に面が乗るので 48 → 60 と増えます）。
///
/// §9.0 は「増加しないケースがあれば、理由を記録して仕様の承認を取ること」と
/// しています。**この判定基準の変更は仕様の承認待ちです**（IMPL-phase1 §7.7）。
void test_monotonic() {
    std::printf(
        "\n  §9.0 の番人: 断片数の深度単調増加"
        "（正準化後で判定。非減少 かつ 深度 0 → 3 で狭義増加）\n");
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();

        // 演算によらず断片の生成は同じですが、選択が変わっても計数が変わらないことを
        // 確かめる意味も込めて 3 演算すべてで見ます。
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            std::size_t prev = 0, first = 0, last = 0;
            std::string raw_line, can_line;
            for (unsigned d = 0; d <= kMaxDepth; ++d) {
                BoolStats st;
                boolean_op(a, b, op, d, &st);
                if (d == 0) {
                    raw_line = std::to_string(st.raw_fragments);
                    can_line = std::to_string(st.fragments);
                    first = st.fragments;
                } else {
                    raw_line += " → " + std::to_string(st.raw_fragments);
                    can_line += " → " + std::to_string(st.fragments);
                    KRI_CHECK_MSG(st.fragments >= prev, std::string("ケース ") + c.id +
                                                            ": 断片数が深度で減った（" +
                                                            std::to_string(prev) + " → " +
                                                            std::to_string(st.fragments) + "）");
                }
                prev = st.fragments;
                last = st.fragments;
            }
            KRI_CHECK_MSG(last > first, std::string("ケース ") + c.id + ": 断片数が深度 0 → " +
                                            std::to_string(kMaxDepth) + " で増えない（" +
                                            std::to_string(first) + " → " + std::to_string(last) +
                                            "）。深度掃引が空回りしています（§9.0）");
            if (op == BoolOp::Union) {
                std::printf("    ケース %-3s %-24s 生 %-26s 正準 %s\n", c.id, c.what,
                            raw_line.c_str(), can_line.c_str());
            }
        }
    }
}

/// 番人自身の検出力。**規律を破った配置を渡せば落ちること。**
///
/// 検査が本当に効いているかは壊してみないと分かりません（SPEC-phase1 §10.5）。
/// ここでは落とさずに、断片数が実際に増えないことだけを確かめます。
void test_guard_detects_undersized() {
    // 座標範囲の隅に小さく置いた立方体 2 個。CP1 で踏んだ配置そのもの。
    const TriMesh a = kritest::grid_scale_box(0, 633);
    const TriMesh b = kritest::grid_scale_box(200, 800);
    KRI_CHECK_MSG(!kritest::size_discipline_ok(a, b),
                  "番人の検出力: この配置はサイズ規律に違反しているはず");
    BoolStats s0, s3;
    boolean_op(a, b, BoolOp::Union, 0, &s0);
    boolean_op(a, b, BoolOp::Union, kMaxDepth, &s3);
    std::printf(
        "\n  番人の検出力（規律違反の配置）: 断片 深度0=%zu 深度%u=%zu"
        "（生 %zu → %zu）\n",
        s0.fragments, kMaxDepth, s3.fragments, s0.raw_fragments, s3.raw_fragments);
    KRI_CHECK_MSG(s0.fragments == s3.fragments,
                  "規律違反の配置なのに断片が増えた。番人の前提が崩れている");
    // **生の数では見逃します。** 面が x=0 に乗るぶん重複だけが増えるためです。
    KRI_CHECK_MSG(s3.raw_fragments > s0.raw_fragments,
                  "生の断片数で判定してよいことになってしまう（この配置では増えるはず）");
}

}  // namespace

int main() {
    test_monotonic();
    test_guard_detects_undersized();
    std::printf("\n");
    return kritest::finish("csg/depth_monotonic");
}

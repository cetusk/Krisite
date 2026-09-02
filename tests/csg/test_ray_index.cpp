// Krisite — レイキャストの 2 次元索引が出力を変えないこと（CP1.5-2）
//
// SPEC-phase5.md §6.3「レイキャストの加速」
//
// **索引は厳密な絞り込みです**（`ray_index.hpp` 冒頭の単調性）。したがって
// **出力は 1 ビットも変わってはいけません。** 幾何や位相ではなく、
// **多角形の平面 ID の列**で比べます（同一版の中の比較なので、それが正しい水準）。
//
// **スープ経路（`csg::boolean`）で見ます。** 二項メッシュ経路（`boolean_op`）は
// **意図的に索引を使いません** — 正解器なので素朴に保ちます
// （`IMPL-phase5.md` §12）。**配線漏れではありません。**
//
// **機構が空回りしていないことも確かめます。** 索引を入れて三角形検査が
// 1 件も減らないなら、それは「通っている」だけで「効いている」ことにはなりません。
// **どのケースで何回減ったかを数え、0 なら落とします。**
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
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

/// 出力の同一性。**幾何ではなくバイトで比べます。**
bool same_output(const PolySoup& x, const PolySoup& y) {
    if (x.polys.size() != y.polys.size()) return false;
    for (std::size_t i = 0; i < x.polys.size(); ++i) {
        const Poly& a = x.polys[i];
        const Poly& b = y.polys[i];
        if (a.frag.support != b.frag.support || a.frag.flipped != b.frag.flipped) return false;
        if (a.src != b.src || a.tag != b.tag) return false;
        if (a.frag.edge != b.frag.edge) return false;
    }
    return true;
}

/// **索引が実際に効いた回数**（検査数が減った実行の数）。0 なら空回りです。
std::size_t effective = 0;
/// 比較した組の数。**式で持ちます**（実測で書くと、比較が減っても PASS になります）。
std::size_t compared = 0;

void run_pair(const kritest::Case& c, const TriMesh& a, const TriMesh& b) {
    std::size_t best_num = 0, best_den = 1;
    const PolySoup sa = from_mesh(a), sb = from_mesh(b);

    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        for (unsigned d = 0; d <= kMaxDepth; ++d) {
            BoolOptions off = kritest::corpus_options(d);
            off.ray_index = false;
            BoolOptions on = kritest::corpus_options(d);
            on.ray_index = true;

            BoolStats s_off, s_on;
            const PolySoup r_off = boolean(sa, sb, op, off, &s_off);
            const PolySoup r_on = boolean(sa, sb, op, on, &s_on);
            ++compared;

            KRI_CHECK_MSG(same_output(r_off, r_on), std::string("ケース ") + c.id + " " +
                                                        op_name(op) + " d" + std::to_string(d) +
                                                        ": 索引の有無で出力が変わった");
            KRI_CHECK_MSG(s_on.ray_tri_tests <= s_off.ray_tri_tests,
                          std::string("ケース ") + c.id + ": 索引で検査数が増えた（" +
                              std::to_string(s_on.ray_tri_tests) + " > " +
                              std::to_string(s_off.ray_tri_tests) + "）");
            if (s_on.ray_tri_tests < s_off.ray_tri_tests) {
                ++effective;
                if (static_cast<std::uint64_t>(s_off.ray_tri_tests) * best_den >
                    static_cast<std::uint64_t>(best_num) * s_on.ray_tri_tests) {
                    best_num = s_off.ray_tri_tests;
                    best_den = s_on.ray_tri_tests ? s_on.ray_tri_tests : 1;
                }
            }
        }
    }
    std::printf("  ケース %-3s 最良の削減 %.1f 倍  %s\n", c.id,
                best_num ? double(best_num) / double(best_den) : 1.0, c.what);
}

void run_case(const kritest::Case& c) {
    run_pair(c, c.make_a(), c.make_b());
}

}  // namespace

int main() {
    std::printf("=== レイキャストの 2 次元索引（CP1.5-2）===\n");
    for (const kritest::Case& c : kritest::corpus()) run_case(c);

    // **期待値は式で持ちます。** 比較の対象が減ったら落ちるように。
    const std::size_t expect = kritest::corpus().size() * 3 * (kMaxDepth + 1);
    KRI_EQ(compared, expect);

    // **空回りの検査。** 索引が 1 件も減らしていないなら、この機構は
    // 「通っている」だけで「効いている」ことになりません。
    KRI_CHECK_MSG(effective > 0, "索引が 1 件も検査数を減らしていない（機構が空回り）");
    std::printf("\n  比較 %zu 組 / 索引が効いた %zu 組（%.0f%%）\n", compared, effective,
                100.0 * double(effective) / double(compared));
    return kritest::finish("ray_index");
}

// Krisite — 自己交差の検査（SPEC-phase3 §5.6 の NSI）
//
// **健全側に倒した検査**なので、確かめるのは 2 方向です。
//
//   自己交差していない入力が **通る**（偽陽性が実用の邪魔をしない）
//   自己交差している入力が **落ちる**（偽陰性がない）
//
// **後者が番人です。** 「常に通る」形になっていたら、この検査は何も言っていません
// （`CLAUDE.md`「機構が空回りしていないことを別に検査する」）。
#include <cstdio>
#include <string>

#include "krisite/mesh/self_intersect.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using krisite::mesh::is_self_intersecting;
using krisite::mesh::SelfIntersectStats;
using krisite::mesh::TriMesh;

namespace {

/// **意図的に自己交差させた入力。** 2 枚の四角形を十字に交差させます。
///
/// **これが落ちなければ、検査は何も検出していません。**
TriMesh crossing_quads() {
    TriMesh m;
    const std::int32_t s = kritest::at(1, 4);
    m.vertices = {
        {-s, -s, 0}, {s, -s, 0}, {s, s, 0}, {-s, s, 0},  // z = 0 の板
        {0, -s, -s}, {0, s, -s}, {0, s, s}, {0, -s, s},  // x = 0 の板（貫く）
    };
    m.triangles = {{0, 1, 2}, {0, 2, 3}, {4, 5, 6}, {4, 6, 7}};
    return m;
}

std::size_t g_clean = 0, g_dirty = 0;

void check(const char* name, const TriMesh& m, bool want) {
    SelfIntersectStats st;
    const bool got = is_self_intersecting(m, &st);
    KRI_CHECK_MSG(got == want, std::string(name) + ": 自己交差の判定が期待と違う（得 " +
                                   (got ? "交差あり" : "交差なし") + " 期待 " +
                                   (want ? "交差あり" : "交差なし") + "）");
    if (got) {
        ++g_dirty;
    } else {
        ++g_clean;
    }
    std::printf(
        "  %-28s %s | 対 %zu（片側で棄却 %zu / 共面 %zu / 共有辺 %zu / 共有点 %zu / 一般 %zu"
        " / **保守的 %zu**）\n",
        name, got ? "**交差あり**" : "交差なし", st.pairs_broad, st.rejected_by_side,
        st.coplanar_exact, st.shared_edge_ok, st.shared_vertex_ok, st.general_exact,
        st.conservative);
}

}  // namespace

int main() {
    std::printf("=== 自己交差の検査（SPEC-phase3 §5.6）===\n");

    // **単純な入力が通ること。** 通らなければ実データを測るまでもありません
    check("立方体（12 三角形）",
          kritest::box(kritest::at(-1, 2), kritest::at(-1, 2), kritest::at(-1, 2),
                       kritest::at(1, 2), kritest::at(1, 2), kritest::at(1, 2)),
          false);
    check("細分した立方体（25）", kritest::cases::case25_a(), false);
    check("細分した立方体（25′）", kritest::cases::case25p_a(), false);
    check("でこぼこの箱（26）", kritest::cases::case26_a(), false);
    check("でこぼこの箱（26′）", kritest::cases::case26p_a(), false);
    check("四面体", kritest::tetra(kritest::at(1, 2)), false);

    // ★ **番人**: 意図的に自己交差させたものが落ちること
    check("**十字に交差した 2 板**", crossing_quads(), true);

    // コーパスの入力すべて（**自己交差しない前提で作ってあります**）
    std::size_t corpus_clean = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        for (const TriMesh& m : {c.make_a(), c.make_b()}) {
            if (!is_self_intersecting(m)) ++corpus_clean;
        }
    }
    const std::size_t want = kritest::corpus().size() * 2;
    std::printf("\n  コーパスの入力 %zu / %zu が「交差なし」\n", corpus_clean, want);
    KRI_CHECK_MSG(corpus_clean == want,
                  "コーパスの入力に自己交差ありと判定されたものがある（偽陽性が多すぎる）" +
                      kritest::pair_msg(want, corpus_clean));

    // **空回りの検査。** 1 件も検出しないなら、この検査は何も言っていません
    KRI_CHECK_MSG(g_dirty > 0, "自己交差を 1 件も検出していない（検査が空回り）");
    return kritest::finish("csg/self_intersect");
}

// Krisite — 平面とセルの閉領域の交差判定（SPEC-phase2.md §2.3）
//
// **この述語は「どのセルをどの平面で切るか」を決めます。** 誤ると継ぎ目に
// T 字接合が出ます（Phase 1 の変異 3 と同型）。
//
// 検査は 3 本立てです。
//
//   (1) 8 隅を回る正解器との一致 — **別経路**。被検体は符号で lo/hi を選んで
//       内積 2 回、正解器は 8 隅すべてを評価して最小・最大を取る
//   (2) **箱の包含に対する単調性** — §2.3 の整合性の証明はここに帰着します。
//       $F \subset \overline{C}$ なら「$F$ と交わる ⇒ $C$ と交わる」
//   (3) 実際の格子での共有面の整合 — (2) の系ですが、**セル面が退化した箱
//       （lo == hi の軸を持つ）でも成り立つ**ことを直接確かめます
//
// **退化を明示的に構成します。** 乱択は「平面が隅をちょうど通る」「面を含む」を
// ほぼ生成しません。対象データでは常態です（CLAUDE.md）。
#include <array>
#include <cstdio>
#include <vector>

#include "krisite/geom/predicates.hpp"
#include "krisite/octree/uniform_grid.hpp"

#include "test_util.hpp"

using krisite::geom::Axis;
using krisite::geom::IPoint;
using krisite::geom::plane_axis_aligned;
using krisite::geom::plane_box_value;
using krisite::geom::plane_crosses_box;
using krisite::geom::plane_from_triangle;
using krisite::geom::PlaneD;
using kritest::Rng;

namespace {

struct Box {
    std::int64_t lo[3];
    std::int64_t hi[3];
};

/// 正解器: **8 隅すべて**を評価して符号の最小・最大を見る。
///
/// 被検体は「N の符号で lo/hi を選ぶ」ので、**選択の論理を共有していません。**
/// 算術そのものの正しさは GMP 差分テスト（`test_plane_box_gmp.cpp`）が受け持ちます。
bool crosses_by_corners(const PlaneD& pl, const Box& b) {
    int lo_sign = +1, hi_sign = -1;
    for (int m = 0; m < 8; ++m) {
        const std::int64_t p[3] = {
            (m & 1) ? b.hi[0] : b.lo[0],
            (m & 2) ? b.hi[1] : b.lo[1],
            (m & 4) ? b.hi[2] : b.lo[2],
        };
        const int s = krisite::arith::sign(plane_box_value(pl, p));
        if (s < lo_sign) lo_sign = s;
        if (s > hi_sign) hi_sign = s;
    }
    return lo_sign <= 0 && hi_sign >= 0;
}

std::int64_t rand_bound(Rng& rng) {
    // セル境界と同じ範囲 [-2^(b-1), +2^(b-1)]。**上限は kCoordMax を超えます**
    const std::int64_t span = -krisite::kCoordMin;
    return static_cast<std::int64_t>(rng.next() % static_cast<std::uint64_t>(2 * span + 1)) - span;
}

Box rand_box(Rng& rng) {
    Box b{};
    for (int t = 0; t < 3; ++t) {
        const std::int64_t u = rand_bound(rng), v = rand_bound(rng);
        b.lo[t] = u < v ? u : v;
        b.hi[t] = u < v ? v : u;
    }
    return b;
}

/// 退化を含む平面の一覧。**軸平行・斜面・法線成分が 0 のもの**を必ず混ぜます。
std::vector<PlaneD> degenerate_planes() {
    std::vector<PlaneD> ps;
    const std::int64_t s = -krisite::kCoordMin;  // 2^(b-1)
    for (Axis ax : {Axis::X, Axis::Y, Axis::Z}) {
        ps.push_back(plane_axis_aligned(ax, 0));
        ps.push_back(plane_axis_aligned(ax, s));   // 座標範囲の上端
        ps.push_back(plane_axis_aligned(ax, -s));  // 下端
        ps.push_back(plane_axis_aligned(ax, s / 2));
    }
    // 斜面。1 成分が 0 のもの（軸に平行だが軸平行でない）も入れる
    const std::int64_t q = s / 4;
    ps.push_back(
        plane_from_triangle(IPoint{0, 0, 0}, IPoint{static_cast<std::int32_t>(q), 0, 0},
                            IPoint{0, static_cast<std::int32_t>(q), static_cast<std::int32_t>(q)}));
    ps.push_back(plane_from_triangle(
        IPoint{0, 0, 0}, IPoint{static_cast<std::int32_t>(q), static_cast<std::int32_t>(q), 0},
        IPoint{0, 0, static_cast<std::int32_t>(q)}));
    ps.push_back(plane_from_triangle(
        IPoint{static_cast<std::int32_t>(-q), static_cast<std::int32_t>(-q), 0},
        IPoint{static_cast<std::int32_t>(q), 0, static_cast<std::int32_t>(q)},
        IPoint{0, static_cast<std::int32_t>(q), static_cast<std::int32_t>(-q)}));
    return ps;
}

// ---- (1) 8 隅の正解器との一致 -----------------------------------------------

void test_matches_corner_oracle() {
    Rng rng(4201);
    const std::vector<PlaneD> deg = degenerate_planes();
    std::size_t crossing = 0, total = 0;
    for (int iter = 0; iter < 40000; ++iter) {
        PlaneD pl;
        if (iter % 3 == 0 && !deg.empty()) {
            pl = deg[static_cast<std::size_t>(iter / 3) % deg.size()];
        } else {
            pl = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                     kritest::rand_point(rng));
        }
        if (krisite::geom::is_degenerate(pl)) continue;  // 退化三角形の平面は対象外
        const Box b = rand_box(rng);
        const bool got = plane_crosses_box(pl, b.lo, b.hi);
        KRI_CHECK_MSG(got == crosses_by_corners(pl, b), "8 隅の正解器と一致しない");
        ++total;
        if (got) ++crossing;
    }
    // **空回り防止**: 交差する場合としない場合の両方が出ていること
    std::printf("    8 隅との一致: %zu 件（うち交差 %zu = %.1f%%）\n", total, crossing,
                100.0 * static_cast<double>(crossing) / static_cast<double>(total));
    KRI_CHECK_MSG(crossing > total / 20, "交差する場合がほとんど出ていない（検査が空回り）");
    KRI_CHECK_MSG(crossing < total - total / 20, "交差しない場合がほとんど出ていない");
}

// ---- (2) 箱の包含に対する単調性 ★ §2.3 の整合性はここに帰着します -------------

void test_monotone_under_inclusion() {
    Rng rng(4202);
    const std::vector<PlaneD> deg = degenerate_planes();
    std::size_t inner_crossing = 0;
    for (int iter = 0; iter < 40000; ++iter) {
        PlaneD pl;
        if (iter % 3 == 0) {
            pl = deg[static_cast<std::size_t>(iter / 3) % deg.size()];
        } else {
            pl = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                     kritest::rand_point(rng));
        }
        if (krisite::geom::is_degenerate(pl)) continue;
        const Box outer = rand_box(rng);
        // outer に含まれる箱を作る（退化した箱 lo == hi も混ぜる）
        Box inner{};
        for (int t = 0; t < 3; ++t) {
            const std::uint64_t span = static_cast<std::uint64_t>(outer.hi[t] - outer.lo[t]) + 1;
            const std::int64_t u = outer.lo[t] + static_cast<std::int64_t>(rng.next() % span);
            const std::int64_t v = outer.lo[t] + static_cast<std::int64_t>(rng.next() % span);
            inner.lo[t] = u < v ? u : v;
            inner.hi[t] = u < v ? v : u;
            if (iter % 5 == 0) inner.hi[t] = inner.lo[t];  // 面・辺・点に退化させる
        }
        if (plane_crosses_box(pl, inner.lo, inner.hi)) {
            ++inner_crossing;
            KRI_CHECK_MSG(plane_crosses_box(pl, outer.lo, outer.hi),
                          "内側の箱と交わるのに外側と交わらない（単調性が破れた）");
        }
    }
    std::printf("    包含に対する単調性: 内側が交差した %zu 件で確認\n", inner_crossing);
    KRI_CHECK_MSG(inner_crossing > 1000, "内側が交差する場合が少なすぎる（検査が空回り）");
}

// ---- (3) 実際の格子での共有面の整合 ------------------------------------------
//
// §2.3 の証明そのものです。隣り合うセルの共有面 F は「lo == hi の軸を持つ箱」で、
// **F と交わる平面は必ず両方のセルと交わる**ことを直接確かめます。

void test_shared_face_consistency() {
    const std::vector<PlaneD> deg = degenerate_planes();
    Rng rng(4203);
    std::size_t checked = 0, face_crossing = 0;
    for (unsigned depth = 1; depth <= 3; ++depth) {
        const krisite::octree::UniformGrid g(depth);
        const std::uint32_t n = g.per_axis();
        for (int iter = 0; iter < 400; ++iter) {
            PlaneD pl;
            if (iter % 3 == 0) {
                pl = deg[static_cast<std::size_t>(iter / 3) % deg.size()];
            } else {
                pl = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                         kritest::rand_point(rng));
            }
            if (krisite::geom::is_degenerate(pl)) continue;
            for (std::uint32_t i = 0; i + 1 < n; ++i) {
                // X 方向に隣り合う 2 セル（他の軸は 0 に固定）
                const std::int64_t c1lo[3] = {g.lo(i), g.lo(0), g.lo(0)};
                const std::int64_t c1hi[3] = {g.hi(i), g.hi(0), g.hi(0)};
                const std::int64_t c2lo[3] = {g.lo(i + 1), g.lo(0), g.lo(0)};
                const std::int64_t c2hi[3] = {g.hi(i + 1), g.hi(0), g.hi(0)};
                // 共有面 F は x = g.hi(i) の退化した箱
                const std::int64_t flo[3] = {g.hi(i), g.lo(0), g.lo(0)};
                const std::int64_t fhi[3] = {g.hi(i), g.hi(0), g.hi(0)};
                ++checked;
                if (!plane_crosses_box(pl, flo, fhi)) continue;
                ++face_crossing;
                KRI_CHECK_MSG(plane_crosses_box(pl, c1lo, c1hi),
                              "共有面と交わるのに片側のセルと交わらない");
                KRI_CHECK_MSG(plane_crosses_box(pl, c2lo, c2hi),
                              "共有面と交わるのに反対側のセルと交わらない");
            }
        }
    }
    std::printf("    共有面の整合: %zu 組（うち面と交差 %zu）\n", checked, face_crossing);
    KRI_CHECK_MSG(face_crossing > 100, "共有面と交わる場合が少なすぎる（検査が空回り）");
}

// ---- (4) 手で構成した退化ケース ----------------------------------------------

void test_explicit_cases() {
    const std::int64_t s = -krisite::kCoordMin;
    const std::int64_t lo[3] = {0, 0, 0};
    const std::int64_t hi[3] = {s / 2, s / 2, s / 2};

    // 面をちょうど含む平面 → 閉領域なので交差
    KRI_CHECK_MSG(plane_crosses_box(plane_axis_aligned(Axis::X, 0), lo, hi),
                  "lo 面を含む平面が交差しない（閉領域で判定していない）");
    KRI_CHECK_MSG(plane_crosses_box(plane_axis_aligned(Axis::X, s / 2), lo, hi),
                  "hi 面を含む平面が交差しない");
    // 1 格子だけ外 → 交差しない
    KRI_CHECK_MSG(!plane_crosses_box(plane_axis_aligned(Axis::X, -1), lo, hi),
                  "lo より 1 小さい平面が交差してしまう");
    KRI_CHECK_MSG(!plane_crosses_box(plane_axis_aligned(Axis::X, s / 2 + 1), lo, hi),
                  "hi より 1 大きい平面が交差してしまう");

    // 隅をちょうど 1 点で通る斜面（x + y + z = 0 は原点 = lo 隅を通る）
    const PlaneD diag = plane_from_triangle(IPoint{0, 0, 0}, IPoint{1, -1, 0}, IPoint{0, 1, -1});
    KRI_CHECK_MSG(plane_crosses_box(diag, lo, hi), "隅を 1 点で通る平面が交差しない");

    // 退化した箱（1 点）でも動く
    const std::int64_t pt[3] = {0, 0, 0};
    KRI_CHECK_MSG(plane_crosses_box(plane_axis_aligned(Axis::X, 0), pt, pt),
                  "1 点の箱で、その点を通る平面が交差しない");
    KRI_CHECK_MSG(!plane_crosses_box(plane_axis_aligned(Axis::X, 1), pt, pt),
                  "1 点の箱で、通らない平面が交差してしまう");

    // 座標範囲の上端 +2^(b-1)（kCoordMax を超える。IPoint では表せない）
    const std::int64_t full_lo[3] = {-s, -s, -s};
    const std::int64_t full_hi[3] = {s, s, s};
    KRI_CHECK_MSG(plane_crosses_box(plane_axis_aligned(Axis::Z, s), full_lo, full_hi),
                  "座標範囲の上端に載る平面が交差しない");
}

}  // namespace

int main() {
    std::printf("\n  plane_crosses_box — SPEC-phase2 §2.3\n");
    test_matches_corner_oracle();
    test_monotone_under_inclusion();
    test_shared_face_consistency();
    test_explicit_cases();
    std::printf("\n");
    return kritest::finish("geom/plane_box");
}

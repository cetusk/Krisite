// Krisite — 退化ケースのテスト
//
// SPEC-phase0.md §8.2
//   ランダムテストは退化ケースをほとんど生成しない。以下を明示的に構成する。
//
//   | ケース                     | 生成方法                                     |
//   |----------------------------|----------------------------------------------|
//   | 共平面                     | 4 点を同一平面上に配置                       |
//   | 共線                       | 3 点を同一直線上に                           |
//   | 重複頂点                   | 同一座標の点を複数                           |
//   | 退化三角形                 | 面積 0                                       |
//   | 構成点が入力頂点に一致     | 3 平面の交点が既存頂点になるよう逆算         |
//   | 構成点同士の一致           | 異なる 3 平面組が同一点を与えるよう構成      |
//   | 格子の端                   | 座標が ±2^(b-1)                              |
//   | 極端に細長い三角形         | 1 点だけ大きくずらす                         |
#include <vector>

#include "test_util.hpp"

using namespace krisite::arith;
using namespace krisite::geom;
using kritest::Rng;

namespace {

constexpr Axis kAxes[3] = {Axis::X, Axis::Y, Axis::Z};

const std::int32_t kMin = static_cast<std::int32_t>(krisite::kCoordMin);
const std::int32_t kMax = static_cast<std::int32_t>(krisite::kCoordMax);

IPoint clampp(std::int64_t x, std::int64_t y, std::int64_t z) {
    auto c = [](std::int64_t v) {
        if (v < krisite::kCoordMin) v = krisite::kCoordMin;
        if (v > krisite::kCoordMax) v = krisite::kCoordMax;
        return static_cast<std::int32_t>(v);
    };
    return IPoint{c(x), c(y), c(z)};
}

// ---- 共平面 -----------------------------------------------------------------

void test_coplanar() {
    Rng rng(101);
    int built = 0;
    for (int iter = 0; iter < 20000 && built < 5000; ++iter) {
        const IPoint a = kritest::rand_small_point(rng, 1 << 10);
        const IPoint b = kritest::rand_small_point(rng, 1 << 10);
        const IPoint c = kritest::rand_small_point(rng, 1 << 10);
        if (orient2d(a, b, c, Axis::X) == 0 && orient2d(a, b, c, Axis::Y) == 0 &&
            orient2d(a, b, c, Axis::Z) == 0) {
            continue;  // 3 点が共線。平面が決まらない
        }
        // d = a + s*(b-a) + t*(c-a) は必ず abc の平面上
        const std::int64_t s = rng.range(-8, 8);
        const std::int64_t t = rng.range(-8, 8);
        const IPoint d =
            clampp(a.x + s * (std::int64_t{b.x} - a.x) + t * (std::int64_t{c.x} - a.x),
                   a.y + s * (std::int64_t{b.y} - a.y) + t * (std::int64_t{c.y} - a.y),
                   a.z + s * (std::int64_t{b.z} - a.z) + t * (std::int64_t{c.z} - a.z));
        // clamp が効いた場合は平面から外れるので、外れていないときだけ検証する
        const PlaneD pl = plane_from_triangle(a, b, c);
        if (side(pl, d) != 0) continue;
        ++built;
        KRI_CHECK(orient3d(a, b, c, d) == 0);
        KRI_CHECK(orient3d(b, a, c, d) == 0);
        KRI_CHECK(orient3d(d, a, b, c) == 0);
        KRI_CHECK(side(pl, to_homogeneous(d)) == 0);
    }
    KRI_CHECK(built > 500);
}

// ---- 共線 -------------------------------------------------------------------

void test_collinear() {
    Rng rng(103);
    for (int iter = 0; iter < 5000; ++iter) {
        const IPoint a = kritest::rand_small_point(rng, 1 << 12);
        const IPoint b = kritest::rand_small_point(rng, 1 << 12);
        const std::int64_t k = rng.range(-16, 16);
        const IPoint c =
            clampp(a.x + k * (std::int64_t{b.x} - a.x), a.y + k * (std::int64_t{b.y} - a.y),
                   a.z + k * (std::int64_t{b.z} - a.z));
        // clamp で直線から外れていないか確認
        const bool on_line = (std::int64_t{c.x} - a.x) == k * (std::int64_t{b.x} - a.x) &&
                             (std::int64_t{c.y} - a.y) == k * (std::int64_t{b.y} - a.y) &&
                             (std::int64_t{c.z} - a.z) == k * (std::int64_t{b.z} - a.z);
        if (!on_line) continue;

        // 共線 → 支持平面は退化（法線が零）
        const PlaneD pl = plane_from_triangle(a, b, c);
        KRI_CHECK(is_degenerate(pl));
        for (Axis ax : kAxes) KRI_CHECK(orient2d(a, b, c, ax) == 0);
        // 共線な 3 点はどんな 4 点目とも共平面
        const IPoint d = kritest::rand_point(rng);
        KRI_CHECK(orient3d(a, b, c, d) == 0);
    }
}

// ---- 重複頂点・退化三角形 ---------------------------------------------------

void test_duplicate_and_zero_area() {
    Rng rng(107);
    for (int iter = 0; iter < 20000; ++iter) {
        const IPoint a = kritest::rand_point(rng);
        const IPoint b = kritest::rand_point(rng);
        const IPoint c = kritest::rand_point(rng);

        // 重複頂点
        KRI_CHECK(orient3d(a, a, b, c) == 0);
        KRI_CHECK(orient3d(a, b, a, c) == 0);
        KRI_CHECK(orient3d(a, b, c, a) == 0);
        KRI_CHECK(orient3d(a, b, b, b) == 0);
        for (Axis ax : kAxes) {
            KRI_CHECK(orient2d(a, b, b, ax) == 0);
            KRI_CHECK(orient2d(a, a, a, ax) == 0);
        }

        // 面積 0 の三角形
        KRI_CHECK(is_degenerate(plane_from_triangle(a, a, b)));
        KRI_CHECK(is_degenerate(plane_from_triangle(a, b, a)));
        KRI_CHECK(is_degenerate(plane_from_triangle(a, a, a)));

        // 退化平面の d も 0（N = 0 なので N・p1 = 0）
        const PlaneD z = plane_from_triangle(a, a, b);
        KRI_CHECK(is_zero(z.d));
        KRI_CHECK(side(z, c) == 0);
    }
}

// ---- 構成点が入力頂点に一致 -------------------------------------------------

/// p を通る平面を作る（p を第 1 頂点にすれば必ず p を含む）。
PlaneD plane_through(Rng& rng, const IPoint& p) {
    return plane_from_triangle(p, kritest::rand_point(rng), kritest::rand_point(rng));
}

void test_construction_hits_input_vertex() {
    Rng rng(109);
    int built = 0;
    for (int iter = 0; iter < 40000 && built < 4000; ++iter) {
        const IPoint p =
            (iter % 3 == 0) ? kritest::rand_extreme_point(rng) : kritest::rand_point(rng);
        const PlaneD a = plane_through(rng, p);
        const PlaneD b = plane_through(rng, p);
        const PlaneD c = plane_through(rng, p);
        if (!kritest::intersects_at_point(a, b, c)) continue;
        ++built;

        const HPointD v = intersect3(a, b, c);
        const HPointD hp = to_homogeneous(p);
        for (Axis ax : kAxes) KRI_CHECK(cmp_h(v, hp, ax) == 0);
        KRI_CHECK(h_equal(v, hp));
        KRI_CHECK(!lex_less(v, hp));
        KRI_CHECK(!lex_less(hp, v));
        KRI_CHECK(side(a, v) == 0);
        KRI_CHECK(side(b, v) == 0);
        KRI_CHECK(side(c, v) == 0);
        // 入力点として見ても平面上
        KRI_CHECK(side(a, p) == 0);
    }
    KRI_CHECK(built > 500);
}

// ---- 構成点同士の一致 -------------------------------------------------------

void test_two_constructions_coincide() {
    Rng rng(113);
    int built = 0;
    for (int iter = 0; iter < 40000 && built < 3000; ++iter) {
        const IPoint p = kritest::rand_point(rng);
        const PlaneD a1 = plane_through(rng, p);
        const PlaneD a2 = plane_through(rng, p);
        const PlaneD a3 = plane_through(rng, p);
        const PlaneD b1 = plane_through(rng, p);
        const PlaneD b2 = plane_through(rng, p);
        const PlaneD b3 = plane_through(rng, p);
        if (!kritest::intersects_at_point(a1, a2, a3)) continue;
        if (!kritest::intersects_at_point(b1, b2, b3)) continue;
        ++built;

        const HPointD u = intersect3(a1, a2, a3);
        const HPointD v = intersect3(b1, b2, b3);
        KRI_CHECK(h_equal(u, v));
        for (Axis ax : kAxes) KRI_CHECK(cmp_h(u, v, ax) == 0);
        // 一方の構成に使った平面は他方の構成点も通る
        KRI_CHECK(side(a1, v) == 0);
        KRI_CHECK(side(b3, u) == 0);
    }
    KRI_CHECK(built > 300);
}

// ---- 格子の端 ---------------------------------------------------------------

void test_grid_extremes() {
    const std::int32_t e[2] = {kMin, kMax};
    std::vector<IPoint> corners;
    for (std::int32_t x : e)
        for (std::int32_t y : e)
            for (std::int32_t z : e) corners.push_back(IPoint{x, y, z});

    // 立方体の 8 隅で orient3d が符号を正しく返すこと（オーバーフローしないこと）
    for (std::size_t i = 0; i < corners.size(); ++i) {
        for (std::size_t j = 0; j < corners.size(); ++j) {
            for (std::size_t k = 0; k < corners.size(); ++k) {
                for (std::size_t l = 0; l < corners.size(); ++l) {
                    const int o = orient3d(corners[i], corners[j], corners[k], corners[l]);
                    KRI_CHECK(o == -orient3d(corners[j], corners[i], corners[k], corners[l]));
                    if (i == j || j == k || k == l || i == k || i == l || j == l) {
                        KRI_CHECK(o == 0);
                    }
                }
            }
        }
    }

    // 端の点を含む平面と構成点
    Rng rng(127);
    int built = 0;
    for (int iter = 0; iter < 40000 && built < 2000; ++iter) {
        const IPoint p = kritest::rand_extreme_point(rng);
        const PlaneD a = plane_from_triangle(p, kritest::rand_extreme_point(rng),
                                             kritest::rand_extreme_point(rng));
        const PlaneD b = plane_from_triangle(p, kritest::rand_extreme_point(rng),
                                             kritest::rand_extreme_point(rng));
        const PlaneD c = plane_from_triangle(p, kritest::rand_extreme_point(rng),
                                             kritest::rand_extreme_point(rng));
        if (!kritest::intersects_at_point(a, b, c)) continue;
        ++built;
        const HPointD v = intersect3(a, b, c);
        KRI_CHECK(side(a, v) == 0);
        KRI_CHECK(side(b, v) == 0);
        KRI_CHECK(side(c, v) == 0);
        KRI_CHECK(h_equal(v, to_homogeneous(p)));
    }
    KRI_CHECK(built > 200);
}

// ---- 極端に細長い三角形 -----------------------------------------------------

void test_sliver_triangles() {
    Rng rng(131);
    int built = 0;
    for (int iter = 0; iter < 20000 && built < 5000; ++iter) {
        // 2 点は隣接、1 点だけ格子の端まで飛ばす
        const IPoint a = kritest::rand_point(rng);
        const IPoint b =
            clampp(std::int64_t{a.x} + rng.range(-1, 1), std::int64_t{a.y} + rng.range(-1, 1),
                   std::int64_t{a.z} + rng.range(-1, 1));
        const IPoint c = kritest::rand_extreme_point(rng);
        const PlaneD pl = plane_from_triangle(a, b, c);
        if (is_degenerate(pl)) continue;
        ++built;
        KRI_CHECK(side(pl, a) == 0);
        KRI_CHECK(side(pl, b) == 0);
        KRI_CHECK(side(pl, c) == 0);
        KRI_CHECK(orient3d(a, b, c, a) == 0);
        // 平面上の点を線形結合で作っても 0
        const std::int64_t s = rng.range(-2, 2);
        const IPoint m =
            clampp(a.x + s * (std::int64_t{b.x} - a.x), a.y + s * (std::int64_t{b.y} - a.y),
                   a.z + s * (std::int64_t{b.z} - a.z));
        if ((std::int64_t{m.x} - a.x) == s * (std::int64_t{b.x} - a.x) &&
            (std::int64_t{m.y} - a.y) == s * (std::int64_t{b.y} - a.y) &&
            (std::int64_t{m.z} - a.z) == s * (std::int64_t{b.z} - a.z)) {
            KRI_CHECK(side(pl, m) == 0);
        }
    }
    KRI_CHECK(built > 1000);
}

}  // namespace

int main() {
    test_coplanar();
    test_collinear();
    test_duplicate_and_zero_area();
    test_construction_hits_input_vertex();
    test_two_constructions_coincide();
    test_grid_extremes();
    test_sliver_triangles();
    return kritest::finish("geom/degenerate");
}

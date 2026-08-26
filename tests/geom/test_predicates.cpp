// Krisite — 述語の性質テスト（正解器なし）
//
// SPEC-phase0.md §8.3
//   - orient3d(a,b,c,d) == -orient3d(b,a,c,d)
//   - 偶置換で符号不変、奇置換で反転（全 24 通り）
//   - side(plane, p) == 0 ⟺ p が plane 上
//   - cmp_h が全順序の公理（反対称律・推移律）を満たす
//   - 同次座標をスカラー倍しても述語の結果が不変
#include <array>
#include <vector>

#include "test_util.hpp"

using namespace krisite::arith;
using namespace krisite::geom;
using kritest::Rng;

namespace {

constexpr Axis kAxes[3] = {Axis::X, Axis::Y, Axis::Z};

// ---- orient3d の置換対称性 --------------------------------------------------

/// 置換の偶奇（+1 / -1）。
int parity(const std::array<int, 4>& p) {
    int inv = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j)
            if (p[i] > p[j]) ++inv;
    return (inv % 2 == 0) ? 1 : -1;
}

void test_orient3d_permutations() {
    Rng rng(11);
    std::array<std::array<int, 4>, 24> perms{};
    {
        std::array<int, 4> p{0, 1, 2, 3};
        int n = 0;
        // 辞書順に 24 通り
        auto next_perm = [](std::array<int, 4>& a) {
            int i = 2;
            while (i >= 0 && a[i] >= a[i + 1]) --i;
            if (i < 0) return false;
            int j = 3;
            while (a[j] <= a[i]) --j;
            std::swap(a[i], a[j]);
            for (int l = i + 1, r = 3; l < r; ++l, --r) std::swap(a[l], a[r]);
            return true;
        };
        do {
            perms[static_cast<std::size_t>(n++)] = p;
        } while (next_perm(p) && n < 24);
        KRI_CHECK(n == 24);
    }

    for (int iter = 0; iter < 5000; ++iter) {
        IPoint q[4];
        // 半分は退化しやすい小さな座標で
        for (int i = 0; i < 4; ++i)
            q[i] = (iter % 2 == 0) ? kritest::rand_point(rng) : kritest::rand_small_point(rng, 3);

        const int base = orient3d(q[0], q[1], q[2], q[3]);
        for (const auto& p : perms) {
            const int got =
                orient3d(q[static_cast<std::size_t>(p[0])], q[static_cast<std::size_t>(p[1])],
                         q[static_cast<std::size_t>(p[2])], q[static_cast<std::size_t>(p[3])]);
            KRI_CHECK(got == parity(p) * base);
        }
        // 明示的に 1 回入れ替え
        KRI_CHECK(orient3d(q[0], q[1], q[2], q[3]) == -orient3d(q[1], q[0], q[2], q[3]));
    }
}

void test_orient2d_permutations() {
    Rng rng(12);
    for (int iter = 0; iter < 20000; ++iter) {
        const IPoint a = kritest::rand_point(rng);
        const IPoint b = kritest::rand_point(rng);
        const IPoint c = kritest::rand_point(rng);
        for (Axis ax : kAxes) {
            const int base = orient2d(a, b, c, ax);
            KRI_CHECK(orient2d(b, a, c, ax) == -base);  // 奇置換
            KRI_CHECK(orient2d(b, c, a, ax) == base);   // 偶置換
            KRI_CHECK(orient2d(c, a, b, ax) == base);   // 偶置換
            KRI_CHECK(orient2d(a, a, c, ax) == 0);      // 重複点
        }
    }
}

// ---- orient3d と side の整合 ------------------------------------------------

void test_orient3d_matches_side() {
    Rng rng(13);
    int checked = 0;
    for (int iter = 0; iter < 20000 && checked < 10000; ++iter) {
        const IPoint a = kritest::rand_point(rng);
        const IPoint b = kritest::rand_point(rng);
        const IPoint c = kritest::rand_point(rng);
        const IPoint d = kritest::rand_point(rng);
        const PlaneD pl = plane_from_triangle(a, b, c);
        if (is_degenerate(pl)) continue;
        ++checked;
        // det[a-d, b-d, c-d] = -N・(d-a) （N = (b-a)x(c-a)）
        KRI_CHECK(orient3d(a, b, c, d) == -side(pl, d));
        // 支持平面上の 3 点は必ず 0
        KRI_CHECK(side(pl, a) == 0);
        KRI_CHECK(side(pl, b) == 0);
        KRI_CHECK(side(pl, c) == 0);
        // 同次点として見ても同じ（w = 1 > 0）
        KRI_CHECK(side(pl, to_homogeneous(d)) == side(pl, d));
        KRI_CHECK(side(pl, to_homogeneous(a)) == 0);
    }
    KRI_CHECK(checked > 1000);
}

// ---- intersect3 の不変条件 --------------------------------------------------

/// 3 平面の交点は、その 3 平面すべての上にある。
void test_intersect3_on_planes() {
    Rng rng(17);
    int checked = 0;
    for (int iter = 0; iter < 40000 && checked < 20000; ++iter) {
        const PlaneD p1 = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                              kritest::rand_point(rng));
        const PlaneD p2 = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                              kritest::rand_point(rng));
        const PlaneD p3 = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                              kritest::rand_point(rng));
        if (!kritest::intersects_at_point(p1, p2, p3)) continue;
        ++checked;
        const HPointD v = intersect3(p1, p2, p3);
        KRI_CHECK(side(p1, v) == 0);
        KRI_CHECK(side(p2, v) == 0);
        KRI_CHECK(side(p3, v) == 0);
        // 平面の順序を変えても同じ点（同次座標としては符号が変わりうる）
        const HPointD v2 = intersect3(p2, p1, p3);
        KRI_CHECK(cmp_h(v, v2, Axis::X) == 0);
        KRI_CHECK(cmp_h(v, v2, Axis::Y) == 0);
        KRI_CHECK(cmp_h(v, v2, Axis::Z) == 0);
    }
    KRI_CHECK(checked > 1000);
}

// ---- 同次座標のスカラー倍不変性 ---------------------------------------------

HPointD scale(const HPointD& h, std::int64_t k) {
    const auto kk = from_i64<1>(k);
    HPointD r{};
    r.x = resize<limbs::kHomoXyz>(mul(h.x, kk));
    r.y = resize<limbs::kHomoXyz>(mul(h.y, kk));
    r.z = resize<limbs::kHomoXyz>(mul(h.z, kk));
    r.w = resize<limbs::kHomoW>(mul(h.w, kk));
    return r;
}

void test_scale_invariance() {
    Rng rng(19);
    int checked = 0;
    for (int iter = 0; iter < 40000 && checked < 8000; ++iter) {
        const PlaneD p1 = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                              kritest::rand_point(rng));
        const PlaneD p2 = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                              kritest::rand_point(rng));
        const PlaneD p3 = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                              kritest::rand_point(rng));
        const PlaneD probe = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                                 kritest::rand_point(rng));
        if (!kritest::intersects_at_point(p1, p2, p3)) continue;
        ++checked;
        const HPointD v = intersect3(p1, p2, p3);

        // |k| <= 1000 なら kHomoXyz の余裕（b=21 で 161 → 192 ビット）に収まる
        std::int64_t k = rng.range(1, 1000);
        if (rng.below(2)) k = -k;

        const HPointD sv = scale(v, k);
        KRI_CHECK(side(probe, sv) == side(probe, v));
        for (Axis ax : kAxes) {
            KRI_CHECK(cmp_h(sv, v, ax) == 0);
        }
        // 別の点との比較順序もスカラー倍で変わらない
        const HPointD other = to_homogeneous(kritest::rand_point(rng));
        for (Axis ax : kAxes) {
            KRI_CHECK(cmp_h(sv, other, ax) == cmp_h(v, other, ax));
            KRI_CHECK(cmp_h(other, sv, ax) == cmp_h(other, v, ax));
        }
    }
    KRI_CHECK(checked > 500);
}

// ---- cmp_h の全順序性 -------------------------------------------------------

void test_cmp_h_total_order() {
    Rng rng(23);
    std::vector<HPointD> pts;
    pts.reserve(64);
    while (pts.size() < 40) {
        const PlaneD p1 = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                              kritest::rand_point(rng));
        const PlaneD p2 = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                              kritest::rand_point(rng));
        const PlaneD p3 = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                              kritest::rand_point(rng));
        if (!kritest::intersects_at_point(p1, p2, p3)) continue;
        pts.push_back(intersect3(p1, p2, p3));
        pts.push_back(to_homogeneous(kritest::rand_small_point(rng, 5)));
        pts.push_back(scale(pts[pts.size() - 2], -3));  // 同じ点の別表現
    }

    for (std::size_t i = 0; i < pts.size(); ++i) {
        // 反射律
        KRI_CHECK(cmp_h_lex(pts[i], pts[i]) == 0);
        for (std::size_t j = 0; j < pts.size(); ++j) {
            // 反対称律
            KRI_CHECK(cmp_h_lex(pts[i], pts[j]) == -cmp_h_lex(pts[j], pts[i]));
            for (Axis ax : kAxes) {
                KRI_CHECK(cmp_h(pts[i], pts[j], ax) == -cmp_h(pts[j], pts[i], ax));
            }
            for (std::size_t k = 0; k < pts.size(); ++k) {
                // 推移律
                if (cmp_h_lex(pts[i], pts[j]) <= 0 && cmp_h_lex(pts[j], pts[k]) <= 0) {
                    KRI_CHECK(cmp_h_lex(pts[i], pts[k]) <= 0);
                }
            }
        }
    }

    // lex_less が狭義全順序であること
    for (std::size_t i = 0; i < pts.size(); ++i) {
        KRI_CHECK(!lex_less(pts[i], pts[i]));
        for (std::size_t j = 0; j < pts.size(); ++j) {
            const bool lt = lex_less(pts[i], pts[j]);
            const bool gt = lex_less(pts[j], pts[i]);
            const bool eq = h_equal(pts[i], pts[j]);
            KRI_CHECK(static_cast<int>(lt) + static_cast<int>(gt) + static_cast<int>(eq) == 1);
        }
    }
}

/// 整数点同士の cmp_h は素朴な比較に一致する。
void test_cmp_h_against_integers() {
    Rng rng(29);
    for (int iter = 0; iter < 50000; ++iter) {
        const IPoint a = kritest::rand_point(rng);
        const IPoint b = kritest::rand_point(rng);
        const HPointD ha = to_homogeneous(a);
        const HPointD hb = to_homogeneous(b);
        auto sgn = [](std::int64_t v) { return (v < 0) ? -1 : (v > 0) ? 1 : 0; };
        KRI_CHECK(cmp_h(ha, hb, Axis::X) == sgn(std::int64_t{a.x} - b.x));
        KRI_CHECK(cmp_h(ha, hb, Axis::Y) == sgn(std::int64_t{a.y} - b.y));
        KRI_CHECK(cmp_h(ha, hb, Axis::Z) == sgn(std::int64_t{a.z} - b.z));

        const bool want_less = (a.x != b.x)   ? (a.x < b.x)
                               : (a.y != b.y) ? (a.y < b.y)
                                              : (a.z < b.z);
        KRI_CHECK(lex_less(ha, hb) == want_less);
    }
}

}  // namespace

int main() {
    test_orient3d_permutations();
    test_orient2d_permutations();
    test_orient3d_matches_side();
    test_intersect3_on_planes();
    test_scale_invariance();
    test_cmp_h_total_order();
    test_cmp_h_against_integers();
    return kritest::finish("geom/predicates");
}

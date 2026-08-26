// Krisite — 平面の同一判定と全順序、軸平行平面の性質テスト
//
// SPEC-phase1.md §3.1（平面 ID）, §3.2（セル面も平面である）, §3.4（向き付けの規約）
//
// plane_same と plane_cmp は別々の式から導いてあるので、**両者が一致すること**が
// 中心的な検査になります。片方だけ間違えれば必ず食い違います。
#include <array>
#include <utility>
#include <vector>

#include "test_util.hpp"

using namespace krisite::arith;
using namespace krisite::geom;
using kritest::Rng;

namespace {

constexpr Axis kAxes[3] = {Axis::X, Axis::Y, Axis::Z};

/// 平面をスカラー倍する（同一の幾何平面を別の係数で表す）。
/// 係数が kNormal / kOffset に収まる範囲でしか倍率を掛けられないので k は小さく。
PlaneD scale_plane(const PlaneD& pl, std::int64_t k) {
    const auto kk = from_i64<1>(k);
    PlaneD r{};
    r.a = resize<limbs::kNormal>(mul(pl.a, kk));
    r.b = resize<limbs::kNormal>(mul(pl.b, kk));
    r.c = resize<limbs::kNormal>(mul(pl.c, kk));
    r.d = resize<limbs::kOffset>(mul(pl.d, kk));
    return r;
}

/// 法線が小さい平面（スカラー倍しても幅に収まるように）。
PlaneD rand_small_plane(Rng& rng) {
    for (;;) {
        const PlaneD pl = plane_from_triangle(kritest::rand_small_point(rng, 64),
                                              kritest::rand_small_point(rng, 64),
                                              kritest::rand_small_point(rng, 64));
        if (!is_degenerate(pl)) return pl;
    }
}

// ---- plane_same と plane_cmp の一致 ------------------------------------------

void test_same_matches_cmp() {
    Rng rng(2026);
    for (int iter = 0; iter < 40000; ++iter) {
        const PlaneD p = rand_small_plane(rng);
        const PlaneD q = (iter % 3 == 0) ? p : rand_small_plane(rng);
        KRI_CHECK(plane_same(p, q) == (plane_cmp(p, q) == 0));
        KRI_CHECK(plane_same(q, p) == (plane_cmp(q, p) == 0));
        // 対称性
        KRI_CHECK(plane_same(p, q) == plane_same(q, p));
        KRI_CHECK(plane_cmp(p, q) == -plane_cmp(q, p));
    }
}

// ---- スカラー倍不変性 --------------------------------------------------------

void test_scale_invariance() {
    Rng rng(2027);
    int checked = 0;
    for (int iter = 0; iter < 40000 && checked < 20000; ++iter) {
        const PlaneD p = rand_small_plane(rng);
        std::int64_t k = rng.range(1, 64);
        if (rng.below(2)) k = -k;  // 符号反転 = 平面の向きを裏返す
        const PlaneD sp = scale_plane(p, k);
        ++checked;

        // 同一の幾何平面。向きが反転していても同一と判定されること
        KRI_CHECK(plane_same(p, sp));
        KRI_CHECK(plane_cmp(p, sp) == 0);

        // 他の平面との順序もスカラー倍で変わらない
        const PlaneD q = rand_small_plane(rng);
        KRI_CHECK(plane_cmp(sp, q) == plane_cmp(p, q));
        KRI_CHECK(plane_same(sp, q) == plane_same(p, q));
    }
    KRI_CHECK(checked > 1000);
}

// ---- 全順序の公理 ------------------------------------------------------------

void test_total_order() {
    Rng rng(2028);
    std::vector<PlaneD> ps;
    ps.reserve(60);
    for (int i = 0; i < 20; ++i) {
        const PlaneD p = rand_small_plane(rng);
        ps.push_back(p);
        ps.push_back(scale_plane(p, 3));   // 同じ平面の別表現
        ps.push_back(scale_plane(p, -2));  // 向きも裏返した別表現
    }
    // 軸平行平面も混ぜる
    for (Axis ax : kAxes) {
        ps.push_back(plane_axis_aligned(ax, 0));
        ps.push_back(plane_axis_aligned(ax, 1 << 10));
    }

    for (std::size_t i = 0; i < ps.size(); ++i) {
        KRI_CHECK(plane_cmp(ps[i], ps[i]) == 0);  // 反射律
        for (std::size_t j = 0; j < ps.size(); ++j) {
            KRI_CHECK(plane_cmp(ps[i], ps[j]) == -plane_cmp(ps[j], ps[i]));  // 反対称律
            for (std::size_t k = 0; k < ps.size(); ++k) {
                if (plane_cmp(ps[i], ps[j]) <= 0 && plane_cmp(ps[j], ps[k]) <= 0) {
                    KRI_CHECK(plane_cmp(ps[i], ps[k]) <= 0);  // 推移律
                }
            }
        }
    }

    // 同値類が plane_same と一致すること
    for (std::size_t i = 0; i < ps.size(); ++i) {
        for (std::size_t j = 0; j < ps.size(); ++j) {
            KRI_CHECK((plane_cmp(ps[i], ps[j]) == 0) == plane_same(ps[i], ps[j]));
        }
    }
}

// ---- 退化平面 ----------------------------------------------------------------

void test_null_plane() {
    Rng rng(2029);
    const IPoint a = kritest::rand_point(rng);
    const PlaneD null_pl = plane_from_triangle(a, a, a);  // 面積 0 → 全成分 0
    KRI_CHECK(is_null(null_pl));
    KRI_CHECK(plane_same(null_pl, null_pl));
    KRI_CHECK(plane_cmp(null_pl, null_pl) == 0);

    // 退化平面は非退化平面と同一ではない（零ベクトルは形式的にはあらゆる平面と比例する）
    for (int iter = 0; iter < 2000; ++iter) {
        const PlaneD p = rand_small_plane(rng);
        KRI_CHECK(!is_null(p));
        KRI_CHECK(!plane_same(null_pl, p));
        KRI_CHECK(!plane_same(p, null_pl));
        KRI_CHECK(plane_cmp(null_pl, p) != 0);
        // plane_same と plane_cmp がここでも一致すること
        KRI_CHECK(plane_same(null_pl, p) == (plane_cmp(null_pl, p) == 0));
    }
}

// ---- 軸平行平面 --------------------------------------------------------------

void test_axis_aligned() {
    Rng rng(2030);
    for (int iter = 0; iter < 20000; ++iter) {
        const std::int64_t c = rng.range(krisite::kCoordMin, -krisite::kCoordMin);
        const Axis ax = kAxes[rng.below(3)];
        const PlaneD pl = plane_axis_aligned(ax, c);
        KRI_CHECK(!is_degenerate(pl));

        // 面上・両側の点で side が期待どおりか（座標範囲内に収まるときだけ）
        const IPoint q = kritest::rand_point(rng);
        const std::int64_t qc = (ax == Axis::X) ? q.x : (ax == Axis::Y) ? q.y : q.z;
        const int want = (qc > c) ? 1 : (qc < c) ? -1 : 0;
        KRI_CHECK(side(pl, q) == want);
        KRI_CHECK(side(pl, to_homogeneous(q)) == want);
    }

    // セル境界の最大値 +2^(b-1) は IPoint の範囲外だが平面としては作れること
    const PlaneD hi = plane_axis_aligned(Axis::X, -krisite::kCoordMin);
    KRI_CHECK(!is_degenerate(hi));
    // 座標範囲内の全点はこの平面より小さい側にある
    Rng rng2(2031);
    for (int i = 0; i < 5000; ++i) {
        KRI_CHECK(side(hi, kritest::rand_extreme_point(rng2)) == -1);
    }

    // 同じ座標なら同一平面、違えば別平面
    for (Axis ax : kAxes) {
        KRI_CHECK(plane_same(plane_axis_aligned(ax, 100), plane_axis_aligned(ax, 100)));
        KRI_CHECK(!plane_same(plane_axis_aligned(ax, 100), plane_axis_aligned(ax, 101)));
    }
    // 軸が違えば別平面
    KRI_CHECK(!plane_same(plane_axis_aligned(Axis::X, 0), plane_axis_aligned(Axis::Y, 0)));
}

/// 立方体の面から作った平面が、対応する軸平行平面と同一と判定されること。
/// ケース 5（面がセル境界と完全一致）で効く性質。
void test_cube_face_matches_grid_plane() {
    const std::int32_t lo = -64, hi = 64;
    // z = hi の面（外向き +Z）を 2 三角形で
    const IPoint p0{lo, lo, hi}, p1{hi, lo, hi}, p2{hi, hi, hi}, p3{lo, hi, hi};
    const PlaneD t0 = plane_from_triangle(p0, p1, p2);
    const PlaneD t1 = plane_from_triangle(p0, p2, p3);
    const PlaneD grid = plane_axis_aligned(Axis::Z, hi);

    KRI_CHECK(plane_same(t0, t1));    // 1 面 2 三角形は同一平面
    KRI_CHECK(plane_same(t0, grid));  // 格子平面とも同一
    KRI_CHECK(plane_cmp(t0, grid) == 0);
    KRI_CHECK(side(t0, p0) == 0);
    KRI_CHECK(side(grid, p0) == 0);
}

// ---- 符号付き体積（向き検査）------------------------------------------------

/// 外向き法線の立方体は正の体積を返す。裏返せば負。
void test_signed_volume() {
    const std::int32_t s = 100;
    const IPoint v[8] = {{0, 0, 0}, {s, 0, 0}, {s, s, 0}, {0, s, 0},
                         {0, 0, s}, {s, 0, s}, {s, s, s}, {0, s, s}};
    // 外向き CCW の 12 三角形
    static const int tri[12][3] = {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
                                   {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                                   {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    constexpr std::size_t LV = limbs::kInputVolume6;
    auto acc = zero<LV>();
    for (const auto& t : tri) {
        acc = add(acc, resize<LV>(tetra_volume6(v[t[0]], v[t[1]], v[t[2]])));
    }
    // 立方体の体積 x6 = 6 * s^3
    const auto want = from_i64<LV>(6ll * s * s * s);
    KRI_CHECK(cmp(acc, want) == 0);
    KRI_CHECK(sign(acc) == 1);

    // 向きを裏返すと符号が反転する
    auto rev = zero<LV>();
    for (const auto& t : tri) {
        rev = add(rev, resize<LV>(tetra_volume6(v[t[0]], v[t[2]], v[t[1]])));
    }
    KRI_CHECK(sign(rev) == -1);
    KRI_CHECK(cmp(rev, neg(want)) == 0);

    // 平行移動で体積は変わらない
    Rng rng(2032);
    for (int iter = 0; iter < 200; ++iter) {
        const std::int32_t dx = static_cast<std::int32_t>(rng.range(-1000, 1000));
        const std::int32_t dy = static_cast<std::int32_t>(rng.range(-1000, 1000));
        const std::int32_t dz = static_cast<std::int32_t>(rng.range(-1000, 1000));
        auto m = zero<LV>();
        for (const auto& t : tri) {
            IPoint q[3];
            for (int k = 0; k < 3; ++k) {
                q[k] = IPoint{v[t[k]].x + dx, v[t[k]].y + dy, v[t[k]].z + dz};
            }
            m = add(m, resize<LV>(tetra_volume6(q[0], q[1], q[2])));
        }
        KRI_CHECK(cmp(m, want) == 0);
    }
}

}  // namespace

int main() {
    test_same_matches_cmp();
    test_scale_invariance();
    test_total_order();
    test_null_plane();
    test_axis_aligned();
    test_cube_face_matches_grid_plane();
    test_signed_volume();
    return kritest::finish("geom/plane_order");
}

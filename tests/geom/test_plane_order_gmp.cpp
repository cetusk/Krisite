// Krisite — 平面の同一判定・全順序・符号付き体積の GMP 差分テスト
//
// SPEC-phase1.md §3.1, §3.2, §3.4
//
// 正解器は被検体と**別経路**で書きます（`docs/ROADMAP.md`「通しての約束」）。
//   - plane_cmp は交差乗算で比較する。正解器は mpq で先頭成分を割って比を作り、直接比較する
//   - plane_same は 6 本の小行列式。正解器は mpz で平面を作り直してから同じ 6 本を計算する
//     （式は同じだが、係数は入力点から独立に再構成したもの）
//
// GMP は LGPL。テスト専用。
#include <cstdlib>

#include "gmp_oracle.hpp"
#include "test_util.hpp"

using namespace krisite::arith;
using namespace krisite::geom;
using kritest::Rng;
using kritest::oracle::to_mpz;

namespace {

constexpr Axis kAxes[3] = {Axis::X, Axis::Y, Axis::Z};

/// 平面 1 枚ぶんの mpz 係数。
struct PlaneZ {
    mpz_t c[4];  // a, b, c, d
    PlaneZ() {
        for (auto& x : c) mpz_init(x);
    }
    ~PlaneZ() {
        for (auto& x : c) mpz_clear(x);
    }
    PlaneZ(const PlaneZ&) = delete;
    PlaneZ& operator=(const PlaneZ&) = delete;
};

/// 入力点から平面を mpz で作り直す（ライブラリの値は使わない）。
void build_plane(PlaneZ& out, const IPoint& p1, const IPoint& p2, const IPoint& p3, mpz_t t0,
                 mpz_t t1) {
    const long u[3] = {static_cast<long>(p2.x) - p1.x, static_cast<long>(p2.y) - p1.y,
                       static_cast<long>(p2.z) - p1.z};
    const long v[3] = {static_cast<long>(p3.x) - p1.x, static_cast<long>(p3.y) - p1.y,
                       static_cast<long>(p3.z) - p1.z};
    for (int i = 0; i < 3; ++i) {
        const int a = (i + 1) % 3, b = (i + 2) % 3;
        mpz_set_si(t0, u[a]);
        mpz_mul_si(t0, t0, v[b]);
        mpz_set_si(t1, u[b]);
        mpz_mul_si(t1, t1, v[a]);
        mpz_sub(out.c[i], t0, t1);
    }
    const long q[3] = {p1.x, p1.y, p1.z};
    mpz_set_ui(out.c[3], 0);
    for (int i = 0; i < 3; ++i) {
        mpz_mul_si(t0, out.c[i], q[i]);
        mpz_add(out.c[3], out.c[3], t0);
    }
    mpz_neg(out.c[3], out.c[3]);  // d = -N・p1（SPEC-phase0 §3.1）
}

bool all_zero(const PlaneZ& p) {
    for (const auto& x : p.c) {
        if (mpz_sgn(x) != 0) return false;
    }
    return true;
}

int lead_of(const PlaneZ& p) {
    for (int i = 0; i < 4; ++i) {
        if (mpz_sgn(p.c[i]) != 0) return i;
    }
    return 4;
}

/// 正解器: 6 本の小行列式で比例判定。
bool same_oracle(const PlaneZ& p, const PlaneZ& q, mpz_t t0, mpz_t t1) {
    const bool pz = all_zero(p), qz = all_zero(q);
    if (pz || qz) return pz && qz;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            mpz_mul(t0, p.c[i], q.c[j]);
            mpz_mul(t1, p.c[j], q.c[i]);
            mpz_sub(t0, t0, t1);
            if (mpz_sgn(t0) != 0) return false;
        }
    }
    return true;
}

/// 正解器: mpq で先頭成分を割って比を作り、直接比較する（交差乗算を使わない）。
int cmp_oracle(const PlaneZ& p, const PlaneZ& q, mpq_t rp, mpq_t rq) {
    const int ip = lead_of(p), iq = lead_of(q);
    if (ip != iq) return (ip < iq) ? -1 : 1;
    if (ip == 4) return 0;
    for (int j = ip + 1; j < 4; ++j) {
        mpq_set_num(rp, p.c[j]);
        mpq_set_den(rp, p.c[ip]);
        mpq_canonicalize(rp);
        mpq_set_num(rq, q.c[j]);
        mpq_set_den(rq, q.c[iq]);
        mpq_canonicalize(rq);
        const int c = mpq_cmp(rp, rq);
        if (c != 0) return (c < 0) ? -1 : 1;
    }
    return 0;
}

int sgn_of(int v) {
    return (v < 0) ? -1 : (v > 0) ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    long iters = (argc > 1) ? std::atol(argv[1]) : 300000L;
    if (iters < 100) iters = 100;

    Rng rng(0xB1ADEu);
    PlaneZ zp, zq;
    mpz_t t0, t1, g;
    mpz_init(t0);
    mpz_init(t1);
    mpz_init(g);
    mpq_t rp, rq;
    mpq_init(rp);
    mpq_init(rq);

    for (long it = 0; it < iters; ++it) {
        // 小さめの座標。スカラー倍しても幅に収まるように
        const bool small = (it % 2 == 0);
        auto pick = [&] {
            return small ? kritest::rand_small_point(rng, 64) : kritest::rand_point(rng);
        };
        IPoint a[3] = {pick(), pick(), pick()};
        IPoint b[3] = {pick(), pick(), pick()};
        // 1/4 は同じ三角形（同一平面の経路を厚く踏む）
        if (it % 4 == 0) {
            b[0] = a[0];
            b[1] = a[1];
            b[2] = a[2];
        }

        const PlaneD p = plane_from_triangle(a[0], a[1], a[2]);
        const PlaneD q = plane_from_triangle(b[0], b[1], b[2]);
        build_plane(zp, a[0], a[1], a[2], t0, t1);
        build_plane(zq, b[0], b[1], b[2], t0, t1);

        // 係数そのものが一致すること（平面構成の差分検査）
        for (int i = 0; i < 4; ++i) {
            mpz_t got;
            mpz_init(got);
            switch (i) {
                case 0:
                    to_mpz(got, p.a);
                    break;
                case 1:
                    to_mpz(got, p.b);
                    break;
                case 2:
                    to_mpz(got, p.c);
                    break;
                default:
                    to_mpz(got, p.d);
                    break;
            }
            KRI_CHECK(mpz_cmp(got, zp.c[i]) == 0);
            mpz_clear(got);
        }

        // 小行列式
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                mpz_mul(t0, zp.c[i], zq.c[j]);
                mpz_mul(t1, zp.c[j], zq.c[i]);
                mpz_sub(t0, t0, t1);
                mpz_t got;
                mpz_init(got);
                to_mpz(got, plane_minor(p, q, i, j));
                KRI_CHECK(mpz_cmp(got, t0) == 0);
                mpz_clear(got);
            }
        }

        // 比例判定と全順序
        KRI_CHECK(plane_same(p, q) == same_oracle(zp, zq, t0, t1));
        KRI_CHECK(plane_cmp(p, q) == cmp_oracle(zp, zq, rp, rq));
        KRI_CHECK(plane_cmp(q, p) == cmp_oracle(zq, zp, rp, rq));
        // 2 つの定義が一致すること（別経路どうしの突き合わせ）
        KRI_CHECK(same_oracle(zp, zq, t0, t1) == (cmp_oracle(zp, zq, rp, rq) == 0));

        // 符号付き体積 x6
        {
            mpz_t m[3][3], det, cof, s0, s1;
            for (auto& row : m) {
                for (auto& e : row) mpz_init(e);
            }
            mpz_init(det);
            mpz_init(cof);
            mpz_init(s0);
            mpz_init(s1);
            const IPoint* tri[3] = {&a[0], &a[1], &a[2]};
            for (int r = 0; r < 3; ++r) {
                mpz_set_si(m[r][0], tri[r]->x);
                mpz_set_si(m[r][1], tri[r]->y);
                mpz_set_si(m[r][2], tri[r]->z);
            }
            mpz_set_ui(det, 0);
            for (int j = 0; j < 3; ++j) {
                const int j1 = (j + 1) % 3, j2 = (j + 2) % 3;
                mpz_mul(s0, m[1][j1], m[2][j2]);
                mpz_mul(s1, m[1][j2], m[2][j1]);
                mpz_sub(cof, s0, s1);
                mpz_mul(cof, cof, m[0][j]);
                mpz_add(det, det, cof);
            }
            mpz_t got;
            mpz_init(got);
            to_mpz(got, tetra_volume6(a[0], a[1], a[2]));
            KRI_CHECK(mpz_cmp(got, det) == 0);
            mpz_clear(got);
            for (auto& row : m) {
                for (auto& e : row) mpz_clear(e);
            }
            mpz_clear(det);
            mpz_clear(cof);
            mpz_clear(s0);
            mpz_clear(s1);
        }

        // 軸平行平面: 係数が (e_axis, -coord) であること
        {
            const Axis ax = kAxes[rng.below(3)];
            const std::int64_t cc = rng.range(krisite::kCoordMin, -krisite::kCoordMin);
            const PlaneD gp = plane_axis_aligned(ax, cc);
            mpz_t got;
            mpz_init(got);
            for (int i = 0; i < 3; ++i) {
                switch (i) {
                    case 0:
                        to_mpz(got, gp.a);
                        break;
                    case 1:
                        to_mpz(got, gp.b);
                        break;
                    default:
                        to_mpz(got, gp.c);
                        break;
                }
                KRI_CHECK(mpz_cmp_si(got, (static_cast<int>(ax) == i) ? 1 : 0) == 0);
            }
            to_mpz(got, gp.d);
            mpz_set_si(t0, static_cast<long>(-cc));
            KRI_CHECK(mpz_cmp(got, t0) == 0);
            mpz_clear(got);
            // 有理点での side が符号どおりか
            const IPoint probe = kritest::rand_point(rng);
            const std::int64_t pc = (ax == Axis::X) ? probe.x : (ax == Axis::Y) ? probe.y : probe.z;
            KRI_CHECK(side(gp, probe) == sgn_of(static_cast<int>((pc > cc) - (pc < cc))));
        }
    }

    mpq_clear(rp);
    mpq_clear(rq);
    mpz_clear(t0);
    mpz_clear(t1);
    mpz_clear(g);
    std::printf("       %ld 件\n", iters);
    return kritest::finish("geom/plane_order_gmp");
}

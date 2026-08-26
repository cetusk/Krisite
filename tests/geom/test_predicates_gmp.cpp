// Krisite — 述語の GMP 差分テスト
//
// SPEC-phase0.md §8.1-2
//   ランダムな IPoint から平面と構成点を作り、各述語の結果を GMP の有理数演算と比較。
//
// 独立性のために、side / cmp_h は **有理数（mpq）で実座標を復元してから** 評価する。
// 固定幅側と同じ「sign(w) を掛ける」変形を使わないので、変形そのものの誤りも検出できる。
//
// 平面の約束は SPEC §3.1 に従い N・x + d = 0（d = -N・p1）。正解器も同じ流儀で書く。
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

/// 作業用の mpz / mpq をループ外に確保しておく。
struct Work {
    mpz_t u[3], v[3], n[3], d, t0, t1, t2;
    mpz_t m[3][3], det, cof;
    mpz_t px[3], pw;
    mpz_t na[3][3], da[3];
    mpq_t q[3], qa, qacc, qt, r0;

    Work() {
        for (int i = 0; i < 3; ++i) {
            mpz_init(u[i]);
            mpz_init(v[i]);
            mpz_init(n[i]);
            mpz_init(px[i]);
            mpq_init(q[i]);
        }
        mpz_init(d);
        mpz_init(t0);
        mpz_init(t1);
        mpz_init(t2);
        mpz_init(det);
        mpz_init(cof);
        mpz_init(pw);
        mpq_init(qa);
        mpq_init(qacc);
        mpq_init(qt);
        mpq_init(r0);
        for (auto& row : m)
            for (auto& e : row) mpz_init(e);
        for (auto& row : na)
            for (auto& e : row) mpz_init(e);
        for (auto& e : da) mpz_init(e);
    }
    ~Work() {
        for (int i = 0; i < 3; ++i) {
            mpz_clear(u[i]);
            mpz_clear(v[i]);
            mpz_clear(n[i]);
            mpz_clear(px[i]);
            mpq_clear(q[i]);
        }
        mpz_clear(d);
        mpz_clear(t0);
        mpz_clear(t1);
        mpz_clear(t2);
        mpz_clear(det);
        mpz_clear(cof);
        mpz_clear(pw);
        mpq_clear(qa);
        mpq_clear(qacc);
        mpq_clear(qt);
        mpq_clear(r0);
        for (auto& row : m)
            for (auto& e : row) mpz_clear(e);
        for (auto& row : na)
            for (auto& e : row) mpz_clear(e);
        for (auto& e : da) mpz_clear(e);
    }
};

/// 3x3 行列式を mpz で計算（余因子展開。ライブラリとは独立に書く）。
void det3_mpz(mpz_ptr out, mpz_t m[3][3], mpz_ptr t0, mpz_ptr t1, mpz_ptr cof) {
    mpz_set_ui(out, 0);
    for (int j = 0; j < 3; ++j) {
        const int j1 = (j + 1) % 3, j2 = (j + 2) % 3;
        mpz_mul(t0, m[1][j1], m[2][j2]);
        mpz_mul(t1, m[1][j2], m[2][j1]);
        mpz_sub(cof, t0, t1);
        mpz_mul(cof, cof, m[0][j]);
        mpz_add(out, out, cof);
    }
}

/// 平面を mpz で作り直す。
void plane_mpz(Work& w, const IPoint& p1, const IPoint& p2, const IPoint& p3) {
    const std::int64_t u[3] = {std::int64_t{p2.x} - p1.x, std::int64_t{p2.y} - p1.y,
                               std::int64_t{p2.z} - p1.z};
    const std::int64_t v[3] = {std::int64_t{p3.x} - p1.x, std::int64_t{p3.y} - p1.y,
                               std::int64_t{p3.z} - p1.z};
    for (int i = 0; i < 3; ++i) {
        mpz_set_si(w.u[i], static_cast<long>(u[i]));
        mpz_set_si(w.v[i], static_cast<long>(v[i]));
    }
    for (int i = 0; i < 3; ++i) {
        const int a = (i + 1) % 3, b = (i + 2) % 3;
        mpz_mul(w.t0, w.u[a], w.v[b]);
        mpz_mul(w.t1, w.u[b], w.v[a]);
        mpz_sub(w.n[i], w.t0, w.t1);
    }
    const std::int64_t c[3] = {p1.x, p1.y, p1.z};
    mpz_set_ui(w.d, 0);
    for (int i = 0; i < 3; ++i) {
        mpz_set_si(w.t0, static_cast<long>(c[i]));
        mpz_mul(w.t0, w.t0, w.n[i]);
        mpz_add(w.d, w.d, w.t0);
    }
    mpz_neg(w.d, w.d);  // 平面は N・x + d = 0、すなわち d = -N・p1（SPEC §3.1）
}

bool same(mpz_srcptr want, const fixed_int<limbs::kNormal>& got) {
    mpz_t g;
    mpz_init(g);
    to_mpz(g, got);
    const bool ok = mpz_cmp(g, want) == 0;
    mpz_clear(g);
    return ok;
}

template <std::size_t N>
bool same_n(mpz_srcptr want, const fixed_int<N>& got) {
    mpz_t g;
    mpz_init(g);
    to_mpz(g, got);
    const bool ok = mpz_cmp(g, want) == 0;
    mpz_clear(g);
    return ok;
}

int sgn_of(int v) {
    return (v < 0) ? -1 : (v > 0) ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    // SPEC §9: 述語 10^6 件以上
    long iters = (argc > 1) ? std::atol(argv[1]) : 1000000L;
    if (iters < 100) iters = 100;

    Rng rng(0x5EED1234u);
    Work w;
    long constructed = 0;

    for (long it = 0; it < iters; ++it) {
        const bool extreme = (it % 3 == 0);
        IPoint p[9];
        for (int i = 0; i < 9; ++i)
            p[i] = extreme ? kritest::rand_extreme_point(rng) : kritest::rand_point(rng);

        // ---- orient3d / orient2d ----
        {
            for (int r = 0; r < 3; ++r) {
                const IPoint& q = p[r];
                mpz_set_si(w.m[r][0], static_cast<long>(std::int64_t{q.x} - p[3].x));
                mpz_set_si(w.m[r][1], static_cast<long>(std::int64_t{q.y} - p[3].y));
                mpz_set_si(w.m[r][2], static_cast<long>(std::int64_t{q.z} - p[3].z));
            }
            det3_mpz(w.det, w.m, w.t0, w.t1, w.cof);
            KRI_CHECK(orient3d(p[0], p[1], p[2], p[3]) == mpz_sgn(w.det));
            KRI_CHECK(same_n(w.det, orient3d_value(p[0], p[1], p[2], p[3])));

            // orient2d は 2x2
            mpz_set_si(w.t0, static_cast<long>((std::int64_t{p[0].x} - p[2].x)));
            mpz_set_si(w.t1, static_cast<long>((std::int64_t{p[1].y} - p[2].y)));
            mpz_mul(w.t0, w.t0, w.t1);
            mpz_set_si(w.t1, static_cast<long>((std::int64_t{p[0].y} - p[2].y)));
            mpz_set_si(w.t2, static_cast<long>((std::int64_t{p[1].x} - p[2].x)));
            mpz_mul(w.t1, w.t1, w.t2);
            mpz_sub(w.det, w.t0, w.t1);
            KRI_CHECK(orient2d(p[0], p[1], p[2], Axis::Z) == mpz_sgn(w.det));
        }

        // ---- 平面の構成 ----
        const PlaneD pl[3] = {plane_from_triangle(p[0], p[1], p[2]),
                              plane_from_triangle(p[3], p[4], p[5]),
                              plane_from_triangle(p[6], p[7], p[8])};
        auto& na = w.na;
        auto& da = w.da;
        for (int k = 0; k < 3; ++k) {
            plane_mpz(w, p[3 * k], p[3 * k + 1], p[3 * k + 2]);
            KRI_CHECK(same(w.n[0], pl[k].a));
            KRI_CHECK(same(w.n[1], pl[k].b));
            KRI_CHECK(same(w.n[2], pl[k].c));
            KRI_CHECK(same_n(w.d, pl[k].d));
            for (int i = 0; i < 3; ++i) mpz_set(na[k][i], w.n[i]);
            mpz_set(da[k], w.d);
        }

        // ---- side(plane, IPoint) を有理数で ----
        {
            mpz_set_ui(w.t2, 0);
            const std::int64_t c[3] = {p[6].x, p[6].y, p[6].z};
            for (int i = 0; i < 3; ++i) {
                mpz_set_si(w.t0, static_cast<long>(c[i]));
                mpz_mul(w.t0, w.t0, na[0][i]);
                mpz_add(w.t2, w.t2, w.t0);
            }
            mpz_add(w.t2, w.t2, da[0]);  // N・p + d
            KRI_CHECK(side(pl[0], p[6]) == mpz_sgn(w.t2));
        }

        // ---- 3 平面の交点 ----
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) mpz_set(w.m[i][j], na[i][j]);
        det3_mpz(w.pw, w.m, w.t0, w.t1, w.cof);

        if (mpz_sgn(w.pw) != 0) {
            ++constructed;
            // 連立は N・X = -d なので、Cramer の置換列は -d（SPEC §3.1 の表の注意点）
            for (int col = 0; col < 3; ++col) {
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) {
                        if (j == col) {
                            mpz_neg(w.m[i][j], da[i]);
                        } else {
                            mpz_set(w.m[i][j], na[i][j]);
                        }
                    }
                det3_mpz(w.px[col], w.m, w.t0, w.t1, w.cof);
            }

            const HPointD hv = intersect3(pl[0], pl[1], pl[2]);
            KRI_CHECK(same_n(w.pw, hv.w));
            KRI_CHECK(same_n(w.px[0], hv.x));
            KRI_CHECK(same_n(w.px[1], hv.y));
            KRI_CHECK(same_n(w.px[2], hv.z));

            // 実座標を有理数で復元
            for (int i = 0; i < 3; ++i) {
                mpq_set_num(w.q[i], w.px[i]);
                mpq_set_den(w.q[i], w.pw);
                mpq_canonicalize(w.q[i]);
            }

            // side(plane, HPoint): sign(N・P - d) を mpq で（sign(w) 変形を使わない）
            for (int k = 0; k < 3; ++k) {
                mpq_set_ui(w.qacc, 0, 1);
                for (int i = 0; i < 3; ++i) {
                    mpq_set_z(w.qa, na[k][i]);
                    mpq_mul(w.qt, w.qa, w.q[i]);
                    mpq_add(w.qacc, w.qacc, w.qt);
                }
                mpq_set_z(w.qa, da[k]);
                mpq_add(w.qacc, w.qacc, w.qa);  // N・P + d
                KRI_CHECK(side(pl[k], hv) == mpq_sgn(w.qacc));
                KRI_CHECK(mpq_sgn(w.qacc) == 0);  // 交点は 3 平面すべての上
            }

            // 交点でない平面に対する side
            const PlaneD probe = plane_from_triangle(p[1], p[4], p[7]);
            plane_mpz(w, p[1], p[4], p[7]);
            mpq_set_ui(w.qacc, 0, 1);
            for (int i = 0; i < 3; ++i) {
                mpq_set_z(w.qa, w.n[i]);
                mpq_mul(w.qt, w.qa, w.q[i]);
                mpq_add(w.qacc, w.qacc, w.qt);
            }
            mpq_set_z(w.qa, w.d);
            mpq_add(w.qacc, w.qacc, w.qa);  // N・P + d
            KRI_CHECK(side(probe, hv) == mpq_sgn(w.qacc));

            // cmp_h: 有理数どうしの比較
            const HPointD hi = to_homogeneous(p[0]);
            const std::int64_t ic[3] = {p[0].x, p[0].y, p[0].z};
            for (int i = 0; i < 3; ++i) {
                mpq_set_si(w.r0, static_cast<long>(ic[i]), 1);
                mpq_canonicalize(w.r0);
                KRI_CHECK(cmp_h(hv, hi, kAxes[i]) == sgn_of(mpq_cmp(w.q[i], w.r0)));
                KRI_CHECK(cmp_h(hi, hv, kAxes[i]) == sgn_of(mpq_cmp(w.r0, w.q[i])));
            }
        }
    }

    std::printf("       %ld 件（うち構成点あり %ld 件）\n", iters, constructed);
    return kritest::finish("geom/predicates_gmp");
}

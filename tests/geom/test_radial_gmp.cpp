// Krisite — radial sort の述語の GMP 差分テスト
//
// `SPEC-phase2.md` §5.1.2.1 の申し送り（辺まわりの二面角の厳密な順序付け）。
// ビット幅の導出は `geom/widths.hpp` の `kRadialDir` / `kRadialDet` / `kRadialDot` /
// `kRadialAlign`。
//
// **正解器は被検体と別経路で書きます。**
//
// 被検体は **2 つの恒等式で簡約した形**を計算します。
//
//   det(d, a×d, b×d) = -|d|^2 det(a, d, b)
//   (a×d)·(b×d)      = (a·b)|d|^2            ← a·d = b·d = 0
//
// **正解器は簡約前の式を mpz でそのまま計算します** — 半平面の方向 m = ±(N×d) を
// 実際に作り、det3(d, m_i, m_j) と m_i·m_j を評価します。
// **恒等式そのものを検査していることになります。** 簡約が誤っていれば落ちます。
//
// **ビット幅の実測もここで行います。** 実測が理論上界を超えたら設計の誤りです。
//
// GMP は LGPL。テスト専用（`KRISITE_BUILD_TESTS_WITH_GMP=ON` でのみビルド）。
#include <cstdio>
#include <cstdlib>

#include "krisite/geom/predicates.hpp"

#include "gmp_oracle.hpp"
#include "test_util.hpp"

using namespace krisite::geom;
using kritest::Rng;
using kritest::oracle::signed_bits;
using kritest::oracle::to_mpz;

namespace {

/// mpz の 3 次元ベクトル（正解器側）。
struct V3 {
    mpz_t x, y, z;
    V3() {
        mpz_init(x);
        mpz_init(y);
        mpz_init(z);
    }
    ~V3() {
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
    }
    V3(const V3&) = delete;
    V3& operator=(const V3&) = delete;
};

/// c = u × v
void ocross(V3* c, const V3& u, const V3& v) {
    mpz_t t1, t2;
    mpz_init(t1);
    mpz_init(t2);
    mpz_mul(t1, u.y, v.z);
    mpz_mul(t2, u.z, v.y);
    mpz_sub(c->x, t1, t2);
    mpz_mul(t1, u.z, v.x);
    mpz_mul(t2, u.x, v.z);
    mpz_sub(c->y, t1, t2);
    mpz_mul(t1, u.x, v.y);
    mpz_mul(t2, u.y, v.x);
    mpz_sub(c->z, t1, t2);
    mpz_clear(t1);
    mpz_clear(t2);
}

/// r = u・v
void odot(mpz_t r, const V3& u, const V3& v) {
    mpz_t t;
    mpz_init(t);
    mpz_mul(r, u.x, v.x);
    mpz_mul(t, u.y, v.y);
    mpz_add(r, r, t);
    mpz_mul(t, u.z, v.z);
    mpz_add(r, r, t);
    mpz_clear(t);
}

/// r = det(a, b, c) = a・(b×c)
void odet3(mpz_t r, const V3& a, const V3& b, const V3& c) {
    V3 bc;
    ocross(&bc, b, c);
    odot(r, a, bc);
}

void neg3(V3* v) {
    mpz_neg(v->x, v->x);
    mpz_neg(v->y, v->y);
    mpz_neg(v->z, v->z);
}

int sgn(const mpz_t x) {
    return mpz_sgn(x);
}

}  // namespace

int main(int argc, char** argv) {
    const long iters = (argc > 1) ? std::strtol(argv[1], nullptr, 10) : 200000;
    Rng rng(0x9e3779b97f4a7c15ull);
    std::printf("\n  radial sort の述語（SPEC-phase2 §5.1.2.1）— GMP 差分 %ld 件\n", iters);

    std::size_t done = 0, deg_skipped = 0;
    std::size_t bits_dir = 0, bits_det = 0, bits_dot = 0, bits_align = 0;
    std::size_t det_zero = 0, dot_zero = 0, align_neg = 0;

    V3 od, omi, omj, oline;
    mpz_t r1, r2, t;
    mpz_init(r1);
    mpz_init(r2);
    mpz_init(t);

    for (long it = 0; it < iters; ++it) {
        // **辺は 2 つの整数点 P, Q で決めます。**
        const IPoint P = kritest::rand_point(rng), Q = kritest::rand_point(rng);
        if (P.x == Q.x && P.y == Q.y && P.z == Q.z) {
            ++deg_skipped;
            continue;
        }
        // **その辺を含む平面を 4 枚作ります**（第 3 点を振る）。
        PlaneD pl[4];
        bool bad = false;
        for (int k = 0; k < 4; ++k) {
            const IPoint R = kritest::rand_point(rng);
            pl[k] = plane_from_triangle(P, Q, R);
            if (is_degenerate(pl[k])) bad = true;
        }
        if (bad) {
            ++deg_skipped;
            continue;
        }

        // ---- (1) 方向 d = N_0 × N_1 ----
        const auto d = radial_dir(pl[0], pl[1]);
        if (krisite::arith::is_zero(d.x) && krisite::arith::is_zero(d.y) &&
            krisite::arith::is_zero(d.z)) {
            ++deg_skipped;  // 2 枚が平行（共面）
            continue;
        }
        to_mpz(od.x, d.x);
        to_mpz(od.y, d.y);
        to_mpz(od.z, d.z);
        for (mpz_srcptr v : {static_cast<mpz_srcptr>(od.x), static_cast<mpz_srcptr>(od.y),
                             static_cast<mpz_srcptr>(od.z)})
            if (signed_bits(v) > bits_dir) bits_dir = signed_bits(v);

        // **正解器: d は辺 PQ と平行か**（別経路。被検体は法線の外積で作っている）
        mpz_set_si(oline.x, static_cast<long>(Q.x) - P.x);
        mpz_set_si(oline.y, static_cast<long>(Q.y) - P.y);
        mpz_set_si(oline.z, static_cast<long>(Q.z) - P.z);
        {
            V3 c;
            ocross(&c, od, oline);
            KRI_CHECK_MSG(mpz_sgn(c.x) == 0 && mpz_sgn(c.y) == 0 && mpz_sgn(c.z) == 0,
                          "radial_dir が辺の方向と平行でない");
        }

        // ---- (2) det と dot を、簡約前の式で照合 ----
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i == j) continue;
                const int ei = ((it + i) & 1) ? 1 : -1;
                const int ej = ((it + j) & 2) ? 1 : -1;
                // **正解器: m = ε(N × d) を実際に作る**
                V3 oni, onj;
                to_mpz(oni.x, pl[i].a);
                to_mpz(oni.y, pl[i].b);
                to_mpz(oni.z, pl[i].c);
                to_mpz(onj.x, pl[j].a);
                to_mpz(onj.y, pl[j].b);
                to_mpz(onj.z, pl[j].c);
                ocross(&omi, oni, od);
                ocross(&omj, onj, od);
                if (ei < 0) neg3(&omi);
                if (ej < 0) neg3(&omj);
                odet3(r1, od, omi, omj);
                odot(r2, omi, omj);

                // **被検体（簡約後）**
                const int sut_det = radial_det(pl[i], d, pl[j]);
                const int sut_dot = radial_dot(pl[i], pl[j]);

                // 恒等式 1: sign(det(d, m_i, m_j)) = -ε_i ε_j sign(det(N_i, d, N_j))
                KRI_CHECK_MSG(sgn(r1) == -ei * ej * sut_det,
                              "det の恒等式が破れた（簡約が誤っている）");
                // 恒等式 2: sign(m_i・m_j) = ε_i ε_j sign(N_i・N_j)
                KRI_CHECK_MSG(sgn(r2) == ei * ej * sut_dot,
                              "内積の恒等式が破れた（簡約が誤っている）");
                if (sut_det == 0) ++det_zero;
                if (sut_dot == 0) ++dot_zero;

                // ビット幅の実測（**簡約後の量**で測る）
                V3 dummy;
                odet3(t, oni, od, onj);
                if (signed_bits(t) > bits_det) bits_det = signed_bits(t);
                odot(t, oni, onj);
                if (signed_bits(t) > bits_dot) bits_dot = signed_bits(t);
            }
        }

        // ---- (3) 向き合わせ ----
        HPointD hu = to_homogeneous(P), hv = to_homogeneous(Q);
        // **同次座標をスカラー倍して w の符号も振る**（実座標は変わらない）
        const std::int64_t ku = static_cast<std::int64_t>(rng.range(-7, 7));
        const std::int64_t kv = static_cast<std::int64_t>(rng.range(-7, 7));
        if (ku != 0 && kv != 0) {
            const auto scale = [](HPointD* h, std::int64_t k) {
                const auto m = krisite::arith::from_i64<1>(k);
                h->x = krisite::arith::resize<limbs::kHomoXyz>(krisite::arith::mul(h->x, m));
                h->y = krisite::arith::resize<limbs::kHomoXyz>(krisite::arith::mul(h->y, m));
                h->z = krisite::arith::resize<limbs::kHomoXyz>(krisite::arith::mul(h->z, m));
                h->w = krisite::arith::resize<limbs::kHomoW>(krisite::arith::mul(h->w, m));
            };
            scale(&hu, ku);
            scale(&hv, kv);
        }
        const int sut_align = radial_align(d, hu, hv);
        // **正解器: 実座標の差 (Q-P) との内積の符号**（同次でなく整数で計算）
        odot(t, od, oline);
        if (signed_bits(t) > bits_align) bits_align = signed_bits(t);
        KRI_CHECK_MSG(sut_align == sgn(t), "radial_align が実座標の向きと合わない");
        if (sut_align < 0) ++align_neg;

        ++done;
    }

    // ---- 退化ケースは明示的に構成します ------------------------------------
    //
    // **乱択は共平面をほぼ生成しません**（上の 20 万件で det = 0 が 0 件）。
    // `CLAUDE.md`「退化ケースを明示的に構成する」に従い、ここで作ります。
    {
        std::size_t built = 0;
        for (long it = 0; it < 4000; ++it) {
            const IPoint P = kritest::rand_point(rng), Q = kritest::rand_point(rng);
            if (P.x == Q.x && P.y == Q.y && P.z == Q.z) continue;
            const IPoint R = kritest::rand_point(rng), Rt = kritest::rand_point(rng);
            const PlaneD pa = plane_from_triangle(P, Q, R);
            const PlaneD pb = plane_from_triangle(P, Q, Rt);
            // **同じ平面・逆向きの法線**（第 1 点と第 2 点を入れ替える）
            const PlaneD pa_rev = plane_from_triangle(Q, P, R);
            if (is_degenerate(pa) || is_degenerate(pb) || is_degenerate(pa_rev)) continue;
            const auto d = radial_dir(pa, pb);
            if (krisite::arith::is_zero(d.x) && krisite::arith::is_zero(d.y) &&
                krisite::arith::is_zero(d.z))
                continue;
            ++built;
            // **同一平面上の半平面どうしは det = 0**
            KRI_CHECK_MSG(radial_det(pa, d, pa) == 0, "同じ平面で det が 0 でない");
            KRI_CHECK_MSG(radial_det(pa, d, pa_rev) == 0, "逆向きの同じ平面で det が 0 でない");
            // **逆向きなので内積は負**
            KRI_CHECK_MSG(radial_dot(pa, pa_rev) < 0, "逆向きの法線の内積が負でない");
            KRI_CHECK_MSG(radial_dot(pa, pa) > 0, "同じ法線の内積が正でない");
            det_zero += 2;
        }
        KRI_CHECK_MSG(built > 1000, "退化ケースの構成が足りない");
        std::printf("    **構成した退化ケース: %zu 件**（同一平面・逆向きの法線）\n", built);
    }
    // ---- 軸平行で内積 0 を作る ----------------------------------------------
    {
        // 辺を z 軸に置き、xz 平面と yz 平面を取る → 法線は (0,±1,0) と (±1,0,0)
        const IPoint P{0, 0, -1000}, Q{0, 0, 1000};
        const PlaneD pxz = plane_from_triangle(P, Q, IPoint{1000, 0, 0});
        const PlaneD pyz = plane_from_triangle(P, Q, IPoint{0, 1000, 0});
        KRI_CHECK_MSG(!is_degenerate(pxz) && !is_degenerate(pyz), "軸平行の構成が退化した");
        const auto d = radial_dir(pxz, pyz);
        KRI_CHECK_MSG(radial_dot(pxz, pyz) == 0, "直交する法線の内積が 0 でない");
        KRI_CHECK_MSG(radial_det(pxz, d, pyz) != 0, "直交する半平面の det が 0 になった");
        ++dot_zero;
        std::printf("    **軸平行で内積 0 を構成しました**\n");
    }
    // ---- 構成点で align の幅を使い切る --------------------------------------
    //
    // **上の本体は整数点（w = 1）しか通していないので、`kRadialAlign` の幅を
    // 使い切りません**（実測 106 / 上界 393）。**構成点で測り直します。**
    {
        std::size_t built = 0;
        mpq_t qa[3], qb[3], dif, acc, dq;
        for (int k = 0; k < 3; ++k) {
            mpq_init(qa[k]);
            mpq_init(qb[k]);
        }
        mpq_init(dif);
        mpq_init(acc);
        mpq_init(dq);
        for (long it = 0; it < 6000; ++it) {
            const IPoint P = kritest::rand_point(rng), Q = kritest::rand_point(rng);
            if (P.x == Q.x && P.y == Q.y && P.z == Q.z) continue;
            const PlaneD pa = plane_from_triangle(P, Q, kritest::rand_point(rng));
            const PlaneD pb = plane_from_triangle(P, Q, kritest::rand_point(rng));
            // **辺を横切る 2 枚**（線上の 2 点を作る）
            const PlaneD c1 = plane_from_triangle(
                kritest::rand_point(rng), kritest::rand_point(rng), kritest::rand_point(rng));
            const PlaneD c2 = plane_from_triangle(
                kritest::rand_point(rng), kritest::rand_point(rng), kritest::rand_point(rng));
            if (is_degenerate(pa) || is_degenerate(pb) || is_degenerate(c1) || is_degenerate(c2))
                continue;
            const auto d = radial_dir(pa, pb);
            if (krisite::arith::is_zero(d.x) && krisite::arith::is_zero(d.y) &&
                krisite::arith::is_zero(d.z))
                continue;
            const HPointD hu = intersect3(pa, pb, c1);
            const HPointD hv = intersect3(pa, pb, c2);
            if (krisite::arith::sign(hu.w) == 0 || krisite::arith::sign(hv.w) == 0) continue;
            if (h_equal(hu, hv)) continue;
            ++built;
            const int sut = radial_align(d, hu, hv);
            // **正解器: 有理数で実座標に直してから差を取る**（除算を使う別経路）
            const auto setq = [](mpq_t* out, const HPointD& h) {
                mpz_t num, den;
                mpz_init(num);
                mpz_init(den);
                to_mpz(den, h.w);
                to_mpz(num, h.x);
                mpq_set_num(out[0], num);
                mpq_set_den(out[0], den);
                to_mpz(num, h.y);
                mpq_set_num(out[1], num);
                mpq_set_den(out[1], den);
                to_mpz(num, h.z);
                mpq_set_num(out[2], num);
                mpq_set_den(out[2], den);
                for (int k = 0; k < 3; ++k) mpq_canonicalize(out[k]);
                mpz_clear(num);
                mpz_clear(den);
            };
            setq(qa, hu);
            setq(qb, hv);
            mpq_set_ui(acc, 0, 1);
            const krisite::arith::fixed_int<limbs::kRadialDir>* dc[3] = {&d.x, &d.y, &d.z};
            for (int k = 0; k < 3; ++k) {
                mpq_sub(dif, qb[k], qa[k]);
                mpz_t dz;
                mpz_init(dz);
                to_mpz(dz, *dc[k]);
                mpq_set_z(dq, dz);
                mpz_clear(dz);
                mpq_mul(dif, dif, dq);
                mpq_add(acc, acc, dif);
            }
            KRI_CHECK_MSG(sut == mpq_sgn(acc), "構成点で radial_align が有理数の答えと違う");
            // 幅の実測（簡約前の e = w_u V_v - w_v V_u と d の内積）
            {
                mpz_t wu, wv, xu, xv, e, s, p;
                mpz_init(wu);
                mpz_init(wv);
                mpz_init(xu);
                mpz_init(xv);
                mpz_init(e);
                mpz_init(s);
                mpz_init(p);
                to_mpz(wu, hu.w);
                to_mpz(wv, hv.w);
                mpz_set_ui(s, 0);
                const krisite::arith::fixed_int<limbs::kHomoXyz>* uu[3] = {&hu.x, &hu.y, &hu.z};
                const krisite::arith::fixed_int<limbs::kHomoXyz>* vv[3] = {&hv.x, &hv.y, &hv.z};
                for (int k = 0; k < 3; ++k) {
                    to_mpz(xu, *uu[k]);
                    to_mpz(xv, *vv[k]);
                    mpz_mul(e, wu, xv);
                    mpz_mul(p, wv, xu);
                    mpz_sub(e, e, p);
                    mpz_t dz;
                    mpz_init(dz);
                    to_mpz(dz, *dc[k]);
                    mpz_mul(e, e, dz);
                    mpz_clear(dz);
                    mpz_add(s, s, e);
                }
                if (signed_bits(s) > bits_align) bits_align = signed_bits(s);
                mpz_clear(wu);
                mpz_clear(wv);
                mpz_clear(xu);
                mpz_clear(xv);
                mpz_clear(e);
                mpz_clear(s);
                mpz_clear(p);
            }
        }
        for (int k = 0; k < 3; ++k) {
            mpq_clear(qa[k]);
            mpq_clear(qb[k]);
        }
        mpq_clear(dif);
        mpq_clear(acc);
        mpq_clear(dq);
        KRI_CHECK_MSG(built > 1000, "構成点の生成が足りない");
        std::printf("    **構成点で照合: %zu 件**（有理数の正解器）\n", built);
    }

    std::printf("    実測 %zu 件 / 退化で除外 %zu\n", done, deg_skipped);
    std::printf("    det が 0（同一平面上の半平面）: %zu / 内積が 0: %zu / d の反転: %zu\n",
                det_zero, dot_zero, align_neg);
    std::printf(
        "    ビット幅（実測 / 理論上界）: dir %zu/%zu, det %zu/%zu, dot %zu/%zu, "
        "align %zu/%zu\n",
        bits_dir, bits::kRadialDir, bits_det, bits::kRadialDet, bits_dot, bits::kRadialDot,
        bits_align, bits::kRadialAlign);

    // **理論上界を実測が超えたら設計の誤り。即座に報告する対象**
    KRI_CHECK_MSG(bits_dir <= bits::kRadialDir, "kRadialDir の実測が理論上界を超えた");
    KRI_CHECK_MSG(bits_det <= bits::kRadialDet, "kRadialDet の実測が理論上界を超えた");
    KRI_CHECK_MSG(bits_dot <= bits::kRadialDot, "kRadialDot の実測が理論上界を超えた");
    KRI_CHECK_MSG(bits_align <= bits::kRadialAlign, "kRadialAlign の実測が理論上界を超えた");
    // **空回り防止**
    KRI_CHECK_MSG(done > static_cast<std::size_t>(iters) / 2, "退化で弾かれすぎている");
    KRI_CHECK_MSG(det_zero > 0, "det が 0 になる場合が 1 件も出ていない");
    KRI_CHECK_MSG(align_neg > done / 20, "d の反転が起きる場合が少なすぎる");

    mpz_clear(r1);
    mpz_clear(r2);
    mpz_clear(t);
    std::printf("\n");
    return kritest::finish("geom/radial_gmp");
}

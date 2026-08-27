// Krisite — 平面とセルの閉領域の交差判定の GMP 差分テスト
//
// SPEC-phase2.md §2.3 / §7
//
// **正解器は被検体と別経路で書きます。** 被検体は
//
//   (1) N の各成分の符号で lo / hi を選び
//   (2) 内積 2 回で N・x + d の最小値と最大値を出し
//   (3) 最小値 <= 0 <= 最大値 を見る
//
// という 3 段です。正解器は **8 隅すべてを mpz で評価して最小・最大を取ります。**
// (1) の選択も (2) の固定幅算術も共有しないので、どちらの誤りも検出できます。
//
// **ビット幅の実測もここで行います**（§7、SPEC-phase0 §8.4 と同じ規律）。
// 実測が理論上界 3b+7 を超えたら設計の誤りなので即座に停止します。
//
// GMP は LGPL。テスト専用（CMake の KRISITE_BUILD_TESTS_WITH_GMP=ON でのみビルド）。
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

struct Work {
    mpz_t a, b, c, d, x, y, z, acc, t;
    Work() {
        mpz_init(a);
        mpz_init(b);
        mpz_init(c);
        mpz_init(d);
        mpz_init(x);
        mpz_init(y);
        mpz_init(z);
        mpz_init(acc);
        mpz_init(t);
    }
    ~Work() {
        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(c);
        mpz_clear(d);
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);
        mpz_clear(acc);
        mpz_clear(t);
    }
};

/// N・x + d を mpz で評価する。
void eval(Work& w, const std::int64_t p[3]) {
    mpz_set_si(w.x, static_cast<long>(p[0]));
    mpz_set_si(w.y, static_cast<long>(p[1]));
    mpz_set_si(w.z, static_cast<long>(p[2]));
    mpz_mul(w.acc, w.a, w.x);
    mpz_mul(w.t, w.b, w.y);
    mpz_add(w.acc, w.acc, w.t);
    mpz_mul(w.t, w.c, w.z);
    mpz_add(w.acc, w.acc, w.t);
    mpz_add(w.acc, w.acc, w.d);
}

std::int64_t rand_bound(Rng& rng) {
    const std::int64_t span = -krisite::kCoordMin;  // 2^(b-1)
    return static_cast<std::int64_t>(rng.next() % static_cast<std::uint64_t>(2 * span + 1)) - span;
}

}  // namespace

int main(int argc, char** argv) {
    const long iters = (argc > 1) ? std::atol(argv[1]) : 200000;
    std::printf("\n  plane_crosses_box の GMP 差分テスト（%ld 件）\n", iters);

    Rng rng(9021);
    Work w;
    std::size_t max_bits = 0, crossing = 0, done = 0;
    // 交差の判定に使った箱の隅の側を数える。**空回り防止**
    std::size_t touching = 0;  // 最小値 == 0 または最大値 == 0（境界に載る）

    for (long it = 0; it < iters; ++it) {
        // **箱を先に作ります。** 平面を先に決めると、境界にちょうど載る場合が
        // 乱択では 1 件も出ません（実測済み）。退化は構成するものです。
        std::int64_t lo[3], hi[3];
        for (int t = 0; t < 3; ++t) {
            const std::int64_t u = rand_bound(rng), v = rand_bound(rng);
            lo[t] = u < v ? u : v;
            hi[t] = u < v ? v : u;
            // 退化した箱（面・辺・点）も作る。乱択はこれをほぼ生成しない
            if (it % 7 == 0) hi[t] = lo[t];
        }

        const Axis kAxes[3] = {Axis::X, Axis::Y, Axis::Z};
        PlaneD pl;
        switch (it % 6) {
            case 0: {
                // **箱の lo 面にちょうど載る**軸平行平面（最小値 = 0）
                const int t = static_cast<int>(it / 6) % 3;
                pl = plane_axis_aligned(kAxes[t], lo[t]);
                break;
            }
            case 1: {
                // **箱の hi 面にちょうど載る**軸平行平面（最大値 = 0）
                const int t = static_cast<int>(it / 6) % 3;
                pl = plane_axis_aligned(kAxes[t], hi[t]);
                break;
            }
            case 2: {
                // 境界の 1 格子だけ外側（交差しないぎりぎり）
                const int t = static_cast<int>(it / 6) % 3;
                const std::int64_t c = (it % 12 < 6) ? lo[t] - 1 : hi[t] + 1;
                if (c < krisite::kCoordMin || c > -krisite::kCoordMin) continue;
                pl = plane_axis_aligned(kAxes[t], c);
                break;
            }
            case 3:
                pl = plane_from_triangle(kritest::rand_extreme_point(rng),
                                         kritest::rand_extreme_point(rng),
                                         kritest::rand_extreme_point(rng));
                break;
            default:
                pl = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                         kritest::rand_point(rng));
                break;
        }
        if (is_degenerate(pl)) continue;

        to_mpz(w.a, pl.a);
        to_mpz(w.b, pl.b);
        to_mpz(w.c, pl.c);
        to_mpz(w.d, pl.d);

        // ---- 正解器: 8 隅すべてを評価 ----
        int lo_sign = +1, hi_sign = -1;
        for (int m = 0; m < 8; ++m) {
            const std::int64_t p[3] = {(m & 1) ? hi[0] : lo[0], (m & 2) ? hi[1] : lo[1],
                                       (m & 4) ? hi[2] : lo[2]};
            eval(w, p);
            const int s = mpz_sgn(w.acc);
            if (s < lo_sign) lo_sign = s;
            if (s > hi_sign) hi_sign = s;
            const std::size_t bits = signed_bits(w.acc);
            if (bits > max_bits) max_bits = bits;
        }
        const bool want = (lo_sign <= 0 && hi_sign >= 0);

        const bool got = plane_crosses_box(pl, lo, hi);
        KRI_CHECK_MSG(got == want, "GMP の 8 隅評価と一致しない");

        // 被検体が使う 2 点のビット幅も実測する（§7 の上界検査）
        {
            const int sa[3] = {krisite::arith::sign(pl.a), krisite::arith::sign(pl.b),
                               krisite::arith::sign(pl.c)};
            std::int64_t pmin[3], pmax[3];
            for (int t = 0; t < 3; ++t) {
                pmin[t] = (sa[t] >= 0) ? lo[t] : hi[t];
                pmax[t] = (sa[t] >= 0) ? hi[t] : lo[t];
            }
            eval(w, pmin);
            const int smin = mpz_sgn(w.acc);
            if (signed_bits(w.acc) > max_bits) max_bits = signed_bits(w.acc);
            eval(w, pmax);
            const int smax = mpz_sgn(w.acc);
            if (signed_bits(w.acc) > max_bits) max_bits = signed_bits(w.acc);
            // **符号による選択が本当に最小・最大を与えているか**（被検体の (1) の検証）
            KRI_CHECK_MSG(smin == lo_sign, "符号で選んだ点が最小値の符号を与えていない");
            KRI_CHECK_MSG(smax == hi_sign, "符号で選んだ点が最大値の符号を与えていない");
            if (smin == 0 || smax == 0) ++touching;
        }

        if (got) ++crossing;
        ++done;
    }

    const std::size_t bound = krisite::geom::bits::kPlaneAabb;
    std::printf("    実測 %zu 件 / 交差 %zu (%.1f%%) / 境界に載る %zu\n", done, crossing,
                100.0 * static_cast<double>(crossing) / static_cast<double>(done), touching);
    std::printf("    ビット幅: 実測 %zu / 理論上界 %zu（%zu リム）\n", max_bits, bound,
                krisite::geom::limbs::kPlaneAabb);

    // §7: **理論上界を実測が超えたら設計の誤り。即座に報告する対象**
    KRI_CHECK_MSG(max_bits <= bound, "実測ビット幅が理論上界 3b+7 を超えた。設計の誤りです");
    // 空回り防止
    KRI_CHECK_MSG(done > static_cast<std::size_t>(iters) / 2, "退化で弾かれすぎている");
    KRI_CHECK_MSG(touching > done / 20,
                  "境界にちょうど載る場合が少なすぎる。乱択では出ないので構成すること");
    KRI_CHECK_MSG(crossing > done / 20 && crossing < done - done / 20,
                  "交差する / しないの片方しか出ていない");

    std::printf("\n");
    return kritest::finish("geom/plane_box_gmp");
}

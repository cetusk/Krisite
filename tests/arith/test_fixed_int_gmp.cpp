// Krisite — fixed_int の GMP 差分テスト
//
// SPEC-phase0.md §8.1-1
//   ランダムな fixed_int × 10^7 組について add / sub / mul / cmp を GMP と比較。
//   ビット幅は境界値（全ビット 1、最小値、0、±1）を重点的に。
//
// GMP は LGPL。テスト専用（KRISITE_BUILD_TESTS_WITH_GMP=ON のときだけビルドされる）。
#include <cstdlib>

#include "gmp_oracle.hpp"
#include "test_util.hpp"

using namespace krisite::arith;
using kritest::Rng;
using kritest::oracle::from_mpz;
using kritest::oracle::signed_bits;
using kritest::oracle::to_mpz;

namespace {

/// 作業用の mpz をループ外で確保しておくための束。
struct Scratch {
    mpz_t za, zb, zr, zt, zpow;
    Scratch() {
        mpz_init(za);
        mpz_init(zb);
        mpz_init(zr);
        mpz_init(zt);
        mpz_init(zpow);
    }
    ~Scratch() {
        mpz_clear(za);
        mpz_clear(zb);
        mpz_clear(zr);
        mpz_clear(zt);
        mpz_clear(zpow);
    }
};

template <std::size_t N>
void check_equal(mpz_srcptr want, const fixed_int<N>& got, const char* op, mpz_srcptr a,
                 mpz_srcptr b) {
    mpz_t g;
    mpz_init(g);
    to_mpz(g, got);
    if (mpz_cmp(g, want) != 0) {
        char* sa = mpz_get_str(nullptr, 10, a);
        char* sb = mpz_get_str(nullptr, 10, b);
        char* sw = mpz_get_str(nullptr, 10, want);
        char* sg = mpz_get_str(nullptr, 10, g);
        ++::kritest::g_checks;
        ::kritest::report(op, __FILE__, __LINE__,
                          std::string("a=") + sa + " b=" + sb + " 期待=" + sw + " 実測=" + sg);
        void (*freefn)(void*, size_t);
        mp_get_memory_functions(nullptr, nullptr, &freefn);
        freefn(sa, std::char_traits<char>::length(sa) + 1);
        freefn(sb, std::char_traits<char>::length(sb) + 1);
        freefn(sw, std::char_traits<char>::length(sw) + 1);
        freefn(sg, std::char_traits<char>::length(sg) + 1);
    } else {
        ++::kritest::g_checks;
    }
    mpz_clear(g);
}

template <std::size_t N>
void run(long iters, std::uint64_t seed) {
    Rng rng(seed);
    Scratch s;
    mpz_ui_pow_ui(s.zpow, 2, static_cast<unsigned long>(64 * N));

    for (long it = 0; it < iters; ++it) {
        const fixed_int<N> a = kritest::rand_biased<N>(rng);
        const fixed_int<N> b = kritest::rand_biased<N>(rng);
        to_mpz(s.za, a);
        to_mpz(s.zb, b);

        // add_widen / sub_widen — 定義上必ず収まる
        mpz_add(s.zr, s.za, s.zb);
        check_equal(s.zr, add_widen(a, b), "add_widen", s.za, s.zb);
        if (kritest::oracle::fits_signed(s.zr, N)) {
            check_equal(s.zr, add(a, b), "add", s.za, s.zb);
        }

        mpz_sub(s.zr, s.za, s.zb);
        check_equal(s.zr, sub_widen(a, b), "sub_widen", s.za, s.zb);
        if (kritest::oracle::fits_signed(s.zr, N)) {
            check_equal(s.zr, sub(a, b), "sub", s.za, s.zb);
        }

        // mul — N+M リムに必ず収まる（SPEC §5.2）
        mpz_mul(s.zr, s.za, s.zb);
        check_equal(s.zr, mul(a, b), "mul", s.za, s.zb);
        KRI_CHECK(kritest::oracle::fits_signed(s.zr, 2 * N));

        // cmp / sign
        const int want_cmp = mpz_cmp(s.za, s.zb);
        const int want_norm = (want_cmp < 0) ? -1 : (want_cmp > 0) ? 1 : 0;
        KRI_CHECK(cmp(a, b) == want_norm);
        KRI_CHECK(sign(a) == mpz_sgn(s.za));

        // neg
        mpz_neg(s.zr, s.za);
        if (kritest::oracle::fits_signed(s.zr, N)) {
            check_equal(s.zr, neg(a), "neg", s.za, s.zb);
        }

        // min_bits（SPEC §8.4 で使う測定器そのものの検証）
        KRI_CHECK(min_bits(a) == signed_bits(s.za));

        // シフト
        const std::size_t sh = static_cast<std::size_t>(rng.below(64 * N + 8));
        mpz_mul_2exp(s.zr, s.za, static_cast<mp_bitcnt_t>(sh));
        mpz_fdiv_r_2exp(s.zt, s.zr, static_cast<mp_bitcnt_t>(64 * N));
        if (mpz_tstbit(s.zt, 64 * N - 1)) mpz_sub(s.zt, s.zt, s.zpow);  // 2 の補数へ
        check_equal(s.zt, shl_bits(a, sh), "shl_bits", s.za, s.zb);

        mpz_fdiv_q_2exp(s.zr, s.za, static_cast<mp_bitcnt_t>(sh));  // 算術右シフト = 床除算
        check_equal(s.zr, shr_bits(a, sh), "shr_bits", s.za, s.zb);
    }
}

}  // namespace

int main(int argc, char** argv) {
    // 既定は SPEC §8.1 の 10^7 組。幅ごとに分割する。
    long total = (argc > 1) ? std::atol(argv[1]) : 10000000L;
    if (total < 4) total = 4;
    const long per = total / 4;

    run<1>(per, 0xC0FFEEull);
    run<2>(per, 0xDECAFull);
    run<3>(per, 0xFACADEull);
    run<4>(per, 0xBADC0DEull);

    std::printf("       %ld 組を検証（幅 1..4 リム）\n", per * 4);
    return kritest::finish("arith/fixed_int_gmp");
}

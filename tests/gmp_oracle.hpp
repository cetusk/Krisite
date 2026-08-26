// Krisite — GMP 正解器（テスト専用）
//
// SPEC-phase0.md §8.1
//
// **ライセンス上の重要な制約**: GMP は LGPL。このヘッダはテストからのみ include し、
// ライブラリ本体（include/krisite/）から参照してはならない。
// ビルドは KRISITE_BUILD_TESTS_WITH_GMP=ON のときだけ行われる。
#ifndef KRISITE_TESTS_GMP_ORACLE_HPP
#define KRISITE_TESTS_GMP_ORACLE_HPP

#include <string>

#include <gmp.h>

#include "krisite/krisite.hpp"

namespace kritest::oracle {

using krisite::arith::fixed_int;

/// RAII な mpz_t。テストコードなので動的確保の禁止は適用されない。
class Z {
public:
    Z() noexcept { mpz_init(v_); }
    explicit Z(long x) noexcept { mpz_init_set_si(v_, x); }
    Z(const Z& o) noexcept { mpz_init_set(v_, o.v_); }
    Z& operator=(const Z& o) noexcept {
        if (this != &o) mpz_set(v_, o.v_);
        return *this;
    }
    ~Z() { mpz_clear(v_); }

    mpz_ptr get() noexcept { return v_; }
    mpz_srcptr get() const noexcept { return v_; }

    std::string str() const {
        char* s = mpz_get_str(nullptr, 10, v_);
        std::string r(s);
        void (*freefn)(void*, size_t);
        mp_get_memory_functions(nullptr, nullptr, &freefn);
        freefn(s, r.size() + 1);
        return r;
    }

private:
    mpz_t v_;
};

/// fixed_int<N> → mpz（2 の補数を符号付き整数として解釈）。
template <std::size_t N>
inline void to_mpz(mpz_ptr out, const fixed_int<N>& x) {
    mpz_import(out, N, -1 /* 最下位ワードが先頭 */, sizeof(std::uint64_t), 0 /* native endian */, 0,
               x.limb.data());
    if (krisite::arith::is_negative(x)) {
        mpz_t m;
        mpz_init(m);
        mpz_ui_pow_ui(m, 2, static_cast<unsigned long>(64 * N));
        mpz_sub(out, out, m);
        mpz_clear(m);
    }
}

template <std::size_t N>
inline Z to_z(const fixed_int<N>& x) {
    Z z;
    to_mpz(z.get(), x);
    return z;
}

/// 値が符号付き 64N ビットに収まるか。
inline bool fits_signed(mpz_srcptr v, std::size_t N) {
    mpz_t lo, hi;
    mpz_init(lo);
    mpz_init(hi);
    mpz_ui_pow_ui(hi, 2, static_cast<unsigned long>(64 * N - 1));
    mpz_neg(lo, hi);
    mpz_sub_ui(hi, hi, 1);
    const bool ok = (mpz_cmp(v, lo) >= 0) && (mpz_cmp(v, hi) <= 0);
    mpz_clear(lo);
    mpz_clear(hi);
    return ok;
}

/// mpz → fixed_int<N>（mod 2^(64N) で 2 の補数に落とす）。
template <std::size_t N>
inline fixed_int<N> from_mpz(mpz_srcptr v) {
    mpz_t t;
    mpz_init(t);
    mpz_fdiv_r_2exp(t, v, static_cast<mp_bitcnt_t>(64 * N));  // 常に非負
    fixed_int<N> r = krisite::arith::zero<N>();
    std::size_t count = 0;
    mpz_export(r.limb.data(), &count, -1, sizeof(std::uint64_t), 0, 0, t);
    mpz_clear(t);
    return r;
}

/// 必要な符号付きビット数（min_bits() の正解器）。
inline std::size_t signed_bits(mpz_srcptr v) {
    if (mpz_sgn(v) == 0) return 1;
    if (mpz_sgn(v) > 0) return mpz_sizeinbase(v, 2) + 1;
    mpz_t t;
    mpz_init(t);
    mpz_neg(t, v);
    mpz_sub_ui(t, t, 1);  // n = bit_length(|v| - 1) + 1
    const std::size_t n = (mpz_sgn(t) == 0) ? 1 : mpz_sizeinbase(t, 2) + 1;
    mpz_clear(t);
    return n;
}

}  // namespace kritest::oracle

#endif  // KRISITE_TESTS_GMP_ORACLE_HPP

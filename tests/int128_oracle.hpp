// Krisite — テスト専用の可搬 128bit 符号付き整数（正解器）
//
// `fixed_int<2>` を突き合わせるための正解器。**64bit 演算だけで書いてある**ので、
// `unsigned __int128` を持たない処理系（MSVC。SPEC-phase0.md §5.3）でも同じ検証が走る。
//
// ライブラリ本体とは独立の実装であることが正解器としての価値なので、
// `krisite::arith` の関数を一切呼ばないこと。乗算も 32bit 部分積を 128bit 加算で
// 積み上げる形にしてあり、本体の「符号なし筆算 + 符号補正」とは別の道筋を通る。
//
// `__int128` がある処理系では i128_selfcheck() で正解器そのものを検証する。
// MSVC ではその検証は走らないが、他のプラットフォームの CI が通っていれば
// 正解器の正しさは担保されている。
#ifndef KRISITE_TESTS_INT128_ORACLE_HPP
#define KRISITE_TESTS_INT128_ORACLE_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "test_util.hpp"

namespace kritest {

/// 符号付き 128bit。2 の補数、lo が下位 64bit。
struct I128 {
    std::uint64_t lo, hi;
};

inline bool i128_is_neg(I128 a) noexcept {
    return (a.hi >> 63) != 0;
}

inline bool i128_is_zero(I128 a) noexcept {
    return a.lo == 0 && a.hi == 0;
}

inline I128 i128_add(I128 a, I128 b) noexcept {
    I128 r;
    r.lo = a.lo + b.lo;
    r.hi = a.hi + b.hi + (r.lo < a.lo ? 1u : 0u);
    return r;
}

inline I128 i128_sub(I128 a, I128 b) noexcept {
    I128 r;
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - (a.lo < b.lo ? 1u : 0u);
    return r;
}

inline I128 i128_negate(I128 a) noexcept {
    return i128_sub(I128{0, 0}, a);
}

/// 符号付き比較。-1 / 0 / +1
inline int i128_cmp(I128 a, I128 b) noexcept {
    const std::int64_t ah = static_cast<std::int64_t>(a.hi);
    const std::int64_t bh = static_cast<std::int64_t>(b.hi);
    if (ah != bh) return ah < bh ? -1 : 1;
    if (a.lo != b.lo) return a.lo < b.lo ? -1 : 1;
    return 0;
}

inline bool i128_eq(I128 a, I128 b) noexcept {
    return a.lo == b.lo && a.hi == b.hi;
}

inline int i128_sign(I128 a) noexcept {
    if (i128_is_neg(a)) return -1;
    return i128_is_zero(a) ? 0 : 1;
}

inline I128 i128_from_i64(std::int64_t v) noexcept {
    return I128{static_cast<std::uint64_t>(v), v < 0 ? ~std::uint64_t{0} : std::uint64_t{0}};
}

inline I128 i128_from_u64(std::uint64_t v) noexcept {
    return I128{v, 0};
}

/// 論理左シフト（mod 2^128）。
inline I128 i128_shl(I128 a, std::size_t s) noexcept {
    if (s >= 128) return I128{0, 0};
    if (s == 0) return a;
    if (s >= 64) return I128{0, a.lo << (s - 64)};
    return I128{a.lo << s, (a.hi << s) | (a.lo >> (64 - s))};
}

/// 算術右シフト。C++20 では符号付きの >> は算術シフトと規定されている。
inline I128 i128_sar(I128 a, std::size_t s) noexcept {
    const std::uint64_t fill = i128_is_neg(a) ? ~std::uint64_t{0} : std::uint64_t{0};
    if (s >= 128) return I128{fill, fill};
    if (s == 0) return a;
    if (s >= 64) {
        const std::int64_t h = static_cast<std::int64_t>(a.hi) >> (s - 64);
        return I128{static_cast<std::uint64_t>(h), fill};
    }
    return I128{(a.lo >> s) | (a.hi << (64 - s)),
                static_cast<std::uint64_t>(static_cast<std::int64_t>(a.hi) >> s)};
}

/// 符号付き 64x64 -> 128。
///
/// 絶対値どうしを 32bit 部分積に分けて 128bit 加算で積み上げ、最後に符号を付ける。
/// |INT64_MIN| = 2^63 は符号なし 64bit に収まるので、最小値も特別扱いは要らない。
inline I128 i128_mul_i64(std::int64_t a, std::int64_t b) noexcept {
    const bool na = a < 0, nb = b < 0;
    const std::uint64_t ua =
        na ? (~static_cast<std::uint64_t>(a) + 1) : static_cast<std::uint64_t>(a);
    const std::uint64_t ub =
        nb ? (~static_cast<std::uint64_t>(b) + 1) : static_cast<std::uint64_t>(b);

    const std::uint64_t a0 = ua & 0xFFFFFFFFull, a1 = ua >> 32;
    const std::uint64_t b0 = ub & 0xFFFFFFFFull, b1 = ub >> 32;

    I128 r{0, 0};
    r = i128_add(r, I128{a0 * b0, 0});
    r = i128_add(r, i128_shl(I128{a0 * b1, 0}, 32));
    r = i128_add(r, i128_shl(I128{a1 * b0, 0}, 32));
    r = i128_add(r, I128{0, a1 * b1});
    return (na != nb) ? i128_negate(r) : r;
}

/// 10 進表記（失敗時の診断用）。
inline std::string i128_str(I128 v) {
    if (i128_is_zero(v)) return "0";
    const bool neg = i128_is_neg(v);
    const I128 m = neg ? i128_negate(v) : v;

    std::uint32_t w[4] = {static_cast<std::uint32_t>(m.lo), static_cast<std::uint32_t>(m.lo >> 32),
                          static_cast<std::uint32_t>(m.hi), static_cast<std::uint32_t>(m.hi >> 32)};
    std::vector<std::uint32_t> groups;  // 下位から 10^9 ごとの桁group
    while (w[0] || w[1] || w[2] || w[3]) {
        std::uint64_t rem = 0;
        for (int i = 3; i >= 0; --i) {
            const std::uint64_t cur = (rem << 32) | w[i];
            w[i] = static_cast<std::uint32_t>(cur / 1000000000ull);
            rem = cur % 1000000000ull;
        }
        groups.push_back(static_cast<std::uint32_t>(rem));
    }

    std::string s = std::to_string(groups.back());
    for (std::size_t i = groups.size() - 1; i-- > 0;) {
        const std::string g = std::to_string(groups[i]);
        s += std::string(9 - g.size(), '0') + g;
    }
    return neg ? "-" + s : s;
}

// ---- fixed_int<2> との相互変換 ----------------------------------------------

inline I128 i128_of(const krisite::arith::fixed_int<2>& x) noexcept {
    return I128{x.limb[0], x.limb[1]};
}

inline krisite::arith::fixed_int<2> fixed_of(I128 v) noexcept {
    krisite::arith::fixed_int<2> r{};
    r.limb[0] = v.lo;
    r.limb[1] = v.hi;
    return r;
}

// ---- 正解器そのものの検証 ---------------------------------------------------

/// `__int128` がある処理系でのみ、I128 の各演算を突き合わせる。
/// MSVC ではこの関数は何もしない（他プラットフォームの CI で担保する）。
inline void i128_selfcheck(int iters = 20000) {
#if defined(__SIZEOF_INT128__)
    using ref = __int128;
    using uref = unsigned __int128;
    auto to_ref = [](I128 v) { return static_cast<ref>((static_cast<uref>(v.hi) << 64) | v.lo); };
    Rng rng(0xA11CE);
    for (int i = 0; i < iters; ++i) {
        const I128 a{rng.next(), rng.next()};
        const I128 b{rng.next(), rng.next()};
        const ref ra = to_ref(a), rb = to_ref(b);

        KRI_CHECK(to_ref(i128_add(a, b)) ==
                  static_cast<ref>(static_cast<uref>(ra) + static_cast<uref>(rb)));
        KRI_CHECK(to_ref(i128_sub(a, b)) ==
                  static_cast<ref>(static_cast<uref>(ra) - static_cast<uref>(rb)));
        KRI_CHECK(to_ref(i128_negate(a)) == static_cast<ref>(~static_cast<uref>(ra) + 1));
        KRI_CHECK(i128_cmp(a, b) == ((ra < rb) ? -1 : (ra > rb) ? 1 : 0));
        KRI_CHECK(i128_sign(a) == ((ra < 0) ? -1 : (ra > 0) ? 1 : 0));

        const std::size_t s = static_cast<std::size_t>(rng.below(128));
        KRI_CHECK(to_ref(i128_shl(a, s)) == static_cast<ref>(static_cast<uref>(ra) << s));
        KRI_CHECK(to_ref(i128_sar(a, s)) == (ra >> s));

        const auto x = static_cast<std::int64_t>(rng.next());
        const auto y = static_cast<std::int64_t>(rng.next());
        KRI_CHECK(to_ref(i128_mul_i64(x, y)) == static_cast<ref>(x) * static_cast<ref>(y));

        // 10 進表記
        {
            ref v = ra;
            std::string want;
            if (v == 0) {
                want = "0";
            } else {
                const bool neg = v < 0;
                uref u = neg ? (~static_cast<uref>(v) + 1) : static_cast<uref>(v);
                while (u != 0) {
                    want.insert(want.begin(), static_cast<char>('0' + static_cast<int>(u % 10)));
                    u /= 10;
                }
                if (neg) want = "-" + want;
            }
            KRI_CHECK(i128_str(a) == want);
        }
    }
#else
    (void)iters;
#endif
}

}  // namespace kritest

#endif  // KRISITE_TESTS_INT128_ORACLE_HPP

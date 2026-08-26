// Krisite — テスト用の最小ヘルパ
//
// 外部テストフレームワークは使わない（新しい依存を増やさないため。CLAUDE.md ライセンス方針）。
// 失敗を数え、main が非零で終了することで ctest に伝える。
#ifndef KRISITE_TESTS_TEST_UTIL_HPP
#define KRISITE_TESTS_TEST_UTIL_HPP

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "krisite/krisite.hpp"

namespace kritest {

// ---- 検証マクロ -------------------------------------------------------------

inline int g_failures = 0;
inline int g_checks = 0;
inline int g_reported = 0;

inline void report(const char* expr, const char* file, int line, const std::string& extra) {
    ++g_failures;
    if (g_reported < 25) {  // 出力が溢れないよう先頭のみ表示
        ++g_reported;
        std::fprintf(stderr, "FAIL %s:%d: %s\n", file, line, expr);
        if (!extra.empty()) std::fprintf(stderr, "     %s\n", extra.c_str());
    }
}

#define KRI_CHECK(cond)                                                \
    do {                                                               \
        ++::kritest::g_checks;                                         \
        if (!(cond)) ::kritest::report(#cond, __FILE__, __LINE__, {}); \
    } while (0)

#define KRI_CHECK_MSG(cond, msg)                                          \
    do {                                                                  \
        ++::kritest::g_checks;                                            \
        if (!(cond)) ::kritest::report(#cond, __FILE__, __LINE__, (msg)); \
    } while (0)

#define KRI_EQ(a, b) KRI_CHECK_MSG((a) == (b), ::kritest::pair_msg((a), (b)))

template <class A, class B>
inline std::string pair_msg(const A& a, const B& b) {
    return "左辺 = " + std::to_string(a) + " / 右辺 = " + std::to_string(b);
}

/// テスト本体の末尾で呼ぶ。
inline int finish(const char* name) {
    if (g_failures == 0) {
        std::printf("[  OK  ] %s (%d 検証)\n", name, g_checks);
        return 0;
    }
    std::printf("[ FAIL ] %s: %d / %d 件失敗\n", name, g_failures, g_checks);
    return 1;
}

// ---- 決定的な乱数 -----------------------------------------------------------
//
// 実装依存を避けるため splitmix64 を自前で持つ。同じ種で常に同じ列。
struct Rng {
    std::uint64_t state;
    explicit Rng(std::uint64_t seed) noexcept : state(seed) {}

    std::uint64_t next() noexcept {
        state += 0x9E3779B97F4A7C15ull;
        std::uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    /// [0, n) の一様乱数（n > 0）。
    std::uint64_t below(std::uint64_t n) noexcept { return next() % n; }
    /// [lo, hi] の一様乱数。
    std::int64_t range(std::int64_t lo, std::int64_t hi) noexcept {
        const std::uint64_t span = static_cast<std::uint64_t>(hi - lo) + 1;
        return lo + static_cast<std::int64_t>(below(span));
    }
};

// ---- fixed_int の生成 -------------------------------------------------------

/// 一様ランダムなリム列。
template <std::size_t N>
inline krisite::arith::fixed_int<N> rand_full(Rng& rng) noexcept {
    krisite::arith::fixed_int<N> x{};
    for (std::size_t i = 0; i < N; ++i) x.limb[i] = rng.next();
    return x;
}

/// 境界値を重点的に混ぜた分布（SPEC §8.1「境界値を重点的に」）。
template <std::size_t N>
inline krisite::arith::fixed_int<N> rand_biased(Rng& rng) noexcept {
    using krisite::arith::fixed_int;
    using krisite::arith::from_i64;
    fixed_int<N> x{};
    switch (rng.below(10)) {
        case 0:
            return krisite::arith::zero<N>();
        case 1:
            return from_i64<N>(1);
        case 2:
            return from_i64<N>(-1);
        case 3:  // 最大値 2^(64N-1) - 1
            for (std::size_t i = 0; i < N; ++i) x.limb[i] = ~std::uint64_t{0};
            x.limb[N - 1] = ~std::uint64_t{0} >> 1;
            return x;
        case 4:  // 最小値 -2^(64N-1)
            for (std::size_t i = 0; i + 1 < N; ++i) x.limb[i] = 0;
            x.limb[N - 1] = std::uint64_t{1} << 63;
            return x;
        case 5:  // 小さい値
            return from_i64<N>(rng.range(-1000, 1000));
        case 6: {  // 上位リムだけ立てる
            for (std::size_t i = 0; i < N; ++i) x.limb[i] = 0;
            x.limb[N - 1] = rng.next();
            return x;
        }
        case 7: {  // 幅をランダムに切り詰める
            x = rand_full<N>(rng);
            const std::size_t s = static_cast<std::size_t>(rng.below(64 * N));
            return krisite::arith::shr_bits(x, s);
        }
        default:
            return rand_full<N>(rng);
    }
}

// ---- 表示 -------------------------------------------------------------------

template <std::size_t N>
inline std::string to_hex(const krisite::arith::fixed_int<N>& x) {
    std::string s = "0x";
    char buf[17];
    for (std::size_t i = N; i-- > 0;) {
        std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(x.limb[i]));
        s += buf;
    }
    return s;
}

}  // namespace kritest

#endif  // KRISITE_TESTS_TEST_UTIL_HPP

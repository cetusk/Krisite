// Krisite — 固定幅符号付き整数
//
// SPEC-phase0.md §5.1
//   N 個の 64bit リムからなる符号付き整数（2 の補数、リトルエンディアン順）。
//   POD であること。動的確保なし。例外なし。
#ifndef KRISITE_ARITH_FIXED_INT_HPP
#define KRISITE_ARITH_FIXED_INT_HPP

#include <array>
#include <cstddef>
#include <cstdint>

#include "krisite/arith/intrinsics.hpp"
#include "krisite/config.hpp"

#if KRISITE_CHECKED_ARITH
#include <cstdio>
#include <cstdlib>
#endif

namespace krisite::arith {

using u64 = std::uint64_t;
using i64 = std::int64_t;

namespace detail {

#if KRISITE_CHECKED_ARITH
/// 検査失敗時の停止。可変な静的状態を持たない（SPEC §5.1 の制約）。
[[noreturn]] inline void check_fail(const char* expr, const char* msg, const char* file,
                                    int line) noexcept {
    std::fprintf(stderr, "krisite: 算術検査に失敗: %s\n  条件: %s\n  位置: %s:%d\n", msg, expr,
                 file, line);
    std::abort();
}
#endif

}  // namespace detail

// KRISITE_CHECKED_ARITH が 1 のときのみ検査する。NDEBUG とは独立（SPEC §3.4）。
#if KRISITE_CHECKED_ARITH
#define KRISITE_CHECK(cond, msg) \
    (static_cast<bool>(cond)     \
         ? void(0)               \
         : ::krisite::arith::detail::check_fail(#cond, (msg), __FILE__, __LINE__))
#else
#define KRISITE_CHECK(cond, msg) void(0)
#endif

/// 64bit リム N 個の符号付き整数。limb[0] が最下位。
///
/// 集成体（aggregate）のまま置く。コンストラクタを持たせないのは、
/// 述語の内部で大量に作られる一時オブジェクトの初期化コストを避けるため（SPEC §5.1）。
template <std::size_t N>
struct fixed_int {
    static_assert(N >= 1, "リム数は 1 以上");

    static constexpr std::size_t kLimbs = N;
    static constexpr std::size_t kBits = 64 * N;

    std::array<u64, N> limb;

    constexpr u64& operator[](std::size_t i) noexcept { return limb[i]; }
    constexpr const u64& operator[](std::size_t i) const noexcept { return limb[i]; }
};

// ---- 生成 -------------------------------------------------------------------

/// すべてのリムが 0。
template <std::size_t N>
constexpr fixed_int<N> zero() noexcept {
    fixed_int<N> r{};
    for (std::size_t i = 0; i < N; ++i) r.limb[i] = 0;
    return r;
}

/// 符号付き 64bit から符号拡張して構築。
template <std::size_t N>
constexpr fixed_int<N> from_i64(i64 v) noexcept {
    const u64 fill = (v < 0) ? ~u64{0} : u64{0};
    fixed_int<N> r{};
    r.limb[0] = static_cast<u64>(v);
    for (std::size_t i = 1; i < N; ++i) r.limb[i] = fill;
    return r;
}

// ---- 符号と大きさ -----------------------------------------------------------

/// 負かどうか。最上位リムの符号ビットのみを見る。
template <std::size_t N>
constexpr bool is_negative(const fixed_int<N>& x) noexcept {
    return (x.limb[N - 1] >> 63) != 0;
}

template <std::size_t N>
constexpr bool is_zero(const fixed_int<N>& x) noexcept {
    u64 acc = 0;
    for (std::size_t i = 0; i < N; ++i) acc |= x.limb[i];
    return acc == 0;
}

/// -1 / 0 / +1。
///
/// SPEC §5.3: 「sign() は最上位リムの符号ビットだけで決まる」。
/// 負の判定は 1 命令で済み、全リム走査は「非負のとき 0 か正かを分ける」ためだけに必要。
/// 述語のホットパスでは符号が非零であることが多く、その場合も走査は短絡できないが、
/// 負側は即座に返るため最上位リムだけで決着する経路が確保されている。
template <std::size_t N>
constexpr int sign(const fixed_int<N>& x) noexcept {
    if (is_negative(x)) return -1;
    return is_zero(x) ? 0 : 1;
}

/// この値を表すのに必要な符号付きビット数（最小の n。-2^(n-1) <= x <= 2^(n-1)-1）。
/// SPEC §8.4（ビット幅の実測）で使う。
template <std::size_t N>
inline std::size_t min_bits(const fixed_int<N>& x) noexcept {
    const u64 fill = is_negative(x) ? ~u64{0} : u64{0};
    std::size_t redundant = 0;
    std::size_t i = N;
    while (i > 0) {
        --i;
        const u64 w = x.limb[i] ^ fill;  // 符号拡張ぶんを 0 に落とす
        if (w == 0) {
            redundant += 64;
            continue;
        }
        redundant += static_cast<std::size_t>(intr::clz64(w));
        break;
    }
    return 64 * N - redundant + 1;
}

// ---- 幅の変換 ---------------------------------------------------------------

/// 符号拡張して広げる（M >= N）。
template <std::size_t M, std::size_t N>
constexpr fixed_int<M> widen(const fixed_int<N>& x) noexcept {
    static_assert(M >= N, "widen は幅を狭められない。narrow を使うこと");
    const u64 fill = is_negative(x) ? ~u64{0} : u64{0};
    fixed_int<M> r{};
    for (std::size_t i = 0; i < N; ++i) r.limb[i] = x.limb[i];
    for (std::size_t i = N; i < M; ++i) r.limb[i] = fill;
    return r;
}

/// 幅を狭める（M <= N）。切り捨てられる上位が符号拡張ぶんであることを検査する。
///
/// SPEC §3 のビット幅解析が正しければ必ず成功する。ここで停止したということは
/// 実測が理論上界を超えたということであり、CLAUDE.md の「即座に報告」対象。
template <std::size_t M, std::size_t N>
inline fixed_int<M> narrow(const fixed_int<N>& x) noexcept {
    static_assert(M <= N, "narrow は幅を広げられない。widen を使うこと");
#if KRISITE_CHECKED_ARITH
    {
        const u64 fill = is_negative(x) ? ~u64{0} : u64{0};
        bool ok = true;
        for (std::size_t i = M; i < N; ++i) ok = ok && (x.limb[i] == fill);
        // 残す最上位リムの符号ビットが、元の符号と一致していること
        ok = ok && (((x.limb[M - 1] >> 63) != 0) == (fill != 0));
        KRISITE_CHECK(ok, "narrow: 値が目標幅に収まらない（SPEC-phase0.md §3 の上界超過）");
    }
#endif
    fixed_int<M> r{};
    for (std::size_t i = 0; i < M; ++i) r.limb[i] = x.limb[i];
    return r;
}

/// M >= N なら widen、M < N なら narrow。
template <std::size_t M, std::size_t N>
inline fixed_int<M> resize(const fixed_int<N>& x) noexcept {
    if constexpr (M >= N) {
        return widen<M>(x);
    } else {
        return narrow<M>(x);
    }
}

}  // namespace krisite::arith

#endif  // KRISITE_ARITH_FIXED_INT_HPP

// Krisite — 固定幅整数の演算
//
// SPEC-phase0.md §5.2
//   積は「幅が増える」ことを型で表現する。これが本設計の要。
//   述語の式を書いた時点で必要なリム数がコンパイル時に確定し、
//   オーバーフローが型レベルで防がれる。
#ifndef KRISITE_ARITH_OPS_HPP
#define KRISITE_ARITH_OPS_HPP

#include <cstddef>

#include "krisite/arith/fixed_int.hpp"

namespace krisite::arith {

namespace detail {

constexpr std::size_t kmax(std::size_t a, std::size_t b) noexcept {
    return a > b ? a : b;
}

/// 2 の補数のまま加算（オーバーフロー検査なし、mod 2^(64N)）。
template <std::size_t N>
inline fixed_int<N> add_raw(const fixed_int<N>& a, const fixed_int<N>& b) noexcept {
    fixed_int<N> r{};
    unsigned char c = 0;
    for (std::size_t i = 0; i < N; ++i) {
        c = intr::addcarry64(c, a.limb[i], b.limb[i], &r.limb[i]);
    }
    return r;
}

/// 2 の補数のまま減算（オーバーフロー検査なし、mod 2^(64N)）。
template <std::size_t N>
inline fixed_int<N> sub_raw(const fixed_int<N>& a, const fixed_int<N>& b) noexcept {
    fixed_int<N> r{};
    unsigned char c = 0;
    for (std::size_t i = 0; i < N; ++i) {
        c = intr::subborrow64(c, a.limb[i], b.limb[i], &r.limb[i]);
    }
    return r;
}

/// r -= x << (64 * off)（mod 2^(64R)）。off + K <= R を仮定。
template <std::size_t R, std::size_t K>
inline void sub_shifted(fixed_int<R>& r, const fixed_int<K>& x, std::size_t off) noexcept {
    unsigned char bor = 0;
    for (std::size_t j = 0; off + j < R; ++j) {
        const u64 xv = (j < K) ? x.limb[j] : u64{0};
        bor = intr::subborrow64(bor, r.limb[off + j], xv, &r.limb[off + j]);
    }
}

}  // namespace detail

// ---- 加減算 -----------------------------------------------------------------

/// 同幅の加算。オーバーフローを検査する（SPEC §5.2 の「検査付き」）。
template <std::size_t N>
inline fixed_int<N> add(const fixed_int<N>& a, const fixed_int<N>& b) noexcept {
    const fixed_int<N> r = detail::add_raw(a, b);
#if KRISITE_CHECKED_ARITH
    const bool an = is_negative(a), bn = is_negative(b), rn = is_negative(r);
    KRISITE_CHECK(!(an == bn && rn != an), "add: 符号付きオーバーフロー");
#endif
    return r;
}

/// 同幅の減算。オーバーフローを検査する。
template <std::size_t N>
inline fixed_int<N> sub(const fixed_int<N>& a, const fixed_int<N>& b) noexcept {
    const fixed_int<N> r = detail::sub_raw(a, b);
#if KRISITE_CHECKED_ARITH
    const bool an = is_negative(a), bn = is_negative(b), rn = is_negative(r);
    KRISITE_CHECK(!(an != bn && rn != an), "sub: 符号付きオーバーフロー");
#endif
    return r;
}

/// 符号反転。x が最小値のときのみオーバーフローする。
template <std::size_t N>
inline fixed_int<N> neg(const fixed_int<N>& x) noexcept {
#if KRISITE_CHECKED_ARITH
    {
        bool is_min = (x.limb[N - 1] == (u64{1} << 63));
        for (std::size_t i = 0; i + 1 < N; ++i) is_min = is_min && (x.limb[i] == 0);
        KRISITE_CHECK(!is_min, "neg: 最小値の符号反転はオーバーフロー");
    }
#endif
    return detail::sub_raw(zero<N>(), x);
}

/// 幅を 1 リム広げて加算。定義上オーバーフローしない。
template <std::size_t N>
inline fixed_int<N + 1> add_widen(const fixed_int<N>& a, const fixed_int<N>& b) noexcept {
    return detail::add_raw(widen<N + 1>(a), widen<N + 1>(b));
}

/// 幅を 1 リム広げて減算。定義上オーバーフローしない。
template <std::size_t N>
inline fixed_int<N + 1> sub_widen(const fixed_int<N>& a, const fixed_int<N>& b) noexcept {
    return detail::sub_raw(widen<N + 1>(a), widen<N + 1>(b));
}

/// 異幅の加算。結果は max(N,M)+1 リム。
template <std::size_t N, std::size_t M>
inline fixed_int<detail::kmax(N, M) + 1> add_mixed(const fixed_int<N>& a,
                                                   const fixed_int<M>& b) noexcept {
    constexpr std::size_t K = detail::kmax(N, M) + 1;
    return detail::add_raw(widen<K>(a), widen<K>(b));
}

/// 異幅の減算。結果は max(N,M)+1 リム。
template <std::size_t N, std::size_t M>
inline fixed_int<detail::kmax(N, M) + 1> sub_mixed(const fixed_int<N>& a,
                                                   const fixed_int<M>& b) noexcept {
    constexpr std::size_t K = detail::kmax(N, M) + 1;
    return detail::sub_raw(widen<K>(a), widen<K>(b));
}

// ---- 乗算 -------------------------------------------------------------------

/// 符号付き乗算。結果は N+M リム。**必ず収まる**（オーバーフローしない）。
///
/// 最悪値は min * min = 2^(64N-1) * 2^(64M-1) = 2^(64(N+M)-2) で、
/// 符号付き 64(N+M) ビットの上限 2^(64(N+M)-1)-1 未満。
///
/// 実装は「符号なし筆算 + 符号補正」。
///   a_u = a + 2^(64N)[a<0], b_u = b + 2^(64M)[b<0] とおくと
///   a*b ≡ a_u*b_u - (a_u << 64M)[b<0] - (b_u << 64N)[a<0]   (mod 2^(64(N+M)))
template <std::size_t N, std::size_t M>
inline fixed_int<N + M> mul(const fixed_int<N>& a, const fixed_int<M>& b) noexcept {
    constexpr std::size_t R = N + M;
    fixed_int<R> r = zero<R>();

    // 符号なし筆算
    for (std::size_t i = 0; i < N; ++i) {
        u64 carry = 0;
        for (std::size_t j = 0; j < M; ++j) {
            u64 lo = 0, hi = 0;
            intr::mul64(a.limb[i], b.limb[j], &lo, &hi);
            const unsigned char c1 = intr::addcarry64(0, r.limb[i + j], lo, &r.limb[i + j]);
            const unsigned char c2 = intr::addcarry64(0, r.limb[i + j], carry, &r.limb[i + j]);
            // 128bit の総和が 2^128-1 を超えないため hi + c1 + c2 は 64bit に収まる
            carry = hi + static_cast<u64>(c1) + static_cast<u64>(c2);
        }
        u64 k = carry;
        for (std::size_t p = i + M; p < R && k != 0; ++p) {
            k = static_cast<u64>(intr::addcarry64(0, r.limb[p], k, &r.limb[p]));
        }
    }

    // 符号補正
    if (is_negative(a)) detail::sub_shifted(r, b, N);
    if (is_negative(b)) detail::sub_shifted(r, a, M);
    return r;
}

// ---- 比較 -------------------------------------------------------------------

/// -1 / 0 / +1。幅が違ってもよい。
template <std::size_t N, std::size_t M>
inline int cmp(const fixed_int<N>& a, const fixed_int<M>& b) noexcept {
    constexpr std::size_t K = detail::kmax(N, M);
    const fixed_int<K> x = widen<K>(a);
    const fixed_int<K> y = widen<K>(b);
    const bool xn = is_negative(x), yn = is_negative(y);
    if (xn != yn) return xn ? -1 : 1;
    // 同符号なら 2 の補数表現の符号なし比較が順序を保つ
    for (std::size_t i = K; i-- > 0;) {
        if (x.limb[i] != y.limb[i]) return (x.limb[i] < y.limb[i]) ? -1 : 1;
    }
    return 0;
}

template <std::size_t N, std::size_t M>
inline bool equal(const fixed_int<N>& a, const fixed_int<M>& b) noexcept {
    return cmp(a, b) == 0;
}

// ---- シフト -----------------------------------------------------------------

/// 論理左シフト（mod 2^(64N)）。s は 0 以上。
template <std::size_t N>
inline fixed_int<N> shl_bits(const fixed_int<N>& x, std::size_t s) noexcept {
    fixed_int<N> r = zero<N>();
    const std::size_t ls = s / 64, bs = s % 64;
    if (ls >= N) return r;
    for (std::size_t i = N; i-- > ls;) {
        const std::size_t src = i - ls;
        u64 v = x.limb[src] << bs;
        if (bs != 0 && src > 0) v |= x.limb[src - 1] >> (64 - bs);
        r.limb[i] = v;
    }
    return r;
}

/// 算術右シフト。
template <std::size_t N>
inline fixed_int<N> shr_bits(const fixed_int<N>& x, std::size_t s) noexcept {
    const u64 fill = is_negative(x) ? ~u64{0} : u64{0};
    fixed_int<N> r{};
    const std::size_t ls = s / 64, bs = s % 64;
    for (std::size_t i = 0; i < N; ++i) {
        const std::size_t src = i + ls;
        u64 v = (src < N) ? x.limb[src] : fill;
        if (bs != 0) {
            const u64 hi = (src + 1 < N) ? x.limb[src + 1] : fill;
            v = (v >> bs) | (hi << (64 - bs));
        }
        r.limb[i] = v;
    }
    return r;
}

// ---- 行列式 -----------------------------------------------------------------
//
// 述語で繰り返し現れる 2x2 / 3x3 の行列式。列ごとに幅が違う場合は、
// 呼び出し側で最大幅に widen してから渡す（正しさ優先。SPEC §0）。

/// 2x2 行列式の 1 行。**`det2` に同じ型の引数を 4 つ並べさせないための器**です。
///
/// > 位置引数が多く型が同じ小関数は、規律ではなく設計で守ること（`CLAUDE.md`）。
///
/// 旧 `det2(a, b, c, d)` は取り違えてもコンパイルが通り、**実際に 2 度間違えました**
/// （`plane_from_triangle` と `plane_from_edge`。どちらも外積の成分）。
/// 行に区切ると、行をまたぐ取り違えは書けなくなります。
template <std::size_t L>
struct row2 {
    fixed_int<L> a, b;
};

/// 3 成分ベクトル。**外積の引数を 6 つ並べさせないための器**です。
template <std::size_t L>
struct vec3 {
    fixed_int<L> x, y, z;
};

/// 行列式
///
///     | r0.a  r0.b |
///     | r1.a  r1.b |  =  r0.a * r1.b - r0.b * r1.a
///
/// |結果| <= 2^(128L-1) なので 2L+1 リムに必ず収まる。
template <std::size_t L>
inline fixed_int<2 * L + 1> det2(const row2<L>& r0, const row2<L>& r1) noexcept {
    return sub_widen(mul(r0.a, r1.b), mul(r0.b, r1.a));
}

/// 外積 $u \times v$。
///
/// **成分ごとに `det2` を書かないでください。** 外積は本プロジェクトで最も頻出する
/// 構成（法線・辺平面・平行判定）で、取り違えの実績が 2 件あります。
template <std::size_t L>
inline vec3<2 * L + 1> cross(const vec3<L>& u, const vec3<L>& v) noexcept {
    return {det2(row2<L>{u.y, u.z}, row2<L>{v.y, v.z}),   // u.y*v.z - u.z*v.y
            det2(row2<L>{u.z, u.x}, row2<L>{v.z, v.x}),   // u.z*v.x - u.x*v.z
            det2(row2<L>{u.x, u.y}, row2<L>{v.x, v.y})};  // u.x*v.y - u.y*v.x
}

/// 3x3 行列式（余因子展開）。
///
///   det = m00*(m11*m22 - m12*m21) - m01*(m10*m22 - m12*m20) + m02*(m10*m21 - m11*m20)
///
/// |det| <= 3 * 2^(64L-1) * 2^(128L-1) = 3*2^(192L-2) < 2^(192L) なので
/// 3L+1 リム（192L+64 ビット）に十分収まる。
template <std::size_t L>
inline fixed_int<3 * L + 1> det3(const fixed_int<L> m[3][3]) noexcept {
    constexpr std::size_t R = 3 * L + 1;
    const fixed_int<R> t0 =
        mul(m[0][0], det2(row2<L>{m[1][1], m[1][2]}, row2<L>{m[2][1], m[2][2]}));
    const fixed_int<R> t1 =
        mul(m[0][1], det2(row2<L>{m[1][0], m[1][2]}, row2<L>{m[2][0], m[2][2]}));
    const fixed_int<R> t2 =
        mul(m[0][2], det2(row2<L>{m[1][0], m[1][1]}, row2<L>{m[2][0], m[2][1]}));
    return add(sub(t0, t1), t2);
}

}  // namespace krisite::arith

#endif  // KRISITE_ARITH_OPS_HPP

// Krisite — 64bit リム演算のプリミティブ
//
// SPEC-phase0.md §5.3
//   - x86-64: _addcarry_u64 / _subborrow_u64 / _mulx_u64
//   - 可搬フォールバック: unsigned __int128
//   - 切替はこのヘッダに閉じ込め、上位からは見えないようにする
//
// このヘッダの関数はすべて noexcept・動的確保なし・グローバル状態なし。
#ifndef KRISITE_ARITH_INTRINSICS_HPP
#define KRISITE_ARITH_INTRINSICS_HPP

#include <cstdint>

// ---- 経路の判定 -------------------------------------------------------------
//
// KRISITE_HAS_INT128 : unsigned __int128 が使えるか
// KRISITE_HAS_X86_ADX: _addcarry_u64 / _subborrow_u64 が使えるか
// KRISITE_HAS_MSVC_64: MSVC の _umul128 / _addcarry_u64 が使えるか

// KRISITE_NO_INT128 / KRISITE_NO_INTRINSICS を定義すると、その経路を無効にできる。
// 手元に無い環境（MSVC ARM64 など）の経路をテストするための逃げ道。
#if defined(__SIZEOF_INT128__) && !defined(KRISITE_NO_INT128)
#define KRISITE_HAS_INT128 1
#else
#define KRISITE_HAS_INT128 0
#endif

#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#if defined(_M_X64)
#define KRISITE_HAS_MSVC_X64 1
#else
#define KRISITE_HAS_MSVC_X64 0
#endif
#else
#define KRISITE_HAS_MSVC_X64 0
#endif

#if (defined(__x86_64__) || defined(_M_X64)) && !defined(KRISITE_NO_INTRINSICS)
#define KRISITE_HAS_X86_64 1
#if defined(__GNUC__) || defined(__clang__)
#include <x86intrin.h>
#endif
#else
#define KRISITE_HAS_X86_64 0
#endif

// SPEC §5.3: 経路は 3 つある。
//   (1) x86-64 の GCC/Clang  : _addcarry_u64 / _subborrow_u64 + __int128
//   (2) MSVC x64             : <intrin.h> の _addcarry_u64 / _subborrow_u64 / _umul128
//                              （MSVC には unsigned __int128 が無いので (3) では不足）
//   (3) 可搬                 : unsigned __int128、または 32bit 分割の筆算
// (3) の最終段（32bit 分割）は MSVC ARM64 のように __int128 も _umul128 も無い環境の
// ための保険。#error で構成が落ちることは無い。

namespace krisite::arith::intr {

using u64 = std::uint64_t;

/// out = a + b + cin、返り値は桁上がり（0 または 1）。
inline unsigned char addcarry64(unsigned char cin, u64 a, u64 b, u64* out) noexcept {
#if KRISITE_HAS_X86_64 && (defined(__GNUC__) || defined(__clang__))
    unsigned long long r = 0;
    const unsigned char c = _addcarry_u64(cin, static_cast<unsigned long long>(a),
                                          static_cast<unsigned long long>(b), &r);
    *out = static_cast<u64>(r);
    return c;
#elif KRISITE_HAS_MSVC_X64
    unsigned __int64 r = 0;
    const unsigned char c = _addcarry_u64(cin, a, b, &r);
    *out = static_cast<u64>(r);
    return c;
#else
    // 可搬経路。__int128 に頼らず 64bit だけで桁上がりを判定する。
    const u64 s = a + b;
    const unsigned char c0 = static_cast<unsigned char>(s < a);
    const u64 t = s + static_cast<u64>(cin);
    const unsigned char c1 = static_cast<unsigned char>(t < s);
    *out = t;
    return static_cast<unsigned char>(c0 | c1);
#endif
}

/// out = a - b - bin、返り値は借り（0 または 1）。
inline unsigned char subborrow64(unsigned char bin, u64 a, u64 b, u64* out) noexcept {
#if KRISITE_HAS_X86_64 && (defined(__GNUC__) || defined(__clang__))
    unsigned long long r = 0;
    const unsigned char c = _subborrow_u64(bin, static_cast<unsigned long long>(a),
                                           static_cast<unsigned long long>(b), &r);
    *out = static_cast<u64>(r);
    return c;
#elif KRISITE_HAS_MSVC_X64
    unsigned __int64 r = 0;
    const unsigned char c = _subborrow_u64(bin, a, b, &r);
    *out = static_cast<u64>(r);
    return c;
#else
    const u64 d = a - b;
    const unsigned char b0 = static_cast<unsigned char>(a < b);
    const u64 t = d - static_cast<u64>(bin);
    const unsigned char b1 = static_cast<unsigned char>(d < static_cast<u64>(bin));
    *out = t;
    return static_cast<unsigned char>(b0 | b1);
#endif
}

/// 符号なし 64x64 -> 128。lo に下位、hi に上位を書く。
inline void mul64(u64 a, u64 b, u64* lo, u64* hi) noexcept {
#if KRISITE_HAS_INT128
    const unsigned __int128 p = static_cast<unsigned __int128>(a) * b;
    *lo = static_cast<u64>(p);
    *hi = static_cast<u64>(p >> 64);
#elif KRISITE_HAS_MSVC_X64
    unsigned __int64 h = 0;
    *lo = static_cast<u64>(_umul128(a, b, &h));
    *hi = static_cast<u64>(h);
#else
    // 32bit 分割の筆算。__int128 も _umul128 も無い環境（MSVC ARM64 等）向け。
    const u64 a0 = a & 0xFFFFFFFFu, a1 = a >> 32;
    const u64 b0 = b & 0xFFFFFFFFu, b1 = b >> 32;
    const u64 p00 = a0 * b0;
    const u64 p01 = a0 * b1;
    const u64 p10 = a1 * b0;
    const u64 p11 = a1 * b1;
    const u64 mid = (p00 >> 32) + (p01 & 0xFFFFFFFFu) + (p10 & 0xFFFFFFFFu);
    *lo = (mid << 32) | (p00 & 0xFFFFFFFFu);
    *hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
#endif
}

/// x の先頭ゼロビット数（x == 0 のとき 64）。
inline int clz64(u64 x) noexcept {
    if (x == 0) return 64;
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(static_cast<unsigned long long>(x));
#elif defined(_MSC_VER)
    unsigned long idx = 0;
    _BitScanReverse64(&idx, x);
    return 63 - static_cast<int>(idx);
#else
    int n = 0;
    for (int i = 63; i >= 0; --i) {
        if ((x >> i) & 1u) break;
        ++n;
    }
    return n;
#endif
}

}  // namespace krisite::arith::intr

#endif  // KRISITE_ARITH_INTRINSICS_HPP

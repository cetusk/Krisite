// Krisite — ビット幅の導出
//
// SPEC-phase0.md §3.1 の表をそのまま constexpr に落としたもの。
// **ここ以外にリム数の数値リテラルを書かないこと**（SPEC §3.4, §6）。
//
// 入力座標が符号付き b ビット（-2^(b-1) <= coord <= 2^(b-1)-1）のとき:
//
//   | 量                    | 式                       | ビット幅 |
//   |-----------------------|--------------------------|---------|
//   | 座標差 p2-p1          |                          | b + 1   |
//   | 法線 N = u x v        | 差の外積                 | 2b + 3  |
//   | オフセット d = -N.p1  | 平面は N.x + d = 0       | 3b + 5  |
//   | w = det(N1,N2,N3)     | 3 平面交点の同次成分     | 6b + 12 |
//   | x,y,z（Cramer）       |                          | 7b + 14 |
//   | side の被符号値       | N.V + d w                | 9b + 20 |
//   | cmp_h の被符号値      | w2 x1 - w1 x2            | 13b + 27|
//
// **Phase 0 で最大幅を要するのは cmp_h**（§3.3）。side ではない。
// b = 21 で 300 ビット = 5 リム、b = 26 で 365 ビット = 6 リム。
//
// なお §3.1 の上界はいずれも厳密最小ではなく、余裕を持たせた保守的な値である。
// 実測が上界を超えないことは tests/arith/test_bitwidth.cpp で確認する（SPEC §8.4）。
// **乱択の実測値でリム数を減らしてはならない**（SPEC §3.4）。
#ifndef KRISITE_GEOM_WIDTHS_HPP
#define KRISITE_GEOM_WIDTHS_HPP

#include <cstddef>

#include "krisite/config.hpp"

namespace krisite::geom {

/// bits ビットを収めるのに必要な 64bit リム数。
constexpr std::size_t limbs_for(std::size_t bits) noexcept {
    return (bits + 63) / 64;
}

/// リム数の最大（行列式で列幅をそろえるときに使う）。
constexpr std::size_t max_limbs(std::size_t a, std::size_t b) noexcept {
    return a > b ? a : b;
}

namespace bits {

inline constexpr std::size_t b = kCoordBits;

inline constexpr std::size_t kCoord = b;             ///< 入力座標
inline constexpr std::size_t kDiff = b + 1;          ///< 座標差
inline constexpr std::size_t kNormal = 2 * b + 3;    ///< 法線成分
inline constexpr std::size_t kOffset = 3 * b + 5;    ///< 平面オフセット d = -N・p1
inline constexpr std::size_t kHomoW = 6 * b + 12;    ///< 構成点の w
inline constexpr std::size_t kHomoXyz = 7 * b + 14;  ///< 構成点の x,y,z
inline constexpr std::size_t kSide = 9 * b + 20;     ///< side() の被符号値
inline constexpr std::size_t kOrient3d = 3 * b + 5;  ///< orient3d の行列式（差分形）
inline constexpr std::size_t kOrient2d = 2 * b + 3;  ///< orient2d の行列式（差分形）

/// cmp_h() の被符号値 w2*x1 - w1*x2。SPEC §3.1「同次点の比較述語」:
///
///   w2*x1 : (6b+12) + (7b+14) = 13b + 26
///   差分   : 13b + 27
///
/// b = 21 で 300 ビット（5 リム）。side() の 209 ビット（4 リム）より広い。
/// すなわち **Phase 0 で最大幅を要求するのは side() ではなく cmp_h()** である（§3.3）。
inline constexpr std::size_t kCmpH = 13 * b + 27;

}  // namespace bits

namespace limbs {

inline constexpr std::size_t kCoord = limbs_for(bits::kCoord);
inline constexpr std::size_t kDiff = limbs_for(bits::kDiff);
inline constexpr std::size_t kNormal = limbs_for(bits::kNormal);
inline constexpr std::size_t kOffset = limbs_for(bits::kOffset);
inline constexpr std::size_t kHomoW = limbs_for(bits::kHomoW);
inline constexpr std::size_t kHomoXyz = limbs_for(bits::kHomoXyz);
inline constexpr std::size_t kSide = limbs_for(bits::kSide);
inline constexpr std::size_t kOrient3d = limbs_for(bits::kOrient3d);
inline constexpr std::size_t kOrient2d = limbs_for(bits::kOrient2d);
inline constexpr std::size_t kCmpH = limbs_for(bits::kCmpH);

/// Phase 0 の述語が要求する最大リム数（SPEC §3.3 の表の最右列）。
/// b = 21 → 5、b = 26 → 6。
inline constexpr std::size_t kMaxPredicate = max_limbs(kSide, kCmpH);

}  // namespace limbs

// 実際の計算経路が §3 の上界を下回っていないことをコンパイル時に確かめる。
static_assert(64 * limbs::kSide >= bits::kSide, "kSide のリム数不足");
static_assert(64 * limbs::kCmpH >= bits::kCmpH, "kCmpH のリム数不足");
static_assert(64 * (limbs::kHomoW + limbs::kHomoXyz + 1) >= bits::kCmpH,
              "cmp_h の計算幅が理論上界を下回っている");

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_WIDTHS_HPP

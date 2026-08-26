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

// ---- Phase 1 で追加した量（SPEC-phase1.md §7）--------------------------------

/// 平面の比例判定に使う 2x2 小行列式 p_i*q_j - p_j*q_i（SPEC-phase1 §3.1）。
///
/// 係数は法線 a,b,c が 2b+3、オフセット d が 3b+5。最大の積は 法線 x オフセットで
/// (2b+3)+(3b+5)-1 = 5b+7、差でさらに +1 して 5b+8。仕様の上界 5b+9 を採る。
inline constexpr std::size_t kPlaneMinor = 5 * b + 9;

/// 平面の全順序に使う交差乗算（SPEC-phase1 §3.1）。
///
/// 実際に到達しうる最大は kPlaneMinor と同じ 5b+8 である。先頭非零成分 i は
/// 法線に限られるため（i が d なら a=b=c=0 で、それは退化三角形の平面）。
/// 仕様の上界 6b+11 は 24 ビット（b=21）保守的だが、安全側なのでそのまま採る。
inline constexpr std::size_t kPlaneOrder = 6 * b + 11;

/// 軸平行平面のオフセット d = -coord（SPEC-phase1 §3.2）。
/// |coord| <= 2^(b-1) なので b+1。kOffset に自明に収まる。
inline constexpr std::size_t kAxisOffset = b + 1;

/// 三角形 1 枚ぶんの符号付き体積 x6 = det(a, b, c)。
/// |det| <= 6*2^(3(b-1)) < 2^(3b) なので 3b+1。
inline constexpr std::size_t kTetraVolume6 = 3 * b + 1;

/// 入力メッシュの三角形数の上限（2 の冪の指数）。符号付き体積の総和幅を決める。
/// SPEC-phase1 §2.1 は「三角形数 100 未満」を想定しているので 2^16 で十分に余裕がある。
inline constexpr std::size_t kMaxTrianglesLog2 = 16;

/// 入力メッシュの符号付き体積 x6 = Σ det(v_i, v_j, v_k)（SPEC-phase1 §3.4 の向き検査）。
///
/// 1 枚ぶんが kTetraVolume6 = 3b+1 ビット、三角形 2^16 個までの総和で 3b+1+16。
///
/// **出力メッシュの体積はこれでは測れません。** 構成点が有理数で共通分母が
/// 三角形数に比例して伸びるため、GMP が要ります（SPEC-phase1 §10.3）。
inline constexpr std::size_t kInputVolume6 = kTetraVolume6 + kMaxTrianglesLog2;

// レイキャスト（軸平行レイ）は新しい述語を必要としません。投影後の内外判定は
// orient2d（2b+3）、前方交差の判定は side(plane, IPoint)（3b+6）と法線成分の符号で
// 決まります。SPEC-phase1 §6.1。

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

// Phase 1（SPEC-phase1.md §7）
inline constexpr std::size_t kPlaneMinor = limbs_for(bits::kPlaneMinor);
inline constexpr std::size_t kPlaneOrder = limbs_for(bits::kPlaneOrder);
inline constexpr std::size_t kInputVolume6 = limbs_for(bits::kInputVolume6);

/// Phase 0 の述語が要求する最大リム数（SPEC §3.3 の表の最右列）。
/// b = 21 → 5、b = 26 → 6。
inline constexpr std::size_t kMaxPredicate = max_limbs(kSide, kCmpH);

}  // namespace limbs

// 平面係数を共通幅で扱うために kOffset >= kNormal が要る。3b+5 > 2b+3 なので常に成立。
static_assert(limbs::kOffset >= limbs::kNormal, "kOffset は kNormal 以上のはず");
static_assert(64 * limbs::kPlaneMinor >= bits::kPlaneMinor, "kPlaneMinor のリム数不足");
static_assert(64 * limbs::kPlaneOrder >= bits::kPlaneOrder, "kPlaneOrder のリム数不足");
static_assert(64 * limbs::kInputVolume6 >= bits::kInputVolume6, "kInputVolume6 のリム数不足");
// 軸平行平面のオフセットは kOffset に収まること
static_assert(bits::kAxisOffset <= bits::kOffset, "軸平行平面のオフセットが kOffset を超える");

// 実際の計算経路が §3 の上界を下回っていないことをコンパイル時に確かめる。
static_assert(64 * limbs::kSide >= bits::kSide, "kSide のリム数不足");
static_assert(64 * limbs::kCmpH >= bits::kCmpH, "kCmpH のリム数不足");
static_assert(64 * (limbs::kHomoW + limbs::kHomoXyz + 1) >= bits::kCmpH,
              "cmp_h の計算幅が理論上界を下回っている");

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_WIDTHS_HPP

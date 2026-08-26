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

/// 同次点を投影した 2D 向き（SPEC-phase1 §6.1 のレイキャスト）。
///
///   O = (b.y-a.y)*(p.z - a.z*p.w) - (b.z-a.z)*(p.y - a.y*p.w)
///   実座標での向き = O / p.w なので、符号は sign(O) * sign(p.w)
///
///   |p.z| < 2^(7b+13)、|a.z*p.w| < 2^(7b+10) → 差は 7b+15 ビット
///   |b.y-a.y| は b+1 ビット → 積が 8b+15、差でさらに +1 して 8b+16
///   安全側に +1 して 8b+17（b=21 で 185 ビット / 3 リム、b=26 で 225 / 4）
///
/// **キャスト元が入力頂点（IPoint）なら既存の orient2d（2b+3）で済みます。**
/// 同次点にフォールバックしたときだけこの幅が要ります。
inline constexpr std::size_t kOrient2dH = 8 * b + 17;

// レイの前方交差の判定は side(plane, ·) と法線成分の符号で決まるので、
// 新しい述語は要りません。キャスト元は ∂B 上に無いことを保証するので
// N・p + d != 0 が保証され、この判定は退化しません（SPEC-phase1 §6.1）。

/// cmp_h() の被符号値 w2*x1 - w1*x2。SPEC §3.1「同次点の比較述語」:
///
///   w2*x1 : (6b+12) + (7b+14) = 13b + 26
///   差分   : 13b + 27
///
/// b = 21 で 300 ビット（5 リム）。side() の 209 ビット（4 リム）より広い。
/// すなわち **Phase 0 で最大幅を要求するのは side() ではなく cmp_h()** である（§3.3）。
inline constexpr std::size_t kCmpH = 13 * b + 27;

/// 断片の 2 頂点の**中点**に対する述語（SPEC-phase1 §6.1 の代表点フォールバック）。
///
/// 分類の代表点として頂点を使う設計は、CP2 のケース 2 で破れました。断片の 4 辺の
/// うち 3 辺が相手の面に載る配置があり、**4 頂点すべてが $\partial B$ 上**になります。
/// 相対内部の点が要るので、対角線の中点を使います。
///
/// $v_i = (X_i : W_i)$ の中点は
/// $$ m = (X_0 W_1 + X_1 W_0 \;:\; 2 W_0 W_1). $$
///
/// **$m$ は構成しません。** `side_value` も `orient2d_h_value` も同次座標について
/// **線形**なので、被符号値 $F$ について
/// $$ F(m) = W_1 F(v_0) + W_0 F(v_1) $$
/// が厳密に成り立ちます。$m$ の $w = 2W_0W_1$ の符号は $\mathrm{sign}(W_0)\mathrm{sign}(W_1)$
/// です。 これで新しい**点の型**を作らずに済み、`kHomoXyz` / `kHomoW` を広げる必要もありません。
///
///   side:       $(6b{+}12) + (9b{+}20) = 15b{+}32$、2 項の和で $15b{+}33$
///   orient2d_h: $(6b{+}12) + (8b{+}17) = 14b{+}29$、2 項の和で $14b{+}30$
///
/// b = 21 で 348 ビット（6 リム）/ 324 ビット（6 リム）。
inline constexpr std::size_t kMidSide = kHomoW + kSide + 1;
inline constexpr std::size_t kMidOrient2dH = kHomoW + kOrient2dH + 1;

/// 断片の 3 頂点の**重心**に対する述語（同上のさらなるフォールバック）。
///
/// 中点は対角線を要求するので**三角形の断片では使えません。** 実際に CP2.5 の
/// ケース 2T で、3 頂点すべてが $\partial B$ 上の三角形断片が出ました。
///
/// $v_i = (X_i : W_i)$ の重心は、$W_i > 0$ に正規化したうえで
/// $$ c = (W_1W_2 X_0 + W_0W_2 X_1 + W_0W_1 X_2 \;:\; 3 W_0W_1W_2). $$
///
/// **凸多角形の 3 頂点が張る三角形の内部は、多角形の相対内部に含まれます。**
/// したがって共線でない 3 頂点を選べば必ず境界から外れます。
///
/// ここでも $c$ は構成せず、線形性から
/// $$ F(c) = W_1W_2 F(v_0) + W_0W_2 F(v_1) + W_0W_1 F(v_2) $$
/// で評価します。符号は $\mathrm{sign}(W_0)\mathrm{sign}(W_1)\mathrm{sign}(W_2)$ 倍。
///
///   side:       $2(6b{+}12) + (9b{+}20) = 21b{+}44$、3 項の和で $21b{+}46$
///   orient2d_h: $2(6b{+}12) + (8b{+}17) = 20b{+}41$、3 項の和で $20b{+}43$
///
/// b = 21 で 487 ビット / 463 ビット（ともに 8 リム）。
inline constexpr std::size_t kTriSide = 2 * kHomoW + kSide + 2;
inline constexpr std::size_t kTriOrient2dH = 2 * kHomoW + kOrient2dH + 2;

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
inline constexpr std::size_t kOrient2dH = limbs_for(bits::kOrient2dH);
inline constexpr std::size_t kMidSide = limbs_for(bits::kMidSide);
inline constexpr std::size_t kMidOrient2dH = limbs_for(bits::kMidOrient2dH);
inline constexpr std::size_t kTriSide = limbs_for(bits::kTriSide);
inline constexpr std::size_t kTriOrient2dH = limbs_for(bits::kTriOrient2dH);

/// Phase 0 の述語が要求する最大リム数（SPEC §3.3 の表の最右列）。
/// b = 21 → 5、b = 26 → 6。
inline constexpr std::size_t kMaxPredicate = max_limbs(kSide, kCmpH);

}  // namespace limbs

// 平面係数を共通幅で扱うために kOffset >= kNormal が要る。3b+5 > 2b+3 なので常に成立。
static_assert(limbs::kOffset >= limbs::kNormal, "kOffset は kNormal 以上のはず");
static_assert(64 * limbs::kPlaneMinor >= bits::kPlaneMinor, "kPlaneMinor のリム数不足");
static_assert(64 * limbs::kPlaneOrder >= bits::kPlaneOrder, "kPlaneOrder のリム数不足");
static_assert(64 * limbs::kInputVolume6 >= bits::kInputVolume6, "kInputVolume6 のリム数不足");
static_assert(64 * limbs::kOrient2dH >= bits::kOrient2dH, "kOrient2dH のリム数不足");
static_assert(64 * limbs::kMidSide >= bits::kMidSide, "kMidSide のリム数不足");
static_assert(64 * limbs::kMidOrient2dH >= bits::kMidOrient2dH, "kMidOrient2dH のリム数不足");
static_assert(64 * limbs::kTriSide >= bits::kTriSide, "kTriSide のリム数不足");
static_assert(64 * limbs::kTriOrient2dH >= bits::kTriOrient2dH, "kTriOrient2dH のリム数不足");
// 軸平行平面のオフセットは kOffset に収まること
static_assert(bits::kAxisOffset <= bits::kOffset, "軸平行平面のオフセットが kOffset を超える");

// 実際の計算経路が §3 の上界を下回っていないことをコンパイル時に確かめる。
static_assert(64 * limbs::kSide >= bits::kSide, "kSide のリム数不足");
static_assert(64 * limbs::kCmpH >= bits::kCmpH, "kCmpH のリム数不足");
static_assert(64 * (limbs::kHomoW + limbs::kHomoXyz + 1) >= bits::kCmpH,
              "cmp_h の計算幅が理論上界を下回っている");

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_WIDTHS_HPP

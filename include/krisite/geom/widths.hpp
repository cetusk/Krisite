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

/// 平面とセルの閉領域の交差判定（SPEC-phase2 §2.3）。
///
/// セル $C$ を平面 $P$ で分割するのは $P \cap \overline{C} \neq \emptyset$ のときだけ。
/// 判定は `N・x + d` の最小値と最大値の符号で行う。$N$ の各成分の符号に応じて
/// `lo` / `hi` を選べば、8 隅を回らずに 2 回の内積で出る。
///
///   セル座標          b+1（上限 2^(b-1) が kCoordMax を超えるため。§3.2）
///   N_i * x_i         (2b+3) + (b+1) = 3b+4
///   N・x（3 項の和）   3b+6
///   N・x + d          **3b+7**
///
/// b = 21 で 70 ビット / 2 リム、b = 26 で 85 ビット / 2 リム。
inline constexpr std::size_t kPlaneAabb = 3 * b + 7;

/// **有理点の 1 軸を、整数の境界と比べる**（`DESIGN-phase5-hotspots.md` §11 の D-1 / D-2）。
///
/// $$\mathrm{sign}(x/w - c) = \mathrm{sign}(w)\cdot\mathrm{sign}(x - c\,w)$$
///
/// **`cmp_h` で代用できますが、高くつきます。** `cmp_h` は同次点どうしの比較なので
/// **乗算 2 回・$13b{+}27$ ビット（5 リム）**です。**相手の $w$ が 1 だと分かっていれば
/// 乗算 1 回で済みます**（実測で 2.34 分の 1。`BENCH.md`）。
///
///   x         kHomoXyz = 7b+14
///   c         b+1（kAxisOffset。セル境界も渡せるように）
///   c*w       (b+1) + (6b+12) = 7b+13
///   x - c*w   max(7b+14, 7b+13) + 1 = **7b+15**
///
/// b = 21 で 162 ビット / **3 リム**、b = 26 で 197 ビット / 4 リム。
inline constexpr std::size_t kAxisIntCmp = 7 * b + 15;

/// 三角形とセルの箱の**厳密な**交差判定（分離軸定理）。
/// `DESIGN-phase5-hotspots.md` §10。**まだ実装していません。導出だけ先に置きます。**
///
/// **いまの割り当ては三角形の AABB とセルの重なりだけを見ています。**
/// **斜めの三角形の AABB は、セルが細かくなるほど過剰になります**
/// （実測: 深度 8 で割り当ての 95.4%、深度 9 で 98.7% が実際には交わらない）。
///
/// **半分を避けるため、すべて 2 倍して整数のまま扱います。**
///
///   セル座標 lo, hi        b+1（上限 2^(b-1) が kCoordMax を超えるため。§3.2）
///   中心 x2  c2 = lo + hi  b+2
///   半径 x2  he = hi - lo  b+2
///   頂点 x2  V = 2v - c2   b+3
///   辺       E = V_i - V_j b+4
///
/// **(1) 辺 × 単位ベクトルの 9 軸。** 軸 A = E × e_k は成分が 2 本だけ非零で |A| <= 2^(b+4)。
///
///   p = A・V（2 項）        (b+4) + (b+3) + 1 = 2b+8
///   r = Σ he|A|（3 項）     (b+2) + (b+4) + 2 = 2b+8
///   → **2b+9**（符号を含む余裕）
///
/// b = 21 で 51 ビット / **1 リム**。b = 26 で 61 ビット / 1 リム。
inline constexpr std::size_t kSatEdgeAxis = 2 * b + 9;

/// **(2) 三角形の法線 1 軸。** N = E0 × E1 で |N| <= 2^(2b+9)。
///
///   p = N・V（3 項）        (2b+9) + (b+3) + 2 = 3b+14
///   r = Σ he|N|（3 項）     (b+2) + (2b+9) + 2 = 3b+13
///   → **3b+14**
///
/// b = 21 で 77 ビット / **2 リム**。b = 26 で 92 ビット / 2 リム。
///
/// > **★ ただし、この軸は既存の `plane_crosses_box`（kPlaneAabb = 3b+7）で代用できます。**
/// > 三角形の支持平面は `PlaneTable` に intern 済みで、係数の幅は 2b+3 に正規化されています。
/// > **法線を作り直すより幅が狭く、GMP 差分テストも既にあります**（`test_plane_box_gmp`）。
/// > **この定数は「法線を自前で作る場合」の上界として置きます。**
inline constexpr std::size_t kSatNormalAxis = 3 * b + 14;

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

/// **中点系・重心系の述語は Phase 3 の段 0 で削除しました**（`SPEC-phase3.md` §2.1）。
///
/// 代表点は「頂点 → 対角線の中点 → 3 頂点の重心」の 3 段で選んでいましたが、
/// 2 段目と 3 段目は**点を組み合わせる**操作で、`kTriSide` が $21b{+}46$
/// （b=21 で 487 ビット / 8 リム）に達していました。
///
/// 段 0 は EMBER §4.4 に倣い、**構成そのものを平面ベースに戻します。**
/// 主経路が軸平行 2 枚 + 支持平面の交点（`kAxisPointXyz`）、
/// 予備経路が辺平面を内側へずらした交点（`kOffsetShifted`）で、
/// **どちらも一般の 3 平面交点の幅に収まります。**
///
/// これで Phase 0 / 1 の最大幅は `kCmpH`（13b+27）に戻ります。

// ---- Phase 3 で追加した量（SPEC-phase3.md §3.1 / §2.1）----------------------

/// **辺平面の法線** $N_e = (p_2 - p_1) \times e_k$（`SPEC-phase3.md` §3.1.2 の案 B）。
///
/// 相手が軸方向の単位ベクトルなので**乗算が起きません。** 成分は座標差そのもの
/// （符号反転を含む）なので、幅は座標差と同じ $b+1$ です。
///
/// $N_e = N_s \times (p_2 - p_1)$（案 A）だと $3b+5$ になり、`side` が $13b+28$、
/// b=26 で 6 リムに増えます。**構成方式の選択が決定的です。**
inline constexpr std::size_t kEdgeNormal = b + 1;

/// **辺平面のオフセット** $d_e = -N_e \cdot p_1$。
///
///   $N_{e,i} \cdot p_{1,i}$ : $(b+1) + b = 2b+1$
///   3 項の和（うち 1 つは 0）  : $2b+2$
inline constexpr std::size_t kEdgeOffset = 2 * b + 2;

/// 辺平面を**内側へずらした**オフセット $d' = d - i\,\sigma$（`SPEC-phase3.md` §2.1.2）。
///
/// $\sigma$ は内側の符号、$i \ge 1$ は小さなオフセット。**`kOffset` と同じリム数に
/// 収める**ので、ずらし幅の上限を $2^{kShiftLog2}$ に制限します。
inline constexpr std::size_t kShiftLog2 = 16;
inline constexpr std::size_t kOffsetShifted = kOffset + 1;

/// 軸平行 2 枚 + 支持平面の交点（`SPEC-phase3.md` §2.1.1 の主経路）。
///
/// 軸平行平面は $N$ が単位ベクトル、$d$ が $b+1$。支持平面は $N_s$ が $2b+3$、
/// $d_s$ が $3b+5$。3x3 行列式の各項は行と列を 1 つずつ選ぶので、
///
///   $w$   : 単位 x 単位 x 法線 = $2b+3$、和で $2b+6$
///   $xyz$ : 最大項は 単位 x 単位 x $d_s$ = $3b+5$、和で $3b+8$
///
/// **一般の構成点（$6b+12$ / $7b+14$）より小さいので、`HPointD` にそのまま入ります。**
/// 新しい点の型は要りません。
inline constexpr std::size_t kAxisPointW = 2 * b + 6;
inline constexpr std::size_t kAxisPointXyz = 3 * b + 8;

/// 2 つの平面の法線の内積の符号（`SPEC-phase3.md` §7 の手順で導出）。
///
/// 面が重なっているとき、**どちら向きに重なっているか**の判定に使います。
///
///   $N_1 \cdot N_2$ の各項 : $(2b+3) + (2b+3) = 4b+6$
///   3 項の和               : $4b+8$
///
/// b = 21 で 92 ビット / 2 リム、b = 26 で 112 ビット / 2 リム。
inline constexpr std::size_t kNormalDot = 4 * b + 8;

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
/// セル境界の座標（b+1 ビット）。`IPoint` では表せない上端を含む（§3.2）
inline constexpr std::size_t kAxisOffset = limbs_for(bits::kAxisOffset);
inline constexpr std::size_t kCmpH = limbs_for(bits::kCmpH);

// Phase 1（SPEC-phase1.md §7）
inline constexpr std::size_t kPlaneMinor = limbs_for(bits::kPlaneMinor);
inline constexpr std::size_t kPlaneOrder = limbs_for(bits::kPlaneOrder);
inline constexpr std::size_t kInputVolume6 = limbs_for(bits::kInputVolume6);
inline constexpr std::size_t kOrient2dH = limbs_for(bits::kOrient2dH);

// Phase 2（SPEC-phase2.md §7）
inline constexpr std::size_t kPlaneAabb = limbs_for(bits::kPlaneAabb);
/// 分離軸定理（`DESIGN-phase5-hotspots.md` §10）。**未実装。導出のみ。**
inline constexpr std::size_t kAxisIntCmp = limbs_for(bits::kAxisIntCmp);
inline constexpr std::size_t kSatEdgeAxis = limbs_for(bits::kSatEdgeAxis);
inline constexpr std::size_t kSatNormalAxis = limbs_for(bits::kSatNormalAxis);

// Phase 3（SPEC-phase3.md §3.1 / §2.1）
inline constexpr std::size_t kEdgeNormal = limbs_for(bits::kEdgeNormal);
inline constexpr std::size_t kEdgeOffset = limbs_for(bits::kEdgeOffset);
inline constexpr std::size_t kAxisPointW = limbs_for(bits::kAxisPointW);
inline constexpr std::size_t kAxisPointXyz = limbs_for(bits::kAxisPointXyz);
inline constexpr std::size_t kNormalDot = limbs_for(bits::kNormalDot);

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
static_assert(64 * limbs::kPlaneAabb >= bits::kPlaneAabb, "kPlaneAabb のリム数不足");
static_assert(64 * limbs::kAxisIntCmp >= bits::kAxisIntCmp, "kAxisIntCmp のリム数不足");
static_assert(bits::kAxisIntCmp >= bits::kHomoXyz + 1, "kAxisIntCmp が x を収められない");
// Phase 3: 辺平面と代表点（SPEC-phase3 §3.1 / §2.1）
static_assert(bits::kEdgeNormal <= bits::kNormal, "辺平面の法線が支持平面の法線を超える");
static_assert(bits::kEdgeOffset <= bits::kOffset, "辺平面のオフセットが kOffset を超える");
static_assert(limbs::kEdgeNormal <= limbs::kNormal && limbs::kEdgeOffset <= limbs::kOffset,
              "辺平面は PlaneD にそのまま収まること");
// ずらしたオフセットが同じリム数に収まること（型を変えずに済む条件）
static_assert(limbs_for(bits::kOffsetShifted) == limbs::kOffset,
              "ずらした d が kOffset のリム数を超える");
static_assert(bits::kShiftLog2 + 1 <= bits::kOffset, "ずらし幅の上限が d の幅を超える");
// 軸平行 2 枚を含む交点が一般の構成点に収まること（新しい点の型を作らない条件）
static_assert(bits::kAxisPointW <= bits::kHomoW && bits::kAxisPointXyz <= bits::kHomoXyz,
              "軸平行を含む交点が HPointD に収まらない");
static_assert(64 * limbs::kNormalDot >= bits::kNormalDot, "kNormalDot のリム数不足");
// 軸平行平面のオフセットは kOffset に収まること
static_assert(bits::kAxisOffset <= bits::kOffset, "軸平行平面のオフセットが kOffset を超える");

// 実際の計算経路が §3 の上界を下回っていないことをコンパイル時に確かめる。
static_assert(64 * limbs::kSide >= bits::kSide, "kSide のリム数不足");
static_assert(64 * limbs::kCmpH >= bits::kCmpH, "kCmpH のリム数不足");
static_assert(64 * (limbs::kHomoW + limbs::kHomoXyz + 1) >= bits::kCmpH,
              "cmp_h の計算幅が理論上界を下回っている");

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_WIDTHS_HPP

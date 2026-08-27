// Krisite — 幾何述語
//
// SPEC-phase0.md §7
//
// すべて動的確保なし・例外なし・グローバル状態なし。返り値は -1 / 0 / +1。
#ifndef KRISITE_GEOM_PREDICATES_HPP
#define KRISITE_GEOM_PREDICATES_HPP

#include "krisite/arith/fixed_int.hpp"
#include "krisite/arith/ops.hpp"
#include "krisite/geom/counters.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/geom/point.hpp"
#include "krisite/geom/widths.hpp"

namespace krisite::geom {

// ---- orient3d ---------------------------------------------------------------

/// 4 点の向き。SPEC §7.2 に従い**座標差の行列式**で書く（中間値のビット幅が小さくなる）。
///
///   det | ax-dx  ay-dy  az-dz |
///       | bx-dx  by-dy  bz-dz |
///       | cx-dx  cy-dy  cz-dz |
///
/// ビット幅: 3(b+1)+2 = 3b+5 → widths.hpp bits::kOrient3d
/// orient3d の被符号値そのもの（SPEC §8.4 のビット幅実測で使う）。
inline arith::fixed_int<3 * limbs::kDiff + 1> orient3d_value(const IPoint& a, const IPoint& b,
                                                             const IPoint& c,
                                                             const IPoint& d) noexcept {
    using arith::fixed_int;
    const fixed_int<limbs::kDiff> m[3][3] = {
        {coord_diff(a.x, d.x), coord_diff(a.y, d.y), coord_diff(a.z, d.z)},
        {coord_diff(b.x, d.x), coord_diff(b.y, d.y), coord_diff(b.z, d.z)},
        {coord_diff(c.x, d.x), coord_diff(c.y, d.y), coord_diff(c.z, d.z)},
    };
    return arith::det3(m);
}

inline int orient3d(const IPoint& a, const IPoint& b, const IPoint& c, const IPoint& d) noexcept {
    return arith::sign(orient3d_value(a, b, c, d));
}

// ---- orient2d ---------------------------------------------------------------

/// 2D の向き。生の 2 座標版。
///
///   det | ax-cx  ay-cy |
///       | bx-cx  by-cy |
///
/// ビット幅: 2(b+1)+1 = 2b+3 → widths.hpp bits::kOrient2d
inline arith::fixed_int<2 * limbs::kDiff + 1> orient2d_value(std::int32_t ax, std::int32_t ay,
                                                             std::int32_t bx, std::int32_t by,
                                                             std::int32_t cx,
                                                             std::int32_t cy) noexcept {
    return arith::det2(coord_diff(ax, cx), coord_diff(ay, cy), coord_diff(bx, cx),
                       coord_diff(by, cy));
}

inline int orient2d(std::int32_t ax, std::int32_t ay, std::int32_t bx, std::int32_t by,
                    std::int32_t cx, std::int32_t cy) noexcept {
    return arith::sign(orient2d_value(ax, ay, bx, by, cx, cy));
}

/// 軸 `along` に沿って投影した 2D の向き。
///
/// 投影後の座標は右手系の巡回順にとる（X を落とす → (y,z)、Y → (z,x)、Z → (x,y)）。
/// この順序により、返る符号は三角形の法線の `along` 成分の符号に一致する。
inline int orient2d(const IPoint& a, const IPoint& b, const IPoint& c, Axis along) noexcept {
    switch (along) {
        case Axis::X:
            return orient2d(a.y, a.z, b.y, b.z, c.y, c.z);
        case Axis::Y:
            return orient2d(a.z, a.x, b.z, b.x, c.z, c.x);
        default:
            return orient2d(a.x, a.y, b.x, b.y, c.x, c.y);
    }
}

// ---- side -------------------------------------------------------------------

/// 入力点の平面に対する側。sign(N・p + d)（平面は N・x + d = 0。SPEC §3.1）。
///
/// ビット幅: N・p が 3b+5、d が 3b+5 なので和は 3b+6。
/// 余裕を見て max(kNormal+kCoord, kOffset)+1 リムで計算する。
inline constexpr std::size_t kSideIPointLimbs =
    max_limbs(limbs::kNormal + limbs::kCoord, limbs::kOffset) + 1;

inline arith::fixed_int<kSideIPointLimbs> side_value(const PlaneD& pl, const IPoint& p) noexcept {
    using namespace arith;
    constexpr std::size_t LP = kSideIPointLimbs;
    auto acc = widen<LP>(mul(pl.a, from_i64<limbs::kCoord>(p.x)));
    acc = add(acc, widen<LP>(mul(pl.b, from_i64<limbs::kCoord>(p.y))));
    acc = add(acc, widen<LP>(mul(pl.c, from_i64<limbs::kCoord>(p.z))));
    acc = add(acc, widen<LP>(pl.d));
    return acc;
}

// ---- 平面とセルの閉領域の交差（SPEC-phase2.md §2.3）--------------------------
//
// **セル $C$ を平面 $P$ で分割するのは $P \cap \overline{C} \neq \emptyset$ のときだけです。**
//
// Phase 1 は両メッシュの全平面ですべてのセルを切っていました（§4.3.1）。継ぎ目に
// T 字接合を出さないためです。絞り込んでも整合するのは、判定を**閉包**で行うからです。
//
//   $C_1$ と $C_2$ が面 $F$ を共有するとき、$P$ が $F$ 上に切断点を生むなら
//   $P$ は $F$ を横切る。$F \subset \overline{C_1} \cap \overline{C_2}$ なので
//   両方の閉包と交わり、**両セルが同じ $P$ で切る。**
//
// **開集合で判定してはいけません。** 面上に載る平面が両側から落ちます
// （SPEC-phase1 §4.2 の閉領域割り当てと同じ規律）。
//
// 「三角形が届くか」で絞ると壊れます。$\mathrm{plane}(T)$ は無限に延びるので、
// $T$ が届かないセルにも切断点を生みます。**それが Phase 1 の変異 3 です。**

/// `N・x + d` をセルの隅 1 点で評価する。ビット幅 3b+7 → widths.hpp bits::kPlaneAabb。
///
/// **`side_value(PlaneD, IPoint)` とは受ける座標の範囲が違います。** セル境界の
/// 上限 `+2^(b-1)` は `kCoordMax` を超えるので `IPoint` に入りません（§3.2）。
inline arith::fixed_int<limbs::kPlaneAabb> plane_box_value(const PlaneD& pl,
                                                           const std::int64_t p[3]) noexcept {
    using namespace arith;
    constexpr std::size_t LB = limbs::kPlaneAabb;
    static_assert(64 * LB >= bits::kPlaneAabb, "kPlaneAabb のリム数が §2.3 の上界を下回っている");
    auto acc = widen<LB>(mul(pl.a, from_i64<limbs::kAxisOffset>(p[0])));
    acc = add(acc, widen<LB>(mul(pl.b, from_i64<limbs::kAxisOffset>(p[1]))));
    acc = add(acc, widen<LB>(mul(pl.c, from_i64<limbs::kAxisOffset>(p[2]))));
    acc = add(acc, widen<LB>(pl.d));
    return acc;
}

/// 平面がセルの**閉領域** `[lo, hi]` と交わるか（SPEC-phase2 §2.3）。
///
/// 8 隅を回る必要はありません。`N` の各成分の符号に応じて `lo` / `hi` を選べば、
/// `N・x + d` の最小値と最大値が**内積 2 回**で出ます。
/// 最小値 $\le 0 \le$ 最大値なら横切ります。
///
/// **成分が 0 の軸はどちらを選んでも同じです**（その軸は寄与しない）。
inline bool plane_crosses_box(const PlaneD& pl, const std::int64_t lo[3],
                              const std::int64_t hi[3]) noexcept {
    const int sa[3] = {arith::sign(pl.a), arith::sign(pl.b), arith::sign(pl.c)};
    std::int64_t pmin[3], pmax[3];
    for (int t = 0; t < 3; ++t) {
        KRISITE_CHECK(lo[t] <= hi[t], "plane_crosses_box: lo > hi");
        pmin[t] = (sa[t] >= 0) ? lo[t] : hi[t];
        pmax[t] = (sa[t] >= 0) ? hi[t] : lo[t];
    }
    return arith::sign(plane_box_value(pl, pmin)) <= 0 &&
           arith::sign(plane_box_value(pl, pmax)) >= 0;
}

inline int side(const PlaneD& pl, const IPoint& p) noexcept {
    KRISITE_COUNT(side_calls);
    return arith::sign(side_value(pl, p));
}

/// **最重要**。構成点の平面に対する側。
///
/// SPEC §7.2:  sign(w) * sign(a*x + b*y + c*z + d*w)
/// 平面 [a,b,c,d] と点 [x,y,z,w] の**単一の 4 次元内積**であり、
/// これが同次 4 元表現（平面 N・x + d = 0）を採る理由そのものである（§3.1）。
///
/// ビット幅: N・V が 9b+19、d*w が 9b+17、総和 9b+20 → §3.1 / widths.hpp bits::kSide
/// side() の被符号値 N・V + d*w（SPEC §8.4 のビット幅実測で使う）。
inline arith::fixed_int<limbs::kSide> side_value(const PlaneD& pl, const HPointD& v) noexcept {
    using namespace arith;
    constexpr std::size_t LS = limbs::kSide;
    static_assert(64 * LS >= bits::kSide, "kSide のリム数が §3.1 の上界を下回っている");

    auto acc = resize<LS>(mul(pl.a, v.x));
    acc = add(acc, resize<LS>(mul(pl.b, v.y)));
    acc = add(acc, resize<LS>(mul(pl.c, v.z)));
    acc = add(acc, resize<LS>(mul(pl.d, v.w)));
    return acc;
}

inline int side(const PlaneD& pl, const HPointD& v) noexcept {
    KRISITE_COUNT(side_calls);
    KRISITE_CHECK(!arith::is_zero(v.w), "side: HPoint の w == 0（不変条件違反）");
    return arith::sign(v.w) * arith::sign(side_value(pl, v));
}

// ---- cmp_h ------------------------------------------------------------------

namespace detail {

/// HPoint の軸成分を取り出す。
inline const arith::fixed_int<limbs::kHomoXyz>& component(const HPointD& h, Axis ax) noexcept {
    switch (ax) {
        case Axis::X:
            return h.x;
        case Axis::Y:
            return h.y;
        default:
            return h.z;
    }
}

}  // namespace detail

/// cmp_h() の被符号値 w2*x1 - w1*x2（SPEC §8.4 のビット幅実測で使う）。
///
/// **Phase 0 で最大幅を要する述語**。ビット幅 13b+27（b=21 で 300 ビット = 5 リム、
/// b=26 で 365 ビット = 6 リム）→ §3.1 / §3.3 / widths.hpp bits::kCmpH
///
/// リム単位の積は kHomoW + kHomoXyz + 1 リムを生むが、§3 の上界に合わせて
/// kCmpH リムへ検査付きで詰める。ここで停止したら上界超過であり、§3.4 の
/// 「即座に報告」対象。
inline arith::fixed_int<limbs::kCmpH> cmp_h_value(const HPointD& h1, const HPointD& h2,
                                                  Axis ax) noexcept {
    static_assert(64 * limbs::kCmpH >= bits::kCmpH, "kCmpH のリム数が §3.1 の上界を下回っている");
    const auto& x1 = detail::component(h1, ax);
    const auto& x2 = detail::component(h2, ax);
    return arith::resize<limbs::kCmpH>(
        arith::sub_widen(arith::mul(h2.w, x1), arith::mul(h1.w, x2)));
}

/// 同次点の軸別比較。x1/w1 と x2/w2 の大小を除算なしで判定する。
///
/// SPEC §7.2:
///   sign(x1/w1 - x2/w2) = sign(w1) * sign(w2) * sign(w2*x1 - w1*x2)
///
/// 早期脱出（SPEC が「最適化として実装すること」と指定しているもの）:
///   - x1 = x2 = 0                    → 0
///   - x1/w1 と x2/w2 の符号が異なる  → 符号から即決
///   - w1 = w2                        → x1 と x2 の比較のみ（sign(w1) 倍が必要）
inline int cmp_h(const HPointD& h1, const HPointD& h2, Axis ax) noexcept {
    using namespace arith;
    const auto& x1 = detail::component(h1, ax);
    const auto& x2 = detail::component(h2, ax);

    KRISITE_CHECK(!is_zero(h1.w) && !is_zero(h2.w), "cmp_h: HPoint の w == 0（不変条件違反）");

    const int sx1 = sign(x1), sx2 = sign(x2);
    if (sx1 == 0 && sx2 == 0) return 0;

    const int sw1 = sign(h1.w), sw2 = sign(h2.w);
    const int s1 = sx1 * sw1;  // x1/w1 の符号
    const int s2 = sx2 * sw2;  // x2/w2 の符号
    if (s1 != s2) return (s1 < s2) ? -1 : 1;

    if (cmp(h1.w, h2.w) == 0) {
        // (x1 - x2) / w1 の符号
        return cmp(x1, x2) * sw1;
    }

    // 一般形。ビット幅は 13b+27 → widths.hpp bits::kCmpH（§3.1 の表には無い量）
    return sign(cmp_h_value(h1, h2, ax)) * sw1 * sw2;
}

/// 頂点テーブル用の全順序。X → Y → Z の辞書式。
inline int cmp_h_lex(const HPointD& h1, const HPointD& h2) noexcept {
    const int cx = cmp_h(h1, h2, Axis::X);
    if (cx != 0) return cx;
    const int cy = cmp_h(h1, h2, Axis::Y);
    if (cy != 0) return cy;
    return cmp_h(h1, h2, Axis::Z);
}

inline bool lex_less(const HPointD& h1, const HPointD& h2) noexcept {
    return cmp_h_lex(h1, h2) < 0;
}

inline bool h_equal(const HPointD& h1, const HPointD& h2) noexcept {
    return cmp_h_lex(h1, h2) == 0;
}

// ---- 平面の同一判定と全順序（SPEC-phase1.md §3.1）----------------------------
//
// GCD による正準化は行いません。除算・剰余・GCD は arith/ に無く、Phase 1 のために
// 算術基盤へ足すのは割に合わないためです（SPEC-phase1 §12 の非目標）。
// 代わりに比例判定と交差乗算で ID を割り当てます。

namespace detail {

/// 平面係数を共通幅で取り出す（0=a, 1=b, 2=c, 3=d）。
/// 法線は kNormal、オフセットは kOffset なので kOffset にそろえる。
inline arith::fixed_int<limbs::kOffset> plane_coeff(const PlaneD& pl, int i) noexcept {
    switch (i) {
        case 0:
            return arith::widen<limbs::kOffset>(pl.a);
        case 1:
            return arith::widen<limbs::kOffset>(pl.b);
        case 2:
            return arith::widen<limbs::kOffset>(pl.c);
        default:
            return pl.d;
    }
}

/// 最初の非零成分の添字。全零なら 4。
inline int plane_lead(const PlaneD& pl) noexcept {
    for (int i = 0; i < 4; ++i) {
        if (!arith::is_zero(plane_coeff(pl, i))) return i;
    }
    return 4;
}

}  // namespace detail

/// 2 平面の 2x2 小行列式 `p_i*q_j - p_j*q_i`。
/// ビット幅 5b+9 → widths.hpp bits::kPlaneMinor
inline arith::fixed_int<limbs::kPlaneMinor> plane_minor(const PlaneD& p, const PlaneD& q, int i,
                                                        int j) noexcept {
    const auto pi = detail::plane_coeff(p, i), pj = detail::plane_coeff(p, j);
    const auto qi = detail::plane_coeff(q, i), qj = detail::plane_coeff(q, j);
    return arith::resize<limbs::kPlaneMinor>(
        arith::sub_widen(arith::mul(pi, qj), arith::mul(pj, qi)));
}

/// 同一の幾何平面か（**向きは問いません**）。SPEC-phase1 §3.1 の比例判定。
///
/// $[a,b,c,d]$ と $[a',b',c',d']$ が比例 ⟺ 2x2 小行列式 6 本がすべて 0。
///
/// 退化平面（4 成分すべて零）は零ベクトルなので形式的にはあらゆる平面と「比例」します。
/// それでは全順序と食い違うので、**両方が退化のときだけ同一**とします。
inline bool plane_same(const PlaneD& p, const PlaneD& q) noexcept {
    const bool pz = is_null(p), qz = is_null(q);
    if (pz || qz) return pz && qz;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            if (!arith::is_zero(plane_minor(p, q, i, j))) return false;
        }
    }
    return true;
}

/// 平面の全順序。-1 / 0 / +1。`PlaneId` 付与のための整列に使います。
///
/// 符号を正規化した（最初の非零成分を正にした）うえで、
///   1. 最初の非零成分の位置で比較
///   2. 同じなら比 `p_j / p_lead` を交差乗算で比較
/// 比は符号正規化で不変なので、係数を書き換える必要はありません。
///
///   sign(p_j/p_i - q_j/q_i) = -sign(p_i*q_j - p_j*q_i) * sign(p_i) * sign(q_i)
///
/// ビット幅 6b+11 → widths.hpp bits::kPlaneOrder（実際に到達するのは 5b+8）
///
/// `plane_cmp(p, q) == 0` と `plane_same(p, q)` は同値です。
/// 別々の式から導いてあるので、テストで突き合わせています。
inline int plane_cmp(const PlaneD& p, const PlaneD& q) noexcept {
    const int ip = detail::plane_lead(p), iq = detail::plane_lead(q);
    if (ip != iq) return (ip < iq) ? -1 : 1;
    if (ip == 4) return 0;  // 両方とも退化平面

    const int sp = arith::sign(detail::plane_coeff(p, ip));
    const int sq = arith::sign(detail::plane_coeff(q, iq));
    for (int j = ip + 1; j < 4; ++j) {
        const int m = arith::sign(plane_minor(p, q, ip, j));
        if (m != 0) return -m * sp * sq;
    }
    return 0;
}

// ---- 同次点を投影した 2D 向き（SPEC-phase1.md §6.1 のレイキャスト）----------

namespace detail {

/// 同次点の成分（軸別）。
inline const arith::fixed_int<limbs::kHomoXyz>& hcomp(const HPointD& p, Axis ax) noexcept {
    return component(p, ax);
}

/// IPoint の成分（軸別）。
inline std::int32_t icomp(const IPoint& p, Axis ax) noexcept {
    return (ax == Axis::X) ? p.x : (ax == Axis::Y) ? p.y : p.z;
}

}  // namespace detail

/// 軸 `along` に沿って投影した平面での `orient2d(a, b, p)` の被符号値。
///
/// 投影後の座標は右手系の巡回順（X を落とす → (y,z)、Y → (z,x)、Z → (x,y)）。
/// 実座標での向きは `O / p.w` なので、符号は `sign(O) * sign(p.w)` になります。
///
/// ビット幅 8b+17 → widths.hpp bits::kOrient2dH
inline arith::fixed_int<limbs::kOrient2dH> orient2d_h_value(const IPoint& a, const IPoint& b,
                                                            const HPointD& p, Axis along) noexcept {
    using namespace arith;
    // 投影後の 2 軸（巡回順）
    const Axis u = (along == Axis::X) ? Axis::Y : (along == Axis::Y) ? Axis::Z : Axis::X;
    const Axis v = (along == Axis::X) ? Axis::Z : (along == Axis::Y) ? Axis::X : Axis::Y;

    // p_c - a_c * p.w
    auto rel = [&](Axis c) {
        return sub_mixed(detail::hcomp(p, c),
                         mul(from_i64<limbs::kCoord>(detail::icomp(a, c)), p.w));
    };
    const auto ru = rel(u), rv = rel(v);
    const auto du = coord_diff(detail::icomp(b, u), detail::icomp(a, u));
    const auto dv = coord_diff(detail::icomp(b, v), detail::icomp(a, v));
    return resize<limbs::kOrient2dH>(sub_widen(mul(du, rv), mul(dv, ru)));
}

/// 軸 `along` に沿って投影した 2D の向き。-1 / 0 / +1。
inline int orient2d_h(const IPoint& a, const IPoint& b, const HPointD& p, Axis along) noexcept {
    KRISITE_CHECK(!arith::is_zero(p.w), "orient2d_h: HPoint の w == 0");
    return arith::sign(orient2d_h_value(a, b, p, along)) * arith::sign(p.w);
}

// ---- 中点に対する述語（SPEC-phase1.md §6.1 の代表点フォールバック）-----------
//
// 中点 m = (X0*W1 + X1*W0 : 2*W0*W1) を**構成せず**に評価します。被符号値 F が
// 同次座標について線形であることから
//
//   F(m) = W1 * F(v0) + W0 * F(v1)
//
// が厳密に成り立ちます。m の w = 2*W0*W1 の符号は sign(W0)*sign(W1) なので、
// 述語の符号は sign(F(m)) * sign(W0) * sign(W1) です。
//
// ビット幅の導出は widths.hpp bits::kMidSide / bits::kMidOrient2dH を参照。

/// side(plane, 中点) の被符号値 W1*side_value(v0) + W0*side_value(v1)。
/// ビット幅 15b+33 → widths.hpp bits::kMidSide
inline arith::fixed_int<limbs::kMidSide> side_value(const PlaneD& pl,
                                                    const HMidPointD& m) noexcept {
    using namespace arith;
    constexpr std::size_t L = limbs::kMidSide;
    static_assert(64 * L >= bits::kMidSide, "kMidSide のリム数が上界を下回っている");
    return add(resize<L>(mul(m.v1.w, side_value(pl, m.v0))),
               resize<L>(mul(m.v0.w, side_value(pl, m.v1))));
}

/// 中点の平面に対する側。中点の w = 2*W0*W1 なので符号は sign(W0)*sign(W1)。
inline int side(const PlaneD& pl, const HMidPointD& m) noexcept {
    KRISITE_COUNT(side_calls);
    KRISITE_CHECK(!arith::is_zero(m.v0.w) && !arith::is_zero(m.v1.w), "side: 中点の端点の w == 0");
    return arith::sign(side_value(pl, m)) * arith::sign(m.v0.w) * arith::sign(m.v1.w);
}

/// orient2d_h(a, b, 中点) の被符号値。ビット幅 14b+30 → widths.hpp bits::kMidOrient2dH
inline arith::fixed_int<limbs::kMidOrient2dH> orient2d_h_value(const IPoint& a, const IPoint& b,
                                                               const HMidPointD& m,
                                                               Axis along) noexcept {
    using namespace arith;
    constexpr std::size_t L = limbs::kMidOrient2dH;
    static_assert(64 * L >= bits::kMidOrient2dH, "kMidOrient2dH のリム数が上界を下回っている");
    return add(resize<L>(mul(m.v1.w, orient2d_h_value(a, b, m.v0, along))),
               resize<L>(mul(m.v0.w, orient2d_h_value(a, b, m.v1, along))));
}

/// 中点を軸 `along` に沿って投影した 2D の向き。
inline int orient2d_h(const IPoint& a, const IPoint& b, const HMidPointD& m, Axis along) noexcept {
    KRISITE_CHECK(!arith::is_zero(m.v0.w) && !arith::is_zero(m.v1.w),
                  "orient2d_h: 中点の端点の w == 0");
    return arith::sign(orient2d_h_value(a, b, m, along)) * arith::sign(m.v0.w) *
           arith::sign(m.v1.w);
}

// ---- 3 頂点の重心に対する述語（三角形の断片用）------------------------------
//
// 中点は対角線を要求するので三角形では使えません。重心も同じ線形性で評価します。
//
//   F(c) = W1*W2*F(v0) + W0*W2*F(v1) + W0*W1*F(v2)
//   c の w = 3*W0*W1*W2 なので符号は sign(W0)*sign(W1)*sign(W2)

/// side(plane, 重心) の被符号値。ビット幅 21b+46 → widths.hpp bits::kTriSide
inline arith::fixed_int<limbs::kTriSide> side_value(const PlaneD& pl,
                                                    const HTriPointD& c) noexcept {
    using namespace arith;
    constexpr std::size_t L = limbs::kTriSide;
    static_assert(64 * L >= bits::kTriSide, "kTriSide のリム数が上界を下回っている");
    auto acc = resize<L>(mul(mul(c.v1.w, c.v2.w), side_value(pl, c.v0)));
    acc = add(acc, resize<L>(mul(mul(c.v0.w, c.v2.w), side_value(pl, c.v1))));
    acc = add(acc, resize<L>(mul(mul(c.v0.w, c.v1.w), side_value(pl, c.v2))));
    return acc;
}

/// 重心の平面に対する側。
inline int side(const PlaneD& pl, const HTriPointD& c) noexcept {
    KRISITE_COUNT(side_calls);
    KRISITE_CHECK(!arith::is_zero(c.v0.w) && !arith::is_zero(c.v1.w) && !arith::is_zero(c.v2.w),
                  "side: 重心の頂点の w == 0");
    return arith::sign(side_value(pl, c)) * arith::sign(c.v0.w) * arith::sign(c.v1.w) *
           arith::sign(c.v2.w);
}

/// orient2d_h(a, b, 重心) の被符号値。ビット幅 20b+43 → widths.hpp bits::kTriOrient2dH
inline arith::fixed_int<limbs::kTriOrient2dH> orient2d_h_value(const IPoint& a, const IPoint& b,
                                                               const HTriPointD& c,
                                                               Axis along) noexcept {
    using namespace arith;
    constexpr std::size_t L = limbs::kTriOrient2dH;
    static_assert(64 * L >= bits::kTriOrient2dH, "kTriOrient2dH のリム数が上界を下回っている");
    auto acc = resize<L>(mul(mul(c.v1.w, c.v2.w), orient2d_h_value(a, b, c.v0, along)));
    acc = add(acc, resize<L>(mul(mul(c.v0.w, c.v2.w), orient2d_h_value(a, b, c.v1, along))));
    acc = add(acc, resize<L>(mul(mul(c.v0.w, c.v1.w), orient2d_h_value(a, b, c.v2, along))));
    return acc;
}

/// 重心を軸 `along` に沿って投影した 2D の向き。
inline int orient2d_h(const IPoint& a, const IPoint& b, const HTriPointD& c, Axis along) noexcept {
    KRISITE_CHECK(!arith::is_zero(c.v0.w) && !arith::is_zero(c.v1.w) && !arith::is_zero(c.v2.w),
                  "orient2d_h: 重心の頂点の w == 0");
    return arith::sign(orient2d_h_value(a, b, c, along)) * arith::sign(c.v0.w) *
           arith::sign(c.v1.w) * arith::sign(c.v2.w);
}

// ---- 入力メッシュの向き検査（SPEC-phase1.md §3.4）----------------------------

/// 三角形 (a,b,c) と原点がなす四面体の符号付き体積 x6 = det(a, b, c)。
///
/// 閉曲面上の全三角形について総和すると、囲む立体の体積 x6 になります。
/// 外向き法線（外から見て CCW）なら正。
///
/// ビット幅: 3b+1 → widths.hpp bits::kInputVolume6（総和ぶんの余裕込み）
inline arith::fixed_int<3 * limbs::kCoord + 1> tetra_volume6(const IPoint& a, const IPoint& b,
                                                             const IPoint& c) noexcept {
    using arith::fixed_int;
    const fixed_int<limbs::kCoord> m[3][3] = {
        {arith::from_i64<limbs::kCoord>(a.x), arith::from_i64<limbs::kCoord>(a.y),
         arith::from_i64<limbs::kCoord>(a.z)},
        {arith::from_i64<limbs::kCoord>(b.x), arith::from_i64<limbs::kCoord>(b.y),
         arith::from_i64<limbs::kCoord>(b.z)},
        {arith::from_i64<limbs::kCoord>(c.x), arith::from_i64<limbs::kCoord>(c.y),
         arith::from_i64<limbs::kCoord>(c.z)},
    };
    return arith::det3(m);
}

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_PREDICATES_HPP

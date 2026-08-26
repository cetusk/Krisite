// Krisite — 幾何述語
//
// SPEC-phase0.md §7
//
// すべて動的確保なし・例外なし・グローバル状態なし。返り値は -1 / 0 / +1。
#ifndef KRISITE_GEOM_PREDICATES_HPP
#define KRISITE_GEOM_PREDICATES_HPP

#include "krisite/arith/fixed_int.hpp"
#include "krisite/arith/ops.hpp"
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

inline int side(const PlaneD& pl, const IPoint& p) noexcept {
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

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_PREDICATES_HPP

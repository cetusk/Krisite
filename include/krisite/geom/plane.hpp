// Krisite — 平面ベース表現と構成
//
// SPEC-phase0.md §6
//
// **符号の約束**（SPEC §3.1）
//   平面は同次 4 元ベクトル [a,b,c,d] として `N・x + d = 0` の形で保持する。
//   すなわち d = -N・p1。
//   この流儀により、平面と点がどちらも同次 4 元ベクトルになり、side が単一の
//   4 次元内積 sign(w) * sign(N・V + d*w) で書ける。
//
//   注意: 3 平面の交点は `N・X = -d` を解くので、Cramer の置換列は **-d** である
//   （§3.1 の表は「1 列を d ベクトルに置換」と書いているが、d = -N・p1 の定義下では
//   置換すべき右辺は -d。intersect3 のコメント参照）。
#ifndef KRISITE_GEOM_PLANE_HPP
#define KRISITE_GEOM_PLANE_HPP

#include "krisite/arith/fixed_int.hpp"
#include "krisite/arith/ops.hpp"
#include "krisite/geom/point.hpp"
#include "krisite/geom/widths.hpp"

namespace krisite::geom {

/// 平面 a*x + b*y + c*z + d = 0（同次 4 元ベクトル [a,b,c,d]）。
template <std::size_t NB, std::size_t DB>
struct Plane {
    arith::fixed_int<NB> a, b, c;  ///< 法線
    arith::fixed_int<DB> d;        ///< オフセット
};

/// Phase 0 の既定の平面型（リム数は §3 から導出）。
using PlaneD = Plane<limbs::kNormal, limbs::kOffset>;

/// 法線が零ベクトル（= 退化三角形から作られた平面）か。
template <std::size_t NB, std::size_t DB>
inline bool is_degenerate(const Plane<NB, DB>& pl) noexcept {
    return arith::is_zero(pl.a) && arith::is_zero(pl.b) && arith::is_zero(pl.c);
}

/// 座標差を fixed_int<kDiff> で返す。|a-b| < 2^b なので b+1 ビットに収まる（§3.1）。
inline arith::fixed_int<limbs::kDiff> coord_diff(std::int32_t a, std::int32_t b) noexcept {
    return arith::from_i64<limbs::kDiff>(static_cast<arith::i64>(a) - static_cast<arith::i64>(b));
}

/// 三角形 (p1,p2,p3) の支持平面。
///
/// N = (p2-p1) x (p3-p1)、d = -N・p1（SPEC §3.1。平面は N・x + d = 0）。
/// 退化三角形（面積 0）では N = 0, d = 0 となる。これは誤りではなく、
/// 呼び出し側が is_degenerate() で判定する。
inline PlaneD plane_from_triangle(const IPoint& p1, const IPoint& p2, const IPoint& p3) noexcept {
    using namespace arith;
    KRISITE_CHECK(in_range(p1) && in_range(p2) && in_range(p3),
                  "plane_from_triangle: 座標が §2 の範囲外");

    const auto ux = coord_diff(p2.x, p1.x);
    const auto uy = coord_diff(p2.y, p1.y);
    const auto uz = coord_diff(p2.z, p1.z);
    const auto vx = coord_diff(p3.x, p1.x);
    const auto vy = coord_diff(p3.y, p1.y);
    const auto vz = coord_diff(p3.z, p1.z);

    // N = u x v。det2(a,b,c,d) = a*d - b*c。ビット幅 2b+3 → §3.1
    PlaneD pl{};
    pl.a = resize<limbs::kNormal>(det2(uy, uz, vy, vz));  // uy*vz - uz*vy
    pl.b = resize<limbs::kNormal>(det2(uz, ux, vz, vx));  // uz*vx - ux*vz
    pl.c = resize<limbs::kNormal>(det2(ux, uy, vx, vy));  // ux*vy - uy*vx

    // d = -N・p1。ビット幅 3b+5 → §3.1
    const auto px = from_i64<limbs::kCoord>(p1.x);
    const auto py = from_i64<limbs::kCoord>(p1.y);
    const auto pz = from_i64<limbs::kCoord>(p1.z);
    auto acc = resize<limbs::kOffset>(mul(pl.a, px));
    acc = add(acc, resize<limbs::kOffset>(mul(pl.b, py)));
    acc = add(acc, resize<limbs::kOffset>(mul(pl.c, pz)));
    pl.d = neg(acc);  // |N・p1| < 2^(3b+4) なので符号反転はオーバーフローしない
    return pl;
}

/// 軸平行平面 `axis = coord`（SPEC-phase1.md §3.2）。
///
/// 法線は +1 の単位ベクトル、d = -coord。平面は `N・x + d = 0` なので
/// `side(plane, p) > 0` は「coord より大きい側」を意味します。
///
/// **セル境界の座標は `IPoint` ではありません。** 深度 k の境界は
/// `-2^(b-1) + m*2^(b-k)`（m = 0..2^k）で、最大値 `+2^(b-1)` は
/// `kCoordMax = 2^(b-1)-1` を超えます。だから `std::int64_t` で受けます。
/// `plane_from_triangle` では作れない平面です。
///
/// ビット幅: d が b+1 で kOffset に自明に収まる → widths.hpp bits::kAxisOffset
inline PlaneD plane_axis_aligned(Axis axis, std::int64_t coord) noexcept {
    KRISITE_CHECK(coord >= kCoordMin && coord <= -kCoordMin,
                  "plane_axis_aligned: 座標が [-2^(b-1), 2^(b-1)] の外");
    PlaneD pl{};
    const auto one = arith::from_i64<limbs::kNormal>(1);
    const auto nil = arith::zero<limbs::kNormal>();
    pl.a = (axis == Axis::X) ? one : nil;
    pl.b = (axis == Axis::Y) ? one : nil;
    pl.c = (axis == Axis::Z) ? one : nil;
    pl.d = arith::from_i64<limbs::kOffset>(-coord);
    return pl;
}

/// 4 成分すべてが零か（退化三角形から作られた平面。比例判定で特別扱いが要る）。
template <std::size_t NB, std::size_t DB>
inline bool is_null(const Plane<NB, DB>& pl) noexcept {
    return is_degenerate(pl) && arith::is_zero(pl.d);
}

/// 3 平面の交点（Cramer）。
///
/// 平面は N・x + d = 0 なので、連立は
///   [a1 b1 c1] [X]   [-d1]
///   [a2 b2 c2] [Y] = [-d2]       （X = x/w など）
///   [a3 b3 c3] [Z]   [-d3]
/// である。w = det(N1,N2,N3)、x = det(1 列目を **-d** に置換)、... とすると
/// X = x/w が成り立つ。ビット幅は w が 6b+12、x,y,z が 7b+14 → §3.1
///
/// SPEC §3.1 の表は「1 列を d ベクトルに置換」と書いているが、d = -N・p1 の定義下では
/// 右辺は -d である。d をそのまま置換すると符号が反転した点が返る。
///
/// 前提: 3 平面が一点で交わること（呼び出し側が保証）。
/// w == 0 は契約違反であり、検査ビルドでは停止する。
inline HPointD intersect3(const PlaneD& p1, const PlaneD& p2, const PlaneD& p3) noexcept {
    using namespace arith;
    // 列によって幅が違うので、行列式は最大幅にそろえて計算する（正しさ優先。SPEC §0）。
    constexpr std::size_t L = max_limbs(limbs::kNormal, limbs::kOffset);

    const fixed_int<L> a[3] = {widen<L>(p1.a), widen<L>(p2.a), widen<L>(p3.a)};
    const fixed_int<L> b[3] = {widen<L>(p1.b), widen<L>(p2.b), widen<L>(p3.b)};
    const fixed_int<L> c[3] = {widen<L>(p1.c), widen<L>(p2.c), widen<L>(p3.c)};
    // 置換列は右辺 -d
    const fixed_int<L> d[3] = {neg(widen<L>(p1.d)), neg(widen<L>(p2.d)), neg(widen<L>(p3.d))};

    HPointD h{};
    {
        const fixed_int<L> m[3][3] = {{a[0], b[0], c[0]}, {a[1], b[1], c[1]}, {a[2], b[2], c[2]}};
        h.w = resize<limbs::kHomoW>(det3(m));
    }
    {
        const fixed_int<L> m[3][3] = {{d[0], b[0], c[0]}, {d[1], b[1], c[1]}, {d[2], b[2], c[2]}};
        h.x = resize<limbs::kHomoXyz>(det3(m));
    }
    {
        const fixed_int<L> m[3][3] = {{a[0], d[0], c[0]}, {a[1], d[1], c[1]}, {a[2], d[2], c[2]}};
        h.y = resize<limbs::kHomoXyz>(det3(m));
    }
    {
        const fixed_int<L> m[3][3] = {{a[0], b[0], d[0]}, {a[1], b[1], d[1]}, {a[2], b[2], d[2]}};
        h.z = resize<limbs::kHomoXyz>(det3(m));
    }
    KRISITE_CHECK(!is_zero(h.w), "intersect3: w == 0（3 平面が一点で交わっていない）");
    return h;
}

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_PLANE_HPP

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
#include "krisite/geom/counters.hpp"
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

/// **辺平面**（`SPEC-phase3.md` §3.1）。辺 $(p_1, p_2)$ の直線を含み、支持平面と一致しない平面。
///
/// 支持平面上の点しか分類しないので、辺平面には自由度があります。**軸方向との外積を
/// 使ってください**（同 §3.1.2 の案 B）。
///
///     N_e = (p2 - p1) x e_k        e_k は軸方向の単位ベクトル
///     d_e = -N_e・p1
///
/// **相手が単位ベクトルなので乗算が起きず、幅が支持平面より小さくなります**
/// （N が b+1、d が 2b+2。`widths.hpp` bits::kEdgeNormal / kEdgeOffset）。
/// $N_s \times (p_2-p_1)$ を使うと $3b+5$ に伸び、`side` が b=26 で 6 リムに増えます。
///
/// **軸の選択は正準に決めます**（条件を満たす最小の k）。実行ごとに違う k を選ぶと
/// 出力が非決定的になります（§3.1.3）。
///
/// 条件は 2 つ。どちらも厳密な整数演算で判定できます。
///
///   1. N != 0                （Δ が e_k と平行だと退化）
///   2. N が N_s と平行でない  （支持平面と一致してはいけない）
///
/// **少なくとも 1 つの k が必ず成立します。** $(\Delta \times e_i) \times (\Delta \times e_j)
/// = (\Delta \cdot (e_i \times e_j))\,\Delta$ なので、$\Delta$ の第 k 成分が 0 のときは
/// i と j が平行になり候補方向は 2 つに減りますが、$N_s$ に平行になり得るのは
/// そのうち高々 1 つです。
///
/// `ns` には支持平面の法線を渡します（条件 2 の判定に使う）。
inline PlaneD plane_from_edge(const IPoint& p1, const IPoint& p2, const PlaneD& support,
                              Axis* chosen = nullptr) noexcept {
    const std::int64_t dx = static_cast<std::int64_t>(p2.x) - p1.x;
    const std::int64_t dy = static_cast<std::int64_t>(p2.y) - p1.y;
    const std::int64_t dz = static_cast<std::int64_t>(p2.z) - p1.z;
    // Δ x e_x = (0, dz, -dy) / Δ x e_y = (-dz, 0, dx) / Δ x e_z = (dy, -dx, 0)
    const std::int64_t cand[3][3] = {{0, dz, -dy}, {-dz, 0, dx}, {dy, -dx, 0}};

    // 条件 2 の判定は **固定幅で行います。** N_e が b+1、N_s が 2b+3 なので
    // 外積の項は 3b+4、差で 3b+5 → kOffset の幅にちょうど収まります。
    // **int64 では溢れます**（b=26 で 82 ビット）。
    using W = arith::fixed_int<limbs::kOffset>;
    const W nsx = arith::resize<limbs::kOffset>(support.a);
    const W nsy = arith::resize<limbs::kOffset>(support.b);
    const W nsz = arith::resize<limbs::kOffset>(support.c);

    int pick = -1;
#if defined(KRISITE_MUTATION_EDGE_AXIS_LAST)
    // SPEC-phase3 §10.5 の変異 14: 軸の選択を非正準にする（最小ではなく最大を選ぶ）。
    //
    // **幾何としては正しい平面が返ります。** 辺の直線を含み支持平面と一致しない、
    // という条件は満たすので、位相も体積も変わりません。壊れるのは**決定性**だけです。
    // 検出するのは `interior` の正準性の検査だけのはずです。
    for (int k = 2; k >= 0; --k) {
#else
    for (int k = 0; k < 3 && pick < 0; ++k) {
#endif
        const std::int64_t* n = cand[k];
        if (n[0] == 0 && n[1] == 0 && n[2] == 0) continue;  // 条件 1: N != 0
        const W nx = arith::from_i64<limbs::kOffset>(n[0]);
        const W ny = arith::from_i64<limbs::kOffset>(n[1]);
        const W nz = arith::from_i64<limbs::kOffset>(n[2]);
        // 条件 2: N x N_s != 0（支持平面と平行でない）
        // det2(a, b, c, d) = a*d - b*c（幅は 2L+1）
        // det2(a, b, c, d) = a*d - b*c。**引数の割り当てを間違えやすい**ので式を併記する。
        const auto cx = arith::det2(ny, nz, nsy, nsz);  // ny*nsz - nz*nsy
        const auto cy = arith::det2(nz, nx, nsz, nsx);  // nz*nsx - nx*nsz
        const auto cz = arith::det2(nx, ny, nsx, nsy);  // nx*nsy - ny*nsx
        if (arith::is_zero(cx) && arith::is_zero(cy) && arith::is_zero(cz)) continue;
        pick = k;
#if defined(KRISITE_MUTATION_EDGE_AXIS_LAST)
        break;
#endif
    }
    KRISITE_CHECK(pick >= 0, "plane_from_edge: 3 通りすべてが退化または支持平面と平行");
    if (chosen != nullptr) *chosen = static_cast<Axis>(pick);

    PlaneD pl{};
    pl.a = arith::from_i64<limbs::kNormal>(cand[pick][0]);
    pl.b = arith::from_i64<limbs::kNormal>(cand[pick][1]);
    pl.c = arith::from_i64<limbs::kNormal>(cand[pick][2]);
    // d = -N・p1。各項は (b+1) + b = 2b+1 なので int64 に収まる（b <= 31）
    auto acc = arith::from_i64<limbs::kOffset>(0);
    acc = arith::add(acc, arith::from_i64<limbs::kOffset>(cand[pick][0] * std::int64_t{p1.x}));
    acc = arith::add(acc, arith::from_i64<limbs::kOffset>(cand[pick][1] * std::int64_t{p1.y}));
    acc = arith::add(acc, arith::from_i64<limbs::kOffset>(cand[pick][2] * std::int64_t{p1.z}));
    pl.d = arith::neg(acc);
    return pl;
}

/// 平面係数を**桁上げ**する（`SPEC-phase3.md` §2.1.2 の刻みを細かくする操作）。
///
/// 平面 $(N, d)$ と $(2^s N, 2^s d)$ は同じ平面です。係数を上げてから $d$ を 1 動かすと、
/// **幾何的な移動量が $2^{-s}$ になります。** `fine = false` なら何もしません。
///
/// 上げ幅は **`kNormal` / `kOffset` の余裕の範囲**に限ります（型を変えないため）。
inline PlaneD plane_scaled_for_shift(const PlaneD& pl, bool fine) noexcept {
    if (!fine) return pl;
    const std::size_t na = arith::min_bits(pl.a), nb = arith::min_bits(pl.b),
                      nc = arith::min_bits(pl.c), nd = arith::min_bits(pl.d);
    const std::size_t nmax = na > nb ? (na > nc ? na : nc) : (nb > nc ? nb : nc);
    if (nmax == 0 || nd == 0) return pl;
    // 余裕の小さいほうに合わせる。ずらし幅（1）ぶんを 1 ビット残す
    const std::size_t room_n = (bits::kNormal > nmax) ? bits::kNormal - nmax : 0;
    const std::size_t room_d = (bits::kOffset > nd + 1) ? bits::kOffset - nd - 1 : 0;
    const std::size_t s = room_n < room_d ? room_n : room_d;
    if (s == 0) return pl;
    PlaneD out{};
    out.a = arith::shl_bits(pl.a, s);
    out.b = arith::shl_bits(pl.b, s);
    out.c = arith::shl_bits(pl.c, s);
    out.d = arith::shl_bits(pl.d, s);
    return out;
}

/// 平面を法線方向に**平行移動**する（`SPEC-phase3.md` §2.1.2 の予備経路）。
///
/// `N・x + d = 0` に対し `d' = d + delta`。**delta > 0 は「side が負の側」を狭めます。**
/// 内側が負の辺平面なら、正の delta で平面が内側へ動きます。
///
/// ずらし幅は $2^{kShiftLog2}$ までに制限します（`widths.hpp` bits::kOffsetShifted）。
inline PlaneD plane_shifted(const PlaneD& pl, std::int64_t delta) noexcept {
    KRISITE_CHECK(delta != 0, "plane_shifted: delta == 0");
    KRISITE_CHECK(delta > -(std::int64_t{1} << bits::kShiftLog2) &&
                      delta < (std::int64_t{1} << bits::kShiftLog2),
                  "plane_shifted: ずらし幅が上限を超える");
    PlaneD out = pl;
    out.d = arith::add(pl.d, arith::from_i64<limbs::kOffset>(delta));
    return out;
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
    KRISITE_COUNT(intersect3_calls);
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

// Krisite — 点の型
//
// SPEC-phase0.md §6
#ifndef KRISITE_GEOM_POINT_HPP
#define KRISITE_GEOM_POINT_HPP

#include <cstdint>

#include "krisite/arith/fixed_int.hpp"
#include "krisite/arith/ops.hpp"
#include "krisite/geom/widths.hpp"

namespace krisite::geom {

/// 入力点（格子上の整数座標）。SPEC §2 の coord。
struct IPoint {
    std::int32_t x, y, z;
};

constexpr bool operator==(const IPoint& a, const IPoint& b) noexcept {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
constexpr bool operator!=(const IPoint& a, const IPoint& b) noexcept {
    return !(a == b);
}

/// 座標が §2 の許容範囲に収まっているか。
constexpr bool in_range(const IPoint& p) noexcept {
    const long long lo = kCoordMin, hi = kCoordMax;
    return p.x >= lo && p.x <= hi && p.y >= lo && p.y <= hi && p.z >= lo && p.z <= hi;
}

/// 同次座標の構成点。実座標は (x/w, y/w, z/w)。w != 0 が不変条件。
template <std::size_t VB, std::size_t WB>
struct HPoint {
    arith::fixed_int<VB> x, y, z;
    arith::fixed_int<WB> w;
};

/// Phase 0 の既定の構成点型（リム数は §3 から導出）。
using HPointD = HPoint<limbs::kHomoXyz, limbs::kHomoW>;

/// 入力点を同次点として見る（w = 1）。
inline HPointD to_homogeneous(const IPoint& p) noexcept {
    HPointD h{};
    h.x = arith::from_i64<limbs::kHomoXyz>(p.x);
    h.y = arith::from_i64<limbs::kHomoXyz>(p.y);
    h.z = arith::from_i64<limbs::kHomoXyz>(p.z);
    h.w = arith::from_i64<limbs::kHomoW>(1);
    return h;
}

/// 軸の指定。
enum class Axis : int { X = 0, Y = 1, Z = 2 };

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_POINT_HPP

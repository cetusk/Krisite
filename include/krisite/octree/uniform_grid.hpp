// Krisite — 固定深度の一様分割
//
// SPEC-phase1.md §2.2（分割深度）, §3.2（セル面も平面である）, §4.2（セル割り当て）, §8
//
// **セル境界は必ず 2 の冪の格子座標に置きます。** §3.2 でセル面を厳密な平面として
// 扱う前提であり、かつ Phase 5 の Morton 線形八分木と同じ構造にするためです。
// ここを外すと八分木を 2 つ持つことになります（SPEC-phase1 §8）。
//
// **深度は実行時パラメータです**（§2.2）。本フェーズの中心的な診断は
// 「深度 0 の出力と深度 2 の出力を同一プロセス内で比較する」ことなので、
// コンパイル時定数にしてはいけません。
#ifndef KRISITE_OCTREE_UNIFORM_GRID_HPP
#define KRISITE_OCTREE_UNIFORM_GRID_HPP

#include <array>
#include <cstdint>

#include "krisite/config.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/geom/point.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace krisite::octree {

/// セルの添字。深度 k のとき各軸 0 .. 2^k - 1。
struct CellIndex {
    std::uint32_t i = 0, j = 0, k = 0;
};

/// 固定深度の一様分割。深度 d で 2^d × 2^d × 2^d セル。
///
/// 深度 0 はセル 1 個（分割なし）で、§2.2 の対照実験に使います。
class UniformGrid {
public:
    explicit UniformGrid(unsigned depth) noexcept : depth_(depth) {
        KRISITE_CHECK(depth + 1 <= kCoordBits, "深度が b-1 を超えている");
        KRISITE_CHECK(depth <= 10, "深度が大きすぎる（セル数が 8^d）");
    }

    unsigned depth() const noexcept { return depth_; }
    std::uint32_t per_axis() const noexcept { return std::uint32_t{1} << depth_; }
    std::size_t cell_count() const noexcept { return std::size_t{1} << (3 * depth_); }

    /// セルの一辺の長さ 2^(b-d)。
    std::int64_t cell_size() const noexcept { return std::int64_t{1} << (kCoordBits - depth_); }

    /// 軸方向の境界座標 `-2^(b-1) + m * 2^(b-d)`（m = 0 .. 2^d）。
    ///
    /// **`IPoint` ではありません。** m = 2^d のとき `+2^(b-1)` になり、
    /// `kCoordMax = 2^(b-1)-1` を超えます（§3.2）。
    std::int64_t bound(std::uint32_t m) const noexcept {
        KRISITE_CHECK(m <= per_axis(), "bound: 境界の添字が範囲外");
        return kCoordMin + static_cast<std::int64_t>(m) * cell_size();
    }

    /// セルの閉領域 [lo, hi]（各軸）。
    std::int64_t lo(std::uint32_t idx) const noexcept { return bound(idx); }
    std::int64_t hi(std::uint32_t idx) const noexcept { return bound(idx + 1); }

    /// セルの 6 面の平面。法線は +1 の単位ベクトルで、d = -境界座標。
    ///
    /// 隣り合うセルは共有面について**同一の平面**を得ます（同じ境界座標から作るため）。
    /// これが §5.1 の「両セルで同じ平面3つ組から生まれる」前提を支えます。
    std::array<geom::PlaneD, 6> cell_planes(const CellIndex& c) const noexcept {
        using geom::Axis;
        using geom::plane_axis_aligned;
        return {
            plane_axis_aligned(Axis::X, lo(c.i)), plane_axis_aligned(Axis::X, hi(c.i)),
            plane_axis_aligned(Axis::Y, lo(c.j)), plane_axis_aligned(Axis::Y, hi(c.j)),
            plane_axis_aligned(Axis::Z, lo(c.k)), plane_axis_aligned(Axis::Z, hi(c.k)),
        };
    }

    /// 深度 d で内部に現れる境界平面（m = 1 .. 2^d - 1）の総数。深度 0 では 0。
    std::size_t interior_plane_count() const noexcept {
        return 3 * static_cast<std::size_t>(per_axis() - 1);
    }

private:
    unsigned depth_;
};

/// 三角形の軸別の境界（AABB）。
struct Aabb {
    std::int64_t lo[3];
    std::int64_t hi[3];
};

inline Aabb triangle_aabb(const geom::IPoint& a, const geom::IPoint& b,
                          const geom::IPoint& c) noexcept {
    const std::int64_t xs[3] = {a.x, b.x, c.x};
    const std::int64_t ys[3] = {a.y, b.y, c.y};
    const std::int64_t zs[3] = {a.z, b.z, c.z};
    Aabb r{};
    for (int t = 0; t < 3; ++t) {
        const std::int64_t* v = (t == 0) ? xs : (t == 1) ? ys : zs;
        r.lo[t] = v[0] < v[1] ? (v[0] < v[2] ? v[0] : v[2]) : (v[1] < v[2] ? v[1] : v[2]);
        r.hi[t] = v[0] > v[1] ? (v[0] > v[2] ? v[0] : v[2]) : (v[1] > v[2] ? v[1] : v[2]);
    }
    return r;
}

/// 三角形をセルに割り当てるか。
///
/// **判定は半開区間 `[lo, hi)` で行います**（SPEC-phase2 §3.4。`adaptive.hpp` に導出）。
/// 開領域にすると、境界に載る三角形が**両側から落ちます**。§10.5 の変異 2 がそれを突きます。
/// 半開はそれとは別で、**境界平面に完全に乗る面を上側のセルだけに入れる**規則です。
///
/// **保守的な判定です。** 三角形の AABB とセルの閉領域が重なるかだけを見るので、
/// 実際には交わらない三角形も割り当てられます。ただし**取りこぼしはしません**
/// （AABB は三角形を含むため）。過剰割り当てはクリップ結果が空になるだけで無害です。
/// 厳密な三角形-箱交差判定（分離軸定理）は Phase 2 の最適化候補です。
inline bool assign_to_cell(const Aabb& tri, const UniformGrid& g, const CellIndex& c) noexcept {
    const std::int64_t clo[3] = {g.lo(c.i), g.lo(c.j), g.lo(c.k)};
    const std::int64_t chi[3] = {g.hi(c.i), g.hi(c.j), g.hi(c.k)};
    for (int t = 0; t < 3; ++t) {
        if (tri.hi[t] < clo[t]) return false;
        if (tri.lo[t] >= chi[t]) return false;  // 半開区間 [lo, hi)（SPEC-phase2 §3.4）
    }
    return true;
}

/// 開領域で判定する版。**本来使ってはいけません。**
/// §10.5 の変異 2 が、これに差し替えるとテストが落ちることを確かめます。
inline bool assign_to_cell_open(const Aabb& tri, const UniformGrid& g,
                                const CellIndex& c) noexcept {
    const std::int64_t clo[3] = {g.lo(c.i), g.lo(c.j), g.lo(c.k)};
    const std::int64_t chi[3] = {g.hi(c.i), g.hi(c.j), g.hi(c.k)};
    for (int t = 0; t < 3; ++t) {
        if (tri.hi[t] <= clo[t]) return false;
        if (tri.lo[t] >= chi[t]) return false;
    }
    return true;
}

}  // namespace krisite::octree

#endif  // KRISITE_OCTREE_UNIFORM_GRID_HPP

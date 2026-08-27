// Krisite — 適応分割
//
// SPEC-phase2.md §3.1（分割の判定）, §3.2（early-out）
//
// **セル $C$ を分割するのは、分割によって仕事が減るときだけです**（§3.1）。
//
//   $C$ が $A$ と $B$ の両方の三角形を含む   交差の可能性がある。分割の候補
//   $C$ が片方しか含まない                   分割不要（§3.2 の early-out）
//   $C$ が含む三角形数が閾値以下             分割の利得より木の管理コストが上回る
//   最大深度に達した                         打ち切り
//
// **閾値と最大深度は実行時パラメータです**（§3.1）。そして
// **固定深度は「常に最大深度まで分割する」という特別な場合**として残します。
// これが §0.1 の正解器（Phase 1 の挙動）であり、**消してはいけない資産**です。
//
// `UniformGrid`（`uniform_grid.hpp`）は深度を 1 つ固定した木で、ここでの
// `SubdivisionPolicy{max_depth, uniform = true}` と同じものを表します。
// **セル境界は同じ 2 の冪の格子座標**なので、両者の葉は完全に一致します。
#ifndef KRISITE_OCTREE_ADAPTIVE_HPP
#define KRISITE_OCTREE_ADAPTIVE_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

#include "krisite/config.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/octree/uniform_grid.hpp"

namespace krisite::octree {

/// 適応分割の葉。**深度ごとに添字の意味が変わる**ので、深度も持ちます。
struct Cell {
    unsigned depth = 0;
    std::uint32_t i = 0, j = 0, k = 0;
};

inline bool operator<(const Cell& a, const Cell& b) noexcept {
    if (a.depth != b.depth) return a.depth < b.depth;
    if (a.i != b.i) return a.i < b.i;
    if (a.j != b.j) return a.j < b.j;
    return a.k < b.k;
}

/// セルの閉領域 [lo, hi]。
struct CellBox {
    std::int64_t lo[3], hi[3];
};

/// 深度 `depth` の軸方向の境界座標 `-2^(b-1) + m * 2^(b-depth)`。
///
/// **`IPoint` ではありません。** m = 2^depth のとき `+2^(b-1)` になり `kCoordMax` を
/// 超えます（`SPEC-phase1.md` §3.2）。
inline std::int64_t cell_bound(unsigned depth, std::uint32_t m) noexcept {
    KRISITE_CHECK(depth + 1 <= kCoordBits, "cell_bound: 深度が b-1 を超えている");
    return kCoordMin + (static_cast<std::int64_t>(m) << (kCoordBits - depth));
}

inline CellBox box_of(const Cell& c) noexcept {
    const std::uint32_t idx[3] = {c.i, c.j, c.k};
    CellBox r{};
    for (int t = 0; t < 3; ++t) {
        r.lo[t] = cell_bound(c.depth, idx[t]);
        r.hi[t] = cell_bound(c.depth, idx[t] + 1);
    }
    return r;
}

/// セルの 6 面の平面。法線は +1 の単位ベクトルで、d = -境界座標。
///
/// **隣り合うセルは共有面について同一の平面を得ます**（同じ境界座標から作るため）。
/// 深さが違っても、粗い側の面は細かい側の面の一部と同一平面です。
inline std::array<geom::PlaneD, 6> cell_planes(const Cell& c) noexcept {
    using geom::Axis;
    using geom::plane_axis_aligned;
    const CellBox b = box_of(c);
    return {
        plane_axis_aligned(Axis::X, b.lo[0]), plane_axis_aligned(Axis::X, b.hi[0]),
        plane_axis_aligned(Axis::Y, b.lo[1]), plane_axis_aligned(Axis::Y, b.hi[1]),
        plane_axis_aligned(Axis::Z, b.lo[2]), plane_axis_aligned(Axis::Z, b.hi[2]),
    };
}

/// **半開区間 `[lo, hi)` での割り当て**（SPEC-phase2 §3.4）。
///
/// > **各セルは各軸で $[\mathrm{lo}, \mathrm{hi})$ を占める。領域は正確に分割される。**
///
/// `SPEC-phase1.md` §4.2 の閉領域割り当てを置き換えたものです。閉領域だと
/// **セル境界平面に完全に乗る面が両側のセルに入り**、適応分割では粗い側の 1 枚と
/// 細かい側の 4 枚が入れ子になって併合されず、**同じ面積が二重に出力されます**（§3.4.1）。
///
/// **§4.2 の目的は維持されます**（§3.4.2）。挙動が変わるのは境界平面に**完全に乗る**面
/// だけで、接するだけの三角形は引き続き両側に入ります。開領域のように
/// **両側から落ちることはありません。**
///
/// | 三角形 | 閉領域 | 半開 |
/// |---|---|---|
/// | $x = m$ をまたぐ | 両セル | 両セル |
/// | $x = m$ に**接する** | 両セル | **両セル** |
/// | $x = m$ に**完全に乗る** | 両セル | **上側のみ** |
///
/// 幾何を失う経路もありません。入力座標は $2^{b-1}-1$ までなので、最上位のセル境界
/// $+2^{b-1}$ に面が乗ることはありません。
inline bool assign_to_cell(const Aabb& tri, const CellBox& c) noexcept {
    for (int t = 0; t < 3; ++t) {
        if (tri.hi[t] < c.lo[t]) return false;
#if defined(KRISITE_MUTATION_CLOSED_CELLS)
        // SPEC-phase2 §9.3 の変異 12: §3.4 の割り当てを閉領域に戻す。
        // **適応分割 × 面がセル境界に乗るケースでしか検出できません**（C が破裂する）。
        if (tri.lo[t] > c.hi[t]) return false;
#else
        if (tri.lo[t] >= c.hi[t]) return false;
#endif
    }
    return true;
}

/// 開領域で判定する版。**本来使ってはいけません**（§10.5 の変異 2）。
inline bool assign_to_cell_open(const Aabb& tri, const CellBox& c) noexcept {
    for (int t = 0; t < 3; ++t) {
        if (tri.hi[t] <= c.lo[t]) return false;
        if (tri.lo[t] >= c.hi[t]) return false;
    }
    return true;
}

/// §3.1 の分割方針。**すべて実行時パラメータです。**
struct SubdivisionPolicy {
    unsigned max_depth = 0;
    /// **固定深度モード。** 真なら常に最大深度まで分割する（Phase 1 の挙動 = 正解器）。
    bool uniform = true;
    /// セルが含む三角形数がこれ以下なら分割しない。0 なら閾値では打ち切らない。
    std::size_t leaf_threshold = 0;
};

/// §3.1 の判定で葉を列挙する。
///
/// `count(cell, &na, &nb)` はそのセルの**閉領域**に割り当たる A / B の三角形数を返す
/// こと。**閉領域であることが要点です**（`SPEC-phase1.md` §4.2）。
///
/// 返す葉は `Cell` の全順序で整列します。**出力の再現性のためです。**
/// 固定深度モードでは `UniformGrid` の三重ループ（i, j, k）と同じ順序になります。
template <class CountFn>
inline std::vector<Cell> build_leaves(const SubdivisionPolicy& p, CountFn count) {
    std::vector<Cell> leaves;
    std::vector<Cell> stack{Cell{0, 0, 0, 0}};
    while (!stack.empty()) {
        const Cell c = stack.back();
        stack.pop_back();

        bool split = false;
        if (c.depth < p.max_depth) {
            if (p.uniform) {
                split = true;  // 固定深度: 常に最大深度まで
            } else {
                std::size_t na = 0, nb = 0;
                count(c, &na, &nb);
                // **両方を含むときだけ分割します**（§3.1）。片方しか無いセルを割っても
                // 交差は生まれず、断片が増えるだけです
                split = (na > 0 && nb > 0) && (na + nb > p.leaf_threshold);
            }
        }
        if (!split) {
            leaves.push_back(c);
            continue;
        }
        for (int t = 0; t < 8; ++t) {
            stack.push_back(Cell{c.depth + 1, c.i * 2 + static_cast<std::uint32_t>(t & 1),
                                 c.j * 2 + static_cast<std::uint32_t>((t >> 1) & 1),
                                 c.k * 2 + static_cast<std::uint32_t>((t >> 2) & 1)});
        }
    }
    std::sort(leaves.begin(), leaves.end());
    return leaves;
}

/// 葉の添字を**最大深度の解像度**に正規化する（`SPEC-phase1.md` §5.4 の局所性）。
///
/// 深さが混ざると添字がそのままでは比べられません。最大深度の格子に写して比べます。
/// **固定深度モードでは恒等写像**なので、Phase 1 の数値がそのまま再現されます。
inline void normalized_index(const Cell& c, unsigned max_depth, std::uint32_t out[3]) noexcept {
    KRISITE_CHECK(c.depth <= max_depth, "normalized_index: 葉の深度が最大深度を超えている");
    const unsigned s = max_depth - c.depth;
    out[0] = c.i << s;
    out[1] = c.j << s;
    out[2] = c.k << s;
}

}  // namespace krisite::octree

#endif  // KRISITE_OCTREE_ADAPTIVE_HPP

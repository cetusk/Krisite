// Krisite — 三角メッシュ
//
// SPEC-phase1.md §2.1（入力）, §3.4（向き付けの規約）
//
// `arith/` の制約（動的確保禁止など）はここには及びません。`std::vector` を使います。
#ifndef KRISITE_MESH_TRI_MESH_HPP
#define KRISITE_MESH_TRI_MESH_HPP

#include <array>
#include <cstdint>
#include <vector>

#include "krisite/arith/ops.hpp"
#include "krisite/geom/point.hpp"
#include "krisite/geom/predicates.hpp"
#include "krisite/geom/widths.hpp"

namespace krisite::mesh {

using VertexId = std::uint32_t;

/// 三角形。頂点添字を反時計回り（外から見て）に持つ。
using Tri = std::array<VertexId, 3>;

/// 入力メッシュ。頂点は整数格子上、三角形は頂点添字。
///
/// **向きの規約**（SPEC-phase1 §3.4）: 外から見て CCW = 外向き法線。
/// これは取り決めではなく、`side(plane, p) > 0` が外側になることから導かれます。
struct TriMesh {
    std::vector<geom::IPoint> vertices;
    std::vector<Tri> triangles;

    bool empty() const noexcept { return triangles.empty(); }
};

/// 符号付き体積 x6 = Σ det(v_i, v_j, v_k)。
///
/// 閉じた向き付き曲面なら、外向き法線のとき正になります。
/// 全頂点が `IPoint` なので厳密な整数で、GMP は要りません。
/// ビット幅: 3b+1+16 → widths.hpp bits::kInputVolume6
inline arith::fixed_int<geom::limbs::kInputVolume6> signed_volume6(const TriMesh& m) noexcept {
    using namespace krisite::geom;
    auto acc = arith::zero<limbs::kInputVolume6>();
    KRISITE_CHECK(m.triangles.size() <= (std::size_t{1} << bits::kMaxTrianglesLog2),
                  "signed_volume6: 三角形数が kMaxTrianglesLog2 の想定を超えている");
    for (const Tri& t : m.triangles) {
        acc = arith::add(acc, arith::resize<limbs::kInputVolume6>(tetra_volume6(
                                  m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]])));
    }
    return acc;
}

/// 向きの規約を満たすか（符号付き体積が正）。SPEC-phase1 §3.4 の入力検査。
///
/// 空メッシュは体積 0 なので false を返します。呼び出し側で区別してください。
inline bool is_outward_oriented(const TriMesh& m) noexcept {
    return arith::sign(signed_volume6(m)) > 0;
}

/// 全頂点が §2 の座標範囲に収まっているか。
inline bool coords_in_range(const TriMesh& m) noexcept {
    for (const geom::IPoint& p : m.vertices) {
        if (!geom::in_range(p)) return false;
    }
    return true;
}

}  // namespace krisite::mesh

#endif  // KRISITE_MESH_TRI_MESH_HPP

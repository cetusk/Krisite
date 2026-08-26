// Krisite — 点が立体の内部にあるかの判定（レイキャスト）
//
// SPEC-phase1.md §6.1
//
// +X 方向のレイを撃ち、入力三角形との交差回数の偶奇で判定します。
//
// **前提: 判定点は立体の境界上にありません。** 呼び出し側が保証します
// （断片の符号ベクトルがすべて非零の頂点を選ぶ）。これにより前方交差の判定
// `N・p + d != 0` が保証され、**記号的摂動は投影の内外判定だけで済みます。**
//
// 記号的摂動: レイ原点を実座標で `(0, ε, ε²)` だけずらしたものとして扱います。
// 投影の向きが 0 になったとき、ε の項 `-(b.z - a.z)`、それも 0 なら ε² の項
// `(b.y - a.y)` で符号を決めます。
#ifndef KRISITE_CSG_RAYCAST_HPP
#define KRISITE_CSG_RAYCAST_HPP

#include "krisite/geom/plane.hpp"
#include "krisite/geom/predicates.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace krisite::csg {

namespace detail {

/// 摂動後の投影向き。0 を返しません。
inline int perturbed_orient(int raw, const geom::IPoint& a, const geom::IPoint& b) noexcept {
    if (raw != 0) return raw;
    // ε の項: -(b.z - a.z)
    if (b.z != a.z) return (a.z > b.z) ? 1 : -1;
    // ε² の項: (b.y - a.y)
    if (b.y != a.y) return (b.y > a.y) ? 1 : -1;
    return 0;  // 投影が退化（三角形の法線 x 成分が 0）。呼び出し側が弾く
}

inline int proj_orient(const geom::IPoint& a, const geom::IPoint& b,
                       const geom::IPoint& p) noexcept {
    return perturbed_orient(geom::orient2d(a.y, a.z, b.y, b.z, p.y, p.z), a, b);
}

inline int proj_orient(const geom::IPoint& a, const geom::IPoint& b,
                       const geom::HPointD& p) noexcept {
    return perturbed_orient(geom::orient2d_h(a, b, p, geom::Axis::X), a, b);
}

/// レイ（+X）が三角形を前方で横切るか。
template <class Point>
bool crosses(const geom::IPoint& a, const geom::IPoint& b, const geom::IPoint& c, const Point& p) {
    const geom::PlaneD pl = geom::plane_from_triangle(a, b, c);
    if (geom::is_degenerate(pl)) return false;
    // 法線の x 成分が 0 ⟺ 投影三角形が退化 ⟺ レイが平面に平行
    const int nx = arith::sign(pl.a);
    if (nx == 0) return false;

    const int o1 = proj_orient(a, b, p);
    const int o2 = proj_orient(b, c, p);
    const int o3 = proj_orient(c, a, p);
    if (o1 == 0 || o2 == 0 || o3 == 0) return false;
    if (o1 != o2 || o2 != o3) return false;  // 投影三角形の外

    // 前方か: t = -(N・p + d)/N.x > 0 ⟺ sign(N・p + d) != sign(N.x)
    const int sp = geom::side(pl, p);
    KRISITE_CHECK(sp != 0, "crosses: 判定点が三角形の平面上にある（呼び出し側の契約違反）");
    return sp != nx;
}

}  // namespace detail

/// 点が閉じた立体 `m` の内部にあるか。
///
/// **`p` が `m` の境界上に無いことを呼び出し側が保証すること。**
template <class Point>
inline bool point_inside(const mesh::TriMesh& m, const Point& p) {
    int hits = 0;
    for (const mesh::Tri& t : m.triangles) {
        if (detail::crosses(m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]], p)) ++hits;
    }
    return (hits % 2) == 1;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_RAYCAST_HPP

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

/// 軸の成分を取り出す。
inline std::int32_t comp(const geom::IPoint& p, geom::Axis a) noexcept {
    return (a == geom::Axis::X) ? p.x : (a == geom::Axis::Y) ? p.y : p.z;
}

/// 投影後の 2 軸（`orient2d_h` と同じ巡回順）。
inline geom::Axis proj_u(geom::Axis along) noexcept {
    return (along == geom::Axis::X)   ? geom::Axis::Y
           : (along == geom::Axis::Y) ? geom::Axis::Z
                                      : geom::Axis::X;
}
inline geom::Axis proj_v(geom::Axis along) noexcept {
    return (along == geom::Axis::X)   ? geom::Axis::Z
           : (along == geom::Axis::Y) ? geom::Axis::X
                                      : geom::Axis::Y;
}

/// 平面の法線の、軸 `along` 方向の成分。
inline int normal_comp_sign(const geom::PlaneD& pl, geom::Axis along) noexcept {
    return (along == geom::Axis::X)   ? arith::sign(pl.a)
           : (along == geom::Axis::Y) ? arith::sign(pl.b)
                                      : arith::sign(pl.c);
}

/// 摂動後の投影向き。0 を返しません。
///
/// レイ原点を実座標で $(\varepsilon_u, \varepsilon_v^2)$（投影後の 2 軸）だけずらした
/// ものとして扱います。投影の向きが 0 になったとき、$\varepsilon$ の項 $-(b_v - a_v)$、
/// それも 0 なら $\varepsilon^2$ の項 $(b_u - a_u)$ で符号を決めます。
///
/// **軸によらず同じ形です**（`orient2d_h` の巡回順に合わせてあります）。
inline int perturbed_orient(int raw, const geom::IPoint& a, const geom::IPoint& b,
                            geom::Axis along) noexcept {
    if (raw != 0) return raw;
    const geom::Axis u = proj_u(along), v = proj_v(along);
    if (comp(b, v) != comp(a, v)) return (comp(a, v) > comp(b, v)) ? 1 : -1;
    if (comp(b, u) != comp(a, u)) return (comp(b, u) > comp(a, u)) ? 1 : -1;
    return 0;  // 投影が退化（法線の along 成分が 0）。呼び出し側が弾く
}

inline int proj_orient(const geom::IPoint& a, const geom::IPoint& b, const geom::IPoint& p,
                       geom::Axis along = geom::Axis::X) noexcept {
    const geom::Axis u = proj_u(along), v = proj_v(along);
    return perturbed_orient(
        geom::orient2d(comp(a, u), comp(a, v), comp(b, u), comp(b, v), comp(p, u), comp(p, v)), a,
        b, along);
}

inline int proj_orient(const geom::IPoint& a, const geom::IPoint& b, const geom::HPointD& p,
                       geom::Axis along = geom::Axis::X) noexcept {
    return perturbed_orient(geom::orient2d_h(a, b, p, along), a, b, along);
}

/// レイ（`along` の正方向）が三角形を前方で横切るか。**平面を渡す版。**
///
/// **支持平面は呼び出し側で一度だけ作ってください**（`SPEC-phase5.md` の CP1.5）。
/// レイキャストは source の全三角形を走査するので、ここで作り直すと
/// **1 レイあたり $O(n)$ 回の平面構成**になります（実測 14.7 億回）。
template <class Point>
bool crosses_with(const geom::PlaneD& pl, const geom::IPoint& a, const geom::IPoint& b,
                  const geom::IPoint& c, const Point& p, geom::Axis along = geom::Axis::X) {
    if (geom::is_degenerate(pl)) return false;
    // 法線の along 成分が 0 ⟺ 投影三角形が退化 ⟺ レイが平面に平行
    const int nx = normal_comp_sign(pl, along);
    if (nx == 0) return false;

    const int o1 = proj_orient(a, b, p, along);
    const int o2 = proj_orient(b, c, p, along);
    const int o3 = proj_orient(c, a, p, along);
    if (o1 == 0 || o2 == 0 || o3 == 0) return false;
    if (o1 != o2 || o2 != o3) return false;  // 投影三角形の外

    // 前方か: t = -(N・p + d)/N_along > 0 ⟺ sign(N・p + d) != sign(N_along)
    const int sp = geom::side(pl, p);
    KRISITE_CHECK(sp != 0, "crosses: 判定点が三角形の平面上にある（呼び出し側の契約違反）");
    return sp != nx;
}

/// 平面を渡さない版（互換）。**ホットパスでは `crosses_with` を使ってください。**
template <class Point>
bool crosses(const geom::IPoint& a, const geom::IPoint& b, const geom::IPoint& c, const Point& p,
             geom::Axis along = geom::Axis::X) {
    return crosses_with(geom::plane_from_triangle(a, b, c), a, b, c, p, along);
}

}  // namespace detail

/// 点が `m` の**境界上**にあるか（三角形の内部・辺・頂点のいずれか）。
///
/// `point_inside` の前提（判定点が境界上に無い）を呼び出し側が確かめるための述語です。
///
/// **「相手の平面上に無い」で代用してはいけません。** 平面上にあることと境界上にある
/// ことは別です。断片の頂点は隣接する切断平面の上に必ず載るので、平面で判定すると
/// 実際には境界から外れた点まで弾いてしまい、代表点が見つからなくなります
/// （CP2 のケース 2 で実際に起きました。IMPL-phase1 §6.6）。
///
/// 射影軸は**法線成分が非零の軸**を選びます。そうすれば射影した三角形が潰れません。
template <class Point>
inline bool point_on_boundary(const mesh::TriMesh& m, const Point& p) {
    for (const mesh::Tri& t : m.triangles) {
        const geom::IPoint& a = m.vertices[t[0]];
        const geom::IPoint& b = m.vertices[t[1]];
        const geom::IPoint& c = m.vertices[t[2]];
        const geom::PlaneD pl = geom::plane_from_triangle(a, b, c);
        if (geom::is_degenerate(pl)) continue;
        if (geom::side(pl, p) != 0) continue;  // 三角形の平面上にすら無い
        const geom::Axis ax = (arith::sign(pl.a) != 0)   ? geom::Axis::X
                              : (arith::sign(pl.b) != 0) ? geom::Axis::Y
                                                         : geom::Axis::Z;
        const int o1 = geom::orient2d_h(a, b, p, ax);
        const int o2 = geom::orient2d_h(b, c, p, ax);
        const int o3 = geom::orient2d_h(c, a, p, ax);
        const bool neg = (o1 < 0) || (o2 < 0) || (o3 < 0);
        const bool pos = (o1 > 0) || (o2 > 0) || (o3 > 0);
        if (!(neg && pos)) return true;  // 内部・辺・頂点のいずれか
    }
    return false;
}

/// **巻き数**（`SPEC-phase3.md` §5.1）。判定点が `m` の表面に載っている場合も扱います。
///
/// 表面に載っている点では巻き数が両側で違うので、**3 つに分けて返します。**
///
///   `w_other`  判定点を含まない面だけから決まる巻き数（表裏で共通）
///   `c_front`  基準平面の法線 $+N$ 側で、載っている面が寄与する分
///   `c_back`   $-N$ 側で寄与する分
///
/// 表側の巻き数は `w_other + c_front`、裏側は `w_other + c_back` です。
///
/// **面に載っている点をレイキャストしてはいけません**（`crosses` の契約違反）。
/// 載っている面を先に除くのが、この関数の役割です。
///
/// 巻き数の寄与は `sign(N_t \cdot x)`（レイは +X 方向）。閉じた向き付き立体なら
/// 内部で 1、外部で 0 になります。**自己交差や入れ子では 2 以上になります**（それが目的）。
/// `planes` を渡すと、**支持平面の作り直しをやめます**（`SPEC-phase5.md` の CP1.5）。
///
/// `planes[j]` は `plane_from_triangle(m.triangles[j] の 3 頂点)` と**同一**でなければ
/// なりません。**`PlaneTable` で intern した平面は使えません** —
/// `intern` は `orientation_differs` を返すので、**表の平面は符号が逆のことがあります。**
///
/// 渡さなければ従来どおり毎回作ります（**出力は同じ**）。
template <class Point>
inline void winding_split(const mesh::TriMesh& m, const Point& p, const geom::PlaneD& ref,
                          int* w_other, int* c_front, int* c_back,
                          const geom::PlaneD* planes = nullptr) {
    *w_other = 0;
    *c_front = 0;
    *c_back = 0;
    int sheet_strict = 0, n_strict = 0, sheet_edge = 0;
    bool any_edge = false;

    // **レイの軸は、基準平面の法線が非零な軸から選びます。**
    // レイがシートと平行だと、どちら側が跨ぐかを決められません（実際に踏みました）。
    const geom::Axis along = (arith::sign(ref.a) != 0)   ? geom::Axis::X
                             : (arith::sign(ref.b) != 0) ? geom::Axis::Y
                                                         : geom::Axis::Z;
    KRISITE_CHECK(detail::normal_comp_sign(ref, along) != 0, "winding_split: 基準平面が退化");

    for (std::size_t j = 0; j < m.triangles.size(); ++j) {
        const mesh::Tri& t = m.triangles[j];
        const geom::IPoint& a = m.vertices[t[0]];
        const geom::IPoint& b = m.vertices[t[1]];
        const geom::IPoint& c = m.vertices[t[2]];
        const geom::PlaneD pl =
            (planes != nullptr) ? planes[j] : geom::plane_from_triangle(a, b, c);
        if (geom::is_degenerate(pl)) continue;

        bool on_face = false, strict = false;
        if (geom::side(pl, p) == 0) {
            const geom::Axis ax = (arith::sign(pl.a) != 0)   ? geom::Axis::X
                                  : (arith::sign(pl.b) != 0) ? geom::Axis::Y
                                                             : geom::Axis::Z;
            const int o1 = geom::orient2d_h(a, b, p, ax);
            const int o2 = geom::orient2d_h(b, c, p, ax);
            const int o3 = geom::orient2d_h(c, a, p, ax);
            const bool neg = (o1 < 0) || (o2 < 0) || (o3 < 0);
            const bool pos = (o1 > 0) || (o2 > 0) || (o3 > 0);
            on_face = !(neg && pos);
            strict = (o1 != 0) && (o2 != 0) && (o3 != 0);
        }
        if (!on_face) {
            // **平面を使い回します。** ここで作り直すと 1 三角形あたり 2 回になります
            if (detail::crosses_with(pl, a, b, c, p, along)) {
                *w_other += detail::normal_comp_sign(pl, along);
            }
            continue;
        }
        // **載っている面（シート）。** レイはここから出るので、跨ぐかどうかは
        // レイの向きと基準法線の関係で決まります（下の分岐）。
        //
        // **辺や頂点に載っていると、同じシートを複数の三角形が共有します。**
        // 内部に載っているものは 1 枚ずつ、境界だけのものは 1 枚として数えます。
        const int contrib = detail::normal_comp_sign(pl, along);
        if (strict) {
            sheet_strict += contrib;
            ++n_strict;
        } else {
            any_edge = true;
            sheet_edge = contrib;
        }
    }

    const int sheet = sheet_strict + ((n_strict == 0 && any_edge) ? sheet_edge : 0);

    // **どちら側がシートを跨ぐか。**
    //
    // 点 $p$ をシートから $\pm N$ に $\varepsilon$ だけ離すと、+X のレイが
    // シートを跨ぐのは「レイが $-N$ 側へ向かうとき」、すなわち $X \cdot N < 0$ のとき。
    //
    //   $X \cdot N < 0$ → 表側（$+N$）のレイがシートを跨ぐ → 表に寄与
    //   $X \cdot N > 0$ → 裏側（$-N$）のレイが跨ぐ → 裏に寄与
    //   $X \cdot N = 0$ → レイがシートと平行 → どちらも跨がない
    const int xn = detail::normal_comp_sign(ref, along);
    if (xn < 0) {
        *c_front = sheet;
    } else if (xn > 0) {
        *c_back = sheet;
    }
}

/// 点が閉じた立体 `m` の内部にあるか。
///
/// **`p` が `m` の境界上に無いことを呼び出し側が保証すること**（`point_on_boundary`）。
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

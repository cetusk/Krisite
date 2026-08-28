// Krisite — 出口（`SPEC-phase3.md` §6）
//
//     to_mesh : PolySoup → TriMesh
//
// **Phase 1 / 2 で作った機構がここに移動します。新規実装ではありません。**
//
//   1. 頂点の併合（平面3つ組 + 値ベースの第2段）   SPEC-phase1 §5.3
//   2. T 頂点の解決（案 D）                        SPEC-phase2 §2.4
//   3. 接触の分裂（辺・頂点）                      SPEC-phase2 §5
//   4. 三角形化                                    SPEC-phase1 §2.4.4
//
// **共平面の再併合（§6.4）と snap rounding（§6.5）は CP2 の範囲外です。**
//
// ---
//
// **縫合はここで独立に実装しています。** `boolean.hpp` の二項実装（§10.1 の正解器）が
// 同じ規則を持っていますが、**同じコードを共有すると正解器になりません**
// （`CLAUDE.md`「正解器は被検体と別経路で書く」）。規則は仕様（§5.3）が正です。
#ifndef KRISITE_CSG_TO_MESH_HPP
#define KRISITE_CSG_TO_MESH_HPP

#include <algorithm>
#include <array>
#include <map>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/tjunction.hpp"
#include "krisite/geom/predicates.hpp"
#include "krisite/mesh/split.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace krisite::csg {

/// 出力メッシュ（構成点 + 三角形）。`boolean.hpp` の `BoolMesh` と同じ形です。
struct SoupMesh {
    std::vector<geom::HPointD> vertices;
    std::vector<mesh::Tri> triangles;
    bool empty() const noexcept { return triangles.empty(); }
};

/// §11 の記録。
struct ToMeshStats {
    std::size_t constructed_points = 0;  ///< 第1段（平面3つ組）で作った点
    std::size_t merged_points = 0;       ///< 第2段（値）の併合後
    std::size_t merged_by_value = 0;     ///< 第2段が併合した数
    TJunctionStats t{};
    mesh::SplitStats split{};
};

struct ToMeshOptions {
    bool split_contacts = true;  ///< §6.3。既定 ON、フラグで無効化可
    bool resolve_t = true;       ///< §6.2 の T 頂点解決
};

/// スープを三角メッシュにする（`SPEC-phase3.md` §6）。
inline SoupMesh to_mesh(const PolySoup& s, const ToMeshOptions& opt = {},
                        ToMeshStats* stats = nullptr) {
    ToMeshStats st;
    SoupMesh out;
    if (s.polys.empty()) {
        if (stats != nullptr) *stats = st;
        return out;
    }

    // ---- 1. 縫合（§5.3）------------------------------------------------------
    //
    // 第1段: 平面3つ組をキーに引く。第2段: 値が厳密に等しい点を併合する。
    // **4 平面以上が 1 点で交わると 3つ組が違っても同じ点になる**ので、第2段が要ります。
    std::map<std::array<PlaneId, 3>, std::uint32_t> by_key;
    std::vector<geom::HPointD> points;
    std::vector<std::vector<std::uint32_t>> raw(s.polys.size());

    for (std::size_t pi = 0; pi < s.polys.size(); ++pi) {
        const Fragment& f = s.polys[pi].frag;
        const std::size_t n = vertex_count(f);
        raw[pi].reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            std::array<PlaneId, 3> k{f.support, f.edge[(i + n - 1) % n], f.edge[i]};
            std::sort(k.begin(), k.end());
            auto it = by_key.find(k);
            if (it == by_key.end()) {
                const auto id = static_cast<std::uint32_t>(points.size());
                points.push_back(fragment_vertex(s.table, f, i));
                it = by_key.emplace(k, id).first;
            }
            raw[pi].push_back(it->second);
        }
    }
    st.constructed_points = points.size();

    std::vector<std::uint32_t> order(points.size());
    for (std::uint32_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
        return geom::lex_less(points[a], points[b]);
    });
    std::vector<std::uint32_t> remap(points.size());
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        const auto id = static_cast<std::uint32_t>(out.vertices.size());
        out.vertices.push_back(points[order[i]]);
        while (j < order.size() && geom::h_equal(points[order[i]], points[order[j]])) {
            remap[order[j]] = id;
            ++j;
        }
        if (j - i > 1) st.merged_by_value += (j - i - 1);
        i = j;
    }
    st.merged_points = out.vertices.size();

    // ---- 2. T 頂点の解決（§6.2）+ 3. 三角形化 --------------------------------
    PlaneVertexIndex index;
    {
        std::vector<PlaneId> sup;
        sup.reserve(s.polys.size());
        for (const Poly& q : s.polys) sup.push_back(q.frag.support);
        std::sort(sup.begin(), sup.end());
        sup.erase(std::unique(sup.begin(), sup.end()), sup.end());
        index.build(s.table, out.vertices, sup);
    }

    std::vector<int> tri_src;
    for (std::size_t pi = 0; pi < s.polys.size(); ++pi) {
        const Fragment& f = s.polys[pi].frag;
        std::vector<std::uint32_t> poly;
        poly.reserve(raw[pi].size());
        for (std::uint32_t v : raw[pi]) poly.push_back(remap[v]);
        if (poly.size() < 3) continue;
        std::vector<PlaneId> edge = f.edge;
        const std::size_t before = out.triangles.size();
        if (opt.resolve_t) {
            const TPolygon tp =
                insert_t_vertices(s.table, out.vertices, index, f.support, edge, poly, &st.t);
            fan_triangulate(tp, out.triangles, &st.t);
        } else {
            TPolygon tp;
            tp.corners = static_cast<std::uint32_t>(poly.size());
            tp.vertex = poly;
            tp.is_corner.assign(poly.size(), 1);
            for (std::uint32_t i = 0; i < poly.size(); ++i) tp.orig.push_back(i);
            fan_triangulate(tp, out.triangles, &st.t);
        }
        tri_src.insert(tri_src.end(), out.triangles.size() - before,
                       static_cast<int>(s.polys[pi].src));
    }

    // ---- 4. 接触の分裂（§6.3）------------------------------------------------
    if (opt.split_contacts && !out.triangles.empty()) {
        std::vector<std::uint32_t> origin;
        out.triangles =
            mesh::split_contacts(out.triangles, out.vertices.size(), &origin, &st.split);
        for (std::uint32_t o : origin) out.vertices.push_back(out.vertices[o]);
    }

    if (stats != nullptr) *stats = st;
    return out;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_TO_MESH_HPP

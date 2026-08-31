// Krisite — 入力メッシュを「平面ごとの多角形」に併合する
//
// SPEC-phase1.md §3.1「立方体の 1 面は 2 三角形ですが、平面は 1 枚です」
//
// **なぜ併合が要るのか。**
//
// 立方体の 1 面を 2 三角形のまま扱うと、その対角辺が問題になります。対角辺は
// 三角形分割の産物であって幾何的な特徴ではなく、**2 平面の交線として表せません**
// （両側の三角形が同一平面なので交線が定まらない）。したがって対角辺上に生まれる
// 交点は「3 平面の交点」の枠に収まらず、§3.3 の頂点表現から外れます。
//
// 併合すれば、多角形の各辺は必ず**別の平面を持つ隣接面**と共有されるので、
// 辺は 2 平面の交線、頂点は 3 平面の交点になります。クリップで生まれる新しい頂点も
// （支持平面, 辺の平面, 切断平面）の 3 つ組で表せて、表現が閉じます。
//
// **Phase 1 の制約**: 併合後の多角形が凸であることを要求します（§4.3 が
// 「凸なので扇でよい」と書いているのはこの前提）。コーパス（立方体・四面体・
// 量子化した回転立方体）はすべて凸面なので成立します。
#ifndef KRISITE_CSG_FACES_HPP
#define KRISITE_CSG_FACES_HPP

#include <algorithm>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "krisite/csg/plane_table.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace krisite::csg {

/// 平面ごとに併合された面（凸多角形）。
///
/// 辺 `edge[i]` は頂点 i と頂点 i+1 を結び、**その辺を定める隣接面の平面**です。
/// 頂点 i = `support` ∩ `edge[(i-1+n)%n]` ∩ `edge[i]`。
struct Face {
    PlaneId support = kNoPlane;
    /// 面の外向き法線が `support` の法線と**逆向き**か。
    bool flipped = false;
    std::vector<PlaneId> edge;
    /// 検証用に、併合前の頂点添字（反時計回り）も持つ。
    std::vector<mesh::VertexId> loop;
    int owner = 0;  ///< 0 = A, 1 = B
};

namespace detail {

inline std::pair<mesh::VertexId, mesh::VertexId> undirected(mesh::VertexId a,
                                                            mesh::VertexId b) noexcept {
    return (a < b) ? std::pair<mesh::VertexId, mesh::VertexId>{a, b}
                   : std::pair<mesh::VertexId, mesh::VertexId>{b, a};
}

}  // namespace detail

/// メッシュを平面ごとの多角形に併合し、`table` に平面を登録する。
///
/// 前提: 閉じた向き付き多様体メッシュ（`mesh::check_topology` で確認しておくこと）。
/// 併合後の各面は単一の閉ループをなす必要があります（穴あきの平面パッチは非対応）。
inline std::vector<Face> build_faces(const mesh::TriMesh& m, int owner, PlaneTable& table) {
    // ---- 三角形ごとの平面 ----
    const std::size_t nt = m.triangles.size();
    std::vector<PlaneRef> tri_plane(nt);
    for (std::size_t i = 0; i < nt; ++i) {
        const auto& t = m.triangles[i];
        const geom::PlaneD pl =
            geom::plane_from_triangle(m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]]);
        KRISITE_CHECK(!geom::is_degenerate(pl), "build_faces: 退化三角形が入力に含まれる");
        tri_plane[i] = table.intern(pl);
    }

    // ---- 有向辺 → 三角形 ----
    std::map<std::pair<mesh::VertexId, mesh::VertexId>, std::size_t> edge_owner;
    for (std::size_t i = 0; i < nt; ++i) {
        const auto& t = m.triangles[i];
        for (int k = 0; k < 3; ++k) {
            edge_owner[{t[k], t[(k + 1) % 3]}] = i;
        }
    }
    // 隣接三角形（無向辺で結ぶ）
    auto neighbor_of = [&](std::size_t tri, int k) -> std::size_t {
        const auto& t = m.triangles[tri];
        auto it = edge_owner.find({t[(k + 1) % 3], t[k]});  // 逆向きの辺を持つ三角形
        KRISITE_CHECK(it != edge_owner.end(), "build_faces: 閉じていないメッシュ");
        return it->second;
    };

    // ---- 同一平面の隣接三角形をグループ化 ----
    std::vector<std::size_t> group(nt, static_cast<std::size_t>(-1));
    std::vector<std::vector<std::size_t>> groups;
    for (std::size_t s = 0; s < nt; ++s) {
        if (group[s] != static_cast<std::size_t>(-1)) continue;
        const std::size_t gid = groups.size();
        groups.emplace_back();
        std::vector<std::size_t> stack{s};
        group[s] = gid;
        while (!stack.empty()) {
            const std::size_t cur = stack.back();
            stack.pop_back();
            groups[gid].push_back(cur);
            for (int k = 0; k < 3; ++k) {
                const std::size_t nb = neighbor_of(cur, k);
                if (group[nb] == static_cast<std::size_t>(-1) &&
                    tri_plane[nb].id == tri_plane[cur].id) {
                    group[nb] = gid;
                    stack.push_back(nb);
                }
            }
        }
    }

    // ---- グループごとに境界ループを取り出す ----
    std::vector<Face> faces;
    faces.reserve(groups.size());
    for (const auto& g : groups) {
        Face f;
        f.owner = owner;
        f.support = tri_plane[g.front()].id;
        f.flipped = tri_plane[g.front()].flipped;

        // 境界辺 = 隣接三角形が別グループにあるもの
        std::map<mesh::VertexId, mesh::VertexId> next;  // 有向辺 a→b
        std::map<mesh::VertexId, PlaneId> edge_plane;   // 辺 a→b を定める隣の平面
        for (std::size_t tri : g) {
            const auto& t = m.triangles[tri];
            for (int k = 0; k < 3; ++k) {
                const std::size_t nb = neighbor_of(tri, k);
                if (group[nb] == group[tri]) continue;  // 内部辺（対角線など）
                const mesh::VertexId a = t[k], b = t[(k + 1) % 3];
                KRISITE_CHECK(next.find(a) == next.end(),
                              "build_faces: 境界ループが単一でない（穴あき平面パッチ？）");
                next[a] = b;
                edge_plane[a] = tri_plane[nb].id;
            }
        }
        KRISITE_CHECK(!next.empty(), "build_faces: 境界辺が無い（メッシュが 1 平面のみ？）");

        // ループを辿る
        const mesh::VertexId start = next.begin()->first;
        mesh::VertexId cur = start;
        do {
            f.loop.push_back(cur);
            f.edge.push_back(edge_plane[cur]);
            cur = next[cur];
        } while (cur != start && f.loop.size() <= next.size());
        KRISITE_CHECK(cur == start && f.loop.size() == next.size(),
                      "build_faces: 境界ループが閉じない");

        // ---- 共線の頂点を落とす（CP1.6）----
        {
            const geom::PlaneD& sp = table.at(f.support);
            const geom::Axis ax = (arith::sign(sp.a) != 0)   ? geom::Axis::X
                                  : (arith::sign(sp.b) != 0) ? geom::Axis::Y
                                                             : geom::Axis::Z;
            bool changed = true;
            while (changed && f.loop.size() > 3) {
                changed = false;
                const std::size_t n = f.loop.size();
                for (std::size_t i = 0; i < n; ++i) {
                    const std::size_t pv = (i + n - 1) % n, nx = (i + 1) % n;
                    if (geom::orient2d(m.vertices[f.loop[pv]], m.vertices[f.loop[i]],
                                       m.vertices[f.loop[nx]], ax) != 0) {
                        continue;
                    }
                    f.edge[pv] = std::min(f.edge[pv], f.edge[i]);
                    f.loop.erase(f.loop.begin() + static_cast<std::ptrdiff_t>(i));
                    f.edge.erase(f.edge.begin() + static_cast<std::ptrdiff_t>(i));
                    changed = true;
                    break;
                }
            }
        }
        faces.push_back(std::move(f));
    }
    return faces;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_FACES_HPP

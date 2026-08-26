// Krisite — ブール演算のパイプライン
//
// SPEC-phase1.md §4（パイプライン）, §4.3（局所 arrangement）, §5（縫合）, §6（分類）
//
//   1. 平面抽出・ID 付与（§3.1）
//   2. 固定深度で八分木を構築。各セルに三角形を割り当てる（§4.2）
//   3. セルごとに局所 arrangement を計算（§4.3）
//   4. 頂点の同一性を解決して大域メッシュに縫合（§5）
//   5. 領域に分割し、内外を分類（§6）
//   6. 演算に応じて選択・向き付けして出力
//
// **手順 4 と 5 の順序に注意。** 縫合を先に済ませてから分類します。
#ifndef KRISITE_CSG_BOOLEAN_HPP
#define KRISITE_CSG_BOOLEAN_HPP

#include <algorithm>
#include <cstdint>
#include <map>
#include <set>
#include <vector>

#include "krisite/csg/faces.hpp"
#include "krisite/csg/fragment.hpp"
#include "krisite/csg/plane_table.hpp"
#include "krisite/csg/raycast.hpp"
#include "krisite/mesh/topology.hpp"
#include "krisite/mesh/tri_mesh.hpp"
#include "krisite/octree/uniform_grid.hpp"

namespace krisite::csg {

enum class BoolOp { Union, Intersection, Difference };

/// 出力メッシュ。頂点は構成点（有理数）なので `mesh::TriMesh` にはできません。
/// 位相検査は三角形の添字だけで行えるので `mesh::check_topology(triangles)` が使えます。
struct BoolMesh {
    std::vector<geom::HPointD> vertices;
    std::vector<mesh::Tri> triangles;
    bool empty() const noexcept { return triangles.empty(); }
};

/// SPEC-phase1 §4.3.3 と §5.4 が要求する計数。
struct BoolStats {
    std::size_t fragments = 0;            ///< 分割後の断片数（§4.3.3）
    std::size_t duplicate_fragments = 0;  ///< 重複割り当てが生んだ重複断片（§5.4）
    std::size_t constructed_points = 0;   ///< 構成点の総数（§5.4 の分母）
    std::size_t merged_by_value = 0;      ///< 第1段が取りこぼし第2段が併合した数（§5.4）
    std::size_t max_planes_at_point = 0;  ///< 1 点に集まる平面の最大枚数（§5.4）
    std::size_t regions = 0;              ///< 相異なる符号ベクトルの数（§6.1）
    std::size_t raycasts = 0;             ///< レイキャスト回数
    std::size_t side_calls = 0;           ///< 参考: side の呼び出し数
    std::size_t intersect3_calls = 0;     ///< 参考: intersect3 の呼び出し数
};

namespace detail {

/// 断片の正準キー（重複検出用）。最小の平面 ID から始まるよう回転する。
inline std::vector<PlaneId> canonical_edges(const std::vector<PlaneId>& e) {
    const std::size_t n = e.size();
    std::size_t best = 0;
    for (std::size_t i = 1; i < n; ++i) {
        if (e[i] < e[best]) best = i;
    }
    std::vector<PlaneId> out;
    out.reserve(n);
    for (std::size_t j = 0; j < n; ++j) out.push_back(e[(best + j) % n]);
    return out;
}

/// 頂点の平面 3 つ組（昇順に正規化）。§5.3 の第1段のキー。
inline std::array<PlaneId, 3> vertex_key(const Fragment& f, std::size_t i) {
    const std::size_t n = f.edge.size();
    std::array<PlaneId, 3> k{f.support, f.edge[(i + n - 1) % n], f.edge[i]};
    std::sort(k.begin(), k.end());
    return k;
}

}  // namespace detail

/// ブール演算。`depth` は八分木の深度（実行時パラメータ。SPEC §2.2）。
inline BoolMesh boolean_op(const mesh::TriMesh& A, const mesh::TriMesh& B, BoolOp op,
                           unsigned depth, BoolStats* stats = nullptr) {
    BoolStats st;

    // ---- 1. 平面抽出・ID 付与 ----
    PlaneTable table;
    const std::vector<Face> faces_a = build_faces(A, 0, table);
    const std::vector<Face> faces_b = build_faces(B, 1, table);
    const std::size_t n_mesh_planes = table.size();
    std::vector<PlaneId> mesh_planes(n_mesh_planes);
    for (std::size_t i = 0; i < n_mesh_planes; ++i) mesh_planes[i] = static_cast<PlaneId>(i);

    // 各面の AABB（セル割り当て用）
    auto face_aabb = [&](const mesh::TriMesh& m, const Face& f) {
        octree::Aabb r{};
        for (int t = 0; t < 3; ++t) {
            r.lo[t] = krisite::kCoordMax;
            r.hi[t] = krisite::kCoordMin;
        }
        for (mesh::VertexId vid : f.loop) {
            const geom::IPoint& p = m.vertices[vid];
            const std::int64_t c[3] = {p.x, p.y, p.z};
            for (int t = 0; t < 3; ++t) {
                r.lo[t] = std::min(r.lo[t], c[t]);
                r.hi[t] = std::max(r.hi[t], c[t]);
            }
        }
        return r;
    };

    // ---- 2 + 3. セルごとに局所 arrangement ----
    const octree::UniformGrid grid(depth);
    std::vector<Fragment> frags;
    std::set<std::pair<PlaneId, std::vector<PlaneId>>> seen;

    const std::uint32_t n = grid.per_axis();
    for (std::uint32_t ci = 0; ci < n; ++ci) {
        for (std::uint32_t cj = 0; cj < n; ++cj) {
            for (std::uint32_t ck = 0; ck < n; ++ck) {
                const octree::CellIndex cell{ci, cj, ck};
                // セル面の平面 ID（保持側つき）: lo 面は +、hi 面は -
                struct CellPlane {
                    PlaneId id;
                    int keep;
                };
                std::vector<CellPlane> cps;
                if (depth > 0) {
                    const auto ps = grid.cell_planes(cell);
                    for (int k = 0; k < 6; ++k) {
                        const PlaneRef r = table.intern(ps[k]);
                        // 平面の代表が裏返っているなら保持側も反転する
                        const int base = (k % 2 == 0) ? +1 : -1;
                        cps.push_back({r.id, r.flipped ? -base : base});
                    }
                }

                for (int which = 0; which < 2; ++which) {
                    const mesh::TriMesh& m = (which == 0) ? A : B;
                    const std::vector<Face>& fs = (which == 0) ? faces_a : faces_b;
                    for (const Face& f : fs) {
                        if (depth > 0 && !octree::assign_to_cell(face_aabb(m, f), grid, cell)) {
                            continue;
                        }
                        Fragment frag = face_to_fragment(f);
                        // セルの 6 面でクリップ
                        bool alive = true;
                        for (const CellPlane& cp : cps) {
                            if (cp.id == frag.support) continue;
                            if (!clip_fragment(table, frag, cp.id, cp.keep)) {
                                alive = false;
                                break;
                            }
                        }
                        if (!alive) continue;

                        // 【両メッシュの全平面】で分割（§4.3.1）
                        std::vector<Fragment> pieces{frag};
                        for (PlaneId q : mesh_planes) {
                            std::vector<Fragment> next;
                            next.reserve(pieces.size());
                            for (const Fragment& p : pieces) {
                                if (q == p.support) {
                                    next.push_back(p);
                                    continue;
                                }
                                const SplitResult r = split_fragment(table, p, q);
                                if (r.has_pos) next.push_back(r.pos);
                                if (r.has_neg) next.push_back(r.neg);
                            }
                            pieces.swap(next);
                        }
                        for (Fragment& p : pieces) {
                            auto key = std::make_pair(p.support, detail::canonical_edges(p.edge));
                            if (!seen.insert(key).second) {
                                ++st.duplicate_fragments;
                                continue;
                            }
                            frags.push_back(std::move(p));
                        }
                    }
                }
            }
        }
    }
    st.fragments = frags.size();

    // ---- 5. 符号ベクトルによる分類（§4.3.2, §6.1）----
    //
    // 相手メッシュの平面に対する符号ベクトルが等しい断片は、同じ凸領域に属するので
    // 内外が一致します。領域ごとに 1 回だけレイキャストします。
    auto other_planes = [&](int owner) {
        std::vector<PlaneId> r;
        for (const Face& f : (owner == 0) ? faces_b : faces_a) {
            if (std::find(r.begin(), r.end(), f.support) == r.end()) r.push_back(f.support);
        }
        std::sort(r.begin(), r.end());
        return r;
    };
    const std::vector<PlaneId> planes_of[2] = {other_planes(0), other_planes(1)};

    std::map<std::pair<int, std::vector<std::int8_t>>, bool> region_inside;
    std::vector<bool> keep(frags.size(), false);

    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        const Fragment& f = frags[fi];
        const std::vector<PlaneId>& qs = planes_of[f.owner];
        std::vector<std::int8_t> sig(qs.size());
        for (std::size_t k = 0; k < qs.size(); ++k) {
            sig[k] = static_cast<std::int8_t>(fragment_sign(table, f, qs[k]));
        }
        const auto key = std::make_pair(f.owner, sig);
        auto it = region_inside.find(key);
        if (it == region_inside.end()) {
            // この領域の代表点を探す: 相手の全平面に対して符号が非零の頂点
            bool decided = false, inside = false;
            for (std::size_t vi = 0; vi < vertex_count(f) && !decided; ++vi) {
                const geom::HPointD v = fragment_vertex(table, f, vi);
                bool ok = true;
                for (PlaneId q : qs) {
                    if (geom::side(table.at(q), v) == 0) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) continue;
                inside = point_inside((f.owner == 0) ? B : A, v);
                decided = true;
                ++st.raycasts;
            }
            KRISITE_CHECK(decided,
                          "boolean_op: 領域の代表点が見つからない（全頂点が相手の平面上）");
            it = region_inside.emplace(key, inside).first;
            ++st.regions;
        }
        const bool inside = it->second;
        switch (op) {
            case BoolOp::Union:
                keep[fi] = !inside;
                break;
            case BoolOp::Intersection:
                keep[fi] = inside;
                break;
            case BoolOp::Difference:
                keep[fi] = (f.owner == 0) ? !inside : inside;
                break;
        }
    }

    // ---- 4. 縫合（§5）----
    //
    // 第1段: 平面3つ組をキーにハッシュ表で引く
    // 第2段: 全構成点を lex_less で整列し、値が厳密に等しいものを併合して再写像する
    std::map<std::array<PlaneId, 3>, std::uint32_t> by_key;
    std::vector<geom::HPointD> points;
    std::vector<std::vector<std::size_t>> point_frag_vertices;

    auto vertex_id = [&](const Fragment& f, std::size_t i) {
        const auto k = detail::vertex_key(f, i);
        auto it = by_key.find(k);
        if (it != by_key.end()) return it->second;
        const auto id = static_cast<std::uint32_t>(points.size());
        points.push_back(fragment_vertex(table, f, i));
        by_key.emplace(k, id);
        return id;
    };

    std::vector<std::vector<std::uint32_t>> polys;
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        if (!keep[fi]) continue;
        const Fragment& f = frags[fi];
        std::vector<std::uint32_t> poly;
        for (std::size_t i = 0; i < vertex_count(f); ++i) poly.push_back(vertex_id(f, i));
        polys.push_back(std::move(poly));
    }
    st.constructed_points = points.size();

    // 第2段: 値ベースの併合
    std::vector<std::uint32_t> order(points.size());
    for (std::uint32_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
        return geom::lex_less(points[a], points[b]);
    });
    std::vector<std::uint32_t> remap(points.size());
    std::vector<geom::HPointD> merged;
#if !defined(KRISITE_MUTATION_NO_STAGE2)
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        const auto id = static_cast<std::uint32_t>(merged.size());
        merged.push_back(points[order[i]]);
        while (j < order.size() && geom::h_equal(points[order[i]], points[order[j]])) {
            remap[order[j]] = id;
            ++j;
        }
        if (j - i > 1) st.merged_by_value += (j - i - 1);
        i = j;
    }
#else
    // §10.5 の変異 1: 第2段を無効化する
    merged = points;
    for (std::uint32_t i = 0; i < remap.size(); ++i) remap[i] = i;
#endif

    // ---- 6. 出力（扇状三角形化）----
    BoolMesh out;
    out.vertices = std::move(merged);
    for (const auto& poly : polys) {
        if (poly.size() < 3) continue;
        for (std::size_t i = 1; i + 1 < poly.size(); ++i) {
            out.triangles.push_back({remap[poly[0]], remap[poly[i]], remap[poly[i + 1]]});
        }
    }

    // ---- §5.4: 1 点に集まる平面の最大枚数（総当たり）----
    for (const geom::HPointD& v : out.vertices) {
        std::size_t cnt = 0;
        for (PlaneId q = 0; q < static_cast<PlaneId>(n_mesh_planes); ++q) {
            if (geom::side(table.at(q), v) == 0) ++cnt;
        }
        st.max_planes_at_point = std::max(st.max_planes_at_point, cnt);
    }

    if (stats) *stats = st;
    return out;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_BOOLEAN_HPP

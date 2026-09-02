// Krisite — 自己交差の検査（`SPEC-phase3.md` §5.6 の NSI）
//
// **NSI は呼び出し側が宣言するフラグ**ですが、**実データで宣言が偽になりました**
// （`IMPL-phase5.md` §25.2）。この検査は「宣言してよいか」を機械的に確かめます。
//
// ## 何を「自己交差」と呼ぶか
//
// 最適化の安全性が要求するのは、正確には
//
// > **三角形 $T_1$ の相対内部を、同じ source の他の三角形 $T_2$ が横切らない**
//
// です。**隣接三角形が共有辺・共有頂点で接するのは自己交差ではありません。**
// 素朴な三角形–三角形交差判定では**全隣接対が引っかかります。**
//
// ## 健全側に倒します
//
// | | |
// |---|---|
// | **偽陽性**（適正なのに交差と報告） | **許容**。NSI を切るだけ |
// | **偽陰性**（交差を見逃す） | **許さない**。誤って NSI を有効にする |
//
// ## 新しい述語は要りません
//
// `orient3d`（$3b+5$）と `side(edge_plane, IPoint)` だけで組めます。
// **ビット幅の再導出も `widths.hpp` の追加も不要です。**
#ifndef KRISITE_MESH_SELF_INTERSECT_HPP
#define KRISITE_MESH_SELF_INTERSECT_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <vector>

#include "krisite/geom/plane.hpp"
#include "krisite/geom/predicates.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace krisite::mesh {

/// §5.6 の検査の記録。**「何件を厳密に判定し、何件を保守的に落としたか」**を分けます。
struct SelfIntersectStats {
    std::size_t pairs_broad = 0;       ///< 広域の絞り込みを通った対
    std::size_t rejected_by_side = 0;  ///< 片側にあるので即棄却
    std::size_t coplanar_exact = 0;    ///< 共面として**厳密に**判定した対
    std::size_t shared_edge_ok = 0;    ///< 辺を共有し、対頂点が厳密に片側
    std::size_t shared_vertex_ok = 0;  ///< 頂点のみ共有し、他 2 点が厳密に同じ側
    std::size_t general_exact = 0;     ///< 頂点を共有せず、辺–三角形で厳密に判定
    std::size_t conservative = 0;      ///< **厳密に判定せず、保守的に交差と報告**
    std::size_t hits = 0;              ///< 自己交差と判定した対
};

namespace detail {

/// 三角形の AABB（`int64`）。
struct TriBox {
    std::int64_t lo[3], hi[3];
};

inline TriBox tri_box(const TriMesh& m, const Tri& t) noexcept {
    TriBox b{{0, 0, 0}, {0, 0, 0}};
    for (int k = 0; k < 3; ++k) {
        const geom::IPoint& p = m.vertices[t[static_cast<std::size_t>(k)]];
        const std::int64_t c[3] = {p.x, p.y, p.z};
        for (int a = 0; a < 3; ++a) {
            if (k == 0 || c[a] < b.lo[a]) b.lo[a] = c[a];
            if (k == 0 || c[a] > b.hi[a]) b.hi[a] = c[a];
        }
    }
    return b;
}

inline bool boxes_overlap(const TriBox& x, const TriBox& y) noexcept {
    for (int a = 0; a < 3; ++a) {
        if (x.hi[a] < y.lo[a] || y.hi[a] < x.lo[a]) return false;
    }
    return true;
}

/// 点 `p`（`t` の平面上）が三角形の**どこにあるか**。辺平面の符号で見ます。
///
/// **向きは第 3 頂点で決めます** — `plane_from_edge` は向きを保証しないので、
/// **規約に頼らず実測します。**
///
///   1  … **相対内部**（厳密に中）
///   0  … 境界の上
///  -1  … 外
inline int point_vs_triangle_coplanar(const TriMesh& m, const Tri& t, const geom::PlaneD& sp,
                                      const geom::IPoint& p) noexcept {
    int on_boundary = 0;
    for (int e = 0; e < 3; ++e) {
        const geom::IPoint& a = m.vertices[t[static_cast<std::size_t>(e)]];
        const geom::IPoint& b = m.vertices[t[static_cast<std::size_t>((e + 1) % 3)]];
        const geom::IPoint& c = m.vertices[t[static_cast<std::size_t>((e + 2) % 3)]];
        const geom::PlaneD ep = geom::plane_from_edge(a, b, sp);
        const int inside = geom::side(ep, c);  // 第 3 頂点の側が「内側」
        if (inside == 0) return 0;             // 退化。**保守的に境界とみなす**
        const int sp2 = geom::side(ep, p) * inside;
        if (sp2 < 0) return -1;
        if (sp2 == 0) on_boundary = 1;
    }
    return on_boundary ? 0 : 1;
}

/// 共面の線分 `pq` が三角形 `t` の**相対内部**と交わるか。
///
/// **辺平面で分離できれば交わりません**（等号を許すので、`t` の辺そのものは分離されます）。
inline bool coplanar_segment_hits_interior(const TriMesh& m, const Tri& t, const geom::PlaneD& sp,
                                           const geom::IPoint& p, const geom::IPoint& q) noexcept {
    for (int e = 0; e < 3; ++e) {
        const geom::IPoint& a = m.vertices[t[static_cast<std::size_t>(e)]];
        const geom::IPoint& b = m.vertices[t[static_cast<std::size_t>((e + 1) % 3)]];
        const geom::IPoint& c = m.vertices[t[static_cast<std::size_t>((e + 2) % 3)]];
        const geom::PlaneD ep = geom::plane_from_edge(a, b, sp);
        const int inside = geom::side(ep, c);
        if (inside == 0) continue;
        if (geom::side(ep, p) * inside <= 0 && geom::side(ep, q) * inside <= 0) return false;
    }
    return true;
}

/// 共面な 2 三角形の**相対内部が重なるか**。分離軸（辺平面）で厳密に判定します。
///
/// **2 次元の凸多角形では、両者の辺法線を試せば分離軸定理は完全**です。
/// 共有辺で接する隣接三角形は、その辺の平面が分離軸になります（等号を許すため）。
inline bool coplanar_overlap(const TriMesh& m, const Tri& t1, const geom::PlaneD& sp1,
                             const Tri& t2) noexcept {
    const auto separates = [&](const Tri& ta, const geom::PlaneD& spa, const Tri& tb) {
        for (int e = 0; e < 3; ++e) {
            const geom::IPoint& a = m.vertices[ta[static_cast<std::size_t>(e)]];
            const geom::IPoint& b = m.vertices[ta[static_cast<std::size_t>((e + 1) % 3)]];
            const geom::IPoint& c = m.vertices[ta[static_cast<std::size_t>((e + 2) % 3)]];
            const geom::PlaneD ep = geom::plane_from_edge(a, b, spa);
            const int inside = geom::side(ep, c);
            if (inside == 0) continue;  // 退化した辺。分離軸にしない
            bool all_out = true;
            for (int k = 0; k < 3 && all_out; ++k) {
                if (geom::side(ep, m.vertices[tb[static_cast<std::size_t>(k)]]) * inside > 0) {
                    all_out = false;
                }
            }
            if (all_out) return true;
        }
        return false;
    };
    if (separates(t1, sp1, t2)) return false;
    const geom::PlaneD sp2 =
        geom::plane_from_triangle(m.vertices[t2[0]], m.vertices[t2[1]], m.vertices[t2[2]]);
    if (geom::is_degenerate(sp2)) return true;  // 退化。**保守的に重なりとみなす**
    return !separates(t2, sp2, t1);
}

/// 同次点 `h` が三角形 `t`（の平面上）の**相対内部**にあるか。
///
/// **交線の端点は有理点**になるので、整数版とは別に要ります。
inline int hpoint_vs_triangle(const TriMesh& m, const Tri& t, const geom::PlaneD& sp,
                              const geom::HPointD& h) noexcept {
    int on_boundary = 0;
    for (int e = 0; e < 3; ++e) {
        const geom::IPoint& a = m.vertices[t[static_cast<std::size_t>(e)]];
        const geom::IPoint& b = m.vertices[t[static_cast<std::size_t>((e + 1) % 3)]];
        const geom::IPoint& c = m.vertices[t[static_cast<std::size_t>((e + 2) % 3)]];
        const geom::PlaneD ep = geom::plane_from_edge(a, b, sp);
        const int inside = geom::side(ep, c);
        if (inside == 0) return 0;
        const int sh = geom::side(ep, h) * inside;
        if (sh < 0) return -1;
        if (sh == 0) on_boundary = 1;
    }
    return on_boundary ? 0 : 1;
}

/// 交線（同次座標の端点 2 つ）が三角形の**相対内部**と交わるか。
///
/// **分離軸で判定します。** ある辺平面について両端が閉じた外側にあれば交わりません。
inline bool hsegment_hits_interior(const TriMesh& m, const Tri& t, const geom::PlaneD& sp,
                                   const geom::HPointD& h0, const geom::HPointD& h1) noexcept {
    for (int e = 0; e < 3; ++e) {
        const geom::IPoint& a = m.vertices[t[static_cast<std::size_t>(e)]];
        const geom::IPoint& b = m.vertices[t[static_cast<std::size_t>((e + 1) % 3)]];
        const geom::IPoint& c = m.vertices[t[static_cast<std::size_t>((e + 2) % 3)]];
        const geom::PlaneD ep = geom::plane_from_edge(a, b, sp);
        const int inside = geom::side(ep, c);
        if (inside == 0) continue;
        if (geom::side(ep, h0) * inside <= 0 && geom::side(ep, h1) * inside <= 0) return false;
    }
    return true;
}

/// **`t2` が `t1` の相対内部と交わるか**（非共面の場合）。
///
/// **辺の交差だけを見ると落とします。** 交線の端点が両方とも `t1` の境界に載る配置が
/// あり、そこでは辺の交差がすべて境界上になります（十字に交差した 2 板で実際に踏みました）。
///
/// **交線そのものを構成して判定します。** 端点は
/// `intersect3(support1, support2, edge_plane)` で厳密に作れます。**新しい述語は要りません。**
inline bool crosses_interior(const TriMesh& m, const Tri& t1, const geom::PlaneD& sp1,
                             const Tri& t2, const geom::PlaneD& sp2, const int s2[3]) {
    geom::HPointD ends[2];
    int ne = 0;
    for (int k = 0; k < 3 && ne < 2; ++k) {
        if (s2[k] == 0) {
            ends[ne++] = geom::to_homogeneous(m.vertices[t2[static_cast<std::size_t>(k)]]);
        }
    }
    for (int k = 0; k < 3 && ne < 2; ++k) {
        const int u = k, v = (k + 1) % 3;
        if (s2[u] == 0 || s2[v] == 0 || s2[u] == s2[v]) continue;
        const geom::IPoint& p = m.vertices[t2[static_cast<std::size_t>(u)]];
        const geom::IPoint& q = m.vertices[t2[static_cast<std::size_t>(v)]];
        const geom::PlaneD ep = geom::plane_from_edge(p, q, sp2);
        ends[ne++] = geom::intersect3(sp1, sp2, ep);
    }
    if (ne == 0) return false;
    if (ne == 1) return hpoint_vs_triangle(m, t1, sp1, ends[0]) > 0;
    return hsegment_hits_interior(m, t1, sp1, ends[0], ends[1]);
}

}  // namespace detail

/// **自己交差しているか**（`SPEC-phase3.md` §5.6）。**健全側に倒します。**
///
/// 真を返したら「NSI を宣言してはいけない」。偽なら宣言してよい。
inline bool is_self_intersecting(const TriMesh& m, SelfIntersectStats* stats = nullptr) {
    SelfIntersectStats st;
    const std::size_t nt = m.triangles.size();
    if (nt < 2) {
        if (stats != nullptr) *stats = st;
        return false;
    }

    std::vector<detail::TriBox> boxes(nt);
    std::vector<geom::PlaneD> planes(nt);
    std::int64_t lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    for (std::size_t i = 0; i < nt; ++i) {
        boxes[i] = detail::tri_box(m, m.triangles[i]);
        const Tri& t = m.triangles[i];
        planes[i] = geom::plane_from_triangle(m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]]);
        for (int a = 0; a < 3; ++a) {
            if (i == 0 || boxes[i].lo[a] < lo[a]) lo[a] = boxes[i].lo[a];
            if (i == 0 || boxes[i].hi[a] > hi[a]) hi[a] = boxes[i].hi[a];
        }
    }

    // ---- 広域の絞り込み（3 次元の一様格子）----
    std::uint32_t r = 1;
    while (static_cast<std::size_t>(r) * r * r < nt && r < 64) ++r;
    std::int64_t step[3];
    for (int a = 0; a < 3; ++a) {
        const std::int64_t span = hi[a] - lo[a];
        step[a] = (span + r - 1) / r;
        if (step[a] < 1) step[a] = 1;
    }
    const auto cell_of = [&](std::int64_t x, int a) -> std::uint32_t {
        if (x <= lo[a]) return 0;
        const std::uint64_t k =
            static_cast<std::uint64_t>(x - lo[a]) / static_cast<std::uint64_t>(step[a]);
        return (k >= r) ? (r - 1) : static_cast<std::uint32_t>(k);
    };
    std::map<std::uint32_t, std::vector<std::uint32_t>> grid;
    for (std::size_t i = 0; i < nt; ++i) {
        for (std::uint32_t z = cell_of(boxes[i].lo[2], 2); z <= cell_of(boxes[i].hi[2], 2); ++z) {
            for (std::uint32_t y = cell_of(boxes[i].lo[1], 1); y <= cell_of(boxes[i].hi[1], 1);
                 ++y) {
                for (std::uint32_t x = cell_of(boxes[i].lo[0], 0); x <= cell_of(boxes[i].hi[0], 0);
                     ++x) {
                    grid[(z * r + y) * r + x].push_back(static_cast<std::uint32_t>(i));
                }
            }
        }
    }

    // ---- 対ごとの厳密な判定 ----
    std::vector<std::uint32_t> cand;
    for (std::size_t i = 0; i < nt; ++i) {
        cand.clear();
        for (std::uint32_t z = cell_of(boxes[i].lo[2], 2); z <= cell_of(boxes[i].hi[2], 2); ++z) {
            for (std::uint32_t y = cell_of(boxes[i].lo[1], 1); y <= cell_of(boxes[i].hi[1], 1);
                 ++y) {
                for (std::uint32_t x = cell_of(boxes[i].lo[0], 0); x <= cell_of(boxes[i].hi[0], 0);
                     ++x) {
                    const auto it = grid.find((z * r + y) * r + x);
                    if (it == grid.end()) continue;
                    for (std::uint32_t j : it->second) {
                        if (j > i) cand.push_back(j);
                    }
                }
            }
        }
        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

        const Tri& t1 = m.triangles[i];
        for (std::uint32_t j : cand) {
            if (!detail::boxes_overlap(boxes[i], boxes[j])) continue;
            ++st.pairs_broad;
            const Tri& t2 = m.triangles[j];
            if (geom::is_degenerate(planes[i]) || geom::is_degenerate(planes[j])) {
                ++st.conservative;
                ++st.hits;
                continue;
            }
            // **「相手の相対内部に入るか」を直接見ます**（場合分けをしません）。
            //
            // 共有辺・共有頂点は**境界どうしの接触**なので、相対内部の判定で
            // 自然に落ちます。**場合分けは不完全でした** — 例えば「頂点のみ共有」で
            // 相手の別の頂点がこちらの平面上にある配置が、立方体で普通に出ます。
            int s2[3];
            int nz = 0, npos = 0, nneg = 0;
            for (int k2 = 0; k2 < 3; ++k2) {
                s2[k2] = geom::side(planes[i], m.vertices[t2[static_cast<std::size_t>(k2)]]);
                if (s2[k2] == 0) {
                    ++nz;
                } else if (s2[k2] > 0) {
                    ++npos;
                } else {
                    ++nneg;
                }
            }
            if (nz == 0 && (npos == 3 || nneg == 3)) {  // 厳密に片側
                ++st.rejected_by_side;
                continue;
            }
            if (nz == 3) {  // 共面 → 厳密に 2D の重なりを見る
                ++st.coplanar_exact;
                if (detail::coplanar_overlap(m, t1, planes[i], t2)) ++st.hits;
                continue;
            }
            // **対称の棄却。** `t1` が `t2` の平面の厳密に片側にあれば交わりません。
            //
            // **これが分離軸の 1 本です。** 交線は `plane(t1) ∩ plane(t2)` なので、
            // `t1` の頂点が `plane(t2)` の片側に揃っていれば、**交線は `t1` を外します。**
            // **落とすと偽陽性が出ます** — 実データの最初の 1 件がこれでした
            // （`IMPL-phase5.md` §28）。三角形の辺法線だけでは分離軸定理は完全になりません。
            int s1[3];
            int nz1 = 0, npos1 = 0, nneg1 = 0;
            for (int k2 = 0; k2 < 3; ++k2) {
                s1[k2] = geom::side(planes[j], m.vertices[t1[static_cast<std::size_t>(k2)]]);
                if (s1[k2] == 0) {
                    ++nz1;
                } else if (s1[k2] > 0) {
                    ++npos1;
                } else {
                    ++nneg1;
                }
            }
            if (nz1 == 0 && (npos1 == 3 || nneg1 == 3)) {
                ++st.rejected_by_side;
                continue;
            }

            // **交線を構成して、相手の相対内部に入るかを見ます。厳密です。**
            ++st.general_exact;
            const bool hit = detail::crosses_interior(m, t1, planes[i], t2, planes[j], s2) ||
                             detail::crosses_interior(m, t2, planes[j], t1, planes[i], s1);
            if (hit) ++st.hits;
        }
        if (st.hits != 0) break;  // **1 件見つかれば十分**（宣言できない）
    }
    if (stats != nullptr) *stats = st;
    return st.hits != 0;
}

}  // namespace krisite::mesh

#endif  // KRISITE_MESH_SELF_INTERSECT_HPP

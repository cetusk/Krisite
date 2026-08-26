// Krisite — ブール演算のパイプライン
//
// SPEC-phase1.md §4（パイプライン）, §4.3（局所 arrangement）, §5（縫合）, §6（分類）
//
//   1. 平面抽出・ID 付与（§3.1）
//   2. 固定深度で八分木を構築。各セルに三角形を割り当てる（§4.2）
//   3. セルごとに局所 arrangement を計算（§4.3）
//   4. 頂点の同一性を解決して大域メッシュに縫合（§5）
//   5. 断片を正準化し、重複割り当てと共平面重複を仕分ける（§4.2, §4.3.2）
//   6. 領域に分割し、内外を分類（§6）
//   7. 演算に応じて選択・向き付けして出力
//
// **手順 4 は 5 より先です。** 断片の同一性は「辺の平面 ID の列」では決まりません。
// 相異なる 2 平面が支持平面と同一の交線を持ち得るためです（例: 支持平面 z=0 に対し
// x=0 と x+z=0 は同じ交線を与える）。したがって
//
//   - 断片の重複判定
//   - 共平面重複の検出（A の断片と B の断片が同一領域を占めるか）
//
// はどちらも**縫合後の大域頂点 ID の集合**で行います。第2段（値ベースの併合）を
// 通した ID は正準なので、平面 ID の別名に影響されません。
#ifndef KRISITE_CSG_BOOLEAN_HPP
#define KRISITE_CSG_BOOLEAN_HPP

#include <algorithm>
#include <cstdint>
#include <map>
#include <vector>
#if defined(KRISITE_DEBUG_REPRESENTATIVE)
#include <cstdio>
#endif

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

/// 断片の分類（SPEC-phase1 §4.3.2, §6.1）。
///
/// 共平面重複は「相手の境界に載っている」ことなので、内外のどちらでもありません。
/// 相手の面との**向きの一致・不一致**で扱いが変わるので 2 値に分けます。
enum class FragClass {
    Outside,           ///< 相手の外部
    Inside,            ///< 相手の内部
    CoplanarSame,      ///< 相手の面と同一領域を占め、外向き法線が同じ向き
    CoplanarOpposite,  ///< 同一領域だが外向き法線が逆向き
};

/// SPEC-phase1 §4.3.3 と §5.4 が要求する計数。
struct BoolStats {
    std::size_t fragments = 0;      ///< 正準化後の断片数（§4.3.3）
    std::size_t raw_fragments = 0;  ///< 正準化前（重複を含む）
    /// 幾何を含むセルの数（§9.0）。**分割が働いているかの直接的な指標**です。
    /// 断片数は間接的で、ケース 5 のように面がセル境界に乗ると動きません。
    std::size_t active_cells = 0;
    std::size_t total_cells = 0;  ///< 8^depth
    /// **併合グループの空間的な広がり**（§5.4、第8版）。
    ///
    /// 並列化で問題になるのは頻度ではなく局所性です。値が一致する構成点は幾何的に
    /// 同一の点なので、同じセルか、その点を共有するセルにしか存在し得ません。
    /// **すべての併合グループが 1 セルとその面・辺・頂点隣接に収まる**なら、
    /// Phase 3 は「セル並列 + 境界併合」で組めます。
    ///
    /// セル添字の差の最大（軸ごとの最大）。0 は同一セル内、1 は隣接まで。
    /// **2 以上が出たら併合の誤り**（本来別の点を同一視した）を疑ってください。
    std::size_t max_merge_span = 0;
    std::size_t merge_groups = 0;              ///< 2 点以上が併合されたグループの数
    std::size_t duplicate_fragments = 0;       ///< 重複割り当てが生んだ重複断片（§5.4）
    std::size_t coplanar_same = 0;             ///< 共平面重複のうち向きが同じ対の数
    std::size_t coplanar_opposite = 0;         ///< 向きが逆の対の数
    std::size_t constructed_points = 0;        ///< 第1段が作った構成点の総数（§5.4 の分母）
    std::size_t merged_points = 0;             ///< 第2段の併合後の点数
    std::size_t merged_by_value = 0;           ///< 第1段が取りこぼし第2段が併合した数（§5.4）
    std::size_t max_planes_at_point = 0;       ///< 1 点に集まる平面の最大枚数（§5.4、セル面込み）
    std::size_t max_mesh_planes_at_point = 0;  ///< 同上、メッシュ平面のみ（対照）
    std::size_t planes_total = 0;              ///< 総当たりの分母（表に載った平面の総数）
    std::size_t mesh_planes = 0;               ///< うちメッシュ由来
    std::size_t regions = 0;                   ///< 相異なる符号ベクトルの数（§6.1）
    std::size_t raycasts = 0;                  ///< レイキャスト回数
    std::size_t midpoint_raycasts = 0;         ///< うち中点へのフォールバック（§6.1）
    std::size_t centroid_raycasts = 0;         ///< うち 3 頂点の重心へのフォールバック
    /// `side` と `intersect3` の呼び出し数（SPEC-phase1 §12）。
    ///
    /// **`KRISITE_COUNT_PREDICATES` を定義したビルドでのみ埋まります。**
    /// 既定ビルドでは 0 のままです（計数のコストを本番に持ち込まないため）。
    std::uint64_t side_calls = 0;
    std::uint64_t intersect3_calls = 0;
};

namespace detail {

/// 頂点の平面 3 つ組（昇順に正規化）。§5.3 の第1段のキー。
inline std::array<PlaneId, 3> vertex_key(const Fragment& f, std::size_t i) {
    const std::size_t n = f.edge.size();
    std::array<PlaneId, 3> k{f.support, f.edge[(i + n - 1) % n], f.edge[i]};
    std::sort(k.begin(), k.end());
    return k;
}

/// 断片の正準キー。支持平面と、縫合後の大域頂点 ID を昇順に並べたもの。
///
/// 凸多角形は頂点集合で一意に定まるので、順序を捨てても領域は復元できます。
using RegionKey = std::pair<PlaneId, std::vector<std::uint32_t>>;

inline RegionKey region_key(PlaneId support, std::vector<std::uint32_t> ids) {
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return RegionKey{support, std::move(ids)};
}

/// 演算と分類から、その断片を出力に採るかを決める（SPEC-phase1 §2.3 の正則化）。
///
/// 共平面重複は「A 側だけを残す」ことで二重出力を避けます。どちらを残すかは
/// 恣意的ですが、向きが問題になる `Difference` では A の向きがそのまま答えなので
/// A に固定するのが自然です。
///
///  | | ∪ | ∩ | A\B |
///  |同方向| A のみ | A のみ | 両方落とす |
///  |逆方向| 両方落とす | 両方落とす | A のみ |
///
/// 同方向の面は「両者の外側が同じ側」なので、∪ と ∩ の境界には残り、A\B では消えます。
/// 逆方向の面は「一方の外側が他方の内側」なので、∪ と ∩ では内部面になって消え、
/// A\B では A の表面として残ります。
inline bool select_fragment(BoolOp op, int owner, FragClass c) noexcept {
    switch (c) {
        case FragClass::Outside:
            return (op == BoolOp::Union)          ? true
                   : (op == BoolOp::Intersection) ? false
                                                  : (owner == 0);
        case FragClass::Inside:
            return (op == BoolOp::Union)          ? false
                   : (op == BoolOp::Intersection) ? true
                                                  : (owner == 1);
        case FragClass::CoplanarSame:
            return owner == 0 && op != BoolOp::Difference;
        case FragClass::CoplanarOpposite:
            return owner == 0 && op == BoolOp::Difference;
    }
    return false;
}

}  // namespace detail

/// ブール演算。`depth` は八分木の深度（実行時パラメータ。SPEC §2.2）。
inline BoolMesh boolean_op(const mesh::TriMesh& A, const mesh::TriMesh& B, BoolOp op,
                           unsigned depth, BoolStats* stats = nullptr) {
    BoolStats st;
#if defined(KRISITE_COUNT_PREDICATES)
    geom::counters::reset();
#endif

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
    //
    // **ここでは重複を落としません。** 同一領域の断片が A と B の両方から出てくるのが
    // 共平面重複そのものなので、辺の平面 ID で潰すと区別がつかなくなります（§4.3.2）。
    const octree::UniformGrid grid(depth);
    std::vector<Fragment> frags;
    // §5.4: 併合の局所性を測るため、断片がどのセルで生まれたかを覚えておく
    std::vector<octree::CellIndex> frag_cell;

    const std::uint32_t n = grid.per_axis();
    for (std::uint32_t ci = 0; ci < n; ++ci) {
        for (std::uint32_t cj = 0; cj < n; ++cj) {
            for (std::uint32_t ck = 0; ck < n; ++ck) {
                const octree::CellIndex cell{ci, cj, ck};
                const std::size_t frags_before = frags.size();
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

                // §10.5 の変異 3: 分割平面を【割り当て集合】に変える。
                // 大域集合を使う理由は §4.3.1（継ぎ目の T 字接合）。
#if defined(KRISITE_MUTATION_LOCAL_PLANES)
                std::vector<PlaneId> local_planes;
                if (depth > 0) {
                    for (int which = 0; which < 2; ++which) {
                        const mesh::TriMesh& lm = (which == 0) ? A : B;
                        for (const Face& f : (which == 0) ? faces_a : faces_b) {
                            if (!octree::assign_to_cell(face_aabb(lm, f), grid, cell)) continue;
                            if (std::find(local_planes.begin(), local_planes.end(), f.support) ==
                                local_planes.end()) {
                                local_planes.push_back(f.support);
                            }
                        }
                    }
                    std::sort(local_planes.begin(), local_planes.end());
                }
                const std::vector<PlaneId>& split_planes = (depth > 0) ? local_planes : mesh_planes;
#else
                const std::vector<PlaneId>& split_planes = mesh_planes;
#endif

                for (int which = 0; which < 2; ++which) {
                    const mesh::TriMesh& m = (which == 0) ? A : B;
                    const std::vector<Face>& fs = (which == 0) ? faces_a : faces_b;
                    for (const Face& f : fs) {
                        // §10.5 の変異 2: 割り当てを開領域に変える。
                        // 閉領域で行う理由は §4.2（共有面上の三角形が両側から落ちる）。
#if defined(KRISITE_MUTATION_OPEN_CELLS)
                        if (depth > 0 &&
                            !octree::assign_to_cell_open(face_aabb(m, f), grid, cell)) {
                            continue;
                        }
#else
                        if (depth > 0 && !octree::assign_to_cell(face_aabb(m, f), grid, cell)) {
                            continue;
                        }
#endif
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
                        for (PlaneId q : split_planes) {
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
                            frags.push_back(std::move(p));
                            frag_cell.push_back(cell);
                        }
                    }
                }
                if (frags.size() != frags_before) ++st.active_cells;
            }
        }
    }
    st.raw_fragments = frags.size();
    st.total_cells = grid.cell_count();

    // ---- 4. 縫合（§5）----
    //
    // 第1段: 平面3つ組をキーに引く
    // 第2段: 全構成点を lex_less で整列し、値が厳密に等しいものを併合して再写像する
    //
    // **選択より先に、全断片について行います。** 断片の正準化（手順 5）が縫合後の
    // 頂点 ID を使うためです。
    std::map<std::array<PlaneId, 3>, std::uint32_t> by_key;
    std::vector<geom::HPointD> points;
    // 構成点ごとの、生成に関わったセル添字の範囲（§5.4 の局所性）
    struct CellRange {
        std::uint32_t lo[3], hi[3];
    };
    std::vector<CellRange> point_cells;

    auto vertex_id = [&](const Fragment& f, std::size_t i, const octree::CellIndex& c) {
        const auto k = detail::vertex_key(f, i);
        const std::uint32_t ci[3] = {c.i, c.j, c.k};
        auto it = by_key.find(k);
        if (it != by_key.end()) {
            CellRange& r = point_cells[it->second];
            for (int t = 0; t < 3; ++t) {
                r.lo[t] = std::min(r.lo[t], ci[t]);
                r.hi[t] = std::max(r.hi[t], ci[t]);
            }
            return it->second;
        }
        const auto id = static_cast<std::uint32_t>(points.size());
        const geom::HPointD v = fragment_vertex(table, f, i);
        KRISITE_CHECK(arith::sign(v.w) != 0,
                      "boolean_op: 構成点の w が 0（3 平面が一点で交わっていない）");
        points.push_back(v);
        point_cells.push_back(CellRange{{ci[0], ci[1], ci[2]}, {ci[0], ci[1], ci[2]}});
        by_key.emplace(k, id);
        return id;
    };

    std::vector<std::vector<std::uint32_t>> raw_polys(frags.size());
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        const Fragment& f = frags[fi];
        raw_polys[fi].reserve(vertex_count(f));
        for (std::size_t i = 0; i < vertex_count(f); ++i) {
            raw_polys[fi].push_back(vertex_id(f, i, frag_cell[fi]));
        }
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
        std::uint32_t lo[3] = {point_cells[order[i]].lo[0], point_cells[order[i]].lo[1],
                               point_cells[order[i]].lo[2]};
        std::uint32_t hi[3] = {point_cells[order[i]].hi[0], point_cells[order[i]].hi[1],
                               point_cells[order[i]].hi[2]};
        while (j < order.size() && geom::h_equal(points[order[i]], points[order[j]])) {
            remap[order[j]] = id;
            for (int t = 0; t < 3; ++t) {
                lo[t] = std::min(lo[t], point_cells[order[j]].lo[t]);
                hi[t] = std::max(hi[t], point_cells[order[j]].hi[t]);
            }
            ++j;
        }
        if (j - i > 1) {
            st.merged_by_value += (j - i - 1);
            ++st.merge_groups;
        }
        // §5.4: 併合グループの空間的な広がり。同一セルなら 0、面・辺・頂点隣接なら 1
        for (int t = 0; t < 3; ++t) {
            st.max_merge_span = std::max(st.max_merge_span, std::size_t{hi[t] - lo[t]});
        }
        i = j;
    }
#else
    // §10.5 の変異 1: 第2段を無効化する
    merged = points;
    for (std::uint32_t i = 0; i < remap.size(); ++i) remap[i] = i;
#endif
    st.merged_points = merged.size();

    // 縫合後の多角形（大域 ID）
    std::vector<std::vector<std::uint32_t>> polys(frags.size());
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        polys[fi].reserve(raw_polys[fi].size());
        for (std::uint32_t v : raw_polys[fi]) polys[fi].push_back(remap[v]);
    }

    // ---- 5. 断片の正準化: 重複割り当てと共平面重複の仕分け ----
    //
    // 同一領域を占める断片をまとめます。
    //   - 同じ owner が複数 → §4.2 の重複割り当て。1 つ残して数える
    //   - 両方の owner が居る → 共平面重複。§4.3.2 の符号 0 の分岐
    std::map<detail::RegionKey, std::array<std::vector<std::size_t>, 2>> groups;
    for (std::size_t fi = 0; fi < frags.size(); ++fi) {
        groups[detail::region_key(frags[fi].support, polys[fi])][frags[fi].owner].push_back(fi);
    }

    std::vector<std::size_t> reps;  // 代表となる断片の添字
    std::vector<FragClass> forced;  // 共平面重複で確定した分類
    std::vector<char> is_coplanar;  // forced が有効か（偽なら手順 6 で決める）
    reps.reserve(groups.size() * 2);
    for (const auto& kv : groups) {
        const auto& by_owner = kv.second;
        const bool both = !by_owner[0].empty() && !by_owner[1].empty();
        FragClass cls = FragClass::Outside;
        if (both) {
            // 外向き法線が同じ向きか。`flipped` は代表平面の法線に対する反転なので、
            // 同じ支持平面どうしなら `flipped` の一致がそのまま向きの一致になります。
            const bool same =
                frags[by_owner[0].front()].flipped == frags[by_owner[1].front()].flipped;
            cls = same ? FragClass::CoplanarSame : FragClass::CoplanarOpposite;
            if (same) {
                ++st.coplanar_same;
            } else {
                ++st.coplanar_opposite;
            }
        }
        for (int o = 0; o < 2; ++o) {
            if (by_owner[o].empty()) continue;
            st.duplicate_fragments += by_owner[o].size() - 1;
            reps.push_back(by_owner[o].front());
            forced.push_back(cls);
            is_coplanar.push_back(both ? 1 : 0);
        }
    }
    st.fragments = reps.size();

    // ---- 6. 符号ベクトルによる分類（§4.3.2, §6.1）----
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
    std::vector<char> keep(reps.size(), 0);

    for (std::size_t ri = 0; ri < reps.size(); ++ri) {
        const Fragment& f = frags[reps[ri]];
        FragClass cls = forced[ri];
        if (!is_coplanar[ri]) {
            const std::vector<PlaneId>& qs = planes_of[f.owner];
            std::vector<std::int8_t> sig(qs.size());
            for (std::size_t k = 0; k < qs.size(); ++k) {
                sig[k] = static_cast<std::int8_t>(fragment_sign(table, f, qs[k]));
            }
            const auto key = std::make_pair(f.owner, sig);
            auto it = region_inside.find(key);
            if (it == region_inside.end()) {
                // 領域の代表点: **相手の境界上に無い**頂点。
                //
                // CP1 では「相手の全平面に対して符号が非零」を条件にしていましたが、
                // これは強すぎました。断片の頂点は隣接する切断平面に必ず載るので、
                // 相手の平面が切断に使われていれば、どの頂点も条件を満たしません。
                // ケース 2 で実際に破れました（4 辺のうち 3 辺が相手の平面上）。
                //
                // 必要なのは「境界上に無いこと」だけです。相手の平面上にあっても、
                // その平面上の面の外側なら境界からは外れています。
                //
                // 頂点の分類が断片全体の分類と一致することは保証されます。$\partial B$ は
                // $B$ の平面の和集合に含まれ、断片は $B$ の平面配置のセルなので、
                // $\partial B$ は断片の**内部を横切りません**。頂点から内部へ引いた線分は
                // 途中で $\partial B$ を跨がないので、頂点が境界上に無ければ同じ側です。
                const mesh::TriMesh& other = (f.owner == 0) ? B : A;
                const std::size_t nv = vertex_count(f);
                std::vector<geom::HPointD> vs(nv);
                for (std::size_t vi = 0; vi < nv; ++vi) vs[vi] = fragment_vertex(table, f, vi);

                bool decided = false, inside = false;
                for (std::size_t vi = 0; vi < nv && !decided; ++vi) {
                    if (point_on_boundary(other, vs[vi])) continue;
                    inside = point_inside(other, vs[vi]);
                    decided = true;
                    ++st.raycasts;
                }
                // 全頂点が境界上のときは対角線の中点に落とす。
                //
                // 凸多角形の対角線の中点は相対内部にあり、相対内部は相手の平面配置の
                // セルの内側なので $\partial B$ から必ず外れます（$\partial B$ は
                // $B$ の平面の和集合に含まれ、セルの内部はどの平面とも交わらない）。
                // したがって頂点数 4 以上なら必ず決まります。
                for (std::size_t i = 0; i < nv && !decided; ++i) {
                    for (std::size_t j = i + 2; j < nv && !decided; ++j) {
                        if (i == 0 && j + 1 == nv) continue;  // 隣接（巡回）
                        const geom::HMidPointD m{vs[i], vs[j]};
                        if (point_on_boundary(other, m)) continue;
                        inside = point_inside(other, m);
                        decided = true;
                        ++st.raycasts;
                        ++st.midpoint_raycasts;
                    }
                }
                // **三角形の断片には対角線がありません。** 3 頂点の重心に落とします。
                // 凸多角形の 3 頂点が張る三角形の内部は多角形の相対内部に含まれるので、
                // 共線でない 3 つ組を選べば必ず境界から外れます。全 3 つ組を試します。
                for (std::size_t i = 0; i < nv && !decided; ++i) {
                    for (std::size_t j = i + 1; j < nv && !decided; ++j) {
                        for (std::size_t k = j + 1; k < nv && !decided; ++k) {
                            const geom::HTriPointD c{vs[i], vs[j], vs[k]};
                            if (point_on_boundary(other, c)) continue;
                            inside = point_inside(other, c);
                            decided = true;
                            ++st.raycasts;
                            ++st.centroid_raycasts;
                        }
                    }
                }
#if defined(KRISITE_DEBUG_REPRESENTATIVE)
                // 代表点が見つからないときに配置を吐く。CP2 で実際に必要になりました。
                // ライブラリ本体の既定ビルドには入りません。
                if (!decided) {
                    std::fprintf(stderr, "代表点なし: owner=%d support=%u 頂点数=%zu 辺=[", f.owner,
                                 f.support, vertex_count(f));
                    for (PlaneId e : f.edge) std::fprintf(stderr, "%u ", e);
                    std::fprintf(stderr, "] qs=[");
                    for (PlaneId q : qs) std::fprintf(stderr, "%u ", q);
                    std::fprintf(stderr, "]\n");
                    for (std::size_t vi = 0; vi < vertex_count(f); ++vi) {
                        const geom::HPointD v = fragment_vertex(table, f, vi);
                        std::fprintf(stderr, "  v%zu:", vi);
                        for (PlaneId q : qs) {
                            std::fprintf(stderr, " %u:%d", q, geom::side(table.at(q), v));
                        }
                        std::fprintf(stderr, "\n");
                    }
                }
#endif
                KRISITE_CHECK(
                    decided,
                    "boolean_op: "
                    "領域の代表点が見つからない（全頂点・全対角線中点・全重心が相手の境界上）");
                it = region_inside.emplace(key, inside).first;
                ++st.regions;
            }
            cls = it->second ? FragClass::Inside : FragClass::Outside;
        }
        keep[ri] = detail::select_fragment(op, f.owner, cls) ? 1 : 0;
    }

    // ---- 7. 出力（扇状三角形化）----
    //
    // 頂点順は「外から見て CCW」= 所有メッシュの外向き（§3.4）。`Difference` で
    // 採用する B の断片は、A\B の境界としては法線が逆になるので**順序を反転**します。
    BoolMesh out;
    out.vertices = std::move(merged);
    for (std::size_t ri = 0; ri < reps.size(); ++ri) {
        if (!keep[ri]) continue;
        std::vector<std::uint32_t> poly = polys[reps[ri]];
        if (poly.size() < 3) continue;
        if (op == BoolOp::Difference && frags[reps[ri]].owner == 1) {
            std::reverse(poly.begin(), poly.end());
        }
        for (std::size_t i = 1; i + 1 < poly.size(); ++i) {
            out.triangles.push_back({poly[0], poly[i], poly[i + 1]});
        }
    }

    // ---- §5.4: 1 点に集まる平面の最大枚数（総当たり）----
    //
    // **セル面も数えます。** §5.2 が壊れる形は {P, F, Q} と {P, F, R}（F はセル面）
    // なので、セル面を外すと肝心の同時交差を数え落とします。
    st.planes_total = table.size();
    st.mesh_planes = n_mesh_planes;
    for (const geom::HPointD& v : out.vertices) {
        std::size_t cnt = 0, mcnt = 0;
        for (PlaneId q = 0; q < static_cast<PlaneId>(table.size()); ++q) {
            if (geom::side(table.at(q), v) == 0) {
                ++cnt;
                if (q < n_mesh_planes) ++mcnt;
            }
        }
        st.max_planes_at_point = std::max(st.max_planes_at_point, cnt);
        st.max_mesh_planes_at_point = std::max(st.max_mesh_planes_at_point, mcnt);
    }

#if defined(KRISITE_COUNT_PREDICATES)
    st.side_calls = geom::counters::side_calls;
    st.intersect3_calls = geom::counters::intersect3_calls;
#endif
    if (stats) *stats = st;
    return out;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_BOOLEAN_HPP

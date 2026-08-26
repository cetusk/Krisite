// Krisite — 位相検査
//
// SPEC-phase1.md §10.1
//
//   - すべての辺がちょうど 2 枚の面に接する（watertight かつ辺多様体）
//   - すべての頂点まわりの面が単一の扇をなす（頂点多様体）
//   - 向きが大域的に整合する
//   - χ = V - E + F = 2(C - g_total)
//   - 出力が空のとき χ は未定義。V = E = F = 0 を代わりに検査する
//
// **この検査器は出力より先に用意します。** 検査器が無いと CP1 の合否を判定できません。
//
// 幾何は一切見ません。**純粋に組合せ的**な検査なので、頂点座標の表現（IPoint か
// 構成点か）に依存しません。§10.2.1 の深度不変性がこの上に載ります。
#ifndef KRISITE_MESH_TOPOLOGY_HPP
#define KRISITE_MESH_TOPOLOGY_HPP

#include <algorithm>
#include <cstdint>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "krisite/mesh/tri_mesh.hpp"

namespace krisite::mesh {

/// 位相検査の結果。
struct TopologyReport {
    bool empty = false;  ///< V = E = F = 0

    std::size_t v = 0;  ///< 参照されている頂点数
    std::size_t e = 0;  ///< 無向辺数
    std::size_t f = 0;  ///< 面数

    bool edge_manifold = false;    ///< すべての無向辺がちょうど 2 面に接する
    bool vertex_manifold = false;  ///< すべての頂点まわりが単一の扇
    bool oriented = false;         ///< 隣接 2 面が共有辺を逆向きに辿る
    bool no_degenerate = false;    ///< 同一頂点を 2 度使う三角形が無い

    std::size_t components = 0;  ///< 連結成分（シェル）数 C
    long long chi = 0;           ///< V - E + F
    long long genus_total = 0;   ///< C - χ/2
    bool chi_even = false;       ///< χ が偶数（そうでなければ g が半整数になり異常）

    /// **非多様体辺の次数**（SPEC-phase1 §10.1 の恒久的診断）。
    ///
    /// $\chi$ は「面が余分」なのか「面が欠けている」のかを区別しません。次数は区別します。
    ///
    ///   次数 1 … その辺の片側に面が無い → **面が欠けている**
    ///   次数 3 以上 … 同じ辺に面が集まりすぎ → **面が余分**、または非多様体な接触
    ///
    /// 落ちたときはこの 2 つを必ず記録してください。デバッグの初手が変わります。
    std::size_t edges_deficient = 0;  ///< 1 面にしか接しない辺の数
    std::size_t edges_excess = 0;     ///< 3 面以上に接する辺の数
    std::size_t max_edge_degree = 0;  ///< 辺の次数の最大

    /// **落ちている頂点のリンク構造**（SPEC-phase1 §9.3.2）。
    ///
    /// 頂点 $v$ に接する三角形 $(v,a,b)$ はリンク上の有向辺 $a \to b$ を与えます。
    /// 辺多様体な閉曲面ならリンクは閉路の直和になり、その本数が**扇の数** $k_v$ です。
    ///
    ///   $k_v = 1$ … 多様体
    ///   $k_v \ge 2$ … 錐が $k_v$ 個出会っている
    ///
    /// **$k_v$ 分裂すると $\chi$ は $k_v - 1$ 増えます**（頂点が $k_v - 1$ 個増え、
    /// 辺と面は変わらない）。`chi_after_split` はその予測値です。
    ///
    /// **扇の数だけでは「ピンチ点」と「円環」を区別できません。** どちらも曲面上では
    /// 閉路 2 本に見えます。区別の手がかりは連結性で、扇が別々の連結成分に属するなら
    /// 錐は分離しており、分裂で多様体化できます。同一成分なら円環の疑いがあり、
    /// **その場合は頂点を複製しても円板になりません**（§9.3.2）。
    std::size_t nonmanifold_vertices = 0;  ///< リンクが単一閉路でない頂点の数
    std::size_t max_vertex_fans = 0;       ///< 頂点まわりの扇の最大数
    long long extra_fans = 0;              ///< Σ(k_v − 1)。分裂で増える χ
    long long chi_after_split = 0;         ///< χ + extra_fans（分裂後の予測 χ）
    /// 扇が**同一の連結成分**に属する頂点の数。**円環の疑いはここに出ます。**
    std::size_t vertices_fans_in_one_component = 0;

    /// SPEC §10.1 のすべてを満たすか。空メッシュは V=E=F=0 だけを見る。
    bool ok() const noexcept {
        if (empty) return v == 0 && e == 0 && f == 0;
        return edge_manifold && vertex_manifold && oriented && no_degenerate && chi_even;
    }
};

namespace detail {

inline std::pair<VertexId, VertexId> undirected(VertexId a, VertexId b) noexcept {
    return (a < b) ? std::pair<VertexId, VertexId>{a, b} : std::pair<VertexId, VertexId>{b, a};
}

/// 素集合（面の連結成分を数える）。
class DisjointSet {
public:
    explicit DisjointSet(std::size_t n) : parent_(n) {
        for (std::size_t i = 0; i < n; ++i) parent_[i] = i;
    }
    std::size_t find(std::size_t x) {
        while (parent_[x] != x) {
            parent_[x] = parent_[parent_[x]];
            x = parent_[x];
        }
        return x;
    }
    void unite(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a != b) parent_[a] = b;
    }

private:
    std::vector<std::size_t> parent_;
};

}  // namespace detail

/// 三角形リストの位相を検査する。頂点座標は見ない。
inline TopologyReport check_topology(const std::vector<Tri>& tris) {
    TopologyReport r;
    r.f = tris.size();
    if (tris.empty()) {
        r.empty = true;
        return r;
    }

    // ---- 退化三角形 ----
    r.no_degenerate = true;
    for (const Tri& t : tris) {
        if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) r.no_degenerate = false;
    }

    // ---- 辺の集計 ----
    // 無向辺 → (接する面数, 有向 (a,b) の出現数, 有向 (b,a) の出現数)
    struct EdgeInfo {
        int count = 0;
        int fwd = 0;  // (u,v) の向き（u < v）
        int bwd = 0;  // (v,u) の向き
        std::size_t face[2] = {0, 0};
    };
    std::map<std::pair<VertexId, VertexId>, EdgeInfo> edges;
    std::set<VertexId> used;

    for (std::size_t fi = 0; fi < tris.size(); ++fi) {
        const Tri& t = tris[fi];
        for (int k = 0; k < 3; ++k) used.insert(t[k]);
        for (int k = 0; k < 3; ++k) {
            const VertexId a = t[k], b = t[(k + 1) % 3];
            EdgeInfo& ei = edges[detail::undirected(a, b)];
            if (ei.count < 2) ei.face[ei.count] = fi;
            ++ei.count;
            if (a < b) {
                ++ei.fwd;
            } else {
                ++ei.bwd;
            }
        }
    }
    r.v = used.size();
    r.e = edges.size();

    // ---- 辺多様体と向きの整合 ----
    r.edge_manifold = true;
    r.oriented = true;
    for (const auto& kv : edges) {
        const EdgeInfo& ei = kv.second;
        if (ei.count != 2) r.edge_manifold = false;
        if (ei.count == 1) ++r.edges_deficient;
        if (ei.count >= 3) ++r.edges_excess;
        r.max_edge_degree = std::max(r.max_edge_degree, static_cast<std::size_t>(ei.count));
        // 向きが整合していれば、共有辺は一方が (u,v)、他方が (v,u)
        if (!(ei.fwd == 1 && ei.bwd == 1)) r.oriented = false;
    }

    // ---- 連結成分（辺を共有する面どうしを併合）----
    detail::DisjointSet ds(tris.size());
    for (const auto& kv : edges) {
        if (kv.second.count == 2) ds.unite(kv.second.face[0], kv.second.face[1]);
    }
    std::set<std::size_t> roots;
    for (std::size_t i = 0; i < tris.size(); ++i) roots.insert(ds.find(i));
    r.components = roots.size();

    // ---- 頂点多様体（頂点まわりのリンクが単一の閉路）----
    //
    // 頂点 v に接する三角形 (v,a,b) は、リンク上の有向辺 a→b を与える。
    // 単一の扇 ⟺ 各頂点の出次数・入次数が 1 で、全体が 1 本の閉路になる。
    {
        // 三角形の添字も持つ。扇がどの連結成分に属するかを見るため（§9.3.2）
        struct Arc {
            VertexId from, to;
            std::size_t tri;
        };
        std::map<VertexId, std::vector<Arc>> link;
        for (std::size_t ti = 0; ti < tris.size(); ++ti) {
            const Tri& t = tris[ti];
            for (int k = 0; k < 3; ++k) {
                const VertexId v = t[k], a = t[(k + 1) % 3], bb = t[(k + 2) % 3];
                link[v].push_back(Arc{a, bb, ti});
            }
        }
        r.vertex_manifold = true;
        for (const auto& kv : link) {
            const auto& arcs = kv.second;
            std::map<VertexId, VertexId> next;
            std::map<VertexId, std::size_t> arc_tri;
            std::map<VertexId, int> indeg;
            bool ok = true;
            for (const Arc& e : arcs) {
                if (!next.emplace(e.from, e.to).second) ok = false;  // 出次数 > 1
                arc_tri.emplace(e.from, e.tri);
                ++indeg[e.to];
            }
            if (!ok) {
                r.vertex_manifold = false;
                ++r.nonmanifold_vertices;
                continue;
            }
            for (const auto& d : indeg) {
                if (d.second != 1) ok = false;
            }
            if (!ok || next.size() != arcs.size()) {
                r.vertex_manifold = false;
                ++r.nonmanifold_vertices;
                continue;
            }
            // リンクは閉路の直和になっている。本数（= 扇の数）を数える
            std::set<VertexId> visited;
            std::size_t fans = 0;
            std::set<std::size_t> fan_components;
            for (const auto& e : next) {
                if (visited.count(e.first)) continue;
                ++fans;
                VertexId cur = e.first;
                std::size_t steps = 0;
                do {
                    visited.insert(cur);
                    fan_components.insert(ds.find(arc_tri[cur]));
                    auto it = next.find(cur);
                    if (it == next.end()) {
                        ok = false;
                        break;
                    }
                    cur = it->second;
                    ++steps;
                } while (cur != e.first && steps <= next.size());
                if (!ok || cur != e.first) break;
            }
            if (!ok || visited.size() != next.size()) {
                r.vertex_manifold = false;
                ++r.nonmanifold_vertices;
                continue;
            }
            r.max_vertex_fans = std::max(r.max_vertex_fans, fans);
            if (fans != 1) {
                r.vertex_manifold = false;
                ++r.nonmanifold_vertices;
                r.extra_fans += static_cast<long long>(fans) - 1;
                // 扇が 1 つの連結成分に収まっている = 錐が分離していない可能性
                if (fan_components.size() == 1) ++r.vertices_fans_in_one_component;
            }
        }
    }

    // ---- オイラー標数と種数 ----
    r.chi = static_cast<long long>(r.v) - static_cast<long long>(r.e) + static_cast<long long>(r.f);
    r.chi_even = (r.chi % 2 == 0);
    // χ = 2(C - g_total) → g_total = C - χ/2
    r.genus_total = static_cast<long long>(r.components) - r.chi / 2;
    // 分裂後の予測 χ（§9.3.2）。k 分裂で頂点が k-1 個増え、辺と面は変わらない
    r.chi_after_split = r.chi + r.extra_fans;
    return r;
}

inline TopologyReport check_topology(const TriMesh& m) {
    return check_topology(m.triangles);
}

}  // namespace krisite::mesh

#endif  // KRISITE_MESH_TOPOLOGY_HPP

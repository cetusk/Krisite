// Krisite — 接触の分裂（非多様体出力の意味論）
//
// SPEC-phase2.md §5
//
// **正則化ブールは一般に多様体出力を保証できません**（`SPEC-phase1.md` §9.3.2）。
// 辺だけ・頂点だけを共有する接触が出力に残ります。Phase 2 はこれを**分裂**させて
// 多様体化します。既定は分裂側です（§5.2。Phase 5 のメッシュ化と GWN が多様体を
// 前提にできると扱いが楽なため）。
//
// **壊れる検査は接触の次元で決まるので、分裂も 2 種類が要ります**（§5.1）。
//
//   辺（次数 4）        `edge_manifold`     辺と、その上の頂点を分裂
//   頂点（k 個の扇）    `vertex_manifold`   頂点を k 個に分裂
//
// ---
//
// ## 実装: 頂点まわりの「扇」で分ける
//
// **両方を 1 つの規則で扱えます。** 頂点 $v$ に接する三角形を、
// **次数 2 の辺だけを辿って**併合します。その同値類が扇です。
//
//   扇が 1 個   多様体な頂点。分裂しない
//   扇が k 個   $k$ 個に複製する
//
// **次数 4 の辺は扇を繋ぎません**（次数 2 でないため）。したがって辺の両端で扇が
// 分かれ、**辺そのものが 2 本に複製されます。** §5.1.1 が言う「辺 → 頂点の順」は、
// この規則では順序の問題が起きません。**扇の計算が辺の次数を見ているからです。**
//
// > **§5.1.2 は「面の連結成分で分ける」と書いています。** ここでは連結性を
// > **頂点のまわりで局所的に**見ています。大域の連結成分より細かく、
// > 「2 枚のシートが別の場所で繋がっている」配置（§5.1.2.1）でも分けられます。
// > **owner では分けません**（§5.1.2）。自己接触で組がシートをまたぎます。
//
// **分裂しても多様体にならない配置は停止させます**（§5.1.2.1）。黙って片方に寄せると
// §5.5.1 の $C$ 不変性で検出はされますが、原因の特定に手間がかかります。
#ifndef KRISITE_MESH_SPLIT_HPP
#define KRISITE_MESH_SPLIT_HPP

#include <chrono>
#include <cstdint>
#include <map>
#include <mutex>
#include <utility>
#include <vector>

#include "krisite/mesh/topology.hpp"
#include "krisite/mesh/tri_mesh.hpp"
#include "krisite/par/thread_pool.hpp"

namespace krisite::mesh {

/// §5.5 の予測と実測（**事前に予測して事後と突き合わせる**）。
///
/// $\Delta\chi = \Delta V - \Delta E$ は $\Delta F = 0$ からの恒等式なので、
/// **そのまま書いても検査になりません。** 分裂の**前に**診断から予測します。
struct SplitStats {
    // ---- 診断から立てた予測 ----
    std::size_t predicted_delta_v = 0;  ///< Σ(扇の数 − 1)
    std::size_t predicted_delta_e = 0;  ///< Σ(辺の次数/2 − 1)。次数 4 の辺は 1 本増える
    long long predicted_delta_chi = 0;  ///< 上の 2 つの差

    // ---- 分裂後の実測 ----
    std::size_t actual_delta_v = 0;
    std::size_t actual_delta_e = 0;
    long long actual_delta_chi = 0;

    std::size_t split_vertices = 0;    ///< 分裂した頂点の数
    std::size_t max_fans = 0;          ///< 扇の個数 k の最大
    std::size_t excess_edges = 0;      ///< 次数 3 以上の辺の本数
    std::size_t excess_endpoints = 0;  ///< その端点の数（重複を除く）
    /// **分裂しても多様体にならなかった配置の数**（§5.1.2.1）。
    /// **radial sort を実装する必要性の判断材料です。** 0 でなければ報告すること。
    std::size_t unresolved = 0;
    /// **early-out で arrangement を省いたセル由来の三角形に接する頂点が分裂した回数**
    /// （SPEC-phase2 §13 の CP5）。
    ///
    /// CP5 の相互作用「early-out × 分裂」が**実際に通ったこと**の指標です。
    /// 省いたセルの断片は分割されないまま分裂に回るので、単独の検査では通りません。
    std::size_t split_from_early_out = 0;

    // ---- 段ごとの時間（`SPEC-phase4.md` §6.2。**逐次部分を見るため**） ----
    double ms_edges = 0.0;  ///< 辺→三角形 / 頂点→三角形の表を作る（逐次）
    double ms_diag = 0.0;   ///< 診断と成分分け（逐次）
    double ms_fan = 0.0;    ///< 扇の計算（**並列**）
    double ms_apply = 0.0;  ///< ID の割り当てと書き戻し（逐次）
    /// §5.5 の検算（`check_topology` を 2 回）。**仕事ではなく検査です。**
    /// 分けて計上しないと、並列化の余地を見誤ります（`IMPL-phase4.md` §2.3）
    double ms_verify = 0.0;
};

namespace detail {

/// 頂点まわりの三角形をまとめる素集合。
class SmallDsu {
public:
    explicit SmallDsu(std::size_t n) : parent_(n) {
        for (std::size_t i = 0; i < n; ++i) parent_[i] = i;
    }
    std::size_t find(std::size_t x) {
        while (parent_[x] != x) x = parent_[x] = parent_[parent_[x]];
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

/// 接触を分裂させ、新しい三角形列を返す（§5.1）。
///
/// 新しく作った頂点の元 ID は `origin` に積みます（座標を複製するのは呼び出し側）。
///
/// `owner` は **§9.3 の変異 9 のためだけ**に受け取ります。既定（`nullptr`）では使いません。
/// **owner で分けてはいけません**（§5.1.2）。$A \setminus B$ のように結果が自分自身に
/// 接触する配置で、組がシートをまたぎます。
inline std::vector<Tri> split_contacts(const std::vector<Tri>& tris, std::size_t vertex_count,
                                       std::vector<std::uint32_t>* origin,
                                       SplitStats* stats = nullptr,
                                       const std::vector<int>* owner = nullptr,
                                       const std::vector<char>* from_early_out = nullptr,
                                       par::ThreadPool* pool = nullptr) {
    using clk = std::chrono::steady_clock;
    const auto ms_since = [](clk::time_point a) {
        return std::chrono::duration<double, std::milli>(clk::now() - a).count();
    };
    auto t_edges = clk::now();

    // ---- 無向辺 → 接する三角形 ----
    std::map<std::pair<VertexId, VertexId>, std::vector<std::size_t>> edge_tris;
    for (std::size_t t = 0; t < tris.size(); ++t) {
        for (int k = 0; k < 3; ++k) {
            edge_tris[detail::undirected(tris[t][k], tris[t][(k + 1) % 3])].push_back(t);
        }
    }

    // ---- 頂点 → 接する三角形 ----
    std::vector<std::vector<std::size_t>> at(vertex_count);
    for (std::size_t t = 0; t < tris.size(); ++t) {
        for (int k = 0; k < 3; ++k) at[tris[t][k]].push_back(t);
    }

    SplitStats st;
    st.ms_edges = ms_since(t_edges);
    auto t_diag = clk::now();
    // 次数 3 以上の辺の診断（§5.5 の予測に使う）
    std::vector<char> endpoint_of_excess(vertex_count, 0);
    for (const auto& kv : edge_tris) {
        if (kv.second.size() <= 2) continue;
        ++st.excess_edges;
        st.predicted_delta_e += kv.second.size() / 2 - 1;
        endpoint_of_excess[kv.first.first] = 1;
        endpoint_of_excess[kv.first.second] = 1;
    }
    for (std::size_t v = 0; v < vertex_count; ++v) {
        st.excess_endpoints += static_cast<std::size_t>(endpoint_of_excess[v]);
    }

    // ---- 次数 3 以上の辺のまわりで、面を【連結成分】で組にする（§5.1.2）★ ----
    //
    // **局所の連結性だけでは足りません。** 接触辺の途中の頂点では、同じ立体の 2 面が
    // その接触辺でしか繋がっていないことがあります（深度 ≥2 の 11b がまさにこれ）。
    // 局所規則だと 4 つの扇に分かれ、**辺が裂けて次数 1 の辺が出ます。実際に踏みました。**
    //
    // **owner で分けてはいけません**（§5.1.2）。$A \setminus B$ のように結果が自分自身に
    // 接触する配置では、各シートが A の面 1 枚と B の面 1 枚から成り、組がシートを
    // またぎます。**正しい基準は面の連結成分です。**
    detail::SmallDsu comp(tris.size());
    for (const auto& kv : edge_tris) {
        if (kv.second.size() == 2) comp.unite(kv.second[0], kv.second[1]);
    }
    // 過剰な辺 → 「同じ組に属する三角形の対」
    std::map<std::pair<VertexId, VertexId>, std::vector<std::vector<std::size_t>>> edge_groups;
    for (const auto& kv : edge_tris) {
        if (kv.second.size() <= 2) continue;
        std::map<std::size_t, std::vector<std::size_t>> by_key;
        for (std::size_t t : kv.second) {
            const std::size_t key =
                (owner != nullptr) ? static_cast<std::size_t>((*owner)[t]) : comp.find(t);
            by_key[key].push_back(t);
        }
        bool ok = by_key.size() == kv.second.size() / 2;
        for (const auto& gk : by_key) {
            if (gk.second.size() != 2) ok = false;
        }
        if (!ok) {
            // **§5.1.2.1: 推測せず停止します。** 2 枚のシートが別の場所で繋がっていると
            // 4 枚が同一成分に入り、連結成分では分けられません。黙って片方に寄せると
            // §5.5.1 の C 不変性で検出はされますが、原因の特定に手間がかかります
            ++st.unresolved;
            continue;
        }
        for (const auto& gk : by_key) edge_groups[kv.first].push_back(gk.second);
    }

    // ---- 頂点ごとに扇を数え、複製の割り当てを決める ----
    //
    // **扇の計算は頂点ごとに独立です**（`SPEC-phase4.md` §3）。読むのは `tris` /
    // `edge_tris` / `edge_groups` だけで、書くのは自分のスロットだけ。
    //
    // **ID の割り当ては逐次で行います。** `next_id++` は共有の状態なので、
    // 並列に配ると番号がスケジューリングに依存します（§4.2）。
    // **頂点の順に配れば、スレッド数に依らず同じ番号になります。**
    std::vector<Tri> out = tris;
    std::vector<std::uint32_t> new_origin;
    auto next_id = static_cast<std::uint32_t>(vertex_count);

    /// 頂点ごとの扇の割り当て（`fan_of[i]` = `at[v][i]` が属する扇の番号）。
    struct FanPlan {
        std::vector<std::uint32_t> fan_of;
        std::size_t fans = 0;
    };
    std::vector<FanPlan> plan(vertex_count);
    st.ms_diag = ms_since(t_diag);

    // **スレッド局所のスクラッチ**（`SPEC-phase4.md` §1.1）。
    //
    // **頂点ごとに `std::map` と `std::vector` を作ってはいけません。** 頂点数が
    // 数万になると、確保だけで並列化の利得が消えます（実測: 並列化しても 1.19 倍。
    // `IMPL-phase4.md` §2.3）。**器を使い回せば確保は 1 スレッドあたり数回**です。
    struct Scratch {
        std::vector<std::size_t> parent;  ///< DSU
        std::vector<std::size_t> root;    ///< 初出順の根
    };
    const unsigned nthreads = (pool != nullptr) ? pool->size() : 1u;
    std::vector<Scratch> tl(nthreads);

    const auto compute = [&](std::size_t v, unsigned tid) {
        const std::vector<std::size_t>& inc = at[v];
        if (inc.size() < 2) return;
#if defined(KRISITE_MUTATION_NO_EDGE_SPLIT)
        // SPEC-phase2 §9.3 の変異 8a: **辺の分裂だけを無効化**（頂点分裂は残す）。
        // **頂点分裂だけでは 11b が直りません**（§5.1）。扇の最大が 1 で適用対象がない
        if (endpoint_of_excess[v]) return;
#endif
#if defined(KRISITE_MUTATION_NO_VERTEX_SPLIT)
        // 変異 8b: **頂点の分裂だけを無効化**（辺分裂は残す）。4T / 4T′ が直らなくなる
        if (!endpoint_of_excess[v]) return;
#endif
        Scratch& sc = tl[tid];
        const std::size_t n = inc.size();
        // **`inc` は小さい**（頂点まわりの三角形数）ので、線形探索が表より速い
        const auto index_of = [&inc, n](std::size_t t) {
            for (std::size_t i = 0; i < n; ++i) {
                if (inc[i] == t) return i;
            }
            return n;
        };
        sc.parent.resize(n);
        for (std::size_t i = 0; i < n; ++i) sc.parent[i] = i;
        const auto find = [&sc](std::size_t x) {
            while (sc.parent[x] != x) x = sc.parent[x] = sc.parent[sc.parent[x]];
            return x;
        };
        const auto unite = [&sc, &find](std::size_t a, std::size_t b) {
            a = find(a);
            b = find(b);
            if (a != b) sc.parent[a] = b;
        };

        for (std::size_t t : inc) {
            for (int k = 0; k < 3; ++k) {
                const VertexId a = tris[t][k], b = tris[t][(k + 1) % 3];
                if (a != v && b != v) continue;  // v に接する辺だけ
                const auto key = detail::undirected(a, b);
                // **`operator[]` を使わないこと。** 無ければ挿入するので、
                // 並列に読むと表を壊します（`const` な `find` で引きます）
                const auto eit = edge_tris.find(key);
                if (eit == edge_tris.end()) continue;
                const std::vector<std::size_t>& sh = eit->second;
                if (sh.size() == 2) {
                    // 多様体な辺は扇を繋ぐ
                    const std::size_t i0 = index_of(sh[0]), i1 = index_of(sh[1]);
                    if (i0 < n && i1 < n) unite(i0, i1);
                } else {
                    // **過剰な辺は「同じ組」の中だけ繋ぐ**（§5.1.2）
                    const auto it = edge_groups.find(key);
                    if (it == edge_groups.end()) continue;  // 分けられなかった辺
                    for (const auto& grp : it->second) {
                        const std::size_t i0 = index_of(grp[0]);
                        if (i0 >= n) continue;
                        for (std::size_t j = 1; j < grp.size(); ++j) {
                            const std::size_t ij = index_of(grp[j]);
                            if (ij < n) unite(i0, ij);
                        }
                    }
                }
            }
        }

        // **扇の番号は `inc` の順で初出順に振ります**（正準な順序）
        sc.root.clear();
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t r = find(i);
            bool seen = false;
            for (std::size_t k = 0; k < sc.root.size() && !seen; ++k) {
                if (sc.root[k] == r) seen = true;
            }
            if (!seen) sc.root.push_back(r);
        }
        // **扇が 1 つなら何も残しません。** ほとんどの頂点がこれで、
        // ここで確保しないことが効きます
        if (sc.root.size() <= 1) return;

        FanPlan& fp = plan[v];
        fp.fans = sc.root.size();
        fp.fan_of.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t r = find(i);
            for (std::size_t k = 0; k < sc.root.size(); ++k) {
                if (sc.root[k] == r) {
                    fp.fan_of[i] = static_cast<std::uint32_t>(k);
                    break;
                }
            }
        }
    };

    /// **1 頂点分の ID を配ります。** 呼ぶ順序が番号を決めるので、
    /// 既定では**頂点の順**に呼びます（§4.2。スレッド数に依らず同じ番号）。
    const auto assign_ids = [&](std::size_t v) {
        const FanPlan& fp = plan[v];
        if (fp.fans <= 1) return;
        const std::vector<std::size_t>& inc = at[v];

        ++st.split_vertices;
        st.max_fans = std::max(st.max_fans, fp.fans);
        // §13 の CP5:「early-out × 分裂」を通った回数
        if (from_early_out != nullptr) {
            for (std::size_t i = 0; i < inc.size(); ++i) {
                if (inc[i] < from_early_out->size() && (*from_early_out)[inc[i]] != 0) {
                    ++st.split_from_early_out;
                    break;
                }
            }
        }
        st.predicted_delta_v += fp.fans - 1;

        std::vector<std::uint32_t> id_of(fp.fans, 0);
        id_of[0] = static_cast<std::uint32_t>(v);  // 最初の扇は元の ID を使う
        for (std::size_t k = 1; k < fp.fans; ++k) {
            id_of[k] = next_id++;
            new_origin.push_back(static_cast<std::uint32_t>(v));
        }
        for (std::size_t i = 0; i < inc.size(); ++i) {
            const std::uint32_t id = id_of[fp.fan_of[i]];
            for (auto& vid : out[inc[i]]) {
                if (vid == v) vid = id;
            }
        }
    };

#if defined(KRISITE_MUTATION_NO_EXIT_BARRIER)
    // 変異 23: **扇の計算と ID の割り当ての間のバリアを外す**（`SPEC-phase4.md` §7.5）。
    //
    // 排他は取るので**競合はありません。** 番号が**完了した順**に配られるだけです。
    // 幾何は同じなので位相も体積も区別せず、**決定性でしか捕まりません。**
    std::mutex mut_id;
    const auto compute_fused = [&](std::size_t v, unsigned tid) {
        compute(v, tid);
        const FanPlan& fp = plan[v];
        if (fp.fans <= 1) return;
        const std::lock_guard<std::mutex> g(mut_id);
        assign_ids(v);
    };
#endif
    auto t_fan = clk::now();
#if defined(KRISITE_MUTATION_NO_EXIT_BARRIER)
    if (pool != nullptr) {
        pool->run(vertex_count, compute_fused);
    } else {
        for (std::size_t v = 0; v < vertex_count; ++v) compute_fused(v, 0u);
    }
    st.ms_fan = ms_since(t_fan);
    auto t_apply = clk::now();
#else
    if (pool != nullptr) {
        pool->run(vertex_count, compute);
    } else {
        for (std::size_t v = 0; v < vertex_count; ++v) compute(v, 0u);
    }
    st.ms_fan = ms_since(t_fan);
    auto t_apply = clk::now();

    // ---- ID の割り当て（**逐次**。頂点の順に配るので決定的）----
    for (std::size_t v = 0; v < vertex_count; ++v) assign_ids(v);
#endif

    st.predicted_delta_chi =
        static_cast<long long>(st.predicted_delta_v) - static_cast<long long>(st.predicted_delta_e);

    st.ms_apply = ms_since(t_apply);
    auto t_verify = clk::now();

    // ---- §5.5 の検算: 予測と実測を突き合わせる ----
    const TopologyReport before = check_topology(tris);
    const TopologyReport after = check_topology(out);
    st.actual_delta_v = after.v - before.v;
    st.actual_delta_e = after.e - before.e;
    st.actual_delta_chi = after.chi - before.chi;
    // **面は増えません**（§5.1.3）。境界の点集合は変わらず、組合せ的な表現だけが変わります
    KRISITE_CHECK(after.f == before.f, "split_contacts: 面が増減した（§5.1.3 に反する）");
    // **分裂しても多様体にならない配置は停止**（§5.1.2.1）。
    // 推測して片方に寄せると §5.5.1 の C 不変性で検出はされますが、原因の特定に手間が
    // かかります。**到達した配置を記録してください**（radial sort の必要性の判断材料）
    if (!before.empty && !(after.edge_manifold && after.vertex_manifold)) ++st.unresolved;

    st.ms_verify = ms_since(t_verify);
    if (origin != nullptr) *origin = std::move(new_origin);
    if (stats != nullptr) *stats = st;
    return out;
}

}  // namespace krisite::mesh

#endif  // KRISITE_MESH_SPLIT_HPP

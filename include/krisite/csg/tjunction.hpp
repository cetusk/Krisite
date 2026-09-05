// Krisite — T 頂点の解決
//
// SPEC-phase2.md §2.4.2（案 D）, §2.4.3（解決手順）, §2.4.4（実装上の規律）
//
// **適応分割では隣接セルの深さが違うので、共有面の上で辺の分かれ方が食い違います。**
// 細かい側は共有面 F 上の断片を内部のセル面平面で切りますが、粗い側は F を 1 枚として
// 扱います。F 上で粗い側の 1 本の辺が細かい側の 2 本の辺に対応する — T 字接合です。
//
// **グリッド平面を両側で共有する道は閉じています**（§2.4.1）。共有面は正の面積を持つので
// 十分細かいレベルのグリッド平面は必ずそこを横切ります。したがって「レベル L(C) 以下を
// 使う」という局所的な規則は、隣接する葉すべてで L が等しいことを要求し、葉の隣接グラフが
// 連結なので L は大域定数に潰れます。
//
// **制約は `Fragment` の辺平面列にあります。** 辺 P の途中に G との交点を入れようとすると
// 辺列が [..., P, G, P, ...] になり、P∩G と G∩P が同一点になって長さ 0 の辺ができます。
// 壊れるのは凸性ではありません（辺の途中に頂点を入れても多角形は弱く凸のままで、
// 扇分割も使えます）。
//
// **この制約は縫合後の多角形（頂点 ID 列）には効きません。** そこでセルごとの arrangement は
// 各セル自身の 6 面だけで行い、**縫合の後に T 頂点を入れます**（案 D）。
//
// **挿入する頂点は既に存在します。** 粗い側が入れる頂点は、細かい側が既に作った頂点
// そのものです。平面3つ組 {P, Q, G} が一致するので SPEC-phase1 §5.3 の第1段で同じ ID に
// なっており、「挿入」は既存 ID を並べ直すだけです。**新しい幾何を作りません。**
#ifndef KRISITE_CSG_TJUNCTION_HPP
#define KRISITE_CSG_TJUNCTION_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

#include "krisite/csg/plane_table.hpp"
#include "krisite/geom/point.hpp"
#include "krisite/geom/predicates.hpp"
#include "krisite/octree/adaptive.hpp"
#include "krisite/par/thread_pool.hpp"

namespace krisite::csg {

/// §2.4.3 / §11 の記録。
struct TJunctionStats {
    std::size_t inserted = 0;      ///< 挿入された T 頂点の数
    std::size_t candidates = 0;    ///< 索引が返した候補の数（絞り込みの効き）
    std::size_t max_per_edge = 0;  ///< 1 本の辺に入った T 頂点の最大数
    /// **残した**面積 0 の三角形の枚数（§2.4.4 (2)）。**捨てません。**
    /// 起点を選べなかった多角形でだけ出ます
    std::size_t degenerate_kept = 0;
    std::size_t apex_fallback = 0;  ///< 起点を選べなかった多角形の数
    std::size_t edges_scanned = 0;  ///< 走査した辺の数（候補数の分母）
    /// **保持された構成点が T 頂点として挿入された回数**（§13 の CP5）。
    ///
    /// CP5 の相互作用「構成点の保持 × T 解決」が**実際に通ったこと**の指標です。
    /// 0 のままなら、その経路を 1 度も踏まずに CP5 が緑になっています。
    std::size_t inserted_from_cache = 0;
};

/// **平面 → その平面上にある大域頂点 ID** の索引（§2.4.3 の手順 1）。
///
/// 多角形の辺 (a,b) は 2 平面 (support, edge) の交線上にあります。したがって候補は
/// **両方の平面に載っている頂点**です。判定は `side(plane, HPoint) == 0`、
/// つまり **Phase 0 からある述語 1 本**で済みます（§2.4.3「新しい述語は要りません」）。
///
/// > **平面3つ組をキーにしてはいけません。** 3つ組は正準ではありません。
/// > **相異なる 2 平面が支持平面と同一の交線を与え得る**ためです
/// > （`IMPL-phase1.md` §2.9。例: 支持平面 z=0 に対し x=0 と x+z=0 はどちらも
/// > 同じ交線を与える）。同じ幾何辺が切る順序によって別の平面で記録されるので、
/// > 平面対で引くと**別名で記録された頂点を取りこぼします。**
/// > 取りこぼすと T 字接合がそのまま残り、次数 1 の辺が出ます。**実際に踏みました。**
namespace detail {

/// スレッド局所に貯めた T 解決の統計を集約する（`SPEC-phase4.md` §1.1）。
///
/// **和と最大を取り違えないこと。** `max_per_edge` を足すと、スレッド数に比例して
/// 増える値になります。
inline void merge_tjunction_stats(TJunctionStats& a, const TJunctionStats& b) {
    a.inserted += b.inserted;
    a.candidates += b.candidates;
    a.degenerate_kept += b.degenerate_kept;
    a.apex_fallback += b.apex_fallback;
    a.edges_scanned += b.edges_scanned;
    a.inserted_from_cache += b.inserted_from_cache;
    a.max_per_edge = std::max(a.max_per_edge, b.max_per_edge);
}

}  // namespace detail

class PlaneVertexIndex {
public:
    /// `planes` の各平面について、その上に載る頂点を集める。
    ///
    /// 計算量は（平面数 × 頂点数）回の `side` です。平面ごとに 1 度だけなので、
    /// 辺ごとに全頂点を走査するより桁で軽くなります。
    /// **平面ごとに独立なので並列にできます**（`SPEC-phase4.md` §3）。
    ///
    /// **結果は平面 ID をキーにした表なので、順序に依存しません**（§4.2）。
    ///
    /// > **計算量は $O(\text{平面数} \times \text{頂点数})$ のままです。** 規模の
    /// > あるコーパスではここが出口の 88% を占めます（`BENCH.md`）。
    /// > **並列化は定数倍しか下げません。** 空間索引で $O(V \log V)$ にするのは
    /// > Phase 5 の課題です。
    void build(const PlaneTable& table, const std::vector<geom::HPointD>& verts,
               const std::vector<PlaneId>& planes, par::ThreadPool* pool = nullptr) {
        std::vector<PlaneId> uniq;
        uniq.reserve(planes.size());
        for (PlaneId p : planes) {
            if (map_.find(p) == map_.end()) uniq.push_back(p);
        }
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());

        std::vector<std::vector<std::uint32_t>> slot(uniq.size());
        const auto work = [&](std::size_t k, unsigned) {
            std::vector<std::uint32_t>& v = slot[k];
            for (std::uint32_t i = 0; i < verts.size(); ++i) {
                if (geom::side(table.at(uniq[k]), verts[i]) == 0) v.push_back(i);
            }
        };
        if (pool != nullptr) {
            pool->run(uniq.size(), work);
        } else {
            for (std::size_t k = 0; k < uniq.size(); ++k) work(k, 0u);
        }
        for (std::size_t k = 0; k < uniq.size(); ++k) map_[uniq[k]] = std::move(slot[k]);
    }

    /// 平面 `p` の上にある頂点。登録していなければ `nullptr`。
    const std::vector<std::uint32_t>* find(PlaneId p) const {
        auto it = map_.find(p);
        return (it == map_.end()) ? nullptr : &it->second;
    }

    std::size_t size() const noexcept { return map_.size(); }

private:
    std::map<PlaneId, std::vector<std::uint32_t>> map_;
};

/// **セルで区切った平面→頂点の索引**（`DESIGN-phase5-hotspots.md` §6.3 の A-3）。
///
/// `PlaneVertexIndex` は $O(\text{平面} \times \text{全頂点})$ で、**実測で出口の 82%、
/// 全体の 75% を占めていました**（`69268`、深度 8 で $4.13 \times 10^{10}$ 回の `side`）。
///
/// **絞り込みの根拠は幾何の事実です。**
///
///   多角形の辺は、その多角形の中にある。
///   多角形は、八分木の 1 つの葉の【閉じた箱】の中にある（幾何がセル境界でクリップされる）。
///   ⟹ 辺の内部に載る頂点は、その葉の閉じた箱の中にある。
///
/// **平面 ID の別名問題の影響を受けません。**「交線を平面 ID の対で索引する」
/// （相異なる 2 平面が同一の交線を与え得る）や「頂点の平面 3 つ組で絞る」
/// （4 枚以上が 1 点で交わる）は**漏れる**ので採っていません（§6.2）。
///
/// **費用**: $\sum_\ell |\text{支持平面}_\ell| \times |V_\ell|$。
/// 実測の削減は **3〜391 倍**で、**出力が占めるセルの数**で決まります（`BENCH.md`）。
///
/// **閉じた箱**であることが要点です。**面にちょうど載る頂点は、隣り合う葉にも属します。**
/// 半開区間で入れると、共有面の上の T 頂点を粗い側で見落とします。
class CellPlaneVertexIndex {
public:
    /// `box[i]` は多角形 `i` が属する葉の箱、`support[i]` はその支持平面。
    ///
    /// **葉は八分木のセルなので、深度と添字で一意に決まります。**
    /// 頂点は「最大深度の格子でどのセルに入るか」を二分探索で求め、
    /// **境界にちょうど載っていれば両側**を候補にします。
    /// **`false` を返したら、この索引は使えません。** 呼び出し側は従来の
    /// `PlaneVertexIndex` に退避してください。
    ///
    /// **`box` が八分木のセルの箱でないときに起こります。** `from_mesh` 直後のスープでは
    /// `Poly::aabb` は**三角形の外接箱**で、セルの箱ではありません
    /// （`boolean` を通ると分類の段でセルの箱に置き換わります）。
    ///
    /// > **退避したことを `groups() == 0` で観測できます。**
    /// > **スープ経路で 0 なら、機構が空回りしています**（`CLAUDE.md`）。
    bool build(const PlaneTable& table, const std::vector<geom::HPointD>& verts,
               const std::vector<octree::Aabb>& box, const std::vector<PlaneId>& support,
               par::ThreadPool* pool = nullptr, std::size_t* locate_tests = nullptr,
               std::size_t* group_tests = nullptr) {
        slot_.assign(box.size(), kNoGroup);
        group_.clear();
        if (box.empty() || verts.empty()) return false;

        // ---- 1. 箱 → 葉（深度と添字）。**箱の辺の長さから深度が決まる** ----
        //
        // **セルの箱でなければ諦めます。** 立方体でない、辺が 2 冪でない、
        // 格子に載っていない、のいずれかで判定できます。
        const auto depth_of = [](std::int64_t side_len) -> unsigned {
            for (unsigned d = 0; d + 1 <= kCoordBits; ++d) {
                if ((std::int64_t{1} << (kCoordBits - d)) == side_len) return d;
            }
            return kNoDepth;
        };
        std::unordered_map<std::uint64_t, std::uint32_t> leaf_id;
        std::vector<std::uint32_t> poly_leaf(box.size());
        unsigned dmax = 0;
        for (std::size_t i = 0; i < box.size(); ++i) {
            const std::int64_t sx = box[i].hi[0] - box[i].lo[0];
            if (sx != box[i].hi[1] - box[i].lo[1] || sx != box[i].hi[2] - box[i].lo[2])
                return false;
            const unsigned d = depth_of(sx);
            if (d == kNoDepth) return false;
            const std::int64_t step = std::int64_t{1} << (kCoordBits - d);
            for (int t = 0; t < 3; ++t) {
                if (((box[i].lo[t] - kCoordMin) % step) != 0) return false;
            }
            const std::uint32_t ix =
                static_cast<std::uint32_t>((box[i].lo[0] - kCoordMin) >> (kCoordBits - d));
            const std::uint32_t iy =
                static_cast<std::uint32_t>((box[i].lo[1] - kCoordMin) >> (kCoordBits - d));
            const std::uint32_t iz =
                static_cast<std::uint32_t>((box[i].lo[2] - kCoordMin) >> (kCoordBits - d));
            dmax = std::max(dmax, d);
            const std::uint64_t k = leaf_key(d, ix, iy, iz);
            auto it = leaf_id.find(k);
            if (it == leaf_id.end()) {
                it = leaf_id.emplace(k, static_cast<std::uint32_t>(leaf_id.size())).first;
            }
            poly_leaf[i] = it->second;
        }

        // ---- 2. 頂点を葉に配る。**閉じた箱**なので、面に載る頂点は複数の葉に入る ----
        std::vector<std::vector<std::uint32_t>> bucket(leaf_id.size());
        for (std::uint32_t v = 0; v < verts.size(); ++v) {
            std::uint32_t lo3[3], hi3[3];
            for (int ax = 0; ax < 3; ++ax) {
                const geom::Axis A = (ax == 0)   ? geom::Axis::X
                                     : (ax == 1) ? geom::Axis::Y
                                                 : geom::Axis::Z;
                // 最大深度の格子で、cell_bound(dmax, m) <= v となる最大の m
                std::uint32_t lo = 0, hi = 1u << dmax;
                while (lo < hi) {
                    const std::uint32_t mid = lo + (hi - lo + 1) / 2;
                    if (locate_tests != nullptr) ++*locate_tests;
                    if (geom::side(geom::plane_axis_aligned(A, octree::cell_bound(dmax, mid)),
                                   verts[v]) >= 0) {
                        lo = mid;
                    } else {
                        hi = mid - 1;
                    }
                }
                const std::uint32_t m = (lo >= (1u << dmax)) ? (1u << dmax) - 1 : lo;
                if (locate_tests != nullptr) ++*locate_tests;
                const bool on_line =
                    geom::side(geom::plane_axis_aligned(A, octree::cell_bound(dmax, m)),
                               verts[v]) == 0;
                lo3[ax] = (on_line && m > 0) ? m - 1 : m;
                hi3[ax] = m;
            }
            for (std::uint32_t i = lo3[0]; i <= hi3[0]; ++i)
                for (std::uint32_t j = lo3[1]; j <= hi3[1]; ++j)
                    for (std::uint32_t k = lo3[2]; k <= hi3[2]; ++k) {
                        // 最大深度のセルから、それを含む葉へ上る（葉は空間を分割する）
                        for (unsigned d = 0; d <= dmax; ++d) {
                            const unsigned sh = dmax - d;
                            auto it = leaf_id.find(leaf_key(d, i >> sh, j >> sh, k >> sh));
                            if (it == leaf_id.end()) continue;
                            std::vector<std::uint32_t>& b = bucket[it->second];
                            if (b.empty() || b.back() != v) b.push_back(v);
                            break;
                        }
                    }
        }

        // ---- 3. (葉, 支持平面) の組ごとに、載っている頂点を集める ----
        std::unordered_map<std::uint64_t, std::uint32_t> gid;
        std::vector<std::pair<std::uint32_t, PlaneId>> gkey;
        for (std::size_t i = 0; i < box.size(); ++i) {
            const std::uint64_t k = (static_cast<std::uint64_t>(poly_leaf[i]) << 32) |
                                    static_cast<std::uint32_t>(support[i]);
            auto it = gid.find(k);
            if (it == gid.end()) {
                it = gid.emplace(k, static_cast<std::uint32_t>(gkey.size())).first;
                gkey.emplace_back(poly_leaf[i], support[i]);
            }
            slot_[i] = it->second;
        }
        group_.assign(gkey.size(), {});
        // 計測用のスロット（`group_tests` が渡されたときだけ確保）
        std::vector<std::size_t> gtests;
        if (group_tests != nullptr) gtests.assign(gkey.size(), 0);
        const auto work = [&](std::size_t g, unsigned) {
            const geom::PlaneD& pl = table.at(gkey[g].second);
            std::vector<std::uint32_t>& out = group_[g];
            for (std::uint32_t v : bucket[gkey[g].first]) {
                if (geom::side(pl, verts[v]) == 0) out.push_back(v);
            }
            // **計測は群ごとのスロットに書きます**（並列区間なので原子操作を避ける。
            // `__atomic_*` は GCC / Clang の組み込みで、**MSVC にありません**）
            if (!gtests.empty()) gtests[g] = bucket[gkey[g].first].size();
        };
        if (pool != nullptr) {
            pool->run(gkey.size(), work);
        } else {
            for (std::size_t g = 0; g < gkey.size(); ++g) work(g, 0u);
        }
        if (group_tests != nullptr) {
            for (std::size_t v : gtests) *group_tests += v;
        }
        return true;
    }

    /// 多角形 `pi` の候補（その葉の箱に入り、その支持平面に載る頂点）。
    const std::vector<std::uint32_t>& candidates(std::size_t pi) const {
        static const std::vector<std::uint32_t> kEmpty;
        return (slot_[pi] == kNoGroup) ? kEmpty : group_[slot_[pi]];
    }

    /// **(葉, 支持平面) の組の数**。0 なら機構が空回りしています。
    std::size_t groups() const noexcept { return group_.size(); }

private:
    static constexpr std::uint32_t kNoGroup = 0xFFFFFFFFu;
    static constexpr unsigned kNoDepth = 0xFFFFFFFFu;
    static std::uint64_t leaf_key(unsigned d, std::uint32_t i, std::uint32_t j,
                                  std::uint32_t k) noexcept {
        return (static_cast<std::uint64_t>(d) << 60) | (static_cast<std::uint64_t>(i) << 40) |
               (static_cast<std::uint64_t>(j) << 20) | static_cast<std::uint64_t>(k);
    }
    std::vector<std::uint32_t> slot_;
    std::vector<std::vector<std::uint32_t>> group_;
};

namespace detail {

/// `a` と `b` が異なる最初の軸と、その向き。すべて同じなら `false`。
inline bool differing_axis(const geom::HPointD& a, const geom::HPointD& b, geom::Axis* out,
                           int* dir) noexcept {
    for (geom::Axis ax : {geom::Axis::X, geom::Axis::Y, geom::Axis::Z}) {
        const int c = geom::cmp_h(a, b, ax);
        if (c != 0) {
            *out = ax;
            *dir = c;
            return true;
        }
    }
    return false;
}

/// **同一直線上にある** `v` が線分 `a`–`b` の内部（両端を含まない）にあるか。
///
/// 共線であることは呼び出し側が保証します（索引が平面対で引くので構造的に成り立つ）。
/// したがって 1 軸の比較だけで決まり、**新しい述語は要りません**（§2.4.3）。
inline bool strictly_between(const geom::HPointD& a, const geom::HPointD& b,
                             const geom::HPointD& v) noexcept {
    geom::Axis ax{};
    int dir = 0;
    if (!differing_axis(a, b, &ax, &dir)) return false;  // a == b（退化した辺）
    return geom::cmp_h(a, v, ax) == dir && geom::cmp_h(v, b, ax) == dir;
}

}  // namespace detail

/// T 頂点を入れたあとの多角形。
///
/// `orig` の意味は `is_corner` で変わります。**扇分割の退化判定に使うのはこれだけで、
/// 幾何を一切見ません**（`fan_triangulate` 参照）。
struct TPolygon {
    std::vector<std::uint32_t> vertex;  ///< 頂点 ID（挿入後の巡回順）
    std::vector<char> is_corner;        ///< 元からあった角か
    /// 角なら**元の頂点添字** i、T 頂点なら**載っている元の辺の添字** j。
    std::vector<std::uint32_t> orig;
    std::uint32_t corners = 0;  ///< 元の角の数
};

/// 多角形の各辺に、その線分の**内部に載る大域頂点をすべて**挿入する（§2.4.3）。
///
/// - `poly[j]` から `poly[j+1]` への辺は平面対 `(support, edge[j])` の交線上にあります
///   （頂点 j = `support ∩ edge[j-1] ∩ edge[j]` なので、辺は `support ∩ edge[j]`）
///
/// **一律に適用してください。** 「隣が持っているから入れる」ではなく
/// **「線分の内部に載る大域頂点はすべて入れる」**という形にすること。片側だけに入れると
/// **T 字接合を直すどころか作ります**（辺 (a,b) が片側で (a,m),(m,b) になり、次数 1 の辺が
/// 3 本できる）。一律なら、同じ辺を共有する多角形は必ず同じ頂点を得ます。
///
/// **1 個ずつ入れてはいけません**（§2.4.3）。生成された部分辺にさらに候補が載る場合を
/// 取りこぼします。ここでは元の線分について候補を**全部集めてから**並べて入れます。
inline TPolygon insert_t_vertices_with(const PlaneTable& table,
                                       const std::vector<geom::HPointD>& verts,
                                       const std::vector<std::uint32_t>& cand_in,
                                       const std::vector<PlaneId>& edge,
                                       const std::vector<std::uint32_t>& poly,
                                       TJunctionStats* stats = nullptr,
                                       const std::vector<char>* from_cache = nullptr) {
    KRISITE_CHECK(poly.size() == edge.size(), "insert_t_vertices: 頂点数と辺数が違う");
    const std::size_t n = poly.size();

    TPolygon out;
    out.corners = static_cast<std::uint32_t>(n);
    out.vertex.reserve(n);
    out.is_corner.reserve(n);
    out.orig.reserve(n);
    if (n < 3) {
        out.vertex = poly;
        out.is_corner.assign(n, 1);
        for (std::uint32_t i = 0; i < n; ++i) out.orig.push_back(i);
        return out;
    }

#if defined(KRISITE_MUTATION_NO_TJUNCTION)
    // SPEC-phase2 §9.3 の変異 10: T 頂点の解決を無効化する。
    //
    // **T 解決は継ぎ目の正しさの単一障害点です**（§2.4.4 (3)）。粗い側に頂点が入らないと
    // 共有面で辺の分かれ方が食い違い、次数 1 の辺が出ます。**ケース 13 で検出されること。**
    out.vertex = poly;
    out.is_corner.assign(n, 1);
    for (std::uint32_t i = 0; i < n; ++i) out.orig.push_back(i);
    (void)table;
    (void)verts;
    (void)cand_in;
    (void)edge;
    (void)stats;
    (void)from_cache;
    return out;
#endif

    for (std::size_t j = 0; j < n; ++j) {
        const std::uint32_t a = poly[j];
        const std::uint32_t b = poly[(j + 1) % n];
        out.vertex.push_back(a);
        out.is_corner.push_back(1);
        out.orig.push_back(static_cast<std::uint32_t>(j));

        if (stats) ++stats->edges_scanned;
        const std::vector<std::uint32_t>* cand = &cand_in;
        if (cand->empty()) continue;
        if (stats) stats->candidates += cand->size();

        // 手順 2: **両方の平面に載っていて**、区間の内部にあるものを【全部】集める。
        //
        // 支持平面の上にあることは索引が保証済み。もう 1 枚を `side` で見れば交線上に
        // あることが確定し、`strictly_between` の前提（共線）が満たされます。
        const geom::PlaneD& qp = table.at(edge[j]);
        std::vector<std::uint32_t> on_edge;
        for (std::uint32_t v : *cand) {
            if (v == a || v == b) continue;
            if (geom::side(qp, verts[v]) != 0) continue;
            if (detail::strictly_between(verts[a], verts[b], verts[v])) on_edge.push_back(v);
        }
        if (on_edge.empty()) continue;

#if defined(KRISITE_MUTATION_TJUNCTION_ONE_AT_A_TIME)
        // SPEC-phase2 §9.3 の変異 11: 1 本の辺に 1 個しか入れない。
        // **1 本の辺に T 頂点が 2 個以上載る配置でしか検出できません**（§8 のケース 13 の要件）。
        on_edge.resize(1);
#endif

        // 手順 3: 線分に沿って整列する。共線なので 1 軸の比較で全順序が決まる
        geom::Axis ax{};
        int dir = 0;
        // **`KRISITE_CHECK` の中で呼ばないこと。** 検査 OFF のビルドで呼び出しごと消え、
        // `ax` / `dir` が未初期化のまま使われます。
        const bool ok = detail::differing_axis(verts[a], verts[b], &ax, &dir);
        KRISITE_CHECK(ok, "insert_t_vertices: 退化した辺に候補が載っている");
        (void)ok;
        std::sort(on_edge.begin(), on_edge.end(), [&](std::uint32_t x, std::uint32_t y) {
            return geom::cmp_h(verts[x], verts[y], ax) == dir;
        });

        for (std::uint32_t v : on_edge) {
            out.vertex.push_back(v);
            out.is_corner.push_back(0);
            out.orig.push_back(static_cast<std::uint32_t>(j));
            // §13 の CP5:「構成点の保持 × T 解決」を通った回数
            if (stats != nullptr && from_cache != nullptr && v < from_cache->size() &&
                (*from_cache)[v] != 0) {
                ++stats->inserted_from_cache;
            }
        }
        if (stats) {
            stats->inserted += on_edge.size();
            stats->max_per_edge = std::max(stats->max_per_edge, on_edge.size());
        }
    }
    return out;
}

/// **旧署名の互換ラッパ**（`PlaneVertexIndex` を引く形）。
///
/// **二項メッシュ経路（`boolean_op`）と `test_tjunction.cpp` が使います。**
/// **二項経路は意図的に素朴なまま**にしています — 正解器は被検体と別経路で書く
/// （`IMPL-phase5.md` §12）。スープ経路の A-3（`CellPlaneVertexIndex`）は
/// **ここを通りません。**
inline TPolygon insert_t_vertices(const PlaneTable& table, const std::vector<geom::HPointD>& verts,
                                  const PlaneVertexIndex& index, PlaneId support,
                                  const std::vector<PlaneId>& edge,
                                  const std::vector<std::uint32_t>& poly,
                                  TJunctionStats* stats = nullptr,
                                  const std::vector<char>* from_cache = nullptr) {
    static const std::vector<std::uint32_t> kEmpty;
    const std::vector<std::uint32_t>* c = index.find(support);
    return insert_t_vertices_with(table, verts, (c == nullptr) ? kEmpty : *c, edge, poly, stats,
                                  from_cache);
}

namespace detail {

/// 出力頂点 `k` が、元の辺 `line` の上に載っているか。**組合せだけで決まります。**
///
/// 元の辺 `L` は元の角 `L` から角 `L+1` へ向かうので、
///   - T 頂点なら「載っている辺の添字が L か」
///   - 角なら「その添字が L か L+1 か」
inline bool on_original_edge(const TPolygon& p, std::size_t k, std::uint32_t line) noexcept {
    const std::uint32_t n = p.corners;
    if (!p.is_corner[k]) return p.orig[k] == line;
    return p.orig[k] == line || p.orig[k] == (line + 1) % n;
}

}  // namespace detail

/// T 頂点を含む凸多角形を扇分割する（§2.4.4 (2)）。
///
/// > **退化三角形は捨てるのではなく、作らないこと。**
///
/// **捨ててはいけません。** 凸多角形の扇分割は $n-2$ 枚が揃って初めて円板になります
/// （$\chi = n - (2n-3) + (n-2) = 1$）。1 枚捨てると、その三角形が持っていた**境界辺が
/// メッシュから消え**、隣の多角形はその辺を持ったままなので**次数 1 の辺が残ります**。
/// **面積は保たれても、組合せ的には保たれません。**
/// 実測: 捨てると §9.1 が 138/574 で落ち、捨てないと 48/574 に減りました。
///
/// **起点を「両隣の元の辺に T 頂点が無い元の角」に取れば、退化は最初から出ません。**
/// 起点が乗る直線はそれに接する 2 本の元の辺だけなので、両隣に T 頂点が無ければ
/// 共線な三つ組は生じません。**そのような角が無い多角形では残します。**
/// 面積 0 でも組合せ的に必要です。枚数は `degenerate_kept` に数えます（§11）。
///
/// **退化の判定に幾何は要りません。** 起点は角なので、起点と共線になり得るのは
/// 起点に接する 2 本の元の辺の上の頂点だけです。両方が同じ辺の上にあるかを見れば足ります。
/// 新しい述語を足さずに済むのが要点です（§2.4.3「新しい述語は要りません」）。
inline void fan_triangulate(const TPolygon& p, std::vector<std::array<std::uint32_t, 3>>& out,
                            TJunctionStats* stats = nullptr) {
    const std::size_t n = p.vertex.size();
    if (n < 3) return;

    // 元の辺ごとに T 頂点が載っているか
    std::vector<char> edge_has_t(p.corners, 0);
    for (std::size_t i = 0; i < n; ++i) {
        if (!p.is_corner[i]) edge_has_t[p.orig[i]] = 1;
    }

    // **起点は「両隣の元の辺に T 頂点が無い元の角」**（§2.4.4 (2)）
    std::size_t apex = n;
    for (std::size_t i = 0; i < n; ++i) {
        if (!p.is_corner[i]) continue;
#if !defined(KRISITE_MUTATION_DROP_DEGENERATE)
        const std::uint32_t oc = p.orig[i];
        const std::uint32_t prev = (oc + p.corners - 1) % p.corners;
        if (edge_has_t[oc] || edge_has_t[prev]) continue;
#endif
        apex = i;
        break;
    }
    if (apex == n) {
        // 選べなかった。**残します**（捨てると境界辺が消えて次数 1 の辺が出る）
        if (stats) ++stats->apex_fallback;
        for (std::size_t i = 0; i < n; ++i) {
            if (p.is_corner[i]) {
                apex = i;
                break;
            }
        }
    }
    KRISITE_CHECK(apex < n && p.is_corner[apex], "fan_triangulate: 元の角が 1 つも無い");

    // 起点に接する元の辺は 2 本。起点の元の頂点添字を c として、辺 c（出る側）と c-1（入る側）
    const std::uint32_t c = p.orig[apex];
    const std::uint32_t lines[2] = {c, (c + p.corners - 1) % p.corners};

    for (std::size_t k = 1; k + 1 < n; ++k) {
        const std::size_t i1 = (apex + k) % n;
        const std::size_t i2 = (apex + k + 1) % n;
        bool degenerate = false;
        for (std::uint32_t line : lines) {
            if (detail::on_original_edge(p, i1, line) && detail::on_original_edge(p, i2, line)) {
                degenerate = true;
                break;
            }
        }
        if (degenerate) {
            if (stats) ++stats->degenerate_kept;
#if defined(KRISITE_MUTATION_DROP_DEGENERATE)
            // SPEC-phase2 §9.3 の変異 13: 作らないのではなく**捨てる**。
            // 境界辺が消えて次数 1 の辺が残ります。**面積は保たれるので体積検査では
            // 捕まりません。**
            continue;
#endif
        }
        out.push_back({p.vertex[apex], p.vertex[i1], p.vertex[i2]});
    }
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_TJUNCTION_HPP

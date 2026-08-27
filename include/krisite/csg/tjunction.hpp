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
#include <utility>
#include <vector>

#include "krisite/csg/plane_table.hpp"
#include "krisite/geom/point.hpp"
#include "krisite/geom/predicates.hpp"

namespace krisite::csg {

/// §2.4.3 / §11 の記録。
struct TJunctionStats {
    std::size_t inserted = 0;            ///< 挿入された T 頂点の数
    std::size_t candidates = 0;          ///< 索引が返した候補の数（絞り込みの効き）
    std::size_t max_per_edge = 0;        ///< 1 本の辺に入った T 頂点の最大数
    std::size_t degenerate_dropped = 0;  ///< 扇分割で捨てた面積 0 の三角形
    std::size_t edges_scanned = 0;       ///< 走査した辺の数（候補数の分母）
};

/// **辺の平面対 → その交線上にある大域頂点 ID** の索引（§2.4.3 の手順 1）。
///
/// 多角形の辺 (a,b) は 2 平面 (P,Q) の交線上にあります。その線上に載り得る頂点は、
/// **平面3つ組が {P,Q} を含むもの**だけです。3つ組を登録すれば 3 つの対が張れます。
///
/// **これが「候補を絞る」ということの全てです。** 線上に載ることは 3つ組から構造的に
/// 決まるので、幾何的な判定は §2.4.3 の手順 2（区間の内部か）だけで済みます。
class EdgeVertexIndex {
public:
    /// 頂点 `vid` の生成 3 つ組を登録する。3 つの対すべてに入ります。
    void add_triple(const std::array<PlaneId, 3>& tri, std::uint32_t vid) {
        add_pair(tri[0], tri[1], vid);
        add_pair(tri[0], tri[2], vid);
        add_pair(tri[1], tri[2], vid);
    }

    /// 平面対 (p,q) の交線上にある頂点。無ければ `nullptr`。
    const std::vector<std::uint32_t>* find(PlaneId p, PlaneId q) const {
        auto it = map_.find(key(p, q));
        return (it == map_.end()) ? nullptr : &it->second;
    }

    std::size_t size() const noexcept { return map_.size(); }

private:
    static std::pair<PlaneId, PlaneId> key(PlaneId p, PlaneId q) noexcept {
        return (p < q) ? std::pair<PlaneId, PlaneId>{p, q} : std::pair<PlaneId, PlaneId>{q, p};
    }

    void add_pair(PlaneId p, PlaneId q, std::uint32_t vid) {
        if (p == q) return;  // 同じ平面が 2 度現れる 3 つ組は線を定めない
        std::vector<std::uint32_t>& v = map_[key(p, q)];
        // 第2段の値併合で複数の 3 つ組が同じ ID に落ちるので、重複を避ける
        if (std::find(v.begin(), v.end(), vid) == v.end()) v.push_back(vid);
    }

    std::map<std::pair<PlaneId, PlaneId>, std::vector<std::uint32_t>> map_;
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
inline TPolygon insert_t_vertices(const std::vector<geom::HPointD>& verts,
                                  const EdgeVertexIndex& index, PlaneId support,
                                  const std::vector<PlaneId>& edge,
                                  const std::vector<std::uint32_t>& poly,
                                  TJunctionStats* stats = nullptr) {
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
    (void)verts;
    (void)index;
    (void)support;
    (void)edge;
    (void)stats;
    return out;
#endif

    for (std::size_t j = 0; j < n; ++j) {
        const std::uint32_t a = poly[j];
        const std::uint32_t b = poly[(j + 1) % n];
        out.vertex.push_back(a);
        out.is_corner.push_back(1);
        out.orig.push_back(static_cast<std::uint32_t>(j));

        if (stats) ++stats->edges_scanned;
        const std::vector<std::uint32_t>* cand = index.find(support, edge[j]);
        if (cand == nullptr) continue;
        if (stats) stats->candidates += cand->size();

        // 手順 2: 区間の内部にあるものを【全部】集める
        std::vector<std::uint32_t> on_edge;
        for (std::uint32_t v : *cand) {
            if (v == a || v == b) continue;
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
        }
        if (stats) {
            stats->inserted += on_edge.size();
            stats->max_per_edge = std::max(stats->max_per_edge, on_edge.size());
        }
    }
    return out;
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
/// **起点は T 頂点でない頂点（元の角）に取ります。** T 頂点を起点にすると共線な三つ組が
/// できるためです。それでも起点の隣接辺に T 頂点があると退化三角形が出ますが、
/// **凸多角形なのでそれらは何も覆いません。捨てても被覆が保たれます。**
///
/// **退化の判定に幾何は要りません。** 起点は角なので、起点と共線になり得るのは
/// 起点に接する 2 本の元の辺の上の頂点だけです。両方が同じ辺の上にあるかを見れば足ります。
/// 新しい述語を足さずに済むのが要点です（§2.4.3「新しい述語は要りません」）。
inline void fan_triangulate(const TPolygon& p, std::vector<std::array<std::uint32_t, 3>>& out,
                            TJunctionStats* stats = nullptr) {
    const std::size_t n = p.vertex.size();
    if (n < 3) return;

    std::size_t apex = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (p.is_corner[i]) {
            apex = i;
            break;
        }
    }
    KRISITE_CHECK(p.is_corner[apex], "fan_triangulate: 元の角が 1 つも無い");

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
            if (stats) ++stats->degenerate_dropped;
            continue;
        }
        out.push_back({p.vertex[apex], p.vertex[i1], p.vertex[i2]});
    }
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_TJUNCTION_HPP

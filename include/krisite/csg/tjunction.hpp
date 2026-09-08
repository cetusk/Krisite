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
    /// **一般解**（§2.4.4 (2)）が作った三角形の数。**0 なら空回りです**
    std::size_t general_used = 0;
    /// 一般解が組めず従来の扇分割に落ちた多角形の数
    std::size_t general_fallback = 0;
    std::size_t edges_scanned = 0;  ///< 走査した辺の数（候補数の分母）
    /// **`side`（平面 × 同次点）を評価した回数**。**照合の本体の演算回数です。**
    ///
    /// > **`candidates` は「候補集合の大きさ × 辺の数」で、
    /// > 実際に述語を評価した回数ではありません**（頂点と一致する候補は
    /// > `side` の前に落ちます）。**効果は演算回数で測る**（`CLAUDE.md`）ので、
    /// > **述語の評価回数を直接数えます。**
    std::size_t side_tests = 0;
    /// **`strictly_between`（`cmp_h` を 2〜6 回）を評価した回数**。
    ///
    /// **`side` が 0 を返した候補だけが到達します。** こちらのほうが 1 回あたり高い
    /// （$13b{+}27$ 対 $9b{+}20$）ので、分けて数えます。
    std::size_t between_tests = 0;
    /// **候補集合を走査した回数**（多角形あたり 1 回なら `polys` と等しくなる）。
    ///
    /// **`edges_scanned` と比べてください。** 辺ごとに走査していれば辺の数と等しく、
    /// 多角形あたり 1 回なら多角形の数と等しくなります。
    std::size_t cand_scans = 0;
    /// **多角形の外接箱の【外】にあった候補の数**（`ToMeshOptions::count_box_reject`）。
    ///
    /// **案 (b)（多角形の箱で候補を落とす）の効果の上限を数えます**（`SPEC-phase5.md` §5.11）。
    /// **箱の外なら、その候補はどの辺の相対内部にも載りません。**
    ///
    /// > **★ これは【効果】であって【費用】ではありません。**
    /// > 落とせる `side` の回数は `box_reject_side_saved` が数えます。
    /// > **前判定そのものの費用は、実装の形を決めてから測ります**
    /// > （`CLAUDE.md`「前判定は、置き換える仕事より安くなければ意味がありません」）。
    ///
    /// **既定では数えません。** 箱を作るのに `cmp_h` を多角形あたり $3(n-1)$ 回、
    /// 判定に候補あたり最大 6 回払うので、**計測の費用は小さくありません。**
    std::size_t cand_outside_box = 0;
    /// **箱で落とせたはずの `side` の評価回数**（落とした候補 × その多角形の辺の数）。
    std::size_t box_reject_side_saved = 0;
    /// **箱の判定に払った `cmp_h` の回数**（前判定の費用の側）。
    std::size_t box_cmp_tests = 0;
    /// **箱の判定にかけた候補の数**（`cand_outside_box` の分母）。
    std::size_t box_cand_tested = 0;
    /// **箱を作るのに払った `cmp_h` の回数**（判定の費用と分けるため）。
    std::size_t box_build_cmp = 0;
    /// **二分探索で飛ばした候補の数**（案 (b2) が実際に落とした数）。
    ///
    /// **0 なら機構が空回りしています**（`CLAUDE.md`）。
    std::size_t cand_skipped_by_range = 0;
    /// **1 軸だけで落とせた候補の数**（軸ごと）。
    ///
    /// **3 軸すべてを見る必要があるかを決める材料です。** 1 軸なら `cmp_h` が 2 回で済みます。
    ///
    /// **そして「候補を軸で整列して区間で絞る」案の効果も、ここから読めます** —
    /// **1 軸の区間の外にある候補は、整列してあれば二分探索で飛ばせます。**
    std::size_t box_reject_axis[3] = {0, 0, 0};
    /// **T 頂点の挿入に費やした CPU 時間**（ミリ秒。`ToMeshOptions::time_stages` で有効）。
    ///
    /// > **★ これは CPU 時間です。壁時計ではありません**（`CLAUDE.md`「CPU 時間と
    /// > 壁時計を混ぜないでください」）。**スレッド局所に貯めて足しているので、
    /// > 合計と比べる前に並列効率で割ってください。**
    ///
    /// **既定では 0 のままです。** 多角形ごとに時計を 2 回読むので、
    /// **計測の費用を本番の経路に持ち込まない**ため旗で囲んでいます（`CLAUDE.md`）。
    double ms_insert_t = 0.0;
    /// **扇分割（`fan_triangulate`）に費やした CPU 時間**（同上）。
    double ms_fan_tri = 0.0;
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
    a.general_used += b.general_used;
    a.general_fallback += b.general_fallback;
    a.edges_scanned += b.edges_scanned;
    a.side_tests += b.side_tests;
    a.between_tests += b.between_tests;
    a.cand_scans += b.cand_scans;
    a.cand_outside_box += b.cand_outside_box;
    a.box_reject_side_saved += b.box_reject_side_saved;
    a.box_cmp_tests += b.box_cmp_tests;
    a.box_cand_tested += b.box_cand_tested;
    a.box_build_cmp += b.box_build_cmp;
    a.cand_skipped_by_range += b.cand_skipped_by_range;
    for (int t = 0; t < 3; ++t) a.box_reject_axis[t] += b.box_reject_axis[t];
    a.inserted_from_cache += b.inserted_from_cache;
    a.ms_insert_t += b.ms_insert_t;
    a.ms_fan_tri += b.ms_fan_tri;
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
               const std::vector<std::uint8_t>& cell_depth, par::ThreadPool* pool = nullptr,
               std::size_t* locate_tests = nullptr, std::size_t* group_tests = nullptr,
               bool sort_candidates = true) {
        slot_.assign(box.size(), kNoGroup);
        group_.clear();
        if (box.empty() || verts.empty()) return false;
        if (cell_depth.size() != box.size()) return false;

        // ---- 1. 箱 → 葉（深度と添字）----
        //
        // **深度は呼び出し側が持っています**（`Poly::cell_depth`。§5.10.6）。
        //
        // > **以前は「箱の辺の長さ」から逆算していました。**
        // > **`Poly::aabb` を「元の多角形の箱 ∩ セル箱」に狭めたので、逆算できません。**
        // > **箱に「割り当ての範囲」と「葉の符号」の 2 つの役目があったのが誤りでした。**
        //
        // **添字は `box.lo` から復元します。** 箱は葉の箱に含まれ、
        // `box.lo` は $[\text{cell.lo}, \text{cell.hi})$ に居るので、床関数で一意です。
        std::unordered_map<std::uint64_t, std::uint32_t> leaf_id;
        std::vector<std::uint32_t> poly_leaf(box.size());
        unsigned dmax = 0;
        for (std::size_t i = 0; i < box.size(); ++i) {
            if (cell_depth[i] == octree::kNoCellDepth)
                return false;  // 深度が不明（`from_mesh` 直後）
            const unsigned d = cell_depth[i];
            if (d + 1 > kCoordBits) return false;
            const std::uint32_t ix =
                static_cast<std::uint32_t>((box[i].lo[0] - kCoordMin) >> (kCoordBits - d));
            const std::uint32_t iy =
                static_cast<std::uint32_t>((box[i].lo[1] - kCoordMin) >> (kCoordBits - d));
            const std::uint32_t iz =
                static_cast<std::uint32_t>((box[i].lo[2] - kCoordMin) >> (kCoordBits - d));
            // **箱が本当にその葉に収まっているか**を確かめます。
            // **収まっていなければ、深度の記録と箱が食い違っています**（退避）。
            for (int t = 0; t < 3; ++t) {
                const std::uint32_t m = (t == 0) ? ix : ((t == 1) ? iy : iz);
                if (box[i].lo[t] < octree::cell_bound(d, m) ||
                    box[i].hi[t] > octree::cell_bound(d, m + 1)) {
                    return false;
                }
            }
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
            // ---- ★ 候補を X 軸で整列する（`DESIGN-phase5-hotspots.md` §13、案 (b2)）----
            //
            // **T 解決の照合は、多角形の X 区間の外にある候補を全部飛ばせます。**
            // **整列してあれば二分探索で飛ばせるので、群ごとに 1 回整列すれば、
            // その群の全多角形で使い回せます**（群あたり 5.6〜13.7 多角形）。
            //
            // **軸は X に固定します**（`SPEC-phase5.md` §5.11）。
            // **3 軸の削減率がほぼ同じ**（83.6 / 83.9 / 83.9、95.8 / 95.2 / 95.8）で、
            // **実行時に選ぶと非決定的になります**（削減率は測らないと分からない）。
            //
            // **同じ X を持つ候補は頂点 ID で決めます。** `cmp_h` は同値に順序を
            // 付けないので、**決めないと入力の並び次第で順序が変わります**
            // （`SPEC-phase4.md` §4.4 と同じ形）。
            if (sort_candidates) {
                std::sort(out.begin(), out.end(), [&verts](std::uint32_t a, std::uint32_t b) {
                    const int c = geom::cmp_h(verts[a], verts[b], geom::Axis::X);
                    return (c != 0) ? (c < 0) : (a < b);
                });
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
        sorted_ = sort_candidates;
        return true;
    }

    /// **候補が X 軸で整列されているか**（案 (b2)）。偽なら二分探索を使えません。
    bool sorted() const noexcept { return sorted_; }

    /// 多角形 `pi` の候補（その葉の箱に入り、その支持平面に載る頂点）。
    const std::vector<std::uint32_t>& candidates(std::size_t pi) const {
        static const std::vector<std::uint32_t> kEmpty;
        return (slot_[pi] == kNoGroup) ? kEmpty : group_[slot_[pi]];
    }

    /// **(葉, 支持平面) の組の数**。0 なら機構が空回りしています。
    std::size_t groups() const noexcept { return group_.size(); }

private:
    bool sorted_ = false;
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
inline TPolygon insert_t_vertices_with(
    const PlaneTable& table, const std::vector<geom::HPointD>& verts,
    const std::vector<std::uint32_t>& cand_in, const std::vector<PlaneId>& edge,
    const std::vector<std::uint32_t>& poly, TJunctionStats* stats = nullptr,
    const std::vector<char>* from_cache = nullptr, bool scan_per_edge = false,
    bool count_box_reject = false, bool sorted_cand = false) {
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

    // ---- 手順 2: 候補集合を【多角形あたり 1 回】走査する（`SPEC-phase5.md` §5.11）----
    //
    // **元は「辺ごとに候補集合を全部走査する」形でした。** 候補集合は
    // (葉, 支持平面) 群で共有されるので**多角形の中で同じ**なのに、
    // **$n$ 角形なら $n$ 回走査し直していました。**
    //
    // **判定は 1 つも変えていません。** 走査の順序を入れ替えただけです。
    // **候補の対（候補 × 辺）は同じ集合を回るので、出力はバイト単位で同一**です。
    //
    // > **候補 $v$ は、高々 1 本の辺の【相対内部】にしか載りません。**
    // > 多角形は（弱く）凸で、`strictly_between` は端点を含まないためです。
    // > **だから当たった時点で打ち切れます。** ただし**当たるのは 1 万分の 1 以下**なので、
    // > 打ち切りが `side` の回数に効く量は小さい（`IMPL-v2.md` §5）。
    //
    // **支持平面の上にあることは索引が保証済み。** もう 1 枚を `side` で見れば
    // 交線上にあることが確定し、`strictly_between` の前提（共線）が満たされます。
    std::vector<std::vector<std::uint32_t>> on_edge(n);

    // ---- 案 (b) の【効果】を数える（`SPEC-phase5.md` §5.11。**計測のみ**）----
    //
    // **多角形の外接箱の外にある候補は、どの辺の相対内部にも載りません。**
    // **落とせる `side` の回数を数えるだけで、実際には落としません**
    // （落とすと出力の同一性を別に示す必要が出ます。**いまは効果だけ知りたい**）。
    //
    // **既定では走りません。** 箱を作る `cmp_h` が多角形あたり $3(n-1)$ 回、
    // 判定が候補あたり最大 6 回。**計測の費用は小さくありません。**
    if (count_box_reject && stats != nullptr && !cand_in.empty()) {
        std::size_t lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
        const geom::Axis kAxis[3] = {geom::Axis::X, geom::Axis::Y, geom::Axis::Z};
        for (std::size_t j = 1; j < n; ++j) {
            for (int t = 0; t < 3; ++t) {
                const auto a3 = static_cast<std::size_t>(t);
                stats->box_cmp_tests += 2;
                stats->box_build_cmp += 2;
                if (geom::cmp_h(verts[poly[j]], verts[poly[lo[a3]]], kAxis[t]) < 0) lo[a3] = j;
                if (geom::cmp_h(verts[poly[j]], verts[poly[hi[a3]]], kAxis[t]) > 0) hi[a3] = j;
            }
        }
        stats->box_cand_tested += cand_in.size();
        for (std::uint32_t v : cand_in) {
            bool outside = false;
            // **軸ごとの内訳を採るので、打ち切りません**（計測なので費用は問いません）。
            // **実装するときは打ち切ります。**
            for (int t = 0; t < 3; ++t) {
                const auto a3 = static_cast<std::size_t>(t);
                stats->box_cmp_tests += 2;
                if (geom::cmp_h(verts[v], verts[poly[lo[a3]]], kAxis[t]) < 0 ||
                    geom::cmp_h(verts[v], verts[poly[hi[a3]]], kAxis[t]) > 0) {
                    outside = true;
                    ++stats->box_reject_axis[a3];
                }
            }
            if (outside) {
                ++stats->cand_outside_box;
                stats->box_reject_side_saved += n;
            }
        }
    }

    // ---- ★ 案 (b2): 多角形の X 区間の外にある候補を、二分探索で飛ばす ----------
    //
    // **候補は群ごとに X 軸で整列済みです**（`CellPlaneVertexIndex::build`）。
    // **辺は多角形の境界上にあるので、辺の区間は多角形の区間に含まれます。**
    // **多角形の区間の外にある候補は、どの辺の相対内部にも載りません。厳密です。**
    //
    // **費用は多角形あたり $2\log \lvert V_g \rvert$ 回の `cmp_h`**（区間の両端）。
    // **辺ごとに絞る形は採りません** — 区間内に残るのは 10〜16 個で、
    // **二分探索より線形走査のほうが安い**からです（`DESIGN` §13）。
    std::size_t c_begin = 0, c_end = cand_in.size();
    if (sorted_cand && !cand_in.empty()) {
        std::size_t lo_i = 0, hi_i = 0;
        for (std::size_t j = 1; j < n; ++j) {
            if (stats) stats->box_cmp_tests += 2;
            if (geom::cmp_h(verts[poly[j]], verts[poly[lo_i]], geom::Axis::X) < 0) lo_i = j;
            if (geom::cmp_h(verts[poly[j]], verts[poly[hi_i]], geom::Axis::X) > 0) hi_i = j;
        }
        // **閉区間 [lo, hi]。** 端と同じ X を持つ候補も残します
        const auto lower =
            std::lower_bound(cand_in.begin(), cand_in.end(), poly[lo_i],
                             [&verts, stats](std::uint32_t a, std::uint32_t key) {
                                 if (stats) ++stats->box_cmp_tests;
                                 return geom::cmp_h(verts[a], verts[key], geom::Axis::X) < 0;
                             });
        const auto upper =
            std::upper_bound(cand_in.begin(), cand_in.end(), poly[hi_i],
                             [&verts, stats](std::uint32_t key, std::uint32_t a) {
                                 if (stats) ++stats->box_cmp_tests;
                                 return geom::cmp_h(verts[key], verts[a], geom::Axis::X) < 0;
                             });
        c_begin = static_cast<std::size_t>(lower - cand_in.begin());
        c_end = static_cast<std::size_t>(upper - cand_in.begin());
        if (stats) stats->cand_skipped_by_range += cand_in.size() - (c_end - c_begin);
    }

    if (!cand_in.empty() && !scan_per_edge) {
        if (stats) ++stats->cand_scans;
        for (std::size_t ci = c_begin; ci < c_end; ++ci) {
            const std::uint32_t v = cand_in[ci];
            for (std::size_t j = 0; j < n; ++j) {
                const std::uint32_t a = poly[j];
                const std::uint32_t b = poly[(j + 1) % n];
                if (v == a || v == b) continue;
                if (stats) ++stats->side_tests;
                if (geom::side(table.at(edge[j]), verts[v]) != 0) continue;
                if (stats) ++stats->between_tests;
                if (detail::strictly_between(verts[a], verts[b], verts[v])) {
                    on_edge[j].push_back(v);
                    break;
                }
            }
        }
    } else if (!cand_in.empty()) {
        // **辺ごとに候補集合を走査する従来の形**（`scan_per_edge`。**比較の基準側**）。
        // **こちらは二分探索を使いません**（基準側なので素朴なまま）。
        //
        // **判定は上とまったく同じです。** 走査の順序だけが違います。
        // **`CLAUDE.md`「機構を追加したら、それを外す経路も用意してください」。**
        // **両方で出力がバイト一致することを `test_tjunction.cpp` が検査します。**
        for (std::size_t j = 0; j < n; ++j) {
            const std::uint32_t a = poly[j];
            const std::uint32_t b = poly[(j + 1) % n];
            if (stats) ++stats->cand_scans;
            const geom::PlaneD& qp = table.at(edge[j]);
            for (std::uint32_t v : cand_in) {
                if (v == a || v == b) continue;
                if (stats) ++stats->side_tests;
                if (geom::side(qp, verts[v]) != 0) continue;
                if (stats) ++stats->between_tests;
                if (detail::strictly_between(verts[a], verts[b], verts[v])) on_edge[j].push_back(v);
            }
        }
    }

    for (std::size_t j = 0; j < n; ++j) {
        const std::uint32_t a = poly[j];
        const std::uint32_t b = poly[(j + 1) % n];
        out.vertex.push_back(a);
        out.is_corner.push_back(1);
        out.orig.push_back(static_cast<std::uint32_t>(j));

        if (stats) ++stats->edges_scanned;
        if (cand_in.empty()) continue;
        // **`candidates` の定義は変えません**（辺ごとの候補数の和）。
        // **変えると過去の記録と比較できなくなります。**
        // **実際に走査した回数は `cand_scans` が別に数えます。**
        if (stats) stats->candidates += cand_in.size();
        if (on_edge[j].empty()) continue;

#if defined(KRISITE_MUTATION_TJUNCTION_ONE_AT_A_TIME)
        // SPEC-phase2 §9.3 の変異 11: 1 本の辺に 1 個しか入れない。
        // **1 本の辺に T 頂点が 2 個以上載る配置でしか検出できません**（§8 のケース 13 の要件）。
        on_edge[j].resize(1);
#endif

        // 手順 3: 線分に沿って整列する。共線なので 1 軸の比較で全順序が決まる
        geom::Axis ax{};
        int dir = 0;
        // **`KRISITE_CHECK` の中で呼ばないこと。** 検査 OFF のビルドで呼び出しごと消え、
        // `ax` / `dir` が未初期化のまま使われます。
        const bool ok = detail::differing_axis(verts[a], verts[b], &ax, &dir);
        KRISITE_CHECK(ok, "insert_t_vertices: 退化した辺に候補が載っている");
        (void)ok;
        std::sort(on_edge[j].begin(), on_edge[j].end(), [&](std::uint32_t x, std::uint32_t y) {
            return geom::cmp_h(verts[x], verts[y], ax) == dir;
        });

        for (std::uint32_t v : on_edge[j]) {
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
            stats->inserted += on_edge[j].size();
            stats->max_per_edge = std::max(stats->max_per_edge, on_edge[j].size());
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
/// **一般解**（`SPEC-phase2.md` §2.4.4 (2) の「一般解」。2026-09-07 に承認・実装）。
///
/// > 元の角だけの多角形を先に三角形化し、T 頂点を持つ元の辺ごとに、
/// > その辺を含む三角形の対頂点へ扇を張る。**退化は 1 枚も出ません。**
///
/// **実装の条件が満たされたので入れました。** 仕様は
/// 「**残した枚数が実際に効いてから**判断してください」としており、
/// **radial sort が解けない直接の原因になりました**（`IMPL-phase5.md` §94 / §95）。
/// コーパスでは 20 枚でしたが、実データでは **1 演算あたり最大 59,465 枚（出力の 0.95%）**です。
///
/// **述語は使いません**（§2.4.4 (2) の禁止。幾何で判定すると $20b+43$ が要る）。
/// **角が一般の位置にあること**（3 つが共線でないこと）だけを前提にします。
/// これは `drop_collinear` が上流で保証します。
///
/// **向きは保たれます。** 元の三角形 $(a,b,o)$ の有向辺 $(a,b)$ を
/// $(a,t_1,o), (t_1,t_2,o), \dots, (t_k,b,o)$ に置き換えるので、
/// **各小三角形が同じ向きの部分辺を持ちます。**
inline bool fan_triangulate_general(const TPolygon& p,
                                    std::vector<std::array<std::uint32_t, 3>>& out) {
    const std::size_t n = p.vertex.size();
    if (p.corners < 3 || n < 3) return false;
    // 角の頂点 ID と、多角形内での位置（元の添字順）
    std::vector<std::uint32_t> corner(p.corners, 0);
    std::vector<std::size_t> cpos(p.corners, 0);
    std::vector<char> seen(p.corners, 0);
    for (std::size_t i = 0; i < n; ++i) {
        if (!p.is_corner[i]) continue;
        const std::uint32_t o = p.orig[i];
        if (o >= p.corners || seen[o]) return false;  // 想定外。従来の経路に任せる
        corner[o] = p.vertex[i];
        cpos[o] = i;
        seen[o] = 1;
    }
    for (std::uint32_t o = 0; o < p.corners; ++o)
        if (!seen[o]) return false;

    // ---- 1. 角だけの扇分割 ----
    std::vector<std::array<std::uint32_t, 3>> tri;
    tri.reserve(n);
    for (std::uint32_t i = 1; i + 1 < p.corners; ++i)
        tri.push_back({corner[0], corner[i], corner[i + 1]});

    // ---- 2. T 頂点を持つ元の辺ごとに、その辺を含む三角形を対頂点から扇に開く ----
    std::vector<std::uint32_t> ts;
    for (std::uint32_t j = 0; j < p.corners; ++j) {
        ts.clear();
        for (std::size_t k = 1; k < n; ++k) {
            const std::size_t i = (cpos[j] + k) % n;
            if (p.is_corner[i]) break;
            if (p.orig[i] != j) return false;  // 並びが想定と違う。従来の経路に任せる
            ts.push_back(p.vertex[i]);
        }
        if (ts.empty()) continue;
        const std::uint32_t a = corner[j], b = corner[(j + 1) % p.corners];
        // **有向辺 (a,b) を持つ三角形はちょうど 1 つ**（扇分割は円板なので境界辺は 1 度だけ）
        std::size_t hit = tri.size();
        int slot = -1;
        for (std::size_t t = 0; t < tri.size() && hit == tri.size(); ++t)
            for (int e = 0; e < 3; ++e)
                if (tri[t][static_cast<std::size_t>(e)] == a &&
                    tri[t][static_cast<std::size_t>((e + 1) % 3)] == b) {
                    hit = t;
                    slot = e;
                    break;
                }
        if (hit == tri.size()) return false;  // 見つからない。従来の経路に任せる
        const std::uint32_t o = tri[hit][static_cast<std::size_t>((slot + 2) % 3)];
        // (a,t1,o), (t1,t2,o), ..., (tk,b,o)
        std::vector<std::array<std::uint32_t, 3>> fan;
        fan.reserve(ts.size() + 1);
        std::uint32_t prev = a;
        for (std::uint32_t t : ts) {
            fan.push_back({prev, t, o});
            prev = t;
        }
        fan.push_back({prev, b, o});
        tri[hit] = fan.front();
        tri.insert(tri.end(), fan.begin() + 1, fan.end());
    }
    // **枚数の検算**: 円板の扇分割は n-2 枚（$\chi = 1$）
    if (tri.size() != n - 2) return false;
    out.insert(out.end(), tri.begin(), tri.end());
    return true;
}

inline void fan_triangulate(const TPolygon& p, std::vector<std::array<std::uint32_t, 3>>& out,
                            TJunctionStats* stats = nullptr, bool general = true) {
    const std::size_t n = p.vertex.size();
    if (n < 3) return;

#if defined(KRISITE_MUTATION_DROP_DEGENERATE)
    // **変異 13 は「従来の扇分割が退化を捨てる」という欠陥を注入するもの**です。
    // **一般解は退化を 1 枚も作らないので、変異の対象が消えます。**
    //
    // **検出器を失わないよう、この変異では一般解を外して評価します**
    // （`CLAUDE.md`「後段で埋める機構は、上流の誤りを覆い隠します。
    // 導入するなら、どの検出器が失われ、どれが残るかを先に書き出してください」）。
    //
    // **失われた検出範囲**: 一般解が既定の経路では、変異 13 は到達しません。
    // **残る検出範囲**: 旗を外した従来の経路（`test_tjunction.cpp` が直接検査）。
    general = false;
#endif
    // **一般解を先に試します**（§2.4.4 (2) の「一般解」。退化を 1 枚も作りません）。
    // **旗で完全に外せます** — 偽なら下の従来の扇分割に落ちます（比較の基準側）。
    if (general) {
        const std::size_t before = out.size();
        if (fan_triangulate_general(p, out)) {
            if (stats) stats->general_used += out.size() - before;
            return;
        }
        out.resize(before);  // 途中まで書いていたら戻す
        if (stats) ++stats->general_fallback;
    }

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

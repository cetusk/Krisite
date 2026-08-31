// Krisite — CP5 全域（SPEC-phase2.md §13）
//
// **CP1〜CP4 で入った機構をすべて有効にして全域を回します。**
//
// 個々の機構は CP1〜CP4 で確かめました。**CP5 が新しく見るのは機構どうしの相互作用**で、
// どのチェックポイントでも単独では通らなかった経路が 3 つあります。
//
//   適応分割 × 分裂          接触辺がセル平面で分割された状態での分裂
//   構成点の保持 × T 解決    保持された点が T 頂点として挿入される経路
//   early-out × 分裂         arrangement を省いたセルの断片が分裂に回る経路
//
// **「全部有効にしたら通った」だけでは足りません。** 3 つの経路を**実際に踏んだこと**を
// 数えて、踏んでいなければ落とします（`CLAUDE.md`「機構を足したら、その機構が空回りして
// いないことを別に検査してください」）。計数は
//
//   TJunctionStats::inserted_from_cache   構成点の保持 × T 解決
//   SplitStats::split_from_early_out      early-out × 分裂
//   適応分割で葉の深さに差があり、かつ分裂が起きた構成の数
//
// ---
//
// ## 分裂 on / off で期待値が反転します（§13 の表）
//
// | | 分裂 OFF | 分裂 ON |
// |---|---|---|
// | §9.3 の除外            | 3 構成       | **0 件** |
// | 値の重複する頂点が無い | 成立         | **成立しない**（分裂が値を複製する） |
// | $\chi$ の偶奇          | 奇数になり得る | 常に偶数 |
//
// **同じテストで両方を見ようとしないこと。** ここは分裂 ON / OFF を**別々に**判定し、
// 両者をまたぐ検査は「分裂で変わってはいけない量」（$C$、$F$、幾何）に限ります。
//
// ## §9.4.2 分裂の不変量は、直接・常時・安く検査する
//
// **体積そのものは GMP が要ります**（構成点は有理数で、三角形ごとに分母が違う。
// `test_volume_gmp.cpp` の注記）。そこで常時走る側は**より強い不変量**で見ます
// （§9.4.2 の分担。体積そのものは `volume_gmp` が GMP ジョブで見ます）。
//
// > **分裂の前後で、頂点の【値】で見た三角形の多重集合が一致すること。**
//
// 分裂は頂点 ID を付け替えるだけで幾何を動かさないので（§5.1.3）、これが一致すれば
// **どんな測度でも一致します。体積の一致はその系です。** 体積そのものの突き合わせは
// `volume_gmp`（CI の GMP ジョブ）が担当します。
#include <algorithm>
#include <array>
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "corpus_expect.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kMaxDepth = 3;
constexpr std::size_t kPasses = kMaxDepth + 2;  ///< 固定深度 0〜3 + 出荷時構成

/// **CP1〜CP4 の機構をすべて有効にした構成。**
///
/// **既定値に依存しません**（§9.4 の CI ジョブで既定が反転するため）。
BoolOptions all_on(unsigned depth, bool adaptive) {
    BoolOptions o;
    o.depth = depth;
    o.cull_planes = true;     // CP1
    o.adaptive = adaptive;    // CP2
    o.leaf_threshold = 0;     // CP2
    o.early_out = true;       // CP2
    o.cache_points = true;    // CP3
    o.split_contacts = true;  // CP4（呼び出し側で上書きする）
    return o;
}

const char* op_name(BoolOp op) {
    switch (op) {
        case BoolOp::Union:
            return "∪";
        case BoolOp::Intersection:
            return "∩";
        default:
            return "\\";
    }
}

/// 頂点の**値**で見た三角形の多重集合（§9.4.2）。
///
/// 値が等しい頂点に同じ番号を振り直し、三角形は**巡回順を保ったまま**最小の番号が
/// 先頭に来るよう回します（向きを落とさないため、並べ替えはしません）。
std::vector<std::array<std::uint32_t, 3>> geometric_key(const BoolMesh& m) {
    const std::size_t n = m.vertices.size();
    std::vector<std::uint32_t> order(n);
    for (std::uint32_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
        return krisite::geom::lex_less(m.vertices[a], m.vertices[b]);
    });
    std::vector<std::uint32_t> canon(n, 0);
    std::uint32_t next = 0;
    for (std::size_t i = 0; i < order.size();) {
        std::size_t j = i;
        while (j < order.size() &&
               krisite::geom::h_equal(m.vertices[order[i]], m.vertices[order[j]])) {
            canon[order[j]] = next;
            ++j;
        }
        ++next;
        i = j;
    }
    std::vector<std::array<std::uint32_t, 3>> out;
    out.reserve(m.triangles.size());
    for (const krisite::mesh::Tri& t : m.triangles) {
        const std::uint32_t c[3] = {canon[t[0]], canon[t[1]], canon[t[2]]};
        int s = 0;
        if (c[1] < c[s]) s = 1;
        if (c[2] < c[s]) s = 2;
        out.push_back({c[s], c[(s + 1) % 3], c[(s + 2) % 3]});
    }
    std::sort(out.begin(), out.end());
    return out;
}

struct Totals {
    std::size_t configs = 0;
    std::size_t split_configs = 0;  ///< 分裂が実際に起きた構成
    // ---- 相互作用の計数（**空回りの番人**）----
    std::size_t adaptive_x_split = 0;   ///< 葉の深さに差があり、かつ分裂が起きた構成
    std::size_t cache_x_tjunction = 0;  ///< 保持された点が T 頂点として入った回数
    /// early-out 由来の三角形に接する頂点の分裂。**retire 済み。0 が正**（§13）
    std::size_t early_out_x_split = 0;
    std::size_t uneven_leaves = 0;  ///< 葉の深さに差が出た構成
    std::size_t early_out_cells = 0;
    std::size_t cache_hits = 0;
    std::size_t t_inserted = 0;
};

Totals g;

void run_case(const kritest::Case& c) {
    const TriMesh a = c.make_a(), b = c.make_b();

    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        for (std::size_t pass = 0; pass < kPasses; ++pass) {
            const bool adaptive = (pass == kPasses - 1);
            const unsigned depth = adaptive ? kMaxDepth : static_cast<unsigned>(pass);
            const std::string tag =
                std::string("ケース ") + c.id + " " + op_name(op) + "（" +
                (adaptive ? std::string("出荷時構成") : ("深度 " + std::to_string(depth))) + "）";

            // 正解器は **Phase 1 の構成**（機構をすべて外した固定深度）です（§0.1）。
            BoolOptions base_off = kritest::phase1_options(depth);
            BoolOptions base_on = base_off;
            base_on.split_contacts = true;
            BoolOptions on_off = all_on(depth, adaptive);
            on_off.split_contacts = false;
            BoolOptions on_on = all_on(depth, adaptive);

            BoolStats s_base_off, s_base_on, s_off, s_on;
            const BoolMesh r_base_off = boolean_op(a, b, op, base_off, &s_base_off);
            const BoolMesh r_base_on = boolean_op(a, b, op, base_on, &s_base_on);
            const BoolMesh r_off = boolean_op(a, b, op, on_off, &s_off);
            const BoolMesh r_on = boolean_op(a, b, op, on_on, &s_on);
            const TopologyReport t_base_off = check_topology(r_base_off.triangles);
            const TopologyReport t_base_on = check_topology(r_base_on.triangles);
            const TopologyReport t_off = check_topology(r_off.triangles);
            const TopologyReport t_on = check_topology(r_on.triangles);
            ++g.configs;

            // ---- (1) 機構を全部入れても答えが変わらないこと（§9.1 の一般化）----
            //
            // **分裂 ON / OFF を別々に突き合わせます。** 期待値が反転する量
            // （除外・χ の偶奇・値の重複）を混ぜないためです（§13 の表）。
            KRI_CHECK_MSG(t_off.components == t_base_off.components,
                          tag + "（分裂 OFF）: C が Phase 1 の構成と違う" +
                              kritest::pair_msg(t_base_off.components, t_off.components));
            KRI_CHECK_MSG(t_off.chi == t_base_off.chi,
                          tag + "（分裂 OFF）: χ が Phase 1 の構成と違う" +
                              kritest::pair_msg(t_base_off.chi, t_off.chi));
            KRI_CHECK_MSG(t_on.components == t_base_on.components,
                          tag + "（分裂 ON）: C が Phase 1 の構成と違う" +
                              kritest::pair_msg(t_base_on.components, t_on.components));
            KRI_CHECK_MSG(t_on.chi == t_base_on.chi,
                          tag + "（分裂 ON）: χ が Phase 1 の構成と違う" +
                              kritest::pair_msg(t_base_on.chi, t_on.chi));

            // ---- (2) 分裂 ON の意味論（§5）----
            //
            // **§5.1.2.2: 対応付けできなかった辺があるときは除外します**（案 A）。
            // 判定は識別子ではなく `unresolved > 0` という機械的な事実です。
            const auto& sp = s_on.split;
            const kritest::Exclusion ex_on = kritest::exclusion_when_split(sp.unresolved, t_on);
            if (ex_on != kritest::Exclusion::None) {
                std::string why;
                KRI_CHECK_MSG(kritest::exclusion_conditions_ok(ex_on, t_on, &why),
                              tag + ": §9.3.1 の適用条件を満たさない（" + why + "）");
                // **裂けていないこと**（§5.1.2.2）
                KRI_CHECK_MSG(t_on.edges_deficient == 0, tag + ": **辺が裂けました**（次数 1 が " +
                                                             std::to_string(t_on.edges_deficient) +
                                                             " 本）");
            } else {
                KRI_CHECK_MSG(t_on.ok(), tag + ": 分裂後も多様体になっていない（§5.3）");
            }
            if (!t_on.empty) KRI_CHECK_MSG(t_on.chi_even, tag + ": 分裂後も χ が奇数（§5.4）");
            KRI_CHECK_MSG(sp.predicted_delta_v == sp.actual_delta_v,
                          tag + ": ΔV の予測と実測が違う" +
                              kritest::pair_msg(sp.predicted_delta_v, sp.actual_delta_v));
            KRI_CHECK_MSG(sp.predicted_delta_e == sp.actual_delta_e,
                          tag + ": ΔE の予測と実測が違う" +
                              kritest::pair_msg(sp.predicted_delta_e, sp.actual_delta_e));
            KRI_CHECK_MSG(sp.predicted_delta_chi == t_on.chi - t_off.chi,
                          tag + ": Δχ が出力の χ の差と一致しない");

            // ---- (3) 分裂で変わってはいけない量（§5.1.2 / §5.1.3）----
            KRI_CHECK_MSG(t_on.components == t_off.components,
                          tag + ": **分裂で C が変わった**（§5.5.1）" +
                              kritest::pair_msg(t_off.components, t_on.components));
            KRI_CHECK_MSG(t_on.f == t_off.f, tag + ": 分裂で面が増減した（§5.1.3）");

            // ---- (4) §9.4.2 分裂の有無で幾何が一致すること（体積の一致はこの系）----
            KRI_CHECK_MSG(geometric_key(r_off) == geometric_key(r_on),
                          tag +
                              ": **分裂で三角形の幾何が変わった**（§9.4.2）。"
                              "分裂は頂点 ID の付け替えだけのはずです");

            // ---- 相互作用の計数 ----
            if (sp.split_vertices > 0) ++g.split_configs;
            if (s_on.leaf_depth_max > s_on.leaf_depth_min) {
                ++g.uneven_leaves;
                if (sp.split_vertices > 0) ++g.adaptive_x_split;
            }
            g.cache_x_tjunction += s_on.t.inserted_from_cache;
            g.early_out_x_split += sp.split_from_early_out;
            g.early_out_cells += s_on.early_out_cells;
            g.cache_hits += s_on.cache_hits;
            g.t_inserted += s_on.t.inserted;
        }
    }
}

/// **空回りの番人。** 3 つの経路を踏んでいなければ、CP5 は何も新しく検証していません。
void check_not_vacuous() {
    // **期待値は実測ではなく式で持ちます**（`CLAUDE.md`）。
    const std::size_t expect = kritest::corpus().size() * 3 * kPasses;
    std::printf("\n  §13 CP5 の記録（%zu ケース × 3 演算 × %zu パス = %zu 構成、期待 %zu）\n",
                kritest::corpus().size(), kPasses, g.configs, expect);
    std::printf("    分裂が起きた構成 %zu / 葉の深さに差が出た構成 %zu\n", g.split_configs,
                g.uneven_leaves);
    std::printf("    T 頂点 %zu（うち保持された点 %zu）/ キャッシュ命中 %zu / early-out セル %zu\n",
                g.t_inserted, g.cache_x_tjunction, g.cache_hits, g.early_out_cells);
    std::printf("\n    **相互作用（どれかが 0 なら CP5 は空回りです）**\n");
    std::printf("      適応分割 × 分裂        %zu 構成\n", g.adaptive_x_split);
    std::printf("      構成点の保持 × T 解決  %zu 回\n", g.cache_x_tjunction);
    std::printf("      early-out × 分裂       %zu 回（retire。0 が正）\n", g.early_out_x_split);

    KRI_CHECK_MSG(g.configs == expect, "構成数が期待と違う。**空回りの疑い**");
    KRI_CHECK_MSG(g.adaptive_x_split > 0,
                  "**適応分割 × 分裂を 1 度も踏んでいない。** 葉の深さに差がある状態で"
                  "接触が分裂する配置がコーパスにありません");
    KRI_CHECK_MSG(g.cache_x_tjunction > 0,
                  "**構成点の保持 × T 解決を 1 度も踏んでいない。** 保持された点が"
                  "T 頂点として挿入される配置がコーパスにありません");
    // **early-out × 分裂 は retire しました**（`SPEC-phase2.md` §13、2026-08-28）。
    //
    // **「作れなかったから外す」ではなく、構造的に到達不能だから外します**
    // （`CLAUDE.md`「番人を retire するときは、証明とセットにしてください」）。
    //
    // > 早期脱出したセル $C$ の**閉じた箱**に相手の三角形が 1 枚も無いなら、
    // > 相手の曲面は $\overline{C}$ に届きません。接触辺（次数 4）は $A$ の 2 枚と
    // > $B$ の 2 枚が辺を共有する形なので、**$C$ 由来の断片が接触に加わるには
    // > 相手の曲面が $\overline{C}$ に届く必要があり、矛盾します。**
    //
    // **Phase 2 では 19 回発火していましたが、それはすべてバグ経由でした。**
    // early-out の存在判定に**割り当て用の半開区間の述語**を使っており、面がセルの
    // 上側境界に乗ると「存在しない」と誤答していました（`IMPL-phase3.md` §7.1）。
    // 閉領域（`overlaps_cell`）で判定するのが正しく、直したら 19 → 0 になりました。
    //
    // ---- ★ 訂正（2026-08-31）: **retire の証明が反証されました** ----
    //
    // 「早期脱出したセルの閉じた箱には**相手の**曲面が無いので到達不能」という証明は、
    // **接触が 2 つの入力の【間】に起きることを暗黙に仮定**していました。
    //
    // **1 つの入力の【中】で接触する場合**（ケース 24 / 24′ の A は、辺だけで接する
    // 2 つの箱）、あるセルで相手が居なくても**自分自身の接触辺が居ます。**
    // したがって early-out したセルの断片に接する頂点が分裂し得ます。**到達可能です。**
    //
    // **したがって「0 が正」には戻しません。記録に留めます**（`CLAUDE.md`
    // 「番人を retire するときは証明とセットで」— 証明が誤りだったので retire を撤回）。
    //
    // **元のバグ（半開区間の述語）を捕まえる力は失われます。** 代わりの検出は
    // `overlaps_cell` の単体テスト（`tests/octree/`）が受け持ちます。
    std::printf("      ※ early-out × 分裂は **記録のみ**（retire の証明が反証。自己接触で到達）\n");
}

}  // namespace

int main() {
    std::printf("\n  CP5 全域 — 機構どうしの相互作用（SPEC-phase2 §13）\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) run_case(c);
    check_not_vacuous();
    std::printf("\n");
    return kritest::finish("csg/cp5");
}

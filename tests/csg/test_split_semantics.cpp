// Krisite — 非多様体出力の意味論（SPEC-phase2.md §5）
//
// **正則化ブールは一般に多様体出力を保証できません**（`SPEC-phase1.md` §9.3.2）。
// Phase 1 は 54 構成のうち 3 件を §9.3 で除外していました。**分裂で 0 件になるはずです**（§5.3）。
//
// ---
//
// ## この検査の背骨は $C$ の不変性です ★
//
// §5.5.1:
//
// > **誤った対応付けは、これでしか捕まりません。**
// > $\Delta V$ と $\Delta E$ は組の分け方に依りません（常に辺 1 本と端点 2 個が増える）。
// > 面も増えません。したがって **$\chi$ も体積も、誤った組を検出しません。**
//
// 誤った組で分けると、本来別のシートだった面が次数 2 の辺で繋がり、**成分が併合されて
// $C$ が減ります。** §9.3 の変異 9（対応付けを owner に戻す）は**この検査でのみ落ちます。**
//
// ## 比較は $(C, \chi)$ で行い、$g$ は使いません（§5.4）
//
// $g = C - \chi/2$ は $\chi$ が偶数のときにしか種数を意味しません。**非多様体出力では
// $\chi$ が奇数になり得ます**（11b の $\cup$ は $\chi=3$）。そこで $g$ を計算すると、
// 整数除算の産物が「種数 1」として表示されます。**Phase 1 で実際に踏んだ誤りです。**
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

struct Totals {
    std::size_t split_vertices = 0;
    std::size_t excess_edges = 0, excess_endpoints = 0;
    std::size_t max_fans = 0;
    std::size_t unresolved = 0;
    std::size_t configs = 0;             ///< 検査した構成の数
    std::size_t nonmanifold_before = 0;  ///< 分裂前に非多様体だった構成
    std::size_t chi_odd_before = 0;      ///< 分裂前に χ が奇数だった構成
    std::size_t excluded_after = 0;      ///< 分裂後に除外が必要だった構成（**0 のはず**）
};

Totals g;

void run_case(const kritest::Case& c) {
    const TriMesh a = c.make_a(), b = c.make_b();
    bool printed = false;

    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        for (unsigned d = 0; d <= kMaxDepth; ++d) {
            BoolOptions off;
            off.depth = d;
            off.split_contacts = false;
            BoolOptions on = off;
            on.split_contacts = true;

            BoolStats s_off, s_on;
            const BoolMesh r_off = boolean_op(a, b, op, off, &s_off);
            const BoolMesh r_on = boolean_op(a, b, op, on, &s_on);
            const TopologyReport t_off = check_topology(r_off.triangles);
            const TopologyReport t_on = check_topology(r_on.triangles);
            ++g.configs;

            const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) + "（深度 " +
                                    std::to_string(d) + "）";

            // ---- §5.5.1 の C 不変性 ★ 対応付けの唯一の番人 ----
            KRI_CHECK_MSG(t_on.components == t_off.components,
                          tag + ": **分裂で C が変わった**（" + std::to_string(t_off.components) +
                              " → " + std::to_string(t_on.components) +
                              "）。§5.1.2 の対応付けが誤っています");
            // ---- §5.1.3: 面は増えない ----
            KRI_CHECK_MSG(t_on.f == t_off.f, tag + ": 面が増減した（§5.1.3）");

            // ---- §5.3: 分裂後は除外 0 件 ----
            //
            // **空メッシュは別扱いです。** `check_topology` は空を V=E=F=0 で見るので、
            // 多様体性や χ の偶奇の旗は立ちません（`ok()` はそれを織り込んでいます）。
            KRI_CHECK_MSG(t_on.ok(), tag +
                                         ": 分裂後も多様体になっていない（§5.3 の帰結が"
                                         "成立していません）");
            if (!t_on.ok()) ++g.excluded_after;
            if (!t_on.empty) {
                // §5.4: **分裂後は χ が偶数**。奇数なら閉曲面ではありません
                KRI_CHECK_MSG(t_on.chi_even, tag + ": 分裂後も χ が奇数（§5.4）");
            }

            // ---- §5.5 の検算: 予測と実測 ----
            const auto& sp = s_on.split;
            KRI_CHECK_MSG(sp.predicted_delta_v == sp.actual_delta_v,
                          tag + ": ΔV の予測 " + std::to_string(sp.predicted_delta_v) + " と実測 " +
                              std::to_string(sp.actual_delta_v) + " が違う");
            KRI_CHECK_MSG(sp.predicted_delta_e == sp.actual_delta_e,
                          tag + ": ΔE の予測 " + std::to_string(sp.predicted_delta_e) + " と実測 " +
                              std::to_string(sp.actual_delta_e) + " が違う");
            KRI_CHECK_MSG(sp.predicted_delta_chi == sp.actual_delta_chi,
                          tag + ": Δχ の予測と実測が違う");
            KRI_CHECK_MSG(sp.predicted_delta_chi == t_on.chi - t_off.chi,
                          tag + ": Δχ が出力の χ の差と一致しない");
            // ---- §5.1.2.1: 分けられない配置に到達したか ----
            KRI_CHECK_MSG(sp.unresolved == 0,
                          tag +
                              ": **分裂しても多様体にならない配置に到達しました**（§5.1.2.1）。"
                              "radial sort が要る配置です。**記録して報告してください**");

            g.split_vertices += sp.split_vertices;
            g.excess_edges += sp.excess_edges;
            g.excess_endpoints += sp.excess_endpoints;
            g.max_fans = std::max(g.max_fans, sp.max_fans);
            g.unresolved += sp.unresolved;
            if (!t_off.ok() && !t_off.empty) ++g.nonmanifold_before;
            if (!t_off.empty && !t_off.chi_even) ++g.chi_odd_before;

            if (sp.split_vertices > 0) {
                if (!printed) {
                    std::printf("\n  ケース %-4s %s\n", c.id, c.what);
                    printed = true;
                }
                std::printf(
                    "    %-2s d%u 次数≥3 の辺 %zu(端点 %zu) 扇の最大 %zu | "
                    "ΔV %zu ΔE %zu Δχ %+lld | χ %lld→%-3lld C %zu→%zu\n",
                    op_name(op), d, sp.excess_edges, sp.excess_endpoints, sp.max_fans,
                    sp.actual_delta_v, sp.actual_delta_e, sp.actual_delta_chi, t_off.chi, t_on.chi,
                    t_off.components, t_on.components);
            }
        }
    }
}

/// **空回りの番人。** 分裂が一度も起きないなら、この検査は何も検証していません。
void check_not_vacuous() {
    std::printf("\n  §5.7 の記録（全ケース × 3 演算 × 深度 0〜%u = %zu 構成）\n", kMaxDepth,
                g.configs);
    std::printf("    分裂した頂点 %zu / 次数≥3 の辺 %zu（端点 %zu）/ 扇の最大 k = %zu\n",
                g.split_vertices, g.excess_edges, g.excess_endpoints, g.max_fans);
    std::printf("    分裂前に非多様体だった構成 %zu / χ が奇数だった構成 %zu\n",
                g.nonmanifold_before, g.chi_odd_before);
    std::printf("    **分裂後に除外が必要だった構成 %zu（0 のはず。§5.3）**\n", g.excluded_after);
    std::printf("    §5.1.2.1 の停止に到達した配置 %zu（radial sort の必要性の判断材料）\n",
                g.unresolved);

    KRI_CHECK_MSG(g.split_vertices > 0, "分裂が一度も起きていない。**空回りです**");
    KRI_CHECK_MSG(g.nonmanifold_before > 0,
                  "分裂前に非多様体になる構成が 1 つも無い。**コーパスが接触を作れていません**");
    KRI_CHECK_MSG(g.chi_odd_before > 0,
                  "分裂前に χ が奇数になる構成が無い。**辺接触（次元 1）を突けていません**");
    KRI_CHECK_MSG(g.excess_edges > 0, "次数 3 以上の辺が 1 本も無い。**辺の分裂を突けていません**");
    KRI_CHECK_MSG(g.max_fans >= 3,
                  "扇の最大が 3 未満。**k=3 の分裂（4T′ の A\\B）を突けていません**");
    KRI_CHECK_MSG(g.excluded_after == 0, "§9.3 の除外が 0 件になっていない（§5.3 の帰結）");
}

}  // namespace

int main() {
    std::printf("\n  接触の分裂 — SPEC-phase2 §5 / CP4\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) run_case(c);
    check_not_vacuous();
    std::printf("\n");
    return kritest::finish("csg/split_semantics");
}

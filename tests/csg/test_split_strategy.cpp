// Krisite — 分割戦略不変性（SPEC-phase2.md §9.1）
//
// > 同一ケース・同一演算について、分割戦略を変えても出力の体積と位相 $(C, \chi)$ が
// > 厳密に一致すること。
//
// **Phase 1 の実装がそのまま正解器です**（§0.1）。外部ライブラリもネットワークも
// 要りません。絞り込みを無効にした側が Phase 1 の挙動そのものです。
//
// ---
//
// **CP1 では一致はもっと強い形で成り立つはずです。**
//
// §2.3 の絞り込みが落とすのは「セルの閉包と交わらない平面」です。断片はセルの
// 閉包に含まれるので、そういう平面では**そもそも切れません**。つまり絞り込みは
// no-op な分割を省いているだけで、**断片の集合は完全に同一**になります。
//
// そこでここでは $(C, \chi)$ より強く **$(V, E, F)$ の一致**まで要求します。
// **CP2（適応分割）では三角形分割が変わるのでこの強い形は使えません。**
// そのときは $(C, \chi)$ と体積に緩めてください。強い形が落ちたら、それは
// 「緩めるべき」ではなく「絞り込みが no-op でないものを落とした」の合図です。
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

/// §11 の記録: 絞り込みの効き。
struct Culling {
    std::size_t slots = 0;  ///< セル × 分割平面の総当たり数
    std::size_t used = 0;   ///< 絞り込んだ後に実際に使った数
    std::size_t max_per_cell = 0;
};

Culling g_cull[kMaxDepth + 1];

/// §2.4.4 (3) の【負の対照】: **固定深度では T 頂点は 1 個も入らないはず。**
///
/// 大域平面集合で切るので、辺の内部に載る大域頂点はその平面での切断点として既に
/// 多角形の角になっています。**入るなら前提が崩れています。**
///
/// **ただし負の対照だけでは何もしないスタブが通ります。** 索引が候補を返していることを
/// 併せて数えます。候補が 0 なら、この検査は何も検証していません。
/// 効くべき場所で効くこと（正の対照）は `tests/csg/test_tjunction.cpp` が持ちます。
struct TControl {
    std::size_t inserted = 0;
    std::size_t candidates = 0;
    std::size_t edges = 0;
    std::size_t dropped = 0;
};

TControl g_t;

void run_case(const kritest::Case& c) {
    const TriMesh a = c.make_a(), b = c.make_b();
    std::printf("\n  ケース %-4s %s\n", c.id, c.what);

    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        for (unsigned d = 0; d <= kMaxDepth; ++d) {
            BoolStats s_off, s_on;
            // 絞り込み無効 = Phase 1 の挙動（正解器）
            const BoolMesh r_off = boolean_op(a, b, op, d, &s_off, false);
            const BoolMesh r_on = boolean_op(a, b, op, d, &s_on, true);
            const TopologyReport t_off = check_topology(r_off.triangles);
            const TopologyReport t_on = check_topology(r_on.triangles);

            const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) + "（深度 " +
                                    std::to_string(d) + "）";

            // ---- §9.1 が要求する形: (C, χ) の一致 ----
            KRI_CHECK_MSG(t_on.components == t_off.components,
                          tag + ": 絞り込みで C が変わった（" + std::to_string(t_off.components) +
                              " → " + std::to_string(t_on.components) + "）");
            KRI_CHECK_MSG(t_on.chi == t_off.chi, tag + ": 絞り込みで χ が変わった（" +
                                                     std::to_string(t_off.chi) + " → " +
                                                     std::to_string(t_on.chi) + "）");
            // **χ の偶奇を併せて記録する**（§5.4）。奇数なら g は種数ではない
            KRI_CHECK_MSG(t_on.chi_even == t_off.chi_even, tag + ": χ の偶奇が変わった");

            // ---- CP1 で成り立つはずの強い形: (V, E, F) まで一致 ----
            KRI_CHECK_MSG(t_on.v == t_off.v && t_on.e == t_off.e && t_on.f == t_off.f,
                          tag + ": 固定深度なのに (V,E,F) が変わった（" + std::to_string(t_off.v) +
                              "/" + std::to_string(t_off.e) + "/" + std::to_string(t_off.f) +
                              " → " + std::to_string(t_on.v) + "/" + std::to_string(t_on.e) + "/" +
                              std::to_string(t_on.f) +
                              "）。§2.3 は no-op な分割だけを省くはずです");
            // 位相の性質そのものも一致すること
            KRI_CHECK_MSG(t_on.edge_manifold == t_off.edge_manifold &&
                              t_on.vertex_manifold == t_off.vertex_manifold &&
                              t_on.oriented == t_off.oriented &&
                              t_on.no_degenerate == t_off.no_degenerate,
                          tag + ": 絞り込みで位相の性質が変わった");
            // 断片数も変わらないはず（no-op を省いただけなので）
            KRI_CHECK_MSG(s_on.fragments == s_off.fragments,
                          tag + ": 断片数が変わった（" + std::to_string(s_off.fragments) + " → " +
                              std::to_string(s_on.fragments) + "）");

            // §2.4.4 (3) の負の対照
            KRI_CHECK_MSG(s_on.t.inserted == 0,
                          tag + ": 固定深度なのに T 頂点が入った（" +
                              std::to_string(s_on.t.inserted) +
                              " 個）。大域平面集合では切断点が既に角になっているはずです");
            KRI_CHECK_MSG(s_off.t.inserted == 0, tag + ": 絞り込み無効側で T 頂点が入った");
            KRI_CHECK_MSG(s_on.t.degenerate_kept == 0, tag + ": T 頂点が無いのに退化三角形が出た");
            g_t.inserted += s_on.t.inserted;
            g_t.candidates += s_on.t.candidates;
            g_t.edges += s_on.t.edges_scanned;
            g_t.dropped += s_on.t.degenerate_kept;

            g_cull[d].slots += s_on.split_plane_slots;
            g_cull[d].used += s_on.split_planes_used;
            g_cull[d].max_per_cell = std::max(g_cull[d].max_per_cell, s_on.max_planes_per_cell);

            if (op == BoolOp::Union) {
                const double pct = s_on.split_plane_slots
                                       ? 100.0 * static_cast<double>(s_on.split_planes_used) /
                                             static_cast<double>(s_on.split_plane_slots)
                                       : 100.0;
                std::printf("    d%u 平面 %zu/%zu (%.1f%%) 最大%zu/セル | C=%zu χ=%-3lld F=%zu\n",
                            d, s_on.split_planes_used, s_on.split_plane_slots, pct,
                            s_on.max_planes_per_cell, t_on.components, t_on.chi, t_on.f);
            }
        }
    }
}

/// §9.0 の番人: **絞り込みが実際に効いていること。**
///
/// 何も落とさないなら、この検査は「常に同じ結果」を確かめているだけで、
/// §2.3 の正しさを一切検証していません。
void check_culling_effective() {
    std::printf("\n  §11 絞り込みの効き（全ケース × 3 演算の合計）\n");
    std::printf("    %-6s %-12s %-12s %-8s %s\n", "深度", "使用", "総当たり", "残存率",
                "最大/セル");
    for (unsigned d = 0; d <= kMaxDepth; ++d) {
        const Culling& c = g_cull[d];
        const double pct =
            c.slots ? 100.0 * static_cast<double>(c.used) / static_cast<double>(c.slots) : 100.0;
        std::printf("    d%-5u %-12zu %-12zu %-7.1f%% %zu\n", d, c.used, c.slots, pct,
                    c.max_per_cell);
    }
    // **深度 0 では落ちません。** 入力点はすべて座標範囲内なので、そこから作った
    // 平面は必ず深度 0 のセル（座標範囲全体）と交わります。
    KRI_CHECK_MSG(g_cull[0].used == g_cull[0].slots,
                  "深度 0 で平面が落ちた。入力から作った平面は座標範囲と必ず交わるはずです");
    // 深度が上がれば落ちること。落ちないなら §2.3 は空回りしています
    for (unsigned d = 1; d <= kMaxDepth; ++d) {
        KRI_CHECK_MSG(g_cull[d].used < g_cull[d].slots,
                      "深度 " + std::to_string(d) +
                          " で 1 枚も落ちていない。§2.3 の絞り込みが空回りしています");
    }
    // 深度が上がるほど残存率は下がるはず（セルが小さくなるので）
    for (unsigned d = 2; d <= kMaxDepth; ++d) {
        const double prev =
            static_cast<double>(g_cull[d - 1].used) / static_cast<double>(g_cull[d - 1].slots);
        const double cur =
            static_cast<double>(g_cull[d].used) / static_cast<double>(g_cull[d].slots);
        KRI_CHECK_MSG(cur <= prev + 1e-12, "深度 " + std::to_string(d) +
                                               " で残存率が上がった（セルは小さくなるので " +
                                               "下がるはずです）");
    }
}

/// §2.4.4 (3) の負の対照が**空虚でない**ことの番人。
void check_t_control_not_vacuous() {
    std::printf("\n  §2.4.4 (3) 負の対照（固定深度）\n");
    std::printf("    挿入 %zu 件 / 走査した辺 %zu 本 / 索引が返した候補 %zu 件 / 退化 %zu 枚\n",
                g_t.inserted, g_t.edges, g_t.candidates, g_t.dropped);
    KRI_CHECK_MSG(g_t.inserted == 0, "固定深度で T 頂点が入った");
    KRI_CHECK_MSG(g_t.edges > 0, "辺を 1 本も走査していない");
    // **ここが要点です。** 候補が 0 なら「挿入 0 件」は索引が空なだけかもしれません。
    KRI_CHECK_MSG(g_t.candidates > 0,
                  "索引が候補を 1 件も返していない。**負の対照が空虚です**"
                  "（何もしないスタブでも通ってしまう状態）");
    // 候補があるのに 1 つも入らない = 判定が効いている、という形になっていること
    KRI_CHECK_MSG(g_t.candidates >= g_t.edges,
                  "辺あたりの候補が 1 未満。索引の引き方が間違っている疑い");
}

}  // namespace

int main() {
    std::printf("\n  分割戦略不変性 — SPEC-phase2 §9.1 / CP1\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) run_case(c);
    check_culling_effective();
    check_t_control_not_vacuous();
    std::printf("\n");
    return kritest::finish("csg/split_strategy");
}

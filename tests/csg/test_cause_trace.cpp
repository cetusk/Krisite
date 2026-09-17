// Krisite — 記録基盤 2 本の限定検証（`SPEC-phase5.md` §5.10.14.149 §3・§4）
//
// **本体（`soup_boolean` など）へは接続しません。** 記録値と保存動作の対照までです。
// **模型の Boolean は回しません。**
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "krisite/csg/diagnostic_trace.hpp"
#include "thingi10k/cause_trace_io.hpp"

namespace diag = krisite::csg::diag;
namespace tio = krithingi::trace_io;

static int g_ok = 0, g_ng = 0;

static void chk(const char* name, bool got, bool want) {
    const bool ok = (got == want);
    (ok ? g_ok : g_ng)++;
    std::printf("  %s %-58s 得た値 %-6s 期待 %s\n", ok ? "ok  " : "**NG**", name,
                got ? "true" : "false", want ? "true" : "false");
}
static void chk_eq(const char* name, long long got, long long want) {
    const bool ok = (got == want);
    (ok ? g_ok : g_ng)++;
    std::printf("  %s %-58s 得た値 %-6lld 期待 %lld\n", ok ? "ok  " : "**NG**", name, got, want);
}

// --- 対照 1: source 別の混在記録（forced と winding_split が同じ領域に共存する）---
static void c1_per_source(diag::TraceBin& bin) {
    diag::RegionRec r;
    r.op = diag::Op::kUnion;
    r.cell = diag::CellKey{3, 1, 2, 3};
    r.region = diag::Ref{0, static_cast<std::uint8_t>(diag::Stage::kClassify), r.cell, 7};
    diag::RegionSourceRec s0;
    s0.source = diag::SourceId{0};
    s0.path = diag::WindingPath::kForced;
    s0.evaluated = true;
    s0.w_front = 1; s0.w_back = 0;
    s0.forced_from = diag::CellKey{3, 1, 2, 2};
    s0.formula = "forced(cell(3,1,2,2)) -> w_front=1, w_back=0";
    diag::RegionSourceRec s1;
    s1.source = diag::SourceId{1};
    s1.path = diag::WindingPath::kWindingSplit;
    s1.evaluated = true;
    s1.w_front = 2; s1.w_back = 1; s1.w_other = 2; s1.c_front = 1; s1.c_back = 0;
    s1.formula = "split(w_other=2, c_front=1, c_back=0) -> w_front=2, w_back=1";
    diag::RegionSourceRec s2;            // **未実施を「値 0」と混ぜない**
    s2.source = diag::SourceId{2};
    r.per_source = {s0, s1, s2};
    r.indicator_front = true;
    r.decision = diag::Decision::kAdopted;
    bin.regions.push_back(r);

    chk("C1 forced と winding_split が同じ領域に共存",
        bin.regions[0].per_source[0].path == diag::WindingPath::kForced
            && bin.regions[0].per_source[1].path == diag::WindingPath::kWindingSplit, true);
    chk("C1 未実施の source を値 0 と区別",
        bin.regions[0].per_source[2].evaluated == false
            && bin.regions[0].per_source[2].path == diag::WindingPath::kNotEvaluated, true);
    chk("C1 供給セルが source ごとに残る",
        bin.regions[0].per_source[0].forced_from == diag::CellKey{3, 1, 2, 2}, true);
}

// --- 対照 2: 接触の隅対応（旧頂点 → 隅 → 同値類 → 新頂点）---
static void c2_corner(diag::TraceBin& bin) {
    diag::VertexSplitRec v;
    v.op = diag::Op::kIsect;
    v.old_vertex = 42;
    v.corners = {
        {diag::OutFaceId{10}, 0, 0, 100},
        {diag::OutFaceId{11}, 2, 0, 100},
        {diag::OutFaceId{12}, 1, 1, 101},
    };
    v.class_count = 2;
    v.resolved = true;
    v.reason = diag::Reason::kContactResolved;
    bin.vertex_splits.push_back(v);

    diag::EdgeContactRec e;
    e.op = diag::Op::kIsect;
    e.edge_u = 42; e.edge_v = 43;
    e.branches = {diag::OutFaceId{10}, diag::OutFaceId{11}, diag::OutFaceId{12}, diag::OutFaceId{13}};
    e.pairs = {{0, 1}, {2, 3}};
    e.resolved = true;
    bin.edge_contacts.push_back(e);

    chk("C2 隅が (元面, 局所隅) で特定できる",
        bin.vertex_splits[0].corners[1].face == diag::OutFaceId{11}
            && bin.vertex_splits[0].corners[1].local_corner == 2, true);
    chk("C2 同じ類の 2 隅が同じ新頂点へ移る",
        bin.vertex_splits[0].corners[0].new_vertex == bin.vertex_splits[0].corners[1].new_vertex, true);
    chk("C2 別の類は別の新頂点",
        bin.vertex_splits[0].corners[2].new_vertex != bin.vertex_splits[0].corners[0].new_vertex, true);
    chk("C2 辺の巡回順と組を別記録で持つ",
        bin.edge_contacts[0].branches.size() == 4 && bin.edge_contacts[0].pairs.size() == 2, true);
}

// --- 対照 3: 結合順と ID（順序を変えても参照が一致する）---
static void c3_merge_order() {
    auto make = [](std::uint32_t local, std::uint32_t k) {
        diag::TraceBin b;
        diag::FragRec f;
        f.cell = diag::CellKey{4, 1, 1, k};
        f.frag = diag::Ref{0, static_cast<std::uint8_t>(diag::Stage::kArrange), f.cell, local};
        f.decision = diag::Decision::kGenerated;
        b.frags.push_back(f);
        return b;
    };
    diag::TraceBin ab = make(0, 5);
    ab.merge(make(0, 9));
    diag::TraceBin ba = make(0, 9);
    ba.merge(make(0, 5));

    // **局所番号は両方 0。セル鍵まで含めて初めて別物と分かります。**
    chk("C3 局所番号だけでは同じに見える", ab.frags[0].frag.local == ba.frags[0].frag.local, true);
    chk("C3 参照全体では結合順で別物と分かる", ab.frags[0].frag == ba.frags[0].frag, false);
    chk("C3 同じ内容は結合順を変えても一致", ab.frags[0].frag == ba.frags[1].frag, true);
    chk_eq("C3 結合で件数が保たれる", (long long)(ab.frags.size() + ba.frags.size()), 4);
}

// --- 対照 4: 全体の上限超過（各器が上限内でも全体で超える）---
static void c4_budget() {
    diag::TraceBudget budget(3);
    diag::CauseTrace t1(&budget), t2(&budget);
    int taken = 0;
    for (int i = 0; i < 2; ++i) taken += t1.room() ? 1 : 0;
    for (int i = 0; i < 2; ++i) taken += t2.room() ? 1 : 0;
    chk_eq("C4 共有予算を超えて取れない", taken, 3);
    chk("C4 超過が固定される", budget.exhausted(), true);
    chk("C4 どちらの器から見ても超過が分かる", t1.exhausted() && t2.exhausted(), true);
    chk("C4 落とした器が特定できる", t2.dropped_here(), true);
    chk("C4 落としていない器と区別できる", t1.dropped_here(), false);
    diag::CauseTrace t3;                 // 予算なし = 上限なし
    chk("C4 予算が無ければ落とさない", t3.room() && !t3.exhausted(), true);
}

// --- 対照 5: 排他的な新規作成と、書出し失敗の固定 ---
static void c5_io(const std::string& dir) {
    const std::string p = dir + "/trace_ok.jsonl";
    tio::TraceWriter w;
    chk("C5 新規なら開ける", w.open(p, "run-1", "250394x45413"), true);

    diag::TraceBin bin;
    c1_per_source(bin);
    c2_corner(bin);
    diag::VerifyRec v;
    v.op = diag::Op::kUnion;
    v.name = "verify_delta";
    v.target = "union";
    v.premise = "PWN";
    v.premise_basis = "**未確認**";
    v.executed = true;
    v.fire_count = 3;
    v.result.text = "-69366264756244193/1";
    v.reached_final = true;
    bin.verifies.push_back(v);
    diag::RemovedFaceRec rm;
    rm.source = diag::SourceId{0};
    rm.src_face = diag::SrcFaceId{17};
    rm.reason = diag::Reason::kSameIndexAfterQuant;   // **共線除去と混ぜない**
    rm.where = "loader.hpp:quantize\tcol\rreturn";     // **制御文字の往復**
    bin.removed.push_back(rm);

    chk("C5 書出しが成功", w.write(bin), true);
    chk("C5 close が成功", w.close(), true);
    chk("C5 失敗が固定されていない", w.failed(), false);

    tio::TraceWriter w2;                 // **同じ道を 2 度開けない（排他的作成）**
    chk("C5 既存の道は開けない", w2.open(p, "run-1", "x"), false);
    chk("C5 失敗が固定される", w2.failed(), true);

    // **既存ファイルの内容が変わっていないこと**（負の対照）
    std::FILE* f = std::fopen(p.c_str(), "rb");
    long size = 0;
    if (f) { std::fseek(f, 0, SEEK_END); size = std::ftell(f); std::fclose(f); }
    chk("C5 既存ファイルが空にされていない", size > 0, true);

    // **開かずに書いたら失敗が固定され、close でも取り消せない**
    tio::TraceWriter w3;
    diag::TraceBin empty;
    chk("C5 開かずに書けば失敗", w3.write(empty), false);
    chk("C5 close は先行失敗を取り消さない", w3.close(), false);
    chk("C5 失敗の理由が残る", !w3.error().empty(), true);
}

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : ".";
    std::printf("=== 記録基盤 2 本の限定検証（診断定義 ON）===\n");
    std::printf("KRI_DIAG_TRACE_ON = %d\n", KRI_DIAG_TRACE_ON);
    chk("C0 診断定義 ON で記録が有効", KRI_DIAG_TRACE_ON == 1, true);

    diag::TraceBin bin;
    std::printf("--- C1 source 別の混在記録 ---\n");
    c1_per_source(bin);
    std::printf("--- C2 接触の隅対応 ---\n");
    c2_corner(bin);
    std::printf("--- C3 結合順と ID ---\n");
    c3_merge_order();
    std::printf("--- C4 全体の上限超過 ---\n");
    c4_budget();
    std::printf("--- C5 排他的作成と書出し失敗 ---\n");
    c5_io(dir);

    std::printf("通過 %d / 失敗 %d\n", g_ok, g_ng);
    std::printf("判定: %s\n", g_ng == 0 ? "限定検証は全件通過" : "**失敗あり**");
    return g_ng == 0 ? 0 : 1;
}

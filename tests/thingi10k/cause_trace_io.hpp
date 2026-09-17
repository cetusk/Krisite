// Krisite — 原因究明用の実行記録の書出し（診断ビルド専用）
//
// `SPEC-phase5.md` §5.10.14.148 §3・§4。
//
// **本体（`include/krisite/`）には入れません。** ファイル I/O をここに閉じ込めます。
//
// ## 書式
//
// - メタデータと判断記録は **版付き JSON Lines**（1 行 1 レコード。`rec` が種別）。
//   **付属の読戻しはこの診断書式専用**で、汎用の来歴管理基盤は作りません。
// - **精密整数・有理数は文字列**。**浮動小数点へ落としません。**
// - **幾何配列は既存のバイナリ書式**（`dump_boolean` の `write_soup` と同じ）。
//   **対応表は別ファイル**にします。
//
// ## 失敗時
//
// **open / write / flush / close をすべて確認**し、1 つでも失敗したら
// **診断未完了**として偽を返します。**自動の再試行・切捨てをしません。**
#ifndef KRISITE_TESTS_CAUSE_TRACE_IO_HPP
#define KRISITE_TESTS_CAUSE_TRACE_IO_HPP

#if defined(KRISITE_DIAG_CAUSE_TRACE)

#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "krisite/csg/diagnostic_trace.hpp"

namespace krithingi::trace_io {

namespace diag = krisite::csg::diag;

/// この診断書式の版。**読戻し側と突き合わせます。**
inline constexpr int kTraceFormatVersion = 1;

// --- 文字列の組み立て（JSON の最小限。**依存を足しません**）------------------
/// **可逆な JSON エスケープ。** **制御文字を空白へ潰しません**
/// （タブ・CR を含む往復が成り立ちます）。
inline std::string esc(const std::string& s) {
    static const char* kHex = "0123456789abcdef";
    std::string o;
    o.reserve(s.size() + 2);
    for (char c : s) {
        const unsigned char u = static_cast<unsigned char>(c);
        switch (c) {
            case '"':  o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\r': o += "\\r"; break;
            case '\t': o += "\\t"; break;
            case '\b': o += "\\b"; break;
            case '\f': o += "\\f"; break;
            default:
                if (u < 0x20) {            // **残りの制御文字は \u00XX で保ちます**
                    o += "\\u00";
                    o += kHex[(u >> 4) & 0xF];
                    o += kHex[u & 0xF];
                } else {
                    o += c;                // UTF-8 のバイトはそのまま
                }
        }
    }
    return o;
}
inline std::string q(const std::string& s) { return "\"" + esc(s) + "\""; }
inline std::string kv(const char* k, const std::string& v) { return q(k) + ":" + v; }
inline std::string num(long long v) { return std::to_string(v); }
inline std::string num(unsigned long long v) { return std::to_string(v); }
inline std::string id(diag::SourceId x) { return x.valid() ? num((unsigned long long)x.v) : "null"; }
template <class T>
inline std::string idv(T x) { return x.valid() ? num((unsigned long long)x.v) : "null"; }
inline std::string boolean(bool b) { return b ? "true" : "false"; }

/// **セル鍵**（Morton を切り詰めません）。
inline std::string cellkey(const diag::CellKey& c) {
    return "{" + kv("depth", num((long long)c.depth)) + "," + kv("i", num((unsigned long long)c.i))
           + "," + kv("j", num((unsigned long long)c.j)) + "," + kv("k", num((unsigned long long)c.k)) + "}";
}
/// **結合順に依存しない参照。** 局所番号だけでは同一性になりません。
inline std::string ref(const diag::Ref& r) {
    return "{" + kv("op", num((long long)r.op)) + "," + kv("stage", num((long long)r.stage)) + ","
           + kv("cell", cellkey(r.cell)) + ","
           + kv("local", r.valid() ? num((unsigned long long)r.local) : std::string("null")) + "}";
}

template <class T, class F>
inline std::string arr(const std::vector<T>& v, F f) {
    std::string o = "[";
    for (std::size_t i = 0; i < v.size(); ++i) { if (i) o += ","; o += f(v[i]); }
    return o + "]";
}

inline const char* op_name(diag::Op o) {
    switch (o) {
        case diag::Op::kUnion: return "union";
        case diag::Op::kIsect: return "isect";
        case diag::Op::kDiffAB: return "diff_ab";
        case diag::Op::kDiffBA: return "diff_ba";
    }
    return "**未定義**";
}
inline const char* stage_name(diag::Stage s) {
    switch (s) {
        case diag::Stage::kPrepare: return "prepare";
        case diag::Stage::kFromMesh: return "from_mesh";
        case diag::Stage::kAssign: return "assign";
        case diag::Stage::kArrange: return "arrange";
        case diag::Stage::kStitch: return "stitch";
        case diag::Stage::kClassify: return "classify";
        case diag::Stage::kEmit: return "emit";
        case diag::Stage::kToMesh: return "to_mesh";
        case diag::Stage::kSplitContacts: return "split_contacts";
        case diag::Stage::kRepair: return "repair";
        case diag::Stage::kVerify: return "verify";
    }
    return "**未定義**";
}
inline const char* state_name(diag::StageState s) {
    switch (s) {
        case diag::StageState::kNotStarted: return "未着手";
        case diag::StageState::kRunning: return "処理中";
        case diag::StageState::kFinished: return "処理終了";
        case diag::StageState::kSaved: return "保存確認済み";
    }
    return "**未定義**";
}
inline const char* decision_name(diag::Decision d) {
    switch (d) {
        case diag::Decision::kGenerated: return "生成";
        case diag::Decision::kExcluded: return "除外";
        case diag::Decision::kAdopted: return "採用";
        case diag::Decision::kDiscarded: return "破棄";
        case diag::Decision::kFlipped: return "反転";
        case diag::Decision::kMergedAway: return "まとめた";
        case diag::Decision::kEmptyClip: return "クリップで空";
    }
    return "**未定義**";
}
inline const char* reason_name(diag::Reason r) {
    switch (r) {
        case diag::Reason::kUnset: return "未設定";
        case diag::Reason::kIndicatorFalse: return "指示関数が偽";
        case diag::Reason::kIndicatorBothSame: return "表裏で指示関数が同じ";
        case diag::Reason::kEarlyOutKnown: return "early-out で確定";
        case diag::Reason::kOutsideCell: return "セルに存在しない";
        case diag::Reason::kClippedEmpty: return "クリップの結果が空";
        case diag::Reason::kDuplicateCoplanar: return "共平面の重複";
        case diag::Reason::kCollinearFace: return "共線（from_mesh）";
        case diag::Reason::kSameIndexAfterQuant: return "量子化で同一索引";
        case diag::Reason::kDegenerateArea: return "面積 0";
        case diag::Reason::kForcedFromCell: return "供給セルから強制";
        case diag::Reason::kWindingSplit: return "巻き数の分裂";
        case diag::Reason::kRepairAccepted: return "修復を採用";
        case diag::Reason::kRepairRejected: return "修復を却下";
        case diag::Reason::kContactResolved: return "接触を解いた";
        case diag::Reason::kContactUnresolved: return "接触が未解決";
    }
    return "**未定義**";
}

// --- 書出し ----------------------------------------------------------------

/// **1 つの記録束を 1 ファイルへ。** 戻り値が偽なら**診断未完了**です。
class TraceWriter {
public:
    /// **排他的な新規作成**（`"wbx"`。C11 / POSIX）。
    ///
    /// **`rb` で存在を見てから `wb` で開く形は採りません** — 確認と作成の間に競合が残り、
    /// **読めない既存ファイル**にも作成してしまいます。
    bool open(const std::string& path, const std::string& run_id, const std::string& pair_id) {
        f_ = std::fopen(path.c_str(), "wbx");
        if (!f_) { fail("排他的に作成できません（既存か作成不能）: " + path); return false; }
        return line("header", kv("format", num((long long)kTraceFormatVersion)) + ","
                             + kv("run", q(run_id)) + "," + kv("pair", q(pair_id)));
    }

    bool write(const diag::TraceBin& b) {
        bool ok = true;
        for (const auto& x : b.stages)
            ok = ok && line("stage", kv("stage", q(stage_name(x.stage))) + ","
                            + kv("state", q(state_name(x.state))) + ","
                            + kv("op", q(op_name(x.op))) + "," + kv("seq", num((unsigned long long)x.seq)));
        for (const auto& x : b.removed)
            ok = ok && line("removed_face", kv("source", id(x.source)) + ","
                            + kv("src_face", idv(x.src_face)) + ","
                            + kv("reason", q(reason_name(x.reason))) + ","
                            + kv("where", q(x.where)));
        for (const auto& x : b.cells)
            ok = ok && line("cell", kv("cell", cellkey(x.cell)) + "," + kv("op", q(op_name(x.op))) + ","
                            + kv("depth", num((long long)x.depth)) + ","
                            + kv("both_sources", boolean(x.both_sources_present)) + ","
                            + kv("early_out", boolean(x.early_out_taken)) + ","
                            + kv("early_out_source", id(x.early_out_source)) + ","
                            + kv("early_out_winding_known", boolean(x.early_out_winding_known)) + ","
                            + kv("early_out_winding", num((long long)x.early_out_winding)) + ","
                            + kv("polys", num((unsigned long long)x.polys_assigned)) + ","
                            + kv("frags", num((unsigned long long)x.frags_after_arrange)));
        for (const auto& x : b.frags)
            ok = ok && line("frag", kv("frag", ref(x.frag)) + "," + kv("op", q(op_name(x.op))) + ","
                            + kv("stage", q(stage_name(x.stage))) + "," + kv("cell", cellkey(x.cell)) + ","
                            + kv("source", id(x.source)) + "," + kv("src_face", idv(x.src_face)) + ","
                            + kv("in_poly", idv(x.in_poly)) + "," + kv("out_poly", idv(x.out_poly)) + ","
                            + kv("parent", ref(x.parent)) + ","
                            + kv("decision", q(decision_name(x.decision))) + ","
                            + kv("reason", q(reason_name(x.reason))) + ","
                            + kv("flipped_before", boolean(x.flipped_before)) + ","
                            + kv("flipped_after", boolean(x.flipped_after)));
        for (const auto& x : b.regions)
            ok = ok && line("region", region_body(x));
        for (const auto& x : b.exits)
            ok = ok && line("exit", kv("op", q(op_name(x.op))) + "," + kv("stage", q(stage_name(x.stage)))
                            + "," + kv("out_face", idv(x.out_face)) + ","
                            + kv("tri_poly", idv(x.tri_poly)) + "," + kv("tri_src", id(x.tri_src)) + ","
                            + kv("tri_tag", num((unsigned long long)x.tri_tag)) + ","
                            + kv("vertex_planes",
                                 arr(x.vertex_planes, [](const diag::ExactValue& v) { return q(v.text); }))
                            + "," + kv("from_edge_split", boolean(x.from_edge_split)) + ","
                            + kv("before", idv(x.before)) + ","
                            + kv("decision", q(decision_name(x.decision))) + ","
                            + kv("reason", q(reason_name(x.reason))));
        for (const auto& x : b.edge_contacts)
            ok = ok && line("edge_contact", kv("op", q(op_name(x.op))) + ","
                            + kv("edge", "[" + num((unsigned long long)x.edge_u) + ","
                                 + num((unsigned long long)x.edge_v) + "]") + ","
                            + kv("branches", arr(x.branches, [](diag::OutFaceId i) { return idv(i); })) + ","
                            + kv("pairs", arr(x.pairs, [](const std::pair<std::uint32_t, std::uint32_t>& p) {
                                   return "[" + num((unsigned long long)p.first) + ","
                                          + num((unsigned long long)p.second) + "]"; })) + ","
                            + kv("resolved", boolean(x.resolved)) + ","
                            + kv("reason", q(reason_name(x.reason))));
        for (const auto& x : b.vertex_splits)
            ok = ok && line("vertex_split", kv("op", q(op_name(x.op))) + ","
                            + kv("old_vertex", num((unsigned long long)x.old_vertex)) + ","
                            + kv("corners", arr(x.corners, [](const diag::CornerRec& c) {
                                   return "{" + kv("face", idv(c.face)) + ","
                                          + kv("local_corner", num((long long)c.local_corner)) + ","
                                          + kv("class", num((unsigned long long)c.klass)) + ","
                                          + kv("new_vertex", c.new_vertex == diag::kNoId
                                                 ? std::string("null")
                                                 : num((unsigned long long)c.new_vertex)) + "}"; })) + ","
                            + kv("class_count", num((unsigned long long)x.class_count)) + ","
                            + kv("resolved", boolean(x.resolved)) + ","
                            + kv("reason", q(reason_name(x.reason))));
        for (const auto& x : b.verifies)
            ok = ok && line("verify", kv("op", q(op_name(x.op))) + "," + kv("name", q(x.name)) + ","
                            + kv("target", q(x.target)) + "," + kv("premise", q(x.premise)) + ","
                            + kv("premise_basis", q(x.premise_basis)) + ","
                            + kv("executed", boolean(x.executed)) + ","
                            + kv("fire_count", num((unsigned long long)x.fire_count)) + ","
                            + kv("result", q(x.result.text)) + ","
                            + kv("reached_final", boolean(x.reached_final)));
        // **落とした件数を必ず残します。** 非零なら診断未完了です。
        // **★ 書出しが失敗していれば、この末尾行も書けているとは限りません。**
        // **読戻し側は「末尾の `dropped` 行が無い」ことを未完了と判定します。**
        ok = line("dropped", kv("count", num((unsigned long long)b.dropped))) && ok;
        return ok && !failed_;
    }

    /// **flush と close を確認**してから真を返します。
    ///
    /// **先行する失敗を、`close` の成功が取り消しません**（`failed_` は固定）。
    bool close() {
        if (!f_) { fail("開いていません"); return false; }
        if (std::fflush(f_) != 0) fail("flush に失敗しました");
        if (std::fclose(f_) != 0) fail("close に失敗しました");
        f_ = nullptr;
        return !failed_;
    }

    /// **一度でも失敗したか**（固定された状態）。
    bool failed() const noexcept { return failed_; }
    const std::string& error() const noexcept { return err_; }

private:
    static const char* path_name(diag::WindingPath p) {
        switch (p) {
            case diag::WindingPath::kNotEvaluated: return "未実施";
            case diag::WindingPath::kDirect: return "direct";
            case diag::WindingPath::kForced: return "forced";
            case diag::WindingPath::kWindingSplit: return "winding_split";
        }
        return "**未定義**";
    }

    static std::string per_source_body(const diag::RegionSourceRec& s) {
        return "{" + kv("source", id(s.source)) + "," + kv("path", q(path_name(s.path))) + ","
               + kv("evaluated", boolean(s.evaluated)) + ","
               + kv("w_front", num((long long)s.w_front)) + ","
               + kv("w_back", num((long long)s.w_back)) + ","
               + kv("forced_from", cellkey(s.forced_from)) + ","
               + kv("w_other", num((long long)s.w_other)) + ","
               + kv("c_front", num((long long)s.c_front)) + ","
               + kv("c_back", num((long long)s.c_back)) + ","
               + kv("formula", q(s.formula)) + "}";
    }

    std::string region_body(const diag::RegionRec& x) const {
        return kv("region", ref(x.region)) + "," + kv("op", q(op_name(x.op))) + ","
               + kv("cell", cellkey(x.cell)) + ","
               + kv("frags", arr(x.frags, [](const diag::Ref& r) { return ref(r); })) + ","
               + kv("representative", ref(x.representative)) + ","
               + kv("support_plane", q(x.support_plane.text)) + ","
               + kv("used_point", boolean(x.used_point)) + "," + kv("point", q(x.point.text)) + ","
               + kv("per_source", arr(x.per_source, per_source_body)) + ","
               + kv("indicator_front", boolean(x.indicator_front)) + ","
               + kv("indicator_back", boolean(x.indicator_back)) + ","
               + kv("decision", q(decision_name(x.decision))) + ","
               + kv("reason", q(reason_name(x.reason))) + "," + kv("out_poly", idv(x.out_poly));
    }

    bool line(const char* rec, const std::string& body) {
        if (!f_) { fail("開いていません"); return false; }
        const std::string s = "{" + kv("rec", q(rec)) + "," + body + "}\n";
        if (std::fwrite(s.data(), 1, s.size(), f_) != s.size()) {
            fail(std::string("書出しに失敗しました: ") + rec);
            return false;
        }
        return true;
    }

    /// **最初の理由を残し、以後は上書きしません。**
    void fail(const std::string& why) {
        failed_ = true;
        if (err_.empty()) err_ = why;
    }

    FILE* f_ = nullptr;
    std::string err_;
    bool failed_ = false;
};

}  // namespace krithingi::trace_io

#endif  // KRISITE_DIAG_CAUSE_TRACE
#endif  // KRISITE_TESTS_CAUSE_TRACE_IO_HPP

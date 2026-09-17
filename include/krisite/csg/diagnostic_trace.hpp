// Krisite — 原因究明用の実行記録（診断ビルド専用）
//
// `SPEC-phase5.md` §5.10.14.148。
//
// ## この文書が満たす制約
//
// - **診断ビルドでだけ有効**（`KRISITE_DIAG_CAUSE_TRACE`）。
//   **通常ビルドでは型も呼出しも消え、Boolean の数値経路・既定値・幾何判断を変えません。**
// - **ファイル I/O と GMP への依存を本体へ入れません。** ここは**器だけ**で、
//   書出しは `tests/thingi10k/cause_trace_io.hpp` が行います。
// - **記録のための再計算をしません。** 判断に使った値をその場で受け取ります。
// - **共有ファイルへワーカーが直接書きません。** 仕事単位ごとの `TraceBin` に貯め、
//   既存の並列処理の**結合後**に `merge` でまとめます。
//
// ## 識別子を分ける（§4）
//
// 実行 / 対 / 演算 / 段 / source / 元面 / セル / 断片 / 出力面 を**別の型**で持ちます。
// **スレッド番号・完了順・アドレスを永続的な同一性に使いません。**
// 内部 ID には、必ず「幾何・元出現へ辿る対応」を併記します。
#ifndef KRISITE_CSG_DIAGNOSTIC_TRACE_HPP
#define KRISITE_CSG_DIAGNOSTIC_TRACE_HPP

#include <cstdint>

#if defined(KRISITE_DIAG_CAUSE_TRACE)

#include <atomic>
#include <cstddef>
#include <iterator>   // std::make_move_iterator（**間接 include に依存しません**）
#include <string>
#include <utility>
#include <vector>

namespace krisite::csg::diag {

/// **未設定を表す番人。** 符号のある量に負を使わないため、専用の最大値を置きます
/// （`CLAUDE.md`「符号のある量に、負を番人の値として使わないでください」）。
inline constexpr std::uint32_t kNoId = 0xFFFFFFFFu;

// ---------------------------------------------------------------------------
// 識別子。**型を分けます**（`CLAUDE.md`「位置引数が多く型が同じ小関数は、
// 規律ではなく設計で守ってください」）。
// ---------------------------------------------------------------------------
#define KRI_DIAG_ID_TYPE(Name)                                                      \
    struct Name {                                                                   \
        std::uint32_t v = kNoId;                                                    \
        constexpr bool valid() const noexcept { return v != kNoId; }                \
        friend constexpr bool operator==(Name a, Name b) noexcept { return a.v == b.v; } \
    }

KRI_DIAG_ID_TYPE(SourceId);   ///< `PolySoup::sources` の添字
KRI_DIAG_ID_TYPE(SrcFaceId);  ///< **元の入力メッシュ**の面番号
KRI_DIAG_ID_TYPE(PolyId);     ///< `PolySoup::polys` の添字
KRI_DIAG_ID_TYPE(OutFaceId);  ///< 出口の三角形番号
#undef KRI_DIAG_ID_TYPE

/// **セルの鍵。** **Morton 符号を 32 bit へ切り詰めません**（衝突するため）。
/// **深度と 3 軸の添字**で持ち、**結合順に依存しません**。
struct CellKey {
    std::uint8_t depth = 0;
    std::uint32_t i = 0, j = 0, k = 0;
    bool valid() const noexcept { return !(depth == 0 && i == 0 && j == 0 && k == 0); }
    friend bool operator==(const CellKey& a, const CellKey& b) noexcept {
        return a.depth == b.depth && a.i == b.i && a.j == b.j && a.k == b.k;
    }
};

/// **仕事単位をまたいで一意な参照。**
///
/// **局所 ID（仕事単位の中のスロット番号）は、単独では同一性になりません。**
/// **演算・段・セル鍵・仕事内番号**を含めて 1 つの参照とし、
/// **結合の順序に依存しない**形にします（§5.10.14.149 §2）。
struct Ref {
    std::uint8_t op = 0;      ///< `Op` の値
    std::uint8_t stage = 0;   ///< `Stage` の値
    CellKey cell{};
    std::uint32_t local = kNoId;  ///< 仕事内番号（**単独では同一性ではありません**）
    bool valid() const noexcept { return local != kNoId; }
    friend bool operator==(const Ref& a, const Ref& b) noexcept {
        return a.op == b.op && a.stage == b.stage && a.cell == b.cell && a.local == b.local;
    }
};

/// 演算。**左右の実引数は別に記録します**（`B\A` で source 番号を戻せるように）。
enum class Op : std::uint8_t { kUnion = 0, kIsect = 1, kDiffAB = 2, kDiffBA = 3 };

/// 段。**段の開始と終了を保存**し、**未到達を「違反 0」にしません**。
enum class Stage : std::uint8_t {
    kPrepare = 0,     ///< 変換・量子化
    kFromMesh,        ///< 入口（スープの生成）
    kAssign,          ///< セルへの割当て
    kArrange,         ///< 配置生成
    kStitch,          ///< 縫合
    kClassify,        ///< 分類
    kEmit,            ///< 出力多角形の採否
    kToMesh,          ///< 三角形化
    kSplitContacts,   ///< 接触分裂
    kRepair,          ///< 修復
    kVerify,          ///< 検算
};

/// 段の状態。**「未着手 / 処理中 / 処理終了 / 保存確認済み」を区別します**（§4）。
enum class StageState : std::uint8_t { kNotStarted = 0, kRunning, kFinished, kSaved };

/// 断片・領域に対して行った判断。**破棄側も残します**（§3 D）。
enum class Decision : std::uint8_t {
    kGenerated = 0,   ///< 生成した
    kExcluded,        ///< その段で除外した（分類へ来なかった）
    kAdopted,         ///< 採用した
    kDiscarded,       ///< 破棄した
    kFlipped,         ///< 向きを反転した
    kMergedAway,      ///< 重複としてまとめた
    kEmptyClip,       ///< クリップで空になった
};

/// 判断の理由。**「なぜ」を文字列の自由記述にしません**（集計と照合ができるように）。
enum class Reason : std::uint8_t {
    kUnset = 0,
    kIndicatorFalse,      ///< 指示関数が偽
    kIndicatorBothSame,   ///< 表裏で指示関数が同じ（現行は先に return する枝）
    kEarlyOutKnown,       ///< early-out で既知の巻き数から確定
    kOutsideCell,         ///< セルに存在しない
    kClippedEmpty,        ///< クリップの結果が空
    kDuplicateCoplanar,   ///< 共平面の重複
    kCollinearFace,       ///< 共線（`from_mesh` の除去）
    kSameIndexAfterQuant, ///< 量子化で同一索引になった（**`kCollinearFace` と混ぜません**）
    kDegenerateArea,      ///< 面積 0
    kForcedFromCell,      ///< 供給セルから強制された値
    kWindingSplit,        ///< 巻き数の分裂で分けた
    kRepairAccepted,      ///< 修復を採用
    kRepairRejected,      ///< 修復を却下
    kContactResolved,     ///< 接触を解いた
    kContactUnresolved,   ///< 接触が未解決
};

// ---------------------------------------------------------------------------
// 記録の単位
// ---------------------------------------------------------------------------

/// **厳密な値**。多倍長・有理数は**文字列**で持ちます（§4）。
/// **浮動小数点へ落としません。**
struct ExactValue {
    std::string text;  ///< 十進表記。有理数は `分子/分母`
};

/// 段の出入り。
struct StageMark {
    Stage stage{};
    StageState state{};
    Op op{};
    std::uint64_t seq = 0;  ///< **段の中の順序**。スレッド番号ではありません
};

/// `from_mesh` が共線として除いた元面（§3 C）。
/// **量子化で同一索引になって落ちた面とは、`reason` で区別します。**
struct RemovedFaceRec {
    SourceId source{};
    SrcFaceId src_face{};
    Reason reason = Reason::kUnset;
    std::string where;  ///< 判定位置（ファイル名と行。呼出し側が渡します）
};

/// セルの判断（§3 D）。**early-out の採否と、使った既知巻き数・source を残します。**
struct CellRec {
    CellKey cell{};
    Op op{};
    std::uint8_t depth = 0;
    bool both_sources_present = false;
    bool early_out_taken = false;
    SourceId early_out_source{};      ///< early-out に使った source
    std::int32_t early_out_winding = 0;
    bool early_out_winding_known = false;  ///< **0 を「未確定」と混ぜません**
    std::uint32_t polys_assigned = 0;
    std::uint32_t frags_after_arrange = 0;
};

/// 断片の 1 判断（§3 D）。**捨てたものにも元の識別子と理由を残します。**
struct FragRec {
    Ref frag{};             ///< **結合順に依存しない参照**
    Op op{};
    Stage stage{};
    CellKey cell{};
    SourceId source{};
    SrcFaceId src_face{};   ///< 元入力の面（**元出現へ辿る対応**）
    PolyId in_poly{};       ///< 入力多角形
    PolyId out_poly{};      ///< 出力多角形（採用時のみ有効）
    Ref parent{};           ///< 分割の親
    Decision decision{};
    Reason reason = Reason::kUnset;
    bool flipped_before = false;
    bool flipped_after = false;
};

/// **source ごとの巻き数の決まり方。**
///
/// **実処理は source ごとに `forced` か `winding_split` かを選びます**
/// （`soup_boolean.hpp:1972` 以降）。**領域に 1 組では混在も各 source の寄与も残りません。**
enum class WindingPath : std::uint8_t {
    kNotEvaluated = 0,  ///< **未実施**（0 を「値 0」と混ぜません）
    kDirect,            ///< その場のレイキャスト
    kForced,            ///< 供給セルから強制
    kWindingSplit,      ///< 巻き数の分裂で分けた
};

/// 領域 × source。**各 source の計算経路・供給セル・寄与を保持します。**
struct RegionSourceRec {
    SourceId source{};
    WindingPath path = WindingPath::kNotEvaluated;
    bool evaluated = false;          ///< **実施状態**（未実施と値 0 を分けます）
    std::int32_t w_front = 0, w_back = 0;
    CellKey forced_from{};           ///< `kForced` のときの供給セル
    std::int32_t w_other = 0;        ///< `kWindingSplit` のとき
    std::int32_t c_front = 0, c_back = 0;
    /// **計算式の接続**。どの値からこの `w_front` / `w_back` を得たかの短い式。
    /// **再計算ではなく、判断に使った値の写しです。**
    std::string formula;
};

/// 分類の領域（§3 D）。**代表点を使わなかった場合は `used_point` を偽にします。**
struct RegionRec {
    Ref region{};
    Op op{};
    CellKey cell{};
    std::vector<Ref> frags;            ///< 構成断片
    Ref representative{};              ///< 選んだ代表
    ExactValue support_plane;          ///< 支持平面の係数（厳密）
    bool used_point = false;           ///< 代表点を使ったか
    ExactValue point;                  ///< 実際の代表点（同次。使わないなら空）
    std::vector<RegionSourceRec> per_source;  ///< **source ごと**（上の理由）
    bool indicator_front = false, indicator_back = false;
    Decision decision{};
    Reason reason = Reason::kUnset;
    PolyId out_poly{};
};

/// 出口（§3 E）。**三角形化 → 接触分裂 → 修復の前後**を、元索引で結びます。
struct ExitRec {
    Op op{};
    Stage stage{};                     ///< kToMesh / kSplitContacts / kRepair
    OutFaceId out_face{};              ///< **元索引**
    PolyId tri_poly{};                 ///< 由来の多角形
    SourceId tri_src{};
    std::uint32_t tri_tag = 0;
    std::vector<ExactValue> vertex_planes;  ///< 頂点の構成平面（厳密）
    bool from_edge_split = false;      ///< 細分点が `edge_split` 由来か
    OutFaceId before{};                ///< 前の段での対応（複製・細分の前）
    Decision decision{};
    Reason reason = Reason::kUnset;
};

/// 接触のうち **辺**（§3 E）。**正常に解けた接触も残します。**
struct EdgeContactRec {
    Op op{};
    std::uint32_t edge_u = 0, edge_v = 0;   ///< **元索引**の辺
    std::vector<OutFaceId> branches;        ///< 過剰辺の巡回順（この順で保存）
    std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;  ///< 組（枝の添字対）
    bool resolved = false;
    Reason reason = Reason::kUnset;
};

/// 隅 1 つ。**`(元面 ID, 局所隅番号)` で特定します**（番号列だけでは復元できません）。
struct CornerRec {
    OutFaceId face{};            ///< 元面 ID
    std::uint8_t local_corner = 0;  ///< 局所隅番号（0/1/2）
    std::uint32_t klass = 0;     ///< 同値類
    std::uint32_t new_vertex = kNoId;  ///< **移った先の新頂点**
};

/// 接触のうち **頂点の隅分割**（§3 E）。**旧頂点 → 同値類 → 新頂点**を対応付けます。
struct VertexSplitRec {
    Op op{};
    std::uint32_t old_vertex = 0;   ///< **元索引**の頂点
    std::vector<CornerRec> corners; ///< その頂点に集まる隅の全部
    std::uint32_t class_count = 0;  ///< 同値類の数（= 分裂後の頂点数）
    bool resolved = false;
    Reason reason = Reason::kUnset;
};

/// 検算 1 件（§3 F）。
/// **「対象 / 必要前提 / 前提の根拠または未確認 / 実施 / 結果 / 最終判定への採用」を分けます。**
struct VerifyRec {
    Op op{};
    std::string name;            ///< 検算の名前（`verify_delta` など）
    std::string target;          ///< 対象
    std::string premise;         ///< 必要前提（PWN / NSI / 入力巻き数 0/1 / 出力位相 は別項目）
    std::string premise_basis;   ///< 根拠、または **`未確認`**
    bool executed = false;       ///< **設定と実施を分けます**
    std::uint32_t fire_count = 0;
    ExactValue result;
    bool reached_final = false;  ///< 最終判定へ採用されたか
};

// ---------------------------------------------------------------------------
// 仕事単位ごとの器と、結合
// ---------------------------------------------------------------------------

/// **1 仕事単位ぶんの記録。** 各ワーカーがこれを持ち、**結合後に 1 度だけ書き出します。**
struct TraceBin {
    std::vector<StageMark> stages;
    std::vector<RemovedFaceRec> removed;
    std::vector<CellRec> cells;
    std::vector<FragRec> frags;
    std::vector<RegionRec> regions;
    std::vector<ExitRec> exits;
    std::vector<EdgeContactRec> edge_contacts;
    std::vector<VertexSplitRec> vertex_splits;
    std::vector<VerifyRec> verifies;
    /// **落とした件数**。**自動の標本化・切捨てをしない**ので、
    /// **非零なら診断未完了**です（§4）。
    std::uint64_t dropped = 0;

    void merge(TraceBin&& o) {
        auto take = [](auto& dst, auto& src) {
            dst.insert(dst.end(), std::make_move_iterator(src.begin()),
                       std::make_move_iterator(src.end()));
            src.clear();
        };
        take(stages, o.stages);
        take(removed, o.removed);
        take(cells, o.cells);
        take(frags, o.frags);
        take(regions, o.regions);
        take(exits, o.exits);
        take(edge_contacts, o.edge_contacts);
        take(vertex_splits, o.vertex_splits);
        take(verifies, o.verifies);
        dropped += o.dropped;
    }

    bool empty() const noexcept {
        return stages.empty() && removed.empty() && cells.empty() && frags.empty()
               && regions.empty() && exits.empty() && edge_contacts.empty()
               && vertex_splits.empty() && verifies.empty();
    }
};

/// **全体の記録予算。** 仕事単位をまたいで**共有**します。
///
/// **各器が自分の上限内でも、全体では超え得ます**（§5.10.14.149 §2）。
/// そこで**予約（`reserve`）で数える**形にし、**超過は固定（sticky）**します。
/// **記録を落としながら対全体を完走する方式は採りません** —
/// 超過したら `exhausted()` が真になり、**呼出し側は安全な合流点で計算を終え、
/// 未完了を保存して後続を起動しません。**
///
/// **記録件数だけを、バイト容量・実メモリー上限の代用にしません** —
/// 容量の上限は駆動の必須引数で別に持ちます。
class TraceBudget {
public:
    explicit TraceBudget(std::uint64_t cap_records = 0) : cap_(cap_records) {}

    /// `n` 件ぶん予約します。**取れなければ超過を固定して偽**を返します。
    bool reserve(std::uint64_t n = 1) noexcept {
        if (cap_ == 0) return true;                 // 0 は「上限なし」
        std::uint64_t cur = used_.load(std::memory_order_relaxed);
        for (;;) {
            if (cur + n > cap_) {
                exhausted_.store(true, std::memory_order_relaxed);
                return false;
            }
            if (used_.compare_exchange_weak(cur, cur + n, std::memory_order_relaxed)) return true;
        }
    }
    bool exhausted() const noexcept { return exhausted_.load(std::memory_order_relaxed); }
    std::uint64_t used() const noexcept { return used_.load(std::memory_order_relaxed); }
    std::uint64_t cap() const noexcept { return cap_; }

private:
    std::uint64_t cap_ = 0;
    std::atomic<std::uint64_t> used_{0};
    std::atomic<bool> exhausted_{false};
};

/// 記録器。**呼出し側はポインタで持ち、`nullptr` なら何もしません。**
///
/// **上限は共有の `TraceBudget` が持ちます。** 器ごとの上限ではありません。
class CauseTrace {
public:
    CauseTrace() = default;
    explicit CauseTrace(TraceBudget* budget) : budget_(budget) {}

    TraceBin& bin() noexcept { return bin_; }
    const TraceBin& bin() const noexcept { return bin_; }
    TraceBudget* budget() const noexcept { return budget_; }

    /// **共有予算から 1 件予約**します。取れなければ `dropped` を数え、偽を返します。
    /// **偽が返ったら、呼出し側は記録を続けず、安全な合流点で停止を要求します。**
    bool room() noexcept {
        if (!budget_ || budget_->reserve(1)) return true;
        ++bin_.dropped;
        return false;
    }
    /// **落としたか**（この器）／**全体で超過したか**（共有）を分けます。
    bool dropped_here() const noexcept { return bin_.dropped != 0; }
    bool exhausted() const noexcept { return budget_ && budget_->exhausted(); }

    /// **結合後も、全体の超過状態を確かめられます。**
    void merge(CauseTrace&& o) { bin_.merge(std::move(o.bin_)); }

private:
    TraceBin bin_;
    TraceBudget* budget_ = nullptr;
};

}  // namespace krisite::csg::diag

// 呼出し側はこのマクロだけを使います。**通常ビルドでは完全に消えます。**
#define KRI_DIAG_TRACE_ON 1
#define KRI_DIAG_PUSH(tr, member, expr)                                   \
    do {                                                                  \
        if ((tr) && (tr)->room()) (tr)->bin().member.push_back(expr);     \
    } while (0)

#else  // KRISITE_DIAG_CAUSE_TRACE

namespace krisite::csg::diag {
/// **通常ビルドの空の器。** 呼出し側は `diag::CauseTrace*` を `nullptr` で渡します。
class CauseTrace;
}  // namespace krisite::csg::diag

#define KRI_DIAG_TRACE_ON 0
#define KRI_DIAG_PUSH(tr, member, expr) \
    do {                                \
    } while (0)

#endif  // KRISITE_DIAG_CAUSE_TRACE

#endif  // KRISITE_CSG_DIAGNOSTIC_TRACE_HPP

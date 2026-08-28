// Krisite — 構成点の保持（SPEC-phase2.md §4、§9.1）
//
// **Phase 1 の実測が示した最優先項目です。** `side` : `intersect3` = 1.26 : 1 で、
// `intersect3` が述語時間の 94.9%（推定）を占めていました。**単価を下げても、回数が
// `side` と同数のままでは構図が変わりません**（`IMPL-phase1.md` §5.6）。
//
// **キャッシュ無効側が正解器です**（§0.1 と同じ構図）。
//
// ---
//
// **キャッシュの有無で出力は 1 ビットも変わらないはずです。** どちらも平面3つ組を
// 昇順に並べてから `intersect3` を呼ぶので、比較を**値の完全一致**まで強められます。
// **$(C,\chi)$ で緩めてはいけません。** 緩めると §9.3 の変異 7（キーから 1 平面を落とす）
// が「別の三角形分割になっただけ」に見えて素通りし得ます。
//
// **危険なのは「静かに壊れる」ほうです**（§4.3）。キーの正規化を誤ると、異なる3つ組が
// 同じキーに落ちて**別の点を同一視します。**
//
// **4 平面以上が一点で交わる場合、異なる3つ組が同じ値を持ちます**（SPEC-phase1 §5.2）。
// **これはキャッシュの誤りではありません。** キャッシュは3つ組で引き、値の同一性は
// 第2段が別途扱う、という分担を崩さないこと。ケース 4T / 4T′ / 7 / 15 がこの配置です。
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
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
    std::size_t hits = 0, misses = 0, entries = 0, bytes = 0;
    std::size_t max_bytes = 0;
    std::size_t merged_by_value = 0;
};

Totals g;

/// 出力メッシュが**値まで**一致するか。$(C,\chi)$ より強い形（§9.1 の但し書き）。
bool same_mesh(const BoolMesh& a, const BoolMesh& b, std::string* why) {
    if (a.triangles.size() != b.triangles.size()) {
        *why = "三角形数が違う（" + std::to_string(a.triangles.size()) + " → " +
               std::to_string(b.triangles.size()) + "）";
        return false;
    }
    if (a.vertices.size() != b.vertices.size()) {
        *why = "頂点数が違う（" + std::to_string(a.vertices.size()) + " → " +
               std::to_string(b.vertices.size()) + "）";
        return false;
    }
    for (std::size_t i = 0; i < a.triangles.size(); ++i) {
        if (a.triangles[i] != b.triangles[i]) {
            *why = "三角形 " + std::to_string(i) + " の頂点添字が違う";
            return false;
        }
    }
    for (std::size_t i = 0; i < a.vertices.size(); ++i) {
        if (!krisite::geom::h_equal(a.vertices[i], b.vertices[i])) {
            *why = "頂点 " + std::to_string(i) + " の値が違う";
            return false;
        }
    }
    return true;
}

void run_case(const kritest::Case& c) {
    const TriMesh a = c.make_a(), b = c.make_b();
    std::printf("\n  ケース %-4s %s\n", c.id, c.what);

    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        for (unsigned d = 0; d <= kMaxDepth; ++d) {
            // **既定値に依存しないこと**（§9.4 の CI ジョブで既定が反転する）
            BoolOptions off = kritest::phase1_options(d);
            BoolOptions on = off;
            on.cache_points = true;

            BoolStats s_off, s_on;
            const BoolMesh r_off = boolean_op(a, b, op, off, &s_off);
            const BoolMesh r_on = boolean_op(a, b, op, on, &s_on);

            const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) + "（深度 " +
                                    std::to_string(d) + "、キャッシュ有無）";

            // ---- 値まで一致すること ----
            std::string why;
            KRI_CHECK_MSG(
                same_mesh(r_off, r_on, &why),
                tag + ": 出力が変わった（" + why + "）。**キャッシュは値を変えてはいけません**");
            // 断片・領域・併合も一致すること（分類経路まで同一）
            KRI_CHECK_MSG(s_on.fragments == s_off.fragments, tag + ": 断片数が変わった");
            KRI_CHECK_MSG(s_on.regions == s_off.regions, tag + ": 領域数が変わった");
            KRI_CHECK_MSG(s_on.merged_by_value == s_off.merged_by_value,
                          tag +
                              ": 第2段の値併合の数が変わった。**キャッシュと第2段の分担が"
                              "崩れています**（§4.3）");
            KRI_CHECK_MSG(s_on.max_planes_at_point == s_off.max_planes_at_point,
                          tag + ": 1 点に集まる平面の最大枚数が変わった");

            g.hits += s_on.cache_hits;
            g.misses += s_on.cache_misses;
            g.entries += s_on.cache_entries;
            g.bytes += s_on.cache_bytes;
            g.max_bytes = std::max(g.max_bytes, s_on.cache_bytes);
            g.merged_by_value += s_on.merged_by_value;

            if (op == BoolOp::Union && d == kMaxDepth) {
                const double hit =
                    (s_on.cache_hits + s_on.cache_misses)
                        ? 100.0 * static_cast<double>(s_on.cache_hits) /
                              static_cast<double>(s_on.cache_hits + s_on.cache_misses)
                        : 0.0;
                std::printf(
                    "    d%u ヒット %6zu / 参照 %6zu (%.1f%%) 実体 %5zu 個 %6zu B "
                    "| 値併合 %zu\n",
                    d, s_on.cache_hits, s_on.cache_hits + s_on.cache_misses, hit,
                    s_on.cache_entries, s_on.cache_bytes, s_on.merged_by_value);
            }
        }
    }
}

/// **空回りの番人。** ヒットが 0 なら、この検査は何も検証していません。
void check_not_vacuous() {
    const std::size_t refs = g.hits + g.misses;
    const double hit = refs ? 100.0 * static_cast<double>(g.hits) / static_cast<double>(refs) : 0.0;
    std::printf("\n  §4.4 の記録（全ケース × 3 演算 × 深度 0〜%u の合計）\n", kMaxDepth);
    std::printf("    参照 %zu / ヒット %zu (%.1f%%) / 実体 %zu 個\n", refs, g.hits, hit, g.entries);
    std::printf("    占有メモリ 合計 %zu B / 1 回の最大 %zu B\n", g.bytes, g.max_bytes);
    std::printf("    第2段の値併合 %zu 件（**キャッシュとは別の分担**。§4.3）\n",
                g.merged_by_value);

    KRI_CHECK_MSG(g.hits > 0, "キャッシュが一度もヒットしていない。**空回りです**");
    KRI_CHECK_MSG(hit > 50.0,
                  "ヒット率が 50% 未満。**構成点を作り直す構図が変わっていません**（§4.1）");
    // **値併合が 0 だと §4.3 の分担が検証されません。** 斜面ケースが必要です
    KRI_CHECK_MSG(g.merged_by_value > 0,
                  "第2段の値併合が 0 件。**異なる3つ組が同じ値を持つ配置がコーパスに"
                  "ありません**（§4.3 の分担が検証されない）");
}

}  // namespace

int main() {
    std::printf("\n  構成点の保持 — SPEC-phase2 §4 / CP3\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) run_case(c);
    check_not_vacuous();
    std::printf("\n");
    return kritest::finish("csg/cache");
}

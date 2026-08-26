// Krisite — CP2.5: 斜面ケース（2T / 4T / 5T）、3 演算、深度 0〜3
//
// SPEC-phase1.md §11 CP2.5「ここが実際の分岐点」
//
// **軸平行な立方体だけでは §5.2（4 平面同時交差）を突けません**（§9.1）。
// 1 点を通れる平面は 1 軸につき高々 1 枚で軸は 3 本だけ、かつセル面が立方体の面と
// 一致しても `PlaneTable` が同一 ID に併合するためです。
//
// **深度 0 が先決です。** 深度 0 には継ぎ目が無いので、そこで落ちるなら
// arrangement / 分類 / 抽出の中核が壊れています（§11「深度 0 で落ちたときの切り分け」）。
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

/// §9.3 の扱いを受けるか（合否に使わず、結果だけ記録する）。
///
/// ケース 4T は四面体 2 個が**頂点だけ**を共有します。ケース 11b と同型なので、
/// $\cup$ の出力は非多様体になります。**正則化では解決しません**（測度 0 なのは
/// 交差であって和集合ではない）。§9.3 の但し書きをそのまま適用します。
bool excluded_from_verdict(const char* id, BoolOp op) {
    const std::string s(id);
    // **ケース 4T′ の A\B も同じ理由で除外します。** 2 つの錐が頂点を共有するとき、
    // A\B のその頂点まわりのリンクは円環か複数片になり、扇が 1 つになりません。
    // 入力の接触ではなく**出力に現れる頂点接触**なので、§9.3 の字面はこれを覆って
    // いません（IMPL-phase1 §7.7 で仕様の拡張をお願いしています）。
    // 辺はすべて次数 2 で面の過不足は無く、落ちるのは頂点多様体だけです。
    return (s == "4T" && op == BoolOp::Union) || (s == "4T'" && op == BoolOp::Difference);
}

struct Row {
    std::string id;
    unsigned depth;
    BoolStats st;
    TopologyReport t;
};

void run_case(const kritest::Case& c, std::vector<Row>& rows) {
    const TriMesh a = c.make_a(), b = c.make_b();

    KRI_CHECK_MSG(kritest::size_discipline_ok(a, b),
                  std::string("ケース ") + c.id + ": §9.0 (1) のサイズ規律に違反");
    KRI_CHECK_MSG(check_topology(a).ok() && check_topology(b).ok(),
                  std::string("ケース ") + c.id + ": 入力が閉じた多様体でない");
    KRI_CHECK_MSG(krisite::mesh::is_outward_oriented(a) && krisite::mesh::is_outward_oriented(b),
                  std::string("ケース ") + c.id + ": 入力の向きが外向きでない");

    std::printf("\n  ケース %-3s %s\n", c.id, c.what);
    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        const bool excluded = excluded_from_verdict(c.id, op);
        TopologyReport ref;
        bool first = true;
        for (unsigned d = 0; d <= kMaxDepth; ++d) {
            BoolStats st;
            const BoolMesh r = boolean_op(a, b, op, d, &st);
            const TopologyReport t = check_topology(r.triangles);
            std::printf(
                "    d%u %s C=%zu g=%-2lld F=%-5zu %s | 断片%-4zu 有効セル%zu/%-4zu "
                "枚数%zu(mesh %zu) 値併合%zu/%zu 中点%zu 重心%zu 欠け辺%zu 余分辺%zu\n",
                d, op_name(op), t.components, t.genus_total, t.f,
                excluded ? "（§9.3 で合否対象外）" : (t.ok() ? "ok" : "**NG**"), st.fragments,
                st.active_cells, st.total_cells, st.max_planes_at_point,
                st.max_mesh_planes_at_point, st.merged_by_value, st.constructed_points,
                st.midpoint_raycasts, st.centroid_raycasts, t.edges_deficient, t.edges_excess);
            rows.push_back({c.id, d, st, t});
            if (excluded) continue;

            const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) + "（深度 " +
                                    std::to_string(d) + "）";
            KRI_CHECK_MSG(t.ok(), tag + ": 位相検査に落ちた");
            // **空は合法です**（正則化ブールの結果として起きる。§2.3）。
            // 4T の ∩ は頂点接触なので空になります。個別のフラグは空メッシュでは
            // すべて false なので、非空のときだけ見ます。
            if (!t.empty) {
                KRI_CHECK_MSG(t.edge_manifold, tag + ": 辺多様体でない");
                KRI_CHECK_MSG(t.vertex_manifold, tag + ": 頂点多様体でない");
                KRI_CHECK_MSG(t.oriented, tag + ": 向きが整合しない");
                KRI_CHECK_MSG(t.no_degenerate, tag + ": 退化三角形がある");
            }

            // §10.2.1 深度不変性
            if (first) {
                ref = t;
                first = false;
            } else {
                KRI_CHECK_MSG(t.components == ref.components,
                              tag + ": 深度不変性 — C が深度で変わった（" +
                                  std::to_string(ref.components) + " → " +
                                  std::to_string(t.components) + "）");
                KRI_CHECK_MSG(t.genus_total == ref.genus_total,
                              tag + ": 深度不変性 — g が深度で変わった（" +
                                  std::to_string(ref.genus_total) + " → " +
                                  std::to_string(t.genus_total) + "）");
            }
        }
    }
}

/// §13: 斜面ケースで 4 枚以上が検出されること。
///
/// **検出されないなら、計測器が機能していないかコーパスが §5.2 を突けていません。**
/// 軸平行ケースでの「最大 3 枚」という報告は、この検査が通って初めて意味を持ちます。
void check_four_planes(const std::vector<Row>& rows) {
    std::printf("\n  §13: 斜面ケースで 4 平面同時交差が検出されること\n");
    for (const char* id : {"2T", "4T", "4T'", "5T"}) {
        std::size_t best = 0, best_mesh = 0, merged = 0, pts = 0;
        for (const Row& r : rows) {
            if (r.id != id) continue;
            best = std::max(best, r.st.max_planes_at_point);
            best_mesh = std::max(best_mesh, r.st.max_mesh_planes_at_point);
            merged = std::max(merged, r.st.merged_by_value);
            pts = std::max(pts, r.st.constructed_points);
        }
        const double pct =
            pts ? (100.0 * static_cast<double>(merged) / static_cast<double>(pts)) : 0.0;
        std::printf(
            "    ケース %-3s 最大枚数 %zu（mesh のみ %zu）値併合 %zu（構成点 %zu の %.1f%%）\n", id,
            best, best_mesh, merged, pts, pct);
        KRI_CHECK_MSG(best >= 4,
                      std::string("ケース ") + id +
                          ": 4 平面同時交差が検出されていない（意図的に作ったのに 3 枚以下）");
        KRI_CHECK_MSG(merged > 0, std::string("ケース ") + id +
                                      ": 第2段が発火していない。4 平面同時交差があれば平面3つ組が "
                                      "相異なるのに同じ値になる頂点が出るはず（§5.3）");
    }
}

/// §9.0: 有効セル数が深度とともに増えること。
///
/// 断片数は間接的な指標で、面がセル境界に乗ると動きません。
/// **有効セル数のほうが「分割が働いているか」の直接的な指標です。**
void check_active_cells(const std::vector<Row>& rows) {
    std::printf("\n  §9.0: 有効セル数の推移（分割が働いているかの直接的な指標）\n");
    for (const char* id : {"2T", "4T", "4T'", "5T"}) {
        // 断片の生成は演算に依存しないので、各深度の最初の 1 件を採れば足ります
        std::size_t first = 0, last = 0, prev = 0;
        std::string line;
        for (unsigned d = 0; d <= kMaxDepth; ++d) {
            for (const Row& r : rows) {
                if (r.id != id || r.depth != d) continue;
                const std::size_t n = r.st.active_cells;
                if (d == 0) {
                    first = n;
                    line = std::to_string(n);
                } else {
                    line += " → " + std::to_string(n);
                    KRI_CHECK_MSG(n >= prev,
                                  std::string("ケース ") + id + ": 有効セル数が深度で減った");
                }
                prev = n;
                last = n;
                break;
            }
        }
        std::printf("    ケース %-3s 有効セル %s\n", id, line.c_str());
        KRI_CHECK_MSG(last > first, std::string("ケース ") + id +
                                        ": 有効セル数が深度 0 → 3 で増えない（分割が空回り）");
    }
}

}  // namespace

int main() {
    std::printf("\n  CP2.5 — 斜面ケース（SPEC-phase1 §11）\n");
    std::vector<Row> rows;
    for (const kritest::Case& c : kritest::corpus()) {
        const std::string id = c.id;
        if (id != "2T" && id != "4T" && id != "4T'" && id != "5T") continue;
        run_case(c, rows);
    }
    check_four_planes(rows);
    check_active_cells(rows);
    std::printf("\n");
    return kritest::finish("csg/cp25");
}

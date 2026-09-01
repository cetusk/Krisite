// Krisite — NSI の宣言が意味論を変えないこと（SPEC-phase3 §5.6 / SPEC-phase5 (c)）
//
// **NSI は三角形分割を変えます。** 切断を省くので断片が減り、多角形が大きくなり、
// 三角形が変わります。**したがってバイトでは守れません**（`IMPL-phase5.md` §20）。
//
// ここで見るのは**位相**です。体積は `test_volume_gmp.cpp` が受け持ちます
// （GMP が要るので別建て。**あちらが無効な構成でも、この検査は走ります**）。
//
// **そして機構が空回りしていないことを検査します。**
//
// > コーパスのサイズ規律（`SPEC-phase1.md` §9.0 (1)）は「両入力が座標範囲の広い部分を
// > 占める」ことを要求するので、**ほとんどのセルに両方の source が居ます。**
// > **単一 source のセルがほぼ生じず、NSI が一度も発火しませんでした。**
//
// ケース 26（離して置き一部だけ重ねる）だけがこの配置を作ります。
// **26 が消えたら、この検査は空回りに戻ります。**
#include <cstdio>
#include <string>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kDepth = 2;

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

struct Run {
    TopologyReport t;
    std::size_t polys = 0;
    std::size_t skipped_cells = 0;
};

Run run_one(const TriMesh& a, const TriMesh& b, BoolOp op, bool nsi) {
    PolySoup sa = from_mesh(a), sb = from_mesh(b);
    if (nsi) {
        sa.nsi.assign(sa.sources.size(), 1);
        sb.nsi.assign(sb.sources.size(), 1);
    }
    BoolOptions o = kritest::corpus_options(kDepth);
    o.local_bsp = true;
    o.split_contacts = true;
    ToMeshOptions tm;
    tm.split_contacts = true;
    BoolStats st;
    const PolySoup soup = boolean(sa, sb, op, o, &st);
    const SoupMesh m = to_mesh(soup, tm);
    Run r;
    r.t = check_topology(m.triangles);
    r.polys = soup.polys.size();
    r.skipped_cells = st.bsp_cells_skipped_nsi;
    return r;
}

std::size_t total_skipped = 0;
std::size_t cases_where_nsi_changed = 0;

}  // namespace

int main() {
    std::printf("=== NSI の宣言（SPEC-phase3 §5.6）===\n");
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        std::size_t skipped = 0;
        bool changed = false;
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            const Run off = run_one(a, b, op, false);
            const Run on = run_one(a, b, op, true);
            const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) + ": ";

            // **宣言しない側で省いてはいけません**（既定は偽）
            KRI_CHECK_MSG(off.skipped_cells == 0, tag + "NSI を宣言していないのに省いた");
            // **位相が一致すること。** 三角形分割は変わってよい
            KRI_CHECK_MSG(off.t.components == on.t.components,
                          tag + "NSI の宣言で連結成分が変わった");
            KRI_CHECK_MSG(off.t.genus_total == on.t.genus_total,
                          tag + "NSI の宣言で種数が変わった");
            KRI_CHECK_MSG(off.t.edge_manifold == on.t.edge_manifold &&
                              off.t.vertex_manifold == on.t.vertex_manifold,
                          tag + "NSI の宣言で多様体性が変わった");
            KRI_CHECK_MSG(off.t.oriented == on.t.oriented, tag + "NSI の宣言で向きが変わった");
            skipped += on.skipped_cells;
            if (off.polys != on.polys) changed = true;
        }
        total_skipped += skipped;
        if (changed) ++cases_where_nsi_changed;
        std::printf("  ケース %-4s 省いたセル %-6zu %s %s\n", c.id, skipped,
                    changed ? "（多角形数が変わった）" : "", c.what);
    }

    // **空回りの検査。** コーパスのサイズ規律が単一 source のセルを排除するので、
    // **ケース 26 が無ければ 0 になります**（実際に 0 でした）。
    KRI_CHECK_MSG(total_skipped > 0,
                  "NSI が 1 セルも省いていない（機構が空回り。ケース 26 が消えた？）");
    KRI_CHECK_MSG(cases_where_nsi_changed > 0,
                  "NSI で多角形数が変わったケースが 1 つも無い（検査が空回り）");
    std::printf("\n  省いたセル 合計 %zu / 多角形数が変わったケース %zu\n", total_skipped,
                cases_where_nsi_changed);
    return kritest::finish("csg/nsi");
}

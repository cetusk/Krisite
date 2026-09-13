// Krisite — **解けずに残る次数 4 の辺を、細分して直す**（案 H）
//
// `DESIGN-phase5-vertex-level.md` §3.4 / §9。`SPEC-phase5.md` §5.10.14.59。
//
// ## この配置は機構から構成しました（実データの最小化ではありません）
//
// §3.4 が確定した機構は次のとおりです。
//
//   辺 e の 2 枚のシートは、両端の頂点まわりで【1 本の閉曲線】になっている（八の字）
//   曲線は節 w を 2 度通るので、剥がしても 1 本の円のまま
//   つまり頂点は割れず、索引付きメッシュでは辺を 2 本に分けられない
//
// **その形を作るには、2 つの楔が辺で接し、かつ【端点の近くで別の道で繋がって】
// いればよい**（馬蹄形の両端が接した形）。
//
//   A       x,y ≥ 0 の箱        z ∈ [0,10]
//   B       x,y ≤ 0 の箱        z ∈ [0,10]      → A と B は z 軸で接する
//   下の橋   両方に体積で重なる   z ∈ [-3,1]      → 下の端点のまわりで繋ぐ
//   上の橋   両方に体積で重なる   z ∈ [9,12]      → 上の端点のまわりで繋ぐ
//
// **残る接触辺は z ∈ [1,9] の 1 本**で、**両端とも扇が 1 個に潰れる**はずです。
//
// ## 検査は【対】で行います
//
// **修復を外した側で従来の失敗が再現しなければ、修復が効いていることを示せません**
// （`CLAUDE.md`「機構を追加したら、それを外す経路も用意してください」）。
#include <cstdint>
#include <cstdio>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"

#include "corpus.hpp"
#include "test_util.hpp"
#include "volume_fp.hpp"

using namespace krisite;
using krisite::mesh::TriMesh;

namespace {

int failures = 0;

void expect(bool ok, const char* what) {
    if (!ok) {
        ++failures;
        std::printf("  **FAIL** %s\n", what);
    }
}

/// 4 つの閉じた箱を 1 枚のメッシュにまとめる（**自己交差します。宣言しません**）。
TriMesh pinch_input(std::int32_t u) {
    TriMesh m;
    const auto append = [&m](const TriMesh& x) {
        const auto off = static_cast<mesh::VertexId>(m.vertices.size());
        for (const auto& v : x.vertices) m.vertices.push_back(v);
        for (const mesh::Tri& t : x.triangles) {
            m.triangles.push_back({static_cast<mesh::VertexId>(t[0] + off),
                                   static_cast<mesh::VertexId>(t[1] + off),
                                   static_cast<mesh::VertexId>(t[2] + off)});
        }
    };
    append(kritest::box(0, 0, 0, 5 * u, 5 * u, 10 * u));                // A
    append(kritest::box(-5 * u, -5 * u, 0, 0, 0, 10 * u));              // B（z 軸で A と接する）
    append(kritest::box(-2 * u, -2 * u, -3 * u, 2 * u, 2 * u, u));      // 下の橋
    append(kritest::box(-2 * u, -2 * u, 9 * u, 2 * u, 2 * u, 12 * u));  // 上の橋
    return m;
}

struct Run {
    std::size_t unresolved_post = 0;
    std::size_t before = 0, edges = 0, collisions = 0, no_planes = 0, no_pair = 0;
    bool edge_manifold = false, vertex_manifold = false;
    std::size_t triangles = 0;
    double vol = 0.0;
    unsigned long long hash = 0;
};

Run go(const TriMesh& m, unsigned depth, bool repair) {
    csg::PolySoup a = csg::from_mesh(m);
    csg::PolySoup b = csg::from_mesh(m);
    csg::BoolOptions o;
    o.depth = depth;
    o.adaptive = true;
    o.cull_planes = true;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    o.threads = 1;
    const csg::PolySoup u = csg::boolean(a, b, csg::BoolOp::Union, o);
    csg::ToMeshOptions tm;
    tm.split_contacts = true;
    tm.verify_split_delta = true;
    tm.repair_unresolved = repair;
    tm.threads = 1;
    csg::ToMeshStats ts;
    const csg::SoupMesh out = csg::to_mesh(u, tm, &ts);
    const mesh::TopologyReport tr = mesh::check_topology(out.triangles);
    Run r;
    r.unresolved_post = ts.split.unresolved_post;
    r.before = ts.split.unresolved_before_repair;
    r.edges = ts.split.repair_edges;
    r.collisions = ts.split.repair_collisions;
    r.no_planes = ts.split.repair_no_planes;
    r.no_pair = ts.split.repair_no_pair;
    r.edge_manifold = tr.edge_manifold;
    r.vertex_manifold = tr.vertex_manifold;
    r.triangles = out.triangles.size();
    r.vol = kritest::volume6_fp(out);
    return r;
}

}  // namespace

int main() {
    std::printf("## 接触辺の細分（案 H）\n\n");
    const std::int32_t u = 1 << 12;
    const TriMesh m = pinch_input(u);
    std::printf("入力: 三角形 %zu 枚（閉じた箱 4 つ。**自己交差します**）\n\n", m.triangles.size());

    std::printf("| 深度 | 修復 | 解けずに残った辺 | 細分した辺 | 衝突 | 平面不足 | 組不足 |");
    std::printf(" 辺多様体 | 頂点多様体 | 三角形 | 6 倍体積 |\n");
    std::printf("|---|---|---:|---:|---:|---:|---:|---|---|---:|---:|\n");

    int fired = 0;
    for (unsigned d : {1u, 2u, 3u}) {
        const Run off = go(m, d, false);
        const Run on = go(m, d, true);
        for (int k = 0; k < 2; ++k) {
            const Run& r = (k == 0) ? off : on;
            std::printf("| %u | %s | %zu | %zu | %zu | %zu | %zu | %s | %s | %zu | %.0f |\n", d,
                        k == 0 ? "外す" : "入れる", r.before, r.edges, r.collisions, r.no_planes,
                        r.no_pair, r.edge_manifold ? "はい" : "**いいえ**",
                        r.vertex_manifold ? "はい" : "**いいえ**", r.triangles, r.vol);
        }
        // **`unresolved_before_repair` は修復の経路でしか埋まりません**
        // （外した側は `unresolved_detail` を集めないため）。**発火の判定は
        // 外した側の `unresolved_post` で行います**
        if (off.unresolved_post > 0) {
            ++fired;
            // **外した側で従来の失敗が再現すること**
            expect(off.unresolved_post > 0, "修復を外すと解けずに残る");
            expect(!off.edge_manifold, "修復を外すと辺多様体でない");
            // **入れた側で直ること**
            expect(on.unresolved_post == 0, "修復を入れると解ける");
            expect(on.edge_manifold && on.vertex_manifold, "修復を入れると多様体になる");
            expect(on.edges == on.before, "残っていた辺をすべて細分した");
            expect(on.before > 0, "修復の前に残っていた辺が記録されている");
            expect(on.collisions == 0 && on.no_planes == 0 && on.no_pair == 0, "諦めた辺が無い");
            // **細分は辺 1 本につき三角形を 4 枚増やします**
            expect(on.triangles == off.triangles + 4 * on.edges, "三角形が 4×細分数だけ増えた");
            // **体積は 1 ビットも変わってはいけません**（細分点は辺の上にある）
            expect(on.vol == off.vol, "**体積が変わっていない**");
        }
    }

    // **空回りの検査**: 1 つの深度でも発火しなければ、この構成は狙いを外しています
    if (fired == 0) {
        std::printf("\n**どの深度でも接触辺が残りませんでした。構成が狙いを外しています**\n");
        ++failures;
    } else {
        std::printf("\n**発火した深度 %d / 3**\n", fired);
    }
    std::printf("\n**不一致 %d 件**\n", failures);
    return failures == 0 ? 0 : 1;
}

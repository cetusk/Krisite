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
#include <string>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"

#include "corpus.hpp"
#include "test_util.hpp"
#include "volume_fp.hpp"

#if defined(KRISITE_TEST_EXACT_VOLUME)
// **厳密な体積は GMP で見ます**（`volume6_fp` は倍精度で、篩にしかなりません）。
// **本体は GMP に依存しません** — このテストのターゲットだけがリンクします
#include "volume_gmp.hpp"
#endif

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
    std::size_t realloc = 0, probes = 0;
    bool edge_manifold = false, vertex_manifold = false;
    std::size_t triangles = 0;
    double vol = 0.0;
    std::string exact_vol;   ///< **GMP のときだけ**。有理数の文字列
    bool restore_ok = true;  ///< 由来から点を復元して一致したか
    std::size_t restored = 0;
};

std::size_t exact_volume_checks = 0;

#if defined(KRISITE_TEST_EXACT_VOLUME)
/// **整数同次座標から有理数の体積を作る**（`volume_gmp.hpp`）。文字列で返して比べます。
std::string exact_volume6(const csg::SoupMesh& m) {
    mpq_t v;
    mpq_init(v);
    kritest::mesh_volume6(v, m);
    char* p = mpq_get_str(nullptr, 10, v);
    std::string r(p);
    void (*fr)(void*, std::size_t) = nullptr;
    mp_get_memory_functions(nullptr, nullptr, &fr);
    fr(p, r.size() + 1);
    mpq_clear(v);
    ++exact_volume_checks;
    return r;
}
#endif

/// **保存した由来から点を作り直し、保存された頂点と厳密に一致するか**（§12.3.3）。
///
/// **符号の一元化だけでは、保存や復元の誤りは防げません。**
bool check_restore(const csg::PolySoup& s, const csg::SoupMesh& m, std::size_t* n) {
    bool ok = true;
    for (std::size_t k = 0; k < m.edge_split.size(); ++k) {
        const csg::EdgeSplitSource& e = m.edge_split[k];
        const geom::PlaneD& p1 = s.table.at(e.p1);
        const geom::PlaneD& p2 = s.table.at(e.p2);
        const geom::PlaneD& p3 = s.table.at(e.p3);
        const geom::PlaneD& p4 = s.table.at(e.p4);
        // 1. 符号が、保存された両端から計算し直した値と一致するか
        const int want = -geom::side(p4, m.vertices[e.v]) * geom::side(p3, m.vertices[e.w]);
        if (want != static_cast<int>(e.sign)) ok = false;
        // 2. 4 枚と符号から作り直した点が、保存された頂点と厳密に一致するか
        const geom::PlaneSum q = geom::plane_sum(p3, p4, e.sign);
        const geom::HPointD again = geom::intersect3(p1, p2, q);
        // 3. 旧 API との一致
        const geom::HPointD old_api =
            geom::edge_interior_point(p1, p2, p3, p4, m.vertices[e.v], m.vertices[e.w]);
        if (!geom::h_equal(again, old_api)) ok = false;
        bool found = false;
        for (std::size_t i = 0; i < m.vertex_split_src.size(); ++i) {
            if (m.vertex_split_src[i] != k) continue;
            found = true;
            if (!geom::h_equal(m.vertices[i], again)) ok = false;
        }
        if (!found) ok = false;
        ++*n;
    }
    return ok;
}

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
    r.realloc = ts.split.repair_vertices_realloc;
    r.probes = ts.split.repair_collide_probes;
    r.restore_ok = check_restore(u, out, &r.restored);
#if defined(KRISITE_TEST_EXACT_VOLUME)
    r.exact_vol = exact_volume6(out);
#endif
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

#if defined(KRISITE_TEST_EXACT_VOLUME)
    std::printf("**体積の比較: GMP の有理数（厳密）**\n\n");
#else
    std::printf("**体積の比較: `volume6_fp` の倍精度（篩）**\n\n");
#endif
    std::printf("| 深度 | 修復 | 残った辺 | 細分 | 衝突 | 平面不足 | 組不足 | 再確保 | 照合 |");
    std::printf(" 辺多様体 | 頂点多様体 | 三角形 | 復元 |\n");
    std::printf("|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---|\n");

    int fired = 0;
    for (unsigned d : {1u, 2u, 3u}) {
        const Run off = go(m, d, false);
        const Run on = go(m, d, true);
        for (int k = 0; k < 2; ++k) {
            const Run& r = (k == 0) ? off : on;
            std::printf(
                "| %u | %s | %zu | %zu | %zu | %zu | %zu | %zu | %zu | %s | %s | %zu |"
                " %s（%zu 件）|\n",
                d, k == 0 ? "外す" : "入れる", r.before, r.edges, r.collisions, r.no_planes,
                r.no_pair, r.realloc, r.probes, r.edge_manifold ? "はい" : "**いいえ**",
                r.vertex_manifold ? "はい" : "**いいえ**", r.triangles,
                r.restore_ok ? "一致" : "**不一致**", r.restored);
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
            // **体積は変わってはいけません**（細分点は辺の上にある）。
            // **倍精度の等値は篩です**（`SPEC-phase5.md` §5.10.14.74）
            expect(on.vol == off.vol, "体積（倍精度の篩）が変わっていない");
#if defined(KRISITE_TEST_EXACT_VOLUME)
            expect(on.exact_vol == off.exact_vol, "**厳密な有理数の体積が変わっていない**");
#endif
            // **保存した由来から点を復元できること**（§12.3.3）
            expect(on.restore_ok, "由来から復元した点が、保存された頂点と厳密に一致する");
            expect(on.restored == on.edges, "細分した辺の数だけ由来が保存されている");
            expect(off.restored == 0, "修復を外すと由来は残らない");
            // **衝突の照合が実際に走ったこと**
            expect(on.probes > 0, "衝突の照合が走った");
#if defined(KRISITE_TEST_EXPECT_REALLOC)
            // **再確保の発火**（`KRISITE_MUTATION_REPAIR_TIGHT_CAPACITY` のとき）
            expect(on.realloc >= 1, "**再確保が発火した**");
#endif
        }
    }

    // **空回りの検査**: 1 つの深度でも発火しなければ、この構成は狙いを外しています
    if (fired == 0) {
        std::printf("\n**どの深度でも接触辺が残りませんでした。構成が狙いを外しています**\n");
        ++failures;
    } else {
        std::printf("\n**発火した深度 %d / 3**\n", fired);
    }
#if defined(KRISITE_TEST_EXACT_VOLUME)
    // **厳密比較が実際に走ったことを確かめます**（走っていなければ検査していないのと同じ）
    if (exact_volume_checks == 0) {
        std::printf("\n**GMP の厳密比較が 1 度も走っていません**\n");
        ++failures;
    } else {
        std::printf("\n**厳密な体積を %zu 回比べました**\n", exact_volume_checks);
    }
#endif
    std::printf("\n**不一致 %d 件**\n", failures);
    return failures == 0 ? 0 : 1;
}

// Krisite — 位相検査器そのもののテスト
//
// SPEC-phase1.md §10.1
//
// **検査器は出力より先に用意します。** そしてこの検査器自身が正しいことを、
// 位相が既知の構成で確かめます。ここが甘いと CP1 の合否判定が信用できません。
//
// とくに重要なのは「壊れた構成をちゃんと落とすこと」です。§10.5 の変異テストは
// この検査器の検出力に乗っているので、緩い検査だと T 字接合が素通りします。
#include <vector>

#include "krisite/mesh/topology.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::mesh;
using kritest::box;
using kritest::combinatorial_torus;
using kritest::concat;
using kritest::cube;
using kritest::flipped;
using kritest::tetra;

namespace {

void expect_closed(const char* what, const TriMesh& m, std::size_t comps, long long genus) {
    const TopologyReport r = check_topology(m);
    KRI_CHECK_MSG(r.ok(), std::string(what) + ": ok() が false");
    KRI_CHECK_MSG(r.edge_manifold, std::string(what) + ": 辺多様体でない");
    KRI_CHECK_MSG(r.vertex_manifold, std::string(what) + ": 頂点多様体でない");
    KRI_CHECK_MSG(r.oriented, std::string(what) + ": 向きが整合しない");
    KRI_CHECK_MSG(r.no_degenerate, std::string(what) + ": 退化三角形がある");
    KRI_CHECK_MSG(r.components == comps, std::string(what) +
                                             ": C = " + std::to_string(r.components) + " 期待 " +
                                             std::to_string(comps));
    KRI_CHECK_MSG(r.genus_total == genus, std::string(what) +
                                              ": g_total = " + std::to_string(r.genus_total) +
                                              " 期待 " + std::to_string(genus));
    // χ = 2(C - g_total)
    KRI_CHECK_MSG(r.chi == 2 * (static_cast<long long>(comps) - genus),
                  std::string(what) + ": χ = " + std::to_string(r.chi));
}

// ---- 正しい構成 --------------------------------------------------------------

void test_valid_shapes() {
    // 立方体: V=8, E=18, F=12 → χ=2, C=1, g=0
    {
        const TriMesh c = cube(0, 0, 0, 100);
        const TopologyReport r = check_topology(c);
        KRI_CHECK(r.v == 8 && r.e == 18 && r.f == 12);
        expect_closed("立方体", c, 1, 0);
    }
    // 四面体: V=4, E=6, F=4 → χ=2
    {
        const TriMesh t = tetra(100);
        const TopologyReport r = check_topology(t);
        KRI_CHECK(r.v == 4 && r.e == 6 && r.f == 4);
        expect_closed("四面体", t, 1, 0);
    }
    // 離れた 2 立方体: C=2, χ=4, g_total=0
    expect_closed("離れた 2 立方体", concat(cube(0, 0, 0, 10), cube(100, 100, 100, 10)), 2, 0);

    // 空洞のある立方体（ケース 10 の A\B の形）: 外殻 + 内殻（内向き）
    // C=2, χ=4, g_total=0
    expect_closed("空洞つき立方体", concat(cube(0, 0, 0, 100), flipped(cube(30, 30, 30, 20))), 2,
                  0);

    // 組合せトーラス: χ=0, C=1, g=1
    {
        const TriMesh t = combinatorial_torus(4, 5);
        const TopologyReport r = check_topology(t);
        KRI_CHECK(r.v == 20 && r.e == 60 && r.f == 40);
        expect_closed("トーラス", t, 1, 1);
    }
    // 2 つのトーラス: C=2, g_total=2, χ=0
    expect_closed("トーラス 2 個", concat(combinatorial_torus(4, 5), combinatorial_torus(3, 6)), 2,
                  2);
}

// ---- 空メッシュ --------------------------------------------------------------

void test_empty() {
    const TriMesh e;
    const TopologyReport r = check_topology(e);
    KRI_CHECK(r.empty);
    KRI_CHECK(r.v == 0 && r.e == 0 && r.f == 0);
    KRI_CHECK(r.ok());  // 空は合法（正則化ブールの結果として起きる。§2.3）
}

// ---- 壊れた構成を落とすこと --------------------------------------------------

void test_detects_breakage() {
    // 1 枚抜く → 穴があき、辺が 1 面にしか接しない
    {
        TriMesh c = cube(0, 0, 0, 100);
        c.triangles.pop_back();
        const TopologyReport r = check_topology(c);
        KRI_CHECK(!r.edge_manifold);
        KRI_CHECK(!r.ok());
    }
    // 1 枚だけ向きを裏返す → 向きが整合しない
    {
        TriMesh c = cube(0, 0, 0, 100);
        const auto tmp = c.triangles[3][1];
        c.triangles[3][1] = c.triangles[3][2];
        c.triangles[3][2] = tmp;
        const TopologyReport r = check_topology(c);
        KRI_CHECK(!r.oriented);
        KRI_CHECK(!r.ok());
        KRI_CHECK(r.edge_manifold);  // 辺の本数は変わらない。向きだけが壊れる
    }
    // 同じ三角形を 2 度入れる → 辺が 4 面に接する
    {
        TriMesh c = cube(0, 0, 0, 100);
        c.triangles.push_back(c.triangles[0]);
        const TopologyReport r = check_topology(c);
        KRI_CHECK(!r.edge_manifold);
        KRI_CHECK(!r.ok());
    }
    // 退化三角形
    {
        TriMesh c = cube(0, 0, 0, 100);
        c.triangles.push_back({1, 1, 2});
        const TopologyReport r = check_topology(c);
        KRI_CHECK(!r.no_degenerate);
        KRI_CHECK(!r.ok());
    }
    // 頂点だけを共有する 2 立方体（ケース 11b の ∪）→ 頂点多様体でない
    {
        const TriMesh m = kritest::two_cubes_sharing_a_vertex(50);
        const TopologyReport r = check_topology(m);
        KRI_CHECK(r.edge_manifold);     // 辺は壊れていない
        KRI_CHECK(r.oriented);          // 向きも整合している
        KRI_CHECK(!r.vertex_manifold);  // **頂点まわりに扇が 2 つできる**
        KRI_CHECK(!r.ok());
    }
}

// ---- T 字接合を検出できること（§10.5 の変異 3 の前提）------------------------
//
// 幾何的に隙間が無くても、片側だけが辺を分割していれば辺多様体が破れます。
// この性質が成り立たないと、§4.3.1 の変異テストが素通りします。
void test_detects_t_junction() {
    // 正方形の面を、片側だけ 2 枚に分割した最小構成
    //   下の面: 1 枚の三角形が辺 (0,1) を持つ
    //   上の面: 辺 (0,m) と (m,1) の 2 枚（m は辺 (0,1) の中点として追加した頂点）
    TriMesh m;
    m.vertices = {{0, 0, 0}, {100, 0, 0}, {50, 100, 0}, {50, 0, 0}, {50, -100, 0}};
    //            0           1            2             3(中点)     4
    m.triangles = {
        {0, 1, 2},  // 分割されていない側
        {0, 3, 4},  // 分割された側（前半）
        {3, 1, 4},  // 分割された側（後半）
    };
    const TopologyReport r = check_topology(m);
    // 辺 (0,1) は 1 面にしか接しない。辺 (0,3) と (3,1) も 1 面ずつ
    KRI_CHECK(!r.edge_manifold);
    KRI_CHECK(!r.ok());
}

// ---- 深度不変性の比較関数（§10.2.1 が使う量）---------------------------------
//
// 位相の同一性は (C, g_total) の組で判定します。空メッシュどうしも一致とみなします。
void test_depth_invariant_key() {
    const TriMesh a = cube(0, 0, 0, 100);
    // 同じ立方体を「分割された」形で作る: 各面を 4 枚に割る代わりに、
    // ここでは検査器が分割に鈍感であることだけ確かめる（頂点数は変わってよい）
    const TopologyReport ra = check_topology(a);
    const TopologyReport rb = check_topology(box(0, 0, 0, 100, 100, 100));
    KRI_CHECK(ra.components == rb.components);
    KRI_CHECK(ra.genus_total == rb.genus_total);

    // 空どうし
    const TopologyReport e1 = check_topology(TriMesh{});
    const TopologyReport e2 = check_topology(TriMesh{});
    KRI_CHECK(e1.empty && e2.empty);
    KRI_CHECK(e1.components == e2.components && e1.genus_total == e2.genus_total);
}

// ---- 符号付き体積と向き検査 --------------------------------------------------

void test_signed_volume() {
    using krisite::arith::sign;
    const std::int32_t s = 100;
    const TriMesh c = cube(0, 0, 0, s);
    KRI_CHECK(is_outward_oriented(c));
    KRI_CHECK(coords_in_range(c));
    KRI_CHECK(sign(signed_volume6(c)) == 1);

    // 裏返すと負
    KRI_CHECK(!is_outward_oriented(flipped(c)));
    KRI_CHECK(sign(signed_volume6(flipped(c))) == -1);

    // 空洞つき: 外殻 100^3 − 内殻 20^3 > 0
    const TriMesh cav = concat(cube(0, 0, 0, 100), flipped(cube(30, 30, 30, 20)));
    KRI_CHECK(is_outward_oriented(cav));
    {
        using namespace krisite::arith;
        constexpr std::size_t LV = krisite::geom::limbs::kInputVolume6;
        const auto want = from_i64<LV>(6ll * (100ll * 100 * 100 - 20ll * 20 * 20));
        KRI_CHECK(cmp(signed_volume6(cav), want) == 0);
    }

    // 空メッシュは体積 0。向き検査は false（呼び出し側で空を区別すること）
    KRI_CHECK(sign(signed_volume6(TriMesh{})) == 0);
    KRI_CHECK(!is_outward_oriented(TriMesh{}));
}

}  // namespace

int main() {
    test_valid_shapes();
    test_empty();
    test_detects_breakage();
    test_detects_t_junction();
    test_depth_invariant_key();
    test_signed_volume();
    return kritest::finish("mesh/topology");
}

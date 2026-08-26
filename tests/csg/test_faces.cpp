// Krisite — 平面 ID 登録と面併合のテスト
//
// SPEC-phase1.md §3.1
//
// 中心的な検査は「頂点が 3 平面の交点として復元できること」です。
// 併合が正しければ、面の頂点 i は support ∩ edge[i-1] ∩ edge[i] で
// 元の入力頂点にぴったり一致するはずです。
#include <set>
#include <vector>

#include "krisite/csg/faces.hpp"
#include "krisite/mesh/topology.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using namespace krisite::geom;
using krisite::mesh::TriMesh;
using kritest::cube;
using kritest::tetra;

namespace {

// ---- PlaneTable --------------------------------------------------------------

void test_plane_table() {
    PlaneTable t;
    const PlaneD px = plane_axis_aligned(Axis::X, 100);
    const PlaneD px2 = plane_axis_aligned(Axis::X, 100);
    const PlaneD py = plane_axis_aligned(Axis::Y, 100);

    const PlaneRef a = t.intern(px);
    KRI_CHECK(a.id == 0 && !a.flipped);
    const PlaneRef b = t.intern(px2);
    KRI_CHECK(b.id == 0 && !b.flipped);  // 同じ平面は同じ ID
    const PlaneRef c = t.intern(py);
    KRI_CHECK(c.id == 1);
    KRI_CHECK(t.size() == 2);

    // 向きを裏返した平面は同じ ID で flipped = true
    PlaneD neg = px;
    neg.a = krisite::arith::neg(px.a);
    neg.b = krisite::arith::neg(px.b);
    neg.c = krisite::arith::neg(px.c);
    neg.d = krisite::arith::neg(px.d);
    const PlaneRef d = t.intern(neg);
    KRI_CHECK(d.id == 0);
    KRI_CHECK_MSG(d.flipped, "符号反転した平面は flipped = true になるはず");
    KRI_CHECK(t.size() == 2);

    // find は登録しない
    KRI_CHECK(t.find(px).id == 0);
    KRI_CHECK(t.find(plane_axis_aligned(Axis::Z, 7)).id == kNoPlane);
    KRI_CHECK(t.size() == 2);
}

/// 多数の平面を入れても ID と全順序が壊れないこと。
void test_plane_table_bulk() {
    PlaneTable t;
    kritest::Rng rng(31);
    std::vector<PlaneD> ps;
    for (int i = 0; i < 300; ++i) {
        PlaneD pl;
        do {
            pl = plane_from_triangle(kritest::rand_small_point(rng, 64),
                                     kritest::rand_small_point(rng, 64),
                                     kritest::rand_small_point(rng, 64));
        } while (is_degenerate(pl));
        ps.push_back(pl);
        t.intern(pl);
    }
    // すべて再登録しても size は増えない
    const std::size_t n = t.size();
    for (const PlaneD& pl : ps) t.intern(pl);
    KRI_CHECK(t.size() == n);
    // 登録した各平面が、同一と判定される代表に写ること
    for (const PlaneD& pl : ps) {
        const PlaneRef r = t.find(pl);
        KRI_CHECK(r.id != kNoPlane);
        KRI_CHECK(plane_same(t.at(r.id), pl));
    }
}

// ---- 面併合 ------------------------------------------------------------------

/// 面の頂点が「support ∩ edge[i-1] ∩ edge[i]」で元の入力頂点に一致すること。
/// **これが併合の正しさの中心的な証拠**です。
void check_vertices_reconstruct(const char* what, const TriMesh& m, const std::vector<Face>& faces,
                                const PlaneTable& table) {
    for (const Face& f : faces) {
        const std::size_t n = f.edge.size();
        KRI_CHECK(n == f.loop.size());
        KRI_CHECK(n >= 3);
        for (std::size_t i = 0; i < n; ++i) {
            const PlaneD& s = table.at(f.support);
            const PlaneD& e0 = table.at(f.edge[(i + n - 1) % n]);
            const PlaneD& e1 = table.at(f.edge[i]);
            if (!kritest::intersects_at_point(s, e0, e1)) {
                KRI_CHECK_MSG(false, std::string(what) + ": 3 平面が一点で交わらない");
                continue;
            }
            const HPointD v = intersect3(s, e0, e1);
            const HPointD want = to_homogeneous(m.vertices[f.loop[i]]);
            KRI_CHECK_MSG(h_equal(v, want),
                          std::string(what) + ": 頂点 " + std::to_string(i) + " が復元できない");
        }
    }
}

void test_cube_faces() {
    PlaneTable table;
    const TriMesh m = cube(-100, -100, -100, 200);
    KRI_CHECK(check_topology(m).ok());
    const std::vector<Face> faces = build_faces(m, 0, table);

    // 立方体: 12 三角形 → 6 面、平面は 6 枚
    KRI_CHECK_MSG(faces.size() == 6,
                  "立方体は 6 面に併合されるはず（実際 " + std::to_string(faces.size()) + "）");
    KRI_CHECK_MSG(table.size() == 6,
                  "立方体の平面は 6 枚（実際 " + std::to_string(table.size()) + "）");
    // 各面は 4 辺
    for (const Face& f : faces) KRI_CHECK(f.edge.size() == 4);
    // 支持平面はすべて相異なる
    std::set<PlaneId> sup;
    for (const Face& f : faces) sup.insert(f.support);
    KRI_CHECK(sup.size() == 6);

    check_vertices_reconstruct("立方体", m, faces, table);
}

void test_tetra_faces() {
    PlaneTable table;
    const TriMesh m = tetra(200);
    KRI_CHECK(check_topology(m).ok());
    const std::vector<Face> faces = build_faces(m, 0, table);

    KRI_CHECK(faces.size() == 4);
    KRI_CHECK(table.size() == 4);
    for (const Face& f : faces) KRI_CHECK(f.edge.size() == 3);
    check_vertices_reconstruct("四面体", m, faces, table);
}

/// 2 つのメッシュを同じ表に入れると、共有する平面は同じ ID になる。
void test_shared_planes_across_meshes() {
    PlaneTable table;
    const TriMesh a = cube(0, 0, 0, 100);
    const TriMesh b = cube(0, 0, 100, 100);  // z = 100 の面で接する
    const std::vector<Face> fa = build_faces(a, 0, table);
    const std::vector<Face> fb = build_faces(b, 1, table);
    KRI_CHECK(fa.size() == 6 && fb.size() == 6);
    // 真上に積むと x と y の側面 4 枚も一致するので、共有は z=100 を含めて 5 枚。
    // 平面は 6 + 6 - 5 = 7 枚
    KRI_CHECK_MSG(table.size() == 7,
                  "積み重ねた立方体の平面は 7 枚（実際 " + std::to_string(table.size()) + "）");

    // z だけで接し、側面は一致しない配置なら共有は z=100 の 1 枚だけ → 11 枚
    PlaneTable t3;
    build_faces(cube(0, 0, 0, 100), 0, t3);
    build_faces(cube(37, 41, 100, 100), 1, t3);
    KRI_CHECK_MSG(t3.size() == 11,
                  "共有が z=100 の 1 枚だけなら 11 枚（実際 " + std::to_string(t3.size()) + "）");

    // 一般位置なら共有なし
    PlaneTable t2;
    build_faces(cube(0, 0, 0, 100), 0, t2);
    build_faces(cube(37, 41, 43, 100), 1, t2);
    KRI_CHECK(t2.size() == 12);
}

/// 面の向きフラグが、元の三角形の外向き法線と整合すること。
void test_orientation_flag() {
    PlaneTable table;
    const TriMesh m = cube(-50, -50, -50, 100);
    const std::vector<Face> faces = build_faces(m, 0, table);
    for (const Face& f : faces) {
        // 面の外向き法線 = flipped ? -N : N。原点は立方体の内部なので
        // 外向き法線に対して原点は負側にあるはず
        const PlaneD& pl = table.at(f.support);
        const int s = side(pl, IPoint{0, 0, 0});
        KRI_CHECK(s != 0);
        // flipped なら符号が反転する
        KRI_CHECK_MSG((f.flipped ? -s : s) == -1,
                      "立方体の内部（原点）は各面の外向き法線に対して負側のはず");
    }
}

}  // namespace

int main() {
    test_plane_table();
    test_plane_table_bulk();
    test_cube_faces();
    test_tetra_faces();
    test_shared_planes_across_meshes();
    test_orientation_flag();
    return kritest::finish("csg/faces");
}

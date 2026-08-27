// Krisite — T 頂点の解決（SPEC-phase2.md §2.4.3 / §2.4.4 (3) の【正の対照】）
//
// **負の対照だけでは、何もしないスタブが通ります**（§2.4.4 (3)）。
// 「固定深度で挿入 0 件」は、候補索引が常に空を返す実装でも `cmp_h` の判定が常に false に
// なる実装でも通ります。そこで**効くべき場所で効くこと**を、手で作った配置で確かめます。
//
// 負の対照（固定深度で挿入 0 件）は `tests/csg/test_split_strategy.cpp` が持ちます。
//
// ここで使う配置は単純です。z = 0 上の正方形 (0,0)-(4,4) で、辺 y = 0 の上に
// x = 1 と x = 3 の 2 点を置きます。**1 本の辺に 2 個**なので、§9.3 の変異 11
// （1 個ずつ挿入する）もここで効きます。
#include <array>
#include <cstdio>
#include <vector>

#include "krisite/csg/tjunction.hpp"

#include "test_util.hpp"

using krisite::csg::EdgeVertexIndex;
using krisite::csg::fan_triangulate;
using krisite::csg::insert_t_vertices;
using krisite::csg::PlaneId;
using krisite::csg::PlaneTable;
using krisite::csg::TJunctionStats;
using krisite::csg::TPolygon;
using krisite::geom::Axis;
using krisite::geom::HPointD;
using krisite::geom::plane_axis_aligned;

namespace {

/// 軸平行平面だけで作る足場。**構成点は必ず `intersect3` で作ります。**
/// 座標を手で書くと、被検体と同じ経路を通らない「別物」を検査してしまいます。
struct Scaffold {
    PlaneTable table;
    std::vector<HPointD> verts;
    EdgeVertexIndex index;

    PlaneId plane(Axis ax, std::int64_t c) { return table.intern(plane_axis_aligned(ax, c)).id; }

    /// 3 平面の交点を頂点として登録し、索引にも入れる。
    std::uint32_t vertex(PlaneId p, PlaneId q, PlaneId r) {
        const auto id = static_cast<std::uint32_t>(verts.size());
        verts.push_back(krisite::geom::intersect3(table.at(p), table.at(q), table.at(r)));
        std::array<PlaneId, 3> tri{p, q, r};
        std::sort(tri.begin(), tri.end());
        index.add_triple(tri, id);
        return id;
    }
};

/// z = 0 上の正方形 (0,0)-(4,4)。頂点 i = support ∩ edge[i-1] ∩ edge[i] の規約に合わせる。
struct Square {
    Scaffold s;
    PlaneId support{}, y0{}, x4{}, y4{}, x0{};
    std::vector<PlaneId> edge;
    std::vector<std::uint32_t> poly;

    Square() {
        support = s.plane(Axis::Z, 0);
        y0 = s.plane(Axis::Y, 0);
        x4 = s.plane(Axis::X, 4);
        y4 = s.plane(Axis::Y, 4);
        x0 = s.plane(Axis::X, 0);
        edge = {y0, x4, y4, x0};
        poly = {
            s.vertex(support, x0, y0),  // (0,0)
            s.vertex(support, y0, x4),  // (4,0)
            s.vertex(support, x4, y4),  // (4,4)
            s.vertex(support, y4, x0),  // (0,4)
        };
    }
};

// ---- 正の対照 --------------------------------------------------------------

void test_inserts_in_order() {
    Square sq;
    // 辺 0（y = 0、(0,0) → (4,0)）の上に 2 点。**順序を逆に登録**して、整列が効くことを見る
    const std::uint32_t t3 = sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, 3));
    const std::uint32_t t1 = sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, 1));

    TJunctionStats st;
    const TPolygon p = insert_t_vertices(sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);

    KRI_CHECK_MSG(st.inserted == 2,
                  "T 頂点が 2 個入るはず（実測 " + std::to_string(st.inserted) + "）");
    KRI_CHECK_MSG(st.max_per_edge == 2, "1 本の辺に 2 個載るはず。**変異 11 の前提**");
    KRI_CHECK(p.vertex.size() == 6);
    // (0,0) → t1 → t3 → (4,0) の順に並ぶこと
    KRI_CHECK_MSG(p.vertex[0] == sq.poly[0] && p.vertex[1] == t1 && p.vertex[2] == t3 &&
                      p.vertex[3] == sq.poly[1],
                  "線分に沿って整列していない（登録順ではなく座標順であること）");
    KRI_CHECK(p.is_corner[0] && !p.is_corner[1] && !p.is_corner[2] && p.is_corner[3]);
    KRI_CHECK(p.corners == 4);
    // 元の角は 4 つとも残る
    std::size_t corners = 0;
    for (char c : p.is_corner) corners += static_cast<std::size_t>(c);
    KRI_CHECK_MSG(corners == 4, "元の角が失われている");
}

/// **区間の外と両端は入れないこと。** ここを緩めると T 字接合を直すどころか作ります。
void test_rejects_outside_and_endpoints() {
    Square sq;
    sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, 5));   // 区間の外
    sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, -2));  // 反対側の外
    // 端点と同一の値を持つ別 ID（第2段が併合し損ねた場合を模す）
    sq.s.vertex(sq.support, sq.y0, sq.x4);

    TJunctionStats st;
    const TPolygon p = insert_t_vertices(sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);
    KRI_CHECK_MSG(st.inserted == 0, "区間の外・端点を入れてしまっている（実測 " +
                                        std::to_string(st.inserted) + " 個）");
    KRI_CHECK(p.vertex.size() == 4);
}

/// **別の辺の上の点を、この辺に入れないこと。** 索引は平面対で引くので構造的に起きない
/// はずですが、対の作り方を間違えると起きます。
void test_other_edge_not_mixed() {
    Square sq;
    const std::uint32_t on_x4 = sq.s.vertex(sq.support, sq.x4, sq.s.plane(Axis::Y, 2));

    TJunctionStats st;
    const TPolygon p = insert_t_vertices(sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);
    KRI_CHECK_MSG(st.inserted == 1, "辺 x=4 の上の点がちょうど 1 個入るはず");
    // 辺 1（x = 4、(4,0) → (4,4)）の内部に入ること
    KRI_CHECK(p.vertex.size() == 5);
    KRI_CHECK_MSG(p.vertex[2] == on_x4, "入る位置が違う（辺 1 の内部のはず）");
    KRI_CHECK(!p.is_corner[2]);
}

// ---- 扇分割（§2.4.4 (2)）---------------------------------------------------

void test_fan_drops_only_degenerate() {
    Square sq;
    sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, 1));
    sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, 3));

    TJunctionStats st;
    const TPolygon p = insert_t_vertices(sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);
    std::vector<std::array<std::uint32_t, 3>> tris;
    fan_triangulate(p, tris, &st);

    // 頂点 6 個の凸多角形だが、起点の隣接辺に T 頂点が 2 個あるので退化が 2 枚出る。
    // **面積は四角形のままなので、正しい三角形は 2 枚。**
    KRI_CHECK_MSG(tris.size() == 2,
                  "三角形が 2 枚のはず（実測 " + std::to_string(tris.size()) + " 枚）");
    KRI_CHECK_MSG(st.degenerate_dropped == 2,
                  "退化を 2 枚捨てるはず（実測 " + std::to_string(st.degenerate_dropped) + "）");
    // 起点は元の角であること
    for (const auto& t : tris) KRI_CHECK(t[0] == sq.poly[0]);
    // 退化していない = 同じ頂点を 2 度使っていない
    for (const auto& t : tris) KRI_CHECK(t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
}

/// **T 頂点が無ければ扇分割は Phase 1 と同一であること。** 負の対照の土台です。
void test_fan_unchanged_without_t() {
    Square sq;
    TJunctionStats st;
    const TPolygon p = insert_t_vertices(sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);
    std::vector<std::array<std::uint32_t, 3>> tris;
    fan_triangulate(p, tris, &st);

    KRI_CHECK(st.inserted == 0 && st.degenerate_dropped == 0);
    KRI_CHECK(tris.size() == 2);
    KRI_CHECK(tris[0][0] == sq.poly[0] && tris[0][1] == sq.poly[1] && tris[0][2] == sq.poly[2]);
    KRI_CHECK(tris[1][0] == sq.poly[0] && tris[1][1] == sq.poly[2] && tris[1][2] == sq.poly[3]);
}

/// 索引の絞り込みが効いていること。**全頂点を返すなら索引の意味がありません。**
void test_index_narrows() {
    Square sq;
    for (int i = 1; i <= 3; ++i) sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, i));
    for (int i = 1; i <= 3; ++i) sq.s.vertex(sq.support, sq.x4, sq.s.plane(Axis::Y, i));

    const auto* on_y0 = sq.s.index.find(sq.support, sq.y0);
    KRI_CHECK(on_y0 != nullptr);
    // 辺 y=0 の線上: 角 2 つ + 追加 3 点 = 5
    KRI_CHECK_MSG(on_y0->size() == 5,
                  "候補数が想定と違う（実測 " + std::to_string(on_y0->size()) + "）");
    KRI_CHECK_MSG(on_y0->size() < sq.s.verts.size(), "索引が全頂点を返している（絞り込み無効）");
}

}  // namespace

int main() {
    std::printf("\n  T 頂点の解決 — SPEC-phase2 §2.4.3 / §2.4.4 (3) の正の対照\n");
    test_inserts_in_order();
    test_rejects_outside_and_endpoints();
    test_other_edge_not_mixed();
    test_fan_drops_only_degenerate();
    test_fan_unchanged_without_t();
    test_index_narrows();
    std::printf("\n");
    return kritest::finish("csg/tjunction");
}

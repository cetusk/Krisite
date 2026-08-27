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

using krisite::csg::fan_triangulate;
using krisite::csg::insert_t_vertices;
using krisite::csg::PlaneId;
using krisite::csg::PlaneTable;
using krisite::csg::PlaneVertexIndex;
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
    PlaneVertexIndex index;
    std::vector<PlaneId> planes;

    PlaneId plane(Axis ax, std::int64_t c) {
        const PlaneId id = table.intern(plane_axis_aligned(ax, c)).id;
        if (std::find(planes.begin(), planes.end(), id) == planes.end()) planes.push_back(id);
        return id;
    }

    /// 3 平面の交点を頂点として登録する。**座標を手で書きません。**
    /// 被検体と同じ経路（`intersect3`）を通らないと、別物を検査することになります。
    std::uint32_t vertex(PlaneId p, PlaneId q, PlaneId r) {
        const auto id = static_cast<std::uint32_t>(verts.size());
        verts.push_back(krisite::geom::intersect3(table.at(p), table.at(q), table.at(r)));
        return id;
    }

    /// **頂点を全部足してから呼ぶこと。** 索引は `side` で張るので、後から足した頂点は
    /// 入りません。
    void build() { index.build(table, verts, planes); }
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

    sq.s.build();
    TJunctionStats st;
    const TPolygon p =
        insert_t_vertices(sq.s.table, sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);

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

    sq.s.build();
    TJunctionStats st;
    const TPolygon p =
        insert_t_vertices(sq.s.table, sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);
    KRI_CHECK_MSG(st.inserted == 0, "区間の外・端点を入れてしまっている（実測 " +
                                        std::to_string(st.inserted) + " 個）");
    KRI_CHECK(p.vertex.size() == 4);
}

/// **別の辺の上の点を、この辺に入れないこと。** 索引は平面対で引くので構造的に起きない
/// はずですが、対の作り方を間違えると起きます。
void test_other_edge_not_mixed() {
    Square sq;
    const std::uint32_t on_x4 = sq.s.vertex(sq.support, sq.x4, sq.s.plane(Axis::Y, 2));

    sq.s.build();
    TJunctionStats st;
    const TPolygon p =
        insert_t_vertices(sq.s.table, sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);
    KRI_CHECK_MSG(st.inserted == 1, "辺 x=4 の上の点がちょうど 1 個入るはず");
    // 辺 1（x = 4、(4,0) → (4,4)）の内部に入ること
    KRI_CHECK(p.vertex.size() == 5);
    KRI_CHECK_MSG(p.vertex[2] == on_x4, "入る位置が違う（辺 1 の内部のはず）");
    KRI_CHECK(!p.is_corner[2]);
}

// ---- 扇分割（§2.4.4 (2)）---------------------------------------------------

/// **退化三角形は捨てるのではなく作らない**（§2.4.4 (2)）。起点を選べば 1 枚も出ません。
void test_fan_avoids_degenerate() {
    Square sq;
    sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, 1));
    sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, 3));

    sq.s.build();
    TJunctionStats st;
    const TPolygon p =
        insert_t_vertices(sq.s.table, sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);
    std::vector<std::array<std::uint32_t, 3>> tris;
    fan_triangulate(p, tris, &st);

    // T 頂点は辺 0 だけに載るので、**元の角 2**（隣接辺は 1 と 2）が起点に選ばれる
    KRI_CHECK_MSG(tris.size() == p.vertex.size() - 2,
                  "扇は n-2 枚そろうこと（実測 " + std::to_string(tris.size()) + " 枚 / n=" +
                      std::to_string(p.vertex.size()) + "）。**1 枚でも欠けると境界辺が消える**");
    KRI_CHECK_MSG(st.degenerate_kept == 0, "起点を選べば退化は出ないはず（実測 " +
                                               std::to_string(st.degenerate_kept) + " 枚）");
    KRI_CHECK_MSG(st.apex_fallback == 0, "起点は選べるはず");
    for (const auto& t : tris) KRI_CHECK_MSG(t[0] == sq.poly[2], "起点が元の角 2 でない");
    for (const auto& t : tris) KRI_CHECK(t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
}

/// **起点を選べない多角形では残します。** 捨てると境界辺が消えて次数 1 の辺が出ます。
void test_fan_fallback_keeps_degenerate() {
    Square sq;
    // 4 辺すべてに T 頂点を 1 個ずつ → どの角も隣接辺に T 頂点を持つ
    sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, 2));
    sq.s.vertex(sq.support, sq.x4, sq.s.plane(Axis::Y, 2));
    sq.s.vertex(sq.support, sq.y4, sq.s.plane(Axis::X, 2));
    sq.s.vertex(sq.support, sq.x0, sq.s.plane(Axis::Y, 2));

    sq.s.build();
    TJunctionStats st;
    const TPolygon p =
        insert_t_vertices(sq.s.table, sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);
    std::vector<std::array<std::uint32_t, 3>> tris;
    fan_triangulate(p, tris, &st);

    KRI_CHECK_MSG(st.inserted == 4,
                  "4 辺に 1 個ずつ入るはず（実測 " + std::to_string(st.inserted) + "）");
    KRI_CHECK_MSG(st.apex_fallback == 1, "起点を選べないはず（どの角も隣に T 頂点がある）");
    // **枚数は n-2 のまま。捨てないので境界辺は全部残ります**
    KRI_CHECK_MSG(tris.size() == p.vertex.size() - 2,
                  "退化を残しても扇は n-2 枚（実測 " + std::to_string(tris.size()) + " 枚）");
    KRI_CHECK_MSG(st.degenerate_kept == 2,
                  "残した退化は 2 枚のはず（実測 " + std::to_string(st.degenerate_kept) + "）");
}

/// **T 頂点が無ければ扇分割は Phase 1 と同一であること。** 負の対照の土台です。
void test_fan_unchanged_without_t() {
    Square sq;
    sq.s.build();
    TJunctionStats st;
    const TPolygon p =
        insert_t_vertices(sq.s.table, sq.s.verts, sq.s.index, sq.support, sq.edge, sq.poly, &st);
    std::vector<std::array<std::uint32_t, 3>> tris;
    fan_triangulate(p, tris, &st);

    KRI_CHECK(st.inserted == 0 && st.degenerate_kept == 0);
    KRI_CHECK(tris.size() == 2);
    KRI_CHECK(tris[0][0] == sq.poly[0] && tris[0][1] == sq.poly[1] && tris[0][2] == sq.poly[2]);
    KRI_CHECK(tris[1][0] == sq.poly[0] && tris[1][1] == sq.poly[2] && tris[1][2] == sq.poly[3]);
}

/// 索引が `side` で張られていること。**平面3つ組で張ってはいけません**（別名問題）。
void test_index_is_geometric() {
    Square sq;
    for (int i = 1; i <= 3; ++i) sq.s.vertex(sq.support, sq.y0, sq.s.plane(Axis::X, i));
    for (int i = 1; i <= 3; ++i) sq.s.vertex(sq.support, sq.x4, sq.s.plane(Axis::Y, i));
    sq.s.build();

    // 支持平面 z=0 の上には**全頂点**が載る（この足場はすべて z=0 上に作っている）
    const auto* on_support = sq.s.index.find(sq.support);
    KRI_CHECK(on_support != nullptr);
    KRI_CHECK_MSG(on_support->size() == sq.s.verts.size(),
                  "支持平面の索引が全頂点を返していない（`side` で張れていない）");

    // 平面 y=0 の上には、その線上の点だけ: 角 2 つ + 追加 3 点 = 5
    const auto* on_y0 = sq.s.index.find(sq.y0);
    KRI_CHECK(on_y0 != nullptr);
    KRI_CHECK_MSG(on_y0->size() == 5,
                  "y=0 上の頂点数が想定と違う（実測 " + std::to_string(on_y0->size()) + "）");
    KRI_CHECK_MSG(on_y0->size() < sq.s.verts.size(), "絞り込みが効いていない");
}

}  // namespace

int main() {
    std::printf("\n  T 頂点の解決 — SPEC-phase2 §2.4.3 / §2.4.4 (3) の正の対照\n");
    test_inserts_in_order();
    test_rejects_outside_and_endpoints();
    test_other_edge_not_mixed();
    test_fan_avoids_degenerate();
    test_fan_fallback_keeps_degenerate();
    test_fan_unchanged_without_t();
    test_index_is_geometric();
    std::printf("\n");
    return kritest::finish("csg/tjunction");
}

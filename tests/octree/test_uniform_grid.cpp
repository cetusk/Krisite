// Krisite — 固定深度一様分割のテスト
//
// SPEC-phase1.md §2.2, §3.2, §4.2
//
// 中心的な検査は 2 つ。
//   1. 隣り合うセルが共有面について**同一の平面**を得ること（§5.1 の前提）
//   2. 閉領域割り当てが、面上で接するだけの三角形を**両方のセルに**渡すこと（§4.2）
#include <set>
#include <vector>

#include "krisite/octree/uniform_grid.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::geom;
using namespace krisite::octree;
using krisite::kCoordMin;
using krisite::mesh::TriMesh;
using kritest::cube;

namespace {

// ---- 境界座標 ----------------------------------------------------------------

void test_bounds() {
    for (unsigned d = 0; d <= 4; ++d) {
        const UniformGrid g(d);
        KRI_CHECK(g.per_axis() == (1u << d));
        KRI_CHECK(g.cell_count() == (std::size_t{1} << (3 * d)));
        // 端が座標範囲の両端に一致すること
        KRI_CHECK(g.bound(0) == kCoordMin);
        KRI_CHECK(g.bound(g.per_axis()) == -kCoordMin);
        // 一様な刻み
        for (std::uint32_t m = 0; m < g.per_axis(); ++m) {
            KRI_CHECK(g.bound(m + 1) - g.bound(m) == g.cell_size());
        }
        // セル 1 個ぶんの幅 x セル数 = 座標範囲の全幅 2^b
        KRI_CHECK(g.cell_size() * g.per_axis() == (std::int64_t{1} << krisite::kCoordBits));
        // 内部境界の枚数
        KRI_CHECK(g.interior_plane_count() == 3u * (g.per_axis() - 1));
    }
    // 深度 0 は分割なし
    KRI_CHECK(UniformGrid(0).cell_count() == 1);
    KRI_CHECK(UniformGrid(0).interior_plane_count() == 0);
}

/// 最大境界 +2^(b-1) は IPoint の範囲外だが、平面としては作れること（§3.2）。
void test_max_bound_outside_ipoint_range() {
    const UniformGrid g(2);
    const std::int64_t top = g.bound(g.per_axis());
    KRI_CHECK(top == -kCoordMin);
    KRI_CHECK(top > krisite::kCoordMax);  // IPoint には置けない
    const PlaneD pl = plane_axis_aligned(Axis::X, top);
    KRI_CHECK(!is_degenerate(pl));
    // 座標範囲内のすべての点はこの平面の負側
    kritest::Rng rng(7);
    for (int i = 0; i < 3000; ++i) {
        KRI_CHECK(side(pl, kritest::rand_extreme_point(rng)) == -1);
    }
}

// ---- 隣接セルが共有面の平面を共有すること（§5.1 の前提）----------------------

void test_shared_face_planes_are_identical() {
    for (unsigned d = 1; d <= 3; ++d) {
        const UniformGrid g(d);
        const std::uint32_t n = g.per_axis();
        for (std::uint32_t i = 0; i + 1 < n; ++i) {
            const CellIndex c0{i, 0, 0}, c1{i + 1, 0, 0};
            const auto p0 = g.cell_planes(c0);
            const auto p1 = g.cell_planes(c1);
            // c0 の +X 面（添字 1）と c1 の -X 面（添字 0）は同じ平面
            KRI_CHECK(plane_same(p0[1], p1[0]));
            KRI_CHECK(plane_cmp(p0[1], p1[0]) == 0);
            // 係数まで完全一致していること（同じ境界座標から作るため）
            KRI_CHECK(krisite::arith::cmp(p0[1].d, p1[0].d) == 0);
        }
        // Y / Z も同様
        for (std::uint32_t j = 0; j + 1 < n; ++j) {
            KRI_CHECK(plane_same(g.cell_planes({0, j, 0})[3], g.cell_planes({0, j + 1, 0})[2]));
        }
        for (std::uint32_t k = 0; k + 1 < n; ++k) {
            KRI_CHECK(plane_same(g.cell_planes({0, 0, k})[5], g.cell_planes({0, 0, k + 1})[4]));
        }
    }
}

/// 異なる境界の平面は別平面であること。
void test_distinct_bounds_are_distinct_planes() {
    const UniformGrid g(3);
    std::vector<PlaneD> ps;
    for (std::uint32_t m = 0; m <= g.per_axis(); ++m) {
        ps.push_back(plane_axis_aligned(Axis::X, g.bound(m)));
    }
    for (std::size_t a = 0; a < ps.size(); ++a) {
        for (std::size_t b2 = 0; b2 < ps.size(); ++b2) {
            KRI_CHECK(plane_same(ps[a], ps[b2]) == (a == b2));
        }
    }
}

// ---- 閉領域割り当て（§4.2）--------------------------------------------------

/// 三角形をどのセルに割り当てるか列挙する。
std::vector<std::size_t> assigned_cells(const UniformGrid& g, const Aabb& box, bool open) {
    std::vector<std::size_t> out;
    const std::uint32_t n = g.per_axis();
    for (std::uint32_t i = 0; i < n; ++i) {
        for (std::uint32_t j = 0; j < n; ++j) {
            for (std::uint32_t k = 0; k < n; ++k) {
                const CellIndex c{i, j, k};
                const bool hit = open ? assign_to_cell_open(box, g, c) : assign_to_cell(box, g, c);
                if (hit) out.push_back((static_cast<std::size_t>(i) * n + j) * n + k);
            }
        }
    }
    return out;
}

/// セル面上に**完全に乗る**三角形は、**半開区間なら上側のセルだけ**に入る（§3.4）。
///
/// **開領域とは別物です。** 開領域だとどのセルにも入らず消えます。片側に落ちるだけでは
/// ありません。境界 x = m の両隣のセルについて、下側は `lo >= chi`、上側は `hi <= clo` で
/// どちらも弾かれるためです。これが §10.5 の変異 2 が突く差で、
/// 予想より深刻な壊れ方をします（幾何が丸ごと欠落する）。
///
/// 半開は「乗る面をどちらか一方に決める」規則で、**幾何は失われません。**
void test_closed_vs_open_assignment() {
    const UniformGrid g(1);  // 2x2x2、境界は x=0（中央）
    const std::int64_t mid = g.bound(1);
    KRI_CHECK(mid == 0);  // b=21 なら -2^20 + 2^20 = 0

    // x = 0 の平面上にぴったり乗る三角形
    const IPoint a{static_cast<std::int32_t>(mid), 10, 10};
    const IPoint b{static_cast<std::int32_t>(mid), 200, 10};
    const IPoint c{static_cast<std::int32_t>(mid), 10, 200};
    const Aabb box = triangle_aabb(a, b, c);

    const auto closed = assigned_cells(g, box, false);
    const auto open = assigned_cells(g, box, true);

    // **半開区間 [lo, hi)**（SPEC-phase2 §3.4）: 境界平面に**完全に乗る**面は
    // **上側のセルだけ**に入る。閉領域だと両側に入り、適応分割で入れ子の重複になる
    std::set<std::uint32_t> ci_half;
    for (std::size_t id : closed) ci_half.insert(static_cast<std::uint32_t>(id / 4));
    KRI_CHECK_MSG(ci_half.size() == 1, "半開なら x 方向は上側のセルだけのはず（§3.4）");
    KRI_CHECK_MSG(*ci_half.begin() == 1, "上側（i=1）に入るはず。lo が境界に一致する側");
    // **開領域とは別物です。** 開領域だとどのセルにも入らず、幾何が丸ごと消えます
    KRI_CHECK_MSG(open.empty(),
                  "開領域だとセル境界上の三角形がどのセルにも入らず消える（変異 2 の狙い）");

    // 境界から 1 格子ずらせば、開領域でもちょうど 1 セルに入る（対照）
    const IPoint a2{static_cast<std::int32_t>(mid) + 1, 10, 10};
    const IPoint b2{static_cast<std::int32_t>(mid) + 1, 200, 10};
    const IPoint c2{static_cast<std::int32_t>(mid) + 1, 10, 200};
    const auto open2 = assigned_cells(g, triangle_aabb(a2, b2, c2), true);
    KRI_CHECK(open2.size() == 1);
}

/// 割り当ては取りこぼさない: 三角形の頂点を含むセルは必ず割り当て対象。
void test_assignment_never_misses() {
    kritest::Rng rng(11);
    for (unsigned d = 0; d <= 2; ++d) {
        const UniformGrid g(d);
        const std::uint32_t n = g.per_axis();
        for (int iter = 0; iter < 2000; ++iter) {
            const IPoint p0 = kritest::rand_point(rng);
            const IPoint p1 = kritest::rand_point(rng);
            const IPoint p2 = kritest::rand_point(rng);
            const Aabb box = triangle_aabb(p0, p1, p2);
            for (const IPoint& p : {p0, p1, p2}) {
                // p を含むセルの添字
                auto idx = [&](std::int64_t v) {
                    std::int64_t t = (v - kCoordMin) / g.cell_size();
                    if (t < 0) t = 0;
                    if (t >= static_cast<std::int64_t>(n)) t = n - 1;
                    return static_cast<std::uint32_t>(t);
                };
                const CellIndex c{idx(p.x), idx(p.y), idx(p.z)};
                KRI_CHECK(assign_to_cell(box, g, c));
            }
        }
    }
}

/// 立方体メッシュ全体を割り当てたとき、深度が上がるとセルあたりの三角形は減る。
void test_cube_assignment_counts() {
    const TriMesh m = cube(-200, -200, -200, 400);
    for (unsigned d = 0; d <= 3; ++d) {
        const UniformGrid g(d);
        std::size_t total = 0, nonempty = 0;
        const std::uint32_t n = g.per_axis();
        for (std::uint32_t i = 0; i < n; ++i) {
            for (std::uint32_t j = 0; j < n; ++j) {
                for (std::uint32_t k = 0; k < n; ++k) {
                    std::size_t cnt = 0;
                    for (const auto& t : m.triangles) {
                        const Aabb box =
                            triangle_aabb(m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]]);
                        if (assign_to_cell(box, g, {i, j, k})) ++cnt;
                    }
                    total += cnt;
                    if (cnt > 0) ++nonempty;
                }
            }
        }
        // 深度 0 は全三角形が唯一のセルに入る
        if (d == 0) {
            KRI_CHECK(total == m.triangles.size());
            KRI_CHECK(nonempty == 1);
        }
        // どの深度でも、すべての三角形は少なくとも 1 セルに入る
        KRI_CHECK(total >= m.triangles.size());
        KRI_CHECK(nonempty >= 1);
    }
}

}  // namespace

int main() {
    test_bounds();
    test_max_bound_outside_ipoint_range();
    test_shared_face_planes_are_identical();
    test_distinct_bounds_are_distinct_planes();
    test_closed_vs_open_assignment();
    test_assignment_never_misses();
    test_cube_assignment_counts();
    return kritest::finish("octree/uniform_grid");
}

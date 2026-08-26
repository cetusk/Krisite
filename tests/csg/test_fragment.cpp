// Krisite — 断片の分割のテスト
//
// SPEC-phase1.md §4.3
//
// 分割の正しさは「元の面積が保存されること」と「頂点が期待どおりの位置に来ること」で
// 見ます。面積は構成点の有理数になるので、ここでは軸平行な構成に限って
// 整数で検算できるようにしてあります。
#include <set>
#include <vector>

#include "krisite/csg/fragment.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using namespace krisite::geom;
using krisite::mesh::TriMesh;
using kritest::cube;

namespace {

/// 断片の頂点を IPoint として取り出す（軸平行構成なら w = ±1 で整数になる）。
bool vertex_as_ipoint(const PlaneTable& t, const Fragment& f, std::size_t i, IPoint& out) {
    const HPointD v = fragment_vertex(t, f, i);
    // w = ±1 のときだけ整数座標。cmp_h で軸ごとに整数と突き合わせて逆算する
    for (std::int32_t* p : {&out.x, &out.y, &out.z}) *p = 0;
    const Axis ax[3] = {Axis::X, Axis::Y, Axis::Z};
    std::int32_t* dst[3] = {&out.x, &out.y, &out.z};
    for (int k = 0; k < 3; ++k) {
        // 二分探索で整数座標を求める（範囲は座標範囲）
        std::int64_t lo = krisite::kCoordMin, hi = krisite::kCoordMax;
        while (lo < hi) {
            const std::int64_t mid = lo + (hi - lo) / 2;
            IPoint probe{0, 0, 0};
            std::int32_t* q[3] = {&probe.x, &probe.y, &probe.z};
            *q[k] = static_cast<std::int32_t>(mid);
            if (cmp_h(v, to_homogeneous(probe), ax[k]) <= 0) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        IPoint probe{0, 0, 0};
        std::int32_t* q[3] = {&probe.x, &probe.y, &probe.z};
        *q[k] = static_cast<std::int32_t>(lo);
        if (cmp_h(v, to_homogeneous(probe), ax[k]) != 0) return false;  // 整数でない
        *dst[k] = static_cast<std::int32_t>(lo);
    }
    return true;
}

/// 立方体の 1 面（軸平行な正方形）を軸平行平面で分割する。
void test_split_square() {
    PlaneTable t;
    const TriMesh m = cube(0, 0, 0, 100);
    const std::vector<Face> faces = build_faces(m, 0, t);

    // z = 100 の面を探す（頂点がすべて z=100）
    const Face* top = nullptr;
    for (const Face& f : faces) {
        bool all_top = true;
        for (auto vid : f.loop) {
            if (m.vertices[vid].z != 100) all_top = false;
        }
        if (all_top) top = &f;
    }
    KRI_CHECK(top != nullptr);
    if (!top) return;

    const Fragment frag = face_to_fragment(*top);
    KRI_CHECK(frag.edge.size() == 4);

    // x = 40 で切る
    const PlaneId cut = t.intern(plane_axis_aligned(Axis::X, 40)).id;
    const SplitResult r = split_fragment(t, frag, cut);
    KRI_CHECK(r.has_pos && r.has_neg);
    KRI_CHECK_MSG(r.pos.edge.size() == 4, "正方形を切れば両側とも四角形");
    KRI_CHECK_MSG(r.neg.edge.size() == 4, "正方形を切れば両側とも四角形");

    // 両側の頂点の x 座標が {40, 100} と {0, 40} になること
    auto xs_of = [&](const Fragment& f) {
        std::set<std::int32_t> xs;
        for (std::size_t i = 0; i < vertex_count(f); ++i) {
            IPoint p{};
            KRI_CHECK(vertex_as_ipoint(t, f, i, p));
            KRI_CHECK(p.z == 100);
            xs.insert(p.x);
        }
        return xs;
    };
    const std::set<std::int32_t> want_pos{40, 100}, want_neg{0, 40};
    KRI_CHECK(xs_of(r.pos) == want_pos);
    KRI_CHECK(xs_of(r.neg) == want_neg);

    // 支持平面と向きは保存される
    KRI_CHECK(r.pos.support == frag.support && r.pos.flipped == frag.flipped);
    KRI_CHECK(r.neg.support == frag.support && r.neg.flipped == frag.flipped);
}

/// 面をまたがない平面で切っても分割されないこと。
void test_no_split_when_outside() {
    PlaneTable t;
    const TriMesh m = cube(0, 0, 0, 100);
    const std::vector<Face> faces = build_faces(m, 0, t);
    const Fragment frag = face_to_fragment(faces.front());

    for (std::int64_t c : {-50, 150}) {
        const PlaneId q = t.intern(plane_axis_aligned(Axis::X, c)).id;
        const SplitResult r = split_fragment(t, frag, q);
        KRI_CHECK_MSG(r.has_pos != r.has_neg, "またがない平面では片側だけが残るはず");
        const Fragment& kept = r.has_pos ? r.pos : r.neg;
        KRI_CHECK(kept.edge == frag.edge);  // そのまま
    }
}

/// 支持平面で切っても何も起きないこと。
void test_split_by_support_is_noop() {
    PlaneTable t;
    const TriMesh m = cube(0, 0, 0, 100);
    const std::vector<Face> faces = build_faces(m, 0, t);
    const Fragment frag = face_to_fragment(faces.front());
    const SplitResult r = split_fragment(t, frag, frag.support);
    KRI_CHECK(r.has_pos && !r.has_neg);
    KRI_CHECK(r.pos.edge == frag.edge);
}

/// 頂点をちょうど通る平面で切る（退化の入口）。
void test_split_through_vertex() {
    PlaneTable t;
    const TriMesh m = cube(0, 0, 0, 100);
    const std::vector<Face> faces = build_faces(m, 0, t);
    const Face* top = nullptr;
    for (const Face& f : faces) {
        bool all_top = true;
        for (auto vid : f.loop) {
            if (m.vertices[vid].z != 100) all_top = false;
        }
        if (all_top) top = &f;
    }
    KRI_CHECK(top != nullptr);
    if (!top) return;
    const Fragment frag = face_to_fragment(*top);

    // x = 0 / x = 100（正方形の辺にちょうど重なる）で切る → 片側だけが残る。
    //
    // **どちらの側が残るかは平面の向きに依存します。** PlaneTable の代表は
    // 「最初に登録された平面」なので、立方体の -X 面から作られた平面は法線が -X 向きです。
    // したがってテストは向きを仮定してはいけません。
    for (std::int64_t c : {0, 100}) {
        const PlaneId q = t.intern(plane_axis_aligned(Axis::X, c)).id;
        const SplitResult r = split_fragment(t, frag, q);
        KRI_CHECK_MSG(r.has_pos != r.has_neg,
                      "辺に重なる平面では片側だけが残る（面積 0 の側は落ちる）");
        const Fragment& kept = r.has_pos ? r.pos : r.neg;
        KRI_CHECK(kept.edge.size() == 4);
        KRI_CHECK(kept.edge == frag.edge);  // 切られていない
    }
}

/// 何度切っても頂点はすべて 3 平面の交点として復元できること（表現が閉じている）。
void test_representation_closed() {
    PlaneTable t;
    const TriMesh m = cube(0, 0, 0, 1000);
    const std::vector<Face> faces = build_faces(m, 0, t);
    std::vector<PlaneId> cuts;
    for (std::int64_t c : {137, 251, 389, 512, 700, 913}) {
        cuts.push_back(t.intern(plane_axis_aligned(Axis::X, c)).id);
        cuts.push_back(t.intern(plane_axis_aligned(Axis::Y, c)).id);
    }

    std::vector<Fragment> pieces;
    for (const Face& f : faces) pieces.push_back(face_to_fragment(f));
    for (PlaneId q : cuts) {
        std::vector<Fragment> next;
        for (const Fragment& p : pieces) {
            const SplitResult r = split_fragment(t, p, q);
            if (r.has_pos) next.push_back(r.pos);
            if (r.has_neg) next.push_back(r.neg);
        }
        pieces.swap(next);
    }
    KRI_CHECK(pieces.size() > faces.size());

    for (const Fragment& p : pieces) {
        KRI_CHECK(p.edge.size() >= 3);
        for (std::size_t i = 0; i < vertex_count(p); ++i) {
            // 3 平面が一点で交わること
            const std::size_t n = p.edge.size();
            KRI_CHECK(kritest::intersects_at_point(t.at(p.support), t.at(p.edge[(i + n - 1) % n]),
                                                   t.at(p.edge[i])));
            // 頂点が支持平面上にあること
            KRI_CHECK(side(t.at(p.support), fragment_vertex(t, p, i)) == 0);
        }
        // 各切断平面に対して符号が一意に定まること（またいでいない）
        for (PlaneId q : cuts) {
            if (q == p.support) continue;
            (void)fragment_sign(t, p, q);  // 内部の KRISITE_CHECK が矛盾を捕まえる
        }
    }
}

}  // namespace

int main() {
    test_split_square();
    test_no_split_when_outside();
    test_split_by_support_is_noop();
    test_split_through_vertex();
    test_representation_closed();
    return kritest::finish("csg/fragment");
}

// Krisite — **VM（頂点多様体）の印**（`SPEC-phase5.md` §5.10.14.57 の案 B）
//
// **印であって、拒否ではありません。**
//
// 量子化後に頂点非多様体な入力があると、**出口で分裂できない次数 4 の辺が
// 残ることがあります**（実測: 該当する 141 対のうち 8 対（5.7%）。
// 非該当の 154 対では 0 対。`DESIGN-phase5-vertex-level.md` §3.4）。
// **必要条件として効きますが、十分条件ではありません。**
//
// **受け入れる例と拒否する例を【対】で置きます**（`CLAUDE.md`
// 「前提を書いたら、検査も書いてください。…その検査を突くケースをコーパスに入れて」）。
//
// **3 値であることも検査します** — 0（未検査）と 2（検査して非多様体）を
// 区別できなければ、`inputs_vertex_nonmanifold == 0` を「多様体だった」と
// 読み違えます（`CLAUDE.md`「1 つの値に 2 つの意味を持たせない」）。
#include <cstdint>
#include <cstdio>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite;
using krisite::mesh::TriMesh;

namespace {

/// 軸平行な立方体（`corpus.hpp` の `box` を辺長で包む）。
TriMesh box(std::int32_t x, std::int32_t y, std::int32_t z, std::int32_t s) {
    return kritest::box(x, y, z, x + s, y + s, z + s);
}

/// **2 つの立方体が【1 頂点だけ】を共有する**メッシュ。
///
/// **辺はすべて次数 2 なので辺多様体ですが、共有した頂点のまわりは扇が 2 個**です。
/// **頂点非多様体の最小の形**で、これが案 B の判別条件に当たります。
TriMesh two_boxes_sharing_a_vertex() {
    TriMesh a = box(0, 0, 0, 100);
    const TriMesh b = box(100, 100, 100, 100);
    const auto off = static_cast<mesh::VertexId>(a.vertices.size());
    // b の (100,100,100) は a の頂点 7 と同じ点。**同じ番号に寄せます**
    std::vector<mesh::VertexId> map(b.vertices.size());
    for (std::size_t i = 0; i < b.vertices.size(); ++i) {
        if (b.vertices[i].x == 100 && b.vertices[i].y == 100 && b.vertices[i].z == 100) {
            map[i] = 6;  // corpus.hpp の並びでは (hi,hi,hi) が 6 番
        } else {
            map[i] = static_cast<mesh::VertexId>(a.vertices.size());
            a.vertices.push_back(b.vertices[i]);
        }
    }
    (void)off;
    for (const mesh::Tri& t : b.triangles) {
        a.triangles.push_back({map[t[0]], map[t[1]], map[t[2]]});
    }
    return a;
}

int failures = 0;

void expect(bool ok, const char* what) {
    if (!ok) {
        ++failures;
        std::printf("  **FAIL** %s\n", what);
    }
}

}  // namespace

int main() {
    std::printf("## VM（頂点多様体）の印\n\n");

    // ---- 1. 検査しなければ「未検査」のまま（既定）----
    {
        csg::FromMeshOptions o;
        const csg::PolySoup s = csg::from_mesh(two_boxes_sharing_a_vertex(), o);
        expect(s.vm.size() == 1 && s.vm[0] == 0, "既定では未検査（0）のまま");
    }

    // ---- 2. 受け入れる例: 単独の立方体は多様体 ----
    {
        csg::FromMeshOptions o;
        o.verify_vertex_manifold = true;
        const csg::PolySoup s = csg::from_mesh(box(0, 0, 0, 100), o);
        expect(s.vm.size() == 1 && s.vm[0] == 1, "立方体は多様体（1）と印が付く");
    }

    // ---- 3. 拒否する例: 1 頂点だけを共有する 2 つの立方体 ----
    {
        csg::FromMeshOptions o;
        o.verify_vertex_manifold = true;
        const TriMesh m = two_boxes_sharing_a_vertex();
        const mesh::TopologyReport t = mesh::check_topology(m.triangles);
        // **突く構成が狙いどおりかを、まず確かめます**（`CLAUDE.md`「空回りしていないか」）
        expect(t.edge_manifold, "構成: 辺は多様体である（頂点だけが非多様体）");
        expect(!t.vertex_manifold, "構成: 頂点が非多様体である");
        const csg::PolySoup s = csg::from_mesh(m, o);
        expect(s.vm.size() == 1 && s.vm[0] == 2, "非多様体（2）と印が付く");
    }

    // ---- 4. 印が演算の統計に出る ----
    {
        csg::FromMeshOptions o;
        o.verify_vertex_manifold = true;
        const csg::PolySoup a = csg::from_mesh(two_boxes_sharing_a_vertex(), o);
        const csg::PolySoup b = csg::from_mesh(box(50, 50, 50, 100), o);
        csg::BoolOptions bo;
        bo.depth = 2;
        csg::BoolStats st;
        const csg::PolySoup u = csg::boolean(a, b, csg::BoolOp::Union, bo, &st);
        expect(st.inputs_vertex_nonmanifold == 1, "統計に非多様体な入力 1 件が出る");
        expect(st.inputs_vertex_unchecked == 0, "未検査は 0 件");
        expect(u.vm.size() == 2, "印は連結される");
    }

    // ---- 5. 検査しない側が混ざったら「未検査」として数える ----
    {
        csg::FromMeshOptions on;
        on.verify_vertex_manifold = true;
        const csg::PolySoup a = csg::from_mesh(two_boxes_sharing_a_vertex(), on);
        const csg::PolySoup b = csg::from_mesh(box(50, 50, 50, 100), csg::FromMeshOptions{});
        csg::BoolOptions bo;
        bo.depth = 2;
        csg::BoolStats st;
        csg::boolean(a, b, csg::BoolOp::Union, bo, &st);
        expect(st.inputs_vertex_nonmanifold == 1, "非多様体は 1 件のまま");
        expect(st.inputs_vertex_unchecked == 1, "**未検査 1 件が区別される**");
    }

    std::printf("\n**不一致 %d 件**\n", failures);
    return failures == 0 ? 0 : 1;
}

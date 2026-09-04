// Krisite — 「内側」の定義（SPEC-phase3 §5.2.0 / §9 のケース 27）
//
// **正しい定義は $\text{内側}_i \iff w_i > 0$ です。**
// 初版の $w \ne 0$ は誤りで、**自己交差した入力では巻き数が負になり、
// 裏返った面に囲まれた領域も「内側」と判定されていました。**
//
// **このケースは既存の検査をすべて素通りします**（§9.1）。
//
//   位相                多様体な出力が返る（立方体が 1 個出るだけ）
//   体積の恒等式        向きの反転は符号を変えるだけで、恒等式は成立する
//   分割戦略不変性      分割の仕方に依らず、同じ誤った出力が出る
//   決定性              誤った出力が決定的に出る
//   二項正解器との一致  **正解器も同じ誤りを持っている**
//
// **だから絶対値で検査します。** 比較する相手がありません。
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite;
using krisite::mesh::TriMesh;

namespace {

csg::BoolOptions opts(unsigned depth, bool adaptive) {
    csg::BoolOptions o;
    o.depth = depth;
    o.adaptive = adaptive;
    o.cull_planes = true;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    return o;
}

/// self-union（$n = 1$）= 同じスープを 2 度使う和。**ケース 19 と同じ形**です。
csg::SoupMesh self_union(const TriMesh& m, const csg::BoolOptions& o, bool legacy) {
    csg::PolySoup a = csg::from_mesh(m), b = csg::from_mesh(m);
    a.indicator.legacy_nonzero_inside = legacy;
    b.indicator.legacy_nonzero_inside = legacy;
    csg::ToMeshOptions tm;
    tm.split_contacts = true;
    return csg::to_mesh(csg::boolean(a, b, csg::BoolOp::Union, o), tm);
}

std::string bytes(const csg::SoupMesh& m) {
    std::string s;
    const auto put = [&s](const void* p, std::size_t n) {
        s.append(static_cast<const char*>(p), n);
    };
    const std::size_t nv = m.vertices.size(), nt = m.triangles.size();
    put(&nv, sizeof nv);
    put(&nt, sizeof nt);
    for (const auto& v : m.vertices) put(&v, sizeof v);
    for (const auto& t : m.triangles) put(&t, sizeof t);
    return s;
}

}  // namespace

int main() {
    std::printf("=== 「内側」の定義（SPEC-phase3 §5.2.0）===\n");

    // **入口の契約**: 向きを反転しても $\partial S = 0$ は保たれます（§4.0）
    for (const auto& c : {std::pair<const char*, TriMesh>{"27a", kritest::cases::case27a()},
                          {"27b", kritest::cases::case27b()},
                          {"27c", kritest::cases::case27c()}}) {
        KRI_CHECK_MSG(mesh::boundary_is_zero(c.second),
                      std::string("ケース ") + c.first + ": 入力が PWN でない（∂S ≠ 0）");
    }
    std::printf("  入口の契約: 27a / 27b / 27c とも ∂S = 0（**向きの反転は壊しません**）\n");

    for (unsigned pass = 0; pass < 4; ++pass) {
        const bool adaptive = (pass == 3);
        const unsigned depth = adaptive ? 5 : pass + 2;
        const csg::BoolOptions o = opts(depth, adaptive);
        const std::string tag =
            std::string("（") + (adaptive ? "適応" : "深度 " + std::to_string(depth)) + "）";

        // ---- 27a: すべての領域が w = -1。**正しくは立体が無い** --------------
        {
            const csg::SoupMesh now = self_union(kritest::cases::case27a(), o, false);
            const csg::SoupMesh old = self_union(kritest::cases::case27a(), o, true);
            KRI_CHECK_MSG(now.triangles.empty(),
                          "ケース 27a" + tag +
                              ": **裏返った立方体が「立体」として出力された**"
                              "（すべての領域が w = -1 なので空になるはず）" +
                              kritest::pair_msg(std::size_t{0}, now.triangles.size()));
            // ★ **旧定義では誤った出力が出ることを、番人として固定します。**
            // これが空になったら、**このケースが何も検査していない**ということです
            KRI_CHECK_MSG(
                !old.triangles.empty(),
                "ケース 27a" + tag + ": **旧定義でも空になった。このケースが空回りしています**");
        }

        // ---- 27b: 一部だけ裏返る。**かたまりの数で違いが出る** ---------------
        {
            const csg::SoupMesh now = self_union(kritest::cases::case27b(), o, false);
            const csg::SoupMesh old = self_union(kritest::cases::case27b(), o, true);
            const mesh::TopologyReport rn = mesh::check_topology(now.triangles);
            const mesh::TopologyReport ro = mesh::check_topology(old.triangles);
            KRI_CHECK_MSG(rn.components == 1, "ケース 27b" + tag + ": かたまりが 1 個でない" +
                                                  kritest::pair_msg(std::size_t{1}, rn.components));
            KRI_CHECK_MSG(ro.components == 2,
                          "ケース 27b" + tag +
                              ": **旧定義でかたまりが 2 個にならない。空回りしています**" +
                              kritest::pair_msg(std::size_t{2}, ro.components));
        }

        // ---- 27c: 対照。**両定義で出力が完全に一致するはず** -----------------
        {
            const csg::SoupMesh now = self_union(kritest::cases::case27c(), o, false);
            const csg::SoupMesh old = self_union(kritest::cases::case27c(), o, true);
            KRI_CHECK_MSG(
                bytes(now) == bytes(old),
                "ケース 27c" + tag + ": **巻き数が負にならない入力なのに、定義で出力が変わった**");
            KRI_CHECK_MSG(!now.triangles.empty(), "ケース 27c" + tag + ": 出力が空（空回り）");
        }
    }
    std::printf("  27a: 正しい定義で空 / 旧定義で非空（**旧定義の誤りを固定**）\n");
    std::printf("  27b: 正しい定義でかたまり 1 個 / 旧定義で 2 個\n");
    std::printf("  27c: **両定義で出力がバイト単位で一致**（対照）\n");
    return kritest::finish("csg/winding");
}

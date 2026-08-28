// Krisite — 入口・中核・出口の分離（SPEC-phase3.md §14 の CP2）
//
//     from_mesh : TriMesh → PolySoup
//     boolean   : PolySoup × PolySoup → PolySoup   ★ CSG について閉じる
//     to_mesh   : PolySoup → TriMesh
//
// **見るのは 4 つです。**
//
//   1. 往復（`to_mesh(from_mesh(m))`）が元の位相を再現すること
//   2. §10.1 二項正解器（`boolean_op`）と $(C, \chi)$ が一致すること
//   3. §10.3 連鎖が丸めを経由しないこと（**中間に TriMesh を作らない**）
//   4. §10.3.1 連鎖でビット幅が伸びないこと ★ 契約の成立条件
//
// **三角形の集合と断片数は一致しません**（§10.1）。面併合をやめたので分割が違います。
#include <cstdio>
#include <string>
#include <type_traits>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kMaxDepth = 3;

/// スープの構成点の最大ビット幅（§10.3.1）。
std::size_t soup_bits(const csg::PolySoup& s) {
    std::size_t mx = 0;
    for (const csg::Poly& q : s.polys) {
        for (std::size_t i = 0; i < csg::vertex_count(q.frag); ++i) {
            const geom::HPointD v = csg::fragment_vertex(s.table, q.frag, i);
            for (std::size_t w : {arith::min_bits(v.x), arith::min_bits(v.y), arith::min_bits(v.z),
                                  arith::min_bits(v.w)}) {
                mx = (w > mx) ? w : mx;
            }
        }
    }
    return mx;
}

const char* op_name(csg::BoolOp op) {
    switch (op) {
        case csg::BoolOp::Union:
            return "∪";
        case csg::BoolOp::Intersection:
            return "∩";
        default:
            return "\\";
    }
}

/// 1. 往復。**入口と出口だけで閉じる検査**なので、中核の誤りが混ざりません。
void test_round_trip() {
    std::size_t n = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        for (int which = 0; which < 2; ++which) {
            const TriMesh m = (which == 0) ? c.make_a() : c.make_b();
            const csg::SoupMesh r = csg::to_mesh(csg::from_mesh(m));
            const TopologyReport t0 = check_topology(m.triangles);
            const TopologyReport t1 = check_topology(r.triangles);
            const std::string tag = std::string("ケース ") + c.id + (which == 0 ? " A" : " B");
            KRI_CHECK_MSG(t1.ok(), tag + ": 往復で多様体でなくなった");
            KRI_CHECK_MSG(
                t0.components == t1.components,
                tag + ": 往復で C が変わった" + kritest::pair_msg(t0.components, t1.components));
            KRI_CHECK_MSG(t0.chi == t1.chi,
                          tag + ": 往復で χ が変わった" + kritest::pair_msg(t0.chi, t1.chi));
            ++n;
        }
    }
    std::printf("    往復 %zu 件\n", n);
}

/// 2. §10.1 二項正解器との一致。**分裂は両方 OFF で比べます**（意味論を揃える）。
void test_matches_binary_oracle() {
    std::size_t n = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (csg::BoolOp op :
             {csg::BoolOp::Union, csg::BoolOp::Intersection, csg::BoolOp::Difference}) {
            for (unsigned d = 0; d <= kMaxDepth; ++d) {
                const csg::BoolOptions o = kritest::phase1_options(d);
                const csg::BoolMesh ref = csg::boolean_op(a, b, op, o);
                csg::ToMeshOptions tm;
                tm.split_contacts = false;
                const csg::SoupMesh got =
                    csg::to_mesh(csg::boolean(csg::from_mesh(a), csg::from_mesh(b), op, o), tm);
                const TopologyReport t0 = check_topology(ref.triangles);
                const TopologyReport t1 = check_topology(got.triangles);
                const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) +
                                        "（深度 " + std::to_string(d) + "）";
                KRI_CHECK_MSG(
                    t0.components == t1.components,
                    tag + ": C が正解器と違う" + kritest::pair_msg(t0.components, t1.components));
                KRI_CHECK_MSG(t0.chi == t1.chi,
                              tag + ": χ が正解器と違う" + kritest::pair_msg(t0.chi, t1.chi));
                ++n;
            }
        }
    }
    std::printf("    正解器との比較 %zu 件\n", n);
}

/// 3. §10.3 連鎖の厳密性 + 4. §10.3.1 ビット幅が伸びないこと。
void test_chain_is_exact() {
    // **型が閉じていること**をコンパイル時に固定する
    static_assert(std::is_same_v<decltype(csg::boolean(std::declval<const csg::PolySoup&>(),
                                                       std::declval<const csg::PolySoup&>(),
                                                       csg::BoolOp::Union, csg::BoolOptions{})),
                                 csg::PolySoup>,
                  "boolean の戻り値が PolySoup でない（型が閉じていない）");

    std::size_t b1 = 0, b2 = 0, b3 = 0, n = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        // 3 つ目は**別のメッシュ**。同じメッシュを 2 度使うと曲面が二重になり、
        // 内外の 1 ビットでは「2 度跨いだ」を表せません（CP3 の WNV が要る）。
        const TriMesh d = kritest::corpus()[1].make_b();
        const csg::BoolOptions o = kritest::phase1_options(1);
        const csg::PolySoup SA = csg::from_mesh(a), SB = csg::from_mesh(b), SD = csg::from_mesh(d);

        const csg::PolySoup s1 = csg::boolean(SA, SB, csg::BoolOp::Union, o);
        const csg::PolySoup s2 = csg::boolean(s1, SD, csg::BoolOp::Difference, o);
        const csg::PolySoup s3 = csg::boolean(s2, SA, csg::BoolOp::Union, o);

        const std::string tag = std::string("ケース ") + c.id;
        // **中間に TriMesh を作っていないこと。** sources は入力そのもので、
        // 段が増えても【増えるだけ】で、丸められた頂点が現れません。
        KRI_CHECK_MSG(s1.source_count() == 2 && s2.source_count() == 3 && s3.source_count() == 4,
                      tag + ": sources の数が段数と合わない");
        for (std::size_t i = 0; i < a.vertices.size(); ++i) {
            KRI_CHECK_MSG(s3.sources[0].vertices[i].x == a.vertices[i].x &&
                              s3.sources[0].vertices[i].y == a.vertices[i].y &&
                              s3.sources[0].vertices[i].z == a.vertices[i].z,
                          tag + ": 連鎖で source の頂点が変わった（丸めが入っている）");
        }
        b1 = std::max(b1, soup_bits(s1));
        b2 = std::max(b2, soup_bits(s2));
        b3 = std::max(b3, soup_bits(s3));

        // **同じメッシュを 2 度使う連鎖**（CP3 の段 1 で開きました）。
        //
        // 内外の 1 ビットでは「同じ曲面を 2 度跨いだ」を表せず、CP2 では 17 / 44 が
        // 食い違っていました。**巻き数なら w = 2 になるだけ**です。
        {
            const csg::PolySoup ub = csg::boolean(csg::boolean(SA, SB, csg::BoolOp::Union, o), SB,
                                                  csg::BoolOp::Difference, o);
            const csg::PolySoup ab = csg::boolean(SA, SB, csg::BoolOp::Difference, o);
            csg::ToMeshOptions t2;
            t2.split_contacts = false;
            const TopologyReport t_ub = check_topology(csg::to_mesh(ub, t2).triangles);
            const TopologyReport t_ab = check_topology(csg::to_mesh(ab, t2).triangles);
            KRI_CHECK_MSG(t_ub.components == t_ab.components && t_ub.chi == t_ab.chi,
                          tag + ": (A∪B)\\B と A\\B が食い違う（同じ曲面を 2 度跨ぐ配置）" +
                              kritest::pair_msg(t_ab.chi, t_ub.chi));
        }

        // 自己整合: (A ∪ B) \ D  ≡  (A \ D) ∪ (B \ D)
        csg::ToMeshOptions tm;
        tm.split_contacts = false;
        const csg::PolySoup rhs =
            csg::boolean(csg::boolean(SA, SD, csg::BoolOp::Difference, o),
                         csg::boolean(SB, SD, csg::BoolOp::Difference, o), csg::BoolOp::Union, o);
        const TopologyReport tl = check_topology(csg::to_mesh(s2, tm).triangles);
        const TopologyReport tr = check_topology(csg::to_mesh(rhs, tm).triangles);
        KRI_CHECK_MSG(
            tl.components == tr.components && tl.chi == tr.chi,
            tag + ": (A∪B)\\D と (A\\D)∪(B\\D) が食い違う" + kritest::pair_msg(tl.chi, tr.chi));
        ++n;
    }

    std::printf("    連鎖 %zu 件。最大ビット幅: 1 段 %zu / 2 段 %zu / 3 段 %zu（上界 %zu）\n", n,
                b1, b2, b3, geom::bits::kHomoXyz);
    // §10.3.1: **連鎖でビット幅は伸びません。** CSG は新しい平面を作らないので、
    // 構成点は入力平面の集合から 3 枚を選んだ交点のままです。
    KRI_CHECK_MSG(b2 == b1 && b3 == b1,
                  "**連鎖でビット幅が伸びた。** どこかで平面を構成しています（§10.3.1）" +
                      kritest::pair_msg(b1, b3));
    KRI_CHECK_MSG(b1 <= geom::bits::kHomoXyz, "構成点が理論上界を超えた");
    KRI_CHECK_MSG(b1 > 0, "構成点が 1 つも無い。**空回りです**");
}

}  // namespace

int main() {
    std::printf("\n  入口・中核・出口の分離 — SPEC-phase3 §14 の CP2\n");
    test_round_trip();
    test_matches_binary_oracle();
    test_chain_is_exact();
    std::printf("\n");
    return kritest::finish("csg/soup");
}

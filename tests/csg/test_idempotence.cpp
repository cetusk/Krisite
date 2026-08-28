// Krisite — §10.2.2 冪等性 と §10.2.3 平面 ID 帰属
//
// SPEC-phase1.md §10.2
//
// **GMP も外部正解器も要りません。** 既定のテスト実行で常に走ります。
//
// §10.2.2 は深度によって検査を変えます。面併合（§3.5）が入力の三角形分割を保存
// しないので、**生の A と比較すると再分割の差で落ちます。**
//
//   深度 0   … A∩A、A∪A が `retriangulate(A)` と三角形単位で一致
//   深度 ≥ 1 … A∩A と A∪A が互いに三角形単位で一致
//
// 深度 ≥ 1 の比較が成立する根拠は、どちらも同じセル分割を通り、どちらも
// 「A の全断片を向きそのままで 1 枚ずつ残す」はずだからです。両者は共平面重複の
// 処理で**逆向きの選択規則**を通るので、規則の食い違いを検出できます。
#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::geom::HPointD;
using krisite::mesh::TriMesh;

namespace {

/// §3.5 の面併合と扇分割だけを通した A（ブール演算を経ていない）。
///
/// **生の入力 A ではありません。** 面併合は入力の三角形分割を保存しないので、
/// 冪等性の比較対象はこちらです（§10.2.2 の但し書き）。
BoolMesh retriangulate(const TriMesh& m) {
    PlaneTable table;
    const std::vector<Face> faces = build_faces(m, 0, table);
    BoolMesh r;
    r.vertices.reserve(m.vertices.size());
    for (const krisite::geom::IPoint& p : m.vertices) {
        r.vertices.push_back(krisite::geom::to_homogeneous(p));
    }
    for (const Face& f : faces) {
        for (std::size_t i = 1; i + 1 < f.loop.size(); ++i) {
            r.triangles.push_back({f.loop[0], f.loop[i], f.loop[i + 1]});
        }
    }
    return r;
}

/// 三角形を「頂点の値」で正準化したもの。順序差と頂点番号の違いを吸収します。
///
/// 巡回順は保ったまま最小の ID を先頭に回すので、**向きは保たれます。**
/// 向きの違いは検出されます。
using CanonTri = std::array<std::uint32_t, 3>;

std::vector<CanonTri> canonical(const BoolMesh& m) {
    // 値で頂点を同一視して番号を振り直す
    std::vector<std::uint32_t> ord(m.vertices.size());
    for (std::uint32_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](std::uint32_t a, std::uint32_t b) {
        return krisite::geom::lex_less(m.vertices[a], m.vertices[b]);
    });
    std::vector<std::uint32_t> id(m.vertices.size(), 0);
    std::uint32_t next = 0;
    for (std::size_t i = 0; i < ord.size();) {
        std::size_t j = i;
        while (j < ord.size() && krisite::geom::h_equal(m.vertices[ord[i]], m.vertices[ord[j]])) {
            id[ord[j]] = next;
            ++j;
        }
        ++next;
        i = j;
    }
    std::vector<CanonTri> out;
    out.reserve(m.triangles.size());
    for (const krisite::mesh::Tri& t : m.triangles) {
        const std::uint32_t v[3] = {id[t[0]], id[t[1]], id[t[2]]};
        const int k = (v[0] <= v[1] && v[0] <= v[2]) ? 0 : (v[1] <= v[2] ? 1 : 2);
        out.push_back(CanonTri{v[k], v[(k + 1) % 3], v[(k + 2) % 3]});
    }
    std::sort(out.begin(), out.end());
    return out;
}

bool same_triangles(const BoolMesh& a, const BoolMesh& b) {
    return canonical(a) == canonical(b);
}

/// §10.2.3: 出力の各三角形が、入力側の平面のいずれかに載っていること。
///
/// **余計な幾何が混入していないことを安く確かめられます。** 三角形の平面を
/// 有理点から作り直す代わりに、「3 頂点すべてが同じ入力平面上にある」ことを見ます。
bool planes_belong_to_input(const BoolMesh& out, const TriMesh& src) {
    PlaneTable table;
    build_faces(src, 0, table);
    for (const krisite::mesh::Tri& t : out.triangles) {
        bool found = false;
        for (PlaneId q = 0; q < static_cast<PlaneId>(table.size()) && !found; ++q) {
            const krisite::geom::PlaneD& pl = table.at(q);
            found = krisite::geom::side(pl, out.vertices[t[0]]) == 0 &&
                    krisite::geom::side(pl, out.vertices[t[1]]) == 0 &&
                    krisite::geom::side(pl, out.vertices[t[2]]) == 0;
        }
        if (!found) return false;
    }
    return true;
}

void run(const kritest::Case& c) {
    const TriMesh a = c.make_a();
    const BoolMesh re = retriangulate(a);
    std::printf("\n  ケース %-4s %s（retriangulate(A) の面数 = %zu）\n", c.id, c.what,
                re.triangles.size());

    for (unsigned d = 0; d <= 3; ++d) {
        BoolStats si, su, sd;
        const BoolMesh inter =
            boolean_op(a, a, BoolOp::Intersection, kritest::corpus_options(d), &si);
        const BoolMesh uni = boolean_op(a, a, BoolOp::Union, kritest::corpus_options(d), &su);
        const BoolMesh diff = boolean_op(a, a, BoolOp::Difference, kritest::corpus_options(d), &sd);
        const std::string tag =
            std::string("ケース ") + c.id + "（深度 " + std::to_string(d) + "）";

        // §10.2.3 分割不変な軽い検査
        KRI_CHECK_MSG(diff.triangles.empty(), tag + ": A\\A が空でない");
        KRI_CHECK_MSG(planes_belong_to_input(inter, a),
                      tag + ": A∩A の面が A の平面に載っていない（余計な幾何の混入）");
        KRI_CHECK_MSG(planes_belong_to_input(uni, a),
                      tag + ": A∪A の面が A の平面に載っていない（余計な幾何の混入）");

        // §10.2.2 冪等性
        bool ok_i = false, ok_u = false;
        if (d == 0) {
            ok_i = same_triangles(inter, re);
            ok_u = same_triangles(uni, re);
            KRI_CHECK_MSG(ok_i, tag + ": A∩A が retriangulate(A) と一致しない");
            KRI_CHECK_MSG(ok_u, tag + ": A∪A が retriangulate(A) と一致しない");
        } else {
            ok_i = ok_u = same_triangles(inter, uni);
            KRI_CHECK_MSG(ok_i, tag +
                                    ": A∩A と A∪A が互いに一致しない（"
                                    "共平面重複の選択規則が食い違っている）");
        }
        std::printf("    d%u  ∩:F=%-5zu ∪:F=%-5zu \\:F=%-3zu  %s\n", d, inter.triangles.size(),
                    uni.triangles.size(), diff.triangles.size(),
                    d == 0 ? (ok_i && ok_u ? "= retriangulate(A)" : "**NG**")
                           : (ok_i ? "∩ = ∪" : "**NG**"));
    }
}

}  // namespace

int main() {
    std::printf("\n  §10.2.2 冪等性 / §10.2.3 平面 ID 帰属\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) run(c);
    std::printf("\n");
    return kritest::finish("csg/idempotence");
}

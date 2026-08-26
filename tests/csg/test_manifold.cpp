// Krisite — 外部正解器 Manifold との突き合わせ（SPEC-phase1.md §10.4）
//
// **正解器は被検体と別経路で書く**（CLAUDE.md「検証」）。§10.3 の体積恒等式は
// 自己整合の検査なので、系統的な誤りが相関して入ると素通りし得ます。Manifold は
// 実装がまったく別なので、**答えのレベルで独立**です。
//
// Manifold は Apache-2.0。**テスト専用で、ライブラリ本体からは絶対に呼びません。**
// 取得は CMake FetchContent で固定タグ（`KRISITE_MANIFOLD_TAG`）。
//
// **比較するのは種数と連結成分数です。** Manifold は浮動小数点なので体積は許容誤差
// つきの比較になりますが、位相量は厳密に一致すべきです（§10.4）。
// 入力座標は整数で $|c| \le 2^{25}$ なので `double` に**誤差なく**載ります。
//
// §9.3 の除外（頂点・辺接触の $\cup$）はここでも適用します。非多様体を出す配置では
// Manifold が別の意味論（多様体化）を採るため、比較しても意味がありません。
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "corpus_expect.hpp"
#include "manifold/manifold.h"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;

namespace {

const char* op_name(BoolOp op) {
    switch (op) {
        case BoolOp::Union:
            return "∪";
        case BoolOp::Intersection:
            return "∩";
        default:
            return "\\";
    }
}

/// 入力メッシュを Manifold に渡す。整数座標なので double に誤差なく載る。
manifold::Manifold to_manifold(const TriMesh& m) {
    manifold::MeshGL64 gl;
    gl.numProp = 3;
    gl.vertProperties.reserve(m.vertices.size() * 3);
    for (const krisite::geom::IPoint& p : m.vertices) {
        gl.vertProperties.push_back(static_cast<double>(p.x));
        gl.vertProperties.push_back(static_cast<double>(p.y));
        gl.vertProperties.push_back(static_cast<double>(p.z));
    }
    gl.triVerts.reserve(m.triangles.size() * 3);
    for (const krisite::mesh::Tri& t : m.triangles) {
        gl.triVerts.push_back(t[0]);
        gl.triVerts.push_back(t[1]);
        gl.triVerts.push_back(t[2]);
    }
    return manifold::Manifold(gl);
}

/// 連結成分数と、成分ごとの種数の総和。
///
/// `Genus()` は単一の立体についてのみ意味を持つので、先に `Decompose()` します。
/// `genus_total = C - χ/2` は非交和で $\sum g_i$ に一致するので、この形で比べられます。
struct Topo {
    std::size_t components = 0;
    long long genus_total = 0;
};

Topo topo_of(const manifold::Manifold& m) {
    Topo r;
    if (m.IsEmpty()) return r;
    const std::vector<manifold::Manifold> parts = m.Decompose();
    r.components = parts.size();
    for (const manifold::Manifold& p : parts) r.genus_total += p.Genus();
    return r;
}

/// §9.3: 孤立した頂点・辺接触は合否判定から除外します。
/// **定義は `corpus_expect.hpp` の 1 箇所だけです。**
bool excluded(const std::string& id, BoolOp op) {
    return kritest::exclusion_of(id, op) != kritest::Exclusion::None;
}

/// 実際に比較した (ケース, 演算, 深度) の組数。
///
/// **「テストが通った」ことと「テストが何かを検証した」ことは別です。**
/// Manifold ジョブは一度、テストが 1 件も登録されないまま緑になりました。
/// 件数を数えて下限を assert しておけば、同じ形の空回りは中身でも検出できます。
std::size_t compared = 0;

void check_case(const kritest::Case& c) {
    const TriMesh A = c.make_a(), B = c.make_b();
    const manifold::Manifold ma = to_manifold(A), mb = to_manifold(B);
    KRI_CHECK_MSG(ma.Status() == manifold::Manifold::Error::NoError,
                  std::string("ケース ") + c.id + ": Manifold が A を受け付けない");
    KRI_CHECK_MSG(mb.Status() == manifold::Manifold::Error::NoError,
                  std::string("ケース ") + c.id + ": Manifold が B を受け付けない");

    std::printf("\n  ケース %-3s %s\n", c.id, c.what);
    for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
        const manifold::OpType mop = (op == BoolOp::Union)          ? manifold::OpType::Add
                                     : (op == BoolOp::Intersection) ? manifold::OpType::Intersect
                                                                    : manifold::OpType::Subtract;
        const Topo want = topo_of(ma.Boolean(mb, mop));

        for (unsigned d = 0; d <= 3; ++d) {
            BoolStats st;
            const BoolMesh r = boolean_op(A, B, op, d, &st);
            const TopologyReport t = check_topology(r.triangles);
            const bool skip = excluded(c.id, op);
            const bool agree =
                (t.components == want.components) && (t.genus_total == want.genus_total);
            std::printf("    d%u %s  Krisite C=%zu g=%-2lld / Manifold C=%zu g=%-2lld  %s\n", d,
                        op_name(op), t.components, t.genus_total, want.components, want.genus_total,
                        skip ? "（§9.3 で合否対象外）" : (agree ? "一致" : "**不一致**"));
            if (skip) continue;
            ++compared;
            KRI_CHECK_MSG(t.components == want.components,
                          std::string("ケース ") + c.id + " " + op_name(op) + "（深度 " +
                              std::to_string(d) + "）: 連結成分数が Manifold と食い違う（" +
                              std::to_string(t.components) + " 対 " +
                              std::to_string(want.components) + "）");
            KRI_CHECK_MSG(t.genus_total == want.genus_total,
                          std::string("ケース ") + c.id + " " + op_name(op) + "（深度 " +
                              std::to_string(d) + "）: 種数が Manifold と食い違う（" +
                              std::to_string(t.genus_total) + " 対 " +
                              std::to_string(want.genus_total) + "）");
        }
    }
}

}  // namespace

int main() {
    std::printf("\n  §10.4 外部正解器（Manifold）との突き合わせ\n");
    KRI_CHECK_MSG(!kritest::corpus().empty(), "コーパスが空");
    for (const kritest::Case& c : kritest::corpus()) check_case(c);

    // 全ケース × 3 演算 × 深度 4 段から §9.3 の除外（各 4 深度）を引いた数
    std::size_t n_excluded = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            if (excluded(c.id, op)) ++n_excluded;
        }
    }
    const std::size_t expect = kritest::corpus().size() * 3 * 4 - n_excluded * 4;
    std::printf("\n  比較した組数: %zu（期待 %zu）\n", compared, expect);
    KRI_CHECK_MSG(compared == expect, "比較した組数が期待と違う（" + std::to_string(compared) +
                                          " 対 " + std::to_string(expect) +
                                          "）。テストが空回りしている疑い");
    return kritest::finish("csg/manifold");
}

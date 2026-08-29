// Krisite — $n$ 項のケース 17 / 18 / 23（`SPEC-phase3.md` §9）
//
// **二項の連鎖とは突き合わせられません。** 二項の `boolean_op` は `BoolMesh` を返し、
// 構成点が有理数なので**入力に戻せない**からです。それが $n$ 項にした理由でもあります。
//
// そこで**別の演算木で同じ立体を作る**のを正解器にします。
//
//     ケース 17   $(A \cup B) \setminus C \;=\; (A \setminus C) \cup (B \setminus C)$
//     ケース 18   $W - T - T - T - T - T \;=\; W - T$
//     ケース 23   $\bigcup_{i=0}^{19} M_i \;=\;$ 解析的に分かる 1 つの箱
//
// **集合の恒等式なので、答えのレベルで独立です**（自己整合の検査ではありません）。
//
// ---
//
// ## 何をもって「同じ立体」とするか
//
// 演算木が違えば**分割も三角形の集合も違います**（それが目的の一部）。そこで
//
// | 検査 | 見るもの |
// |---|---|
// | $(C, \chi)$ | 位相が同じか |
// | **$X \setminus Y$ と $Y \setminus X$ が空** | **点集合として同じか** |
// | $X$ が空でない | **空回りの番人**（両方空なら差も空で素通りする） |
//
// **差が両方空であることは、正則化ブールの意味で $X = Y$ と同値です。**
// これは「同じ機械で確かめている」ことになりますが、**演算木が違う**ので
// 同じ間違いをする経路にはなっていません。体積の突き合わせは `volume_gmp` の担当です。
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kMaxDepth = 3;
constexpr std::size_t kPasses = kMaxDepth + 2;  ///< 固定深度 0〜3 + 適応分割

BoolOptions all_on(unsigned depth, bool adaptive) {
    BoolOptions o;
    o.depth = depth;
    o.cull_planes = true;
    o.adaptive = adaptive;
    o.leaf_threshold = 0;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    return o;
}

std::size_t g_checks = 0;
std::size_t g_empty_diffs = 0;

/// $X$ と $Y$ が正則化ブールの意味で等しいこと。
void check_same(const PolySoup& x, const PolySoup& y, const BoolOptions& o,
                const std::string& tag) {
    ToMeshOptions tm;
    tm.split_contacts = true;
    const TopologyReport rx = check_topology(to_mesh(x, tm).triangles);
    const TopologyReport ry = check_topology(to_mesh(y, tm).triangles);

    // **空回りの番人。** 両方空なら差も空になり、以下の検査は素通りします
    KRI_CHECK_MSG(!rx.empty, tag + ": **左辺が空**。恒等式の検査が素通りします");
    KRI_CHECK_MSG(rx.ok(), tag + ": 左辺が多様体でない");
    KRI_CHECK_MSG(ry.ok(), tag + ": 右辺が多様体でない");
    KRI_CHECK_MSG(rx.components == ry.components,
                  tag + ": C が違う" + kritest::pair_msg(rx.components, ry.components));
    KRI_CHECK_MSG(rx.chi == ry.chi, tag + ": χ が違う" + kritest::pair_msg(rx.chi, ry.chi));
    // **意図した形になっているかをログに残します。** 「通った」だけでは、
    // コーパスが狙った退化を作れているかが分かりません
    if (g_checks == 0 || tag.find("深度 0") != std::string::npos) {
        std::printf("      %s: C=%zu χ=%lld（種数 %lld）三角形 %zu\n", tag.c_str(), rx.components,
                    rx.chi, (2 * static_cast<long long>(rx.components) - rx.chi) / 2, rx.f);
    }

    // **点集合として同じか。** 位相が同じでも別の立体でありえます
    for (int dir = 0; dir < 2; ++dir) {
        const PolySoup& p = (dir == 0) ? x : y;
        const PolySoup& q = (dir == 0) ? y : x;
        const SoupMesh d = to_mesh(boolean(p, q, BoolOp::Difference, o), tm);
        KRI_CHECK_MSG(d.triangles.empty(),
                      tag + (dir == 0 ? ": 左 \\ 右 が空でない" : ": 右 \\ 左 が空でない") +
                          "（三角形 " + std::to_string(d.triangles.size()) + " 枚）");
        ++g_empty_diffs;
    }
    ++g_checks;
}

/// ケース 17: 分配法則。
void run_distributive(const kritest::NaryCase& c, const BoolOptions& o, const std::string& tag) {
    const std::vector<TriMesh> m = c.make();
    KRI_CHECK_MSG(m.size() == 3, tag + ": メッシュが 3 枚でない");
    const PolySoup a = from_mesh(m[0]), b = from_mesh(m[1]), d = from_mesh(m[2]);

    // 左辺: (A ∪ B) \ C
    const PolySoup lhs = boolean(boolean(a, b, BoolOp::Union, o), d, BoolOp::Difference, o);
    // 右辺: (A \ C) ∪ (B \ C)。**演算木の形が違います**
    const PolySoup rhs = boolean(boolean(a, d, BoolOp::Difference, o),
                                 boolean(b, d, BoolOp::Difference, o), BoolOp::Union, o);

    KRI_CHECK_MSG(lhs.source_count() == 3, tag + ": 左辺の source 数が 3 でない");
    // **右辺は C を 2 度使うので source が 4 つになります。** それでも同じ立体です
    // （同じ曲面を 2 度跨ぐ配置。内外の 1 ビットでは表せません）
    KRI_CHECK_MSG(rhs.source_count() == 4, tag + ": 右辺の source 数が 4 でない");
    check_same(lhs, rhs, o, tag);
}

/// ケース 18: 同じ形状を繰り返し引く。
void run_repeated(const kritest::NaryCase& c, const BoolOptions& o, const std::string& tag) {
    const std::vector<TriMesh> m = c.make();
    KRI_CHECK_MSG(m.size() == 2, tag + ": メッシュが 2 枚でない");
    const PolySoup w = from_mesh(m[0]);

    // 左辺: W - T - T - ... - T（**毎回 from_mesh するので source が増えます**）
    PolySoup lhs = w;
    for (int i = 0; i < c.repeat; ++i) {
        lhs = boolean(lhs, from_mesh(m[1]), BoolOp::Difference, o);
    }
    const PolySoup rhs = boolean(w, from_mesh(m[1]), BoolOp::Difference, o);

    KRI_CHECK_MSG(lhs.source_count() == static_cast<std::size_t>(c.repeat) + 1,
                  tag + ": 左辺の source 数が repeat + 1 でない");
    KRI_CHECK_MSG(rhs.source_count() == 2, tag + ": 右辺の source 数が 2 でない");
    check_same(lhs, rhs, o, tag);
}

/// ケース 23: 多数の和が、解析的に分かる立体と一致すること。
void run_union_of_many(const kritest::NaryCase& c, const BoolOptions& o, const std::string& tag) {
    const std::vector<TriMesh> m = c.make();
    KRI_CHECK_MSG(m.size() >= 2, tag + ": メッシュが 2 枚未満");
    KRI_CHECK_MSG(c.expect != nullptr, tag + ": 期待値が無い");

    PolySoup lhs = from_mesh(m[0]);
    for (std::size_t i = 1; i < m.size(); ++i) {
        lhs = boolean(lhs, from_mesh(m[i]), BoolOp::Union, o);
    }
    KRI_CHECK_MSG(lhs.source_count() == m.size(), tag + ": source 数がメッシュ数と合わない");

    // **期待値は解析的に分かる立体です。** 恒等式は自己整合の検査なので、
    // 答えのレベルで独立な対照を置きます（`SPEC-phase1.md` §10.3 と同じ構図）。
    check_same(lhs, from_mesh(c.expect()), o, tag);
}

/// ケース 19 / 20a: **1 つのメッシュとして与えた 2 つの閉曲面**が、
/// 別々に与えて和を取った結果と一致すること。
///
/// **入力は PWN です**（各閉曲面が個別に $\partial = 0$）。多様体である必要は
/// ありません。**受け入れられることと、正しく解けることの両方**を見ます。
void run_self_union(const kritest::NaryCase& c, const BoolOptions& o, const std::string& tag) {
    const std::vector<TriMesh> m = c.make();
    KRI_CHECK_MSG(m.size() == 3, tag + ": メッシュが 3 枚でない");
    KRI_CHECK_MSG(krisite::mesh::boundary_is_zero(m[0]), tag + ": 入力が PWN でない（∂S ≠ 0）");

    // self-union: 同じスープを 2 度使う和 = 「$w \ne 0$」の領域
    const PolySoup self = from_mesh(m[0]);
    const PolySoup lhs = boolean(self, from_mesh(m[0]), BoolOp::Union, o);
    const PolySoup rhs = boolean(from_mesh(m[1]), from_mesh(m[2]), BoolOp::Union, o);
    check_same(lhs, rhs, o, tag);
}

void run_case(const kritest::NaryCase& c) {
    // §9.0 (1) のサイズ規律
    KRI_CHECK_MSG(kritest::size_discipline_ok(c.make()),
                  std::string("ケース ") + c.id + ": サイズ規律を満たしていません（§9.0）");
    for (std::size_t pass = 0; pass < kPasses; ++pass) {
        const bool adaptive = (pass == kPasses - 1);
        const unsigned depth = adaptive ? kMaxDepth : static_cast<unsigned>(pass);
        const std::string tag = std::string("ケース ") + c.id + "（" +
                                (adaptive ? "適応" : "深度 " + std::to_string(depth)) + "）";
        const BoolOptions o = all_on(depth, adaptive);
        switch (c.kind) {
            case kritest::NaryCase::Kind::Distributive:
                run_distributive(c, o, tag);
                break;
            case kritest::NaryCase::Kind::RepeatedDifference:
                run_repeated(c, o, tag);
                break;
            case kritest::NaryCase::Kind::UnionOfMany:
                run_union_of_many(c, o, tag);
                break;
            case kritest::NaryCase::Kind::SelfUnion:
                run_self_union(c, o, tag);
                break;
        }
    }
}

/// **入口の契約の検査**（`SPEC-phase3.md` §4.0）。
///
/// **受け入れる例と拒否する例を対で置きます**（`CLAUDE.md`「前提を書いたら、検査も
/// 書いてください。…受け入れる例と拒否する例の**対**で」）。片方だけでは
///
/// - 受け入れだけ → 「どこまで受け入れるか」が分からない
/// - 拒否だけ → 「非多様体を本当に受け入れられるか」が分からない
///
/// `from_mesh` は表明で止まりますが、**表明は検査ビルドでしか効きません。**
/// 契約の窓口は `mesh::boundary_is_zero` のほうで、ここではその述語を直接見ます。
void test_pwn_contract() {
    struct Item {
        const char* id;
        const char* what;
        TriMesh mesh;
        bool want_pwn;
        bool want_edge_manifold;
    };
    const std::vector<Item> items = {
        {"20a", "辺を共有した 2 つの閉曲面", kritest::cases::case20a()[0], true, false},
        {"20b", "1 辺に 3 枚", kritest::cases::case20b(), false, false},
        {"19", "自己交差した閉曲面", kritest::cases::case19()[0], true, true},
    };
    for (const Item& it : items) {
        const bool pwn = krisite::mesh::boundary_is_zero(it.mesh);
        const TopologyReport r = check_topology(it.mesh.triangles);
        KRI_CHECK_MSG(pwn == it.want_pwn, std::string("ケース ") + it.id + "（" + it.what +
                                              "）: PWN の判定が想定と違う");
        KRI_CHECK_MSG(r.edge_manifold == it.want_edge_manifold,
                      std::string("ケース ") + it.id + ": 辺多様体の判定が想定と違う");
        std::printf("    ケース %-3s %-24s PWN=%d 辺多様体=%d\n", it.id, it.what, (int)pwn,
                    (int)r.edge_manifold);
    }
    // **非多様体でも PWN なら受け入れる**ことが要点です。両者が一致するなら、
    // `edge_manifold` で代用できてしまい、§4.0 の検査を足す意味がありません
    KRI_CHECK_MSG(krisite::mesh::boundary_is_zero(kritest::cases::case20a()[0]) &&
                      !check_topology(kritest::cases::case20a()[0].triangles).edge_manifold,
                  "**PWN かつ非多様体の例が無い。** §4.0 の検査が `edge_manifold` と"
                  "区別できていません");
}

}  // namespace

int main() {
    std::printf("\n  n 項のケース 17 / 18 / 19 / 20 / 23（SPEC-phase3 §9）\n");
    KRI_CHECK_MSG(!kritest::nary_corpus().empty(), "n 項のコーパスが空");
    test_pwn_contract();
    for (const kritest::NaryCase& c : kritest::nary_corpus()) {
        std::printf("    ケース %-3s %s\n", c.id, c.what);
        run_case(c);
    }
    // **期待値は式で持たせます**（`CLAUDE.md`）
    const std::size_t want = kritest::nary_corpus().size() * kPasses;
    KRI_CHECK_MSG(g_checks == want,
                  "恒等式の検査数が式と合わない" + kritest::pair_msg(want, g_checks));
    KRI_CHECK_MSG(g_empty_diffs == want * 2, "差の検査数が検査数の 2 倍でない");
    std::printf("    恒等式 %zu 件（両方向の差 %zu 件が空）\n", g_checks, g_empty_diffs);
    std::printf("\n");
    return kritest::finish("csg/nary");
}

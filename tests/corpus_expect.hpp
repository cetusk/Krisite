// Krisite — コーパスの期待値と §9.3 の除外表
//
// **除外の定義を 1 箇所に集めます。** CP2.5 / CP3 / Manifold の 3 つのテストが同じ
// 判断をする必要があり、別々に書くと必ずずれます。
//
// SPEC-phase1.md §9.2（ケース 10 の期待値）, §9.3（除外）, §9.3.1（除外の適用条件）
#ifndef KRISITE_TESTS_CORPUS_EXPECT_HPP
#define KRISITE_TESTS_CORPUS_EXPECT_HPP

#include <string>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"

namespace kritest {

/// §9.3 の除外の種類。
///
/// **仕様は「頂点多様体の検査を外す」と書いていますが、それでは足りません**
/// （`IMPL-phase1.md` §7.9 で仕様の訂正をお願いしています）。
///
/// | 接触 | 壊れる検査 | 例 |
/// |---|---|---|
/// | 頂点接触 | `vertex_manifold`（扇が $k$ 個できる） | 4T の $\cup$、4T′ の $A\setminus B$ |
/// | **辺接触** | **`edge_manifold`**（辺に 4 面が接する） | **11b の $\cup$** |
///
/// §9.3.1 の適用条件は「すべての辺の次数が 2」を課していますが、**辺接触では
/// その条件自体が成り立ちません。** 条件が想定しているのは頂点接触だけです。
enum class Exclusion {
    None,
    VertexContact,  ///< `vertex_manifold` のみ外す。辺の次数は 2 のままを要求
    EdgeContact,    ///< `edge_manifold` のみ外す。次数 4 の辺だけを許す
};

/// §9.3 の除外表。**ここだけが定義です。**
inline Exclusion exclusion_of(const std::string& id, krisite::csg::BoolOp op) {
    using krisite::csg::BoolOp;
    if (id == "4T" && op == BoolOp::Union) return Exclusion::VertexContact;
    if (id == "4T'" && op == BoolOp::Difference) return Exclusion::VertexContact;
    if (id == "11b" && op == BoolOp::Union) return Exclusion::EdgeContact;
    return Exclusion::None;
}

/// §9.3.1 の適用条件を満たすか。**除外を広げすぎないための番人です。**
///
/// 除外するのは 1 つの検査だけで、**他はすべて拘束のまま**であることを確かめます。
/// 満たさないなら、それは除外すべき退化ではなくバグです。
inline bool exclusion_conditions_ok(Exclusion e, const krisite::mesh::TopologyReport& t,
                                    std::string* why) {
    if (!t.no_degenerate) {
        *why = "退化三角形がある";
        return false;
    }
    if (t.edges_deficient != 0) {
        *why = "面が欠けている辺がある（次数 1）";
        return false;
    }
    if (e == Exclusion::VertexContact) {
        if (!t.oriented) {
            *why = "向きが整合しない";
            return false;
        }
        if (t.edges_excess != 0) {
            *why = "辺の次数が 2 でない（頂点接触なら面の過不足は無いはず）";
            return false;
        }
        if (t.vertex_manifold) {
            *why = "頂点多様体が通っている（除外が不要）";
            return false;
        }
        return true;
    }
    if (e == Exclusion::EdgeContact) {
        // **`oriented` は使えません。** 次数 2 を前提にした検査なので、辺に 4 枚の面が
        // 接すると向きが正しくても必ず false になります。次数 2 の辺に限った検査と、
        // 接触辺で両向きの本数が釣り合っていることを見ます。
        if (!t.oriented_on_manifold_edges) {
            *why = "次数 2 の辺で向きが整合しない";
            return false;
        }
        if (!t.excess_edges_balanced) {
            *why = "接触辺で両向きの本数が釣り合わない";
            return false;
        }
        if (t.edges_excess == 0) {
            *why = "次数 3 以上の辺が無い（除外が不要）";
            return false;
        }
        if (t.max_edge_degree != 4) {
            *why = "接触辺の次数が 4 でない（2 つの立体の辺接触なら 4 のはず）";
            return false;
        }
        // **頂点多様体は落ちます。** 接触辺の上の頂点では 4 枚の面が出会うので、
        // リンクが単純閉路の直和になりません。要求すべきなのは
        // 「非多様体な頂点がすべて接触辺の上にあること」です。
        if (t.nonmanifold_vertices_off_excess != 0) {
            *why = "接触辺から離れた場所に非多様体な頂点がある";
            return false;
        }
        return true;
    }
    *why = "除外していない";
    return false;
}

/// 期待する連結成分数と種数。分からないケースは `known = false`。
struct ExpectedTopo {
    bool known = false;
    bool empty = false;
    std::size_t components = 0;
    long long genus_total = 0;
};

/// §9.2 と §9.1 が明示している期待値。
///
/// **導けるものだけを書きます。** 分からないものを推測で書くと、実装の誤りを
/// 期待値のほうに合わせてしまいます。
inline ExpectedTopo expected_topo(const std::string& id, krisite::csg::BoolOp op) {
    using krisite::csg::BoolOp;
    // ケース 10: B が A の内部に完全に含まれる（§9.2 の表）
    if (id == "10") {
        if (op == BoolOp::Union) return {true, false, 1, 0};         // A
        if (op == BoolOp::Intersection) return {true, false, 1, 0};  // B。**空ではない**
        return {true, false, 2, 0};  // A の内部に空洞。外殻 + 内殻で χ = 4
    }
    // ケース 8: A ∩ A = A、A \ A = 空（§9.1）
    if (id == "8") {
        if (op == BoolOp::Difference) return {true, true, 0, 0};
        return {true, false, 1, 0};
    }
    // ケース 11a / 11b: 接触のみなので ∩ は空、A \ B = A（§9.1、§9.3）
    if (id == "11a" || id == "11b") {
        if (op == BoolOp::Intersection) return {true, true, 0, 0};
        if (op == BoolOp::Difference) return {true, false, 1, 0};
        if (id == "11a") return {true, false, 1, 0};  // 接触面が内部面になって消える
        return {};                                    // 11b の ∪ は非多様体（§9.3）
    }
    // ケース 4T: 頂点接触なので ∩ は空、A \ B = A
    if (id == "4T") {
        if (op == BoolOp::Intersection) return {true, true, 0, 0};
        if (op == BoolOp::Difference) return {true, false, 1, 0};
        return {};
    }
    // 交わらない 2 立方体
    if (id == "D") {
        if (op == BoolOp::Union) return {true, false, 2, 0};
        if (op == BoolOp::Intersection) return {true, true, 0, 0};
        return {true, false, 1, 0};
    }
    return {};
}

}  // namespace kritest

#endif  // KRISITE_TESTS_CORPUS_EXPECT_HPP

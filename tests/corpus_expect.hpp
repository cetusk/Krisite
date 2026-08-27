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

/// §9.3 の除外 — **接触の次元**（SPEC-phase1 §9.3.0、第11版）。
///
/// | 次元 | 壊れる検査 | 例 |
/// |---|---|---|
/// | **0（頂点）** | `vertex_manifold`（リンクが $k$ 個の扇） | 4T ∪、4T′ $A\setminus B$ |
/// | **1（辺）** | `edge_manifold`（辺に 4 枚の面）。深度 ≥1 では `vertex_manifold` も | 11b ∪ |
///
/// **`oriented` は除外しません**（§9.3.0.1）。検査を「各辺で $\#(u,v)=\#(v,u)$」に
/// 一般化したので、接触辺でも正しく機能します。
///
/// 次元は**§13 の「除外件数を接触の次元別に記録」のためのラベル**です。
/// 適用条件（`exclusion_conditions_ok`）は次元に依りません。
enum class Exclusion {
    None,
    VertexContact,  ///< 接触の次元 0
    EdgeContact,    ///< 接触の次元 1
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
/// **接触の種類に依らない形で書きます**（SPEC-phase1 §9.3.1、第11版）。
/// 「どの辺が接触辺か」を独立に知る手段が無いので、**辺次数の性質だけ**で判定します。
///
/// | 条件 | 根拠 |
/// |---|---|
/// | すべての辺の次数が偶数 | 奇数次数は曲面に境界がある証拠。**次数 1 を含みます** |
/// | 次数の最大が 4 以下 | 立体 2 個の辺接触なら 2+2 で 4。**5 以上は別の異常** |
/// | 各辺で $\#(u,v) = \#(v,u)$ | 一般化した `oriented`（§9.3.0.1）。**除外しません** |
/// | 非多様体な頂点がすべて説明できる | 次数 3 以上の辺に接するか $k \ge 2$ 個の扇 |
///
/// §10.3 の体積検査（恒等式・解析値）も条件のひとつですが、GMP が要るので
/// `csg/test_volume_gmp.cpp` が別に受け持ちます。**あちらは除外を適用しません。**
inline bool exclusion_conditions_ok(Exclusion e, const krisite::mesh::TopologyReport& t,
                                    std::string* why) {
    if (e == Exclusion::None) {
        *why = "除外していない";
        return false;
    }
    if (!t.no_degenerate) {
        *why = "退化三角形がある";
        return false;
    }
    if (t.edges_odd_degree != 0) {
        *why = "次数が奇数の辺がある（曲面に境界がある = 面の過不足）";
        return false;
    }
    if (t.max_edge_degree > 4) {
        *why = "辺の次数の最大が " + std::to_string(t.max_edge_degree) +
               "（5 以上は接触では説明できない異常）";
        return false;
    }
    if (!t.oriented) {
        *why = "向きが整合しない（各辺で #(u,v) = #(v,u) が崩れている）";
        return false;
    }
    if (t.nonmanifold_vertices_unexplained != 0) {
        *why = "接触で説明できない非多様体頂点がある";
        return false;
    }
    // **除外が実際に必要であることも確かめます。** 通っているのに除外していたら、
    // 除外表のほうが古くなっています。
    if (t.edge_manifold && t.vertex_manifold) {
        *why = "多様体検査が通っている（除外が不要）";
        return false;
    }
    return true;
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
    // SPEC-phase2 §8 の追加ケース。**導けるものは書きます**（§10.3 の但し書きと同じ規律）。
    //
    // 13 / 15: B の各部品が A の面をまたぐので、∪ は A に「こぶ」が付いた球、
    //          ∩ は部品の内側だけが残って部品数ぶんの成分、A∖B は A に「くぼみ」。
    //          くぼみは空洞ではないので成分は増えません。
    if (id == "13" || id == "15") {
        const std::size_t parts = (id == "13") ? 4 : 2;
        if (op == BoolOp::Union) return {true, false, 1, 0};
        if (op == BoolOp::Intersection) return {true, false, parts, 0};
        return {true, false, 1, 0};
    }
    // 14: 離れているのでケース D と同じ形
    if (id == "14") {
        if (op == BoolOp::Union) return {true, false, 2, 0};
        if (op == BoolOp::Intersection) return {true, true, 0, 0};
        return {true, false, 1, 0};
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

// Krisite — 4 出力の厳密体積整合性（テスト専用。LGPL の GMP に依存）
//
// **`DESIGN-phase5-vertex-level.md` §22.16 の式 1 を実装します。**
//
// $$ |A \cup B| = |A \setminus B| + |B \setminus A| + |A \cap B| $$
//
// **4 つの出力だけで閉じている**ので、生入力の巻き数が 0/1 であることを仮定しません。
//
// > **★ これは独立した正解器ではありません**（§22.16）。
// > **4 出力が全部空でも成立し、共通の欠落や誤差の相殺も検出しません。**
// > **記録は「4 出力の厳密体積整合性が一致／不一致」と書いてください。**
// > **一致を、入力に対する正しさ全体の証拠にしないこと。**
//
// **ライブラリ本体からは絶対に呼びません**（`CLAUDE.md` のライセンス方針）。
#ifndef KRISITE_TESTS_GMP_IDENTITY_HPP
#define KRISITE_TESTS_GMP_IDENTITY_HPP

#include <gmp.h>

#include "volume_gmp.hpp"

namespace kritest {

/// 4 出力の 6 倍体積を入れる器。**添字の意味を型で固定できないので、定数で示します。**
enum : int { kOpUnion = 0, kOpIsect = 1, kOpDiffAB = 2, kOpDiffBA = 3, kOpCount = 4 };

/// 式 1 の残差 $|A\cup B| - (|A\setminus B| + |B\setminus A| + |A\cap B|)$ を
/// **6 倍体積のまま** `res` に入れ、0 なら真を返します。
///
/// **`v6` は 4 つの出力の 6 倍体積**（`kOp*` の順）。**呼び出し側が用意します。**
inline bool identity1_residual(mpq_ptr res, const mpq_t v6[kOpCount]) {
    mpq_t rhs;
    mpq_init(rhs);
    mpq_set(rhs, v6[kOpDiffAB]);
    mpq_add(rhs, rhs, v6[kOpDiffBA]);
    mpq_add(rhs, rhs, v6[kOpIsect]);
    mpq_sub(res, v6[kOpUnion], rhs);
    mpq_clear(rhs);
    return mpq_sgn(res) == 0;
}

/// 式 2 の残差 $(|A\cup B| + |A\cap B|) - (|A| + |B|)$。
///
/// > **★ 前提つきです**（§22.16）。**$w_A, w_B$ が 0 か 1 しか取らないこと。**
/// > **前提を独立に確かめる手段はまだありません。**
/// > **確かめられないうちは、破れを欠陥の根拠にしないでください。**
inline bool identity2_residual(mpq_ptr res, const mpq_t v6[kOpCount], mpq_srcptr va6,
                               mpq_srcptr vb6) {
    mpq_t lhs, rhs;
    mpq_init(lhs);
    mpq_init(rhs);
    mpq_set(lhs, v6[kOpUnion]);
    mpq_add(lhs, lhs, v6[kOpIsect]);
    mpq_set(rhs, va6);
    mpq_add(rhs, rhs, vb6);
    mpq_sub(res, lhs, rhs);
    mpq_clear(lhs);
    mpq_clear(rhs);
    return mpq_sgn(res) == 0;
}

}  // namespace kritest

#endif  // KRISITE_TESTS_GMP_IDENTITY_HPP

// Krisite — 4 出力の厳密体積整合性の【検査器】を、合成で検定する
//
// **`DESIGN-phase5-vertex-level.md` §22.17。**
//
// **実データに掛ける前に、検査器が破れを検出できることを確かめます**
// （`CLAUDE.md`「判定器を先に検定した形が要点」）。
//
// 構成: $A = [0,4]^3$、$B = [2,6] \times [0,4] \times [0,4]$
//
// | 判定対象 | 体積 | 6 倍体積 |
// |---|---:|---:|
// | ∪ | 96 | **576** |
// | ∩ / A∖B / B∖A | 各 32 | **各 192** |
// | N1 の残差（∩ を ∪ に差し替え） | −64 | **−384** |
// | N2 の残差（A∖B を空に） | ＋32 | **＋192** |
//
// **関係式だけに頼りません。各出力の解析値とも照合します。**
#include <cstdio>
#include <cstring>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"

#include "corpus.hpp"
#include "gmp_identity.hpp"

using namespace krisite;
using namespace krisite::csg;

namespace {

int failures = 0;

/// `mpq` が整数 `want` と一致するか。**丸めも許容差もありません。**
void expect_q(mpq_srcptr got, long want, const char* what) {
    mpq_t w;
    mpq_init(w);
    mpq_set_si(w, want, 1);
    if (mpq_equal(got, w) == 0) {
        char* s = mpq_get_str(nullptr, 10, got);
        std::printf("**不一致** %s: 期待 %ld、実測 %s\n", what, want, s);
        void (*freefunc)(void*, std::size_t);
        mp_get_memory_functions(nullptr, nullptr, &freefunc);
        freefunc(s, std::strlen(s) + 1);
        ++failures;
    } else {
        std::printf("  %s = %ld  一致\n", what, want);
    }
    mpq_clear(w);
}

}  // namespace

int main() {
    std::printf("# 4 出力の厳密体積整合性 — 合成の検定\n\n");
    std::printf("| 設定 | 値 |\n|---|---|\n");
    std::printf("| A | [0,4]^3 |\n| B | [2,6]x[0,4]x[0,4] |\n");
    std::printf("| 演算 | 4（∪ / ∩ / A∖B / B∖A）|\n\n");

    const mesh::TriMesh a = kritest::box(0, 0, 0, 4, 4, 4);
    const mesh::TriMesh b = kritest::box(2, 0, 0, 6, 4, 4);
    const PolySoup A = from_mesh(a), B = from_mesh(b);

    BoolOptions o;
    ToMeshOptions tm;
    tm.split_contacts = true;

    // **4 演算**。B∖A は引数を入れ替えて Difference を呼びます
    const SoupMesh out[kritest::kOpCount] = {
        to_mesh(boolean(A, B, BoolOp::Union, o), tm),
        to_mesh(boolean(A, B, BoolOp::Intersection, o), tm),
        to_mesh(boolean(A, B, BoolOp::Difference, o), tm),
        to_mesh(boolean(B, A, BoolOp::Difference, o), tm),
    };
    const char* name[kritest::kOpCount] = {"∪", "∩", "A∖B", "B∖A"};
    const long want6[kritest::kOpCount] = {576, 192, 192, 192};

    mpq_t v6[kritest::kOpCount];
    for (int k = 0; k < kritest::kOpCount; ++k) {
        mpq_init(v6[k]);
        kritest::mesh_volume6(v6[k], out[k]);
    }

    // ---- 1. 各出力の解析値との照合（関係式だけに頼らない）----
    std::printf("## 個別の 6 倍体積\n\n");
    for (int k = 0; k < kritest::kOpCount; ++k) {
        char buf[64];
        std::snprintf(buf, sizeof buf, "%s の 6 倍体積", name[k]);
        expect_q(v6[k], want6[k], buf);
    }

    // ---- 2. 式 1 の整合性 ----
    std::printf("\n## 式 1（4 出力の整合性）\n\n");
    mpq_t res;
    mpq_init(res);
    const bool ok1 = kritest::identity1_residual(res, v6);
    expect_q(res, 0, "正例の残差");
    if (!ok1) std::printf("**正例で整合しませんでした**\n");

    // ---- 3. 負例 N1: ∩ の出力を ∪ の出力に差し替える ----
    std::printf("\n## 負例\n\n");
    {
        mpq_t bad[kritest::kOpCount];
        for (int k = 0; k < kritest::kOpCount; ++k) {
            mpq_init(bad[k]);
            mpq_set(bad[k], v6[k]);
        }
        kritest::mesh_volume6(bad[kritest::kOpIsect], out[kritest::kOpUnion]);  // ★ 差し替え
        const bool consistent = kritest::identity1_residual(res, bad);
        expect_q(res, -384, "N1 の残差（∩ → ∪）");
        if (consistent) {
            std::printf("**N1 を検出できませんでした**\n");
            ++failures;
        }
        for (auto& q : bad) mpq_clear(q);
    }

    // ---- 4. 負例 N2: A∖B の出力を空にする ----
    {
        mpq_t bad[kritest::kOpCount];
        for (int k = 0; k < kritest::kOpCount; ++k) {
            mpq_init(bad[k]);
            mpq_set(bad[k], v6[k]);
        }
        const SoupMesh empty;
        kritest::mesh_volume6(bad[kritest::kOpDiffAB], empty);  // ★ 空に差し替え
        const bool consistent = kritest::identity1_residual(res, bad);
        expect_q(res, 192, "N2 の残差（A∖B を空に）");
        if (consistent) {
            std::printf("**N2 を検出できませんでした**\n");
            ++failures;
        }
        for (auto& q : bad) mpq_clear(q);
    }

    // ---- 5. 式 2（前提つき）。この構成は巻き数が 0/1 なので成り立ちます ----
    std::printf("\n## 式 2（前提つき。この構成は巻き数が 0/1）\n\n");
    {
        mpq_t va, vb;
        mpq_init(va);
        mpq_init(vb);
        kritest::input_volume6(va, a);
        kritest::input_volume6(vb, b);
        expect_q(va, 384, "入力 A の 6 倍体積");
        expect_q(vb, 384, "入力 B の 6 倍体積");
        kritest::identity2_residual(res, v6, va, vb);
        expect_q(res, 0, "式 2 の残差");
        mpq_clear(va);
        mpq_clear(vb);
    }

    mpq_clear(res);
    for (auto& q : v6) mpq_clear(q);
    std::printf("\n**不一致 %d 件**\n", failures);
    return failures == 0 ? 0 : 1;
}

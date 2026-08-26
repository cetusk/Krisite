// Krisite — 体積の恒等式（SPEC-phase1.md §10.3）
//
// **位相検査だけでは「正しい立体か」は分かりません。** 辺が 2 面に接し向きが整合して
// いても、選んだ断片が間違っていれば別の立体になります。$\chi$ は面が余分なのか
// 欠けているのかも区別しません。体積の恒等式はそこを直接突きます。
//
//   |A ∪ B| + |A ∩ B| = |A| + |B|
//   |A \ B|           = |A| - |A ∩ B|
//   |A ∪ B|           = |A| + |B| - |A ∩ B|
//
// **出力の体積は GMP が要ります。** 構成点は有理数で、共通分母が三角形数に比例して
// 伸びるためです（`widths.hpp` bits::kInputVolume6 の注記）。入力側は整数座標なので
// `mesh::signed_volume6` で厳密に出ます。**両者が一致することも検査します**
// （別経路どうしの突き合わせ。CLAUDE.md「正解器は被検体と別経路で書く」）。
//
// GMP は LGPL。テスト専用（`KRISITE_BUILD_TESTS_WITH_GMP=ON` のときだけビルドされる）。
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "gmp_oracle.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using krisite::mesh::TriMesh;
using kritest::oracle::to_mpz;

namespace {

/// 閉じた向き付き三角形メッシュの符号付き体積 x6 = Σ det(v0, v1, v2)。
///
/// 頂点は同次座標 $(X : W)$ なので、三角形 1 枚の寄与は
/// $$ \frac{\det(X_0, X_1, X_2)}{W_0 W_1 W_2}. $$
/// 分子は整数行列式、分母は 3 つの $w$ の積です。**三角形ごとに分母が違う**ので、
/// `mpq` で足し合わせます（ここが GMP を要する理由）。
void bool_mesh_volume6(mpq_ptr out, const BoolMesh& m) {
    mpq_set_ui(out, 0, 1);
    if (m.triangles.empty()) return;

    mpz_t x[3][3], w[3], num, den, t0, t1, cof;
    for (auto& row : x) {
        for (auto& e : row) mpz_init(e);
    }
    for (auto& e : w) mpz_init(e);
    mpz_init(num);
    mpz_init(den);
    mpz_init(t0);
    mpz_init(t1);
    mpz_init(cof);
    mpq_t term;
    mpq_init(term);

    for (const krisite::mesh::Tri& tri : m.triangles) {
        for (int r = 0; r < 3; ++r) {
            const krisite::geom::HPointD& v = m.vertices[tri[r]];
            to_mpz(x[r][0], v.x);
            to_mpz(x[r][1], v.y);
            to_mpz(x[r][2], v.z);
            to_mpz(w[r], v.w);
        }
        // 分子: 3x3 行列式（余因子展開）
        mpz_set_ui(num, 0);
        for (int j = 0; j < 3; ++j) {
            const int j1 = (j + 1) % 3, j2 = (j + 2) % 3;
            mpz_mul(t0, x[1][j1], x[2][j2]);
            mpz_mul(t1, x[1][j2], x[2][j1]);
            mpz_sub(cof, t0, t1);
            mpz_mul(cof, cof, x[0][j]);
            mpz_add(num, num, cof);
        }
        mpz_mul(den, w[0], w[1]);
        mpz_mul(den, den, w[2]);
        KRI_CHECK_MSG(mpz_sgn(den) != 0, "体積: 構成点の w が 0");
        mpq_set_num(term, num);
        mpq_set_den(term, den);
        mpq_canonicalize(term);
        mpq_add(out, out, term);
    }

    mpq_clear(term);
    for (auto& row : x) {
        for (auto& e : row) mpz_clear(e);
    }
    for (auto& e : w) mpz_clear(e);
    mpz_clear(num);
    mpz_clear(den);
    mpz_clear(t0);
    mpz_clear(t1);
    mpz_clear(cof);
}

/// 入力メッシュの体積 x6（整数座標なので固定幅整数で出る）を mpq にする。
void input_volume6(mpq_ptr out, const TriMesh& m) {
    mpz_t z;
    mpz_init(z);
    to_mpz(z, krisite::mesh::signed_volume6(m));
    mpq_set_z(out, z);
    mpz_clear(z);
}

/// 軸平行な直方体の体積 x6 を mpz で（int64 では溢れる: 6*(2^21)^3 > 2^63）。
void box_volume6(mpz_ptr out, const std::int64_t ext[3]) {
    mpz_set_si(out, 6);
    for (int t = 0; t < 3; ++t) mpz_mul_si(out, out, static_cast<long>(ext[t]));
}

std::string to_str(mpq_srcptr q) {
    char* s = mpq_get_str(nullptr, 10, q);
    std::string r(s);
    void (*freefn)(void*, std::size_t);
    mp_get_memory_functions(nullptr, nullptr, &freefn);
    freefn(s, r.size() + 1);
    return r;
}

void check_case(const kritest::Case& c) {
    const TriMesh A = c.make_a(), B = c.make_b();

    mpq_t va, vb, vu, vi, vd, vu0, lhs, rhs;
    for (mpq_ptr p : {va, vb, vu, vi, vd, vu0, lhs, rhs}) mpq_init(p);

    input_volume6(va, A);
    input_volume6(vb, B);
    std::printf("\n  ケース %-3s %s\n", c.id, c.what);

    // §10.3: 解析的な期待体積との比較。**恒等式だけでは足りません。**
    // 恒等式は自己整合の検査なので、系統的な誤りが相関して入ると素通りし得ます。
    mpq_t want_i;
    mpq_init(want_i);
    bool has_analytic = false;
    if (c.box_pair) {
        const kritest::Aabb64 ba = kritest::aabb_of(A), bb = kritest::aabb_of(B);
        std::int64_t ea[3], eb[3], eo[3];
        for (int t = 0; t < 3; ++t) {
            ea[t] = ba.hi[t] - ba.lo[t];
            eb[t] = bb.hi[t] - bb.lo[t];
        }
        kritest::overlap_extent(ba, bb, eo);
        mpz_t z;
        mpz_init(z);
        box_volume6(z, ea);
        mpq_set_z(lhs, z);
        KRI_CHECK_MSG(mpq_equal(lhs, va) != 0,
                      std::string("ケース ") + c.id + ": |A| が解析値と食い違う");
        box_volume6(z, eb);
        mpq_set_z(lhs, z);
        KRI_CHECK_MSG(mpq_equal(lhs, vb) != 0,
                      std::string("ケース ") + c.id + ": |B| が解析値と食い違う");
        box_volume6(z, eo);
        mpq_set_z(want_i, z);
        mpz_clear(z);
        has_analytic = true;
        std::printf("    解析値: |A|=%s |B|=%s |A∩B|=%s（x6）\n", to_str(va).c_str(),
                    to_str(vb).c_str(), to_str(want_i).c_str());
    }

    for (unsigned d = 0; d <= 3; ++d) {
        BoolStats st;
        const BoolMesh u = boolean_op(A, B, BoolOp::Union, d, &st);
        const BoolMesh i = boolean_op(A, B, BoolOp::Intersection, d, &st);
        const BoolMesh df = boolean_op(A, B, BoolOp::Difference, d, &st);
        bool_mesh_volume6(vu, u);
        bool_mesh_volume6(vi, i);
        bool_mesh_volume6(vd, df);

        // |A ∪ B| + |A ∩ B| = |A| + |B|
        mpq_add(lhs, vu, vi);
        mpq_add(rhs, va, vb);
        const bool ok1 = mpq_equal(lhs, rhs) != 0;

        // |A \ B| = |A| - |A ∩ B|
        mpq_sub(rhs, va, vi);
        const bool ok2 = mpq_equal(vd, rhs) != 0;

        std::printf("    深度%u  ∪+∩=A+B %s   A\\B=A-∩ %s", d, ok1 ? "OK" : "**NG**",
                    ok2 ? "OK" : "**NG**");
        if (!ok1 || !ok2) {
            std::printf("\n      |A|=%s |B|=%s |∪|=%s |∩|=%s |\\|=%s", to_str(va).c_str(),
                        to_str(vb).c_str(), to_str(vu).c_str(), to_str(vi).c_str(),
                        to_str(vd).c_str());
        }
        std::printf("\n");

        KRI_CHECK_MSG(ok1, std::string("ケース ") + c.id + "（深度 " + std::to_string(d) +
                               "）: |A∪B| + |A∩B| != |A| + |B|");
        KRI_CHECK_MSG(ok2, std::string("ケース ") + c.id + "（深度 " + std::to_string(d) +
                               "）: |A\\B| != |A| - |A∩B|");

        // 解析値との比較（導出できるケースのみ）
        if (has_analytic) {
            KRI_CHECK_MSG(mpq_equal(vi, want_i) != 0,
                          std::string("ケース ") + c.id + "（深度 " + std::to_string(d) +
                              "）: |A∩B| が解析値と食い違う（得 " + to_str(vi) + " 期待 " +
                              to_str(want_i) + "）");
            mpq_add(rhs, va, vb);
            mpq_sub(rhs, rhs, want_i);
            KRI_CHECK_MSG(mpq_equal(vu, rhs) != 0, std::string("ケース ") + c.id + "（深度 " +
                                                       std::to_string(d) +
                                                       "）: |A∪B| が解析値と食い違う");
            mpq_sub(rhs, va, want_i);
            KRI_CHECK_MSG(mpq_equal(vd, rhs) != 0, std::string("ケース ") + c.id + "（深度 " +
                                                       std::to_string(d) +
                                                       "）: |A\\B| が解析値と食い違う");
        }

        // 深度不変性の体積版: 深度によらず同じ値であること（§10.2.1 の補強）
        if (d == 0) {
            mpq_set(vu0, vu);  // 深度 0 の ∪ を基準に取っておく（lhs は毎回上書きされる）
        } else {
            KRI_CHECK_MSG(mpq_equal(vu0, vu) != 0, std::string("ケース ") + c.id +
                                                       ": ∪ の体積が深度で変わった（深度 " +
                                                       std::to_string(d) + "）");
        }
    }

    mpq_clear(want_i);
    for (mpq_ptr p : {va, vb, vu, vi, vd, vu0, lhs, rhs}) mpq_clear(p);
}

}  // namespace

int main() {
    std::printf("\n  §10.3 体積の恒等式（GMP。位相では捕まらない誤りを突く）\n");
    for (const kritest::Case& c : kritest::corpus()) check_case(c);
    std::printf("\n");
    return kritest::finish("csg/volume_gmp");
}

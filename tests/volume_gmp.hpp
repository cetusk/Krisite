// Krisite — GMP による厳密な体積（テスト専用。LGPL の GMP に依存）
//
// **`SPEC-phase5.md` §3.0 の分担**: 浮動小数点は篩、**GMP は小標本と、
// 篩に引っかかったもの**。全件に掛けると意味がありません。
//
// **ライブラリ本体からは絶対に呼びません**（`CLAUDE.md` のライセンス方針）。
#ifndef KRISITE_TESTS_VOLUME_GMP_HPP
#define KRISITE_TESTS_VOLUME_GMP_HPP

#include <string>

#include <gmp.h>

#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/tri_mesh.hpp"

#include "gmp_oracle.hpp"
#include "test_util.hpp"

namespace kritest {

using kritest::oracle::to_mpz;

/// 閉じた向き付き三角形メッシュの符号付き体積 x6 = Σ det(v0, v1, v2)。
///
/// 頂点は同次座標 $(X : W)$ なので、三角形 1 枚の寄与は
/// $$ \frac{\det(X_0, X_1, X_2)}{W_0 W_1 W_2}. $$
/// 分子は整数行列式、分母は 3 つの $w$ の積です。**三角形ごとに分母が違う**ので、
/// `mpq` で足し合わせます（ここが GMP を要する理由）。
/// **`BoolMesh` と `SoupMesh` の両方**で使えるようにテンプレートにしています
/// （どちらも `vertices` / `triangles` を持つだけの型）。CP4 で**スープ経路の体積**も
/// 見るために要ります。
template <class Mesh>
inline void mesh_volume6(mpq_ptr out, const Mesh& m) {
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
inline void input_volume6(mpq_ptr out, const krisite::mesh::TriMesh& m) {
    mpz_t z;
    mpz_init(z);
    to_mpz(z, krisite::mesh::signed_volume6(m));
    mpq_set_z(out, z);
    mpz_clear(z);
}

/// 軸平行な直方体の体積 x6 を mpz で（int64 では溢れる: 6*(2^21)^3 > 2^63）。
inline void box_volume6(mpz_ptr out, const std::int64_t ext[3]) {
    mpz_set_si(out, 6);
    for (int t = 0; t < 3; ++t) mpz_mul_si(out, out, static_cast<long>(ext[t]));
}

inline std::string to_str(mpq_srcptr q) {
    char* s = mpq_get_str(nullptr, 10, q);
    std::string r(s);
    void (*freefunc)(void*, std::size_t) = nullptr;
    mp_get_memory_functions(nullptr, nullptr, &freefunc);
    freefunc(s, r.size() + 1);
    return r;
}

}  // namespace kritest

#endif  // KRISITE_TESTS_VOLUME_GMP_HPP

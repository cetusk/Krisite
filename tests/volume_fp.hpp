// Krisite — 浮動小数点の体積（**篩**。SPEC-phase5 §3.0）
//
// **これは検査ではなく篩です。** 目的は「GMP に回すものを絞る」ことで、
// **厳密性はここでは求めません。**
//
// 構成点は有理数なので、浮動小数点では丸めが乗ります。**閾値は実測から決めます**
// （`CLAUDE.md`「閾値は実測から決めてください。数字を発明しないこと」）。
// 測り方は `tests/csg/fp_volume_calib.cpp`、決めた値は `kIdentityTol`。
#ifndef KRISITE_TESTS_VOLUME_FP_HPP
#define KRISITE_TESTS_VOLUME_FP_HPP

#include <cmath>
#include <cstddef>

#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace kritest {

/// 多倍長整数を double に落とす。**上位リムから畳みます。**
///
/// 幅は最大でも $13b+27$（$b=26$ で 365 ビット）なので、
/// **double の指数範囲（約 1024 ビット）には収まります。** 落ちるのは仮数だけです。
template <std::size_t N>
inline double to_double(const krisite::arith::fixed_int<N>& v) {
    if (krisite::arith::sign(v) < 0) return -to_double(krisite::arith::neg(v));
    double r = 0.0;
    for (std::size_t i = N; i-- > 0;) {
        r = r * 18446744073709551616.0 + static_cast<double>(v[i]);
    }
    return r;
}

/// 出力メッシュの符号付き体積 × 6（浮動小数点）。
inline double volume6_fp(const krisite::csg::SoupMesh& m) {
    double acc = 0.0;
    for (const krisite::mesh::Tri& t : m.triangles) {
        double p[3][3];
        for (int k = 0; k < 3; ++k) {
            const krisite::geom::HPointD& v = m.vertices[t[static_cast<std::size_t>(k)]];
            const double w = to_double(v.w);
            p[k][0] = to_double(v.x) / w;
            p[k][1] = to_double(v.y) / w;
            p[k][2] = to_double(v.z) / w;
        }
        acc += p[0][0] * (p[1][1] * p[2][2] - p[1][2] * p[2][1]) -
               p[0][1] * (p[1][0] * p[2][2] - p[1][2] * p[2][0]) +
               p[0][2] * (p[1][0] * p[2][1] - p[1][1] * p[2][0]);
    }
    return acc;
}

/// 入力メッシュの符号付き体積 × 6。**こちらは厳密**（`signed_volume6`）を double に落とすだけ。
inline double volume6_fp(const krisite::mesh::TriMesh& m) {
    return to_double(krisite::mesh::signed_volume6(m));
}

/// 恒等式 $|A \cup B| + |A \cap B| = |A| + |B|$ の相対誤差。
inline double identity_error(double vu, double vi, double va, double vb) {
    const double scale = std::max(std::fabs(va), std::fabs(vb));
    if (scale == 0.0) return 0.0;
    return std::fabs((vu + vi) - (va + vb)) / scale;
}

/// 恒等式 $|A \setminus B| = |A| - |A \cap B|$ の相対誤差。
inline double difference_error(double vd, double va, double vi) {
    const double scale = std::max(std::fabs(va), 1.0);
    return std::fabs(vd - (va - vi)) / scale;
}

/// **篩の閾値。実測から導きました**（`CLAUDE.md`「数字を発明しないこと」）。
///
/// | 標本 | $P$ | 恒等式の相対誤差 |
/// |---|---|---|
/// | コーパス 112 標本（**GMP で正しさが確認済み**） | $\le 2{,}000$ | 最大 **6.8e-15** |
/// | 実データ 16 対 | 5万〜68万 | 最大 **1.6e-13** |
///
/// **$P$ は良い説明変数ではありません**（$R^2 = 0.32$）。
/// 正規化しても幅が縮まらないので（$\varepsilon/\sqrt{P}$ で 57 倍、$\varepsilon/P$ で 31 倍）、
/// **16 点からスケーリング則を作りません。** 平坦な閾値に、外挿ぶんの余裕を足します。
///
/// 外挿: $P = 10^7$（100 倍）でも $\varepsilon/\sqrt{P}$ 則で 6e-13、$\varepsilon/P$ 則で 4e-12。
/// **1e-9 はその 250 倍以上**の余裕があります。
///
/// **★ 感度の限界を理解して使ってください。** 格子 1 単位ぶんの四面体が
/// 欠けても、相対誤差は $10^{-18}$ 程度で**丸めの床より下**です。
/// **この篩は「大きな誤り」しか捕まえません。** 小さな誤りは
/// **小標本の GMP** が受け持ちます（§3.0 の分担）。
inline constexpr double kIdentityTol = 1e-9;

}  // namespace kritest

#endif  // KRISITE_TESTS_VOLUME_FP_HPP

// Krisite — 浮動小数点の体積恒等式の【校正】（SPEC-phase5 §3.0）
//
// **正しいと分かっている入力（コーパス）で相対誤差の分布を測ります。**
// **閾値を発明しないため**です（`CLAUDE.md`）。
//
// コーパスは GMP の体積検査を通っているので、**ここで出る誤差は丸めだけ**です。
#include <algorithm>
#include <cstdio>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"

#include "corpus.hpp"
#include "test_util.hpp"
#include "volume_fp.hpp"

using namespace krisite::csg;
using krisite::mesh::TriMesh;

int main() {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    std::printf("=== 浮動小数点の体積恒等式の校正（b=%d）===\n", KRISITE_COORD_BITS);
    std::printf("| ケース | 深度 | 恒等式の相対誤差 | 差の相対誤差 |\n|---|---:|---:|---:|\n");
    std::vector<double> errs, derrs;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        const double va = kritest::volume6_fp(a), vb = kritest::volume6_fp(b);
        for (unsigned d = 0; d <= 3; ++d) {
            BoolOptions o = kritest::corpus_options(d);
            o.local_bsp = true;
            o.split_contacts = true;
            ToMeshOptions tm;
            tm.split_contacts = true;
            double v[3];
            int k = 0;
            for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
                v[k++] =
                    kritest::volume6_fp(to_mesh(boolean(from_mesh(a), from_mesh(b), op, o), tm));
            }
            const double e = kritest::identity_error(v[0], v[1], va, vb);
            const double de = kritest::difference_error(v[2], va, v[1]);
            errs.push_back(e);
            derrs.push_back(de);
            std::printf("| %-4s | %u | %.3e | %.3e |\n", c.id, d, e, de);
        }
    }
    std::sort(errs.begin(), errs.end());
    std::sort(derrs.begin(), derrs.end());
    const auto q = [](const std::vector<double>& v, double p) {
        return v[static_cast<std::size_t>(p * static_cast<double>(v.size() - 1))];
    };
    std::printf(
        "\n## 分布（%zu 標本）\n\n| 量 | 中央 | 90% | 99% | **最大** "
        "|\n|---|---:|---:|---:|---:|\n",
        errs.size());
    std::printf("| 恒等式 | %.3e | %.3e | %.3e | **%.3e** |\n", q(errs, 0.5), q(errs, 0.9),
                q(errs, 0.99), errs.back());
    std::printf("| 差 | %.3e | %.3e | %.3e | **%.3e** |\n", q(derrs, 0.5), q(derrs, 0.9),
                q(derrs, 0.99), derrs.back());
    // **回帰の番人。** 丸めの床が上がったら、篩の前提が崩れます。
    // **コーパスの実測（6.83e-15）の 100 倍**を上限に置きます。
    KRI_CHECK_MSG(errs.back() < 1e-12, "コーパスでの恒等式の相対誤差が想定を超えた");
    KRI_CHECK_MSG(derrs.back() < 1e-12, "コーパスでの差の相対誤差が想定を超えた");
    // **篩の閾値は、コーパスの床より十分上にあること**（さもなくば全件が引っかかる）
    KRI_CHECK_MSG(kritest::kIdentityTol > errs.back() * 100.0,
                  "篩の閾値がコーパスの丸めの床に近すぎる");
    return kritest::finish("csg/fp_volume_calib");
}

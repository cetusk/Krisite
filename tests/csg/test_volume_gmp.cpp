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
// （別経路どうしの突き合わせ。`docs/ROADMAP.md`「正解器は被検体と別経路で書く」）。
//
// GMP は LGPL。テスト専用（`KRISITE_BUILD_TESTS_WITH_GMP=ON` のときだけビルドされる）。
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"

#include "corpus.hpp"
#include "gmp_oracle.hpp"
#include "test_util.hpp"
#include "volume_gmp.hpp"

using namespace krisite::csg;
using krisite::mesh::TriMesh;
using kritest::box_volume6;
using kritest::input_volume6;
using kritest::mesh_volume6;
using kritest::oracle::to_mpz;

namespace {

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

    // 深度 0〜3 は固定深度、**4 番目は適応分割 + early-out**（SPEC-phase2 §3）。
    //
    // **体積は「断片の選択の誤り」を捕まえます**（SPEC-phase1 §10.2 の但し書き）。
    // early-out はセルの隅 1 点で分類を決めるので、**誤ればここに出ます。**
    // 位相検査は「割れているか」しか見ないので、この 2 つは役割が違います。
    for (unsigned d = 0; d <= 4; ++d) {
        BoolStats st;
        BoolOptions opt;
        opt.depth = (d <= 3) ? d : 3;
        opt.adaptive = (d == 4);
        opt.early_out = (d == 4);
        // **4 番目は最適化を全部入れた構成**（適応分割 + early-out + 構成点の保持）。
        // 出荷時の構成なので、ここが体積で守られていることに意味があります。
        opt.cache_points = (d == 4);
        // **既定に依存させないこと**（§9.4 の CI ジョブで既定が反転する）
        opt.split_contacts = true;
        const BoolMesh u = boolean_op(A, B, BoolOp::Union, opt, &st);
        const BoolMesh i = boolean_op(A, B, BoolOp::Intersection, opt, &st);
        const BoolMesh df = boolean_op(A, B, BoolOp::Difference, opt, &st);
        mesh_volume6(vu, u);
        mesh_volume6(vi, i);
        mesh_volume6(vd, df);

        // |A ∪ B| + |A ∩ B| = |A| + |B|
        mpq_add(lhs, vu, vi);
        mpq_add(rhs, va, vb);
        const bool ok1 = mpq_equal(lhs, rhs) != 0;

        // |A \ B| = |A| - |A ∩ B|
        mpq_sub(rhs, va, vi);
        const bool ok2 = mpq_equal(vd, rhs) != 0;

        std::printf("    %-8s ∪+∩=A+B %s   A\\B=A-∩ %s",
                    (d <= 3) ? ("深度" + std::to_string(d)).c_str() : "適応+全部",
                    ok1 ? "OK" : "**NG**", ok2 ? "OK" : "**NG**");
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

        // ---- §9.4.2: **分裂の有無で体積が一致すること** ----
        //
        // 「分裂は体積を変えない」は §5.1.3 の**定理**です。変異 8 / 9 の素通りで
        // 間接的に固定するより、**不変量そのものを直接検査する**ほうが確かです。
        // 分裂の誤りは体積では捕まりませんが、**分裂が幾何を動かしていないこと**は
        // ここでしか言えません。
        {
            BoolOptions ns = opt;
            ns.split_contacts = false;
            const BoolMesh u2 = boolean_op(A, B, BoolOp::Union, ns, &st);
            const BoolMesh i2 = boolean_op(A, B, BoolOp::Intersection, ns, &st);
            const BoolMesh d2 = boolean_op(A, B, BoolOp::Difference, ns, &st);
            mpq_t nu, ni, nd;
            mpq_init(nu);
            mpq_init(ni);
            mpq_init(nd);
            mesh_volume6(nu, u2);
            mesh_volume6(ni, i2);
            mesh_volume6(nd, d2);
            const bool same =
                mpq_equal(nu, vu) != 0 && mpq_equal(ni, vi) != 0 && mpq_equal(nd, vd) != 0;
            KRI_CHECK_MSG(same, std::string("ケース ") + c.id + "（深度 " + std::to_string(d) +
                                    "）: **分裂の有無で体積が変わった**（§9.4.2 / §5.1.3）");
            for (mpq_ptr q : {nu, ni, nd}) mpq_clear(q);
        }

        // ---- CP4: **局所 BSP ↔ 過剰分割の体積一致**（`SPEC-phase3.md` §10.1）----
        //
        // 位相検査は「割れているか」しか見ません。**断片の選択の誤り**は体積でしか
        // 出ないので、分割方式を替えたときはここを通す必要があります。
        //
        // **(b) はスープ経路そのものの検査です。** $(C, \chi)$ は面の向きを反転しても
        // 変わらないので、ここを見ていなかったために**出力の向きの誤りが CP2 から
        // CP3 まで残りました**（`IMPL-phase3.md` §6.1、変異 15）。
        //
        // **両側を明示的に指定します。** 既定に依存させると、既定が変わったときに
        // 自分自身との比較になります（`CLAUDE.md`）。
        {
            mpq_t sv[3][3];
            for (auto& row : sv) {
                for (auto& q : row) mpq_init(q);
            }
            // **3 通り**: 0 = 過剰分割 / 1 = 局所 BSP / 2 = 局所 BSP + NSI。
            //
            // **NSI は三角形分割を変えるので、バイトでは守れません**
            // （`IMPL-phase5.md` §20）。**体積が残る唯一の不変量**です。
            // 過剰分割は NSI をまったく使わないので、**最も独立な対照**になります。
            for (int bsp = 0; bsp < 3; ++bsp) {
                BoolOptions so = opt;
                so.local_bsp = (bsp >= 1);
                so.split_contacts = true;
                krisite::csg::ToMeshOptions tmo;
                tmo.split_contacts = true;
                int k = 0;
                for (BoolOp bop : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
                    krisite::csg::PolySoup sa = krisite::csg::from_mesh(A);
                    krisite::csg::PolySoup sb = krisite::csg::from_mesh(B);
                    if (bsp == 2) {
                        // **コーパスの入力は自己交差しません**（§9.0 の前提）。
                        sa.nsi.assign(sa.sources.size(), 1);
                        sb.nsi.assign(sb.sources.size(), 1);
                    }
                    const krisite::csg::SoupMesh sm =
                        krisite::csg::to_mesh(krisite::csg::boolean(sa, sb, bop, so), tmo);
                    mesh_volume6(sv[bsp][k++], sm);
                }
            }
            // (a) 分割方式によらず同じ体積
            for (int k = 0; k < 3; ++k) {
                KRI_CHECK_MSG(mpq_equal(sv[0][k], sv[1][k]) != 0,
                              std::string("ケース ") + c.id + "（深度 " + std::to_string(d) +
                                  "）: **局所 BSP と過剰分割で体積が違う**（§5.4 / §10.1）");
                KRI_CHECK_MSG(mpq_equal(sv[1][k], sv[2][k]) != 0,
                              std::string("ケース ") + c.id + "（深度 " + std::to_string(d) +
                                  "）: **NSI の宣言で体積が変わった**"
                                  "（`SPEC-phase3.md` §5.6。切ってよい前提が破れている）");
            }
            // (b) 二項正解器と同じ体積か
            {
                mpq_ptr ref[3] = {vu, vi, vd};
                const char* nm[3] = {"∪", "∩", "\\"};
                for (int k = 0; k < 3; ++k) {
                    KRI_CHECK_MSG(mpq_equal(sv[1][k], ref[k]) != 0,
                                  std::string("ケース ") + c.id + "（深度 " + std::to_string(d) +
                                      "）: スープの " + nm[k] + " の体積が二項正解器と違う（得 " +
                                      to_str(sv[1][k]) + " 期待 " + to_str(ref[k]) + "）");
                }
            }
            for (auto& row : sv) {
                for (auto& q : row) mpq_clear(q);
            }
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

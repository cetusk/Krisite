// Krisite — 辺平面の構成と、断片の相対内部の点（SPEC-phase3.md §3.1 / §2.1）
//
// **Phase 3 の段 0 と CP1 の検査です。**
//
// 段 0 は代表点の構成を平面ベースに戻します。Phase 1 / 2 の「頂点 → 対角線の中点 →
// 3 頂点の重心」は 2 段目以降が点を組み合わせる操作で、被符号値が 21b+46
// （b=21 で 487 ビット / 8 リム）に達していました。
//
// ここで見るのは 3 つです。
//
//   1. 辺平面が**辺の直線を含み、支持平面と一致しない**こと（§3.1.2）
//   2. 軸の選択が**正準**であること（同じ入力で同じ平面。§3.1.3）
//   3. 代表点が**相対内部**にあること、そして**両方の経路が実際に使われる**こと
//
// **3 の後半が空回りの番人です。** 予備経路（角のオフセット）が一度も走らないなら、
// その経路は検証されていません。
#include <cstdio>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/interior.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;
using namespace krisite::geom;
using krisite::mesh::TriMesh;

namespace {

/// 軸 `j` が棄却されるべきか（N = 0 か、支持平面の法線と平行）。
///
/// **実装と同じ判定を、テスト側で独立に書きます。** 実装の条件式をそのまま呼ぶと
/// 「同じ変形を使う正解器」になります（`CLAUDE.md`）。
bool rejected_reason_holds(const IPoint& p1, const IPoint& p2, const PlaneD& sp, int j) {
    const std::int64_t d[3] = {static_cast<std::int64_t>(p2.x) - p1.x,
                               static_cast<std::int64_t>(p2.y) - p1.y,
                               static_cast<std::int64_t>(p2.z) - p1.z};
    const std::int64_t cand[3][3] = {{0, d[2], -d[1]}, {-d[2], 0, d[0]}, {d[1], -d[0], 0}};
    const std::int64_t* n = cand[j];
    if (n[0] == 0 && n[1] == 0 && n[2] == 0) return true;  // 条件 1 で棄却
    // 条件 2: N x N_s == 0 か。**支持平面の法線を整数に落とせる範囲で比較します。**
    // 三角形の法線は 2b+3 ビットなので、b <= 21 なら int64 に収まります（45 ビット）。
    // 積は最大 45 + 22 = 67 ビットで int64 を超えるため、__int128 相当の広い型で見ます。
    using W = krisite::arith::fixed_int<limbs::kOffset>;
    const W nx = krisite::arith::from_i64<limbs::kOffset>(n[0]);
    const W ny = krisite::arith::from_i64<limbs::kOffset>(n[1]);
    const W nz = krisite::arith::from_i64<limbs::kOffset>(n[2]);
    const W sx = krisite::arith::resize<limbs::kOffset>(sp.a);
    const W sy = krisite::arith::resize<limbs::kOffset>(sp.b);
    const W sz = krisite::arith::resize<limbs::kOffset>(sp.c);
    // **正解器なので `arith::cross` は使いません。** 被検体と同じ関数を通すと、
    // 同じ間違いをして両方が通ります（`CLAUDE.md`「正解器は被検体と別経路で書く」）。
    // ここでは外積の成分を `det2` で直接書き下します。
    using krisite::arith::row2;
    const auto cx = krisite::arith::det2(row2<limbs::kOffset>{ny, nz},
                                         row2<limbs::kOffset>{sy, sz});  // ny*sz - nz*sy
    const auto cy = krisite::arith::det2(row2<limbs::kOffset>{nz, nx},
                                         row2<limbs::kOffset>{sz, sx});  // nz*sx - nx*sz
    const auto cz = krisite::arith::det2(row2<limbs::kOffset>{nx, ny},
                                         row2<limbs::kOffset>{sx, sy});  // nx*sy - ny*sx
    return krisite::arith::is_zero(cx) && krisite::arith::is_zero(cy) &&
           krisite::arith::is_zero(cz);
}

/// 三角形 1 枚から支持平面と 3 枚の辺平面を作り、断片に組み立てる。
Fragment make_triangle_fragment(PlaneTable& t, const IPoint& a, const IPoint& b, const IPoint& c) {
    const PlaneD sp = plane_from_triangle(a, b, c);
    Fragment f;
    f.support = t.intern(sp).id;
    f.owner = 0;
    // 頂点 i = support ∩ edge[i-1] ∩ edge[i] なので、辺の順序を合わせる
    f.edge = {t.intern(plane_from_edge(c, a, sp)).id, t.intern(plane_from_edge(a, b, sp)).id,
              t.intern(plane_from_edge(b, c, sp)).id};
    return f;
}

void test_edge_plane_properties() {
    struct Case {
        const char* what;
        IPoint a, b, c;
    };
    const Case cases[] = {
        {"軸平行（z=0 の三角形）", {0, 0, 0}, {100, 0, 0}, {0, 100, 0}},
        {"斜面", {1, 2, 3}, {70, -50, 110}, {-40, 90, 20}},
        {"1 成分が 0 の辺（Δ_z = 0）", {0, 0, 5}, {60, 40, 5}, {10, -30, 90}},
        {"軸に平行な辺（Δ が e_x と平行）", {-20, 7, 7}, {50, 7, 7}, {3, 60, -11}},
        {"細長い三角形", {0, 0, 0}, {1000, 1, 0}, {500, 2, 3}},
    };
    for (const Case& cs : cases) {
        const PlaneD sp = plane_from_triangle(cs.a, cs.b, cs.c);
        const IPoint* v[3] = {&cs.a, &cs.b, &cs.c};
        for (int i = 0; i < 3; ++i) {
            const IPoint& p1 = *v[i];
            const IPoint& p2 = *v[(i + 1) % 3];
            const IPoint& p3 = *v[(i + 2) % 3];
            Axis k{};
            const PlaneD e = plane_from_edge(p1, p2, sp, &k);
            const std::string tag = std::string(cs.what) + " 辺 " + std::to_string(i);
            // 1. 辺の両端を含む
            KRI_CHECK_MSG(side(e, p1) == 0, tag + ": 端点 p1 が辺平面上にない");
            KRI_CHECK_MSG(side(e, p2) == 0, tag + ": 端点 p2 が辺平面上にない");
            // 2. 支持平面と一致しない（第 3 点が載っていないことで確認）
            KRI_CHECK_MSG(side(e, p3) != 0, tag + ": 第 3 点が辺平面上にある（支持平面と一致）");
            // 3. 正準（同じ入力で同じ軸・同じ係数）
            Axis k2{};
            const PlaneD e2 = plane_from_edge(p1, p2, sp, &k2);
            KRI_CHECK_MSG(k == k2 && plane_cmp(e, e2) == 0, tag + ": 軸の選択が非正準");
            // 4. **条件を満たす最小の k** であること（正準性の定義そのもの）。
            //    小さい k が棄却された理由は「N = 0」か「N が N_s と平行」のどちらか
            //    でなければなりません。**片方でも成り立たないなら選択は最小ではない。**
            for (int j = 0; j < static_cast<int>(k); ++j) {
                KRI_CHECK_MSG(rejected_reason_holds(p1, p2, sp, j),
                              tag + ": より小さい軸 " + std::to_string(j) +
                                  " が使えるのに選ばれていない（正準でない）");
            }
            // 5. 辺平面のビット幅が支持平面を超えない（§3.1.2 の主張）
            KRI_CHECK_MSG(krisite::arith::min_bits(e.a) <= bits::kEdgeNormal &&
                              krisite::arith::min_bits(e.b) <= bits::kEdgeNormal &&
                              krisite::arith::min_bits(e.c) <= bits::kEdgeNormal,
                          tag + ": 辺平面の法線が b+1 を超えた");
            KRI_CHECK_MSG(krisite::arith::min_bits(e.d) <= bits::kEdgeOffset,
                          tag + ": 辺平面のオフセットが 2b+2 を超えた");
        }
    }
}

void test_interior_point_is_inside() {
    struct Case {
        const char* what;
        IPoint a, b, c;
    };
    const Case cases[] = {
        {"軸平行", {0, 0, 0}, {100, 0, 0}, {0, 100, 0}},
        {"斜面", {1, 2, 3}, {70, -50, 110}, {-40, 90, 20}},
        {"細長い（幅 1）", {0, 0, 0}, {1000, 1, 0}, {500, 2, 0}},
        {"小さい（辺 2）", {0, 0, 0}, {2, 0, 0}, {0, 2, 0}},
        {"原点から離れた大きな座標",
         {100000, -90000, 80000},
         {100050, -90000, 80000},
         {100000, -89950, 80000}},
    };
    for (const Case& cs : cases) {
        PlaneTable t;
        const Fragment f = make_triangle_fragment(t, cs.a, cs.b, cs.c);
        InteriorStats st;
        const HPointD x = interior_point(t, f, nullptr, &st);
        const std::string tag = cs.what;
        // 支持平面の上にある
        KRI_CHECK_MSG(side(t.at(f.support), x) == 0, tag + ": 代表点が支持平面上にない");
        // すべての辺平面に対して**厳密に内側**（符号が非零で、第 3 頂点と同じ側）
        const std::size_t n = vertex_count(f);
        for (std::size_t k = 0; k < n; ++k) {
            const int s = side(t.at(f.edge[k]), x);
            KRI_CHECK_MSG(s != 0, tag + ": 代表点が辺平面上にある（相対内部でない）");
        }
        KRI_CHECK_MSG(st.axis_line + st.corner_offset == 1, tag + ": どちらの経路も決まらない");
    }
}

/// **空回りの番人。** コーパス全体で両方の経路が使われることを固定する。
void test_both_paths_are_exercised() {
    std::size_t axis = 0, corner = 0, failed = 0, configs = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            for (unsigned d = 0; d <= 3; ++d) {
                BoolStats st;
                boolean_op(a, b, op, kritest::phase1_options(d), &st);
                axis += st.interior.axis_line;
                corner += st.interior.corner_offset;
                failed += st.interior.axis_failed;
                ++configs;
            }
        }
    }
    std::printf("\n  代表点の経路（%zu 構成、SPEC-phase3 §2.1）\n", configs);
    std::printf("    主経路（軸平行直線） %zu / 予備経路（角のオフセット） %zu\n", axis, corner);
    std::printf("    主経路が外れた回数 %zu\n", failed);

    KRI_CHECK_MSG(axis > 0, "主経路が一度も決まっていない。**空回りです**");
    KRI_CHECK_MSG(corner > 0,
                  "**予備経路が一度も走っていない。** 角のオフセットが検証されていません");
}

/// **★ float ヒント `approx` が符号と桁を保つこと**（`DESIGN-phase5-hotspots.md` §20.3）。
///
/// > **初版は「最上位リムを符号付き、それ以外を符号なしとして積む」形でした。**
/// > **【すべての負の値】で $2^{128}$ を返します。**
/// > $-1$（$N$ リム）は `{~0, ..., ~0}` で、最上位は $-1$ ですが次のリムは $2^{64}-1$。
/// > **double の仮数は 53 ビットなので $2^{64}-1$ は $2^{64}$ に丸められ、
/// > $-1 \cdot 2^{64} + 2^{64} = 0$ になって符号が消えます。**
///
/// **厳密性は無傷でした**（ヒントは候補を出すだけで、判定は厳密演算）。
/// **だから正しさの検査では 3 フェーズにわたって見つかりませんでした。**
/// **主経路の成功率が 2.5〜25.9% に落ちていただけです**（修正後 70.5〜86.2%）。
///
/// > **`CLAUDE.md`「機構を足したら、それが空回りしていないことを別に検査してください」。**
/// > **`test_both_paths_are_exercised` は「両方の経路が使われる」ことしか見ておらず、
/// > 【どちらが主経路か】を見ていませんでした。**
void test_approx_sign_and_magnitude() {
    std::size_t checks = 0;
    const long long vals[] = {0, 1, -1, 2, -2, 1000, -1000, 1048575, -1048575, 33554431, -33554431};
    for (long long v : vals) {
        // **リム数を変えても同じ**であること（キャンセルはリム数に依存しました）
        const double a2 = krisite::csg::detail::approx(krisite::arith::from_i64<2>(v));
        const double a4 = krisite::csg::detail::approx(krisite::arith::from_i64<4>(v));
        const double a6 = krisite::csg::detail::approx(krisite::arith::from_i64<6>(v));
        const auto want = static_cast<double>(v);
        KRI_CHECK_MSG(a2 == want, "approx(2 リム) が真の値と違う" + kritest::pair_msg(want, a2));
        KRI_CHECK_MSG(a4 == want, "approx(4 リム) が真の値と違う" + kritest::pair_msg(want, a4));
        KRI_CHECK_MSG(a6 == want, "approx(6 リム) が真の値と違う" + kritest::pair_msg(want, a6));
        checks += 3;
    }
    // **仮数に収まらない大きさ**では、符号と相対誤差で見ます。
    {
        const auto big = krisite::arith::mul(krisite::arith::from_i64<4>(-1234567890123456789LL),
                                             krisite::arith::from_i64<4>(1000000007LL));
        const double got = krisite::csg::detail::approx(big);
        KRI_CHECK_MSG(got < 0, "approx: 大きい負の値で符号が消えた");
        const double want = -1234567890123456789.0 * 1000000007.0;
        const double rel = (got - want) / want;
        KRI_CHECK_MSG(rel > -1e-12 && rel < 1e-12,
                      "approx: 大きい値の相対誤差が大きすぎる " + std::to_string(rel));
        checks += 2;
    }
    std::printf("  approx の符号と桁: %zu 件\n", checks);
}

}  // namespace

int main() {
    std::printf("\n  辺平面と代表点 — SPEC-phase3 §3.1 / §2.1（CP1）\n");
    test_edge_plane_properties();
    test_interior_point_is_inside();
    test_both_paths_are_exercised();
    test_approx_sign_and_magnitude();
    std::printf("\n");
    return kritest::finish("csg/interior");
}

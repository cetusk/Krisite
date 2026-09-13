// Krisite — **辺の内部の点を平面で構成する**述語の GMP 差分テスト
//
// `DESIGN-phase5-vertex-level.md` §7（案 H）。ビット幅の導出は
// `geom/widths.hpp` の `kSumNormal` / `kSumOffset` / `kHomoW` / `kHomoXyz`。
//
// **被検体**: `edge_interior_point(P1,P2,P3,P4,v,w)`。
// $Q = P_3 \pm P_4$ を作り、$P_1 \cap P_2 \cap Q$ を返す。
//
// **正解器は別経路です** — **有理数（mpq）で実座標に直して**、
// 点が線分 $vw$ の内部にあるかを $m = v + t(w-v)$、$0 < t < 1$ で確かめます。
// 被検体は同次整数のまま符号だけで判定するので、**経路が違います。**
//
// **退化した配置を明示的に構成します**（`CLAUDE.md`「乱択は共平面・共線・重複頂点を
// ほぼ生成しません」「対称性の高い入力はバグを構造的に隠します」）。
//
// **ビット幅の実測もここで行います。** 実測が理論上界を超えたら設計の誤りです。
//
// GMP は LGPL。テスト専用（`KRISITE_BUILD_TESTS_WITH_GMP=ON` でのみビルド）。
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "krisite/geom/predicates.hpp"

#include "gmp_oracle.hpp"
#include "test_util.hpp"

using namespace krisite;
using namespace krisite::geom;
using kritest::Rng;
using kritest::oracle::signed_bits;
using kritest::oracle::to_mpz;

namespace {

/// 3 点から平面を作る（`plane_from_triangle` と同じ向きの規約）。
PlaneD pl(long long x1, long long y1, long long z1, long long x2, long long y2, long long z2,
          long long x3, long long y3, long long z3) {
    const IPoint a{static_cast<std::int32_t>(x1), static_cast<std::int32_t>(y1),
                   static_cast<std::int32_t>(z1)};
    const IPoint b{static_cast<std::int32_t>(x2), static_cast<std::int32_t>(y2),
                   static_cast<std::int32_t>(z2)};
    const IPoint c{static_cast<std::int32_t>(x3), static_cast<std::int32_t>(y3),
                   static_cast<std::int32_t>(z3)};
    return plane_from_triangle(a, b, c);
}

/// 同次点を mpq の 3 成分へ（**除算する別経路**）。
struct Q3 {
    mpq_t x, y, z;
    Q3() {
        mpq_init(x);
        mpq_init(y);
        mpq_init(z);
    }
    ~Q3() {
        mpq_clear(x);
        mpq_clear(y);
        mpq_clear(z);
    }
    Q3(const Q3&) = delete;
    Q3& operator=(const Q3&) = delete;
};

void to_q3(Q3* out, const HPointD& h) {
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    to_mpz(d, h.w);
    to_mpz(n, h.x);
    mpq_set_num(out->x, n);
    mpq_set_den(out->x, d);
    mpq_canonicalize(out->x);
    to_mpz(n, h.y);
    mpq_set_num(out->y, n);
    mpq_set_den(out->y, d);
    mpq_canonicalize(out->y);
    to_mpz(n, h.z);
    mpq_set_num(out->z, n);
    mpq_set_den(out->z, d);
    mpq_canonicalize(out->z);
    mpz_clear(n);
    mpz_clear(d);
}

std::size_t max_bits_w = 0, max_bits_xyz = 0, max_bits_qn = 0, max_bits_qd = 0;

/// 固定幅整数の**有効なビット幅**（符号ビットを含む）を mpz 経由で測る。
template <std::size_t N>
std::size_t bits_of(const krisite::arith::fixed_int<N>& x) {
    mpz_t z;
    mpz_init(z);
    to_mpz(z, x);
    const std::size_t r = signed_bits(z);
    mpz_clear(z);
    return r;
}
std::size_t checked = 0, failures = 0;

void fail(const char* what, const char* tag) {
    ++failures;
    std::printf("  **FAIL** %s（%s）\n", what, tag);
}

/// 1 件を検査する。**正解器は mpq で線分の媒介変数 t を求める。**
void check(const char* tag, const PlaneD& p1, const PlaneD& p2, const PlaneD& p3, const PlaneD& p4,
           const HPointD& v, const HPointD& w) {
    // 前提が成り立っているかを先に確かめる（**空回りしていないこと**）
    if (side(p1, v) != 0 || side(p1, w) != 0 || side(p2, v) != 0 || side(p2, w) != 0) {
        fail("前提: p1 / p2 が両端を通っていない", tag);
        return;
    }
    if (side(p3, v) != 0 || side(p3, w) == 0) {
        fail("前提: p3 が v を通り w を通らない、になっていない", tag);
        return;
    }
    if (side(p4, w) != 0 || side(p4, v) == 0) {
        fail("前提: p4 が w を通り v を通らない、になっていない", tag);
        return;
    }
    const int s4v = side(p4, v), s3w = side(p3, w);
    const PlaneSum q = plane_sum(p3, p4, -s4v * s3w);
    const HPointD m = edge_interior_point(p1, p2, p3, p4, v, w);
    ++checked;

    max_bits_w = std::max(max_bits_w, bits_of(m.w));
    max_bits_xyz = std::max(max_bits_xyz, bits_of(m.x));
    max_bits_xyz = std::max(max_bits_xyz, bits_of(m.y));
    max_bits_xyz = std::max(max_bits_xyz, bits_of(m.z));
    max_bits_qn = std::max(max_bits_qn, bits_of(q.a));
    max_bits_qn = std::max(max_bits_qn, bits_of(q.b));
    max_bits_qn = std::max(max_bits_qn, bits_of(q.c));
    max_bits_qd = std::max(max_bits_qd, bits_of(q.d));

    // ---- 正解器: mpq で t を求める ----
    Q3 qv, qw, qm;
    to_q3(&qv, v);
    to_q3(&qw, w);
    to_q3(&qm, m);
    mpq_t num, den, t, tk, tmp;
    mpq_init(num);
    mpq_init(den);
    mpq_init(t);
    mpq_init(tk);
    mpq_init(tmp);
    bool have_t = false, on_line = true;
    mpq_ptr vv[3] = {qv.x, qv.y, qv.z};
    mpq_ptr ww[3] = {qw.x, qw.y, qw.z};
    mpq_ptr mm[3] = {qm.x, qm.y, qm.z};
    for (int k = 0; k < 3; ++k) {
        mpq_sub(den, ww[k], vv[k]);
        mpq_sub(num, mm[k], vv[k]);
        if (mpq_sgn(den) == 0) {
            if (mpq_sgn(num) != 0) on_line = false;
            continue;
        }
        mpq_div(tk, num, den);
        if (!have_t) {
            mpq_set(t, tk);
            have_t = true;
        } else if (mpq_cmp(t, tk) != 0) {
            on_line = false;
        }
    }
    mpq_set_ui(tmp, 1, 1);
    const bool inside = have_t && on_line && mpq_sgn(t) > 0 && mpq_cmp(t, tmp) < 0;
    if (!inside) fail("構成した点が辺の内部にない", tag);
    mpq_clear(num);
    mpq_clear(den);
    mpq_clear(t);
    mpq_clear(tk);
    mpq_clear(tmp);
}

/// 4 平面から $v = P_1\cap P_2\cap P_3$、$w = P_1\cap P_2\cap P_4$ を作って検査する。
void run(const char* tag, const PlaneD& p1, const PlaneD& p2, const PlaneD& p3, const PlaneD& p4) {
    const HPointD v = intersect3(p1, p2, p3);
    const HPointD w = intersect3(p1, p2, p4);
    check(tag, p1, p2, p3, p4, v, w);
}

}  // namespace

int main(int argc, char** argv) {
    const long n_rand = (argc > 1) ? std::atol(argv[1]) : 200000;
    std::printf("## 辺の内部の点（平面による構成）の GMP 差分テスト\n\n");
    std::printf("| 設定 | 値 |\n|---|---|\n");
    std::printf("| b（座標ビット） | %d |\n", KRISITE_COORD_BITS);
    std::printf("| 乱択の件数 | %ld |\n\n", n_rand);

    const long long M = (1LL << (KRISITE_COORD_BITS - 1)) - 1;

    // ---- 1. 退化した配置を【明示的に】構成する ----
    //
    // **乱択は共平面・共線・軸平行・極値をほぼ生成しません**（`CLAUDE.md`）。
    // **そして対称な配置だけでは退化の一部しか作れない**ので、斜めも入れます。

    // 1a. 軸平行（最も対称。x 軸に沿った辺）
    run("軸平行", pl(0, 0, 0, 1000, 0, 0, 0, 0, 1000),  // z=0 ではなく y=0 平面
        pl(0, 0, 0, 1000, 0, 0, 0, 1000, 0),            // z=0 平面
        pl(0, 0, 0, 0, 1000, 0, 0, 0, 1000),            // x=0 平面
        pl(500, 0, 0, 500, 1000, 0, 500, 0, 1000));     // x=500 平面

    // 1b. **極値**（座標の上限付近。幅が最大になる）
    run("極値", pl(-M, -M, -M, M, -M, -M, -M, -M, M), pl(-M, -M, -M, M, -M, -M, M, M, -M),
        pl(-M, -M, -M, -M, M, -M, -M, -M, M), pl(M, -M, -M, M, M, -M, M, -M, M));

    // 1c. **斜め**（軸に整列しない。非対称）
    run("斜め", pl(1, 2, 3, 101, 47, 5, 7, 113, 61), pl(1, 2, 3, 101, 47, 5, 89, 13, 211),
        pl(1, 2, 3, 3, 5, 7, 11, 13, 17), pl(101, 47, 5, 103, 51, 11, 111, 59, 23));

    // 1d. **辺が非常に短い**（両端が近い。実データの残存辺がこの形だった）
    run("短い辺", pl(0, 0, 0, 1000, 0, 0, 0, 0, 1000), pl(0, 0, 0, 1000, 0, 0, 0, 1000, 0),
        pl(0, 0, 0, 0, 1000, 0, 0, 0, 1000), pl(1, 0, 0, 1, 1000, 0, 1, 0, 1000));

    // 1e. **P3 と P4 が平行**（和が 0 ベクトルになりかける。符号の選択が効く）
    run("P3 と P4 が平行", pl(0, 0, 0, 1000, 0, 0, 0, 0, 1000), pl(0, 0, 0, 1000, 0, 0, 0, 1000, 0),
        pl(0, 0, 0, 0, 1000, 0, 0, 0, 1000),
        pl(700, 0, 0, 700, 0, 1000, 700, 1000, 0));  // 向きが逆の x=700

    // 1f. **P3 と P4 が大きく傾いている**（和の法線が長くなる）
    run("大きく傾いた P3 / P4", pl(0, 0, 0, M, 0, 0, 0, 0, M), pl(0, 0, 0, M, 0, 0, 0, M, 0),
        pl(0, 0, 0, 1, M, 0, 1, 0, M), pl(M, 0, 0, M - 1, M, 0, M - 1, 0, M));

    const std::size_t explicit_cases = checked;
    std::printf("**明示的に構成した退化配置 %zu 件**\n\n", explicit_cases);

    // ---- 2. 乱択（対称でない配置を広く踏む）----
    Rng rng(0x51D3ull);
    const auto rc = [&]() { return static_cast<long long>(rng.range(-M, M)); };
    long skipped = 0;
    for (long i = 0; i < n_rand; ++i) {
        const PlaneD p1 = pl(rc(), rc(), rc(), rc(), rc(), rc(), rc(), rc(), rc());
        const PlaneD p2 = pl(rc(), rc(), rc(), rc(), rc(), rc(), rc(), rc(), rc());
        const PlaneD p3 = pl(rc(), rc(), rc(), rc(), rc(), rc(), rc(), rc(), rc());
        const PlaneD p4 = pl(rc(), rc(), rc(), rc(), rc(), rc(), rc(), rc(), rc());
        if (is_degenerate(p1) || is_degenerate(p2) || is_degenerate(p3) || is_degenerate(p4)) {
            ++skipped;
            continue;
        }
        // 3 平面が一点で交わらない組は飛ばす（`intersect3` の前提）
        const auto d12 = radial_dir(p1, p2);
        if (arith::is_zero(d12.x) && arith::is_zero(d12.y) && arith::is_zero(d12.z)) {
            ++skipped;
            continue;
        }
        const HPointD v = intersect3(p1, p2, p3);
        const HPointD w = intersect3(p1, p2, p4);
        if (arith::is_zero(v.w) || arith::is_zero(w.w)) {
            ++skipped;
            continue;
        }
        if (side(p3, w) == 0 || side(p4, v) == 0) {
            ++skipped;  // 退化（p3 が w も通る等）。前提が成り立たない
            continue;
        }
        check("乱択", p1, p2, p3, p4, v, w);
    }

    std::printf("| 項目 | 値 |\n|---|---:|\n");
    std::printf("| 検査した件数 | %zu |\n", checked);
    std::printf("| 前提が成り立たず飛ばした件数 | %ld |\n", skipped);
    std::printf("| **不一致** | **%zu** |\n\n", failures);

    const long long b = KRISITE_COORD_BITS;
    struct Row {
        const char* name;
        std::size_t got;
        long long bound;
    } rows[] = {
        {"Q の法線成分", max_bits_qn, 2 * b + 4},
        {"Q のオフセット", max_bits_qd, 3 * b + 6},
        {"構成点の w", max_bits_w, 6 * b + 13},
        {"構成点の x,y,z", max_bits_xyz, 7 * b + 15},
    };
    std::printf("| 量 | 実測の最大 | 理論上界 | 判定 |\n|---|---:|---:|---|\n");
    for (const Row& r : rows) {
        const bool over = static_cast<long long>(r.got) > r.bound;
        if (over) ++failures;
        std::printf("| %s | %zu | %lld | %s |\n", r.name, r.got, r.bound,
                    over ? "**上界を超えた**" : "上界の中");
    }

    // ---- 3. **対照: 符号の選択を外すと、点が辺の外に出ること** ----
    //
    // **「不一致 0」が「検査が空回りしている」と区別できなければ使えません**
    // （`CLAUDE.md`「機構が空回りしていないことを別に検査してください」）。
    // **符号を常に +1 にすると、辺の外へ出る件が必ずあるはず**です。
    std::size_t control_outside = 0, control_total = 0;
    {
        Rng r2(0x9A17ull);
        const auto rc2 = [&]() { return static_cast<long long>(r2.range(-M, M)); };
        for (long i = 0; i < 2000; ++i) {
            const PlaneD p1 = pl(rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2());
            const PlaneD p2 = pl(rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2());
            const PlaneD p3 = pl(rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2());
            const PlaneD p4 = pl(rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2(), rc2());
            if (is_degenerate(p1) || is_degenerate(p2) || is_degenerate(p3) || is_degenerate(p4)) {
                continue;
            }
            const auto d12 = radial_dir(p1, p2);
            if (arith::is_zero(d12.x) && arith::is_zero(d12.y) && arith::is_zero(d12.z)) continue;
            const HPointD v = intersect3(p1, p2, p3);
            const HPointD w = intersect3(p1, p2, p4);
            if (side(p3, w) == 0 || side(p4, v) == 0) continue;
            ++control_total;
            // **符号を選ばない**（常に +1）
            const PlaneSum q = plane_sum(p3, p4, 1);
            if (side(q, v) * side(q, w) >= 0) ++control_outside;
        }
    }
    std::printf("\n| 対照（符号を選ばない） | 値 |\n|---|---:|\n");
    std::printf("| 試した件数 | %zu |\n", control_total);
    std::printf("| **辺を分けられなかった件数** | **%zu** |\n", control_outside);
    if (control_outside == 0) {
        std::printf("\n**対照が 1 件も外れません。検査が空回りしています**\n");
        ++failures;
    }

    // **空回りの検査**: 退化配置が 1 件も通っていなければ、検査していないのと同じ
    if (explicit_cases < 6) {
        std::printf("\n**明示的な退化配置が %zu 件しか通っていません**\n", explicit_cases);
        ++failures;
    }
    std::printf("\n**不一致 %zu 件**\n", failures);
    return failures == 0 ? 0 : 1;
}

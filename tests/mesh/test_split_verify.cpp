// Krisite — 分裂の検証（増分計算）が、従来経路と同じ答えを返すこと
//
// `SPEC-phase2.md` §5.5 の検算と、§5.1.2.1 の事後の多様体性は、もともと
// `mesh::check_topology` を 2 回呼んで求めていました。**実測で出口の大半を占めます**
// （`HANDOVER.md` §4.3）。**同値な増分計算に置き換えました**（`IMPL-v2.md` §2）。
//
// ---
//
// ## このテストが要る理由 ★
//
// **`test_cp5n.cpp` が、コーパス 420 構成で両経路の一致を検査しています。**
// **しかし 420 構成すべてで「両方とも真」でした**（事後の非多様体 0 件）。
//
// > **一致は、偽を返す側を一度も踏んでいなければ何も言っていません。**
// > `CLAUDE.md`「番人が数えるべきは『機構が動いたか』ではなく、
// > **『その変異が観測可能になる条件が揃ったか』**です」。
//
// **そこで機会を作る入力をここに置きます。** コーパスではなく**テスト専用**です
// （`SPEC-phase4.md` §7.3。**退化を狙う入力とは性格が違います** — こちらは
// 「検証が偽を返す条件」を作るためのもので、幾何としては壊れた入力を含みます）。
//
// `split_contacts` は**三角形の索引だけ**を受け取るので、手で三角形を並べれば
// 幾何を作らずに条件を揃えられます。
//
// ## 押さえる分岐
//
//   (a) 辺の次数が 2 でない        `edge_manifold = false`
//   (b) 辺の向きが揃っていない      `vertex_manifold = false`（リンクの出次数が 2 になる）
//   (c) 退化三角形がある            前提が破れ、**正解器に落ちる**（`verify_fallback`）
//   (d) 多様体な入力                両方とも真（真を返す側の対照）
//
// **(a)(b)(c) では `check_topology` と答えが一致することを直接比べます。**
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/mesh/split.hpp"
#include "krisite/mesh/topology.hpp"

#include "test_util.hpp"

using krisite::mesh::check_topology;
using krisite::mesh::SplitOptions;
using krisite::mesh::SplitStats;
using krisite::mesh::TopologyReport;
using krisite::mesh::Tri;

namespace {

std::size_t g_false_manifold = 0;  ///< 事後の検査が偽になった件数（**空回りの番人**）
std::size_t g_fallback = 0;        ///< 前提の破れで正解器に落ちた件数

/// 増分計算と従来経路で `split_contacts` を回し、答えが一致することを確かめる。
///
/// **比べるのは、事後の検査が使う 2 つの量だけ**です（`unresolved_post` と
/// $\Delta V$ / $\Delta E$）。**それがこの置き換えの範囲そのもの**です。
void compare(const std::string& tag, const std::vector<Tri>& tris, std::size_t vertex_count) {
    SplitOptions fast;
    fast.verify_manifold = true;
    fast.verify_delta = true;
    fast.verify_naive = false;
    SplitOptions naive = fast;
    naive.verify_naive = true;

    SplitStats sf{}, sn{};
    std::vector<std::uint32_t> of, on;
    const std::vector<Tri> rf =
        krisite::mesh::split_contacts(tris, vertex_count, &of, &sf, nullptr, nullptr, nullptr, fast);
    const std::vector<Tri> rn = krisite::mesh::split_contacts(tris, vertex_count, &on, &sn, nullptr,
                                                             nullptr, nullptr, naive);

    // **出力は経路に依りません**（検証は出力に触りません）
    KRI_CHECK_MSG(rf == rn, tag + ": 検証の経路で出力が変わった");
    KRI_CHECK_MSG(of == on, tag + ": 検証の経路で複製された頂点が変わった");

    KRI_CHECK_MSG(sf.unresolved_post == sn.unresolved_post,
                  tag + ": 事後の非多様体の判定が違う" +
                      kritest::pair_msg(sf.unresolved_post, sn.unresolved_post));
    KRI_CHECK_MSG(sf.unresolved == sn.unresolved,
                  tag + ": unresolved が違う" + kritest::pair_msg(sf.unresolved, sn.unresolved));
    KRI_CHECK_MSG(sf.actual_delta_v == sn.actual_delta_v,
                  tag + ": ΔV が違う" + kritest::pair_msg(sf.actual_delta_v, sn.actual_delta_v));
    KRI_CHECK_MSG(sf.actual_delta_e == sn.actual_delta_e,
                  tag + ": ΔE が違う" + kritest::pair_msg(sf.actual_delta_e, sn.actual_delta_e));
    KRI_CHECK_MSG(sf.actual_delta_chi == sn.actual_delta_chi,
                  tag + ": Δχ が違う" +
                      kritest::pair_msg(sf.actual_delta_chi, sn.actual_delta_chi));

    // **従来経路の答えを、`check_topology` から直接も確かめます。**
    // 正解器の正解器です（`sn` は同じ関数の中で `check_topology` を呼んでいるので、
    // ここが独立した突き合わせになります）
    const TopologyReport before = check_topology(tris);
    const TopologyReport after = check_topology(rn);
    const bool ok_manifold = after.edge_manifold && after.vertex_manifold;
    const std::size_t expect_post = (!tris.empty() && !ok_manifold) ? 1u : 0u;
    if (sf.verify_fallback == 0) {
        KRI_CHECK_MSG(sf.unresolved_post == expect_post,
                      tag + ": 増分計算が check_topology と食い違う" +
                          kritest::pair_msg(sf.unresolved_post, expect_post));
        KRI_CHECK_MSG(sf.actual_delta_v == after.v - before.v,
                      tag + ": ΔV が check_topology と食い違う" +
                          kritest::pair_msg(sf.actual_delta_v, after.v - before.v));
        KRI_CHECK_MSG(sf.actual_delta_e == after.e - before.e,
                      tag + ": ΔE が check_topology と食い違う" +
                          kritest::pair_msg(sf.actual_delta_e, after.e - before.e));
    }

    g_false_manifold += (sf.unresolved_post != 0) ? 1u : 0u;
    g_fallback += (sf.verify_fallback != 0) ? 1u : 0u;
    std::printf("    %-34s 事後の非多様体 %zu / 退避 %zu / ΔV %zu / ΔE %zu\n", tag.c_str(),
                sf.unresolved_post, sf.verify_fallback, sf.actual_delta_v, sf.actual_delta_e);
}

/// 四面体（頂点 `base`..`base+3`）。**外向きに向き付けます。**
void tetra(std::vector<Tri>& out, std::uint32_t base) {
    const std::uint32_t a = base, b = base + 1, c = base + 2, d = base + 3;
    out.push_back({a, c, b});
    out.push_back({a, b, d});
    out.push_back({b, c, d});
    out.push_back({c, a, d});
}

}  // namespace

int main() {
    std::printf("\n  分裂の検証（増分計算） — 従来経路との一致（IMPL-v2 §2）\n");

    // ---- (d) 対照: 多様体な入力。両方とも真 ---------------------------------
    {
        std::vector<Tri> t;
        tetra(t, 0);
        compare("(d) 四面体 1 個（多様体）", t, 4);
    }

    // ---- (a) 辺の次数が 4 のまま残る -----------------------------------------
    //
    // **辺 (0,1) に 4 枚が接し、その 4 枚が「別の場所」で環に繋がっています。**
    // 4 枚が同じ連結成分に入るので組が作れず（§5.1.2.1）、`RadialGeom` を渡して
    // いないので radial sort も働かず、**分裂させずに次数 4 のまま残ります**（§5.1.2.2）。
    //
    // **`IMPL-phase5.md` §98 の残件と同じ側の分岐です**（出力の辺の次数が 2 でない）。
    {
        // **幾何としては成り立ちません。`split_contacts` は索引しか見ません。**
        std::vector<Tri> t = {
            {0, 1, 2}, {1, 0, 3}, {0, 1, 4}, {1, 0, 5},
            // 4 枚を 1 つの連結成分に繋ぐ「別の場所」（次数 2 の辺で 4 枚を環にする）
            {2, 3, 6}, {3, 4, 6}, {4, 5, 6}, {5, 2, 6},
        };
        compare("(a) 次数 4 の辺が残る", t, 7);
    }

    // ---- (b) 辺の向きが揃っていない ------------------------------------------
    //
    // **1 枚だけ裏返します。** 辺の次数は 2 のままですが、$\#(u,v) \ne \#(v,u)$ に
    // なるので、その頂点でリンクの出次数が 2 になり `vertex_manifold` が偽になります。
    {
        std::vector<Tri> t;
        tetra(t, 0);
        std::swap(t[0][1], t[0][2]);  // 1 枚だけ向きを反転
        compare("(b) 向きが揃っていない", t, 4);
    }

    // ---- (c) 退化三角形。前提が破れて正解器に落ちる ---------------------------
    {
        std::vector<Tri> t;
        tetra(t, 0);
        t.push_back({0, 0, 1});  // 同じ頂点を 2 度使う三角形
        compare("(c) 退化三角形（前提が破れる）", t, 4);
    }

    // ---- 空回りの番人 --------------------------------------------------------
    //
    // **偽を返す側を踏んでいなければ、一致は何も言っていません**（冒頭の注記）。
    KRI_CHECK_MSG(g_false_manifold >= 2,
                  "**検証が偽を返す側を踏んでいません。** (a)(b) が条件を作れていません");
    KRI_CHECK_MSG(g_fallback == 1,
                  "**前提の破れ（退化三角形）で正解器に落ちる経路を踏んでいません**");
    std::printf("    偽を返した %zu 件 / 正解器へ退避 %zu 件\n", g_false_manifold, g_fallback);

    return kritest::finish("mesh/split_verify");
}

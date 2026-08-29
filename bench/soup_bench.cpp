// Krisite — スープ経路の計数と時間（`SPEC-phase3.md` §11）
//
// **性能ではなく設計判断のための基準線です。** 完了条件ではありません。
//
// §11 が要求する分担で測ります。
//
//     from_mesh   入口。凸分割の片数・平面数・時間
//     boolean     中核。時間・断片数・early-out 発火率
//     to_mesh     出口。時間・併合された頂点数・T 頂点数
//
// > **EMBER の 1.6 ms はスープを出すまでの時間です。** 縫合・T 解決・三角形化まで
// > 含めた時間と比べると、**同じものを測っていません**（§11）。だから分けて出します。
//
// **単発の差に原因を当てないこと。** 実行ごとに ±15% 程度動きます（`CLAUDE.md`）。
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"

#include "corpus.hpp"

using namespace krisite;
using namespace krisite::csg;

namespace {

using Clock = std::chrono::steady_clock;

double ms_since(Clock::time_point t0) {
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(Clock::now() - t0)
        .count();
}

BoolOptions all_on(unsigned depth, bool adaptive) {
    BoolOptions o;
    o.depth = depth;
    o.cull_planes = true;
    o.adaptive = adaptive;
    o.leaf_threshold = 0;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    return o;
}

/// 構成点の最大ビット幅（`test_soup.cpp` と同じ数え方）。
std::size_t soup_bits(const PolySoup& s) {
    std::size_t w = 0;
    for (const Poly& q : s.polys) {
        for (std::size_t i = 0; i < vertex_count(q.frag); ++i) {
            const geom::HPointD v = fragment_vertex(s.table, q.frag, i);
            for (std::size_t b : {arith::min_bits(v.x), arith::min_bits(v.y), arith::min_bits(v.z),
                                  arith::min_bits(v.w)}) {
                w = std::max(w, b);
            }
        }
    }
    return w;
}

/// 入口・中核・出口を分けて計測する（§11 の ★）。
void stage_timing() {
    std::printf("\n## 入口・中核・出口の内訳（全 %zu ケース × 3 演算、適応分割 + 全機構）\n\n",
                kritest::corpus().size());
    double t_from = 0, t_bool = 0, t_mesh = 0;
    std::size_t polys_in = 0, planes_in = 0, frags = 0, tris = 0;
    std::size_t merged = 0, t_inserted = 0, raycasts = 0, regions = 0, eo_cells = 0, cells = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const mesh::TriMesh a = c.make_a(), b = c.make_b();
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            auto t0 = Clock::now();
            const PolySoup sa = from_mesh(a), sb = from_mesh(b);
            t_from += ms_since(t0);
            polys_in += sa.polys.size() + sb.polys.size();
            planes_in += sa.table.size() + sb.table.size();

            BoolStats st{};
            t0 = Clock::now();
            const PolySoup r = boolean(sa, sb, op, all_on(3, true), &st);
            t_bool += ms_since(t0);
            frags += st.fragments;
            raycasts += st.raycasts;
            regions += st.regions;
            eo_cells += st.early_out_cells;
            cells += st.total_cells;

            ToMeshOptions tm;
            tm.split_contacts = true;
            ToMeshStats ts{};
            t0 = Clock::now();
            const SoupMesh m = to_mesh(r, tm, &ts);
            t_mesh += ms_since(t0);
            tris += m.triangles.size();
            merged += ts.merged_by_value;
            t_inserted += ts.t.inserted;
        }
    }
    const double total = t_from + t_bool + t_mesh;
    std::printf("| 区分 | 時間 (ms) | 割合 | 主な計数 |\n");
    std::printf("|---|---:|---:|---|\n");
    std::printf("| `from_mesh` | %.1f | %.1f%% | 多角形 %zu、平面 %zu |\n", t_from,
                100.0 * t_from / total, polys_in, planes_in);
    std::printf(
        "| **`boolean`（中核）** | **%.1f** | **%.1f%%** | 断片 %zu、領域 %zu、"
        "レイキャスト %zu |\n",
        t_bool, 100.0 * t_bool / total, frags, regions, raycasts);
    std::printf("| `to_mesh` | %.1f | %.1f%% | 三角形 %zu、値で併合 %zu、T 頂点 %zu |\n", t_mesh,
                100.0 * t_mesh / total, tris, merged, t_inserted);
    std::printf("| 合計 | %.1f | 100.0%% | セル %zu（うち early-out %zu = %.1f%%）|\n", total,
                cells, eo_cells,
                100.0 * static_cast<double>(eo_cells) / static_cast<double>(cells));
}

/// 連鎖の段数に対する断片数・領域数・分類コストの推移（§11）。
///
/// **CP2 の分類コストは sources 数に比例するはずでした。** 大域レイキャストのまま
/// なので、ここが item 2（セグメントトレース、Phase 5 へ降格）の必要性の判断材料です。
void chain_scaling() {
    std::printf("\n## 連鎖の段数に対する推移（適応分割 + 全機構）\n\n");
    // **`st.regions` と `st.raycasts` は同じ地点で加算されるので、比は常に 1.0 です。**
    // 「一致しようがない指標」は検査になりません。分類コストの段数依存は
    // **断片あたりのレイキャスト**で見ます。
    std::printf("| 段 | source 数 | 断片 | 分類キー数 | 断片あたり | 最大ビット幅 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|\n");
    const BoolOptions o = all_on(3, true);
    for (int stage = 1; stage <= 3; ++stage) {
        std::size_t frags = 0, regions = 0, raycasts = 0, srcs = 0, bits = 0;
        for (const kritest::Case& c : kritest::corpus()) {
            const mesh::TriMesh a = c.make_a(), b = c.make_b();
            const mesh::TriMesh d = kritest::corpus()[1].make_b();
            const PolySoup sa = from_mesh(a), sb = from_mesh(b), sd = from_mesh(d);
            BoolStats st{};
            PolySoup r = boolean(sa, sb, BoolOp::Union, o, &st);
            if (stage >= 2) r = boolean(r, sd, BoolOp::Difference, o, &st);
            if (stage >= 3) r = boolean(r, sa, BoolOp::Union, o, &st);
            frags += st.fragments;
            regions += st.regions;
            raycasts += st.raycasts;
            srcs = r.source_count();
            bits = std::max(bits, soup_bits(r));
        }
        std::printf("| %d | %zu | %zu | %zu | %.3f | **%zu** |\n", stage, srcs, frags, regions,
                    frags == 0 ? 0.0 : static_cast<double>(raycasts) / static_cast<double>(frags),
                    bits);
    }
}

/// 局所 BSP と過剰分割の断片数の比（§5.4 の効果）。
void bsp_ratio() {
    std::printf("\n## 局所 BSP と過剰分割（適応分割 + 全機構）\n\n");
    std::size_t raw_o = 0, raw_b = 0, out_o = 0, out_b = 0, used = 0, skipped = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const mesh::TriMesh a = c.make_a(), b = c.make_b();
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            BoolOptions oo = all_on(3, true);
            oo.local_bsp = false;
            BoolOptions ob = all_on(3, true);
            ob.local_bsp = true;
            BoolStats so{}, sb{};
            boolean(from_mesh(a), from_mesh(b), op, oo, &so);
            boolean(from_mesh(a), from_mesh(b), op, ob, &sb);
            raw_o += so.raw_fragments;
            raw_b += sb.raw_fragments;
            out_o += so.fragments;
            out_b += sb.fragments;
            used += sb.bsp_cuts_used;
            skipped += sb.bsp_cuts_skipped;
        }
    }
    std::printf("| | 過剰分割 | 局所 BSP | 比 |\n|---|---:|---:|---:|\n");
    std::printf("| 生の断片 | %zu | **%zu** | **%.1f%%** |\n", raw_o, raw_b,
                100.0 * static_cast<double>(raw_b) / static_cast<double>(raw_o));
    std::printf("| 正準化後 | %zu | **%zu** | **%.1f%%** |\n", out_o, out_b,
                100.0 * static_cast<double>(out_b) / static_cast<double>(out_o));
    std::printf("\n切断候補 %zu（切った %zu / 省いた %zu = %.1f%% を省略）\n", used + skipped, used,
                skipped,
                100.0 * static_cast<double>(skipped) / static_cast<double>(used + skipped));
}

}  // namespace

int main() {
    std::printf("# Krisite — スープ経路の計数と時間（SPEC-phase3 §11）\n");
    std::printf("\n**単発の差に原因を当てないこと。実行ごとに ±15%% 程度動きます。**\n");
    stage_timing();
    chain_scaling();
    bsp_ratio();
    std::printf("\n");
    return 0;
}

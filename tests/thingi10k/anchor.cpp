// Krisite — 清浄な絶対時間の錨（`PERF.md` #2′ / CP1.9）
//
// **絶対時間は、機械が空いていることを確認してからでないと使えません**
// （`IMPL-phase5.md` §26）。この道具は**錨を採るためだけ**のものです。
//
// **錨は 3 点採ります**（小・基準・大）。1 点では絶対スケールしか与えませんが、
// **3 点あれば比の模型そのものを検定できます。**
//
//   T2 / T1  対  ΣW_b^(2) / ΣW_b^(1)   が合えば模型が使える
//
// **構成を必ず併記すること。** 錨は 1 組しかないので、条件が曖昧だと後で使えません。
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"
#include "krisite/par/thread_pool.hpp"

#include "thingi10k/loader.hpp"

using namespace krisite;

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    const int reps = (argc > 1) ? std::atoi(argv[1]) : 5;
    const unsigned nthreads = (std::getenv("KRI_THREADS") != nullptr)
                                  ? static_cast<unsigned>(std::atoi(std::getenv("KRI_THREADS")))
                                  : 16u;

    // ---- 構成の宣言（**錨に必須**）----
    std::printf("# 錨（清浄な絶対時間）\n\n");
    std::printf("| 項目 | 値 |\n|---|---|\n");
    std::printf("| b（座標ビット幅） | **%d** |\n", KRISITE_COORD_BITS);
    std::printf("| スレッド | **%u** |\n", nthreads);
    std::printf("| **算術の検査**（`KRISITE_CHECKED_ARITH`） | **%s** |\n",
                (KRISITE_CHECKED_ARITH != 0) ? "**ON**" : "**OFF**");
    std::printf("| **§5.5 の検算**（`check_topology` × 2） | **常に ON**（旗が無い） |\n");
    std::printf("| NSI の宣言 | **OFF**（既定） |\n");
    std::printf("| レイキャストの索引 | ON（既定） |\n");
    std::printf("| 適応分割 / 深度 | ON / 6 |\n");
    std::printf("| 演算 | ∪ ∩ ＼ の 3 つ + `to_mesh`（`check_one` と同じ） |\n");
    std::printf("| 反復 | **%d 回。最小値を採用**（競合下では最小値） |\n", reps);
    std::printf(
        "| **汚染の判定** | **最大/最小 > %.2f なら棄却**（内側の load average は"
        "外部の負荷を見ないので、**ばらつきで判定します**。`IMPL-phase5.md` §26） |\n",
        1.10);
    std::printf("| 最適化 | `-O2` |\n\n");

    std::printf(
        "| 対 | 入力 n | 秒（最小） | 秒（中央） | 秒（最大） | P | 領域 | レイ |"
        " 三角形検査 | ΣP_葉 | 葉あたり最大 | arrange | classify | stitch | 判定 | 最大/最小 |\n");
    std::printf(
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|\n");

    for (int a = 2; a + 1 < argc; a += 2) {
        const auto qa =
            krithingi::quantize(krithingi::load_kmesh(argv[a]), krithingi::make_transform(1000));
        const auto qb = krithingi::quantize(krithingi::load_kmesh(argv[a + 1]),
                                            krithingi::make_transform(1001));
        const std::size_t n = qa.mesh.triangles.size() + qb.mesh.triangles.size();
        if (qa.mesh.triangles.empty() || qb.mesh.triangles.empty()) {
            std::printf("| %s | **読み込み失敗** |\n", argv[a]);
            continue;
        }
        const csg::PolySoup A = csg::from_mesh(qa.mesh), B = csg::from_mesh(qb.mesh);
        par::ThreadPool pool(nthreads);

        std::vector<double> ts;
        csg::BoolStats best{};
        double best_ms_arr = 0, best_ms_cls = 0, best_ms_sti = 0;
        std::size_t P = 0;
        for (int r = 0; r < reps; ++r) {
            csg::BoolOptions o;
            o.depth = 6;
            o.adaptive = true;
            o.cull_planes = true;
            o.early_out = true;
            o.cache_points = true;
            o.local_bsp = true;
            o.split_contacts = true;
            o.threads = nthreads;
            o.pool = &pool;
            csg::ToMeshOptions tm;
            tm.split_contacts = true;
            tm.threads = nthreads;
            tm.pool = &pool;

            csg::BoolStats acc{};
            std::size_t polys = 0;
            const auto t0 = std::chrono::steady_clock::now();
            for (csg::BoolOp op :
                 {csg::BoolOp::Union, csg::BoolOp::Intersection, csg::BoolOp::Difference}) {
                csg::BoolStats bs;
                const csg::PolySoup soup = csg::boolean(A, B, op, o, &bs);
                csg::to_mesh(soup, tm);
                polys += soup.polys.size();
                acc.regions += bs.regions;
                acc.raycasts += bs.raycasts;
                acc.ray_tri_tests += bs.ray_tri_tests;
                acc.leaf_input_total += bs.leaf_input_total;
                acc.leaf_input_max = std::max(acc.leaf_input_max, bs.leaf_input_max);
                acc.ms_arrange += bs.ms_arrange;
                acc.ms_classify += bs.ms_classify;
                acc.ms_stitch += bs.ms_stitch;
            }
            const double dt =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            ts.push_back(dt);
            if (r == 0 || dt < *std::min_element(ts.begin(), ts.end() - 1)) {
                best = acc;
                best_ms_arr = acc.ms_arrange;
                best_ms_cls = acc.ms_classify;
                best_ms_sti = acc.ms_stitch;
                P = polys;
            }
        }
        std::sort(ts.begin(), ts.end());
        std::printf(
            "| `%s` | %zu | **%.3f** | %.3f | %.3f | %zu | %zu | %zu | %zu | %zu | %zu |"
            " %.0f | %.0f | %.0f | %s | %.3f |\n",
            argv[a] + 22, n, ts.front(), ts[ts.size() / 2], ts.back(), P, best.regions,
            best.raycasts, best.ray_tri_tests, best.leaf_input_total, best.leaf_input_max,
            best_ms_arr, best_ms_cls, best_ms_sti,
            (ts.back() / ts.front() <= 1.10) ? "**清浄**" : "**棄却（汚染）**",
            ts.back() / ts.front());
    }
    return 0;
}

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
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"
#include "krisite/par/thread_pool.hpp"

#include "thingi10k/loader.hpp"

using namespace krisite;

namespace {

/// **`cp1.txt` の添字から変換を引きます**（`IMPL-phase5.md` §33 / §40）。
///
/// **`cp1.cpp` は `prepare(raw, 1000 + i)` を使います**（$i$ は一覧での添字）。
/// **固定値を使うと別の入力を測ることになります。** §33 で `edge_source.cpp` を
/// この形にしましたが、**このツールの監査を怠って同じ誤りが残っていました。**
///
/// **手で変換を書けない形にするのが対処**です。呼び出し側は `"AxB"` を渡すだけ。
std::map<std::string, std::size_t> kri_cp1_index(const char* list) {
    std::map<std::string, std::size_t> index;
    std::ifstream f(list);
    std::string id, nf;
    std::size_t i = 0;
    while (f >> id >> nf) index[id] = i++;
    return index;
}

}  // namespace

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
    // ---- ★ 定常状態で測る（`PERF.md` §1.5 / §5.4）------------------------------
    //
    // **「冷えた 1 回目」は最速値であって代表値ではありません。**
    // CP2 / CP3 は何時間も連続で回るので、**定常（熱い）が代表的**です。
    // そして**冷やすのは待つしかない**のに対し、**温めるのは走らせればよい。**
    //
    // **ドリフト検査を測定に埋め込みます**（アイドルより確実。`/proc/loadavg` が
    // 使えなかったときと同じ形）。**先頭の対を最後にもう一度測り、差で判定します。**
    const double warm_sec =
        (std::getenv("KRI_WARMUP_SEC") != nullptr) ? std::atof(std::getenv("KRI_WARMUP_SEC")) : 0.0;
    // **判定は 3%。** 判定したい差は 1.70 対 1.87 = 10% で、ドリフト δ は
    // 係数の幅を最大 $(1+\delta)$ 倍に膨らませます。**信号の 1/3 以下**が要ります。
    // **発明した数字ではありません。**
    const double kDriftTol = 0.03;
    std::printf("| **ウォームアップ** | **%.0f 秒**（先頭の対を空回し。定常に入れるため） |\n",
                warm_sec);
    std::printf("| **ドリフト検査** | **先頭の対を最後にもう一度測る。差 > %.0f%% なら棄却** |\n",
                kDriftTol * 100);
    std::printf(
        "| **測定順** | **大きさ順にしない**（後の対ほど熱く、残留ドリフトが規模と"
        "相関して $W_b$ の指数に直接効きます） |\n");
    std::printf("| 最適化 | `-O2` |\n\n");

    std::printf(
        "| 対 | 入力 n | 秒（最小） | P | **Σ多角形²** | 葉 | **BSP 枠** | **BSP 使用** |"
        " **辺/多角形** | **平面/多角形** | arrange |"
        " **収集 %%** | **存在 %%** | **準備 %%** | **断片 %%** | **共平面 %%** |"
        " classify | stitch | 判定 |\n");
    std::printf(
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
        "---:|---:|---:|---|\n");

    const auto index = kri_cp1_index("data/thingi10k/cp1.txt");
    // ドリフト検査用。**先頭の対の最小値**を覚えておき、最後の測定と比べます
    std::string drift_key;
    double drift_first = 0.0, drift_last = 0.0;
    bool warmed = false;

    // **対は `"AxB"` で渡します。** ファイル名を渡す形だと変換を手で書くことになり、
    // §33 / §40 の取り違えが再発します
    for (int a = 2; a < argc; ++a) {
        const std::string key = argv[a];
        const std::size_t xp = key.find('x');
        if (xp == std::string::npos || index.count(key.substr(0, xp)) == 0 ||
            index.count(key.substr(xp + 1)) == 0) {
            std::printf("| `%s` | **一覧に無い** |\n", key.c_str());
            continue;
        }
        const std::string ida = key.substr(0, xp), idb = key.substr(xp + 1);
        const auto qa =
            krithingi::quantize(krithingi::load_kmesh("data/thingi10k/kmesh/" + ida + ".kmesh"),
                                krithingi::make_transform(1000 + index.at(ida)));
        const auto qb =
            krithingi::quantize(krithingi::load_kmesh("data/thingi10k/kmesh/" + idb + ".kmesh"),
                                krithingi::make_transform(1000 + index.at(idb)));
        const std::size_t n = qa.mesh.triangles.size() + qb.mesh.triangles.size();
        if (qa.mesh.triangles.empty() || qb.mesh.triangles.empty()) {
            std::printf("| `%s` | **読み込み失敗** |\n", key.c_str());
            continue;
        }
        const csg::PolySoup A = csg::from_mesh(qa.mesh), B = csg::from_mesh(qb.mesh);
        par::ThreadPool pool(nthreads);

        // ---- ウォームアップ（**先頭の対だけ。結果は捨てます**）------------------
        if (!warmed) {
            warmed = true;
            const auto w0 = std::chrono::steady_clock::now();
            std::size_t nrun = 0;
            while (std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count() <
                   warm_sec) {
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
                for (csg::BoolOp op :
                     {csg::BoolOp::Union, csg::BoolOp::Intersection, csg::BoolOp::Difference}) {
                    csg::to_mesh(csg::boolean(A, B, op, o), tm);
                }
                ++nrun;
            }
            if (warm_sec > 0) {
                std::printf("<!-- ウォームアップ %.0f 秒 / %zu 巡（`%s`）-->\n", warm_sec, nrun,
                            key.c_str());
                std::fflush(stdout);
            }
        }

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
                acc.leaf_input_sq += bs.leaf_input_sq;
                acc.leaf_poly_sq += bs.leaf_poly_sq;
                // **局所 BSP の切断数。** 合成が実データを代表しているかの切り分けに要る
                // （`PERF.md` §1.1 の完了条件。**時間ではなく仕事の量で比べます**）
                acc.leaf_nonempty += bs.leaf_nonempty;
                acc.leaf_single_src += bs.leaf_single_src;
                acc.bsp_cut_slots += bs.bsp_cut_slots;
                acc.bsp_cuts_used += bs.bsp_cuts_used;
                // **無次元群**（`PERF.md` §1.8）。合成と「単位の中身」を揃えるため
                acc.fragments += bs.fragments;
                acc.frag_edges_total += bs.frag_edges_total;
                acc.frag_edges_count += bs.frag_edges_count;
                acc.frag_edges_max = std::max(acc.frag_edges_max, bs.frag_edges_max);
                acc.leaf_planes_total += bs.leaf_planes_total;
                acc.ms_arr_gather += bs.ms_arr_gather;
                acc.ms_arr_present += bs.ms_arr_present;
                acc.ms_arr_prep += bs.ms_arr_prep;
                acc.ms_arr_frag += bs.ms_arr_frag;
                acc.ms_arr_coplanar += bs.ms_arr_coplanar;
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
        // **先頭の対の最小値**を基準に、**最後の測定**と比べます
        if (drift_key.empty()) {
            drift_key = key;
            drift_first = ts.front();
        }
        if (key == drift_key) drift_last = ts.front();
        // **段別の割合**（`PERF.md` §1.1）
        const double atot = best.ms_arr_gather + best.ms_arr_present + best.ms_arr_prep +
                            best.ms_arr_frag + best.ms_arr_coplanar;
        const double aq = atot > 0 ? 100.0 / atot : 0.0;
        std::printf(
            "| `%s` | %zu | **%.3f** | %zu | **%zu** | %zu | %zu | %zu | **%.2f** | **%.2f** |"
            " %.0f | **%.1f** | **%.1f** | **%.1f** | **%.1f** | **%.1f** | %.0f | %.0f | %s |\n",
            key.c_str(), n, ts.front(), P, best.leaf_poly_sq, best.leaf_nonempty,
            best.bsp_cut_slots, best.bsp_cuts_used,
            best.frag_edges_count ? double(best.frag_edges_total) / double(best.frag_edges_count)
                                  : 0.0,
            best.leaf_input_total ? double(best.leaf_planes_total) / double(best.leaf_input_total)
                                  : 0.0,
            best_ms_arr, best.ms_arr_gather * aq, best.ms_arr_present * aq, best.ms_arr_prep * aq,
            best.ms_arr_frag * aq, best.ms_arr_coplanar * aq, best_ms_cls, best_ms_sti,
            (ts.back() / ts.front() <= 1.10) ? "**清浄**" : "**棄却（汚染）**");
    }
    // ---- ★ ドリフトの判定（測定自身に定常を判定させます）------------------------
    if (!drift_key.empty() && drift_last > 0 && drift_last != drift_first) {
        const double d = std::abs(drift_last - drift_first) / std::min(drift_first, drift_last);
        std::printf(
            "\n**ドリフト検査**（`%s`）: 最初 **%.3f s** / 最後 **%.3f s** / "
            "**差 %.1f%%** → **%s**\n",
            drift_key.c_str(), drift_first, drift_last, d * 100,
            (d <= kDriftTol) ? "**定常。採用**" : "**棄却。ウォームアップを延ばして再実行**");
    } else if (!drift_key.empty()) {
        std::printf("\n**ドリフト検査は未実施**（先頭の対を最後にもう一度渡してください）\n");
    }
    return 0;
}

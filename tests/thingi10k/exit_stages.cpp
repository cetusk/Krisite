// Krisite — 出口（`to_mesh`）の段の内訳を、実データで測る
//
// `SPEC-phase5.md` §6「性能（CP4 以降）」の「出口の内訳」の記録項目そのものです。
// **`HANDOVER.md` §9 の 6「三角形化の内訳を測る」**に対応します。
//
// ---
//
// ## 測り方の規律（`CLAUDE.md`）
//
// **第一に、設定と対象を最初に出します。** 何を何件回すかを、回し始める前に。
//
// **第二に、比は同一実行の中で対にします。** 時間帯をまたぐと 1.5〜1.6 倍動くので、
// **旧経路と新経路を同じプロセスの中で交互に回します。**
//
// **第三に、中核（`boolean`）は 1 度だけ回し、その結果を使い回します。**
// 比較対象を出口だけに限るためです。
//
// **第四に、出力が 1 ビットも変わらないことを確かめます。** 検証の経路を変えただけ
// なので、変わってはいけません。
//
// ## 使い方
//
//     exit_stages <model-id> [depth=6] [threads=16] [reps=2]
//     exit_stages <idA>x<idB> [depth=6] [threads=16] [reps=2]
//
// 引数が 1 つの ID なら self-union、`AxB` なら 2 項の 3 演算。
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"

#include "loader.hpp"

using namespace krisite;

namespace {

/// 出口の 1 回分。**時間は段ごとに、構造は演算回数で。**
struct Run {
    double total = 0;
    csg::ToMeshStats st{};
    unsigned long long hash = 0;
};

/// **出力のハッシュ**（バイト一致の判定。三角形と頂点の両方を混ぜます）。
unsigned long long hash_mesh(const csg::SoupMesh& m) {
    unsigned long long h = 1469598103934665603ull;
    const auto mix = [&h](unsigned long long x) {
        h ^= x;
        h *= 1099511628211ull;
    };
    mix(m.triangles.size());
    mix(m.vertices.size());
    for (const mesh::Tri& t : m.triangles) {
        for (int k = 0; k < 3; ++k) mix(t[static_cast<std::size_t>(k)]);
    }
    return h;
}

Run run_once(const csg::PolySoup& soup, const csg::ToMeshOptions& tm) {
    Run r;
    const auto t0 = std::chrono::steady_clock::now();
    const csg::SoupMesh m = csg::to_mesh(soup, tm, &r.st);
    r.total = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    r.hash = hash_mesh(m);
    return r;
}

void print_breakdown(const char* tag, const Run& r) {
    const csg::ToMeshStats& s = r.st;
    std::printf("| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |\n", tag, r.total, s.ms_construct,
                s.ms_merge, s.ms_index, s.ms_tri, s.ms_split);
    // **内部で整合するはずの量を並べます**（`CLAUDE.md`）。
    // **5 段の和と `ms_total` が一致しなければ、計時されていない段があります。**
    // **`ms_total` と壁時計の差は、`to_mesh` の【外】にあります。**
    const double sum = s.ms_construct + s.ms_merge + s.ms_index + s.ms_tri + s.ms_split;
    std::printf(
        "|   └ 検算（ms） | 5 段の和 %.1f | `to_mesh` の全体 %.1f | 差 %+.1f | "
        "壁時計 %.1f | 外側の差 %+.1f | |\n",
        sum, s.ms_total, s.ms_total - sum, r.total * 1000.0, r.total * 1000.0 - s.ms_total);
    const mesh::SplitStats& p = s.split;
    std::printf(
        "|   └ 分裂の内訳（ms） | 辺表 %.1f | 診断 %.1f | 扇 %.1f | 書戻し %.1f | "
        "**検証 %.1f** | |\n",
        p.ms_edges, p.ms_diag, p.ms_fan, p.ms_apply, p.ms_verify);
}

}  // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    if (argc < 2) {
        std::fprintf(stderr, "使い方: exit_stages <id | idAxidB> [depth] [threads] [reps]\n");
        return 2;
    }
    const std::string target = argv[1];
    const unsigned depth = (argc > 2) ? static_cast<unsigned>(std::atoi(argv[2])) : 6;
    const unsigned nthreads = (argc > 3) ? static_cast<unsigned>(std::atoi(argv[3])) : 16;
    const int reps = (argc > 4) ? std::atoi(argv[4]) : 2;

    const std::size_t x = target.find('x');
    const bool pair = (x != std::string::npos);
    const std::string ida = pair ? target.substr(0, x) : target;
    const std::string idb = pair ? target.substr(x + 1) : target;

    // **★ 設定と対象を最初に出します**（`CLAUDE.md`「測定の出力に、設定と
    // 実際に測った対象を必ず書いてください」「対象の数も設定の一部です」）
    std::printf("\n## 出口の段の内訳（`SPEC-phase5.md` §6）\n\n");
    std::printf("| 設定 | 値 |\n|---|---|\n");
    std::printf("| 対象 | `%s`（%s） |\n", target.c_str(), pair ? "2 項 3 演算" : "self-union");
    std::printf("| 深度 | %u |\n", depth);
    std::printf("| スレッド | %u |\n", nthreads);
    std::printf("| 反復 | %d（**最小値を採ります**） |\n", reps);
    std::printf("| b（座標ビット） | %d |\n", KRISITE_COORD_BITS);
    std::printf("| 検算（§5.5） | **ON**（`SPEC-phase5.md` §3.2） |\n\n");

    par::ThreadPool pool(nthreads);
    const auto qa = krithingi::quantize(
        krithingi::load_kmesh("data/thingi10k/kmesh/" + ida + ".kmesh"),
        krithingi::make_transform(1000 + static_cast<unsigned>(std::atoi(ida.c_str())) % 997));
    const auto qb =
        pair ? krithingi::quantize(krithingi::load_kmesh("data/thingi10k/kmesh/" + idb + ".kmesh"),
                                   krithingi::make_transform(
                                       2000 + static_cast<unsigned>(std::atoi(idb.c_str())) % 997))
             : qa;

    csg::BoolOptions o;
    o.depth = depth;
    o.adaptive = true;
    o.cull_planes = true;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    o.threads = nthreads;
    o.pool = &pool;

    const csg::BoolOp ops[3] = {csg::BoolOp::Union, csg::BoolOp::Intersection,
                                csg::BoolOp::Difference};
    const char* op_name[3] = {"∪", "∩", "＼"};
    const int nops = pair ? 3 : 1;
    std::printf("**入力**: `%s` %zu 三角形 / `%s` %zu 三角形\n\n", ida.c_str(),
                qa.mesh.triangles.size(), idb.c_str(), qb.mesh.triangles.size());

    for (int oi = 0; oi < nops; ++oi) {
        const csg::PolySoup A = csg::from_mesh(qa.mesh);
        const csg::PolySoup B = csg::from_mesh(qb.mesh);
        csg::BoolStats bs;
        const auto tc0 = std::chrono::steady_clock::now();
        const csg::PolySoup soup = csg::boolean(A, B, ops[oi], o, &bs);
        const double core =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count();

        std::printf("### %s — 中核 %.3f s（多角形 %zu）\n\n", op_name[oi], core, soup.polys.size());

        csg::ToMeshOptions fast;
        fast.threads = nthreads;
        fast.pool = &pool;
        fast.verify_split_delta = true;  // §3.2
        // **段の内訳を採ります**（既定は偽。計測の費用を本番に持ち込まないため）
        fast.time_stages = true;
        // **候補の区間の絞り込み（案 (b2)）を外した側**（`DESIGN` §13）。**比較の基準**
        csg::ToMeshOptions unsorted = fast;
        unsorted.sort_candidates = false;
        csg::ToMeshOptions naive = fast;
        naive.verify_split_naive = true;
        // **T 解決の走査順の比較**（`SPEC-phase5.md` §5.11）。**判定は同じ。順序だけ**
        csg::ToMeshOptions per_edge = fast;
        per_edge.scan_per_edge = true;

        // **交互に回します**（`CLAUDE.md`「測定の順序が対象の大きさと相関しない
        // ようにしてください」の同じ理由。1 回目だけが冷えている効果を分散させる）
        Run best_fast, best_naive, best_pe;
        for (int i = 0; i < reps; ++i) {
            const Run f = run_once(soup, fast);
            const Run n = run_once(soup, naive);
            const Run e = run_once(soup, per_edge);
            if (i == 0 || f.total < best_fast.total) best_fast = f;
            if (i == 0 || n.total < best_naive.total) best_naive = n;
            if (i == 0 || e.total < best_pe.total) best_pe = e;
        }

        std::printf(
            "| 経路 | 合計 s | 構成点 | 値で併合 | 平面索引 | 扇（T+三角形化） | "
            "接触の分裂 |\n|---|---:|---:|---:|---:|---:|---:|\n");
        print_breakdown("**増分計算（既定）**", best_fast);
        print_breakdown("従来（check_topology×2）", best_naive);
        std::printf("\n");

        // **出力は 1 ビットも変わってはいけません**
        std::printf("| 出力のハッシュ | %s |\n",
                    best_fast.hash == best_naive.hash ? "**一致**" : "**★ 食い違い（重大）**");
        std::printf("| 検証の時間の比 | %.2f 倍 |\n",
                    best_fast.st.split.ms_verify > 0
                        ? best_naive.st.split.ms_verify / best_fast.st.split.ms_verify
                        : 0.0);
        std::printf("| 出口の合計の比 | %.2f 倍 |\n",
                    best_fast.total > 0 ? best_naive.total / best_fast.total : 0.0);
        std::printf("| `unresolved`（増分 / 従来） | %zu / %zu |\n", best_fast.st.split.unresolved,
                    best_naive.st.split.unresolved);
        std::printf("| うち事後の非多様体 | %zu / %zu |\n", best_fast.st.split.unresolved_post,
                    best_naive.st.split.unresolved_post);
        std::printf("| ΔV（増分 / 従来） | %zu / %zu |\n", best_fast.st.split.actual_delta_v,
                    best_naive.st.split.actual_delta_v);
        std::printf("| ΔE（増分 / 従来） | %zu / %zu |\n", best_fast.st.split.actual_delta_e,
                    best_naive.st.split.actual_delta_e);
        std::printf("| 前提の破れで退避 | %zu |\n", best_fast.st.split.verify_fallback);
        std::printf("| `degenerate_kept` | %zu |\n", best_fast.st.t.degenerate_kept);
        // **A-3 の群の規模**（案「軸で整列して二分探索」の整列の費用を見積もるため）。
        // **整列は群ごとに 1 回で、その群の全多角形で使い回せます。**
        std::printf("| A-3 の群 | %zu 群 / 多角形 %zu（群あたり %.1f 多角形） |\n",
                    best_fast.st.cell_index_groups, soup.polys.size(),
                    best_fast.st.cell_index_groups > 0
                        ? static_cast<double>(soup.polys.size()) / best_fast.st.cell_index_groups
                        : 0.0);
        std::printf("| 三角形 / 頂点 | %zu / %zu |\n", best_fast.st.t.general_used,
                    best_fast.st.merged_points);
        // **`ms_tri` の中身**（`HANDOVER.md` §9 の 6「三角形化の内訳を測る」）。
        //
        // **★ これは CPU 時間です。** 壁時計（`ms_tri`）と直接比べないこと
        // （`CLAUDE.md`「CPU 時間と壁時計を混ぜないでください」）。
        // **比率としてなら読めます** — 同じ並列区間の中の 2 段だからです。
        {
            const double a = best_fast.st.t.ms_insert_t, b = best_fast.st.t.ms_fan_tri;
            std::printf(
                "| **扇の内訳（CPU 時間 ms）** | T 解決 %.1f / 三角形化 %.1f"
                "（%.0f%% : %.0f%%） |\n",
                a, b, (a + b > 0) ? 100.0 * a / (a + b) : 0.0,
                (a + b > 0) ? 100.0 * b / (a + b) : 0.0);
            std::printf(
                "| 壁時計の `ms_tri` と、CPU 時間の和 | %.1f ms / %.1f ms"
                "（並列効率 %.2f） |\n",
                best_fast.st.ms_tri, a + b,
                best_fast.st.ms_tri > 0 ? (a + b) / best_fast.st.ms_tri : 0.0);
            std::printf(
                "| T 解決の演算回数 | 走査した辺 %zu / 索引の候補 %zu / "
                "挿入 %zu / 1 辺の最大 %zu |\n",
                best_fast.st.t.edges_scanned, best_fast.st.t.candidates, best_fast.st.t.inserted,
                best_fast.st.t.max_per_edge);
            // **★ 述語の評価回数**（`CLAUDE.md`「効果は演算回数で測ってください」）。
            // **候補数ではなく、これが照合の本体の費用です。**
            std::printf(
                "| **述語の評価回数** | `side` %zu / `strictly_between` %zu / "
                "候補集合の走査 %zu |\n",
                best_fast.st.t.side_tests, best_fast.st.t.between_tests, best_fast.st.t.cand_scans);
            std::printf("| 一般解 | 作った三角形 %zu / 従来へ落ちた多角形 %zu |\n",
                        best_fast.st.t.general_used, best_fast.st.t.general_fallback);
            // **走査順の新旧を、同じ実行の中で比べます**（§5.11）
            const auto& pe = best_pe.st.t;
            std::printf(
                "| **走査順（辺ごと → 多角形あたり 1 回）** | 走査 %zu → %zu（%.2f 倍）"
                " / `side` %zu → %zu（%.3f 倍） |\n",
                pe.cand_scans, best_fast.st.t.cand_scans,
                best_fast.st.t.cand_scans > 0
                    ? static_cast<double>(pe.cand_scans) / best_fast.st.t.cand_scans
                    : 0.0,
                pe.side_tests, best_fast.st.t.side_tests,
                best_fast.st.t.side_tests > 0
                    ? static_cast<double>(pe.side_tests) / best_fast.st.t.side_tests
                    : 0.0);
            std::printf(
                "| 同、扇の CPU 時間と出口の壁時計 | %.1f → %.1f ms / %.3f → %.3f s "
                "（%s） |\n",
                pe.ms_insert_t + pe.ms_fan_tri, a + b, best_pe.total, best_fast.total,
                best_pe.hash == best_fast.hash ? "**ハッシュ一致**" : "**★ 食い違い（重大）**");
        }
        // ---- 案 (b2) の効果（`DESIGN` §13）--------------------------------------
        {
            const Run bx = run_once(soup, unsorted);
            const auto& t = bx.st.t;
            std::printf(
                "| **案 (b2) 候補の区間の絞り込み** | 飛ばした %zu / `side` %zu → %zu"
                "（**%.2f 倍**） |\n",
                best_fast.st.t.cand_skipped_by_range, t.side_tests, best_fast.st.t.side_tests,
                best_fast.st.t.side_tests > 0
                    ? static_cast<double>(t.side_tests) / best_fast.st.t.side_tests
                    : 0.0);
            std::printf(
                "| 同、扇の CPU 時間と出口の壁時計 | %.1f → %.1f ms / %.3f → %.3f s"
                "（%s） |\n",
                t.ms_insert_t + t.ms_fan_tri,
                best_fast.st.t.ms_insert_t + best_fast.st.t.ms_fan_tri, bx.total, best_fast.total,
                bx.hash == best_fast.hash ? "**ハッシュ一致**" : "**★ 食い違い（重大）**");
            std::printf("| 同、挿入した T 頂点 | %zu → %zu（**1 個も変わってはいけない**） |\n",
                        t.inserted, best_fast.st.t.inserted);
            // **整列の費用は索引の段（`ms_index`）に乗ります。** 同一実行の中で比べます
            std::printf("| **整列の費用**（索引の段） | %.1f → %.1f ms（差 %+.1f ms） |\n",
                        bx.st.ms_index, best_fast.st.ms_index,
                        best_fast.st.ms_index - bx.st.ms_index);
            std::printf("| 同、扇の段（壁時計） | %.1f → %.1f ms |\n", bx.st.ms_tri,
                        best_fast.st.ms_tri);
            std::printf(
                "| 二分探索に払う `cmp_h` | %zu（多角形あたり %.1f 回） |\n",
                best_fast.st.t.box_cmp_tests,
                best_fast.st.cell_index_groups > 0
                    ? static_cast<double>(best_fast.st.t.box_cmp_tests) / best_fast.st.t.cand_scans
                    : 0.0);
        }
        std::printf("\n");
    }
    return 0;
}

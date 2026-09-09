// Krisite — 構成点キャッシュの実装の A/B を、40 対で測る（`SPEC-phase5.md` §5.10.12.4）
//
// ## なぜ 40 対か
//
// **`std::map` → 開番地法の効果は、1 対 2 模型でしか確かめていません**
// （CPU 1.30 / 1.34 倍）。**分布が見えません。**
//
// ## ★ 条件 1 — 同一実行の中で A/B にします
//
// **40 対を 2 回回すと、時間帯の交絡（`BENCH.md` の 1.5〜1.6 倍）が入ります。**
// **旗を演算ごとに切り替えて、同じプロセスで交互に回します。**
//
// **そして先に回す側を演算ごとに入れ替えます。**
// **片方を必ず先に回すと、機械の温まりが一方に偏ります。**
//
// ## ★ 条件 2 — 対ごとの比の【分布】を出します
//
// **中央値だけでなく最小と最大。** **1.0 倍の対があるかどうかが知りたいことです。**
// **効かない対は「探索の回数が少ない」で説明できるはずです。説明できなければ別の要因です。**
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/par/thread_pool.hpp"

#include "thingi10k/loader.hpp"

using namespace krisite;

namespace {

unsigned long long hash_soup(const csg::PolySoup& s) {
    unsigned long long h = 1469598103934665603ull;
    const auto mix = [&h](unsigned long long v) {
        h ^= v;
        h *= 1099511628211ull;
    };
    mix(s.polys.size());
    for (const csg::Poly& q : s.polys) {
        mix(q.frag.support);
        mix(q.frag.flipped ? 1u : 0u);
        mix(q.src);
        mix(q.tag);
        for (csg::PlaneId e : q.frag.edge) mix(e);
    }
    return h;
}

struct PairResult {
    std::string key;
    double cpu_map = 0, cpu_hash = 0;
    double wall_map = 0, wall_hash = 0;
    std::size_t lookups = 0, hits = 0;
    std::size_t tris = 0;
    bool same = true;
};

}  // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    const std::string base = (argc > 1) ? argv[1] : "data/thingi10k/ember_o3";
    const unsigned depth = (argc > 2) ? static_cast<unsigned>(std::atoi(argv[2])) : 6;
    const unsigned nthreads = (argc > 3) ? static_cast<unsigned>(std::atoi(argv[3])) : 8;

    // ---- 一覧を読む（**添字がそのまま変換の種**。`thingi_cp1.cpp` と同じ規約）------
    std::vector<std::string> ids;
    {
        std::ifstream f(base + ".txt");
        std::string id;
        long long faces = 0;
        while (f >> id >> faces) ids.push_back(id);
    }
    std::map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (index.find(ids[i]) == index.end()) index.emplace(ids[i], i);
    }
    std::vector<std::pair<std::string, std::string>> pairs;
    {
        std::ifstream f(base + "_only.txt");
        std::string line;
        while (std::getline(f, line)) {
            const std::size_t x = line.find('x');
            if (x == std::string::npos) continue;
            pairs.emplace_back(line.substr(0, x), line.substr(x + 1));
        }
    }

    // **★ 設定と対象の数を先に出します**（`CLAUDE.md`）
    std::printf("\n## 構成点キャッシュの A/B — 40 対（`SPEC-phase5.md` §5.10.12.4）\n\n");
    std::printf("| 設定 | 値 |\n|---|---|\n");
    std::printf("| 一覧 | `%s.txt`（%zu 模型）/ `%s_only.txt` |\n", base.c_str(), ids.size(),
                base.c_str());
    std::printf("| 深度 | %u（適応） |\n", depth);
    std::printf("| スレッド | %u |\n", nthreads);
    std::printf("| b | %d |\n", KRISITE_COORD_BITS);
    std::printf("| 演算 | 3（∪ / ∩ / ＼） |\n");
    std::printf("| A/B | **同一プロセス。演算ごとに先攻を入れ替える** |\n");
    std::printf("| 測るもの | 中核（`boolean`）のみ。出口は含みません |\n\n");
    std::printf("**これから %zu 対を回します**（%zu 演算 × 2 実装 = %zu 回）\n\n", pairs.size(),
                pairs.size() * 3, pairs.size() * 6);
    if (pairs.empty()) {
        std::fprintf(stderr, "対が 0 件です。一覧を確かめてください\n");
        return 2;
    }

    par::ThreadPool pool(nthreads);
    csg::BoolOptions base_o;
    base_o.depth = depth;
    base_o.adaptive = true;
    base_o.cull_planes = true;
    base_o.early_out = true;
    base_o.cache_points = true;
    base_o.local_bsp = true;
    base_o.split_contacts = true;
    base_o.threads = nthreads;
    base_o.pool = &pool;

    std::vector<PairResult> res;
    std::size_t done = 0;
    for (const auto& pr : pairs) {
        const auto ia = index.find(pr.first), ib = index.find(pr.second);
        if (ia == index.end() || ib == index.end()) {
            std::printf("  一覧に無いので飛ばす %sx%s\n", pr.first.c_str(), pr.second.c_str());
            continue;
        }
        const auto qa = krithingi::quantize(
            krithingi::load_kmesh("data/thingi10k/kmesh/" + pr.first + ".kmesh"),
            krithingi::make_transform(1000 + ia->second));
        const auto qb = krithingi::quantize(
            krithingi::load_kmesh("data/thingi10k/kmesh/" + pr.second + ".kmesh"),
            krithingi::make_transform(1000 + ib->second));
        if (qa.mesh.triangles.empty() || qb.mesh.triangles.empty()) {
            std::printf("  読めないので飛ばす %sx%s\n", pr.first.c_str(), pr.second.c_str());
            continue;
        }
        PairResult r;
        r.key = pr.first + "x" + pr.second;
        r.tris = qa.mesh.triangles.size() + qb.mesh.triangles.size();
        const csg::PolySoup A = csg::from_mesh(qa.mesh), B = csg::from_mesh(qb.mesh);
        const csg::BoolOp ops[3] = {csg::BoolOp::Union, csg::BoolOp::Intersection,
                                    csg::BoolOp::Difference};
        std::printf("  → 開始 %s（入力 %zu 三角形）\n", r.key.c_str(), r.tris);
        std::fflush(stdout);
        for (int oi = 0; oi < 3; ++oi) {
            unsigned long long h[2] = {0, 0};
            // **★ 先攻を演算ごとに入れ替えます**（機械の温まりが一方に偏らないように）
            const int first = oi % 2;
            for (int t = 0; t < 2; ++t) {
                const int use_map = (t == 0) ? first : (1 - first);
                csg::BoolOptions o = base_o;
                o.point_cache_map = (use_map != 0);
                csg::BoolStats st;
                const auto w0 = std::chrono::steady_clock::now();
                const std::clock_t c0 = std::clock();
                const csg::PolySoup s = csg::boolean(A, B, ops[oi], o, &st);
                const double cs = static_cast<double>(std::clock() - c0) / CLOCKS_PER_SEC;
                const double ws =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
                h[use_map] = hash_soup(s);
                if (use_map) {
                    r.cpu_map += cs;
                    r.wall_map += ws;
                } else {
                    r.cpu_hash += cs;
                    r.wall_hash += ws;
                    r.lookups += st.cache_hits + st.cache_misses;
                    r.hits += st.cache_hits;
                }
            }
            if (h[0] != h[1]) r.same = false;
        }
        res.push_back(r);
        ++done;
        std::printf("    完了 %s: CPU %.3f → %.3f s（**%.2f 倍**）/ 探索 %zu / %s\n", r.key.c_str(),
                    r.cpu_map, r.cpu_hash, r.cpu_hash == 0 ? 0.0 : r.cpu_map / r.cpu_hash,
                    r.lookups, r.same ? "出力一致" : "**★ 食い違い（重大）**");
        std::fflush(stdout);
    }

    // ---- 分布 ---------------------------------------------------------------
    std::printf("\n### 結果（%zu 対）\n\n", res.size());
    std::printf(
        "| 対 | 入力三角形 | CPU `std::map` | CPU 開番地 | **比** | 壁時計の比 | "
        "探索 | 命中率 | 出力 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|---:|---|\n");
    std::vector<double> ratios;
    std::size_t mismatch = 0;
    for (const PairResult& r : res) {
        const double ratio = r.cpu_hash == 0 ? 0.0 : r.cpu_map / r.cpu_hash;
        const double wratio = r.wall_hash == 0 ? 0.0 : r.wall_map / r.wall_hash;
        ratios.push_back(ratio);
        if (!r.same) ++mismatch;
        std::printf("| `%s` | %zu | %.3f s | %.3f s | **%.2f** | %.2f | %zu | %.1f%% | %s |\n",
                    r.key.c_str(), r.tris, r.cpu_map, r.cpu_hash, ratio, wratio, r.lookups,
                    r.lookups == 0
                        ? 0.0
                        : 100.0 * static_cast<double>(r.hits) / static_cast<double>(r.lookups),
                    r.same ? "一致" : "**食い違い**");
    }
    if (!ratios.empty()) {
        std::vector<double> s = ratios;
        std::sort(s.begin(), s.end());
        double logsum = 0;
        for (double v : s) logsum += std::log(v);
        std::printf("\n| 統計 | 値 |\n|---|---:|\n");
        std::printf("| 対の数 | %zu |\n", s.size());
        std::printf("| **最小** | **%.2f 倍** |\n", s.front());
        std::printf("| 第 1 四分位 | %.2f 倍 |\n", s[s.size() / 4]);
        std::printf("| **中央値** | **%.2f 倍** |\n", s[s.size() / 2]);
        std::printf("| 第 3 四分位 | %.2f 倍 |\n", s[s.size() * 3 / 4]);
        std::printf("| **最大** | **%.2f 倍** |\n", s.back());
        std::printf("| 幾何平均 | %.2f 倍 |\n", std::exp(logsum / static_cast<double>(s.size())));
        std::printf("| **出力の食い違い** | **%zu 対**（0 でなければ重大） |\n", mismatch);
    }
    return mismatch == 0 ? 0 : 1;
}

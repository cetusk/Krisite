// Krisite — 出力の次数 3 以上の辺で、面がどの source から来たか（`SPEC-phase5.md` (d″)）
//
// **相関ではなく機構を見ます。** 失敗 7 対は「自己交差する模型を含む」ことが
// $p = 4.6\times10^{-7}$ で示されましたが、**それは相関です。**
//
//   4 枚が A/A または B/B  →  **自己接触**。機構の理解が正しい
//   A–B 混在               →  **2 立体の接触**。理解が誤っている
//
// ## ★ 変換は `cp1.txt` の添字から決めます
//
// **`cp1.cpp` は `prepare(raw, 1000 + i)` を使います**（$i$ は一覧での添字）。
// **固定値を使うと別の入力を測ることになります** — 実際に踏みました
// （`IMPL-phase5.md` §33）。**一覧を読んで添字を引くのが唯一の正しい方法です。**
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
using namespace krisite::csg;

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    const std::string list = "data/thingi10k/cp1.txt";
    std::map<std::string, std::size_t> index;
    {
        std::ifstream f(list);
        std::string id, nf;
        std::size_t i = 0;
        while (f >> id >> nf) index[id] = i++;
    }
    // **出力はスレッド数に依りません**（`SPEC-phase4.md` §7.1。
    // 実データでも逐次 / 4 / 16 が完全一致することを確認済み。`IMPL-phase5.md` §33.2）。
    // **速さのためだけに使います。**
    const unsigned nthreads = (std::getenv("KRI_THREADS") != nullptr)
                                  ? static_cast<unsigned>(std::atoi(std::getenv("KRI_THREADS")))
                                  : 16u;
    par::ThreadPool pool(nthreads);

    std::printf(
        "| 対 | 演算 | 次数>=3 の辺 | **A/A** | **B/B** | **A-B 混在** | 余分辺 | 最大次数 |\n");
    std::printf("|---|---|---:|---:|---:|---:|---:|---:|\n");

    for (int k = 1; k < argc; ++k) {
        const std::string key = argv[k];
        const std::size_t x = key.find('x');
        const std::string ida = key.substr(0, x), idb = key.substr(x + 1);
        if (index.count(ida) == 0 || index.count(idb) == 0) {
            std::printf("| %s | **一覧に無い** |\n", key.c_str());
            continue;
        }
        // **添字から変換を作ります**（`cp1.cpp` と同じ）
        const auto qa =
            krithingi::quantize(krithingi::load_kmesh("data/thingi10k/kmesh/" + ida + ".kmesh"),
                                krithingi::make_transform(1000 + index[ida]));
        const auto qb =
            krithingi::quantize(krithingi::load_kmesh("data/thingi10k/kmesh/" + idb + ".kmesh"),
                                krithingi::make_transform(1000 + index[idb]));
        const PolySoup A = from_mesh(qa.mesh), B = from_mesh(qb.mesh);
        const char* nm[3] = {"∪", "∩", "＼"};
        int oi = 0;
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            BoolOptions o;
            o.depth = 6;
            o.adaptive = true;
            o.cull_planes = true;
            o.early_out = true;
            o.cache_points = true;
            o.local_bsp = true;
            o.split_contacts = true;
            o.threads = nthreads;
            o.pool = &pool;
            ToMeshOptions tm;
            tm.split_contacts = true;
            tm.threads = nthreads;
            tm.pool = &pool;
            const SoupMesh m = to_mesh(boolean(A, B, op, o), tm);
            std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<int>> e;
            for (std::size_t i = 0; i < m.triangles.size(); ++i) {
                const mesh::Tri& t = m.triangles[i];
                for (int j = 0; j < 3; ++j) {
                    auto u = t[static_cast<std::size_t>(j)];
                    auto v = t[static_cast<std::size_t>((j + 1) % 3)];
                    if (v < u) std::swap(u, v);
                    e[{u, v}].push_back(i < m.tri_src.size() ? m.tri_src[i] : -1);
                }
            }
            std::size_t deg = 0, aa = 0, bb = 0, mix = 0;
            for (const auto& kv : e) {
                if (kv.second.size() < 3) continue;
                ++deg;
                bool h0 = false, h1 = false;
                for (int s : kv.second) {
                    if (s == 0) {
                        h0 = true;
                    } else if (s == 1) {
                        h1 = true;
                    }
                }
                if (h0 && h1) {
                    ++mix;
                } else if (h0) {
                    ++aa;
                } else {
                    ++bb;
                }
            }
            const auto tp = mesh::check_topology(m.triangles);
            std::printf("| %s | %s | %zu | %zu | %zu | **%zu** | %zu | %zu |\n",
                        (oi == 0 ? key.c_str() : ""), nm[oi], deg, aa, bb, mix, tp.edges_excess,
                        tp.max_edge_degree);
            ++oi;
        }
    }
    return 0;
}

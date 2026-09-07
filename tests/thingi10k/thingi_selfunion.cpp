// Krisite — 自己交差する模型の self-union（`SPEC-phase5.md` §2.9.1.5）
//
// **★ この実行ファイルは `build/tests/*` とは別物です**（実データのドライバ）。
//
// **94 模型の検査は「CP2 で失敗した 47 対に現れる模型」を対象にしており、
// 失敗した対から作った集合なので self-union 一般の検査ではありませんでした**
// （`IMPL-phase5.md` §84.4）。ここでは **CP1 の 1,000 模型のうち、量子化後に
// 自己交差する 70 模型**を、性質で選んで回します。
//
// **失敗の有無だけでなく、`135071` と同じ形かを判定できる項目を採ります**（§85.4）。
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/geom/predicates.hpp"
#include "krisite/mesh/self_intersect.hpp"
#include "krisite/mesh/topology.hpp"

#include "thingi10k/loader.hpp"

using namespace krisite;

namespace {

/// 非多様体の構造（§85.4）。**`135071` と同じ形かを判定できるように採ります。**
struct Shape {
    std::size_t excess_edges = 0;  ///< 次数 3 以上の無向辺
    std::size_t verts = 0;         ///< その辺に関与する頂点
    std::size_t on_input = 0;      ///< うち入力の頂点そのもの（残りは構成点）
    std::size_t components = 0;    ///< 過剰辺だけで作った部分グラフの連結成分
    std::size_t deg1 = 0;          ///< 端点（次数 1）。**0 なら環、非零なら鎖**
    std::size_t max_deg = 0;       ///< 部分グラフの頂点次数の最大
};

Shape analyze(const csg::SoupMesh& m, const std::vector<geom::HPointD>& sorted_input) {
    Shape s;
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::size_t> deg;
    for (const auto& t : m.triangles)
        for (int e = 0; e < 3; ++e) {
            std::uint32_t x = t[e], y = t[(e + 1) % 3];
            if (x > y) std::swap(x, y);
            ++deg[{x, y}];
        }
    std::map<std::uint32_t, std::vector<std::uint32_t>> adj;
    for (const auto& kv : deg) {
        if (kv.second <= 2) continue;
        ++s.excess_edges;
        adj[kv.first.first].push_back(kv.first.second);
        adj[kv.first.second].push_back(kv.first.first);
    }
    s.verts = adj.size();
    for (const auto& kv : adj) {
        s.max_deg = std::max(s.max_deg, kv.second.size());
        if (kv.second.size() == 1) ++s.deg1;
        if (std::binary_search(sorted_input.begin(), sorted_input.end(), m.vertices[kv.first],
                               geom::lex_less))
            ++s.on_input;
    }
    // 連結成分
    std::set<std::uint32_t> seen;
    for (const auto& kv : adj) {
        if (seen.count(kv.first)) continue;
        ++s.components;
        std::vector<std::uint32_t> st{kv.first};
        seen.insert(kv.first);
        while (!st.empty()) {
            const std::uint32_t v = st.back();
            st.pop_back();
            for (std::uint32_t w : adj[v])
                if (seen.insert(w).second) st.push_back(w);
        }
    }
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    const std::string list = (argc > 1) ? argv[1] : "data/thingi10k/selfint_b21.txt";
    const unsigned depth = (argc > 2) ? static_cast<unsigned>(std::atoi(argv[2])) : 6;
    const unsigned nthreads = (argc > 3) ? static_cast<unsigned>(std::atoi(argv[3])) : 8;
    const std::string out_path = "data/thingi10k/selfunion_results.txt";

    // **量子化の変換は CP1 と同じ**（`cp1.txt` の並び順が種）。**違えると別の入力になります。**
    std::map<std::string, std::size_t> idx;
    {
        std::ifstream f("data/thingi10k/cp1.txt");
        std::string a, b;
        std::size_t i = 0;
        while (f >> a >> b) idx[a] = i++;
    }
    std::vector<std::string> ids;
    {
        std::ifstream f(list);
        std::string id;
        while (f >> id)
            if (!id.empty() && id[0] != '#') ids.push_back(id);
    }
    // 済みの模型は飛ばす（再開できるように）
    std::set<std::string> done;
    {
        std::ifstream f(out_path);
        std::string line;
        while (std::getline(f, line))
            if (!line.empty()) done.insert(line.substr(0, line.find(' ')));
    }
    std::size_t planned = 0;
    for (const std::string& id : ids)
        if (!done.count(id)) ++planned;

    // **★ 設定と対象の数を必ず出します**（`IMPL-phase5.md` §82）
    std::printf("\n## self-union（自己交差する模型）\n\n");
    std::printf("| 設定 | 値 |\n|---|---|\n");
    std::printf("| 一覧 | `%s` |\n", list.c_str());
    std::printf("| 記録先 | `%s` |\n", out_path.c_str());
    std::printf("| 深度 | %u |\n", depth);
    std::printf("| スレッド | %u |\n", nthreads);
    std::printf("| b（座標ビット） | %d |\n", KRISITE_COORD_BITS);
    std::printf("| **NSI** | **宣言しない**（対象が自己交差しているので宣言できません） |\n");
    std::printf("| 接触の分裂 | 有効 |\n\n");
    std::printf("**これから %zu 模型を回します**（一覧 %zu / 既済 %zu）\n\n", planned, ids.size(),
                done.size());

    par::ThreadPool pool(nthreads);
    std::ofstream out(out_path, std::ios::app);
    const auto t_all = std::chrono::steady_clock::now();
    std::size_t ok = 0, bad = 0, n = 0;
    for (const std::string& id : ids) {
        if (done.count(id)) continue;
        ++n;
        const auto q =
            krithingi::quantize(krithingi::load_kmesh("data/thingi10k/kmesh/" + id + ".kmesh"),
                                krithingi::make_transform(1000 + idx.at(id)));
        std::printf("  [%zu/%zu] %s（三角形 %zu）… ", n, planned, id.c_str(),
                    q.mesh.triangles.size());
        std::fflush(stdout);
        const auto t0 = std::chrono::steady_clock::now();
        // self-union = **同じスープを 2 度使う和**（`test_winding.cpp` の `self_union` と同形）
        csg::PolySoup a = csg::from_mesh(q.mesh), b = csg::from_mesh(q.mesh);
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
        csg::BoolStats st;
        const csg::PolySoup u = csg::boolean(a, b, csg::BoolOp::Union, o, &st);
        csg::ToMeshStats ts;
        const csg::SoupMesh m = csg::to_mesh(u, {}, &ts);
        const double sec =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const mesh::TopologyReport r = mesh::check_topology(m.triangles);
        std::vector<geom::HPointD> inv;
        inv.reserve(q.mesh.vertices.size());
        for (const auto& v : q.mesh.vertices) inv.push_back(geom::to_homogeneous(v));
        std::sort(inv.begin(), inv.end(), geom::lex_less);
        const Shape sh = analyze(m, inv);
        const bool good = r.edge_manifold && sh.excess_edges == 0;
        if (good)
            ++ok;
        else
            ++bad;
        std::printf("%s（%.1f s、三角形 %zu、**過剰辺 %zu**、unresolved %zu）\n",
                    good ? "ok" : "**FAIL**", sec, m.triangles.size(), sh.excess_edges,
                    ts.split.unresolved);
        out << id << ' ' << (good ? "ok" : "FAIL") << ' ' << q.mesh.triangles.size() << ' '
            << m.triangles.size() << ' ' << sec << ' ' << sh.excess_edges << ' '
            << r.edges_deficient << ' ' << r.max_edge_degree << ' ' << ts.split.unresolved << ' '
            << sh.verts << ' ' << sh.on_input << ' ' << sh.components << ' ' << sh.deg1 << ' '
            << sh.max_deg << ' ' << r.components << ' ' << r.chi << ' ' << st.regions_negative_w
            << ' ' << st.regions_w_ge2 << '\n';
        out.flush();
    }
    const double all =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_all).count();
    std::printf("\n**模型 %zu / ok %zu / FAIL %zu / %.1f s**\n", n, ok, bad, all);
    return 0;
}

// Krisite — Phase 5: 量子化のモデル単位の統計（`SPEC-phase5.md` §1.5.0 / §6）
//
// **ブール演算を回しません。** 量子化と位相検査だけです（全 1,000 件で数分）。
//
// **狙いは 2 つ。**
//
//   1. 失敗の分類。**併合された頂点は、離れた 2 点を同じ格子点に落とすので、
//      三角形を退化させずにメッシュを自己接触させます。** これは既知の失敗
//      （ケース 24）と同じ配置です。落ちた三角形はその副産物にすぎません
//   2. **$b$ の設計の根拠**（§1.0）。$b$ を上げれば併合は減るはずです。
//      減り方を実データで測れば、上界計算だけで支えてきた選択に根拠が付きます
//
// **「入力は多様体」は選択の時点の性質です。** メタデータで多様体と分類された
// モデルが、**量子化後も多様体とは限りません。** その差を数字にします。
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "krisite/mesh/self_intersect.hpp"
#include "krisite/mesh/topology.hpp"

#include "thingi10k/loader.hpp"

using namespace krisite;

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    const std::string list = (argc > 1) ? argv[1] : "data/thingi10k/cp1.txt";

    std::vector<std::string> ids;
    std::vector<std::size_t> nf;
    {
        std::ifstream f(list);
        std::string id;
        std::size_t n = 0;
        while (f >> id >> n) {
            ids.push_back(id);
            nf.push_back(n);
        }
    }
    std::printf("# 量子化の統計（b=%d、%zu 件）\n", KRISITE_COORD_BITS, ids.size());
    std::printf(
        "id 元面数 量子化後 落ち 併合 ∂S=0 辺多様体 頂点多様体 成分 χ 余分辺 不足辺"
        " 自己交差\n");

    std::size_t n_merged = 0, n_dropped = 0, n_nonmanifold = 0, n_boundary = 0, n_si = 0;
    std::size_t merged_total = 0;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        const krithingi::RawMesh raw =
            krithingi::load_kmesh("data/thingi10k/kmesh/" + ids[i] + ".kmesh");
        // **seed は cp1.cpp と同じ**（1000 + 添字）。違う変換だと別の話になります
        const krithingi::Quantized q =
            krithingi::quantize(raw, krithingi::make_transform(1000 + i));
        const mesh::TopologyReport t = mesh::check_topology(q.mesh.triangles);
        // **自己交差**（`SPEC-phase3.md` §5.6）。**変換が cp1.cpp と同じ（1000 + 添字）で
        // なければ意味がありません** — 別の変換では量子化の結果が変わります
        const bool si = mesh::is_self_intersecting(q.mesh);
        const bool bz = mesh::boundary_is_zero(q.mesh);
        const bool nm = !(t.edge_manifold && t.vertex_manifold);
        if (q.merged_vertices > 0) ++n_merged;
        if (q.dropped_degenerate > 0) ++n_dropped;
        if (nm) ++n_nonmanifold;
        if (si) ++n_si;
        if (!bz) ++n_boundary;
        merged_total += q.merged_vertices;
        std::printf("%s %zu %zu %zu %zu %d %d %d %zu %lld %zu %zu %d\n", ids[i].c_str(), nf[i],
                    q.mesh.triangles.size(), q.dropped_degenerate, q.merged_vertices, (int)bz,
                    (int)t.edge_manifold, (int)t.vertex_manifold, t.components, t.chi,
                    t.edges_excess, t.edges_deficient, (int)si);
        if ((i + 1) % 200 == 0) std::fprintf(stderr, "  %zu / %zu\n", i + 1, ids.size());
    }
    std::printf("\n# 集計（b=%d）\n", KRISITE_COORD_BITS);
    std::printf("# 件数 %zu / 併合ありのモデル %zu / 落ちありのモデル %zu\n", ids.size(), n_merged,
                n_dropped);
    std::printf("# 併合された頂点（延べ） %zu\n", merged_total);
    std::printf("# **量子化後に非多様体** %zu / ∂S≠0 %zu / **自己交差 %zu**\n", n_nonmanifold,
                n_boundary, n_si);
    return 0;
}

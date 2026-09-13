// Krisite — **入力そのものを `split_contacts` に掛ける**（`IMPL-phase5.md` §98 の再検討）
//
// **狙い**: 出口に残った次数 4 の辺が、**入力から受け継いだ形かどうか**を分ける。
//
// self-union で同じ本数が再現したので（8/8）、次の問いは
// 「**入力の非多様体な辺のうち、どれが解けないのか**」です。
// ブール演算を回さず、量子化した入力を直接 `split_contacts` に渡します。
//
// **量子化の変換の種は一覧での並び順で決まります**（`thingi_cp1.cpp` と同じ規則）。
// **種の一覧を間違えると、別の入力を測ることになります。**
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "krisite/geom/plane.hpp"
#include "krisite/mesh/split.hpp"
#include "krisite/mesh/topology.hpp"

#include "thingi10k/loader.hpp"

using namespace krisite;

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    const std::string list = (argc > 1) ? argv[1] : "data/thingi10k/selfint_b21.txt";
    const std::string seed_list = (argc > 2) ? argv[2] : "data/thingi10k/cp1.txt";

    std::map<std::string, std::size_t> idx;
    {
        std::ifstream f(seed_list);
        std::string a, b;
        std::size_t i = 0;
        while (f >> a >> b) idx[a] = i++;
    }
    std::vector<std::string> ids;
    {
        std::ifstream f(list);
        std::string id;
        while (f >> id) {
            if (!id.empty() && id[0] != '#') ids.push_back(id);
        }
    }
    std::printf("\n## 入力を直接 `split_contacts` に掛ける\n\n");
    std::printf("| 設定 | 値 |\n|---|---|\n");
    std::printf("| 一覧 | `%s` |\n", list.c_str());
    std::printf("| 種の一覧 | `%s` |\n", seed_list.c_str());
    std::printf("| b | %d |\n\n", KRISITE_COORD_BITS);
    std::printf("**これから %zu 模型を回します**\n\n", ids.size());
    std::printf(
        "| 模型 | 三角形 | 過剰辺 | radial 試行 | radial 解決 | 組めず | 分裂後に非多様体 |");
    std::printf(" 解けずに残った辺 |\n|---|---:|---:|---:|---:|---:|---:|---:|\n");

    for (const std::string& id : ids) {
        if (idx.count(id) == 0) {
            std::printf("| %s | **種の一覧に無い** | | | | | | |\n", id.c_str());
            continue;
        }
        const krithingi::RawMesh raw =
            krithingi::load_kmesh("data/thingi10k/kmesh/" + id + ".kmesh");
        const krithingi::Quantized q =
            krithingi::quantize(raw, krithingi::make_transform(1000 + idx.at(id)));
        const mesh::TopologyReport t = mesh::check_topology(q.mesh.triangles);
        mesh::SplitStats st;
        mesh::SplitOptions sopt;
        sopt.diag_unresolved = true;
        sopt.verify_delta = true;
        // **radial sort には幾何が要ります**（平面の法線と頂点座標）
        std::vector<geom::PlaneD> nrm(q.mesh.triangles.size());
        for (std::size_t i = 0; i < q.mesh.triangles.size(); ++i) {
            const mesh::Tri& tr = q.mesh.triangles[i];
            nrm[i] = geom::plane_from_triangle(q.mesh.vertices[tr[0]], q.mesh.vertices[tr[1]],
                                               q.mesh.vertices[tr[2]]);
        }
        std::vector<geom::HPointD> hv(q.mesh.vertices.size());
        for (std::size_t i = 0; i < q.mesh.vertices.size(); ++i) {
            hv[i] = geom::to_homogeneous(q.mesh.vertices[i]);
        }
        mesh::RadialGeom rg;
        rg.normal = &nrm;
        rg.vertices = &hv;
        std::vector<std::uint32_t> origin;
        const std::vector<mesh::Tri> out =
            mesh::split_contacts(q.mesh.triangles, q.mesh.vertices.size(), &origin, &st, nullptr,
                                 nullptr, nullptr, sopt, &rg);
        std::printf("| %s | %zu | %zu | %zu | %zu | %zu | %zu | **%zu** |\n", id.c_str(),
                    q.mesh.triangles.size(), t.edges_excess, st.radial_attempted,
                    st.radial_resolved, st.unsplit_edges, st.unresolved_post,
                    st.unresolved_detail.size());
        for (const auto& ue : st.unresolved_detail) {
            std::printf("|  | 辺 %u-%u 次数 %zu | %u: 三角形 %zu・過剰辺 %zu・扇 %zu→%zu |", ue.a,
                        ue.b, ue.degree, ue.a, ue.inc_a, ue.excess_a, ue.fans2_a, ue.fans_a);
            std::printf(" %u: 三角形 %zu・過剰辺 %zu・扇 %zu→%zu | | | |\n", ue.b, ue.inc_b,
                        ue.excess_b, ue.fans2_b, ue.fans_b);
        }
        (void)out;
    }
    return 0;
}

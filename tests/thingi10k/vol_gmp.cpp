// Krisite — 実データの小標本で体積恒等式を【厳密に】確かめる（SPEC-phase5 §3.0）
//
// **浮動小数点の篩は「大きな誤り」しか捕まえません**（格子 1 単位の四面体は
// 相対誤差 $10^{-18}$ で丸めの床より下）。**小さな誤りはここが受け持ちます。**
//
// **小標本に限ります。** 全件に掛けると篩の意味がなくなります。
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"

#include "thingi10k/loader.hpp"
#include "volume_gmp.hpp"

using namespace krisite;
using namespace krisite::csg;

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    const bool nsi_on = (std::getenv("KRI_NSI") != nullptr);
    std::printf("# GMP による体積恒等式（b=%d、NSI %s）\n\n", KRISITE_COORD_BITS,
                nsi_on ? "ON" : "OFF");
    std::printf("| 対 | ∪+∩ = A+B | A\\B = A-∩ | 位相 |\n|---|---|---|---|\n");
    int bad = 0, n = 0;
    for (int a = 1; a + 1 < argc; a += 2) {
        const auto qa =
            krithingi::quantize(krithingi::load_kmesh(argv[a]), krithingi::make_transform(1000));
        const auto qb = krithingi::quantize(krithingi::load_kmesh(argv[a + 1]),
                                            krithingi::make_transform(1001));
        PolySoup A = from_mesh(qa.mesh), B = from_mesh(qb.mesh);
        if (nsi_on) {
            A.nsi.assign(A.sources.size(), 1);
            B.nsi.assign(B.sources.size(), 1);
        }
        BoolOptions o;
        o.depth = 6;
        o.adaptive = true;
        o.cull_planes = true;
        o.early_out = true;
        o.cache_points = true;
        o.local_bsp = true;
        o.split_contacts = true;
        ToMeshOptions tm;
        tm.split_contacts = true;

        mpq_t va, vb, v[3], lhs, rhs;
        mpq_init(va);
        mpq_init(vb);
        mpq_init(lhs);
        mpq_init(rhs);
        for (auto& q : v) mpq_init(q);
        kritest::input_volume6(va, qa.mesh);
        kritest::input_volume6(vb, qb.mesh);
        bool topo = true;
        int k = 0;
        for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
            const SoupMesh m = to_mesh(boolean(A, B, op, o), tm);
            const mesh::TopologyReport t = mesh::check_topology(m.triangles);
            if (!t.empty && !(t.edge_manifold && t.vertex_manifold && t.oriented)) topo = false;
            kritest::mesh_volume6(v[k++], m);
        }
        mpq_add(lhs, v[0], v[1]);
        mpq_add(rhs, va, vb);
        const bool id_ok = mpq_equal(lhs, rhs) != 0;
        mpq_sub(lhs, va, v[1]);
        const bool df_ok = mpq_equal(lhs, v[2]) != 0;
        ++n;
        if (!id_ok || !df_ok) ++bad;
        std::printf("| %s | %s | %s | %s |\n", argv[a] + 22, id_ok ? "**一致**" : "**不一致**",
                    df_ok ? "**一致**" : "**不一致**", topo ? "ok" : "**NG**");
        mpq_clear(va);
        mpq_clear(vb);
        mpq_clear(lhs);
        mpq_clear(rhs);
        for (auto& q : v) mpq_clear(q);
    }
    std::printf("\n**%d 対中 %d 対で恒等式が破れました**\n", n, bad);
    return bad == 0 ? 0 : 1;
}

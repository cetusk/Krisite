// ---------------------------------------------------------------------------
// O3（切る三角形が触れない断片は切らない）の条件 1 を突く構成（SPEC-phase3 §5.4、SPEC-phase5
// §5.10.14.28）
//
//   P1, P2  同じ支持平面（z = 0）に載り、一部が重なる（A の底面 x∈[0,40]、B の底面 x∈[20,60]）
//   T       P1 の重ならない側（x∈[6,14]）にだけ箱が触れる三角形（B の四面体の面。平面 x + 4z = 30）
//
// 平面 x + 4z = 30 は z = 0 で x = 30（重なりの中）を通るので、
//   O3 なし: P1 も P2 もこの平面で切られ、重なり [20,40] の断片の頂点集合は一致する
//   O3 あり: P2 は T の箱に触れないので切られず、P1 だけ切られる → 頂点集合が食い違う可能性
// これを旗の ON / OFF で確かめる。壊れていれば位相か体積に出る。
// ---------------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"

#include "corpus.hpp"
#include "test_util.hpp"
#include "volume_fp.hpp"

using namespace krisite;

namespace {

unsigned long long hash_mesh(const csg::SoupMesh& m) {
    unsigned long long h = 1469598103934665603ull;
    const auto mix = [&h](unsigned long long v) {
        h ^= v;
        h *= 1099511628211ull;
    };
    for (const auto& v : m.vertices) {
        for (int k = 0; k < 4; ++k) {
            const auto& c = k == 0 ? v.x : k == 1 ? v.y : k == 2 ? v.z : v.w;
            unsigned char buf[sizeof(c)];
            std::memcpy(buf, &c, sizeof(c));
            for (std::size_t l = 0; l < sizeof(c); ++l) mix(buf[l]);
        }
    }
    for (const auto& t : m.triangles) {
        mix(t[0]);
        mix(t[1]);
        mix(t[2]);
    }
    return h;
}

void append(mesh::TriMesh& dst, const mesh::TriMesh& src) {
    const auto off = static_cast<std::uint32_t>(dst.vertices.size());
    dst.vertices.insert(dst.vertices.end(), src.vertices.begin(), src.vertices.end());
    for (const mesh::Tri& t : src.triangles)
        dst.triangles.push_back({t[0] + off, t[1] + off, t[2] + off});
}

/// 外向きの四面体（符号付き体積が正になるように面の向きを決める）
mesh::TriMesh tetra(geom::IPoint a, geom::IPoint b, geom::IPoint c, geom::IPoint d) {
    mesh::TriMesh m;
    m.vertices = {a, b, c, d};
    const auto vol = [&](const geom::IPoint& p, const geom::IPoint& q, const geom::IPoint& r,
                         const geom::IPoint& s) {
        const long double qx = q.x - p.x, qy = q.y - p.y, qz = q.z - p.z;
        const long double rx = r.x - p.x, ry = r.y - p.y, rz = r.z - p.z;
        const long double sx = s.x - p.x, sy = s.y - p.y, sz = s.z - p.z;
        return qx * (ry * sz - rz * sy) - qy * (rx * sz - rz * sx) + qz * (rx * sy - ry * sx);
    };
    if (vol(a, b, c, d) > 0) {
        m.triangles = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}};
    } else {
        m.triangles = {{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {0, 2, 3}};
    }
    return m;
}

}  // namespace

int main() {
    // 単位 u = 2^(b-7)。座標はすべて正なので、深度 1 までは 1 つの葉に全部入る
    const std::int32_t u = kritest::at(1, 64);
    const auto P = [u](int x, int y, int z) { return geom::IPoint{x * u, y * u, z * u}; };
    const mesh::TriMesh A = kritest::box(0, 0, 0, 40 * u, 20 * u, 20 * u);
    mesh::TriMesh B = kritest::box(20 * u, 0, 0, 60 * u, 20 * u, 20 * u);
    // **T は支持平面 z = 0 を跨ぐこと**（`cuts_for`
    // は支持平面に触れる三角形の平面だけを切断集合に入れる）。 面 (6,19,0),(12,9,2),(10,17,-1)
    // は平面 x + y + 2z = 25 に載り、z = 0 との交線 x + y = 25 は P1 の重ならない側（x <
    // 20）と重なり（x∈[20,25]）の両方を通る。箱 x∈[6,12] は P2（x ≥ 20）に触れない。
    append(B, tetra(P(6, 19, 0), P(12, 9, 2), P(10, 17, -1), P(9, 14, -3)));

    std::printf(
        "| 深度 | 演算 | O3 | 領域 | 出力多角形 | 三角形 | 位相 | 6 倍体積 | ハッシュ | "
        "飛ばした（箱 / 厳密） | 実際に分けた "
        "|\n|---|---|---|---:|---:|---:|---|---:|---|---|---:|\n");
    int broken = 0, idle = 0;
    for (unsigned d = 0; d <= 3; ++d) {
        for (csg::BoolOp op :
             {csg::BoolOp::Union, csg::BoolOp::Intersection, csg::BoolOp::Difference}) {
            unsigned long long h[3] = {0, 0, 0};
            double vol[3] = {0, 0, 0};
            bool ok[3] = {false, false, false};
            for (int mode : {0, 1, 2}) {
                csg::BoolOptions o = kritest::corpus_options(d);
                o.bsp_skip_disjoint = mode;
                csg::BoolStats st;
                const csg::PolySoup s =
                    csg::boolean(csg::from_mesh(A), csg::from_mesh(B), op, o, &st);
                csg::ToMeshOptions tm;
                tm.split_contacts = true;
                const csg::SoupMesh m = csg::to_mesh(s, tm);
                const mesh::TopologyReport r = mesh::check_topology(m.triangles);
                h[mode] = hash_mesh(m);
                vol[mode] = kritest::volume6_fp(m);
                ok[mode] = r.ok();
                std::printf(
                    "| %u | %s | %d | %zu | %zu | %zu | %s | %.6g | `%016llx` | %zu / %zu | %zu "
                    "|\n",
                    d,
                    op == csg::BoolOp::Union          ? "∪"
                    : op == csg::BoolOp::Intersection ? "∩"
                                                      : "＼",
                    mode, st.regions, s.polys.size(), m.triangles.size(),
                    r.ok() ? "ok" : "**壊れ**", vol[mode], h[mode], st.bsp_skip_box,
                    st.bsp_skip_exact, st.bsp_split_actual);
                if (mode != 0 && st.bsp_skip_box + st.bsp_skip_exact == 0) ++idle;
            }
            for (int mode : {1, 2}) {
                const double sc = std::max(std::abs(vol[0]), 1.0);
                const bool same_vol = std::abs(vol[mode] - vol[0]) / sc < 1e-9;
                if (!ok[mode] || !same_vol) ++broken;
            }
        }
    }
    std::printf(
        "\n**壊れた（位相が崩れた、または体積が従来と違う）組: %d**（深度 4 × 演算 3 × 旗 2 = 24 "
        "組のうち）\n",
        broken);
    std::printf("**機構が発火しなかった組（飛ばした数 0）: %d**（空回りなら構成の誤り）\n", idle);
    // **空回りは失敗**（機構が発火していなければ、この構成は何も検査していない）
    return (broken == 0 && idle == 0) ? 0 : 1;
}

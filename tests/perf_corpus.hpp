// Krisite — 性能測定用のコーパス（`SPEC-phase4.md` §8.2）
//
// **正しさのコーパスとは目的が違います。混ぜないでください。**
//
// | | 狙い | 規模 |
// |---|---|---|
// | `corpus.hpp` | **退化**（共平面、共線、接触、非対称） | 小さい（全 44 メッシュで 488 三角形） |
// | ここ | **規模**（スケーリングの測定） | 数万〜十数万三角形 |
//
// **既存のコーパスではスケーリングが測れません。** 1 回のブール演算あたりの仕事が
// 1 ms 程度しかなく、**スレッドのディスパッチが並列区間の 1/3 を占めます**
// （実測。`IMPL-phase4.md` §1.2）。
//
// **期待値は持ちません**（§8.2）。既存の恒等式検査が通れば十分です。
#ifndef KRISITE_TESTS_PERF_CORPUS_HPP
#define KRISITE_TESTS_PERF_CORPUS_HPP

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "krisite/config.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace kriperf {

using krisite::geom::IPoint;
using krisite::mesh::Tri;
using krisite::mesh::TriMesh;

/// 緯度経度で切った球（閉多様体）。
///
/// **極は扇、中間は四角形を 2 分割**します。半径を大きく取れば量子化で退化しません
/// （最短辺は $r \sin(\pi/n_{lat})$）。
inline TriMesh sphere(std::int64_t r, int nlat, int nlon, std::int64_t cx, std::int64_t cy,
                      std::int64_t cz) {
    TriMesh m;
    const double pi = 3.14159265358979323846;
    const auto q = [](double v) { return static_cast<std::int32_t>(std::llround(v)); };
    m.vertices.push_back({q(cx), q(cy), q(cz + r)});  // 北極 = 0
    for (int i = 1; i < nlat; ++i) {
        const double phi = pi * static_cast<double>(i) / nlat;
        for (int j = 0; j < nlon; ++j) {
            const double th = 2.0 * pi * static_cast<double>(j) / nlon;
            m.vertices.push_back({q(cx + r * std::sin(phi) * std::cos(th)),
                                  q(cy + r * std::sin(phi) * std::sin(th)),
                                  q(cz + r * std::cos(phi))});
        }
    }
    m.vertices.push_back({q(cx), q(cy), q(cz - r)});  // 南極 = 最後
    const auto ring = [nlon](int i, int j) {
        return static_cast<std::uint32_t>(1 + (i - 1) * nlon + (j % nlon));
    };
    const auto south = static_cast<std::uint32_t>(m.vertices.size() - 1);
    for (int j = 0; j < nlon; ++j) m.triangles.push_back({0u, ring(1, j), ring(1, j + 1)});
    for (int i = 1; i + 1 < nlat; ++i) {
        for (int j = 0; j < nlon; ++j) {
            m.triangles.push_back({ring(i, j), ring(i + 1, j), ring(i, j + 1)});
            m.triangles.push_back({ring(i, j + 1), ring(i + 1, j), ring(i + 1, j + 1)});
        }
    }
    for (int j = 0; j < nlon; ++j)
        m.triangles.push_back({south, ring(nlat - 1, j + 1), ring(nlat - 1, j)});
    return m;
}

/// 軸平行な直方体（外向き法線）。
inline TriMesh box(std::int64_t lo0, std::int64_t lo1, std::int64_t lo2, std::int64_t hi0,
                   std::int64_t hi1, std::int64_t hi2) {
    const auto c = [](std::int64_t v) { return static_cast<std::int32_t>(v); };
    TriMesh m;
    m.vertices = {{c(lo0), c(lo1), c(lo2)}, {c(hi0), c(lo1), c(lo2)}, {c(hi0), c(hi1), c(lo2)},
                  {c(lo0), c(hi1), c(lo2)}, {c(lo0), c(lo1), c(hi2)}, {c(hi0), c(lo1), c(hi2)},
                  {c(hi0), c(hi1), c(hi2)}, {c(lo0), c(hi1), c(hi2)}};
    const std::uint32_t f[12][3] = {{0, 3, 2}, {0, 2, 1}, {4, 5, 6}, {4, 6, 7},
                                    {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                                    {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    for (const auto& t : f) m.triangles.push_back({t[0], t[1], t[2]});
    return m;
}

inline TriMesh concat(const TriMesh& a, const TriMesh& b) {
    TriMesh r = a;
    const auto off = static_cast<std::uint32_t>(a.vertices.size());
    r.vertices.insert(r.vertices.end(), b.vertices.begin(), b.vertices.end());
    for (const Tri& t : b.triangles) r.triangles.push_back({t[0] + off, t[1] + off, t[2] + off});
    return r;
}

/// $n^3$ 個の立方体の格子（1 つのメッシュ）。**PWN です**（各立方体が閉じている）。
inline TriMesh cube_grid(int n, std::int64_t pitch, std::int64_t half, std::int64_t origin) {
    TriMesh m;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                const std::int64_t x = origin + pitch * i, y = origin + pitch * j,
                                   z = origin + pitch * k;
                m = concat(m, box(x - half, y - half, z - half, x + half, y + half, z + half));
            }
        }
    }
    return m;
}

struct PerfCase {
    std::string id;
    std::string what;
    TriMesh a, b;
};

/// **一様な分布**と**偏った分布**の 2 件（`SPEC-phase4.md` §8.2）。
///
/// `scale` で規模を変えられます。**測定は規模を変えて 2 点以上取ること** —
/// 1 点だけでは、飽和しているのか伸びているのかが分かりません。
inline std::vector<PerfCase> perf_corpus(int scale) {
    const std::int64_t r = std::int64_t{1} << (krisite::kCoordBits - 2);  // 2^(b-2)
    std::vector<PerfCase> v;
    v.push_back({"sphere2", "細分した球 x 2（一様）",
                 sphere(r, 8 * scale, 16 * scale, -r / 3, 0, 0),
                 sphere(r, 8 * scale, 16 * scale, r / 3, r / 7, 0)});
    const int n = 2 * scale;
    v.push_back({"grid", "立方体の格子 x 球（偏り）",
                 cube_grid(n, (2 * r) / n, r / (2 * n), -r + (2 * r) / (2 * n)),
                 sphere(r, 6 * scale, 12 * scale, 0, 0, 0)});
    return v;
}

}  // namespace kriperf

#endif  // KRISITE_TESTS_PERF_CORPUS_HPP

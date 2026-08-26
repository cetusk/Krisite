// Krisite — テストコーパスのメッシュ構築
//
// SPEC-phase1.md §9。ファイル入出力は非目標なので、テストコードで直接構築します。
//
// 向きの規約は「外から見て CCW」= 外向き法線（SPEC-phase1 §3.4）。
#ifndef KRISITE_TESTS_CORPUS_HPP
#define KRISITE_TESTS_CORPUS_HPP

#include <cstdint>
#include <vector>

#include "krisite/mesh/tri_mesh.hpp"

namespace kritest {

using krisite::geom::IPoint;
using krisite::mesh::Tri;
using krisite::mesh::TriMesh;

/// 軸平行な直方体。`lo` <= `hi` を各軸で満たすこと。外向き法線。
///
/// 頂点番号:
///   0=(lo,lo,lo) 1=(hi,lo,lo) 2=(hi,hi,lo) 3=(lo,hi,lo)
///   4=(lo,lo,hi) 5=(hi,lo,hi) 6=(hi,hi,hi) 7=(lo,hi,hi)
inline TriMesh box(std::int32_t lox, std::int32_t loy, std::int32_t loz, std::int32_t hix,
                   std::int32_t hiy, std::int32_t hiz) {
    TriMesh m;
    m.vertices = {
        {lox, loy, loz}, {hix, loy, loz}, {hix, hiy, loz}, {lox, hiy, loz},
        {lox, loy, hiz}, {hix, loy, hiz}, {hix, hiy, hiz}, {lox, hiy, hiz},
    };
    // 各面を外から見て CCW に。-Z 面は (0,3,2),(0,2,1) の向き
    m.triangles = {
        {0, 2, 1}, {0, 3, 2},  // -Z
        {4, 5, 6}, {4, 6, 7},  // +Z
        {0, 1, 5}, {0, 5, 4},  // -Y
        {1, 2, 6}, {1, 6, 5},  // +X
        {2, 3, 7}, {2, 7, 6},  // +Y
        {3, 0, 4}, {3, 4, 7},  // -X
    };
    return m;
}

/// 一辺 `s` の立方体を `(ox,oy,oz)` に置く。
inline TriMesh cube(std::int32_t ox, std::int32_t oy, std::int32_t oz, std::int32_t s) {
    return box(ox, oy, oz, ox + s, oy + s, oz + s);
}

/// 八分木は**座標範囲全体**を分割します（SPEC-phase1 §3.2）。
///
/// したがって入力が座標範囲の隅に小さく置かれていると、深度を上げてもセル境界が
/// 幾何をまたがず、**深度が効きません。** 実データでは量子化がモデルのバウンディング
/// ボックスを格子全体に写すので、テストの入力も座標範囲に見合った大きさにします。
///
/// b = 21 のとき深度 2 の内部境界は -2^19, 0, +2^19。
inline TriMesh grid_scale_box(std::int64_t lo, std::int64_t hi) {
    const auto c = [](std::int64_t v) { return static_cast<std::int32_t>(v); };
    return box(c(lo), c(lo), c(lo), c(hi), c(hi), c(hi));
}

/// 三角形の向きを全部裏返す（内向き法線のシェルを作る）。
inline TriMesh flipped(const TriMesh& m) {
    TriMesh r = m;
    for (Tri& t : r.triangles) {
        const auto tmp = t[1];
        t[1] = t[2];
        t[2] = tmp;
    }
    return r;
}

/// 2 つのメッシュを、頂点を共有せずに連結する（別シェルとして並べる）。
inline TriMesh concat(const TriMesh& a, const TriMesh& b) {
    TriMesh r = a;
    const auto off = static_cast<std::uint32_t>(a.vertices.size());
    r.vertices.insert(r.vertices.end(), b.vertices.begin(), b.vertices.end());
    for (const Tri& t : b.triangles) r.triangles.push_back({t[0] + off, t[1] + off, t[2] + off});
    return r;
}

/// 正四面体風（格子上なので正確な正四面体ではない）。外向き法線。
inline TriMesh tetra(std::int32_t s) {
    TriMesh m;
    m.vertices = {{0, 0, 0}, {s, 0, 0}, {0, s, 0}, {0, 0, s}};
    m.triangles = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    return m;
}

/// **組合せ的な**トーラス（n x m のグリッドを周期境界で貼る）。
///
/// `check_topology` は幾何を見ないので、種数の計算を検証する目的なら
/// 頂点座標は何でも構いません。V = nm, E = 3nm, F = 2nm → χ = 0, g = 1。
inline TriMesh combinatorial_torus(std::uint32_t n, std::uint32_t m) {
    TriMesh r;
    auto vid = [n, m](std::uint32_t i, std::uint32_t j) { return (i % n) * m + (j % m); };
    for (std::uint32_t i = 0; i < n; ++i) {
        for (std::uint32_t j = 0; j < m; ++j) {
            r.vertices.push_back(IPoint{static_cast<std::int32_t>(i), static_cast<std::int32_t>(j),
                                        static_cast<std::int32_t>((i * 7 + j * 13) % 11)});
        }
    }
    for (std::uint32_t i = 0; i < n; ++i) {
        for (std::uint32_t j = 0; j < m; ++j) {
            const auto a = vid(i, j), b = vid(i + 1, j), c = vid(i + 1, j + 1), d = vid(i, j + 1);
            r.triangles.push_back({a, b, c});
            r.triangles.push_back({a, c, d});
        }
    }
    return r;
}

/// 2 つの立方体が 1 頂点だけを共有する構成（SPEC-phase1 §9.2 のケース 11b）。
///
/// **和集合は非多様体になります。** 共有頂点まわりに扇が 2 つできるので、
/// `check_topology` の頂点多様体検査が落ちるはずです。
inline TriMesh two_cubes_sharing_a_vertex(std::int32_t s) {
    const TriMesh a = cube(0, 0, 0, s);  // 頂点 6 = (s,s,s)
    const TriMesh b = cube(s, s, s, s);  // 頂点 0 = (s,s,s)
    TriMesh r = concat(a, b);
    // b 側の頂点 0（連結後の添字 8）を a 側の頂点 6 に同一視する
    const std::uint32_t dup = static_cast<std::uint32_t>(a.vertices.size()) + 0;
    for (Tri& t : r.triangles) {
        for (auto& v : t) {
            if (v == dup) v = 6;
        }
    }
    return r;
}

}  // namespace kritest

#endif  // KRISITE_TESTS_CORPUS_HPP

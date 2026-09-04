// Krisite — テストコーパスのメッシュ構築
//
// SPEC-phase1.md §9。ファイル入出力は非目標なので、テストコードで直接構築します。
//
// 向きの規約は「外から見て CCW」= 外向き法線（SPEC-phase1 §3.4）。
#ifndef KRISITE_TESTS_CORPUS_HPP
#define KRISITE_TESTS_CORPUS_HPP

#include <array>
#include <cstdint>
#include <map>
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

/// **座標を絶対値で書かないこと**（SPEC-phase1 §9.0）。
///
/// `b` は CMake オプションで変わり、CI は b=21 と b=26 の両方を回します。絶対値で
/// 書くと b=26 側で入力が座標範囲の 1/32 になり、深度掃引が空回りします。
/// §9.0 の罠が **b 依存で再発する**わけです。
///
/// そこで位置は $2^{b-1}$ を 1 とみなした**比**で書きます。分母は 2 の冪に限ること。
/// そうすれば b が変わっても格子点のまま、セル境界との相対関係も保たれます。
///
/// b = 21 なら `at(-1, 1)` = -2^20 = `kCoordMin`、`at(1, 2)` = 2^19。
constexpr std::int32_t at(int num, int den) noexcept {
    return static_cast<std::int32_t>((-static_cast<std::int64_t>(krisite::kCoordMin) * num) / den);
}

/// 深度 `d` のセル境界の間隔。`at()` で書いた値がこの倍数なら境界に乗ります。
constexpr std::int64_t cell_size_at(unsigned d) noexcept {
    return std::int64_t{1} << (krisite::kCoordBits - d);
}

/// 比で指定する立方体。分子 2 つと共通の分母を渡す。
inline TriMesh ratio_box(int lo_num, int hi_num, int den) {
    return box(at(lo_num, den), at(lo_num, den), at(lo_num, den), at(hi_num, den), at(hi_num, den),
               at(hi_num, den));
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

/// 2 つの立方体が 1 **辺**だけを共有する構成（SPEC-phase1 §9.3、接触の次元 1）。
///
/// **和集合は非多様体になります。** 共有辺には 4 枚の面が接するので、
/// `check_topology` の辺多様体検査が落ちます。
///
/// **向きの整合は落ちません。** §9.3.0.1 の一般形（各辺で $\#(u,v)=\#(v,u)$）なら
/// 次数 4 の接触辺でも $2 = 2$ になるためです。
inline TriMesh two_cubes_sharing_an_edge(std::int32_t s) {
    const TriMesh a = cube(0, 0, 0, s);  // 頂点 2 = (s,s,0)、頂点 6 = (s,s,s)
    const TriMesh b = cube(s, s, 0, s);  // 頂点 0 = (s,s,0)、頂点 4 = (s,s,s)
    TriMesh r = concat(a, b);
    const auto off = static_cast<std::uint32_t>(a.vertices.size());
    for (Tri& t : r.triangles) {
        for (auto& v : t) {
            if (v == off + 0) v = 2;
            if (v == off + 4) v = 6;
        }
    }
    return r;
}

// ---- SPEC-phase1 §9.1 のケース ---------------------------------------------
//
// 追加するたびに `kCorpus` に足してください。§9.0 の番人（断片数の深度単調増加）と
// CP の掃引はこの表を回ります。**表に載せないケースは自動検査から漏れます。**

/// コーパスの 1 ケース。
struct Case {
    const char* id;    ///< §9.1 の番号
    const char* what;  ///< 何を突くか
    TriMesh (*make_a)();
    TriMesh (*make_b)();
    /// 両入力が軸平行な直方体か。**真なら解析的な期待体積を導けます**（§10.3）。
    ///
    /// 恒等式は自己整合の検査なので、系統的な誤りが相関して入ると素通りし得ます。
    /// 期待値との比較は**答えのレベルで独立**です。
    bool box_pair = false;
};

/// 軸平行境界箱。
struct Aabb64 {
    std::int64_t lo[3], hi[3];
};

inline Aabb64 aabb_of(const TriMesh& m) {
    Aabb64 r{{krisite::kCoordMax, krisite::kCoordMax, krisite::kCoordMax},
             {krisite::kCoordMin, krisite::kCoordMin, krisite::kCoordMin}};
    for (const IPoint& p : m.vertices) {
        const std::int64_t c[3] = {p.x, p.y, p.z};
        for (int t = 0; t < 3; ++t) {
            r.lo[t] = (c[t] < r.lo[t]) ? c[t] : r.lo[t];
            r.hi[t] = (c[t] > r.hi[t]) ? c[t] : r.hi[t];
        }
    }
    return r;
}

/// 2 つの箱の交差の各軸の長さ（交わらなければ 0）。
inline void overlap_extent(const Aabb64& a, const Aabb64& b, std::int64_t out[3]) {
    for (int t = 0; t < 3; ++t) {
        const std::int64_t lo = (a.lo[t] > b.lo[t]) ? a.lo[t] : b.lo[t];
        const std::int64_t hi = (a.hi[t] < b.hi[t]) ? a.hi[t] : b.hi[t];
        out[t] = (hi > lo) ? (hi - lo) : 0;
    }
}

namespace cases {

/// ケース 1: 一般位置。どの深度（0〜3）のセル境界にも面が乗らない。
///
/// 分子 9, 39, 27, 57（`kCoordMin` からの格子単位／$2^{b-6}$）はいずれも 8 の倍数で
/// ないので、深度 3 の境界間隔 $2^{b-3}$ = 8 単位に乗りません。
inline TriMesh case1_a() {
    return ratio_box(-23, 7, 32);
}
inline TriMesh case1_b() {
    return ratio_box(-5, 25, 32);
}

/// ケース 2: 面が完全共平面。$z = -2^{b-1}$ の面を共有し、**外向き法線は同じ向き**。
inline TriMesh case2_a() {
    return ratio_box(-1, 0, 1);
}
inline TriMesh case2_b() {
    return box(at(-1, 2), at(-1, 2), at(-1, 1), at(1, 2), at(1, 2), at(1, 2));
}

/// ケース 5: 面がセル境界と完全一致 ★ CP2 の主戦場。
///
/// A の面は $\{-2^{b-1}, 0\}$ = 深度 1 以降の境界、B の面は $\pm 2^{b-2}$ = 深度 2 以降の
/// 境界にちょうど乗ります。§5.2 の 4 平面同時交差と §4.2 の重複割り当てを同時に突きます。
inline TriMesh case5_a() {
    return ratio_box(-1, 0, 1);
}
inline TriMesh case5_b() {
    return ratio_box(-1, 1, 2);
}

/// ケース 6: ケース 5 を 1 格子だけずらした対照。
///
/// **ずらしただけで挙動が変わるなら危険**、を検出するための対です。
inline TriMesh case6_a() {
    return box(at(-1, 1) + 1, at(-1, 1) + 1, at(-1, 1) + 1, 1, 1, 1);
}
inline TriMesh case6_b() {
    return box(at(-1, 2) + 1, at(-1, 2) + 1, at(-1, 2) + 1, at(1, 2) + 1, at(1, 2) + 1,
               at(1, 2) + 1);
}

/// ケース 5T: セル角を**斜めに**通る面を持たせたケース 5 の派生。★
///
/// **軸平行な立方体だけのコーパスでは、1 点に集まる平面が 3 枚を超えられません。**
/// 1 軸あたり平面は高々 1 枚しか点を通れず（平行な 2 平面は交わらない）、
/// 軸は 3 本しかないためです。セル面が立方体の面と一致しても、両者は同一の幾何平面
/// なので `PlaneTable` で同じ ID に併合され、枚数は増えません。
///
/// つまり §5.2 の 4 平面同時交差は、**斜めの平面がなければ原理的に起きません。**
/// そこで四面体の斜面 $x+y+z=0$ を原点（深度 1 以降のセル角、かつ A の角）に通します。
/// 原点を通る平面は $x=0$, $y=0$, $z=0$, $x+y+z=0$ の 4 枚になります。
///
/// **この対は §5.4 の計測が 4 枚を検出できることの証拠でもあります。**
/// 検出できない計測器で「最大 3 枚」と報告しても何も言っていません。
inline TriMesh case5t_a() {
    return ratio_box(-1, 0, 1);
}
inline TriMesh case5t_b() {
    TriMesh m;
    // a, b, c は x+y+z = 0 上。d はその外
    m.vertices = {
        {at(1, 2), at(-1, 2), 0},
        {0, at(1, 2), at(-1, 2)},
        {at(-1, 2), 0, at(1, 2)},
        {at(-1, 2), at(-1, 2), at(-1, 2)},
    };
    m.triangles = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
    if (!krisite::mesh::is_outward_oriented(m)) m = flipped(m);
    return m;
}

/// 平面 $x+y+z=0$ 上に底面を持ち、負側に頂点がある四面体。
///
/// 底面の 3 点は成分和が 0 なのでこの平面上にあります。`(dx,dy,dz)` の和も 0 に
/// すること（平行移動しても底面が同じ平面に載るように）。
inline TriMesh slanted_tetra(std::int32_t dx, std::int32_t dy, std::int32_t dz) {
    TriMesh m;
    m.vertices = {
        {at(1, 2) + dx, at(-1, 2) + dy, 0 + dz},
        {0 + dx, at(1, 2) + dy, at(-1, 2) + dz},
        {at(-1, 2) + dx, 0 + dy, at(1, 2) + dz},
        {at(-1, 2) + dx, at(-1, 2) + dy, at(-1, 2) + dz},
    };
    m.triangles = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
    if (!krisite::mesh::is_outward_oriented(m)) m = flipped(m);
    return m;
}

/// ケース 2T: 共平面接触を**斜面**で（SPEC-phase1 §9.1）。
///
/// ケース 2 の斜面版です。2 つの四面体が平面 $x+y+z=0$ を共有し、底面が部分的に
/// 重なります。どちらの本体も同じ側にあるので**同方向**の共平面重複になります。
/// 平行移動は $(1,-1,0)\cdot 2^{b-3}$ で、成分和が 0 なので底面は同じ平面に留まります。
inline TriMesh case2t_a() {
    return slanted_tetra(0, 0, 0);
}
inline TriMesh case2t_b() {
    return slanted_tetra(at(1, 4), at(-1, 4), 0);
}

/// ケース 4T: 四面体 2 個が**頂点だけ**を共有（SPEC-phase1 §9.1）。**頂点に 6 平面**。
///
/// 原点に集まる面は各四面体につき 3 枚。**軸平行にしてはいけません。**
/// 軸平行だと 3 枚が $x{=}0, y{=}0, z{=}0$ になり、両者が `PlaneTable` で併合されて
/// 3 枚に潰れます（ケース 4 が 3 枚止まりなのと同じ理由）。そこで斜めに置きます。
///
/// **和集合は非多様体です**（ケース 11b と同型。§9.3）。合否には使いません。
inline TriMesh corner_tetra(int sgn) {
    const auto s = [sgn](int num, int den) {
        const std::int32_t v = at(num, den);
        // 正側は kCoordMax = 2^(b-1)-1 を超えないよう 1 だけ内側に寄せる
        return static_cast<std::int32_t>(sgn > 0 ? v - 1 : -v);
    };
    TriMesh m;
    m.vertices = {
        {0, 0, 0},
        {s(1, 1), s(1, 4), s(1, 8)},
        {s(1, 8), s(1, 1), s(1, 4)},
        {s(1, 4), s(1, 8), s(1, 1)},
    };
    m.triangles = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    if (!krisite::mesh::is_outward_oriented(m)) m = flipped(m);
    return m;
}
inline TriMesh case4t_a() {
    return corner_tetra(-1);
}
inline TriMesh case4t_b() {
    return corner_tetra(+1);
}

/// ケース 4T′: 頂点を共有し、**体積も重なる**四面体 2 個（SPEC-phase1 §9.1）。
///
/// **4T の $\cup$ は §9.3 で除外されるので、$\cup$ だけが高多重度の同時交差を
/// 通りません。** 4T′ がその穴を塞ぎます。
///
/// 4T は両者が別々の八分儀にあるため、共有頂点まわりに扇が 2 つできて非多様体に
/// なります。4T′ は両者を同じ側に置いて重ねるので、**扇が 1 つにつながり多様体に
/// なります。** 共有頂点に集まる平面は 6 枚のままです。
inline TriMesh cyclic_tetra(int n1, int n2, int n3, int den1, int den2, int den3) {
    TriMesh m;
    m.vertices = {
        {0, 0, 0},
        {at(n1, den1), at(n2, den2), at(n3, den3)},
        {at(n3, den3), at(n1, den1), at(n2, den2)},
        {at(n2, den2), at(n3, den3), at(n1, den1)},
    };
    m.triangles = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    if (!krisite::mesh::is_outward_oriented(m)) m = flipped(m);
    return m;
}
inline TriMesh case4tp_a() {
    return cyclic_tetra(-1, -1, -1, 1, 4, 8);
}
inline TriMesh case4tp_b() {
    return cyclic_tetra(-1, -1, -1, 2, 1, 16);
}

/// メッシュを平行移動する。
inline TriMesh translated(const TriMesh& m, std::int32_t dx, std::int32_t dy, std::int32_t dz) {
    TriMesh r = m;
    for (IPoint& p : r.vertices) {
        p = IPoint{p.x + dx, p.y + dy, p.z + dz};
    }
    return r;
}

/// ケース 3: 立方体 2 個、**辺が一致**（共線）。
///
/// $x=-2^{b-1}$, $y=-2^{b-1}$ の稜線を両者が共有し、体積も重なります。
/// 共線の退化を突きます。
inline TriMesh case3_a() {
    return ratio_box(-1, 0, 1);
}
inline TriMesh case3_b() {
    return box(at(-1, 1), at(-1, 1), at(-1, 2), at(1, 2), at(1, 2), at(1, 2));
}

/// ケース 4: 立方体 2 個、**頂点が一致**（点の一致）。
///
/// $(0, 0, -2^{b-1})$ が両者の頂点になります。**軸平行なので 3 枚止まりです**
/// （§9.1 の注記）。頂点が一致するなら、その 3 平面も一致して併合されるためです。
inline TriMesh case4_a() {
    return ratio_box(-1, 0, 1);
}
inline TriMesh case4_b() {
    return box(at(-1, 2), at(-1, 2), at(-1, 1), 0, 0, at(1, 2));
}

/// ケース 7: 格子に量子化した**回転立方体**。
///
/// $M = \begin{pmatrix} 1&2&2 \\ 2&1&-2 \\ -2&2&-1 \end{pmatrix}$ は
/// $M^{\mathsf T} M = 9I$、$\det M = 27$ の相似変換です。
///
/// **整数座標を整数座標に写すので、量子化による退化が入りません。** 丸めると辺が
/// わずかに非平面になり、面併合（§3.5）の前提が壊れます。この行列ならその心配が
/// ありません。得られる立方体は一辺 $3L$ で、**面はどれも軸に平行ではありません。**
///
/// $\det M > 0$ なので向きは保たれ、`box` の三角形リストをそのまま使えます。
inline TriMesh rotated_cube(std::int32_t l, std::int32_t ox, std::int32_t oy, std::int32_t oz) {
    static const std::int64_t kM[3][3] = {{1, 2, 2}, {2, 1, -2}, {-2, 2, -1}};
    TriMesh m = box(0, 0, 0, l, l, l);
    for (IPoint& p : m.vertices) {
        const std::int64_t v[3] = {p.x, p.y, p.z};
        std::int64_t r[3];
        for (int i = 0; i < 3; ++i) {
            r[i] = kM[i][0] * v[0] + kM[i][1] * v[1] + kM[i][2] * v[2];
        }
        p = IPoint{static_cast<std::int32_t>(r[0] + ox), static_cast<std::int32_t>(r[1] + oy),
                   static_cast<std::int32_t>(r[2] + oz)};
    }
    return m;
}
inline TriMesh case7_a() {
    // **相手に含まれないようずらすこと。** 含まれてしまうと A\B が空になり、
    // ケース 10 と同じものを測ることになります。
    return rotated_cube(at(1, 8), at(1, 8), 0, 0);
}
inline TriMesh case7_b() {
    return ratio_box(-1, 1, 2);
}

/// ケース 9: 四面体 2 個（一般位置）。鋭角。三角形数が最小。
///
/// 頂点も平面も共有しません。**斜面だけで構成されるので、`clip_edges` の
/// 三角形分岐（2.3.1 で直した経路）を最も濃く踏みます。**
inline TriMesh case9_a() {
    return cyclic_tetra(-1, -1, -1, 1, 4, 8);
}
inline TriMesh case9_b() {
    // **平行移動の向きに注意。** 元の四面体は $-2^{b-1}$ に接しているので、
    // 負方向へずらすと座標範囲を突き抜けます（`coords_in_range` が落ちます）。
    return translated(cyclic_tetra(-1, -1, -1, 2, 1, 16), at(1, 4), at(1, 8), at(1, 8));
}

/// ケース 10: 一方が他方を**完全に含む**。
///
/// **境界どうしの交差曲線が空**です。$A \setminus B$ は 2 シェル（外殻は外向き、
/// 内殻は内向き）になり、$\chi = 4$ になる唯一のケースです（§9.2）。
inline TriMesh case10_a() {
    return ratio_box(-1, 1, 2);
}
inline TriMesh case10_b() {
    // **B は A の内部に完全に含まれること。** 面を共有すると「完全に含む」ではなく
    // なり、$A \setminus B$ が 2 シェルにならないので §9.2 の期待値を検査できません。
    return ratio_box(-1, 1, 4);
}

/// ケース 11a: **面**接触のみ。
///
/// $x = 0$ で接するだけで体積は重なりません。正則化（§2.3）により
/// $\cap$ は空、$\cup$ は接触面が内部面になって消えます。**逆方向の共平面重複**です。
inline TriMesh case11a_a() {
    return ratio_box(-1, 0, 1);
}
inline TriMesh case11a_b() {
    return box(0, at(-1, 2), at(-1, 2), at(1, 2), at(1, 2), at(1, 2));
}

/// ケース 11b: **辺**接触のみ。
///
/// $x=0, y=0$ の稜線だけを共有します。**$\cup$ の出力は非多様体になる**ので
/// §9.3 で合否から除外します（辺まわりに 4 枚の面が接する）。
inline TriMesh case11b_a() {
    return ratio_box(-1, 0, 1);
}
inline TriMesh case11b_b() {
    return box(0, 0, at(-1, 1), at(1, 2), at(1, 2), 0);
}

/// ケース 12: セル境界を**またぐ大きな立方体**。継ぎ目の本数を増やします。
inline TriMesh case12_a() {
    return ratio_box(-3, 3, 4);
}
inline TriMesh case12_b() {
    // A からはみ出させること。含まれるとケース 10 と同じものを測ることになります。
    return box(at(-1, 4), at(-1, 4), at(-1, 4), at(1, 1) - 1, at(1, 1) - 1, at(1, 1) - 1);
}

/// ケース 8: 同一の立方体 2 個。全 6 平面を共有し、すべて同方向。
inline TriMesh case8() {
    return ratio_box(-1, 1, 2);
}

/// 補助: 交わらない 2 立方体（z 方向にだけ離す）。
///
/// x と y では原点をまたぐので、深度 1 でも分割が効きます。**両方を原点から離すと
/// 深度 1 で断片が増えず、§9.0 の番人に落ちます。**
inline TriMesh disjoint_a() {
    return box(at(-15, 16), at(-15, 16), at(-15, 16), at(5, 16), at(5, 16), at(-9, 16));
}
inline TriMesh disjoint_b() {
    return box(at(-5, 16), at(-5, 16), at(9, 16), at(15, 16), at(15, 16), at(15, 16));
}

/// **面が規則的に細分された立方体**（`SPEC-phase5.md` CP1.6）。
///
/// 1 面を $k \times k$ の格子に分けます。**併合すると 1 枚の正方形になりますが、
/// 各辺に中点が残ります。** 中点は角ではないので、その両側の辺は同じ直線上にあり、
/// **共線の頂点を落とさないと `support ∩ edge[i-1] ∩ edge[i]` が退化します。**
///
/// **CP1.6 で実際に落ちた形です**（48 三角形、深度 0）。`polysoup.hpp` には
/// `drop_collinear` がありましたが、`faces.hpp` に無く、
/// **二項メッシュ経路だけが落ちていました。**
///
/// **Thingi10K の 336 対では出ません。** あちらはスープ経路だからです。
/// **CAD 由来のメッシュでは、平坦面の規則的な分割は常態です。**
inline TriMesh tess_box(int lo_num, int hi_num, int den, int k) {
    TriMesh m;
    std::map<std::array<int, 3>, std::uint32_t> id;
    const auto coord = [&](int i) { return at(lo_num + (hi_num - lo_num) * i / k, den); };
    const auto vid = [&](int x, int y, int z) {
        const std::array<int, 3> key{x, y, z};
        const auto it = id.find(key);
        if (it != id.end()) return it->second;
        const std::uint32_t n = static_cast<std::uint32_t>(m.vertices.size());
        m.vertices.push_back(IPoint{coord(x), coord(y), coord(z)});
        id.emplace(key, n);
        return n;
    };
    // 面ごとに (原点, u 方向, v 方向)。u × v が外向きになる順に並べてあります
    const int F[6][9] = {
        {0, 0, 1, 1, 0, 0, 0, 1, 0},  // z = k
        {0, 0, 0, 0, 1, 0, 1, 0, 0},  // z = 0
        {1, 0, 0, 0, 1, 0, 0, 0, 1},  // x = k
        {0, 0, 0, 0, 0, 1, 0, 1, 0},  // x = 0
        {0, 1, 0, 0, 0, 1, 1, 0, 0},  // y = k
        {0, 0, 0, 1, 0, 0, 0, 0, 1},  // y = 0
    };
    for (const auto& f : F) {
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                const auto vat = [&](int a, int b) {
                    return vid(f[0] * k + f[3] * a + f[6] * b, f[1] * k + f[4] * a + f[7] * b,
                               f[2] * k + f[5] * a + f[8] * b);
                };
                const std::uint32_t p00 = vat(i, j), p10 = vat(i + 1, j), p11 = vat(i + 1, j + 1),
                                    p01 = vat(i, j + 1);
                m.triangles.push_back({p00, p10, p11});
                m.triangles.push_back({p00, p11, p01});
            }
        }
    }
    return m;
}

/// **面をでこぼこにした立方体**（`SPEC-phase5.md` (c)）。
///
/// **細分しただけでは平面が 6 枚のままです。** 局所 BSP の切断は平面ごとなので、
/// **三角形を増やしても切断は増えません。** 面の内部の格子点をその面の法線方向に
/// ずらすと、**三角形ごとに別の平面**になり、切断が意味を持ちます。
///
/// ずらす量は箱の 1.4% 以下（$7/512$ 対 $14/16$）なので、自己交差しません。
/// **比で書くので $b$ が変わっても格子点のまま**です。
inline TriMesh bumpy_box(int lo_num, int hi_num, int den, int k) {
    TriMesh m;
    std::map<std::array<int, 3>, std::uint32_t> id;
    const auto coord = [&](int i) { return lo_num + (hi_num - lo_num) * i / k; };
    const auto vid = [&](int x, int y, int z) {
        const std::array<int, 3> key{x, y, z};
        const auto it = id.find(key);
        if (it != id.end()) return it->second;
        // **面の内部の格子点だけ**、その面の外向き法線方向にずらす
        const int nb = (x == 0) + (x == k) + (y == 0) + (y == k) + (z == 0) + (z == k);
        std::int32_t d[3] = {0, 0, 0};
        if (nb == 1) {
            const unsigned h =
                static_cast<unsigned>(x * 73856093 ^ y * 19349663 ^ z * 83492791) % 7u + 1u;
            const std::int32_t off = at(static_cast<int>(h), den * 32);
            if (x == 0) {
                d[0] = -off;
            } else if (x == k) {
                d[0] = off;
            } else if (y == 0) {
                d[1] = -off;
            } else if (y == k) {
                d[1] = off;
            } else if (z == 0) {
                d[2] = -off;
            } else {
                d[2] = off;
            }
        }
        const std::uint32_t n = static_cast<std::uint32_t>(m.vertices.size());
        m.vertices.push_back(IPoint{static_cast<std::int32_t>(at(coord(x), den) + d[0]),
                                    static_cast<std::int32_t>(at(coord(y), den) + d[1]),
                                    static_cast<std::int32_t>(at(coord(z), den) + d[2])});
        id.emplace(key, n);
        return n;
    };
    const int F[6][9] = {
        {0, 0, 1, 1, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 1, 0, 1, 0, 0}, {1, 0, 0, 0, 1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 1, 0, 1, 0}, {0, 1, 0, 0, 0, 1, 1, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 1},
    };
    for (const auto& f : F) {
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                const auto vat = [&](int a, int b) {
                    return vid(f[0] * k + f[3] * a + f[6] * b, f[1] * k + f[4] * a + f[7] * b,
                               f[2] * k + f[5] * a + f[8] * b);
                };
                const std::uint32_t p00 = vat(i, j), p10 = vat(i + 1, j), p11 = vat(i + 1, j + 1),
                                    p01 = vat(i, j + 1);
                m.triangles.push_back({p00, p10, p11});
                m.triangles.push_back({p00, p11, p01});
            }
        }
    }
    return m;
}

/// ケース 26: **2 つの立体を離して置き、一部だけを重ねる**（`SPEC-phase5.md` (c)）。
///
/// **§9.0 (1) のサイズ規律の例外**（対の AABB は規律を満たしますが、
/// **個々の立体は座標範囲の一部しか占めません**）。規律の目的は
/// 「退化を偶然にしない」ことで、**配置の多様性を狭めることではありません。**
///
/// **狙いは単一 source のセルを作ること。** 実データでは 61% のセルがそうでしたが、
/// **コーパスは 1 つも作れていませんでした**（`IMPL-phase5.md` §20.4）。
/// NSI の最適化を CI で守るには、この配置が要ります。
inline TriMesh case26_a() {
    return bumpy_box(-15, -1, 16, 2);
}
inline TriMesh case26_b() {
    return bumpy_box(-3, 11, 16, 2);
}

/// ケース 26′: **26 の対照。** 立体の大きさと重なりの量を変えます。
/// **番人が 1 ケース依存だと、そのケースが消えた瞬間に強度が 0 になります。**
///
/// **細分は $k=2$ に留めること。** でこぼこの箱は三角形ごとに別の平面を持つので、
/// **NSI を切ったときの局所 BSP が $O(n_\text{葉}^2)$ で効きます。**
/// $k=3$（108 三角形）にしたら `interior` が 3.6 → 33.8 秒、`soup` が 17 → 78 秒に
/// なり、**Windows の CI で `parallel` の並列判定が負荷に負けました**
/// （`IMPL-phase5.md` §24）。**コーパスの費用は CI の予算です。**
inline TriMesh case26p_a() {
    return bumpy_box(-13, -1, 16, 2);
}
inline TriMesh case26p_b() {
    return bumpy_box(-5, 7, 16, 2);
}

/// **向きを逆にした立方体を足す**（`SPEC-phase3.md` §9 のケース 27）。
///
/// **面を反転すると法線が内向きになり、巻き数が負になります。**
/// **$\partial S = 0$ は保たれます** — 各辺の向きも一緒に反転するので、
/// 隣接する面との打ち消しが成立したままです。
inline void add_oriented_box(TriMesh& m, std::int32_t lo, std::int32_t hi, bool flip) {
    const auto o = static_cast<std::uint32_t>(m.vertices.size());
    for (int k = 0; k < 8; ++k)
        m.vertices.push_back({(k & 1) ? hi : lo, (k & 2) ? hi : lo, (k & 4) ? hi : lo});
    static const std::uint32_t t[12][3] = {{0, 3, 1}, {0, 2, 3}, {4, 5, 7}, {4, 7, 6},
                                           {0, 1, 5}, {0, 5, 4}, {1, 3, 7}, {1, 7, 5},
                                           {3, 2, 6}, {3, 6, 7}, {2, 0, 4}, {2, 4, 6}};
    for (const auto& x : t) {
        if (flip) {
            m.triangles.push_back({o + x[0], o + x[2], o + x[1]});
        } else {
            m.triangles.push_back({o + x[0], o + x[1], o + x[2]});
        }
    }
}

/// **ケース 27a: 逆向きの立方体 1 個**（12 三角形。`SPEC-phase3.md` §9.1）。
///
/// **すべての領域が $w = -1$ なので、正しくは「立体は無い」**です。
/// **旧定義（$w \ne 0$）では、裏返った立方体を「立体」として出力します。**
///
/// **既存の検査をすべて素通りします** — 多様体な出力が返り、体積の恒等式は
/// 符号が変わるだけで成立し、**正解器も同じ誤りを持っています。**
inline TriMesh case27a() {
    TriMesh m;
    add_oriented_box(m, at(-1, 2), 0, true);
    return m;
}

/// **ケース 27b: 大きい立方体 + 外に置いた逆向きの小さい立方体**（24 三角形）。
///
/// **一部だけが裏返る形**で、実データに近いと考えられます。
/// **かたまりの数で違いが出ます**（旧定義で 2、正しい定義で 1）。
inline TriMesh case27b() {
    TriMesh m;
    add_oriented_box(m, at(-1, 2), 0, false);
    add_oriented_box(m, at(1, 8), at(3, 8), true);
    return m;
}

/// **ケース 27c: 27b と同じ配置で、どちらも正しい向き**（対照）。
///
/// **巻き数が負にならないので、$w \ne 0$ と $w > 0$ で出力が完全に一致するはず**です。
/// **定義を変えたときに、自己交差しない入力の挙動が変わっていないことを守ります。**
inline TriMesh case27c() {
    TriMesh m;
    add_oriented_box(m, at(-1, 2), 0, false);
    add_oriented_box(m, at(1, 8), at(3, 8), false);
    return m;
}

/// ケース 25: 細分された立方体の対（$k=2$ で 48 三角形。**最小の再現**）。
inline TriMesh case25_a() {
    return tess_box(-12, 4, 16, 2);
}
inline TriMesh case25_b() {
    return tess_box(-4, 12, 16, 2);
}
/// ケース 25': 細分の粗さを変えても同じ立体になること（$k=4$）。
inline TriMesh case25p_a() {
    return tess_box(-12, 4, 16, 4);
}
inline TriMesh case25p_b() {
    return tess_box(-4, 12, 16, 4);
}

// ---- SPEC-phase2 §8 の追加ケース -------------------------------------------
//
// **座標の単位は $H/256$（$H = 2^{b-1}$）にそろえます。** 深度 0〜3 のセル境界は
// $-H + m \cdot H/4$、この単位で言えば $-256 + 64m$ です。追加ケースの座標は
// **どれも 64 の倍数（からのずれ）に乗せていません。** 乗せるとケース 5（面がセル境界と
// 一致）と同じ退化を同時に測ることになり、何を突いているのか分からなくなります。

/// ケース 13: **密度が偏る対**（§8）。**§2.4 の本丸。**
///
/// A は粗い大立方体、B は A の +X 面をまたぐ小さな立方体 4 個です。適応分割では
/// **B の近傍だけが深く分割され、それ以外は浅いまま**になるので、隣接セルの深さが
/// 2 段以上違う配置ができます。
///
/// **B が A の面をまたぐことが要点です。** 内部に浮かせると、B の三角形しか含まない
/// セルは §3.1 の判定（A と B の両方を含むか）で分割されず、深さの差ができません。
inline TriMesh case13_a() {
    return box(at(-193, 256), at(-193, 256), at(-193, 256), at(63, 256), at(63, 256), at(63, 256));
}
inline TriMesh case13_b() {
    TriMesh r;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            const int dy = -1 + 12 * i, dz = -1 + 12 * j;  // 一辺 8、隙間 4（互いに素）
            r = concat(r, box(at(59, 256), at(dy, 256), at(dz, 256), at(67, 256), at(dy + 8, 256),
                              at(dz + 8, 256)));
        }
    }
    return r;
}

/// ケース 14: 一方が極端に小さく、離れている（§8）。**early-out が支配的になる配置。**
///
/// B は深度 3 のセル 1 個の内部に収まる大きさで、A から離れています。
/// 相手が居ないセルが大半になるので、§3.2 の early-out の発火率がここで出ます。
inline TriMesh case14_a() {
    return box(at(-225, 256), at(-225, 256), at(-225, 256), at(31, 256), at(31, 256), at(31, 256));
}
inline TriMesh case14_b() {
    return box(at(101, 256), at(101, 256), at(101, 256), at(103, 256), at(103, 256), at(103, 256));
}

/// ケース 15: **細かい構造が斜面を持つ**（§8）。§2.4 × 4 平面同時交差。
///
/// ケース 13 の小立方体を四面体に替えたものです。斜面 $x+y+z=\mathrm{const}$ を持つので、
/// **軸平行では原理的に届かない 4 枚以上の同時交差**（`SPEC-phase1.md` §9.1）を
/// 深さの差と同時に突きます。
inline TriMesh case15_a() {
    return case13_a();
}
inline TriMesh case15_b() {
    // **1 個目の斜面をセル角に通します。** 四面体 `tetra(s)` を $(o_x,o_y,o_z)$ に置くと
    // 斜面は $x+y+z = s+o_x+o_y+o_z$。$(56,-1,-1)$ と $s=10$ で $x+y+z=64$ になり、
    // **セル角 $(64,0,0)$（$H/256$ 単位。64 は深度 3、0 は深度 1 の境界）を通ります。**
    // その点では斜面・$x{=}64$・$y{=}0$・$z{=}0$ の **4 平面が交わります**（§8 のケース 15）。
    //
    // **軸平行だけでは原理的に 3 枚を超えられません**（`SPEC-phase1.md` §9.1）。
    // ここを外すと、このケースは「斜面がある」だけで 4 平面同時交差を突けません。
    TriMesh r = translated(tetra(at(10, 256)), at(56, 256), at(-1, 256), at(-1, 256));
    // 2 個目は一般位置（$x+y+z=96$ はセル角を通らない）。対照として置きます
    return concat(r, translated(tetra(at(10, 256)), at(56, 256), at(15, 256), at(15, 256)));
}

/// 四角柱（断面は任意の凸四角形）。`px`/`py` は $+z$ から見て反時計回りに並べること。
/// 外から見て CCW（外向き法線）。
inline TriMesh prism4(const std::int32_t px[4], const std::int32_t py[4], std::int32_t zlo,
                      std::int32_t zhi) {
    TriMesh m;
    for (int i = 0; i < 4; ++i) m.vertices.push_back(IPoint{px[i], py[i], zlo});
    for (int i = 0; i < 4; ++i) m.vertices.push_back(IPoint{px[i], py[i], zhi});
    m.triangles.push_back({0, 2, 1});  // 底面（法線 -z）
    m.triangles.push_back({0, 3, 2});
    m.triangles.push_back({4, 5, 6});  // 上面（法線 +z）
    m.triangles.push_back({4, 6, 7});
    for (std::uint32_t i = 0; i < 4; ++i) {
        const std::uint32_t j = (i + 1) % 4;
        m.triangles.push_back({i, j, j + 4});
        m.triangles.push_back({i, j + 4, i + 4});
    }
    return m;
}

/// ケース 16: **自己接触**（SPEC-phase2 §8）。**§5.1.2 の対応付けの唯一の検証手段。**
///
/// 平板 $A$ から菱形プリズム $B$ を抜きます。$B$ の 4 頂点が $A$ の側面の中点に
/// **ちょうど載る**ので、$A \setminus B$ が 4 隅に分かれ、**隣り合う隅が縦の辺で接します。**
///
/// **13 / 14 / 15 は自己接触を作りません。** これらだけでは §5.1.2 の対応付けが
/// **owner でも連結成分でも同じ結果になり、素通りで通ります**（§8.1）。
/// $A \setminus B$ の各シートは $A$ の面 1 枚と $B$ の面 1 枚から成るので、
/// **owner で分けると組がシートをまたぎます。**
///
/// 菱形の辺は $x+y=\mathrm{const}$ なので**斜面**を含み、§5.2 も同時に突きます。
/// 中心を $(7,7)$ に置いてあるのは、深度 0〜3 のセル境界（$H/256$ 単位で $-256+64m$）に
/// 乗せないためです。
inline TriMesh case16_a() {
    return box(at(-121, 256), at(-121, 256), at(-9, 256), at(135, 256), at(135, 256), at(9, 256));
}
inline TriMesh case16_b() {
    // 菱形の頂点は A の側面の中点。$|x-7| + |y-7| = 128$
    const std::int32_t px[4] = {at(135, 256), at(7, 256), at(-121, 256), at(7, 256)};
    const std::int32_t py[4] = {at(7, 256), at(135, 256), at(7, 256), at(-121, 256)};
    return prism4(px, py, at(-130, 256), at(130, 256));
}

/// 断面が任意の**星型**単純多角形の角柱（頂点 0 からの扇で三角形化）。
inline TriMesh star_prism(const std::vector<std::int32_t>& px, const std::vector<std::int32_t>& py,
                          std::int32_t zlo, std::int32_t zhi) {
    const std::size_t n = px.size();
    TriMesh m;
    for (std::size_t i = 0; i < n; ++i) m.vertices.push_back({px[i], py[i], zlo});
    for (std::size_t i = 0; i < n; ++i) m.vertices.push_back({px[i], py[i], zhi});
    const auto u = [](std::size_t i) { return static_cast<std::uint32_t>(i); };
    const auto d = [n](std::size_t i) { return static_cast<std::uint32_t>(i + n); };
    for (std::size_t i = 1; i + 1 < n; ++i) {
        m.triangles.push_back({u(0), u(i + 1), u(i)});  // 底面（法線 -z）
        m.triangles.push_back({d(0), d(i), d(i + 1)});  // 上面（法線 +z）
    }
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        m.triangles.push_back({u(i), u(j), d(j)});
        m.triangles.push_back({u(i), d(j), d(i)});
    }
    return m;
}

/// ケース 21 の A: **L 字柱**（非凸の面）。
///
/// **凸分割が 1 枚に併合できない**面を持ちます。半空間の交わりは凸多角形しか
/// 表せないので、凸性の判定が壊れると**凸包**になり、点集合が変わります。
inline TriMesh case21_a() {
    const std::vector<std::int32_t> px = {at(-121, 256), at(135, 256), at(135, 256),
                                          at(7, 256),    at(7, 256),   at(-121, 256)};
    const std::vector<std::int32_t> py = {at(-121, 256), at(-121, 256), at(7, 256),
                                          at(7, 256),    at(135, 256),  at(135, 256)};
    return star_prism(px, py, at(-128, 256), at(128, 256));
}

/// ケース 21 の B: **穴あき平板**（外枠と内枠を側面でつないだ角柱）。
///
/// 上下面が**穴あきの領域**なので、凸分割の併合が穴を塞いではいけません。
/// 外枠を反時計回り、内枠を時計回りに置き、4 本の橋で 1 つの帯にしています。
inline TriMesh case21_b() {
    TriMesh m;
    const std::int32_t o[4][2] = {{at(-135, 256), at(-135, 256)},
                                  {at(121, 256), at(-135, 256)},
                                  {at(121, 256), at(121, 256)},
                                  {at(-135, 256), at(121, 256)}};
    const std::int32_t i4[4][2] = {{at(-57, 256), at(-57, 256)},
                                   {at(43, 256), at(-57, 256)},
                                   {at(43, 256), at(43, 256)},
                                   {at(-57, 256), at(43, 256)}};
    const std::int32_t zlo = at(-40, 256), zhi = at(40, 256);
    // 頂点: 外枠 lo(0-3) / 内枠 lo(4-7) / 外枠 hi(8-11) / 内枠 hi(12-15)
    for (const auto& p : o) m.vertices.push_back({p[0], p[1], zlo});
    for (const auto& p : i4) m.vertices.push_back({p[0], p[1], zlo});
    for (const auto& p : o) m.vertices.push_back({p[0], p[1], zhi});
    for (const auto& p : i4) m.vertices.push_back({p[0], p[1], zhi});
    const auto ol = [](int k) { return static_cast<std::uint32_t>(k % 4); };
    const auto il = [](int k) { return static_cast<std::uint32_t>(4 + k % 4); };
    const auto oh = [](int k) { return static_cast<std::uint32_t>(8 + k % 4); };
    const auto ih = [](int k) { return static_cast<std::uint32_t>(12 + k % 4); };
    for (int k = 0; k < 4; ++k) {
        // 底面（法線 -z）: 外枠と内枠のあいだの帯
        m.triangles.push_back({ol(k), il(k), ol(k + 1)});
        m.triangles.push_back({ol(k + 1), il(k), il(k + 1)});
        // 上面（法線 +z）
        m.triangles.push_back({oh(k), oh(k + 1), ih(k)});
        m.triangles.push_back({oh(k + 1), ih(k + 1), ih(k)});
        // 外側面
        m.triangles.push_back({ol(k), ol(k + 1), oh(k + 1)});
        m.triangles.push_back({ol(k), oh(k + 1), oh(k)});
        // 内側面（法線は内向きの穴を向く = 立体から見て外向き）
        m.triangles.push_back({il(k), ih(k), ih(k + 1)});
        m.triangles.push_back({il(k), ih(k + 1), il(k + 1)});
    }
    return m;
}

/// ケース 24: **連結成分では対応付けられない自己接触**（`SPEC-phase2.md` §5.1.2.2）。
///
/// **Phase 5 の CP1 が実データで踏んだ配置を最小化したもの**（6,292 三角形 → 36）。
///
/// 次数 4 の辺に集まる 4 枚が**同じ連結成分**に属すると、面の連結成分では組に
/// 分けられません。§5.1.2 の基準が破れる**唯一の**配置です。
///
/// ```
///   A = 辺 (x=0, y=0) だけで接する 2 つの箱      ← ここが次数 4 の辺になる
///   B = 上から両方に被さる箱                      ← これが 1 成分に繋ぐ
/// ```
///
/// **橋が無ければ通ります。** 対照は下の 3 つで、落ちるのは C = 1 のときだけです。
///
/// | 構成 | 出力の成分数 | 対応付け |
/// |---|---:|---|
/// | 11b（別々の 2 立体が辺接触） | 2 | できる |
/// | 24 の A × 離れた小箱 | 3 | できる |
/// | **24（橋あり）** | **1** | **できない** |
///
/// **実データでは 306 対のうち 3 対で発火しました。** コーパスでは 0 件だった経路です。
inline TriMesh case24_a() {
    return concat(box(at(-1, 2), at(-1, 2), at(-1, 2), 0, 0, at(1, 2)),
                  box(0, 0, at(-1, 2), at(1, 2), at(1, 2), at(1, 2)));
}
inline TriMesh case24_b() {
    return box(at(-1, 2), at(-1, 2), at(1, 4), at(1, 2), at(1, 2), at(1, 2));
}

/// ケース 24′: **24 の対照。** 橋を渡さないので成分が分かれ、対応付けできます。
///
/// **2 つの基準が一致する配置だけでテストしない**（`CLAUDE.md`）の実践です。
/// 24 だけを入れると「落ちること」は分かっても「何が原因か」が分かりません。
inline TriMesh case24p_b() {
    return box(at(-1, 2), at(-1, 2), at(5, 8), at(-1, 4), at(-1, 4), at(3, 4));
}

}  // namespace cases

/// §9.1 の実装済みケース。CP3 に向けて増やしていきます。
inline const std::vector<Case>& corpus() {
    static const std::vector<Case> k = {
        {"1", "一般位置", cases::case1_a, cases::case1_b, true},
        {"2", "面が完全共平面", cases::case2_a, cases::case2_b, true},
        {"3", "辺が一致（共線）", cases::case3_a, cases::case3_b, true},
        {"4", "頂点が一致", cases::case4_a, cases::case4_b, true},
        {"5", "面がセル境界と一致", cases::case5_a, cases::case5_b, true},
        {"7", "量子化した回転立方体", cases::case7_a, cases::case7_b},
        {"9", "四面体 2 個（一般位置）", cases::case9_a, cases::case9_b},
        {"10", "一方が他方を完全に含む", cases::case10_a, cases::case10_b, true},
        {"11a", "面接触のみ", cases::case11a_a, cases::case11a_b, true},
        {"11b", "辺接触のみ", cases::case11b_a, cases::case11b_b, true},
        {"12", "セル境界をまたぐ大立方体", cases::case12_a, cases::case12_b, true},
        {"2T", "共平面接触を斜面で", cases::case2t_a, cases::case2t_b},
        {"4T", "四面体2個が頂点を共有", cases::case4t_a, cases::case4t_b},
        {"4T'", "頂点共有 + 体積も重なる", cases::case4tp_a, cases::case4tp_b},
        {"5T", "セル角を斜面が通る", cases::case5t_a, cases::case5t_b},
        {"6", "セル境界から 1 格子ずれ", cases::case6_a, cases::case6_b, true},
        {"8", "同一の立方体", cases::case8, cases::case8, true},
        {"D", "交わらない 2 立方体", cases::disjoint_a, cases::disjoint_b, true},
        // SPEC-phase2 §8。**13 が §2.4 の本丸。** B が複数の箱・四面体なので
        // `box_pair`（単一の AABB 対を前提にした解析体積）は使えません
        {"13", "密度が偏る対", cases::case13_a, cases::case13_b},
        {"14", "極端に小さく離れている", cases::case14_a, cases::case14_b, true},
        {"15", "細かい構造が斜面を持つ", cases::case15_a, cases::case15_b},
        // **§5.1.2 の対応付けはこのケースでしか検証できません**（§8.1）
        {"16", "自己接触（平板 − 菱形プリズム）", cases::case16_a, cases::case16_b},
        // **Phase 5 の CP1 が実データで踏んだ配置**（`SPEC-phase2.md` §5.1.2.2）。
        // 次数 4 の辺の 4 枚が同じ連結成分に入るので、成分では分けられません
        {"24", "連結成分で分けられない自己接触", cases::case24_a, cases::case24_b},
        {"24'", "24 の対照（橋なし。成分が分かれる）", cases::case24_a, cases::case24p_b},
        // **Phase 5 の CP1.6 が踏んだ配置**（`SPEC-phase5.md` CP1.6）。
        // 面が規則的に細分されており、併合すると辺の中点が共線頂点として残ります。
        // **落とさないと `intersect3` が退化します**（深度 0 で落ちていました）
        {"25", "面が 2x2 に細分された立方体（共線頂点）", cases::case25_a, cases::case25_b, true},
        {"25'", "25 を 4x4 に細分（粗さを変えても同じ立体）", cases::case25p_a, cases::case25p_b,
         true},
        // **Phase 5 の (c)**。離して置き、一部だけ重ねる。**単一 source のセルを作る**
        // 唯一のケースで、NSI の最適化（`SPEC-phase3.md` §5.6）を CI で守れます。
        // **§9.0 (1) のサイズ規律の例外**（対の AABB は満たす）
        {"26", "離して置き一部だけ重ねる（単一 source のセル）", cases::case26_a, cases::case26_b},
        // **26 の対照（規模と重なり方を変える）。** 番人が 1 ケース依存だと、
        // そのケースが消えた瞬間に強度が 0 になります
        {"26'", "26 の対照（大きさと重なりが違う）", cases::case26p_a, cases::case26p_b},
    };
    return k;
}

// ---- $n$ 項のケース（`SPEC-phase3.md` §9 のケース 17 / 18）--------------------
//
// **二項の連鎖とは突き合わせられません。** 二項の `boolean_op` は `BoolMesh` を返し、
// 構成点が有理数なので**入力に戻せない**からです（それが $n$ 項にした理由でもあります）。
//
// そこで**別の演算木で同じ立体を作る**のを正解器にします。集合の恒等式なので、
// **答えのレベルで独立**です（自己整合の検査ではありません）。

namespace cases {

/// ケース 17: $n=3$。**共平面で接する 2 箱を、斜めの角柱が貫く。**
///
/// A と B は $x$ で重なり、$y$ / $z$ の 4 面が**完全に共平面**です（§9.1 の退化）。
/// C は菱形断面の角柱で、**斜面**を持ちます（軸平行だけでは 4 平面同時交差を
/// 突けません。`SPEC-phase1.md` §9.1）。
inline std::vector<TriMesh> case17() {
    const TriMesh a =
        box(at(-121, 256), at(-121, 256), at(-121, 256), at(7, 256), at(135, 256), at(135, 256));
    const TriMesh b =
        box(at(-7, 256), at(-121, 256), at(-121, 256), at(135, 256), at(135, 256), at(135, 256));
    // 菱形 $|x-7| + |y-7| = 64$。**中心を原点からずらす**（対称性はバグを隠します）
    const std::int32_t px[4] = {at(71, 256), at(7, 256), at(-57, 256), at(7, 256)};
    const std::int32_t py[4] = {at(7, 256), at(71, 256), at(7, 256), at(-57, 256)};
    const TriMesh c = prism4(px, py, at(-140, 256), at(140, 256));
    return {a, b, c};
}

/// ケース 18: **同一形状を 5 回引く**（$W - T - T - T - T - T$）。
///
/// **内外の 1 ビットでは表せません。** 同じ曲面を何度も跨ぐので、真偽値版では
/// 2 回目以降の減算が「既に外側」と衝突します（`IMPL-phase3.md` §2.4）。
/// **巻き数なら $w_T$ が 5 つ並ぶだけ**です。
///
/// T は W を $z$ 方向に貫くので、結果は**種数 1**（菱形の穴が開いた箱）になります。
inline std::vector<TriMesh> case18() {
    const TriMesh w =
        box(at(-121, 256), at(-121, 256), at(-121, 256), at(135, 256), at(135, 256), at(135, 256));
    const std::int32_t px[4] = {at(97, 256), at(7, 256), at(-83, 256), at(7, 256)};
    const std::int32_t py[4] = {at(7, 256), at(97, 256), at(7, 256), at(-83, 256)};
    const TriMesh t = prism4(px, py, at(-140, 256), at(140, 256));
    return {w, t};
}

/// 同一座標の頂点を溶接する。**連結しただけでは組合せ的に多様体のまま**なので、
/// 非多様体な入力を作るには添字を共有させる必要があります。
inline TriMesh weld(const TriMesh& m) {
    std::map<std::array<std::int32_t, 3>, std::uint32_t> id;
    TriMesh r;
    std::vector<std::uint32_t> remap(m.vertices.size());
    for (std::size_t i = 0; i < m.vertices.size(); ++i) {
        const IPoint& p = m.vertices[i];
        const std::array<std::int32_t, 3> k{p.x, p.y, p.z};
        const auto it = id.find(k);
        if (it != id.end()) {
            remap[i] = it->second;
        } else {
            remap[i] = static_cast<std::uint32_t>(r.vertices.size());
            id.emplace(k, remap[i]);
            r.vertices.push_back(p);
        }
    }
    for (const Tri& t : m.triangles) r.triangles.push_back({remap[t[0]], remap[t[1]], remap[t[2]]});
    return r;
}

/// ケース 19: **自己交差したメッシュ**（重なる 2 つの箱を 1 つのメッシュとして与える）。
///
/// **$\partial S = 0$ なので PWN です**（各箱が個別に閉じている）。self-union が
/// 修復として機能すること、すなわち**別々に与えて和を取った結果と一致すること**を突きます。
inline std::vector<TriMesh> case19() {
    const TriMesh a =
        box(at(-121, 256), at(-121, 256), at(-121, 256), at(39, 256), at(39, 256), at(39, 256));
    const TriMesh b =
        box(at(-25, 256), at(-25, 256), at(-25, 256), at(135, 256), at(135, 256), at(135, 256));
    return {concat(a, b), a, b};
}

/// ケース 20a: **PWN な非多様体入力**（辺を共有して貼り合わせた 2 つの閉曲面）。
///
/// 共有辺には 4 枚が集まりますが、**各箱が個別に $\partial = 0$ なので和も 0** です。
/// **`edge_manifold` は落ちますが PWN です。** 受け入れられることを突きます。
inline std::vector<TriMesh> case20a() {
    // **$z$ 方向の辺を共有します**（$x=7$, $y=7$ の鉛直線）。
    // 角ではなく辺で接するので、共有辺には 4 枚が集まります。
    const TriMesh a =
        box(at(-121, 256), at(-121, 256), at(-121, 256), at(7, 256), at(7, 256), at(135, 256));
    const TriMesh b =
        box(at(7, 256), at(7, 256), at(-121, 256), at(135, 256), at(135, 256), at(135, 256));
    // **溶接します。** 連結しただけでは頂点が重複し、組合せ的には多様体のままです
    return {weld(concat(a, b)), a, b};
}

/// ケース 20b: **PWN でない入力**（1 辺に 3 枚）。**拒否できることを突きます。**
///
/// 箱に三角形を 1 枚足すと、その辺の向き付き係数の和が $\pm1$ になり $\partial S \ne 0$。
/// **奇数個の $\pm1$ の和は 0 になりません。**
inline TriMesh case20b() {
    TriMesh m =
        box(at(-121, 256), at(-121, 256), at(-121, 256), at(135, 256), at(135, 256), at(135, 256));
    const auto v = static_cast<std::uint32_t>(m.vertices.size());
    m.vertices.push_back({at(7, 256), at(7, 256), at(200, 256)});
    m.triangles.push_back({m.triangles[0][0], m.triangles[0][1], v});
    return m;
}

/// ケース 23: **20 個の立方体が同じ上下面（と側面）を共有**（EMBER Fig. 25）。
///
/// **共平面重複が 20 枚**になります。`SPEC-phase2.md` §4.3.2 の選択規則は
/// 「同方向 / 逆方向 × 3 演算」の表で**二項**なので、この配置を通したことがありません。
/// **`SPEC-phase3.md` §5.4.1（全順序）の要否がここで決まります。**
///
/// $x$ 方向に 32 きざみで幅 40 の箱を 20 個並べます（**隣と 8 だけ重なる**）。
/// $y$ / $z$ は 20 個すべてで同一なので、4 枚の平面に 20 枚のシートが載ります。
///
/// **この配置は意図的に軸平行です。** 突きたいのは共平面重複の枚数であって
/// 4 平面同時交差ではありません（そちらはケース 17 / 2T / 15 が持っています）。
inline std::vector<TriMesh> case23() {
    std::vector<TriMesh> v;
    v.reserve(20);
    for (int i = 0; i < 20; ++i) {
        const int x0 = -300 + 32 * i;
        v.push_back(box(at(x0, 640), at(-300, 640), at(-300, 640), at(x0 + 40, 640), at(340, 640),
                        at(340, 640)));
    }
    return v;
}

/// ケース 23 の期待値: 20 個の和は 1 つの箱になります。
inline TriMesh case23_expect() {
    return box(at(-300, 640), at(-300, 640), at(-300, 640), at(-300 + 32 * 19 + 40, 640),
               at(340, 640), at(340, 640));
}

}  // namespace cases

/// $n$ 項のケース。**恒等式とセットで持ちます。**
struct NaryCase {
    enum class Kind {
        /// $(A \cup B) \setminus C = (A \setminus C) \cup (B \setminus C)$
        Distributive,
        /// $W - T - \dots - T = W - T$（`repeat` 回引く）
        RepeatedDifference,
        /// $\bigcup_i M_i = E$（`expect` が返す立体）
        UnionOfMany,
        /// $\mathrm{selfunion}(M_0) = M_1 \cup M_2$（`make` は $\{M_0, M_1, M_2\}$）
        SelfUnion,
    };
    const char* id;
    const char* what;
    Kind kind;
    std::vector<TriMesh> (*make)();
    int repeat = 0;
    /// `UnionOfMany` のときの期待値。**解析的に分かる立体を置くこと**
    /// （恒等式は自己整合の検査なので、答えのレベルで独立な対照が要ります）。
    TriMesh (*expect)() = nullptr;
};

inline const std::vector<NaryCase>& nary_corpus() {
    static const std::vector<NaryCase> k = {
        {"17", "n=3。共平面で接する 2 箱を斜めの角柱が貫く", NaryCase::Kind::Distributive,
         cases::case17, 0},
        {"18", "同一形状を 5 回引く（中間結果を経由しない）", NaryCase::Kind::RepeatedDifference,
         cases::case18, 5},
        {"19", "自己交差メッシュの self-union（修復として機能すること）", NaryCase::Kind::SelfUnion,
         cases::case19, 0},
        {"20a", "PWN な非多様体入力（辺を共有した 2 つの閉曲面）", NaryCase::Kind::SelfUnion,
         cases::case20a, 0},
        {"23", "20 個の立方体が同じ上下面を共有（共平面重複 20 枚）", NaryCase::Kind::UnionOfMany,
         cases::case23, 0, cases::case23_expect},
    };
    return k;
}

/// **スープ経路だけで回すケース**（`SPEC-phase3.md` §9 のケース 21）。
///
/// > **二項実装は非凸な面を扱えません。** `build_faces` が共平面の三角形を面に
/// > まとめ、**扇で三角形化**するので、非凸な面では別の立体になります
/// > （穴あきの面では境界ループが 2 本になり、そもそも構築が止まります）。
/// >
/// > **§4.1 の凸分割は性能の最適化ではなく、非凸な面を扱う唯一の経路です。**
/// > だから主コーパス（二項正解器と突き合わせるもの）には入れられません。
inline const std::vector<Case>& soup_only_corpus() {
    static const std::vector<Case> k = {
        {"21", "非凸の面（L 字柱 と 穴あき平板）", cases::case21_a, cases::case21_b},
    };
    return k;
}

/// §9.0 (1) のサイズ規律: 両入力の AABB の和が各軸で座標範囲の半分以上あるか。
inline bool size_discipline_ok(const TriMesh& a, const TriMesh& b) {
    std::int64_t lo[3] = {krisite::kCoordMax, krisite::kCoordMax, krisite::kCoordMax};
    std::int64_t hi[3] = {krisite::kCoordMin, krisite::kCoordMin, krisite::kCoordMin};
    for (const TriMesh* m : {&a, &b}) {
        for (const IPoint& p : m->vertices) {
            const std::int64_t c[3] = {p.x, p.y, p.z};
            for (int t = 0; t < 3; ++t) {
                lo[t] = (c[t] < lo[t]) ? c[t] : lo[t];
                hi[t] = (c[t] > hi[t]) ? c[t] : hi[t];
            }
        }
    }
    const std::int64_t need = -static_cast<std::int64_t>(krisite::kCoordMin);  // 2^(b-1)
    for (int t = 0; t < 3; ++t) {
        if (hi[t] - lo[t] < need) return false;
    }
    return true;
}

/// §9.0 (1) のサイズ規律の $n$ 項版。**全メッシュの AABB の和**で見ます。
inline bool size_discipline_ok(const std::vector<TriMesh>& ms) {
    std::int64_t lo[3] = {krisite::kCoordMax, krisite::kCoordMax, krisite::kCoordMax};
    std::int64_t hi[3] = {krisite::kCoordMin, krisite::kCoordMin, krisite::kCoordMin};
    for (const TriMesh& m : ms) {
        for (const IPoint& p : m.vertices) {
            const std::int64_t c[3] = {p.x, p.y, p.z};
            for (int t = 0; t < 3; ++t) {
                lo[t] = (c[t] < lo[t]) ? c[t] : lo[t];
                hi[t] = (c[t] > hi[t]) ? c[t] : hi[t];
            }
        }
    }
    const std::int64_t need = -static_cast<std::int64_t>(krisite::kCoordMin);  // 2^(b-1)
    for (int t = 0; t < 3; ++t) {
        if (hi[t] - lo[t] < need) return false;
    }
    return true;
}

}  // namespace kritest

#endif  // KRISITE_TESTS_CORPUS_HPP

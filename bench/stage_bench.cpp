// Krisite — 段別の費用係数 $c_b$ を合成入力で測る（`PERF.md` §1.1）
//
// **費用模型は $T = \sum_b W_b \times c_b$ です。** arrange の $W_b$ は
// $\sum_\ell P_\ell^2$ と確定しました（`IMPL-phase5.md` §31）が、**$c_b$ は
// 1 つも確定していません。** CP1.9 の前提がこれです。
//
// ## なぜ合成が要るか（`PERF.md` §1.2）
//
// 実データでは 3 つの量が**同時に動きます。**
//
//     Σ_ℓ P_ℓ²      葉あたりの仕事の総量
//     葉の数         ディスパッチのオーバーヘッド
//     分布の偏り     裾が支配する量
//
// **どれが効いているかを実データでは分けられません。** 合成でしか分けられない、
// というのが §1.2 の要点です。そこで**列を独立に振ります。**
//
//     つまみ 1  Σ_ℓ P_ℓ² を振る（葉の数を固定）      → c_b の値
//     つまみ 2  葉の数を振る（Σ を固定）             → ディスパッチが混じっていないか
//     つまみ 3  分布の偏りを振る（Σ と葉数を固定）   → 中央の 19% が説明できるか
//
// ## 完了条件（`PERF.md` §1.1）
//
// **合成が実データを代表しているか**を確かめるまで、$c_b$ は使えません。
//
//     arrange(錨)   対   c_b^合成 × Σ_ℓ P_ℓ²(錨)
//
//     合う     合成が実データを代表している。c_b が使える
//     合わない 作業集合・分布・メモリ配置が実データと違う。**c_b は使えません**
//
// > 合成で綺麗な数字が出ても、実データに当てはまらなければ意味がありません
// > （`CLAUDE.md`「合成で測った係数は、実測に当てはめて確かめる」）。
//
// ## 合成入力の作り方
//
// **葉 1 つに 1 つの塊**を置きます。深度 `d` のセルを `L` 個選び、各セルの中に
//
//     A … 一列に並べた `t` 個の小立方体（互いに離れている）
//     B … 半個ぶんずらした `t` 個の小立方体（A と食い違って重なる）
//
// を入れます。**多角形数は $P_\ell = 24t$**（三角形化なので 12 三角形 × 2t）。
//
// > **なぜ立方体の列か。** 頂点は格子点でなければならず、退化した三角形を
// > 作ってはいけません。多角形（櫛形など）を使うと**凹多角形の三角形化**が要り、
// > そこで退化を踏みます。**離れた立方体の合併なら PWN であることが自明**で、
// > 多角形数を `t` で厳密に決められます。
//
// **塊はセルに完全に収まります。** はみ出すと葉が分かれ、つまみ 2 と 3 が
// 混ざります（`leaf_nonempty` を出力に出しているのはそのためです）。
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/mesh/tri_mesh.hpp"

using namespace krisite;
using namespace krisite::csg;

namespace {

using Clock = std::chrono::steady_clock;

/// 深度 `d` のセルの一辺。
std::int64_t cell_size(unsigned d) noexcept {
    return std::int64_t{1} << (kCoordBits - d);
}

/// `(lo..hi)` の軸平行立方体を `m` に足す（**閉じているので合併しても PWN**）。
void add_box(mesh::TriMesh& m, std::int64_t lox, std::int64_t loy, std::int64_t loz,
             std::int64_t hix, std::int64_t hiy, std::int64_t hiz) {
    const std::uint32_t o = static_cast<std::uint32_t>(m.vertices.size());
    const auto v = [](std::int64_t x, std::int64_t y, std::int64_t z) {
        return geom::IPoint{static_cast<std::int32_t>(x), static_cast<std::int32_t>(y),
                            static_cast<std::int32_t>(z)};
    };
    m.vertices.push_back(v(lox, loy, loz));
    m.vertices.push_back(v(hix, loy, loz));
    m.vertices.push_back(v(hix, hiy, loz));
    m.vertices.push_back(v(lox, hiy, loz));
    m.vertices.push_back(v(lox, loy, hiz));
    m.vertices.push_back(v(hix, loy, hiz));
    m.vertices.push_back(v(hix, hiy, hiz));
    m.vertices.push_back(v(lox, hiy, hiz));
    static const std::uint32_t tri[12][3] = {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
                                             {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                                             {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    for (const auto& t : tri) m.triangles.push_back({o + t[0], o + t[1], o + t[2]});
}

/// 塊を 1 つ置く。`t` 個の立方体を x 方向に並べる。`shift` で B をずらす。
///
/// セル `(ci,cj,ck)` の内側に**余白を残して**収めます。
void add_cluster(mesh::TriMesh& m, unsigned d, int ci, int cj, int ck, int t, bool shift) {
    const std::int64_t cs = cell_size(d);
    const std::int64_t x0 = std::int64_t{kCoordMin} + std::int64_t{ci} * cs;
    const std::int64_t y0 = std::int64_t{kCoordMin} + std::int64_t{cj} * cs;
    const std::int64_t z0 = std::int64_t{kCoordMin} + std::int64_t{ck} * cs;
    // **セルの内側 1/2 に収めます。** 端に寄せると隣のセルに漏れます
    const std::int64_t span = cs / 2;
    const std::int64_t pitch = span / t;        // 立方体 1 個ぶんの間隔
    const std::int64_t side = (pitch * 3) / 4;  // 立方体の一辺（間に隙間を残す）
    const std::int64_t bx = x0 + cs / 4;
    const std::int64_t by = y0 + cs / 4;
    const std::int64_t bz = z0 + cs / 4;
    for (int k = 0; k < t; ++k) {
        // **B は半ピッチずらして A と食い違わせます。** 重なりが無いと
        // 交差の仕事がゼロになり、arrange の費用を測れません
        const std::int64_t sx = bx + std::int64_t{k} * pitch + (shift ? pitch / 2 : 0);
        add_box(m, sx, by, bz, sx + side, by + side, bz + side);
    }
}

struct Config {
    const char* knob;  ///< どのつまみを振ったか
    unsigned d;        ///< 塊を置く深度
    int L;             ///< 塊の数（= 空でない葉の数の期待値）
    int t;             ///< 塊あたりの立方体数（$P_\ell = 24t$）
    int t2;            ///< 偏りを作るときの第 2 の値（0 なら一様）
    int n_big;         ///< `t2` を使う塊の数
};

/// 深度 `d` の格子から `L` 個のセルを選ぶ（**離して置きます**）。
std::vector<std::array<int, 3>> pick_cells(unsigned d, int L) {
    const int n = 1 << d;
    std::vector<std::array<int, 3>> out;
    for (int i = 0; i < n && (int)out.size() < L; ++i) {
        for (int j = 0; j < n && (int)out.size() < L; ++j) {
            for (int k = 0; k < n && (int)out.size() < L; ++k) out.push_back({i, j, k});
        }
    }
    return out;
}

struct Run {
    double ms_arrange = 0, ms_classify = 0, ms_stitch = 0;
    std::size_t poly_sq = 0, input_sq = 0, leaves = 0, polys = 0;
    /// **局所 BSP の切断数。** 時間ではなく**仕事の量**で実データと比べるため
    std::size_t bsp_slots = 0, bsp_used = 0, single = 0;
    /// **無次元群**（`PERF.md` §1.8）。実データと揃っているかを見ます
    std::size_t frag_edges = 0, frag_edges_max = 0, leaf_planes = 0, leaf_in = 0, frags = 0;
    std::size_t frag_n = 0;
};

Run measure(const Config& c, unsigned threads, par::ThreadPool* pool, int reps) {
    const auto cells = pick_cells(c.d, c.L);
    mesh::TriMesh ma, mb;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const int t = (c.t2 > 0 && (int)i < c.n_big) ? c.t2 : c.t;
        add_cluster(ma, c.d, cells[i][0], cells[i][1], cells[i][2], t, false);
        add_cluster(mb, c.d, cells[i][0], cells[i][1], cells[i][2], t, true);
    }
    const PolySoup A = from_mesh(ma), B = from_mesh(mb);
    BoolOptions o;
    o.depth = c.d;  // **深度を塊の深度に合わせます。** 深く割ると葉数が増えます
    o.adaptive = true;
    o.local_bsp = true;
    o.ray_index = true;
    o.threads = threads;
    o.pool = pool;
    Run best;
    // **最小値を取ります**（`CLAUDE.md`「単発の測定値の差に原因を当てない」）
    for (int r = 0; r < reps; ++r) {
        BoolStats bs;
        const PolySoup out = boolean(A, B, BoolOp::Union, o, &bs);
        if (r == 0 || bs.ms_arrange < best.ms_arrange) {
            best.ms_arrange = bs.ms_arrange;
            best.ms_classify = bs.ms_classify;
            best.ms_stitch = bs.ms_stitch;
        }
        best.poly_sq = bs.leaf_poly_sq;
        best.input_sq = bs.leaf_input_sq;
        best.leaves = bs.leaf_nonempty;
        best.bsp_slots = bs.bsp_cut_slots;
        best.bsp_used = bs.bsp_cuts_used;
        best.single = bs.leaf_single_src;
        best.frag_edges = bs.frag_edges_total;
        best.frag_edges_max = bs.frag_edges_max;
        best.leaf_planes = bs.leaf_planes_total;
        best.leaf_in = bs.leaf_input_total;
        best.frags = bs.fragments;
        best.frag_n = bs.frag_edges_count;
        best.polys = out.polys.size();
    }
    return best;
}

}  // namespace

int main(int argc, char** argv) {
    const int reps = (argc > 1) ? std::atoi(argv[1]) : 3;
    const unsigned threads = (argc > 2) ? static_cast<unsigned>(std::atoi(argv[2])) : 1;
    const char* only = (argc > 3) ? argv[3] : nullptr;  // つまみの絞り込み（"深" など）
    // **設定を出力に書きます**（`IMPL-phase5.md` §29 の再発防止）。
    // 書かないと、あとで「どの構成で測ったか」が分かりません
    std::printf("=== 段別費用係数の合成ベンチ（PERF.md §1.1）===\n");
    std::printf("b=%d / 反復 %d（最小値）/ スレッド %u / 索引 ON / 局所 BSP ON / NSI 宣言なし\n\n",
                KRISITE_COORD_BITS, reps, threads);

    par::ThreadPool pool(threads);

    // **$\sum_\ell P_\ell^2 = L \times (24t)^2$**（偏らせるときは $576\sum_\ell t_\ell^2$）。
    // つまみ 2 と 3 では **9,437,184 に厳密に揃えて**あります。
    // **「ほぼ一定」では列が独立になりません。**
    const Config cfgs[] = {
        // つまみ 1: Σ P² を振る（葉の数は 8 で固定）。**c_b の値そのもの**
        {"Σ", 2, 8, 4, 0, 0},
        {"Σ", 2, 8, 8, 0, 0},
        {"Σ", 2, 8, 16, 0, 0},
        {"Σ", 2, 8, 32, 0, 0},
        {"Σ", 2, 8, 64, 0, 0},
        // つまみ 2: 葉の数を振る（Σ を厳密に固定。$t \propto 1/\sqrt{L}$）。
        // **ディスパッチのオーバーヘッドが c_b に混じっていないか**
        {"葉", 2, 1, 128, 0, 0},
        {"葉", 2, 4, 64, 0, 0},
        {"葉", 2, 16, 32, 0, 0},
        {"葉", 2, 64, 16, 0, 0},
        // つまみ 3: 偏りを振る（葉数 16 と Σ を厳密に固定）。**中央の 19%**
        //   $t_2^2 n + t^2 (16-n) = 16384$ を満たす整数解を並べています
        {"偏", 2, 16, 32, 0, 0},    // 一様。最大/平均 = 1.00
        {"偏", 2, 16, 28, 68, 1},   // 最大/平均 = 2.23
        {"偏", 2, 16, 24, 88, 1},   // 最大/平均 = 3.14
        {"偏", 2, 16, 16, 112, 1},  // 最大/平均 = 5.09
        {"偏", 2, 16, 8, 88, 2},    // 2 つ大きい（最大/平均 = 4.89）
        // つまみ 2 の続き: **実データの葉数まで上げる**（Σ は同じ 9,437,184）。
        //
        // > **`ThreadPool::kDefaultMinItems = 64`。** 葉が 64 未満だと arrange は
        // > **逐次に落ちます。** 深度 2 では葉が足りず、16 スレッドでも
        // > **並列経路を一度も踏みませんでした。** 錨は葉が数千あります。
        //
        // **完了条件（錨の再現）はこの段で判定します。**
        {"深", 3, 256, 8, 0, 0},
        {"深", 4, 1024, 4, 0, 0},
        {"深", 4, 4096, 2, 0, 0},
    };

    std::printf(
        "| つまみ | 葉 L | t | 大 t₂×n | 空でない葉 | **最大/平均** | 多角形 |"
        " **Σ P_ℓ²** | **BSP 枠** | **BSP 使用** | arrange ms |"
        " **辺/多角形** | **最大辺** | **平面/多角形** | **枠/ΣP²** |"
        " **c_arrange (ns/単位)** |\n");
    std::printf("|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for (const Config& c : cfgs) {
        // **つまみを絞れるようにします。** 足した段だけを回すため
        if (only != nullptr && std::string(c.knob) != only) continue;
        const Run r = measure(c, threads, &pool, reps);
        // **偏りを数字で出します。** 「偏らせた」と言うだけでは列が独立か分かりません
        const double mean =
            (c.t2 > 0) ? double(c.t2 * c.n_big + c.t * (c.L - c.n_big)) / double(c.L) : double(c.t);
        const double skew = (c.t2 > 0 ? double(std::max(c.t2, c.t)) : double(c.t)) / mean;
        char big[32] = "-";
        if (c.t2 > 0) std::snprintf(big, sizeof(big), "%d×%d", c.t2, c.n_big);
        const double u = r.poly_sq ? 1e6 / double(r.poly_sq) : 0.0;  // ms → ns/単位
        std::printf(
            "| %s | %d | %d | %s | %zu | %.2f | %zu | **%zu** | %zu | %zu | %.2f |"
            " **%.2f** | **%zu** | **%.2f** | **%.4f** | **%.3f** |\n",
            c.knob, c.L, c.t, big, r.leaves, skew, r.polys, r.poly_sq, r.bsp_slots, r.bsp_used,
            r.ms_arrange, r.frag_n ? double(r.frag_edges) / double(r.frag_n) : 0.0,
            r.frag_edges_max, r.leaf_in ? double(r.leaf_planes) / double(r.leaf_in) : 0.0,
            r.poly_sq ? double(r.bsp_slots) / double(r.poly_sq) : 0.0, r.ms_arrange * u);
        std::fflush(stdout);
    }
    return 0;
}

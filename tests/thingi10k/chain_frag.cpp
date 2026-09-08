// Krisite — **連鎖で断片数がどう増えるか**（`SPEC-phase5.md` §5.10.5.7、依頼 3）
//
// ## 何を測るか
//
// **`(A ∪ B) \ D` の 1 段目と 2 段目の断片数を、単発の `A \ D` と比べます。**
//
// **`CONTRACTS.md` は「連鎖でビット幅が伸びない」と書いていますが、
// 断片数については何も言っていません。**
//
// ## なぜ測るか
//
// `DESIGN-phase5-hotspots.md` §14.8.3 が**推測**として次を書きました。
//
// > 1 段目で境界に沿って切られた 2 つの断片が、2 段目で違うセルに入る。
//
// **これが正しければ、連鎖するほど断片が増えます。**
// **推測のままにせず、機構まで降りるための計器です**（`.claude/rules/deduction.md`）。
//
// ## 出力する量
//
// | 量 | 意味 |
// |---|---|
// | 入力多角形 | その段が受け取った `PolySoup::polys` の数（**膨張率の分母**） |
// | `raw_fragments` | 正準化前の断片数（重複を含む） |
// | `fragments` | 出力多角形数（= 領域の数） |
// | `region_cross_cell` | 同じ `region_key` の断片が複数のセルに分かれた群 |
// | `..._same_edges` | **そのうち、辺平面集合が完全に一致するもの**（同じ多角形） |
// | `..._axis_support` | **そのうち、支持平面が軸平行**（= セル境界面） |
//
// **最後の 2 つが機構の切り分けです。**
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/par/thread_pool.hpp"

#include "thingi10k/loader.hpp"
#include "volume_fp.hpp"

using namespace krisite;

namespace {

/// **出力の同一性を版をまたいで比べるためのハッシュ**（`SPEC-phase5.md` §5.10.6 の検査）。
///
/// **同一版の中の比較ではなく、コードを変える前後で比べるので、
/// 値そのものを出力に書きます**（`CLAUDE.md`「版をまたいだ比較」）。
unsigned long long hash_mesh(const csg::SoupMesh& m) {
    unsigned long long h = 1469598103934665603ull;
    const auto mix = [&h](unsigned long long v) {
        h ^= v;
        h *= 1099511628211ull;
    };
    mix(m.triangles.size());
    mix(m.vertices.size());
    for (const mesh::Tri& t : m.triangles) {
        for (int k = 0; k < 3; ++k) mix(t[k]);
    }
    for (const geom::HPointD& v : m.vertices) {
        // **同次座標の全リムを混ぜます**（値そのものを比べるため）
        for (std::size_t i = 0; i < v.x.kLimbs; ++i) mix(v.x[i]);
        for (std::size_t i = 0; i < v.y.kLimbs; ++i) mix(v.y[i]);
        for (std::size_t i = 0; i < v.z.kLimbs; ++i) mix(v.z[i]);
        for (std::size_t i = 0; i < v.w.kLimbs; ++i) mix(v.w[i]);
    }
    return h;
}

struct Row {
    std::string name;
    std::size_t in_polys = 0;
    csg::BoolStats st;
};

std::vector<Row> g_rows;

void run_stage(const char* name, const csg::PolySoup& X, const csg::PolySoup& Y, csg::BoolOp op,
               const csg::BoolOptions& o, csg::PolySoup* out) {
    Row r;
    r.name = name;
    r.in_polys = X.polys.size() + Y.polys.size();
    const csg::PolySoup s = csg::boolean(X, Y, op, o, &r.st);
    g_rows.push_back(r);
    if (out != nullptr) *out = s;
}

void print_rows() {
    std::printf("\n| 段 | 入力多角形 | `raw_fragments` | `fragments` | 膨張率 |\n");
    std::printf("|---|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        std::printf("| %s | %zu | %zu | %zu | **%.3f 倍** |\n", r.name.c_str(), r.in_polys,
                    r.st.raw_fragments, r.st.fragments,
                    r.in_polys == 0
                        ? 0.0
                        : static_cast<double>(r.st.fragments) / static_cast<double>(r.in_polys));
    }
    // ---- ★ 到達可能性解析の機会（`SPEC-phase5.md` §5.10.8.4）------------------
    //
    // **「捨てられる葉が 0 なら、案 A を実装する意味がありません。」**
    // **葉の数だけでなく $\sum P_\ell$ と $\sum P_\ell^2$ を見ます**
    // （局所 BSP が $O(P_\ell^2)$ なので、後者が本当の効きです）。
    std::printf(
        "\n| 段 | 葉（非空） | 確定あり | **定値** | $\\sum P$ 三角形 | $\\sum P^2$ 三角形 | "
        "$\\sum P$ 多角形 | $\\sum P^2$ 多角形 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const auto pc = [](std::size_t a, std::size_t b) {
            return b == 0 ? 0.0 : 100.0 * static_cast<double>(a) / static_cast<double>(b);
        };
        std::printf(
            "| %s | %zu | %zu | **%zu（%.1f%%）** | %.1f%% | **%.1f%%** | %.1f%% | "
            "**%.1f%%** |\n",
            r.name.c_str(), r.st.leaf_nonempty, r.st.eo_forced_leaves, r.st.eo_const_leaves,
            pc(r.st.eo_const_leaves, r.st.leaf_nonempty),
            pc(r.st.eo_const_input, r.st.leaf_input_total),
            pc(r.st.eo_const_input_sq, r.st.leaf_input_sq),
            pc(r.st.eo_const_polys, r.st.leaf_poly_total),
            pc(r.st.eo_const_poly_sq, r.st.leaf_poly_sq));
    }
    // **★ 代理ではなく、実際に省ける仕事**（`CLAUDE.md`「数えている量が費用の代理か」）。
    std::printf(
        "\n| 段 | 定値の葉の `bsp_cut_slots` | 全体 | **割合** | 定値の葉の断片 | 全体 | "
        "**割合** |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const auto pc = [](std::size_t a, std::size_t b) {
            return b == 0 ? 0.0 : 100.0 * static_cast<double>(a) / static_cast<double>(b);
        };
        std::printf("| %s | %zu | %zu | **%.1f%%** | %zu | %zu | **%.1f%%** |\n", r.name.c_str(),
                    r.st.eo_const_bsp_slots, r.st.bsp_cut_slots,
                    pc(r.st.eo_const_bsp_slots, r.st.bsp_cut_slots), r.st.eo_const_frags,
                    r.st.frag_edges_count, pc(r.st.eo_const_frags, r.st.frag_edges_count));
    }

    // ---- ★ 代表点の構成（`SPEC-phase5.md` §5.10.9.3）--------------------------
    //
    // **EMBER §4.4 の fast-path は、代表点の構成として既に実装されています**
    // （`interior.hpp` の主経路 = float ヒント + 軸平行直線）。
    // **各段が何回試され、何回成功しているかを数えます。**
    std::printf(
        "\n| 段 | 代表点の構成 | **主経路の成功** | 主経路の失敗 | 予備の成功 | "
        "予備の試行 | `side`（代表点） | `side`/構成 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const csg::InteriorStats& it = r.st.interior;
        const std::size_t tot = it.axis_line + it.corner_offset;
        std::printf(
            "| %s | %zu | **%zu（%.1f%%）** | %zu | %zu | %zu | %zu | %.1f |\n", r.name.c_str(),
            tot, it.axis_line,
            tot == 0 ? 0.0 : 100.0 * static_cast<double>(it.axis_line) / static_cast<double>(tot),
            it.axis_failed, it.corner_offset, it.corner_tries, it.side_tests,
            tot == 0 ? 0.0 : static_cast<double>(it.side_tests) / static_cast<double>(tot));
    }

    // **★ 主経路の失敗の内訳**（§20.2。「外れた」だけでは機構が決まりません）。
    std::printf(
        "\n| 段 | 主経路の失敗 | 範囲外 | **内部でない** | 退化 | 辺の平均 | "
        "**範囲外の最大 $|c|/\\mathrm{max}$（‰）** |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const csg::InteriorStats& it = r.st.interior;
        std::printf("| %s | %zu | %zu | **%zu** | %zu | %.2f | **%zu** |\n", r.name.c_str(),
                    it.axis_failed, it.axis_out_of_range, it.axis_outside, it.axis_degenerate,
                    r.st.frag_edges_count == 0 ? 0.0
                                               : static_cast<double>(r.st.frag_edges_total) /
                                                     static_cast<double>(r.st.frag_edges_count),
                    it.axis_range_max_permille);
    }

    // ---- ★ 段ごとの `side`（`KRISITE_COUNT_PREDICATES` のときだけ）------------
    if (g_rows[0].st.side_calls_arrange + g_rows[0].st.side_calls_classify > 0) {
        std::printf(
            "\n| 段 | `side` arrange | `side` 分類 | うち代表点 | **代表点の割合（全体）** "
            "| `intersect3` arrange | `intersect3` 分類 |\n");
        std::printf("|---|---:|---:|---:|---:|---:|---:|\n");
        for (const Row& r : g_rows) {
            const std::uint64_t tot = r.st.side_calls_arrange + r.st.side_calls_classify;
            std::printf("| %s | %llu | %llu | %zu | **%.2f%%** | %llu | %llu |\n", r.name.c_str(),
                        static_cast<unsigned long long>(r.st.side_calls_arrange),
                        static_cast<unsigned long long>(r.st.side_calls_classify),
                        r.st.interior.side_tests,
                        tot == 0 ? 0.0
                                 : 100.0 * static_cast<double>(r.st.interior.side_tests) /
                                       static_cast<double>(tot),
                        static_cast<unsigned long long>(r.st.intersect3_arrange),
                        static_cast<unsigned long long>(r.st.intersect3_classify));
        }
    }

    // ---- ★★ `side` の被符号値の【実際の】幅（`SPEC-phase5.md` §5.10.10 の案 E）------
    //
    // **測るのは最終的な被符号値の幅です**（中間結果の最大ではありません）。
    // **上界と並べて出します** — 張り付いていれば案 E は成立しません。
    if (g_rows[0].st.side_calls_arrange + g_rows[0].st.side_calls_classify > 0) {
        std::printf("\n**上界**: `bits::kSide` = %zu ビット（$9b+20$、$b$ = %d）→ %zu リム\n",
                    static_cast<std::size_t>(geom::bits::kSide), KRISITE_COORD_BITS,
                    static_cast<std::size_t>(geom::limbs::kSide));
        std::printf("\n| 段 | 演算 | `side` | **≤64** | **≤128** | **≤192** | **>192** | 最大 |\n");
        std::printf("|---|---|---:|---:|---:|---:|---:|---:|\n");
        for (const Row& r : g_rows) {
            const auto row = [&](const char* stage, std::uint64_t n, std::uint64_t a,
                                 std::uint64_t b, std::uint64_t c, std::uint64_t d) {
                const double t = static_cast<double>(n);
                std::printf(
                    "| %s | %s | %llu | **%.1f%%** | **%.1f%%** | **%.1f%%** | "
                    "**%.1f%%** | %llu |\n",
                    r.name.c_str(), stage, static_cast<unsigned long long>(n),
                    n == 0 ? 0.0 : 100.0 * static_cast<double>(a) / t,
                    n == 0 ? 0.0 : 100.0 * static_cast<double>(b) / t,
                    n == 0 ? 0.0 : 100.0 * static_cast<double>(c) / t,
                    n == 0 ? 0.0 : 100.0 * static_cast<double>(d) / t,
                    static_cast<unsigned long long>(r.st.side_wmax));
            };
            row("arrange",
                r.st.side_w64_arrange + r.st.side_w128_arrange + r.st.side_w192_arrange +
                    r.st.side_wmore_arrange,
                r.st.side_w64_arrange, r.st.side_w128_arrange, r.st.side_w192_arrange,
                r.st.side_wmore_arrange);
            row("分類",
                r.st.side_w64_classify + r.st.side_w128_classify + r.st.side_w192_classify +
                    r.st.side_wmore_classify,
                r.st.side_w64_classify, r.st.side_w128_classify, r.st.side_w192_classify,
                r.st.side_wmore_classify);
        }
    }

    // ---- ★ 分類の費用（依頼 2 の材料）--------------------------------------
    //
    // **参照点の伝播で置き換えたいのは、この「領域ごとの大域レイキャスト」です。**
    // **局所トレースの費用の上限は「葉の中の多角形数」なので、両方を並べます。**
    std::printf(
        "\n| 段 | 領域（レイ） | 隅のレイ | 三角形検査 | /レイ | 葉あたり多角形 | 分類の割合 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const double total = r.st.ms_prepare + r.st.ms_arrange + r.st.ms_stitch + r.st.ms_classify;
        std::printf("| %s | %zu | %zu | %zu | %.1f | %.1f | %.1f%% |\n", r.name.c_str(),
                    r.st.raycasts, r.st.early_out_raycasts, r.st.ray_tri_tests,
                    r.st.raycasts == 0 ? 0.0
                                       : static_cast<double>(r.st.ray_tri_tests) /
                                             static_cast<double>(r.st.raycasts),
                    r.st.leaf_nonempty == 0 ? 0.0
                                            : static_cast<double>(r.st.raw_fragments) /
                                                  static_cast<double>(r.st.leaf_nonempty),
                    total == 0 ? 0.0 : 100.0 * r.st.ms_classify / total);
    }

    // **★ 機構の候補を切り分けるための量**（§15.2）。
    // **葉の数と深さが段で変わるなら、増分はセル境界の切断で説明できます。**
    std::printf(
        "\n| 段 | 葉（非空） | 深度 min/max | `bsp_cuts_used` | `active_cells` | **箱が狭まった** "
        "|\n");
    std::printf("|---|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        std::printf("| %s | %zu | %u / %u | %zu | %zu | %zu |\n", r.name.c_str(),
                    r.st.leaf_nonempty, r.st.leaf_depth_min, r.st.leaf_depth_max,
                    r.st.bsp_cuts_used, r.st.active_cells, r.st.out_aabb_narrowed);
    }
    std::printf(
        "\n| 段 | 群 | 食い違い | **セルまたぎ** | うち辺平面が同一 | うち支持平面が軸平行 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        std::printf("| %s | %zu | %zu | **%zu** | %zu | %zu |\n", r.name.c_str(),
                    r.st.region_cmp_groups, r.st.region_cmp_mismatch, r.st.region_cross_cell,
                    r.st.region_cross_cell_same_edges, r.st.region_cross_cell_axis_support);
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::setvbuf(stdout, nullptr, _IOLBF, 4096);
    if (argc < 4) {
        std::fprintf(stderr, "使い方: chain_frag <idA> <idB> <idD> [depth] [threads]\n");
        return 2;
    }
    const std::string ida = argv[1], idb = argv[2], idd = argv[3];
    const unsigned depth = (argc > 4) ? static_cast<unsigned>(std::atoi(argv[4])) : 5;
    const unsigned nthreads = (argc > 5) ? static_cast<unsigned>(std::atoi(argv[5])) : 8;
    // **★ 適応分割を切る経路**（§15.3 の予言を試すため）。
    // **固定深度なら段 1 と段 2 の格子が一致します。**
    const bool adaptive = (argc > 6) ? (std::atoi(argv[6]) != 0) : true;

    // **★ 設定と対象を最初に出します**（`CLAUDE.md`）
    std::printf("\n## 連鎖の断片数（`SPEC-phase5.md` §5.10.5.7 / 依頼 3）\n\n");
    std::printf("| 設定 | 値 |\n|---|---|\n");
    std::printf("| 対象 | A=`%s` B=`%s` D=`%s` |\n", ida.c_str(), idb.c_str(), idd.c_str());
    std::printf("| 深度 | %u（%s） |\n", depth, adaptive ? "**適応**" : "**固定**");
    std::printf("| スレッド | %u |\n", nthreads);
    std::printf("| b（座標ビット） | %d |\n", KRISITE_COORD_BITS);
    std::printf("| 走らせる段の数 | **3**（単発 1 段 + 連鎖 2 段） |\n\n");

    auto load = [](const std::string& id, unsigned seed_base) {
        return krithingi::quantize(
            krithingi::load_kmesh("data/thingi10k/kmesh/" + id + ".kmesh"),
            krithingi::make_transform(seed_base +
                                      static_cast<unsigned>(std::atoi(id.c_str())) % 997));
    };
    const auto qa = load(ida, 1000), qb = load(idb, 2000), qd = load(idd, 3000);
    if (qa.mesh.triangles.empty() || qb.mesh.triangles.empty() || qd.mesh.triangles.empty()) {
        std::fprintf(stderr, "模型を読めませんでした\n");
        return 3;
    }
    std::printf("**入力**: `%s` %zu / `%s` %zu / `%s` %zu 三角形\n", ida.c_str(),
                qa.mesh.triangles.size(), idb.c_str(), qb.mesh.triangles.size(), idd.c_str(),
                qd.mesh.triangles.size());

    // **★ 量子化後の座標範囲を出します**（`DESIGN` §20.2 の切り分け）。
    // **`interior_point` の主経路は、重心が `(kCoordMin, kCoordMax)` の【開区間】に
    // 無いと弾かれます。** 端に寄っていれば、そこが失敗の理由です。
    {
        const auto span = [](const mesh::TriMesh& m) {
            std::int64_t lo = krisite::kCoordMax, hi = krisite::kCoordMin;
            for (const geom::IPoint& v : m.vertices) {
                const std::int64_t c[3] = {v.x, v.y, v.z};
                for (int k = 0; k < 3; ++k) {
                    lo = std::min(lo, c[k]);
                    hi = std::max(hi, c[k]);
                }
            }
            return std::pair<std::int64_t, std::int64_t>{lo, hi};
        };
        const auto sa = span(qa.mesh), sd = span(qd.mesh);
        std::printf("\n**座標範囲**: `kCoordMax` = %lld / A = [%lld, %lld] / D = [%lld, %lld]\n",
                    static_cast<long long>(krisite::kCoordMax), static_cast<long long>(sa.first),
                    static_cast<long long>(sa.second), static_cast<long long>(sd.first),
                    static_cast<long long>(sd.second));
    }

    par::ThreadPool pool(nthreads);
    csg::BoolOptions o;
    o.depth = depth;
    o.adaptive = adaptive;
    o.cull_planes = true;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    o.threads = nthreads;
    o.pool = &pool;
    // **突き合わせだけ**（出力は変えません）。セルまたぎの計数はこの旗の下です。
    o.verify_region_key = true;

    const csg::PolySoup A = csg::from_mesh(qa.mesh);
    const csg::PolySoup B = csg::from_mesh(qb.mesh);
    const csg::PolySoup D = csg::from_mesh(qd.mesh);

    // ---- 単発（対照）--------------------------------------------------------
    csg::PolySoup u1;
    run_stage("単発 `A＼D`", A, D, csg::BoolOp::Difference, o, &u1);
    // ---- ★ 冪等な連鎖（決め手の対照。§15.1 の予測 2）------------------------
    //
    // **`(A＼D)＼D` の幾何は `A＼D` と同一のはず**（差は冪等）。
    // **にもかかわらず断片が増えるなら、増分は【幾何】ではなく【表現】から来ています。**
    csg::PolySoup u2;
    run_stage("**冪等** `(A＼D)＼D`", u1, D, csg::BoolOp::Difference, o, &u2);
    // ---- 連鎖 ---------------------------------------------------------------
    csg::PolySoup s1;
    run_stage("連鎖 1 段目 `A∪B`", A, B, csg::BoolOp::Union, o, &s1);
    run_stage("連鎖 2 段目 `(A∪B)＼D`", s1, D, csg::BoolOp::Difference, o, nullptr);

    print_rows();

    // ---- ★ 冪等の対照は【幾何が同じか】まで見ます ---------------------------
    //
    // **断片数だけでは「表現の問題」と「幾何が変わった」を区別できません。**
    {
        csg::ToMeshOptions tm;
        tm.split_contacts = true;
        const csg::SoupMesh m1 = csg::to_mesh(u1, tm);
        const csg::SoupMesh m2 = csg::to_mesh(u2, tm);
        const mesh::TopologyReport r1 = mesh::check_topology(m1.triangles);
        const mesh::TopologyReport r2 = mesh::check_topology(m2.triangles);
        std::printf("\n### 冪等の対照 — 幾何が同じか\n\n");
        std::printf(
            "| | 三角形 | 頂点 | 成分 $C$ | $\\chi$ | 多様体 |\n|---|---:|---:|---:|---:|---|\n");
        std::printf("| `A＼D` | %zu | %zu | %zu | %lld | %s |\n", m1.triangles.size(),
                    m1.vertices.size(), r1.components, static_cast<long long>(r1.chi),
                    r1.ok() ? "はい" : "**いいえ**");
        std::printf("| `(A＼D)＼D` | %zu | %zu | %zu | %lld | %s |\n", m2.triangles.size(),
                    m2.vertices.size(), r2.components, static_cast<long long>(r2.chi),
                    r2.ok() ? "はい" : "**いいえ**");
        std::printf("\n**★ ハッシュ**: `A＼D` = `%016llx` / `(A＼D)＼D` = `%016llx`\n",
                    hash_mesh(m1), hash_mesh(m2));
        // **★ バイト一致しないときに「意味論は同じ」を示す側**（`SPEC-phase5.md` §5.10.6）。
        //
        // **これは篩です。厳密な検査ではありません**（`volume_fp.hpp` の注記）。
        // **格子 1 単位ぶんの欠損は丸めの床より下で、捕まりません。**
        {
            const double v1 = kritest::volume6_fp(m1), v2 = kritest::volume6_fp(m2);
            const double scale = std::max(std::fabs(v1), 1.0);
            const double err = std::fabs(v2 - v1) / scale;
            std::printf(
                "**体積（篩）**: `A＼D` = %.10e / `(A＼D)＼D` = %.10e / "
                "相対差 **%.3e**（閾値 %.0e → %s）\n",
                v1, v2, err, kritest::kIdentityTol,
                err <= kritest::kIdentityTol ? "**通過**" : "**★ 超過（重大）**");
        }
        const bool same_geom = (r1.components == r2.components) && (r1.chi == r2.chi);
        std::printf(
            "\n**判定**: 位相は %s。断片数は %s。\n", same_geom ? "**一致**" : "**不一致（重大）**",
            g_rows[1].st.fragments > g_rows[0].st.fragments ? "**増えた**" : "増えていない");
    }

    // **★ 対照の要点**: 2 段目の入力は 1 段目の【出力】です。
    // **単発の A（`from_mesh` が作った tight な外接箱）と、
    // 1 段目の出力（外接箱が【セル箱】）を比べてください**（`polysoup.hpp` の `Poly::aabb`）。
    std::printf(
        "\n**注**: 2 段目の入力多角形 %zu 枚のうち、`soup_boolean` の出力は "
        "外接箱が**セル箱**です（`Poly::aabb` は「保守的」）。\n",
        s1.polys.size());
    return 0;
}

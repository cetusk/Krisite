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
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <map>
#include <new>
#include <string>
#include <thread>
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

// ---- ★★ 記憶の確保を数える（`SPEC-phase5.md` §5.10.11）------------------------
//
// **中核の CPU 時間の 86% が述語でないことが分かったので、種類別に分けます。**
//
// > **プロファイラ（`perf` / `valgrind`）はこの環境にありません。**
// > **計数で組みます。** 大域の `operator new` を置き換えて回数と量を数え、
// > **単価は同じ機械で測ります。**
//
// **これは駆動（テスト側）の実装で、ライブラリ本体には触れていません。**
// **`std::map` の節点も個別に確保されるので、この計数に含まれます。**
namespace kricount {

std::atomic<std::uint64_t> alloc_count{0};
std::atomic<std::uint64_t> alloc_bytes{0};
/// **★ 発生箇所（印）ごとの計数**（§5.10.12.4。残る 93.4% を刻む）。
constexpr int kTags = 20;
std::atomic<std::uint64_t> by_tag[kTags] = {};
std::atomic<std::uint64_t> bytes_by_tag[kTags] = {};
const char* const kTagName[kTags] = {"その他（印なし）",
                                     "前処理",
                                     "葉の列挙",
                                     "arrange: 設定",
                                     "arrange: 断片の生成（frag の複製 + pieces）",
                                     "arrange: 共平面",
                                     "arrange: 出力へ",
                                     "縫合: 構成点",
                                     "縫合: 仕分け（ids / packed）",
                                     "分類: 準備",
                                     "分類: 代表点",
                                     "分類: レイ",
                                     "分類: 出力の多角形",
                                     "early-out の隅",
                                     "arrange: 切断ループ（next）",
                                     "arrange: split の中（新 edge）",
                                     "縫合: 仕分け（map の節点 + 値）",
                                     "",
                                     "",
                                     ""};
/// **数えるのは中核の中だけ**（駆動自身の確保を混ぜないため）。
std::atomic<bool> enabled{false};

}  // namespace kricount

void* operator new(std::size_t n) {
    if (kricount::enabled.load(std::memory_order_relaxed)) {
        kricount::alloc_count.fetch_add(1, std::memory_order_relaxed);
        kricount::alloc_bytes.fetch_add(n, std::memory_order_relaxed);
#if defined(KRISITE_COUNT_PREDICATES)
        const int tag = krisite::geom::counters::alloc_tag;
        if (tag >= 0 && tag < kricount::kTags) {
            kricount::by_tag[tag].fetch_add(1, std::memory_order_relaxed);
            kricount::bytes_by_tag[tag].fetch_add(n, std::memory_order_relaxed);
        }
#endif
    }
    void* p = std::malloc(n == 0 ? 1 : n);
    if (p == nullptr) throw std::bad_alloc();
    return p;
}
void operator delete(void* p) noexcept {
    std::free(p);
}
void operator delete(void* p, std::size_t) noexcept {
    std::free(p);
}

namespace {

/// **出力の同一性を版をまたいで比べるためのハッシュ**（`SPEC-phase5.md` §5.10.6 の検査）。
///
/// **同一版の中の比較ではなく、コードを変える前後で比べるので、
/// 値そのものを出力に書きます**（`CLAUDE.md`「版をまたいだ比較」）。
/// **順序を除いた鍵**（`tests/csg/test_soup.cpp` の `geometric_key` と同じ形）。
/// 頂点を値で正準化し、三角形を最小の頂点から始めて整列します。
std::vector<std::array<std::uint32_t, 3>> geometric_key(const csg::SoupMesh& m) {
    std::vector<std::uint32_t> ord(m.vertices.size());
    for (std::uint32_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](std::uint32_t a, std::uint32_t b) {
        return geom::lex_less(m.vertices[a], m.vertices[b]);
    });
    std::vector<std::uint32_t> canon(m.vertices.size(), 0);
    std::uint32_t next = 0;
    for (std::size_t i = 0; i < ord.size();) {
        std::size_t j = i;
        while (j < ord.size() && geom::h_equal(m.vertices[ord[i]], m.vertices[ord[j]])) {
            canon[ord[j]] = next;
            ++j;
        }
        ++next;
        i = j;
    }
    std::vector<std::array<std::uint32_t, 3>> out;
    out.reserve(m.triangles.size());
    for (const mesh::Tri& t : m.triangles) {
        const std::uint32_t c[3] = {canon[t[0]], canon[t[1]], canon[t[2]]};
        int s = 0;
        if (c[1] < c[s]) s = 1;
        if (c[2] < c[s]) s = 2;
        out.push_back({c[s], c[(s + 1) % 3], c[(s + 2) % 3]});
    }
    std::sort(out.begin(), out.end());
    return out;
}

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
    /// **★ 見積もりの分母**（`SPEC-phase5.md` §5.10.10）。
    ///
    /// **`side` が中核に占める割合を出すには、中核の【CPU 時間】が要ります。**
    /// **壁時計では、スレッド数で割られた値と `side` の CPU 時間を比べることになります。**
    double wall_s = 0.0, cpu_s = 0.0;
    /// **記憶の確保**（`operator new` の回数と量）。
    std::uint64_t allocs = 0, alloc_bytes = 0;
    /// **逐次部分も含む述語の計数**（`cmp_h` は縫合の整列で効きます）。
    std::uint64_t cmp_h = 0, side_ip = 0;
    /// **発生箇所ごとの確保**（§5.10.12.4。残る 93.4% を刻む）。
    std::uint64_t tag_count[20] = {}, tag_bytes[20] = {};
    /// **断片の生成の内訳**（§5.10.12.4。`edge` の small-array の見積もり）。
    std::uint64_t split_calls = 0, split_early = 0, split_both = 0, make_calls = 0, make_ok = 0,
                  clip_calls = 0;
    std::uint64_t edge_hist[10] = {};
};

std::vector<Row> g_rows;
/// **並列効率を出すためのスレッド数**（`print_rows` から見えるように保持します）。
unsigned g_threads = 1;

void run_stage(const char* name, const csg::PolySoup& X, const csg::PolySoup& Y, csg::BoolOp op,
               const csg::BoolOptions& o, csg::PolySoup* out) {
    Row r;
    r.name = name;
    r.in_polys = X.polys.size() + Y.polys.size();
    const std::uint64_t a0 = kricount::alloc_count.load(std::memory_order_relaxed);
    const std::uint64_t b0 = kricount::alloc_bytes.load(std::memory_order_relaxed);
#if defined(KRISITE_COUNT_ALLOC) && defined(KRISITE_COUNT_PREDICATES)
    std::uint64_t t0c[20], t0b[20];
    for (int k = 0; k < 20; ++k) {
        t0c[k] = kricount::by_tag[k].load(std::memory_order_relaxed);
        t0b[k] = kricount::bytes_by_tag[k].load(std::memory_order_relaxed);
    }
#endif
#if defined(KRISITE_COUNT_PREDICATES)
    const std::uint64_t h0 = geom::counters::cmp_h_calls.load(std::memory_order_relaxed);
    const std::uint64_t i0 = geom::counters::side_ipoint_calls.load(std::memory_order_relaxed);
    const auto ld = [](const std::atomic<std::uint64_t>& a) {
        return a.load(std::memory_order_relaxed);
    };
    const std::uint64_t f0[6] = {
        ld(geom::counters::frag_split_calls), ld(geom::counters::frag_split_early),
        ld(geom::counters::frag_split_both),  ld(geom::counters::frag_make_calls),
        ld(geom::counters::frag_make_ok),     ld(geom::counters::frag_clip_calls)};
    std::uint64_t e0[10];
    for (int k = 0; k < 10; ++k) e0[k] = ld(geom::counters::frag_edge_hist[k]);
#endif
    kricount::enabled.store(true, std::memory_order_relaxed);
    const auto t0 = std::chrono::steady_clock::now();
    const std::clock_t c0 = std::clock();
    const csg::PolySoup s = csg::boolean(X, Y, op, o, &r.st);
    r.cpu_s = static_cast<double>(std::clock() - c0) / CLOCKS_PER_SEC;
    r.wall_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    kricount::enabled.store(false, std::memory_order_relaxed);
#if defined(KRISITE_COUNT_PREDICATES)
    r.cmp_h = geom::counters::cmp_h_calls.load(std::memory_order_relaxed) - h0;
    r.side_ip = geom::counters::side_ipoint_calls.load(std::memory_order_relaxed) - i0;
    r.split_calls = ld(geom::counters::frag_split_calls) - f0[0];
    r.split_early = ld(geom::counters::frag_split_early) - f0[1];
    r.split_both = ld(geom::counters::frag_split_both) - f0[2];
    r.make_calls = ld(geom::counters::frag_make_calls) - f0[3];
    r.make_ok = ld(geom::counters::frag_make_ok) - f0[4];
    r.clip_calls = ld(geom::counters::frag_clip_calls) - f0[5];
    for (int k = 0; k < 10; ++k) r.edge_hist[k] = ld(geom::counters::frag_edge_hist[k]) - e0[k];
#endif
    r.allocs = kricount::alloc_count.load(std::memory_order_relaxed) - a0;
    r.alloc_bytes = kricount::alloc_bytes.load(std::memory_order_relaxed) - b0;
#if defined(KRISITE_COUNT_ALLOC) && defined(KRISITE_COUNT_PREDICATES)
    for (int k = 0; k < 20; ++k) {
        r.tag_count[k] = kricount::by_tag[k].load(std::memory_order_relaxed) - t0c[k];
        r.tag_bytes[k] = kricount::bytes_by_tag[k].load(std::memory_order_relaxed) - t0b[k];
    }
#endif
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

    // ---- ★★ E2 の判定が【選ぶ】リム数（実際の幅と並べます）--------------------
    if (g_rows[0].st.side_disp1 + g_rows[0].st.side_disp2 + g_rows[0].st.side_disp3 +
            g_rows[0].st.side_disp4 >
        0) {
        std::printf("\n| 段 | 判定の対象 | **1 リム** | **2 リム** | **3 リム** | **4 リム** |\n");
        std::printf("|---|---:|---:|---:|---:|---:|\n");
        for (const Row& r : g_rows) {
            const double t = static_cast<double>(r.st.side_disp1 + r.st.side_disp2 +
                                                 r.st.side_disp3 + r.st.side_disp4);
            std::printf("| %s | %.0f | **%.1f%%** | **%.1f%%** | **%.1f%%** | **%.1f%%** |\n",
                        r.name.c_str(), t, t == 0 ? 0.0 : 100.0 * (double)r.st.side_disp1 / t,
                        t == 0 ? 0.0 : 100.0 * (double)r.st.side_disp2 / t,
                        t == 0 ? 0.0 : 100.0 * (double)r.st.side_disp3 / t,
                        t == 0 ? 0.0 : 100.0 * (double)r.st.side_disp4 / t);
        }
    }

    // ---- ★★ 見積もり: リム乗算の回数（時間を測る前に出します）------------------
    if (g_rows[0].st.side_mul_now > 0) {
        std::printf(
            "\n| 段 | `side` | いまのリム乗算 | **E2 のリム乗算** | **比** | "
            "**判定が読むリム** | /呼び出し |\n");
        std::printf("|---|---:|---:|---:|---:|---:|---:|\n");
        for (const Row& r : g_rows) {
            const std::uint64_t calls =
                r.st.side_disp1 + r.st.side_disp2 + r.st.side_disp3 + r.st.side_disp4;
            std::printf("| %s | %llu | %llu | %llu | **%.2f 倍** | %llu | **%.1f** |\n",
                        r.name.c_str(), static_cast<unsigned long long>(calls),
                        static_cast<unsigned long long>(r.st.side_mul_now),
                        static_cast<unsigned long long>(r.st.side_mul_e2),
                        r.st.side_mul_e2 == 0 ? 0.0
                                              : static_cast<double>(r.st.side_mul_now) /
                                                    static_cast<double>(r.st.side_mul_e2),
                        static_cast<unsigned long long>(r.st.side_disp_limbreads),
                        calls == 0 ? 0.0
                                   : static_cast<double>(r.st.side_disp_limbreads) /
                                         static_cast<double>(calls));
        }
    }

    // ---- ★★ 見積もりの分母: 中核の CPU 時間と、`side` が占める割合 --------------
    //
    // **`side` の単価は `BENCH.md` の 7.80 ns（b=21、4 リム）を使います。**
    // **これは【別の測定】なので、由来は「推測」です**（`.claude/rules/provenance.md`）。
    std::printf(
        "\n| 段 | 壁時計 | **CPU 時間** | 並列効率 | `side` の CPU（7.80 ns/回） | "
        "**中核に占める割合** |\n");
    std::printf("|---|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const std::uint64_t calls =
            r.st.side_disp1 + r.st.side_disp2 + r.st.side_disp3 + r.st.side_disp4;
        const double side_s = static_cast<double>(calls) * 7.80e-9;
        std::printf("| %s | %.3f s | **%.3f s** | %.2f | %.3f s | **%.1f%%** |\n", r.name.c_str(),
                    r.wall_s, r.cpu_s,
                    r.wall_s == 0.0 ? 0.0 : r.cpu_s / (r.wall_s * static_cast<double>(g_threads)),
                    side_s, r.cpu_s == 0.0 ? 0.0 : 100.0 * side_s / r.cpu_s);
    }

    // ---- ★★★ `std::map` の探索の単価（同じ機械・同じ規模で測ります）--------------
    //
    // **`PointCache` は `std::map<array<PlaneId,3>, Entry>` で、
    // 実測でエントリ 70 万・78 MB。** 木の探索は log n 段のポインタ追跡で、
    // **この大きさはキャッシュに載りません。**
    double map_ns = 0.0;
    {
        struct Entry96 {
            std::uint64_t pad[12];
        };  // `HPointD` に近い大きさ（96 バイト）
        std::map<std::array<std::uint32_t, 3>, Entry96> probe;
        // **★ `PointCache` はスレッド局所です**（`tl_cache[tid]`）。
        // **エントリ数は「全体 ÷ スレッド数」なので、木は小さくなります。**
        const std::size_t kN = 88000;
        std::uint64_t rs = 12345;
        const auto rnd = [&rs] {
            rs = rs * 6364136223846793005ull + 1442695040888963407ull;
            return static_cast<std::uint32_t>(rs >> 33);
        };
        std::vector<std::array<std::uint32_t, 3>> keys;
        keys.reserve(kN);
        for (std::size_t i = 0; i < kN; ++i) {
            std::array<std::uint32_t, 3> k{rnd() % 400000, rnd() % 400000, rnd() % 400000};
            std::sort(k.begin(), k.end());
            keys.push_back(k);
            probe.emplace(k, Entry96{});
        }
        const std::size_t kQ = 5000000;
        std::size_t found = 0;
        const std::clock_t m0 = std::clock();
        for (std::size_t i = 0; i < kQ; ++i) {
            if (probe.find(keys[(i * 2654435761u) % keys.size()]) != probe.end()) ++found;
        }
        const double ns_rand = 1e9 * static_cast<double>(std::clock() - m0) / CLOCKS_PER_SEC / kQ;
        // **(b) 局所性あり**（64 個の作業集合）= **下限**。
        // **`split_fragment` は同じ断片の頂点を続けて引くので、実際は両者の間です。**
        const std::clock_t m1 = std::clock();
        for (std::size_t i = 0; i < kQ; ++i) {
            if (probe.find(keys[(i % 64) + 1000]) != probe.end()) ++found;
        }
        const double ns_local = 1e9 * static_cast<double>(std::clock() - m1) / CLOCKS_PER_SEC / kQ;
        map_ns = 0.5 * (ns_rand + ns_local);
        std::printf(
            "\n| 引き方 | ns/回 |\n|---|---:|\n| **局所性なし**（上限） | **%.1f** |\n"
            "| **局所性あり**（下限） | **%.1f** |\n| 中間（以下で使います） | %.1f |\n",
            ns_rand, ns_local, map_ns);
        std::printf(
            "\n**`std::map` の探索の単価**（%zu 節点、96 バイトの値、この機械で実測）: "
            "**%.1f ns**（命中 %zu / %zu）\n",
            probe.size(), map_ns, found, kQ);
    }

    // ---- ★★★ 種類別の内訳（`SPEC-phase5.md` §5.10.11）------------------------
    //
    // **段ごとではなく「何に時間を使っているか」で分けます。**
    // **単価は `BENCH.md`（同じ機械の別の測定。由来は「推測」）と、
    // この実行で測った確保の単価を使います。**
    //
    // **和が 100% になることを確かめます**（`CLAUDE.md`）。**残りは「未計上」に出ます。**
    {
        // **確保の単価をこの機械で測ります**（`std::map` の節点に近い 64 バイトで）。
        const std::size_t kProbe = 4000000;
        std::vector<void*> keep(1024, nullptr);
        const std::clock_t pc0 = std::clock();
        for (std::size_t i = 0; i < kProbe; ++i) {
            void* q = std::malloc(64);
            std::free(keep[i & 1023]);
            keep[i & 1023] = q;
        }
        const double alloc_ns1 =
            1e9 * static_cast<double>(std::clock() - pc0) / CLOCKS_PER_SEC / kProbe;
        for (void* q : keep) std::free(q);
        // **★ 8 スレッドでも測ります。** 中核は並列なので確保器の競合が入ります。
        // **単一スレッドの温まった値は【下限】です。**
        double alloc_ns8 = alloc_ns1;
        {
            const std::size_t per = kProbe / g_threads;
            std::vector<std::thread> th;
            const std::clock_t q0 = std::clock();
            for (unsigned t = 0; t < g_threads; ++t) {
                th.emplace_back([per] {
                    std::vector<void*> k(1024, nullptr);
                    for (std::size_t i = 0; i < per; ++i) {
                        void* q = std::malloc(64);
                        std::free(k[i & 1023]);
                        k[i & 1023] = q;
                    }
                    for (void* q : k) std::free(q);
                });
            }
            for (auto& t : th) t.join();
            alloc_ns8 = 1e9 * static_cast<double>(std::clock() - q0) / CLOCKS_PER_SEC /
                        static_cast<double>(per * g_threads);
        }
        const double alloc_ns = alloc_ns8;
        std::printf(
            "\n**確保の単価**（この機械で実測。64 バイトの確保 + 解放）: "
            "**1 スレッド %.1f ns / %u スレッド %.1f ns**（後者を使います）\n",
            alloc_ns1, g_threads, alloc_ns8);

        std::printf(
            "\n| 段 | 中核の CPU | 述語 | 確保 | **点キャッシュ** | **未計上** | 述語 %% | "
            "確保 %% | **点キャッシュ %%** | **未計上 %%** |\n");
        std::printf("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
        for (const Row& r : g_rows) {
            const std::uint64_t sh =
                r.st.side_disp1 + r.st.side_disp2 + r.st.side_disp3 + r.st.side_disp4;
            const double pred_s =
                static_cast<double>(sh) * 7.80e-9 + static_cast<double>(r.side_ip) * 2.35e-9 +
                static_cast<double>(r.st.intersect3_arrange + r.st.intersect3_classify) *
                    261.76e-9 +
                static_cast<double>(r.cmp_h) * 11.25e-9;
            const double alloc_s = static_cast<double>(r.allocs) * alloc_ns * 1e-9;
            const double map_s =
                static_cast<double>(r.st.cache_hits + r.st.cache_misses) * map_ns * 1e-9;
            const double rest = r.cpu_s - pred_s - alloc_s - map_s;
            std::printf(
                "| %s | %.3f s | %.3f s | %.3f s | **%.3f s** | **%.3f s** | %.1f%% | "
                "%.1f%% | **%.1f%%** | **%.1f%%** |\n",
                r.name.c_str(), r.cpu_s, pred_s, alloc_s, map_s, rest,
                r.cpu_s == 0 ? 0.0 : 100.0 * pred_s / r.cpu_s,
                r.cpu_s == 0 ? 0.0 : 100.0 * alloc_s / r.cpu_s,
                r.cpu_s == 0 ? 0.0 : 100.0 * map_s / r.cpu_s,
                r.cpu_s == 0 ? 0.0 : 100.0 * rest / r.cpu_s);
        }
        std::printf(
            "\n| 段 | `side`(点) | `side`(整数) | `intersect3` | `cmp_h` | "
            "**確保の回数** | 確保の量 |\n");
        std::printf("|---|---:|---:|---:|---:|---:|---:|\n");
        for (const Row& r : g_rows) {
            const std::uint64_t sh =
                r.st.side_disp1 + r.st.side_disp2 + r.st.side_disp3 + r.st.side_disp4;
            std::printf(
                "| %s | %llu | %llu | %llu | %llu | **%llu** | %.1f MB |\n", r.name.c_str(),
                static_cast<unsigned long long>(sh), static_cast<unsigned long long>(r.side_ip),
                static_cast<unsigned long long>(r.st.intersect3_arrange + r.st.intersect3_classify),
                static_cast<unsigned long long>(r.cmp_h), static_cast<unsigned long long>(r.allocs),
                static_cast<double>(r.alloc_bytes) / 1048576.0);
        }
    }

    // ---- ★★ 段の内訳（CPU 時間ではなく壁時計。既にある計器）--------------------
    //
    // **種類別で 80% が未計上だったので、【どの段に】あるかを見ます。**
    std::printf("\n| 段 | 前処理 | 葉の列挙 | **arrange** | **縫合** | 分類 | 合計 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const double t =
            r.st.ms_prepare + r.st.ms_leaves + r.st.ms_arrange + r.st.ms_stitch + r.st.ms_classify;
        std::printf("| %s | %.1f%% | %.1f%% | **%.1f%%** | **%.1f%%** | %.1f%% | %.3f s |\n",
                    r.name.c_str(), t == 0 ? 0.0 : 100.0 * r.st.ms_prepare / t,
                    t == 0 ? 0.0 : 100.0 * r.st.ms_leaves / t,
                    t == 0 ? 0.0 : 100.0 * r.st.ms_arrange / t,
                    t == 0 ? 0.0 : 100.0 * r.st.ms_stitch / t,
                    t == 0 ? 0.0 : 100.0 * r.st.ms_classify / t, t / 1000.0);
    }
    std::printf(
        "\n| 段 | arrange の内訳: 収集 | 存在判定 | 準備 | **断片の生成** | 共平面 | "
        "縫合（葉ごと） |\n");
    std::printf("|---|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const double t = r.st.ms_arr_gather + r.st.ms_arr_present + r.st.ms_arr_prep +
                         r.st.ms_arr_frag + r.st.ms_arr_coplanar + r.st.ms_arr_stitch;
        std::printf("| %s | %.1f%% | %.1f%% | %.1f%% | **%.1f%%** | %.1f%% | %.1f%% |\n",
                    r.name.c_str(), t == 0 ? 0.0 : 100.0 * r.st.ms_arr_gather / t,
                    t == 0 ? 0.0 : 100.0 * r.st.ms_arr_present / t,
                    t == 0 ? 0.0 : 100.0 * r.st.ms_arr_prep / t,
                    t == 0 ? 0.0 : 100.0 * r.st.ms_arr_frag / t,
                    t == 0 ? 0.0 : 100.0 * r.st.ms_arr_coplanar / t,
                    t == 0 ? 0.0 : 100.0 * r.st.ms_arr_stitch / t);
    }

    // ---- ★ 分類の内訳（`SPEC-phase5.md` §5.10.14.7。まず刻むだけ）------------------------
    //
    // **CPU（全スレッドの和）で 4
    // 区分。壁時計は「並列の領域ループ」と「逐次（順序の構築と結合）」。**
    std::printf(
        "\n| 段 | 分類の内訳（CPU）: 準備 | 代表点 | **レイ** | 判定 + 出力 | CPU の和 | "
        "並列部の壁時計 | 逐次部の壁時計 | 並列効率 | 領域 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const double c = r.st.ms_cl_prep + r.st.ms_cl_rep + r.st.ms_cl_ray + r.st.ms_cl_out;
        const double seq = r.st.ms_classify - r.st.ms_cl_par_wall;
        std::printf(
            "| %s | %.1f%% | %.1f%% | **%.1f%%** | %.1f%% | %.3f s | %.3f s | %.3f s | %.2f | %zu "
            "|\n",
            r.name.c_str(), c == 0 ? 0.0 : 100.0 * r.st.ms_cl_prep / c,
            c == 0 ? 0.0 : 100.0 * r.st.ms_cl_rep / c, c == 0 ? 0.0 : 100.0 * r.st.ms_cl_ray / c,
            c == 0 ? 0.0 : 100.0 * r.st.ms_cl_out / c, c / 1000.0, r.st.ms_cl_par_wall / 1000.0,
            seq / 1000.0,
            r.st.ms_cl_par_wall == 0 ? 0.0 : c / (r.st.ms_cl_par_wall * (double)g_threads),
            r.st.regions);
    }

    // ---- ★ レイの候補の漏斗（§5.10.14.7。索引が返す候補のうち何割を使うか）----------------
    std::printf(
        "\n| 段 | レイ | 索引の候補 / レイ | 投影 AABB を通過 / レイ | 前方 / レイ | "
        "**厳密判定に回った / レイ** | "
        "寄与 / レイ | 安い前判定の回数 / レイ |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const double n = r.st.raycasts == 0 ? 1.0 : (double)r.st.raycasts;
        std::printf("| %s | %zu | %.1f | %.2f | %.2f | **%.2f** | %.2f | %.1f |\n", r.name.c_str(),
                    r.st.raycasts, r.st.ray_tri_tests / n, r.st.ray_tri_aabb / n,
                    r.st.ray_tri_fwd / n, r.st.ray_tri_kept / n, r.st.ray_tri_hits / n,
                    r.st.ray_cheap_tests / n);
    }

    // ---- ★ 索引の粒度（`SPEC-phase5.md` §5.10.14.11。候補はどの段から来るか）------------
    //
    // **段 0 が最も細かい。粗い段の項目は、その粗いセルに落ちるすべてのレイが見ます。**
    {
        std::printf(
            "\n| 段（演算） | 候補の総数 | 索引の段ごとの候補（細 → 粗）と割合 "
            "|\n|---|---:|---|\n");
        for (const Row& r : g_rows) {
            std::size_t tot = 0;
            for (int l = 0; l < 12; ++l) tot += r.st.ray_cand_level[l];
            std::printf("| %s | %zu | ", r.name.c_str(), tot);
            for (std::size_t l = 0; l < r.st.ray_levels_max && l < 12; ++l) {
                std::printf("%s%zu（%.1f%%）", l ? " / " : "", r.st.ray_cand_level[l],
                            tot == 0 ? 0.0 : 100.0 * (double)r.st.ray_cand_level[l] / (double)tot);
            }
            std::printf(" |\n");
        }
    }

    // ---- ★ 縫合の内訳と、並列化の前提の確認（`SPEC-phase5.md` §5.10.13）-----------
    std::printf(
        "\n| 段 | 縫合の内訳（壁時計）: 構成点 + 3 つ組の表 | 整列 | 番号付け | **仕分け** | "
        "縫合の合計 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const double t =
            r.st.ms_st_points + r.st.ms_st_sort + r.st.ms_st_remap + r.st.ms_st_regions;
        std::printf(
            "| %s | %.1f%% | %.1f%% | %.1f%% | **%.1f%%** | %.3f s（段の計 %.3f s） |\n",
            r.name.c_str(), t == 0 ? 0.0 : 100.0 * r.st.ms_st_points / t,
            t == 0 ? 0.0 : 100.0 * r.st.ms_st_sort / t, t == 0 ? 0.0 : 100.0 * r.st.ms_st_remap / t,
            t == 0 ? 0.0 : 100.0 * r.st.ms_st_regions / t, t / 1000.0, r.st.ms_stitch / 1000.0);
    }
    if (g_rows[0].st.leaf_adj_pairs > 0) {
        std::printf(
            "\n| 段 | 構成点 | 葉をまたぐ点 | うち深度が混ざる | **境界上の点** | 葉 | 深度 | "
            "接する葉の対 | うち深度が違う | 1 葉あたり最大 | 併合群 | 最大広がり | "
            "またぐのに境界に無い |\n");
        std::printf("|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|\n");
        for (const Row& r : g_rows) {
            std::printf(
                "| %s | %zu | **%zu**（%.2f%%） | **%zu** | **%zu**（%.2f%%） | %zu | %u〜%u | %zu "
                "| "
                "%zu | %zu | %zu | %zu | **%zu** |\n",
                r.name.c_str(), r.st.constructed_points, r.st.pt_multi_leaf,
                r.st.constructed_points == 0
                    ? 0.0
                    : 100.0 * (double)r.st.pt_multi_leaf / (double)r.st.constructed_points,
                r.st.pt_mixed_depth, r.st.pt_on_boundary,
                r.st.constructed_points == 0
                    ? 0.0
                    : 100.0 * (double)r.st.pt_on_boundary / (double)r.st.constructed_points,
                r.st.leaf_nonempty, r.st.leaf_depth_min, r.st.leaf_depth_max, r.st.leaf_adj_pairs,
                r.st.leaf_adj_mixed, r.st.leaf_adj_max, r.st.merge_groups, r.st.max_merge_span,
                r.st.pt_multi_not_boundary);
        }
    }

    // ---- ★★★ 構成点キャッシュの探索（`split_fragment` が頂点ごとに呼びます）--------
    //
    // **`split_fragment` は `s[i] = side(qp, fragment_vertex(t, f, i, cache))` です。**
    // **`fragment_vertex` は `PointCache`（`std::map<array<PlaneId,3>, HPointD>`）の探索。**
    // **`side` の回数と同じだけ探索が走ります。**
    std::printf("\n| 段 | 命中 | 失敗 | **探索の合計** | 命中率 | 登録 | 量 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const std::size_t look = r.st.cache_hits + r.st.cache_misses;
        std::printf("| %s | %zu | %zu | **%zu** | %.1f%% | %zu | %.1f MB |\n", r.name.c_str(),
                    r.st.cache_hits, r.st.cache_misses, look,
                    look == 0
                        ? 0.0
                        : 100.0 * static_cast<double>(r.st.cache_hits) / static_cast<double>(look),
                    r.st.cache_entries, static_cast<double>(r.st.cache_bytes) / 1048576.0);
    }

    // ---- ★★★ 確保の発生箇所ごとの内訳（§5.10.12.4。残る 93.4% を刻む）--------------
    //
    // **回数だけでなく【CPU 時間あたりの密度】も出します**（仕様側の条件）。
    // **和が 100% になることを、印なし（その他）を含めて確かめます。**
    if (g_rows[0].allocs > 0) {
        for (std::size_t ri = 0; ri < g_rows.size(); ++ri) {
            const Row& r = g_rows[ri];
            std::uint64_t sum = 0;
            for (int k = 0; k < 20; ++k) sum += r.tag_count[k];
            if (sum == 0) continue;
            std::printf("\n#### %s — 確保の発生箇所（全 %llu 回、CPU %.3f 秒）\n\n", r.name.c_str(),
                        (unsigned long long)r.allocs, r.cpu_s);
            std::printf("| 発生箇所 | 回数 | **割合** | 量 | **密度（回 / CPU 秒）** |\n");
            std::printf("|---|---:|---:|---:|---:|\n");
            for (int k = 0; k < 20; ++k) {
                if (r.tag_count[k] == 0 && k != 0) continue;
                std::printf("| %s | %llu | **%.1f%%** | %.1f MB | %.0f |\n", kricount::kTagName[k],
                            (unsigned long long)r.tag_count[k],
                            100.0 * (double)r.tag_count[k] / (double)r.allocs,
                            (double)r.tag_bytes[k] / 1048576.0,
                            r.cpu_s == 0 ? 0.0 : (double)r.tag_count[k] / r.cpu_s);
            }
            std::printf("| **和** | **%llu** | **%.1f%%** | | |\n", (unsigned long long)sum,
                        100.0 * (double)sum / (double)r.allocs);
        }
    }

    // ---- ★★ 断片の生成の内訳と、`edge` の small-array の見積もり（§5.10.12.4）------
    if (g_rows[0].split_calls > 0) {
        std::printf(
            "\n| 段 | `split` 呼び出し | 早期 return | 2 つに切った | `make` | 新 `edge` | "
            "`clip` | **`split` 由来の確保** | 全確保 | **割合** |\n");
        std::printf("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
        for (const Row& r : g_rows) {
            // `s` は呼び出しごと、`keep` は make ごと、新 edge は make_ok、複製は早期 return ごと
            const std::uint64_t due = r.split_calls + r.make_calls + r.make_ok + r.split_early;
            std::printf(
                "| %s | %llu | %llu | %llu | %llu | %llu | %llu | **%llu** | %llu | "
                "**%.1f%%** |\n",
                r.name.c_str(), (unsigned long long)r.split_calls,
                (unsigned long long)r.split_early, (unsigned long long)r.split_both,
                (unsigned long long)r.make_calls, (unsigned long long)r.make_ok,
                (unsigned long long)r.clip_calls, (unsigned long long)due,
                (unsigned long long)r.allocs,
                r.allocs == 0 ? 0.0 : 100.0 * (double)due / (double)r.allocs);
        }
        std::printf("\n| 段 | 3 | 4 | 5 | 6 | 7 | 8 | ≥9 | **≤4 累積** | **≤8 累積** | 最大 |\n");
        std::printf("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
        for (const Row& r : g_rows) {
            std::uint64_t tot = 0;
            for (int k = 0; k < 10; ++k) tot += r.edge_hist[k];
            const double t = tot == 0 ? 1.0 : (double)tot;
            const double c4 = (r.edge_hist[3] + r.edge_hist[4]) / t;
            double c8 = 0;
            for (int k = 3; k <= 8; ++k) c8 += r.edge_hist[k];
            c8 /= t;
            std::printf(
                "| %s | %.1f%% | %.1f%% | %.1f%% | %.1f%% | %.1f%% | %.1f%% | %.2f%% | "
                "**%.1f%%** | **%.2f%%** | %zu |\n",
                r.name.c_str(), 100.0 * r.edge_hist[3] / t, 100.0 * r.edge_hist[4] / t,
                100.0 * r.edge_hist[5] / t, 100.0 * r.edge_hist[6] / t, 100.0 * r.edge_hist[7] / t,
                100.0 * r.edge_hist[8] / t, 100.0 * r.edge_hist[9] / t, 100.0 * c4, 100.0 * c8,
                r.st.frag_edges_max);
        }
    }

    // ---- ★ 分類の費用（依頼 2 の材料）--------------------------------------
    //
    // **参照点の伝播で置き換えたいのは、この「領域ごとの大域レイキャスト」です。**
    // **局所トレースの費用の上限は「葉の中の多角形数」なので、両方を並べます。**
    std::printf(
        "\n| 段 | 領域（レイ） | 隅のレイ | 三角形検査 | /レイ | 葉あたり多角形 | "
        "**中核の秒** | 分類の割合 |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|---:|\n");
    for (const Row& r : g_rows) {
        const double total = r.st.ms_prepare + r.st.ms_arrange + r.st.ms_stitch + r.st.ms_classify;
        std::printf("| %s | %zu | %zu | %zu | %.1f | %.1f | **%.3f** | %.1f%% |\n", r.name.c_str(),
                    r.st.raycasts, r.st.early_out_raycasts, r.st.ray_tri_tests,
                    r.st.raycasts == 0 ? 0.0
                                       : static_cast<double>(r.st.ray_tri_tests) /
                                             static_cast<double>(r.st.raycasts),
                    r.st.leaf_nonempty == 0 ? 0.0
                                            : static_cast<double>(r.st.raw_fragments) /
                                                  static_cast<double>(r.st.leaf_nonempty),
                    total / 1000.0, total == 0 ? 0.0 : 100.0 * r.st.ms_classify / total);
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
    // **★ 鍵の突き合わせ（`verify_region_key`）は既定で切ります**（第 7 引数で 1 にすると入る）。
    // **本番では偽なので、確保の内訳を測るときに入れると仕分けの段が水増しされます**
    // （実際に踏みました。`old_regions` の分が 14.5% に乗っていました）。
    const bool verify_rk = (argc > 7) ? (std::atoi(argv[7]) != 0) : false;
    // **★ 器の使い回し（G1〜G3）のビット集合**（第 8 引数。既定 7 = 本番と同じ）。
    // **確保の回数を項目ごとに測るには 0 / 1 / 2 / 4 / 7 で回します。**
    const unsigned alloc_reuse = (argc > 8) ? static_cast<unsigned>(std::atoi(argv[8])) : 7u;

    // **★ 設定と対象を最初に出します**（`CLAUDE.md`）
    std::printf("\n## 連鎖の断片数（`SPEC-phase5.md` §5.10.5.7 / 依頼 3）\n\n");
    std::printf("| 設定 | 値 |\n|---|---|\n");
    std::printf("| 対象 | A=`%s` B=`%s` D=`%s` |\n", ida.c_str(), idb.c_str(), idd.c_str());
    std::printf("| 深度 | %u（%s） |\n", depth, adaptive ? "**適応**" : "**固定**");
    std::printf("| 鍵の突き合わせ | %s |\n",
                verify_rk ? "**入れる**（本番は偽）" : "切る（本番と同じ）");
    std::printf("| スレッド | %u |\n", nthreads);
    std::printf("| 器の使い回し（alloc_reuse） | %u（G1=%d G2=%d G3=%d） |\n", alloc_reuse,
                (alloc_reuse & 1u) != 0, (alloc_reuse & 2u) != 0, (alloc_reuse & 4u) != 0);
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

    g_threads = nthreads;
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
    o.verify_region_key = verify_rk;
    o.alloc_reuse = alloc_reuse;
    o.measure_classify = true;   // 分類の内訳（§5.10.14.7）
    o.record_ray_levels = true;  // 索引の粒度（§5.10.14.11）
    // **縫合の前提の確認**（第 9 引数。既定 0。$O(L^2)$ の葉の対の数え上げを含む）
    o.measure_stitch = (argc > 9) && (std::atoi(argv[9]) != 0);

    const csg::PolySoup A = csg::from_mesh(qa.mesh);
    const csg::PolySoup B = csg::from_mesh(qb.mesh);
    const csg::PolySoup D = csg::from_mesh(qd.mesh);

    // ---- ★ 索引の段ごとの項目数（§5.10.14.11。§11.3 の形。模型・軸ごと）------------------
    {
        // **模型ごと・軸ごとの段の項目数**（§11.3 の形。三角形は自分の大きさに合った段に入る）
        std::printf(
            "\n| 模型 | 三角形 | 軸 | 段の数 | 段ごとの項目数（細 → 粗） | 最多の段 "
            "|\n|---|---:|---|---:|---|---:|\n");
        const struct {
            const char* name;
            const mesh::TriMesh* m;
        } models[3] = {{ida.c_str(), &qa.mesh}, {idb.c_str(), &qb.mesh}, {idd.c_str(), &qd.mesh}};
        for (const auto& mm : models) {
            for (int ax = 0; ax < 3; ++ax) {
                csg::RayIndex ix;
                ix.build(*mm.m, static_cast<geom::Axis>(ax));
                std::size_t best = 0;
                std::printf("| `%s` | %zu | %c | %zu | ", mm.name, mm.m->triangles.size(),
                            "XYZ"[ax], ix.levels());
                for (std::size_t l = 0; l < ix.levels(); ++l) {
                    if (ix.items_at(l) > ix.items_at(best)) best = l;
                    std::printf("%s%zu", l ? " / " : "", ix.items_at(l));
                }
                std::printf(" | **%zu** |\n", best);
            }
        }
    }

    // ---- ★★ 構成点キャッシュの A/B（`SPEC-phase5.md` §5.10.12）------------------
    //
    // **同一実行の中で `std::map` と開番地法を比べます。**
    // **時間帯の交絡（`BENCH.md` の 1.5〜1.6 倍）を受けません。**
    {
        std::printf("\n### 構成点キャッシュの A/B（同一実行）\n\n");
        std::printf("| 実装 | 壁時計 | CPU 時間 | 探索 | 命中率 | 出力ハッシュ |\n");
        std::printf("|---|---:|---:|---:|---:|---|\n");
        unsigned long long h[2] = {0, 0};
        double cpu[2] = {0, 0};
        for (int use_map = 1; use_map >= 0; --use_map) {
            csg::BoolOptions ab = o;
            ab.point_cache_map = (use_map != 0);
            csg::BoolStats st;
            const auto t0 = std::chrono::steady_clock::now();
            const std::clock_t c0 = std::clock();
            const csg::PolySoup s2 = csg::boolean(A, D, csg::BoolOp::Difference, ab, &st);
            const double cs = static_cast<double>(std::clock() - c0) / CLOCKS_PER_SEC;
            const double ws =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            csg::ToMeshOptions tm2;
            tm2.split_contacts = true;
            h[use_map] = hash_mesh(csg::to_mesh(s2, tm2));
            cpu[use_map] = cs;
            const std::size_t look = st.cache_hits + st.cache_misses;
            std::printf(
                "| %s | %.3f s | **%.3f s** | %zu | %.1f%% | `%016llx` |\n",
                use_map ? "`std::map`（従来）" : "**開番地法（既定）**", ws, cs, look,
                look == 0 ? 0.0
                          : 100.0 * static_cast<double>(st.cache_hits) / static_cast<double>(look),
                h[use_map]);
        }
        std::printf("\n**出力**: %s / **CPU 時間の比**: **%.2f 倍**\n",
                    h[0] == h[1] ? "**バイト一致**" : "**★ 食い違い（重大）**",
                    cpu[0] == 0 ? 0.0 : cpu[1] / cpu[0]);
    }

    // ---- ★★ 縫合の A/B（`SPEC-phase5.md` §5.10.13.3。葉ごと + 境界の類だけ大域）--------
    //
    // **従来（大域の表と整列、逐次）と新（葉ごと + 境界だけ大域）を同一実行で比べます。**
    // **領域の順序が変わるので出力のバイトは違って当然で、一致は【順序を除いた鍵】で見ます。**
    // **決定性はスレッド数 1 と 8 のバイト一致で見ます。**
    {
        std::printf("\n### 縫合の A/B（同一実行）\n\n");
        std::printf(
            "| 実装 | スレッド | 壁時計 | CPU 時間 | 縫合（壁） | 葉の中の縫合（CPU） | 境界の類 / "
            "類 | "
            "出力ハッシュ |\n|---|---:|---:|---:|---:|---:|---:|---|\n");
        std::vector<std::array<std::uint32_t, 3>> gk[3];
        unsigned long long h[3] = {0, 0, 0};
        double wall[3] = {0, 0, 0}, stw[3] = {0, 0, 0};
        for (int cfg = 0; cfg < 3; ++cfg) {
            csg::BoolOptions ab = o;
            ab.stitch_parallel = (cfg != 0);
            if (cfg == 2) {
                ab.threads = 1;
                ab.pool = nullptr;
            }
            csg::BoolStats st;
            const auto t0 = std::chrono::steady_clock::now();
            const std::clock_t c0 = std::clock();
            const csg::PolySoup s2 = csg::boolean(A, D, csg::BoolOp::Difference, ab, &st);
            const double cs = static_cast<double>(std::clock() - c0) / CLOCKS_PER_SEC;
            const double ws =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            csg::ToMeshOptions tm2;
            tm2.split_contacts = true;
            const csg::SoupMesh m = csg::to_mesh(s2, tm2);
            h[cfg] = hash_mesh(m);
            gk[cfg] = geometric_key(m);
            wall[cfg] = ws;
            stw[cfg] = st.ms_stitch / 1000.0;
            std::printf(
                "| %s | %u | %.3f s | %.3f s | **%.3f s** | %.3f s | %zu / %zu | `%016llx` |\n",
                cfg == 0 ? "従来（大域の表と整列）" : "**葉ごと + 境界だけ大域**",
                cfg == 2 ? 1u : nthreads, ws, cs, st.ms_stitch / 1000.0, st.ms_arr_stitch / 1000.0,
                st.stitch_boundary_classes, st.stitch_classes, h[cfg]);
        }
        std::printf(
            "\n**順序を除いた鍵**: %s / **スレッド数 1 と 8**: %s / **バイト**: %s / "
            "**縫合（壁）の比（従来 ÷ 新）**: **%.2f 倍** / 全体（壁）: %.2f 倍\n",
            gk[0] == gk[1] ? "**一致**" : "**★ 食い違い（重大）**",
            h[1] == h[2] ? "**バイト一致**" : "**★ 食い違い（決定性が壊れている）**",
            h[0] == h[1] ? "一致（順序も同じ）" : "違う（領域の順序が変わる。想定どおり）",
            stw[1] == 0 ? 0.0 : stw[0] / stw[1], wall[1] == 0 ? 0.0 : wall[0] / wall[1]);
    }

    // ---- ★★ 断片の切断の A/B（`SPEC-phase5.md` §5.10.12.4。早期 return を移動に）------
    //
    // **同一実行の中で、従来（`SplitResult` を値で返す）と移動版を比べます。**
    {
        std::printf("\n### 断片の切断の A/B（同一実行）\n\n");
        std::printf(
            "| 実装 | 壁時計 | CPU 時間 | 確保 | 出力ハッシュ |\n|---|---:|---:|---:|---|\n");
        unsigned long long h[2] = {0, 0};
        double cpu[2] = {0, 0};
        std::uint64_t al[2] = {0, 0};
        for (int legacy = 1; legacy >= 0; --legacy) {
            csg::BoolOptions ab = o;
            ab.split_legacy = (legacy != 0);
            csg::BoolStats st;
            const std::uint64_t a0 = kricount::alloc_count.load(std::memory_order_relaxed);
            kricount::enabled.store(true, std::memory_order_relaxed);
            const auto t0 = std::chrono::steady_clock::now();
            const std::clock_t c0 = std::clock();
            const csg::PolySoup s2 = csg::boolean(A, D, csg::BoolOp::Difference, ab, &st);
            const double cs = static_cast<double>(std::clock() - c0) / CLOCKS_PER_SEC;
            const double ws =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            kricount::enabled.store(false, std::memory_order_relaxed);
            al[legacy] = kricount::alloc_count.load(std::memory_order_relaxed) - a0;
            csg::ToMeshOptions tm2;
            tm2.split_contacts = true;
            h[legacy] = hash_mesh(csg::to_mesh(s2, tm2));
            cpu[legacy] = cs;
            std::printf("| %s | %.3f s | **%.3f s** | %llu | `%016llx` |\n",
                        legacy ? "従来（`SplitResult` を値で）" : "**移動（既定）**", ws, cs,
                        static_cast<unsigned long long>(al[legacy]), h[legacy]);
        }
        std::printf("\n**出力**: %s / **CPU 時間の比**: **%.2f 倍** / 確保: %llu → %llu\n",
                    h[0] == h[1] ? "**バイト一致**" : "**★ 食い違い（重大）**",
                    cpu[0] == 0 ? 0.0 : cpu[1] / cpu[0], (unsigned long long)al[1],
                    (unsigned long long)al[0]);
    }

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
    csg::PolySoup s1, s2;
    run_stage("連鎖 1 段目 `A∪B`", A, B, csg::BoolOp::Union, o, &s1);
    run_stage("連鎖 2 段目 `(A∪B)＼D`", s1, D, csg::BoolOp::Difference, o, &s2);

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
        // **連鎖の 2 段もバイト一致で見ます**（G1 の「出力は不変のはず」を 4 段すべてで確かめる）
        std::printf("**★ ハッシュ（連鎖）**: `A∪B` = `%016llx` / `(A∪B)＼D` = `%016llx`\n",
                    hash_mesh(csg::to_mesh(s1, tm)), hash_mesh(csg::to_mesh(s2, tm)));
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

// Krisite — side : intersect3 の呼び出し比率（SPEC-phase1.md §12）
//
// **最適化はしません。比率だけを測ります。**
//
// IMPL-phase0 §7 の問い: `intersect3` は 262 ns で `side`（11 ns）より 24 倍重いが、
// 「構成点は一度作って何度も `side` にかける」なら支配的とは限らない。
// **呼び出し比率を知らずに重い関数を最適化しても無駄になります。**
//
// この実行ファイルだけ `KRISITE_COUNT_PREDICATES` を定義してビルドします。
// 既定ビルドには計数のコストを持ち込みません。
#include <chrono>
#include <cstdio>

#include "krisite/csg/boolean.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite::csg;

namespace {

const char* op_name(BoolOp op) {
    return op == BoolOp::Union ? "∪" : op == BoolOp::Intersection ? "∩" : "\\";
}

}  // namespace

int main() {
#if !defined(KRISITE_COUNT_PREDICATES)
    std::printf("\n  KRISITE_COUNT_PREDICATES が定義されていません。計数できません。\n");
    KRI_CHECK_MSG(false, "計測用の実行ファイルは KRISITE_COUNT_PREDICATES 付きでビルドすること");
    return kritest::finish("csg/call_ratio");
#else
    std::printf("\n  §12 side : intersect3 の呼び出し比率\n");
    std::printf("    %-5s %-4s %-3s %-12s %-12s %-8s %s\n", "ケース", "深度", "演算", "side",
                "intersect3", "比", "断片");

    std::uint64_t tot_side = 0, tot_i3 = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const auto a = c.make_a(), b = c.make_b();
        for (unsigned d = 0; d <= 3; ++d) {
            for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
                BoolStats st;
                // **既定値に依存しないこと**（§9.4 の CI ジョブで既定が反転する）
                boolean_op(a, b, op, kritest::phase1_options(d), &st);
                KRI_CHECK_MSG(st.side_calls > 0,
                              "side が 1 回も呼ばれていない（計数が効いていない）");
                KRI_CHECK_MSG(st.intersect3_calls > 0, "intersect3 が 1 回も呼ばれていない");
                tot_side += st.side_calls;
                tot_i3 += st.intersect3_calls;
                if (d == 3 && op == BoolOp::Union) {
                    std::printf("    %-6s %-5u %-4s %-13llu %-13llu %-9.1f %zu\n", c.id, d,
                                op_name(op), (unsigned long long)st.side_calls,
                                (unsigned long long)st.intersect3_calls,
                                static_cast<double>(st.side_calls) /
                                    static_cast<double>(st.intersect3_calls),
                                st.fragments);
                }
            }
        }
    }

    const double ratio = static_cast<double>(tot_side) / static_cast<double>(tot_i3);
    std::printf("\n    合計: side = %llu, intersect3 = %llu, 比 = %.2f : 1\n",
                (unsigned long long)tot_side, (unsigned long long)tot_i3, ratio);

    // Phase 0 の実測（pred_bench）: side 11.25 ns / intersect3 262 ns
    const double t_side = 11.25, t_i3 = 262.0;
    const double share =
        (static_cast<double>(tot_i3) * t_i3) /
        (static_cast<double>(tot_i3) * t_i3 + static_cast<double>(tot_side) * t_side) * 100.0;
    std::printf("    Phase 0 の実測時間で重み付けすると intersect3 が全体の %.1f%%\n", share);
    std::printf("    （side 11.25 ns / intersect3 262 ns。BENCH.md）\n");

    // **比率が 1 未満なら、intersect3 のほうが呼ばれていることになります。**
    // その場合 IMPL-phase0 §7 の前提（構成点を一度作って何度も side にかける）が
    // 崩れているので、報告対象です。
    KRI_CHECK_MSG(ratio > 1.0,
                  "side の呼び出しが intersect3 より少ない。IMPL-phase0 §7 の前提が崩れています");

    // ---- §4 の構成点の保持で比がどう変わるか -------------------------------
    //
    // **Phase 1 の前提はここで作り直されます。** 単価ではなく回数を減らす手です。
    std::uint64_t c_side = 0, c_i3 = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const auto a = c.make_a(), b = c.make_b();
        for (unsigned d = 0; d <= 3; ++d) {
            for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
                BoolStats st;
                BoolOptions opt = kritest::phase1_options(d);
                opt.cache_points = true;
                boolean_op(a, b, op, opt, &st);
                c_side += st.side_calls;
                c_i3 += st.intersect3_calls;
            }
        }
    }
    const double cratio = static_cast<double>(c_side) / static_cast<double>(c_i3);
    const double cshare =
        (static_cast<double>(c_i3) * t_i3) /
        (static_cast<double>(c_i3) * t_i3 + static_cast<double>(c_side) * t_side) * 100.0;
    std::printf("\n  §4 構成点の保持あり\n");
    std::printf("    合計: side = %llu, intersect3 = %llu, 比 = %.2f : 1\n",
                (unsigned long long)c_side, (unsigned long long)c_i3, cratio);
    std::printf("    intersect3 の呼び出しは %.1f%% に減り、時間の占有は %.1f%% → %.1f%%\n",
                100.0 * static_cast<double>(c_i3) / static_cast<double>(tot_i3), share, cshare);

    // **回数が減っていること。** 減らないならメモ化が効いていません
    KRI_CHECK_MSG(c_i3 < tot_i3, "キャッシュを入れても intersect3 の回数が減っていない");
    KRI_CHECK_MSG(cratio > ratio, "比が改善していない。§4.1 の構図が変わっていません");

    // ---- §11 の【時間】★ Phase 1 で測っていないもの -------------------------
    //
    // **ばらつきに注意してください**（`BENCH.md` は同一条件でも ±15% 動くと記録）。
    // **1 割前後の差に意味を持たせないこと。** ここでは合否に使いません。
    auto measure = [&](bool cache, bool adaptive, bool early) {
        const auto t0 = std::chrono::steady_clock::now();
        for (const kritest::Case& c : kritest::corpus()) {
            const auto a = c.make_a(), b = c.make_b();
            for (BoolOp op : {BoolOp::Union, BoolOp::Intersection, BoolOp::Difference}) {
                BoolOptions opt = kritest::phase1_options(3);
                opt.cache_points = cache;
                opt.adaptive = adaptive;
                opt.early_out = early;
                boolean_op(a, b, op, opt, nullptr);
            }
        }
        const auto t1 = std::chrono::steady_clock::now();
        return std::chrono::duration<double, std::milli>(t1 - t0).count();
    };
    std::printf("\n  §11 実時間（全 18 ケース × 3 演算 × 深度 3。**合否には使いません**）\n");
    std::printf("    %-34s %8.1f ms\n", "固定深度", measure(false, false, false));
    std::printf("    %-34s %8.1f ms\n", "固定深度 + 構成点の保持", measure(true, false, false));
    std::printf("    %-34s %8.1f ms\n", "適応分割", measure(false, true, false));
    std::printf("    %-34s %8.1f ms\n", "適応分割 + early-out", measure(false, true, true));
    std::printf("    %-34s %8.1f ms\n", "適応分割 + early-out + 保持", measure(true, true, true));
    std::printf("    **同一条件でも ±15%% 動きます。1 割前後の差に意味を持たせないこと。**\n");

    std::printf("\n");
    return kritest::finish("csg/call_ratio");
#endif
}

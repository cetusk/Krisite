// Krisite — side : intersect3 の呼び出し比率（SPEC-phase1.md §12）
//
// **最適化はしません。比率だけを測ります。**
//
// IMPL-phase0 §7 の問い: `intersect3` は 262 ns で `side`（11 ns）より 24 倍重いが、
// 「構成点は一度作って何度も `side` にかける」なら支配的とは限らない。
// **呼び出し比率を知らずに重い関数を最適化しても無駄になります**（CLAUDE.md「計測」）。
//
// この実行ファイルだけ `KRISITE_COUNT_PREDICATES` を定義してビルドします。
// 既定ビルドには計数のコストを持ち込みません。
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
                boolean_op(a, b, op, d, &st);
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
    std::printf("\n");
    return kritest::finish("csg/call_ratio");
#endif
}

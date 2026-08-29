// Krisite — 並列実行の最小の器（`SPEC-phase4.md` §2）
//
// **中央同期キューを 1 本のアトミックな添字で表します。**
//
// EMBER §4.5.4 の work-stealing は「各分割ステップで一方を再帰し、もう一方をキューへ。
// 待機中のスレッドがキューから次の部分問題を取る」形です。**Krisite は葉を先に
// 列挙する**（`octree::build_leaves`）ので、タスクは最初から全部そろっています。
// キューは固定長の配列で足り、取り出しはアトミックな添字の加算 1 回です。
//
// > **したがって burn-in がありません**（EMBER §5.3）。再帰でタスクを生む構造だと
// > 完全な並列性に達するまで時間がかかりますが、こちらは最初から全並列です。
// > **EMBER §6 の「静的な事前クリップ」は、Krisite では構造として既に入っています。**
//
// **決定性はこの器が持ちません。呼び出し側の責任です**（`SPEC-phase4.md` §4.2）。
// 各タスクは**自分の添字のスロットにだけ書き**、結合は**添字順**に行ってください。
#ifndef KRISITE_PAR_PARALLEL_FOR_HPP
#define KRISITE_PAR_PARALLEL_FOR_HPP

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <thread>
#include <vector>

namespace krisite::par {

/// このマシンで使えるスレッド数（取得できなければ 1）。
inline unsigned hardware_threads() noexcept {
    const unsigned n = std::thread::hardware_concurrency();
    return n == 0 ? 1u : n;
}

/// `fn(i, tid)` を $i = 0 \dots n-1$ について実行する。**順序は保証しません。**
///
/// `tid` は $0 \dots t-1$ で、**スレッド局所の器を引くための添字**です
/// （`SPEC-phase4.md` §1.1「可変な部分はスレッド局所に」）。
///
/// **`threads <= 1` ならスレッドを作りません。** 逐次経路がそのまま残るので、
/// §7.1 の「逐次実装との一致」が同一プロセスで比較できます。
template <class F>
inline void parallel_for(std::size_t n, unsigned threads, F&& fn) {
    if (threads <= 1 || n <= 1) {
        for (std::size_t i = 0; i < n; ++i) fn(i, 0u);
        return;
    }
    const auto t = static_cast<unsigned>(std::min<std::size_t>(threads, n));
    std::atomic<std::size_t> next{0};
    std::vector<std::thread> pool;
    pool.reserve(t);
    for (unsigned k = 0; k < t; ++k) {
        pool.emplace_back([&fn, &next, n, k] {
            for (;;) {
                // **`relaxed` で足ります。** 同期はタスク間ではなく join で取ります
                const std::size_t i = next.fetch_add(1, std::memory_order_relaxed);
                if (i >= n) break;
                fn(i, k);
            }
        });
    }
    for (std::thread& th : pool) th.join();
}

}  // namespace krisite::par

#endif  // KRISITE_PAR_PARALLEL_FOR_HPP

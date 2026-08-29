// Krisite — 永続スレッドプール（`SPEC-phase4.md` §2）
//
// **呼び出しのたびにスレッドを作ってはいけません。**
//
// `parallel_for` は 1 回ごとに生成・join するので、**タスクが小さいと生成コストが
// 支配します。** 実測: 8 スレッド × 132 回で **26.2 ms**（中核 64.7 ms の 40%）。
// **並列化した結果、並列化のコストが最大の項目になっていました。**
//
// ここでは待機しているスレッドを持ち回します。`run` は
//
//   1. 仕事を置いて世代を進める
//   2. **呼び出し元も 1 本として働く**（`tid = 0`）
//   3. **走り出したワーカーだけの完了を待つ**
//
// **「全員の完了」を待ってはいけません。** タスクが小さいと呼び出し元が全部さらって
// しまい、**何もしなかったワーカーが起きるのを待つだけ**になります。
// 実測: 出口を 8 スレッドで回すと `real 31.9s` に対し `user 1.8s`。**94% が待ち**でした。
//
// **決定性はこの器が持ちません。呼び出し側の責任です**（§4.2）。
// 各タスクは自分の添字のスロットにだけ書き、結合は添字順に行うこと。
//
// **`run` は入れ子にできません。** 所有スレッドからだけ呼んでください。
//
// ---
//
// ## ディスパッチの下限（`SPEC-phase4.md` §6.3）★
//
// **項目数が少ない段は、並列にすると遅くなります。**
//
//   索引       0.9 ms → 2.9 ms（**3 倍遅い**）
//   分裂・扇   2.8 ms → 3.8 ms（1.4 倍遅い）
//
// **1 項目あたりの仕事がディスパッチのコストを下回る**ためです。
// そこで `min_items` を下回る呼び出しは**その場で逐次実行**します。
//
// **項目数だけでは段を分けられません。** 1 項目あたりの仕事が段ごとに違うからです。
//
//   索引       8 項目で **損**（1 項目 1.7 µs、1 回 0.014 ms）
//   三角形化  72 項目で **得**（1 項目 2.2 µs、1 回 0.155 ms）
//
// 効くかどうかを決めるのは **1 回あたりの仕事とディスパッチのコストの比**です。
// そこで**下限は呼び出しごとに渡せる**ようにし、既定はプールが持ちます。
//
// 損益分岐は $\text{ディスパッチ} / (1 \text{項目あたりの仕事})$ で、
// ディスパッチは実測 **0.030 ms/回**（索引が 0.9 → 2.9 ms、66 回）。
//
// **これは「遅くならない」保証であって、最適な閾値ではありません。**
// 精密なチューニングは Phase 5 です（`IMPL-phase4.md` §4）。
#ifndef KRISITE_PAR_THREAD_POOL_HPP
#define KRISITE_PAR_THREAD_POOL_HPP

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace krisite::par {

class ThreadPool {
public:
    /// **ディスパッチの下限の既定値**（`SPEC-phase4.md` §6.3）。
    ///
    /// 項目数がこれ未満なら逐次で回します。**実測から決めた値**で、最適値では
    /// ありません（`IMPL-phase4.md` §4.1）。
    static constexpr std::size_t kDefaultMinItems = 64;

    /// `threads` 本で動く（0 か 1 なら**スレッドを作らない**）。
    explicit ThreadPool(unsigned threads) : n_(threads == 0 ? 1u : threads) {
        if (n_ <= 1) return;
        workers_.reserve(n_ - 1);
        for (unsigned k = 1; k < n_; ++k) workers_.emplace_back([this, k] { worker(k); });
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    ~ThreadPool() {
        if (n_ <= 1) return;
        {
            const std::lock_guard<std::mutex> g(m_);
            stop_ = true;
            ++gen_;
        }
        cv_.notify_all();
        for (std::thread& t : workers_) t.join();
    }

    unsigned size() const noexcept { return n_; }

    /// ディスパッチの下限（**実行時パラメータ**。§6.3）。
    ///
    /// **0 を渡すと下限そのものを無効化します**（呼び出しごとの上書きも効きません）。
    /// 検査で「必ず並列の経路を通す」ために要ります。
    /// **下限は性能のための機構なので、正しさの検査では外してください**
    /// （`IMPL-phase4.md` §4.2）。
    std::size_t min_items() const noexcept { return min_items_; }
    void set_min_items(std::size_t n) noexcept { min_items_ = n; }

    /// `fn(i, tid)` を $i = 0 \dots n-1$ について実行し、全部終わるまで待つ。
    ///
    /// **順序は保証しません。** `tid` は $0 \dots \text{size}()-1$ で、
    /// スレッド局所の器を引くための添字です。
    /// `min_items` を渡すと、この呼び出しだけ下限を上書きします（§6.3）。
    /// **1 項目あたりの仕事が軽い段では、既定より大きくしてください。**
    template <class F>
    void run(std::size_t n, F&& fn, std::size_t min_items = 0) {
        // **プールが 0 なら下限は無効**。呼び出しごとの上書きより強い
        const std::size_t floor =
            (min_items_ == 0) ? 0 : ((min_items != 0) ? min_items : min_items_);
        // **下限を下回る呼び出しは、その場で回します**（§6.3）
        if (n_ <= 1 || n < floor) {
            for (std::size_t i = 0; i < n; ++i) fn(i, 0u);
            return;
        }
        {
            const std::lock_guard<std::mutex> g(m_);
            job_ = std::function<void(std::size_t, unsigned)>(std::forward<F>(fn));
            total_ = n;
            next_.store(0, std::memory_order_relaxed);
            closed_ = false;
            ++gen_;
        }
        cv_.notify_all();
        drain(0);  // **呼び出し元も 1 本として働きます**
        // ここを抜けた時点で**全項目が誰かに取られています**（`next_ >= total_`）。
        // あとは**実際に走っているワーカー**が終わるのを待てば十分です。
        {
            std::unique_lock<std::mutex> lk(m_);
            cv_done_.wait(lk, [this] { return running_ == 0; });
            // **同じ危険区間で閉じます。** 遅れて起きたワーカーはこれを見て何もしません
            closed_ = true;
        }
    }

private:
    void drain(unsigned tid) {
        for (;;) {
            // **`relaxed` で足ります。** 同期は世代と完了数のほうで取ります
            const std::size_t i = next_.fetch_add(1, std::memory_order_relaxed);
            if (i >= total_) break;
            job_(i, tid);
        }
    }

    void worker(unsigned tid) {
        std::size_t seen = 0;
        for (;;) {
            std::unique_lock<std::mutex> lk(m_);
            cv_.wait(lk, [this, &seen] { return gen_ != seen; });
            seen = gen_;
            if (stop_) return;
            // **閉じた世代には触れません。** `job_` を読むのも危険です
            if (closed_) continue;
            ++running_;
            lk.unlock();
            drain(tid);
            lk.lock();
            --running_;
            lk.unlock();
            cv_done_.notify_one();
        }
    }

    unsigned n_;
    std::vector<std::thread> workers_;
    std::mutex m_;
    std::condition_variable cv_, cv_done_;
    std::size_t gen_ = 0;
    bool stop_ = false;
    std::function<void(std::size_t, unsigned)> job_;
    std::size_t total_ = 0;
    std::atomic<std::size_t> next_{0};
    /// **走り出したワーカーの数。** 呼び出し元はこれが 0 になるまで待ちます
    unsigned running_ = 0;
    /// 世代を閉じた印。**遅れて起きたワーカーへの合図**です
    bool closed_ = true;
    std::size_t min_items_ = kDefaultMinItems;
};

}  // namespace krisite::par

#endif  // KRISITE_PAR_THREAD_POOL_HPP

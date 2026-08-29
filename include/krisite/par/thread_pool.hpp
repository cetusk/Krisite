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

    /// `fn(i, tid)` を $i = 0 \dots n-1$ について実行し、全部終わるまで待つ。
    ///
    /// **順序は保証しません。** `tid` は $0 \dots \text{size}()-1$ で、
    /// スレッド局所の器を引くための添字です。
    template <class F>
    void run(std::size_t n, F&& fn) {
        if (n_ <= 1 || n <= 1) {
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
};

}  // namespace krisite::par

#endif  // KRISITE_PAR_THREAD_POOL_HPP

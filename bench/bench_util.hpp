// Krisite — ベンチマーク用の最小ヘルパ
#ifndef KRISITE_BENCH_UTIL_HPP
#define KRISITE_BENCH_UTIL_HPP

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

namespace kribench {

/// 最適化で消されないように値を握りつぶす。
template <class T>
inline void sink(const T& v) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "r,m"(v) : "memory");
#else
    volatile const T* p = &v;
    (void)p;
#endif
}

struct Result {
    std::string name;
    double ns_per_op;
    double mops;
};

inline std::vector<Result> g_results;

/// fn() を n 回まわして 1 回あたりの時間を測る。3 回計測して最小を採る。
template <class F>
inline void run(const char* name, long n, F&& fn) {
    double best = 1e300;
    for (int rep = 0; rep < 3; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        fn(n);
        const auto t1 = std::chrono::steady_clock::now();
        const double ns =
            std::chrono::duration_cast<std::chrono::duration<double, std::nano>>(t1 - t0).count();
        best = (ns / static_cast<double>(n) < best) ? ns / static_cast<double>(n) : best;
    }
    g_results.push_back(Result{name, best, 1000.0 / best});
    std::printf("  %-34s %10.2f ns/op %10.2f Mops/s\n", name, best, 1000.0 / best);
}

inline void print_markdown_table() {
    std::printf("\n<!-- docs/BENCH.md に貼り付ける -->\n\n");
    std::printf("| 演算 | ns/op | Mops/s |\n|---|---:|---:|\n");
    for (const auto& r : g_results) {
        std::printf("| `%s` | %.2f | %.1f |\n", r.name.c_str(), r.ns_per_op, r.mops);
    }
    std::printf("\n");
}

}  // namespace kribench

#endif  // KRISITE_BENCH_UTIL_HPP

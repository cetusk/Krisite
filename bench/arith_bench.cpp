// Krisite — fixed_int のスループット計測
//
// SPEC-phase0.md §10「まず正しさ。最適化は基準線を取ってから」に従い、
// ここでは基準線を取るだけで最適化はしない。
#include <vector>

#include "bench_util.hpp"
#include "test_util.hpp"

using namespace krisite::arith;
using kritest::Rng;

namespace {

constexpr std::size_t kPool = 4096;  // L1 に収まる程度

template <std::size_t N>
std::vector<fixed_int<N>> make_pool(std::uint64_t seed) {
    Rng rng(seed);
    std::vector<fixed_int<N>> v(kPool);
    for (auto& x : v) x = kritest::rand_full<N>(rng);
    return v;
}

template <std::size_t N>
void bench_width(const char* label) {
    const auto a = make_pool<N>(1), b = make_pool<N>(2);
    char name[96];

    std::snprintf(name, sizeof(name), "add<%s>", label);
    kribench::run(name, 5000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(detail::add_raw(a[k], b[k]));
        }
    });

    std::snprintf(name, sizeof(name), "mul<%s>", label);
    kribench::run(name, 2000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(mul(a[k], b[k]));
        }
    });

    std::snprintf(name, sizeof(name), "cmp<%s>", label);
    kribench::run(name, 5000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(cmp(a[k], b[k]));
        }
    });

    std::snprintf(name, sizeof(name), "sign<%s>", label);
    kribench::run(name, 10000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(sign(a[k]));
        }
    });
}

}  // namespace

int main() {
    std::printf("\nKrisite arith bench (b = %zu, リム数はビット幅から導出)\n\n",
                krisite::kCoordBits);
    bench_width<1>("1 リム / 64bit");
    bench_width<2>("2 リム / 128bit");
    bench_width<4>("4 リム / 256bit");
    kribench::print_markdown_table();
    return 0;
}

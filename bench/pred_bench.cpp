// Krisite — 述語のスループット計測
//
// SPEC-phase0.md §9:
//   「side() のスループットをベンチで計測し、数値を docs/BENCH.md に記録」
#include <array>
#include <vector>

#include "bench_util.hpp"
#include "test_util.hpp"

using namespace krisite::arith;
using namespace krisite::geom;
using kritest::Rng;

namespace {

constexpr std::size_t kPool = 2048;

}  // namespace

int main() {
    Rng rng(20260826);

    std::vector<IPoint> pts(kPool);
    for (auto& p : pts) p = kritest::rand_point(rng);

    std::vector<PlaneD> planes(kPool);
    for (auto& pl : planes) {
        do {
            pl = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                     kritest::rand_point(rng));
        } while (is_degenerate(pl));
    }

    // 3 平面の添字は必ず相異なるものを選ぶ。i, 7i+1, 13i+5 のような写像は
    // i = 170 などで 2 つが一致し、その組は絶対に交点を持たない（無限ループになる）。
    auto triple_of = [](std::size_t i) {
        return std::array<std::size_t, 3>{i & (kPool - 1), (i + 1) & (kPool - 1),
                                          (i + 2) & (kPool - 1)};
    };

    std::vector<HPointD> hpts;
    hpts.reserve(kPool);
    while (hpts.size() < kPool) {
        const auto t = triple_of(hpts.size());
        if (!kritest::intersects_at_point(planes[t[0]], planes[t[1]], planes[t[2]])) {
            planes[t[0]] = plane_from_triangle(kritest::rand_point(rng), kritest::rand_point(rng),
                                               kritest::rand_point(rng));
            continue;
        }
        hpts.push_back(intersect3(planes[t[0]], planes[t[1]], planes[t[2]]));
    }

    std::printf("\nKrisite predicate bench (b = %zu, side は %zu ビット / %zu リムで評価)\n\n",
                krisite::kCoordBits, bits::kSide, limbs::kSide);

    kribench::run("side(plane, HPoint)  ★最重要", 5000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(side(planes[k], hpts[k]));
        }
    });

    kribench::run("side(plane, IPoint)", 5000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(side(planes[k], pts[k]));
        }
    });

    kribench::run("orient3d(IPoint x4)", 5000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(orient3d(pts[k], pts[(k + 1) & (kPool - 1)], pts[(k + 2) & (kPool - 1)],
                                    pts[(k + 3) & (kPool - 1)]));
        }
    });

    kribench::run("orient2d(IPoint x3, Z)", 10000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(
                orient2d(pts[k], pts[(k + 1) & (kPool - 1)], pts[(k + 2) & (kPool - 1)], Axis::Z));
        }
    });

    kribench::run("cmp_h(HPoint x2, X)", 5000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(cmp_h(hpts[k], hpts[(k + 1) & (kPool - 1)], Axis::X));
        }
    });

    kribench::run("lex_less(HPoint x2)", 2000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(lex_less(hpts[k], hpts[(k + 1) & (kPool - 1)]));
        }
    });

    kribench::run("plane_from_triangle", 2000000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i) & (kPool - 1);
            kribench::sink(plane_from_triangle(pts[k], pts[(k + 1) & (kPool - 1)],
                                               pts[(k + 2) & (kPool - 1)]));
        }
    });

    // intersect3 は w != 0 が契約なので、成立する組だけを先に選んでおく
    std::vector<std::array<std::size_t, 3>> triples;
    triples.reserve(kPool);
    for (std::size_t k = 0; k < kPool; ++k) {
        const auto t = triple_of(k);
        if (kritest::intersects_at_point(planes[t[0]], planes[t[1]], planes[t[2]])) {
            triples.push_back(t);
        }
    }
    const std::size_t ntri = triples.size();
    kribench::run("intersect3(Plane x3)", 300000, [&](long n) {
        for (long i = 0; i < n; ++i) {
            const auto& t = triples[static_cast<std::size_t>(i) % ntri];
            kribench::sink(intersect3(planes[t[0]], planes[t[1]], planes[t[2]]));
        }
    });

    kribench::print_markdown_table();
    return 0;
}

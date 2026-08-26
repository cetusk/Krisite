// Krisite — ビット幅の実測
//
// SPEC-phase0.md §8.4
//   §3 の解析が正しいことを実測で確認する。各演算の実際の使用ビット幅の最大値を記録し、
//   (1) 理論上界を超えないこと、(2) 上界が過大でないこと、を検証する。
//
// 「過大でない」の基準は SPEC に無いので、**リム数が理論値と一致すること**を条件とする。
// すなわち実測値が上界より小さくても、必要リム数が変わらなければ設計上は問題ない。
#include <cstdio>
#include <initializer_list>

#include "test_util.hpp"

using namespace krisite::arith;
using namespace krisite::geom;
using kritest::Rng;

namespace {

struct Gauge {
    const char* name;
    std::size_t theory;    ///< §3 の上界（ビット）
    std::size_t measured;  ///< 実測の最大（ビット）
    std::size_t samples;
};

template <std::size_t N>
void feed(Gauge& g, const fixed_int<N>& x) {
    const std::size_t n = min_bits(x);
    if (n > g.measured) g.measured = n;
    ++g.samples;
}

}  // namespace

int main() {
    constexpr std::size_t b = krisite::kCoordBits;

    Gauge g_diff{"座標差 (b+1)", bits::kDiff, 0, 0};
    Gauge g_normal{"法線 N (2b+3)", bits::kNormal, 0, 0};
    Gauge g_offset{"オフセット d (3b+5)", bits::kOffset, 0, 0};
    Gauge g_w{"構成点 w (6b+12)", bits::kHomoW, 0, 0};
    Gauge g_xyz{"構成点 x,y,z (7b+14)", bits::kHomoXyz, 0, 0};
    Gauge g_side{"side(plane,HPoint) (9b+20)", bits::kSide, 0, 0};
    Gauge g_sidei{"side(plane,IPoint) (3b+6)", 3 * b + 6, 0, 0};
    Gauge g_o3d{"orient3d (3b+5)", bits::kOrient3d, 0, 0};
    Gauge g_o2d{"orient2d (2b+3)", bits::kOrient2d, 0, 0};
    Gauge g_cmph{"cmp_h (13b+27)", bits::kCmpH, 0, 0};
    // Phase 1（SPEC-phase1.md §7）
    Gauge g_pmin{"平面の小行列式 (5b+9)", bits::kPlaneMinor, 0, 0};
    Gauge g_vol{"四面体の体積x6 (3b+1)", bits::kTetraVolume6, 0, 0};
    Gauge g_o2dh{"orient2d_h (8b+17)", bits::kOrient2dH, 0, 0};
    Gauge g_msid{"side(plane,中点) (15b+33)", bits::kMidSide, 0, 0};
    Gauge g_mo2dh{"orient2d_h(中点) (14b+30)", bits::kMidOrient2dH, 0, 0};
    Gauge g_tsid{"side(plane,重心) (21b+46)", bits::kTriSide, 0, 0};
    Gauge g_to2dh{"orient2d_h(重心) (20b+43)", bits::kTriOrient2dH, 0, 0};

    Rng rng(20260826);

    // 極端な座標を厚めに引くことで上界に近づける
    auto pick = [&rng](int mode) {
        return (mode == 0) ? kritest::rand_extreme_point(rng) : kritest::rand_point(rng);
    };

    for (int iter = 0; iter < 400000; ++iter) {
        const int mode = (iter % 4 == 0) ? 1 : 0;  // 3/4 は格子の端寄り
        const IPoint p[9] = {pick(mode), pick(mode), pick(mode), pick(mode), pick(mode),
                             pick(mode), pick(mode), pick(mode), pick(mode)};

        feed(g_diff, coord_diff(p[0].x, p[1].x));
        feed(g_diff, coord_diff(p[0].y, p[1].y));
        feed(g_diff, coord_diff(p[0].z, p[1].z));

        const PlaneD pl0 = plane_from_triangle(p[0], p[1], p[2]);
        const PlaneD pl1 = plane_from_triangle(p[3], p[4], p[5]);
        const PlaneD pl2 = plane_from_triangle(p[6], p[7], p[8]);
        for (const PlaneD* pp : {&pl0, &pl1, &pl2}) {
            feed(g_normal, pp->a);
            feed(g_normal, pp->b);
            feed(g_normal, pp->c);
            feed(g_offset, pp->d);
        }

        feed(g_o3d, orient3d_value(p[0], p[1], p[2], p[3]));
        feed(g_o2d, orient2d_value(p[0].x, p[0].y, p[1].x, p[1].y, p[2].x, p[2].y));
        feed(g_sidei, side_value(pl0, p[6]));

        if (!kritest::intersects_at_point(pl0, pl1, pl2)) continue;
        const HPointD v = intersect3(pl0, pl1, pl2);
        feed(g_w, v.w);
        feed(g_xyz, v.x);
        feed(g_xyz, v.y);
        feed(g_xyz, v.z);
        // side は「その平面上に無い」構成点で測る（構成に使った平面では常に 0 になる）
        const PlaneD probe = plane_from_triangle(p[0], p[4], p[8]);
        feed(g_side, side_value(probe, v));

        // cmp_h は 2 つの構成点が要る。もう 1 点を作る
        const PlaneD pl3 = plane_from_triangle(p[1], p[4], p[7]);
        if (!kritest::intersects_at_point(pl0, pl1, pl3)) continue;
        const HPointD u = intersect3(pl0, pl1, pl3);
        feed(g_cmph, cmp_h_value(v, u, Axis::X));
        feed(g_cmph, cmp_h_value(v, u, Axis::Y));
        feed(g_cmph, cmp_h_value(v, u, Axis::Z));

        // Phase 1: 平面の小行列式と符号付き体積
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) feed(g_pmin, plane_minor(pl0, pl1, i, j));
        }
        feed(g_vol, tetra_volume6(p[0], p[1], p[2]));
        // レイキャストの投影向き（同次点版）
        feed(g_o2dh, orient2d_h_value(p[0], p[1], v, Axis::X));
        feed(g_o2dh, orient2d_h_value(p[2], p[3], u, Axis::Y));
        // §6.1 の代表点フォールバック（2 構成点の中点）
        const HMidPointD mid{v, u};
        feed(g_msid, side_value(probe, mid));
        feed(g_mo2dh, orient2d_h_value(p[0], p[1], mid, Axis::X));
        // 3 頂点の重心（三角形の断片用のフォールバック）
        const PlaneD pl4 = plane_from_triangle(p[2], p[5], p[8]);
        if (!kritest::intersects_at_point(pl0, pl2, pl4)) continue;
        const HTriPointD tri3{v, u, intersect3(pl0, pl2, pl4)};
        feed(g_tsid, side_value(probe, tri3));
        feed(g_to2dh, orient2d_h_value(p[0], p[1], tri3, Axis::X));
    }

    Gauge* all[] = {&g_diff,  &g_normal, &g_offset, &g_w,    &g_xyz,  &g_side,
                    &g_sidei, &g_o3d,    &g_o2d,    &g_cmph, &g_pmin, &g_vol,
                    &g_o2dh,  &g_msid,   &g_mo2dh,  &g_tsid, &g_to2dh};

    std::printf("\n  b = %zu のビット幅実測（SPEC-phase0.md §8.4）\n", b);
    std::printf("  実測値は乱択で到達した下限であり、真の最大値ではない。\n");
    std::printf(
        "  検証するのは「理論上界を超えないこと」（硬）と「上界が過大でないこと」（参考）。\n\n");
    std::printf("  %-30s %6s %6s %5s %8s %8s %s\n", "量", "理論", "実測", "余裕", "理論リム",
                "実測リム", "標本数");
    int oversized = 0;
    for (const Gauge* g : all) {
        const std::size_t tl = limbs_for(g->theory);
        const std::size_t ml = limbs_for(g->measured);
        std::printf("  %-30s %6zu %6zu %5zu %8zu %8zu %zu%s\n", g->name, g->theory, g->measured,
                    g->theory - g->measured, tl, ml, g->samples,
                    (ml == tl) ? "" : "  ← 乱択が届かず");

        // (1) 理論上界を超えないこと — 超えたら CLAUDE.md「即座に報告」対象
        KRI_CHECK_MSG(g->measured <= g->theory,
                      std::string(g->name) + ": 実測 " + std::to_string(g->measured) +
                          " ビットが理論上界 " + std::to_string(g->theory) + " ビットを超過");
        // (2) 上界が過大でないこと — 乱択では真の最大に届かないので参考値にとどめる。
        //     リム数が 2 以上ずれたら設計を見直す必要があるので、そこだけ硬い検査にする。
        KRI_CHECK_MSG(ml + 1 >= tl, std::string(g->name) + ": 理論 " + std::to_string(tl) +
                                        " リムに対し実測 " + std::to_string(ml) +
                                        " リム。上界が明らかに過大");
        if (ml != tl) ++oversized;
        KRI_CHECK(g->samples > 1000);
    }
    if (oversized > 0) {
        std::printf(
            "\n  注: SPEC §3.4「乱択の実測値でリム数を減らさないこと」。乱択は真の最大に\n"
            "      届かないので、この列がリム境界を下回っていても型は変更しない。減らす判断は\n"
            "      上界を再導出したときだけ行う。\n");
    }
    std::printf(
        "\n  注: 入力メッシュの体積は三角形ごとの総和なので、蓄積後の幅は\n"
        "      %zu ビット（%zu リム、三角形 2^%zu 枚まで）を確保してある。\n"
        "      総和の実測はメッシュのテスト側で行う。\n",
        bits::kInputVolume6, limbs::kInputVolume6, bits::kMaxTrianglesLog2);
    std::printf("\n  Phase 0 の述語が要求する最大リム数（§3.3）: %zu\n", limbs::kMaxPredicate);
    KRI_CHECK(limbs::kMaxPredicate == limbs_for(bits::kCmpH));
    std::printf("\n");

    return kritest::finish("arith/bitwidth");
}

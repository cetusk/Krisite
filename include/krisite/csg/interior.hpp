// Krisite — 断片の相対内部の点（代表点）
//
// SPEC-phase3.md §2.1（段 0）。**構成を平面ベースに戻します。**
//
// Phase 1 / 2 は「頂点 → 対角線の中点 → 3 頂点の重心」の 3 段で代表点を選んでいました。
// 2 段目以降は**点を組み合わせる**操作で、被符号値が 21b+46（b=21 で 487 ビット）に
// 達し、`widths.hpp` に中点系・重心系の定数を 4 本抱えていました。
//
// EMBER §4.4 に倣い、代表点そのものを**平面の交点**として構成します。
//
//   主経路   多角形の重心を【浮動小数点】で求めて最近傍の整数 c に丸め、
//            支持平面の法線に最も近い軸に沿った軸平行直線と支持平面の交点を取る。
//            **判定はすべて厳密です。** float は候補を出すだけで、正しさに関与しません
//   予備経路 角を 1 つ取り、その 2 枚の辺平面を【内側へずらして】交点を取る
//
// **どちらも一般の 3 平面交点なので、`side` は 9b+20 のままです。**
// 主経路は軸平行 2 枚を含むのでさらに小さくなります（widths.hpp bits::kAxisPointXyz）。
#ifndef KRISITE_CSG_INTERIOR_HPP
#define KRISITE_CSG_INTERIOR_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "krisite/csg/fragment.hpp"
#include "krisite/csg/plane_table.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/geom/point.hpp"
#include "krisite/geom/predicates.hpp"

namespace krisite::csg {

/// §11 の記録用。**どちらの経路で決まったか**を数えます。
struct InteriorStats {
    std::size_t axis_line = 0;      ///< 主経路で決まった回数
    std::size_t corner_offset = 0;  ///< 予備経路で決まった回数
    std::size_t axis_failed = 0;    ///< 主経路が外れた回数
    std::size_t corner_tries = 0;   ///< 予備経路で試した角 x 倍率の総数
};

namespace detail {

/// 固定幅整数の**近似値**（float ヒント専用。正しさには関与しません）。
///
/// 符号反転を使わずに済ませます。最上位リムを符号付き、それ以外を符号なしとして
/// 桁を積むだけなので、最小値（-2^(64N-1)）でも破綻しません。
template <std::size_t N>
inline double approx(const arith::fixed_int<N>& x) noexcept {
    double r = static_cast<double>(static_cast<std::int64_t>(x.limb[N - 1]));
    for (std::size_t i = N - 1; i-- > 0;) {
        r = r * 18446744073709551616.0 + static_cast<double>(x.limb[i]);
    }
    return r;
}

/// 辺平面 `e` について、断片の**内側の符号**を返す（頂点のうち載っていないものの符号）。
///
/// 凸多角形なので、辺の上に無い頂点はすべて同じ側にあります。
inline int inward_sign(const PlaneTable& t, const std::vector<geom::HPointD>& vs, PlaneId e) {
    for (const geom::HPointD& v : vs) {
        const int s = geom::side(t.at(e), v);
        if (s != 0) return s;
    }
    return 0;  // 全頂点が載っている = 退化した断片
}

}  // namespace detail

/// 断片の**相対内部**の点を構成する（SPEC-phase3 §2.1）。
///
/// 返す点は、支持平面上にあり、**すべての辺平面に対して内側（符号が非零で内側と一致）**
/// です。相対内部の点は相手の平面配置のセルの内部にあるので、$\partial B$ には
/// 決して載りません（`SPEC-phase1.md` §6.1 の論証）。
/// `variant` で**別の内部点**を求めます（`SPEC-phase3.md` §3.3 / §5.5）。
///
/// セグメントトレースの経路は、内部点が軸に整列していると退化しやすくなります
/// （軸平行な入力では、内部点を通る軸平行線が他の面の辺をちょうど通る）。
/// **経路が退化したら別の内部点で試す**ための引数です。
///
///   variant = 0        主経路（float ヒント + 軸平行直線）
///   variant >= 1       予備経路（角ごと・刻みごと。順に別の角を使う）
inline geom::HPointD interior_point(PlaneTable& t, const Fragment& f, PointCache* cache = nullptr,
                                    InteriorStats* st = nullptr,
                                    std::array<PlaneId, 3>* planes = nullptr, unsigned variant = 0,
                                    geom::IPoint* anchor = nullptr, geom::Axis* anchor_axis = nullptr) {
    const std::size_t n = vertex_count(f);
    KRISITE_CHECK(n >= 3, "interior_point: 頂点が 3 未満");

    std::vector<geom::HPointD> vs(n);
    for (std::size_t i = 0; i < n; ++i) vs[i] = fragment_vertex(t, f, i, cache);

    std::vector<int> inward(n);
    for (std::size_t k = 0; k < n; ++k) inward[k] = detail::inward_sign(t, vs, f.edge[k]);

    // 候補が相対内部にあるか。**辺平面すべてに対して内側で、かつ載っていないこと。**
    auto strictly_inside = [&](const geom::HPointD& x) {
        for (std::size_t k = 0; k < n; ++k) {
            if (inward[k] == 0) return false;  // 退化した断片
            if (geom::side(t.at(f.edge[k]), x) != inward[k]) return false;
        }
        return true;
    };

    const geom::PlaneD& sp = t.at(f.support);

    // ---- 主経路: float ヒント + 軸平行直線（EMBER §4.4）------------------------
    if (variant == 0) {
        double cx = 0.0, cy = 0.0, cz = 0.0;
        for (const geom::HPointD& v : vs) {
            const double w = detail::approx(v.w);
            cx += detail::approx(v.x) / w;
            cy += detail::approx(v.y) / w;
            cz += detail::approx(v.z) / w;
        }
        const double inv = 1.0 / static_cast<double>(n);
        const double cd[3] = {cx * inv, cy * inv, cz * inv};

        // 支持平面の法線に最も近い軸を選ぶ（**厳密に**比較する。float は座標だけ）
        const arith::fixed_int<geom::limbs::kNormal> nc[3] = {sp.a, sp.b, sp.c};
        std::size_t axis = 0;
        for (std::size_t i = 1; i < 3; ++i) {
            const auto ai = arith::is_negative(nc[i]) ? arith::neg(nc[i]) : nc[i];
            const auto ab = arith::is_negative(nc[axis]) ? arith::neg(nc[axis]) : nc[axis];
            if (arith::cmp(ai, ab) > 0) axis = i;
        }

        std::int64_t ic[3];
        bool in_range = true;
        for (std::size_t i = 0; i < 3; ++i) {
            const double r = (cd[i] >= 0.0) ? (cd[i] + 0.5) : (cd[i] - 0.5);
            if (!(r > static_cast<double>(kCoordMin) && r < static_cast<double>(kCoordMax))) {
                in_range = false;
                break;
            }
            ic[i] = static_cast<std::int64_t>(r);
        }
        if (in_range) {
            // 軸平行直線 = 残り 2 軸の軸平行平面の交わり
            const geom::Axis other[3][2] = {{geom::Axis::Y, geom::Axis::Z},
                                            {geom::Axis::X, geom::Axis::Z},
                                            {geom::Axis::X, geom::Axis::Y}};
            const geom::PlaneD a0 =
                geom::plane_axis_aligned(other[axis][0], ic[static_cast<int>(other[axis][0])]);
            const geom::PlaneD a1 =
                geom::plane_axis_aligned(other[axis][1], ic[static_cast<int>(other[axis][1])]);
            // w = ±N_s[axis] != 0（最大成分を選んだので非零）
            // 交わらない組は作れない（軸の選び方から起きないはずだが、安全側に）
            if (!geom::planes_meet_at_point(a0, a1, sp)) return vs[0];
            const geom::HPointD x = geom::intersect3(a0, a1, sp);
            if (strictly_inside(x)) {
                if (st != nullptr) ++st->axis_line;
                if (planes != nullptr) {
                    // **トレースの始点に要ります**（§3.3 の経路は平面 3 つ組で作る）。
                    // 軸平行平面を表に登録します。幅は kAxisOffset <= kOffset で収まります。
                    *planes = {t.intern(a0).id, t.intern(a1).id, f.support};
                }
                // **主経路の整数アンカー**（§3.3.0）。この点と x は軸平行な線で結べます。
                if (anchor != nullptr) {
                    *anchor = geom::IPoint{static_cast<std::int32_t>(ic[0]),
                                           static_cast<std::int32_t>(ic[1]),
                                           static_cast<std::int32_t>(ic[2])};
                }
                if (anchor_axis != nullptr) *anchor_axis = static_cast<geom::Axis>(axis);
                return x;
            }
        }
        if (st != nullptr) ++st->axis_failed;
    }

    // ---- 予備経路: 角の 2 平面を内側へずらす（EMBER §4.4 の第 2 の方法）--------
    //
    // 辺平面を内側へ動かすと、交点はその角の近傍の**内部**に入ります。
    // 動かす量は d の +-1 が最小ですが、平面係数をスケールアップすれば
    // 幾何的な移動量を細かくできます（EMBER §4.4）。**幅の余裕の範囲で行います。**
    // variant >= 1 なら、その番号ぶんだけ後ろの角から試す（別の点になる）
    const std::size_t start = (variant == 0) ? 0 : ((variant - 1) % n);
    for (std::size_t ii = 0; ii < n; ++ii) {
        const std::size_t i = (start + ii) % n;
        const PlaneId ea = f.edge[(i + n - 1) % n];
        const PlaneId eb = f.edge[i];
        const std::size_t ka = (i + n - 1) % n, kb = i;
        if (inward[ka] == 0 || inward[kb] == 0) continue;

        for (int step = 0; step < 2; ++step) {
            const geom::PlaneD pa = geom::plane_scaled_for_shift(t.at(ea), step == 1);
            const geom::PlaneD pb = geom::plane_scaled_for_shift(t.at(eb), step == 1);
            // 内側が負なら d を増やすと平面が内側へ動く（side = N・x + d・w）
            const geom::PlaneD sa = geom::plane_shifted(pa, -inward[ka]);
            const geom::PlaneD sb = geom::plane_shifted(pb, -inward[kb]);
            if (st != nullptr) ++st->corner_tries;
            // ずらしても平面の向きは変わらないので、元の 3 枚が一点で交わるなら
            // ずらした 3 枚も交わります（法線の行列式は同じ）。角の 2 辺は
            // 支持平面上で交わるので退化しません。
            if (!geom::planes_meet_at_point(sp, sa, sb)) continue;
            const geom::HPointD x = geom::intersect3(sp, sa, sb);
            if (strictly_inside(x)) {
                if (st != nullptr) ++st->corner_offset;
                if (planes != nullptr) {
                    *planes = {f.support, t.intern(sa).id, t.intern(sb).id};
                }
                return x;
            }
        }
    }

    KRISITE_CHECK(false,
                  "interior_point: 相対内部の点を構成できない"
                  "（主経路も予備経路も外れた。SPEC-phase3 §2.1）");
    return vs[0];
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_INTERIOR_HPP

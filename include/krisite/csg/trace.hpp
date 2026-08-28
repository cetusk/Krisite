// Krisite — セグメントトレースによる巻き数の伝播（SPEC-phase3.md §3.3 / §5.5）
//
// **レイキャストを置き換えます。**
//
// レイキャストは分類を**大域的な問題**にします（無限遠まで全多角形と交差判定する）。
// 既知の巻き数を持つ**参照点**があれば、そこからの経路が跨いだ面だけを数えれば済み、
// 分類が局所問題になります（EMBER §3.4）。**Phase 4 の並列化の前提**です。
//
// ---
//
// ## 2 つの同次点を結ぶ線分は作れません（EMBER §3.2）
//
// 差を取ると分母が $x_4 x'_4$ になり、幅が倍になります。**しかし経路なら作れます。**
//
//   $x = (s_0, s_1, s_2)$、$x_{ref} = (r_0, r_1, r_2)$ のとき
//   $x_1 = (s_0, s_1, r_2)$ は $x$ と平面 $s_0, s_1$ を共有する → **同じ線の上**
//   $x_2 = (s_0, r_1, r_2)$ も同様
//
// **平面を 1 枚ずつ置き換えていく**操作です。順序は任意なので、経路が退化したら
// 別の順序を試せます（§3.3）。
//
// ## ビット幅は増えません ★
//
// 経路の点も交点も、**入力平面から 3 枚選んだ交点**です。判定はすべて `side`。
// **新しい述語は要りません**（`SPEC-phase3.md` §7 の「セグメント–多角形の交差」と
// 「セグメントトレースで使う述語」は、既存の `intersect3` と `side` に落ちます）。
#ifndef KRISITE_CSG_TRACE_HPP
#define KRISITE_CSG_TRACE_HPP

#include <array>
#include <vector>

#include "krisite/csg/fragment.hpp"
#include "krisite/csg/plane_table.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/geom/predicates.hpp"
#if defined(KRISITE_DEBUG_TRACE)
#include <cstdio>
#endif

namespace krisite::csg {

/// 平面 3 つ組で表した点。**座標を持ちません。**
using PPoint = std::array<PlaneId, 3>;

/// トレースの対象。多角形 1 枚と、その巻き数への寄与。
struct TracePoly {
    PlaneId support = kNoPlane;
    std::vector<PlaneId> edge;
    std::vector<std::int8_t> inward;  ///< 辺平面ごとの内側の符号（前もって求めておく）
    std::uint32_t comp = 0;           ///< どの WNV 成分か
    int orient = +1;                  ///< 外向き法線が支持平面の法線と同じなら +1
};

/// トレースの結果。`ok` が偽なら**経路が退化**しているので、別の経路を試します。
struct TraceResult {
    bool ok = false;
    std::vector<std::int32_t> dw;  ///< 始点から終点への巻き数の変化
};

namespace detail {

/// 線分 $(l_0, l_1)$ 上の点 $p$ が、端平面 $e$ について「$q$ と同じ側」か。
inline bool same_side(const PlaneTable& t, PlaneId e, const geom::HPointD& p,
                      const geom::HPointD& q, bool* degenerate) {
    const int sp = geom::side(t.at(e), p);
    const int sq = geom::side(t.at(e), q);
    if (sp == 0 || sq == 0) {
        *degenerate = true;
        return false;
    }
    return sp == sq;
}

}  // namespace detail

/// 線分 1 本ぶんのトレース。$(l_0, l_1)$ が線、$a$ から $b$ へ進む。
///
/// $a$ は端平面 `ea` の上に、$b$ は `eb` の上にあります。
///
/// **跨いだ向きは端点の符号から決まります。** 差を取る必要はありません
/// （$\mathrm{side}$ の符号が反対なら跨いでいて、終点側の符号が進行方向）。
inline bool trace_segment(const PlaneTable& t, const std::vector<TracePoly>& polys, PlaneId l0,
                          PlaneId l1, PlaneId ea, PlaneId eb, const geom::HPointD& a,
                          const geom::HPointD& b, std::vector<std::int32_t>* w,
                          bool a_on_surface = false, bool b_on_surface = false) {
    for (const TracePoly& q : polys) {
        const geom::PlaneD& sp = t.at(q.support);
        const int sa = geom::side(sp, a);
        const int sb = geom::side(sp, b);
        if (sa == 0 && sb == 0) {
#if defined(KRISITE_DEBUG_TRACE)
            std::fprintf(stderr, "  退化: 線分が支持平面 %u に乗っている\n", q.support);
#endif
            return false;  // 線分が支持平面に乗っている → 退化
        }
        if (sa == 0 || sb == 0) {
            // **端点が支持平面の上。** 分類する面の上に始点を置くので、これは常態です。
            //
            //   その端点が「面の上」だと分かっている（始点）→ 跨がない。表裏は
            //   呼び出し側が `delta` で扱う
            //   そうでない中間点が面に載っている → 経路が面をかすめている → 退化
            const bool known = (sa == 0) ? a_on_surface : b_on_surface;
            if (known) continue;
            const geom::HPointD& on = (sa == 0) ? a : b;
            bool inside_poly = true;
            for (std::size_t k = 0; k < q.edge.size(); ++k) {
                if (geom::side(t.at(q.edge[k]), on) != q.inward[k]) {
                    inside_poly = false;
                    break;
                }
            }
            if (inside_poly) {
#if defined(KRISITE_DEBUG_TRACE)
                std::fprintf(stderr, "  退化: 端点が多角形 %u の上（known=%d）\n", q.support,
                             (int)known);
#endif
                return false;  // 多角形の上を通っている → 退化
            }
            continue;                       // 平面上だが多角形の外 → 跨がない
        }
        if (sa == sb) continue;  // 跨いでいない

        // 交点。線が支持平面と平行なら跨げない（w = 0 は起きない: sa != sb なので）
        const geom::HPointD p = geom::intersect3(t.at(l0), t.at(l1), sp);

        // 多角形の内部か（**厳密に内側**。辺の上なら経路が退化）
        bool inside = true;
        for (std::size_t k = 0; k < q.edge.size(); ++k) {
            const int s = geom::side(t.at(q.edge[k]), p);
            if (s == 0) {
#if defined(KRISITE_DEBUG_TRACE)
                std::fprintf(stderr, "  退化: 交点が辺平面 %u の上\n", q.edge[k]);
#endif
                return false;
            }
            if (s != q.inward[k]) {
                inside = false;
                break;
            }
        }
        if (!inside) continue;

        // 線分の内側か（両端平面について、それぞれ反対の端点と同じ側）
        bool degenerate = false;
        const bool between = detail::same_side(t, ea, p, b, &degenerate) &&
                             detail::same_side(t, eb, p, a, &degenerate);
        if (degenerate) {
#if defined(KRISITE_DEBUG_TRACE)
            std::fprintf(stderr, "  退化: 交点が端平面の上\n");
#endif
            return false;
        }
        if (!between) continue;

        // **跨いだ向き。** 終点が支持平面の正側なら、法線の側へ抜けた。
        // 外向き法線が +N（orient = +1）なら、それは立体から【出た】ので -1。
        const int dir = (sb > 0) ? -q.orient : q.orient;
        (*w)[q.comp] += dir;
    }
    return true;
}

/// 経路（最大 3 セグメント）で $x$ から $x_{ref}$ へトレースする（§3.3）。
///
/// $x$ を定義する平面を 1 枚ずつ $x_{ref}$ のものに置き換えます。**置き換える順序は
/// 任意**なので、退化したら別の順序を試します（6 通り）。
inline TraceResult trace_path(const PlaneTable& t, const std::vector<TracePoly>& polys,
                              const PPoint& x, const PPoint& xref,
                              const std::vector<std::int32_t>& w_ref, PointCache* cache = nullptr) {
    TraceResult r;
    static constexpr int kOrder[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                                         {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (const auto& ord : kOrder) {
        // 中間点: x の平面を ord の順に xref の平面へ置き換える
        PPoint cur = x;
        std::vector<std::int32_t> w = w_ref;
        bool ok = true;
        // **参照点から出発して x へ向かう**ほうが自然だが、置き換えは x 側から
        // 定義されるので、経路を作ってから逆向きに足す
        std::vector<PPoint> pts{x};
        for (int k = 0; k < 3; ++k) {
            cur[static_cast<std::size_t>(ord[k])] = xref[static_cast<std::size_t>(ord[k])];
            pts.push_back(cur);
        }
        // pts.back() == xref のはず（3 枚とも置き換えた）
        std::vector<geom::HPointD> vals;
        vals.reserve(pts.size());
        for (const PPoint& p : pts) {
            // **交わるかを先に確かめます**（`intersect3` は交わらない組を表明で弾く）
            if (!geom::planes_meet_at_point(t.at(p[0]), t.at(p[1]), t.at(p[2]))) {
#if defined(KRISITE_DEBUG_TRACE)
                std::fprintf(stderr, "  退化: 経路の点が作れない（平面が一点で交わらない）\n");
#endif
                ok = false;
                break;
            }
            const geom::HPointD v =
                (cache != nullptr) ? cache->get(t, p[0], p[1], p[2])
                                   : geom::intersect3(t.at(p[0]), t.at(p[1]), t.at(p[2]));
            if (arith::is_zero(v.w)) {
#if defined(KRISITE_DEBUG_TRACE)
                std::fprintf(stderr, "  退化: 経路の点が作れない（w=0）\n");
#endif
                ok = false;
                break;
            }
            vals.push_back(v);
        }
        if (!ok) continue;

        // xref から x へ向かって足す（w は xref のもの）
        for (int k = 2; k >= 0 && ok; --k) {
            // セグメント pts[k] → pts[k+1] は、共有する 2 平面が線
            const std::size_t j = static_cast<std::size_t>(ord[k]);
            std::array<PlaneId, 2> line{};
            std::size_t n = 0;
            for (std::size_t m = 0; m < 3; ++m) {
                if (m != j) line[n++] = pts[static_cast<std::size_t>(k)][m];
            }
            // 端平面: 始点側は pts[k][j]、終点側は pts[k+1][j]
            const PlaneId ea = pts[static_cast<std::size_t>(k)][j];
            const PlaneId eb = pts[static_cast<std::size_t>(k) + 1][j];
            if (ea == eb) continue;  // 置き換えても同じ平面 → 動いていない
            // 経路の始点（k == 0 の a 側）は**分類する面の上**にあります
            ok = trace_segment(t, polys, line[0], line[1], ea, eb,
                               vals[static_cast<std::size_t>(k)],
                               vals[static_cast<std::size_t>(k) + 1], &w, k == 0, false);
        }
        if (!ok) continue;
        r.ok = true;
        r.dw = std::move(w);
        return r;
    }
    return r;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_TRACE_HPP

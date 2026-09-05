// Krisite — 点が立体の内部にあるかの判定（レイキャスト）
//
// SPEC-phase1.md §6.1
//
// +X 方向のレイを撃ち、入力三角形との交差回数の偶奇で判定します。
//
// **前提: 判定点は立体の境界上にありません。** 呼び出し側が保証します
// （断片の符号ベクトルがすべて非零の頂点を選ぶ）。これにより前方交差の判定
// `N・p + d != 0` が保証され、**記号的摂動は投影の内外判定だけで済みます。**
//
// 記号的摂動: レイ原点を実座標で `(0, ε, ε²)` だけずらしたものとして扱います。
// 投影の向きが 0 になったとき、ε の項 `-(b.z - a.z)`、それも 0 なら ε² の項
// `(b.y - a.y)` で符号を決めます。
#ifndef KRISITE_CSG_RAYCAST_HPP
#define KRISITE_CSG_RAYCAST_HPP

#include "krisite/csg/ray_index.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/geom/predicates.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace krisite::csg {

namespace detail {

// `comp` / `proj_u` / `proj_v` は `ray_index.hpp` にあります（索引と共有）。

/// 平面の法線の、軸 `along` 方向の成分。
inline int normal_comp_sign(const geom::PlaneD& pl, geom::Axis along) noexcept {
    return (along == geom::Axis::X)   ? arith::sign(pl.a)
           : (along == geom::Axis::Y) ? arith::sign(pl.b)
                                      : arith::sign(pl.c);
}

/// 摂動後の投影向き。0 を返しません。
///
/// レイ原点を実座標で $(\varepsilon_u, \varepsilon_v^2)$（投影後の 2 軸）だけずらした
/// ものとして扱います。投影の向きが 0 になったとき、$\varepsilon$ の項 $-(b_v - a_v)$、
/// それも 0 なら $\varepsilon^2$ の項 $(b_u - a_u)$ で符号を決めます。
///
/// **軸によらず同じ形です**（`orient2d_h` の巡回順に合わせてあります）。
inline int perturbed_orient(int raw, const geom::IPoint& a, const geom::IPoint& b,
                            geom::Axis along) noexcept {
    if (raw != 0) return raw;
    const geom::Axis u = proj_u(along), v = proj_v(along);
    if (comp(b, v) != comp(a, v)) return (comp(a, v) > comp(b, v)) ? 1 : -1;
    if (comp(b, u) != comp(a, u)) return (comp(b, u) > comp(a, u)) ? 1 : -1;
    return 0;  // 投影が退化（法線の along 成分が 0）。呼び出し側が弾く
}

inline int proj_orient(const geom::IPoint& a, const geom::IPoint& b, const geom::IPoint& p,
                       geom::Axis along = geom::Axis::X) noexcept {
    const geom::Axis u = proj_u(along), v = proj_v(along);
    return perturbed_orient(
        geom::orient2d(comp(a, u), comp(a, v), comp(b, u), comp(b, v), comp(p, u), comp(p, v)), a,
        b, along);
}

inline int proj_orient(const geom::IPoint& a, const geom::IPoint& b, const geom::HPointD& p,
                       geom::Axis along = geom::Axis::X) noexcept {
    return perturbed_orient(geom::orient2d_h(a, b, p, along), a, b, along);
}

/// レイ（`along` の正方向）が三角形を前方で横切るか。**平面を渡す版。**
///
/// **支持平面は呼び出し側で一度だけ作ってください**（`SPEC-phase5.md` の CP1.5）。
/// レイキャストは source の全三角形を走査するので、ここで作り直すと
/// **1 レイあたり $O(n)$ 回の平面構成**になります（実測 14.7 億回）。
template <class Point>
bool crosses_with(const geom::PlaneD& pl, const geom::IPoint& a, const geom::IPoint& b,
                  const geom::IPoint& c, const Point& p, geom::Axis along = geom::Axis::X) {
    if (geom::is_degenerate(pl)) return false;
    // 法線の along 成分が 0 ⟺ 投影三角形が退化 ⟺ レイが平面に平行
    const int nx = normal_comp_sign(pl, along);
    if (nx == 0) return false;

    const int o1 = proj_orient(a, b, p, along);
    const int o2 = proj_orient(b, c, p, along);
    const int o3 = proj_orient(c, a, p, along);
    if (o1 == 0 || o2 == 0 || o3 == 0) return false;
    if (o1 != o2 || o2 != o3) return false;  // 投影三角形の外

    // 前方か: t = -(N・p + d)/N_along > 0 ⟺ sign(N・p + d) != sign(N_along)
    const int sp = geom::side(pl, p);
    KRISITE_CHECK(sp != 0, "crosses: 判定点が三角形の平面上にある（呼び出し側の契約違反）");
    return sp != nx;
}

/// 平面を渡さない版（互換）。**ホットパスでは `crosses_with` を使ってください。**
template <class Point>
bool crosses(const geom::IPoint& a, const geom::IPoint& b, const geom::IPoint& c, const Point& p,
             geom::Axis along = geom::Axis::X) {
    return crosses_with(geom::plane_from_triangle(a, b, c), a, b, c, p, along);
}

}  // namespace detail

/// 点が `m` の**境界上**にあるか（三角形の内部・辺・頂点のいずれか）。
///
/// `point_inside` の前提（判定点が境界上に無い）を呼び出し側が確かめるための述語です。
///
/// **「相手の平面上に無い」で代用してはいけません。** 平面上にあることと境界上にある
/// ことは別です。断片の頂点は隣接する切断平面の上に必ず載るので、平面で判定すると
/// 実際には境界から外れた点まで弾いてしまい、代表点が見つからなくなります
/// （CP2 のケース 2 で実際に起きました。IMPL-phase1 §6.6）。
///
/// 射影軸は**法線成分が非零の軸**を選びます。そうすれば射影した三角形が潰れません。
template <class Point>
inline bool point_on_boundary(const mesh::TriMesh& m, const Point& p) {
    for (const mesh::Tri& t : m.triangles) {
        const geom::IPoint& a = m.vertices[t[0]];
        const geom::IPoint& b = m.vertices[t[1]];
        const geom::IPoint& c = m.vertices[t[2]];
        const geom::PlaneD pl = geom::plane_from_triangle(a, b, c);
        if (geom::is_degenerate(pl)) continue;
        if (geom::side(pl, p) != 0) continue;  // 三角形の平面上にすら無い
        const geom::Axis ax = (arith::sign(pl.a) != 0)   ? geom::Axis::X
                              : (arith::sign(pl.b) != 0) ? geom::Axis::Y
                                                         : geom::Axis::Z;
        const int o1 = geom::orient2d_h(a, b, p, ax);
        const int o2 = geom::orient2d_h(b, c, p, ax);
        const int o3 = geom::orient2d_h(c, a, p, ax);
        const bool neg = (o1 < 0) || (o2 < 0) || (o3 < 0);
        const bool pos = (o1 > 0) || (o2 > 0) || (o3 > 0);
        if (!(neg && pos)) return true;  // 内部・辺・頂点のいずれか
    }
    return false;
}

/// `winding_split` の補助データ。**すべて任意**で、渡さなければ従来どおり動きます。
///
/// **位置引数を並べずに構造で区切ります**（`CLAUDE.md`「位置引数が多く型が同じ
/// 小関数は、規律ではなく設計で守ってください」）。
struct RaySupport {
    /// 支持平面。`planes[j]` は `plane_from_triangle(m.triangles[j] の 3 頂点)` と
    /// **同一**でなければなりません。**`PlaneTable` で intern した平面は使えません** —
    /// `intern` は `orientation_differs` を返すので、**表の平面は符号が逆のことがあります。**
    const geom::PlaneD* planes = nullptr;
    /// 軸ごとの 2 次元索引（`Axis` の値で引きます）。**どの軸を使うかは
    /// 基準平面から中で決まる**ので、3 軸ぶん渡します。
    const RayIndex* index[3] = {nullptr, nullptr, nullptr};
    /// 実際に走査した三角形の数を足し込みます（計測用。任意）。
    std::size_t* tested = nullptr;
    /// **前判定を通った候補の数**（`prefilter` が偽なら候補と同数。計測用。任意）。
    /// **`tested` との比が D の絞り込みの効きです。**
    std::size_t* kept = nullptr;
    /// **そのうち実際に寄与した三角形の数**（計測用。任意）。
    ///
    /// **`tested` との比が、索引の絞り込みの効きです**
    /// （`DESIGN-phase5-hotspots.md` §11）。**A-3 と同じ形の問い**で、
    /// A-3 では「索引が返す候補の 99.9% が捨てられている」でした。
    std::size_t* hits = nullptr;
    /// **投影した AABB に判定点が入る候補の数**（計測用。任意）。
    ///
    /// **索引はセルの単位でしか絞れません。** 三角形が段 $l$ の 4 セルを占めても、
    /// **その三角形自身の投影 AABB はもっと小さい**ことがあります。
    /// **この差が、安い前判定で落とせる分です**（`DESIGN-phase5-hotspots.md` §11）。
    std::size_t* aabb_pass = nullptr;
    /// **★ D-2 → D-1 の前判定を有効にする**（`DESIGN-phase5-hotspots.md` §11）。
    ///
    /// **索引はセルの単位でしか絞れません。** 三角形が段 $l$ の 4 セルを占めても、
    /// **その三角形自身の投影 AABB はもっと小さい**ことがあります。
    /// **実測でレイあたりの候補 210〜734 のうち、寄与するのは 2.5〜2.8 だけ**でした。
    ///
    ///   **D-2**  レイの前方にあるか（`along` の最大座標。**比較 1 回**）
    ///   **D-1**  投影した AABB に判定点が入るか（**比較 2.3〜2.7 回**、早期打ち切り）
    ///
    /// **どちらも厳密な絞り込みです。** 落とすのは「必ず寄与しないもの」だけ。
    ///
    /// > **★ 境界は【残す】側です。** 判定点が AABB の面にちょうど載る場合、
    /// > あるいは三角形の `along` 最大座標にちょうど一致する場合は、落としません。
    /// > **落とすと格子線の上の交点を見逃します**（`SPEC-phase1.md` §9.3 と同じ形）。
    bool prefilter = false;
    /// **安い前判定が実際に走った回数**（早期打ち切りあり。計測用）。
    /// **費用は「回数 × 単価」で出るので、時間のばらつきに依りません。**
    std::size_t* cheap_tests = nullptr;
    /// **レイの前方にある候補の数**（D-1 を掛けない、D-2 単独の効き）。
    std::size_t* fwd_only = nullptr;
    /// **投影 AABB を通り、かつレイの前方にある候補の数**（計測用。任意）。
    ///
    /// **D-1 の【後で】D-2 が何を落とすかを順に数えます。**
    /// 独立とは限らないので、順に適用した実測が要ります。
    std::size_t* fwd_pass = nullptr;
    /// **段ごとの候補数**（計測用。任意。`kRayIndexMaxLevels` 個の配列）。
    ///
    /// **大きい三角形が粗い段に入っているか**を見ます
    /// （`SPEC-phase3.md` §5.5.0 の「未調査」）。
    std::size_t* per_level = nullptr;
};

/// **巻き数**（`SPEC-phase3.md` §5.1）。判定点が `m` の表面に載っている場合も扱います。
///
/// 表面に載っている点では巻き数が両側で違うので、**3 つに分けて返します。**
///
///   `w_other`  判定点を含まない面だけから決まる巻き数（表裏で共通）
///   `c_front`  基準平面の法線 $+N$ 側で、載っている面が寄与する分
///   `c_back`   $-N$ 側で寄与する分
///
/// 表側の巻き数は `w_other + c_front`、裏側は `w_other + c_back` です。
///
/// **面に載っている点をレイキャストしてはいけません**（`crosses` の契約違反）。
/// 載っている面を先に除くのが、この関数の役割です。
///
/// 巻き数の寄与は `sign(N_t \cdot x)`（レイは +X 方向）。閉じた向き付き立体なら
/// 内部で 1、外部で 0 になります。**自己交差や入れ子では 2 以上になります**（それが目的）。
///
/// `sup`（`RaySupport`）で支持平面と 2 次元索引を渡せます（`SPEC-phase5.md` の CP1.5）。
/// **渡さなくても答えは同じ**です。速さだけが変わります。
template <class Point>
inline void winding_split(const mesh::TriMesh& m, const Point& p, const geom::PlaneD& ref,
                          int* w_other, int* c_front, int* c_back, const RaySupport& sup = {}) {
    *w_other = 0;
    *c_front = 0;
    *c_back = 0;
    int sheet_strict = 0, n_strict = 0, sheet_edge = 0;
    bool any_edge = false;

    // **レイの軸は、基準平面の法線が非零な軸から選びます。**
    // レイがシートと平行だと、どちら側が跨ぐかを決められません（実際に踏みました）。
    const geom::Axis along = (arith::sign(ref.a) != 0)   ? geom::Axis::X
                             : (arith::sign(ref.b) != 0) ? geom::Axis::Y
                                                         : geom::Axis::Z;
    KRISITE_CHECK(detail::normal_comp_sign(ref, along) != 0, "winding_split: 基準平面が退化");

    // **候補の絞り込み**（`SPEC-phase5.md` §6.3）。索引があれば、判定点の
    // (u, v) セルにある三角形だけを見ます。**落とすものはありません** —
    // 根拠は `ray_index.hpp` 冒頭の単調性の議論です。
    const RayIndex* idx = sup.index[static_cast<int>(along)];
    if (idx != nullptr && !idx->ready()) idx = nullptr;

    auto visit = [&](std::size_t j) {
        const mesh::Tri& t = m.triangles[j];
        const geom::IPoint& a = m.vertices[t[0]];
        const geom::IPoint& b = m.vertices[t[1]];
        const geom::IPoint& c = m.vertices[t[2]];
        const geom::PlaneD pl =
            (sup.planes != nullptr) ? sup.planes[j] : geom::plane_from_triangle(a, b, c);
        if (geom::is_degenerate(pl)) return;

        // ---- ★ D-2 → D-1 の前判定（`prefilter` が真のときだけ）--------------------
        //
        // **順序は D-2 が先です。** 1 回の比較で 42〜45% を落としてから、
        // D-1 の 2.3〜2.7 回を残りに掛けるほうが安い（`BENCH.md` の実測）。
        if (sup.prefilter) {
            // **D-2**: 三角形が丸ごとレイの手前にあるなら、前方では交わりません。
            // **境界は残します**（`>` であって `>=` ではない）。
            std::int64_t amax = static_cast<std::int64_t>(detail::comp(a, along));
            {
                const std::int64_t b1 = static_cast<std::int64_t>(detail::comp(b, along));
                const std::int64_t c1 = static_cast<std::int64_t>(detail::comp(c, along));
                if (b1 > amax) amax = b1;
                if (c1 > amax) amax = c1;
            }
            if (geom::cmp_axis_int(p, amax, along) > 0) return;
            // **D-1**: 投影した AABB の外なら交わりません。**境界は残します。**
            const geom::Axis uu = detail::proj_u(along), vv = detail::proj_v(along);
            std::int64_t ulo = static_cast<std::int64_t>(detail::comp(a, uu)), uhi = ulo;
            std::int64_t vlo = static_cast<std::int64_t>(detail::comp(a, vv)), vhi = vlo;
            for (const geom::IPoint* qq : {&b, &c}) {
                const std::int64_t cu = static_cast<std::int64_t>(detail::comp(*qq, uu));
                const std::int64_t cv = static_cast<std::int64_t>(detail::comp(*qq, vv));
                if (cu < ulo) ulo = cu;
                if (cu > uhi) uhi = cu;
                if (cv < vlo) vlo = cv;
                if (cv > vhi) vhi = cv;
            }
            if (geom::cmp_axis_int(p, ulo, uu) < 0) return;
            if (geom::cmp_axis_int(p, uhi, uu) > 0) return;
            if (geom::cmp_axis_int(p, vlo, vv) < 0) return;
            if (geom::cmp_axis_int(p, vhi, vv) > 0) return;
        }

        if (sup.kept != nullptr) ++*sup.kept;

        if (sup.fwd_only != nullptr) {
            // **D-2 単独**（D-1 を掛けない）
            std::int64_t amax = static_cast<std::int64_t>(detail::comp(a, along));
            const std::int64_t ab2 = static_cast<std::int64_t>(detail::comp(b, along));
            const std::int64_t ac2 = static_cast<std::int64_t>(detail::comp(c, along));
            if (ab2 > amax) amax = ab2;
            if (ac2 > amax) amax = ac2;
            if (geom::side(geom::plane_axis_aligned(along, amax), p) <= 0) ++*sup.fwd_only;
        }
        if (sup.aabb_pass != nullptr) {
            // **投影した AABB に判定点が入るか**（計測。実装では安い比較で書けます）
            //
            // **早期打ち切りありで、実際に走った比較の回数を数えます。**
            // 費用は「回数 × 1 回の単価」で出るので、時間のばらつきに依りません。
            const geom::Axis uu = detail::proj_u(along), vv = detail::proj_v(along);
            std::int64_t ulo = detail::comp(a, uu), uhi = ulo;
            std::int64_t vlo = detail::comp(a, vv), vhi = vlo;
            for (const geom::IPoint* q : {&b, &c}) {
                const std::int64_t cu = detail::comp(*q, uu), cv = detail::comp(*q, vv);
                ulo = std::min(ulo, cu);
                uhi = std::max(uhi, cu);
                vlo = std::min(vlo, cv);
                vhi = std::max(vhi, cv);
            }
            bool in = true;
            {
                const geom::Axis ax4[4] = {uu, uu, vv, vv};
                const std::int64_t bd4[4] = {ulo, uhi, vlo, vhi};
                const int want[4] = {+1, -1, +1, -1};
                for (int k = 0; k < 4; ++k) {
                    if (sup.cheap_tests != nullptr) ++*sup.cheap_tests;
                    const int sg = geom::side(geom::plane_axis_aligned(ax4[k], bd4[k]), p);
                    if (want[k] > 0 ? (sg < 0) : (sg > 0)) {
                        in = false;
                        break;
                    }
                }
            }
            if (in) {
                ++*sup.aabb_pass;
                // **D-2**: レイは `along` の正方向へ進むので、三角形の `along` 最大座標が
                // 判定点より手前なら、前方では交わりません。**1 回の比較で決まります。**
                if (sup.fwd_pass != nullptr) {
                    std::int64_t amax = static_cast<std::int64_t>(detail::comp(a, along));
                    const std::int64_t ab = static_cast<std::int64_t>(detail::comp(b, along));
                    const std::int64_t ac = static_cast<std::int64_t>(detail::comp(c, along));
                    if (ab > amax) amax = ab;
                    if (ac > amax) amax = ac;
                    if (geom::side(geom::plane_axis_aligned(along, amax), p) <= 0) {
                        ++*sup.fwd_pass;
                    }
                }
            }
        }
        bool on_face = false, strict = false;
        if (geom::side(pl, p) == 0) {
            const geom::Axis ax = (arith::sign(pl.a) != 0)   ? geom::Axis::X
                                  : (arith::sign(pl.b) != 0) ? geom::Axis::Y
                                                             : geom::Axis::Z;
            const int o1 = geom::orient2d_h(a, b, p, ax);
            const int o2 = geom::orient2d_h(b, c, p, ax);
            const int o3 = geom::orient2d_h(c, a, p, ax);
            const bool neg = (o1 < 0) || (o2 < 0) || (o3 < 0);
            const bool pos = (o1 > 0) || (o2 > 0) || (o3 > 0);
            on_face = !(neg && pos);
            strict = (o1 != 0) && (o2 != 0) && (o3 != 0);
        }
        if (!on_face) {
            // **平面を使い回します。** ここで作り直すと 1 三角形あたり 2 回になります
            if (detail::crosses_with(pl, a, b, c, p, along)) {
                *w_other += detail::normal_comp_sign(pl, along);
                if (sup.hits != nullptr) ++*sup.hits;
            }
            return;
        }
        if (sup.hits != nullptr) ++*sup.hits;  // 載っている面も「寄与した」に数えます
        // **載っている面（シート）。** レイはここから出るので、跨ぐかどうかは
        // レイの向きと基準法線の関係で決まります（下の分岐）。
        //
        // **辺や頂点に載っていると、同じシートを複数の三角形が共有します。**
        // 内部に載っているものは 1 枚ずつ、境界だけのものは 1 枚として数えます。
        const int contrib = detail::normal_comp_sign(pl, along);
        if (strict) {
            sheet_strict += contrib;
            ++n_strict;
        } else {
            any_edge = true;
            sheet_edge = contrib;
        }
    };

    std::size_t visited = 0;
    if (idx == nullptr) {
        for (std::size_t j = 0; j < m.triangles.size(); ++j) visit(j);
        visited = m.triangles.size();
    } else {
        // **昇順にマージします。** `sheet_edge` は「最後に見た 1 枚」で決まるので、
        // **全走査と同じ順序でなければ出力が変わり得ます。**
        //
        // 段は互いに素なので、段ごとの区間を k 路マージすれば昇順になります
        // （段は 11 以下。`uint32` の比較 11 回は、述語 1 回より十分安い）。
        const std::uint32_t* cur[kRayIndexMaxLevels];
        const std::uint32_t* end[kRayIndexMaxLevels];
        const std::size_t nl = idx->candidates_levels();
        idx->candidates(p, cur, end);
        for (;;) {
            std::size_t best = nl;
            std::uint32_t bv = 0;
            for (std::size_t l = 0; l < nl; ++l) {
                if (cur[l] == end[l]) continue;
                if (best == nl || *cur[l] < bv) {
                    best = l;
                    bv = *cur[l];
                }
            }
            if (best == nl) break;
            ++cur[best];
            visit(bv);
            ++visited;
            if (sup.per_level != nullptr) ++sup.per_level[best];
        }
    }
    if (sup.tested != nullptr) *sup.tested += visited;

    const int sheet = sheet_strict + ((n_strict == 0 && any_edge) ? sheet_edge : 0);

    // **どちら側がシートを跨ぐか。**
    //
    // 点 $p$ をシートから $\pm N$ に $\varepsilon$ だけ離すと、+X のレイが
    // シートを跨ぐのは「レイが $-N$ 側へ向かうとき」、すなわち $X \cdot N < 0$ のとき。
    //
    //   $X \cdot N < 0$ → 表側（$+N$）のレイがシートを跨ぐ → 表に寄与
    //   $X \cdot N > 0$ → 裏側（$-N$）のレイが跨ぐ → 裏に寄与
    //   $X \cdot N = 0$ → レイがシートと平行 → どちらも跨がない
    const int xn = detail::normal_comp_sign(ref, along);
    if (xn < 0) {
        *c_front = sheet;
    } else if (xn > 0) {
        *c_back = sheet;
    }
}

/// 点が閉じた立体 `m` の内部にあるか。
///
/// **`p` が `m` の境界上に無いことを呼び出し側が保証すること**（`point_on_boundary`）。
template <class Point>
inline bool point_inside(const mesh::TriMesh& m, const Point& p, const RaySupport& sup = {}) {
    int hits = 0;
    std::size_t visited = 0;
    // **レイは +X 固定です。** 索引も X 軸のものだけを使います。
    //
    // ここは偶奇しか見ないので**順序に依存しません**が、`winding_split` と
    // 同じ道具を使うために昇順のマージをそのまま通します。
    const RayIndex* idx = sup.index[static_cast<int>(geom::Axis::X)];
    if (idx != nullptr && !idx->ready()) idx = nullptr;
    const auto visit = [&](std::size_t j) {
        const mesh::Tri& t = m.triangles[j];
        const geom::PlaneD pl =
            (sup.planes != nullptr)
                ? sup.planes[j]
                : geom::plane_from_triangle(m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]]);
        if (detail::crosses_with(pl, m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]], p)) {
            ++hits;
        }
    };
    if (idx == nullptr) {
        for (std::size_t j = 0; j < m.triangles.size(); ++j) visit(j);
        visited = m.triangles.size();
    } else {
        const std::uint32_t* cur[kRayIndexMaxLevels];
        const std::uint32_t* end[kRayIndexMaxLevels];
        const std::size_t nl = idx->candidates_levels();
        idx->candidates(p, cur, end);
        for (;;) {
            std::size_t best = nl;
            std::uint32_t bv = 0;
            for (std::size_t l = 0; l < nl; ++l) {
                if (cur[l] == end[l]) continue;
                if (best == nl || *cur[l] < bv) {
                    best = l;
                    bv = *cur[l];
                }
            }
            if (best == nl) break;
            ++cur[best];
            visit(bv);
            ++visited;
        }
    }
    if (sup.tested != nullptr) *sup.tested += visited;
    return (hits % 2) == 1;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_RAYCAST_HPP

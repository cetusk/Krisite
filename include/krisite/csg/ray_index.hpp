// Krisite — レイキャストの候補を絞る 2 次元階層格子
//
// SPEC-phase5.md §6.3「レイキャストの加速」
//
// **レイは軸平行です。** `winding_split` は基準平面の法線が非零な軸 `along` を選び、
// その正方向にレイを飛ばします。したがってレイは
//
//   「投影面 (u, v) 上の 1 点を通る、`along` に平行な直線」
//
// です。三角形がこのレイと交わり得るのは、**その三角形の (u, v) 投影の AABB が
// 判定点の (u, v) を含むときに限ります。** 一般のレイ–八分木の走査は要りません。
//
// **正しさの根拠は単調性です。** セル番号を返す `cell_of` は座標について単調
// 非減少なので
//
//   p_u ∈ [u_min(t), u_max(t)]  ⟹  cell_of(p_u) ∈ [cell_of(u_min(t)), cell_of(u_max(t))]
//
// が成り立ちます。三角形を「自分の AABB が覆うセル」すべてに入れておけば、
// **判定点のセルを見るだけで、交わり得る三角形をひとつも落としません。**
// これは保守的な絞り込みではなく **厳密な絞り込み**です（誤差も許容も入りません）。
//
// **なぜ一様格子ではなく階層か。** 実データでは三角形の大きさが桁で違います。
// 一様格子に入れると、大きい三角形 1 枚が $R^2$ 個の項目を生みます
// （実測: 2,048 三角形の `474825` で **757,343 項目 = 1 枚あたり 370**）。
// 上限を設けて溢れたものを毎回走査する形にすると、**その模型では 88.9% が
// 溢れて全走査に退化しました。** どちらも成り立ちません。
//
// **階層にすると、三角形は「自分の大きさに合った段」に 4 セル以下で入ります。**
// 項目数は段の合計で $4n$ 以下に収まり、問い合わせは段ごとに 1 セルずつ見るだけです。
//
// 段 $\ell$ のセルは段 0 の $2^\ell$ 倍なので、**位置決めは段 0 の二分探索 1 回で足ります**
// （$\lfloor q/(a 2^\ell)\rfloor = \lfloor \lfloor q/a \rfloor / 2^\ell \rfloor$）。
//
// **判定点は同次座標です。** 位置決めは既存の `cmp_h` で行います
// （w = 1 の同次点として格子線を渡す）。**新しい述語は導入していません。**
#ifndef KRISITE_CSG_RAY_INDEX_HPP
#define KRISITE_CSG_RAY_INDEX_HPP

#include <cstdint>
#include <vector>

#include "krisite/geom/predicates.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace krisite::csg {

namespace detail {

/// 軸の成分を取り出す。
inline std::int32_t comp(const geom::IPoint& p, geom::Axis a) noexcept {
    return (a == geom::Axis::X) ? p.x : (a == geom::Axis::Y) ? p.y : p.z;
}

/// 投影後の 2 軸（`orient2d_h` と同じ巡回順）。
inline geom::Axis proj_u(geom::Axis along) noexcept {
    return (along == geom::Axis::X)   ? geom::Axis::Y
           : (along == geom::Axis::Y) ? geom::Axis::Z
                                      : geom::Axis::X;
}
inline geom::Axis proj_v(geom::Axis along) noexcept {
    return (along == geom::Axis::X)   ? geom::Axis::Z
           : (along == geom::Axis::Y) ? geom::Axis::X
                                      : geom::Axis::Y;
}

/// 判定点の軸成分と整数 `g` の大小。`sign(p_ax - g)` を返します。
inline int cmp_coord(const geom::IPoint& p, std::int64_t g, geom::Axis ax) noexcept {
    const std::int64_t x = comp(p, ax);
    return (x < g) ? -1 : (x > g) ? 1 : 0;
}

/// 同次点の版。**既存の `cmp_h` を使います**（w = 1 の同次点との比較）。
///
/// `cmp_h` は指定した軸の成分と `w` しか読まないので、他の成分は 0 のままで
/// 構いません。**1 レイあたり $O(\log R)$ 回**しか呼ばれないので、
/// `cmp_h`（13b+27 ビット）の重さは効きません。
inline int cmp_coord(const geom::HPointD& p, std::int64_t g, geom::Axis ax) noexcept {
    geom::HPointD q{};
    q.w = arith::from_i64<geom::limbs::kHomoW>(1);
    (ax == geom::Axis::X   ? q.x
     : ax == geom::Axis::Y ? q.y
                           : q.z) = arith::from_i64<geom::limbs::kHomoXyz>(g);
    return geom::cmp_h(p, q, ax);
}

}  // namespace detail

/// `along` に垂直な投影面に張った 2 次元階層格子。
///
/// **1 つの軸につき 1 つ必要です。** `along` は基準平面ごとに変わるので、
/// 呼び出し側は 3 軸ぶん作ります（構築は $O(n)$ なので無視できます）。
class RayIndex {
public:
    /// 三角形を入れる段は「覆うセルがこれ以下になる最小の段」です。
    ///
    /// **各軸 2 セル以下 = 2×2。** どの段でも 1 枚あたり 4 項目以下なので、
    /// **全項目数は $4n$ 以下**に収まります（段は互いに素）。

    void build(const mesh::TriMesh& m, geom::Axis along) {
        along_ = along;
        u_ = detail::proj_u(along);
        v_ = detail::proj_v(along);
        levels_.clear();
        n_tri_ = m.triangles.size();
        if (n_tri_ == 0) {
            res0_ = 0;
            return;
        }

        std::int64_t umin = INT64_MAX, umax = INT64_MIN, vmin = INT64_MAX, vmax = INT64_MIN;
        for (const geom::IPoint& p : m.vertices) {
            const std::int64_t cu = detail::comp(p, u_), cv = detail::comp(p, v_);
            umin = (cu < umin) ? cu : umin;
            umax = (cu > umax) ? cu : umax;
            vmin = (cv < vmin) ? cv : vmin;
            vmax = (cv > vmax) ? cv : vmax;
        }
        u0_ = umin;
        v0_ = vmin;

        // **最下段は 1 セルあたり三角形 1 枚を狙って $R = \lceil\sqrt{n}\rceil$。**
        std::uint32_t r = 1;
        while (static_cast<std::size_t>(r) * r < n_tri_ && r < 1024) ++r;
        res0_ = r;
        du0_ = span_step(umax - umin, r);
        dv0_ = span_step(vmax - vmin, r);

        // 段の数: 解像度が 1 になるまで
        std::uint32_t n_lv = 1;
        while ((res0_ >> n_lv) != 0) ++n_lv;
        levels_.resize(n_lv);
        for (std::uint32_t l = 0; l < n_lv; ++l) {
            levels_[l].res = static_cast<std::uint32_t>((res0_ + (1u << l) - 1) >> l);
            if (levels_[l].res == 0) levels_[l].res = 1;
            levels_[l].offsets.assign(static_cast<std::size_t>(levels_[l].res) * levels_[l].res + 1,
                                      0);
        }

        // 各三角形の段と、段 0 でのセル範囲
        std::vector<std::uint32_t> lv(n_tri_), lu(n_tri_), hu(n_tri_), lvv(n_tri_), hv(n_tri_);
        for (std::size_t j = 0; j < n_tri_; ++j) {
            std::int64_t e[4];
            tri_extent(m, j, e);
            const std::uint32_t a = cell0(e[0], u0_, du0_), b = cell0(e[1], u0_, du0_);
            const std::uint32_t c = cell0(e[2], v0_, dv0_), d = cell0(e[3], v0_, dv0_);
            lu[j] = a;
            hu[j] = b;
            lvv[j] = c;
            hv[j] = d;
            // **各軸で 2 セル以下**になる最小の段。積で 4 以下にすると
            // 1×4 のような並びが通ってしまい、**四隅だけでは中が抜けます。**
            std::uint32_t l = 0;
            for (; l + 1 < n_lv; ++l) {
                if (((b >> l) - (a >> l)) <= 1 && ((d >> l) - (c >> l)) <= 1) break;
            }
            lv[j] = l;
            ++levels_[l].offsets[cell_at(l, a >> l, c >> l) + 1];
            if ((b >> l) != (a >> l)) ++levels_[l].offsets[cell_at(l, b >> l, c >> l) + 1];
            if ((d >> l) != (c >> l)) ++levels_[l].offsets[cell_at(l, a >> l, d >> l) + 1];
            if ((b >> l) != (a >> l) && (d >> l) != (c >> l))
                ++levels_[l].offsets[cell_at(l, b >> l, d >> l) + 1];
        }
        for (Level& L : levels_) {
            for (std::size_t i = 0; i + 1 < L.offsets.size(); ++i) L.offsets[i + 1] += L.offsets[i];
            L.items.resize(L.offsets.back());
            L.fill.assign(L.offsets.begin(), L.offsets.end() - 1);
        }
        for (std::size_t j = 0; j < n_tri_; ++j) {
            const std::uint32_t l = lv[j];
            const std::uint32_t a = lu[j] >> l, b = hu[j] >> l, c = lvv[j] >> l, d = hv[j] >> l;
            Level& L = levels_[l];
            const auto put = [&](std::uint32_t ku, std::uint32_t kv) {
                L.items[L.fill[cell_at(l, ku, kv)]++] = static_cast<std::uint32_t>(j);
            };
            put(a, c);
            if (b != a) put(b, c);
            if (d != c) put(a, d);
            if (b != a && d != c) put(b, d);
        }
        for (Level& L : levels_) L.fill.clear();
    }

    bool ready() const noexcept { return res0_ != 0; }
    std::size_t levels() const noexcept { return levels_.size(); }
    /// `candidates` が埋める区間の数（段の数）。
    std::size_t candidates_levels() const noexcept { return levels_.size(); }
    std::size_t items() const noexcept {
        std::size_t n = 0;
        for (const Level& L : levels_) n += L.items.size();
        return n;
    }
    /// 段 `l` の項目数（計測用。`DESIGN-phase5-hotspots.md` §11）。
    /// **大きい三角形は粗い段に入ります。**段ごとの分布が偏りを示します。
    std::size_t items_at(std::size_t l) const noexcept { return levels_[l].items.size(); }
    /// 段 `l` の解像度（1 辺のセル数）。
    std::uint32_t res_at(std::size_t l) const noexcept { return levels_[l].res; }

    /// 判定点 `p` の候補（**三角形番号の昇順**）を `out` に集めます。
    ///
    /// **昇順であることが要ります。** `winding_split` は接触面の寄与を
    /// 「最後に見た 1 枚」で決める分岐を持つので、**走査順が変わると出力が
    /// 変わり得ます**（全走査と同じ順序であることが不変性の条件）。
    ///
    /// 段は互いに素なので、段ごとの区間を k 路マージするだけで昇順になります。
    template <class Point>
    void candidates(const Point& p, const std::uint32_t** first, const std::uint32_t** last) const {
        const std::uint32_t ku0 = locate(p, u_, u0_, du0_);
        const std::uint32_t kv0 = locate(p, v_, v0_, dv0_);
        for (std::uint32_t l = 0; l < levels_.size(); ++l) {
            const Level& L = levels_[l];
            const std::size_t c = cell_at(l, ku0 >> l, kv0 >> l);
            first[l] = L.items.data() + L.offsets[c];
            last[l] = L.items.data() + L.offsets[c + 1];
        }
    }

private:
    struct Level {
        std::uint32_t res = 1;
        std::vector<std::uint32_t> offsets, items, fill;
    };

    static std::int64_t span_step(std::int64_t span, std::uint32_t r) noexcept {
        const std::int64_t s = (span + r - 1) / r;
        return (s < 1) ? 1 : s;
    }
    std::size_t cell_at(std::uint32_t l, std::uint32_t ku, std::uint32_t kv) const noexcept {
        return static_cast<std::size_t>(kv) * levels_[l].res + ku;
    }
    /// **単調非減少**。これが絞り込みの正しさの根拠です。
    std::uint32_t cell0(std::int64_t x, std::int64_t origin, std::int64_t step) const noexcept {
        if (x <= origin) return 0;
        const std::uint64_t k =
            static_cast<std::uint64_t>(x - origin) / static_cast<std::uint64_t>(step);
        return (k >= res0_) ? (res0_ - 1) : static_cast<std::uint32_t>(k);
    }
    /// 有理点の版。`cell0` と同じ値を二分探索で求めます（**1 軸あたり 1 回**）。
    template <class Point>
    std::uint32_t locate(const Point& p, geom::Axis ax, std::int64_t origin,
                         std::int64_t step) const {
        std::uint32_t lo = 0, hi = res0_ - 1;
        while (lo < hi) {
            const std::uint32_t mid = lo + (hi - lo + 1) / 2;
            if (detail::cmp_coord(p, origin + static_cast<std::int64_t>(mid) * step, ax) >= 0)
                lo = mid;
            else
                hi = mid - 1;
        }
        return lo;
    }
    void tri_extent(const mesh::TriMesh& m, std::size_t j, std::int64_t* out) const noexcept {
        const mesh::Tri& t = m.triangles[j];
        std::int64_t umin = INT64_MAX, umax = INT64_MIN, vmin = INT64_MAX, vmax = INT64_MIN;
        for (std::size_t k = 0; k < 3; ++k) {
            const geom::IPoint& p = m.vertices[t[k]];
            const std::int64_t cu = detail::comp(p, u_), cv = detail::comp(p, v_);
            umin = (cu < umin) ? cu : umin;
            umax = (cu > umax) ? cu : umax;
            vmin = (cv < vmin) ? cv : vmin;
            vmax = (cv > vmax) ? cv : vmax;
        }
        out[0] = umin;
        out[1] = umax;
        out[2] = vmin;
        out[3] = vmax;
    }

    geom::Axis along_ = geom::Axis::X, u_ = geom::Axis::Y, v_ = geom::Axis::Z;
    std::int64_t u0_ = 0, v0_ = 0, du0_ = 1, dv0_ = 1;
    std::uint32_t res0_ = 0;
    std::size_t n_tri_ = 0;
    std::vector<Level> levels_;
};

/// 段の数の上限（`candidates` に渡す配列の大きさ）。$R \le 1024$ なので 11 で足ります。
inline constexpr std::size_t kRayIndexMaxLevels = 12;

}  // namespace krisite::csg

#endif  // KRISITE_CSG_RAY_INDEX_HPP

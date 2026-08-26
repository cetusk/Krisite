// Krisite — 平面 ID の登録
//
// SPEC-phase1.md §3.1
//
// **GCD による正準化は行いません。** 除算・剰余・GCD は `arith/` に無く、Phase 1 の
// ために足すのは非目標です（§12）。代わりに `plane_cmp` の全順序で整列し、
// 隣接併合で ID を割り当てます。
//
// 符号正規化で平面から向きの情報が消えるので、**登録時に「表か裏か」を返します**。
// 三角形側がこのフラグを持ちます（§3.1）。
#ifndef KRISITE_CSG_PLANE_TABLE_HPP
#define KRISITE_CSG_PLANE_TABLE_HPP

#include <algorithm>
#include <cstdint>
#include <vector>

#include "krisite/geom/plane.hpp"
#include "krisite/geom/predicates.hpp"

namespace krisite::csg {

using PlaneId = std::uint32_t;
inline constexpr PlaneId kNoPlane = static_cast<PlaneId>(-1);

/// 平面の登録結果。
struct PlaneRef {
    PlaneId id = kNoPlane;
    /// 登録した平面の法線が、表に格納された代表の法線と**逆向き**か。
    bool flipped = false;
};

/// 幾何平面 → `PlaneId` の対応表。同一平面（向きを問わない）は同一 ID になります。
class PlaneTable {
public:
    /// 平面を登録して ID を得る。既に等価な平面があればその ID を返す。
    PlaneRef intern(const geom::PlaneD& pl) {
        KRISITE_CHECK(!geom::is_null(pl), "PlaneTable: 退化平面は登録できない");
        // order_ は plane_cmp 昇順の ID 列。二分探索する
        auto it = std::lower_bound(order_.begin(), order_.end(), pl,
                                   [this](PlaneId id, const geom::PlaneD& q) {
                                       return geom::plane_cmp(planes_[id], q) < 0;
                                   });
        if (it != order_.end() && geom::plane_cmp(planes_[*it], pl) == 0) {
            return PlaneRef{*it, orientation_differs(planes_[*it], pl)};
        }
        const auto id = static_cast<PlaneId>(planes_.size());
        planes_.push_back(pl);
        order_.insert(it, id);
        return PlaneRef{id, false};
    }

    /// 既に登録済みの平面を探す。無ければ `kNoPlane`。
    PlaneRef find(const geom::PlaneD& pl) const {
        auto it = std::lower_bound(order_.begin(), order_.end(), pl,
                                   [this](PlaneId id, const geom::PlaneD& q) {
                                       return geom::plane_cmp(planes_[id], q) < 0;
                                   });
        if (it != order_.end() && geom::plane_cmp(planes_[*it], pl) == 0) {
            return PlaneRef{*it, orientation_differs(planes_[*it], pl)};
        }
        return PlaneRef{};
    }

    const geom::PlaneD& at(PlaneId id) const noexcept { return planes_[id]; }
    std::size_t size() const noexcept { return planes_.size(); }

private:
    /// 比例する 2 平面の向きが逆か。先頭非零成分の符号で決まる。
    static bool orientation_differs(const geom::PlaneD& rep, const geom::PlaneD& pl) noexcept {
        for (int i = 0; i < 4; ++i) {
            const int sr = arith::sign(geom::detail::plane_coeff(rep, i));
            if (sr != 0) {
                const int sp = arith::sign(geom::detail::plane_coeff(pl, i));
                return sr != sp;
            }
        }
        return false;
    }

    std::vector<geom::PlaneD> planes_;
    std::vector<PlaneId> order_;
};

}  // namespace krisite::csg

#endif  // KRISITE_CSG_PLANE_TABLE_HPP

// Krisite — 構成点の保持（メモ化）
//
// SPEC-phase2.md §4
//
// **Phase 1 の実測が示した最優先項目です。** `side` : `intersect3` の呼び出し比が
// **1.26 : 1** で、`intersect3` が述語時間の 94.9%（推定）を占めていました
// （`BENCH.md` §12、`IMPL-phase1.md` §5.6）。
//
// 原因は `split_fragment` が
//
//     s[i] = geom::side(qp, fragment_vertex(t, f, i));
//
// の形で、**述語を評価するたびに構成点を作り直している**ことです。
// **単価を下げても、回数が `side` と同数のままでは構図が変わりません。**
//
// **構成点は正規化した平面3つ組で一意に決まります**（§3.3、§5.1）。これをキーに
// メモ化します。
//
// ---
//
// **キャッシュはグローバルに持ちません**（§4.2）。`STYLE.md` の算術コードの制約
// （グローバル変数・可変な静的変数の禁止）と、Phase 3 の並列化のためです。
// **呼び出し側が持ち、明示的に引き回します。**
//
// **4 平面以上が一点で交わる場合、異なる3つ組が同じ値を持ちます**（`SPEC-phase1.md` §5.2）。
// **これはキャッシュの誤りではありません。** キャッシュは3つ組で引き、値の同一性は
// §5.3 の第2段が別途扱う、という分担を崩さないこと（§4.3）。
#ifndef KRISITE_CSG_POINT_CACHE_HPP
#define KRISITE_CSG_POINT_CACHE_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>

#include "krisite/csg/plane_table.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/geom/point.hpp"

namespace krisite::csg {

/// 平面3つ組 → 構成点のメモ化（§4.2）。
///
/// **キーは昇順に正規化した3つ組です。** `intersect3` は引数の順序で符号が変わりますが
/// （行列式の置換）、$V$ と $-V$ は同じ射影点で、すべての述語が
/// `sign(w)` を掛ける形になっているため**符号は結果に影響しません。**
/// それでも**常に昇順で評価**します。そうすればキャッシュの有無で出力が 1 ビットも
/// 変わらず、§9.1 の比較が最も強い形（値の完全一致）で成立します。
class PointCache {
public:
    using Key = std::array<PlaneId, 3>;

    /// 3 平面の交点。**キャッシュに無ければ計算して覚えます。**
    const geom::HPointD& get(const PlaneTable& table, PlaneId a, PlaneId b, PlaneId c) {
        Key k{a, b, c};
        std::sort(k.begin(), k.end());
#if defined(KRISITE_MUTATION_CACHE_KEY_DROP)
        // SPEC-phase2 §9.3 の変異 7: キーから 1 平面を落とす。
        // **異なる3つ組が同じキーに落ち、別の点を同一視します。**
        // 静かに壊れる種類の誤りなので、位相・体積の両方で突く必要があります。
        k[2] = k[0];
#endif
        auto it = map_.find(k);
        if (it != map_.end()) {
            ++hits_;
            return it->second;
        }
        ++misses_;
        return map_.emplace(k, geom::intersect3(table.at(k[0]), table.at(k[1]), table.at(k[2])))
            .first->second;
    }

    std::size_t hits() const noexcept { return hits_; }
    std::size_t misses() const noexcept { return misses_; }
    std::size_t entries() const noexcept { return map_.size(); }

    /// 占有メモリの概算（§4.4 の記録）。**`HPoint` は 7 リム規模**（`BENCH.md`）なので
    /// 構成点が増えると効きます。ノードの管理領域は含みません。
    std::size_t bytes() const noexcept {
        return map_.size() * (sizeof(Key) + sizeof(geom::HPointD));
    }

private:
    std::map<Key, geom::HPointD> map_;
    std::size_t hits_ = 0;
    std::size_t misses_ = 0;
};

}  // namespace krisite::csg

#endif  // KRISITE_CSG_POINT_CACHE_HPP

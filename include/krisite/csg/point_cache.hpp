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
#include <cstdint>
#include <map>
#include <memory>
#include <utility>
#include <vector>

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

    PointCache() = default;
    /// **`use_map` が真なら `std::map` を使います**（`SPEC-phase5.md` §5.10.12 の正解器）。
    ///
    /// **旗の ON / OFF で出力がバイト一致することを検査するために残しています。**
    /// **`CLAUDE.md`「従来の経路を旗で残し、両者が一致することを検査してください。
    /// 残さないと、検査を弱めていないことを自分で示せません」。**
    explicit PointCache(bool use_map) : use_map_(use_map) {}

    /// 3 平面の交点。**キャッシュに無ければ計算して覚えます。**
    ///
    /// > **★ 返す参照は、その後の `get` で無効になりません**（§5.10.12）。
    /// > **開番地法の表は伸びますが、点そのものは【塊に分けて】確保しており、
    /// > 一度置いた場所から動きません。**
    const geom::HPointD& get(const PlaneTable& table, PlaneId a, PlaneId b, PlaneId c) {
        Key k{a, b, c};
        std::sort(k.begin(), k.end());
#if defined(KRISITE_MUTATION_CACHE_KEY_DROP)
        // SPEC-phase2 §9.3 の変異 7: キーから 1 平面を落とす。
        // **異なる3つ組が同じキーに落ち、別の点を同一視します。**
        // 静かに壊れる種類の誤りなので、位相・体積の両方で突く必要があります。
        k[2] = k[0];
#endif
        if (use_map_) {
            auto it = map_.find(k);
            if (it != map_.end()) {
                ++hits_;
                ++it->second.hits;
                return it->second.point;
            }
            ++misses_;
            return map_
                .emplace(k,
                         Entry{geom::intersect3(table.at(k[0]), table.at(k[1]), table.at(k[2])), 0})
                .first->second.point;
        }
        // ---- 開番地法（線形探索）------------------------------------------------
        if (buckets_.empty()) rehash(64);
        const std::size_t slot = probe(k);
        if (buckets_[slot].idx != kEmpty) {
            ++hits_;
            Entry& e = at(buckets_[slot].idx);
            ++e.hits;
            return e.point;
        }
        ++misses_;
        const auto idx = static_cast<std::uint32_t>(count_);
        push(Entry{geom::intersect3(table.at(k[0]), table.at(k[1]), table.at(k[2])), 0});
        buckets_[slot].key = k;
        buckets_[slot].idx = idx;
        // **占有率が 0.7 を超えたら倍にします。** 線形探索は詰まると急に遅くなります。
        if ((count_ + 1) * 10 > buckets_.size() * 7) rehash(buckets_.size() * 2);
        return at(idx).point;
    }

    /// **この3つ組が 1 度でも「覚えていた側」から返されたか**（SPEC-phase2 §13 の CP5）。
    ///
    /// CP5 の「構成点の保持 × T 解決」は、**保持された点がそのまま T 頂点として
    /// 挿入される経路**です。両方の機構が同じ実行で動いただけでは足りないので、
    /// **同じ点が両方を通ったこと**を数えられるようにします。
    bool served_from_cache(const Key& k) const noexcept {
        if (use_map_) {
            auto it = map_.find(k);
            return it != map_.end() && it->second.hits > 0;
        }
        if (buckets_.empty()) return false;
        const std::size_t slot = probe(k);
        return buckets_[slot].idx != kEmpty && at(buckets_[slot].idx).hits > 0;
    }

    std::size_t hits() const noexcept { return hits_; }
    std::size_t misses() const noexcept { return misses_; }
    std::size_t entries() const noexcept { return use_map_ ? map_.size() : count_; }

    /// 占有メモリの概算（§4.4 の記録）。**`HPoint` は 7 リム規模**（`BENCH.md`）なので
    /// 構成点が増えると効きます。ノードの管理領域は含みません。
    std::size_t bytes() const noexcept { return entries() * (sizeof(Key) + sizeof(Entry)); }

private:
    /// 点と、**その3つ組が覚えていた側から返された回数**。
    struct Entry {
        geom::HPointD point;
        std::uint32_t hits;
    };

    /// **表の枠**。鍵を枠に置くので、探索の比較で点まで辿りません（16 バイト）。
    struct Bucket {
        Key key{};
        std::uint32_t idx = kEmptyInit;
    };

    static constexpr std::uint32_t kEmpty = 0xFFFFFFFFu;
    static constexpr std::uint32_t kEmptyInit = 0xFFFFFFFFu;
    /// **点は塊に分けて確保します。** 表が伸びても、返した参照が無効になりません。
    static constexpr std::size_t kChunk = 512;

    static std::uint64_t hash_key(const Key& k) noexcept {
        std::uint64_t h = (static_cast<std::uint64_t>(k[0]) << 32) ^ k[1];
        h *= 0x9E3779B97F4A7C15ull;
        h ^= h >> 29;
        h += k[2];
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 32;
        return h;
    }

    /// 鍵の枠を探す。**見つからなければ空き枠を返します。**
    std::size_t probe(const Key& k) const noexcept {
        const std::size_t mask = buckets_.size() - 1;
        std::size_t i = static_cast<std::size_t>(hash_key(k)) & mask;
        while (buckets_[i].idx != kEmpty && buckets_[i].key != k) i = (i + 1) & mask;
        return i;
    }

    Entry& at(std::size_t i) noexcept { return chunks_[i / kChunk][i % kChunk]; }
    const Entry& at(std::size_t i) const noexcept { return chunks_[i / kChunk][i % kChunk]; }

    void push(Entry&& e) {
        if (count_ % kChunk == 0) chunks_.push_back(std::make_unique<Entry[]>(kChunk));
        chunks_[count_ / kChunk][count_ % kChunk] = std::move(e);
        ++count_;
    }

    void rehash(std::size_t cap) {
        std::vector<Bucket> nb(cap);
        for (const Bucket& old : buckets_) {
            if (old.idx == kEmpty) continue;
            const std::size_t mask = cap - 1;
            std::size_t i = static_cast<std::size_t>(hash_key(old.key)) & mask;
            while (nb[i].idx != kEmpty) i = (i + 1) & mask;
            nb[i] = old;
        }
        buckets_.swap(nb);
    }

    bool use_map_ = false;
    std::vector<Bucket> buckets_;
    std::vector<std::unique_ptr<Entry[]>> chunks_;
    std::size_t count_ = 0;
    std::map<Key, Entry> map_;  ///< `use_map_` のときだけ使う正解器
    std::size_t hits_ = 0;
    std::size_t misses_ = 0;
};

}  // namespace krisite::csg

#endif  // KRISITE_CSG_POINT_CACHE_HPP

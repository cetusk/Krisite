// Krisite — 断片（凸多角形）と平面による分割
//
// SPEC-phase1.md §4.3
//
// 断片は支持平面と辺の平面だけで表します。頂点は
//   頂点 i = support ∩ edge[(i-1+n)%n] ∩ edge[i]
// で復元します。**クリップで生まれる頂点も 3 平面の交点になる**ので、
// 表現がクリップについて閉じます。
//
// 交差線分ではなく**平面全体**で分割します（§4.3）。凸多角形を半平面で切ると
// 凸のままなので、制約付き三角形分割が不要になり扇状三角形化で済みます。
#ifndef KRISITE_CSG_FRAGMENT_HPP
#define KRISITE_CSG_FRAGMENT_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include "krisite/csg/faces.hpp"
#include "krisite/csg/plane_table.hpp"
#include "krisite/csg/point_cache.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/geom/predicates.hpp"

namespace krisite::csg {

/// 凸多角形の断片。
struct Fragment {
    PlaneId support = kNoPlane;
    bool flipped = false;  ///< 外向き法線が support の法線と逆か
    std::vector<PlaneId> edge;
    int owner = 0;  ///< 0 = A, 1 = B
    /// **切断の符号列**（`DESIGN-phase5-hotspots.md` §14）。
    ///
    /// **共平面重複の仕分け（`region_key`）の鍵です。**
    /// **これがあるので、中核は頂点 ID を作らずに済みます**（縫合が要らない）。
    ///
    /// **鍵は `(セル, 支持平面, この符号列)`。** セルが要る理由は §14.7.2 —
    /// **第 2 段の切断（共平面揃え）がセルごとに決まる**ので、
    /// **符号だけでは違うセルの領域が混ざります**（実測で「粗すぎ」が全件）。
    ///
    /// **1 切断 = 1 ビット。** 実測で断片あたり 13〜117 ビット（最大 328）。
    std::vector<std::uint64_t> cutbits;
    /// **符号列のビット数。** `cutbits` の長さだけでは区別できません
    /// （65 ビットと 128 ビットは、どちらも 2 ワード）
    std::uint32_t ncuts = 0;
#if defined(KRISITE_EXPERIMENT_REGION_HIST)
    /// **切断の履歴**（`RESEARCH-perf.md` §S3.5 の実験。**実験用です**）。
    ///
    /// **問い**: 共平面重複の仕分け（`region_key`）を、
    /// **頂点 ID ではなく切断の符号列**で行えるか。
    /// できれば**中核の縫合（実測で 29.2%）が消えます。**
    ///
    /// > **`edge` だけでは足りません。** `split_fragment` は
    /// > **pos 側にも neg 側にも切断平面を `edge` に足す**ので、
    /// > **切断の両側が同じ辺平面の集合を持ちます。符号が要ります。**
    ///
    /// **コンパイル時の旗で囲んであります。** 既定では
    /// **`Fragment` にこのメンバは存在せず、メモリは 1 バイトも増えません。**
    /// **実験が終わったら、鍵として残すか消すかを別に判断します。**
    std::vector<std::pair<PlaneId, std::int8_t>> hist;
#endif
};

/// 頂点数 = 辺数。
inline std::size_t vertex_count(const Fragment& f) noexcept {
    return f.edge.size();
}

/// 頂点 i = support ∩ edge[(i-1+n)%n] ∩ edge[i]。
///
/// **`cache` を渡すとメモ化されます**（SPEC-phase2 §4）。渡さなければ Phase 1 の挙動
/// そのもので、**それが §9.1 の正解器**です。
///
/// **キャッシュの有無で戻り値は 1 ビットも変わりません。** どちらも平面3つ組を昇順に
/// 並べてから `intersect3` を呼ぶためです（`point_cache.hpp` の但し書き）。
inline geom::HPointD fragment_vertex(const PlaneTable& t, const Fragment& f, std::size_t i,
                                     PointCache* cache = nullptr) {
    const std::size_t n = f.edge.size();
    const PlaneId a = f.support, b = f.edge[(i + n - 1) % n], c = f.edge[i];
    if (cache != nullptr) return cache->get(t, a, b, c);
    std::array<PlaneId, 3> k{a, b, c};
    std::sort(k.begin(), k.end());
    return geom::intersect3(t.at(k[0]), t.at(k[1]), t.at(k[2]));
}

namespace detail {

/// **切断の符号を 1 ビット足す**（`DESIGN-phase5-hotspots.md` §14）。
///
/// **「切られなかった」ことも履歴の一部です**（§14.1）。
/// 全頂点が片側にある場合も符号を記録しないと、
/// **相手の多角形が切られたときに、同じ領域が違う鍵になります。**
inline void push_cut(Fragment& f, bool positive) noexcept {
    const std::uint32_t i = f.ncuts++;
    if ((i >> 6) >= f.cutbits.size()) f.cutbits.push_back(0);
    if (positive) f.cutbits[i >> 6] |= (std::uint64_t{1} << (i & 63));
}

}  // namespace detail

/// 断片を平面 `q` で分割した結果。存在しない側は `edge` が空。
struct SplitResult {
    Fragment pos;  ///< side(q, ·) >= 0 の側
    Fragment neg;  ///< side(q, ·) <= 0 の側
    bool has_pos = false;
    bool has_neg = false;
};

namespace detail {

/// 保持する側（k = +1 か -1）の辺列を作る。`s[i]` は頂点 i の side。
///
/// 凸多角形なので、保持される辺は巡回的に連続します。その連続区間の末尾に
/// 切断平面 `q` を 1 枚だけ足せば閉じます。
///
/// **開始位置は「境界が切断線から戻ってくる辺」です。**
///
/// 以前は「直前が非保持である最初の保持辺」で探していましたが、**全辺が保持される
/// 場合に開始点が見つからず、その側を丸ごと落としていました。** 三角形を 1 頂点だけ
/// 外側に切るとこれが起きます。符号が $(+,+,-)$ のとき、3 辺はいずれも
/// 「端点の少なくとも一方が真に内側」を満たすためです。残るべき四角形が消えます。
///
/// **軸平行な立方体を軸平行な平面で切るかぎり断片は常に長方形なので、この配置は
/// 一度も現れません。** 斜面を持つ入力（四面体、回転立方体）で初めて出ます。
///
/// そこで開始位置を、**脱出辺**（内側から切断線へ出る辺。凸多角形では一意）の
/// 次にある最初の保持辺として求めます。全辺が保持される場合も含めて一様に決まります。
inline std::vector<PlaneId> clip_edges(const std::vector<PlaneId>& edge, const std::vector<int>& s,
                                       int k, PlaneId q) {
    const std::size_t n = edge.size();
    std::vector<char> keep(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const int a = s[i] * k, b = s[(i + 1) % n] * k;
        // 辺の保持部分が正の長さを持つのは、端点の少なくとも一方が真に内側のとき
        keep[i] = (a > 0 || b > 0) ? 1 : 0;
    }
    // 脱出辺: 始点が真に内側で、終点が内側でない辺。凸多角形では高々 1 本
    std::size_t exit_edge = n;
    for (std::size_t i = 0; i < n; ++i) {
        if (s[i] * k > 0 && s[(i + 1) % n] * k <= 0) {
            exit_edge = i;
            break;
        }
    }
    if (exit_edge == n) return {};  // 真に内側の頂点が無い

    // 脱出辺の次から巡回して、最初の保持辺を開始位置にする
    std::size_t start = n;
    for (std::size_t j = 1; j <= n; ++j) {
        const std::size_t i = (exit_edge + j) % n;
        if (keep[i]) {
            start = i;
            break;
        }
    }
    if (start == n) return {};  // 保持なし（全部落ちた）

    std::vector<PlaneId> out;
    for (std::size_t j = 0; j < n; ++j) {
        const std::size_t i = (start + j) % n;
        if (!keep[i]) break;
        out.push_back(edge[i]);
    }
    out.push_back(q);  // 保持区間の末尾（= 脱出辺の直後）で切断平面が閉じる
    return out;
}

}  // namespace detail

/// 断片を平面 `q`（ID）で分割する。`q == f.support` のときは分割しません。
inline SplitResult split_fragment(const PlaneTable& t, const Fragment& f, PlaneId q,
                                  PointCache* cache = nullptr) {
    SplitResult r;
    if (q == f.support) {
        r.pos = f;
        r.has_pos = true;
        return r;
    }
    const std::size_t n = f.edge.size();
    const geom::PlaneD& qp = t.at(q);
    std::vector<int> s(n);
    bool any_pos = false, any_neg = false;
    // **切断の履歴**（実験用。上記）。**`q == f.support` の場合は追記しません** —
    // その分岐は上で早期に返しており、**共平面の 2 つの多角形で同じ扱い**なので、
    // 鍵の比較には影響しません。
    for (std::size_t i = 0; i < n; ++i) {
        s[i] = geom::side(qp, fragment_vertex(t, f, i, cache));
        if (s[i] > 0) any_pos = true;
        if (s[i] < 0) any_neg = true;
    }
    if (!any_neg) {  // すべて >= 0
        r.pos = f;
        r.has_pos = true;
        detail::push_cut(r.pos, true);
#if defined(KRISITE_EXPERIMENT_REGION_HIST)
        // **★ 切断が効かなくても符号は記録します。**
        //
        // **相手の多角形はこの平面で切られるかもしれません。**
        // そのとき相手の pos 側は `(q,+1)` を持つので、
        // **こちらが持たないと、同じ領域が違う鍵になります。**
        r.pos.hist.emplace_back(q, static_cast<std::int8_t>(1));
#endif
        return r;
    }
    if (!any_pos) {  // すべて <= 0
        r.neg = f;
        r.has_neg = true;
        detail::push_cut(r.neg, false);
#if defined(KRISITE_EXPERIMENT_REGION_HIST)
        r.neg.hist.emplace_back(q, static_cast<std::int8_t>(-1));
#endif
        return r;
    }

    auto make = [&](int k, Fragment& out) {
        std::vector<PlaneId> e = detail::clip_edges(f.edge, s, k, q);
        if (e.size() < 3) return false;
        out.support = f.support;
        out.flipped = f.flipped;
        out.owner = f.owner;
        out.edge = std::move(e);
        out.cutbits = f.cutbits;
        out.ncuts = f.ncuts;
        detail::push_cut(out, k > 0);
#if defined(KRISITE_EXPERIMENT_REGION_HIST)
        out.hist = f.hist;
        out.hist.emplace_back(q, static_cast<std::int8_t>(k));
#endif
        return true;
    };
    r.has_pos = make(+1, r.pos);
    r.has_neg = make(-1, r.neg);
    return r;
}

/// 断片を半平面 `side(q, ·) * k >= 0` にクリップする（片側だけ残す）。
inline bool clip_fragment(const PlaneTable& t, Fragment& f, PlaneId q, int k,
                          PointCache* cache = nullptr) {
    const SplitResult r = split_fragment(t, f, q, cache);
    if (k > 0) {
        if (!r.has_pos) return false;
        f = r.pos;
    } else {
        if (!r.has_neg) return false;
        f = r.neg;
    }
    return true;
}

/// 断片の平面 `q` に対する符号。頂点の符号の総意をとる（0 でない値が一意に定まる）。
///
/// 分割後の断片は `q` の片側に収まっているので、0 でない符号は高々 1 種類です。
/// すべて 0 なら断片は `q` 上に載っています（`q == support` のとき）。
inline int fragment_sign(const PlaneTable& t, const Fragment& f, PlaneId q,
                         PointCache* cache = nullptr) {
    const std::size_t n = f.edge.size();
    int acc = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const int s = geom::side(t.at(q), fragment_vertex(t, f, i, cache));
        if (s != 0) {
            KRISITE_CHECK(acc == 0 || acc == s, "fragment_sign: 断片が平面をまたいでいる");
            acc = s;
        }
    }
    return acc;
}

/// 面をそのまま断片にする。
inline Fragment face_to_fragment(const Face& f) {
    Fragment r;
    r.support = f.support;
    r.flipped = f.flipped;
    r.owner = f.owner;
    r.edge = f.edge;
    return r;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_FRAGMENT_HPP

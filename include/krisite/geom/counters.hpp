// Krisite — 述語の呼び出し回数の計測（測定専用。既定ビルドでは無効）
//
// SPEC-phase1.md §12: 「`side` と `intersect3` の呼び出し比率を計測してください。
// 最適化はしませんが、比率は Phase 2 の判断材料です」
//
// IMPL-phase0 §7 が問うているのは「`intersect3` が 262 ns で突出して重いが、
// 構成点は一度作って何度も `side` にかけるはずなので、支配的とは限らない」でした。
// **比率を知らずに重い関数を最適化しても無駄になります。**
//
// **既定ビルドでは何も起きません。** `KRISITE_COUNT_PREDICATES` を定義した
// ときだけ計数します。計測専用の実行ファイルでのみ定義してください。
//
// **可変な静的変数を置くことについて。** 本プロジェクトが禁じているのは
// `include/krisite/arith/` 配下です（並列化の前提を守るため）。ここは `geom/` で、
// かつ既定ビルドには存在せず、`thread_local` なのでスレッド間で干渉しません。
// 代案として述語にカウンタを引き回す案がありましたが、`side` は `split_fragment` や
// `crosses` の内側から呼ばれるため、呼び出し経路すべてに引数を足すことになり、
// **測定のために本番の API を歪める**ことになるので採りませんでした。
#ifndef KRISITE_GEOM_COUNTERS_HPP
#define KRISITE_GEOM_COUNTERS_HPP

#include <cstdint>

namespace krisite::geom {

#if defined(KRISITE_COUNT_PREDICATES)

namespace counters {

inline thread_local std::uint64_t side_calls = 0;
inline thread_local std::uint64_t intersect3_calls = 0;

// ---- ★ `side` の被符号値の【実際の】ビット幅（`SPEC-phase5.md` §5.10.10 の案 E）----
//
// **測るのは「計算過程で現れる最大の中間結果の幅」ではなく、
// 【最終的な被符号値の幅】です。** 前者は上界（$9b+20$）で決まりますが、
// **後者は入力の分布で決まります。**
//
// **区分は 64 / 128 / 192 / それ以上の 4 つです。** リム数の段が 64 ビットごとなので、
// **192 ビット（3 リム）で済むなら 4 リムから 1 リム減ります**
// （Nehring-Wirxel の Table 1 で 192b → 256b が 103 → 142 サイクル、1.38 倍）。
inline thread_local std::uint64_t side_w64 = 0;
inline thread_local std::uint64_t side_w128 = 0;
inline thread_local std::uint64_t side_w192 = 0;
inline thread_local std::uint64_t side_wmore = 0;
/// **観測した最大幅**（上界と並べるために要ります。張り付いていれば案 E は不成立）。
inline thread_local std::uint64_t side_wmax = 0;

inline void reset() noexcept {
    side_calls = 0;
    intersect3_calls = 0;
    side_w64 = 0;
    side_w128 = 0;
    side_w192 = 0;
    side_wmore = 0;
    side_wmax = 0;
}

}  // namespace counters

#define KRISITE_COUNT(which) (++::krisite::geom::counters::which)

#else

#define KRISITE_COUNT(which) ((void)0)

#endif

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_COUNTERS_HPP

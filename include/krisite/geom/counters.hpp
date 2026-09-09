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

#include <atomic>
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

// ---- ★ E2 の判定が【選ぶ】リム数（`SPEC-phase5.md` §5.10.10）--------------------
//
// **§21.4 は「被符号値の実際の幅」を測りました。それは E1（計算した後で見る）の値です。**
// **E2 はオペランドの幅から積の上界をその場で決めるので、
// 【判定が選ぶリム数】は必ずそれ以上になります。**
//
// **両方を並べないと、判定がどれだけ保守的かが分かりません。**
inline thread_local std::uint64_t side_disp1 = 0;
inline thread_local std::uint64_t side_disp2 = 0;
inline thread_local std::uint64_t side_disp3 = 0;
inline thread_local std::uint64_t side_disp4 = 0;

// ---- ★ 見積もりのための「リム乗算の回数」------------------------------------
//
// **`arith::mul` は筆算なので、費用は `N × M`（リム乗算の回数）に比例します**
// （`ops.hpp` の 137〜151 行）。**時間を測る前に、演算回数で見積もりを出すために数えます。**
//
//   `side_mul_now`  いまの費用（**型の**リム数の積の和）
//   `side_mul_e2`   E2 の費用（**実際に使っている**リム数の積の和）
//
// **比が、乗算の仕事がどれだけ減るかの見積もりです。**
inline thread_local std::uint64_t side_mul_now = 0;
inline thread_local std::uint64_t side_mul_e2 = 0;
/// **判定そのものの費用**（`used_limbs` が読んだリムの数）。
///
/// **`CLAUDE.md`「前判定は、置き換える仕事より安くなければ意味がありません」。**
/// **削減する乗算の数と、この数を並べないと見積もりになりません。**
inline thread_local std::uint64_t side_disp_limbreads = 0;

// ---- ★ 述語の内訳（`SPEC-phase5.md` §5.10.11 の「種類別」）--------------------
//
// **`side` と `intersect3` だけでは、述語の合計が出ません。**
// **`cmp_h` は縫合の整列（全構成点を `lex_less` で並べる）で効きます。**
///
/// > **★ この 2 つだけ `std::atomic` です。**
/// > **縫合（`lex_less` による全構成点の整列）は【逐次部分】で、
/// > 並列区間の外にあります。** 段ごとの差分では拾えません。
/// > **計測専用のビルドにしか存在しないので、競合の費用は本番に出ません。**
inline std::atomic<std::uint64_t> cmp_h_calls{0};
inline std::atomic<std::uint64_t> side_ipoint_calls{0};

// ---- ★ 断片の生成の内訳（`SPEC-phase5.md` §5.10.12.4。`edge` の small-array の見積もり）--
//
// **`split_fragment` の中で作られる `std::vector` は 3 種あります**
// （符号の作業配列 `s`、`clip_edges` の `keep`、新しい `edge`）。
// **確保の回数を「断片の生成」に帰着させるために、それぞれ数えます。**
// **原子的にしてあるのは、駆動が `boolean()` の前後で差分を取るためです**（計測ビルドのみ）。
inline std::atomic<std::uint64_t> frag_split_calls{
    0};  ///< `split_fragment` の呼び出し（`s` の確保）
inline std::atomic<std::uint64_t> frag_split_early{0};  ///< 早期 return（`edge` の複製 1 回）
inline std::atomic<std::uint64_t> frag_split_both{0};   ///< 2 つに切った（`make` 2 回）
inline std::atomic<std::uint64_t> frag_make_calls{0};   ///< `make`（`keep` の確保）
inline std::atomic<std::uint64_t> frag_make_ok{0};      ///< 新しい `edge` を作った
inline std::atomic<std::uint64_t> frag_clip_calls{0};   ///< `clip_fragment` の呼び出し
/// **作った `edge` の要素数の分布**（0..8 と、9 以上）。**区分は small-array の枠に合わせます。**
inline std::atomic<std::uint64_t> frag_edge_hist[10] = {};

inline void reset() noexcept {
    side_calls = 0;
    intersect3_calls = 0;
    side_w64 = 0;
    side_w128 = 0;
    side_w192 = 0;
    side_wmore = 0;
    side_wmax = 0;
    side_disp1 = 0;
    side_disp2 = 0;
    side_disp3 = 0;
    side_disp4 = 0;
    side_mul_now = 0;
    side_mul_e2 = 0;
    side_disp_limbreads = 0;
    cmp_h_calls.store(0);
    side_ipoint_calls.store(0);
    frag_split_calls.store(0);
    frag_split_early.store(0);
    frag_split_both.store(0);
    frag_make_calls.store(0);
    frag_make_ok.store(0);
    frag_clip_calls.store(0);
    for (auto& h : frag_edge_hist) h.store(0);
}

}  // namespace counters

#define KRISITE_COUNT(which) (++::krisite::geom::counters::which)
/// **原子的な計数**（`cmp_h_calls` / `side_ipoint_calls`）。
#define KRISITE_COUNT_ATOMIC(which) \
    (::krisite::geom::counters::which.fetch_add(1, ::std::memory_order_relaxed))

#else

#define KRISITE_COUNT(which) ((void)0)
#define KRISITE_COUNT_ATOMIC(which) ((void)0)

#endif

}  // namespace krisite::geom

#endif  // KRISITE_GEOM_COUNTERS_HPP

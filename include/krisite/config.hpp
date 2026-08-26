// Krisite — 共通コンパイル時設定
//
// SPEC-phase0.md §2（座標モデル）, §3.4（実装上の要求）
#ifndef KRISITE_CONFIG_HPP
#define KRISITE_CONFIG_HPP

#include <cstddef>

// 入力座標のビット幅 b。CMake オプション KRISITE_COORD_BITS から与えられる。
#ifndef KRISITE_COORD_BITS
#define KRISITE_COORD_BITS 21
#endif

// 全演算のオーバーフロー検査。CMake オプション KRISITE_CHECKED_ARITH、
// または Debug ビルドで自動的に 1 になる。
#ifndef KRISITE_CHECKED_ARITH
#define KRISITE_CHECKED_ARITH 0
#endif

namespace krisite {

/// 入力座標のビット幅 b（符号ビット込み）。
///
/// SPEC-phase0.md §3.1 の表は「座標差 p2-p1 が b+1 ビット」から始まる。
/// これが成り立つのは座標そのものが符号付き b ビット、すなわち
/// `-2^(b-1) <= coord <= 2^(b-1) - 1` のときである。
/// b = 21 のとき区間は [-2^20, 2^20-1] で、符号なし 21 ビット（64bit Morton の
/// 1 軸ぶん）を原点中心にずらしたものに一致する。SPEC §2 の「|coord| < 2^b」は
/// 概略の記述であり、ビット幅解析の基準は §3.1 の表とする。
inline constexpr std::size_t kCoordBits = static_cast<std::size_t>(KRISITE_COORD_BITS);

static_assert(kCoordBits >= 4, "b が小さすぎます");
static_assert(kCoordBits <= 31, "IPoint は int32 座標。b <= 31 であること");

/// 座標の許容範囲（両端含む）。
inline constexpr long long kCoordMin = -(1LL << (kCoordBits - 1));
inline constexpr long long kCoordMax = (1LL << (kCoordBits - 1)) - 1;

}  // namespace krisite

#endif  // KRISITE_CONFIG_HPP

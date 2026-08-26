# Krisite — コーディング規約

CLAUDE.md の「コーディング規約」で「一度決めたら記録し、以降は機械的に従う」とされた
選択をここに固定します。**議論はここで終わり、以降は `.clang-format` に従ってください。**

## 決定事項

| 対象 | 決定 | 例 |
|---|---|---|
| 型名 | `snake_case`（標準ライブラリ風） | `fixed_int`, `mpz_t` に合わせる |
| 幾何の型名 | `PascalCase` | `IPoint`, `HPoint`, `Plane`, `Axis` |
| 関数・変数 | `snake_case` | `plane_from_triangle`, `min_bits` |
| 定数 | `kCamelCase` | `kCoordBits`, `kHomoXyz`, `kSide` |
| マクロ | `UPPER_SNAKE` + `KRISITE_` 接頭辞 | `KRISITE_CHECKED_ARITH`, `KRISITE_CHECK` |
| 名前空間 | `snake_case` | `krisite::arith`, `krisite::geom` |
| テンプレート引数 | `PascalCase` または 1 文字大文字 | `N`, `M`, `L`, `NB`, `DB` |

型名が 2 系統あるのは意図的です。`arith` 層は標準ライブラリの数値型と並べて書かれるので
`snake_case`、`geom` 層は独自の幾何概念なので `PascalCase` とします。
**層をまたいで混ぜないこと。**

## 書式

`.clang-format` に定義済み。要点だけ:

- インデント 4、タブ禁止
- 1 行 100 桁
- 中括弧は同じ行（Attach）
- `#include` は「自ヘッダ → C++ 標準 → 外部 → プロジェクト」の順、各群内でアルファベット順

## コメント

- 論文の式を実装した箇所には**必ず出典**（`EMBER §4.2 式(7)` のように節と式番号まで）
- ビット幅を伴う式には `SPEC-phase0.md §3.1` への参照
- 数式や論文からの引用は原文の記法を保つ。それ以外の説明文は日本語

```cpp
// SPEC-phase0.md §7.2: 構成点の平面に対する側
//   sign(w) * sign(N・V - d*w)
// ビット幅: 9b+20 → widths.hpp bits::kSide
```

## 算術コードの禁止事項（`include/krisite/arith/` 配下）

CLAUDE.md より。**例外なし**。

- 動的メモリ確保（`new`, `malloc`, `std::vector`, `std::string`）
- 例外送出（すべて `noexcept`）
- グローバル変数・可変な静的変数
- 仮想関数
- リム数の実行時変更（テンプレートパラメータで固定）
- **リム数の数値リテラル**。`geom/widths.hpp` の `constexpr` から導出すること

## テストコードは対象外

`tests/` と `bench/` には上記の禁止事項は適用しません（`std::vector` も GMP も使ってよい）。
ただし GMP はライセンス上、`KRISITE_BUILD_TESTS_WITH_GMP=ON` のときだけリンクされる
`*_gmp.cpp` からのみ参照してください。

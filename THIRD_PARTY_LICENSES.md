# Third-Party Licenses / 第三者コンポーネントのライセンス

**配布物は外部依存を一切持ちません。**

Krisite のライブラリ本体は `include/krisite/` 配下のヘッダのみで構成され、
C++20 標準ライブラリ以外に依存しません。したがって **Krisite を利用する側が
第三者ライセンスの義務を負うことはありません。**

以下に挙げるのは **テストとベンチマークのためだけに使う**コンポーネントです。
いずれも既定では無効で、明示的に CMake オプションを立てたときにのみ取得・
リンクされます。**配布物には含まれません。**

---

## 一覧

| コンポーネント | ライセンス | 用途 | 既定 | 配布物への混入 |
|---|---|---|---|---|
| [GMP](https://gmplib.org/) | LGPL-3.0 / GPL-2.0 (dual) | 差分テストの正解器 | **OFF** | なし |
| [Manifold](https://github.com/elalish/manifold) | Apache-2.0 | 位相の外部正解器 | **OFF** | なし |

Shewchuk の robust predicates（パブリックドメイン）は**参照した文献**であって、
コードの取り込みはありません（§「再実装について」）。

---

## GMP

- **用途**: `tests/` の差分テストで、固定幅整数演算と幾何述語の正解値を出すため。
  出力メッシュの体積（有理数）の厳密計算にも使う
- **有効化**: `-DKRISITE_BUILD_TESTS_WITH_GMP=ON`（既定 **OFF**）
- **リンク先**: テスト実行ファイルのみ。`krisite` ターゲットには一切リンクしない
- **取得**: システムにインストールされたものを `find_path` / `find_library` で探す。
  ソースの取り込み（vendoring）はしない

**LGPL のコードを本体に混ぜないことを、文書ではなく機構で保証しています。**
CI の GMP ジョブは、テストを外した構成（`KRISITE_BUILD_TESTS=OFF`）で成果物を作り、
`ldd` に `gmp` が現れたら失敗します。

## Manifold

- **用途**: `tests/csg/test_manifold.cpp` で、ブール演算の出力の連結成分数と種数を
  独立な実装と突き合わせるため
- **有効化**: `-DKRISITE_BUILD_TESTS_WITH_MANIFOLD=ON`（既定 **OFF**）
- **リンク先**: `manifold_oracle` テスト実行ファイルのみ
- **取得**: CMake `FetchContent` で**特定コミットに固定**して取得
  （`KRISITE_MANIFOLD_COMMIT`）。vendoring はしない
- **推移的依存**: `MANIFOLD_CROSS_SECTION` / `MANIFOLD_PAR` / `MANIFOLD_DOWNLOADS`
  などをすべて OFF にして、追加の依存を引かない

**ライセンスの確認も機構に落としています。** CMake の構成時に取得先の `LICENSE` を
実際に読み、Apache License でなければ `FATAL_ERROR` で停止します。README の記述は
信用しません。CI は取得コミットの SHA と `LICENSE` の冒頭をログに残します。

---

## 再実装について

Krisite のアルゴリズムは**論文からの再実装**です。著作権は表現を保護し、
アルゴリズムそのものは保護しません。

**GPL / LGPL のコードは参照・引用・移植していません。** 具体的には、
CGAL、`MarcoAttene/*`、`gcherchi/*`、`mangoleaves/OpenMeshCraft`、MeshLab / VCGlib、
LASlib のいずれからもコードを取り込んでいません。

参考にした文献は各 `docs/SPEC-phase<N>.md` の参考文献表にあります。
論文 PDF はリポジトリに含めていません。

---

## Krisite 自身のライセンス

MIT License。`LICENSE` を参照してください。

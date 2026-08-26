<p align="center">
  <img src="assets/krisite-logo.svg#gh-light-mode-only" alt="Krisite" height="72">
  <img src="assets/krisite-logo-dark.svg#gh-dark-mode-only" alt="Krisite" height="72">
</p>

<p align="center"><em>Krisite — exact, plane-based geometry for point clouds and meshes</em></p>

---

`Krisite` は 3D データを扱う C++20 ヘッダオンリーライブラリです。
最終目標は点群圧縮・メッシュ化・厳密ブール演算の統合で、ブール演算の中核として
**EMBER**（Trettner, Nehring-Wirxel, Kobbelt, SIGGRAPH 2022）の再現実装を目指します。

名前の由来は **kris**（短剣、かつ crystal の語幹 κρύσταλλος）+ **-ite**（鉱物の接尾辞）。
平面ベース表現を採る本ライブラリにとって、結晶が平面で囲まれた立体であることは
単なる比喩ではありません。

## 設計の核

頂点を明示的な座標ではなく **3 枚の平面の交点** として表現し、すべての幾何述語を
**固定幅の多倍長整数**で厳密に評価します。浮動小数点は境界（ファイル入出力）でのみ扱います。

- 動的メモリ確保が発生しない（並列化と相性が良い）
- 浮動小数点フィルタとそのフォールバック経路が不要
- 退化構成（共平面・共線・重複頂点）でも分岐が増えない

平面も点も同次 4 元ベクトルとして持ちます。平面は `[a, b, c, d]` で `N·x + d = 0`
（`d = -N·p₁`）、構成点は `[x, y, z, w]`（実座標は `V/w`）。この流儀により、
最重要述語 `side` が **単一の 4 次元内積** `sign(w)·sign(a·x + b·y + c·z + d·w)` になります。

## 現在の状況

**Phase 0（固定幅厳密整数演算と平面ベース述語）を実装中です。**
仕様は [`docs/SPEC-phase0.md`](docs/SPEC-phase0.md)。

| Phase | 内容 | 状態 |
|---|---|---|
| **0** | 固定幅厳密整数 + 平面ベース述語 | **実装中** |
| 1 | 出力抽出の最小検証（立方体2個、固定深度分割、単スレッド） | 未着手 |
| 2 | 適応的再帰分割 + early-out 判定 | 未着手 |
| 3 | work-stealing 並列化、継ぎ目の整合性検証 | 未着手 |
| 4 | Thingi10K 全件検証 | 未着手 |
| 5+ | 点群コーデック、GWN、メッシュ化 | 未着手 |

Phase 1 が最大の難所です（Blender が EMBER 実装を断念したのがこの段階）。
Phase 0 は Phase 1 の可否判定を早く行うための土台です。

## ドキュメント

| ファイル | 内容 |
|---|---|
| [`docs/SPEC-phase0.md`](docs/SPEC-phase0.md) | Phase 0 の仕様。ビット幅解析、述語一覧、テスト要件 |
| [`docs/STYLE.md`](docs/STYLE.md) | コーディング規約（命名、書式、算術コードの制約） |
| [`CLAUDE.md`](CLAUDE.md) | 開発方針。**ライセンス方針を含む** |
| [`tools/fix-spec.py`](tools/fix-spec.py) | 仕様書の編集漏れを機械的に直すスクリプト |

## ライセンス

**MIT。** この制約は他のすべてに優先します。

ライブラリ本体は外部依存を一切持ちません。GMP（LGPL）は差分テストの正解器として
`KRISITE_BUILD_TESTS_WITH_GMP=ON` のときにテストバイナリからのみリンクされ、
配布物には含まれません。GPL/LGPL のコード（CGAL、Indirect_Predicates、
OpenMeshCraft、VCGlib 等）は参照・引用・移植していません。方針の詳細は
[`CLAUDE.md`](CLAUDE.md) を参照してください。

## 参考文献

| 略称 | 文献 |
|---|---|
| **EMBER** | Trettner, Nehring-Wirxel, Kobbelt. *EMBER: Exact Mesh Booleans via Efficient & Robust Local Arrangements.* ACM TOG 41(4), SIGGRAPH 2022. |
| **OEBSP** | Nehring-Wirxel, Trettner, Kobbelt. *Fast Exact Booleans for Iterated CSG using Octree-Embedded BSPs.* CAD 135, 2021. |
| **Levy24** | Bruno Lévy. *Exact predicates, exact constructions and combinatorics for mesh CSG.* arXiv:2405.12949. |
| **Shewchuk97** | Shewchuk. *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates.* DCG 18(3), 1997. |

実装はすべて論文からの再実装です。

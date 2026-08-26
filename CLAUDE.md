# CLAUDE.md — krisite プロジェクト

## 言語

**このプロジェクトでは日本語でやりとりします。**

- 応答・説明・質問・提案はすべて日本語
- コミットメッセージ、Issue、PR 説明も日本語
- **ただし以下は英語**: 識別子、型名、関数名、ファイル名、`docs/` 配下の見出し以外の技術用語
- コードコメントは日本語可。ただし数式や論文からの引用は原文の記法を保つこと

---

## プロジェクト概要

`Krisite` は、3D データを扱う C++20 ライブラリです。最終目標は点群圧縮・メッシュ化・
厳密ブール演算の統合ですが、**現在は Phase 0 のみ実装中**です。

| 項目 | 値 |
|---|---|
| 名前空間 | `krisite::`（短縮別名 `namespace kri = krisite;` を提供してよい） |
| ヘッダ | `include/krisite/...` |
| CMake 接頭辞 | `KRISITE_` |
| ファイル拡張子 | `.kris`（Phase 5 以降） |
| マジックバイト | ASCII `"KRIS"` + `uint32` バージョン |

名前の由来は **kris**（短剣、かつ crystal の語幹）+ **-ite**（鉱物の接尾辞）です。
削ることと結晶をひとつの語に込めています。README のタグラインは下記を既定とします。

```
Krisite — exact, plane-based geometry for point clouds and meshes
```

### 現在のフェーズ

**Phase 0: 固定幅厳密整数演算と平面ベース述語**

仕様は `docs/SPEC-phase0.md` にあります。**作業前に必ず全文を読んでください。**

Phase 0 の目的は、後続すべてが依存する算術基盤を作ることです。ここが不正確だと
上流が静かに壊れます。**性能より正しさを優先してください。**

### 全体ロードマップ（参考。Phase 0 では着手しない）

| Phase | 内容 |
|---|---|
| **0** | **固定幅厳密整数 + 平面ベース述語** ← 現在地 |
| 1 | 出力抽出の最小検証（立方体2個、固定深度分割、単スレッド） |
| 2 | 適応的再帰分割 + early-out 判定 |
| 3 | work-stealing 並列化、継ぎ目の整合性検証 |
| 4 | Thingi10K 全件検証 |
| 5+ | 点群コーデック、GWN、メッシュ化 |

Phase 1 が最大の難所です（Blender が EMBER 実装を断念したのがこの段階）。
Phase 0 は Phase 1 の可否判定を早く行うための土台です。

---

## ライセンス方針（最重要）

**このプロジェクトは MIT ライセンスで公開します。この制約は他のすべてに優先します。**

### 絶対に禁止

以下のコードを **参照・引用・移植・コピー** しないでください。ライセンス汚染は
後から取り返しがつきません。

| 対象 | ライセンス |
|---|---|
| CGAL | GPL |
| `MarcoAttene/Indirect_Predicates` | LGPL-3.0 |
| `MarcoAttene/TMesh_Kernel`, `CDT`, その他 Attene 氏のリポジトリ | LGPL / GPL |
| `gcherchi/*`（本体は MIT だが LGPL 依存を含む） | 実質 LGPL |
| `mangoleaves/OpenMeshCraft` | GPL-3.0 |
| MeshLab / VCGlib | GPL |
| LASlib（LAStools 同梱版） | LGPL |

**判断に迷ったら実装せず、必ず確認を求めてください。**

### 条件付きで許可

| 対象 | 条件 |
|---|---|
| GMP / MPFR | **テストのみ**。`KRISITE_BUILD_TESTS_WITH_GMP=ON`（既定 OFF）でのみ有効。ライブラリ本体に絶対にリンクしない |
| Geogram (BSD-3) | **コードは使わない**。Levy24 論文は読んでよい |
| Manifold (Apache-2.0) | Phase 1 以降の差分テスト用正解器としてのみ。本体から呼ばない |

### 制限なく許可

- **論文からの再実装**。著作権は表現を保護し、アルゴリズムは保護しません
- Shewchuk の predicates（パブリックドメイン）
- 自分で導出した数式・アルゴリズム

### 依存を追加するとき

**新しい外部依存を入れる前に必ず確認を求めてください。** 勝手に追加しないこと。
確認時には以下を報告してください。

1. ライセンス名（LICENSE ファイルを実際に読んで確認したもの。README の記述を信用しない）
2. その依存の推移的依存とそのライセンス
3. ヘッダオンリーか、静的リンクか、動的リンクか
4. 自前実装した場合の工数見積もり

---

## ビルドとコマンド

```bash
# 構成（開発時の既定）
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DKRISITE_CHECKED_ARITH=ON \
  -DKRISITE_BUILD_TESTS_WITH_GMP=ON

# ビルド
cmake --build build

# テスト
ctest --test-dir build --output-on-failure

# 単一テストを実行
ctest --test-dir build -R fixed_int --output-on-failure

# ベンチマーク
cmake -B build-rel -G Ninja -DCMAKE_BUILD_TYPE=Release -DKRISITE_CHECKED_ARITH=OFF
cmake --build build-rel
./build-rel/bench/pred_bench
```

### CMake オプション

| オプション | 既定 | 意味 |
|---|---|---|
| `KRISITE_COORD_BITS` | 21 | 入力座標のビット幅 `b` |
| `KRISITE_CHECKED_ARITH` | Debug で ON | 全演算のオーバーフロー検査 |
| `KRISITE_BUILD_TESTS` | ON | テストのビルド |
| `KRISITE_BUILD_TESTS_WITH_GMP` | **OFF** | GMP 差分テスト（LGPL、テスト専用） |
| `KRISITE_BUILD_BENCH` | OFF | ベンチマークのビルド |

---

## コーディング規約

### C++ 標準

C++20。コンパイラは GCC 13+ / Clang 16+ / MSVC 2022+。

### 算術コードの制約（`include/krisite/arith/` 配下）

以下は **例外なく守ってください**。並列化と性能の前提が崩れます。

- **動的メモリ確保を行わない**（`new`, `malloc`, `std::vector`, `std::string` 禁止）
- **例外を投げない**（すべて `noexcept`）
- **グローバル変数・可変な静的変数を持たない**（スレッド安全性のため）
- **仮想関数を使わない**
- リム数はテンプレートパラメータで持ち、実行時に変えない

### 型でビット幅を表現する

**これが本プロジェクトの設計の要です。**

```cpp
// 良い: 幅の増加が型に現れる。オーバーフローがコンパイル時に防がれる
template <std::size_t N, std::size_t M>
fixed_int<N + M> mul(const fixed_int<N>&, const fixed_int<M>&) noexcept;

// 悪い: 実行時にしか壊れない
fixed_int<4> mul(const fixed_int<4>&, const fixed_int<4>&);
```

リム数を数値リテラルで書かないこと。`constexpr` 関数から導出してください。

### 命名

- 型: `PascalCase`（`FixedInt` ではなく標準ライブラリ風に `fixed_int` でもよい。統一されていればどちらでも）
- 関数・変数: `snake_case`
- 定数: `kCamelCase` または `UPPER_SNAKE`。統一すること
- 名前空間: `krisite::arith`, `krisite::geom`

一度決めたら `.clang-format` と `docs/STYLE.md` に記録し、以降は機械的に従ってください。

### コメント

論文の式を実装した箇所には、**必ず出典を書いてください。**

```cpp
// EMBER §4.2 式(7): 構成点の平面に対する側
// sign(w) * sign(a*x + b*y + c*z + d*w)
// ビット幅: 9b+20 → SPEC-phase0.md §3.1 参照
```

数ヶ月後に読み返したとき、論文のどこを見ればいいか分かることが重要です。

---

## 作業の進め方

### テストを先に書く

**述語を実装する前に、その述語の GMP 差分テストを書いてください。**
正解器なしに述語を書くと、間違いに気づけません。

### 退化ケースを軽視しない

ランダムテストは共平面・共線・重複頂点をほぼ生成しません。
`docs/SPEC-phase0.md` §8.2 の表にあるケースを **明示的に構成**してください。
このプロジェクトが対象とする実データ（CAD、量子化された点群）では、
退化構成は例外ではなく常態です。

### 早めに止まって確認を求める場面

以下に該当したら、進めずに報告してください。

- 仕様書に書かれていない設計判断が必要になったとき
- 新しい述語を追加する必要が生じたとき（ビット幅の再解析が必要）
- §3 の理論上界を実測が超えたとき（**これは重大です。即座に報告**）
- 外部依存を追加したくなったとき
- 仕様書の記述が実装上おかしいと気づいたとき（仕様のバグかもしれません）

**「たぶんこうだろう」で進めないでください。** 算術基盤の誤りは上流で発見が困難です。

### やらないこと

`docs/SPEC-phase0.md` §10 の非目標を勝手に実装しないでください。特に:

- SIMD 化・浮動小数点フィルタ → **ベンチで基準線を取ってから判断します**
- 八分木、メッシュ構造、ファイル入出力 → Phase 1 以降
- 「ついでに」のリファクタリング → 指示された範囲に留めてください

---

## Git 運用

- ブランチ: `phase0/<topic>`（例: `phase0/fixed-int-mul`）
- コミットは論理単位で小さく。「テスト追加」と「実装」は分ける
- コミットメッセージは日本語。1 行目は 50 字以内の要約
- **テストが通らない状態でコミットしない**
- `main` に直接コミットしない

```
fixed_int の乗算を実装

x86-64 では _mulx_u64、それ以外は __int128 にフォールバック。
GMP 差分テスト 10^7 件でパスを確認。

SPEC-phase0.md §5.2
```

---

## 参考文献

`docs/SPEC-phase0.md` §1 の表を参照してください。
論文 PDF は `docs/papers/` に置いてください（**Git には含めない**。`.gitignore` 済み）。

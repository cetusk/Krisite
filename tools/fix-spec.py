#!/usr/bin/env python3
"""docs/SPEC-phase0.md に残っている編集漏れ 4 点を修正する。

2026-08 の改訂（平面を N・x + d = 0 に変更、cmp_h のビット幅追加）に伴う取りこぼしを
機械的に直します。実装側（include/krisite/）はすでに修正後の内容に従っています。

修正内容:
  1. §3.1 冒頭の前提文が旧版（|coord| < 2^b）のまま → §2 の符号付き b ビットに合わせる
  2. §3.1 構成点の表「1 列を d ベクトルに置換」 → 正しくは -d
  3. §3.1 に置換列の符号についての注記を追加（2 の理由）
  4. §6 のコード片のコメントが旧流儀「平面 N·p = d」のまま
  5. §8.2 退化ケース表「格子の端 | 座標が ±(2^b − 1)」 → §2 の両端に合わせる

使い方:
  python tools/fix-spec.py            差分を表示して適用（.bak を残す）
  python tools/fix-spec.py --check    差分を表示するだけ（書き換えない）
  python tools/fix-spec.py --no-backup

安全策:
  - すべての置換箇所が「ちょうど 1 箇所」見つかることを **書き換える前に** 確認します。
    ひとつでも見つからなければ何も書き換えずに終了します（全か無か）。
  - 適用済みの箇所は自動で読み飛ばすので、二重に実行しても壊れません。
  - CRLF / BOM は元のファイルの形式をそのまま保ちます。
"""

import argparse
import difflib
import os
import shutil
import sys

SPEC_RELPATH = os.path.join("docs", "SPEC-phase0.md")

# 注記本文。前後に空行が入るよう、先頭と末尾に改行を持たせている。
CRAMER_NOTE = (
    "\n"
    r"> **置換列の符号に注意**: 平面を $N \cdot x + d = 0$ と置いたので、3 平面の交点が" "\n"
    r"> 満たす連立は $N \cdot X = -d$ です。したがって Cramer で置き換える列は $d$ ではなく" "\n"
    r"> $-d$ になります。$d$ をそのまま置換すると、符号が反転した点 $-V$ が返ります。" "\n"
    ">\n"
    r"> 等価な流儀として $w = -\det(N_1, N_2, N_3)$ と定義して置換列を $d$ のままにする手も" "\n"
    r"> ありますが、本書は $w = \det(N_1, N_2, N_3)$ を採ります。" "\n"
)

# `·` は U+00B7、`−` は U+2212。原文の文字をそのまま使うこと。
PLANE_STRUCT_OLD = (
    "// 平面 N·p = d\n"
    "template <std::size_t NB, std::size_t DB>\n"
    "struct Plane {\n"
    "    arith::fixed_int<NB> a, b, c;   // 法線\n"
    "    arith::fixed_int<DB> d;         // オフセット\n"
    "};"
)
PLANE_STRUCT_NEW = (
    "// 平面 N·p + d = 0"
    "（同次 4 元ベクトル [a,b,c,d]。d = -N·p1 → §3.1）\n"
    "template <std::size_t NB, std::size_t DB>\n"
    "struct Plane {\n"
    "    arith::fixed_int<NB> a, b, c;   // 法線\n"
    "    arith::fixed_int<DB> d;         // オフセット d = -N·p1\n"
    "};"
)


class Edit:
    """1 箇所の置換。`marker` が本文にあれば適用済みとみなす。"""

    def __init__(self, num, title, old, new, marker):
        self.num = num
        self.title = title
        self.old = old
        self.new = new
        self.marker = marker


EDITS = [
    Edit(
        1,
        "§3.1 の前提文を §2 の座標範囲に合わせる",
        old=(
            "入力座標の絶対値が $2^b$ 未満とする。"
            "以下、すべて符号付きのビット数"
            "（符号ビット込み）。"
        ),
        new=(
            r"入力座標は **符号付き $b$ ビット**とする（§2。$-2^{b-1} \le \mathrm{coord} < 2^{b-1}$）。"
            "\n"
            "以下、すべて符号付きのビット数（符号ビット込み）。"
        ),
        marker="入力座標は **符号付き $b$ ビット**とする",
    ),
    Edit(
        2,
        "§3.1 構成点の表: Cramer の置換列を -d に直す",
        old="| $x, y, z$ | Cramer（1 列を $d$ ベクトルに置換） | $7b + 14$ |",
        new="| $x, y, z$ | Cramer（1 列を $-d$ ベクトルに置換） | $7b + 14$ |",
        marker="Cramer（1 列を $-d$ ベクトルに置換）",
    ),
    Edit(
        3,
        "§3.1: 置換列の符号についての注記を追加",
        old="\n**主要述語** — 構成点 $V$ が平面 $(N, d)$ のどちら側か：",
        new=CRAMER_NOTE + "\n**主要述語** — 構成点 $V$ が平面 $(N, d)$ のどちら側か：",
        marker="**置換列の符号に注意**",
    ),
    Edit(
        4,
        "§6: Plane のコメントを新しい流儀に直す",
        old=PLANE_STRUCT_OLD,
        new=PLANE_STRUCT_NEW,
        marker="// 平面 N·p + d = 0",
    ),
    Edit(
        5,
        "§8.2: 「格子の端」を §2 の両端に直す",
        old="| 格子の端 | 座標が ±(2^b − 1) |",
        new="| 格子の端 | 座標が -2^(b-1) と 2^(b-1) − 1（§2 の両端） |",
        marker="座標が -2^(b-1) と 2^(b-1)",
    ),
]


def find_spec(explicit):
    """SPEC のパスを決める。スクリプトの位置からリポジトリ根を推定する。"""
    if explicit:
        return explicit
    here = os.path.dirname(os.path.abspath(__file__))
    for base in (os.path.dirname(here), os.getcwd()):
        cand = os.path.join(base, SPEC_RELPATH)
        if os.path.isfile(cand):
            return cand
    return os.path.join(os.path.dirname(here), SPEC_RELPATH)


def read_spec(path):
    """(本文, BOMの有無, 改行がCRLFか) を返す。本文の改行は LF に正規化する。"""
    with open(path, "rb") as f:
        raw = f.read()
    has_bom = raw.startswith(b"\xef\xbb\xbf")
    if has_bom:
        raw = raw[3:]
    text = raw.decode("utf-8")
    is_crlf = "\r\n" in text
    if is_crlf:
        text = text.replace("\r\n", "\n")
    return text, has_bom, is_crlf


def write_spec(path, text, has_bom, is_crlf):
    if is_crlf:
        text = text.replace("\n", "\r\n")
    raw = text.encode("utf-8")
    if has_bom:
        raw = b"\xef\xbb\xbf" + raw
    with open(path, "wb") as f:
        f.write(raw)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    ap = argparse.ArgumentParser(
        description="docs/SPEC-phase0.md の編集漏れを修正する",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--check", action="store_true", help="差分を表示するだけで書き換えない")
    ap.add_argument("--no-backup", action="store_true", help=".bak を作らない")
    ap.add_argument("--spec", default=None, help="SPEC-phase0.md のパス（既定は自動検出）")
    args = ap.parse_args()

    path = find_spec(args.spec)
    if not os.path.isfile(path):
        print("[ERROR] SPEC が見つかりません: %s" % path)
        print("        --spec でパスを指定してください。")
        return 1

    print("対象: %s" % path)
    original, has_bom, is_crlf = read_spec(path)
    print("形式: UTF-8%s / 改行 %s" % (" (BOM)" if has_bom else "", "CRLF" if is_crlf else "LF"))
    print("")

    # --- 書き換える前に全箇所を検査する（全か無か）---
    #
    # marker の判定を先に行うこと。注記の追加のような「挿入」の編集では、適用後も
    # old（挿入位置のアンカー）が残るため、old の有無だけで判断すると二重に挿入される。
    pending, already, problems = [], [], []
    for e in EDITS:
        if e.marker in original:
            already.append(e)
            continue
        n = original.count(e.old)
        if n == 1:
            pending.append(e)
        elif n == 0:
            problems.append("  [%d] 該当箇所が見つかりません: %s" % (e.num, e.title))
        else:
            problems.append("  [%d] 該当箇所が %d 箇所あります（1 箇所のはず）: %s" % (e.num, n, e.title))

    for e in already:
        print("  [%d] 適用済み  %s" % (e.num, e.title))
    for e in pending:
        print("  [%d] 適用する  %s" % (e.num, e.title))

    if problems:
        print("")
        print("[ERROR] 想定と違う箇所があるため、**何も書き換えずに** 終了します。")
        for p in problems:
            print(p)
        print("")
        print("        SPEC が本スクリプトの想定より新しい可能性があります。")
        print("        手で直すか、スクリプトの EDITS を更新してください。")
        return 1

    if not pending:
        print("")
        print("すべて適用済みです。変更はありません。")
        return 0

    updated = original
    for e in pending:
        updated = updated.replace(e.old, e.new, 1)

    print("")
    print("--- 差分 " + "-" * 62)
    diff = difflib.unified_diff(
        original.splitlines(keepends=True),
        updated.splitlines(keepends=True),
        fromfile=SPEC_RELPATH + " (before)",
        tofile=SPEC_RELPATH + " (after)",
        n=2,
    )
    for line in diff:
        sys.stdout.write(line if line.endswith("\n") else line + "\n")
    print("-" * 70)
    print("")

    if args.check:
        print("--check のため書き換えていません。")
        return 0

    if not args.no_backup:
        backup = path + ".bak"
        shutil.copyfile(path, backup)
        print("バックアップ: %s" % backup)

    write_spec(path, updated, has_bom, is_crlf)
    print("%d 箇所を修正しました。" % len(pending))
    print("")
    print("確認してください:")
    print("  - §3.1 の前提文と構成点の表（置換列が -d になっているか）")
    print("  - §6 の Plane のコメント")
    print("  - §8.2 の「格子の端」")
    return 0


if __name__ == "__main__":
    sys.exit(main())

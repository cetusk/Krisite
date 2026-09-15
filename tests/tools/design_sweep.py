#!/usr/bin/env python3
"""docs/DESIGN-phase5-vertex-level.md §43 の機械的な突き合わせ。

10 巡のレビューで、確実な誤りの過半が「前の巡の修正が、同じ記号を使う
他の箇所に波及していない」形だった。記号ごとに全出現を並べ、旧い書き方が
残っていたら落ちる。

**この番人自身が空回りしないよう、変異試験を同梱している**（--selftest）。
各検査について「わざと旧い記述に戻したら落ちる」ことを、その場で確かめる。
検査を足したら、対応する変異も足すこと。足さなければ selftest が落ちる。

使い方:
    python3 tests/tools/design_sweep.py              # 突き合わせ
    python3 tests/tools/design_sweep.py --selftest   # 変異試験（番人の検出力）
"""
import io, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DOC = os.path.join(HERE, "..", "..", "docs", "DESIGN-phase5-vertex-level.md")


def section43(text):
    return text[text.index("## 43. "):]


# --- 検査 ---------------------------------------------------------------
# forbid: §43 の中に【あってはならない】正規表現（旧い記述）
# mutate: その検査を発火させる変異（置換前, 置換後）。selftest が使う
CHECKS = [
    dict(name="A5 の対象は S",
         forbid=r"A1〜A5|A5[^\n]*\*\*T のみ\*\*|\| \*\*A5\*\* \| \*\*T\*\*",
         mutate=("| A5 | **下の表のとおり**（**S のみ**） |",
                 "| A5 | **下の表のとおり**（**T のみ**） |")),
    dict(name="|S| は模型ごと",
         forbid=r"\\lvert S\\rvert|\\binom\{\\lvert S\\rvert\}",
         mutate=(r"\binom{\lvert S_A\rvert}{2}+\binom{\lvert S_B\rvert}{2}",
                 r"2\binom{\lvert S\rvert}{2}")),
    dict(name="k_max は k_0 相対",
         forbid=r"k_\{\\max\} = 40(?!\D*\+)",
         mutate=(r"上限 $k_{\max} = k_0 + 40$", r"上限 $k_{\max} = 40$")),
    dict(name="(1) は H′ で書く",
         forbid=r"\*\*\(1\)\*\* \| \*\*前提成立・H 成立",
         mutate=("| **(1)** | **前提成立・H′ 成立・一致。**",
                 "| **(1)** | **前提成立・H 成立・一致。**")),
    dict(name="決定表は P-静/P-動 で書く（循環を作らない）",
         forbid=r"\| \*\*2\*\* \| \*\*前提通過が破れた",
         mutate=("| **2** | **P-静が破れた**", "| **2** | **前提通過が破れた**")),
    dict(name="共線面で必ず A1 が破れるとは書かない",
         forbid=r"共線面があると、S は必ず A1 を破ります(?!\*\*$)",
         mutate=("> **★ 共線面を S から除くと、A1 が破れることがあります**",
                 "> **★ 共線面があると、S は必ず A1 を破ります**")),
    dict(name="R_A+R_B=R を「P3 そのもの」と書かない",
         forbid=r"これは P3 そのものです(?!」と書いたのは誤り)",  # 取り下げの引用は除く
         mutate=("**上の式と【代数的に同値】です**", "**これは P3 そのものです**")),
    dict(name="経緯の節への依存を作らない",
         forbid=r"（§(?:2[2-9]|3[0-4]|3[5-9]|4[0-2])(?:\.\d+)? を見|詳しくは §(?:3[5-9]|4[0-2])",
         mutate=("**列は同ディレクトリ `README.md` の 34 列**",
                 "**列は詳しくは §38.1 を見てください**")),
]

# 式で持つ期待値（実測した数を書かない。文書から導ける関係だけを検査する）
def expectations(b):
    out = []
    budgets = [int(x) for x in re.findall(r"\| (\d+) 秒 \|", b)]
    out.append(("段の予算の合計が本文の合計と一致",
                sum(budgets) == 630 and "**合計 630 秒**" in b))
    out.append(("A5 の対の上限 = 2*C(1496,2)",
                "2{,}236{,}520" in b and 2 * (1496 * 1495 // 2) == 2236520))
    out.append(("代表点 = 2*(|S_A|+|S_B|)、共線面 0 枚で 5,984",
                "5,984" in b and 2 * (1496 + 1496) == 5984))
    out.append(("合成対照 = S 7 + I 7 + U 8 = 22",
                "**22 構成**" in b and 7 + 7 + 8 == 22))
    out.append(("自己検査の数: S/I/U の表の行数",
                len(re.findall(r"^\| S[1-7] \|", b, re.M)) == 7 and
                len(re.findall(r"^\| I[1-7] \|", b, re.M)) == 7 and
                len(re.findall(r"^\| U[1-8] \|", b, re.M)) == 8))
    out.append(("判定は (1)(2-a)(2-b)(3)(4a)(4b)(4c) の 7 種",
                all(("**(%s)**" % k) in b for k in
                    ["1", "2-a", "2-b", "3", "4a", "4b", "4c"])))
    return out


def tables_ok(b):
    L = b.split("\n"); bad = 0; n = 0; i = 0
    while i < len(L):
        if L[i].startswith("|"):
            blk = []; j = i
            while j < len(L) and L[j].startswith("|"):
                blk.append(L[j]); j += 1
            n += 1
            sep = len(blk) > 1 and re.match(r"^\|[\s:\-|]+\|$", blk[1].strip())
            cols = [len(re.sub(r"\\\|", "", r).split("|")) for r in blk]
            if not sep or len(set(cols)) != 1:
                print("  ★表が異常: %s" % blk[0][:70]); bad += 1
            i = j
        else:
            i += 1
    return n, bad


def run(text, quiet=False):
    b = section43(text); bad = 0
    n, tb = tables_ok(b); bad += tb
    if not quiet:
        print("表 %d 個 / 異常 %d" % (n, tb))
    for c in CHECKS:
        hits = [m for m in re.finditer(c["forbid"], b)]
        if hits:
            bad += len(hits)
            if not quiet:
                for h in hits:
                    print("  ★%s: %r" % (c["name"], b[h.start():h.start() + 60]))
        elif not quiet:
            print("  ok  %s" % c["name"])
    for name, ok in expectations(b):
        if not ok:
            bad += 1
            if not quiet:
                print("  ★期待値が合いません: %s" % name)
        elif not quiet:
            print("  ok  %s" % name)
    return bad


def selftest(text):
    """各検査について、変異を入れたら落ちることを確かめる。"""
    base = run(text, quiet=True)
    if base:
        print("★ 変異前に既に %d 件落ちています。先にそちらを直してください。" % base)
        return 1
    fail = 0
    for c in CHECKS:
        a, m = c["mutate"]
        if a not in text:
            print("★ %s: 変異の適用先 %r が見つかりません" % (c["name"], a[:40]))
            fail += 1; continue
        if run(text.replace(a, m, 1), quiet=True) == 0:
            print("★ %s: 変異を入れても落ちません（空回り）" % c["name"])
            fail += 1
        else:
            print("  ok  %s — 変異を検出" % c["name"])
    # 表の検査も変異させる
    k = text.index("## 43. ")   # §43 の中の区切り行を消す
    mut = text[:k] + text[k:].replace("|---|---|\n", "", 1)
    if run(mut, quiet=True) == 0:
        print("★ 表の検査: 区切り行を消しても落ちません"); fail += 1
    else:
        print("  ok  表の検査 — 区切り行の削除を検出")
    print("変異 %d 件 / 検出できなかったもの %d 件" % (len(CHECKS) + 1, fail))
    return fail


if __name__ == "__main__":
    t = io.open(DOC, encoding="utf-8").read()
    sys.exit(selftest(t) if "--selftest" in sys.argv else run(t))

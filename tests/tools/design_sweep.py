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
    """§43 の本文だけを返す。§44 が足されたらそこで切る（12 巡目の指摘）。"""
    b = text[text.index("## 43. "):]
    m = re.search(r"^## 4[4-9]\.", b[10:], re.M)
    return b[: 10 + m.start()] if m else b


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
         mutate=("| **(1)** | **前提成立・H′ 成立・E1 が一致。**",
                 "| **(1)** | **前提成立・H 成立・E1 が一致。**")),
    dict(name="決定表は P-静/P-動 で書く（循環を作らない）",
         forbid=r"\| \*\*2\*\* \| \*\*前提通過が破れた",
         mutate=("| **2** | **P-静が破れた**", "| **2** | **前提通過が破れた**")),
    dict(name="共線面で必ず A1 が破れるとは書かない",
         forbid=r"共線面があると、S は必ず A1 を破ります(?!\*\*$)",
         mutate=("> **★ 共線面を S から除くと、A1 が破れることがあります**",
                 "> **★ 共線面があると、S は必ず A1 を破ります**")),
    dict(name="R_A+R_B=R を「P3 そのもの」と書かない",
         forbid=r"これは P3 そのものです(?!」と(?:書いたのは誤り|いう名付け))",  # 取り下げの引用は除く
         mutate=("**$R_A+R_B=R$ と【必要十分】で同値**", "**これは P3 そのものです**")),
    dict(name="取り下げた記述が復活していない（P3′ の強さ・仮定）",
         forbid=r"P3′ は P3 より弱い|P3′ を仮定したうえで",
         mutate=("> **★ P3 と P3′ は【比較不能】です**（どちらも他方を導きません。§43.3）。",
                 "> **P3′ は P3 より弱い条件です**（§43.3）。")),
    dict(name="E3 の右辺が全箇所で -R_B",
         forbid=r"内側成分\}_j\) = -R(?!_B)",
         mutate=(r"$\sum_j S(\text{内側成分}_j) = -R_B$。**B のみ**",
                 r"$\sum_j S(\text{内側成分}_j) = -R$。**B のみ**")),
    dict(name="頂点数・面数は N_v / N_f（n_f は法線）",
         forbid=r"\$n_v\$|\$n_f\$ は本体|\\times 3n_[vf]",
         mutate=("先頭の $N_v$ / $N_f$ は本体です", "先頭の $n_v$ / $n_f$ は本体です")),
    dict(name="経緯の節への依存を作らない",
         forbid=r"（§(?:2[2-9]|3[0-4]|3[5-9]|4[0-2])(?:\.\d+)? を見|詳しくは §(?:3[5-9]|4[0-2])",
         mutate=("**列は同ディレクトリ `README.md` の 34 列**",
                 "**列は詳しくは §38.1 を見てください**")),
]

def num(b, pat, cast=int):
    """文書から数を読む。見つからなければ None（＝期待値が落ちる）。"""
    m = re.search(pat, b)
    return cast(m.group(1).replace(",", "").replace("{", "").replace("}", "")) if m else None


# 式で持つ期待値。**数はすべて文書から読み、スクリプトに直書きしない。**
# 直書きすると、文書だけを書き換えたときに検出できない（10 巡目の指摘）。
def expectations(b):
    out = []
    budgets = [int(x) for x in re.findall(r"\| (\d+) 秒 \|", b)]
    total = num(b, r"\*\*合計 ([\d,]+) 秒\*\*")
    D = num(b, r"\*\*\$D = ([\d,]+)\$ 秒\*\*")
    out.append(("段の予算の合計 = 本文の合計 かつ <= D",
                bool(budgets) and total is not None and sum(budgets) == total
                and D is not None and total <= D))

    F = num(b, r"面数はどちらも ([\d,]+)")                      # 1,496
    pairs = num(b, r"\$([\d,{}]+)\$ 対で、これが上限")                 # 2,236,520
    reps = num(b, r"共線面が 0 枚なら ([\d,]+) 個")              # 5,984
    out.append(("A5 の対の上限 = 2*C(F,2)",
                None not in (F, pairs) and pairs == 2 * (F * (F - 1) // 2)))
    out.append(("代表点 = 2*(F+F)", None not in (F, reps) and reps == 4 * F))

    nS = len(re.findall(r"^\| S\d+ \|", b, re.M))
    nI = len(re.findall(r"^\| I\d+ \|", b, re.M))
    nU = len(re.findall(r"^\| U\d+ \|", b, re.M))
    nc = num(b, r"\*\*(\d+) 構成\*\*")
    n14 = num(b, r"（(\d+) 構成）")
    out.append(("構成数 = S+I+U", nc is not None and nc == nS + nI + nU))
    out.append(("自己検査の構成数 = S+I", n14 is not None and n14 == nS + nI))
    # 行数は「表の行数」と「本文の構成数」の一致で見る（直書きしない）
    n7 = num(b, r"\*\*(\d+) 構成\*\*")
    out.append(("S/I/U の行数が 1 以上で、構成数と整合",
                min(nS, nI, nU) >= 1 and n7 is not None and nS + nI + nU == n7))

    bits_xyz = num(b, r"`kHomoXyz = (\d)b\+15`")
    bits_w = num(b, r"`kHomoW = (\d)b\+13`")
    bb = num(b, r"\$b = (\d+)\$")
    tot_bits = num(b, r"= (\d+)\$ ビット")
    if None not in (bits_xyz, bits_w, bb, tot_bits):
        x, w = bits_xyz * bb + 15, bits_w * bb + 13
        out.append(("ビット幅: 3*(7b+15)+(6b+13) = 文書の 625",
                    3 * x + w == tot_bits and ("**%d ビット**" % x) in b
                    and ("**%d ビット**" % w) in b))
    else:
        out.append(("ビット幅の式が読めません", False))

    by79 = num(b, r"真のビット詰め (\d+) バイト")
    by81 = num(b, r"バイト境界丸め (\d+) バイト")
    if None in (tot_bits, by79, by81):
        out.append(("バイト数の行が読めません（削除されたか、書式が変わった）", False))
    if None not in (tot_bits, by79, bits_xyz, bits_w, bb):
        x, w = bits_xyz * bb + 15, bits_w * bb + 13
        out.append(("79 = ceil(625/8)", by79 == -(-tot_bits // 8)))
        out.append(("81 = 3*ceil(162/8)+ceil(139/8)",
                    by81 == 3 * (-(-x // 8)) + (-(-w // 8))))

    n261 = num(b, r"\*\*(\d+)\*\* \| \*\*実測\*\*（\d+ 行を数えた）")
    nrows = num(b, r"。合計 ([\d,]+) 行）")
    npair = num(b, r"相異なる【対】\*\* \| \*\*(\d+)\*\*")
    nmark = num(b, r"の印がある行\*\* \| \*\*(\d+)\*\*")
    res = [os.path.join(HERE, "..", "..", "data", "thingi10k", f)
           for f in ("cp2b_results.txt", "cp3_results.txt")]
    if all(os.path.exists(f) for f in res) and None not in (n261, nrows, npair, nmark):
        hit, keys, mark = 0, set(), 0
        tot = 0
        for f in res:
            for ln in io.open(f, encoding="utf-8"):
                c = ln.split()
                if len(c) < 31:
                    continue
                tot += 1
                try:
                    ve, de = float(c[29]), float(c[30])
                except ValueError:
                    continue
                if ve > 1e-9 or de > 1e-9:
                    hit += 1; keys.add(c[0])
                    if "体積の篩" in ln:
                        mark += 1
        out.append(("§43.12 の 261 / 256 / 257 / 590 が実データと一致",
                    (hit, len(keys), mark, tot) == (n261, npair, nmark, nrows)))
    else:
        out.append(("§43.12 の数を実データと照合できません", False))
    n3op = num(b, r"実測 ([\d,]+)\*\* —")
    n4op = num(b, r"\| 候補 261 実行[^|]*\| ([\d,]+) \|")
    out.append(("4 演算の面数 = 3 演算 * 4/3",
                None not in (n3op, n4op) and n4op * 3 == n3op * 4))

    # --- 15 巡目: 値の検算。実データと導出で確かめる ---
    dnum = num(b, r"分子 ([\d,]+) 桁 / 分母")
    dden = num(b, r"分母 ([\d,]+) 桁")
    ev = os.path.join(HERE, "..", "..", "docs", "evidence", "gmp_diag_r1",
                      "cp3_gmp_results.txt")
    if os.path.exists(ev) and None not in (dnum, dden):
        row = [l.split() for l in io.open(ev, encoding="utf-8")
               if l.startswith("250394x45413")]
        if row and "/" in row[0][28]:
            n_, d_ = row[0][28].split("/")
            out.append(("列 29 の分子・分母の桁数が実データと一致",
                        (len(n_), len(d_)) == (dnum, dden)))
        else:
            out.append(("列 29 を読めません", False))
    else:
        out.append(("桁数の検算ができません", False))

    seedA = num(b, r"添字 4232 → 種 (\d+)")
    lst = os.path.join(HERE, "..", "..", "data", "thingi10k", "cp3.txt")
    if os.path.exists(lst) and seedA is not None:
        ids = [l.split()[0] for l in io.open(lst, encoding="utf-8") if l.split()]
        out.append(("種 = 1000 + 模型の添字 が実データと一致",
                    "250394" in ids and 1000 + ids.index("250394") == seedA))
    else:
        out.append(("種を実データと照合できません", False))

    a1 = num(b, r"\$(\d+)\+48=(\d+)\$")
    m7 = re.search(r"\$(\d+)\+(\d+)=(\d+)\$", b)
    out.append(("S7 の S_内 の和が合う",
                m7 is not None and int(m7.group(1)) + int(m7.group(2)) == int(m7.group(3))))

    # E3 の反例: 3 殻（外・中が外向き、内が内向き）から R_B と ΣS(内側) を再計算
    rb = num(b, r"\$R_B = (-?\d+)\$")
    si = num(b, r"\$\\sum_j S\(\\text\{内側成分\}_j\) = (\d+)\$")
    v = [6 * 10 ** 3, 6 * 6 ** 3, 6 * 2 ** 3]          # [0,10]^3 / [2,8]^3 / [4,6]^3
    sg = [1, 1, -1]
    reg = [(sum(sg[: k + 1]), v[k] - (v[k + 1] if k + 1 < 3 else 0)) for k in range(3)]
    RB = sum(x for w, x in reg if w > 0) - sum(w * x for w, x in reg)
    out.append(("E3 の反例の R_B と ΣS(内側) が再計算と一致",
                None not in (rb, si) and rb == RB and si == v[1]))

    Dv = num(b, r"\*\*\$D = ([\d,]+)\$ 秒\*\*")
    worst = num(b, r"最悪の総時間 \$D\+30 = ([\d,]+)\$ 秒")
    reap = num(b, r"回収上限 (\d+) 秒")
    out.append(("最悪の総時間 = D + 回収上限",
                None not in (Dv, worst, reap) and Dv + reap == worst))

    f5 = num(b, r"\| 5 対（\*\*診断で回した対\*\*） \| ([\d,]+) \|")
    mb60 = num(b, r"\| ([\d,]+) \| \*\*351 MB\*\*") if False else None
    if f5 is not None:
        tot = 0
        for f in ("cp2b_gmp_results.txt", "cp3_gmp_results.txt"):
            q = os.path.join(HERE, "..", "..", "docs", "evidence", "gmp_diag_r1", f)
            if os.path.exists(q):
                for l in io.open(q, encoding="utf-8"):
                    c = l.split()
                    if len(c) > 17:
                        tot += sum(int(c[k]) for k in range(14, 18))
        out.append(("5 対の面数が実データと一致", tot == f5))
        m351 = re.search(r"\| 5 対[^|]*\| [\d,]+ \| \*\*(\d+) MB\*\* \| \*\*(\d+) MB\*\*", b)
        out.append(("351 / 307 MB が 60F / 52.5F と一致",
                    m351 is not None
                    and round(60 * f5 / 1e6) == int(m351.group(1))
                    and round(52.5 * f5 / 1e6) == int(m351.group(2))))
    else:
        out.append(("5 対の面数を読めません", False))

    out.append(("判定は (1)(2-a)(2-b)(3)(4a)(4b)(4c) の 7 種が【意味の表】にある",
                all(re.search(r"^\| \*\*\(%s\)\*\* \|" % re.escape(k), b, re.M)
                    for k in ["1", "2-a", "2-b", "3", "4a", "4b", "4c"])))

    rows = [int(x) for x in re.findall(r"^\| \*\*(\d+)\*\* \| \*\*", b, re.M)]
    out.append(("決定表の順序が 1..N で連続し、最終行が (1)",
                rows == list(range(1, len(rows) + 1)) and len(rows) >= 8))

    # 方向の表の中だけを見る（§43 全体だと I1〜I7 の座標と衝突する。16 巡目の指摘）
    tbl = re.search(r"\| # \| 方向 \|.*?\n\n", b, re.S)
    dirs = re.findall(r"\$\((\d),(\d),(\d)\)\$", tbl.group(0)) if tbl else []
    ndir = num(b, r"レイの (\d+) 方向")
    out.append(("レイの方向が、本文の本数と一致し、相異なり、零でない",
                ndir is not None and len(dirs) == ndir and len(set(dirs)) == ndir
                and ("0", "0", "0") not in set(dirs)))

    # --- 保存値: 文書の数が、実データと一致するか（12 巡目の指摘） ---
    SA = num(b, r"\| \$S\(A\)\$ \| \*\*33\*\* \| `(-?\d+)`")
    SB = num(b, r"\| \$S\(B\)\$ \| \*\*34\*\* \| `(-?\d+)`")
    R = num(b, r"\| \$R\$ \| \*\*28\*\* \| `(-?\d+)`")
    ev = os.path.join(HERE, "..", "..", "docs", "evidence", "gmp_diag_r1",
                      "cp3_gmp_results.txt")
    if os.path.exists(ev) and None not in (SA, SB, R):
        row = [l.split() for l in io.open(ev, encoding="utf-8")
               if l.startswith("250394x45413")]
        if row:
            c = row[0]
            # 列番号も文書から読む（スクリプトに直書きしない。13 巡目の指摘）
            cSA = num(b, r"\| \$S\(A\)\$ \| \*\*(\d+)\*\*")
            cSB = num(b, r"\| \$S\(B\)\$ \| \*\*(\d+)\*\*")
            cR = num(b, r"\| \$R\$ \| \*\*(\d+)\*\*")
            ok = None not in (cSA, cSB, cR) and (
                int(c[cSA - 1]), int(c[cSB - 1]), int(c[cR - 1])) == (SA, SB, R)
            out.append(("保存値 S(A)/S(B)/R が cp3_gmp_results.txt と一致（列番号も文書から）",
                        ok))
        else:
            out.append(("保存値の行が見つかりません", False))
    else:
        out.append(("保存値を読めません（式が変わったか、証拠が無い）", False))

    # --- 予言値の内部整合: 殻間 + 2*内側 = S(B)、内側 = -R ---
    shell = num(b, r"\| 殻間 \| 1 \| `(\d+)`")
    inner = num(b, r"\| 内側 \| 2 \| `(\d+)`")
    out.append(("予言値: 殻間+2*内側 = S(B) かつ 内側 = -R",
                None not in (shell, inner, SB, R)
                and shell + 2 * inner == SB and inner == -R))

    # --- 式の形（判断を左右するので、文字列として固定する） ---
    # 式は正規表現ではなく【部分文字列】で固定する。
    # 正規表現はエスケープを 1 段誤ると黙って常に不一致／常に一致になり、
    # 番人が空回りする（12 巡目に実際に起きた）。
    forms = [
        ("R の一次定義", r"(\lvert A\cup B\rvert+\lvert A\cap B\rvert) - (S(A)+S(B))"),
        ("P3′ の右辺", r"\lvert A\cup B\rvert + \lvert A\cap B\rvert \;=\; \mu(\mathcal A)+\mu(\mathcal B)"),
        ("E1", r"**E1** | **$R_A + R_B = R$**"),
        ("E2 の符号と閾値", r"R_X = -\sum_{\text{領域}:\,w_X\ge2}(w_X-1)\,\mathrm{vol} \;+\; \sum_{\text{領域}:\,w_X\le-1}\lvert w_X\rvert\,\mathrm{vol}"),
        ("E3 の右辺は -R_B（定義の行）", r"**E3** | **$\sum_j S(\text{内側成分}_j) = -R_B$**"),
        ("E4", r"**E4** | **$\sum_j S_j = S(X)$**"),
        ("領域の体積の式", r"\mathrm{vol}(\text{領域}_j) = \lvert S_j \rvert - \sum"),
        ("k_0 の式", r"k_0 = \lceil \log_2 \max_i \lvert (n_f)_i \rvert \rceil"),
        ("A5 の対の数", r"\binom{\lvert S_A\rvert}{2}+\binom{\lvert S_B\rvert}{2}"),
        ("|S_A| の式", r"\lvert S_A \rvert = 1496 - c_A"),
    ]
    for nm, pat in forms:
        out.append(("式が変わっていない: " + nm, pat in b))

    # --- 決定表の判定値（最終行が (1)、各行の判定が定義済みの 7 種） ---
    rows2 = re.findall(r"^\| \*\*(\d+)\*\* \|.*\| \*\*\(([^)]*)\)\*\* \|$", b, re.M)
    # 判定値の集合は【意味の表】から読む。直書きすると、判定を 1 つ増やす
    # という正しい変更で偽陽性になる（15 巡目の指摘）
    known = set(re.findall(r"^\| \*\*\(([^)]+)\)\*\* \|", b, re.M))
    out.append(("決定表の最終行が (1)", bool(rows2) and rows2[-1][1] == "1"))
    out.append(("決定表の判定値がすべて定義済みの 7 種",
                bool(rows2) and all(v in known for _, v in rows2)))
    out.append(("(1) は最終行にだけ現れる",
                bool(rows2) and [v for _, v in rows2].count("1") == 1))
    return out


def first_break(P, F):
    """A0-1 → A0-2 → A0-3 → A1(T) → A2(T) → A3(T) → A4(T) → A1(S) → … の順に
    評価し、最初に破れた検査の名前を返す（§43.4 の評価順）。

    ★ 文言ではなく【幾何】を見る。17 巡目まで、合成対照の構成が
    「期待する停止箇所で本当に止まるか」に検出器が 1 つも無かった。
    """
    from collections import Counter

    def cross(a, c, d):
        u = [c[i] - a[i] for i in range(3)]
        v = [d[i] - a[i] for i in range(3)]
        return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0])

    if any(max(f) >= len(P) or min(f) < 0 for f in F):
        return "A0-1"
    if any(len(set(f)) < 3 for f in F):
        return "A0-2"
    col = [f for f in F if cross(P[f[0]], P[f[1]], P[f[2]]) == (0, 0, 0)]
    for tag, G in (("T", F), ("S", [f for f in F if f not in col])):
        ud, dr = Counter(), Counter()
        for f in G:
            for k in range(3):
                x, y = f[k], f[(k + 1) % 3]
                ud[frozenset((x, y))] += 1
                dr[(x, y)] += 1
        if set(ud.values()) != {2}:
            return "A1(%s)" % tag
        if set(dr.values()) != {1}:
            return "A2(%s)" % tag
        for v in {x for f in G for x in f}:
            es = [tuple(x for x in f if x != v) for f in G if v in f]
            adj = {}
            for e in es:
                adj.setdefault(e[0], []).append(e[1])
                adj.setdefault(e[1], []).append(e[0])
            if not adj or any(len(t) != 2 for t in adj.values()):
                return "A3(%s)" % tag
            seen = {next(iter(adj))}
            st = [next(iter(adj))]
            while st:                       # 単一の閉路か（連結か）
                x = st.pop()
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y); st.append(y)
            if len(seen) != len(adj):
                return "A3(%s)" % tag
        if len({tuple(sorted(f)) for f in G}) != len(G):
            return "A4(%s)" % tag
    return "A5(S) 以降"


def structure(b):
    """「中核の語・節を直して、参照元を洗っていない」形を機械的に見る。

    ★ 何を見て、何を見ないかを先に書く。
      見る  : 参照の実在（§43.N / §43.N(x) / 順序 N）、小節の連続、
              **参照の固定**（「この主張はこの節を指すはず」を表で宣言）、
              **記号の取り残し**（一般化した記号が古いまま残っていないか）
      見ない: 「説明の中身が参照先と合っているか」。
              13〜14 巡目に語の重なりで測ろうとしたが、**正しい言い換えで落ち、
              実際の欠陥は素通り**した。当てにならない番人は流されるので採らない。
              この範囲は独立レビューが担う。
    """
    out = []
    heads, subs = {}, {}
    cur = None
    for ln in b.split("\n"):
        m = re.match(r"^### (43\.\d+) (.+)$", ln)
        if m:
            cur = m.group(1); heads[cur] = m.group(2); subs[cur] = {}
        m = re.match(r"^#### \(([a-z])\) (.+)$", ln)
        if m and cur:
            subs[cur][m.group(1)] = m.group(2)
    rows = set(re.findall(r"^\| \*\*(\d+)\*\* \| \*\*", b, re.M))

    refs = set(re.findall(r"§(43\.\d+)", b))
    out.append(("§43.N への参照がすべて実在する",
                all(r in heads for r in refs),
                sorted(r for r in refs if r not in heads)))
    srefs = re.findall(r"§(43\.\d+)\(([a-z])\)", b)
    bad = [f"§{n}({x})" for n, x in srefs if x not in subs.get(n, {})]
    out.append(("§43.N(x) への参照がすべて実在する", not bad, bad))
    gaps = [n for n, d in subs.items() if d and
            sorted(d) != [chr(ord("a") + k) for k in range(len(d))]]
    out.append(("小節が (a) から連続している", not gaps, gaps))
    miss = sorted(x for x in set(re.findall(r"順序 (\d+)", b)) if x not in rows)
    out.append(("「順序 N」で送る先の行が実在する", not miss, miss))

    secs0 = {}
    cur = None
    for ln in b.split("\n"):
        m = re.match(r"^### (43\.\d+) ", ln)
        if m:
            cur = m.group(1); secs0[cur] = []
        if cur:
            secs0[cur].append(ln)
    secs0 = {k: "\n".join(v) for k, v in secs0.items()}

    # ★ 参照の固定: 「この語を含む行は、この節を指す」。改番の取り残しを捕まえる。
    #   行を特定する語と、要求する参照先を【宣言】する。
    PINS = [
        ("採用領域が", "§43.5(g)"),
        ("領域の列挙がそれに依る", "§43.5(f)"),
        ("$\\mathrm{vol}$ と領域の定義は", "§43.5(f)"),
        ("H′ の定義は", "§43.5(e)"),
        ("「内側成分」の定義は", "§43.5(e)"),
        ("その段が書くと決めたファイル", "§43.8"),
    ]
    bad = []
    for key, want in PINS:
        hits = [ln for ln in b.split("\n") if key in ln]
        if not hits:
            # ★ 鍵が消えたら落とす。消えたまま黙ると、言い換え 1 回で pin が死ぬ
            #   （15 巡目の指摘。MUST / LEFTOVER と同じ「閉じて壊れる」向きに揃える）
            bad.append("固定の鍵「%s」が §43 から消えました（言い換えたなら表も直す）" % key)
            continue
        for ln in hits:
            if want not in ln:
                bad.append("「%s」を含む行が %s を指していません: %s"
                           % (key, want, ln.strip()[:60]))
    out.append(("参照の固定（改番の取り残し）", not bad, bad[:4]))

    # ★ 必須の項目: 一度書くと決めたものが、消えていないか。
    #   参照元の主張（「§43.13 に項目として立ててあります」など）が空振りしないよう、
    #   実際に在るべき語を宣言する。
    MUST = [
        # --- 15 巡目: 今回書いた本文に網を掛ける。追記と同時に網を広げる ---
        ("43.10", "に【クランプ】する", "経路 P 手順 7 のクランプ"),
        ("43.10", r"$(p_i[k] - \mathrm{mid}[k]) \times", "経路 P 手順 6 の中心合わせ"),
        ("43.10", "3 軸の最大", "経路 P 手順 3 の ext"),
        ("43.10", "| 8 | **同一格子点の併合**", "経路 P 手順 8"),
        ("43.10", "KMSH", "原本 .kmesh の書式"),
        ("43.10", "（最近接・偶数優先）", "nearbyint の丸めモード"),
        ("43.10", "リトルエンディアン固定", "バイト順"),
        ("43.10", "Python では `round()`（偶数優先）を使います", "Python 側の丸め"),
        ("43.10", "$r[0..8]$ の 9 行、続いて $\\text{shift}[0..2]$ の 3 行", "hex の並び"),
        ("43.10", "`int32` $\\times 3N_v$", "quantized の頂点の型"),
        ("43.10", "最小の元の添字", "併合の代表の選び方"),
        ("43.9", "3 辺が 1 枚に減って", "I1 が破る辺の本数"),
        ("43.8", "$> 300 \\times 0.8$", "段 6 の起動判定の余裕"),
        ("43.8", "0x2545F4914F6CDD1D", "段 4 の乱数の定数"),
        ("43.8", "種は 1", "段 4 の種（0 は不動点）"),
        ("43.4", "A1(T)→A2(T)", "A1〜A4 を T→S の順に掛けること"),
        ("43.4", "A1 → A2 → A3 → A4 → A5", "検査の評価順"),
        ("43.13", "共線面があったときに判定を続ける機構", "共線面の未設計項目"),
        ("43.13", "経路 P で変換の生成までバイト一致", "経路 P の未確認項目"),
        ("43.13", "同一性 1・2 に残る libm 依存", "libm 依存の未確認項目"),
        ("43.7", "$S_j = 0$", "順序 7 の S_j=0 の受け皿"),
        ("43.7", "$R_A = 0$", "順序 7 の R_A=0 の自己検査"),
        ("43.8", "transform_A.hex", "経路 P へ渡す変換の受け渡しファイル"),
        ("43.3", "P3′", "P3′ の定義"),
        ("43.5", "包含", "包含の森の手順"),
    ]
    bad = [d for n, key, d in MUST if key not in secs0.get(n, "")]
    out.append(("必ず在るべき項目が消えていない", not bad, bad))

    # ★ 記号の取り残し: 一般化した記号が、古いまま残っていないか。
    #   (節, その節にあってはならない正規表現, 説明)
    LEFTOVER = [
        # --- 15 巡目: 誤って書いた値が戻っていないか ---
        ("43.3", r"分子が 7,782 桁(?!」と書いたのは誤り)", "列 29 の桁数の読み違い"),
        ("43.10", r"geom/widths\.hpp`。\$b\$ は", "kCoordMax の所在（config.hpp が正）"),
        ("43.8", r"\*\*種は 0\*\*|同じ生成器。種 0", "段 4 の種 0（不動点。出力が永久に 0）"),
        ("43.5", r"これが E2 の実体です(?!」と書いたのは言い過ぎ)", "検算を E2 と同一視"),
        ("43.5", r"\\mathrm\{measure\}\(w_B>0\) = ", "measure の手順が B 専用"),
        ("43.5", r"\\cdot\\mathrm\{vol\}\(\\text\{領域\}_j\) = S\(B\)", "検算が B 専用"),
        ("43.3", r"\*\*E2\*\* \| \*\*\$R_B =", "E2 が B 専用"),
        ("43.3", r"\*\*E4\*\* \| \*\*\$\\sum_j S_j = S\(B\)\$", "E4 が B 専用"),
    ]
    bad = [d for n, pat, d in LEFTOVER if re.search(pat, secs0.get(n, ""))]
    out.append(("記号の取り残し（一般化の波及漏れ）", not bad, bad))

    # ★ 性質で見る（文言を鍵にしない。16 巡目の指摘。
    #   文言を鍵にすると「項目を正しく直す」と「項目を消す」が区別できない）
    s9 = secs0.get("43.9", "")
    irows = re.findall(r"^\| (I\d+) \|[^|]*\| ([^|]*)\|", s9, re.M)
    bad = []
    for name, stop in irows:
        if not re.search(r"A(?:0-3|[1-5])\s*[（(]\s*(?:T|S)\s*(?:側)?\s*[）)]", stop):
            bad.append("%s の期待する停止箇所に面集合の併記がありません: %s"
                       % (name, stop.strip()[:40]))
    have = {n for n, _ in irows}
    for name in have:
        if not (re.search(r"\*\*%s\*\*: " % name, s9)
                or re.search(r"^\| \*\*%s\*\* \|" % name, s9, re.M)):
            bad.append("%s の構成が書かれていません" % name)
    out.append(("I 表の各行に、面集合を併記した停止箇所と、構成がある",
                bool(irows) and not bad, bad[:4]))

    # ★ コードの引用を、実ファイルと照合する（16 巡目の指摘。
    #   旧い誤記を禁じるだけでは、新しい別の誤った行番号が通る）
    # 範囲で引いているものは、範囲の【両端】も照合する
    CITE_RANGE = [
        ("tests/thingi10k/loader.hpp", 44, 65, "load_kmesh", "}"),
        ("tests/thingi10k/loader.hpp", 173, 178, "auto next", "};"),
        ("include/krisite/config.hpp", 36, 37, "kCoordMin", "kCoordMax"),
    ]
    CITE = [
        ("include/krisite/config.hpp", 36, "kCoordMin"),
        ("include/krisite/config.hpp", 37, "kCoordMax"),
        ("include/krisite/geom/widths.hpp", 70, "kHomoW"),
        ("include/krisite/geom/widths.hpp", 79, "kHomoXyz"),
        ("include/krisite/csg/polysoup.hpp", 116, "w[n.src] > 0"),
        ("tests/thingi10k/loader.hpp", 44, "load_kmesh"),
        ("tests/thingi10k/loader.hpp", 158, "a == b"),
        ("tests/thingi10k/loader.hpp", 173, "auto next"),
        ("tests/thingi10k/thingi_cp1.cpp", 1690, "prepare(raw, 1000 + i)"),
        ("CMakeLists.txt", 10, "KRISITE_COORD_BITS"),

    ]
    bad = []
    for path, line, token in CITE:
        base = os.path.basename(path)
        # §43 がその行番号で引いているか（basename でも full path でも）
        # 引用が【無い】ことは欠陥ではない（正当に消すことがある。17 巡目の指摘）。
        # 在るなら、その行が実コードと合っているかだけを見る。
        cited = False
        for m in re.finditer(r"`(?:[\w/.]*/)?%s:(\d+)(?:-(\d+))?" % re.escape(base), b):
            lo = int(m.group(1)); hi = int(m.group(2)) if m.group(2) else lo
            if lo <= line <= hi:
                cited = True; break
        if not cited:
            continue
        f = os.path.join(HERE, "..", "..", *path.split("/"))
        if not os.path.exists(f):
            bad.append("%s がありません" % path); continue
        L2 = io.open(f, encoding="utf-8").read().split("\n")
        if line > len(L2) or token not in L2[line - 1]:
            bad.append("%s:%d に %r がありません（実コードとずれています）"
                       % (base, line, token))
    for path, lo, hi, tlo, thi in CITE_RANGE:
        base = os.path.basename(path)
        found = [(int(x), int(y)) for x, y in
                 re.findall(r"`(?:[\w/.]*/)?%s:(\d+)-(\d+)`" % re.escape(base), b)]
        if (lo, hi) not in found:
            bad.append("%s:%d-%d の範囲の引用がありません（文書: %s）"
                       % (base, lo, hi, found))
            continue
        f = os.path.join(HERE, "..", "..", *path.split("/"))
        L2 = io.open(f, encoding="utf-8").read().split("\n") if os.path.exists(f) else []
        if hi > len(L2) or tlo not in L2[lo - 1] or thi not in L2[hi - 1]:
            bad.append("%s:%d-%d の範囲の両端が実コードとずれています" % (base, lo, hi))
    # 宣言していない範囲の引用があれば落とす（同じ引用が文書に複数あるとき、
    # 片方だけ書き換わったのを捕まえる。16 巡目の指摘）
    # 範囲の引用は、宣言してあれば両端を照合する。宣言に無い範囲は、
    # 【その範囲の両端が実コードにあるか】を直接見る（宣言表の保守を強いない）。
    declared = {(os.path.basename(pp), l, h) for pp, l, h, _, _ in CITE_RANGE}
    known = {os.path.basename(pp): pp for pp, *_ in CITE_RANGE}
    known.update({os.path.basename(pp): pp for pp, *_ in CITE})
    for base, l, h in re.findall(r"`(?:[\w/.]*/)?([\w.]+\.(?:hpp|cpp|py|txt)):(\d+)-(\d+)`", b):
        if (base, int(l), int(h)) in declared or base not in known:
            continue
        f = os.path.join(HERE, "..", "..", *known[base].split("/"))
        n = len(io.open(f, encoding="utf-8").read().split("\n")) if os.path.exists(f) else 0
        if int(h) > n or int(l) > int(h):
            bad.append("範囲の引用が実ファイルの外です: %s:%s-%s（全 %d 行）"
                       % (base, l, h, n))
    out.append(("コードの引用が実ファイルと一致する", not bad, bad[:4]))

    # ★ 文書の【値】を、実コードの値と照合する（16 巡目の指摘。
    #   行番号だけでなく、引用した定数そのものを突き合わせる）
    SRC = [
        ("fill", r"\$\\mathrm\{fill\} = ([\d.]+)\$",
         "tests/thingi10k/loader.hpp", r"double fill = ([\d.]+)"),
        ("kCoordMin の式", r"\\mathtt\{kCoordMin\} = -2\^\{(b-1|b)\}",
         "include/krisite/config.hpp", r"kCoordMin = -\(1LL << \(kCoordBits - (1)\)\)"),
        ("xorshift の右シフト 1", r"シフト (\d+) / \d+ / \d+",
         "tests/thingi10k/loader.hpp", r"seed \^= seed >> (\d+);"),
        ("xorshift の左シフト", r"シフト \d+ / (\d+) / \d+",
         "tests/thingi10k/loader.hpp", r"seed \^= seed << (\d+);"),
        ("xorshift の右シフト 2", r"シフト \d+ / \d+ / (\d+)",
         "tests/thingi10k/loader.hpp", r"seed >> 12;\s*\n.*\n\s*seed \^= seed >> (\d+);"),
        ("b の既定", r"\$b = (\d+)\$", "CMakeLists.txt", r"KRISITE_COORD_BITS (\d+)"),
    ]
    bad = []
    for name, dpat, path, spat in SRC:
        dm = re.search(dpat, b)
        f = os.path.join(HERE, "..", "..", *path.split("/"))
        sm = re.search(spat, io.open(f, encoding="utf-8").read()) if os.path.exists(f) else None
        if dm is None or sm is None:
            bad.append("%s: 文書側=%s コード側=%s（読めません）"
                       % (name, dm and dm.group(1), sm and sm.group(1)))
        else:
            dv = dm.group(1).replace("b-1", "1").replace("b", "0")
            if dv != sm.group(1):
                bad.append("%s: 文書 %r ≠ コード %r" % (name, dm.group(1), sm.group(1)))
    out.append(("文書の値が実コードと一致する", not bad, bad[:4]))

    # ★ 合成対照の構成が、期待する停止箇所で本当に止まるか（17 巡目の指摘）
    #   文書の表から頂点と面を読み、評価順のとおりに計算して突き合わせる。
    bad = []
    rows = re.findall(r"^\| \*\*(I\d+)\*\* \|[^|]*\| `([^`]*)` \| `([^`]*)` \|", b, re.M)
    if not rows:
        bad.append("I の構成の表（頂点と面）が読めません")
    for name, vs, fs in rows:
        try:
            P = [tuple(int(x) for x in t.strip("()").split(","))
                 for t in re.findall(r"\([-\d, ]+\)", vs)]
            F = [tuple(int(x) for x in t.strip("()").split(","))
                 for t in re.findall(r"\([-\d, ]+\)", fs)]
        except ValueError:
            bad.append("%s の頂点・面を読めません" % name); continue
        if name == "I1":                     # 基準の立方体と頂点を共有しないので
            got = first_break(P, F)          # 追加ぶんだけで判定が決まる
        else:
            got = first_break(P, F)
        # 期待する停止箇所は I 表（(b)）の 3 列目。行の中の最後の A?(T|S) を取る
        if name == "I1":                 # I1 は「共線の面を 1 枚含む」入力
            from collections import Counter as _C
            def _cr(a, c, d):
                u = [c[i] - a[i] for i in range(3)]; v = [d[i] - a[i] for i in range(3)]
                return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
            nc = sum(1 for f in F if _cr(P[f[0]], P[f[1]], P[f[2]]) == (0, 0, 0))
            if nc != 1:
                bad.append("I1 の共線面が %d 枚です（1 枚であるべき）" % nc)
        m2 = re.search(r"^\| %s \|[^|]*\|([^|]*)\|" % name, b, re.M)
        cand = re.findall(r"(A[\w-]+)\s*[（(]\s*(T|S)\s*(?:側)?\s*[）)]",
                          m2.group(1)) if m2 else []
        want = "%s(%s)" % cand[-1] if cand else None   # 表記を正規化して比べる
        if want is None:
            bad.append("%s の期待する停止箇所を読めません" % name)
        elif got != want:
            bad.append("%s は %s で止まります（期待は %s）" % (name, got, want))
    out.append(("I の構成が、期待する停止箇所で止まる（幾何を計算）", not bad, bad[:4]))

    # U 表: 方向を引数に取る関数（winding_by_ray / 退化判定）の回数は 16
    su = secs0.get("43.9", "")
    bad = []
    prev = ""
    for name, rest in re.findall(r"^\| (U\d+) \|(.*)$", su, re.M):
        cells = [c.strip() for c in rest.split("|")]
        if not cells:
            continue
        fn = prev if cells[0].startswith("同上") else cells[0]   # 「同上」を解決する
        prev = fn
        if len(cells) >= 4 and "winding_by_ray" in fn:
            if cells[-2] not in ("16", "**16**"):
                bad.append("%s の呼び出し回数が %r（方向を引数に取るなら 16）"
                           % (name, cells[-2]))
    out.append(("U 表の呼び出し回数が、方向の本数と整合する", not bad, bad[:4]))

    # レイ: 本文の本数・表の本数・「すべて評価する」の本数が一致する
    n1 = re.search(r"レイの (\d+) 方向", b)
    n2 = re.search(r"(\d+) 方向を【すべて】評価", b)
    tbl = re.search(r"\| # \| 方向 \|.*?\n\n", b, re.S)
    ndirs = len(re.findall(r"\$\((\d),(\d),(\d)\)\$", tbl.group(0))) if tbl else 0
    out.append(("レイの本数が、本文・表・「すべて評価」で一致する",
                None not in (n1, n2) and int(n1.group(1)) == int(n2.group(1)) == ndirs,
                [n1 and n1.group(1), n2 and n2.group(1), ndirs]))

    # 段 4 の標本: 各 N 対 × 2 = 合計
    m = re.search(r"各 ([\d,]+) 対、合わせて ([\d,]+) 対", b)
    out.append(("段 4 の標本数が内部で整合する",
                m is not None
                and int(m.group(1).replace(",", "")) * 2 == int(m.group(2).replace(",", "")),
                m.groups() if m else None))
    return out


def tables_ok(b, quiet=False):
    # 引用ブロック（"> "）の中の表も対象にする（11 巡目）。
    # ただし、引用の内外がまたがる表は Markdown として壊れるので、
    # 剥がす【前】に引用の深さが表の中で揃っていることも見る（14 巡目）。
    raw = b.split("\n")
    L = [re.sub(r"^> ?", "", x) for x in raw]
    depth = [1 if x.startswith(">") else 0 for x in raw]
    bad = 0; n = 0; i = 0
    while i < len(L):
        if L[i].startswith("|"):
            blk = []; j = i
            while j < len(L) and L[j].startswith("|"):
                blk.append(L[j]); j += 1
            n += 1
            if len(set(depth[i:j])) != 1:
                if not quiet:
                    print("  ★表が引用の内外にまたがる: %s" % blk[0][:60])
                bad += 1
            sep = len(blk) > 1 and re.match(r"^\|[\s:\-|]+\|$", blk[1].strip())
            cols = [len(re.sub(r"\\\|", "", r).split("|")) for r in blk]
            if not sep or len(set(cols)) != 1:
                if not quiet:
                    print("  ★表が異常: %s" % blk[0][:70])
                bad += 1
            i = j
        else:
            i += 1
    return n, bad


def run(text, quiet=False, collect=None):
    """collect に list を渡すと、発火した検査の名前を積む（主検出器の記録）。"""
    b = section43(text); bad = 0
    n, tb = tables_ok(b, quiet); bad += tb
    if tb and collect is not None:
        collect.append("表")
    if not quiet:
        print("表 %d 個 / 異常 %d" % (n, tb))
    for c in CHECKS:
        hits = [m for m in re.finditer(c["forbid"], b)]
        if hits:
            bad += len(hits)
            if collect is not None:
                collect.append("検査:" + c["name"])
            if not quiet:
                for h in hits:
                    print("  ★%s: %r" % (c["name"], b[h.start():h.start() + 60]))
        elif not quiet:
            print("  ok  %s" % c["name"])
    for name, ok, detail in structure(b):
        if not ok:
            bad += 1
            if collect is not None:
                collect.append("構造:" + name)
            if not quiet:
                print("  ★構造: %s → %s" % (name, detail))
        elif not quiet:
            print("  ok  構造: %s" % name)
    for name, ok in expectations(b):
        if not ok:
            bad += 1
            if collect is not None:
                collect.append("期待値:" + name)
            if not quiet:
                print("  ★期待値が合いません: %s" % name)
        elif not quiet:
            print("  ok  %s" % name)
    return bad


def mutate(text, a, m):
    """§43 の【中】だけを変異させる。外（経緯の節）を書き換えても
    番人は落ちないので、外を変異させると空回りを見逃す（11 巡目の指摘）。"""
    k = text.index("## 43. ")
    if a not in text[k:]:
        return None
    return text[:k] + text[k:].replace(a, m, 1)


# 期待値の変異ごとの【主検出器】。巻き添えを検査済みに数えないために宣言する。
EXP_MAIN = {
    "面数 1,496 → 1,500": "A5 の対の上限 = 2*C(F,2)",
    "A5 の対の上限だけを変える": "A5 の対の上限 = 2*C(F,2)",
    "代表点の数だけを変える": "代表点 = 2*(F+F)",
    "段 6 の予算だけ増やす": "段の予算の合計 = 本文の合計 かつ <= D",
    "構成数 22 → 21": "構成数 = S+I+U",
    "14 構成 → 13 構成": "自己検査の構成数 = S+I",
    "S1 の行を消す": "S/I/U の行数が 1 以上で、構成数と整合",
    "162 ビット → 160 ビット": "ビット幅: 3*(7b+15)+(6b+13) = 文書の 625",
    "79 バイト → 78 バイト": "79 = ceil(625/8)",
    "81 バイト → 80 バイト": "81 = 3*ceil(162/8)+ceil(139/8)",
    "バイト数の行を消す": "バイト数の行が読めません（削除されたか、書式が変わった）",
    "保存値 S(A) の末尾を変える": "保存値 S(A)/S(B)/R が cp3_gmp_results.txt と一致（列番号も文書から）",
    # 列番号を壊すと「読めません」側が先に落ちるので、そちらを主検出器にする
    "保存値の列番号を壊す": "保存値を読めません（式が変わったか、証拠が無い）",
    "予言値 殻間の末尾を変える": "予言値: 殻間+2*内側 = S(B) かつ 内側 = -R",
    "§43.12 の 261 を 262 に": "§43.12 の 261 / 256 / 257 / 590 が実データと一致",
    "§43.12 の 257 を 258 に": "§43.12 の 261 / 256 / 257 / 590 が実データと一致",
    "判定 (4b) の意味行を消す": "判定は (1)(2-a)(2-b)(3)(4a)(4b)(4c) の 7 種が【意味の表】にある",
    "決定表の 1 行を消す": "決定表の順序が 1..N で連続し、最終行が (1)",
    "決定表の最終行を (4c) に": "決定表の最終行が (1)",
    "決定表に未定義の判定を入れる": "決定表の判定値がすべて定義済みの 7 種",
    "(1) を途中の行にも置く": "(1) は最終行にだけ現れる",
    "レイ方向を重複させる": "レイの方向が、本文の本数と一致し、相異なり、零でない",
    "4 演算の面数を 4/3 でなくする": "4 演算の面数 = 3 演算 * 4/3",
    "R の一次定義の符号を反転": "式が変わっていない: R の一次定義",
    "P3′ の右辺を差に": "式が変わっていない: P3′ の右辺",
    "E1 の式を変える": "式が変わっていない: E1",
    "E2 の符号を + に": "式が変わっていない: E2 の符号と閾値",
    "E2 の閾値を w>=1 に": "式が変わっていない: E2 の符号と閾値",
    "E3 の右辺を -R に戻す": "式が変わっていない: E3 の右辺は -R_B（定義の行）",
    "E4 の式を変える": "式が変わっていない: E4",
    "領域の体積の式を和に": "式が変わっていない: 領域の体積の式",
    "k_0 を log10 に": "式が変わっていない: k_0 の式",
    "A5 の対の数を積に": "式が変わっていない: A5 の対の数",
    "|S_A| の符号を反転": "式が変わっていない: |S_A| の式",
    "存在しない §43.N を参照": "§43.N への参照がすべて実在する",
    "存在しない順序へ送る": "「順序 N」で送る先の行が実在する",
    "存在しない小節を参照": "§43.N(x) への参照がすべて実在する",
    "小節を飛ばす（(f) を (h) に）": "小節が (a) から連続している",
    "I1 の共線面の向きを戻す": "I の構成が、期待する停止箇所で止まる（幾何を計算）",
    "I1 の座標を壊す": "I の構成が、期待する停止箇所で止まる（幾何を計算）",
    "I4 の面の向きを反転": "I の構成が、期待する停止箇所で止まる（幾何を計算）",
    "U8 の回数を 1 に戻す": "U 表の呼び出し回数が、方向の本数と整合する",
    "fill を実コードと違える": "文書の値が実コードと一致する",
    "コードの引用の行番号を壊す": "コードの引用が実ファイルと一致する",
    "実ファイルの外を引く": "コードの引用が実ファイルと一致する",
    "I 表から面集合の併記を落とす": "I 表の各行に、面集合を併記した停止箇所と、構成がある",
    "I の構成を 1 つ消す": "I 表の各行に、面集合を併記した停止箇所と、構成がある",
    "レイの本数を食い違わせる": "レイの本数が、本文・表・「すべて評価」で一致する",
    "段 4 の標本数を食い違わせる": "段 4 の標本数が内部で整合する",
    "桁数を実データと違える": "列 29 の分子・分母の桁数が実データと一致",
    "種を実データと違える": "種 = 1000 + 模型の添字 が実データと一致",
    "S7 の和を違える": "S7 の S_内 の和が合う",
    "E3 反例の R_B を違える": "E3 の反例の R_B と ΣS(内側) が再計算と一致",
    "D と最悪の総時間を食い違わせる": "最悪の総時間 = D + 回収上限",
    "5 対の面数を違える": "5 対の面数が実データと一致",
    "351 MB を違える": "351 / 307 MB が 60F / 52.5F と一致",
    "手順 6 の mid を消す": "必ず在るべき項目が消えていない",
    "手順 8 を消す": "必ず在るべき項目が消えていない",
    "評価順を入れ替える": "必ず在るべき項目が消えていない",
    "I4 の停止箇所から併記を落とす": "I 表の各行に、面集合を併記した停止箇所と、構成がある",
    "丸めモードの注を消す": "必ず在るべき項目が消えていない",
    "KMSH の書式を消す": "必ず在るべき項目が消えていない",
    "経路 P のクランプを消す": "必ず在るべき項目が消えていない",
    "段 4 の種を 0 に戻す": "記号の取り残し（一般化の波及漏れ）",
    "A1〜A4 の T→S の順を消す": "必ず在るべき項目が消えていない",
    "桁数の誤りを戻す": "記号の取り残し（一般化の波及漏れ）",
    "必須項目を消す（§43.13 の共線面）": "必ず在るべき項目が消えていない",
    "必須項目を消す（順序 7 の S_j=0）": "必ず在るべき項目が消えていない",
    "参照の固定を破る（採用領域の参照を (f) に戻す）": "参照の固定（改番の取り残し）",
    "記号の取り残しを作る（E2 を B 専用に戻す）": "記号の取り残し（一般化の波及漏れ）",
    "measure の手順を B 専用に戻す": "記号の取り残し（一般化の波及漏れ）",
    "存在しない順序へ送る": "「順序 N」で送る先の行が実在する",
}


def selftest(text):
    """各検査について、変異を入れたら【その検査だけが】落ちることを確かめる。"""
    base = run(text, quiet=True)
    if base:
        print("★ 変異前に既に %d 件落ちています。先にそちらを直してください。" % base)
        print("   （このまま変異試験に進むと、網が黙って縮みます）")
        return 1
    fail = 0
    covered = set()
    for c in CHECKS:
        a, m = c["mutate"]
        mut = mutate(text, a, m)
        if mut is None:
            print("★ %s: 変異の適用先 %r が §43 に見つかりません" % (c["name"], a[:40]))
            fail += 1; continue
        got = []
        if run(mut, quiet=True, collect=got) == 0:
            print("★ %s: 変異を入れても落ちません（空回り）" % c["name"])
            fail += 1
        elif "検査:" + c["name"] not in got:
            # 意図した検出器が発火しなかった（別の検査が拾っただけ）。
            # CLAUDE.md「変異ごとに主検出器も記録してください」
            print("★ %s: 主検出器が発火しません → %s" % (c["name"], got))
            fail += 1
        else:
            covered.add("検査:" + c["name"])   # 巻き添えは数えない（EXP_MUT と揃える）
            print("  ok  %-30s — 主検出器が発火（%s）" % (c["name"], "／".join(got)))
    # 期待値にも変異を当てる（11 巡目の指摘。CHECKS だけでは 6 件が無検査だった）
    EXP_MUT = [
        ("面数 1,496 → 1,500", "面数はどちらも 1,496", "面数はどちらも 1,500"),
        ("段 6 の予算だけ増やす", "| 300 秒 |", "| 400 秒 |"),
        ("構成数 22 → 21", "**22 構成**", "**21 構成**"),
        ("162 ビット → 160 ビット", "**162 ビット**", "**160 ビット**"),
        ("79 バイト → 78 バイト", "真のビット詰め 79 バイト", "真のビット詰め 78 バイト"),
        ("判定 (4b) の意味行を消す", "| **(4b)** | 判定不能（**こちらの未対応**） |\n", ""),
        ("決定表の 1 行を消す", "| **8** | **$B$ に内側成分が 0 個**", "| **99** | **$B$ に内側成分が 0 個**"),
        ("レイ方向を重複させる", "$(1,1,2)$", "$(1,1,1)$"),
        ("4 演算の面数を 4/3 でなくする", "| 369,156,480 |", "| 369,156,481 |"),
        ("A5 の対の上限だけを変える", "$2{,}236{,}520$ 対で、これが上限",
                                        "$2{,}236{,}521$ 対で、これが上限"),
        ("代表点の数だけを変える", "共線面が 0 枚なら 5,984 個", "共線面が 0 枚なら 5,985 個"),
        # --- 13 巡目: 構造（改番・参照の取り残し）ぶん ---
        ("存在しない §43.N を参照", "§43.5(f)", "§43.99"),
        ("存在しない小節を参照", "§43.5(g)", "§43.5(z)"),
        ("小節を飛ばす（(f) を (h) に）", "#### (f) 領域の列挙", "#### (h) 領域の列挙"),
        ("I1 の共線面の向きを戻す", "`(0,2,1) (0,1,3) (1,2,3) (2,0,3)`",
                                     "`(0,1,2) (0,1,3) (1,2,3) (2,0,3)`"),
        ("I1 の座標を壊す", "(25,10,0)`", "(25,0,0)`"),
        ("I4 の面の向きを反転", "`(1,2,3) (0,3,2)", "`(1,3,2) (0,3,2)"),
        ("U8 の回数を 1 に戻す", "**16 本とも無効と判定し、有効 0 本を返す** | 16 |",
                                 "**有効 0 本**を返す | 1 |"),
        ("fill を実コードと違える", "$\\mathrm{fill} = 0.6$", "$\\mathrm{fill} = 0.5$"),
        ("コードの引用の行番号を壊す", "config.hpp:36-37", "config.hpp:40-41"),
        ("実ファイルの外を引く", "`loader.hpp:44-65`", "`loader.hpp:44-99999`"),
        ("I 表から面集合の併記を落とす", "| I3 | 1 辺に 3 枚の面 | **A1(T)** |", "| I3 | 1 辺に 3 枚の面 | **A1** |"),
        ("I の構成を 1 つ消す", "> **I5**: 基準から面を 1 枚取り除く。", ""),
        ("レイの本数を食い違わせる", "レイの 16 方向", "レイの 20 方向"),
        ("段 4 の標本数を食い違わせる", "各 5,000 対、合わせて 10,000 対", "各 5,000 対、合わせて 12,000 対"),
        ("桁数を実データと違える", "**分子 3,899 桁 / 分母 3,882 桁**", "**分子 3,900 桁 / 分母 3,882 桁**"),
        ("種を実データと違える", "添字 4232 → 種 5232", "添字 4232 → 種 5233"),
        ("S7 の和を違える", "$1296+48=1344$", "$1296+48=1444$"),
        ("E3 反例の R_B を違える", "$R_B = -1248$", "$R_B = -1348$"),
        ("D と最悪の総時間を食い違わせる", "**$D = 900$ 秒**", "**$D = 800$ 秒**"),
        ("5 対の面数を違える", "| 5,851,180 |", "| 5,851,181 |"),
        ("351 MB を違える", "**351 MB**", "**451 MB**"),
        ("手順 6 の mid を消す", r"$(p_i[k] - \mathrm{mid}[k]) \times", r"$(p_i[k]) \times"),
        ("手順 8 を消す", "| 8 | **同一格子点の併合**（下） |", ""),
        ("評価順を入れ替える", "A1 → A2 → A3 → A4 → A5", "A3 → A2 → A1 → A4 → A5"),
        ("I4 の停止箇所から併記を落とす", "| I4 | 頂点で 2 つの錐が接する | **A3(T)**", "| I4 | 頂点で 2 つの錐が接する | **A3**"),
        ("丸めモードの注を消す", "（最近接・偶数優先）", ""),
        ("KMSH の書式を消す", 'マジック `"KMSH"`（4 バイト）', "マジック（省略）"),
        ("経路 P のクランプを消す", "に【クランプ】する", "にする"),
        ("段 4 の種を 0 に戻す", "> **種は 1**", "> **種は 0**"),
        ("A1〜A4 の T→S の順を消す", "**A1(T)→A2(T)→A3(T)→A4(T)→A1(S)→…**", "**順は任意**"),
        ("桁数の誤りを戻す", "**分子 3,899 桁 / 分母 3,882 桁**", "**分子が 7,782 桁**"),
        ("必須項目を消す（§43.13 の共線面）",
         "- **共線面があったときに判定を続ける機構**（**未設計**。", "- **（消した）**（"),
        ("必須項目を消す（順序 7 の S_j=0）",
         "**／ どれかの成分で $S_j = 0$**", "**／ （消した）**"),
        ("参照の固定を破る（採用領域の参照を (f) に戻す）",
         "採用領域が $w>0$ だからです**（§43.5(g)）", "採用領域が $w>0$ だからです**（§43.5(f)）"),
        ("記号の取り残しを作る（E2 を B 専用に戻す）",
         "| **E2** | **$R_X = -", "| **E2** | **$R_B = -"),
        ("measure の手順を B 専用に戻す",
         r"**$\mathrm{measure}(w_X>0) = \sum", r"**$\mathrm{measure}(w_B>0) = \sum"),
        ("存在しない順序へ送る", "（§43.7 の順序 6）", "（§43.7 の順序 99）"),
        ("§43.12 の 261 を 262 に", "| **条件を満たす【行】** | **261** |",
                                    "| **条件を満たす【行】** | **262** |"),
        ("§43.12 の 257 を 258 に", "の印がある行** | **257** |", "の印がある行** | **258** |"),
        ("保存値の列番号を壊す", "| $S(A)$ | **33** |", "| $S(A)$ | **35** |"),
        ("バイト数の行を消す", "**頂点: リム表現 96 バイト / バイト境界丸め 81 バイト / 真のビット詰め 79 バイト**",
                               "**頂点の大きさは省略します**"),
        # --- 12 巡目の指摘で足した検査ぶん ---
        ("14 構成 → 13 構成", "（14 構成）", "（13 構成）"),
        ("S1 の行を消す", "| S1 | **同方向の入れ子**", "| Sx | **同方向の入れ子**"),
        ("81 バイト → 80 バイト", "バイト境界丸め 81 バイト", "バイト境界丸め 80 バイト"),
        ("保存値 S(A) の末尾を変える", "`47157422889525707`", "`47157422889525708`"),
        ("予言値 殻間の末尾を変える", "`222847097096154239`", "`222847097096154238`"),
        ("R の一次定義の符号を反転", r"- (S(A)+S(B))$$", r"+ (S(A)+S(B))$$"),
        ("P3′ の右辺を差に", r"\mu(\mathcal A)+\mu(\mathcal B)$ |", r"\mu(\mathcal A)-\mu(\mathcal B)$ |"),
        ("E1 の式を変える", "**E1** | **$R_A + R_B = R$**", "**E1** | **$R_A - R_B = R$**"),
        ("E2 の符号を + に", r"R_X = -\sum_{\text{領域}:\,w_X\ge2}", r"R_X = +\sum_{\text{領域}:\,w_X\ge2}"),
        ("E2 の閾値を w>=1 に", r"\text{領域}:\,w_X\ge2}(w_X-1)", r"\text{領域}:\,w_X\ge1}(w_X-1)"),
        ("E3 の右辺を -R に戻す", r"**E3** | **$\sum_j S(\text{内側成分}_j) = -R_B$**", r"**E3** | **$\sum_j S(\text{内側成分}_j) = -R$**"),
        ("E4 の式を変える", r"**E4** | **$\sum_j S_j = S(X)$**", r"**E4** | **$\sum_j S_j = S(A)$**"),
        ("領域の体積の式を和に", r"\mathrm{vol}(\text{領域}_j) = \lvert S_j \rvert - \sum",
                                r"\mathrm{vol}(\text{領域}_j) = \lvert S_j \rvert + \sum"),
        ("k_0 を log10 に", r"k_0 = \lceil \log_2 \max_i", r"k_0 = \lceil \log_{10} \max_i"),
        ("A5 の対の数を積に", r"\binom{\lvert S_A\rvert}{2}+\binom{\lvert S_B\rvert}{2}",
                             r"\binom{\lvert S_A\rvert}{2}\cdot\binom{\lvert S_B\rvert}{2}"),
        ("|S_A| の符号を反転", r"\lvert S_A \rvert = 1496 - c_A", r"\lvert S_A \rvert = 1496 + c_A"),
        ("決定表の最終行を (4c) に", "| **10** | **上のどれにも当たらない**（＝ 内側成分が 1 個以上あり、E1 が厳密一致し、自己検査も通った） | **(1)** |",
                                     "| **10** | **上のどれにも当たらない**（＝ 内側成分が 1 個以上あり、E1 が厳密一致し、自己検査も通った） | **(4c)** |"),
        ("決定表に未定義の判定を入れる", "| **(4a)** |\n| **5**", "| **(5z)** |\n| **5**"),
        ("(1) を途中の行にも置く", "| **8** | **$B$ に内側成分が 0 個**（**H′ が成立しない**） | **(2-a)** |",
                                   "| **8** | **$B$ に内側成分が 0 個**（**H′ が成立しない**） | **(1)** |"),
    ]
    for name, a, m in EXP_MUT:
        mut = mutate(text, a, m)
        if mut is None:
            print("★ 期待値の変異 %s: 適用先 %r が §43 に見つかりません" % (name, a[:36]))
            fail += 1; continue
        got = []
        if run(mut, quiet=True, collect=got) == 0:
            print("★ 期待値の変異 %s: 落ちません（空回り）" % name); fail += 1
        else:
            # 主検出器を【宣言】し、それが発火したときだけ検査済みにする
            # （13 巡目の指摘 H3。巻き添えを「検査済み」に数えない）
            main = EXP_MAIN.get(name)
            if main is None:
                print("★ 期待値の変異 %s: 主検出器が宣言されていません" % name); fail += 1
            elif not any(g.endswith(main) for g in got):
                print("★ 期待値の変異 %s: 主検出器 %r が発火しません → %s"
                      % (name, main, got)); fail += 1
            else:
                covered.update(g for g in got if g.endswith(main))
                print("  ok  期待値の変異 %-28s — 主検出器が発火（%s）" % (name, main))

    # 表の検査も変異させる
    k = text.index("## 43. ")   # §43 の中の区切り行を消す
    mut = text[:k] + text[k:].replace("|---|---|\n", "", 1)
    if run(mut, quiet=True) == 0:
        print("★ 表の検査: 区切り行を消しても落ちません"); fail += 1
    else:
        print("  ok  表の検査 — 区切り行の削除を検出")
    # 変異が 1 つも当たっていない検査を落とす（12 巡目の指摘。
    # docstring が「足さなければ selftest が落ちる」と書いていたのに未実装だった）
    all_names = (["検査:" + c["name"] for c in CHECKS]
                 + ["構造:" + n for n, _, _ in structure(section43(text))]
                 + ["期待値:" + n for n, _ in expectations(section43(text))])
    dead = sorted(set(EXP_MAIN) - {n for n, _, _ in EXP_MUT})
    if dead:
        print("★ EXP_MAIN に、対応する変異が無い鍵が %d 件あります:" % len(dead))
        for d in dead:
            print("     %s" % d)
        fail += len(dead)
    naked = [n for n in all_names if n not in covered]
    if naked:
        print("★ 変異が登録されていない検査が %d 件あります（空回りか判別できません）:"
              % len(naked))
        for n in naked:
            print("     %s" % n)
        fail += len(naked)
    print("変異 %d 件 / 検出できなかったもの・無検査 %d 件"
          % (len(CHECKS) + len(EXP_MUT) + 1, fail))
    return fail


if __name__ == "__main__":
    t = io.open(DOC, encoding="utf-8").read()
    sys.exit(selftest(t) if "--selftest" in sys.argv else run(t))

#!/usr/bin/env python3
"""合成入力による、幾何検査の部品の試験。

**この試験の入力と期待値は、このファイルが原本です。** 文書と二重に持ちません。
期待値は、箱の解析値か、手で確認できる局所構造（辺に接する面の枚数など）から与えます。

**通過しても「実入力の巻き数・包含・残差検査が完了した」ことにはなりません。**
**部品が合成入力で期待どおり動いたことまでです。**

実行枠: 計算＋終了猶予 90 秒 / 仮想アドレス空間 4 GiB。
**90 秒で実際に打ち切ります**（`signal.alarm`。戻ってから測るだけでは、
止まらない入力で永久に戻りません）。
"""
import os, resource, signal, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inspect_geom as g

BUDGET_S = 90.0
AS_BYTES = 4 * 1024**3

# --- 合成入力（コードが原本） --------------------------------------------
def cube(lo, hi, s=+1):
    return g.box((lo, lo, lo), (hi, hi, hi), s)


def open_surface():
    """開口: 立方体から面を 1 枚取り除く"""
    P, F = cube(0, 10)
    return P, F[1:]


def duplicate_face():
    """重複: 立方体に、同じ 3 頂点・同じ向きの面をもう 1 枚足す"""
    P, F = cube(0, 10)
    return P, F + [F[0]]


def degenerate():
    """退化: 立方体と、頂点を共有しない「潰れた四面体」（共線面 1 枚）"""
    P, F = cube(0, 10)
    n = len(P)
    P = P + [(20, 0, 0), (30, 0, 0), (25, 0, 0), (25, 10, 0)]
    F = F + [(n+0, n+2, n+1), (n+0, n+1, n+3), (n+1, n+2, n+3), (n+2, n+0, n+3)]
    return P, F


def nonmanifold():
    """非多様体: 原点を共有する 2 つの四面体"""
    P = [(0, 0, 0), (10, 0, 0), (0, 10, 0), (0, 0, 10),
         (-10, 0, 0), (0, -10, 0), (0, 0, -10)]
    F = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1),
         (4, 5, 6), (0, 6, 5), (0, 4, 6), (0, 5, 4)]
    return P, F


# --- 期待値（解析値・局所構造） ------------------------------------------
BOXES = [
    ("単一",            [((0, 10), +1)],                       6000,  6000,     0, 0),
    ("同方向の入れ子",   [((0, 10), +1), ((2, 8), +1)],         7296,  6000, -1296, 1),
    ("反転空洞",        [((0, 10), +1), ((2, 8), -1)],         4704,  4704,     0, 0),
    ("離れた 2 個",      [((0, 10), +1), ((20, 30), +1)],      12000, 12000,     0, 0),
]
A5_BOXES = [
    ("接触", (0, 0, 0), (10, 10, 10), (5, 0, 10), (15, 10, 20), "接触あり・交差なし"),
    ("交差", (0, 0, 0), (10, 10, 10), (5, 5, 5), (15, 15, 15),  "交差あり"),
]
BAD = [
    ("開口",     open_surface,  "A1", 3),     # 取り除いた面の 3 辺が 1 枚に減る
    ("重複",     duplicate_face, "A1", 3),    # 重ねた面の 3 辺が 3 枚になる
    ("退化",     degenerate,    "A0-3", 1),   # 共線面がちょうど 1 枚
    ("非多様体", nonmanifold,   "A3", 1),     # 原点のリンクが 2 閉路
]


def run():
    ok = fail = 0

    def chk(name, got, want):
        nonlocal ok, fail
        good = got == want
        ok, fail = ok + good, fail + (not good)
        print("  %-4s %-34s 得た値 %-22s 期待 %s"
              % ("ok" if good else "★NG", name, got, want))

    print("=== 1. 単純な箱 — 体積・領域・残差・内側成分 ===")
    for name, spec, sb, mu, r, ninner in BOXES:
        P, F = g.join(*[cube(a, b, s) for (a, b), s in spec])
        a0 = g.check_a0(P, F)
        chk(name + ": A0〜A4",
            (a0[0], a0[1], a0[2], g.check_a1(F), g.check_a2(F),
             g.check_a3(F), g.check_a4(F)), ([], [], [], [], [], [], []))
        chk(name + ": A5 の不正", g.check_a5(P, F)[:2], (0, 0))
        a = g.analyze(P, F)
        chk(name + ": S(B)", a["S_total"], sb)
        chk(name + ": measure(w>0)", a["measure"], mu)
        chk(name + ": 残差 R", a["R"], r)
        chk(name + ": 内側成分の個数", len(a["inner"]), ninner)
        if ninner:
            chk(name + ": ΣS(内側) = -R", a["S_inner"], -r)

    print("\n=== 2. A5 が破れる箱（接触・交差） ===")
    for name, l1, h1, l2, h2, note in A5_BOXES:
        P, F = g.join(g.box(l1, h1, +1), g.box(l2, h2, +1))
        chk(name + ": A1〜A4 は通る",
            (g.check_a1(F), g.check_a2(F), g.check_a3(F), g.check_a4(F)),
            ([], [], [], []))
        nx, nt, bad = g.check_a5(P, F)
        chk(name + ": 頂点の共有は 0", len(set(P[:8]) & set(P[8:])), 0)
        if name == "接触":
            chk(name + ": 交差は 0 件", nx, 0)
            chk(name + ": 接触が 1 件以上", nt > 0, True)
        else:
            chk(name + ": 交差が 1 件以上", nx > 0, True)
        print("       （内訳: 交差 %d 件 / 接触 %d 件 — %s）" % (nx, nt, note))

    print("\n=== 3. 小さな不正入力 — 期待する検査で止まるか ===")
    for name, make, where, n in BAD:
        P, F = make()
        oor, dup, col = g.check_a0(P, F)
        a1, a3 = g.check_a1(F), g.check_a3(F)
        if where == "A0-3":
            chk(name + ": 共線面の枚数", len(col), n)
            S = [f for f in F if f not in col]
            chk(name + ": T は A1 を通る", g.check_a1(F), [])
            chk(name + ": S の A1 が破れる辺の本数", len(g.check_a1(S)), 3)
        elif where == "A1":
            chk(name + ": 共線面なし", len(col), 0)
            chk(name + ": A1 が破れる辺の本数", len(a1), n)
        elif where == "A3":
            chk(name + ": A1・A2 は通る", (g.check_a1(F), g.check_a2(F)), ([], []))
            chk(name + ": A3 が破れる頂点の数", len(a3), n)
            chk(name + ": その頂点のリンクの閉路数", a3[0][1] if a3 else None, 2)

    print("\n=== 4. 深い入れ子（深さ 3）— 親の選び方が効く範囲 ===")
    # 期待値は手で出したもの。領域は [外殻と中殻の間 / 中殻と内殻の間 / 内殻の中]
    DEEP = [
        ("すべて +1", [+1, +1, +1], [1, 2, 3], 7344, 6000, -1344, 2),
        ("+1 / -1 / -1", [+1, -1, -1], [1, 0, -1], 4656, 4704, 48, 0),
    ]
    for name, sg, w_want, sb, mu, r, ninner in DEEP:
        P, F = g.join(*[cube(a, b, s) for (a, b), s
                        in zip([(0, 10), (2, 8), (4, 6)], sg)])
        a = g.analyze(P, F)
        chk(name + ": 親（最も内側の包含者）", a["parent"], [None, 0, 1])
        chk(name + ": 領域ごとの巻き数", a["region_w"], w_want)
        chk(name + ": 領域ごとの体積", a["region_vol"], [4704, 1248, 48])
        chk(name + ": S(B)", a["S_total"], sb)
        chk(name + ": measure(w>0)", a["measure"], mu)
        chk(name + ": 残差 R", a["R"], r)
        chk(name + ": 内側成分の個数", len(a["inner"]), ninner)

    print("\n=== 5. 巻き数の直接の試験（符号・退化・方向の取り直し） ===")
    P, F = cube(0, 10, +1)
    chk("外向きの箱の内部は +1", g.winding(P, F, (5, 5, 5))[0], 1)
    chk("外部は 0", g.winding(P, F, (20, 5, 5))[0], 0)
    Pm, Fm = cube(0, 10, -1)
    chk("内向きの箱の内部は -1", g.winding(Pm, Fm, (5, 5, 5))[0], -1)
    chk("面の上は測れない", g.winding(P, F, (0, 5, 5))[0], None)
    chk("辺の上は測れない", g.winding(P, F, (0, 0, 5))[0], None)
    chk("頂点の上は測れない", g.winding(P, F, (0, 0, 0))[0], None)
    # (5,3,3) からの +x は x=10 面の対角線に当たる → 先頭の方向が退化し、取り直す
    w, di, nv = g.winding(P, F, (5, 3, 3))
    chk("対角線に当たる点でも測れる", w, 1)
    chk("★ 方向の取り直しが【効く】（先頭以外を使う）", di != 0, True)
    chk("有効な方向は 2 本以上、かつ全方向ではない", 2 <= nv < len(g.DIRS), True)

    print("\n=== 6. 共面で面積を共有しない接触（s=0） ===")
    # 辺だけで触れる 2 箱。頂点は 1 つも共有しない（z をずらしてある）
    P, F = g.join(g.box((0, 0, 0), (10, 10, 10), +1),
                  g.box((10, 10, 3), (20, 20, 13), +1))
    chk("頂点の共有は 0", len(set(P[:8]) & set(P[8:])), 0)
    chk("A1〜A4 は通る", (g.check_a1(F), g.check_a2(F), g.check_a3(F),
                        g.check_a4(F)), ([], [], [], []))
    nx, nt, bad = g.check_a5(P, F)
    chk("交差は 0 件", nx, 0)
    # ★ 「共面・s=0・面積 0・交わりは空でない」対だけを取り出して判定を見る。
    #   nt > 0 だけだと、共面でない接触で通ってしまい、この経路を検査できない。
    cop0 = []
    for i in range(len(F)):
        for j in range(i + 1, len(F)):
            n1, n2 = g.normal(P, F[i]), g.normal(P, F[j])
            if g.cross(n1, n2) != (0, 0, 0):
                continue
            if g.dot(n1, P[F[j][0]]) != g.dot(n1, P[F[i][0]]):
                continue
            if set(F[i]) & set(F[j]):
                continue
            ax = max(range(3), key=lambda k: abs(n1[k]))
            ij = [k for k in range(3) if k != ax]
            poly = g._inter2([P[k] for k in F[i]], [P[k] for k in F[j]], ij)
            if poly and g._area2(poly, ij) == 0:
                cop0.append((i, j))
    chk("共面・s=0・面積 0・交わり非空の対がある", len(cop0) > 0, True)
    chk("★ その対をすべて【不正】と判定する",
        all(not g.pair_a5(P, F[i], F[j])[1] for i, j in cop0), True)
    print("       （内訳: 交差 %d 件 / 接触 %d 件、うち共面・面積 0 が %d 件）"
          % (nx, nt, len(cop0)))

    print("\n=== 7. 方向が食い違う入力（開いた曲面） ===")
    # 開いた曲面では、抜ける方向と抜けない方向で巻き数が違う → 測れないと返すべき
    Po, Fo = open_surface()
    w, di, nv = g.winding(Po, Fo, (5, 5, 5))
    chk("開いた曲面の内側では測れない", w, None)
    chk("有効な方向は 2 本以上ある（0 本ではなく、食い違いで None）", nv >= 2, True)
    chk("★ 方向に正負の両方がある", min(min(d) for d in g.DIRS) < 0, True)

    return ok, fail


def _on_alarm(sig, frm):
    raise TimeoutError("実行枠 %.0f 秒を超えました" % BUDGET_S)


if __name__ == "__main__":
    resource.setrlimit(resource.RLIMIT_AS, (AS_BYTES, AS_BYTES))
    signal.signal(signal.SIGALRM, _on_alarm)
    signal.alarm(int(BUDGET_S))
    t0 = time.monotonic()
    try:
        ok, fail = run()
    except TimeoutError as e:
        print("\n★ %s（打ち切りました）" % e)
        sys.exit(2)
    finally:
        signal.alarm(0)
    dt = time.monotonic() - t0
    print("\n通過 %d / 失敗 %d   経過 %.2f 秒（枠 %.0f 秒で打ち切り）"
          % (ok, fail, dt, BUDGET_S))
    sys.exit(1 if fail else 0)

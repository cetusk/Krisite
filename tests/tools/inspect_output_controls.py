#!/usr/bin/env python3
r"""案 I の段 I-0 — 合成対照（**壊した出力を落とすこと**を確かめます）。

**§44.8 の表に書いた正負の対照だけを回します。** 実データは読みません。
"""
import os, sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inspect_geom as ig
import inspect_output as io_

RES = []


def chk(name, got, want, log):
    ok = (got == want)
    RES.append(ok)
    log("  %s %-56s 得た値 %-30s 期待 %s" % ("ok  " if ok else "★NG", name, got, want))


def e1_violations(F):
    cnt = {}
    for f in F:
        for k in range(3):
            cnt[(f[k], f[(k + 1) % 3])] = cnt.get((f[k], f[(k + 1) % 3]), 0) + 1
    return [e for e in cnt if cnt[e] != cnt.get((e[1], e[0]), 0)]


def run(log):
    # ---- E0 ----
    C = ig.box((0, 0, 0), (10, 10, 10), +1)
    chk("E0 正: 正しい箱に違反なし", [f for f in C[1] if len(set(f)) < 3], [], log)
    chk("E0 負: 3 頂点が相異でない面を数える",
        len([f for f in C[1] + [(0, 0, 1)] if len(set(f)) < 3]), 1, log)

    # ---- E1 ----
    chk("E1 正: 閉じた箱は違反 0", len(e1_violations(C[1])), 0, log)
    chk("E1 負: 面を 1 枚落とすと違反が出る", len(e1_violations(C[1][:-1])) > 0, True, log)
    rev = list(C[1]); rev[0] = (rev[0][0], rev[0][2], rev[0][1])
    chk("E1 負: 1 枚だけ向きを反転すると違反が出る", len(e1_violations(rev)) > 0, True, log)

    # ---- 退化面 ----
    P2 = [(0, 0, 0), (10, 0, 0), (20, 0, 0), (0, 10, 0)]
    chk("退化面 正: 非退化な面は 0 枚",
        sum(1 for f in [(0, 1, 3)] if io_.degenerate(_h(P2), f)), 0, log)
    chk("退化面 負: 共線の 3 点は 1 枚",
        sum(1 for f in [(0, 1, 2)] if io_.degenerate(_h(P2), f)), 1, log)

    # ---- E4 の分類 ----
    P = [(0, 0, 0), (10, 0, 0), (0, 10, 0), (10, 10, 0),
         (2, 2, -5), (2, 2, 5), (8, 3, 0), (1, 1, 0), (9, 1, 0), (1, 9, 0),
         (5, 0, 0), (5, 0, 10), (9, 4, 5)]
    chk("E4 正: 辺を共有する 2 枚（共面）", io_.classify_e4(P, (0, 1, 2), (1, 2, 3)),
        "正規（共有する単体と一致）", log)
    chk("E4 負: 同じ 3 頂点の 2 枚", io_.classify_e4(P, (0, 1, 2), (0, 2, 1)),
        "不正（正の面積を共有）", log)
    chk("E4 負: 共面で重なる 2 枚", io_.classify_e4(P, (0, 1, 2), (7, 8, 9)),
        "不正（正の面積を共有）", log)
    chk("E4 負: 非共面で貫く 2 枚", io_.classify_e4(P, (0, 1, 2), (4, 5, 6)),
        "不正（真の交差）", log)
    chk("E4 未判定: 一方の頂点が他方の辺の内部に触れる（**非退化どうし**）",
        io_.classify_e4(P, (0, 1, 2), (10, 11, 12)), "未判定（3 分類に当たりません）", log)

    # ---- C1（合成の 2 立体。**離れた 2 箱なので、$\\cup$ の正解が手で作れます**）----
    log("  --- C1（合成。入力は離れた 2 箱。直線は両方の内部を通ります）---")
    A = ig.box((0, 0, 0), (10, 10, 10), +1)
    B = ig.box((20, 0, 0), (30, 10, 10), +1)
    # ★ y = z にすると、箱の面の三角形分割の対角線にちょうど載ります（辺に当たる）。
    #   **非退化な直線を選ぶには、対角線を外す必要があります。**
    o = (Fr(-100), Fr(11, 2), Fr(13, 3))
    d = (1, 0, 0)

    def ev(mesh):
        return io_.events(mesh[0], mesh[1], o, d)

    ea, eb = ev(A), ev(B)
    chk("C1 前提: 入力の直線が非退化", (ea is not None) and (eb is not None), True, log)

    def score(X):
        """(C1-c の区間数, C1-b の区間数)。退化なら None。"""
        ex = ev(X)
        if ex is None or ea is None or eb is None:
            return None
        ts = sorted({t for t, _ in ex + ea + eb})
        c = b_ = 0
        for k in range(len(ts) + 1):
            m = (ts[0] - 1 if k == 0 else
                 (ts[-1] + 1 if k == len(ts) else (ts[k - 1] + ts[k]) / 2))
            wx = sum(sg for t, sg in ex if t < m)
            wa = sum(sg for t, sg in ea if t < m)
            wb = sum(sg for t, sg in eb if t < m)
            want = 1 if (wa > 0 or wb > 0) else 0
            if (1 if wx > 0 else 0) != want:
                c += 1
            elif wx != want:
                b_ += 1
        return (c, b_)

    # 正: 離れた 2 箱の ∪ は、2 箱を並べたもの
    chk("C1 正: 正しい ∪ → (C1-c, C1-b) = (0, 0)", score(ig.join(A, B)), (0, 0), log)
    # 負（C1-c）: B の領域を丸ごと欠く（**E1 は満たしたまま**）
    sa = score(A)
    chk("C1 負: B の領域を欠く → C1-c が出る", sa is not None and sa[0] > 0, True, log)
    chk("C1 負: B の領域を欠いても E1 は通る", len(e1_violations(A[1])), 0, log)
    # 負（C1-b）: A の内部に余分な閉殻を足す（**[w>0] は同じで w=2**）
    extra = ig.join(ig.join(A, B), ig.box((2, 2, 2), (8, 8, 8), +1))
    se = score(extra)
    chk("C1 負: 余分な内側の閉殻 → C1-c は 0、C1-b が出る",
        None if se is None else (se[0], se[1] > 0), (0, True), log)
    chk("C1 負: 余分な内側の閉殻でも E1 は通る", len(e1_violations(extra[1])), 0, log)
    # 負（未評価）: 面を 1 枚落とす → E1 が破れるので C1 は未評価
    broken = (ig.join(A, B)[0], ig.join(A, B)[1][:-1])
    chk("C1 未評価: 面を 1 枚落とすと E1 が破れる（**C1-c にはしません**）",
        len(e1_violations(broken[1])) > 0, True, log)

    # ★ 判定を終了値へ届けます。**ログに書くだけでは駆動が成功と読みます。**
    #   （この 3 行を落としていて、`run(log)` が `None` を返していました。
    #    **失敗しても終了値 0 になる形**でした。）
    ng = RES.count(False)
    log("通過 %d / 失敗 %d" % (RES.count(True), ng))
    log("判定: %s" % ("合成対照は全件通過" if ng == 0 else "★ 失敗あり。後続を起動しません"))
    return 0 if ng == 0 else 1

def _h(P):
    """整数点を同次座標（w=1）に。"""
    return [(p[0], p[1], p[2], 1) for p in P]


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    f = open(os.path.join(out, "stage_controls.txt"), "w")

    def log(s):
        print(s)
        sys.stdout.flush()
        f.write(s + "\n")
        f.flush()

    log("=== 段 controls ===")
    try:
        rc = run(log)
    finally:
        f.close()
    sys.exit(rc)

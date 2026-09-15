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


# ★ 対照は【本体の実装】を呼びます。**写しを検査すると、本体が壊れたままになります。**
e1_violations = io_.e1_violations


def run(log):
    # ---- E0（**本体の `e0_violations` を呼びます**）----
    C = ig.box((0, 0, 0), (10, 10, 10), +1)
    HV = _h(C[0])
    chk("E0 正: 正しい箱に違反なし", io_.e0_violations(HV, C[1]), ([], [], []), log)
    chk("E0 負: 3 頂点が相異でない面", io_.e0_violations(HV, C[1] + [(0, 0, 1)])[1], [12], log)
    chk("E0 負: 範囲外の添字", io_.e0_violations(HV, C[1] + [(0, 1, 99)])[0], [12], log)
    bad_w = list(HV); bad_w[0] = (HV[0][0], HV[0][1], HV[0][2], 0)
    chk("E0 負: 同次座標の分母が 0", io_.e0_violations(bad_w, C[1])[2], [0], log)

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
    # ★ **非共面で、座標では辺を共有するが添字は別**（併合漏れの形）。
    #   **添字で判定すると「未判定」へ流れます。座標で判定すれば「正規」です。**
    Pm = [(0, 0, 0), (10, 0, 0), (0, 10, 0),          # 0..2
          (0, 0, 0), (10, 0, 0), (0, 0, 10)]          # 3..5（3,4 は 0,1 と同じ座標）
    chk("E4 正: 非共面・**座標で**辺を共有・添字は別 → 正規",
        io_.classify_e4(Pm, (0, 1, 2), (3, 4, 5)), "正規（共有する単体と一致）", log)
    chk("E4 正: 非共面・添字でも辺を共有 → 正規",
        io_.classify_e4([(0, 0, 0), (10, 0, 0), (0, 10, 0), (0, 0, 10)],
                        (0, 1, 2), (0, 1, 3)), "正規（共有する単体と一致）", log)

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
        """`(C1-c の区間数, C1-b の区間数)`。退化なら `None`。**本体の `compare_line` を呼びます。**"""
        ex = ev(X)
        if ex is None or ea is None or eb is None:
            return None
        r = io_.compare_line(ex, ea, eb, "union")
        return None if r is None else (r[0], r[1])

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
    # ---- `read_soup` の往復（**書式の書き手と読み手が一致するか**）----
    import struct, tempfile
    nx, nw = 3, 3
    VV = [(5, -7, 9, 1), (-(1 << 100), 3, -4, 2)]
    FF = [(0, 1, 0)]
    with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as tf:
        tf.write(struct.pack("<II", len(VV), len(FF)))
        for v in VV:
            for comp, n in ((v[0], nx), (v[1], nx), (v[2], nx), (v[3], nw)):
                u = comp & ((1 << (64 * n)) - 1)
                tf.write(struct.pack("<%dQ" % n, *[(u >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
                                                   for i in range(n)]))
        tf.write(struct.pack("<%dI" % (3 * len(FF)), *[c for t in FF for c in t]))
        tmp = tf.name
    got, gotf, _ = io_.read_soup(tmp, nx, nw)
    os.unlink(tmp)
    chk("read_soup 往復: 負の値と 100 ビット超も戻る", (got, gotf), (VV, FF), log)

    # ---- `c2_sample` の正負（**本体を呼びます**）----
    with tempfile.TemporaryDirectory() as td:
        ap = os.path.join(td, "a.bin")
        with open(ap, "wb") as f:
            f.write(struct.pack("<II", 3, 1))
            f.write(struct.pack("<9i", 0, 0, 0, 10, 0, 0, 0, 10, 0))
            f.write(struct.pack("<3I", 0, 1, 2))
        bp = os.path.join(td, "b.bin")
        with open(bp, "wb") as f:
            f.write(struct.pack("<II", 3, 1))
            f.write(struct.pack("<9i", 100, 0, 0, 110, 0, 0, 100, 10, 0))
            f.write(struct.pack("<3I", 0, 1, 2))

        class _A:
            in_a, in_b = ap, bp

        P2 = [(Fr(1), Fr(1), Fr(0)), (Fr(4), Fr(1), Fr(0)), (Fr(1), Fr(4), Fr(0)),
              (Fr(1), Fr(1), Fr(1)), (Fr(4), Fr(1), Fr(1)), (Fr(1), Fr(4), Fr(1)),
              (Fr(-5), Fr(-5), Fr(0)), (Fr(20), Fr(-5), Fr(0)), (Fr(-5), Fr(20), Fr(0))]
        chk("C2 正: 入力面を細分した断片",
            io_.c2_sample(P2, (0, 1, 2), _A).startswith("単一の入力三角形に収まる"), True, log)
        chk("C2 負: 別の平面に載せた面 → **違反 C2**",
            io_.c2_sample(P2, (3, 4, 5), _A).startswith("**違反 C2**"), True, log)
        chk("C2 負: はみ出した面 → 未判定（支持平面は一致）",
            io_.c2_sample(P2, (6, 7, 8), _A).startswith("未判定（支持平面は一致"), True, log)

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

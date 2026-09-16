#!/usr/bin/env python3
r"""段 J-0 — G1・G2・G4 の合成対照（**本体の実装を呼びます。写しを検査しません**）。

**この子から別のプロセスを起動しません**（`SPEC-phase5.md` §5.10.14.111）。
**実データは読みません。** G2 は合成の出力一式を作り、`stage_measure` を**直接呼びます**。
"""
import os, struct, sys, tempfile
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inspect_geom as ig
import inspect_output as io_

RES = []


def chk(name, got, want, log):
    ok = (got == want)
    RES.append(ok)
    log("  %s %-52s 得た値 %-26s 期待 %s" % ("ok  " if ok else "★NG", name, got, want))


def _h(P, w=1):
    """整数点を同次座標に（`w` を掛けた表現。**同じ点を別の `w` で表せます**）。"""
    return [(p[0] * w, p[1] * w, p[2] * w, w) for p in P]


# ===== G1 =====
def g1(log):
    log("--- G1 体積の既知値対照（**本体の `volume6_rational` を呼びます**）---")
    U = ig.box((0, 0, 0), (1, 1, 1), +1)
    C = ig.box((0, 0, 0), (10, 10, 10), +1)
    CI = ig.box((0, 0, 0), (10, 10, 10), -1)
    chk("G1-1 単位立方体, w=1", io_.volume6_rational(_h(U[0]), U[1]), 6, log)
    chk("G1-2 辺 10 の立方体", io_.volume6_rational(_h(C[0]), C[1]), 6000, log)
    chk("G1-3 向きを反転", io_.volume6_rational(_h(CI[0]), CI[1]), -6000, log)
    chk("G1-4 同じ点を w=7 で表す（値が変わらない）",
        io_.volume6_rational(_h(C[0], 7), C[1]), 6000, log)
    D = ig.join(ig.box((0, 0, 0), (12, 12, 12), +1), ig.box((40, 0, 0), (46, 6, 6), +1))
    chk("G1-5 2 成分の和（加法性のみ。位置関係は検査しません）",
        io_.volume6_rational(_h(D[0]), D[1]), 6 * (12**3 + 6**3), log)
    # 一般の整数頂点の四面体。**原点も 0 成分も含めません**
    # （**含めると `det3` の第 3 項が消え、変異が打ち消しで素通りします** — 実測）。
    # **期待値は被検査関数から作りません。** Sarrus の展開を別に書きます。
    T = [(1, 2, 3), (4, 1, 2), (2, 5, 1), (3, 2, 6)]
    TF = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]

    def det_edges(p, q, r, s_):
        a = tuple(q[k] - p[k] for k in range(3))
        b = tuple(r[k] - p[k] for k in range(3))
        c = tuple(s_[k] - p[k] for k in range(3))
        return (a[0] * b[1] * c[2] + a[1] * b[2] * c[0] + a[2] * b[0] * c[1]
                - a[2] * b[1] * c[0] - a[0] * b[2] * c[1] - a[1] * b[0] * c[2])

    det = det_edges(*T)
    chk("G1-6 一般の四面体（原点も 0 成分も含まない）",
        io_.volume6_rational(_h(T), TF), det, log)

    # 変異 2 件。**主検出器を指定します。**
    import inspect_output as m
    orig_det3, orig_pt = m.det3, m.pt
    try:
        m.det3 = lambda a, b, c: (a[0] * (b[1] * c[2] - b[2] * c[1])
                                  - a[1] * (b[0] * c[2] - b[2] * c[0]))   # 第 3 項を落とす
        chk("G1 変異 1（det3 の第 3 項を落とす）→ 主検出器 G1-6 が落ちる",
            m.volume6_rational(_h(T), TF) != det, True, log)
    finally:
        m.det3 = orig_det3
    try:
        m.pt = lambda v: (Fr(v[0]), Fr(v[1]), Fr(v[2]))                   # 分母を落とす
        chk("G1 変異 2（pt の分母を落とす）→ 主検出器 G1-4 が落ちる",
            m.volume6_rational(_h(C[0], 7), C[1]) != 6000, True, log)
    finally:
        m.pt = orig_pt


# ===== 合成の出力一式（§44.8 の書式）=====
def write_soup_py(path, V, F, nx=3, nw=3):
    with open(path, "wb") as f:
        f.write(struct.pack("<II", len(V), len(F)))
        for v in V:
            for comp, n in ((v[0], nx), (v[1], nx), (v[2], nx), (v[3], nw)):
                u = comp & ((1 << (64 * n)) - 1)
                f.write(struct.pack("<%dQ" % n,
                                    *[(u >> (64 * i)) & 0xFFFFFFFFFFFFFFFF for i in range(n)]))
        f.write(struct.pack("<%dI" % (3 * len(F)), *[c for t in F for c in t]))


def write_quant(path, P, F):
    with open(path, "wb") as f:
        f.write(struct.pack("<II", len(P), len(F)))
        f.write(struct.pack("<%di" % (3 * len(P)), *[c for p in P for c in p]))
        f.write(struct.pack("<%dI" % (3 * len(F)), *[c for t in F for c in t]))


class _Args:
    pass


def make_case(td, break_e1):
    r"""合成の出力一式・入力・保存物を作り、`_Args` を返す。

    **立体は $10^5$ 規模**にします（`c1_sample` は原点 $\pm200{,}000$・方向 $\pm99$ で
    直線を引き、**A と B の両方に当たること**を条件にするため。§45.3）。
    """
    A = ig.box((-60000, -60000, -60000), (60000, 60000, 60000), +1)
    B = ig.box((60001, -60000, -60000), (180000, 60000, 60000), +1)
    U = ig.join(A, B)                      # 離れた 2 箱なので ∪ の正解が手で作れます
    outs = {"union": U, "isect": A, "diff_ab": A, "diff_ba": B}
    d = os.path.join(td, "run")
    os.makedirs(d, exist_ok=True)
    for op, m in outs.items():
        V, F = _h(m[0]), list(m[1])
        if break_e1 and op == "union":
            F = F[:-1]                     # 面を 1 枚落とす → E1 が破れる
        write_soup_py(os.path.join(d, "b_out_%s.bin" % op), V, F)
    with open(os.path.join(d, "b_run_meta.txt"), "w") as f:
        f.write("kHomoXyz=3\nkHomoW=3\n")
    write_quant(os.path.join(td, "a.bin"), A[0], A[1])
    write_quant(os.path.join(td, "b.bin"), B[0], B[1])
    sv = os.path.join(td, "saved.txt")
    vol = io_.volume6_rational(_h(U[0]), U[1])
    cols = ["k"] + ["0"] * 27 + [str(vol)] + ["0"] * 6
    with open(sv, "w") as f:
        f.write(" ".join(cols) + "\n")
    a = _Args()
    a.out, a.read_from, a.sums = td, d, os.path.join(td, "none")
    a.in_a, a.in_b, a.saved, a.key = os.path.join(td, "a.bin"), os.path.join(td, "b.bin"), sv, "k"
    return a


def run_case(a, log):
    lines = []
    io_.stage_measure(a, lines.append)
    for x in lines:
        log("    | " + x)
    return lines


# ===== G2 =====
def g2(log):
    log("--- G2 $E1$ 違反 → C1 未評価（**`stage_measure` を直接呼びます**）---")
    for tag, brk in (("正", False), ("負（∪ から面を 1 枚落とす）", True)):
        with tempfile.TemporaryDirectory() as td:
            a = make_case(td, brk)
            # 照合は合成なので通りません。**幾何処理まで進めるため、照合を差し替えます。**
            ov = io_.verify_inputs
            io_.verify_inputs = lambda _a, _l: True
            try:
                lines = run_case(a, log)
            finally:
                io_.verify_inputs = ov
            txt = "\n".join(lines)
            chk("G2 %s: E0 が通る" % tag, "★ E0 が破れました" not in txt, True, log)
            if brk:
                chk("G2 %s: C1 が未評価" % tag, "C1: **未評価**" in txt, True, log)
                chk("G2 %s: C1-c / C1-b が立たない" % tag,
                    ("`C1-c`" in txt and "個 /" in txt) is False or "採用" not in txt, True, log)
                chk("G2 %s: 分岐到達を記録" % tag,
                    "**分岐到達**: E1 違反 → C1 未評価" in txt, True, log)
            else:
                chk("G2 %s: C1 が採用される" % tag, "C1: 試行" in txt, True, log)
                chk("G2 %s: 不一致 0" % tag, "`C1-c` 0 個 / `C1-b` 0 個" in txt, True, log)
                chk("G2 %s: 分岐到達を記録" % tag,
                    "**分岐到達**: E1 通過 → C1 評価" in txt, True, log)


# ===== G4 =====
def g4(log):
    log("--- G4 区間の端点・巻き数まで固定（**離れた 2 箱、手で選んだ直線**）---")
    A = ig.box((0, 0, 0), (10, 10, 10), +1)
    B = ig.box((20, 0, 0), (30, 10, 10), +1)
    o, d = (Fr(-100), Fr(11, 2), Fr(13, 3)), (1, 0, 0)
    ea, eb = io_.events(A[0], A[1], o, d), io_.events(B[0], B[1], o, d)
    ex = io_.events(*ig.join(A, B), o, d)
    r = io_.compare_line(ex, ea, eb, "union")
    chk("G4 正しい ∪ の (C1-c, C1-b)", (r[0], r[1]), (0, 0), log)
    chk("G4 区間数", r[2], 5, log)
    got = [(x[0], x[1], x[3], x[4], x[5], x[6], x[7]) for x in r[3]]
    want = [("-inf", "100", 0, 0, 0, 0, "ok"),
            ("100", "110", 1, 1, 0, 1, "ok"),
            ("110", "120", 0, 0, 0, 0, "ok"),
            ("120", "130", 1, 0, 1, 1, "ok"),
            ("130", "+inf", 0, 0, 0, 0, "ok")]
    chk("G4 区間の端点と巻き数まで一致", got, want, log)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    f = open(os.path.join(out, "stage_J0_controls.txt"), "w")

    def log(s):
        print(s)
        sys.stdout.flush()
        f.write(s + "\n")
        f.flush()

    log("=== 段 J-0（G1 / G2 / G4 の合成対照）===")
    try:
        g1(log)
        g2(log)
        g4(log)
        ng = RES.count(False)
        log("通過 %d / 失敗 %d" % (RES.count(True), ng))
        log("判定: %s" % ("合成対照は全件通過" if ng == 0 else "★ 失敗あり。後続を起動しません"))
    finally:
        f.close()
    sys.exit(0 if RES.count(False) == 0 else 1)

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
            if tag == "正":
                log("    ★ 合成の保存物の列 29 は `volume6_rational` の出力をそのまま"
                    "書いたものです。**ここでの「列 29 と一致」は恒等式**で、"
                    "**G1 の回帰の根拠は J-2（実データ）です。**")
            if brk:
                chk("G2 %s: C1 が未評価" % tag, "C1: **未評価**" in txt, True, log)
                chk("G2 %s: C1 の集計が無い（段の集計 None）" % tag,
                    "段の集計 None" in txt, True, log)
                chk("G2 %s: G4 未解消として記録される" % tag,
                    "**G4 未解消**" in txt, True, log)
                chk("G2 %s: 分岐到達を記録" % tag,
                    "**分岐到達**: E1 違反 → C1 未評価" in txt, True, log)
            else:
                chk("G2 %s: C1 が採用される" % tag, "C1: 試行" in txt, True, log)
                chk("G2 %s: 不一致 0" % tag, "`C1-c` 0 個 / `C1-b` 0 個" in txt, True, log)
                chk("G2 %s: 分岐到達を記録" % tag,
                    "**分岐到達**: E1 通過 → C1 評価" in txt, True, log)


# ===== G4 =====
def g4(log):
    log("--- G4 既定の 4 構成（**代表点も期待値に含めます**）と、読戻しの取りこぼし ---")
    A = ig.box((0, 0, 0), (10, 10, 10), +1)
    B = ig.box((20, 0, 0), (30, 10, 10), +1)
    o, d = (Fr(-100), Fr(11, 2), Fr(13, 3)), (1, 0, 0)
    ea, eb = io_.events(A[0], A[1], o, d), io_.events(B[0], B[1], o, d)

    def rows_of(mesh):
        ex = io_.events(mesh[0], mesh[1], o, d)
        r = io_.compare_line(ex, ea, eb, "union")
        return r

    # --- 構成 1: 正 ---
    r = rows_of(ig.join(A, B))
    chk("G4-1 正: (C1-c, C1-b)", (r[0], r[1]), (0, 0), log)
    chk("G4-1 正: 区間数", r[2], 5, log)
    chk("G4-1 正: 全 8 列（**代表点を含む**）",
        r[3],
        [("-inf", "100", "99", 0, 0, 0, 0, "ok"),
         ("100", "110", "105", 1, 1, 0, 1, "ok"),
         ("110", "120", "115", 0, 0, 0, 0, "ok"),
         ("120", "130", "125", 1, 0, 1, 1, "ok"),
         ("130", "+inf", "131", 0, 0, 0, 0, "ok")], log)

    # --- 構成 2: C1-c（B の領域を丸ごと欠く。**E1 は満たしたまま**）---
    r = rows_of(A)
    chk("G4-2 C1-c: (C1-c, C1-b)", (r[0], r[1]), (1, 0), log)
    chk("G4-2 C1-c: 全 8 列",
        r[3],
        [("-inf", "100", "99", 0, 0, 0, 0, "ok"),
         ("100", "110", "105", 1, 1, 0, 1, "ok"),
         ("110", "120", "115", 0, 0, 0, 0, "ok"),
         ("120", "130", "125", 0, 0, 1, 1, "C1-c"),
         ("130", "+inf", "131", 0, 0, 0, 0, "ok")], log)

    # --- 構成 3: C1-b（A の内部に余分な閉殻。[w>0] は同じで w=2）---
    r = rows_of(ig.join(ig.join(A, B), ig.box((2, 2, 2), (8, 8, 8), +1)))
    chk("G4-3 C1-b: (C1-c, C1-b)", (r[0], r[1]), (0, 1), log)
    chk("G4-3 C1-b: 区間数", r[2], 7, log)
    chk("G4-3 C1-b: C1-b の行（代表点と巻き数）",
        [x for x in r[3] if x[7] == "C1-b"],
        [("102", "108", "105", 2, 1, 0, 1, "C1-b")], log)

    # --- 構成 4: 未評価（E1 が破れる。一覧は見出しだけ）---
    broken = (ig.join(A, B)[0], ig.join(A, B)[1][:-1])
    chk("G4-4 未評価: E1 が破れる", len(io_.e1_violations(broken[1])) > 0, True, log)
    with tempfile.TemporaryDirectory() as td:
        q = os.path.join(td, "iv.txt")
        io_.write_intervals(q, [])                       # 未評価のときの書き出し
        chk("G4-4 未評価: 一覧は見出しだけ（読戻しは 0 行）", io_.read_intervals(q), [], log)
        chk("G4-4 未評価: ファイルは空でない",
            os.path.getsize(q) > 0, True, log)

    # --- 読戻しが【正常行の欠落】を検出すること ---
    rows = rows_of(ig.join(A, B))[3]
    with tempfile.TemporaryDirectory() as td:
        q = os.path.join(td, "iv.txt")
        io_.write_intervals(q, rows)
        back = io_.read_intervals(q)
        want = [tuple(str(x) for x in r) for r in rows]
        chk("G4-5 往復: 行数と各行が一致", (len(back) == len(want), back == want), (True, True), log)
        # **正常行を 1 行落とす**。件数（C1-c / C1-b）は 0 対 0 のままです。
        io_.write_intervals(q, rows[:2] + rows[3:])
        back2 = io_.read_intervals(q)
        c2 = (sum(1 for t in back2 if t[-1] == "C1-c"), sum(1 for t in back2 if t[-1] == "C1-b"))
        chk("G4-5 欠落: 件数の比較では検出できない（0 対 0 のまま）", c2, (0, 0), log)
        chk("G4-5 欠落: **行数の比較で検出できる**", len(back2) == len(want), False, log)
        chk("G4-5 欠落: **各行の比較でも検出できる**", back2 == want, False, log)


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
        # **件数は式で持ちます**（実測した数を書かない。`CLAUDE.md`）。
        n_g1 = 6 + 2                    # 対照 6 + 変異 2
        n_g2 = (1 + 3) + (1 + 4)        # 正: E0 + 3 件 / 負: E0 + 4 件
        n_g4 = (3 + 2 + 3 + 3) + 4     # 正 3 / C1-c 2 / C1-b 3 / 未評価 3 / 往復と欠落 4
        want_n = n_g1 + n_g2 + n_g4
        log("対照の件数 %d（期待 %d = G1 %d + G2 %d + G4 %d）"
            % (len(RES), want_n, n_g1, n_g2, n_g4))
        if len(RES) != want_n:
            RES.append(False)
            ng += 1
            log("  ★NG 件数が合いません。**呼ばれなかった対照があります。**")
        log("通過 %d / 失敗 %d" % (RES.count(True), ng))
        log("判定: %s" % ("合成対照は全件通過" if ng == 0 else "★ 失敗あり。後続を起動しません"))
    finally:
        f.close()
    sys.exit(0 if RES.count(False) == 0 else 1)

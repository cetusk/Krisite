#!/usr/bin/env python3
r"""③ 限定検査の合成対照 — `REVIEW-phase5-correctness-literature.md` §26.4。

**`inspect_index3` の内部関数 `run_manifest` を直接呼びます。**
**この子から別のプロセスを起動しません。実データは読みません。**

15 の対照群を、**静的なケース表**に列挙します。**予定ケースと実行ケースの鍵集合を照合**し、
件数の合計だけで呼出し漏れを見逃さない形にします。
"""
import builtins, hashlib, io, os, shutil, struct, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inspect_geom as ig
import inspect_output as io_
import inspect_index3 as i3

RES = []
RAN = []


def chk(name, got, want, log):
    ok = (got == want)
    RES.append(ok)
    log("  %s %-56s 得た値 %-30s 期待 %s" % ("ok  " if ok else "★NG", name, got, want))


# ---- 合成の書き出し ----------------------------------------------------
def _limbs(v, n):
    m = (1 << (64 * n)) - 1
    v &= m
    return [(v >> (64 * i)) & 0xFFFFFFFFFFFFFFFF for i in range(n)]


def write_soup(path, V4, F, nx=3, nw=3):
    """`read_soup` が読む書式（§44.8）。`V4` は同次 4 成分の整数。"""
    out = bytearray(struct.pack("<II", len(V4), len(F)))
    for (X, Y, Z, W) in V4:
        for val, n in ((X, nx), (Y, nx), (Z, nx), (W, nw)):
            for w in _limbs(val, n):
                out += struct.pack("<Q", w)
    for f in F:
        out += struct.pack("<III", f[0], f[1], f[2])
    b = bytes(out)
    with open(path, "wb") as fp:
        fp.write(b)
    return b


def homo(P, w=1):
    """整数 3 点を同次 4 成分へ（`w` 倍の表現。**同じ点を別の `w` で表せます**）。"""
    return [(p[0] * w, p[1] * w, p[2] * w, w) for p in P]


def two_boxes():
    """離れた 2 箱（正常例）。"""
    return ig.join(ig.box((0, 0, 0), (10, 10, 10), +1),
                   ig.box((20, 0, 0), (30, 10, 10), +1))


def merge_index(P, F, pairs):
    """`pairs` の `(b, a)` について、面索引の `b` を `a` へ読み替えます。

    **座標では併合しません。** 頂点配列はそのまま残し、**索引だけ**を同一視します。
    """
    m = dict(pairs)
    return P, [tuple(m.get(x, x) for x in f) for f in F]


def find(P, q):
    return P.index(tuple(q))


# ---- 合成マニフェストの組み立て ----------------------------------------
def make_dir(td, soups, meta_text="kHomoXyz=3\nkHomoW=3\n", limbs=(3, 3),
             sha_override=None, meta_missing=False, truncate=None):
    r"""`soups` は `{op: (V4, F)}`。戻り値は `man`。

    **期待ハッシュは、書いたバイト列から先に固定します**（§26.4）。
    `sha_override` は**同一性不成立の対照だけ**で使います。
    """
    d = os.path.join(td, "run")
    os.makedirs(d, exist_ok=True)
    out = {}
    for op in i3.OPS:
        V4, F = soups[op]
        p = os.path.join(d, "b_out_%s.bin" % op)
        b = write_soup(p, V4, F)
        if truncate and truncate[0] == op:
            b = b[:truncate[1]]
            with open(p, "wb") as fp:
                fp.write(b)
        out[op] = (p, hashlib.sha256(b).hexdigest())
    if sha_override:
        op, h = sha_override
        out[op] = (out[op][0], h)
    mp = os.path.join(d, "b_run_meta.txt")
    if meta_missing:
        mp = os.path.join(d, "no_such_meta.txt")
        msha = "0" * 64
    else:
        mb = meta_text.encode()
        with open(mp, "wb") as fp:
            fp.write(mb)
        msha = hashlib.sha256(mb).hexdigest()
    return {"meta": (mp, msha), "limbs": limbs, "out": out}


def normal_soups():
    P, F = two_boxes()
    return dict((op, (homo(P), list(F))) for op in i3.OPS)


def soups_with(bad_op, V4, F):
    s = normal_soups()
    s[bad_op] = (V4, F)
    return s


# ---- 実行 --------------------------------------------------------------
def run_case(key, man, td, log, patches=None, out_exists=False):
    """`run_manifest` を 1 回呼び、`(rc, 本文, 呼出し回数)` を返します。"""
    RAN.append(key)
    calls = {"degenerate": 0, "A1": 0, "A2": 0, "table": 0}
    ov = {"d": io_.degenerate, "a1": ig.check_a1, "a2": ig.check_a2,
          "t": i3.build_link_table, "open": builtins.open}

    def wrap(nm, fn):
        def g(*a, **k):
            calls[nm] += 1
            return fn(*a, **k)
        return g

    io_.degenerate = wrap("degenerate", ov["d"])
    ig.check_a1 = wrap("A1", ov["a1"])
    ig.check_a2 = wrap("A2", ov["a2"])
    i3.build_link_table = wrap("table", ov["t"])
    if patches:
        for k, v in patches.items():
            if k == "table":
                i3.build_link_table = wrap("table", v(ov["t"]))
            elif k == "open":
                builtins.open = v(ov["open"])
    lines = []
    outdir = os.path.join(td, "evi")
    if out_exists:
        os.makedirs(outdir, exist_ok=True)
    try:
        rc = i3.run_manifest(man, outdir, lines.append, {"synthetic": key})
    finally:
        io_.degenerate, ig.check_a1, ig.check_a2 = ov["d"], ov["a1"], ov["a2"]
        i3.build_link_table, builtins.open = ov["t"], ov["open"]
    for x in lines:
        log("    | " + x)
    return rc, "\n".join(lines), calls


def final(txt):
    xs = [x for x in txt.split("\n") if x.startswith("判定: ")]
    return xs[-1] if xs else "**判定行なし**"


# ---- 15 の対照群 -------------------------------------------------------
def g01_normal(log):
    P, F = two_boxes()
    Pt, Ft = ig.join(ig.box((0, 0, 0), (10, 10, 10), +1),
                     ig.box((10, 10, 10), (20, 20, 20), +1))
    for key, soups in (("01a-離れた2箱", normal_soups()),
                       ("01b-点接触だが別索引", dict((op, (homo(Pt), list(Ft)))
                                                     for op in i3.OPS))):
        with tempfile.TemporaryDirectory() as td:
            rc, txt, c = run_case(key, make_dir(td, soups), td, log)
            chk("%s: 戻り値" % key, rc, 0, log)
            chk("%s: 最終判定行" % key, final(txt).startswith("判定: ③ 限定検査 成功"), True, log)
            chk("%s: link まで到達" % key, c["table"], 4, log)
            # **旧方式との一致**（小対照に限る比較基準。実データでは二重実行しません）
            V4, FF = soups["union"]
            chk("%s: 旧 check_a3 と一致" % key, ig.check_a3(FF), [], log)


def g02_index(log):
    P, F = two_boxes()
    for key, FF in (("02a-範囲外索引", [(0, 1, 999)] + list(F[1:])),
                    ("02b-同一索引の反復", [(0, 0, 1)] + list(F[1:]))):
        with tempfile.TemporaryDirectory() as td:
            rc, txt, c = run_case(key, make_dir(td, soups_with("union", homo(P), FF)), td, log)
            chk("%s: 戻り値" % key, rc, 1, log)
            chk("%s: 理由が union/e0" % key, final(txt), "判定: ★ 失敗（理由=union/e0）", log)
            chk("%s: 後段の呼出し 0" % key,
                (c["degenerate"], c["A1"], c["A2"], c["table"]), (0, 0, 0, 0), log)
            chk("%s: 後続の出力は未着手" % key, "**未着手**（前の出力で停止）" in txt, True, log)


def g03_w0(log):
    P, F = two_boxes()
    V4 = homo(P)
    V4[0] = (V4[0][0], V4[0][1], V4[0][2], 0)
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("03-同次分母 0", make_dir(td, soups_with("union", V4, list(F))),
                              td, log)
        chk("03: 戻り値", rc, 1, log)
        chk("03: 理由が union/e0", final(txt), "判定: ★ 失敗（理由=union/e0）", log)
        chk("03: degenerate 未呼出し", c["degenerate"], 0, log)
    # **関数の対照**: `degenerate` は `None` を返す。**成功と読みません。**
    RAN.append("03b-degenerate が None")
    chk("03b: W=0 の面で degenerate が None", io_.degenerate(V4, F[0]), None, log)


def g04_collinear(log):
    P = [(0, 0, 0), (1, 0, 0), (2, 0, 0)]
    F = [(0, 1, 2)]
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("04-共線", make_dir(td, soups_with("union", homo(P), F)), td, log)
        chk("04: 戻り値", rc, 1, log)
        chk("04: 理由が union/degenerate", final(txt),
            "判定: ★ 失敗（理由=union/degenerate）", log)
        chk("04: ③ は未評価（A1/A2/表を呼ばない）", (c["A1"], c["A2"], c["table"]), (0, 0, 0), log)
        chk("04: 表示が ② の一部", "**② の一部**" in txt, True, log)


def g05_scale(log):
    P, F = two_boxes()
    Pc = [(0, 0, 0), (1, 0, 0), (2, 0, 0)]
    for key, base, w, want_rc, want_final in (
            ("05a-正常 ×7", (P, list(F)), 7, 0, None),
            ("05b-正常 ×(-3)", (P, list(F)), -3, 0, None),
            ("05c-共線 ×5", (Pc, [(0, 1, 2)]), 5, 1,
             "判定: ★ 失敗（理由=union/degenerate）")):
        with tempfile.TemporaryDirectory() as td:
            soups = soups_with("union", homo(base[0], w), base[1])
            rc, txt, c = run_case(key, make_dir(td, soups), td, log)
            chk("%s: 戻り値（同次表現に依らない）" % key, rc, want_rc, log)
            if want_final:
                chk("%s: 理由" % key, final(txt), want_final, log)


def g06_boundary(log):
    P, F = two_boxes()
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("06-面 1 枚欠落",
                              make_dir(td, soups_with("union", homo(P), list(F[:-1]))), td, log)
        chk("06: 戻り値", rc, 1, log)
        chk("06: 理由が union/A1", final(txt), "判定: ★ 失敗（理由=union/A1）", log)
        chk("06: A2・表は未評価", (c["A2"], c["table"]), (0, 0), log)
        chk("06: **E1 も破れる**（A1 固有の対照ではない）",
            io_.e1_violations(list(F[:-1])) != [], True, log)


def g07_degree4(log):
    P, F = ig.join(ig.box((0, 0, 0), (1, 1, 1), +1),
                   ig.box((1, 1, 0), (2, 2, 1), +1))
    n = 8
    pairs = [(n + find(P[n:], (1, 1, 0)), find(P[:n], (1, 1, 0))),
             (n + find(P[n:], (1, 1, 1)), find(P[:n], (1, 1, 1)))]
    P2, F2 = merge_index(P, F, pairs)
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("07-次数 4・収支 0",
                              make_dir(td, soups_with("union", homo(P2), F2)), td, log)
        chk("07: 戻り値", rc, 1, log)
        chk("07: 理由が union/A1", final(txt), "判定: ★ 失敗（理由=union/A1）", log)
        chk("07: **E1 は通る**（収支 0）", io_.e1_violations(F2), [], log)
        chk("07: A3 は未評価（表を呼ばない）", c["table"], 0, log)


def g08_orient(log):
    P, F = two_boxes()
    F2 = [(F[0][0], F[0][2], F[0][1])] + list(F[1:])
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("08-面 1 枚の向き反転",
                              make_dir(td, soups_with("union", homo(P), F2)), td, log)
        chk("08: 戻り値", rc, 1, log)
        chk("08: 理由が union/A2", final(txt), "判定: ★ 失敗（理由=union/A2）", log)
        chk("08: A1 は通る", ig.check_a1(F2), [], log)
        chk("08: A3 は未評価（表を呼ばない）", c["table"], 0, log)


def g09_two_circles(log):
    P, F = ig.join(ig.box((0, 0, 0), (1, 1, 1), +1),
                   ig.box((1, 1, 1), (2, 2, 2), +1))
    n = 8
    pairs = [(n + find(P[n:], (1, 1, 1)), find(P[:n], (1, 1, 1)))]
    P2, F2 = merge_index(P, F, pairs)
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("09-二つのリンク円",
                              make_dir(td, soups_with("union", homo(P2), F2)), td, log)
        chk("09: 戻り値", rc, 1, log)
        chk("09: 理由が union/link", final(txt), "判定: ★ 失敗（理由=union/link）", log)
        chk("09: A1・A2 は通る", (ig.check_a1(F2), ig.check_a2(F2)), ([], []), log)
        chk("09: 旧 check_a3 も 2 円を報告",
            [t for t in ig.check_a3(F2) if t[1] == 2] != [], True, log)
        chk("09: 別索引のままなら 1 円（正常例と区別）", ig.check_a3(F), [], log)


def g10_separate(log):
    P = [(0, 0, 0), (10, 0, 0), (0, 10, 0)]
    F = [(0, 1, 2), (0, 2, 1)]
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("10-正逆 2 面（③ 成立・② 不成立）",
                              make_dir(td, dict((op, (homo(P), list(F))) for op in i3.OPS)),
                              td, log)
        chk("10: 戻り値（③ は成立）", rc, 0, log)
        chk("10: A4（重複面）は ③ に混ぜない", ig.check_a4(F) != [], True, log)
        chk("10: 多重リンク辺を保持（旧方式も 1 円）", ig.check_a3(F), [], log)


def g11_empty(log):
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("11-対象なし（空の面集合）",
                              make_dir(td, dict((op, ([], [])) for op in i3.OPS)), td, log)
        chk("11: 戻り値（全称条件は成立）", rc, 0, log)
        chk("11: 未評価ではない（link 通過）", "link=通過" in txt or "すべて通過" in txt, True, log)


def g12_table(log):
    P, F = two_boxes()

    def drop(fn):
        def g(FF):
            L = fn(FF)
            for v in sorted(L):
                if L[v]:
                    L[v] = L[v][1:]
                    break
            return L
        return g

    def dupl(fn):
        def g(FF):
            L = fn(FF)
            for v in sorted(L):
                if L[v]:
                    L[v] = L[v] + [L[v][0]]
                    break
            return L
        return g

    def misplace(fn):
        def g(FF):
            L = fn(FF)
            ks = sorted(L)
            a, b = ks[0], ks[1]
            r = L[a][0]
            L[a] = L[a][1:]
            L[b] = L[b] + [r]              # **総数は保つが所属が違う**
            return L
        return g

    for key, pf in (("12a-登録欠落", drop), ("12b-重複登録", dupl), ("12c-誤所属", misplace)):
        with tempfile.TemporaryDirectory() as td:
            rc, txt, c = run_case(key, make_dir(td, normal_soups()), td, log,
                                  patches={"table": pf})
            chk("%s: 戻り値" % key, rc, 1, log)
            chk("%s: 理由が union/link(不整合)" % key, final(txt),
                "判定: ★ 失敗（理由=union/link(不整合)/E:不整合）", log)
            chk("%s: **A3 の反例へ誤帰属しない**" % key,
                "**③ の反例にしません**" in txt, True, log)


def g13_entry(log):
    P, F = two_boxes()
    cases = [("13a-期待ハッシュ不一致", dict(sha_override=("union", "0" * 64)), "union/I"),
             ("13b-meta 欠落", dict(meta_missing=True), "I:meta"),
             ("13c-不正なリム設定", dict(meta_text="kHomoXyz=4\nkHomoW=3\n"), "I:meta"),
             ("13d-長さ不一致", dict(truncate=("union", 20)), "union/I")]
    for key, kw, want in cases:
        with tempfile.TemporaryDirectory() as td:
            man = make_dir(td, normal_soups(), **kw)
            rc, txt, c = run_case(key, man, td, log)
            chk("%s: 戻り値" % key, rc, 1, log)
            chk("%s: 理由" % key, final(txt), "判定: ★ 失敗（理由=%s）" % want, log)
            chk("%s: 幾何の後段を呼ばない" % key,
                (c["degenerate"], c["A1"], c["A2"], c["table"]), (0, 0, 0, 0), log)


def g14_save(log):
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("14a-保存先が既存", make_dir(td, normal_soups()), td, log,
                              out_exists=True)
        chk("14a: 戻り値", rc, 1, log)
        chk("14a: 理由は S だけ", final(txt), "判定: ★ 失敗（理由=S）", log)
        chk("14a: 幾何は全部通っている", "/ D∧T 通過" in txt, True, log)

    def short_read(orig):
        def g(path, *a, **k):
            f = orig(path, *a, **k)
            if isinstance(path, str) and path.endswith("index3_result.txt") \
                    and (not a or "r" in str(a[0])) and "w" not in str(a[0] if a else "r"):
                s = f.read()
                f.close()
                return io.StringIO("\n".join(s.splitlines()[:-1]) + "\n")
            return f
        return g

    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("14b-読戻しの内容欠落", make_dir(td, normal_soups()), td, log,
                              patches={"open": short_read})
        chk("14b: 戻り値", rc, 1, log)
        chk("14b: 理由は S だけ", final(txt), "判定: ★ 失敗（理由=S）", log)
        chk("14b: 幾何欠陥として表示しない", "/ D∧T 通過" in txt, True, log)


def g15_four(log):
    P, F = two_boxes()
    bad = [(0, 1, 999)] + list(F[1:])
    for i, op in enumerate(i3.OPS):
        key = "15%d-%s で失敗" % (i, op)
        with tempfile.TemporaryDirectory() as td:
            rc, txt, c = run_case(key, make_dir(td, soups_with(op, homo(P), bad)), td, log)
            chk("%s: 戻り値" % key, rc, 1, log)
            chk("%s: 失敗名と段" % key, final(txt), "判定: ★ 失敗（理由=%s/e0）" % op, log)
            chk("%s: 先行する出力は link まで到達" % key, c["table"], i, log)
            chk("%s: 後続は未着手" % key, txt.count("**未着手**（前の出力で停止）"),
                len(i3.OPS) - 1 - i, log)
    with tempfile.TemporaryDirectory() as td:
        rc, txt, c = run_case("15z-全 4 件正常", make_dir(td, normal_soups()), td, log)
        chk("15z: 戻り値", rc, 0, log)
        chk("15z: 4 出力すべて link まで到達", c["table"], 4, log)


GROUPS = [g01_normal, g02_index, g03_w0, g04_collinear, g05_scale, g06_boundary,
          g07_degree4, g08_orient, g09_two_circles, g10_separate, g11_empty,
          g12_table, g13_entry, g14_save, g15_four]

# **静的なケース表**（予定）。実行した鍵集合と照合します。
PLANNED = [
    "01a-離れた2箱", "01b-点接触だが別索引",
    "02a-範囲外索引", "02b-同一索引の反復",
    "03-同次分母 0", "03b-degenerate が None",
    "04-共線",
    "05a-正常 ×7", "05b-正常 ×(-3)", "05c-共線 ×5",
    "06-面 1 枚欠落",
    "07-次数 4・収支 0",
    "08-面 1 枚の向き反転",
    "09-二つのリンク円",
    "10-正逆 2 面（③ 成立・② 不成立）",
    "11-対象なし（空の面集合）",
    "12a-登録欠落", "12b-重複登録", "12c-誤所属",
    "13a-期待ハッシュ不一致", "13b-meta 欠落", "13c-不正なリム設定", "13d-長さ不一致",
    "14a-保存先が既存", "14b-読戻しの内容欠落",
    "150-union で失敗", "151-isect で失敗", "152-diff_ab で失敗", "153-diff_ba で失敗",
    "15z-全 4 件正常",
]


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    f = open(os.path.join(out, "stage_index3_controls.txt"), "w")

    def log(s):
        print(s)
        sys.stdout.flush()
        f.write(s + "\n")
        f.flush()

    log("=== ③ 限定検査の合成対照（REVIEW §26.4 の 15 群）===")
    log("予定ケース %d 件 / 対照群 %d 件" % (len(PLANNED), len(GROUPS)))
    try:
        for g in GROUPS:
            log("--- %s ---" % g.__name__)
            g(log)
        ng = RES.count(False)
        log("実行ケース %d 件 / 予定 %d 件" % (len(RAN), len(PLANNED)))
        miss = sorted(set(PLANNED) - set(RAN))
        extra = sorted(set(RAN) - set(PLANNED))
        if miss or extra:
            RES.append(False)
            ng += 1
            log("  ★NG 鍵集合が違います。未実行 %s / 予定外 %s" % (miss, extra))
        else:
            log("  ok   予定ケースと実行ケースの鍵集合が一致")
        if len(RAN) != len(set(RAN)):
            RES.append(False)
            ng += 1
            log("  ★NG 同じ鍵が 2 度実行されています")
        log("判定の件数 %d / 通過 %d / 失敗 %d" % (len(RES), RES.count(True), ng))
        log("判定: %s" % ("合成対照は全件通過" if ng == 0 else "★ 失敗あり"))
        sys.exit(0 if ng == 0 else 1)
    finally:
        f.close()

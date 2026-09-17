#!/usr/bin/env python3
r"""③（索引接続）の限定検査 — `REVIEW-phase5-correctness-literature.md` §26。

**この子から別のプロセスを起動しません。** 実データは**固定の 4 出力だけ**を選べます。
**入力 A/B は読みません**（③ は索引と保存頂点だけで閉じます）。

成功条件（§26.1）:

    rc = 0  ⟺  I ∧ E ∧ S ∧ ⋀_o ( D_o ∧ T_o )

- `I`  対象 4 ファイルと meta の同一性・書式・リム設定・長さ
- `D_o` 参照範囲 / 3 索引相異 / 同次分母非零 / **全使用面の非共線**（**非共線は ② の一部**）
- `T_o` **③**: 全使用索引辺の次数 2・向き収支 0、全使用頂点のリンク 1 円
- `E`  4 出力の全必要判定を完了し、**接続面表の不変条件**も成立
- `S`  結果表・反例・識別情報の保存と、その読戻し・チェックサム記録に成功
"""
import argparse, hashlib, os, struct, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inspect_geom as ig
import inspect_output as io_

OPS = ("union", "isect", "diff_ab", "diff_ba")

# ---- 基準マニフェスト（承認済み束から**転記**。検査対象から作りません。§26.2）----
REAL_ROOT = "data/evidence/output_20260916"
REAL_SUMS_PATH = REAL_ROOT + "/SHA256SUMS"
REAL_SUMS_SHA = "3c6e163794f3f5aa520305a758ffa29815a2202934724aa7b05165935755be70"
REAL_META = REAL_ROOT + "/run/b_run_meta.txt"
REAL_META_SHA = "96dc84ad71d45848e879b9cfd918a98b0465ae0659a99406262120ab6ada83a3"
REAL_OUT = {
    "union":   (REAL_ROOT + "/run/b_out_union.bin",
                "a06d320901bfffcbb93653ca4f0eb21708e9aa5b624adcb9f3395193534fe82a"),
    "isect":   (REAL_ROOT + "/run/b_out_isect.bin",
                "44a7bc3a2523dc8f988b0ace121e2aa06f29a149706379937a7628f0457ca7a6"),
    "diff_ab": (REAL_ROOT + "/run/b_out_diff_ab.bin",
                "3fbe2bd704669243186934cf1235551e00d71f8aa812e3bb6ea7a4ccf6674b1b"),
    "diff_ba": (REAL_ROOT + "/run/b_out_diff_ba.bin",
                "1876020a4e715c6c9aa6f5d318523398be0c9ac480c6a762c92ddeb13abbab4c"),
}
REAL_LIMBS = (3, 3)          # kHomoXyz / kHomoW。**meta を解析して一致を確かめます**

# 段の名前（固定順）。**未着手と未評価を分けます。**
STAGES = ("I", "e0", "degenerate", "A1", "A2", "link")


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def read_exact(path):
    with open(path, "rb") as f:
        b = f.read()
    return b


def parse_meta(text):
    d = {}
    for line in text.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


# ---- 接続面表（§26.3）--------------------------------------------------
def build_link_table(F):
    r"""元面の各隅を 1 回ずつ登録します。**重複を消しません。**

    $L[a]\ni(f,0,b,c)$, $L[b]\ni(f,1,c,a)$, $L[c]\ni(f,2,a,b)$。
    """
    L = {}
    for fi, f in enumerate(F):
        for pos in range(3):
            v = f[pos]
            x = f[(pos + 1) % 3]
            y = f[(pos + 2) % 3]
            L.setdefault(v, []).append((fi, pos, x, y))
    return L


def table_invariants(F, L):
    """表の不変条件。戻り値: `(成立か, 理由, 計数の辞書)`。**③ の反例にしません。**"""
    used = set()
    for f in F:
        used.update(f)
    bad = []
    seen = set()
    n_rec = 0
    for v, recs in L.items():
        for (fi, pos, x, y) in recs:
            n_rec += 1
            if not (0 <= fi < len(F)) or pos not in (0, 1, 2):
                bad.append(("範囲外", v, fi, pos))
                continue
            f = F[fi]
            if f[pos] != v or f[(pos + 1) % 3] != x or f[(pos + 2) % 3] != y:
                bad.append(("元面と不一致", v, fi, pos))
            if (fi, pos) in seen:
                bad.append(("重複", v, fi, pos))
            seen.add((fi, pos))
    cnt = {"使用頂点": len(used), "表の鍵": len(L), "登録": n_rec, "期待登録": 3 * len(F)}
    if set(L.keys()) != used:
        bad.append(("鍵集合が使用頂点集合と違う", len(set(L.keys()) ^ used), None, None))
    if len(seen) != 3 * len(F):
        bad.append(("隅の被覆が全単射でない", len(seen), 3 * len(F), None))
    return (not bad), bad, cnt


def link_state(recs):
    """1 頂点のリンク。戻り値: `(単一閉路か, 次数の違反, 成分数)`。"""
    deg = {}
    adj = {}
    for (_fi, _pos, x, y) in recs:
        deg[x] = deg.get(x, 0) + 1
        deg[y] = deg.get(y, 0) + 1
        adj.setdefault(x, []).append(y)
        adj.setdefault(y, []).append(x)
    badd = sorted((k, v) for k, v in deg.items() if v != 2)
    if badd:
        return False, badd, None
    seen, n = set(), 0
    for s in adj:
        if s in seen:
            continue
        n += 1
        st = [s]
        seen.add(s)
        while st:
            u = st.pop()
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    st.append(w)
    return (n == 1), [], n


# ---- 1 出力の評価 ------------------------------------------------------
def check_one(name, path, want_sha, nx, nw, log):
    r"""戻り値: `(段ごとの状態, 反例, 計数)`。

    状態は **`通過` / `違反` / `未評価` / `不整合`** の 4 つ。
    **初めて不成立・内部不整合を得た段で後続を止め、後続は `未評価`** にします。
    """
    st = dict((s, "未着手") for s in STAGES)
    ex = {}
    cnt = {}

    # --- I ---
    try:
        raw = read_exact(path)
    except Exception as e:
        st["I"] = "違反"
        ex["I"] = "読取り例外: %r" % (e,)
        return st, ex, cnt
    got = sha256_bytes(raw)
    cnt["バイト長"] = len(raw)
    if got != want_sha:
        st["I"] = "違反"
        ex["I"] = "ハッシュ不一致: 得た %s / 期待 %s" % (got, want_sha)
        return st, ex, cnt
    try:
        V, F, sha2 = io_.read_soup(path, nx, nw)
    except Exception as e:
        st["I"] = "違反"
        ex["I"] = "read_soup が失敗: %r" % (e,)
        return st, ex, cnt
    if sha2 != want_sha:
        st["I"] = "不整合"
        ex["I"] = "read_soup の読取りバイト列が期待と違います: %s" % sha2
        return st, ex, cnt
    st["I"] = "通過"
    cnt["頂点(保存)"] = len(V)
    cnt["面"] = len(F)
    log("[%s] I 通過 / 保存頂点 %d / 面 %d / sha %s" % (name, len(V), len(F), got[:16]))

    # --- e0（参照範囲・3 索引相異・同次分母）---
    oor, dup, w0 = io_.e0_violations(V, F)
    cnt["e0 判定"] = 2 * len(F) + len(V)
    if oor or dup or w0:
        st["e0"] = "違反"
        ex["e0"] = {"参照範囲": oor[:1], "3索引相異": dup[:1], "分母0": w0[:1],
                    "件数": (len(oor), len(dup), len(w0))}
        for s in ("degenerate", "A1", "A2", "link"):
            st[s] = "未評価"
        log("[%s] ★ e0 違反（範囲 %d / 相異 %d / 分母 %d）。後続は未評価"
            % (name, len(oor), len(dup), len(w0)))
        return st, ex, cnt
    st["e0"] = "通過"

    # --- 非共線（**② の一部。③ には混ぜません**）---
    used_f = list(range(len(F)))
    deg_bad, deg_none = [], []
    for i in used_f:
        r = io_.degenerate(V, F[i])
        if r is None:
            deg_none.append(i)
        elif r:
            deg_bad.append(i)
    cnt["degenerate 判定"] = len(used_f)
    if deg_none:
        st["degenerate"] = "不整合"
        ex["degenerate"] = "e0 を通ったのに None: 面 %s（%d 件）" % (deg_none[:1], len(deg_none))
        for s in ("A1", "A2", "link"):
            st[s] = "未評価"
        log("[%s] ★ degenerate が None（前段との内部不整合）。後続は未評価" % name)
        return st, ex, cnt
    if deg_bad:
        st["degenerate"] = "違反"
        ex["degenerate"] = "共線の面 %s（%d 件。**② の一部**）" % (deg_bad[:1], len(deg_bad))
        for s in ("A1", "A2", "link"):
            st[s] = "未評価"
        log("[%s] ★ 非共線の補助条件が破れました（%d 枚。**② の一部**）。③ は未評価"
            % (name, len(deg_bad)))
        return st, ex, cnt
    st["degenerate"] = "通過"

    # --- A1（無向辺の次数 2）---
    a1 = ig.check_a1(F)
    cnt["A1 判定"] = 3 * len(F)
    if a1:
        st["A1"] = "違反"
        ex["A1"] = "辺 %s（%d 件）" % (a1[:1], len(a1))
        for s in ("A2", "link"):
            st[s] = "未評価"
        log("[%s] ★ A1 違反 %d 件。A2・link は未評価" % (name, len(a1)))
        return st, ex, cnt
    st["A1"] = "通過"

    # --- A2（有向辺がちょうど 1 回 = 向き収支 0）---
    a2 = ig.check_a2(F)
    cnt["A2 判定"] = 3 * len(F)
    if a2:
        st["A2"] = "違反"
        ex["A2"] = "有向辺 %s（%d 件）" % (a2[:1], len(a2))
        st["link"] = "未評価"
        log("[%s] ★ A2 違反 %d 件。link は未評価" % (name, len(a2)))
        return st, ex, cnt
    st["A2"] = "通過"

    # --- 接続面表 → リンク ---
    L = build_link_table(F)
    ok_inv, bad_inv, c2 = table_invariants(F, L)
    cnt.update(c2)
    if not ok_inv:
        st["link"] = "不整合"
        ex["link"] = "表の不変条件: %s（%d 件）" % (bad_inv[:1], len(bad_inv))
        log("[%s] ★ 接続面表の不変条件が破れました（%d 件）。**③ の反例にしません**"
            % (name, len(bad_inv)))
        return st, ex, cnt
    bad_v = []
    for v in sorted(L):
        ok, badd, ncomp = link_state(L[v])
        if not ok:
            bad_v.append((v, badd[:1], ncomp))
    cnt["link 判定(頂点)"] = len(L)
    if bad_v:
        st["link"] = "違反"
        ex["link"] = "頂点 %s（%d 件）" % (bad_v[:1], len(bad_v))
        log("[%s] ★ A3 違反 %d 頂点" % (name, len(bad_v)))
        return st, ex, cnt
    st["link"] = "通過"
    log("[%s] e0・非共線・A1・A2・link すべて通過（使用頂点 %d / 登録 %d）"
        % (name, c2["使用頂点"], c2["登録"]))
    return st, ex, cnt


# ---- 段の本体 ----------------------------------------------------------
def run_manifest(man, outdir, log, deps):
    r"""`man` は `{"meta": (path, sha), "limbs": (nx, nw), "out": {op: (path, sha)}}`。

    戻り値は終了値。**保存・読戻しまで終えてから成功を出します。**
    """
    t_all = time.process_time()
    log("=== ③ 限定検査（REVIEW §26）===")
    log("対象 %d 出力（固定順 %s）" % (len(OPS), "/".join(OPS)))
    for k, v in sorted(deps.items()):
        log("依存の版 %s = %s" % (k, v))

    # --- I: meta ---
    I_ok, I_why = True, ""
    nx = nw = None
    try:
        mraw = read_exact(man["meta"][0])
        msha = sha256_bytes(mraw)
        if msha != man["meta"][1]:
            I_ok, I_why = False, "meta のハッシュ不一致: %s" % msha
        else:
            md = parse_meta(mraw.decode("utf-8", "replace"))
            if "kHomoXyz" not in md or "kHomoW" not in md:
                I_ok, I_why = False, "meta に kHomoXyz/kHomoW がありません"
            else:
                try:
                    nx, nw = int(md["kHomoXyz"]), int(md["kHomoW"])
                except ValueError:
                    I_ok, I_why = False, "リム数が整数ではありません"
                if I_ok and not (nx > 0 and nw > 0):
                    I_ok, I_why = False, "リム数が正ではありません: %r/%r" % (nx, nw)
                if I_ok and (nx, nw) != man["limbs"]:
                    I_ok, I_why = False, "リム設定が固定対象と違います: %d/%d" % (nx, nw)
    except Exception as e:
        I_ok, I_why = False, "meta の読取り例外: %r" % (e,)
    log("I(meta): %s%s" % ("通過" if I_ok else "★ 不成立", "" if I_ok else " — " + I_why))

    res = {}
    reasons = []
    if not I_ok:
        reasons.append("I:meta")
        for op in OPS:
            res[op] = (dict((s, "未着手") for s in STAGES), {"I": I_why}, {})
    else:
        stop = False
        for op in OPS:
            if stop:
                res[op] = (dict((s, "未着手") for s in STAGES), {}, {})
                log("[%s] **未着手**（前の出力で停止）" % op)
                continue
            p, h = man["out"][op]
            st, ex, cnt = check_one(op, p, h, nx, nw, log)
            res[op] = (st, ex, cnt)
            bad = [s for s in STAGES if st[s] in ("違反", "不整合")]
            if bad:
                s0 = bad[0]
                reasons.append("%s/%s%s" % (op, s0, "(不整合)" if st[s0] == "不整合" else ""))
                stop = True

    # --- E: 全必要判定の完了 ---
    E_ok = I_ok and all(res[op][0]["link"] == "通過" for op in OPS)
    # **不整合は E を偽にします**（違反 0 に化けさせないため）
    E_bad = [op for op in OPS if any(res[op][0][s] == "不整合" for s in STAGES)]
    if E_bad:
        E_ok = False
        if "E:不整合" not in reasons:
            reasons.append("E:不整合")
    D_T = all(res[op][0]["link"] == "通過" for op in OPS)
    log("E(全評価の完了): %s" % ("通過" if E_ok else "★ 不成立"))

    # --- S: 保存と読戻し ---
    S_ok, S_why = True, ""
    try:
        os.makedirs(outdir, exist_ok=False)
    except Exception as e:
        S_ok, S_why = False, "保存先を作れません（既存か作成不能）: %r" % (e,)
    lines = []
    if S_ok:
        lines.append("# ③ 限定検査の結果（保存前の判定内訳）")
        lines.append("manifest_meta\t%s\t%s" % (man["meta"][0], man["meta"][1]))
        for k in sorted(deps):
            lines.append("dep\t%s\t%s" % (k, deps[k]))
        lines.append("I\t%s\t%s" % ("通過" if I_ok else "不成立", I_why))
        for op in OPS:
            st, ex, cnt = res[op]
            lines.append("out\t%s\t%s\t%s" % (op, man["out"][op][0] if I_ok else "-",
                                              man["out"][op][1]))
            lines.append("state\t%s\t%s" % (op, "\t".join("%s=%s" % (s, st[s]) for s in STAGES)))
            lines.append("count\t%s\t%s" % (op, "\t".join("%s=%s" % (k, cnt[k])
                                                          for k in sorted(cnt))))
            for k in sorted(ex):
                lines.append("evidence\t%s\t%s\t%s" % (op, k, ex[k]))
        lines.append("E\t%s" % ("通過" if E_ok else "不成立"))
        rpath = os.path.join(outdir, "index3_result.txt")
        try:
            with open(rpath, "w") as f:
                f.write("\n".join(lines) + "\n")
            with open(rpath) as f:
                back = f.read().splitlines()
            if back != lines:
                S_ok, S_why = False, "読戻しが書いた内容と違います（行数 %d / %d）" % (
                    len(back), len(lines))
        except Exception as e:
            S_ok, S_why = False, "結果の保存・読戻しに失敗: %r" % (e,)
    if S_ok:
        try:
            spath = os.path.join(outdir, "SHA256SUMS")
            ent = ["%s  %s" % (sha256_bytes(read_exact(rpath)), "index3_result.txt")]
            with open(spath, "w") as f:
                f.write("\n".join(ent) + "\n")
            with open(spath) as f:
                if f.read().splitlines() != ent:
                    S_ok, S_why = False, "チェックサム一覧の読戻しが違います"
        except Exception as e:
            S_ok, S_why = False, "チェックサム一覧の保存に失敗: %r" % (e,)
    if not S_ok:
        sys.stderr.write("保存不成立: %s\n" % S_why)
        reasons.append("S")
    log("S(保存と読戻し): %s%s" % ("通過" if S_ok else "★ 不成立",
                                   "" if S_ok else " — " + S_why))

    rc = 0 if (I_ok and E_ok and S_ok and D_T) else 1
    log("内訳: I %s / E %s / S %s / D∧T %s"
        % (*(("通過" if x else "★ 不成立") for x in (I_ok, E_ok, S_ok, D_T)),))
    log("CPU %.3f 秒" % (time.process_time() - t_all))
    if rc == 0:
        log("判定: ③ 限定検査 成功（**固定保存物の ③ と補助条件の成立であって、"
            "② 全体・①・④ の成立ではありません**）")
    else:
        log("判定: ★ 失敗（理由=%s）" % "/".join(reasons))
    return rc


def real_manifest():
    return {"meta": (REAL_META, REAL_META_SHA), "limbs": REAL_LIMBS,
            "out": dict(REAL_OUT)}


def dep_versions(extra=()):
    d = {}
    here = os.path.dirname(os.path.abspath(__file__))
    for n in ("inspect_geom.py", "inspect_output.py", "inspect_index3.py") + tuple(extra):
        p = os.path.join(here, n)
        try:
            d[n] = sha256_bytes(read_exact(p))[:16]
        except Exception:
            d[n] = "**読めません**"
    return d


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--real", action="store_true",
                   help="固定の保存 4 出力を検査します（**これ以外の実データは選べません**）")
    p.add_argument("--out", required=True, help="**新規**の保存先")
    a = p.parse_args(argv)
    if not a.real:
        sys.stderr.write("--real が要ります（合成は対照プログラムから内部関数を呼びます）\n")
        return 1

    def log(s):
        print(s)
        sys.stdout.flush()
    return run_manifest(real_manifest(), a.out, log, dep_versions())


if __name__ == "__main__":
    sys.exit(main())

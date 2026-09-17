#!/usr/bin/env python3
r"""記録基盤 2 本の限定検証 — コンパイルと合成対照（`SPEC-phase5.md` §5.10.14.149 §4）。

**この子から別のプロセスを起動しません**（コンパイラの起動だけは行います）。
**本体へは接続しません。模型の Boolean も回しません。**

段:
  1. **診断定義 ON** で `test_cause_trace` をコンパイル・実行
  2. **診断定義 OFF** で `test_cause_trace_off` をコンパイル・実行（**依存漏れの検査**）
  3. 書かれた JSON Lines の **往復**（制御文字を含む）と、**末尾 `dropped` 行**の確認
"""
import json, os, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
CXX = os.environ.get("KRI_CXX", "zigcxx")
RES = []


def chk(name, got, want):
    ok = (got == want)
    RES.append(ok)
    print("  %s %-58s 得た値 %-28s 期待 %s" % ("ok  " if ok else "**NG**", name, got, want))
    sys.stdout.flush()


def compile_one(src, out, defs):
    cmd = [CXX, "-std=c++20", "-O0", "-g0",
           "-I", os.path.join(ROOT, "include"), "-I", os.path.join(ROOT, "tests"),
           os.path.join(ROOT, src), "-o", out] + ["-D" + d for d in defs]
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr), cmd


def main(outdir):
    print("=== 記録基盤 2 本の限定検証 ===")
    print("コンパイラ %s / 出力先 %s" % (CXX, outdir))
    os.makedirs(outdir, exist_ok=False)     # **新規先だけ**

    # --- 段 1: 診断定義 ON ---
    print("--- 段 1: 診断定義 ON のコンパイルと実行 ---")
    exe_on = os.path.join(outdir, "test_cause_trace")
    rc, log, cmd = compile_one("tests/csg/test_cause_trace.cpp", exe_on, ["KRISITE_DIAG_CAUSE_TRACE"])
    with open(os.path.join(outdir, "compile_on.log"), "w") as f:
        f.write(" ".join(cmd) + "\n\n" + log)
    chk("段 1: ON のコンパイルが通る", rc, 0)
    if rc != 0:
        print(log[-2000:])
        return 1
    run_dir = os.path.join(outdir, "run_on")
    os.makedirs(run_dir, exist_ok=False)
    p = subprocess.run([exe_on, run_dir], capture_output=True, text=True)
    with open(os.path.join(outdir, "run_on.log"), "w") as f:
        f.write(p.stdout + p.stderr)
    for line in p.stdout.splitlines():
        print("    | " + line)
    chk("段 1: ON の実行が全件通過", p.returncode, 0)

    # --- 段 2: 診断定義 OFF（依存漏れの検査）---
    print("--- 段 2: 診断定義 OFF の単独 include ---")
    exe_off = os.path.join(outdir, "test_cause_trace_off")
    rc2, log2, cmd2 = compile_one("tests/csg/test_cause_trace_off.cpp", exe_off, [])
    with open(os.path.join(outdir, "compile_off.log"), "w") as f:
        f.write(" ".join(cmd2) + "\n\n" + log2)
    chk("段 2: OFF のコンパイルが通る", rc2, 0)
    if rc2 != 0:
        print(log2[-2000:])
        return 1
    p2 = subprocess.run([exe_off], capture_output=True, text=True)
    with open(os.path.join(outdir, "run_off.log"), "w") as f:
        f.write(p2.stdout + p2.stderr)
    for line in p2.stdout.splitlines():
        print("    | " + line)
    chk("段 2: OFF の実行が通過", p2.returncode, 0)
    # **OFF のバイナリに記録の文字列が残っていないこと**（診断が漏れていない直接の検査）
    with open(exe_off, "rb") as f:
        blob = f.read()
    chk("段 2: OFF のバイナリに記録の語が無い",
        (b"winding_split" in blob) or (b"vertex_split" in blob), False)

    # --- 段 3: JSON の往復 ---
    print("--- 段 3: JSON Lines の往復 ---")
    jp = os.path.join(run_dir, "trace_ok.jsonl")
    chk("段 3: 記録ファイルがある", os.path.exists(jp), True)
    with open(jp, encoding="utf-8") as f:
        raw = f.read().splitlines()
    recs = []
    bad = None
    for i, line in enumerate(raw):
        try:
            recs.append(json.loads(line))
        except Exception as e:
            bad = (i, repr(e))
            break
    chk("段 3: 全行が JSON として読める", bad, None)
    if bad:
        return 1
    kinds = [r.get("rec") for r in recs]
    chk("段 3: 先頭が header", kinds[0] if kinds else None, "header")
    chk("段 3: **末尾が dropped**（末尾欠落は未完了）", kinds[-1] if kinds else None, "dropped")
    chk("段 3: 落とした件数 0", [r for r in recs if r["rec"] == "dropped"][0]["count"], 0)

    rm = [r for r in recs if r["rec"] == "removed_face"]
    chk("段 3: 制御文字が可逆に往復する（タブ・CR）",
        rm[0]["where"] if rm else None, "loader.hpp:quantize\tcol\rreturn")
    chk("段 3: 量子化の同一索引を共線除去と区別",
        rm[0]["reason"] if rm else None, "量子化で同一索引")

    reg = [r for r in recs if r["rec"] == "region"]
    ps = reg[0]["per_source"] if reg else []
    chk("段 3: source 別の記録が 3 件", len(ps), 3)
    chk("段 3: forced と winding_split が混在",
        sorted({x["path"] for x in ps}), ["forced", "winding_split", "未実施"])
    chk("段 3: 供給セルが Morton 切詰めでない鍵",
        ps[0]["forced_from"] if ps else None, {"depth": 3, "i": 1, "j": 2, "k": 2})
    chk("段 3: 未実施の source が値 0 と区別される",
        [x["evaluated"] for x in ps], [True, True, False])

    vs = [r for r in recs if r["rec"] == "vertex_split"]
    cor = vs[0]["corners"] if vs else []
    chk("段 3: 隅が (元面, 局所隅) で特定できる",
        [(c["face"], c["local_corner"]) for c in cor], [(10, 0), (11, 2), (12, 1)])
    chk("段 3: 旧頂点 → 類 → 新頂点が辿れる",
        [(c["class"], c["new_vertex"]) for c in cor], [(0, 100), (0, 100), (1, 101)])

    ver = [r for r in recs if r["rec"] == "verify"]
    chk("段 3: 検算の前提の根拠が **未確認** のまま残る",
        ver[0]["premise_basis"] if ver else None, "**未確認**")
    chk("段 3: 有理数が文字列で残る",
        ver[0]["result"] if ver else None, "-69366264756244193/1")

    ng = RES.count(False)
    print("判定の件数 %d / 通過 %d / 失敗 %d" % (len(RES), RES.count(True), ng))
    print("判定: %s" % ("限定検証は全件通過" if ng == 0 else "**失敗あり**"))
    return 0 if ng == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))

#!/usr/bin/env python3
r"""G1・G2・G4 の限定試験を、段ごとに監督して回します（`SPEC-phase5.md` §5.10.14.111）。

**J-0 と J-2 を直接、直列に起動します。** **その子から別のプロセスを起動しません。**
**既存の監督（`run_gmp_diag.py` の `spawn` / `supervise`）を再利用し、改修しません。**
**G3・J-1・ビルド・被検体の再実行・案 II・対象拡大は対象外です。**

  * 段の枠 … J-0 120 秒 / J-2 120 秒（**終了猶予込み**）
  * 全体の絶対期限 … `--deadline`（既定 300 秒）
  * `TERM` → `KILL` … `--grace`（既定 10 秒）。**枠の内側**
  * 回収の上限 … `--reap`（既定 30 秒）。**枠の外**
  * 仮想アドレス空間 … `--as-gib`（既定 4 GiB）。**単一ワーカー**

**保存済みの証拠は読むだけです。** `--read-from` の束は**削除の対象にしません**。
**`--out` が保護対象と同じ実体を指したら、削除処理の【前】に止めます。**
"""
import argparse, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_gmp_diag as sup

STAGES = [("J-0", 120), ("J-2", 120)]
ARTIFACTS = {
    "J-0": ["stage_J0_controls.txt"],
    "J-2": ["stage_measure.txt", "stage_measure_c1_intervals.txt"],
}
# **変更しないと宣言した証拠束**（`--out` がこれらと同じか、含む／含まれるなら起動しません）
PROTECTED = [
    "data/evidence/output_20260916",
    "data/evidence/nesting_20260915",
    "data/evidence/rev_oos_20260915",
    "data/logs/inspect/20260915-205108",
    "docs/evidence",
]


def _inside(a, b):
    a, b = os.path.realpath(a), os.path.realpath(b)
    return a == b or a.startswith(b + os.sep) or b.startswith(a + os.sep)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="")
    p.add_argument("--read-from", default="data/evidence/output_20260916/run")
    p.add_argument("--sums", default="data/evidence/output_20260916/SHA256SUMS")
    p.add_argument("--deadline", type=float, default=300.0)
    p.add_argument("--grace", type=float, default=10.0)
    p.add_argument("--reap", type=float, default=30.0)
    p.add_argument("--as-gib", type=float, default=4.0)
    p.add_argument("--in-a", default="data/logs/inspect/20260915-205108/C_A_quantized.bin")
    p.add_argument("--in-b", default="data/logs/inspect/20260915-205108/C_B_quantized.bin")
    p.add_argument("--saved", default="docs/evidence/gmp_diag_r1/cp3_gmp_results.txt")
    p.add_argument("--key", default="250394x45413")
    a = p.parse_args(argv)

    for nm, v in (("--deadline", a.deadline), ("--grace", a.grace),
                  ("--reap", a.reap), ("--as-gib", a.as_gib)):
        if not (v == v) or v in (float("inf"), float("-inf")) or v <= 0:
            p.error("%s は有限の正の数であること: %r" % (nm, v))
    if a.grace >= a.deadline:
        p.error("--grace (%g) は --deadline (%g) より小さいこと" % (a.grace, a.deadline))

    out = a.out or os.path.join("data", "logs", "gapfix", time.strftime("%Y%m%d-%H%M%S"))
    # ★ 保護対象との重なりを、【削除処理の前に】検査します。
    for q in PROTECTED + [a.read_from]:
        if os.path.exists(q) and _inside(out, q):
            print("★ --out (%s) が保護対象 (%s) と重なります。起動しません。" % (out, q))
            return 7
    # ★ 既に存在したら起動しません（秒単位の ID は非衝突の保証ではありません）。
    try:
        os.makedirs(out, exist_ok=False)
    except FileExistsError:
        print("★ 保存先が既に存在します: %s。起動しません。" % out)
        return 7
    as_bytes = int(a.as_gib * 1024**3)
    here = os.path.dirname(os.path.abspath(__file__))

    t0 = time.monotonic()
    d_abs = t0 + a.deadline
    print("これから %d 段を回します: %s" % (len(STAGES), " → ".join(n for n, _ in STAGES)))
    print("  読み先（**読むだけ**） %s / 照合表 %s" % (a.read_from, a.sums))
    print("  入力 %s, %s" % (a.in_a, a.in_b))
    print("  保存物 %s / 鍵 %s" % (a.saved, a.key))
    print("  段の枠 %s / 合計 %g 秒 / 絶対期限 %g 秒 / 猶予 %g 秒（枠の内側）/ "
          "回収 %g 秒（枠の外）/ 単一ワーカー %d バイト"
          % ({n: b for n, b in STAGES}, sum(b for _, b in STAGES), a.deadline,
             a.grace, a.reap, as_bytes))
    inj = {k: os.environ[k] for k in os.environ if k.startswith("KRI_DIAG_TEST_")}
    print("  試験専用の注入: %s" % (inj if inj else "無し"))
    print("  書き先 %s（新規に作りました）" % out)
    print("  **G3・J-1・ビルド・被検体の再実行・案 II・対象拡大は対象外です。**")
    sys.stdout.flush()

    cmds = {
        "J-0": [sys.executable, os.path.join(here, "inspect_gapfix_controls.py"), out],
        "J-2": [sys.executable, os.path.join(here, "inspect_output.py"),
                "--stage", "measure", "--out", out, "--read-from", a.read_from,
                "--sums", a.sums, "--in-a", a.in_a, "--in-b", a.in_b,
                "--saved", a.saved, "--key", a.key],
    }
    state = {n: "未起動" for n, _ in STAGES}

    def report(rc):
        print("実行状態: %s" % state)
        notrun = [n for n, _ in STAGES if state[n] == "未起動"]
        print("未起動の段: %s" % (notrun if notrun else "無し"))
        print("**未測定・未評価・不一致は成功にしません。**")
        print("出力: %s" % out)
        return rc

    for name, budget in STAGES:
        now = time.monotonic()
        if d_abs - now <= a.grace:
            print("[%s] ★ 残り %.1f 秒では終了猶予 %g 秒を確保できません。起動しません。"
                  % (name, d_abs - now, a.grace))
            return report(4)
        t_kill = min(now + budget, d_abs)
        t_term = max(now, t_kill - a.grace)
        log = os.path.join(out, "child_%s.log" % name)
        print("[%s] 起動: %s" % (name, " ".join(cmds[name])))
        print("[%s]   TERM は %.1f 秒後 / KILL は %.1f 秒後 / 回収の上限 %g 秒"
              % (name, t_term - now, t_kill - now, a.reap))
        sys.stdout.flush()
        pid = sup.spawn(cmds[name], dict(os.environ), log, as_bytes)
        code, rss, st = sup.supervise(pid, t_term, t_kill, hard_cap=a.reap)
        dt = time.monotonic() - now
        try:
            sys.stdout.write(open(log, encoding="utf-8", errors="replace").read())
        except OSError:
            pass
        print("[%s] 実行状態 %s / 終了値 %s / ピーク RSS %s / 経過 %.2f 秒"
              % (name, st, "不明" if code is None else code,
                 "不明" if rss is None else "%d KiB" % rss, dt))
        sys.stdout.flush()
        if st == sup.UNREAPED:
            state[name] = "回収不能"
            print("[%s] ★ 回収できていません（pid=%d）。後続を起動しません。" % (name, pid))
            return report(2)
        if st == sup.CUT:
            state[name] = "打ち切り"
            print("[%s] ★ 打ち切り、または期限後の回収。**未完了として区切ります。**" % name)
            return report(3)
        if code == 127:
            state[name] = "起動失敗"
            return report(6)
        if code != 0:
            state[name] = "失敗"
            print("[%s] ★ 子が終了値 %s で終わりました。後続を起動しません。" % (name, code))
            return report(1)
        state[name] = "完了"
        miss = [f for f in ARTIFACTS[name]
                if not os.path.exists(os.path.join(out, f))
                or os.path.getsize(os.path.join(out, f)) == 0]
        if miss:
            print("[%s] ★ 成果物が無いか空です: %s" % (name, ", ".join(miss)))
            return report(5)
        # ★ 一覧は見出しだけを「非空」と数えません。
        ip = os.path.join(out, "stage_measure_c1_intervals.txt")
        if name == "J-2" and sum(1 for L in open(ip) if not L.startswith("#")) < 2:
            print("[J-2] ★ 区間の一覧が見出しだけです。")
            return report(5)

    print("全 %d 段を回しました。経過 %.2f 秒（絶対期限 %g 秒）"
          % (len(STAGES), time.monotonic() - t0, a.deadline))
    return report(0)


if __name__ == "__main__":
    sys.exit(main())

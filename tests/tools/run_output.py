#!/usr/bin/env python3
r"""案 I（`DESIGN-phase5-vertex-level.md` §44.8）を、段ごとに監督して回します。

**既存の監督（`run_gmp_diag.py` の `spawn` / `supervise`）を再利用します。改修しません。**

  * 段の枠 … I-0 90 / I-1 30 / I-2 120 / I-3 30 / I-4 120 秒（**合計 390 秒。終了猶予込み**）
  * 全体の絶対期限 … `--deadline`（既定 480 秒）
  * `TERM` → `KILL` … `--grace`（既定 10 秒）。**枠の内側**
  * 回収の枠 … `--reap`（既定 30 秒）。**枠の外**
  * 仮想アドレス空間 … `--as-gib`（既定 4 GiB）。**単一ワーカー**（同時に走る子は 1 つ）
  * 被検体のスレッド数 … `--threads`（既定 8。**保存された実行と同じ**。単一ワーカーとは別）

**停止条件**（§44.7 (D)。**自動延長・自動再試行・対象拡大なし**）::

    実行状態が 起動失敗 / 打ち切り / 回収不能 / 失敗
    成果物が無いか空
    入力の同一性が破れた（I-1）
    E0 が破れた（I-4 の中で止まります）
    合成対照が落ちた（I-0）
    絶対期限までに終了猶予を確保できない

**★ 錨の不一致（I-3）は停止条件ではありません。** 単価は錨と無関係です。
**案 I が全段通っても、案 II へ自動では進みません。**
"""
import argparse, hashlib, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_gmp_diag as sup

STAGES = [("I-0", 90), ("I-1", 30), ("I-2", 120), ("I-3", 30), ("I-4", 120)]
ARTIFACTS = {
    "I-0": ["stage_controls.txt"],
    "I-1": ["stage_identity.txt"],
    "I-2": ["b_out_union.bin", "b_out_isect.bin", "b_out_diff_ab.bin", "b_out_diff_ba.bin",
            "b_run_meta.txt"],
    "I-3": ["stage_anchor.txt"],
    "I-4": ["stage_measure.txt"],
}
# **錨の不一致で止めない段**（§44.7 (D)）
NO_STOP_ON_NONZERO = set()


def sha256_of(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="")
    p.add_argument("--deadline", type=float, default=480.0)
    p.add_argument("--grace", type=float, default=10.0)
    p.add_argument("--reap", type=float, default=30.0)
    p.add_argument("--as-gib", type=float, default=4.0)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--dumper", default="build/tests/dump_boolean")
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
    if a.threads < 1:
        p.error("--threads は 1 以上であること")

    out = a.out or os.path.join("data", "logs", "output", time.strftime("%Y%m%d-%H%M%S"))
    os.makedirs(out, exist_ok=True)
    # ★ 書く前に、前の実行の成果物と子のログを消します。
    stale = 0
    for name, _ in STAGES:
        for f in ARTIFACTS[name] + ["child_%s.log" % name]:
            q = os.path.join(out, f)
            if os.path.exists(q):
                os.unlink(q)
                stale += 1
    as_bytes = int(a.as_gib * 1024**3)
    here = os.path.dirname(os.path.abspath(__file__))

    sa, sb = sha256_of(a.in_a), sha256_of(a.in_b)
    t0 = time.monotonic()
    d_abs = t0 + a.deadline
    # ★ 設定と対象の数を、最初に全部出します。
    print("これから %d 段を回します: %s" % (len(STAGES), " → ".join(n for n, _ in STAGES)))
    print("  対象: 固定入力対 %s。**4 演算。合計 76,592 三角形（保存値）**" % a.key)
    print("  入力 A %s（sha256 %s）" % (a.in_a, sa))
    print("  入力 B %s（sha256 %s）" % (a.in_b, sb))
    print("  **再量子化しません。この同じファイルを被検体へ渡します。**")
    print("  段の枠 %s / 合計 %g 秒 / 絶対期限 %g 秒 / 猶予 %g 秒（枠の内側）/ "
          "回収 %g 秒（枠の外）/ 単一ワーカー %d バイト / 被検体 %d スレッド"
          % ({n: b for n, b in STAGES}, sum(b for _, b in STAGES), a.deadline,
             a.grace, a.reap, as_bytes, a.threads))
    inj = {k: os.environ[k] for k in
           ("KRI_DIAG_TEST_UNREAPED_AT", "KRI_DIAG_TEST_SPAWN_DELAY",
            "KRI_DIAG_TEST_CHILD_DELAY", "KRI_DIAG_TEST_SPAWN_LOG")
           if k in os.environ}
    print("  試験専用の注入: %s" % (inj if inj else "無し"))
    print("  出力 %s（前の実行の成果物 %d 個を消しました）" % (out, stale))
    sys.stdout.flush()

    py = sys.executable
    cmds = {
        "I-0": [py, os.path.join(here, "inspect_output_controls.py"), out],
        "I-1": [py, os.path.join(here, "inspect_output.py"), "--stage", "identity",
                "--out", out, "--in-a", a.in_a, "--in-b", a.in_b],
        "I-2": [a.dumper, a.in_a, sa, a.in_b, sb, os.path.join(out, "b"),
                "6", str(a.threads), "1", "1", "0", "0", "16", "2", "1", "1"],
        "I-3": [py, os.path.join(here, "inspect_output.py"), "--stage", "anchor",
                "--out", out, "--saved", a.saved, "--key", a.key],
        "I-4": [py, os.path.join(here, "inspect_output.py"), "--stage", "measure",
                "--out", out, "--in-a", a.in_a, "--in-b", a.in_b,
                "--saved", a.saved, "--key", a.key],
    }

    state = {n: "未起動" for n, _ in STAGES}

    def report(rc):
        """★ **停止したときも (A) 実行状態をそのまま並べます**（§44.7 (E)）。"""
        print("実行状態: %s" % state)
        notrun = [n for n, _ in STAGES if state[n] == "未起動"]
        print("未起動の段: %s" % (notrun if notrun else "無し"))
        print("**案 II へは自動で進みません。**")
        print("出力: %s" % out)
        return rc

    for name, budget in STAGES:
        now = time.monotonic()
        # ★ **この枝は、いまの枠の算術では到達しません**（自己検定で確かめました）。
        #   段が「完了」するのは子が `t_term` より前に終わったときで、
        #   `t_term <= d_abs - grace` なので、完了した時刻 T は必ず `T < d_abs - grace`、
        #   すなわち `d_abs - T > grace`。**次の段は必ず起動できます。**
        #   最初の段も `grace < deadline` を引数検査で強制しているので当たりません。
        #   **残すのは、枠の算術を将来変えたときの備えです。空回りしていることを明記します。**
        if d_abs - now <= a.grace:
            print("[%s] ★ 残り %.1f 秒では終了猶予 %g 秒を確保できません。起動しません。"
                  % (name, d_abs - now, a.grace))
            return report(4)
        t_kill = min(now + budget, d_abs)
        t_term = max(now, t_kill - a.grace)
        log = os.path.join(out, "child_%s.log" % name)
        print("[%s] 起動: %s" % (name, " ".join(cmds[name])))
        print("[%s]   TERM は %.1f 秒後 / KILL は %.1f 秒後 / 回収の枠 %g 秒"
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
            print("[%s] ★ 回収できていません（pid=%d、群 %d）。後続を起動しません。"
                  % (name, pid, pid))
            return report(2)
        if st == sup.CUT:
            state[name] = "打ち切り"
            print("[%s] ★ 打ち切り、または期限後の回収。**未完了として区切ります。**" % name)
            return report(3)
        if code == 127:
            state[name] = "起動失敗"
            print("[%s] ★ 起動できなかった可能性があります（終了値 127）。後続を起動しません。" % name)
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
            print("[%s] ★ 成果物が無いか空です: %s。後続を起動しません。" % (name, ", ".join(miss)))
            return report(5)

    print("全 %d 段を回しました。経過 %.2f 秒（絶対期限 %g 秒）"
          % (len(STAGES), time.monotonic() - t0, a.deadline))
    print("**案 I は E2〜E5・C1・C2 の全数判定をせず、$Q$ の認定もしません。**")
    return report(0)


if __name__ == "__main__":
    sys.exit(main())

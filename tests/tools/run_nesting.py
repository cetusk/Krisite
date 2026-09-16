#!/usr/bin/env python3
"""保存配列 2 本に限定した検証を、**段ごとに監督して**回します。

**既存の監督（`run_gmp_diag.py` の `spawn` / `supervise`）を再利用します。**
監督は改修しません。

  * 段の枠 … 合成 90 秒 / A5 確認 240 秒 / 包含・算術 30 秒（**終了猶予込み**）
  * 全体の絶対期限 … `--deadline`（既定 360 秒）
  * `TERM` → `KILL` … `--grace`（既定 10 秒）。**枠の内側**
  * 回収の枠 … `--reap`（既定 30 秒）。**枠の外**
  * 仮想アドレス空間 … `--as-gib`（既定 4 GiB）。**単一ワーカー**

**後続を起動しない条件**（**自動延長も自動再試行も対象拡大もしません**）::

    回収不能 / 打ち切り / 子の終了値が 0 でない / 成果物が無いか空
    絶対期限までに終了猶予を確保できない

**ビルドも再量子化もしません。** 受理済みの保存配列を読むだけです。
"""
import argparse, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_gmp_diag as sup

STAGES = [("controls", 90), ("a5", 240), ("nesting", 30)]
ARTIFACTS = {"controls": ["stage_controls.txt"],
             "a5": ["stage_a5.txt"],
             "nesting": ["stage_nesting.txt"]}


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="")
    p.add_argument("--deadline", type=float, default=360.0)
    p.add_argument("--grace", type=float, default=10.0)
    p.add_argument("--reap", type=float, default=30.0)
    p.add_argument("--as-gib", type=float, default=4.0)
    p.add_argument("--arrays", default="data/logs/inspect/20260915-205108")
    p.add_argument("--saved", default="docs/evidence/gmp_diag_r1/cp3_gmp_results.txt")
    p.add_argument("--key", default="250394x45413")
    a = p.parse_args(argv)

    for nm, v in (("--deadline", a.deadline), ("--grace", a.grace),
                  ("--reap", a.reap), ("--as-gib", a.as_gib)):
        if not (v == v) or v in (float("inf"), float("-inf")) or v <= 0:
            p.error("%s は有限の正の数であること: %r" % (nm, v))
    if a.grace >= a.deadline:
        p.error("--grace (%g) は --deadline (%g) より小さいこと" % (a.grace, a.deadline))

    out = a.out or os.path.join("data", "logs", "nesting", time.strftime("%Y%m%d-%H%M%S"))
    os.makedirs(out, exist_ok=True)
    # ★ 成果物だけでなく**子のログも**消します。
    #   `spawn` の子は exec の前に落ちることがあり（`os._exit(127)`）、
    #   その場合ログを開かないので、**前の実行のログを表示してしまいます。**
    stale = 0
    for name, _ in STAGES:
        for f in ARTIFACTS[name] + ["child_%s.log" % name]:
            q = os.path.join(out, f)
            if os.path.exists(q):
                os.unlink(q)
                stale += 1
    as_bytes = int(a.as_gib * 1024**3)

    t0 = time.monotonic()
    d_abs = t0 + a.deadline
    print("これから %d 段を回します: %s" % (len(STAGES), " → ".join(n for n, _ in STAGES)))
    print("  保存配列 %s（**再量子化しません。ビルドもしません**）" % a.arrays)
    print("  保存物 %s / 鍵 %s" % (a.saved, a.key))
    print("  段の枠 %s / 合計 %g 秒 / 絶対期限 %g 秒 / 猶予 %g 秒（枠の内側）/ "
          "回収 %g 秒（枠の外）/ 単一ワーカー %d バイト"
          % ({n: b for n, b in STAGES}, sum(b for _, b in STAGES), a.deadline,
             a.grace, a.reap, as_bytes))
    inj = {k: os.environ[k] for k in
           ("KRI_DIAG_TEST_UNREAPED_AT", "KRI_DIAG_TEST_SPAWN_DELAY",
            "KRI_DIAG_TEST_CHILD_DELAY", "KRI_DIAG_TEST_SPAWN_LOG")
           if k in os.environ}
    print("  試験専用の注入: %s" % (inj if inj else "無し"))
    print("  出力 %s（前の実行の成果物 %d 個を消しました）" % (out, stale))
    sys.stdout.flush()

    here = os.path.dirname(os.path.abspath(__file__))
    for name, budget in STAGES:
        now = time.monotonic()
        if d_abs - now <= a.grace:
            print("[%s] ★ 残り %.1f 秒では終了猶予 %g 秒を確保できません。起動しません。"
                  % (name, d_abs - now, a.grace))
            return 4
        t_kill = min(now + budget, d_abs)
        t_term = max(now, t_kill - a.grace)
        log = os.path.join(out, "child_%s.log" % name)
        cmd = [sys.executable, os.path.join(here, "inspect_nesting.py"),
               "--stage", name, "--out", out, "--arrays", a.arrays,
               "--saved", a.saved, "--key", a.key]
        print("[%s] 起動: %s" % (name, " ".join(cmd)))
        print("[%s]   TERM は %.1f 秒後 / KILL は %.1f 秒後 / 回収の枠 %g 秒"
              % (name, t_term - now, t_kill - now, a.reap))
        sys.stdout.flush()

        pid = sup.spawn(cmd, dict(os.environ), log, as_bytes)
        code, rss, state = sup.supervise(pid, t_term, t_kill, hard_cap=a.reap)
        dt = time.monotonic() - now
        try:
            sys.stdout.write(open(log, encoding="utf-8", errors="replace").read())
        except OSError:
            pass
        print("[%s] 状態 %s / 終了値 %s / ピーク RSS %s / 経過 %.2f 秒"
              % (name, state, "不明" if code is None else code,
                 "不明" if rss is None else "%d KiB" % rss, dt))
        sys.stdout.flush()

        if state == sup.UNREAPED:
            print("[%s] ★ 回収できていません（pid=%d、群 %d）。後続を起動しません。"
                  % (name, pid, pid))
            return 2
        if state == sup.CUT:
            print("[%s] ★ 打ち切り、または期限後の回収。**未完了として区切ります。**"
                  "後続を起動しません。" % name)
            return 3
        if code != 0:
            print("[%s] ★ 子が終了値 %s で終わりました。後続を起動しません。" % (name, code))
            return 1
        miss = [f for f in ARTIFACTS[name]
                if not os.path.exists(os.path.join(out, f))
                or os.path.getsize(os.path.join(out, f)) == 0]
        if miss:
            print("[%s] ★ 成果物が無いか空です: %s。後続を起動しません。" % (name, ", ".join(miss)))
            return 5

    print("全 %d 段を回しました。経過 %.2f 秒（絶対期限 %g 秒）"
          % (len(STAGES), time.monotonic() - t0, a.deadline))
    print("出力: %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""GMP 整合性診断の投入（DESIGN-phase5-vertex-level.md §22.26 / §23.11）。

**予算と上限を守るのはこの層の仕事です。** 駆動は **1 対ずつ子プロセス**で回します。

  * **全体の期限 D は【全計画に 1 つ】**（既定 1200 秒 = 20 分）。
    **CP をまたいでも作り直しません。** 準備・終了猶予もその内側です。
  * **各対の期限** = ``min(per_pair, 残り − 猶予)``。**0 以下なら起動しません。**
  * **D − 猶予** までに ``TERM``、**D** までに ``KILL``。
  * **1 ワーカープロセスの仮想アドレス空間**を ``RLIMIT_AS`` で制限
    （**RSS の制限ではありません**）。
  * **ピーク RSS は `os.wait4` で【その子だけ】の値を採ります**
    （``RUSAGE_CHILDREN`` は先に終わった子の影響を受けます）。
  * **完了は「終了値 0」では数えません。** **その対の行が結果に在り、
    実施済みで判定が完了していること**を確かめます。

**OS による終了・回収の遅延があれば、超過として報告します。**
"""
import argparse
import os
import resource
import signal
import sys
import time

FIXED_COLS = 34


def parse_args(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True, help="駆動（thingi_cp1_gmp）")
    ap.add_argument("--target", action="append", required=True, metavar="一覧:計画",
                    help="`模型の一覧:計画の一覧` を CP ごとに繰り返し指定します")
    ap.add_argument("--args", required=True, help="駆動へ渡す引数列（第 2 引数以降）")
    ap.add_argument("--deadline", type=float, default=1200.0, help="全計画に 1 つの期限（秒）")
    ap.add_argument("--grace", type=float, default=10.0, help="終了猶予（秒。期限の内側）")
    ap.add_argument("--per-pair", type=float, default=300.0, help="1 対の上限（秒）")
    ap.add_argument("--as-gib", type=float, default=8.0, help="仮想アドレス空間の上限（GiB）")
    ap.add_argument("--logdir", required=True, help="対ごとの出力を残す場所")
    ap.add_argument("--run-id", default=None, help="実行の識別子（ログ名に入ります）")
    ap.add_argument("--dry-run", action="store_true", help="起動せず、計画だけ出す")
    return ap.parse_args(argv)


def preexec(as_bytes):
    def f():
        resource.setrlimit(resource.RLIMIT_AS, (as_bytes, as_bytes))
        os.setsid()  # **子を独立した群にして、群ごと止められるようにします**
    return f


def spawn(cmd, env, log_path, as_bytes):
    """fork + exec します。**`os.wait4` を使うため `subprocess` は通しません。**"""
    pid = os.fork()
    if pid == 0:  # 子
        try:
            fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(fd, 1)
            os.dup2(fd, 2)
            preexec(as_bytes)()
            os.execvpe(cmd[0], cmd, env)
        except BaseException:
            os._exit(127)
    return pid


def run_one(a, key, only_path, budget, log_path):
    """1 対を回します。戻り値: (状態, 終了値, 秒, その子のピーク RSS[KiB])。"""
    env = dict(os.environ, KRI_GMP_ONLY=only_path)
    cmd = [os.path.abspath(a.bin)] + [a.list_of_key] + a.args.split()
    t0 = time.monotonic()
    pid = spawn(cmd, env, log_path, int(a.as_gib * (1 << 30)))
    status, rc, rss = "ok", None, 0
    deadline = t0 + budget
    while True:
        wpid, st, ru = os.wait4(pid, os.WNOHANG)
        if wpid == pid:
            rc = os.waitstatus_to_exitcode(st) if hasattr(os, "waitstatus_to_exitcode") else (
                os.WEXITSTATUS(st) if os.WIFEXITED(st) else -os.WTERMSIG(st))
            rss = ru.ru_maxrss
            break
        if time.monotonic() >= deadline:
            status = "打ち切り"
            try:
                os.killpg(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            t_kill = time.monotonic() + a.grace
            while time.monotonic() < t_kill:
                wpid, st, ru = os.wait4(pid, os.WNOHANG)
                if wpid == pid:
                    rc, rss = -1, ru.ru_maxrss
                    break
                time.sleep(0.05)
            else:
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                _, st, ru = os.wait4(pid, 0)
                rc, rss = -9, ru.ru_maxrss
            break
        time.sleep(0.05)
    dt = time.monotonic() - t0
    if status == "ok" and rc != 0:
        status = "失敗"
    return status, rc, dt, rss


def row_complete(res_path, key):
    """**その対の行が在り、実施済みで判定が完了しているか**を確かめます。"""
    try:
        with open(res_path, errors="replace") as f:
            for line in f:
                t = line.split()
                if not t or t[0] != key:
                    continue
                return (len(t) >= FIXED_COLS and t[1] == "ok" and t[2] == "1"
                        and t[3] == "id1=ok" and all(t[18 + k] == "15" for k in range(4))
                        and all(t[22 + k] == "0" for k in range(4)))
    except OSError:
        return False
    return False


def main(argv=None):
    a = parse_args(argv)
    run_id = a.run_id or time.strftime("%Y%m%d_%H%M%S")
    plan = []  # (list_path, plan_path, key)
    for t in a.target:
        if ":" not in t:
            print(f"**--target は `一覧:計画` の形です**: {t}")
            return 2
        lp, pp = t.rsplit(":", 1)
        for f in (a.bin, lp, pp):
            if not os.path.exists(f):
                print(f"**ありません**: {f}")
                return 2
        keys = [ln.split()[0] for ln in open(pp) if ln.strip()]
        if not keys:
            print(f"**計画が空です**: {pp}")
            return 2
        for k in keys:
            plan.append((lp, pp, k))
    os.makedirs(a.logdir, exist_ok=True)

    print(f"**これから {len(plan)} 対を回します**（全計画に 1 つの期限 {a.deadline:.0f} 秒、"
          f"各対 {a.per_pair:.0f} 秒、猶予 {a.grace:.0f} 秒、"
          f"仮想アドレス空間 {a.as_gib} GiB、実行 {run_id}）")
    for lp, pp, k in plan:
        print(f"  対象: {k}  （一覧 {lp}）")
    if a.dry_run:
        print("**計画だけ出しました。起動していません。**")
        return 0

    t_end = time.monotonic() + a.deadline
    done, skipped, bad = 0, 0, 0
    for lp, pp, k in plan:
        base = os.path.splitext(lp)[0]
        res_path = base + "_gmp_results.txt"
        stem = os.path.basename(base)
        if row_complete(res_path, k):
            print(f"  {k}: **済み**（検証した再利用。起動しません）")
            done += 1
            continue
        remain = t_end - time.monotonic()
        budget = min(a.per_pair, remain - a.grace)
        if budget <= 0:
            print(f"  {k}: **予算切れのため起動しません**（残り {remain:.1f} 秒）")
            skipped += 1
            continue
        only_path = os.path.join(a.logdir, f"{stem}_{run_id}_one.txt")
        with open(only_path, "w") as f:
            f.write(k + "\n")
        log_path = os.path.join(a.logdir, f"{stem}_{run_id}_{k.replace('/', '_')}.log")
        a.list_of_key = lp
        env_plan = os.environ.copy()
        env_plan["KRI_GMP_PLAN"] = pp
        os.environ["KRI_GMP_PLAN"] = pp
        status, rc, dt, rss = run_one(a, k, only_path, budget, log_path)
        okrow = row_complete(res_path, k)
        print(f"  {k}: {status}（終了値 {rc}、{dt:.1f} 秒、ピーク RSS {rss / 1024:.0f} MiB、"
              f"期限 {budget:.1f} 秒、行 {'完全' if okrow else '**不完全**'}）  ログ: {log_path}")
        if status != "ok" or not okrow:
            bad += 1
            try:
                for ln in open(log_path, errors="replace").read().splitlines()[-5:]:
                    print(f"      | {ln}")
            except OSError:
                pass
        else:
            done += 1
    over = time.monotonic() - t_end
    print(f"\n**済み {done} / 失敗・打ち切り {bad} / 未起動 {skipped}**")
    if over > 0:
        print(f"**★ 期限を {over:.1f} 秒超過しました**（OS の終了・回収の遅延）")
    return 0 if (bad == 0 and skipped == 0) else 1


if __name__ == "__main__":
    sys.exit(main())

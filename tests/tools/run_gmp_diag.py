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
import math
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
    ap.add_argument("--synth-check", default=None, metavar="実行ファイル",
                    help="**合成検定**（検査器の正例・負例）。**同じ期限の内側で、対より先に回します**")
    ap.add_argument("--dry-run", action="store_true", help="起動せず、計画だけ出す")
    return ap.parse_args(argv)


def spawn(cmd, env, log_path, as_bytes):
    """fork + exec します。**`os.wait4` を使うため `subprocess` は通しません。**

    **★ 子は【いちばん最初に】`setsid()` します**（§29）。
    **ログを開いてから群を作ると、その間に親が `killpg` を撃っても届きません。**

    **試験専用の遅延**（既定では読まれても 0。本番の経路に影響しません）:

    * ``KRI_DIAG_TEST_SPAWN_DELAY`` … 親が **起動の前**に待つ秒数（§28）
    * ``KRI_DIAG_TEST_CHILD_DELAY`` … 子が **群を作る前**に待つ秒数（§29）
    """
    d = os.environ.get("KRI_DIAG_TEST_SPAWN_DELAY")
    if d:
        time.sleep(float(d))
    pid = os.fork()
    if pid == 0:  # 子
        try:
            cd = os.environ.get("KRI_DIAG_TEST_CHILD_DELAY")
            if cd:
                time.sleep(float(cd))
            os.setsid()          # ★ まずここ
            fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(fd, 1)
            os.dup2(fd, 2)
            resource.setrlimit(resource.RLIMIT_AS, (as_bytes, as_bytes))
            os.execvpe(cmd[0], cmd, env)
        except BaseException:
            os._exit(127)
    return pid


def signal_child(pid, sig):
    """**群へ送り、群がまだ無ければ【子へ直接】送ります**（§29）。

    **送れたかを返します。** **送信の失敗を、成功と同じに扱わないため**です。
    """
    try:
        os.killpg(pid, sig)
        return True
    except (ProcessLookupError, PermissionError):
        pass
    try:
        os.kill(pid, sig)
        return True
    except ProcessLookupError:
        return False


def exit_code(st):
    """待機状態を終了値に直します（**シグナルは負で返します**）。"""
    if hasattr(os, "waitstatus_to_exitcode"):
        try:
            return os.waitstatus_to_exitcode(st)
        except ValueError:
            pass
    return os.WEXITSTATUS(st) if os.WIFEXITED(st) else -os.WTERMSIG(st)


def supervise(pid, t_term, t_kill, hard_cap=30.0):
    """**絶対の期限で子を待ちます。** 戻り値: (終了値, ピーク RSS[KiB], 期限を守れなかったか)。

    * ``t_term`` … **遅くともこの時刻に `TERM`**（全体の期限 − 猶予 で上から抑えます）
    * ``t_kill`` … **遅くともこの時刻に `KILL`**（全体の期限で上から抑えます）

    **★ 3 つを分けます**（§29）。

    1. **終了値**（子が何を返したか）
    2. **期限を守れたか**（**回収の時刻が `t_term` を過ぎていれば、守れていません**）
    3. **停止を送ったか**

    **「終了値 0 で回収できた」を、そのまま成功に渡しません。**
    **監督の開始が遅れた場合、子が期限を過ぎて自然終了していることがあります。**

    **★ 無期限の待機はしません。** `KILL` は届くまで毎周回送り直し、
    それでも回収できなければ ``hard_cap`` 秒で諦めて報告します。
    """
    sent_term = False
    while True:
        wpid, st, ru = os.wait4(pid, os.WNOHANG)
        if wpid == pid:
            # **回収できた時刻で、期限を守れたかを判定します**
            late = time.monotonic() > t_term
            return exit_code(st), ru.ru_maxrss, (sent_term or late)
        now = time.monotonic()
        if not sent_term and now >= t_term:
            sent_term = True
            if not signal_child(pid, signal.SIGTERM):
                # **送れませんでした**（群がまだ無い / もう居ない）。
                # **次の周回で回収を試み、届かなければ KILL へ進みます。**
                pass
        if sent_term and now >= t_kill:
            signal_child(pid, signal.SIGKILL)  # **届くまで毎周回送り直します**
            if now >= t_kill + hard_cap:
                print(f"**子 {pid} を {hard_cap:.0f} 秒かけても回収できませんでした**")
                return None, 0, True
        time.sleep(0.01)


def deadlines(a, t_end, budget):
    """**起動の【前】に、絶対の TERM / KILL 時刻を決めます**（§28）。

    **どちらも全体の期限 D で上から抑えます。**
    """
    now = time.monotonic()
    t_term = min(now + budget, t_end - a.grace)
    t_kill = min(t_term + a.grace, t_end)
    return t_term, t_kill


def run_one(a, key, only_path, budget, log_path, t_end):
    """1 対を回します。戻り値: (状態, 終了値, 秒, その子のピーク RSS[KiB])。"""
    env = dict(os.environ, KRI_GMP_ONLY=only_path)
    cmd = [os.path.abspath(a.bin), a.list_of_key] + a.args.split()
    t_term, t_kill = deadlines(a, t_end, budget)   # ★ spawn の【前】に決めます
    t0 = time.monotonic()
    pid = spawn(cmd, env, log_path, int(a.as_gib * (1 << 30)))
    sd = os.environ.get("KRI_DIAG_TEST_SUPERVISE_DELAY")   # 試験専用（§29）
    if sd:
        time.sleep(float(sd))
    rc, rss, cut = supervise(pid, t_term, t_kill)
    dt = time.monotonic() - t0
    status = "打ち切り" if cut else ("ok" if rc == 0 else "失敗")
    return status, rc, dt, rss


def check_only(a, lp, pp, only_path, log_path, budget, t_end):
    """**駆動に照合だけさせます**（量子化もブール演算もしません）。

    **★ 再利用してよいかを、この層では決めません**（§25.1）。
    **meta・入力・バイナリ・設定・行の検査は、駆動の 1 か所だけが持ちます。**

    **★ 期限つきで監督します**（§27）。**打ち切ったら `None`** を返します。

    戻り値: 0 = すべて済み（**検証した再利用**）/ 3 = 回す対がある /
            `None` = 打ち切り / それ以外 = 拒否
    """
    env = dict(os.environ, KRI_GMP_ONLY=only_path, KRI_GMP_PLAN=pp, KRI_GMP_CHECK_ONLY="1")
    cmd = [os.path.abspath(a.bin), lp] + a.args.split()
    t_term, t_kill = deadlines(a, t_end, budget)   # ★ spawn の【前】に決めます
    pid = spawn(cmd, env, log_path, int(a.as_gib * (1 << 30)))
    rc, _, cut = supervise(pid, t_term, t_kill)
    return None if cut else rc


def main(argv=None):
    a = parse_args(argv)
    # **★ 期限は、準備を始める【前】から数えます**（§25.3）
    t_begin = time.monotonic()
    for name, v in (("--deadline", a.deadline), ("--grace", a.grace),
                    ("--per-pair", a.per_pair), ("--as-gib", a.as_gib)):
        if not math.isfinite(v) or v <= 0:
            print(f"**{name} が不正です**: {v}")
            return 2
    if a.grace >= a.deadline:
        print(f"**--grace は --deadline より小さくしてください**: {a.grace} >= {a.deadline}")
        return 2
    t_end = t_begin + a.deadline
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

    # ---- 合成検定（**同じ期限の内側**。§26）------------------------------------
    #
    # **検査器が破れを検出できることを、実データの前に確かめます**
    # （`CLAUDE.md`「判定器を先に検定した形が要点」）。
    # **通らなければ、対は 1 つも起動しません。**
    if a.synth_check:
        if not os.path.exists(a.synth_check):
            print(f"**合成検定の実行ファイルがありません**: {a.synth_check}")
            return 2
        remain = t_end - time.monotonic()
        budget = min(a.per_pair, remain - a.grace)
        if budget <= 0:
            print(f"**合成検定の予算がありません**（残り {remain:.1f} 秒）")
            return 1
        slog = os.path.join(a.logdir, f"synth_{run_id}.log")
        t_term, t_kill = deadlines(a, t_end, budget)   # ★ spawn の【前】に決めます
        pid = spawn([os.path.abspath(a.synth_check)], dict(os.environ), slog,
                    int(a.as_gib * (1 << 30)))
        t0 = time.monotonic()
        rc, _, cut = supervise(pid, t_term, t_kill)
        if cut:
            rc = None
        print(f"  合成検定: {'通過' if rc == 0 else ('**打ち切り**' if rc is None else '**不通過**')}"
              f"（終了値 {rc}、{time.monotonic() - t0:.1f} 秒）  ログ: {slog}")
        if rc != 0:
            tail(slog)
            print("**合成検定が通らないので、対を 1 つも起動しません。**")
            return 1

    done, skipped, bad, stopped = 0, 0, 0, False
    for lp, pp, k in plan:
        stem = os.path.basename(os.path.splitext(lp)[0])
        if stopped:
            print(f"  {k}: **停止条件に当たったので起動しません**")
            skipped += 1
            continue
        # ---- 0. 予算（**照合も予算の内側です**。§25.3）----
        remain = t_end - time.monotonic()
        budget = min(a.per_pair, remain - a.grace)
        if budget <= 0:
            print(f"  {k}: **予算切れのため起動しません**（残り {remain:.1f} 秒）")
            skipped += 1
            stopped = True
            continue
        only_path = os.path.join(a.logdir, f"{stem}_{run_id}_one.txt")
        with open(only_path, "w") as f:
            f.write(k + "\n")
        # ---- 1. 照合（駆動に聞きます。量子化もブール演算もしません）----
        clog = os.path.join(a.logdir, f"{stem}_{run_id}_{k}_check.log")
        rc = check_only(a, lp, pp, only_path, clog, budget, t_end)
        if rc is None:
            print(f"  {k}: **照合が期限で打ち切られました**。**未検証を成功に数えません。**"
                  f"  ログ: {clog}")
            bad += 1
            stopped = True
            continue
        if rc == 0:
            print(f"  {k}: **済み**（駆動が照合した再利用。起動しません）  ログ: {clog}")
            done += 1
            continue
        if rc != 3:
            print(f"  {k}: **照合に失敗しました**（終了値 {rc}）。**再計算しません。**  ログ: {clog}")
            tail(clog)
            bad += 1
            stopped = True          # **拒否は停止条件です**（§25.3）
            continue
        # ---- 2. 本計算（照合に使った時間も引いた残りで）----
        remain = t_end - time.monotonic()
        budget = min(a.per_pair, remain - a.grace)
        if budget <= 0:
            print(f"  {k}: **照合の後で予算切れになりました**（残り {remain:.1f} 秒）")
            skipped += 1
            stopped = True
            continue
        log_path = os.path.join(a.logdir, f"{stem}_{run_id}_{k}.log")
        a.list_of_key = lp
        os.environ["KRI_GMP_PLAN"] = pp
        status, rc, dt, rss = run_one(a, k, only_path, budget, log_path, t_end)
        # ---- 3. 完了は、もう一度【駆動に照合させて】数えます ----
        vlog = os.path.join(a.logdir, f"{stem}_{run_id}_{k}_verify.log")
        vremain = t_end - time.monotonic()
        vbudget = min(a.per_pair, vremain - a.grace)
        if vbudget <= 0:
            # **照合できないまま「済み」にしません**（§27）
            print(f"  {k}: **事後の照合に予算が残っていません**（残り {vremain:.1f} 秒）。"
                  f"**未検証を成功に数えません。**")
            bad += 1
            stopped = True
            continue
        vrc = check_only(a, lp, pp, only_path, vlog, vbudget, t_end)
        okrow = (vrc == 0)
        print(f"  {k}: {status}（終了値 {rc}、{dt:.1f} 秒、ピーク RSS {rss / 1024:.0f} MiB、"
              f"期限 {budget:.1f} 秒、照合 {'通過' if okrow else '**不通過**'}）  ログ: {log_path}")
        if status != "ok" or not okrow:
            bad += 1
            tail(log_path)
            if not okrow:
                tail(vlog)
            stopped = True          # **失敗したら次へ進みません**（§25.3）
        else:
            done += 1
    over = time.monotonic() - t_end
    print(f"\n**済み {done} / 失敗・打ち切り {bad} / 未起動 {skipped}**")
    if over > 0:
        print(f"**★ 期限を {over:.1f} 秒超過しました**（OS の終了・回収の遅延）")
    return 0 if (bad == 0 and skipped == 0) else 1


def tail(path, n=5):
    """失敗したときに、その場でログの末尾を出します。"""
    try:
        for ln in open(path, errors="replace").read().splitlines()[-n:]:
            print(f"      | {ln}")
    except OSError:
        pass


if __name__ == "__main__":
    sys.exit(main())

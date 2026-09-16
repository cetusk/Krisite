#!/usr/bin/env python3
"""停止・回収の**対照**（外部監督が働くことを確かめるための、最小の子）。

**これは被検体ではありません。** 監督（`run_inspect.py` / `run_gmp_diag.py`）が
「止まらない子を止められるか」「止めた子を回収できたか」を確かめるための道具です。

  * `never` … `SIGTERM` を無視して回り続ける。**`KILL` でしか止まりません**
  * `term`  … `SIGTERM` を受けたら終了値 3 で終わる
  * `alive` … 指定した **argv の要素**を持つプロセスが**走っているか**を `/proc` から調べる

**`alive` が示せること・示せないこと**（**主張をここまでに限定してください**）:

  * 示せるのは「**その argv 要素を持つプロセスが走っているか**」までです。
  * **回収漏れ（ゾンビ）は、この照合では見つかりません。**
    回収されていない子は `/proc/<pid>/cmdline` が**空**なので、argv と突き合わせられません。
    **`回収不能` が意味するのはまさに回収されていない子なので、
    「一致 0 件」を「回収漏れが無い」の根拠にしないでください。**
    回収できたことの根拠は、監督の `wait4` が終了値を返したことのほうです。
    そのため `alive` は、**照合できなかったゾンビの件数と一覧を別に出します**。
  * `ps -ef | grep <語>` / `pgrep -f <語>` は使いません。
    **あの形は grep 自身と呼び出し側シェルの命令行に語が現れるため、常に偽陽性を出します**
    （実際に踏みました）。照合は `cmdline` を `\0` で切った**要素との完全一致**で行い、
    **自分自身と自分の祖先を除きます**。
  * **監督自身の命令行にも語が現れます**（`... -- python3 stopctl.py never`）。
    祖先でなければ一致として出ます。**偽陽性（安全側）**ですが、
    一致の行に命令行を全部出すので読み手が区別できます。
  * **読めなかった `/proc` の項目は件数を出します。** 「一致 0 件」は
    「**見られた範囲で 0 件**」です。

終了値: **0** = 一致 0 件かつゾンビ 0 件 / **1** = 一致あり / **2** = 一致 0 件だがゾンビあり（照合不能）
"""
import os, signal, sys, time


def run_never():
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    sys.stderr.write("止まらない子: TERM を無視します\n")
    sys.stderr.flush()
    while True:
        time.sleep(0.05)


def run_term():
    def on_term(sig, frm):
        sys.stderr.write("TERM を受けたので終了値 3 で終わります\n")
        sys.stderr.flush()
        sys.exit(3)
    signal.signal(signal.SIGTERM, on_term)
    sys.stderr.write("止まる子: TERM で終了値 3\n")
    sys.stderr.flush()
    while True:
        time.sleep(0.05)


def _stat(pid):
    """`(comm, state, ppid)` を返します。**コマンド名に空白や `)` が入り得る**ので、
    最後の `)` で切ります。読めなければ `None`。"""
    try:
        with open("/proc/%d/stat" % pid) as f:
            s = f.read()
        head, tail = s.split("(", 1)[1].rsplit(")", 1)
        fields = tail.split()
        return head, fields[0], int(fields[1])
    except Exception:
        return None


def run_alive(target):
    me = os.getpid()
    anc, p = set(), me
    while p > 1:
        anc.add(p)
        st = _stat(p)
        p = st[2] if st else 0
        if p in anc:
            break

    hit, zombie, unread = [], [], 0
    for e in os.listdir("/proc"):
        if not e.isdigit():
            continue
        pid = int(e)
        if pid in anc:
            continue
        try:
            with open("/proc/%d/cmdline" % pid, "rb") as f:
                argv = [a for a in f.read().split(b"\0") if a]
        except (FileNotFoundError, ProcessLookupError):
            continue          # 走査の途中で終わったもの。取りこぼしではありません
        except OSError:
            unread += 1
            continue
        if argv:
            if any(a.decode("utf-8", "replace") == target for a in argv):
                hit.append((pid, b" ".join(argv).decode("utf-8", "replace")))
            continue
        # cmdline が空 = カーネルスレッド、または【回収されていない子（ゾンビ）】。
        # ゾンビは argv と突き合わせられないので、target に関わらず全件を出します。
        st = _stat(pid)
        if st and st[1] == "Z":
            zombie.append((pid, st[0], st[2]))

    for pid, c in hit:
        print("★ 走行中 pid=%d  %s" % (pid, c))
    for pid, comm, ppid in zombie:
        print("★ ゾンビ pid=%d  comm=%s  親 pid=%d （cmdline が空なので照合できません）"
              % (pid, comm, ppid))
    print("一致 %d 件 / ゾンビ %d 件 / 読めなかった項目 %d 件"
          "（自分と祖先 %d 個を除外。照合は argv の要素の完全一致）"
          % (len(hit), len(zombie), unread, len(anc)))
    if hit:
        return 1
    return 2 if zombie else 0


def main(argv):
    if len(argv) < 2 or argv[1] not in ("never", "term", "alive"):
        sys.stderr.write("使い方: stopctl.py never|term|alive [照合する argv の要素]\n")
        return 64
    mode = argv[1]
    if mode == "never":
        run_never()
    elif mode == "term":
        run_term()
    else:
        if len(argv) < 3:
            sys.stderr.write("alive には照合する argv の要素が要ります\n")
            return 64
        return run_alive(argv[2])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

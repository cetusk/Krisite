#!/bin/bash
# 負の対照の判定器（`SPEC-phase5.md` §5.10.14.74 の 1）。
#
#   bash tests/tools/negctl_check.sh <ログ> <期待する読み出し箇所（例: to_mesh.hpp）>
#
# **「ログのどこかに 2 つの語がある」では合格にしません。**
# **問題の【読み出しのスタック】だけを見ます** — 解放側・確保側のスタックに
# その名前があっても受理しません。
#
# ASan の出力は次の順です。
#
#   ==N==ERROR: AddressSanitizer: heap-use-after-free on address ...
#   READ of size 8 at 0x... thread T0        ← ここから
#       #0 ... in ... <file>:<line>
#       #1 ...
#   0x... is located ... inside of ...       ← ここまで
#   freed by thread T0 here:                 ← 以降は見ない
#   previously allocated by thread T0 here:
set -u
log="${1:?ログを渡してください}"
want="${2:-to_mesh.hpp}"
[ -f "$log" ] || { echo "ログがありません: $log" >&2; exit 2; }

awk -v want="$want" '
    /ERROR: AddressSanitizer: heap-use-after-free/ { err = 1; next }
    err && /^[[:space:]]*(READ|WRITE) of size/     { acc = 1; next }
    acc && (/^[[:space:]]*0x[0-9a-fA-F]+ is located/ || /^[[:space:]]*freed by/ ||
            /^[[:space:]]*previously allocated/ || /^SUMMARY:/ || /^$/) { acc = 0 }
    acc && index($0, want) > 0 { found = 1 }
    END { exit found ? 0 : 1 }
' "$log"

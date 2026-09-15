#!/usr/bin/env bash
# 投入スクリプトの自己検査（DESIGN-phase5-vertex-level.md §22.26）
#
# **模擬のバイナリだけを使います。** 実データも本物の駆動も動かしません。
set -u
R="python3 tests/tools/run_gmp_diag.py"
OK=0; NG=0
ok() { OK=$((OK+1)); printf '  OK   %s\n' "$1"; }
ng() { NG=$((NG+1)); printf '  **NG** %s\n' "$1"; }
chk() { if [ "$2" = "$3" ]; then ok "$1（$2）"; else ng "$1: 期待 $3、実測 $2"; fi; }

T=$(mktemp -d)
# 模擬の駆動: 環境変数で振る舞いを変えます
cat > "$T/mock" <<'M'
#!/usr/bin/env bash
case "${MOCK:-ok}" in
  ok)    exit 0 ;;
  slow)  sleep 30; exit 0 ;;
  fail)  exit 1 ;;
  hog)   python3 -c "a=bytearray(64*1024*1024*1024)"; exit $? ;;   # 64 GiB を確保しにいく
esac
M
chmod +x "$T/mock"
printf 'a\nb\nc\n' > "$T/only.txt"
: > "$T/list.txt"

echo "# 投入スクリプトの自己検査"
echo
echo "## 1. 計画だけ（起動しない）"
out=$($R --bin "$T/mock" --list "$T/list.txt" --only "$T/only.txt" --args "0" --dry-run); rc=$?
chk "終了値" "$rc" 0
chk "件数を先に出す" "$(echo "$out" | grep -c 'これから 3 対を回します')" 1
chk "起動していないと言う" "$(echo "$out" | grep -c '起動していません')" 1

echo
echo "## 2. 正常"
out=$(MOCK=ok $R --bin "$T/mock" --list "$T/list.txt" --only "$T/only.txt" --args "0"); rc=$?
chk "終了値" "$rc" 0
chk "済み 3" "$(echo "$out" | grep -c '済み 3 / 失敗・打ち切り 0 / 未起動 0')" 1

echo
echo "## 3. 失敗が終了値に届く"
out=$(MOCK=fail $R --bin "$T/mock" --list "$T/list.txt" --only "$T/only.txt" --args "0"); rc=$?
chk "終了値" "$rc" 1
chk "失敗 3" "$(echo "$out" | grep -c '失敗・打ち切り 3')" 1

echo
echo "## 4. 対ごとの期限で打ち切る"
t0=$(date +%s)
out=$(MOCK=slow $R --bin "$T/mock" --list "$T/list.txt" --only "$T/only.txt" --args "0" \
      --deadline 60 --grace 2 --per-pair 3); rc=$?
t1=$(date +%s)
chk "終了値" "$rc" 1
chk "3 対とも打ち切り" "$(echo "$out" | grep -cE '^  [a-c]: 打ち切り')" 3
chk "全体が 20 秒以内に戻る（3 対 x 3 秒 + 猶予）" "$([ $((t1-t0)) -le 20 ] && echo はい || echo いいえ)" はい

echo
echo "## 5. 全体の期限で【起動しない】"
out=$(MOCK=slow $R --bin "$T/mock" --list "$T/list.txt" --only "$T/only.txt" --args "0" \
      --deadline 8 --grace 2 --per-pair 5); rc=$?
chk "終了値" "$rc" 1
n=$(echo "$out" | grep -c '予算切れのため起動しません')
chk "予算切れで起動しない対がある（1 件以上）" "$([ "$n" -ge 1 ] && echo はい || echo いいえ)" はい

echo
echo "## 6. 仮想アドレス空間の上限"
printf 'a\n' > "$T/one.txt"
out=$(MOCK=hog $R --bin "$T/mock" --list "$T/list.txt" --only "$T/one.txt" --args "0" \
      --deadline 60 --per-pair 30 --as-gib 1); rc=$?
chk "終了値（確保に失敗して異常終了）" "$rc" 1
chk "上限が効いた証拠" "$(echo "$out" | grep -c 'MemoryError')" 1
chk "ピーク RSS を出す" "$(echo "$out" | grep -c 'ピーク RSS')" 1

echo
echo "## 7. 入力が無い"
rc=$($R --bin "$T/mock" --list "$T/list.txt" --only "$T/nope.txt" --args "0" > /dev/null 2>&1; echo $?)
chk "終了値" "$rc" 2
: > "$T/empty.txt"
rc=$($R --bin "$T/mock" --list "$T/list.txt" --only "$T/empty.txt" --args "0" > /dev/null 2>&1; echo $?)
chk "空の一覧 → 終了値" "$rc" 2

echo
printf '**OK %d / NG %d**\n' "$OK" "$NG"
[ "$NG" -eq 0 ]

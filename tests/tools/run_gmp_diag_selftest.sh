#!/usr/bin/env bash
# 投入の層の自己検査（DESIGN-phase5-vertex-level.md §23.11）
#
# **模擬のバイナリだけを使います。** 実データも本物の駆動も動かしません。
# **実駆動とつないだ試験は `gmp_diag_integration.sh`** です（役割が違います）。
set -u
R="python3 tests/tools/run_gmp_diag.py"
OK=0; NG=0
ok() { OK=$((OK+1)); printf '  OK   %s\n' "$1"; }
ng() { NG=$((NG+1)); printf '  **NG** %s\n' "$1"; }
chk() { if [ "$2" = "$3" ]; then ok "$1（$2）"; else ng "$1: 期待 $3、実測 $2"; fi; }

T=$(mktemp -d); mkdir -p "$T/log"
cat > "$T/mock" <<'M'
#!/usr/bin/env bash
# 模擬の駆動。**結果の行は自分で書きます**（投入の層が「行の完全性」を見るため）
key=$(head -1 "${KRI_GMP_ONLY:-/dev/null}" 2>/dev/null | awk '{print $1}')
base="${1%.txt}"
row() {  # 34 列の完全な行
  printf '%s ok 1 id1=ok id2=ok 12 12 0.1 %016x %016x 1 1 0 0 4 4 4 4 15 15 15 15 0 0 0 0 0 0 1 1 1 1 1 1\n' \
    "$key" 1 2 >> "${base}_gmp_results.txt"
}
# **照合専用**（KRI_GMP_CHECK_ONLY=1）。0 = 済み / 3 = 回す対がある / 2 = 拒否
if [ "${KRI_GMP_CHECK_ONLY:-0}" = "1" ]; then
  case "${MOCK_CHECK:-auto}" in
    reject) exit 2 ;;
    done)   exit 0 ;;
    slow)   sleep 0.7; exit 0 ;;          # ★ 仕様担当が再現した条件
    hang)   sleep 300; exit 0 ;;          # 止まったまま
    *) grep -q "^$key " "${base}_gmp_results.txt" 2>/dev/null && exit 0 || exit 3 ;;
  esac
fi
case "${MOCK:-ok}" in
  ok)    row; exit 0 ;;
  slow)  sleep 30; exit 0 ;;
  fail)  exit 1 ;;
  norow) exit 0 ;;                       # 終了値 0 なのに行を書かない
  hog)   python3 -c "a=bytearray(64*1024*1024*1024)"; exit $? ;;
esac
M
chmod +x "$T/mock"
printf 'a1xa2\nb1xb2\nc1xc2\n' > "$T/plan.txt"
: > "$T/list.txt"
TGT="$T/list.txt:$T/plan.txt"
clean() { rm -f "$T/list_gmp_results.txt"; }

echo "# 投入の層の自己検査"
echo
echo "## 1. 計画だけ（起動しない）"
out=$($R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" --dry-run); rc=$?
chk "終了値" "$rc" 0
chk "件数を先に出す" "$(echo "$out" | grep -c 'これから 3 対を回します')" 1
chk "起動していないと言う" "$(echo "$out" | grep -c '起動していません')" 1

echo
echo "## 2. 正常（行が完全なら済み）"
clean; out=$(MOCK=ok $R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" --run-id s2); rc=$?
chk "終了値" "$rc" 0
chk "済み 3" "$(echo "$out" | grep -c '済み 3 / 失敗・打ち切り 0 / 未起動 0')" 1

echo
echo "## 3. ★ 終了値 0 でも、行が無ければ済みにしない"
clean; out=$(MOCK=norow $R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" --run-id s3); rc=$?
chk "終了値" "$rc" 1
chk "照合が不通過（1 件目で停止）" "$(echo "$out" | grep -c '照合 \*\*不通過\*\*')" 1

echo
echo "## 4. 失敗が終了値に届く"
clean; out=$(MOCK=fail $R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" --run-id s4); rc=$?
chk "終了値" "$rc" 1
chk "1 件目で停止（失敗 1 / 未起動 2）" "$(echo "$out" | grep -c '失敗・打ち切り 1 / 未起動 2')" 1

echo
echo "## 5. 対ごとの期限で打ち切る"
clean; t0=$(date +%s)
out=$(MOCK=slow $R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" --run-id s5 \
      --deadline 60 --grace 2 --per-pair 3); rc=$?
t1=$(date +%s)
chk "終了値" "$rc" 1
chk "1 件目で打ち切り、以降は起動しない" "$(echo "$out" | grep -cE '^  [abc][0-9]x[abc][0-9]: 打ち切り')" 1
chk "全体が 20 秒以内に戻る" "$([ $((t1-t0)) -le 20 ] && echo はい || echo いいえ)" はい

echo
echo "## 6. 全計画に 1 つの期限（起動しない対が出る）"
# **照合も予算の内側**なので、予算が無ければ 1 対目から起動しません（決定的）
clean; out=$(MOCK=ok $R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" --run-id s6 \
      --deadline 2.001 --grace 2 --per-pair 5); rc=$?
chk "終了値" "$rc" 1
chk "1 対目で予算切れ（起動前か照合後）" "$(echo "$out" | grep -c '予算切れ')" 1
chk "以降は停止条件で起動しない" "$(echo "$out" | grep -c '停止条件に当たったので起動しません')" 2
chk "行を作っていない" "$([ -f "$T/list_gmp_results.txt" ] && echo あり || echo なし)" なし

echo
echo "## 7. CP を 2 つ渡しても期限は 1 つ"
printf 'd1xd2\n' > "$T/plan2.txt"; : > "$T/list2.txt"
clean; out=$(MOCK=ok $R --bin "$T/mock" --target "$TGT" --target "$T/list2.txt:$T/plan2.txt" \
      --args "0" --logdir "$T/log" --run-id s7 --dry-run)
chk "4 対（2 つの CP）" "$(echo "$out" | grep -c 'これから 4 対を回します')" 1
chk "期限は 1 つと表示" "$(echo "$out" | grep -c '全計画に 1 つの期限')" 1

echo
echo "## 8. 仮想アドレス空間の上限"
printf 'a1xa2\n' > "$T/one.txt"
clean; out=$(MOCK=hog $R --bin "$T/mock" --target "$T/list.txt:$T/one.txt" --args "0" \
      --logdir "$T/log" --run-id s8 --deadline 60 --per-pair 30 --as-gib 1); rc=$?
chk "終了値" "$rc" 1
chk "上限が効いた証拠" "$(echo "$out" | grep -c 'MemoryError')" 1
chk "ピーク RSS を出す" "$(echo "$out" | grep -c 'ピーク RSS')" 1

echo
echo "## 9. 入力が無い / 形が違う"
rc=$($R --bin "$T/mock" --target "$T/list.txt:$T/nope.txt" --args "0" --logdir "$T/log" > /dev/null 2>&1; echo $?)
chk "計画が無い → 終了値" "$rc" 2
: > "$T/empty.txt"
rc=$($R --bin "$T/mock" --target "$T/list.txt:$T/empty.txt" --args "0" --logdir "$T/log" > /dev/null 2>&1; echo $?)
chk "空の計画 → 終了値" "$rc" 2
rc=$($R --bin "$T/mock" --target "$T/list.txt" --args "0" --logdir "$T/log" > /dev/null 2>&1; echo $?)
chk "target の形が違う → 終了値" "$rc" 2


echo "## 10. 不正な予算の指定"
for bad in "--deadline -5" "--deadline nan" "--per-pair 0" "--as-gib -1" "--grace 1200"; do
  rc=$($R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" $bad > /dev/null 2>&1; echo $?)
  chk "$bad → 終了値" "$rc" 2
done


echo "## 11. 合成検定を期限の内側で回す"
printf '#!/usr/bin/env bash\nexit ${SYNTH_RC:-0}\n' > "$T/synth"; chmod +x "$T/synth"
clean; out=$(MOCK=ok SYNTH_RC=0 $R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" \
      --run-id s11a --synth-check "$T/synth"); rc=$?
chk "通過 → 終了値" "$rc" 0
chk "合成検定を先に回す" "$(echo "$out" | grep -c '合成検定: 通過')" 1
chk "対も回る" "$(echo "$out" | grep -c '済み 3')" 1

clean; out=$(MOCK=ok SYNTH_RC=1 $R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" \
      --run-id s11b --synth-check "$T/synth"); rc=$?
chk "不通過 → 終了値" "$rc" 1
chk "対を 1 つも起動しない" "$(echo "$out" | grep -c '対を 1 つも起動しません')" 1
chk "結果を作っていない" "$([ -f "$T/list_gmp_results.txt" ] && echo あり || echo なし)" なし

rc=$($R --bin "$T/mock" --target "$TGT" --args "0" --logdir "$T/log" --run-id s11c \
     --synth-check "$T/nope" > /dev/null 2>&1; echo $?)
chk "合成検定が無い → 終了値" "$rc" 2


echo "## 12. ★ 照合の子が期限で止まること（仕様担当の再現）"
# 照合だけ 0.7 秒待つ模擬。全体 0.3 秒・猶予 0.05 秒・対ごと 0.2 秒
clean; t0=$(date +%s%N)
out=$(MOCK=ok MOCK_CHECK=slow $R --bin "$T/mock" --target "$TGT" --args "0" \
      --logdir "$T/log" --run-id s12 --deadline 0.3 --grace 0.05 --per-pair 0.2 2>&1); rc=$?
t1=$(date +%s%N); ms=$(( (t1-t0)/1000000 ))
chk "終了値（成功にしない）" "$rc" 1
chk "済み 0" "$(echo "$out" | grep -c '済み 0')" 1
chk "未検証を成功に数えないと言う" "$(echo "$out" | grep -c '未検証を成功に数えません')" 1
chk "全体が 1.5 秒以内に戻る（無期限待機でない）" "$([ "$ms" -le 1500 ] && echo はい || echo いいえ)" はい
ok "実測 ${ms} ミリ秒"

echo
echo "## 13. 照合が止まったまま（KILL まで）"
clean; t0=$(date +%s%N)
out=$(MOCK=ok MOCK_CHECK=hang $R --bin "$T/mock" --target "$TGT" --args "0" \
      --logdir "$T/log" --run-id s13 --deadline 2 --grace 0.3 --per-pair 0.5 2>&1); rc=$?
t1=$(date +%s%N); ms=$(( (t1-t0)/1000000 ))
chk "終了値" "$rc" 1
chk "照合が打ち切られたと言う" "$(echo "$out" | grep -c '照合が期限で打ち切られました')" 1
chk "以降を起動しない" "$(echo "$out" | grep -c '停止条件に当たったので起動しません')" 2
chk "全体が 3 秒以内に戻る" "$([ "$ms" -le 3000 ] && echo はい || echo いいえ)" はい
ok "実測 ${ms} ミリ秒"

echo
echo "## 14. 事後の照合に予算が残らない"
# 本計算が対ごとの期限をすべて使い切ると、事後の照合に予算が残りません
clean; out=$(MOCK=slow MOCK_CHECK=auto $R --bin "$T/mock" --target "$TGT" --args "0" \
      --logdir "$T/log" --run-id s14 --deadline 0.9 --grace 0.3 --per-pair 0.5 2>&1); rc=$?
chk "終了値" "$rc" 1
chk "済みにしない" "$(echo "$out" | grep -c '済み 0')" 1
# **照合まで届けば「不通過」、予算が尽きていれば「未検証を成功に数えません」。**
# **どちらでも、済みには数えません。**
n=$(echo "$out" | grep -cE '照合 \*\*不通過\*\*|未検証を成功に数えません')
chk "未検証・不通過のどちらかで止まる" "$([ "$n" -ge 1 ] && echo はい || echo いいえ)" はい
chk "以降を起動しない" "$(echo "$out" | grep -c '停止条件に当たったので起動しません')" 2

echo
printf '**OK %d / NG %d**\n' "$OK" "$NG"
[ "$NG" -eq 0 ]

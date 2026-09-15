#!/usr/bin/env bash
# GMP 整合性診断の【駆動そのもの】の回帰試験（DESIGN-phase5-vertex-level.md §22.25）
#
# **合成の入力だけを使います。** 実データを指定して「拒否されるはず」を試す形は採りません
# （拒否されなければ実データ計算が走ってしまうため）。
#
# 使い方: bash tests/tools/gmp_diag_selftest.sh [駆動] [変異つき駆動]
set -u
BIN="${1:-build/tests/thingi_cp1_gmp}"
MUT="${2:-build/tests/thingi_cp1_gmp_mut}"
ARGS="0 6 2 0 0 1 1 0 0 16 2 1 0 1 1"   # 索引の突き合わせは 0（診断は 1 巡）
OK=0; NG=0
ok()   { OK=$((OK+1)); printf '  OK   %s\n' "$1"; }
ng()   { NG=$((NG+1)); printf '  **NG** %s\n' "$1"; }
chk()  { if [ "$2" = "$3" ]; then ok "$1（$2）"; else ng "$1: 期待 $3、実測 $2"; fi; }

[ -x "$BIN" ] || { echo "駆動がありません: $BIN"; exit 9; }
[ -x "$MUT" ] || { echo "変異つき駆動がありません: $MUT"; exit 9; }

# 合成の場を作ります（毎回新しく）
setup() {
    T=$(mktemp -d); mkdir -p "$T/kmesh"
    python3 tests/tools/mk_synth_kmesh.py "$T/kmesh" > "$T/synth.txt" || return 1
    printf '1001x1002 12 12\n' > "$T/synth_pairs.txt"
    printf '1001x1002\n' > "$T/synth_gmp_only.txt"
}
run()  { timeout 120 "$1" "$T/synth.txt" $ARGS > "$T/run.log" 2>&1; echo $?; }

echo "# GMP 整合性診断 — 駆動の回帰試験"
echo

echo "## 1. 正常な経路"
setup || exit 9
rc=$(run "$BIN")
chk "終了値" "$rc" 0
chk "結果の行数" "$(wc -l < "$T/synth_gmp_results.txt" 2>/dev/null || echo 0)" 1
chk "固定部の列数（34 以上）" "$(awk '{print (NF>=34)?"はい":"いいえ"}' "$T/synth_gmp_results.txt")" はい
chk "実施の列" "$(awk '{print $3}' "$T/synth_gmp_results.txt")" 1
chk "式 1" "$(awk '{print $4}' "$T/synth_gmp_results.txt")" id1=ok
chk "状態" "$(awk '{print $2}' "$T/synth_gmp_results.txt")" ok
chk "meta が作られた" "$([ -s "$T/synth_gmp.meta" ] && echo はい || echo いいえ)" はい
chk "対の一覧を上書きしない" "$(cat "$T/synth_pairs.txt")" "1001x1002 12 12"
SAVED="$T"

echo
echo "## 2. 再開（済みの対を飛ばす）"
rc=$(run "$BIN")
chk "終了値" "$rc" 0
chk "行が増えない" "$(wc -l < "$T/synth_gmp_results.txt")" 1
chk "これから 0 対" "$(grep -c 'これから 0 対を回します' "$T/run.log")" 1

echo
echo "## 3. 式 1 の不一致（変異）"
setup || exit 9
rc=$(KRI_DIAG_MUTATE=id1 timeout 120 "$MUT" "$T/synth.txt" $ARGS > "$T/run.log" 2>&1; echo $?)
chk "終了値（失敗が届く）" "$rc" 1
chk "状態" "$(awk '{print $2}' "$T/synth_gmp_results.txt")" FAIL
chk "式 1 の列" "$(awk '{print $4}' "$T/synth_gmp_results.txt")" id1=ng
chk "集計に届く" "$(grep -c '失敗\*\* 1' "$T/run.log")" 1

echo
echo "## 4. 4 演算目の位相が不良（変異）"
setup || exit 9
rc=$(KRI_DIAG_MUTATE=topo4 timeout 120 "$MUT" "$T/synth.txt" $ARGS > "$T/run.log" 2>&1; echo $?)
chk "終了値" "$rc" 1
chk "状態" "$(awk '{print $2}' "$T/synth_gmp_results.txt")" FAIL
chk "理由に位相" "$(grep -c '4演算の位相が不良' "$T/run.log")" 1

echo
echo "## 5. 途中で切れた行"
setup || exit 9
run "$BIN" > /dev/null
cut -d' ' -f1-4 "$T/synth_gmp_results.txt" > "$T/cut.txt" && mv "$T/cut.txt" "$T/synth_gmp_results.txt"
rc=$(run "$BIN")
chk "終了値（止まる）" "$rc" 2
chk "再開に使えないと言う" "$(grep -c '再開に使えません' "$T/run.log")" 1
chk "計算に入らない" "$(grep -c '開始 1001x1002' "$T/run.log")" 0

echo
echo "## 6. meta が違う"
setup || exit 9
run "$BIN" > /dev/null
rc=$(timeout 120 "$BIN" "$T/synth.txt" 0 6 4 0 0 1 1 0 0 16 2 1 0 1 1 > "$T/run.log" 2>&1; echo $?)
chk "終了値" "$rc" 2
chk "meta の不一致と言う" "$(grep -c 'meta が一致しません' "$T/run.log")" 1

echo
echo "## 7a. 結果があるのに meta が無い（出所不明）"
setup || exit 9
mkdir "$T/synth_gmp_results.txt"      # ★ 「在る」が meta は無い
rc=$(run "$BIN")
chk "終了値" "$rc" 2
chk "出所不明と言う" "$(grep -c '結果があるのに meta がありません' "$T/run.log")" 1
chk "量子化に入らない" "$(grep -c '量子化 ' "$T/run.log")" 0
rmdir "$T/synth_gmp_results.txt"

echo
echo "## 7b. 結果を書けない（meta はある）"
setup || exit 9
run "$BIN" > /dev/null                # meta と結果を作る
rm "$T/synth_gmp_results.txt"         # 結果だけ消す
chmod a-w "$T"                        # ★ 作れなくする
rc=$(run "$BIN")
chmod u+w "$T"
chk "終了値" "$rc" 2
chk "開けないと言う" "$(grep -c '結果ファイルを開けません' "$T/run.log")" 1

echo
echo "## 8. 不正な引数・設定の矛盾"
setup || exit 9
rc=$(timeout 60 "$BIN" "$T/synth.txt" 0 6 2 0 0 1 1 0 0 16 2 1 0 1 x > "$T/run.log" 2>&1; echo $?)
chk "第 16 引数が不正 → 終了値" "$rc" 2
chk "量子化に入らない" "$(grep -c '量子化' "$T/run.log")" 0
rc=$(timeout 60 "$BIN" "$T/synth.txt" 0 6 2 1 0 1 1 0 0 16 2 1 0 1 1 > "$T/run.log" 2>&1; echo $?)
chk "診断と索引の突き合わせ → 終了値" "$rc" 2
chk "量子化に入らない" "$(grep -c '量子化' "$T/run.log")" 0

echo
echo "## 9. 対象が足りない"
setup || exit 9
: > "$T/synth_gmp_only.txt"
rc=$(run "$BIN"); chk "空の一覧 → 終了値" "$rc" 2
setup || exit 9
printf '9999x8888\n' > "$T/synth_gmp_only.txt"
rc=$(run "$BIN"); chk "対の一覧に無い → 終了値" "$rc" 2
setup || exit 9
rm "$T/synth_pairs.txt"
rc=$(run "$BIN"); chk "対の一覧が無い → 終了値" "$rc" 2

echo
echo "## 10. 対の一覧と量子化が食い違う"
setup || exit 9
printf '1001x1002 99 12\n' > "$T/synth_pairs.txt"
rc=$(run "$BIN")
chk "終了値" "$rc" 2
chk "食い違いを言う" "$(grep -c '対の一覧と量子化が一致しません' "$T/run.log")" 1
chk "ブール演算に入らない" "$(grep -c '開始 1001x1002' "$T/run.log")" 0

echo
printf '**OK %d / NG %d**\n' "$OK" "$NG"
[ "$NG" -eq 0 ]

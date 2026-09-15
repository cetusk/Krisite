#!/usr/bin/env bash
# 投入の層と実駆動を【つないだ】回帰試験（DESIGN-phase5-vertex-level.md §23.12）
#
# **1 対の試験と、meta を扱わない模擬バイナリの試験では、
# 「同じ CP の 2 対目が meta 不一致で止まる」形を検出できませんでした。**
#
# **合成の入力だけを使います。**
set -u
BIN="${1:-build/tests/thingi_cp1_gmp}"
MUT="${2:-build/tests/thingi_cp1_gmp_mut}"
R="python3 tests/tools/run_gmp_diag.py"
ARGS="0 6 2 0 0 1 1 0 0 16 2 1 0 1 1"
OK=0; NG=0
ok() { OK=$((OK+1)); printf '  OK   %s\n' "$1"; }
ng() { NG=$((NG+1)); printf '  **NG** %s\n' "$1"; }
chk() { if [ "$2" = "$3" ]; then ok "$1（$2）"; else ng "$1: 期待 $3、実測 $2"; fi; }

# 3 対ぶんの合成（4 模型）。CP を 2 つ作ります
setup() {
    T=$(mktemp -d); mkdir -p "$T/kmesh" "$T/log"
    python3 - "$T/kmesh" <<'PY'
import sys
sys.path.insert(0, "tests/tools")
from mk_synth_kmesh import cube, write
out = sys.argv[1]
for i, (o, s) in enumerate([((0,0,0),4), ((2,0,0),4), ((0,0,0),6), ((3,1,0),5)]):
    write(f"{out}/{2001+i}.kmesh", *cube(*o, s))
PY
    printf '2001 12\n2002 12\n' > "$T/cpA.txt"
    printf '2003 12\n2004 12\n' > "$T/cpB.txt"
    # 対の一覧（駆動が読むだけ。三角形数は量子化後の値＝12）
    printf '2001x2002 12 12\n' > "$T/cpA_pairs.txt"
    printf '2003x2004 12 12\n' > "$T/cpB_pairs.txt"
    printf '2001x2002\n' > "$T/cpA_gmp_only.txt"
    printf '2003x2004\n' > "$T/cpB_gmp_only.txt"
}

echo "# 投入の層と実駆動をつないだ回帰試験"
echo
echo "## 1. 同じ CP で 2 対（meta が固定できること）"
setup
# CP A に 2 対目を足します（模型 2003 を CP A の一覧に加える）
printf '2001 12\n2002 12\n2003 12\n' > "$T/cpA.txt"
printf '2001x2002 12 12\n2001x2003 12 12\n' > "$T/cpA_pairs.txt"
printf '2001x2002\n2001x2003\n' > "$T/cpA_gmp_only.txt"
out=$($R --bin "$BIN" --target "$T/cpA.txt:$T/cpA_gmp_only.txt" --args "$ARGS" \
        --logdir "$T/log" --run-id t1 --deadline 300 --per-pair 120 2>&1); rc=$?
echo "$out" | sed 's/^/    /' | tail -6
chk "終了値" "$rc" 0
chk "2 対とも済み" "$(echo "$out" | grep -c '済み 2 / 失敗・打ち切り 0 / 未起動 0')" 1
chk "結果が 2 行" "$(wc -l < "$T/cpA_gmp_results.txt")" 2
chk "meta は 1 つ" "$(ls "$T"/cpA_gmp.meta | wc -l)" 1

echo
echo "## 2. 再開（検証した再利用。起動しない）"
out=$($R --bin "$BIN" --target "$T/cpA.txt:$T/cpA_gmp_only.txt" --args "$ARGS" \
        --logdir "$T/log" --run-id t2 --deadline 300 --per-pair 120 2>&1); rc=$?
chk "終了値" "$rc" 0
chk "2 対とも再利用" "$(echo "$out" | grep -c '検証した再利用')" 2
chk "行が増えない" "$(wc -l < "$T/cpA_gmp_results.txt")" 2

echo
echo "## 3. CP を 2 つまたいで、期限は 1 つ"
setup
out=$($R --bin "$BIN" --target "$T/cpA.txt:$T/cpA_gmp_only.txt" \
        --target "$T/cpB.txt:$T/cpB_gmp_only.txt" --args "$ARGS" \
        --logdir "$T/log" --run-id t3 --deadline 300 --per-pair 120 2>&1); rc=$?
chk "終了値" "$rc" 0
chk "2 対（CP をまたぐ）" "$(echo "$out" | grep -c 'これから 2 対を回します')" 1
chk "CP A の結果" "$(wc -l < "$T/cpA_gmp_results.txt")" 1
chk "CP B の結果" "$(wc -l < "$T/cpB_gmp_results.txt")" 1
chk "meta が CP ごとに" "$(ls "$T"/cpA_gmp.meta "$T"/cpB_gmp.meta | wc -l)" 2

echo
echo "## 4. unresolved だけを正にする変異"
setup
out=$(KRI_DIAG_MUTATE=unres4 $R --bin "$MUT" --target "$T/cpA.txt:$T/cpA_gmp_only.txt" \
        --args "$ARGS" --logdir "$T/log" --run-id t4 --deadline 300 --per-pair 120 2>&1); rc=$?
chk "終了値（失敗が届く）" "$rc" 1
chk "状態" "$(awk '{print $2}' "$T/cpA_gmp_results.txt")" FAIL
chk "理由に unresolved" "$(grep -c 'unresolvedが残る' "$T"/log/*t4*.log)" 1
chk "行は不完全と判定" "$(echo "$out" | grep -c '行 \*\*不完全\*\*')" 1

echo
echo "## 5. 計画に無い対は回さない"
setup
printf '2001x2002\n' > "$T/cpA_gmp_only.txt"
printf '9999x8888\n' > "$T/other.txt"
out=$($R --bin "$BIN" --target "$T/cpA.txt:$T/other.txt" --args "$ARGS" \
        --logdir "$T/log" --run-id t5 --deadline 120 --per-pair 60 2>&1); rc=$?
chk "終了値" "$rc" 1
chk "対の一覧に無いと言う" "$(grep -c '対の一覧に無い対です' "$T"/log/*t5*.log)" 1

echo
echo "## 6. ピーク RSS が【その子だけ】の値であること"
setup
# 2 対目は 1 対目より小さい模型。RUSAGE_CHILDREN なら 1 対目に引きずられます
printf '2001 12\n2002 12\n2003 12\n' > "$T/cpA.txt"
printf '2001x2002 12 12\n2001x2003 12 12\n' > "$T/cpA_pairs.txt"
printf '2001x2002\n2001x2003\n' > "$T/cpA_gmp_only.txt"
out=$($R --bin "$BIN" --target "$T/cpA.txt:$T/cpA_gmp_only.txt" --args "$ARGS" \
        --logdir "$T/log" --run-id t6 --deadline 300 --per-pair 120 2>&1)
r1=$(echo "$out" | grep '2001x2002:' | sed 's/.*ピーク RSS \([0-9]*\) MiB.*/\1/')
r2=$(echo "$out" | grep '2001x2003:' | sed 's/.*ピーク RSS \([0-9]*\) MiB.*/\1/')
ok "1 対目 ${r1} MiB / 2 対目 ${r2} MiB（**同じ子の値でないこと**を見ます）"
chk "2 つの値が独立に出ている" "$([ -n "$r1" ] && [ -n "$r2" ] && echo はい || echo いいえ)" はい

echo
printf '**OK %d / NG %d**\n' "$OK" "$NG"
[ "$NG" -eq 0 ]

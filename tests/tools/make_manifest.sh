#!/bin/bash
# 投入の基準（manifest）を作る（`DESIGN-phase5-vertex-level.md` §14.2）。
#
#   bash tests/tools/make_manifest.sh <バイナリ> [出力先] [基準の並び...]
#
# **ファイルを移動しません。** 候補は `KRI_BIN` で名指しして投入します。
# **旧バイナリと候補をそのまま保存でき、移動後のパスの食い違いも起きません。**
#
# **一時ファイルに作り、検査してから確定します。**
# **途中で失敗したら、そこで止まります**（`|| true` で握り潰しません）。
# **既存の manifest は上書きしません**（`--force` を付けたときだけ置き換えます）。
set -u

fail() { printf '**中止**: %s\n' "$1" >&2; exit 2; }

FORCE=0
if [ "${1:-}" = "--force" ]; then FORCE=1; shift; fi
BIN="${1:?投入するバイナリを渡してください}"
OUT="${2:-tests/tools/run_cp23.manifest}"
shift 2 2>/dev/null || shift $# 
BASES="${*:-cp2b cp3}"

[ -x "$BIN" ] || fail "実行ファイルがありません: $BIN"
[ "$FORCE" = 1 ] || [ ! -e "$OUT" ] || fail "既に存在します（置き換えるなら --force）: $OUT"

# **列数の取得は、終了状態を先に見ます**
COLERR="$(mktemp)"; TMP="$(mktemp)"
trap 'rm -f "$COLERR" "$TMP"' EXIT
cols="$("$BIN" --cols 2>"$COLERR")"
rc=$?
if [ "$rc" != 0 ]; then
    printf '  --cols の stderr:\n' >&2; head -5 "$COLERR" >&2
    fail "--cols が異常終了しました（exit ${rc}）: $BIN"
fi
case "${cols:-x}" in ''|*[!0-9]*) fail "--cols が数を返しません（${cols}）: $BIN" ;; esac
[ "$cols" -ge 7 ] && [ "$cols" -le 4096 ] || fail "--cols が範囲の外です（${cols}）"

bin_sha="$(sha256sum "$BIN" | cut -d' ' -f1)" || fail "指紋を作れません: $BIN"
printf 'bin=%s\n' "$bin_sha" >> "$TMP" || fail "一時ファイルに書けません"
for b in $BASES; do
    only="data/thingi10k/${b}_only.txt"
    [ -f "$only" ] || fail "標本の一覧がありません: $only"
    n="$(grep -c . "$only")" || fail "件数を数えられません: $only"
    [ "$n" -gt 0 ] || fail "標本の一覧が空です: $only"
    sort "$only" | uniq -d | grep -q . && fail "標本にキーの重複があります: $only"
    keys="$(sort "$only" | sha256sum | cut -d' ' -f1)" || fail "ハッシュを作れません: $only"
    printf '%s keys=%s n=%s cols=%s\n' "$b" "$keys" "$n" "$cols" >> "$TMP" \
        || fail "一時ファイルに書けません"
done

# **確定の前に、作ったものを読み直して検査します**
grep -q '^bin=[0-9a-f]\{64\}$' "$TMP" || fail "作った基準の bin= が不正です"
for b in $BASES; do
    line="$(awk -v b="$b" '$1==b{print}' "$TMP")"
    [ -n "$line" ] || fail "作った基準に ${b} の行がありません"
    echo "$line" | grep -qE "^${b} keys=[0-9a-f]{64} n=[0-9]+ cols=[0-9]+$" \
        || fail "作った基準の ${b} の行が不正です: $line"
done

mv "$TMP" "$OUT" || fail "確定できません: $OUT"
trap 'rm -f "$COLERR"' EXIT
chmod 0644 "$OUT"
printf '作りました: %s\n' "$OUT"
cat "$OUT"

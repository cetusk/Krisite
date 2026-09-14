# -*- coding: utf-8 -*-
"""CP3 の層化標本を選ぶ。

**規模で 5 層に切り、層の中では【性質の被覆】を見て選びます。**

規模で切る理由: 時間の見積もりが立ち、CP2 と同じ軸で読めるため。
性質を見る理由: CP3 を回す目的が「CP2 に無い性質の模型で新しい失敗が出るか」なので、
**層の中で偏ると目的を外します**（`CLAUDE.md`「除外や例外は識別子ではなく性質で」は
検査対象の選び方にも当てはまります）。

**等間隔だけで選ぶと、性質が偏っていても気づけません。**
"""
import sys, collections

def load_props(path):
    p = {}
    for ln in open(path, encoding='utf-8'):
        if ln.startswith('#') or ln.startswith('id'):
            continue
        f = ln.split()
        if len(f) < 13:
            continue
        p[f[0]] = dict(em=int(f[6]), vm=int(f[7]), comp=int(f[8]), si=int(f[12]),
                       drop=int(f[3]), merge=int(f[4]))
    return p

def cls_of(props, a, b):
    """対の性質のクラス。**両端の性質を合わせて 1 つの符号にします。**"""
    def one(k):
        d = props.get(k)
        if d is None:
            return 'unknown'
        return ('nm' if (d['em'] == 0 or d['vm'] == 0) else 'mf') + \
               ('+si' if d['si'] else '') + ('+mc' if d['comp'] > 1 else '')
    return tuple(sorted((one(a), one(b))))

def main(pairs_path, quant_path, out_path, per=59, strata=5):
    props = load_props(quant_path)
    rows = []
    for ln in open(pairs_path, encoding='utf-8'):
        f = ln.split()
        a, b = f[0].split('x')
        rows.append((int(f[1]) + int(f[2]), f[0], cls_of(props, a, b)))
    rows.sort()
    n = len(rows)
    print(f"母集団 {n} 対 / 性質のクラス {len(set(c for _,_,c in rows))} 種")
    all_cnt = collections.Counter(c for _, _, c in rows)
    picked = []
    for i in range(strata):
        seg = rows[round(i*n/strata):round((i+1)*n/strata)]
        # **まず各クラスから 1 対ずつ**（被覆）。残りを等間隔で埋める
        seen, first = set(), []
        for r in seg:
            if r[2] not in seen:
                seen.add(r[2])
                first.append(r)
        take = first[:per]
        if len(take) < per:
            rest = [r for r in seg if r not in take]
            step = max(1, len(rest) / (per - len(take)))
            take += [rest[min(int(round(j*step)), len(rest)-1)] for j in range(per - len(take))]
        # 重複を除いて per 件にそろえる
        uniq, keys = [], set()
        for r in take:
            if r[1] not in keys:
                keys.add(r[1]); uniq.append(r)
        j = 0
        while len(uniq) < per and j < len(seg):
            if seg[j][1] not in keys:
                keys.add(seg[j][1]); uniq.append(seg[j])
            j += 1
        picked.append(sorted(uniq))
        cov = collections.Counter(c for _, _, c in uniq)
        seg_cls = collections.Counter(c for _, _, c in seg)
        print(f"  層 {i+1}: {seg[0][0]:,}–{seg[-1][0]:,} / 母集団 {len(seg)} 対・{len(seg_cls)} クラス"
              f" → 標本 {len(uniq)} 対・**{len(cov)} クラス被覆**")
    flat = [r for s in picked for r in s]
    with open(out_path, 'w', encoding='utf-8') as f:
        for _, key, _ in flat:
            f.write(key + "\n")
    got = collections.Counter(c for _, _, c in flat)
    print(f"\n選んだ {len(flat)} 対 / 被覆したクラス {len(got)} / 全クラス {len(all_cnt)}")
    miss = [c for c in all_cnt if c not in got]
    if miss:
        print(f"**被覆できなかったクラス {len(miss)} 種**（母集団での件数の合計 "
              f"{sum(all_cnt[c] for c in miss)} 対）")
        for c in sorted(miss, key=lambda c: -all_cnt[c])[:5]:
            print(f"    {c} : {all_cnt[c]} 対")
    else:
        print("**すべてのクラスを被覆しました**")

main(*sys.argv[1:4], per=int(sys.argv[4]) if len(sys.argv) > 4 else 59)

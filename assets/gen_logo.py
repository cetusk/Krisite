import re, os
D = os.path.dirname(os.path.abspath(__file__))
mark = open(os.path.join(D,'krisite-mark.svg')).read()
CRYSTAL = re.search(r'<g>\n(.*?)\n  </g>', mark, re.S).group(1)

FS, LS, X0 = 62, 3, 205
W = {'K':722,'r':389,'\u0131':278,'s':556,'t':333,'e':556}
word = "Kr\u0131s\u0131te"
xs, x = [], X0
for ch in word:
    xs.append(round(x,1)); x += W[ch]*FS/1000 + LS
XSTR = " ".join(str(v) for v in xs)

BASE = 100.0
CAP  = BASE - 0.717*FS
i1 = xs[2] + W['\u0131']*FS/2000
i2 = xs[4] + W['\u0131']*FS/2000

# --- K の縦画（|）の実寸に合わせる ---
LSB   = 0.076*FS          # Helvetica Bold 'K' の左サイドベアリング
STEMW = 0.142*FS          # 縦画の幅
SL    = xs[0] + LSB       # 縦画 左端
SR    = SL + STEMW        # 縦画 右端
XL, XR = SL-2.0, SR+2.0   # 余剰は2重クリップで落とす

SLOPE = 0.95
A, B, C, TOP = 66.0, 75.0, 84.0, 38.0     # x=SL 基準
def y_at(Y, x): return Y + (SL-x)*SLOPE
def band(y_top, y_bot):
    return (f"{XL:.1f},{y_at(y_top,XL):.1f} {XR:.1f},{y_at(y_top,XR):.1f} "
            f"{XR:.1f},{y_at(y_bot,XR):.1f} {XL:.1f},{y_at(y_bot,XL):.1f}")

LEG, YBOT = 10.0, 62.5     # i の三角パッチ（下げた）
def tri(c):
    x0 = c-LEG/2
    return f"{x0:.1f},{YBOT:.1f} {x0:.1f},{YBOT-LEG:.1f} {x0+LEG:.1f},{YBOT:.1f}"

FONT = "'Helvetica Neue', Helvetica, Arial, sans-serif"

# --- タグライン：幅を実測して viewBox を自動調整 ---
TAG = "exact, plane-based geometry for point clouds and meshes"
TFS, TLS = 15.5, 1.5
TW = {' ':278,',':278,'-':333,'a':556,'b':556,'c':500,'d':556,'e':556,'f':278,'g':556,
      'h':556,'i':222,'l':222,'m':833,'n':556,'o':556,'p':556,'r':333,'s':500,'t':278,
      'u':556,'x':500,'y':500}
TX  = xs[0] + 3
TWID = sum(TW[c]*TFS/1000 + TLS for c in TAG)
VIEW_W = max(TX + TWID, xs[-1] + W[word[-1]]*FS/1000) + 22

def build(word_fill, tag_fill):
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {VIEW_W:.0f} 200" width="{VIEW_W:.0f}" height="200" role="img" aria-label="Krisite">
  <title>Krisite</title>
  <style>
    .w {{ font-family: {FONT}; font-size: {FS}px; font-weight: 700; }}
  </style>

  <clipPath id="k-clip">
    <text class="w" x="{xs[0]}" y="{BASE}">K</text>
  </clipPath>
  <!-- 縦画の帯だけを通す。腕の付け根へのはみ出しをここで断つ -->
  <clipPath id="stem-clip">
    <rect x="{xs[0]-4:.1f}" y="0" width="{SR-xs[0]+4:.1f}" height="200"/>
  </clipPath>

  <g>
{CRYSTAL}
  </g>

  <text class="w" fill="{word_fill}" x="{XSTR}" y="{BASE}">{word}</text>

  <!-- K の左上先端：縦画の幅ちょうどに収めた3段の断片。境界は腕の「/」と平行 -->
  <g clip-path="url(#k-clip)">
    <g clip-path="url(#stem-clip)">
      <polygon points="{band(A, TOP)}" fill="#5EEAD4"/>
      <polygon points="{band(B, A)}"   fill="#2DD4BF"/>
      <polygon points="{band(C, B)}"   fill="#0D9488"/>
    </g>
  </g>

  <polygon points="{tri(i1)}" fill="#2DD4BF"/>
  <polygon points="{tri(i2)}" fill="#2DD4BF"/>

  <text fill="{tag_fill}" x="{TX:.0f}" y="128"
        font-family="{FONT}" font-size="{TFS}" letter-spacing="{TLS}">{TAG}</text>
</svg>
"""

open(os.path.join(D,'krisite-logo.svg'),'w').write(build("#2E3741", "#6E7A87"))
open(os.path.join(D,'krisite-logo-dark.svg'),'w').write(build("#E9EEF3", "#8C99A6"))

import xml.etree.ElementTree as ET
for p in [os.path.join(D,f) for f in ['krisite-logo.svg','krisite-logo-dark.svg','krisite-mark.svg']]:
    ET.parse(p); print(p.split('/')[-1], "XML OK")
print("K 縦画: x = %.1f 〜 %.1f (幅 %.1f)" % (SL, SR, STEMW))
print("断片 塗り範囲: x = %.1f 〜 %.1f" % (XL, XR))
print("三角パッチ: 一辺 %.0f, 下端 y=%.1f (x-ハイト上端 y=%.1f)" % (LEG, YBOT, BASE-0.52*FS))
print("タグライン: x=%.0f 幅%.0f 右端%.0f -> viewBox 幅 %.0f" % (TX, TWID, TX+TWID, VIEW_W))

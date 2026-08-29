// Krisite — 凸多角形スープ（中間表現）
//
// SPEC-phase3.md §1「確定した契約」
//
//     from_mesh : TriMesh(整数) → PolySoup      入口
//     boolean   : PolySoup × … → PolySoup       ★ CSG について閉じる
//     to_mesh   : PolySoup → TriMesh            出口
//
// **中間表現を平面ベースの凸多角形に保つので、連鎖で丸めが入りません。**
// $(A \cup B) \setminus C$ を計算するのに中間結果を整数へ落とす必要がなくなります。
//
// ---
//
// ## 分類の台は【生成 0 の整数メッシュ】です
//
// スープは、自分が**どの入力メッシュから作られたか**（`sources`）と、
// **各入力の内外からどう決まるか**（`indicator`）を持ちます。
//
// 断片の分類は「相手のスープの内側か」ではなく「各入力メッシュの内側か」で行い、
// その真偽ベクトルに指示関数を適用します。**入力メッシュは整数座標なので、
// 何段連鎖してもレイキャストが厳密に行えます。**
//
// > **これは WNV（`SPEC-phase3.md` §5.1）の真偽値版です。** WNV は $\mathbb{Z}^n$ で
// > 巻き数を数えますが、閉多様体入力なら内外の 1 ビットで足ります。**CP3 で
// > $\mathbb{Z}^n$ に一般化します。** 指示関数（§5.2）はそのまま使えます。
#ifndef KRISITE_CSG_POLYSOUP_HPP
#define KRISITE_CSG_POLYSOUP_HPP

#include <cstdint>
#include <vector>

#include "krisite/csg/fragment.hpp"
#include "krisite/csg/plane_table.hpp"
#include "krisite/geom/plane.hpp"
#include "krisite/mesh/tri_mesh.hpp"
#include "krisite/octree/uniform_grid.hpp"

namespace krisite::csg {

/// スープの 1 枚。**凸多角形**（支持平面 + 辺平面）と、由来と外接箱。
struct Poly {
    Fragment frag;          ///< 支持平面・辺平面・向き（`owner` は使わない）
    octree::Aabb aabb{};    ///< **保守的な**外接箱（セル割り当てに使う）
    std::uint32_t src = 0;  ///< どの入力メッシュの面か（`sources` の添字）
    std::uint32_t tag = 0;  ///< 由来のタグ（§4.3。元の多角形 ID）
};

/// 指示関数（`SPEC-phase3.md` §5.2）。**表ではなく【関数】として持ちます。**
///
/// > 真理値表を materialize しないこと。`sources` が $n$ 個なら $2^n$ 項になり、
/// > 工具を 20 個引くミリングでは 100 万項です。そして **CP3 で $\mathbb{Z}^n$ に
/// > 一般化すると、表そのものが作れません**（定義域が無限）。
///
/// 合成は関数合成、評価は呼び出しです。式木を平らな配列で持ち、**最後の要素が根**。
/// 子の添字は必ず自分より小さいので、前から 1 回走査すれば評価できます。
///
/// **定義域は「各 source の状態」です。** CP2 は内外の 0/1 ですが、CP3 で巻き数
/// （$\mathbb{Z}^n$）に変えても「**非零なら内側**」の規約は変わりません（§5.1）。
struct Indicator {
    enum class Kind : std::uint8_t { Source, Not, And, Or };
    struct Node {
        Kind kind = Kind::Source;
        std::uint32_t src = 0;  ///< Kind::Source のとき、sources の添字
        std::uint32_t a = 0, b = 0;
    };
    std::vector<Node> nodes;

    bool eval(const std::vector<std::int32_t>& w) const {
        KRISITE_CHECK(!nodes.empty(), "Indicator: 空の式");
        std::vector<char> v(nodes.size(), 0);
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            const Node& n = nodes[i];
            switch (n.kind) {
                case Kind::Source:
                    KRISITE_CHECK(n.src < w.size(), "Indicator: source の添字が範囲外");
                    v[i] = (w[n.src] != 0) ? 1 : 0;  // **非零なら内側**（§5.1）
                    break;
                case Kind::Not:
                    v[i] = v[n.a] ? 0 : 1;
                    break;
                case Kind::And:
                    v[i] = (v[n.a] && v[n.b]) ? 1 : 0;
                    break;
                case Kind::Or:
                    v[i] = (v[n.a] || v[n.b]) ? 1 : 0;
                    break;
            }
        }
        return v.back() != 0;
    }
};

/// 「source `i` の内側」だけを見る指示関数。
inline Indicator indicator_source(std::uint32_t src) {
    Indicator f;
    f.nodes.push_back({Indicator::Kind::Source, src, 0, 0});
    return f;
}

/// 演算の種類（`boolean.hpp` の `BoolOp` に依存しないよう、ここで持つ）。
enum class Compose { Union, Intersection, Difference };

/// **合成は関数合成です。** `y` の source 添字と子の添字をずらして連結し、根を足すだけ。
inline Indicator compose(const Indicator& x, const Indicator& y, std::uint32_t src_offset,
                         Compose how) {
    Indicator f = x;
    const auto base = static_cast<std::uint32_t>(f.nodes.size());
    const std::uint32_t root_x = base - 1;
    for (const Indicator::Node& n : y.nodes) {
        Indicator::Node m = n;
        if (m.kind == Indicator::Kind::Source) {
            m.src += src_offset;
        } else {
            m.a += base;
            if (m.kind == Indicator::Kind::And || m.kind == Indicator::Kind::Or) m.b += base;
        }
        f.nodes.push_back(m);
    }
    const auto root_y = static_cast<std::uint32_t>(f.nodes.size() - 1);
    switch (how) {
        case Compose::Union:
            f.nodes.push_back({Indicator::Kind::Or, 0, root_x, root_y});
            break;
        case Compose::Intersection:
            f.nodes.push_back({Indicator::Kind::And, 0, root_x, root_y});
            break;
        case Compose::Difference: {
            f.nodes.push_back({Indicator::Kind::Not, 0, root_y, 0});
            const auto not_y = static_cast<std::uint32_t>(f.nodes.size() - 1);
            f.nodes.push_back({Indicator::Kind::And, 0, root_x, not_y});
            break;
        }
    }
    return f;
}

/// 凸多角形スープ。**CSG について閉じた中間表現。**
struct PolySoup {
    PlaneTable table;
    std::vector<Poly> polys;
    /// **生成 0 の整数メッシュ。** 分類の台であり、連鎖しても増えるだけで丸められない
    std::vector<mesh::TriMesh> sources;
    /// 指示関数（§5.2）。**表ではなく関数**
    Indicator indicator;

    std::size_t source_count() const noexcept { return sources.size(); }
};

/// ビット列を「各 source の状態」に開く。`own` の面上の断片について、
/// その source の内外を `own_inside` で与えます（表裏の 2 通りを作るため）。
inline std::vector<std::int32_t> to_membership(std::uint32_t bits, std::size_t n, std::uint32_t own,
                                               bool own_inside) {
    std::vector<std::int32_t> w(n, 0);
    for (std::size_t i = 0; i < n; ++i) w[i] = ((bits >> i) & 1u) ? 1 : 0;
    w[own] = own_inside ? 1 : 0;
    return w;
}

/// 入力メッシュ 1 枚を、そのまま「内側 = そのメッシュの内側」のスープにする。
///
/// **面併合は行いません**（`SPEC-phase3.md` §3.1.4 で「必須」から「任意の最適化」へ）。
/// 辺平面は §3.1 の構成（軸方向との外積）で作るので、隣接情報が要りません。
/// **PWN な polygon soup や自己交差する入力を受け入れる前提**がこれで整います。
/// `from_mesh` の設定（`SPEC-phase3.md` §4.1）。
struct FromMeshOptions {
    /// 併合の上限。**0 なら三角形化**（CP1〜CP5 の挙動 = §4.1.1 の正解器）。
    ///
    /// > 凸分割は「三角形化してから、凸性が保たれる限り貪欲に併合する」形なので、
    /// > **上限を 0 にすれば三角形化**になります。コードは一本化されています。
    ///
    /// **切り替えは実行時です**（§4.1.2）。コンパイル時定数だと同一プロセスで
    /// 比較できません。
    std::size_t max_merges = 0;
};

/// 上限なし（凸分割を最後まで行う）。
inline constexpr std::size_t kMergeAll = static_cast<std::size_t>(-1);

namespace detail {

/// 支持平面の法線の絶対値が最大の軸。**投影して 2D で扱うため**の軸です。
inline geom::Axis dominant_axis(const geom::PlaneD& pl) {
    const std::size_t bx = arith::min_bits(pl.a), by = arith::min_bits(pl.b),
                      bz = arith::min_bits(pl.c);
    // **幅で選びます。** 値そのものは固定幅整数なので、比較には `cmp_abs` が要ります。
    // 幅は同値のとき区別しませんが、**非零であればどれを選んでも投影は退化しません**
    // （法線の非零成分に対応する軸へ潰すと面積が保たれるため）。
    if (bx >= by && bx >= bz && !arith::is_zero(pl.a)) return geom::Axis::X;
    if (by >= bz && !arith::is_zero(pl.b)) return geom::Axis::Y;
    if (!arith::is_zero(pl.c)) return geom::Axis::Z;
    return arith::is_zero(pl.a) ? geom::Axis::Y : geom::Axis::X;
}

/// 頂点ループが凸か（曲がりの符号が一貫しているか）。**整数演算だけ**で決まります。
///
/// **入口では頂点が整数座標なので、通常の手法が使えます**（`SPEC-phase3.md` §4.1）。
/// 同次座標で三角形分割ができないのは中間表現の話です。
inline bool loop_is_convex(const mesh::TriMesh& m, const std::vector<std::uint32_t>& loop,
                           geom::Axis ax) {
    const std::size_t n = loop.size();
    if (n < 3) return false;
#if defined(KRISITE_MUTATION_MERGE_IGNORE_CONVEXITY)
    // 変異: **凸性を確かめない**（`SPEC-phase3.md` §4.1.2）。
    //
    // `Fragment` は半空間の交わりなので、非凸なループを渡すと**凸包**になります。
    // **位相は保たれることがあります**（$\chi$ も $C$ も変わらない）。
    // 「三角形化と凸分割で両方向の差が空」でしか出ません。
    (void)ax;
    return true;
#else
    int seen = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const int o = geom::orient2d(m.vertices[loop[(i + n - 1) % n]], m.vertices[loop[i]],
                                     m.vertices[loop[(i + 1) % n]], ax);
        if (o == 0) continue;  // 共線は凸性を壊しません
        if (seen == 0) {
            seen = o;
        } else if (seen != o) {
            return false;
        }
    }
    return seen != 0;
#endif
}

/// 共線の頂点を落とす。**落とさないと辺平面が重複します。**
///
/// 頂点 $i$ は `support ∩ edge[i-1] ∩ edge[i]` なので、連続する 2 辺が同じ平面に
/// なると `intersect3` が退化します。**共線の頂点は角ではないので落とせます。**
///
/// > 隣の多角形がそこに角を持っていれば T 頂点になりますが、
/// > **案 D の T 頂点解決（`SPEC-phase2.md` §2.4.3）が縫合後に埋めます。**
inline void drop_collinear(const mesh::TriMesh& m, std::vector<std::uint32_t>& loop,
                           geom::Axis ax) {
    bool changed = true;
    while (changed && loop.size() > 3) {
        changed = false;
        const std::size_t n = loop.size();
        for (std::size_t i = 0; i < n; ++i) {
            if (geom::orient2d(m.vertices[loop[(i + n - 1) % n]], m.vertices[loop[i]],
                               m.vertices[loop[(i + 1) % n]], ax) == 0) {
                loop.erase(loop.begin() + static_cast<std::ptrdiff_t>(i));
                changed = true;
                break;
            }
        }
    }
}

/// 有向辺 `(u, v)` がループのどこにあるか。無ければ `npos`。
inline std::size_t find_edge(const std::vector<std::uint32_t>& loop, std::uint32_t u,
                             std::uint32_t v) {
    const std::size_t n = loop.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (loop[i] == u && loop[(i + 1) % n] == v) return i;
    }
    return static_cast<std::size_t>(-1);
}

/// ループ `a` と `b` を共有辺で貼り合わせる。`a` は `(u,v)`、`b` は `(v,u)` を持つこと。
inline std::vector<std::uint32_t> splice(const std::vector<std::uint32_t>& a, std::size_t ia,
                                         const std::vector<std::uint32_t>& b, std::size_t ib) {
    const std::size_t na = a.size(), nb = b.size();
    std::vector<std::uint32_t> r;
    r.reserve(na + nb - 2);
    // `a` を $v$ から一周して $u$ まで（$v = a[i_a{+}1]$、$u = a[i_a]$）
    for (std::size_t k = 0; k < na; ++k) r.push_back(a[(ia + 1 + k) % na]);
    // `b` の $u$ と $v$ の**間**だけ（$u = b[i_b{+}1]$、$v = b[i_b]$）
    for (std::size_t k = 0; k + 2 < nb; ++k) r.push_back(b[(ib + 2 + k) % nb]);
    KRISITE_CHECK(r.size() == na + nb - 2, "splice: 頂点数が合わない");
    return r;
}

}  // namespace detail

inline PolySoup from_mesh(const mesh::TriMesh& m, const FromMeshOptions& opt = {}) {
    PolySoup s;
    s.sources.push_back(m);
    s.indicator = indicator_source(0);  // 「source 0 の内側」

    // ---- 1. 三角形を「支持平面 + 頂点ループ」にする ---------------------------
    //
    // **併合の上限が 0 ならここで終わり**、三角形化そのものになります（§4.1）。
    struct Piece {
        PlaneId support = kNoPlane;
        bool flipped = false;
        std::vector<std::uint32_t> loop;  ///< `m.vertices` への添字
        std::uint32_t tag = 0;            ///< **最小の三角形添字**（正準）
        bool alive = true;
    };
    std::vector<Piece> pieces;
    pieces.reserve(m.triangles.size());
    for (std::size_t ti = 0; ti < m.triangles.size(); ++ti) {
        const mesh::Tri& t = m.triangles[ti];
        const geom::PlaneD sp =
            geom::plane_from_triangle(m.vertices[t[0]], m.vertices[t[1]], m.vertices[t[2]]);
        if (geom::is_degenerate(sp)) continue;  // 面積 0 の三角形は入口で落とす
        const PlaneRef ref = s.table.intern(sp);
        Piece p;
        p.support = ref.id;
        p.flipped = ref.flipped;
        p.loop = {t[0], t[1], t[2]};
        p.tag = static_cast<std::uint32_t>(ti);
        pieces.push_back(std::move(p));
    }

    // ---- 2. 凸性が保たれる限り貪欲に併合する（§4.1）--------------------------
    //
    // **走査順を正準に固定します**（§4.1.2）。貪欲な併合は順序に依存するので、
    // 決めておかないと同じ入力で違う結果が出ます。ここでは
    //
    //   (1) 支持平面が同じで **向きも同じ**（`flipped` が一致）
    //   (2) 有向辺 $(u,v)$ と $(v,u)$ を共有する
    //   (3) 貼り合わせたループが**凸**
    //
    // を満たす組を、**添字の小さい順**に併合します。
    // **`flipped` を条件に入れるのを忘れないこと。** 裏表が逆の面を貼ると、
    // 巻き順の基準が混ざります。
    std::size_t merges = 0;
    if (opt.max_merges > 0) {
        bool changed = true;
        while (changed && merges < opt.max_merges) {
            changed = false;
            for (std::size_t i = 0; i < pieces.size() && merges < opt.max_merges; ++i) {
                if (!pieces[i].alive) continue;
                const geom::Axis ax = detail::dominant_axis(s.table.at(pieces[i].support));
                for (std::size_t j = i + 1; j < pieces.size(); ++j) {
                    if (!pieces[j].alive) continue;
                    if (pieces[j].support != pieces[i].support) continue;
                    if (pieces[j].flipped != pieces[i].flipped) continue;
                    // 共有する有向辺を探す
                    const std::size_t na = pieces[i].loop.size();
                    std::size_t ia = static_cast<std::size_t>(-1), ib = ia;
                    for (std::size_t k = 0; k < na; ++k) {
                        const std::uint32_t u = pieces[i].loop[k];
                        const std::uint32_t v = pieces[i].loop[(k + 1) % na];
                        const std::size_t f = detail::find_edge(pieces[j].loop, v, u);
                        if (f != static_cast<std::size_t>(-1)) {
                            ia = k;
                            ib = f;
                            break;
                        }
                    }
                    if (ia == static_cast<std::size_t>(-1)) continue;
                    std::vector<std::uint32_t> merged =
                        detail::splice(pieces[i].loop, ia, pieces[j].loop, ib);
                    if (!detail::loop_is_convex(m, merged, ax)) continue;
                    detail::drop_collinear(m, merged, ax);
                    if (merged.size() < 3 || !detail::loop_is_convex(m, merged, ax)) continue;
                    pieces[i].loop = std::move(merged);
                    pieces[i].tag = std::min(pieces[i].tag, pieces[j].tag);
                    pieces[j].alive = false;
                    ++merges;
                    changed = true;
                    break;
                }
            }
        }
    }

    // ---- 3. 断片に落とす ------------------------------------------------------
    s.polys.reserve(pieces.size());
    for (const Piece& p : pieces) {
        if (!p.alive) continue;
        // **値で受けること。** `intern` が平面表を伸ばすと参照が無効になります
        // （実際に踏みました。支持平面と平行な候補が棄却されなくなり、
        // 辺平面が支持平面と一致して `intersect3` が w=0 で落ちます）。
        const geom::PlaneD sp = s.table.at(p.support);
        Poly q;
        q.frag.support = p.support;
        q.frag.flipped = p.flipped;
        q.frag.owner = 0;
        q.src = 0;
        q.tag = p.tag;
        // 頂点 i = support ∩ edge[i-1] ∩ edge[i] なので、辺 (v_{i-1}, v_i) の平面を
        // edge[i] に置く（`fragment.hpp` の規約）。
        const std::size_t n = p.loop.size();
        q.frag.edge.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            const geom::IPoint& a = m.vertices[p.loop[(i + n - 1) % n]];
            const geom::IPoint& b = m.vertices[p.loop[i]];
            q.frag.edge.push_back(s.table.intern(geom::plane_from_edge(a, b, sp)).id);
        }
        for (int k = 0; k < 3; ++k) {
            q.aabb.lo[k] = krisite::kCoordMax;
            q.aabb.hi[k] = krisite::kCoordMin;
        }
        for (std::uint32_t vi : p.loop) {
            const geom::IPoint& v = m.vertices[vi];
            const std::int64_t c[3] = {v.x, v.y, v.z};
            for (int k = 0; k < 3; ++k) {
                q.aabb.lo[k] = std::min(q.aabb.lo[k], c[k]);
                q.aabb.hi[k] = std::max(q.aabb.hi[k], c[k]);
            }
        }
        s.polys.push_back(std::move(q));
    }
    return s;
}

}  // namespace krisite::csg

#endif  // KRISITE_CSG_POLYSOUP_HPP

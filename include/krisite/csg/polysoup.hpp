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
inline PolySoup from_mesh(const mesh::TriMesh& m) {
    PolySoup s;
    s.sources.push_back(m);
    s.indicator = indicator_source(0);  // 「source 0 の内側」
    s.polys.reserve(m.triangles.size());
    for (std::size_t ti = 0; ti < m.triangles.size(); ++ti) {
        const mesh::Tri& t = m.triangles[ti];
        const geom::IPoint& p0 = m.vertices[t[0]];
        const geom::IPoint& p1 = m.vertices[t[1]];
        const geom::IPoint& p2 = m.vertices[t[2]];
        const geom::PlaneD sp = geom::plane_from_triangle(p0, p1, p2);
        if (geom::is_degenerate(sp)) continue;  // 面積 0 の三角形は入口で落とす
        const PlaneRef ref = s.table.intern(sp);

        Poly q;
        q.frag.support = ref.id;
        q.frag.flipped = ref.flipped;
        q.frag.owner = 0;
        q.src = 0;
        q.tag = static_cast<std::uint32_t>(ti);
        // 頂点 i = support ∩ edge[i-1] ∩ edge[i] なので、辺 (p_{i-1}, p_i) の平面を
        // edge[i] に置く（`fragment.hpp` の規約）。
        const geom::IPoint* vp[3] = {&p0, &p1, &p2};
        for (int i = 0; i < 3; ++i) {
            const geom::IPoint& a = *vp[(i + 2) % 3];
            const geom::IPoint& b = *vp[i];
            q.frag.edge.push_back(s.table.intern(geom::plane_from_edge(a, b, sp)).id);
        }
        for (int k = 0; k < 3; ++k) {
            q.aabb.lo[k] = krisite::kCoordMax;
            q.aabb.hi[k] = krisite::kCoordMin;
        }
        for (const geom::IPoint* p : vp) {
            const std::int64_t c[3] = {p->x, p->y, p->z};
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

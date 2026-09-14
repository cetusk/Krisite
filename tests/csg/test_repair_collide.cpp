// Krisite — **修復の拒否分岐**（衝突・平面不足・組不足）
//
// `DESIGN-phase5-vertex-level.md` §12.4、`SPEC-phase5.md` §5.10.14.74 の 2。
//
// ## これは【人工の入力による単体試験】です
//
// **公開経路（ブール演算）から「既存頂点との衝突」「別の辺の細分点との衝突」に
// 到達させる幾何は構成していません。** ここで確かめるのは
// **「その状態に置かれたとき、狙った分岐が発火し、メッシュを書き換えないこと」**です。
//
// ## 不変性の範囲
//
// **「その候補の処理開始から拒否までのメッシュ内容」だけ**を主張します。
// **関数入口の `vertex_split_src.assign`・他候補の成功・診断カウンタは含みません。**
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "krisite/csg/to_mesh.hpp"

#include "test_util.hpp"

using namespace krisite;

namespace {

int failures = 0;

void expect(bool ok, const char* what) {
    if (!ok) {
        ++failures;
        std::printf("  **FAIL** %s\n", what);
    }
}

geom::IPoint P(std::int32_t x, std::int32_t y, std::int32_t z) {
    return geom::IPoint{x, y, z};
}

/// 辺 $(V,W)$ = z 軸の $[0,10]$ に 4 枚が接する、最小の配置。
///
/// **`repair_unresolved_edges` の前提をすべて満たします** —
/// 4 枚あり、組が 2 つあり、両端に載る平面 2 枚（平行でない）と
/// 片側だけに載る平面が 1 枚ずつある。
struct Fixture {
    csg::PolySoup soup;
    csg::SoupMesh mesh;
    csg::ToMeshStats st;
    csg::PlaneId p1{}, p2{}, p3{}, p4{};
    std::uint32_t a = 0, b = 0;    ///< 1 本目の辺
    std::uint32_t a2 = 0, b2 = 0;  ///< 2 本目の辺（同じ位置の複製）
    geom::HPointD mid{};           ///< 期待される細分点
};

/// `extra` が真なら、同じ位置の 2 本目の辺（複製された頂点）も作ります。
Fixture make(bool extra) {
    Fixture f;
    const geom::PlaneD pl1 = geom::plane_from_triangle(P(0, 0, 0), P(1, 0, 0), P(0, 0, 1));
    const geom::PlaneD pl2 = geom::plane_from_triangle(P(0, 0, 0), P(0, 0, 1), P(0, 1, 0));
    const geom::PlaneD pl3 = geom::plane_from_triangle(P(0, 0, 0), P(1, 0, 0), P(0, 1, 0));
    const geom::PlaneD pl4 = geom::plane_from_triangle(P(0, 0, 10), P(1, 0, 10), P(0, 1, 10));
    f.p1 = f.soup.table.intern(pl1).id;
    f.p2 = f.soup.table.intern(pl2).id;
    f.p3 = f.soup.table.intern(pl3).id;
    f.p4 = f.soup.table.intern(pl4).id;
    f.mid = geom::edge_interior_point(pl1, pl2, pl3, pl4, geom::to_homogeneous(P(0, 0, 0)),
                                      geom::to_homogeneous(P(0, 0, 10)));

    // 頂点は **lex_less の昇順**に置きます（前半の二分探索の前提）
    std::vector<geom::IPoint> pts = {P(0, 0, 0), P(0, 0, 10), P(1, 0, 0),
                                     P(0, 1, 0), P(-1, 0, 0), P(0, -1, 0)};
    if (extra) {
        pts.push_back(P(0, 0, 0));   // a2（同じ位置の複製）
        pts.push_back(P(0, 0, 10));  // b2
    }
    std::vector<geom::HPointD> hv;
    hv.reserve(pts.size());
    for (const geom::IPoint& p : pts) hv.push_back(geom::to_homogeneous(p));
    std::vector<std::uint32_t> ord(hv.size());
    for (std::uint32_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(), [&](std::uint32_t x, std::uint32_t y) {
        return geom::lex_less(hv[x], hv[y]);
    });
    std::vector<std::uint32_t> pos(hv.size());
    for (std::uint32_t i = 0; i < ord.size(); ++i) pos[ord[i]] = i;
    for (std::uint32_t i : ord) f.mesh.vertices.push_back(hv[i]);
    f.a = pos[0];
    f.b = pos[1];
    if (extra) {
        f.a2 = pos[6];
        f.b2 = pos[7];
    }
    const std::uint32_t q1 = pos[2], q2 = pos[3], q3 = pos[4], q4 = pos[5];

    f.mesh.vertex_key.assign(f.mesh.vertices.size(), {csg::kNoPlane, csg::kNoPlane, csg::kNoPlane});
    f.mesh.vertex_key[f.a] = {f.p1, f.p2, f.p3};
    f.mesh.vertex_key[f.b] = {f.p1, f.p2, f.p4};
    if (extra) {
        f.mesh.vertex_key[f.a2] = {f.p1, f.p2, f.p3};
        f.mesh.vertex_key[f.b2] = {f.p1, f.p2, f.p4};
    }
    f.mesh.vertex_merged.assign(f.mesh.vertices.size(), 1);

    const auto add_poly = [&](csg::PlaneId sup) {
        csg::Poly q;
        q.frag.support = sup;
        f.soup.polys.push_back(q);
        return static_cast<std::uint32_t>(f.soup.polys.size() - 1);
    };
    const auto add_tri = [&](std::uint32_t x, std::uint32_t y, std::uint32_t z, csg::PlaneId sup) {
        f.mesh.triangles.push_back({x, y, z});
        f.mesh.tri_poly.push_back(add_poly(sup));
        f.mesh.tri_src.push_back(0);
        f.mesh.tri_tag.push_back(0);
    };
    add_tri(f.a, f.b, q1, f.p1);
    add_tri(f.a, f.b, q2, f.p2);
    add_tri(f.a, f.b, q3, f.p1);
    add_tri(f.a, f.b, q4, f.p2);
    if (extra) {
        add_tri(f.a2, f.b2, q1, f.p1);
        add_tri(f.a2, f.b2, q2, f.p2);
        add_tri(f.a2, f.b2, q3, f.p1);
        add_tri(f.a2, f.b2, q4, f.p2);
    }

    f.st.merged_points = f.mesh.vertices.size();
    mesh::SplitStats::UnresolvedEdge ue{};
    ue.out_a = f.a;
    ue.out_b = f.b;
    ue.degree = 4;
    ue.pair_groups = 2;
    ue.pair_tris[0] = 0;
    ue.pair_tris[1] = 2;
    ue.pair_tris[2] = 1;
    ue.pair_tris[3] = 3;
    f.st.split.unresolved_detail.push_back(ue);
    if (extra) {
        mesh::SplitStats::UnresolvedEdge u2{};
        u2.out_a = f.a2;
        u2.out_b = f.b2;
        u2.degree = 4;
        u2.pair_groups = 2;
        u2.pair_tris[0] = 4;
        u2.pair_tris[1] = 6;
        u2.pair_tris[2] = 5;
        u2.pair_tris[3] = 7;
        f.st.split.unresolved_detail.push_back(u2);
    }
    return f;
}

/// **その候補の処理で、メッシュの内容が変わっていないこと。**
///
/// **長さだけでなく、補助配列の【内容】も比べます。**
struct Snap {
    std::vector<geom::HPointD> verts;
    std::vector<mesh::Tri> tris;
    std::vector<std::array<csg::PlaneId, 3>> keys;
    std::vector<std::uint32_t> merged, tpoly, vsrc, ttag;
    std::vector<int> tsrc;
    std::vector<csg::EdgeSplitSource> esplit;
};

Snap snap(const csg::SoupMesh& m) {
    Snap s;
    s.verts = m.vertices;
    s.tris = m.triangles;
    s.keys = m.vertex_key;
    s.merged = m.vertex_merged;
    s.tpoly = m.tri_poly;
    s.vsrc = m.vertex_split_src;
    s.ttag = m.tri_tag;
    s.tsrc = m.tri_src;
    s.esplit = m.edge_split;
    return s;
}

bool same_edge_split(const csg::EdgeSplitSource& a, const csg::EdgeSplitSource& b) {
    return a.v == b.v && a.w == b.w && a.p1 == b.p1 && a.p2 == b.p2 && a.p3 == b.p3 &&
           a.p4 == b.p4 && a.sign == b.sign;
}

/// **すべての配列を、長さと内容で比べます。**
///
/// **頂点は `h_equal`（幾何としての等値）**で比べます。**同次座標の表現までは見ません**
/// — スカラー倍が違っても同じ点なら一致と扱います。他の配列は値そのものを比べます。
bool same_content(const Snap& s, const csg::SoupMesh& m) {
    if (s.verts.size() != m.vertices.size() || s.tris.size() != m.triangles.size()) return false;
    if (s.keys.size() != m.vertex_key.size() || s.esplit.size() != m.edge_split.size()) {
        return false;
    }
    for (std::size_t i = 0; i < s.verts.size(); ++i) {
        if (!geom::h_equal(s.verts[i], m.vertices[i])) return false;
    }
    for (std::size_t i = 0; i < s.tris.size(); ++i) {
        if (s.tris[i] != m.triangles[i]) return false;
    }
    for (std::size_t i = 0; i < s.keys.size(); ++i) {
        if (s.keys[i] != m.vertex_key[i]) return false;
    }
    if (s.merged != m.vertex_merged || s.tpoly != m.tri_poly) return false;
    if (s.vsrc != m.vertex_split_src || s.ttag != m.tri_tag || s.tsrc != m.tri_src) return false;
    for (std::size_t i = 0; i < s.esplit.size(); ++i) {
        if (!same_edge_split(s.esplit[i], m.edge_split[i])) return false;
    }
    return true;
}

/// **成功した候補の影響と、拒否した候補の影響を分ける**ための対照。
bool same_mesh(const csg::SoupMesh& a, const csg::SoupMesh& b) {
    return same_content(snap(a), b);
}

}  // namespace

int main() {
    std::printf("## 修復の拒否分岐（人工の入力。公開経路からの到達は未確認）\n\n");
    std::printf("| # | 構成 | 狙う分岐 | 細分 | 衝突 | 平面不足 | 組不足 | 照合 | 内容不変 |\n");
    std::printf("|---|---|---|---:|---:|---:|---:|---:|---|\n");

    // ---- 1. 衝突なし ----
    {
        Fixture f = make(false);
        csg::detail::repair_unresolved_edges(f.soup, f.mesh, f.st);
        const auto& sp = f.st.split;
        std::printf("| 1 | 衝突なし | — | %zu | %zu | %zu | %zu | %zu | — |\n", sp.repair_edges,
                    sp.repair_collisions, sp.repair_no_planes, sp.repair_no_pair,
                    sp.repair_collide_probes);
        expect(sp.repair_edges == 1, "1: 細分した");
        expect(sp.repair_collisions == 0 && sp.repair_no_planes == 0 && sp.repair_no_pair == 0,
               "1: 拒否していない");
        // ---- 2. 同一辺の意図した 2 頂点は、衝突として拒否しない ----
        expect(f.mesh.vertices.size() == 8, "2: 頂点が 2 個増えた（同じ位置の 2 複製）");
        expect(geom::h_equal(f.mesh.vertices[6], f.mesh.vertices[7]), "2: 2 複製は同じ位置");
        expect(geom::h_equal(f.mesh.vertices[6], f.mid), "2: 位置が期待どおり");
        std::printf("| 2 | 同一辺の 2 複製 | 拒否しない | %zu | %zu | — | — | — | — |\n",
                    sp.repair_edges, sp.repair_collisions);
    }

    // ---- 3. 既存頂点との衝突 ----
    {
        Fixture f = make(false);
        // **前半（整列済み）に細分点と同じ位置の頂点を混ぜます。**
        // 前半は lex_less の昇順でなければならないので、入れ直して並べ替えます
        f.mesh.vertices.push_back(f.mid);
        std::vector<geom::HPointD>& v = f.mesh.vertices;
        std::stable_sort(v.begin(), v.end(), [](const geom::HPointD& x, const geom::HPointD& y) {
            return geom::lex_less(x, y);
        });
        f.st.merged_points = v.size();
        f.mesh.vertex_key.resize(v.size(), {csg::kNoPlane, csg::kNoPlane, csg::kNoPlane});
        f.mesh.vertex_merged.resize(v.size(), 1);
        // **並べ替えたので、辺の添字を取り直します**（前提を壊さないため）
        for (std::uint32_t i = 0; i < v.size(); ++i) {
            if (geom::h_equal(v[i], geom::to_homogeneous(P(0, 0, 0)))) f.a = i;
            if (geom::h_equal(v[i], geom::to_homogeneous(P(0, 0, 10)))) f.b = i;
        }
        for (mesh::Tri& t : f.mesh.triangles) {
            for (int k = 0; k < 3; ++k) {
                // 添字は並べ替えで動くので、位置で引き直します
                (void)k;
            }
        }
        // 三角形は作り直します（並べ替えで添字が変わったため）
        f.mesh.triangles.clear();
        f.mesh.tri_poly.clear();
        f.mesh.tri_src.clear();
        f.mesh.tri_tag.clear();
        const auto idx = [&](geom::IPoint p) {
            const geom::HPointD h = geom::to_homogeneous(p);
            for (std::uint32_t i = 0; i < v.size(); ++i) {
                if (geom::h_equal(v[i], h)) return i;
            }
            return std::uint32_t{0};
        };
        const std::uint32_t q[4] = {idx(P(1, 0, 0)), idx(P(0, 1, 0)), idx(P(-1, 0, 0)),
                                    idx(P(0, -1, 0))};
        const csg::PlaneId sup[4] = {f.p1, f.p2, f.p1, f.p2};
        for (int i = 0; i < 4; ++i) {
            f.mesh.triangles.push_back({f.a, f.b, q[i]});
            csg::Poly pq;
            pq.frag.support = sup[i];
            f.soup.polys.push_back(pq);
            f.mesh.tri_poly.push_back(static_cast<std::uint32_t>(f.soup.polys.size() - 1));
            f.mesh.tri_src.push_back(0);
            f.mesh.tri_tag.push_back(0);
        }
        f.mesh.vertex_key[f.a] = {f.p1, f.p2, f.p3};
        f.mesh.vertex_key[f.b] = {f.p1, f.p2, f.p4};
        f.st.split.unresolved_detail[0].out_a = f.a;
        f.st.split.unresolved_detail[0].out_b = f.b;
        f.st.split.unresolved_detail[0].pair_tris[0] = 0;
        f.st.split.unresolved_detail[0].pair_tris[1] = 2;
        f.st.split.unresolved_detail[0].pair_tris[2] = 1;
        f.st.split.unresolved_detail[0].pair_tris[3] = 3;

        // **候補の処理が始まる直前の内容**（入口の assign は含めない）
        f.mesh.vertex_split_src.assign(f.mesh.vertices.size(), csg::kNoEdgeSplit);
        const Snap before = snap(f.mesh);
        csg::detail::repair_unresolved_edges(f.soup, f.mesh, f.st);
        const auto& sp = f.st.split;
        const bool inv = same_content(before, f.mesh);
        std::printf(
            "| 3 | 既存頂点と一致 | `repair_collisions` | %zu | %zu | %zu | %zu | %zu |"
            " %s |\n",
            sp.repair_edges, sp.repair_collisions, sp.repair_no_planes, sp.repair_no_pair,
            sp.repair_collide_probes, inv ? "はい" : "**いいえ**");
        expect(sp.repair_collisions == 1, "3: **衝突として数えた**");
        expect(sp.repair_edges == 0, "3: 細分していない");
        expect(inv, "3: その候補でメッシュの内容が変わっていない");
    }

    // ---- 4. 別の辺の細分点との衝突 ----
    //
    // **対照**: 同じ初期状態から **1 本目だけ**処理した結果と比べます。
    // **これで「成功した候補の影響」と「拒否した候補の影響」を分けられます。**
    {
        Fixture ctl = make(true);
        ctl.st.split.unresolved_detail.resize(1);
        csg::detail::repair_unresolved_edges(ctl.soup, ctl.mesh, ctl.st);

        Fixture f = make(true);
        csg::detail::repair_unresolved_edges(f.soup, f.mesh, f.st);
        const auto& sp = f.st.split;
        const bool same = same_mesh(ctl.mesh, f.mesh);
        std::printf(
            "| 4 | 別の辺の細分点と一致 | `repair_collisions` | %zu | %zu | %zu | %zu |"
            " %zu | %s |\n",
            sp.repair_edges, sp.repair_collisions, sp.repair_no_planes, sp.repair_no_pair,
            sp.repair_collide_probes, same ? "はい" : "**いいえ**");
        expect(sp.repair_edges == 1, "4: 1 本目だけ細分した");
        expect(sp.repair_collisions == 1, "4: **2 本目を衝突として数えた**");
        expect(f.mesh.edge_split.size() == 1, "4: 由来は 1 件だけ");
        expect(same,
               "4: **2 本目は、1 本目だけ処理した対照と一致**（頂点は h_equal、他は値の一致）");
    }

    std::printf("\n**不一致 %d 件**\n", failures);
    return failures == 0 ? 0 : 1;
}

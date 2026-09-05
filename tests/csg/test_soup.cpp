// Krisite — 入口・中核・出口の分離（SPEC-phase3.md §14 の CP2）
//
//     from_mesh : TriMesh → PolySoup
//     boolean   : PolySoup × PolySoup → PolySoup   ★ CSG について閉じる
//     to_mesh   : PolySoup → TriMesh
//
// **見るのは 4 つです。**
//
//   1. 往復（`to_mesh(from_mesh(m))`）が元の位相を再現すること
//   2. §10.1 二項正解器（`boolean_op`）と $(C, \chi)$ が一致すること
//   3. §10.3 連鎖が丸めを経由しないこと（**中間に TriMesh を作らない**）
//   4. §10.3.1 連鎖でビット幅が伸びないこと ★ 契約の成立条件
//
// **三角形の集合と断片数は一致しません**（§10.1）。面併合をやめたので分割が違います。
#include <algorithm>
#include <array>
#include <cstdio>
#include <string>
#include <type_traits>
#include <vector>

#include "krisite/csg/boolean.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/topology.hpp"

#include "corpus.hpp"
#include "test_util.hpp"

using namespace krisite;
using krisite::mesh::check_topology;
using krisite::mesh::TopologyReport;
using krisite::mesh::TriMesh;

namespace {

constexpr unsigned kMaxDepth = 3;

/// スープの構成点の最大ビット幅（§10.3.1）。
std::size_t soup_bits(const csg::PolySoup& s) {
    std::size_t mx = 0;
    for (const csg::Poly& q : s.polys) {
        for (std::size_t i = 0; i < csg::vertex_count(q.frag); ++i) {
            const geom::HPointD v = csg::fragment_vertex(s.table, q.frag, i);
            for (std::size_t w : {arith::min_bits(v.x), arith::min_bits(v.y), arith::min_bits(v.z),
                                  arith::min_bits(v.w)}) {
                mx = (w > mx) ? w : mx;
            }
        }
    }
    return mx;
}

const char* op_name(csg::BoolOp op) {
    switch (op) {
        case csg::BoolOp::Union:
            return "∪";
        case csg::BoolOp::Intersection:
            return "∩";
        default:
            return "\\";
    }
}

/// 頂点の**値**で見た三角形の多重集合（順序非依存の比較。`SPEC-phase2.md` §9.4.2）。
std::vector<std::array<std::uint32_t, 3>> geometric_key(const csg::SoupMesh& m) {
    std::vector<std::uint32_t> ord(m.vertices.size());
    for (std::uint32_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](std::uint32_t a, std::uint32_t b) {
        return geom::lex_less(m.vertices[a], m.vertices[b]);
    });
    std::vector<std::uint32_t> canon(m.vertices.size(), 0);
    std::uint32_t next = 0;
    for (std::size_t i = 0; i < ord.size();) {
        std::size_t j = i;
        while (j < ord.size() && geom::h_equal(m.vertices[ord[i]], m.vertices[ord[j]])) {
            canon[ord[j]] = next;
            ++j;
        }
        ++next;
        i = j;
    }
    std::vector<std::array<std::uint32_t, 3>> out;
    out.reserve(m.triangles.size());
    for (const mesh::Tri& t : m.triangles) {
        const std::uint32_t c[3] = {canon[t[0]], canon[t[1]], canon[t[2]]};
        int s = 0;
        if (c[1] < c[s]) s = 1;
        if (c[2] < c[s]) s = 2;
        out.push_back({c[s], c[(s + 1) % 3], c[(s + 2) % 3]});
    }
    std::sort(out.begin(), out.end());
    return out;
}

/// **並列化の前提の検査**（`SPEC-phase3.md` §14 の CP3 の判定）。
///
/// > 判定基準は「葉の分類が、可変な共有状態に触れずに書けるか」です。
///
/// 分類が可変な共有状態に依存していれば、**領域を回す順序を変えると結果が変わります。**
/// 依存していなければ、順序を変えても幾何は同じです。
///
/// **これは「書ける」ことの証拠であって、証明ではありません。** 依存があっても
/// たまたま同じ結果になる可能性は残ります。ただし、依存が入ったときに落ちる網には
/// なります（`CLAUDE.md`「機構を足したら、その機構が空回りしていないことを別に検査」）。
void test_classification_is_order_independent() {
    std::size_t n = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (csg::BoolOp op :
             {csg::BoolOp::Union, csg::BoolOp::Intersection, csg::BoolOp::Difference}) {
            for (unsigned d = 0; d <= kMaxDepth; ++d) {
                csg::BoolOptions o1 = kritest::phase1_options(d);
                csg::BoolOptions o2 = o1;
                o2.reverse_regions = true;
                csg::ToMeshOptions tm;
                tm.split_contacts = false;
                const csg::SoupMesh m1 =
                    csg::to_mesh(csg::boolean(csg::from_mesh(a), csg::from_mesh(b), op, o1), tm);
                const csg::SoupMesh m2 =
                    csg::to_mesh(csg::boolean(csg::from_mesh(a), csg::from_mesh(b), op, o2), tm);
                const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) +
                                        "（深度 " + std::to_string(d) + "）";
                KRI_CHECK_MSG(geometric_key(m1) == geometric_key(m2),
                              tag +
                                  ": **領域を回す順序で結果が変わった。** 分類が可変な"
                                  "共有状態に依存しています（並列化の前提が崩れます）");
                ++n;
            }
        }
    }
    std::printf("    順序非依存 %zu 件（並列化の前提。§14 の CP3 の判定）\n", n);
}

/// **CP4: 局所 BSP ↔ 過剰分割**（`SPEC-phase3.md` §5.4 / §10.1 / §14 の CP4）。
///
/// **一致するのは $(C, \chi)$ と体積で、三角形の集合と断片数は一致しません。**
/// 断片数が減ることが目的なので、**減っていること自体も検査します。**
///
/// 体積は GMP を要するので `csg/test_volume_gmp.cpp` の担当です。ここでは位相と
/// 断片数、そして**判定が空回りしていないこと**を見ます。
void test_local_bsp_vs_over_subdivision() {
    std::size_t n = 0, differ = 0;
    std::size_t raw_over = 0, raw_bsp = 0, out_over = 0, out_bsp = 0;
    std::size_t cuts_used = 0, cuts_skipped = 0, slots = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (csg::BoolOp op :
             {csg::BoolOp::Union, csg::BoolOp::Intersection, csg::BoolOp::Difference}) {
            for (unsigned d = 0; d <= kMaxDepth; ++d) {
                // **両側を明示します。** 既定に依存させると、既定が変わったときに
                // 自分自身との比較になります（`CLAUDE.md`）。
                csg::BoolOptions oo = kritest::phase1_options(d);
                oo.local_bsp = false;
                csg::BoolOptions ob = oo;
                ob.local_bsp = true;
                csg::ToMeshOptions tm;
                tm.split_contacts = false;
                csg::BoolStats so{}, sb{};
                const csg::SoupMesh mo = csg::to_mesh(
                    csg::boolean(csg::from_mesh(a), csg::from_mesh(b), op, oo, &so), tm);
                const csg::SoupMesh mb = csg::to_mesh(
                    csg::boolean(csg::from_mesh(a), csg::from_mesh(b), op, ob, &sb), tm);
                const TopologyReport t0 = check_topology(mo.triangles);
                const TopologyReport t1 = check_topology(mb.triangles);
                const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) +
                                        "（深度 " + std::to_string(d) + "）";
                KRI_CHECK_MSG(t0.components == t1.components,
                              tag + ": 局所 BSP で C が変わった" +
                                  kritest::pair_msg(t0.components, t1.components));
                KRI_CHECK_MSG(t0.chi == t1.chi, tag + ": 局所 BSP で χ が変わった" +
                                                    kritest::pair_msg(t0.chi, t1.chi));
                KRI_CHECK_MSG(t1.empty || t1.oriented,
                              tag + ": 局所 BSP の出力の向きが整合していない");
                KRI_CHECK_MSG(sb.raw_fragments <= so.raw_fragments,
                              tag + ": **局所 BSP のほうが生の断片が多い**" +
                                  kritest::pair_msg(so.raw_fragments, sb.raw_fragments));
                if (geometric_key(mo) != geometric_key(mb)) ++differ;
                raw_over += so.raw_fragments;
                raw_bsp += sb.raw_fragments;
                out_over += so.fragments;
                out_bsp += sb.fragments;
                slots += sb.bsp_cut_slots;
                cuts_used += sb.bsp_cuts_used;
                cuts_skipped += sb.bsp_cuts_skipped;
                KRI_CHECK_MSG(so.bsp_cut_slots == 0,
                              tag + ": 過剰分割の側で局所 BSP の判定が走っている");
                ++n;
            }
        }
    }
    // **空回りの検査。** どちらかが 0 なら判定が一方に倒れていて、機構が働いていません。
    KRI_CHECK_MSG(cuts_used > 0, "局所 BSP が 1 枚も切っていない（判定が空回り）");
    KRI_CHECK_MSG(cuts_skipped > 0,
                  "局所 BSP が 1 枚も省いていない（過剰分割と同じ。判定が空回り）");
    KRI_CHECK_MSG(slots == cuts_used + cuts_skipped, "候補の総数が使用 + 省略と合わない");
    // **三角形の集合は一致しないはず**（それが目的）。すべて一致したら置き換えが no-op。
    KRI_CHECK_MSG(differ > 0, "**局所 BSP と過剰分割で出力が 1 件も変わらない**（no-op）");
    KRI_CHECK_MSG(raw_bsp < raw_over, "**生の断片が減っていない**（§14 の CP4 の目的）");
    std::printf("    局所 BSP %zu 件: 生の断片 %zu → %zu（%.1f%%）, 正準化後 %zu → %zu（%.1f%%）\n",
                n, raw_over, raw_bsp,
                100.0 * static_cast<double>(raw_bsp) / static_cast<double>(raw_over), out_over,
                out_bsp, 100.0 * static_cast<double>(out_bsp) / static_cast<double>(out_over));
    std::printf("    切断候補 %zu（切った %zu / 省いた %zu）, 出力が変わった %zu / %zu 件\n", slots,
                cuts_used, cuts_skipped, differ, n);
}

/// 1. 往復。**入口と出口だけで閉じる検査**なので、中核の誤りが混ざりません。
void test_round_trip() {
    std::size_t n = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        for (int which = 0; which < 2; ++which) {
            const TriMesh m = (which == 0) ? c.make_a() : c.make_b();
            const csg::SoupMesh r = csg::to_mesh(csg::from_mesh(m));
            const TopologyReport t0 = check_topology(m.triangles);
            const TopologyReport t1 = check_topology(r.triangles);
            const std::string tag = std::string("ケース ") + c.id + (which == 0 ? " A" : " B");
            KRI_CHECK_MSG(t1.ok(), tag + ": 往復で多様体でなくなった");
            KRI_CHECK_MSG(
                t0.components == t1.components,
                tag + ": 往復で C が変わった" + kritest::pair_msg(t0.components, t1.components));
            KRI_CHECK_MSG(t0.chi == t1.chi,
                          tag + ": 往復で χ が変わった" + kritest::pair_msg(t0.chi, t1.chi));
            ++n;
        }
    }
    std::printf("    往復 %zu 件\n", n);
}

/// 2. §10.1 二項正解器との一致。**分裂は両方 OFF で比べます**（意味論を揃える）。
void test_matches_binary_oracle() {
    std::size_t n = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (csg::BoolOp op :
             {csg::BoolOp::Union, csg::BoolOp::Intersection, csg::BoolOp::Difference}) {
            for (unsigned d = 0; d <= kMaxDepth; ++d) {
                const csg::BoolOptions o = kritest::phase1_options(d);
                const csg::BoolMesh ref = csg::boolean_op(a, b, op, o);
                csg::ToMeshOptions tm;
                tm.split_contacts = false;
                const csg::SoupMesh got =
                    csg::to_mesh(csg::boolean(csg::from_mesh(a), csg::from_mesh(b), op, o), tm);
                const TopologyReport t0 = check_topology(ref.triangles);
                const TopologyReport t1 = check_topology(got.triangles);
                const std::string tag = std::string("ケース ") + c.id + " " + op_name(op) +
                                        "（深度 " + std::to_string(d) + "）";
                KRI_CHECK_MSG(
                    t0.components == t1.components,
                    tag + ": C が正解器と違う" + kritest::pair_msg(t0.components, t1.components));
                KRI_CHECK_MSG(t0.chi == t1.chi,
                              tag + ": χ が正解器と違う" + kritest::pair_msg(t0.chi, t1.chi));
                // **(C, χ) は面の向きを反転しても変わりません。**
                // 向きの誤りはここでしか出ないので、別に固定します（変異 15）。
                // 接触の分裂を切っているので `ok()` は使えません（次数 4 の辺が残る）。
                // `oriented` は各辺で #(u,v) = #(v,u) を見るので次数 4 でも有効です。
                // **空の出力（交わらない立体の ∩ など）では向きは定義されません。**
                KRI_CHECK_MSG(t1.empty || t1.oriented, tag + ": **出力の向きが整合していない**");
                ++n;
            }
        }
    }
    std::printf("    正解器との比較 %zu 件\n", n);
}

/// 3. §10.3 連鎖の厳密性 + 4. §10.3.1 ビット幅が伸びないこと。
void test_chain_is_exact() {
    // **型が閉じていること**をコンパイル時に固定する
    static_assert(std::is_same_v<decltype(csg::boolean(std::declval<const csg::PolySoup&>(),
                                                       std::declval<const csg::PolySoup&>(),
                                                       csg::BoolOp::Union, csg::BoolOptions{})),
                                 csg::PolySoup>,
                  "boolean の戻り値が PolySoup でない（型が閉じていない）");

    std::size_t b1 = 0, b2 = 0, b3 = 0, n = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        // 3 つ目は**別のメッシュ**。同じメッシュを 2 度使うと曲面が二重になり、
        // 内外の 1 ビットでは「2 度跨いだ」を表せません（CP3 の WNV が要る）。
        const TriMesh d = kritest::corpus()[1].make_b();
        const csg::BoolOptions o = kritest::phase1_options(1);
        const csg::PolySoup SA = csg::from_mesh(a), SB = csg::from_mesh(b), SD = csg::from_mesh(d);

        const csg::PolySoup s1 = csg::boolean(SA, SB, csg::BoolOp::Union, o);
        const csg::PolySoup s2 = csg::boolean(s1, SD, csg::BoolOp::Difference, o);
        const csg::PolySoup s3 = csg::boolean(s2, SA, csg::BoolOp::Union, o);

        const std::string tag = std::string("ケース ") + c.id;
        // **中間に TriMesh を作っていないこと。** sources は入力そのもので、
        // 段が増えても【増えるだけ】で、丸められた頂点が現れません。
        KRI_CHECK_MSG(s1.source_count() == 2 && s2.source_count() == 3 && s3.source_count() == 4,
                      tag + ": sources の数が段数と合わない");
        for (std::size_t i = 0; i < a.vertices.size(); ++i) {
            KRI_CHECK_MSG(s3.sources[0].vertices[i].x == a.vertices[i].x &&
                              s3.sources[0].vertices[i].y == a.vertices[i].y &&
                              s3.sources[0].vertices[i].z == a.vertices[i].z,
                          tag + ": 連鎖で source の頂点が変わった（丸めが入っている）");
        }
        b1 = std::max(b1, soup_bits(s1));
        b2 = std::max(b2, soup_bits(s2));
        b3 = std::max(b3, soup_bits(s3));

        // **同じメッシュを 2 度使う連鎖**（CP3 の段 1 で開きました）。
        //
        // 内外の 1 ビットでは「同じ曲面を 2 度跨いだ」を表せず、CP2 では 17 / 44 が
        // 食い違っていました。**巻き数なら w = 2 になるだけ**です。
        {
            const csg::PolySoup ub = csg::boolean(csg::boolean(SA, SB, csg::BoolOp::Union, o), SB,
                                                  csg::BoolOp::Difference, o);
            const csg::PolySoup ab = csg::boolean(SA, SB, csg::BoolOp::Difference, o);
            csg::ToMeshOptions t2;
            t2.split_contacts = false;
            const TopologyReport t_ub = check_topology(csg::to_mesh(ub, t2).triangles);
            const TopologyReport t_ab = check_topology(csg::to_mesh(ab, t2).triangles);
            KRI_CHECK_MSG(t_ub.components == t_ab.components && t_ub.chi == t_ab.chi,
                          tag + ": (A∪B)\\B と A\\B が食い違う（同じ曲面を 2 度跨ぐ配置）" +
                              kritest::pair_msg(t_ab.chi, t_ub.chi));
        }

        // 自己整合: (A ∪ B) \ D  ≡  (A \ D) ∪ (B \ D)
        csg::ToMeshOptions tm;
        tm.split_contacts = false;
        const csg::PolySoup rhs =
            csg::boolean(csg::boolean(SA, SD, csg::BoolOp::Difference, o),
                         csg::boolean(SB, SD, csg::BoolOp::Difference, o), csg::BoolOp::Union, o);
        const TopologyReport tl = check_topology(csg::to_mesh(s2, tm).triangles);
        const TopologyReport tr = check_topology(csg::to_mesh(rhs, tm).triangles);
        KRI_CHECK_MSG(
            tl.components == tr.components && tl.chi == tr.chi,
            tag + ": (A∪B)\\D と (A\\D)∪(B\\D) が食い違う" + kritest::pair_msg(tl.chi, tr.chi));
        ++n;
    }

    std::printf("    連鎖 %zu 件。最大ビット幅: 1 段 %zu / 2 段 %zu / 3 段 %zu（上界 %zu）\n", n,
                b1, b2, b3, geom::bits::kHomoXyz);
    // §10.3.1: **連鎖でビット幅は伸びません。** CSG は新しい平面を作らないので、
    // 構成点は入力平面の集合から 3 枚を選んだ交点のままです。
    KRI_CHECK_MSG(b2 == b1 && b3 == b1,
                  "**連鎖でビット幅が伸びた。** どこかで平面を構成しています（§10.3.1）" +
                      kritest::pair_msg(b1, b3));
    KRI_CHECK_MSG(b1 <= geom::bits::kHomoXyz, "構成点が理論上界を超えた");
    KRI_CHECK_MSG(b1 > 0, "構成点が 1 つも無い。**空回りです**");
}

/// **A-3（セルで区切った T 字接合の索引）**（`DESIGN-phase5-hotspots.md` §6.3）。
///
/// **見るのは 3 つです。**
///
///   1. **旗の ON / OFF で出力がバイト単位で一致すること**（厳密な絞り込みなので）
///   2. **挿入された T 頂点が 1 個も変わらないこと**（絞り込みが厳密であることの直接の検査）
///   3. ★ **機構が実際に発火したこと**（`cell_index_groups > 0`）
///
/// **3 が要ります。** 箱がセルの箱でなければ従来の索引に退避する形なので、
/// **退避したまま緑になり得ます**（`CLAUDE.md`「足した機構が実際に発火したことを、
/// テスト自身に確かめさせること」）。
///
/// **`from_mesh` 直後のスープでは `Poly::aabb` は三角形の外接箱**なので、
/// **必ず `boolean` を通した出力で確かめます。**
void test_cell_index() {
    std::printf("  A-3: セルで区切った T 字接合の索引\n");
    std::size_t n = 0, fired = 0, groups_total = 0;
    for (const kritest::Case& c : kritest::corpus()) {
        const TriMesh a = c.make_a(), b = c.make_b();
        for (unsigned d = 0; d <= kMaxDepth; ++d) {
            const csg::BoolOptions o = kritest::phase1_options(d);
            for (csg::BoolOp op : {csg::BoolOp::Union, csg::BoolOp::Intersection,
                                   csg::BoolOp::Difference}) {
                const csg::PolySoup s =
                    csg::boolean(csg::from_mesh(a), csg::from_mesh(b), op, o);
                if (s.polys.empty()) continue;
                csg::ToMeshOptions on, off;
                on.cell_index = true;
                off.cell_index = false;
                csg::ToMeshStats st_on, st_off;
                const csg::SoupMesh m_on = csg::to_mesh(s, on, &st_on);
                const csg::SoupMesh m_off = csg::to_mesh(s, off, &st_off);
                const std::string tag =
                    std::string(c.id) + " 深度 " + std::to_string(d) + " 演算 " +
                    std::to_string(static_cast<int>(op));
                ++n;
                if (st_on.cell_index_groups > 0) ++fired;
                groups_total += st_on.cell_index_groups;
                // 1. バイト一致（頂点も三角形も同じ並びで同じ値）
                KRI_CHECK_MSG(m_on.vertices.size() == m_off.vertices.size(),
                              tag + ": A-3 で頂点数が変わった" +
                                  kritest::pair_msg(m_on.vertices.size(), m_off.vertices.size()));
                KRI_CHECK_MSG(m_on.triangles == m_off.triangles,
                              tag + ": A-3 で三角形の並びが変わった");
                bool same = m_on.vertices.size() == m_off.vertices.size();
                for (std::size_t i = 0; same && i < m_on.vertices.size(); ++i) {
                    same = geom::cmp_h_lex(m_on.vertices[i], m_off.vertices[i]) == 0;
                }
                KRI_CHECK_MSG(same, tag + ": A-3 で頂点の値が変わった");
                // 2. T 頂点が 1 個も変わらない（**絞り込みが厳密であることの直接の検査**）
                KRI_CHECK_MSG(st_on.t.inserted == st_off.t.inserted,
                              tag + ": A-3 で T 頂点の数が変わった。**絞り込みが漏れています**" +
                                  kritest::pair_msg(st_on.t.inserted, st_off.t.inserted));
                // 絞り込みが効いていること（候補が増えることはない）
                KRI_CHECK_MSG(st_on.t.candidates <= st_off.t.candidates,
                              tag + ": A-3 で候補が増えた" +
                                  kritest::pair_msg(st_on.t.candidates, st_off.t.candidates));
            }
        }
    }
    std::printf("    %zu 件。**機構が発火したのは %zu 件**（(葉,支持平面) の組 計 %zu）\n", n,
                fired, groups_total);
    // 3. ★ **空回りの検査。** 退避したまま緑になっていないか。
    KRI_CHECK_MSG(n > 0, "A-3: 比較が 1 件も回っていない。**空回りです**");
    KRI_CHECK_MSG(fired == n,
                  "**A-3 が発火していない件がある。** `boolean` の出力なら箱はセルの箱の"
                  "はずで、退避は起きません" +
                      kritest::pair_msg(fired, n));
}

}  // namespace

int main() {
    std::printf("\n  入口・中核・出口の分離 — SPEC-phase3 §14 の CP2 / CP3 / CP4\n");
    test_round_trip();
    test_matches_binary_oracle();
    test_chain_is_exact();
    test_classification_is_order_independent();
    test_local_bsp_vs_over_subdivision();
    test_cell_index();
    std::printf("\n");
    return kritest::finish("csg/soup");
}

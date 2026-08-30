// Krisite — Thingi10K の読み込みと量子化（`SPEC-phase5.md` §2）
//
// **これはテスト側のツールです。** ライブラリ本体は整数座標のメッシュしか受けません。
// 実データは float なので、**量子化はここで行います。**
//
// ---
//
// ## 量子化は入力の性質を変えます ★
//
// `∂S` は組合せ的な量なので**量子化で変わりません**（三角形の並びが同じなら同じ）。
// **変わるのは幾何です。**
//
//   面積 0 の三角形   量子化で 3 頂点が同一格子点に落ちると生まれる
//   自己交差         もとの float では交わっていなくても、動かせば交わり得る
//
// **したがって「Thingi10K のメタデータで自己交差なし」は
// 「量子化後に自己交差なし」を意味しません。** 分類は**量子化後の性質**で行います
// （`SPEC-phase5.md` §1。Krisite が実際に受け取るのは量子化後のメッシュです）。
#ifndef KRISITE_TESTS_THINGI10K_LOADER_HPP
#define KRISITE_TESTS_THINGI10K_LOADER_HPP

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <array>
#include <map>
#include <string>
#include <vector>

#include "krisite/config.hpp"
#include "krisite/mesh/tri_mesh.hpp"

namespace krithingi {

/// float64 の生メッシュ（`fetch.py` が書く `.kmesh`）。
struct RawMesh {
    std::vector<double> v;                 ///< 3 * nv
    std::vector<std::uint32_t> f;          ///< 3 * nf
    std::size_t nv = 0, nf = 0;
    bool ok = false;
};

/// `.kmesh` を読む。**リトルエンディアン固定**（`fetch.py` と対）。
inline RawMesh load_kmesh(const std::string& path) {
    RawMesh m;
    std::FILE* fp = std::fopen(path.c_str(), "rb");
    if (fp == nullptr) return m;
    char magic[4];
    std::uint32_t ver = 0, nv = 0, nf = 0;
    if (std::fread(magic, 1, 4, fp) != 4 || std::string(magic, 4) != "KMSH" ||
        std::fread(&ver, 4, 1, fp) != 1 || ver != 1 || std::fread(&nv, 4, 1, fp) != 1 ||
        std::fread(&nf, 4, 1, fp) != 1) {
        std::fclose(fp);
        return m;
    }
    m.nv = nv;
    m.nf = nf;
    m.v.resize(std::size_t{nv} * 3);
    m.f.resize(std::size_t{nf} * 3);
    const bool r1 = std::fread(m.v.data(), 8, m.v.size(), fp) == m.v.size();
    const bool r2 = std::fread(m.f.data(), 4, m.f.size(), fp) == m.f.size();
    std::fclose(fp);
    m.ok = r1 && r2;
    return m;
}

/// 剛体変換 + 一様スケール（`SPEC-phase5.md` §2.1 の「ランダムな変換」）。
///
/// **記録して再現できることが要件です。** 生成は `make_transform` が seed から行います。
/// **正規化のあとに掛けます。** 先に AABB を座標範囲へ合わせてから、
/// `shift`（正規化後の半径に対する割合）だけずらします。
/// **こうしないと 2 つの立体が完全に重なります**（§2.1 は「大きく重なる」であって
/// 「一致する」ではありません）。
struct Transform {
    double r[9];      ///< 回転（行優先）
    double shift[3];  ///< 正規化後のずらし量（半径に対する割合）
};

/// 量子化の結果。**分類に必要なものを全部返します**（`SPEC-phase5.md` §1）。
struct Quantized {
    krisite::mesh::TriMesh mesh;
    std::size_t dropped_degenerate = 0;  ///< 量子化で面積 0 になり落とした三角形
    std::size_t merged_vertices = 0;     ///< 同じ格子点に落ちて併合された頂点
    bool out_of_range = false;           ///< 座標が §2 の範囲を外れた（起きてはいけない）
};

/// `[lo, hi]` に丸める。**四捨五入**（切り捨てだと原点側に偏ります）。
inline std::int32_t quantize_one(double x) {
    const double lo = static_cast<double>(krisite::kCoordMin);
    const double hi = static_cast<double>(krisite::kCoordMax);
    const double r = std::nearbyint(x);
    return static_cast<std::int32_t>(r < lo ? lo : (r > hi ? hi : r));
}

/// 変換して量子化する。
///
/// **AABB を座標範囲いっぱいに伸ばします。** 入力のスケールはモデルごとにばらばらで、
/// そのまま丸めると分解能を捨てます（`SPEC-phase1.md` §2 の量子化と同じ考え方）。
///
/// **面積 0 の三角形は落とします。** `from_mesh` も落としますが、
/// **落とした数を数えるのはここです**（§1 の内訳に要る）。
inline Quantized quantize(const RawMesh& raw, const Transform& tr, double fill = 0.6) {
    Quantized q;
    if (!raw.ok || raw.nv == 0 || raw.nf == 0) return q;

    // ---- 1. 変換 ----
    std::vector<double> p(raw.v.size());
    for (std::size_t i = 0; i < raw.nv; ++i) {
        const double x = raw.v[3 * i], y = raw.v[3 * i + 1], z = raw.v[3 * i + 2];
        for (int k = 0; k < 3; ++k) {
            p[3 * i + k] = tr.r[3 * k] * x + tr.r[3 * k + 1] * y + tr.r[3 * k + 2] * z;
        }
    }

    // ---- 2. AABB を座標範囲に合わせる ----
    double lo[3] = {p[0], p[1], p[2]}, hi[3] = {p[0], p[1], p[2]};
    for (std::size_t i = 1; i < raw.nv; ++i) {
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::fmin(lo[k], p[3 * i + k]);
            hi[k] = std::fmax(hi[k], p[3 * i + k]);
        }
    }
    double ext = 0;
    for (int k = 0; k < 3; ++k) ext = std::fmax(ext, hi[k] - lo[k]);
    if (!(ext > 0)) return q;  // 退化した入力
    const double range = static_cast<double>(krisite::kCoordMax) * fill;
    const double scale = range / ext;

    // ---- 3. 量子化 + 同一格子点の併合 ----
    std::map<std::array<std::int32_t, 3>, std::uint32_t> at;
    std::vector<std::uint32_t> remap(raw.nv);
    for (std::size_t i = 0; i < raw.nv; ++i) {
        std::array<std::int32_t, 3> c{};
        for (int k = 0; k < 3; ++k) {
            const double mid = 0.5 * (lo[k] + hi[k]);
            c[static_cast<std::size_t>(k)] =
                quantize_one((p[3 * i + k] - mid) * scale + tr.shift[k] * range);
        }
        const auto it = at.find(c);
        if (it == at.end()) {
            const auto id = static_cast<std::uint32_t>(q.mesh.vertices.size());
            q.mesh.vertices.push_back({c[0], c[1], c[2]});
            at.emplace(c, id);
            remap[i] = id;
        } else {
            remap[i] = it->second;
            ++q.merged_vertices;
        }
    }
    for (const auto& v : q.mesh.vertices) {
        if (!krisite::geom::in_range(v)) q.out_of_range = true;
    }

    // ---- 4. 三角形。**頂点が潰れたものは落とします** ----
    for (std::size_t i = 0; i < raw.nf; ++i) {
        const std::uint32_t a = remap[raw.f[3 * i]], b = remap[raw.f[3 * i + 1]],
                            c = remap[raw.f[3 * i + 2]];
        if (a == b || b == c || c == a) {
            ++q.dropped_degenerate;
            continue;
        }
        q.mesh.triangles.push_back({a, b, c});
    }
    return q;
}

/// seed から変換を作る。**同じ seed なら同じ変換**（§2.1 の再現性）。
///
/// 回転は四元数から。**軸と角を独立に振ると分布が偏る**ので、
/// 単位四元数を正規分布から作ります（Marsaglia）。
inline Transform make_transform(std::uint64_t seed) {
    // xorshift64*。**標準の乱数エンジンは実装差があるので使いません**（再現性のため）
    auto next = [&seed]() {
        seed ^= seed >> 12;
        seed ^= seed << 25;
        seed ^= seed >> 27;
        return (seed * 0x2545F4914F6CDD1DULL) >> 11;
    };
    auto uniform = [&next]() { return static_cast<double>(next()) / 9007199254740992.0; };
    auto normal = [&uniform]() {
        // Box-Muller。**0 を避ける**
        const double u1 = std::fmax(uniform(), 1e-12), u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    };
    double qt[4];
    double n = 0;
    for (double& x : qt) {
        x = normal();
        n += x * x;
    }
    n = std::sqrt(n);
    if (!(n > 0)) {
        qt[0] = 1;
        qt[1] = qt[2] = qt[3] = 0;
        n = 1;
    }
    for (double& x : qt) x /= n;
    const double w = qt[0], x = qt[1], y = qt[2], z = qt[3];

    Transform t{};
    t.r[0] = 1 - 2 * (y * y + z * z);
    t.r[1] = 2 * (x * y - z * w);
    t.r[2] = 2 * (x * z + y * w);
    t.r[3] = 2 * (x * y + z * w);
    t.r[4] = 1 - 2 * (x * x + z * z);
    t.r[5] = 2 * (y * z - x * w);
    t.r[6] = 2 * (x * z - y * w);
    t.r[7] = 2 * (y * z + x * w);
    t.r[8] = 1 - 2 * (x * x + y * y);
    // **ずらし量は半径の 35% まで。** `fill = 0.6` と合わせて 0.6 * 1.35 = 0.81 なので
    // 座標範囲に収まります。**大きく重なりつつ、はみ出る部分も残ります**
    for (double& d : t.shift) d = (uniform() * 2.0 - 1.0) * 0.35;
    return t;
}

}  // namespace krithingi

#endif  // KRISITE_TESTS_THINGI10K_LOADER_HPP

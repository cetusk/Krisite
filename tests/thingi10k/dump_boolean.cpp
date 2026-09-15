// Krisite — 4 つのブール出力を厳密に書き出す（案 I。`DESIGN-phase5-vertex-level.md` §44.8）
//
// **これは診断の道具です。** 被検体（`include/krisite/`）を**呼びます**が、
// **正解器ではありません。** 出力を保存するだけで、判定はしません。
//
// **★ 受理済みの量子化配列を【読み】、SHA256 を照合し、【その同じ配列を被検体へ渡します】。**
// **再量子化しません**（§44.6-1）。**別のプログラムが配列を再現しても、
// この駆動が実際に使った配列が同じことにはならないため**です。
//
// **GMP には依存しません**（`hash_mesh` 相当の計算は自前で行います。§44.8）。
//
// 使い方:
//   dump_boolean <A.bin> <A の sha256> <B.bin> <B の sha256> <出力の接頭辞> \
//                <depth> <threads> <nsi_mode> <verify_delta> \
//                <fine_k> <fine_budget> <fine_mb> <skip_disjoint> <skip_boxside> <repair>
//
// **保存された実行の設定**（`cp3_gmp.meta` の `args` から。§44.8）:
//   depth=6 threads=8 nsi_mode=1 verify_delta=1 fine_k=0 fine_budget=0
//   fine_mb=16 skip_disjoint=2 skip_boxside=1 repair=1
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "krisite/csg/polysoup.hpp"
#include "krisite/csg/soup_boolean.hpp"
#include "krisite/csg/to_mesh.hpp"
#include "krisite/mesh/tri_mesh.hpp"
#include "krisite/par/thread_pool.hpp"

#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__)
#error "この道具はリトルエンディアン固定です（書式の取り決め。§44.8）"
#endif

namespace kri = krisite;
namespace csg = krisite::csg;
namespace geom = krisite::geom;
namespace mesh = krisite::mesh;
namespace par = krisite::par;

namespace {

// ---- SHA-256（FIPS 180-4。アルゴリズムなので著作権の対象外。自前実装）----
struct Sha256 {
    std::uint32_t h[8] = {0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
                          0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
    std::uint8_t buf[64]{};
    std::size_t n = 0;
    std::uint64_t total = 0;

    static std::uint32_t ror(std::uint32_t x, int r) { return (x >> r) | (x << (32 - r)); }

    void block(const std::uint8_t* p) {
        static const std::uint32_t k[64] = {
            0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
            0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
            0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
            0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
            0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
            0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
            0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
            0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
            0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
            0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
            0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
        std::uint32_t w[64];
        for (int i = 0; i < 16; ++i) {
            w[i] = (std::uint32_t(p[4 * i]) << 24) | (std::uint32_t(p[4 * i + 1]) << 16) |
                   (std::uint32_t(p[4 * i + 2]) << 8) | std::uint32_t(p[4 * i + 3]);
        }
        for (int i = 16; i < 64; ++i) {
            const std::uint32_t s0 = ror(w[i - 15], 7) ^ ror(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const std::uint32_t s1 = ror(w[i - 2], 17) ^ ror(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        std::uint32_t a = h[0], b = h[1], c = h[2], d = h[3];
        std::uint32_t e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            const std::uint32_t S1 = ror(e, 6) ^ ror(e, 11) ^ ror(e, 25);
            const std::uint32_t ch = (e & f) ^ (~e & g);
            const std::uint32_t t1 = hh + S1 + ch + k[i] + w[i];
            const std::uint32_t S0 = ror(a, 2) ^ ror(a, 13) ^ ror(a, 22);
            const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const std::uint32_t t2 = S0 + maj;
            hh = g; g = f; f = e; e = d + t1;
            d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d;
        h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }

    void update(const void* data, std::size_t len) {
        const auto* p = static_cast<const std::uint8_t*>(data);
        total += len;
        while (len > 0) {
            const std::size_t take = (64 - n < len) ? (64 - n) : len;
            std::memcpy(buf + n, p, take);
            n += take; p += take; len -= take;
            if (n == 64) { block(buf); n = 0; }
        }
    }

    std::string hex() {
        const std::uint64_t bits = total * 8;
        std::uint8_t pad = 0x80;
        update(&pad, 1);
        pad = 0;
        while (n != 56) update(&pad, 1);
        std::uint8_t len8[8];
        for (int i = 0; i < 8; ++i) len8[i] = std::uint8_t(bits >> (56 - 8 * i));
        update(len8, 8);  // `bits` は先に取ってあるので、`total` の増加は結果に効きません
        char out[65];
        for (int i = 0; i < 8; ++i) std::snprintf(out + 8 * i, 9, "%08x", h[i]);
        return std::string(out, 64);
    }
};

std::string sha256_hex(const std::vector<std::uint8_t>& b) {
    Sha256 s;
    s.update(b.data(), b.size());
    return s.hex();
}

// ---- 量子化済み配列の読み込み（`dump_quantized` が書いた書式）----
struct Loaded {
    mesh::TriMesh mesh;
    std::string sha;
    bool ok = false;
};

Loaded read_quantized(const std::string& path) {
    Loaded r;
    std::FILE* fp = std::fopen(path.c_str(), "rb");
    if (fp == nullptr) return r;
    std::fseek(fp, 0, SEEK_END);
    const long sz = std::ftell(fp);
    std::fseek(fp, 0, SEEK_SET);
    std::vector<std::uint8_t> b(static_cast<std::size_t>(sz < 0 ? 0 : sz));
    const bool rd = !b.empty() && std::fread(b.data(), 1, b.size(), fp) == b.size();
    std::fclose(fp);
    if (!rd || b.size() < 8) return r;
    std::uint32_t nv = 0, nf = 0;
    std::memcpy(&nv, b.data(), 4);
    std::memcpy(&nf, b.data() + 4, 4);
    const std::size_t need = 8 + std::size_t{nv} * 12 + std::size_t{nf} * 12;
    if (b.size() != need) return r;
    r.mesh.vertices.resize(nv);
    for (std::uint32_t i = 0; i < nv; ++i) {
        std::int32_t c[3];
        std::memcpy(c, b.data() + 8 + std::size_t{i} * 12, 12);
        r.mesh.vertices[i] = {c[0], c[1], c[2]};
    }
    r.mesh.triangles.resize(nf);
    const std::size_t off = 8 + std::size_t{nv} * 12;
    for (std::uint32_t i = 0; i < nf; ++i) {
        std::uint32_t t[3];
        std::memcpy(t, b.data() + off + std::size_t{i} * 12, 12);
        r.mesh.triangles[i] = {t[0], t[1], t[2]};
    }
    r.sha = sha256_hex(b);
    r.ok = true;
    return r;
}

// ---- 出力の書き出し（§44.8 の書式）----
//
// リトルエンディアン固定。`uint32 Nv` / `uint32 Nf` /
// 頂点ごとに x の kHomoXyz リム、y の kHomoXyz リム、z の kHomoXyz リム、
// w の kHomoW リム（いずれも `uint64` を little-endian で） / `uint32` x 3Nf（面）。
// **マジックも版も余白も付けません。**
bool write_soup(const std::string& path, const csg::SoupMesh& m) {
    std::FILE* fp = std::fopen(path.c_str(), "wb");
    if (fp == nullptr) return false;
    const auto nv = static_cast<std::uint32_t>(m.vertices.size());
    const auto nf = static_cast<std::uint32_t>(m.triangles.size());
    bool ok = std::fwrite(&nv, 4, 1, fp) == 1 && std::fwrite(&nf, 4, 1, fp) == 1;
    std::vector<std::uint64_t> row(3 * geom::limbs::kHomoXyz + geom::limbs::kHomoW);
    for (const geom::HPointD& v : m.vertices) {
        std::size_t k = 0;
        for (std::size_t l = 0; l < geom::limbs::kHomoXyz; ++l) row[k++] = v.x[l];
        for (std::size_t l = 0; l < geom::limbs::kHomoXyz; ++l) row[k++] = v.y[l];
        for (std::size_t l = 0; l < geom::limbs::kHomoXyz; ++l) row[k++] = v.z[l];
        for (std::size_t l = 0; l < geom::limbs::kHomoW; ++l) row[k++] = v.w[l];
        ok = ok && std::fwrite(row.data(), 8, row.size(), fp) == row.size();
    }
    std::vector<std::uint32_t> fs(std::size_t{nf} * 3);
    for (std::size_t i = 0; i < nf; ++i) {
        fs[3 * i] = m.triangles[i][0];
        fs[3 * i + 1] = m.triangles[i][1];
        fs[3 * i + 2] = m.triangles[i][2];
    }
    ok = ok && (fs.empty() || std::fwrite(fs.data(), 4, fs.size(), fp) == fs.size());
    return (std::fclose(fp) == 0) && ok;
}

/// `thingi_cp1.cpp:285-300` の `hash_mesh` と**同じ計算**。**GMP は呼びません。**
unsigned long long hash_mesh(const csg::SoupMesh& m) {
    unsigned long long h = 1469598103934665603ull;
    const auto mix = [&h](unsigned long long v) {
        h ^= v;
        h *= 1099511628211ull;
    };
    mix(m.vertices.size());
    mix(m.triangles.size());
    for (const geom::HPointD& v : m.vertices) {
        for (std::size_t l = 0; l < geom::limbs::kHomoXyz; ++l) {
            mix(v.x[l]);
            mix(v.y[l]);
            mix(v.z[l]);
        }
        for (std::size_t l = 0; l < geom::limbs::kHomoW; ++l) mix(v.w[l]);
    }
    for (const mesh::Tri& t : m.triangles) {
        mix(t[0]);
        mix(t[1]);
        mix(t[2]);
    }
    return h;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 16) {
        std::fprintf(stderr,
                     "使い方: dump_boolean <A.bin> <A の sha256> <B.bin> <B の sha256> "
                     "<出力の接頭辞> <depth> <threads> <nsi_mode> <verify_delta> "
                     "<fine_k> <fine_budget> <fine_mb> <skip_disjoint> <skip_boxside> <repair>\n");
        return 2;
    }
    const std::string pa = argv[1], sa = argv[2], pb = argv[3], sb = argv[4], prefix = argv[5];
    const unsigned depth = static_cast<unsigned>(std::atoi(argv[6]));
    const unsigned nthreads = static_cast<unsigned>(std::atoi(argv[7]));
    const int nsi_mode = std::atoi(argv[8]);
    const bool verify_delta = std::atoi(argv[9]) != 0;
    const std::size_t fine_k = std::strtoul(argv[10], nullptr, 10);
    const std::size_t fine_budget = std::strtoul(argv[11], nullptr, 10);
    const std::size_t fine_mb = std::strtoul(argv[12], nullptr, 10);
    const int skip_disjoint = std::atoi(argv[13]);
    const bool skip_boxside = std::atoi(argv[14]) != 0;
    const bool repair = std::atoi(argv[15]) != 0;

    // ---- 設定と対象を先に全部出す（`CLAUDE.md`）----
    std::printf("dump_boolean\n");
    std::printf("  入力 A %s（期待 sha256 %s）\n", pa.c_str(), sa.c_str());
    std::printf("  入力 B %s（期待 sha256 %s）\n", pb.c_str(), sb.c_str());
    std::printf("  出力の接頭辞 %s\n", prefix.c_str());
    std::printf("  設定: depth=%u threads=%u nsi_mode=%d verify_delta=%d fine_k=%zu "
                "fine_budget=%zu fine_mb=%zu skip_disjoint=%d skip_boxside=%d repair=%d\n",
                depth, nthreads, nsi_mode, verify_delta ? 1 : 0, fine_k, fine_budget, fine_mb,
                skip_disjoint, skip_boxside ? 1 : 0, repair ? 1 : 0);
    std::printf("  b = %zu / kHomoXyz = %zu リム / kHomoW = %zu リム\n", kri::kCoordBits,
                geom::limbs::kHomoXyz, geom::limbs::kHomoW);
    std::printf("  **再量子化しません。読んだ配列をそのまま渡します。**\n");
    std::fflush(stdout);

    // ---- 入力を読み、SHA256 を照合する（**違えば起動しません**）----
    const Loaded A0 = read_quantized(pa);
    const Loaded B0 = read_quantized(pb);
    if (!A0.ok || !B0.ok) {
        std::fprintf(stderr, "入力を読めません\n");
        return 3;
    }
    std::printf("  A 読込: 頂点 %zu / 面 %zu / sha256 %s / %s\n", A0.mesh.vertices.size(),
                A0.mesh.triangles.size(), A0.sha.c_str(), A0.sha == sa ? "一致" : "★ 不一致");
    std::printf("  B 読込: 頂点 %zu / 面 %zu / sha256 %s / %s\n", B0.mesh.vertices.size(),
                B0.mesh.triangles.size(), B0.sha.c_str(), B0.sha == sb ? "一致" : "★ 不一致");
    if (A0.sha != sa || B0.sha != sb) {
        std::fprintf(stderr, "入力のハッシュが違います。起動しません。\n");
        return 4;
    }

    // ---- 被検体（`thingi_cp1.cpp` の `check_one` と同じ組み立て）----
    csg::FromMeshOptions fm;
    fm.verify_nsi = (nsi_mode == 1);
    csg::PolySoup A = csg::from_mesh(A0.mesh, fm), B = csg::from_mesh(B0.mesh, fm);
    if (nsi_mode == 2) {
        A.nsi.assign(A.sources.size(), 1);
        B.nsi.assign(B.sources.size(), 1);
    }

    par::ThreadPool pool(nthreads);
    csg::BoolOptions o;
    o.measure_classify = true;
    o.measure_frag = true;
    o.record_ray_levels = true;
    o.ray_index_fine_cells = fine_k;
    o.ray_index_fine_budget = fine_budget;
    o.ray_index_fine_bytes = fine_mb << 20;
    o.bsp_skip_disjoint = skip_disjoint;
    o.bsp_skip_boxside = skip_boxside;
    o.depth = depth;
    o.adaptive = true;
    o.leaf_threshold = 0;
    o.cull_planes = true;
    o.early_out = true;
    o.cache_points = true;
    o.local_bsp = true;
    o.split_contacts = true;
    o.threads = nthreads;
    o.pool = &pool;

    csg::ToMeshOptions tm;
    tm.split_contacts = true;
    tm.threads = pool.size();
    tm.pool = &pool;
    tm.verify_split_delta = verify_delta;
    tm.diag_unresolved = false;
    tm.repair_unresolved = repair;

    const char* name[4] = {"union", "isect", "diff_ab", "diff_ba"};
    csg::SoupMesh out[4];
    const auto t0 = std::chrono::steady_clock::now();
    int k = 0;
    for (csg::BoolOp op :
         {csg::BoolOp::Union, csg::BoolOp::Intersection, csg::BoolOp::Difference}) {
        const csg::PolySoup s = csg::boolean(A, B, op, o, nullptr);
        out[k] = csg::to_mesh(s, tm, nullptr);
        std::printf("  %s: 三角形 %zu / 頂点 %zu\n", name[k], out[k].triangles.size(),
                    out[k].vertices.size());
        std::fflush(stdout);
        ++k;
    }
    {
        const csg::PolySoup s4 = csg::boolean(B, A, csg::BoolOp::Difference, o, nullptr);
        out[3] = csg::to_mesh(s4, tm, nullptr);
        std::printf("  %s: 三角形 %zu / 頂点 %zu\n", name[3], out[3].triangles.size(),
                    out[3].vertices.size());
    }
    const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // ---- 4 演算の合成ハッシュ（`thingi_cp1.cpp:1180` と同じ式）----
    const unsigned long long hash4 = hash_mesh(out[0]) ^ (hash_mesh(out[1]) * 3) ^
                                     (hash_mesh(out[2]) * 7) ^ (hash_mesh(out[3]) * 11);

    for (int i = 0; i < 4; ++i) {
        const std::string path = prefix + "_out_" + name[i] + ".bin";
        if (!write_soup(path, out[i])) {
            std::fprintf(stderr, "書き込みに失敗しました: %s\n", path.c_str());
            return 5;
        }
    }

    // ---- 演算別の識別情報 ----
    {
        const std::string path = prefix + "_run_meta.txt";
        std::FILE* fp = std::fopen(path.c_str(), "w");
        if (fp == nullptr) {
            std::fprintf(stderr, "書けません: %s\n", path.c_str());
            return 5;
        }
        std::fprintf(fp, "b=%zu\nkHomoXyz=%zu\nkHomoW=%zu\n", kri::kCoordBits,
                     geom::limbs::kHomoXyz, geom::limbs::kHomoW);
        std::fprintf(fp, "in_a=%s\nin_a_sha256=%s\nin_b=%s\nin_b_sha256=%s\n", pa.c_str(),
                     A0.sha.c_str(), pb.c_str(), B0.sha.c_str());
        std::fprintf(fp,
                     "depth=%u\nthreads=%u\nnsi_mode=%d\nverify_delta=%d\nfine_k=%zu\n"
                     "fine_budget=%zu\nfine_mb=%zu\nskip_disjoint=%d\nskip_boxside=%d\nrepair=%d\n",
                     depth, nthreads, nsi_mode, verify_delta ? 1 : 0, fine_k, fine_budget, fine_mb,
                     skip_disjoint, skip_boxside ? 1 : 0, repair ? 1 : 0);
        std::fprintf(fp, "requantized=0\n");
        for (int i = 0; i < 4; ++i) {
            std::fprintf(fp, "op%d=%s tri=%zu vert=%zu hash=%llu\n", i, name[i],
                         out[i].triangles.size(), out[i].vertices.size(), hash_mesh(out[i]));
        }
        std::fprintf(fp, "hash4=%016llx\n", hash4);
        std::fprintf(fp, "seconds=%.6f\n", secs);
        if (std::fclose(fp) != 0) {
            std::fprintf(stderr, "書き込みに失敗しました: %s\n", path.c_str());
            return 5;
        }
    }

    std::printf("  4 演算ハッシュ %016llx / 経過 %.3f 秒\n", hash4, secs);
    std::printf("  面数 %zu %zu %zu %zu\n", out[0].triangles.size(), out[1].triangles.size(),
                out[2].triangles.size(), out[3].triangles.size());
    return 0;
}

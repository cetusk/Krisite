// Krisite — 量子化済みの入力を書き出す（経路 C。`DESIGN-phase5-vertex-level.md` §43.10）
//
// **これは診断の道具です。ブール演算も分類器も述語も呼びません。**
// 既存のローダー（`loader.hpp`）の `load_kmesh` / `make_transform` / `quantize` を
// そのまま通し、**結果と、使った変換の 12 個の `double` を書き出します。**
//
// **変換の【生成】は経路 P では再現しません**（`std::log` / `std::cos` / `std::sqrt` を
// 通るため、libm の実装差で最下位ビットが違い得ます）。**だからここが出す変換を
// 経路 P に渡し、経路 P は【量子化と併合だけ】を独立に実装します**（§43.10）。
//
// 使い方:
//   dump_quantized <入力 .kmesh> <seed> <出力の接頭辞>
// 書くもの:
//   <接頭辞>_transform.hex   r[0..8] と shift[0..2] を 1 行 1 個の %a（C99 の 16 進浮動小数点）
//   <接頭辞>_quantized.bin   uint32 Nv / uint32 Nf / int32 x 3Nv / uint32 x 3Nf
//                            **マジックも版も余白も付けません。リトルエンディアン固定。**
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "loader.hpp"

#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__)
#error "この道具はリトルエンディアン固定です（書式の取り決め。§43.10）"
#endif

namespace {

/// **書けたかを必ず確かめます。** 上流の失敗を、次の段の前提にしないため。
bool write_all(std::FILE* fp, const void* p, std::size_t n) {
    return n == 0 || std::fwrite(p, 1, n, fp) == n;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "使い方: dump_quantized <入力 .kmesh> <seed> <出力の接頭辞>\n");
        return 2;
    }
    const std::string in = argv[1];
    const std::uint64_t seed = std::strtoull(argv[2], nullptr, 10);
    const std::string prefix = argv[3];

    const krithingi::RawMesh raw = krithingi::load_kmesh(in);
    if (!raw.ok) {
        std::fprintf(stderr, "読み込みに失敗しました: %s\n", in.c_str());
        return 3;
    }

    const krithingi::Transform tr = krithingi::make_transform(seed);

    // ---- 変換の 12 個の double（経路 P へ渡す）----
    {
        const std::string path = prefix + "_transform.hex";
        std::FILE* fp = std::fopen(path.c_str(), "w");
        if (fp == nullptr) {
            std::fprintf(stderr, "書けません: %s\n", path.c_str());
            return 4;
        }
        bool ok = true;
        for (int i = 0; i < 9; ++i) ok = ok && std::fprintf(fp, "%a\n", tr.r[i]) > 0;
        for (int i = 0; i < 3; ++i) ok = ok && std::fprintf(fp, "%a\n", tr.shift[i]) > 0;
        ok = std::fclose(fp) == 0 && ok;
        if (!ok) {
            std::fprintf(stderr, "書き込みに失敗しました: %s\n", path.c_str());
            return 4;
        }
    }

    const krithingi::Quantized q = krithingi::quantize(raw, tr);

    // ---- 量子化済みの配列 ----
    {
        const std::string path = prefix + "_quantized.bin";
        std::FILE* fp = std::fopen(path.c_str(), "wb");
        if (fp == nullptr) {
            std::fprintf(stderr, "書けません: %s\n", path.c_str());
            return 5;
        }
        const auto nv = static_cast<std::uint32_t>(q.mesh.vertices.size());
        const auto nf = static_cast<std::uint32_t>(q.mesh.triangles.size());
        std::vector<std::int32_t> vs(std::size_t{nv} * 3);
        for (std::size_t i = 0; i < nv; ++i) {
            vs[3 * i] = q.mesh.vertices[i].x;
            vs[3 * i + 1] = q.mesh.vertices[i].y;
            vs[3 * i + 2] = q.mesh.vertices[i].z;
        }
        std::vector<std::uint32_t> fs(std::size_t{nf} * 3);
        for (std::size_t i = 0; i < nf; ++i) {
            fs[3 * i] = q.mesh.triangles[i][0];
            fs[3 * i + 1] = q.mesh.triangles[i][1];
            fs[3 * i + 2] = q.mesh.triangles[i][2];
        }
        bool ok = write_all(fp, &nv, 4) && write_all(fp, &nf, 4) &&
                  write_all(fp, vs.data(), vs.size() * 4) &&
                  write_all(fp, fs.data(), fs.size() * 4);
        ok = std::fclose(fp) == 0 && ok;
        if (!ok) {
            std::fprintf(stderr, "書き込みに失敗しました: %s\n", path.c_str());
            return 5;
        }
    }

    // ---- 設定と、実際に扱った対象を出力に書く ----
    std::printf("入力 %s / seed %llu / 出力の接頭辞 %s\n", in.c_str(),
                static_cast<unsigned long long>(seed), prefix.c_str());
    std::printf("原本 頂点 %zu / 面 %zu\n", raw.nv, raw.nf);
    std::printf("量子化後 頂点 %zu / 面 %zu / 併合した頂点 %zu / 落とした退化面 %zu / 範囲外 %d\n",
                q.mesh.vertices.size(), q.mesh.triangles.size(), q.merged_vertices,
                q.dropped_degenerate, q.out_of_range ? 1 : 0);
    std::printf("b = %zu / kCoordMin = %lld / kCoordMax = %lld / fill = 0.6\n", krisite::kCoordBits,
                static_cast<long long>(krisite::kCoordMin),
                static_cast<long long>(krisite::kCoordMax));
    return 0;
}

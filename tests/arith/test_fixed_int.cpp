// Krisite — fixed_int の性質テスト（正解器なし）
//
// SPEC-phase0.md §8.3 に相当する算術側の不変条件。
// GMP を必要としないので既定のビルドで常に実行される。
#include <cstdint>

#include "test_util.hpp"

using namespace krisite::arith;
using kritest::Rng;

namespace {

using i128 = __int128;
using u128 = unsigned __int128;

fixed_int<2> from_i128(i128 v) noexcept {
    fixed_int<2> r{};
    const u128 u = static_cast<u128>(v);
    r.limb[0] = static_cast<std::uint64_t>(u);
    r.limb[1] = static_cast<std::uint64_t>(u >> 64);
    return r;
}

i128 to_i128(const fixed_int<2>& x) noexcept {
    return static_cast<i128>((static_cast<u128>(x.limb[1]) << 64) | x.limb[0]);
}

std::string i128_str(i128 v) {
    if (v == 0) return "0";
    const bool neg = v < 0;
    u128 u = neg ? (~static_cast<u128>(v) + 1) : static_cast<u128>(v);
    std::string s;
    while (u != 0) {
        s.insert(s.begin(), static_cast<char>('0' + static_cast<int>(u % 10)));
        u /= 10;
    }
    return neg ? "-" + s : s;
}

// ---- 基本 -------------------------------------------------------------------

void test_basics() {
    KRI_CHECK(sign(zero<4>()) == 0);
    KRI_CHECK(is_zero(zero<4>()));
    KRI_CHECK(sign(from_i64<4>(1)) == 1);
    KRI_CHECK(sign(from_i64<4>(-1)) == -1);
    KRI_CHECK(is_negative(from_i64<3>(-1)));
    KRI_CHECK(!is_negative(from_i64<3>(0)));

    // min_bits: -2^(n-1) <= x <= 2^(n-1)-1 を満たす最小の n
    KRI_CHECK(min_bits(from_i64<2>(0)) == 1);
    KRI_CHECK(min_bits(from_i64<2>(-1)) == 1);
    KRI_CHECK(min_bits(from_i64<2>(1)) == 2);
    KRI_CHECK(min_bits(from_i64<2>(-2)) == 2);
    KRI_CHECK(min_bits(from_i64<2>(2)) == 3);
    KRI_CHECK(min_bits(from_i64<2>(127)) == 8);
    KRI_CHECK(min_bits(from_i64<2>(-128)) == 8);
    KRI_CHECK(min_bits(from_i64<2>(128)) == 9);

    // 符号拡張して幅を変えても値は変わらない
    Rng rng(1);
    for (int i = 0; i < 10000; ++i) {
        const fixed_int<2> a = kritest::rand_biased<2>(rng);
        KRI_CHECK(cmp(a, widen<5>(a)) == 0);
        KRI_CHECK(cmp(a, resize<2>(widen<5>(a))) == 0);
        KRI_CHECK(min_bits(a) == min_bits(widen<5>(a)));
    }
}

// ---- __int128 との突き合わせ ------------------------------------------------

void test_against_int128() {
    Rng rng(20240826);
    for (int iter = 0; iter < 200000; ++iter) {
        const fixed_int<2> a = kritest::rand_biased<2>(rng);
        const fixed_int<2> b = kritest::rand_biased<2>(rng);
        const i128 ia = to_i128(a), ib = to_i128(b);

        // cmp
        const int want_cmp = (ia < ib) ? -1 : (ia > ib) ? 1 : 0;
        KRI_CHECK_MSG(cmp(a, b) == want_cmp, i128_str(ia) + " vs " + i128_str(ib));

        // sign
        const int want_sign = (ia < 0) ? -1 : (ia > 0) ? 1 : 0;
        KRI_CHECK(sign(a) == want_sign);

        // add / sub は 3 リムに広げれば必ず収まる
        const i128 hi_a = ia >> 1, hi_b = ib >> 1;  // オーバーフローしない範囲に落とす
        const fixed_int<2> ha = from_i128(hi_a), hb = from_i128(hi_b);
        KRI_CHECK_MSG(to_i128(add(ha, hb)) == hi_a + hi_b, i128_str(hi_a) + " + " + i128_str(hi_b));
        KRI_CHECK_MSG(to_i128(sub(ha, hb)) == hi_a - hi_b, i128_str(hi_a) + " - " + i128_str(hi_b));

        // add_widen / sub_widen は入力の全域でオーバーフローしない
        KRI_CHECK(cmp(sub_mixed(add_widen(a, b), b), a) == 0);
        KRI_CHECK(cmp(add_mixed(sub_widen(a, b), b), a) == 0);

        // シフト
        const std::size_t s = static_cast<std::size_t>(rng.below(128));
        KRI_CHECK_MSG(to_i128(shl_bits(a, s)) == static_cast<i128>(static_cast<u128>(ia) << s),
                      i128_str(ia) + " << " + std::to_string(s));
        KRI_CHECK_MSG(to_i128(shr_bits(a, s)) == (ia >> s),
                      i128_str(ia) + " >> " + std::to_string(s));
    }
}

void test_mul_against_int128() {
    Rng rng(7);
    for (int iter = 0; iter < 200000; ++iter) {
        const fixed_int<1> a = kritest::rand_biased<1>(rng);
        const fixed_int<1> b = kritest::rand_biased<1>(rng);
        const auto ia = static_cast<std::int64_t>(a.limb[0]);
        const auto ib = static_cast<std::int64_t>(b.limb[0]);
        const i128 want = static_cast<i128>(ia) * static_cast<i128>(ib);
        const fixed_int<2> got = mul(a, b);
        KRI_CHECK_MSG(to_i128(got) == want, std::to_string(ia) + " * " + std::to_string(ib) +
                                                " -> " + i128_str(to_i128(got)) + " 期待 " +
                                                i128_str(want));
    }

    // 最小値 x 最小値 が N+M リムにちょうど収まること（SPEC §5.2 の最悪ケース）
    fixed_int<1> mn{};
    mn.limb[0] = std::uint64_t{1} << 63;
    const fixed_int<2> sq = mul(mn, mn);
    // (-2^63)^2 = 2^126
    KRI_CHECK(sq.limb[0] == 0);
    KRI_CHECK(sq.limb[1] == (std::uint64_t{1} << 62));
    KRI_CHECK(sign(sq) == 1);
}

// ---- 代数的な恒等式 ---------------------------------------------------------

void test_identities() {
    Rng rng(99);
    for (int iter = 0; iter < 50000; ++iter) {
        const fixed_int<3> a = kritest::rand_biased<3>(rng);
        const fixed_int<3> b = kritest::rand_biased<3>(rng);
        const fixed_int<3> c = kritest::rand_biased<3>(rng);

        // (a + b) - b == a （幅を広げて必ず成立させる）
        const fixed_int<4> ab = add_widen(a, b);
        KRI_CHECK(cmp(sub_mixed(ab, b), a) == 0);

        // 交換則
        KRI_CHECK(cmp(mul(a, b), mul(b, a)) == 0);

        // 分配則: a*(b+c) == a*b + a*c
        const fixed_int<3 + 4> lhs = mul(a, add_widen(b, c));
        const fixed_int<3 + 3 + 1> rhs = add_mixed(mul(a, b), mul(a, c));
        KRI_CHECK(cmp(lhs, rhs) == 0);

        // 結合則
        KRI_CHECK(cmp(mul(mul(a, b), c), mul(a, mul(b, c))) == 0);

        // 符号
        KRI_CHECK(sign(mul(a, b)) == sign(a) * sign(b));

        // neg(neg(a)) == a （最小値以外）
        if (min_bits(a) < 64 * 3) {
            KRI_CHECK(cmp(neg(neg(a)), a) == 0);
            KRI_CHECK(sign(neg(a)) == -sign(a));
        }

        // 反対称律・推移律（cmp の全順序性）
        KRI_CHECK(cmp(a, b) == -cmp(b, a));
        if (cmp(a, b) <= 0 && cmp(b, c) <= 0) KRI_CHECK(cmp(a, c) <= 0);

        // min_bits は乗算で高々和になる
        KRI_CHECK(min_bits(mul(a, b)) <= min_bits(a) + min_bits(b));
    }
}

// ---- 行列式 -----------------------------------------------------------------

void test_determinants() {
    Rng rng(1234);
    for (int iter = 0; iter < 20000; ++iter) {
        fixed_int<1> m[3][3];
        for (auto& row : m)
            for (auto& e : row) e = from_i64<1>(rng.range(-100000, 100000));

        // 2 行が同じなら 0
        {
            fixed_int<1> q[3][3];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) q[i][j] = m[i][j];
            q[1][0] = q[0][0];
            q[1][1] = q[0][1];
            q[1][2] = q[0][2];
            KRI_CHECK(sign(det3(q)) == 0);
        }

        // 行の入れ替えで符号反転
        {
            fixed_int<1> q[3][3];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) q[i][j] = m[i][j];
            for (int j = 0; j < 3; ++j) {
                q[0][j] = m[1][j];
                q[1][j] = m[0][j];
            }
            KRI_CHECK(sign(det3(q)) == -sign(det3(m)));
        }

        // det2 の反対称性
        const auto a = from_i64<1>(rng.range(-100000, 100000));
        const auto b = from_i64<1>(rng.range(-100000, 100000));
        const auto c = from_i64<1>(rng.range(-100000, 100000));
        const auto d = from_i64<1>(rng.range(-100000, 100000));
        KRI_CHECK(sign(det2(a, b, c, d)) == -sign(det2(b, a, d, c)));
        KRI_CHECK(sign(det2(a, b, a, b)) == 0);
    }
}

}  // namespace

int main() {
    test_basics();
    test_against_int128();
    test_mul_against_int128();
    test_identities();
    test_determinants();
    return kritest::finish("arith/fixed_int");
}

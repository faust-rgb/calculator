// ============================================================================
// 大整数乘法算法实现
// ============================================================================
//
// 本文件实现高效的大整数乘法算法：
// - 朴素乘法 O(n^2)
// - Karatsuba 算法 O(n^1.585)
// - Toom-Cook 3 算法 O(n^1.465)
// - NTT（数论变换）算法 O(n log n)
// ============================================================================

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

#include "precise_decimal.h"

#include <algorithm>
#include <mutex>

namespace precise {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;  ///< 大整数基数（10^9）
constexpr int kBaseDigits = 9;          ///< 每个 chunk 的十进制位数

constexpr std::size_t KARATSUBA_THRESHOLD = 128;
constexpr std::size_t NTT_THRESHOLD = 2048;

// ============================================================================
// 辅助函数
// ============================================================================

BigIntData bigint_low(const BigIntData& v, std::size_t n) {
    if (n >= v.size()) return v;
    return BigIntData(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(n));
}

BigIntData bigint_high(const BigIntData& v, std::size_t n) {
    if (n >= v.size()) return {0};
    return BigIntData(v.begin() + static_cast<std::ptrdiff_t>(n), v.end());
}

void trim_bigint(BigIntData& v) {
    while (v.size() > 1 && v.back() == 0) v.pop_back();
}

BigIntData multiply_bigint_by_uint32(const BigIntData& v, uint32_t n) {
    if (n == 0 || (v.size() == 1 && v[0] == 0)) return {0};
    if (n == 1) return v;
    BigIntData res;
    res.resize(v.size() + 1);
    size_t new_size = mul_uint32_raw(v.data(), v.size(), n, res.data());
    res.resize(new_size);
    return res;
}

BigIntData shift_bigint(const BigIntData& v, std::size_t n) {
    if (v.empty() || (v.size() == 1 && v[0] == 0)) return {0};
    if (n == 0) return v;
    BigIntData res(n, 0);
    res.insert(res.end(), v.begin(), v.end());
    return res;
}

// 前向声明
BigIntData multiply_bigint_karatsuba(const BigIntData& lhs, const BigIntData& rhs);
BigIntData multiply_bigint_toom3(const BigIntData& lhs, const BigIntData& rhs);
BigIntData multiply_bigint_ntt(const BigIntData& lhs, const BigIntData& rhs);

// ============================================================================
// 朴素乘法
// ============================================================================

BigIntData multiply_bigint_naive(const BigIntData& lhs, const BigIntData& rhs) {
    if ((lhs.size() == 1 && lhs[0] == 0) || (rhs.size() == 1 && rhs[0] == 0)) return {0};

    BigIntData res(lhs.size() + rhs.size(), 0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i] == 0) continue;

        uint64_t carry = 0;
        uint64_t l_val = lhs[i];
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            uint64_t cur = res[i + j] + l_val * rhs[j] + carry;
            res[i + j] = static_cast<uint32_t>(cur % kBase);
            carry = cur / kBase;
        }
        if (carry) {
            std::size_t k = i + rhs.size();
            while (carry && k < res.size()) {
                uint64_t cur = res[k] + carry;
                res[k] = static_cast<uint32_t>(cur % kBase);
                carry = cur / kBase;
                k++;
            }
        }
    }

    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

// ============================================================================
// NTT (数论变换) 乘法实现
// ============================================================================

namespace ntt {

constexpr uint32_t P1 = 998244353;
constexpr uint32_t P2 = 1004535809;
constexpr uint32_t P3 = 469762049;

template <uint32_t P>
struct NTTConfig {
    static constexpr uint32_t mod = P;
    static constexpr uint32_t g = 3;

    static uint32_t power(uint32_t a, uint32_t b) {
        uint32_t res = 1;
        a %= mod;
        while (b > 0) {
            if (b & 1) res = static_cast<uint32_t>((static_cast<uint64_t>(res) * a) % mod);
            a = static_cast<uint32_t>((static_cast<uint64_t>(a) * a) % mod);
            b >>= 1;
        }
        return res;
    }

    static uint32_t inv(uint32_t n) {
        return power(n, mod - 2);
    }

    static void transform(uint32_t* a, std::size_t n, bool invert) {
        for (std::size_t i = 1, j = 0; i < n; i++) {
            std::size_t bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) std::swap(a[i], a[j]);
        }

        for (std::size_t len = 2; len <= n; len <<= 1) {
            uint32_t wlen = power(g, (mod - 1) / static_cast<uint32_t>(len));
            if (invert) wlen = inv(wlen);
            for (std::size_t i = 0; i < n; i += len) {
                uint32_t w = 1;
                for (std::size_t j = 0; j < len / 2; j++) {
                    uint32_t u = a[i + j];
                    uint32_t v = static_cast<uint32_t>((static_cast<uint64_t>(a[i + j + len / 2]) * w) % mod);
                    a[i + j] = (u + v < mod) ? (u + v) : (u + v - mod);
                    a[i + j + len / 2] = (u >= v) ? (u - v) : (u + mod - v);
                    w = static_cast<uint32_t>((static_cast<uint64_t>(w) * wlen) % mod);
                }
            }
        }

        if (invert) {
            uint32_t n_inv = inv(static_cast<uint32_t>(n));
            for (std::size_t i = 0; i < n; ++i) a[i] = static_cast<uint32_t>((static_cast<uint64_t>(a[i]) * n_inv) % mod);
        }
    }
};

BigIntData crt(const uint32_t* r1, const uint32_t* r2, const uint32_t* r3, std::size_t n, uint32_t base) {
    BigIntData res;
    res.reserve(n + 2);

    const uint64_t M12 = static_cast<uint64_t>(P1) * P2;
    const uint32_t invP1P2 = NTTConfig<P2>::inv(P1);
    const uint32_t invM12P3 = NTTConfig<P3>::inv(static_cast<uint32_t>(M12 % P3));

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
    unsigned __int128 carry = 0;
    for (std::size_t i = 0; i < n; ++i) {
        uint32_t k1 = static_cast<uint32_t>((static_cast<uint64_t>(r2[i] + P2 - r1[i]) * invP1P2) % P2);
        uint64_t r12 = r1[i] + static_cast<uint64_t>(k1) * P1;

        uint32_t r12_mod_p3 = static_cast<uint32_t>(r12 % P3);
        uint32_t k2 = static_cast<uint32_t>((static_cast<uint64_t>(r3[i] + P3 - r12_mod_p3) * invM12P3) % P3);

        unsigned __int128 val = carry + r12 + (unsigned __int128)k2 * M12;
        res.push_back(static_cast<uint32_t>(val % base));
        carry = val / base;
    }
#pragma GCC diagnostic pop

    while (carry > 0) {
        res.push_back(static_cast<uint32_t>(carry % base));
        carry /= base;
    }

    return res;
}

} // namespace ntt

BigIntData multiply_bigint_ntt(const BigIntData& lhs, const BigIntData& rhs) {
    std::size_t n = 1;
    while (n < lhs.size() + rhs.size()) n <<= 1;

    static thread_local BigIntArena arena;
    const std::size_t required_arena_words = n * 6 + lhs.size() + rhs.size() + 32;
    if (arena.capacity() < required_arena_words) {
        arena = BigIntArena(required_arena_words);
    }
    arena.reset();

    uint32_t* fa1 = arena.allocate(n, 0);
    uint32_t* fa2 = arena.allocate(n, 0);
    uint32_t* fa3 = arena.allocate(n, 0);
    uint32_t* fb1 = arena.allocate(n, 0);
    uint32_t* fb2 = arena.allocate(n, 0);
    uint32_t* fb3 = arena.allocate(n, 0);

    // 复制并取模
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        fa1[i] = lhs[i] % ntt::P1;
        fa2[i] = lhs[i] % ntt::P2;
        fa3[i] = lhs[i] % ntt::P3;
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        fb1[i] = rhs[i] % ntt::P1;
        fb2[i] = rhs[i] % ntt::P2;
        fb3[i] = rhs[i] % ntt::P3;
    }

    ntt::NTTConfig<ntt::P1>::transform(fa1, n, false);
    ntt::NTTConfig<ntt::P1>::transform(fb1, n, false);
    ntt::NTTConfig<ntt::P2>::transform(fa2, n, false);
    ntt::NTTConfig<ntt::P2>::transform(fb2, n, false);
    ntt::NTTConfig<ntt::P3>::transform(fa3, n, false);
    ntt::NTTConfig<ntt::P3>::transform(fb3, n, false);

    for (std::size_t i = 0; i < n; ++i) {
        fa1[i] = static_cast<uint32_t>((static_cast<uint64_t>(fa1[i]) * fb1[i]) % ntt::P1);
        fa2[i] = static_cast<uint32_t>((static_cast<uint64_t>(fa2[i]) * fb2[i]) % ntt::P2);
        fa3[i] = static_cast<uint32_t>((static_cast<uint64_t>(fa3[i]) * fb3[i]) % ntt::P3);
    }

    ntt::NTTConfig<ntt::P1>::transform(fa1, n, true);
    ntt::NTTConfig<ntt::P2>::transform(fa2, n, true);
    ntt::NTTConfig<ntt::P3>::transform(fa3, n, true);

    BigIntData result = ntt::crt(fa1, fa2, fa3, n, 1000000000);
    while (result.size() > 1 && result.back() == 0) result.pop_back();
    return result;
}

// ============================================================================
// Karatsuba 乘法
// ============================================================================

BigIntData multiply_bigint_karatsuba(const BigIntData& lhs, const BigIntData& rhs) {
    std::size_t n = std::max(lhs.size(), rhs.size());

    if (n <= KARATSUBA_THRESHOLD) {
        return multiply_bigint_naive(lhs, rhs);
    }

    std::size_t m = n / 2;

    BigIntData a0 = bigint_low(lhs, m);
    BigIntData a1 = bigint_high(lhs, m);
    BigIntData b0 = bigint_low(rhs, m);
    BigIntData b1 = bigint_high(rhs, m);

    trim_bigint(a0);
    trim_bigint(a1);
    trim_bigint(b0);
    trim_bigint(b1);

    BigIntData z0 = multiply_bigint_karatsuba(a0, b0);
    BigIntData z2 = multiply_bigint_karatsuba(a1, b1);

    BigIntData a0_plus_a1 = add_bigint(a0, a1);
    BigIntData b0_plus_b1 = add_bigint(b0, b1);
    BigIntData z1_temp = multiply_bigint_karatsuba(a0_plus_a1, b0_plus_b1);
    BigIntData z1 = subtract_bigint(subtract_bigint(z1_temp, z0), z2);

    BigIntData result = add_bigint(
        add_bigint(shift_bigint(z2, 2 * m), shift_bigint(z1, m)),
        z0
    );

    trim_bigint(result);
    return result;
}

// ============================================================================
// Toom-Cook 3 乘法
// ============================================================================

BigIntData divide_bigint_by_2(BigIntData v) {
    if (v.empty() || (v.size() == 1 && v[0] == 0)) return {0};

    uint64_t rem = 0;
    for (int i = static_cast<int>(v.size()) - 1; i >= 0; --i) {
        uint64_t cur = rem * kBase + v[i];
        v[i] = static_cast<uint32_t>(cur / 2);
        rem = cur % 2;
    }
    while (v.size() > 1 && v.back() == 0) v.pop_back();
    return v;
}

BigIntData divide_bigint_by_3(BigIntData v) {
    if (v.empty() || (v.size() == 1 && v[0] == 0)) return {0};

    uint64_t rem = 0;
    for (int i = static_cast<int>(v.size()) - 1; i >= 0; --i) {
        uint64_t cur = rem * kBase + v[i];
        v[i] = static_cast<uint32_t>(cur / 3);
        rem = cur % 3;
    }
    while (v.size() > 1 && v.back() == 0) v.pop_back();
    return v;
}

struct SignedBigInt {
    BigIntData data;
    bool negative;

    SignedBigInt() : data{0}, negative(false) {}
    SignedBigInt(const BigIntData& d, bool neg = false) : data(d), negative(neg) {}
    SignedBigInt(BigIntData&& d, bool neg = false) : data(std::move(d)), negative(neg) {}

    bool is_zero() const {
        return data.empty() || (data.size() == 1 && data[0] == 0);
    }

    void normalize() {
        while (data.size() > 1 && data.back() == 0) data.pop_back();
        if (is_zero()) negative = false;
    }

    SignedBigInt operator-() const {
        SignedBigInt result = *this;
        if (!is_zero()) result.negative = !result.negative;
        return result;
    }

    static SignedBigInt add(const SignedBigInt& a, const SignedBigInt& b);
    static SignedBigInt sub(const SignedBigInt& a, const SignedBigInt& b);
};

SignedBigInt SignedBigInt::add(const SignedBigInt& a, const SignedBigInt& b) {
    if (a.negative == b.negative) {
        SignedBigInt result(add_bigint(a.data, b.data), a.negative);
        result.normalize();
        return result;
    }

    int cmp = compare_bigint(a.data, b.data);
    if (cmp == 0) {
        return SignedBigInt{0, false};
    } else if (cmp > 0) {
        SignedBigInt result(subtract_bigint(a.data, b.data), a.negative);
        result.normalize();
        return result;
    } else {
        SignedBigInt result(subtract_bigint(b.data, a.data), b.negative);
        result.normalize();
        return result;
    }
}

SignedBigInt SignedBigInt::sub(const SignedBigInt& a, const SignedBigInt& b) {
    return add(a, -b);
}

BigIntData multiply_bigint_toom3(const BigIntData& lhs, const BigIntData& rhs) {
    std::size_t n = std::max(lhs.size(), rhs.size());

    if (n <= 256) {
        return multiply_bigint_karatsuba(lhs, rhs);
    }

    std::size_t m = (n + 2) / 3;

    BigIntData x0 = bigint_low(lhs, m);
    BigIntData x1 = bigint_low(bigint_high(lhs, m), m);
    BigIntData x2 = bigint_high(lhs, 2 * m);

    BigIntData y0 = bigint_low(rhs, m);
    BigIntData y1 = bigint_low(bigint_high(rhs, m), m);
    BigIntData y2 = bigint_high(rhs, 2 * m);

    trim_bigint(x0); trim_bigint(x1); trim_bigint(x2);
    trim_bigint(y0); trim_bigint(y1); trim_bigint(y2);

    BigIntData x0_plus_x2 = add_bigint(x0, x2);
    BigIntData y0_plus_y2 = add_bigint(y0, y2);

    BigIntData m0 = multiply_bigint_toom3(x0, y0);
    BigIntData m_inf = multiply_bigint_toom3(x2, y2);

    BigIntData p1 = add_bigint(x0_plus_x2, x1);
    BigIntData q1 = add_bigint(y0_plus_y2, y1);
    BigIntData m1 = multiply_bigint_toom3(p1, q1);

    SignedBigInt s_p_neg1, s_q_neg1;

    int cmp_x1 = compare_bigint(x0_plus_x2, x1);
    if (cmp_x1 >= 0) {
        s_p_neg1 = SignedBigInt(subtract_bigint(x0_plus_x2, x1), false);
    } else {
        s_p_neg1 = SignedBigInt(subtract_bigint(x1, x0_plus_x2), true);
    }

    int cmp_y1 = compare_bigint(y0_plus_y2, y1);
    if (cmp_y1 >= 0) {
        s_q_neg1 = SignedBigInt(subtract_bigint(y0_plus_y2, y1), false);
    } else {
        s_q_neg1 = SignedBigInt(subtract_bigint(y1, y0_plus_y2), true);
    }

    bool m_neg1_neg = s_p_neg1.negative != s_q_neg1.negative;
    BigIntData m_neg1 = multiply_bigint_toom3(s_p_neg1.data, s_q_neg1.data);

    BigIntData x1_2 = multiply_bigint_by_uint32(x1, 2);
    BigIntData x2_4 = multiply_bigint_by_uint32(x2, 4);
    BigIntData y1_2 = multiply_bigint_by_uint32(y1, 2);
    BigIntData y2_4 = multiply_bigint_by_uint32(y2, 4);

    BigIntData p2 = add_bigint(add_bigint(x0, x1_2), x2_4);
    BigIntData q2 = add_bigint(add_bigint(y0, y1_2), y2_4);
    BigIntData m2 = multiply_bigint_toom3(p2, q2);

    SignedBigInt s_m0(m0, false);
    SignedBigInt s_m1(m1, false);
    SignedBigInt s_m_neg1(m_neg1, m_neg1_neg);
    SignedBigInt s_m2(m2, false);
    SignedBigInt s_m_inf(m_inf, false);

    SignedBigInt s_c0 = s_m0;
    SignedBigInt s_c4 = s_m_inf;

    // 3. r0 = (P(1) - P(-1)) / 2 = c1 + c3
    SignedBigInt s_r0 = SignedBigInt::sub(s_m1, s_m_neg1);
    s_r0.data = divide_bigint_by_2(std::move(s_r0.data));
    s_r0.normalize();

    // 4. r1 = (P(1) + P(-1)) / 2 = c0 + c2 + c4
    SignedBigInt s_r1 = SignedBigInt::add(s_m1, s_m_neg1);
    s_r1.data = divide_bigint_by_2(std::move(s_r1.data));
    s_r1.normalize();

    // 5. c2 = r1 - c0 - c4
    SignedBigInt s_c2 = SignedBigInt::sub(SignedBigInt::sub(s_r1, s_c0), s_c4);
    s_c2.normalize();

    // 6. r2 = (P(2) - P(-1)) / 3 = c1 + c2 + 3c3 + 5c4
    SignedBigInt s_r2 = SignedBigInt::sub(s_m2, s_m_neg1);
    s_r2.data = divide_bigint_by_3(std::move(s_r2.data));
    s_r2.normalize();

    // 7. 2c3 = r2 - c2 - 5c4 - r0
    SignedBigInt s_c4_5;
    s_c4_5.data = multiply_bigint_by_uint32(s_c4.data, 5);
    s_c4_5.negative = s_c4.negative;

    SignedBigInt s_2c3 = SignedBigInt::sub(SignedBigInt::sub(SignedBigInt::sub(s_r2, s_c2), s_c4_5), s_r0);
    s_2c3.normalize();
    
    // 8. c3 = (2c3) / 2
    SignedBigInt s_c3 = s_2c3;
    s_c3.data = divide_bigint_by_2(std::move(s_c3.data));
    s_c3.normalize();

    // 9. c1 = r0 - c3
    SignedBigInt s_c1 = SignedBigInt::sub(s_r0, s_c3);
    s_c1.normalize();

    if (s_c0.negative || s_c1.negative || s_c2.negative || s_c3.negative || s_c4.negative) {
        return multiply_bigint_karatsuba(lhs, rhs);
    }

    BigIntData result = s_c0.data;
    result = add_bigint(result, shift_bigint(s_c1.data, m));
    result = add_bigint(result, shift_bigint(s_c2.data, 2 * m));
    result = add_bigint(result, shift_bigint(s_c3.data, 3 * m));
    result = add_bigint(result, shift_bigint(s_c4.data, 4 * m));

    trim_bigint(result);
    return result;
}

// ============================================================================
// 主乘法函数（自动选择算法 - 优化版）
// ============================================================================

BigIntData multiply_bigint(const BigIntData& lhs, const BigIntData& rhs) {
    if ((lhs.size() == 1 && lhs[0] == 0) || (rhs.size() == 1 && rhs[0] == 0)) return {0};

    std::size_t max_size = std::max(lhs.size(), rhs.size());
    std::size_t min_size = std::min(lhs.size(), rhs.size());

    // 1. 小规模：朴素乘法
    //if (min_size <= 32 || max_size <= KARATSUBA_THRESHOLD) {
    //    return multiply_bigint_naive(lhs, rhs);
    //}

    // 2. 中等规模：零分配 Karatsuba（减少递归中的内存分配开销）
    if (max_size <= 512) {
        //return multiply_bigint_zero_allocation(lhs, rhs);
        return multiply_bigint_karatsuba(lhs, rhs);
    }

    // 3. 中大规模：Toom-Cook 3
    if (max_size <= 4096) {
        return multiply_bigint_toom3(lhs, rhs);
    }

    // 4. 特大规模：优化的 NTT
    return multiply_bigint_ntt_optimized(lhs, rhs);
}


} // namespace precise

#pragma GCC diagnostic pop

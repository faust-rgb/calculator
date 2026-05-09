// ============================================================================
// 精确小数实现 (优化版：基数 10^9)
// ============================================================================
//
// 本文件实现了 PreciseDecimal 类型的所有功能，包括：
// - 基于大整数的精确算术运算
// - NTT（数论变换）和 Karatsuba 高效乘法算法
// - Knuth 算法 D 用于高效大整数除法
// - 高精度数学函数（三角函数、指数对数等）
//
// 性能优化：
// - 小规模乘法使用朴素算法
// - 中等规模使用 Karatsuba 算法 (O(n^1.585))
// - 大规模使用 NTT 算法 (O(n log n))
// ============================================================================

// __int128 是 GCC 扩展，用于高精度整数运算
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

#include "precise_decimal.h"
#include "rational.h"
#include "types/stored_value.h"
#include "core/common/calculator_exceptions.h"
#include "math/mymath.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;  ///< 大整数基数（10^9）
constexpr int kBaseDigits = 9;          ///< 每个 chunk 的十进制位数
int g_default_scale = 40;               ///< 默认精度（小数点后位数）
constexpr long double kDisplayZeroEps = 1e-16;      ///< 显示零的阈值
constexpr long double kDisplayIntegerEps = 1e-9;    ///< 显示整数的阈值

// 前向声明辅助函数
PreciseDecimal one();
PreciseDecimal two();
PreciseDecimal half();
PreciseDecimal decimal_from_uint(uint32_t value);
PreciseDecimal scale_epsilon(int extra_digits = 4);
int series_iteration_limit(int minimum);

/**
 * @brief 规范化显示用的浮点数
 *
 * 对非常小的值返回零，对接近整数的值返回整数。
 *
 * @param value 输入值
 * @return 规范化后的值
 */
long double normalize_display_decimal(long double value) {
    if (mymath::abs(value) < kDisplayZeroEps) return 0.0L;
    if (mymath::abs(value) > kDisplayIntegerEps &&
        mymath::abs(value - mymath::round(value)) < kDisplayIntegerEps) {
        return mymath::round(value);
    }
    return value;
}

/**
 * @brief 格式化浮点数为字符串
 *
 * @param value 浮点数值
 * @param precision 显示精度（默认 12 位）
 * @return 格式化后的字符串
 */
std::string format_decimal(long double value, int precision = 12) {
    value = normalize_display_decimal(value);
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    return out.str();
}

/**
 * @brief 计算 a / b，返回商和余数 (均为 BigInt)
 *
 * 这是 Knuth 算法 D 的声明，在后面实现。
 */
void div_bigint(const std::vector<uint32_t>& num,
                const std::vector<uint32_t>& den,
                std::vector<uint32_t>* quotient,
                std::vector<uint32_t>* remainder);

// ============================================================================
// 大整数基础操作函数
// ============================================================================

/**
 * @brief 将整数字符串转换为基数为 10^9 的大整数数组
 *
 * 字符串表示从高位到低位，而结果数组是小端序（低位在前）。
 *
 * @param s 整数字符串，如 "123456789012345"
 * @return 大整数数组，小端序
 *
 * @example "1234567890123" -> [890123, 234567, 1]
 */
std::vector<uint32_t> string_to_bigint(const std::string& s) {
    if (s.empty() || s == "0") return {0};
    std::vector<uint32_t> res;
    for (int i = static_cast<int>(s.size()); i > 0; i -= kBaseDigits) {
        if (i < kBaseDigits) {
            res.push_back(static_cast<uint32_t>(std::stoul(s.substr(0, static_cast<std::size_t>(i)))));
        } else {
            res.push_back(static_cast<uint32_t>(std::stoul(s.substr(static_cast<std::size_t>(i - kBaseDigits), kBaseDigits))));
        }
    }
    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

/**
 * @brief 将基数为 10^9 的大整数数组转换为字符串
 *
 * 数组是小端序（低位在前），输出字符串是高位在前。
 *
 * @param v 大整数数组，小端序
 * @return 整数字符串
 *
 * @example [890123, 234567, 1] -> "1234567890123"
 */
std::string bigint_to_string(const std::vector<uint32_t>& v) {
    if (v.empty() || (v.size() == 1 && v[0] == 0)) return "0";
    std::ostringstream oss;
    oss << v.back();
    for (int i = static_cast<int>(v.size()) - 2; i >= 0; --i) {
        oss << std::setfill('0') << std::setw(kBaseDigits) << v[static_cast<std::size_t>(i)];
    }
    return oss.str();
}

/**
 * @brief 比较两个大整数的大小
 *
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return -1 表示 lhs < rhs，0 表示相等，1 表示 lhs > rhs
 */
int compare_bigint(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs) {
    if (lhs.size() != rhs.size()) return lhs.size() < rhs.size() ? -1 : 1;
    for (int i = static_cast<int>(lhs.size()) - 1; i >= 0; --i) {
        if (lhs[static_cast<std::size_t>(i)] != rhs[static_cast<std::size_t>(i)]) {
            return lhs[static_cast<std::size_t>(i)] < rhs[static_cast<std::size_t>(i)] ? -1 : 1;
        }
    }
    return 0;
}

/**
 * @brief 大整数加法
 *
 * 实现两个大整数的精确加法，自动处理进位。
 *
 * @param lhs 左操作数（小端序）
 * @param rhs 右操作数（小端序）
 * @return 和（小端序）
 */
std::vector<uint32_t> add_bigint(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs) {
    std::vector<uint32_t> res;
    res.reserve(std::max(lhs.size(), rhs.size()) + 1);
    uint32_t carry = 0;
    for (std::size_t i = 0; i < std::max(lhs.size(), rhs.size()) || carry; ++i) {
        uint64_t sum = carry + (i < lhs.size() ? lhs[i] : 0) + (i < rhs.size() ? rhs[i] : 0);
        res.push_back(static_cast<uint32_t>(sum % kBase));
        carry = static_cast<uint32_t>(sum / kBase);
    }
    return res;
}

/**
 * @brief 大整数减法（要求 lhs >= rhs）
 *
 * 实现两个大整数的减法，自动处理借位。
 * 注意：此函数假设 lhs >= rhs，否则结果不正确。
 *
 * @param lhs 被减数（小端序）
 * @param rhs 减数（小端序）
 * @return 差（小端序）
 */
std::vector<uint32_t> subtract_bigint(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs) {
    std::vector<uint32_t> res = lhs;
    int32_t borrow = 0;
    for (std::size_t i = 0; i < rhs.size() || borrow; ++i) {
        int64_t sub = static_cast<int64_t>(res[i]) - borrow - (i < rhs.size() ? rhs[i] : 0);
        if (sub < 0) {
            sub += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res[i] = static_cast<uint32_t>(sub);
    }
    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

// ============================================================================
// 乘法算法阈值配置
// ============================================================================

/**
 * @brief Karatsuba 乘法阈值
 *
 * 当数字长度小于此值时使用朴素乘法，否则使用 Karatsuba 算法。
 * 32 是经验值，在此规模下 Karatsuba 的额外开销不划算。
 */
constexpr std::size_t KARATSUBA_THRESHOLD = 32;

/**
 * @brief NTT 乘法阈值
 *
 * 当数字长度大于此值时使用 NTT（数论变换）算法。
 * 256 是经验值，在此规模以上 NTT 的 O(n log n) 优势显现。
 */
constexpr std::size_t NTT_THRESHOLD = 256;

// ============================================================================
// NTT (数论变换) 乘法实现
// ============================================================================
//
// NTT 是 FFT 在模数域的变体，用于高效计算多项式乘法。
// 这里使用三个不同的质数进行 NTT，然后用中国剩余定理（CRT）合并结果。
// ============================================================================

namespace ntt {

/**
 * @struct NTTConfig
 * @brief NTT 配置模板类
 *
 * 封装 NTT 变换所需的配置参数和辅助函数。
 * @tparam P 模数（必须是形如 k*2^n+1 的质数）
 */
template <uint32_t P>
struct NTTConfig {
    static constexpr uint32_t mod = P;  ///< 模数
    static constexpr uint32_t g = 3;    ///< 原根（对所选模数有效）

    /**
     * @brief 快速幂运算（模意义下）
     * @param a 底数
     * @param b 指数
     * @return a^b mod P
     */
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

    /**
     * @brief 求模逆元
     * @param n 要求逆元的数
     * @return n^(-1) mod P
     */
    static uint32_t inv(uint32_t n) {
        return power(n, mod - 2);
    }

    /**
     * @brief 执行 NTT 变换或逆变换
     *
     * 这是最核心的 NTT 实现，使用蝶形运算。
     *
     * @param a 待变换的数组（原地修改）
     * @param invert true 表示逆变换，false 表示正变换
     */
    static void transform(std::vector<uint32_t>& a, bool invert) {
        std::size_t n = a.size();
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
            for (uint32_t& x : a) x = static_cast<uint32_t>((static_cast<uint64_t>(x) * n_inv) % mod);
        }
    }
};

// ============================================================================
// NTT 使用的质数常量
// ============================================================================
//
// 三个适用于 NTT 的质数，都形如 k*2^23+1，支持最大 2^23 长度的变换。
// 这些质数足够大，使得 CRT 合并后的结果能正确表示乘积。

constexpr uint32_t P1 = 998244353;     ///< 第一个质数
constexpr uint32_t P2 = 1004535809;    ///< 第二个质数
constexpr uint32_t P3 = 469762049;     ///< 第三个质数

/**
 * @brief 中国剩余定理（CRT）合并三个 NTT 结果
 *
 * 将三个不同模数下的 NTT 结果合并为最终结果。
 * 使用 __int128_t 处理中间大数运算。
 *
 * @param r1 第一个质数模下的结果
 * @param r2 第二个质数模下的结果
 * @param r3 第三个质数模下的结果
 * @param base 输出基数（这里是 10^9）
 * @return 合并后的大整数数组
 */
std::vector<uint32_t> crt(const std::vector<uint32_t>& r1,
                          const std::vector<uint32_t>& r2,
                          const std::vector<uint32_t>& r3,
                          uint32_t base) {
    std::size_t n = r1.size();
    std::vector<uint32_t> res;
    res.reserve(n + 2);

    const uint64_t M12 = static_cast<uint64_t>(P1) * P2;
    const uint32_t invP1P2 = NTTConfig<P2>::inv(P1);
    const uint32_t invM12P3 = NTTConfig<P3>::inv(static_cast<uint32_t>(M12 % P3));

    unsigned __int128 carry = 0;
    for (std::size_t i = 0; i < n; ++i) {
        // r = r1 (mod P1)
        // r = r1 + P1 * k1 (mod P2) => k1 = (r2 - r1) * invP1 (mod P2)
        uint32_t k1 = static_cast<uint32_t>((static_cast<uint64_t>(r2[i] + P2 - r1[i]) * invP1P2) % P2);
        uint64_t r12 = r1[i] + static_cast<uint64_t>(k1) * P1;

        // r = r12 + M12 * k2 (mod P3) => k2 = (r3 - r12) * invM12 (mod P3)
        uint32_t r12_mod_p3 = static_cast<uint32_t>(r12 % P3);
        uint32_t k2 = static_cast<uint32_t>((static_cast<uint64_t>(r3[i] + P3 - r12_mod_p3) * invM12P3) % P3);
        
        unsigned __int128 val = carry + r12 + (unsigned __int128)k2 * M12;
        res.push_back(static_cast<uint32_t>(val % base));
        carry = val / base;
    }
    
    while (carry > 0) {
        res.push_back(static_cast<uint32_t>(carry % base));
        carry /= base;
    }
    
    return res;
}

} // namespace ntt

/**
 * @brief 大整数乘法 - NTT 算法
 *
 * 使用数论变换（NTT）实现高效的大整数乘法。
 * 时间复杂度 O(N log N)，适合超大规模计算。
 *
 * 算法流程：
 * 1. 将两个大整数视为多项式
 * 2. 对三个不同质数执行 NTT 变换
 * 3. 在变换域进行点乘
 * 4. 逆变换并用 CRT 合并结果
 *
 * @param lhs 左操作数（小端序）
 * @param rhs 右操作数（小端序）
 * @return 乘积（小端序）
 */
std::vector<uint32_t> multiply_bigint_ntt(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs) {
    std::size_t n = 1;
    while (n < lhs.size() + rhs.size()) n <<= 1;

    std::vector<uint32_t> fa1(n, 0), fb1(n, 0);
    std::vector<uint32_t> fa2(n, 0), fb2(n, 0);
    std::vector<uint32_t> fa3(n, 0), fb3(n, 0);

    for (std::size_t i = 0; i < lhs.size(); ++i) {
        fa1[i] = fa2[i] = fa3[i] = lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        fb1[i] = fb2[i] = fb3[i] = rhs[i];
    }

    ntt::NTTConfig<ntt::P1>::transform(fa1, false);
    ntt::NTTConfig<ntt::P1>::transform(fb1, false);
    ntt::NTTConfig<ntt::P2>::transform(fa2, false);
    ntt::NTTConfig<ntt::P2>::transform(fb2, false);
    ntt::NTTConfig<ntt::P3>::transform(fa3, false);
    ntt::NTTConfig<ntt::P3>::transform(fb3, false);

    for (std::size_t i = 0; i < n; ++i) {
        fa1[i] = static_cast<uint32_t>((static_cast<uint64_t>(fa1[i]) * fb1[i]) % ntt::P1);
        fa2[i] = static_cast<uint32_t>((static_cast<uint64_t>(fa2[i]) * fb2[i]) % ntt::P2);
        fa3[i] = static_cast<uint32_t>((static_cast<uint64_t>(fa3[i]) * fb3[i]) % ntt::P3);
    }

    ntt::NTTConfig<ntt::P1>::transform(fa1, true);
    ntt::NTTConfig<ntt::P2>::transform(fa2, true);
    ntt::NTTConfig<ntt::P3>::transform(fa3, true);

    std::vector<uint32_t> result = ntt::crt(fa1, fa2, fa3, 1000000000); // kBase = 10^9
    while (result.size() > 1 && result.back() == 0) result.pop_back();
    return result;
}

/**
 * @brief 大整数乘法 - Karatsuba 算法
 *
 * 使用分治策略的 Karatsuba 乘法算法。
 * 时间复杂度 O(n^1.585)，比朴素 O(n^2) 更快。
 *
 * 算法原理：
 * 将 x * y 分解为：
 *   (a1 * B^m + a0) * (b1 * B^m + b0)
 * = a1*b1 * B^(2m) + [(a0+a1)(b0+b1) - a0*b0 - a1*b1] * B^m + a0*b0
 * 只需要 3 次递归乘法而非 4 次。
 *
 * @param lhs 左操作数（小端序）
 * @param rhs 右操作数（小端序）
 * @return 乘积（小端序）
 */
std::vector<uint32_t> multiply_bigint_karatsuba(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs);

/**
 * @brief 朴素大整数乘法（用于小规模数据）
 *
 * 直接按照手工乘法的方式计算，时间复杂度 O(n^2)。
 * 虽然效率低，但在小规模数据下开销最小。
 *
 * @param lhs 左操作数（小端序）
 * @param rhs 右操作数（小端序）
 * @return 乘积（小端序）
 */
std::vector<uint32_t> multiply_bigint_naive(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs) {
    if ((lhs.size() == 1 && lhs[0] == 0) || (rhs.size() == 1 && rhs[0] == 0)) return {0};
    std::vector<uint32_t> res(lhs.size() + rhs.size(), 0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        uint32_t carry = 0;
        for (std::size_t j = 0; j < rhs.size() || carry; ++j) {
            uint64_t cur = res[i + j] +
                           static_cast<uint64_t>(lhs[i]) * (j < rhs.size() ? rhs[j] : 0) + carry;
            res[i + j] = static_cast<uint32_t>(cur % kBase);
            carry = static_cast<uint32_t>(cur / kBase);
        }
    }
    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

/**
 * @brief 大整数左移位（相当于乘以 base^n）
 *
 * 在大整数数组前面添加 n 个零，相当于数值乘以 base^n。
 *
 * @param v 大整数数组
 * @param n 移位量
 * @return 移位后的数组
 */
std::vector<uint32_t> shift_bigint(const std::vector<uint32_t>& v, std::size_t n) {
    if (v.empty() || (v.size() == 1 && v[0] == 0)) return {0};
    if (n == 0) return v;
    std::vector<uint32_t> res(n, 0);
    res.insert(res.end(), v.begin(), v.end());
    return res;
}

/**
 * @brief 大整数截取低位 n 个 chunk
 *
 * @param v 大整数数组
 * @param n 要截取的 chunk 数量
 * @return 截取后的数组
 */
std::vector<uint32_t> bigint_low(const std::vector<uint32_t>& v, std::size_t n) {
    if (n >= v.size()) return v;
    return std::vector<uint32_t>(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(n));
}

/**
 * @brief 大整数截取高位（从第 n 个 chunk 开始）
 *
 * @param v 大整数数组
 * @param n 起始 chunk 索引
 * @return 截取后的数组
 */
std::vector<uint32_t> bigint_high(const std::vector<uint32_t>& v, std::size_t n) {
    if (n >= v.size()) return {0};
    return std::vector<uint32_t>(v.begin() + static_cast<std::ptrdiff_t>(n), v.end());
}

/**
 * @brief 移除大整数前导零
 *
 * 原地修改，移除数组末尾的所有零（保留至少一个元素）。
 *
 * @param v 大整数数组（原地修改）
 */
void trim_bigint(std::vector<uint32_t>& v) {
    while (v.size() > 1 && v.back() == 0) v.pop_back();
}

/**
 * @brief 大整数乘法 - Karatsuba 算法实现
 *
 * 递归实现 Karatsuba 乘法，小规模时退化为朴素乘法。
 */
std::vector<uint32_t> multiply_bigint_karatsuba(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs) {
    std::size_t n = std::max(lhs.size(), rhs.size());

    // 基准情况：小规模使用朴素乘法
    if (n <= KARATSUBA_THRESHOLD) {
        return multiply_bigint_naive(lhs, rhs);
    }

    std::size_t m = n / 2;

    // 分解: lhs = a1 * B^m + a0, rhs = b1 * B^m + b0
    std::vector<uint32_t> a0 = bigint_low(lhs, m);
    std::vector<uint32_t> a1 = bigint_high(lhs, m);
    std::vector<uint32_t> b0 = bigint_low(rhs, m);
    std::vector<uint32_t> b1 = bigint_high(rhs, m);

    trim_bigint(a0);
    trim_bigint(a1);
    trim_bigint(b0);
    trim_bigint(b1);

    // 递归计算三个乘积
    std::vector<uint32_t> z0 = multiply_bigint_karatsuba(a0, b0);  // a0 * b0
    std::vector<uint32_t> z2 = multiply_bigint_karatsuba(a1, b1);  // a1 * b1

    // z1 = (a0 + a1) * (b0 + b1) - z0 - z2
    std::vector<uint32_t> a0_plus_a1 = add_bigint(a0, a1);
    std::vector<uint32_t> b0_plus_b1 = add_bigint(b0, b1);
    std::vector<uint32_t> z1_temp = multiply_bigint_karatsuba(a0_plus_a1, b0_plus_b1);
    std::vector<uint32_t> z1 = subtract_bigint(subtract_bigint(z1_temp, z0), z2);

    // 结果 = z2 * B^(2m) + z1 * B^m + z0
    std::vector<uint32_t> result = add_bigint(
        add_bigint(shift_bigint(z2, 2 * m), shift_bigint(z1, m)),
        z0
    );

    trim_bigint(result);
    return result;
}

/**
 * @brief 大整数乘法（自动选择最优算法）
 *
 * 根据输入规模自动选择最适合的乘法算法：
 * - 小规模（<= 32 chunks）：朴素乘法
 * - 中等规模（<= 256 chunks）：Karatsuba 算法
 * - 大规模（> 256 chunks）：NTT 算法
 *
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 乘积
 */
std::vector<uint32_t> multiply_bigint(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs) {
    // 零值快速返回
    if ((lhs.size() == 1 && lhs[0] == 0) || (rhs.size() == 1 && rhs[0] == 0)) return {0};

    std::size_t max_size = std::max(lhs.size(), rhs.size());
    
    // 小规模直接使用朴素乘法
    if (max_size <= KARATSUBA_THRESHOLD) {
        return multiply_bigint_naive(lhs, rhs);
    }

    // 中等规模使用 Karatsuba
    if (max_size <= NTT_THRESHOLD) {
        return multiply_bigint_karatsuba(lhs, rhs);
    }

    // 大规模使用 NTT
    return multiply_bigint_ntt(lhs, rhs);
}

/**
 * @brief 大整数乘以一个小整数
 *
 * 当乘数是单个 uint32_t 时使用此函数，效率更高。
 *
 * @param v 大整数数组
 * @param n 小整数乘数
 * @return 乘积
 */
std::vector<uint32_t> multiply_bigint_by_uint32(const std::vector<uint32_t>& v, uint32_t n) {
    if (n == 0 || (v.size() == 1 && v[0] == 0)) return {0};
    if (n == 1) return v;
    std::vector<uint32_t> res;
    res.reserve(v.size() + 1);
    uint32_t carry = 0;
    for (uint32_t x : v) {
        uint64_t cur = static_cast<uint64_t>(x) * n + carry;
        res.push_back(static_cast<uint32_t>(cur % kBase));
        carry = static_cast<uint32_t>(cur / kBase);
    }
    if (carry) res.push_back(carry);
    return res;
}

/**
 * @brief Knuth 算法 D：高效大整数除法
 *
 * 实现经典的 Knuth 除法算法，时间复杂度 O(n*m)。
 * 相比逐位试商的 O(n*m*10) 有显著提升。
 *
 * 算法步骤：
 * 1. D1: 归一化 - 使除数最高位 >= base/2
 * 2. D2-D7: 主循环，计算每一位商
 * 3. D8: 反归一化得到最终余数
 *
 * @param num 被除数
 * @param den 除数
 * @param quotient 输出：商
 * @param remainder 输出：余数
 * @throws std::runtime_error 如果除数为零
 */
void div_bigint(const std::vector<uint32_t>& num,
                const std::vector<uint32_t>& den,
                std::vector<uint32_t>* quotient,
                std::vector<uint32_t>* remainder) {
    if (den.empty() || (den.size() == 1 && den[0] == 0)) {
        throw std::runtime_error("division by zero in div_bigint");
    }

    // 处理被除数为零
    if (num.empty() || (num.size() == 1 && num[0] == 0)) {
        *quotient = {0};
        *remainder = {0};
        return;
    }

    // 确定实际长度（去除前导零）
    std::size_t n = num.size();
    while (n > 1 && num[n - 1] == 0) --n;
    std::size_t m = den.size();
    while (m > 1 && den[m - 1] == 0) --m;

    // 除数大于被除数
    if (m > n || (m == n && compare_bigint(
            std::vector<uint32_t>(num.begin(), num.begin() + static_cast<std::ptrdiff_t>(n)),
            std::vector<uint32_t>(den.begin(), den.begin() + static_cast<std::ptrdiff_t>(m))) < 0)) {
        *quotient = {0};
        *remainder = std::vector<uint32_t>(num.begin(), num.begin() + static_cast<std::ptrdiff_t>(n));
        while (remainder->size() > 1 && remainder->back() == 0) remainder->pop_back();
        return;
    }

    // 特殊情况：除数为单 chunk
    if (m == 1) {
        std::vector<uint32_t> q(n, 0);
        uint64_t rem = 0;
        uint32_t d = den[0];
        for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
            uint64_t cur = rem * kBase + num[i];
            q[i] = static_cast<uint32_t>(cur / d);
            rem = cur % d;
        }
        while (q.size() > 1 && q.back() == 0) q.pop_back();
        *quotient = std::move(q);
        *remainder = {static_cast<uint32_t>(rem)};
        return;
    }

    // Knuth 算法 D
    // D1: 归一化 - 找到 d 使得 den[m-1] * d >= kBase/2
    uint32_t d = 1;
    uint64_t top_den = den[m - 1];
    while (top_den < kBase / 2) {
        top_den *= 2;
        d *= 2;
    }

    // 归一化被除数和除数
    std::vector<uint32_t> u(n + 1, 0);  // 被除数，多一位防止溢出
    std::vector<uint32_t> v(m, 0);       // 除数

    uint64_t carry = 0;
    for (std::size_t i = 0; i < n; ++i) {
        uint64_t cur = static_cast<uint64_t>(num[i]) * d + carry;
        u[i] = static_cast<uint32_t>(cur % kBase);
        carry = cur / kBase;
    }
    u[n] = static_cast<uint32_t>(carry);

    carry = 0;
    for (std::size_t i = 0; i < m; ++i) {
        uint64_t cur = static_cast<uint64_t>(den[i]) * d + carry;
        v[i] = static_cast<uint32_t>(cur % kBase);
        carry = cur / kBase;
    }

    // D2-D7: 主循环
    std::vector<uint32_t> q(n - m + 1, 0);

    for (int j = static_cast<int>(n - m); j >= 0; --j) {
        // D3: 计算试商 q_hat
        uint64_t numerator_val = static_cast<uint64_t>(u[j + m]) * kBase + u[j + m - 1];
        uint64_t q_hat = numerator_val / v[m - 1];
        uint64_t r_hat = numerator_val % v[m - 1];

        // 修正 q_hat
        while (q_hat >= kBase ||
               (q_hat * v[m - 2] > r_hat * kBase + (j + static_cast<int>(m) - 2 >= 0 ? u[j + m - 2] : 0))) {
            --q_hat;
            r_hat += v[m - 1];
            if (r_hat >= kBase) break;
        }

        // D4: 乘法并减法
        int64_t borrow = 0;
        for (std::size_t i = 0; i < m; ++i) {
            uint64_t product = q_hat * v[i];
            int64_t diff = static_cast<int64_t>(u[j + i]) - static_cast<int64_t>(product % kBase) - borrow;
            borrow = static_cast<int64_t>(product / kBase);
            if (diff < 0) {
                diff += kBase;
                ++borrow;
            }
            u[j + i] = static_cast<uint32_t>(diff);
        }
        int64_t diff = static_cast<int64_t>(u[j + m]) - borrow;
        u[j + m] = static_cast<uint32_t>(diff);

        // D5: 测试余数
        q[j] = static_cast<uint32_t>(q_hat);
        if (diff < 0) {
            // D6: 加回
            --q[j];
            uint64_t c = 0;
            for (std::size_t i = 0; i < m; ++i) {
                uint64_t sum = static_cast<uint64_t>(u[j + i]) + v[i] + c;
                u[j + i] = static_cast<uint32_t>(sum % kBase);
                c = sum / kBase;
            }
            u[j + m] += static_cast<uint32_t>(c);
        }
    }

    // D8: 反归一化得到余数
    while (q.size() > 1 && q.back() == 0) q.pop_back();
    *quotient = std::move(q);

    // 计算余数
    std::vector<uint32_t> rem(m, 0);
    uint64_t r = 0;
    for (int i = static_cast<int>(m) - 1; i >= 0; --i) {
        uint64_t cur = r * kBase + u[i];
        rem[i] = static_cast<uint32_t>(cur / d);
        r = cur % d;
    }
    while (rem.size() > 1 && rem.back() == 0) rem.pop_back();
    *remainder = std::move(rem);
}

/**
 * @brief 大整数乘以 10^n
 *
 * 通过移动小数点位置来实现，分两步：
 * 1. 移动完整的 chunk（乘以 base^k）
 * 2. 移动剩余的位数（乘以 10^r）
 *
 * @param v 大整数数组
 * @param n 10 的幂次
 * @return 乘积
 */
std::vector<uint32_t> multiply_bigint_by_power_of_10(std::vector<uint32_t> v, int n) {
    if (n <= 0) return v;
    if (v.size() == 1 && v[0] == 0) return v;

    int chunk_shift = n / kBaseDigits;
    int digit_shift = n % kBaseDigits;

    if (digit_shift > 0) {
        uint32_t multiplier = 1;
        for (int i = 0; i < digit_shift; ++i) multiplier *= 10;
        v = multiply_bigint_by_uint32(v, multiplier);
    }

    if (chunk_shift > 0) {
        v.insert(v.begin(), static_cast<std::size_t>(chunk_shift), 0);
    }
    return v;
}

/**
 * @brief 对齐两个 PreciseDecimal 的 scale
 *
 * 通过扩展较小 scale 的数字来对齐两个数的小数点位置。
 * 例如：1.23 (scale=2) 和 0.1 (scale=1) 变为 1.23 和 0.10
 *
 * @param lhs 第一个数（可能被修改）
 * @param rhs 第二个数（可能被修改）
 */
void align_precise_scales(PreciseDecimal* lhs, PreciseDecimal* rhs) {
    if (lhs->scale == rhs->scale) return;
    if (lhs->scale < rhs->scale) {
        lhs->data = multiply_bigint_by_power_of_10(lhs->data, rhs->scale - lhs->scale);
        lhs->scale = rhs->scale;
    } else {
        rhs->data = multiply_bigint_by_power_of_10(rhs->data, lhs->scale - rhs->scale);
        rhs->scale = lhs->scale;
    }
}

// ============================================================================
// 高精度数学常量和辅助函数
// ============================================================================

/// 返回高精度整数 1
PreciseDecimal one() { return PreciseDecimal(1LL); }

/// 返回高精度整数 2
PreciseDecimal two() { return PreciseDecimal(2LL); }

/// 返回高精度小数 0.5
PreciseDecimal half() { return PreciseDecimal("0.5"); }

/// 从 uint32_t 构造高精度小数
PreciseDecimal decimal_from_uint(uint32_t value) {
    return PreciseDecimal(static_cast<long long>(value));
}

/**
 * @brief 生成精度相关的 epsilon 值
 *
 * 用于级数展开的收敛判断。
 *
 * @param extra_digits 额外精度位数
 * @return 10^(-(default_scale + extra_digits))
 */
PreciseDecimal scale_epsilon(int extra_digits) {
    return PreciseDecimal("1e-" + std::to_string(PrecisionContext::get_default_scale() + extra_digits));
}

/**
 * @brief 计算级数迭代上限
 *
 * @param minimum 最小迭代次数
 * @return 迭代次数上限
 */
int series_iteration_limit(int minimum) {
    return std::max(minimum, PrecisionContext::get_default_scale() * 2 + 20);
}

/**
 * @brief 将值规约到正模数范围内
 *
 * 计算 value mod modulus，保证结果在 [0, modulus) 范围内。
 *
 * @param value 输入值
 * @param modulus 模数（必须为正）
 * @return 规约后的值
 */
PreciseDecimal reduce_mod_positive(const PreciseDecimal& value, const PreciseDecimal& modulus) {
    PreciseDecimal q = precise::floor(value / modulus);
    PreciseDecimal r = value - q * modulus;
    if (r < PreciseDecimal(0LL)) r += modulus;
    return r;
}

/**
 * @brief 规约三角函数参数到 [-π, π] 范围
 *
 * 通过模 2π 运算将参数规约到主值区间，提高三角函数计算的效率和精度。
 *
 * @param x 输入弧度值
 * @return 规约后的弧度值，在 [-π, π] 范围内
 */
PreciseDecimal reduce_trig_argument(const PreciseDecimal& x) {
    // 使用预计算的高精度常量，避免动态计算引入误差
    const PreciseDecimal p = precise::pi();
    const PreciseDecimal two_pi_val = precise::two_pi();
    PreciseDecimal r = reduce_mod_positive(x, two_pi_val);
    if (r > p) r -= two_pi_val;
    return r;
}

} // namespace

int PrecisionContext::get_default_scale() { return g_default_scale; }
void PrecisionContext::set_default_scale(int s) { g_default_scale = s; }

// ============================================================================
// PreciseDecimal 方法实现
// ============================================================================

PreciseDecimal::PreciseDecimal(long long value) {
    negative = value < 0;
    data = string_to_bigint(std::to_string(negative ? -value : value));
    scale = 0;
}

PreciseDecimal::PreciseDecimal(long double value) {
    *this = from_decimal_literal(format_decimal(value, 15));
}

PreciseDecimal::PreciseDecimal(const std::string& token) {
    *this = from_decimal_literal(token);
}

void PreciseDecimal::normalize() {
    if (data.empty()) {
        data = {0};
        scale = 0;
        negative = false;
        return;
    }

    // 去除末尾的零（缩减 scale）
    while (scale > 0 && !is_zero()) {
        // 检查最后一个 chunk 是否被 10 整除
        if (data[0] % 10 == 0) {
            // 对所有 chunks 进行除以 10 的操作
            uint32_t rem = 0;
            for (int i = static_cast<int>(data.size()) - 1; i >= 0; --i) {
                uint64_t current = data[i] + static_cast<uint64_t>(rem) * kBase;
                data[i] = static_cast<uint32_t>(current / 10);
                rem = static_cast<uint32_t>(current % 10);
            }
            scale--;
            // 如果最高位变 0，移除它
            while (data.size() > 1 && data.back() == 0) data.pop_back();
        } else {
            break;
        }
    }

    if (is_zero()) {
        negative = false;
        scale = 0;
    }
}

bool PreciseDecimal::is_zero() const {
    return data.empty() || (data.size() == 1 && data[0] == 0);
}

PreciseDecimal PreciseDecimal::operator-() const {
    PreciseDecimal res = *this;
    if (!res.is_zero()) res.negative = !res.negative;
    return res;
}

PreciseDecimal& PreciseDecimal::operator+=(const PreciseDecimal& rhs) {
    *this = add_precise_decimal(*this, rhs);
    return *this;
}

PreciseDecimal& PreciseDecimal::operator-=(const PreciseDecimal& rhs) {
    *this = subtract_precise_decimal(*this, rhs);
    return *this;
}

PreciseDecimal& PreciseDecimal::operator*=(const PreciseDecimal& rhs) {
    *this = multiply_precise_decimal(*this, rhs);
    return *this;
}

PreciseDecimal& PreciseDecimal::operator/=(const PreciseDecimal& rhs) {
    *this = divide_precise_decimal(*this, rhs);
    return *this;
}

bool PreciseDecimal::operator==(const PreciseDecimal& rhs) const {
    return compare_precise_decimal(*this, rhs) == 0;
}

bool PreciseDecimal::operator<(const PreciseDecimal& rhs) const {
    return compare_precise_decimal(*this, rhs) < 0;
}

std::string PreciseDecimal::to_string() const {
    if (is_zero()) return "0";

    std::string s = bigint_to_string(data);
    std::string res;
    if (negative) res.push_back('-');

    if (scale <= 0) {
        res += s;
        if (scale < 0) res.append(static_cast<std::size_t>(-scale), '0');
        // 限制最多输出50位
        if (res.size() > 50) {
            bool is_neg = !res.empty() && res[0] == '-';
            if (is_neg) {
                res = res.substr(0, 51);  // 负号 + 50位数字
            } else {
                res = res.substr(0, 50);
            }
        }
        return res;
    }

    if (static_cast<int>(s.size()) <= scale) {
        res += "0.";
        res.append(static_cast<std::size_t>(scale - static_cast<int>(s.size())), '0');
        res += s;
    } else {
        std::size_t dot_pos = s.size() - static_cast<std::size_t>(scale);
        res += s.substr(0, dot_pos);
        res.push_back('.');
        res += s.substr(dot_pos);
    }

    // 限制最多输出50位（不含小数点和负号）
    std::size_t sign_len = (negative ? 1 : 0);
    std::size_t dot_pos = res.find('.');
    if (dot_pos != std::string::npos) {
        // 小数形式：统计整数部分和小数部分的数字总数
        std::size_t int_digits = dot_pos - sign_len;
        std::size_t frac_digits = res.size() - dot_pos - 1;
        std::size_t total_digits = int_digits + frac_digits;
        if (total_digits > 50) {
            // 保留前50位数字
            std::size_t keep_frac = 50 - int_digits;
            if (keep_frac > 0) {
                res = res.substr(0, dot_pos + 1 + keep_frac);
            } else {
                // 整数部分已经超过50位，截断整数部分
                res = res.substr(0, sign_len + 50);
            }
        }
    } else {
        // 整数形式
        if (res.size() - sign_len > 50) {
            res = res.substr(0, sign_len + 50);
        }
    }

    return res;
}

long double PreciseDecimal::to_double() const {
    if (is_zero()) return 0.0L;

    // 找出第一个非零 chunk 的索引
    int start_chunk = static_cast<int>(data.size()) - 1;
    while (start_chunk >= 0 && data[start_chunk] == 0) start_chunk--;

    if (start_chunk < 0) return 0.0L;

    // 取最高位的 4 个 chunk (约 36 位十进制)，足以填满 long double 的 64 位尾数
    long double res = 0;
    int end_chunk = std::max(0, start_chunk - 3);
    for (int i = start_chunk; i >= end_chunk; --i) {
        res = res * static_cast<long double>(kBase) + static_cast<long double>(data[i]);
    }

    // 计算指数偏移: (end_chunk * 9) - scale
    long double exponent = static_cast<long double>(end_chunk) * 9.0L - static_cast<long double>(scale);

    if (exponent != 0.0L) {
        res *= mymath::pow(10.0L, exponent);
    }

    if (!mymath::isfinite(static_cast<long double>(res))) {
        return negative ? -mymath::infinity() : mymath::infinity();
    }

    return static_cast<long double>(negative ? -res : res);
}

PreciseDecimal PreciseDecimal::from_digits(std::string raw_digits, int raw_scale, bool is_negative) {
    PreciseDecimal res;
    res.data = string_to_bigint(raw_digits);
    res.scale = raw_scale;
    res.negative = is_negative;
    res.normalize();
    return res;
}

PreciseDecimal PreciseDecimal::from_integer_string(const std::string& integer_text, bool is_negative) {
    return from_digits(integer_text, 0, is_negative);
}

PreciseDecimal PreciseDecimal::from_decimal_literal(const std::string& token) {
    std::string text;
    text.reserve(token.size());
    for (char ch : token) {
        if (!std::isspace(static_cast<unsigned char>(ch))) text.push_back(ch);
    }
    if (text.empty()) return {};

    bool is_negative = false;
    if (text.front() == '+' || text.front() == '-') {
        is_negative = text.front() == '-';
        text.erase(text.begin());
    }
    if (text.empty()) return {};

    std::string significand = text;
    int exponent_adjust = 0;
    const std::size_t exponent_pos = text.find_first_of("eE");
    if (exponent_pos != std::string::npos) {
        significand = text.substr(0, exponent_pos);
        exponent_adjust = std::stoi(text.substr(exponent_pos + 1));
    }

    const std::size_t dot_pos = significand.find('.');
    if (dot_pos == std::string::npos) {
        return from_digits(significand, -exponent_adjust, is_negative);
    }

    std::string digits_only = significand.substr(0, dot_pos);
    digits_only += significand.substr(dot_pos + 1);
    return from_digits(
        digits_only,
        static_cast<int>(significand.size() - dot_pos - 1) - exponent_adjust,
        is_negative);
}

// ============================================================================
// PreciseDecimal 算术运算
// ============================================================================

PreciseDecimal operator+(PreciseDecimal lhs, const PreciseDecimal& rhs) { lhs += rhs; return lhs; }
PreciseDecimal operator+(PreciseDecimal lhs, long double rhs) { lhs += PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator+(long double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) + rhs; }
PreciseDecimal operator+(PreciseDecimal lhs, double rhs) { lhs += PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator+(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) + rhs; }
PreciseDecimal operator+(PreciseDecimal lhs, int rhs) { lhs += PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator+(int lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) + rhs; }

PreciseDecimal operator-(PreciseDecimal lhs, const PreciseDecimal& rhs) { lhs -= rhs; return lhs; }
PreciseDecimal operator-(PreciseDecimal lhs, long double rhs) { lhs -= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator-(long double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) - rhs; }
PreciseDecimal operator-(PreciseDecimal lhs, double rhs) { lhs -= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator-(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) - rhs; }
PreciseDecimal operator-(PreciseDecimal lhs, int rhs) { lhs -= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator-(int lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) - rhs; }

PreciseDecimal operator*(PreciseDecimal lhs, const PreciseDecimal& rhs) { lhs *= rhs; return lhs; }
PreciseDecimal operator*(PreciseDecimal lhs, long double rhs) { lhs *= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator*(long double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) * rhs; }
PreciseDecimal operator*(PreciseDecimal lhs, double rhs) { lhs *= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator*(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) * rhs; }
PreciseDecimal operator*(PreciseDecimal lhs, int rhs) { lhs *= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator*(int lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) * rhs; }

PreciseDecimal operator/(PreciseDecimal lhs, const PreciseDecimal& rhs) { lhs /= rhs; return lhs; }
PreciseDecimal operator/(PreciseDecimal lhs, long double rhs) { lhs /= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator/(long double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) / rhs; }
PreciseDecimal operator/(PreciseDecimal lhs, double rhs) { lhs /= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator/(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) / rhs; }
PreciseDecimal operator/(PreciseDecimal lhs, int rhs) { lhs /= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator/(int lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) / rhs; }

bool operator==(const PreciseDecimal& lhs, long double rhs) { return lhs == PreciseDecimal(rhs); }
bool operator==(long double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) == rhs; }
bool operator!=(const PreciseDecimal& lhs, long double rhs) { return !(lhs == rhs); }
bool operator!=(long double lhs, const PreciseDecimal& rhs) { return !(lhs == rhs); }
bool operator<(const PreciseDecimal& lhs, long double rhs) { return lhs < PreciseDecimal(rhs); }
bool operator<(long double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) < rhs; }
bool operator>(const PreciseDecimal& lhs, long double rhs) { return lhs > PreciseDecimal(rhs); }
bool operator>(long double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) > rhs; }
bool operator<=(const PreciseDecimal& lhs, long double rhs) { return lhs <= PreciseDecimal(rhs); }
bool operator<=(long double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) <= rhs; }
bool operator>=(const PreciseDecimal& lhs, long double rhs) { return lhs >= PreciseDecimal(rhs); }
bool operator>=(long double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) >= rhs; }

bool operator==(const PreciseDecimal& lhs, double rhs) { return lhs == PreciseDecimal(rhs); }
bool operator==(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) == rhs; }
bool operator!=(const PreciseDecimal& lhs, double rhs) { return !(lhs == rhs); }
bool operator!=(double lhs, const PreciseDecimal& rhs) { return !(lhs == rhs); }
bool operator<(const PreciseDecimal& lhs, double rhs) { return lhs < PreciseDecimal(rhs); }
bool operator<(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) < rhs; }
bool operator>(const PreciseDecimal& lhs, double rhs) { return lhs > PreciseDecimal(rhs); }
bool operator>(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) > rhs; }
bool operator<=(const PreciseDecimal& lhs, double rhs) { return lhs <= PreciseDecimal(rhs); }
bool operator<=(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) <= rhs; }
bool operator>=(const PreciseDecimal& lhs, double rhs) { return lhs >= PreciseDecimal(rhs); }
bool operator>=(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) >= rhs; }

bool operator==(const PreciseDecimal& lhs, int rhs) { return lhs == PreciseDecimal(rhs); }
bool operator==(int lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) == rhs; }
bool operator!=(const PreciseDecimal& lhs, int rhs) { return !(lhs == rhs); }
bool operator!=(int lhs, const PreciseDecimal& rhs) { return !(lhs == rhs); }
bool operator<(const PreciseDecimal& lhs, int rhs) { return lhs < PreciseDecimal(rhs); }
bool operator<(int lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) < rhs; }
bool operator>(const PreciseDecimal& lhs, int rhs) { return lhs > PreciseDecimal(rhs); }
bool operator>(int lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) > rhs; }
bool operator<=(const PreciseDecimal& lhs, int rhs) { return lhs <= PreciseDecimal(rhs); }
bool operator<=(int lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) <= rhs; }
bool operator>=(const PreciseDecimal& lhs, int rhs) { return lhs >= PreciseDecimal(rhs); }
bool operator>=(int lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) >= rhs; }

std::ostream& operator<<(std::ostream& os, const PreciseDecimal& val) {
    os << val.to_string();
    return os;
}

namespace precise {

PreciseDecimal abs(const PreciseDecimal& val) {
    PreciseDecimal res = val;
    res.negative = false;
    return res;
}

PreciseDecimal sqrt(const PreciseDecimal& val) {
    if (val.is_zero()) return PreciseDecimal(0LL);
    if (val.negative) throw PreciseDecimalUnsupported("sqrt of negative number");

    // Newton 迭代: x_{n+1} = 0.5 * (x_n + val / x_n)
    // 初始猜测使用 long double 版 sqrt
    PreciseDecimal x(mymath::sqrt(val.to_double()));
    const PreciseDecimal one_half = half();

    const int iterations = std::max(12, PrecisionContext::get_default_scale() / 8 + 8);
    for (int i = 0; i < iterations; ++i) {
        x = one_half * (x + val / x);
    }
    return x;
}

PreciseDecimal pow(const PreciseDecimal& base, long long exp) {
    return pow_precise_decimal(base, exp);
}

PreciseDecimal pow(const PreciseDecimal& base, const PreciseDecimal& exp) {
    if (exp.is_zero()) return one();
    if (base.is_zero()) {
        if (exp.negative) throw std::domain_error("0^negative is undefined");
        return PreciseDecimal(0LL);
    }
    if (base == one()) return one();
    if (base.negative) {
        // 对于负数底数，只有整数指数才有意义
        if (exp.scale == 0) {
            return pow(base, static_cast<long long>(exp.to_double()));
        }
        throw std::domain_error("negative base with non-integer exponent");
    }
    // 使用 exp(ln(base) * exp) 计算
    return precise::exp(precise::ln(base) * exp);
}

PreciseDecimal floor(const PreciseDecimal& val) {
    if (val.scale <= 0) return val;
    std::string s = bigint_to_string(val.data);
    if (static_cast<int>(s.size()) <= val.scale) {
        return val.negative ? PreciseDecimal(-1LL) : PreciseDecimal(0LL);
    }
    std::string integer_part = s.substr(0, s.size() - static_cast<std::size_t>(val.scale));
    PreciseDecimal res = PreciseDecimal::from_integer_string(integer_part, val.negative);
    if (val.negative) res -= PreciseDecimal(1LL);
    return res;
}

PreciseDecimal ceil(const PreciseDecimal& val) {
    if (val.scale <= 0) return val;
    PreciseDecimal f = floor(val);
    if (f == val) return val;
    return f + PreciseDecimal(1LL);
}

PreciseDecimal round(const PreciseDecimal& val) {
    if (val.negative) return precise::ceil(val - half());
    return precise::floor(val + half());
}

PreciseDecimal pi() {
    return PreciseDecimal(
        "3.141592653589793238462643383279502884197169399375105820974944592307816406286"
        "208998628034825342117067982148086513282306647093844609550582231725359408128481"
        "117450284102701938521105559644622948954930381964428810975665933446128475648233");
}

PreciseDecimal two_pi() {
    // 预计算的高精度 2π 常量，避免动态计算引入误差
    return PreciseDecimal(
        "6.283185307179586476925286766559005768394338798750211641949889184615632812572"
        "417997256069650684234135964296173026564613294187689219101164434507187816256962"
        "234900568205403877042211119289245897909860763928857621951331866892256951296466");
}

PreciseDecimal half_pi() {
    // 预计算的高精度 π/2 常量
    return PreciseDecimal(
        "1.5707963267948966192313216916397514420985846996875529104874722961539082031431"
        "0449931401741267105853399107404325664115332354692223077529111586285977040642405"
        "587251420513509692605527798223114744974651909822144054878329667323064237824116");
}

PreciseDecimal e() {
    return PreciseDecimal(
        "2.718281828459045235360287471352662497757247093699959574966967627724076630353"
        "547594571382178525166427427466391932003059921817413596629043572900334295260595"
        "630738132328627943490763233829880753195251019011573834187930702154089149934884");
}

PreciseDecimal exp(const PreciseDecimal& x) {
    if (x.is_zero()) return one();
    if (x.negative) return one() / exp(precise::abs(x));
    int k = 0;
    PreciseDecimal r = x;
    const PreciseDecimal threshold("0.01");
    while (r > threshold) {
        r /= two();
        k++;
    }
    PreciseDecimal term = one();
    PreciseDecimal sum = one();
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(80);
    for (int i = 1; i < limit; ++i) {
        term = term * r / decimal_from_uint(static_cast<uint32_t>(i));
        sum += term;
        if (precise::abs(term) < epsilon) break;
    }
    for (int i = 0; i < k; ++i) sum = sum * sum;
    return sum;
}

PreciseDecimal ln(const PreciseDecimal& x) {
    if (x <= PreciseDecimal(0LL)) throw std::domain_error("ln of non-positive number");
    if (x == one()) return PreciseDecimal(0LL);

    // 对于接近 1 的值，使用专门的 Taylor 级数展开
    // ln(1 + ε) = ε - ε²/2 + ε³/3 - ε⁴/4 + ...
    // 当 ε 很小时，这个级数收敛很快
    PreciseDecimal diff = x - one();
    int scale = PrecisionContext::get_default_scale();

    // 如果 |x - 1| < 1e-10，使用 Taylor 级数直接计算
    if (precise::abs(diff) < PreciseDecimal("1e-10")) {
        PreciseDecimal term = diff;
        PreciseDecimal sum = diff;
        PreciseDecimal diff2 = diff * diff;
        const PreciseDecimal epsilon = scale_epsilon();
        const int limit = series_iteration_limit(100);

        for (int i = 2; i < limit; ++i) {
            term = -term * diff / decimal_from_uint(static_cast<uint32_t>(i));
            sum += term;
            if (precise::abs(term) < epsilon) break;
        }
        return sum;
    }

    // 对于其他值，使用 Newton 迭代
    // 使用更高精度的初始猜测：通过 to_long_double() 获取更多精度
    PreciseDecimal y(mymath::Scalar(mymath::ln(x.to_long_double())));
    const int iterations = std::max(16, scale / 6 + 10);
    for (int i = 0; i < iterations; ++i) y = y + x / exp(y) - one();
    return y;
}

PreciseDecimal log10(const PreciseDecimal& x) {
    static const PreciseDecimal ln10 = ln(PreciseDecimal(10LL));
    return ln(x) / ln10;
}

PreciseDecimal sin(const PreciseDecimal& x) {
    if (x.is_zero()) return PreciseDecimal(0LL);
    PreciseDecimal r = reduce_trig_argument(x);

    PreciseDecimal term = r;
    PreciseDecimal sum = r;
    PreciseDecimal r2 = r * r;
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(80);
    for (int i = 1; i < limit; ++i) {
        term = -term * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i) * (2 * i + 1)));
        sum += term;
        if (precise::abs(term) < epsilon) break;
    }
    return sum;
}

PreciseDecimal cos(const PreciseDecimal& x) {
    if (x.is_zero()) return one();
    PreciseDecimal r = reduce_trig_argument(x);

    PreciseDecimal term = one();
    PreciseDecimal sum = one();
    PreciseDecimal r2 = r * r;
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(80);
    for (int i = 1; i < limit; ++i) {
        term = -term * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i - 1) * (2 * i)));
        sum += term;
        if (precise::abs(term) < epsilon) break;
    }
    return sum;
}

PreciseDecimal tan(const PreciseDecimal& x) {
    // 检查是否接近 pi/2 + k*pi
    // 使用较宽松的容差 (1e-12) 以捕捉从 long double 转换来的近似 pi/2
    PreciseDecimal p = pi();
    PreciseDecimal half_p = p / two();
    PreciseDecimal r = x - precise::floor(x / p) * p; // mod pi
    if (precise::abs(r - half_p) < PreciseDecimal("1e-12")) {
        throw std::domain_error("tan undefined near singularity");
    }
    PreciseDecimal s = sin(x);
    PreciseDecimal c = cos(x);
    if (c.is_zero()) throw std::domain_error("tan undefined (cos is zero)");
    return s / c;
}

PreciseDecimal asin(const PreciseDecimal& x) {
    if (precise::abs(x) > one()) throw std::domain_error("asin out of range");
    if (x == one()) return pi() / two();
    if (x == -one()) return -pi() / two();
    return atan(x / precise::sqrt(one() - x * x));
}

PreciseDecimal acos(const PreciseDecimal& x) {
    return pi() / two() - asin(x);
}

PreciseDecimal atan(const PreciseDecimal& x) {
    if (precise::abs(x) > one()) {
        if (x.negative) return -pi() / two() - atan(one() / x);
        return pi() / two() - atan(one() / x);
    }
    PreciseDecimal y(mymath::atan(x.to_double()));
    const int iterations = std::max(12, PrecisionContext::get_default_scale() / 8 + 8);
    for (int i = 0; i < iterations; ++i) {
        PreciseDecimal ty = tan(y);
        y = y - (ty - x) / (one() + ty * ty);
    }
    return y;
}

PreciseDecimal sinh(const PreciseDecimal& x) {
    PreciseDecimal ex = exp(x);
    return (ex - one() / ex) / two();
}

PreciseDecimal cosh(const PreciseDecimal& x) {
    PreciseDecimal ex = exp(x);
    return (ex + one() / ex) / two();
}

PreciseDecimal tanh(const PreciseDecimal& x) {
    PreciseDecimal ex = exp(x);
    PreciseDecimal einv = one() / ex;
    return (ex - einv) / (ex + einv);
}

PreciseDecimal asinh(const PreciseDecimal& x) {
    return ln(x + precise::sqrt(x * x + one()));
}

PreciseDecimal acosh(const PreciseDecimal& x) {
    if (x < one()) throw std::domain_error("acosh out of range");
    return ln(x + precise::sqrt(x * x - one()));
}

PreciseDecimal atanh(const PreciseDecimal& x) {
    if (precise::abs(x) >= one()) throw std::domain_error("atanh out of range");
    return ln((one() + x) / (one() - x)) / two();
}

} // namespace precise

PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (lhs.negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return subtract_precise_decimal(lhs, rhs_flipped);
    }

    PreciseDecimal l = lhs;
    PreciseDecimal r = rhs;
    align_precise_scales(&l, &r);

    PreciseDecimal res;
    res.data = add_bigint(l.data, r.data);
    res.scale = l.scale;
    res.negative = lhs.negative;
    res.normalize();
    return res;
}

PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (lhs.negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return add_precise_decimal(lhs, rhs_flipped);
    }

    PreciseDecimal l = lhs;
    PreciseDecimal r = rhs;
    align_precise_scales(&l, &r);

    int cmp = compare_bigint(l.data, r.data);
    if (cmp == 0) return {};

    PreciseDecimal res;
    res.scale = l.scale;
    if (cmp > 0) {
        res.data = subtract_bigint(l.data, r.data);
        res.negative = lhs.negative;
    } else {
        res.data = subtract_bigint(r.data, l.data);
        res.negative = !lhs.negative;
    }
    res.normalize();
    return res;
}

PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    PreciseDecimal res;
    res.data = multiply_bigint(lhs.data, rhs.data);
    res.scale = lhs.scale + rhs.scale;
    res.negative = lhs.negative != rhs.negative;
    res.normalize();
    return res;
}

PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (rhs.is_zero()) throw std::runtime_error("division by zero");
    if (lhs.is_zero()) return {};

    int target_scale = PrecisionContext::get_default_scale();
    int numerator_shift = target_scale + rhs.scale - lhs.scale;
    std::vector<uint32_t> numerator = lhs.data;
    if (numerator_shift >= 0) {
        numerator = multiply_bigint_by_power_of_10(numerator, numerator_shift);
    } else {
        std::vector<uint32_t> divisor = multiply_bigint_by_power_of_10({1}, -numerator_shift);
        std::vector<uint32_t> truncated, ignored_remainder;
        div_bigint(numerator, divisor, &truncated, &ignored_remainder);
        numerator = truncated;
    }
    std::vector<uint32_t> denominator = rhs.data;

    std::vector<uint32_t> q, r;
    div_bigint(numerator, denominator, &q, &r);

    PreciseDecimal res;
    res.data = q;
    res.scale = target_scale;
    res.negative = lhs.negative != rhs.negative;
    res.normalize();
    return res;
}

int compare_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (lhs.is_zero() && rhs.is_zero()) return 0;
    if (lhs.negative != rhs.negative) return lhs.negative ? -1 : 1;

    PreciseDecimal l = lhs;
    PreciseDecimal r = rhs;
    align_precise_scales(&l, &r);

    int res = compare_bigint(l.data, r.data);
    return lhs.negative ? -res : res;
}

PreciseDecimal pow_precise_decimal(const PreciseDecimal& base, long long exponent) {
    if (exponent < 0) throw PreciseDecimalUnsupported("negative exponent not supported for precise decimal");
    if (exponent == 0) return PreciseDecimal::from_integer_string("1", false);
    if (exponent == 1) return base;

    // 特殊情况优化
    if (base.is_zero()) return PreciseDecimal(0LL);
    if (base == one()) return one();
    if (base == two()) {
        // 2^n 可以通过位移快速计算
        std::vector<uint32_t> result = {1};
        for (long long i = 0; i < exponent; ++i) {
            result = multiply_bigint_by_uint32(result, 2);
        }
        PreciseDecimal res;
        res.data = result;
        res.scale = 0;
        res.negative = false;
        return res;
    }

    // 快速幂算法 (二分幂)
    PreciseDecimal res = PreciseDecimal::from_integer_string("1", false);
    PreciseDecimal b = base;
    long long exp = exponent;

    while (exp > 0) {
        if (exp & 1) {
            res = multiply_precise_decimal(res, b);
        }
        b = multiply_precise_decimal(b, b);
        exp >>= 1;
    }
    return res;
}

// ============================================================================
// 辅助函数
// ============================================================================

std::string rational_to_precise_decimal_text(const Rational& value) {
    PreciseDecimal numerator = PreciseDecimal::from_integer_string(
        std::to_string(value.numerator < 0 ? -value.numerator : value.numerator),
        value.numerator < 0);
    const PreciseDecimal denominator =
        PreciseDecimal::from_integer_string(std::to_string(value.denominator), false);
    return divide_precise_decimal(numerator, denominator).to_string();
}

std::string stored_value_precise_decimal_text(const StoredValue& value) {
    if (value.exact) return rational_to_precise_decimal_text(value.rational);
    if (value.has_precise_decimal_text) return value.precise_decimal_text;
    return format_decimal(normalize_display_decimal(value.decimal.to_long_double()), 15);
}

#pragma GCC diagnostic pop

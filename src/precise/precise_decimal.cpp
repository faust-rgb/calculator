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
#include <climits>
#include <cmath>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

// ============================================================================
// BigIntData 实现 (SSO)
// ============================================================================

BigIntData::BigIntData() : ptr(sso), size_(1), capacity_(SSO_CAP) { sso[0] = 0; }
BigIntData::BigIntData(size_t n, uint32_t val) : ptr(sso), size_(static_cast<uint32_t>(n)), capacity_(SSO_CAP) {
    if (n > SSO_CAP) {
        ptr = new uint32_t[n];
        capacity_ = static_cast<uint32_t>(n);
    }
    std::fill(ptr, ptr + size_, val);
}
BigIntData::BigIntData(const BigIntData& other) : ptr(sso), size_(other.size_), capacity_(SSO_CAP) {
    if (other.size_ > SSO_CAP) {
        ptr = new uint32_t[other.size_];
        capacity_ = other.size_;
    }
    std::copy(other.ptr, other.ptr + other.size_, ptr);
}
BigIntData::BigIntData(BigIntData&& other) noexcept : ptr(sso), size_(other.size_), capacity_(SSO_CAP) {
    if (other.ptr != other.sso) {
        ptr = other.ptr;
        capacity_ = other.capacity_;
        other.ptr = other.sso;
        other.capacity_ = SSO_CAP;
    } else {
        std::copy(other.sso, other.sso + other.size_, sso);
    }
    other.size_ = 0;
}
BigIntData::BigIntData(std::initializer_list<uint32_t> list) : ptr(sso), size_(static_cast<uint32_t>(list.size())), capacity_(SSO_CAP) {
    if (list.size() > SSO_CAP) {
        ptr = new uint32_t[list.size()];
        capacity_ = static_cast<uint32_t>(list.size());
    }
    std::copy(list.begin(), list.end(), ptr);
}
BigIntData::BigIntData(const uint32_t* first, const uint32_t* last) : ptr(sso), size_(0), capacity_(SSO_CAP) {
    size_t n = last - first;
    if (n > SSO_CAP) {
        ptr = new uint32_t[n];
        capacity_ = static_cast<uint32_t>(n);
    }
    size_ = static_cast<uint32_t>(n);
    std::copy(first, last, ptr);
}
BigIntData::~BigIntData() {
    if (ptr != sso) delete[] ptr;
}
BigIntData& BigIntData::operator=(const BigIntData& other) {
    if (this == &other) return *this;
    if (other.size_ > capacity_) {
        if (ptr != sso) delete[] ptr;
        ptr = new uint32_t[other.size_];
        capacity_ = other.size_;
    }
    size_ = other.size_;
    std::copy(other.ptr, other.ptr + other.size_, ptr);
    return *this;
}
BigIntData& BigIntData::operator=(BigIntData&& other) noexcept {
    if (this == &other) return *this;
    if (ptr != sso) delete[] ptr;
    if (other.ptr != other.sso) {
        ptr = other.ptr;
        capacity_ = other.capacity_;
        other.ptr = other.sso;
        other.capacity_ = SSO_CAP;
    } else {
        ptr = sso;
        capacity_ = SSO_CAP;
        std::copy(other.sso, other.sso + other.size_, sso);
    }
    size_ = other.size_;
    other.size_ = 0;
    return *this;
}
void BigIntData::push_back(uint32_t val) {
    if (size_ == capacity_) reserve(capacity_ * 2);
    ptr[size_++] = val;
}
void BigIntData::reserve(size_t n) {
    if (n <= capacity_) return;
    uint32_t* new_ptr = new uint32_t[n];
    std::copy(ptr, ptr + size_, new_ptr);
    if (ptr != sso) delete[] ptr;
    ptr = new_ptr;
    capacity_ = static_cast<uint32_t>(n);
}
void BigIntData::resize(size_t n, uint32_t val) {
    if (n > capacity_) reserve(n);
    if (n > size_) std::fill(ptr + size_, ptr + n, val);
    size_ = static_cast<uint32_t>(n);
}
void BigIntData::erase(const uint32_t* first, const uint32_t* last) {
    size_t count = last - first;
    size_t start_idx = first - ptr;
    if (count == 0) return;
    std::move(ptr + start_idx + count, ptr + size_, ptr + start_idx);
    size_ -= static_cast<uint32_t>(count);
}
void BigIntData::insert(uint32_t* pos, const uint32_t* first, const uint32_t* last) {
    size_t count = last - first;
    size_t offset = pos - ptr;
    if (size_ + count > capacity_) reserve(std::max(size_ + count, (size_t)capacity_ * 2));
    std::move_backward(ptr + offset, ptr + size_, ptr + size_ + count);
    std::copy(first, last, ptr + offset);
    size_ += static_cast<uint32_t>(count);
}
void BigIntData::insert(uint32_t* pos, size_t count, uint32_t val) {
    size_t offset = pos - ptr;
    if (size_ + count > capacity_) reserve(std::max(size_ + count, (size_t)capacity_ * 2));
    std::move_backward(ptr + offset, ptr + size_, ptr + size_ + count);
    std::fill(ptr + offset, ptr + offset + count, val);
    size_ += static_cast<uint32_t>(count);
}

namespace {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;  ///< 大整数基数（10^9）
constexpr int kBaseDigits = 9;          ///< 每个 chunk 的十进制位数
thread_local int g_default_scale = 20;  ///< 默认精度（小数点后位数）
constexpr long double kDisplayZeroEps = 1e-16;      ///< 显示零的阈值
constexpr long double kDisplayIntegerEps = 1e-9;    ///< 显示整数的阈值

thread_local bool g_suppress_normalization = false;

/**
 * @struct NormalizationSuppressor
 * @brief 自动压制规范化的 RAII 助手
 */
struct NormalizationSuppressor {
    bool old_state;
    NormalizationSuppressor() {
        old_state = g_suppress_normalization;
        g_suppress_normalization = true;
    }
    ~NormalizationSuppressor() {
        g_suppress_normalization = old_state;
    }
};

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
    char buf[128];
    // Use %.*Lg to format long double with dynamic precision
    int len = snprintf(buf, sizeof(buf), "%.*Lg", precision, value);
    if (len > 0 && len < static_cast<int>(sizeof(buf))) {
        return std::string(buf, len);
    }
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    return out.str();
}

/**
 * @brief 原始指针算术优化 (Raw Pointer Primitives)
 * 这些函数直接操作原始指针，避免 BigIntData 封装开销
 */

inline size_t add_raw(const uint32_t* a, size_t an, const uint32_t* b, size_t bn, uint32_t* res) {
    uint32_t carry = 0;
    size_t i = 0;
    size_t min_n = std::min(an, bn);
    for (; i < min_n; ++i) {
        uint64_t sum = static_cast<uint64_t>(a[i]) + b[i] + carry;
        res[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
    }
    const uint32_t* longer = (an > bn) ? a : b;
    size_t max_n = std::max(an, bn);
    for (; i < max_n; ++i) {
        uint64_t sum = static_cast<uint64_t>(longer[i]) + carry;
        res[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
    }
    if (carry) {
        res[i++] = carry;
    }
    return i;
}

inline size_t sub_raw(const uint32_t* a, size_t an, const uint32_t* b, size_t bn, uint32_t* res) {
    int32_t borrow = 0;
    size_t i = 0;
    for (; i < bn; ++i) {
        int64_t sub = static_cast<int64_t>(a[i]) - borrow - b[i];
        if (sub < 0) {
            sub += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res[i] = static_cast<uint32_t>(sub);
    }
    for (; i < an; ++i) {
        int64_t sub = static_cast<int64_t>(a[i]) - borrow;
        if (sub < 0) {
            sub += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res[i] = static_cast<uint32_t>(sub);
    }
    while (i > 1 && res[i - 1] == 0) i--;
    return i;
}

/**
 * @brief 原地加法：a += b，复用 a 的内存
 * @return 新的长度
 */
inline size_t add_inplace_raw(uint32_t* a, size_t an, const uint32_t* b, size_t bn, size_t capacity) {
    uint32_t carry = 0;
    size_t i = 0;
    size_t min_n = std::min(an, bn);
    for (; i < min_n; ++i) {
        uint64_t sum = static_cast<uint64_t>(a[i]) + b[i] + carry;
        a[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
    }
    // 处理 b 剩余部分（当 bn > an 时）
    for (; i < bn && i < capacity; ++i) {
        uint64_t sum = static_cast<uint64_t>(b[i]) + carry;
        a[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
    }
    // 处理 a 剩余部分（当 an > bn 时，需要继续处理进位）
    while (i < an && i < capacity) {
        uint64_t sum = static_cast<uint64_t>(a[i]) + carry;
        a[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
        ++i;
    }
    // 如果还有进位且容量足够
    while (carry && i < capacity) {
        a[i++] = carry % kBase;
        carry = carry / kBase;
    }
    return std::max(an, std::max(bn, i));
}

/**
 * @brief 原地减法：a -= b（要求 a >= b）
 * @return 新的长度
 */
inline size_t sub_inplace_raw(uint32_t* a, size_t an, const uint32_t* b, size_t bn) {
    int32_t borrow = 0;
    size_t i = 0;
    for (; i < bn; ++i) {
        int64_t sub = static_cast<int64_t>(a[i]) - borrow - b[i];
        if (sub < 0) {
            sub += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        a[i] = static_cast<uint32_t>(sub);
    }
    for (; i < an; ++i) {
        int64_t sub = static_cast<int64_t>(a[i]) - borrow;
        if (sub < 0) {
            sub += kBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        a[i] = static_cast<uint32_t>(sub);
    }
    // 去除前导零
    while (an > 1 && a[an - 1] == 0) --an;
    return an;
}

/**
 * @brief 原始指针比较
 * @return -1 表示 a < b，0 表示相等，1 表示 a > b
 */
inline int compare_raw(const uint32_t* a, size_t an, const uint32_t* b, size_t bn) {
    if (an != bn) return an < bn ? -1 : 1;
    for (size_t i = an; i > 0; --i) {
        if (a[i - 1] != b[i - 1]) {
            return a[i - 1] < b[i - 1] ? -1 : 1;
        }
    }
    return 0;
}

/**
 * @brief 原始指针乘以小整数
 * @return 新的长度
 */
inline size_t mul_uint32_raw(const uint32_t* a, size_t an, uint32_t b, uint32_t* res) {
    if (b == 0) { res[0] = 0; return 1; }
    uint32_t carry = 0;
    size_t i = 0;
    for (; i < an; ++i) {
        uint64_t cur = static_cast<uint64_t>(a[i]) * b + carry;
        res[i] = static_cast<uint32_t>(cur % kBase);
        carry = static_cast<uint32_t>(cur / kBase);
    }
    if (carry) res[i++] = carry;
    return i;
}

/**
 * @brief 原地乘以小整数
 */
inline size_t mul_uint32_inplace_raw(uint32_t* a, size_t an, uint32_t b, size_t capacity) {
    if (b == 0) { a[0] = 0; return 1; }
    uint32_t carry = 0;
    size_t i = 0;
    for (; i < an; ++i) {
        uint64_t cur = static_cast<uint64_t>(a[i]) * b + carry;
        a[i] = static_cast<uint32_t>(cur % kBase);
        carry = static_cast<uint32_t>(cur / kBase);
    }
    while (carry && i < capacity) {
        a[i++] = carry % kBase;
        carry /= kBase;
    }
    return i;
}

/**
 * @brief 计算 a / b，返回商和余数 (均为 BigInt)
 *
 * 这是 Knuth 算法 D 的声明，在后面实现。
 */
void div_bigint(const BigIntData& num,
                const BigIntData& den,
                BigIntData* quotient,
                BigIntData* remainder);

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
BigIntData string_to_bigint(const std::string& s) {
    if (s.empty() || s == "0") return {0};
    BigIntData res;
    res.clear();
    for (int i = static_cast<int>(s.size()); i > 0; i -= kBaseDigits) {
        if (i < kBaseDigits) {
            res.push_back(static_cast<uint32_t>(std::stoul(s.substr(0, static_cast<std::size_t>(i)))));
        } else {
            res.push_back(static_cast<uint32_t>(std::stoul(s.substr(static_cast<std::size_t>(i - kBaseDigits), kBaseDigits))));
        }
    }
    while (res.size() > 1 && res.back() == 0) res.pop_back();
    if (res.empty()) res.push_back(0);
    return res;
}

/**
 * @brief 将基数为 10^9 的大整数数组转换为字符串
 *
 * 数组是小端序（低位在前），输出字符串是高位在前。
 */
std::string bigint_to_string(const BigIntData& v) {
    if (v.empty() || (v.size() == 1 && v[0] == 0)) return "0";
    
    std::string s;
    s.reserve(v.size() * kBaseDigits);
    
    s.append(std::to_string(v.back()));
    for (int i = static_cast<int>(v.size()) - 2; i >= 0; --i) {
        std::string chunk = std::to_string(v[i]);
        if (chunk.size() < kBaseDigits) {
            s.append(kBaseDigits - chunk.size(), '0');
        }
        s.append(chunk);
    }
    return s;
}

/**
 * @brief 比较两个大整数的大小
 *
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return -1 表示 lhs < rhs，0 表示相等，1 表示 lhs > rhs
 */
int compare_bigint(const BigIntData& lhs, const BigIntData& rhs) {
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
BigIntData add_bigint(const BigIntData& lhs, const BigIntData& rhs) {
    BigIntData res;
    res.resize(std::max(lhs.size(), rhs.size()) + 1);
    size_t new_size = add_raw(lhs.data(), lhs.size(), rhs.data(), rhs.size(), res.data());
    res.resize(new_size);
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
BigIntData subtract_bigint(const BigIntData& lhs, const BigIntData& rhs) {
    BigIntData res;
    res.resize(lhs.size());
    size_t new_size = sub_raw(lhs.data(), lhs.size(), rhs.data(), rhs.size(), res.data());
    res.resize(new_size);
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
constexpr std::size_t KARATSUBA_THRESHOLD = 128;

/**
 * @brief NTT 乘法阈值
 *
 * 当数字长度大于此值时使用 NTT（数论变换）算法。
 * 2048 是经验值，在此规模以上 NTT 的 O(n log n) 优势显现。
 */
constexpr std::size_t NTT_THRESHOLD = 2048;

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
     * @param n 数组长度（必须是 2 的幂）
     * @param invert true 表示逆变换，false 表示正变换
     */
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
 * @param n 数组长度
 * @param base 输出基数（这里是 10^9）
 * @return 合并后的大整数数组
 */
BigIntData crt(const uint32_t* r1,
               const uint32_t* r2,
               const uint32_t* r3,
               std::size_t n,
               uint32_t base) {
    BigIntData res;
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

    std::copy(lhs.begin(), lhs.end(), fa1);
    std::copy(lhs.begin(), lhs.end(), fa2);
    std::copy(lhs.begin(), lhs.end(), fa3);
    std::copy(rhs.begin(), rhs.end(), fb1);
    std::copy(rhs.begin(), rhs.end(), fb2);
    std::copy(rhs.begin(), rhs.end(), fb3);

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

    BigIntData result = ntt::crt(fa1, fa2, fa3, n, 1000000000); // kBase = 10^9
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
BigIntData multiply_bigint_karatsuba(const BigIntData& lhs, const BigIntData& rhs);

/**
 * @brief 大整数乘法 - Toom-Cook 3 算法
 *
 * 将数字分成 3 部分进行分治乘法。
 * 时间复杂度 O(n^1.465)，比 Karatsuba 更快。
 * 适合中等规模（256-4096 chunks）的乘法。
 */
BigIntData multiply_bigint_toom3(const BigIntData& lhs, const BigIntData& rhs);

/**
 * @brief 大整数乘以一个小整数
 *
 * 当乘数是单个 uint32_t 时使用此函数，效率更高。
 * 使用原始指针 API 优化。
 *
 * @param v 大整数数组
 * @param n 小整数乘数
 * @return 乘积
 */
BigIntData multiply_bigint_by_uint32(const BigIntData& v, uint32_t n);

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
// SIMD 优化的朴素乘法（循环展开）
// ============================================================================

/**
 * @brief SIMD 优化的朴素乘法（使用循环展开）
 *
 * 利用循环展开并行处理多个乘法运算。
 * 适合中等规模（32-256 chunks）的乘法。
 */
BigIntData multiply_bigint_naive_simd(const BigIntData& lhs, const BigIntData& rhs) {
    if ((lhs.size() == 1 && lhs[0] == 0) || (rhs.size() == 1 && rhs[0] == 0)) return {0};

    const std::size_t lhs_size = lhs.size();
    const std::size_t rhs_size = rhs.size();
    BigIntData res(lhs_size + rhs_size, 0);

    const uint32_t* __restrict__ lhs_ptr = lhs.data();
    const uint32_t* __restrict__ rhs_ptr = rhs.data();
    uint32_t* __restrict__ res_ptr = res.data();

    // 使用循环展开优化
    for (std::size_t i = 0; i < lhs_size; ++i) {
        const uint32_t l_val = lhs_ptr[i];
        if (l_val == 0) continue;

        uint64_t carry = 0;
        std::size_t j = 0;

        // 主循环：每次处理 4 个元素
        const std::size_t rhs_size_4 = rhs_size - (rhs_size % 4);
        for (; j < rhs_size_4; j += 4) {
            uint64_t cur0 = res_ptr[i + j] + static_cast<uint64_t>(l_val) * rhs_ptr[j] + carry;
            uint64_t cur1 = res_ptr[i + j + 1] + static_cast<uint64_t>(l_val) * rhs_ptr[j + 1] + (cur0 / kBase);
            uint64_t cur2 = res_ptr[i + j + 2] + static_cast<uint64_t>(l_val) * rhs_ptr[j + 2] + (cur1 / kBase);
            uint64_t cur3 = res_ptr[i + j + 3] + static_cast<uint64_t>(l_val) * rhs_ptr[j + 3] + (cur2 / kBase);

            res_ptr[i + j] = static_cast<uint32_t>(cur0 % kBase);
            res_ptr[i + j + 1] = static_cast<uint32_t>(cur1 % kBase);
            res_ptr[i + j + 2] = static_cast<uint32_t>(cur2 % kBase);
            res_ptr[i + j + 3] = static_cast<uint32_t>(cur3 % kBase);

            carry = cur3 / kBase;
        }

        // 处理剩余元素
        for (; j < rhs_size; ++j) {
            uint64_t cur = res_ptr[i + j] + static_cast<uint64_t>(l_val) * rhs_ptr[j] + carry;
            res_ptr[i + j] = static_cast<uint32_t>(cur % kBase);
            carry = cur / kBase;
        }

        // 处理进位
        if (carry) {
            std::size_t k = i + rhs_size;
            while (carry && k < res.size()) {
                uint64_t cur = res_ptr[k] + carry;
                res_ptr[k] = static_cast<uint32_t>(cur % kBase);
                carry = cur / kBase;
                ++k;
            }
        }
    }

    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

/**
 * @brief 大整数乘以一个小整数
 *
 * 当乘数是单个 uint32_t 时使用此函数，效率更高。
 * 使用原始指针 API 优化。
 *
 * @param v 大整数数组
 * @param n 小整数乘数
 * @return 乘积
 */
BigIntData multiply_bigint_by_uint32(const BigIntData& v, uint32_t n) {
    if (n == 0 || (v.size() == 1 && v[0] == 0)) return {0};
    if (n == 1) return v;
    BigIntData res;
    res.resize(v.size() + 1);
    size_t new_size = mul_uint32_raw(v.data(), v.size(), n, res.data());
    res.resize(new_size);
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
BigIntData shift_bigint(const BigIntData& v, std::size_t n) {
    if (v.empty() || (v.size() == 1 && v[0] == 0)) return {0};
    if (n == 0) return v;
    BigIntData res(n, 0);
    res.insert(res.end(), v.begin(), v.end());
    return res;
}

/**
 * @brief 统计大整数末尾的零的个数
 */
int count_trailing_zeros(const BigIntData& v) {
    if (v.empty() || (v.size() == 1 && v[0] == 0)) return 0;
    int count = 0;
    std::size_t i = 0;
    while (i < v.size() && v[i] == 0) {
        count += kBaseDigits;
        i++;
    }
    if (i < v.size()) {
        uint32_t val = v[i];
        while (val > 0 && val % 10 == 0) {
            count++;
            val /= 10;
        }
    }
    return count;
}

/**
 * @brief 大整数除以 10^n (原地修改或返回新向量)
 */
BigIntData divide_bigint_by_pow10(BigIntData v, int n) {
    if (n <= 0) return v;
    int chunk_shift = n / kBaseDigits;
    int digit_shift = n % kBaseDigits;

    if (chunk_shift >= static_cast<int>(v.size())) return {0};

    if (chunk_shift > 0) {
        v.erase(v.begin(), v.begin() + chunk_shift);
    }

    if (digit_shift > 0) {
        uint32_t divisor = 1;
        for (int i = 0; i < digit_shift; ++i) divisor *= 10;

        uint32_t rem = 0;
        for (int i = static_cast<int>(v.size()) - 1; i >= 0; --i) {
            uint64_t current = v[i] + static_cast<uint64_t>(rem) * kBase;
            v[i] = static_cast<uint32_t>(current / divisor);
            rem = static_cast<uint32_t>(current % divisor);
        }
    }
    while (v.size() > 1 && v.back() == 0) v.pop_back();
    if (v.empty()) v = {0};
    return v;
}

/**
 * @brief 检查大整数是否能被 10^n 整除
 */
bool is_bigint_multiple_of_pow10(const BigIntData& v, int n) {
    if (n <= 0) return true;
    return count_trailing_zeros(v) >= n;
}

/**
 * @brief 大整数截取低位 n 个 chunk
 *
 * @param v 大整数数组
 * @param n 要截取的 chunk 数量
 * @return 截取后的数组
 */
BigIntData bigint_low(const BigIntData& v, std::size_t n) {
    if (n >= v.size()) return v;
    return BigIntData(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(n));
}

/**
 * @brief 大整数截取高位（从第 n 个 chunk 开始）
 *
 * @param v 大整数数组
 * @param n 起始 chunk 索引
 * @return 截取后的数组
 */
BigIntData bigint_high(const BigIntData& v, std::size_t n) {
    if (n >= v.size()) return {0};
    return BigIntData(v.begin() + static_cast<std::ptrdiff_t>(n), v.end());
}

/**
 * @brief 移除大整数前导零
 *
 * 原地修改，移除数组末尾的所有零（保留至少一个元素）。
 *
 * @param v 大整数数组（原地修改）
 */
void trim_bigint(BigIntData& v) {
    while (v.size() > 1 && v.back() == 0) v.pop_back();
}

/**
 * @brief 大整数乘法 - Karatsuba 算法实现
 *
 * 递归实现 Karatsuba 乘法，小规模时退化为朴素乘法。
 */
BigIntData multiply_bigint_karatsuba(const BigIntData& lhs, const BigIntData& rhs) {
    std::size_t n = std::max(lhs.size(), rhs.size());

    // 基准情况：小规模使用朴素乘法
    if (n <= KARATSUBA_THRESHOLD) {
        return multiply_bigint_naive(lhs, rhs);
    }

    std::size_t m = n / 2;

    // 分解: lhs = a1 * B^m + a0, rhs = b1 * B^m + b0
    BigIntData a0 = bigint_low(lhs, m);
    BigIntData a1 = bigint_high(lhs, m);
    BigIntData b0 = bigint_low(rhs, m);
    BigIntData b1 = bigint_high(rhs, m);

    trim_bigint(a0);
    trim_bigint(a1);
    trim_bigint(b0);
    trim_bigint(b1);

    // 递归计算三个乘积
    BigIntData z0 = multiply_bigint_karatsuba(a0, b0);  // a0 * b0
    BigIntData z2 = multiply_bigint_karatsuba(a1, b1);  // a1 * b1

    // z1 = (a0 + a1) * (b0 + b1) - z0 - z2
    BigIntData a0_plus_a1 = add_bigint(a0, a1);
    BigIntData b0_plus_b1 = add_bigint(b0, b1);
    BigIntData z1_temp = multiply_bigint_karatsuba(a0_plus_a1, b0_plus_b1);
    BigIntData z1 = subtract_bigint(subtract_bigint(z1_temp, z0), z2);

    // 结果 = z2 * B^(2m) + z1 * B^m + z0
    BigIntData result = add_bigint(
        add_bigint(shift_bigint(z2, 2 * m), shift_bigint(z1, m)),
        z0
    );

    trim_bigint(result);
    return result;
}

/**
 * @brief 大整数乘法 - Karatsuba 算法（并行版本）
 *
 * 使用 OpenMP 并行计算三个子乘积。
 * 适合大规模（> 512 chunks）的乘法。
 */
BigIntData multiply_bigint_karatsuba_parallel(const BigIntData& lhs, const BigIntData& rhs) {
    std::size_t n = std::max(lhs.size(), rhs.size());

    // 基准情况：小规模使用普通 Karatsuba
    if (n <= 512) {
        return multiply_bigint_karatsuba(lhs, rhs);
    }

    std::size_t m = n / 2;

    // 分解
    BigIntData a0 = bigint_low(lhs, m);
    BigIntData a1 = bigint_high(lhs, m);
    BigIntData b0 = bigint_low(rhs, m);
    BigIntData b1 = bigint_high(rhs, m);

    trim_bigint(a0);
    trim_bigint(a1);
    trim_bigint(b0);
    trim_bigint(b1);

    BigIntData z0, z1_temp, z2;
    BigIntData a0_plus_a1 = add_bigint(a0, a1);
    BigIntData b0_plus_b1 = add_bigint(b0, b1);

    // 并行计算三个乘积
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            z0 = multiply_bigint_karatsuba(a0, b0);
        }
        #pragma omp section
        {
            z2 = multiply_bigint_karatsuba(a1, b1);
        }
        #pragma omp section
        {
            z1_temp = multiply_bigint_karatsuba(a0_plus_a1, b0_plus_b1);
        }
    }

    // z1 = z1_temp - z0 - z2
    BigIntData z1 = subtract_bigint(subtract_bigint(z1_temp, z0), z2);

    // 结果 = z2 * B^(2m) + z1 * B^m + z0
    BigIntData result = add_bigint(
        add_bigint(shift_bigint(z2, 2 * m), shift_bigint(z1, m)),
        z0
    );

    trim_bigint(result);
    return result;
}

// ============================================================================
// Toom-Cook 3 乘法算法
// ============================================================================

/**
 * @brief 大整数除以 3（用于 Toom-Cook）
 */
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

/**
 * @brief 大整数乘法 - Toom-Cook 3 算法
 *
 * 将数字分成 3 部分进行分治乘法。
 * 时间复杂度 O(n^1.465)，比 Karatsuba 更快。
 * 适合中等规模（256-4096 chunks）的乘法。
 *
 * 算法原理：
 * 将 x, y 分成 3 部分：
 *   x = x2 * B^(2m) + x1 * B^m + x0
 *   y = y2 * B^(2m) + y1 * B^m + y0
 *
 * 计算 5 个点值（在 0, 1, -1, 2, -2 处求值）：
 *   m(0) = x0 * y0
 *   m(1) = (x0+x1+x2) * (y0+y1+y2)
 *   m(-1) = (x0-x1+x2) * (y0-y1+y2)
 *   m(2) = (x0+2x1+4x2) * (y0+2y1+4y2)
 *   m(-2) = (x0-2x1+4x2) * (y0-2y1+4y2)
 *   m(∞) = x2 * y2
 *
 * 然后通过插值得到结果。
 */
BigIntData multiply_bigint_toom3(const BigIntData& lhs, const BigIntData& rhs) {
    std::size_t n = std::max(lhs.size(), rhs.size());

    // 基准情况：小规模使用 Karatsuba
    if (n <= 256) {
        return multiply_bigint_karatsuba(lhs, rhs);
    }

    // 分割点
    std::size_t m = (n + 2) / 3;

    // 分解为 3 部分
    BigIntData x0 = bigint_low(lhs, m);
    BigIntData x1 = bigint_low(bigint_high(lhs, m), m);
    BigIntData x2 = bigint_high(lhs, 2 * m);

    BigIntData y0 = bigint_low(rhs, m);
    BigIntData y1 = bigint_low(bigint_high(rhs, m), m);
    BigIntData y2 = bigint_high(rhs, 2 * m);

    trim_bigint(x0); trim_bigint(x1); trim_bigint(x2);
    trim_bigint(y0); trim_bigint(y1); trim_bigint(y2);

    // 求值点计算
    // p(0) = x0
    // p(1) = x0 + x1 + x2
    // p(-1) = x0 - x1 + x2
    // p(2) = x0 + 2*x1 + 4*x2
    // p(∞) = x2

    BigIntData x0_plus_x2 = add_bigint(x0, x2);
    BigIntData y0_plus_y2 = add_bigint(y0, y2);

    // m0 = x0 * y0
    BigIntData m0 = multiply_bigint_toom3(x0, y0);

    // m_inf = x2 * y2
    BigIntData m_inf = multiply_bigint_toom3(x2, y2);

    // m1 = (x0 + x1 + x2) * (y0 + y1 + y2)
    BigIntData p1 = add_bigint(x0_plus_x2, x1);
    BigIntData q1 = add_bigint(y0_plus_y2, y1);
    BigIntData m1 = multiply_bigint_toom3(p1, q1);

    // m_neg1 = (x0 - x1 + x2) * (y0 - y1 + y2)
    BigIntData p_neg1, q_neg1;
    bool p_neg1_neg = false, q_neg1_neg = false;

    int cmp_x1 = compare_bigint(x0_plus_x2, x1);
    if (cmp_x1 >= 0) {
        p_neg1 = subtract_bigint(x0_plus_x2, x1);
    } else {
        p_neg1 = subtract_bigint(x1, x0_plus_x2);
        p_neg1_neg = true;
    }

    int cmp_y1 = compare_bigint(y0_plus_y2, y1);
    if (cmp_y1 >= 0) {
        q_neg1 = subtract_bigint(y0_plus_y2, y1);
    } else {
        q_neg1 = subtract_bigint(y1, y0_plus_y2);
        q_neg1_neg = true;
    }

    BigIntData m_neg1 = multiply_bigint_toom3(p_neg1, q_neg1);
    bool m_neg1_neg = p_neg1_neg != q_neg1_neg;

    // m2 = (x0 + 2*x1 + 4*x2) * (y0 + 2*y1 + 4*y2)
    BigIntData x1_2 = multiply_bigint_by_uint32(x1, 2);
    BigIntData x2_4 = multiply_bigint_by_uint32(x2, 4);
    BigIntData y1_2 = multiply_bigint_by_uint32(y1, 2);
    BigIntData y2_4 = multiply_bigint_by_uint32(y2, 4);

    BigIntData p2 = add_bigint(add_bigint(x0, x1_2), x2_4);
    BigIntData q2 = add_bigint(add_bigint(y0, y1_2), y2_4);
    BigIntData m2 = multiply_bigint_toom3(p2, q2);

    // 插值计算结果系数
    // c0 = m0
    // c4 = m_inf
    // c3 = (m2 - m1) / 3 - (m_neg1 - m0) / 2 - 2*m_inf
    // c1 = (m1 - m_neg1) / 2
    // c2 = m1 - m0 - m_inf - c1

    BigIntData c0 = m0;
    BigIntData c4 = m_inf;

    // c1 = (m1 - m_neg1) / 2
    BigIntData c1;
    if (m_neg1_neg) {
        c1 = divide_bigint_by_pow10(add_bigint(m1, m_neg1), 1);
    } else {
        BigIntData diff = subtract_bigint(m1, m_neg1);
        c1 = divide_bigint_by_pow10(diff, 1);
    }

    // c2 = m1 - m0 - m_inf - c1
    BigIntData c2 = subtract_bigint(subtract_bigint(subtract_bigint(m1, m0), m_inf), c1);

    // c3 = (m2 - m1) / 3 - (m_neg1 - m0) / 2 - 2*m_inf
    BigIntData m2_minus_m1 = subtract_bigint(m2, m1);
    BigIntData m2_minus_m1_div3 = divide_bigint_by_3(m2_minus_m1);

    BigIntData m_neg1_minus_m0;
    if (m_neg1_neg) {
        m_neg1_minus_m0 = add_bigint(m_neg1, m0);
    } else {
        m_neg1_minus_m0 = subtract_bigint(m_neg1, m0);
    }
    BigIntData half_diff = divide_bigint_by_pow10(m_neg1_minus_m0, 1);

    BigIntData two_m_inf = multiply_bigint_by_uint32(m_inf, 2);
    BigIntData c3 = subtract_bigint(subtract_bigint(m2_minus_m1_div3, half_diff), two_m_inf);

    // 结果 = c0 + c1*B^m + c2*B^(2m) + c3*B^(3m) + c4*B^(4m)
    BigIntData result = c0;
    result = add_bigint(result, shift_bigint(c1, m));
    result = add_bigint(result, shift_bigint(c2, 2 * m));
    result = add_bigint(result, shift_bigint(c3, 3 * m));
    result = add_bigint(result, shift_bigint(c4, 4 * m));

    trim_bigint(result);
    return result;
}

/**
 * @brief 大整数乘法（自动选择最优算法）
 *
 * 根据输入规模自动选择最适合的乘法算法：
 * - 小规模（<= 128 chunks）：朴素乘法
 * - 中等规模（<= 512 chunks）：Karatsuba 算法
 * - 中大规模（<= 4096 chunks）：Toom-Cook 3 算法
 * - 大规模（> 4096 chunks）：NTT 算法
 *
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 乘积
 */
BigIntData multiply_bigint(const BigIntData& lhs, const BigIntData& rhs) {
    // 零值快速返回
    if ((lhs.size() == 1 && lhs[0] == 0) || (rhs.size() == 1 && rhs[0] == 0)) return {0};

    std::size_t max_size = std::max(lhs.size(), rhs.size());
    std::size_t min_size = std::min(lhs.size(), rhs.size());

    // 如果其中一个操作数很小（长乘短），或者两者规模都很小，直接使用朴素乘法
    // 朴素乘法 O(N*M) 在 M 很小时比 Karatsuba/NTT 更快
    if (min_size <= 32 || max_size <= KARATSUBA_THRESHOLD) {
        return multiply_bigint_naive(lhs, rhs);
    }

    // 中等规模使用 Karatsuba
    if (max_size <= 512) {
        return multiply_bigint_karatsuba(lhs, rhs);
    }

    // Toom-Cook interpolation uses signed intermediate values; until that path
    // is fully sign-aware, keep medium-size products on the stable Karatsuba path.
    if (max_size <= 4096) {
        return multiply_bigint_karatsuba(lhs, rhs);
    }

    // 大规模使用 NTT
    return multiply_bigint_ntt(lhs, rhs);
}

// ============================================================================
// Newton-Raphson 除法（大规模除法优化）
// ============================================================================

/**
 * @brief 计算大整数的近似倒数（Newton-Raphson 方法）
 *
 * 使用迭代公式：x_{n+1} = x_n * (2 - d * x_n)
 * 收敛速度为二次收敛。
 *
 * @param d 除数
 * @param precision 所需精度（位数）
 * @return 1/d 的近似值（放大 2^precision 倍）
 */
BigIntData reciprocal_newton(const BigIntData& d, std::size_t precision) {
    if (d.empty() || (d.size() == 1 && d[0] == 0)) {
        throw std::runtime_error("division by zero in reciprocal_newton");
    }

    // 初始近似：使用浮点数
    std::size_t d_bits = d.size() * 32;  // 估计位数

    // 初始猜测：x0 = 1/d 的高位部分
    // 使用 64 位浮点倒数作为起点
    uint64_t d_high = 0;
    if (d.size() >= 2) {
        d_high = (static_cast<uint64_t>(d[d.size() - 1]) << 32) | d[d.size() - 2];
    } else {
        d_high = d[0];
    }

    // 初始近似：x0 ≈ 2^64 / d_high
    uint64_t x0 = (static_cast<uint64_t>(1) << 63) / (d_high >> 1);

    BigIntData x;
    x.push_back(static_cast<uint32_t>(x0));
    x.push_back(static_cast<uint32_t>(x0 >> 32));

    // 迭代直到达到所需精度
    std::size_t current_precision = 64;
    BigIntData two = {0, 0, 1};  // 2^64

    while (current_precision < precision) {
        // x = x * (2 - d * x) = 2*x - d*x*x
        BigIntData d_times_x = multiply_bigint(d, x);
        BigIntData x_squared = multiply_bigint(x, x);
        BigIntData d_times_x_squared = multiply_bigint(d, x_squared);

        // 2*x (左移)
        BigIntData two_x;
        two_x.reserve(x.size() + 1);
        uint64_t carry = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            uint64_t val = static_cast<uint64_t>(x[i]) * 2 + carry;
            two_x.push_back(static_cast<uint32_t>(val));
            carry = val >> 32;
        }
        if (carry) two_x.push_back(static_cast<uint32_t>(carry));

        // x = two_x - d_times_x_squared（需要调整精度）
        if (compare_bigint(two_x, d_times_x_squared) >= 0) {
            x = subtract_bigint(two_x, d_times_x_squared);
        } else {
            // 结果为负，需要调整
            x = subtract_bigint(d_times_x_squared, two_x);
        }

        current_precision *= 2;
    }

    // 调整到所需精度
    if (x.size() > (precision + 31) / 32) {
        x.resize((precision + 31) / 32);
    }

    return x;
}

/**
 * @brief Barrett 约减器
 *
 * 预计算除数的倒数近似值，用于快速重复除法。
 * 将除法转换为乘法，复杂度从 O(n*m) 降为 O(n*log(m))。
 */
struct BarrettReducer {
    BigIntData divisor;       ///< 除数
    BigIntData reciprocal;    ///< floor(2^(2n) / divisor)
    std::size_t shift;        ///< 2n，其中 n 是除数的位数

    /**
     * @brief 构造 Barrett 约减器
     * @param d 除数
     */
    explicit BarrettReducer(const BigIntData& d) : divisor(d) {
        if (d.empty() || (d.size() == 1 && d[0] == 0)) {
            throw std::runtime_error("BarrettReducer: divisor cannot be zero");
        }

        // 计算位数
        std::size_t n = divisor.size();

        // 计算 2^(2n) / divisor
        shift = 2 * n * 32;  // 每个 chunk 32 位

        // 2^(2n) 表示为 BigIntData
        BigIntData two_pow_2n(shift / 32 + 1, 0);
        two_pow_2n[shift / 32] = static_cast<uint32_t>(1) << (shift % 32);
        if (shift % 32 == 0) {
            two_pow_2n.resize(shift / 32);
            two_pow_2n.push_back(1);
        }

        // 计算 reciprocal = 2^(2n) / divisor
        BigIntData rem;
        div_bigint(two_pow_2n, divisor, &reciprocal, &rem);
    }

    /**
     * @brief 快速约减：计算 x mod divisor
     * @param x 被约减数
     * @return x mod divisor
     */
    BigIntData reduce(const BigIntData& x) const {
        if (compare_bigint(x, divisor) < 0) {
            return x;
        }

        // q = floor(x / divisor) ≈ floor(x * reciprocal / 2^shift)
        BigIntData x_times_r = multiply_bigint(x, reciprocal);

        // 右移 shift 位
        std::size_t limb_shift = shift / 32;
        if (limb_shift >= x_times_r.size()) {
            return x;  // x < divisor
        }

        BigIntData q(x_times_r.begin() + static_cast<std::ptrdiff_t>(limb_shift), x_times_r.end());

        // r = x - q * divisor
        BigIntData q_times_d = multiply_bigint(q, divisor);
        BigIntData r;
        if (compare_bigint(x, q_times_d) >= 0) {
            r = subtract_bigint(x, q_times_d);
        } else {
            // 修正：q 偏大
            q = subtract_bigint(q, BigIntData{1});
            q_times_d = multiply_bigint(q, divisor);
            r = subtract_bigint(x, q_times_d);
        }

        // 最终修正
        while (compare_bigint(r, divisor) >= 0) {
            r = subtract_bigint(r, divisor);
        }

        return r;
    }

    /**
     * @brief 快速除法：计算 x / divisor
     * @param x 被除数
     * @return 商
     */
    BigIntData divide(const BigIntData& x) const {
        if (compare_bigint(x, divisor) < 0) {
            return {0};
        }

        // q ≈ floor(x * reciprocal / 2^shift)
        BigIntData x_times_r = multiply_bigint(x, reciprocal);

        std::size_t limb_shift = shift / 32;
        if (limb_shift >= x_times_r.size()) {
            return {0};
        }

        BigIntData q(x_times_r.begin() + static_cast<std::ptrdiff_t>(limb_shift), x_times_r.end());

        // 修正 q
        BigIntData q_times_d = multiply_bigint(q, divisor);
        BigIntData r;
        if (compare_bigint(x, q_times_d) >= 0) {
            r = subtract_bigint(x, q_times_d);
        } else {
            q = subtract_bigint(q, BigIntData{1});
            q_times_d = multiply_bigint(q, divisor);
        }

        // 检查是否需要进一步修正
        while (compare_bigint(multiply_bigint(add_bigint(q, BigIntData{1}), divisor), x) <= 0) {
            q = add_bigint(q, BigIntData{1});
        }

        return q;
    }
};

/**
 * @brief 使用 Newton-Raphson 方法的大整数除法
 *
 * 适合大规模除法（除数位数 > 1000）。
 * 时间复杂度 O(n log n)。
 *
 * @param num 被除数
 * @param den 除数
 * @param quotient 输出：商
 * @param remainder 输出：余数
 */
void div_bigint_newton(const BigIntData& num,
                       const BigIntData& den,
                       BigIntData* quotient,
                       BigIntData* remainder) {
    if (den.empty() || (den.size() == 1 && den[0] == 0)) {
        throw std::runtime_error("division by zero in div_bigint_newton");
    }

    // 小规模使用 Knuth 算法
    if (den.size() < 64 || num.size() < 64) {
        div_bigint(num, den, quotient, remainder);
        return;
    }

    // 使用 Barrett 约减
    BarrettReducer reducer(den);
    *quotient = reducer.divide(num);
    *remainder = reducer.reduce(num);
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
void div_bigint(const BigIntData& num,
                const BigIntData& den,
                BigIntData* quotient,
                BigIntData* remainder) {
    std::size_t m = den.size();
    while (m > 0 && den[m - 1] == 0) --m;
    if (m == 0) {
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
    // 除数大于被除数
    if (m > n || (m == n && compare_bigint(
            BigIntData(num.begin(), num.begin() + static_cast<std::ptrdiff_t>(n)),
            BigIntData(den.begin(), den.begin() + static_cast<std::ptrdiff_t>(m))) < 0)) {
        *quotient = {0};
        *remainder = BigIntData(num.begin(), num.begin() + static_cast<std::ptrdiff_t>(n));
        while (remainder->size() > 1 && remainder->back() == 0) remainder->pop_back();
        return;
    }

    // 特殊情况：除数为单 chunk
    if (m == 1) {
        BigIntData q(n, 0);
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
    uint32_t d = static_cast<uint32_t>(kBase / (static_cast<uint64_t>(den[m - 1]) + 1));
    if (d == 0) d = 1;

    // 归一化被除数和除数
    BigIntData u(n + 1, 0);  // 被除数，多一位防止溢出
    BigIntData v(m, 0);       // 除数

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
    BigIntData q(n - m + 1, 0);

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
    BigIntData rem(m, 0);
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
BigIntData multiply_bigint_by_power_of_10(BigIntData v, int n) {
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
        v.insert(v.begin(), static_cast<std::size_t>(chunk_shift), 0u);
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

long double long_double_power_of_10(int exponent) {
    long double result = 1.0L;
    long double factor = 10.0L;
    int n = exponent < 0 ? -exponent : exponent;

    while (n > 0) {
        if (n % 2 == 1) result *= factor;
        factor *= factor;
        n /= 2;
    }

    return exponent < 0 ? 1.0L / result : result;
}

PreciseDecimal decimal_power_of_10(int exponent) {
    PreciseDecimal result(1LL);
    if (exponent >= 0) {
        result.data = multiply_bigint_by_power_of_10(result.data, exponent);
    } else {
        result.scale = -exponent;
    }
    result.normalize();
    return result;
}

PreciseDecimal sqrt_initial_guess(const PreciseDecimal& val) {
    std::string digits = bigint_to_string(val.data);
    int decimal_exponent = static_cast<int>(digits.size()) - val.scale - 1;
    int mantissa_digits = std::min<int>(18, digits.size());
    std::string mantissa_text = digits.substr(0, static_cast<std::size_t>(mantissa_digits));

    PreciseDecimal mantissa = PreciseDecimal::from_digits(mantissa_text, mantissa_digits - 1, false);
    if (decimal_exponent % 2 != 0) {
        mantissa *= PreciseDecimal(10LL);
        --decimal_exponent;
    }

    PreciseDecimal x = mantissa > one() ? mantissa : one();
    const PreciseDecimal one_half = half();
    for (int i = 0; i < 8; ++i) {
        x = one_half * (x + mantissa / x);
    }

    return x * decimal_power_of_10(decimal_exponent / 2);
}

PreciseDecimal ln_near_one(const PreciseDecimal& x) {
    const PreciseDecimal z = (x - one()) / (x + one());
    const PreciseDecimal z2 = z * z;
    
    NormalizationSuppressor suppressor;

    PreciseDecimal term = z;
    PreciseDecimal sum = z;
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(200);

    for (int n = 3; n < limit * 2; n += 2) {
        term *= z2;
        PreciseDecimal add = term / decimal_from_uint(static_cast<uint32_t>(n));
        sum += add;
        if (precise::abs(add) < epsilon) break;
    }

    g_suppress_normalization = false;
    sum.normalize();
    return two() * sum;
}

bool try_integral_exponent(const PreciseDecimal& value, long long* exponent) {
    PreciseDecimal rounded = precise::round(value);
    if (precise::abs(value - rounded) > scale_epsilon(2)) return false;
    *exponent = static_cast<long long>(rounded.to_double());
    return true;
}

bool try_thirds_exponent(const PreciseDecimal& value, long long* thirds) {
    PreciseDecimal scaled = value * PreciseDecimal(3LL);
    PreciseDecimal rounded = precise::round(scaled);
    if (precise::abs(scaled - rounded) > PreciseDecimal("1e-25")) return false;
    *thirds = static_cast<long long>(rounded.to_double());
    return true;
}

PreciseDecimal cbrt_precise_decimal(const PreciseDecimal& val) {
    if (val.is_zero()) return PreciseDecimal(0LL);
    PreciseDecimal x = sqrt_initial_guess(precise::abs(val));
    if (x.is_zero()) x = one();
    const PreciseDecimal abs_val = precise::abs(val);
    for (int i = 0; i < std::max(14, PrecisionContext::get_default_scale() / 6 + 8); ++i) {
        x = (two() * x + abs_val / (x * x)) / PreciseDecimal(3LL);
    }
    return val.negative ? -x : x;
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

PreciseDecimal reduce_sin_argument(const PreciseDecimal& x, bool* negate) {
    PreciseDecimal r = reduce_trig_argument(x);
    const PreciseDecimal p = precise::pi();
    const PreciseDecimal half_p = precise::half_pi();
    *negate = false;
    if (r > half_p) {
        r = p - r;
    } else if (r < -half_p) {
        r = -p - r;
    }
    return r;
}

PreciseDecimal reduce_cos_argument(const PreciseDecimal& x, bool* negate) {
    PreciseDecimal r = reduce_trig_argument(x);
    const PreciseDecimal p = precise::pi();
    const PreciseDecimal half_p = precise::half_pi();
    *negate = false;
    if (r > half_p) {
        r = p - r;
        *negate = true;
    } else if (r < -half_p) {
        r = -p - r;
        *negate = true;
    }
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
    if (value == LLONG_MIN) {
        // Special case: LLONG_MIN = -9223372036854775808
        // -LLONG_MIN overflows, so handle it directly
        data = string_to_bigint("9223372036854775808");
    } else {
        data = string_to_bigint(std::to_string(negative ? -value : value));
    }
    scale = 0;
}

PreciseDecimal::PreciseDecimal(long double value) {
    if (mymath::isnan(value)) {
        *this = PreciseDecimal::nan();
        return;
    }
    if (mymath::isinf(value)) {
        *this = value < 0.0L ? PreciseDecimal::neg_infinity() : PreciseDecimal::infinity();
        return;
    }
    // Use higher precision for small values to avoid truncation to zero
    int precision = 40;
    if (mymath::abs(value) > 0.0L && mymath::abs(value) < 1e-10L) {
        // For very small values, use scientific notation with full precision
        std::ostringstream out;
        out << std::scientific << std::setprecision(35) << value;
        *this = from_decimal_literal(out.str());
    } else {
        *this = from_decimal_literal(format_decimal(value, precision));
    }
}

PreciseDecimal::PreciseDecimal(const std::string& token) {
    *this = from_decimal_literal(token);
}

void PreciseDecimal::normalize() {
    if (is_inf || is_nan || g_suppress_normalization) return;

    if (data.empty() || (data.size() == 1 && data[0] == 0)) {
        data = {0};
        scale = 0;
        negative = false;
        return;
    }

    // 去除末尾的零（缩减 scale）
    if (scale > 0) {
        int zeros = count_trailing_zeros(data);
        int to_remove = std::min(zeros, scale);
        if (to_remove > 0) {
            data = divide_bigint_by_pow10(std::move(data), to_remove);
            scale -= to_remove;
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

/**
 * @brief 原地加法优化
 * 当 this 的 capacity 足够时，直接在原内存上操作，避免重新分配
 */
PreciseDecimal& PreciseDecimal::operator+=(const PreciseDecimal& rhs) {
    // Handle infinity
    if (is_inf && rhs.is_inf) {
        if (negative != rhs.negative) {
            *this = PreciseDecimal::nan();
        }
        return *this;
    }
    if (is_inf) return *this;
    if (rhs.is_inf) {
        *this = rhs;
        return *this;
    }

    // Handle NaN
    if (is_nan || rhs.is_nan) {
        *this = PreciseDecimal::nan();
        return *this;
    }

    if (negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return *this -= rhs_flipped;
    }

    // 对齐 scale
    PreciseDecimal r = rhs;
    if (scale < r.scale) {
        data = multiply_bigint_by_power_of_10(std::move(data), r.scale - scale);
        scale = r.scale;
    } else if (scale > r.scale) {
        r.data = multiply_bigint_by_power_of_10(std::move(r.data), scale - r.scale);
    }

    // 原地加法：复用 this 的内存
    size_t max_len = std::max(data.size(), r.data.size()) + 1;
    if (data.capacity() >= max_len) {
        // 容量足够时也要先扩展逻辑长度，否则后续 resize 会把写入到容量区的高位清零。
        const size_t old_size = data.size();
        data.resize(max_len, 0);
        size_t new_size = add_inplace_raw(data.data(), old_size, r.data.data(), r.data.size(), data.size());
        data.resize(new_size);
    } else {
        // 需要扩容
        BigIntData new_data(max_len);
        size_t new_size = add_raw(data.data(), data.size(), r.data.data(), r.data.size(), new_data.data());
        new_data.resize(new_size);
        data = std::move(new_data);
    }

    normalize();
    return *this;
}

/**
 * @brief 原地减法优化
 */
PreciseDecimal& PreciseDecimal::operator-=(const PreciseDecimal& rhs) {
    // Handle infinity
    if (is_inf && rhs.is_inf) {
        if (negative == rhs.negative) {
            *this = PreciseDecimal::nan();
        }
        return *this;
    }
    if (is_inf) return *this;
    if (rhs.is_inf) {
        *this = rhs;
        negative = !negative;
        return *this;
    }

    // Handle NaN
    if (is_nan || rhs.is_nan) {
        *this = PreciseDecimal::nan();
        return *this;
    }

    if (negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return *this += rhs_flipped;
    }

    // 对齐 scale
    PreciseDecimal r = rhs;
    if (scale < r.scale) {
        data = multiply_bigint_by_power_of_10(std::move(data), r.scale - scale);
        scale = r.scale;
    } else if (scale > r.scale) {
        r.data = multiply_bigint_by_power_of_10(std::move(r.data), scale - r.scale);
    }

    // 比较大小
    int cmp = compare_raw(data.data(), data.size(), r.data.data(), r.data.size());
    if (cmp == 0) {
        data = {0};
        scale = 0;
        negative = false;
        return *this;
    }

    if (cmp > 0) {
        // this > rhs，原地减法
        size_t new_size = sub_inplace_raw(data.data(), data.size(), r.data.data(), r.data.size());
        data.resize(new_size);
    } else {
        // this < rhs，需要交换
        BigIntData new_data(r.data.size());
        size_t new_size = sub_raw(r.data.data(), r.data.size(), data.data(), data.size(), new_data.data());
        new_data.resize(new_size);
        data = std::move(new_data);
        negative = !negative;
    }

    normalize();
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
    if (is_nan || rhs.is_nan) return false;
    if (is_inf || rhs.is_inf) {
        return is_inf && rhs.is_inf && negative == rhs.negative;
    }
    return compare_precise_decimal(*this, rhs) == 0;
}

bool PreciseDecimal::operator<(const PreciseDecimal& rhs) const {
    if (is_nan || rhs.is_nan) return false;
    if (is_inf || rhs.is_inf) {
        if (is_inf && rhs.is_inf) return negative && !rhs.negative;
        if (is_inf) return negative;
        return !rhs.negative;
    }
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
        res *= long_double_power_of_10(static_cast<int>(exponent));
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
    PreciseDecimal x = sqrt_initial_guess(val);
    const PreciseDecimal one_half = half();
    const PreciseDecimal epsilon = scale_epsilon();

    NormalizationSuppressor suppressor; // Suppress normalization in the loop

    const int iterations = std::max(12, PrecisionContext::get_default_scale() / 8 + 8);
    for (int i = 0; i < iterations; ++i) {
        PreciseDecimal next = one_half * (x + val / x);
        if (precise::abs(next - x) <= epsilon * std::max(one(), precise::abs(next))) {
            g_suppress_normalization = false; // Restore to normalize the result
            next.normalize();
            return next;
        }
        x = next;
    }
    g_suppress_normalization = false;
    x.normalize();
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

    long long integer_exp = 0;
    if (try_integral_exponent(exp, &integer_exp)) {
        if (integer_exp < 0) return one() / pow(base, -integer_exp);
        return pow(base, integer_exp);
    }

    long long thirds_exp = 0;
    const bool is_thirds = try_thirds_exponent(exp, &thirds_exp);

    if (base.negative) {
        if (is_thirds && thirds_exp == 1) {
            return cbrt_precise_decimal(base);
        }
        if (is_thirds && thirds_exp == 2) {
            PreciseDecimal root = cbrt_precise_decimal(base);
            return root * root;
        }
        if (is_thirds && thirds_exp == -1) {
            return one() / cbrt_precise_decimal(base);
        }
        throw std::domain_error("negative base with non-integer exponent");
    }

    if (exp == half()) return precise::sqrt(base);
    if (is_thirds && thirds_exp == 1) {
        return cbrt_precise_decimal(base);
    }
    if (is_thirds && thirds_exp == 2) {
        PreciseDecimal root = cbrt_precise_decimal(base);
        return root * root;
    }
    if (is_thirds && thirds_exp == -1) {
        return one() / cbrt_precise_decimal(base);
    }

    // 使用 exp(ln(base) * exp) 计算
    return precise::exp(precise::ln(base) * exp);
}

PreciseDecimal floor(const PreciseDecimal& val) {
    if (val.scale <= 0) return val;

    PreciseDecimal res;
    res.data = divide_bigint_by_pow10(val.data, val.scale);
    res.scale = 0;
    res.negative = val.negative;
    res.normalize();

    if (val.negative && !is_bigint_multiple_of_pow10(val.data, val.scale)) {
        res -= PreciseDecimal(1LL);
    }
    return res;
}

PreciseDecimal ceil(const PreciseDecimal& val) {
    if (val.scale <= 0) return val;
    if (is_bigint_multiple_of_pow10(val.data, val.scale)) {
        PreciseDecimal res = val;
        res.normalize();
        return res;
    }
    PreciseDecimal f = precise::floor(val);
    return f + PreciseDecimal(1LL);
}

PreciseDecimal round(const PreciseDecimal& val) {
    if (val.is_zero()) return val;
    PreciseDecimal half_val = half();
    if (val.negative) return precise::ceil(val - half_val);
    return precise::floor(val + half_val);
}

PreciseDecimal pi() {
    static const PreciseDecimal p(
        "3.141592653589793238462643383279502884197169399375105820974944592307816406286"
        "208998628034825342117067982148086513282306647093844609550582231725359408128481"
        "117450284102701938521105559644622948954930381964428810975665933446128475648233");
    return p;
}

PreciseDecimal two_pi() {
    static const PreciseDecimal tp(
        "6.283185307179586476925286766559005768394338798750211641949889184615632812572"
        "417997256069650684234135964296173026564613294187689219101164434507187816256962"
        "234900568205403877042211119289245897909860763928857621951331866892256951296466");
    return tp;
}

PreciseDecimal half_pi() {
    static const PreciseDecimal hp(
        "1.5707963267948966192313216916397514420985846996875529104874722961539082031431"
        "0449931401741267105853399107404325664115332354692223077529111586285977040642405"
        "587251420513509692605527798223114744974651909822144054878329667323064237824116");
    return hp;
}

PreciseDecimal e() {
    static const PreciseDecimal val_e(
        "2.718281828459045235360287471352662497757247093699959574966967627724076630353"
        "547594571382178525166427427466391932003059921817413596629043572900334295260595"
        "630738132328627943490763233829880753195251019011573834187930702154089149934884");
    return val_e;
}

// ============================================================================
// 优化的数学函数：exp2, log2
// ============================================================================

/**
 * @brief 预计算的 ln(2) 常量
 */
PreciseDecimal ln2() {
    static const PreciseDecimal val_ln2(
        "0.693147180559945309417232121458176568075500134360255254120680009493393621969"
        "694715605863326996418687542001481020570685733685520235758130557032670751635075"
        "96519478973694723683769222563258791851371881657850144027689067671439205263563");
    return val_ln2;
}

/**
 * @brief 预计算的 log2(e) 常量
 */
PreciseDecimal log2e() {
    static const PreciseDecimal val_log2e(
        "1.442695040888963407359924681001892137426645954152985934135449406931109219181"
        "118578976889241724881935745812784855883514987960813830091768091530425806977874"
        "64647182847503653289639175625393204426815406783430527636057096275114863783363");
    return val_log2e;
}

/**
 * @brief 计算 2^x（优化的指数函数）
 *
 * 使用分解方法：2^x = 2^int_part * 2^frac_part
 * - 2^int_part 通过位移快速计算
 * - 2^frac_part 使用泰勒级数
 *
 * @param x 指数
 * @return 2^x
 */
PreciseDecimal exp2(const PreciseDecimal& x) {
    if (x.is_zero()) return one();

    ScopedPrecision guard(8);

    // 分离整数部分和小数部分
    PreciseDecimal int_part = floor(x);
    PreciseDecimal frac_part = x - int_part;

    // 处理负指数
    if (int_part.negative) {
        return one() / exp2(-x);
    }

    // 计算 2^int_part（通过位移）
    long long n = static_cast<long long>(int_part.to_double());
    PreciseDecimal int_result = one();
    if (n > 0) {
        // 2^n 可以通过重复平方快速计算
        BigIntData power_of_2 = {1};
        for (long long i = 0; i < n; ++i) {
            power_of_2 = multiply_bigint_by_uint32(power_of_2, 2);
        }
        int_result.data = power_of_2;
        int_result.scale = 0;
        int_result.normalize();
    }

    // 计算 2^frac_part（使用泰勒级数）
    // 2^y = e^(y * ln(2)) = sum((y * ln(2))^k / k!)
    PreciseDecimal y = frac_part;
    PreciseDecimal ln2_val = ln2();

    // 使用泰勒级数展开
    NormalizationSuppressor suppressor;

    PreciseDecimal term = one();
    PreciseDecimal sum = one();
    PreciseDecimal y_ln2 = y * ln2_val;
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(150);

    for (int i = 1; i < limit; ++i) {
        term = term * y_ln2 / decimal_from_uint(static_cast<uint32_t>(i));
        sum += term;
        if (precise::abs(term) < epsilon) break;
    }

    g_suppress_normalization = false;
    sum.normalize();

    return int_result * sum;
}

/**
 * @brief 计算 log2(x)（优化的对数函数）
 *
 * log2(x) = ln(x) / ln(2)
 *
 * @param x 输入值（必须为正数）
 * @return log2(x)
 */
PreciseDecimal log2(const PreciseDecimal& x) {
    if (x <= PreciseDecimal(0LL)) throw std::domain_error("log2 of non-positive number");
    if (x == one()) return PreciseDecimal(0LL);
    if (x == two()) return one();

    // log2(x) = ln(x) / ln(2)
    return ln(x) / ln2();
}

/**
 * @brief 快速幂运算（优化版）
 *
 * 对于整数指数，使用二分幂算法。
 * 对于分数指数，使用 exp(ln(base) * exp)。
 *
 * @param base 底数
 * @param exp 指数
 * @return base^exp
 */
PreciseDecimal pow_fast(const PreciseDecimal& base, const PreciseDecimal& exponent) {
    if (exponent.is_zero()) return one();
    if (base.is_zero()) {
        if (exponent.negative) throw std::domain_error("0^negative is undefined");
        return PreciseDecimal(0LL);
    }
    if (base == one()) return one();
    if (base == two()) return exp2(exponent);

    // 检查是否为整数指数
    long long int_exp = 0;
    if (try_integral_exponent(exponent, &int_exp)) {
        if (int_exp < 0) return one() / pow_fast(base, -int_exp);
        return pow(base, int_exp);
    }

    // 检查是否为 1/2（平方根）
    if (exponent == half()) return precise::sqrt(base);

    // 检查是否为 1/3（立方根）
    long long thirds_exp = 0;
    if (try_thirds_exponent(exponent, &thirds_exp)) {
        if (thirds_exp == 1) return cbrt_precise_decimal(base);
        if (thirds_exp == 2) {
            PreciseDecimal root = cbrt_precise_decimal(base);
            return root * root;
        }
        if (thirds_exp == -1) return one() / cbrt_precise_decimal(base);
    }

    // 一般情况：使用 exp(ln(base) * exponent)
    if (base.negative) {
        throw std::domain_error("negative base with non-integer exponent");
    }

    return precise::exp(ln(base) * exponent);
}

PreciseDecimal exp(const PreciseDecimal& x) {
    if (x.is_zero()) return one();

    const long double x_ld = static_cast<long double>(x);
    if (x_ld >= mymath::kLnDoubleMax) return PreciseDecimal::infinity();
    if (x_ld <= mymath::kLnDoubleDenormMin) return PreciseDecimal(0LL);

    ScopedPrecision guard(8); // Add 8 guard digits

    if (x.negative) {
        PreciseDecimal abs_x = precise::abs(x);
        return one() / exp(abs_x);
    }

    int k = 0;
    PreciseDecimal r = x;
    const PreciseDecimal threshold("0.01");
    while (r > threshold) {
        r /= two();
        k++;
    }

    NormalizationSuppressor suppressor;

    PreciseDecimal term = one();
    PreciseDecimal sum = one();
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(150); // Increased iteration limit to support high precision
    for (int i = 1; i < limit; ++i) {
        term = term * r / decimal_from_uint(static_cast<uint32_t>(i));
        sum += term;
        if (precise::abs(term) < epsilon) break;
    }
    for (int i = 0; i < k; ++i) sum = sum * sum;
    
    g_suppress_normalization = false;
    sum.normalize();
    return sum;
}

PreciseDecimal ln(const PreciseDecimal& x) {
    if (x <= PreciseDecimal(0LL)) throw std::domain_error("ln of non-positive number");
    if (x == one()) return PreciseDecimal(0LL);

    ScopedPrecision guard(8); // Add 8 guard digits

    PreciseDecimal reduced = x;
    int sqrt_count = 0;
    const PreciseDecimal lower("0.75");
    const PreciseDecimal upper("1.5");
    while ((reduced < lower || reduced > upper) && sqrt_count < 32) {
        reduced = precise::sqrt(reduced);
        ++sqrt_count;
    }

    PreciseDecimal result = ln_near_one(reduced);
    for (int i = 0; i < sqrt_count; ++i) result *= two();
    return result;
}

PreciseDecimal log10(const PreciseDecimal& x) {
    static const PreciseDecimal ln10 = ln(PreciseDecimal(10LL));
    return ln(x) / ln10;
}

PreciseDecimal sin(const PreciseDecimal& x) {
    if (x.is_zero()) return PreciseDecimal(0LL);
    
    ScopedPrecision guard(8);
    bool negate = false;
    PreciseDecimal r = reduce_sin_argument(x, &negate);

    NormalizationSuppressor suppressor;

    PreciseDecimal term = r;
    PreciseDecimal sum = r;
    PreciseDecimal r2 = r * r;
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(150);
    for (int i = 1; i < limit; ++i) {
        term = -term * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i) * (2 * i + 1)));
        sum += term;
        if (precise::abs(term) < epsilon) break;
    }
    
    g_suppress_normalization = false;
    sum.normalize();
    return negate ? -sum : sum;
}

PreciseDecimal cos(const PreciseDecimal& x) {
    if (x.is_zero()) return one();
    
    ScopedPrecision guard(8);
    bool negate = false;
    PreciseDecimal r = reduce_cos_argument(x, &negate);

    NormalizationSuppressor suppressor;

    PreciseDecimal term = one();
    PreciseDecimal sum = one();
    PreciseDecimal r2 = r * r;
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(150);
    for (int i = 1; i < limit; ++i) {
        term = -term * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i - 1) * (2 * i)));
        sum += term;
        if (precise::abs(term) < epsilon) break;
    }
    
    g_suppress_normalization = false;
    sum.normalize();
    return negate ? -sum : sum;
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
    if (x.is_zero()) return PreciseDecimal(0LL);
    
    ScopedPrecision guard(8); // Add guard digits
    
    if (x == one()) return pi() / PreciseDecimal(4LL);
    if (x == -one()) return -pi() / PreciseDecimal(4LL);
    if (precise::abs(x) > one()) {
        if (x.negative) return -pi() / two() - atan(one() / x);
        return pi() / two() - atan(one() / x);
    }

    NormalizationSuppressor suppressor;

    PreciseDecimal term = x;
    PreciseDecimal sum = x;
    const PreciseDecimal x2 = x * x;
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(200);
    for (int i = 1; i < limit; ++i) {
        term = -term * x2;
        PreciseDecimal add = term / decimal_from_uint(static_cast<uint32_t>(2 * i + 1));
        sum += add;
        if (precise::abs(add) < epsilon) break;
    }
    
    g_suppress_normalization = false;
    sum.normalize();
    return sum;
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
    // Handle infinity
    if (lhs.is_inf && rhs.is_inf) {
        if (lhs.negative != rhs.negative) {
            return PreciseDecimal::nan(); // inf + (-inf) = nan
        }
        return lhs; // inf + inf = inf (with same sign)
    }
    if (lhs.is_inf) return lhs;
    if (rhs.is_inf) return rhs;

    // Handle NaN
    if (lhs.is_nan || rhs.is_nan) return PreciseDecimal::nan();

    if (lhs.negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return subtract_precise_decimal(lhs, rhs_flipped);
    }

    // 优化：尝试原地修改，减少内存分配
    PreciseDecimal l = lhs;
    PreciseDecimal r = rhs;
    align_precise_scales(&l, &r);

    PreciseDecimal res;
    res.scale = l.scale;
    res.negative = lhs.negative;

    // 使用原始指针 API 进行加法
    size_t max_len = std::max(l.data.size(), r.data.size()) + 1;
    res.data.resize(max_len);
    size_t new_size = add_raw(l.data.data(), l.data.size(), r.data.data(), r.data.size(), res.data.data());
    res.data.resize(new_size);
    res.normalize();
    return res;
}

PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    // Handle infinity
    if (lhs.is_inf && rhs.is_inf) {
        if (lhs.negative == rhs.negative) {
            return PreciseDecimal::nan(); // inf - inf = nan
        }
        return lhs; // inf - (-inf) = inf
    }
    if (lhs.is_inf) return lhs;
    if (rhs.is_inf) {
        PreciseDecimal res = rhs;
        res.negative = !res.negative; // -(-inf) = inf
        return res;
    }

    // Handle NaN
    if (lhs.is_nan || rhs.is_nan) return PreciseDecimal::nan();

    if (lhs.negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return add_precise_decimal(lhs, rhs_flipped);
    }

    PreciseDecimal l = lhs;
    PreciseDecimal r = rhs;
    align_precise_scales(&l, &r);

    // 使用原始指针比较
    int cmp = compare_raw(l.data.data(), l.data.size(), r.data.data(), r.data.size());
    if (cmp == 0) return {};

    PreciseDecimal res;
    res.scale = l.scale;
    res.data.resize(std::max(l.data.size(), r.data.size()));

    if (cmp > 0) {
        size_t new_size = sub_raw(l.data.data(), l.data.size(), r.data.data(), r.data.size(), res.data.data());
        res.data.resize(new_size);
        res.negative = lhs.negative;
    } else {
        size_t new_size = sub_raw(r.data.data(), r.data.size(), l.data.data(), l.data.size(), res.data.data());
        res.data.resize(new_size);
        res.negative = !lhs.negative;
    }
    res.normalize();
    return res;
}

PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    // Handle infinity
    if (lhs.is_inf || rhs.is_inf) {
        if (lhs.is_zero() || rhs.is_zero()) {
            return PreciseDecimal::nan(); // 0 * inf = nan
        }
        PreciseDecimal res = lhs.is_inf ? lhs : rhs;
        res.negative = lhs.negative != rhs.negative;
        return res;
    }

    // Handle NaN
    if (lhs.is_nan || rhs.is_nan) return PreciseDecimal::nan();

    PreciseDecimal res;
    res.data = multiply_bigint(lhs.data, rhs.data);
    res.scale = lhs.scale + rhs.scale;
    res.negative = lhs.negative != rhs.negative;
    res.normalize();
    return res;
}

PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    // Handle infinity
    if (lhs.is_inf && rhs.is_inf) {
        return PreciseDecimal::nan(); // inf / inf = nan
    }
    if (lhs.is_inf) {
        PreciseDecimal res = lhs;
        res.negative = lhs.negative != rhs.negative;
        return res;
    }
    if (rhs.is_inf) {
        return PreciseDecimal(0LL); // finite / inf = 0
    }

    // Handle NaN
    if (lhs.is_nan || rhs.is_nan) return PreciseDecimal::nan();

    if (rhs.is_zero()) throw std::runtime_error("division by zero");
    if (lhs.is_zero()) return {};

    int target_scale = PrecisionContext::get_default_scale();
    int numerator_shift = target_scale + rhs.scale - lhs.scale;
    BigIntData numerator = lhs.data;
    if (numerator_shift >= 0) {
        numerator = multiply_bigint_by_power_of_10(numerator, numerator_shift);
    } else {
        BigIntData divisor = multiply_bigint_by_power_of_10({1}, -numerator_shift);
        BigIntData truncated, ignored_remainder;
        div_bigint(numerator, divisor, &truncated, &ignored_remainder);
        numerator = truncated;
    }
    BigIntData denominator = rhs.data;

    BigIntData q, r;
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

    int res = 0;
    if (lhs.scale == rhs.scale) {
        res = compare_bigint(lhs.data, rhs.data);
    } else if (lhs.scale < rhs.scale) {
        BigIntData shifted_lhs = multiply_bigint_by_power_of_10(lhs.data, rhs.scale - lhs.scale);
        res = compare_bigint(shifted_lhs, rhs.data);
    } else {
        BigIntData shifted_rhs = multiply_bigint_by_power_of_10(rhs.data, lhs.scale - rhs.scale);
        res = compare_bigint(lhs.data, shifted_rhs);
    }

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
        BigIntData result = {1};
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

// ============================================================================
// BinaryBigInt 成员函数实现
// ============================================================================

// 二进制大整数加法（原始指针版本）
inline size_t binary_add_raw(const uint32_t* a, size_t an, const uint32_t* b, size_t bn, uint32_t* res) {
    uint64_t carry = 0;
    size_t i = 0;
    size_t min_n = std::min(an, bn);
    for (; i < min_n; ++i) {
        uint64_t sum = static_cast<uint64_t>(a[i]) + b[i] + carry;
        res[i] = static_cast<uint32_t>(sum);
        carry = sum >> 32;
    }
    const uint32_t* longer = (an > bn) ? a : b;
    size_t max_n = std::max(an, bn);
    for (; i < max_n; ++i) {
        uint64_t sum = static_cast<uint64_t>(longer[i]) + carry;
        res[i] = static_cast<uint32_t>(sum);
        carry = sum >> 32;
    }
    if (carry) {
        res[i++] = static_cast<uint32_t>(carry);
    }
    return i;
}

// 二进制大整数减法（原始指针版本，要求 a >= b）
inline size_t binary_sub_raw(const uint32_t* a, size_t an, const uint32_t* b, size_t bn, uint32_t* res) {
    uint64_t borrow = 0;
    size_t i = 0;
    for (; i < bn; ++i) {
        uint64_t diff = static_cast<uint64_t>(a[i]) - b[i] - borrow;
        res[i] = static_cast<uint32_t>(diff);
        borrow = (diff >> 63) & 1;
    }
    for (; i < an; ++i) {
        uint64_t diff = static_cast<uint64_t>(a[i]) - borrow;
        res[i] = static_cast<uint32_t>(diff);
        borrow = (diff >> 63) & 1;
    }
    while (i > 1 && res[i - 1] == 0) --i;
    return i;
}

// 二进制大整数比较
inline int binary_compare_raw(const uint32_t* a, size_t an, const uint32_t* b, size_t bn) {
    if (an != bn) return an < bn ? -1 : 1;
    for (size_t i = an; i > 0; --i) {
        if (a[i - 1] != b[i - 1]) {
            return a[i - 1] < b[i - 1] ? -1 : 1;
        }
    }
    return 0;
}

// ============================================================================
// 进制转换优化：预计算幂次
// ============================================================================

/**
 * @brief 预计算的 10^9 幂次表（用于十进制转二进制）
 *
 * powers_of_base[i] = (10^9)^i 的二进制表示
 * 使用快速幂计算，避免重复计算。
 */
class PowersOfBaseCache {
public:
    static constexpr std::size_t MAX_POWER = 64;  // 支持最大 64 * 9 = 576 位十进制数

    // 10^9 的二进制表示（预计算）
    // 10^9 = 1000000000 = 0x3B9ACA00
    static const BigIntData& get(std::size_t i) {
        if (i < cache_.size() && !cache_[i].empty()) {
            return cache_[i];
        }
        ensure_cached(i);
        return cache_[i];
    }

    static void ensure_cached(std::size_t max_index) {
        if (max_index < cache_.size()) return;

        std::lock_guard<std::mutex> lock(mutex_);
        if (max_index < cache_.size()) return;

        // 扩展缓存
        while (cache_.size() <= max_index) {
            std::size_t idx = cache_.size();
            if (idx == 0) {
                cache_.push_back({1});
            } else if (idx == 1) {
                // 10^9 = 0x3B9ACA00
                cache_.push_back({0x3B9ACA00u});
            } else {
                // 使用快速幂：powers[idx] = powers[idx/2] * powers[idx/2] * (idx 奇数 ? base : 1)
                BigIntData result = multiply_bigint(cache_[idx / 2], cache_[idx / 2]);
                if (idx % 2 == 1) {
                    result = multiply_bigint(result, cache_[1]);
                }
                cache_.push_back(std::move(result));
            }
        }
    }

private:
    static std::vector<BigIntData> cache_;
    static std::mutex mutex_;
};

std::vector<BigIntData> PowersOfBaseCache::cache_;
std::mutex PowersOfBaseCache::mutex_;

/**
 * @brief 预计算的 2^32 幂次表（用于二进制转十进制）
 */
class PowersOf2_32Cache {
public:
    static constexpr std::size_t MAX_POWER = 64;

    static const BigIntData& get(std::size_t i) {
        if (i < cache_.size() && !cache_[i].empty()) {
            return cache_[i];
        }
        ensure_cached(i);
        return cache_[i];
    }

    static void ensure_cached(std::size_t max_index) {
        if (max_index < cache_.size()) return;

        std::lock_guard<std::mutex> lock(mutex_);
        if (max_index < cache_.size()) return;

        while (cache_.size() <= max_index) {
            std::size_t idx = cache_.size();
            if (idx == 0) {
                cache_.push_back({1});
            } else {
                // 2^32 = [0, 1]（两个 chunk）
                // powers[idx] = powers[idx-1] * 2^32
                BigIntData prev = cache_[idx - 1];
                BigIntData result;
                result.reserve(prev.size() + 1);
                result.push_back(0);  // 左移一个 chunk
                for (const auto& val : prev) {
                    result.push_back(val);
                }
                cache_.push_back(std::move(result));
            }
        }
    }

private:
    static std::vector<BigIntData> cache_;
    static std::mutex mutex_;
};

std::vector<BigIntData> PowersOf2_32Cache::cache_;
std::mutex PowersOf2_32Cache::mutex_;

// 十进制基数转二进制基数（优化版：使用预计算幂次）
BigIntData decimal_to_binary(const BigIntData& decimal) {
    if (decimal.empty() || (decimal.size() == 1 && decimal[0] == 0)) {
        return {0};
    }

    // 预计算所需的幂次
    PowersOfBaseCache::ensure_cached(decimal.size());

    // 使用快速累加：value = sum(decimal[i] * (10^9)^i)
    // 采用分治策略累加
    if (decimal.size() <= 8) {
        // 小规模：直接累加
        BigIntData result = {decimal[0]};
        for (size_t i = 1; i < decimal.size(); ++i) {
            BigIntData term = multiply_bigint_by_uint32(PowersOfBaseCache::get(i), decimal[i]);
            result = add_bigint(result, term);
        }
        return result;
    }

    // 大规模：分治累加
    // 将 decimal 分成两半，分别计算后合并
    std::size_t mid = decimal.size() / 2;

    BigIntData low(decimal.begin(), decimal.begin() + static_cast<std::ptrdiff_t>(mid));
    BigIntData high(decimal.begin() + static_cast<std::ptrdiff_t>(mid), decimal.end());

    BigIntData low_result = decimal_to_binary(low);
    BigIntData high_result = decimal_to_binary(high);

    // high_result 需要乘以 (10^9)^mid
    BigIntData power_mid = PowersOfBaseCache::get(mid);
    high_result = multiply_bigint(high_result, power_mid);

    return add_bigint(low_result, high_result);
}

// 二进制基数转十进制基数（优化版：使用预计算幂次）
BigIntData binary_to_decimal(const BigIntData& binary) {
    if (binary.empty() || (binary.size() == 1 && binary[0] == 0)) {
        return {0};
    }

    // 预计算所需的幂次
    PowersOf2_32Cache::ensure_cached(binary.size());

    // 使用分治累加
    if (binary.size() <= 8) {
        // 小规模：直接累加
        BigIntData result = {binary[0]};
        for (size_t i = 1; i < binary.size(); ++i) {
            BigIntData term = multiply_bigint_by_uint32(PowersOf2_32Cache::get(i), binary[i]);
            result = add_bigint(result, term);
        }
        return result;
    }

    // 大规模：分治累加
    std::size_t mid = binary.size() / 2;

    BigIntData low(binary.begin(), binary.begin() + static_cast<std::ptrdiff_t>(mid));
    BigIntData high(binary.begin() + static_cast<std::ptrdiff_t>(mid), binary.end());

    BigIntData low_result = binary_to_decimal(low);
    BigIntData high_result = binary_to_decimal(high);

    // high_result 需要乘以 (2^32)^mid
    BigIntData power_mid = PowersOf2_32Cache::get(mid);
    high_result = multiply_bigint(high_result, power_mid);

    return add_bigint(low_result, high_result);
}

// 二进制大整数乘法（朴素算法，利用位移优化）
BigIntData binary_mul_naive(const BigIntData& a, const BigIntData& b) {
    if ((a.size() == 1 && a[0] == 0) || (b.size() == 1 && b[0] == 0)) {
        return {0};
    }

    BigIntData res(a.size() + b.size(), 0);

    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == 0) continue;

        uint64_t carry = 0;
        uint64_t a_val = a[i];
        for (size_t j = 0; j < b.size(); ++j) {
            // 关键优化：使用位移代替除法/模运算
            uint64_t product = a_val * b[j];
            uint64_t sum = static_cast<uint64_t>(res[i + j]) + (product & 0xFFFFFFFFULL) + carry;
            res[i + j] = static_cast<uint32_t>(sum);
            carry = (sum >> 32) + (product >> 32);
        }
        if (carry) {
            size_t k = i + b.size();
            while (carry && k < res.size()) {
                uint64_t sum = static_cast<uint64_t>(res[k]) + carry;
                res[k] = static_cast<uint32_t>(sum);
                carry = sum >> 32;
                ++k;
            }
            if (carry) {
                res.push_back(static_cast<uint32_t>(carry));
            }
        }
    }

    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

// 二进制大整数乘以小整数
BigIntData binary_mul_uint32(const BigIntData& a, uint32_t b) {
    if (b == 0 || (a.size() == 1 && a[0] == 0)) {
        return {0};
    }
    if (b == 1) return a;

    BigIntData res;
    res.reserve(a.size() + 1);
    uint64_t carry = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        uint64_t product = static_cast<uint64_t>(a[i]) * b + carry;
        res.push_back(static_cast<uint32_t>(product));
        carry = product >> 32;
    }
    if (carry) res.push_back(static_cast<uint32_t>(carry));
    return res;
}

// 二进制大整数除法（优化版本）
void binary_div(const BigIntData& num, const BigIntData& den,
                BigIntData* quotient, BigIntData* remainder) {
    // 检查除数是否为 2 的幂次
    bool is_power_of_two = false;
    int shift_amount = 0;
    if (den.size() == 1) {
        uint32_t d = den[0];
        if (d > 0 && (d & (d - 1)) == 0) {
            is_power_of_two = true;
            while (d > 1) {
                d >>= 1;
                ++shift_amount;
            }
        }
    } else if (den.size() == 2 && den[1] == 1 && den[0] == 0) {
        // 2^32
        is_power_of_two = true;
        shift_amount = 32;
    }

    if (is_power_of_two) {
        // 使用位移优化
        if (shift_amount == 0) {
            *quotient = num;
            *remainder = {0};
            return;
        }

        if (shift_amount < 32) {
            // 除以 2^k，k < 32
            *quotient = BigIntData(num.size(), 0);
            *remainder = {num[0] & ((1U << shift_amount) - 1)};
            for (size_t i = 0; i < num.size(); ++i) {
                (*quotient)[i] = num[i] >> shift_amount;
                if (i + 1 < num.size()) {
                    (*quotient)[i] |= num[i + 1] << (32 - shift_amount);
                }
            }
        } else {
            // 除以 2^32 或更大
            size_t limb_shift = shift_amount / 32;
            int bit_shift = shift_amount % 32;

            if (limb_shift >= num.size()) {
                *quotient = {0};
                *remainder = num;
            } else {
                *quotient = BigIntData(num.size() - limb_shift, 0);
                for (size_t i = 0; i < num.size() - limb_shift; ++i) {
                    (*quotient)[i] = num[i + limb_shift] >> bit_shift;
                    if (i + limb_shift + 1 < num.size() && bit_shift > 0) {
                        (*quotient)[i] |= num[i + limb_shift + 1] << (32 - bit_shift);
                    }
                }
                *remainder = BigIntData(num.begin(), num.begin() + static_cast<std::ptrdiff_t>(std::min(limb_shift, num.size())));
            }
        }
        while (quotient->size() > 1 && quotient->back() == 0) quotient->pop_back();
        while (remainder->size() > 1 && remainder->back() == 0) remainder->pop_back();
        if (quotient->empty()) quotient->push_back(0);
        if (remainder->empty()) remainder->push_back(0);
        return;
    }

    // 一般情况：使用 Knuth 算法 D
    div_bigint(num, den, quotient, remainder);
}

BinaryBigInt::BinaryBigInt(uint64_t val) {
    if (val == 0) {
        limbs = {0};
    } else {
        limbs.push_back(static_cast<uint32_t>(val & 0xFFFFFFFFULL));
        if (val > 0xFFFFFFFFULL) {
            limbs.push_back(static_cast<uint32_t>(val >> 32));
        }
    }
}

BinaryBigInt::BinaryBigInt(const BigIntData& decimal_data) {
    limbs = decimal_to_binary(decimal_data);
}

BinaryBigInt::BinaryBigInt(BigIntData&& decimal_data) {
    limbs = decimal_to_binary(decimal_data);
}

BigIntData BinaryBigInt::to_decimal() const {
    return binary_to_decimal(limbs);
}

size_t BinaryBigInt::effective_size() const {
    size_t size = limbs.size();
    while (size > 1 && limbs[size - 1] == 0) --size;
    return size;
}

void BinaryBigInt::normalize() {
    while (limbs.size() > 1 && limbs.back() == 0) {
        limbs.pop_back();
    }
    if (limbs.empty()) limbs.push_back(0);
}

BinaryBigInt BinaryBigInt::add(const BinaryBigInt& a, const BinaryBigInt& b) {
    BinaryBigInt result;
    result.limbs.resize(std::max(a.limbs.size(), b.limbs.size()) + 1);
    size_t new_size = binary_add_raw(a.limbs.data(), a.limbs.size(),
                                      b.limbs.data(), b.limbs.size(),
                                      result.limbs.data());
    result.limbs.resize(new_size);
    return result;
}

BinaryBigInt BinaryBigInt::sub(const BinaryBigInt& a, const BinaryBigInt& b) {
    BinaryBigInt result;
    result.limbs.resize(a.limbs.size());
    size_t new_size = binary_sub_raw(a.limbs.data(), a.limbs.size(),
                                      b.limbs.data(), b.limbs.size(),
                                      result.limbs.data());
    result.limbs.resize(new_size);
    return result;
}

BinaryBigInt BinaryBigInt::mul(const BinaryBigInt& a, const BinaryBigInt& b) {
    BinaryBigInt result;
    result.limbs = binary_mul_naive(a.limbs, b.limbs);
    return result;
}

BinaryBigInt BinaryBigInt::mul_uint32(const BinaryBigInt& a, uint32_t b) {
    BinaryBigInt result;
    result.limbs = binary_mul_uint32(a.limbs, b);
    return result;
}

void BinaryBigInt::div(const BinaryBigInt& num, const BinaryBigInt& den,
                       BinaryBigInt* quotient, BinaryBigInt* remainder) {
    binary_div(num.limbs, den.limbs, &quotient->limbs, &remainder->limbs);
}

int BinaryBigInt::compare(const BinaryBigInt& a, const BinaryBigInt& b) {
    return binary_compare_raw(a.limbs.data(), a.limbs.size(),
                               b.limbs.data(), b.limbs.size());
}

BinaryBigInt BinaryBigInt::operator<<(int bits) const {
    if (bits <= 0) return *this;
    if (is_zero()) return *this;

    size_t limb_shift = bits / 32;
    int bit_shift = bits % 32;

    BinaryBigInt result;
    result.limbs.resize(limbs.size() + limb_shift + 1, 0);

    if (bit_shift == 0) {
        // 只有 limb 位移
        for (size_t i = 0; i < limbs.size(); ++i) {
            result.limbs[i + limb_shift] = limbs[i];
        }
    } else {
        // 需要位位移
        uint32_t carry = 0;
        for (size_t i = 0; i < limbs.size(); ++i) {
            uint64_t val = (static_cast<uint64_t>(limbs[i]) << bit_shift) | carry;
            result.limbs[i + limb_shift] = static_cast<uint32_t>(val);
            carry = static_cast<uint32_t>(val >> 32);
        }
        if (carry) {
            result.limbs[limbs.size() + limb_shift] = carry;
        }
    }

    result.normalize();
    return result;
}

BinaryBigInt BinaryBigInt::operator>>(int bits) const {
    if (bits <= 0) return *this;
    if (is_zero()) return *this;

    size_t limb_shift = bits / 32;
    int bit_shift = bits % 32;

    if (limb_shift >= limbs.size()) {
        return BinaryBigInt(0);
    }

    BinaryBigInt result;
    result.limbs.resize(limbs.size() - limb_shift, 0);

    if (bit_shift == 0) {
        for (size_t i = 0; i < result.limbs.size(); ++i) {
            result.limbs[i] = limbs[i + limb_shift];
        }
    } else {
        for (size_t i = 0; i < result.limbs.size(); ++i) {
            result.limbs[i] = limbs[i + limb_shift] >> bit_shift;
            if (i + limb_shift + 1 < limbs.size()) {
                result.limbs[i] |= limbs[i + limb_shift + 1] << (32 - bit_shift);
            }
        }
    }

    result.normalize();
    return result;
}

#pragma GCC diagnostic pop

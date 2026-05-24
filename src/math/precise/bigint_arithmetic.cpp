// ============================================================================
// 大整数基础算术操作
// ============================================================================
//
// 本文件实现大整数的基础操作，包括：
// - 字符串与大整数转换
// - 加法、减法
// - 比较操作
// - 原始指针算术优化
// ============================================================================

#include "precise_decimal.h"
#include "app/default_precision.h"
#include "math/mymath.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>

namespace precise {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;  ///< 大整数基数（10^9）
constexpr int kBaseDigits = 9;          ///< 每个 chunk 的十进制位数

// 使用 app 命名空间中的统一精度设置
using app::get_default_scale;
using app::set_default_scale;

// 显示阈值 - 基于 default_scale 动态计算
inline long double display_zero_eps() {
    return app::display_zero_threshold();
}

inline long double display_integer_eps() {
    return app::display_integer_threshold();
}

thread_local bool g_suppress_normalization = false;

// ============================================================================
// 原始指针算术优化 (Raw Pointer Primitives)
// ============================================================================
// 使用 SIMD 优化版本，在无 AVX2 时回退到标量实现

size_t add_raw(const uint32_t* a, size_t an, const uint32_t* b, size_t bn, uint32_t* res) {
#ifdef __AVX2__
    return add_raw_optimized(a, an, b, bn, res);
#else
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
#endif
}

size_t sub_raw(const uint32_t* a, size_t an, const uint32_t* b, size_t bn, uint32_t* res) {
#ifdef __AVX2__
    return sub_raw_optimized(a, an, b, bn, res);
#else
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
#endif
}

size_t add_inplace_raw(uint32_t* a, size_t an, const uint32_t* b, size_t bn, size_t capacity) {
    uint32_t carry = 0;
    size_t i = 0;
    size_t min_n = std::min(an, bn);
    for (; i < min_n; ++i) {
        uint64_t sum = static_cast<uint64_t>(a[i]) + b[i] + carry;
        a[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
    }
    for (; i < bn && i < capacity; ++i) {
        uint64_t sum = static_cast<uint64_t>(b[i]) + carry;
        a[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
    }
    while (i < an && i < capacity) {
        uint64_t sum = static_cast<uint64_t>(a[i]) + carry;
        a[i] = static_cast<uint32_t>(sum % kBase);
        carry = static_cast<uint32_t>(sum / kBase);
        ++i;
    }
    while (carry && i < capacity) {
        a[i++] = carry % kBase;
        carry = carry / kBase;
    }
    return std::max(an, std::max(bn, i));
}

size_t sub_inplace_raw(uint32_t* a, size_t an, const uint32_t* b, size_t bn) {
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
    while (an > 1 && a[an - 1] == 0) --an;
    return an;
}

int compare_raw(const uint32_t* a, size_t an, const uint32_t* b, size_t bn) {
#ifdef __AVX2__
    return compare_raw_optimized(a, an, b, bn);
#else
    if (an != bn) return an < bn ? -1 : 1;
    for (size_t i = an; i > 0; --i) {
        if (a[i - 1] != b[i - 1]) {
            return a[i - 1] < b[i - 1] ? -1 : 1;
        }
    }
    return 0;
#endif
}

size_t mul_uint32_raw(const uint32_t* a, size_t an, uint32_t b, uint32_t* res) {
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

size_t mul_uint32_inplace_raw(uint32_t* a, size_t an, uint32_t b, size_t capacity) {
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

// ============================================================================
// 大整数基础操作函数
// ============================================================================

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

int compare_bigint(const BigIntData& lhs, const BigIntData& rhs) {
#ifdef __AVX2__
    return compare_raw_optimized(lhs.data(), lhs.size(), rhs.data(), rhs.size());
#else
    if (lhs.size() != rhs.size()) return lhs.size() < rhs.size() ? -1 : 1;
    for (int i = static_cast<int>(lhs.size()) - 1; i >= 0; --i) {
        if (lhs[static_cast<std::size_t>(i)] != rhs[static_cast<std::size_t>(i)]) {
            return lhs[static_cast<std::size_t>(i)] < rhs[static_cast<std::size_t>(i)] ? -1 : 1;
        }
    }
    return 0;
#endif
}

BigIntData add_bigint(const BigIntData& lhs, const BigIntData& rhs) {
    BigIntData res;
    res.resize(std::max(lhs.size(), rhs.size()) + 1);
    size_t new_size = add_raw(lhs.data(), lhs.size(), rhs.data(), rhs.size(), res.data());
    res.resize(new_size);
    return res;
}

BigIntData subtract_bigint(const BigIntData& lhs, const BigIntData& rhs) {
    BigIntData res;
    res.resize(lhs.size());
    size_t new_size = sub_raw(lhs.data(), lhs.size(), rhs.data(), rhs.size(), res.data());
    res.resize(new_size);
    return res;
}

// ============================================================================
// 辅助函数
// ============================================================================

long double normalize_display_decimal(long double value) {
    if (mymath::abs(value) < display_zero_eps()) return 0.0L;
    if (mymath::abs(value) > display_integer_eps() &&
        mymath::abs(value - mymath::round(value)) < display_integer_eps()) {
        return mymath::round(value);
    }
    return value;
}

std::string format_decimal(long double value, int precision) {
    value = normalize_display_decimal(value);
    char buf[128];
    int len = snprintf(buf, sizeof(buf), "%.*Lg", precision, value);
    if (len > 0 && len < static_cast<int>(sizeof(buf))) {
        return std::string(buf, len);
    }
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    return out.str();
}

// ============================================================================
// 全局状态访问
// ============================================================================

bool& get_suppress_normalization() {
    return g_suppress_normalization;
}

} // namespace precise
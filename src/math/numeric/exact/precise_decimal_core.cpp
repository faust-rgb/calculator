// ============================================================================
// PreciseDecimal 核心方法实现
// ============================================================================
//
// 本文件实现 PreciseDecimal 类的核心方法：
// - 构造函数
// - 基本运算符
// - 字符串转换
// - 比较操作
// ============================================================================

#include "precise_decimal.h"
#include "math/numeric/rational/rational.h"
#include "core/value/stored_value.h"
#include "core/common/calculator_exceptions.h"
#include "core/services/format_utils.h"
#include "math/mymath.h"
#include "math/runtime/precision/default_precision.h"

#include <algorithm>
#include <cctype>
#include <climits>
#include <iomanip>
#include <sstream>
#include <string>

namespace precise {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;
constexpr int kBaseDigits = 9;

using app::get_default_scale;
using app::set_default_scale;

inline long double display_zero_eps() {
    return app::display_zero_threshold();
}

inline long double display_integer_eps() {
    return app::display_integer_threshold();
}

long double long_double_power_of_10(int exponent) {
    return std::pow(10.0L, static_cast<long double>(exponent));
}

// ============================================================================
// NormalizationSuppressor 实现
// ============================================================================

NormalizationSuppressor::NormalizationSuppressor()
    : old_state_(g_suppress_normalization) {
    g_suppress_normalization = true;
}

NormalizationSuppressor::~NormalizationSuppressor() {
    g_suppress_normalization = old_state_;
}

ScopedNormalizationEnable::ScopedNormalizationEnable()
    : old_state_(g_suppress_normalization) {
    g_suppress_normalization = false;
}

ScopedNormalizationEnable::~ScopedNormalizationEnable() {
    g_suppress_normalization = old_state_;
}

// ============================================================================
// PreciseDecimal 构造函数
// ============================================================================

PreciseDecimal::PreciseDecimal(long long value) {
    negative = value < 0;
    if (value == LLONG_MIN) {
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
    int precision = 40;
    if (mymath::abs(value) > 0.0L && mymath::abs(value) < display_zero_eps()) {
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

// ============================================================================
// PreciseDecimal 基本方法
// ============================================================================

void PreciseDecimal::normalize() {
    if (is_inf || is_nan || g_suppress_normalization) return;

    if (data.empty() || (data.size() == 1 && data[0] == 0)) {
        data = {0};
        scale = 0;
        negative = false;
        return;
    }

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
    if (is_inf || is_nan) return false;
    return data.empty() || (data.size() == 1 && data[0] == 0);
}

bool PreciseDecimal::is_integer_value() const {
    if (is_inf || is_nan) return false;
    return scale <= 0;
}

PreciseDecimal PreciseDecimal::operator-() const {
    PreciseDecimal res = *this;
    if (res.is_inf) {
        // 无穷大翻转负号，避免被空的 data 误判为零导致符号无法翻转
        res.negative = !res.negative;
        return res;
    }
    if (!res.is_zero()) res.negative = !res.negative;
    return res;
}

// ============================================================================
// PreciseDecimal 算术运算符
// ============================================================================

PreciseDecimal& PreciseDecimal::operator+=(const PreciseDecimal& rhs) {
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

    if (is_nan || rhs.is_nan) {
        *this = PreciseDecimal::nan();
        return *this;
    }

    if (negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return *this -= rhs_flipped;
    }

    PreciseDecimal r = rhs;
    if (scale < r.scale) {
        data = multiply_bigint_by_power_of_10(std::move(data), r.scale - scale);
        scale = r.scale;
    } else if (scale > r.scale) {
        r.data = multiply_bigint_by_power_of_10(std::move(r.data), scale - r.scale);
    }

    size_t max_len = std::max(data.size(), r.data.size()) + 1;
    if (data.capacity() >= max_len) {
        const size_t old_size = data.size();
        data.resize(max_len, 0);
        size_t new_size = add_inplace_raw(data.data(), old_size, r.data.data(), r.data.size(), data.size());
        data.resize(new_size);
    } else {
        BigIntData new_data(max_len);
        size_t new_size = add_raw(data.data(), data.size(), r.data.data(), r.data.size(), new_data.data());
        new_data.resize(new_size);
        data = std::move(new_data);
    }

    normalize();
    return *this;
}

PreciseDecimal& PreciseDecimal::operator-=(const PreciseDecimal& rhs) {
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

    if (is_nan || rhs.is_nan) {
        *this = PreciseDecimal::nan();
        return *this;
    }

    if (negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return *this += rhs_flipped;
    }

    PreciseDecimal r = rhs;
    if (scale < r.scale) {
        data = multiply_bigint_by_power_of_10(std::move(data), r.scale - scale);
        scale = r.scale;
    } else if (scale > r.scale) {
        r.data = multiply_bigint_by_power_of_10(std::move(r.data), scale - r.scale);
    }

    int cmp = compare_raw(data.data(), data.size(), r.data.data(), r.data.size());
    if (cmp == 0) {
        data = {0};
        scale = 0;
        negative = false;
        return *this;
    }

    if (cmp > 0) {
        size_t new_size = sub_inplace_raw(data.data(), data.size(), r.data.data(), r.data.size());
        data.resize(new_size);
    } else {
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

// ============================================================================
// 延迟规范化接口
// ============================================================================

PreciseDecimal& PreciseDecimal::add_without_normalize(const PreciseDecimal& rhs) {
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

    if (is_nan || rhs.is_nan) {
        *this = PreciseDecimal::nan();
        return *this;
    }

    if (negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        PreciseDecimal r = rhs_flipped;
        if (scale < r.scale) {
            data = multiply_bigint_by_power_of_10(std::move(data), r.scale - scale);
            scale = r.scale;
        } else if (scale > r.scale) {
            r.data = multiply_bigint_by_power_of_10(std::move(r.data), scale - r.scale);
        }
        int cmp = compare_raw(data.data(), data.size(), r.data.data(), r.data.size());
        if (cmp == 0) {
            data = {0};
            scale = 0;
            negative = false;
            return *this;
        }
        if (cmp > 0) {
            size_t new_size = sub_inplace_raw(data.data(), data.size(), r.data.data(), r.data.size());
            data.resize(new_size);
        } else {
            BigIntData new_data(r.data.size());
            size_t new_size = sub_raw(r.data.data(), r.data.size(), data.data(), data.size(), new_data.data());
            new_data.resize(new_size);
            data = std::move(new_data);
            negative = !negative;
        }
        return *this;
    }

    PreciseDecimal r = rhs;
    if (scale < r.scale) {
        data = multiply_bigint_by_power_of_10(std::move(data), r.scale - scale);
        scale = r.scale;
    } else if (scale > r.scale) {
        r.data = multiply_bigint_by_power_of_10(std::move(r.data), scale - r.scale);
    }

    size_t max_len = std::max(data.size(), r.data.size()) + 1;
    if (data.capacity() >= max_len) {
        const size_t old_size = data.size();
        data.resize(max_len, 0);
        size_t new_size = add_inplace_raw(data.data(), old_size, r.data.data(), r.data.size(), data.size());
        data.resize(new_size);
    } else {
        BigIntData new_data(max_len);
        size_t new_size = add_raw(data.data(), data.size(), r.data.data(), r.data.size(), new_data.data());
        new_data.resize(new_size);
        data = std::move(new_data);
    }

    return *this;
}

PreciseDecimal& PreciseDecimal::mul_without_normalize(const PreciseDecimal& rhs) {
    PreciseDecimal product = multiply_precise_decimal(*this, rhs);
    NormalizationSuppressor suppressor;
    *this = product;
    return *this;
}

// ============================================================================
// 比较操作符
// ============================================================================

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

// ============================================================================
// 字符串转换
// ============================================================================

std::string PreciseDecimal::to_string(int precision) const {
    if (is_nan) return "nan";
    if (is_inf) return negative ? "-inf" : "inf";
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

    // 如果指定了精度，截断小数部分
    if (precision >= 0) {
        std::size_t dot_pos = res.find('.', res[0] == '-' ? 1 : 0);
        if (dot_pos != std::string::npos) {
            std::size_t max_len = dot_pos + 1 + static_cast<std::size_t>(precision);
            if (res.size() > max_len) {
                res.resize(max_len);
            }
            // 移除末尾的零
            while (!res.empty() && res.back() == '0') res.pop_back();
            if (!res.empty() && res.back() == '.') res.pop_back();
        }
    }

    return res;
}

long double PreciseDecimal::to_double() const {
    if (is_nan) return mymath::quiet_nan();
    if (is_inf) return negative ? -mymath::infinity() : mymath::infinity();
    if (is_zero()) return 0.0L;

    try {
        std::string s = to_string();
        if (s.empty()) return 0.0L;
        return std::stold(s);
    } catch (...) {
        return 0.0L;
    }
}

// ============================================================================
// 静态工厂方法
// ============================================================================

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
// 二元运算符
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

// ============================================================================
// DelayedNormalizationScope 实现
// ============================================================================

DelayedNormalizationScope::DelayedNormalizationScope()
    : old_state_(g_suppress_normalization) {
    g_suppress_normalization = true;
}

DelayedNormalizationScope::~DelayedNormalizationScope() {
    g_suppress_normalization = old_state_;
}

// ============================================================================
// 流输出运算符
// ============================================================================

namespace {

std::string trim_trailing_zeros(std::string text) {
    std::size_t dot = text.find('.');
    if (dot == std::string::npos) return text;
    while (!text.empty() && text.back() == '0') text.pop_back();
    if (!text.empty() && text.back() == '.') text.pop_back();
    if (text == "-0") return "0";
    return text.empty() ? "0" : text;
}

std::string format_decimal_truncated(const PreciseDecimal& value, int precision) {
    if (value.is_zero()) return "0";

    const bool negative = value < PreciseDecimal(0.0L);
    PreciseDecimal abs_value = negative ? -value : value;

    std::string text = abs_value.to_string();
    std::string digits;
    digits.reserve(text.size());
    std::size_t dot_digits = std::string::npos;
    for (char ch : text) {
        if (ch == '.') {
            dot_digits = digits.size();
        } else if (std::isdigit(static_cast<unsigned char>(ch))) {
            digits.push_back(ch);
        }
    }
    if (digits.empty()) return "0";
    if (dot_digits == std::string::npos) dot_digits = digits.size();

    std::size_t first_significant = 0;
    while (first_significant < digits.size() && digits[first_significant] == '0') {
        ++first_significant;
    }
    if (first_significant == digits.size()) return "0";

    const std::size_t keep_until =
        std::min(digits.size(), first_significant + static_cast<std::size_t>(precision));
    const std::string kept = digits.substr(0, keep_until);

    std::string result;
    if (negative) result.push_back('-');
    if (dot_digits == 0) {
        result += "0.";
        result += kept;
    } else if (kept.size() <= dot_digits) {
        result += kept;
    } else {
        result += kept.substr(0, dot_digits);
        result.push_back('.');
        result += kept.substr(dot_digits);
    }
    return trim_trailing_zeros(result);
}

int process_display_precision() {
    return PrecisionContext::get_default_scale();
}

} // namespace

std::ostream& operator<<(std::ostream& os, const PreciseDecimal& val) {
    int precision = os.precision() > 0 ? static_cast<int>(os.precision())
                                       : process_display_precision();
    os << format_decimal_truncated(val, precision);
    return os;
}

// ============================================================================
// multiply_add_without_normalize 实现
// ============================================================================

PreciseDecimal& PreciseDecimal::multiply_add_without_normalize(const PreciseDecimal& a, const PreciseDecimal& b) {
    PreciseDecimal product = multiply_precise_decimal(a, b);
    NormalizationSuppressor suppressor;
    return add_without_normalize(product);
}

} // namespace precise

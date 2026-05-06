// ============================================================================
// 精确小数实现 (优化版：基数 10^9)
// ============================================================================

#include "precise_decimal.h"
#include "rational.h"
#include "types/stored_value.h"
#include "core/calculator_exceptions.h"
#include "math/mymath.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr uint32_t kBase = 1000000000;
constexpr int kBaseDigits = 9;
int g_default_scale = 40;
constexpr double kDisplayZeroEps = 1e-16;
constexpr double kDisplayIntegerEps = 1e-9;

double normalize_display_decimal(double value) {
    if (mymath::abs(value) < kDisplayZeroEps) return 0.0;
    if (mymath::abs(value) > kDisplayIntegerEps &&
        mymath::abs(value - mymath::round(value)) < kDisplayIntegerEps) {
        return mymath::round(value);
    }
    return value;
}

std::string format_decimal(double value, int precision = 12) {
    value = normalize_display_decimal(value);
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    return out.str();
}

/**
 * @brief 计算 a / b，返回商和余数 (均为 BigInt)
 */
void div_bigint(const std::vector<uint32_t>& num, 
                const std::vector<uint32_t>& den,
                std::vector<uint32_t>* quotient,
                std::vector<uint32_t>* remainder);

/**
 * @brief 将整数字符串转换为基数为 10^9 的向量
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
 * @brief 将基数为 10^9 的向量转换为字符串
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
 * @brief 大整数减法 (要求 lhs >= rhs)
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

/**
 * @brief 大整数乘法
 */
std::vector<uint32_t> multiply_bigint(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs) {
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
 * @brief 大整数乘以一个小整数
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

void div_bigint(const std::vector<uint32_t>& num, 
                const std::vector<uint32_t>& den,
                std::vector<uint32_t>* quotient,
                std::vector<uint32_t>* remainder) {
    if (den.empty() || (den.size() == 1 && den[0] == 0)) {
        throw std::runtime_error("division by zero in div_bigint");
    }

    *remainder = {0};
    quotient->clear();
    
    std::string num_s = bigint_to_string(num);
    std::string den_s = bigint_to_string(den);
    std::string q_s;
    std::string rem_s = "0";

    for (char d : num_s) {
        if (rem_s == "0") rem_s = d;
        else rem_s.push_back(d);
        
        std::size_t first = rem_s.find_first_not_of('0');
        if (first == std::string::npos) rem_s = "0";
        else if (first > 0) rem_s = rem_s.substr(first);

        int q = 0;
        if (rem_s.size() > den_s.size() || (rem_s.size() == den_s.size() && rem_s >= den_s)) {
            for (int candidate = 9; candidate >= 1; --candidate) {
                std::vector<uint32_t> p_bigint = multiply_bigint_by_uint32(den, static_cast<uint32_t>(candidate));
                std::string p_s = bigint_to_string(p_bigint);
                if (p_s.size() < rem_s.size() || (p_s.size() == rem_s.size() && p_s <= rem_s)) {
                    q = candidate;
                    std::vector<uint32_t> rem_bigint = subtract_bigint(string_to_bigint(rem_s), p_bigint);
                    rem_s = bigint_to_string(rem_bigint);
                    break;
                }
            }
        }
        q_s.push_back(static_cast<char>('0' + q));
    }
    
    *quotient = string_to_bigint(q_s);
}

/**
 * @brief 大整数乘以 10^n
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

} // namespace

int PrecisionContext::get_default_scale() { return g_default_scale; }
void PrecisionContext::set_default_scale(int s) { g_default_scale = s; }

// ============================================================================
// PreciseDecimal 方法实现
// ============================================================================

// ============================================================================
// PreciseDecimal 方法实现
// ============================================================================

PreciseDecimal::PreciseDecimal(long long value) {
    negative = value < 0;
    data = string_to_bigint(std::to_string(negative ? -value : value));
    scale = 0;
}

PreciseDecimal::PreciseDecimal(double value) {
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
            uint32_t remainder = 0;
            for (int i = static_cast<int>(data.size()) - 1; i >= 0; --i) {
                uint64_t current = data[i] + static_cast<uint64_t>(remainder) * kBase;
                data[i] = static_cast<uint32_t>(current / 10);
                remainder = static_cast<uint32_t>(current % 10);
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

double PreciseDecimal::to_double() const {
    if (is_zero()) return 0.0;
    
    long double res = 0;
    long double p = 1.0L;
    for (uint32_t x : data) {
        res += static_cast<long double>(x) * p;
        p *= static_cast<long double>(kBase);
    }

    if (scale > 0) {
        res /= mymath::pow(10.0L, static_cast<long double>(scale));
    } else if (scale < 0) {
        res *= mymath::pow(10.0L, static_cast<long double>(-scale));
    }

    return static_cast<double>(negative ? -res : res);
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
    std::string significand = token;
    int exponent_adjust = 0;
    const std::size_t exponent_pos = token.find_first_of("eE");
    if (exponent_pos != std::string::npos) {
        significand = token.substr(0, exponent_pos);
        exponent_adjust = std::stoi(token.substr(exponent_pos + 1));
    }

    const std::size_t dot_pos = significand.find('.');
    if (dot_pos == std::string::npos) {
        return from_digits(significand, -exponent_adjust, false);
    }

    std::string digits_only = significand.substr(0, dot_pos);
    digits_only += significand.substr(dot_pos + 1);
    return from_digits(
        digits_only,
        static_cast<int>(significand.size() - dot_pos - 1) - exponent_adjust,
        false);
}

// ============================================================================
// PreciseDecimal 算术运算
// ============================================================================

PreciseDecimal operator+(PreciseDecimal lhs, const PreciseDecimal& rhs) { lhs += rhs; return lhs; }
PreciseDecimal operator+(PreciseDecimal lhs, double rhs) { lhs += PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator+(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) + rhs; }

PreciseDecimal operator-(PreciseDecimal lhs, const PreciseDecimal& rhs) { lhs -= rhs; return lhs; }
PreciseDecimal operator-(PreciseDecimal lhs, double rhs) { lhs -= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator-(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) - rhs; }

PreciseDecimal operator*(PreciseDecimal lhs, const PreciseDecimal& rhs) { lhs *= rhs; return lhs; }
PreciseDecimal operator*(PreciseDecimal lhs, double rhs) { lhs *= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator*(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) * rhs; }

PreciseDecimal operator/(PreciseDecimal lhs, const PreciseDecimal& rhs) { lhs /= rhs; return lhs; }
PreciseDecimal operator/(PreciseDecimal lhs, double rhs) { lhs /= PreciseDecimal(rhs); return lhs; }
PreciseDecimal operator/(double lhs, const PreciseDecimal& rhs) { return PreciseDecimal(lhs) / rhs; }

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
    // 初始猜测使用 double 版 sqrt
    PreciseDecimal x(mymath::sqrt(val.to_double()));
    const PreciseDecimal half(0.5);
    
    // 迭代 8 次通常足以达到 40+ 位精度 (Newton 迭代是二阶收敛)
    for (int i = 0; i < 8; ++i) {
        x = half * (x + val / x);
    }
    return x;
}

PreciseDecimal pow(const PreciseDecimal& base, long long exp) {
    return pow_precise_decimal(base, exp);
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
    std::vector<uint32_t> numerator = multiply_bigint_by_power_of_10(lhs.data, target_scale);
    std::vector<uint32_t> denominator = rhs.data;
    
    std::vector<uint32_t> q, r;
    div_bigint(numerator, denominator, &q, &r);

    PreciseDecimal res;
    res.data = q;
    res.scale = lhs.scale + target_scale - rhs.scale;
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

    PreciseDecimal res = PreciseDecimal::from_integer_string("1", false);
    PreciseDecimal b = base;
    while (exponent > 0) {
        if (exponent % 2 == 1) res = multiply_precise_decimal(res, b);
        b = multiply_precise_decimal(b, b);
        exponent /= 2;
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
    return format_decimal(normalize_display_decimal(value.decimal));
}

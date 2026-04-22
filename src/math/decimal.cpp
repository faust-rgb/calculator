/**
 * @file decimal.cpp
 * @brief 任意精度十进制数实现
 */

#include "decimal.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace decimal {

// 每个存储单元存储的十进制位数
constexpr int32_t kDigitsPerWord = 4;
constexpr uint32_t kWordBase = 10000;  // 10^4

// 默认最大精度（小数位数）
int32_t Decimal::max_precision_ = 50;

// ============================================================================
// 构造函数
// ============================================================================

Decimal::Decimal() : exponent_(0), sign_(1) {
    mantissa_.push_back(0);
}

Decimal::Decimal(int64_t value) : exponent_(0) {
    if (value < 0) {
        sign_ = -1;
        value = -value;
    } else {
        sign_ = 1;
    }

    if (value == 0) {
        mantissa_.push_back(0);
        return;
    }

    // 将值分解为 4 位一组的十进制数字
    uint64_t v = static_cast<uint64_t>(value);
    while (v > 0) {
        mantissa_.push_back(static_cast<uint16_t>(v % kWordBase));
        v /= kWordBase;
    }
}

Decimal::Decimal(const std::string& value) : exponent_(0), sign_(1) {
    std::string s = value;

    // 处理符号
    if (!s.empty() && s[0] == '-') {
        sign_ = -1;
        s = s.substr(1);
    } else if (!s.empty() && s[0] == '+') {
        s = s.substr(1);
    }

    // 处理科学计数法（如 1e-10, 1.5e+5）
    size_t e_pos = s.find_first_of("eE");
    int32_t exp_adjust = 0;
    if (e_pos != std::string::npos) {
        std::string exp_str = s.substr(e_pos + 1);
        s = s.substr(0, e_pos);
        // 解析指数
        bool exp_neg = false;
        size_t exp_start = 0;
        if (!exp_str.empty() && exp_str[0] == '-') {
            exp_neg = true;
            exp_start = 1;
        } else if (!exp_str.empty() && exp_str[0] == '+') {
            exp_start = 1;
        }
        for (size_t i = exp_start; i < exp_str.size(); ++i) {
            if (std::isdigit(static_cast<unsigned char>(exp_str[i]))) {
                exp_adjust = exp_adjust * 10 + (exp_str[i] - '0');
            }
        }
        if (exp_neg) {
            exp_adjust = -exp_adjust;
        }
    }

    // 找到小数点
    size_t dot_pos = s.find('.');
    std::string int_part, frac_part;
    if (dot_pos == std::string::npos) {
        int_part = s;
        frac_part = "";
    } else {
        int_part = s.substr(0, dot_pos);
        frac_part = s.substr(dot_pos + 1);
    }

    // 移除前导零
    size_t int_start = 0;
    while (int_start < int_part.size() && int_part[int_start] == '0') {
        ++int_start;
    }
    int_part = int_part.substr(int_start);

    // 移除后导零
    while (!frac_part.empty() && frac_part.back() == '0') {
        frac_part.pop_back();
    }

    // 计算指数调整
    exponent_ = -static_cast<int32_t>(frac_part.size()) + exp_adjust;

    // 合并整数部分和小数部分
    std::string all_digits = int_part + frac_part;

    if (all_digits.empty()) {
        mantissa_.push_back(0);
        sign_ = 1;
        exponent_ = 0;
        return;
    }

    // 将数字字符串从后向前每 4 位分组
    for (int i = static_cast<int>(all_digits.size()); i > 0; i -= kDigitsPerWord) {
        int start = std::max(0, i - kDigitsPerWord);
        int len = i - start;
        uint16_t word = 0;
        for (int j = 0; j < len; ++j) {
            char c = all_digits[start + j];
            if (!std::isdigit(static_cast<unsigned char>(c))) {
                throw std::invalid_argument("Invalid decimal string: " + value);
            }
            word = word * 10 + (c - '0');
        }
        // 如果不是完整的 4 位，需要补零
        for (int j = len; j < kDigitsPerWord; ++j) {
            word *= 10;
        }
        mantissa_.push_back(word);
    }

    normalize();
}

Decimal::Decimal(double value) : exponent_(0) {
    if (value < 0) {
        sign_ = -1;
        value = -value;
    } else {
        sign_ = 1;
    }

    if (std::isnan(value) || std::isinf(value)) {
        throw std::invalid_argument("Cannot convert NaN or Inf to Decimal");
    }

    if (value == 0.0) {
        mantissa_.push_back(0);
        return;
    }

    // 使用 ostringstream 进行精确转换
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    *this = Decimal(oss.str());
}

Decimal::Decimal(std::vector<uint16_t> mantissa, int32_t exponent, int8_t sign)
    : mantissa_(std::move(mantissa)), exponent_(exponent), sign_(sign) {
    normalize();
}

// ============================================================================
// 内部辅助函数
// ============================================================================

void Decimal::normalize() {
    // 移除尾数前导零
    while (mantissa_.size() > 1 && mantissa_.back() == 0) {
        mantissa_.pop_back();
    }

    // 零的标准化表示
    if (mantissa_.size() == 1 && mantissa_[0] == 0) {
        sign_ = 1;
        exponent_ = 0;
        return;
    }

    // 处理尾数末尾的零（可以合并到指数中）
    // 例如：123000 * 10^-3 可以变成 123 * 10^0
    while (!mantissa_.empty() && mantissa_[0] == 0) {
        // 检查是否可以整体除以 10
        bool can_divide = true;
        for (size_t i = 0; i < mantissa_.size(); ++i) {
            if (mantissa_[i] % 10 != 0) {
                can_divide = false;
                break;
            }
        }
        if (!can_divide) {
            break;
        }

        // 整体除以 10
        uint32_t carry = 0;
        for (int i = static_cast<int>(mantissa_.size()) - 1; i >= 0; --i) {
            uint32_t val = mantissa_[i] + carry * kWordBase;
            mantissa_[i] = static_cast<uint16_t>(val / 10);
            carry = val % 10;
        }
        ++exponent_;
    }

    // 再次移除前导零
    while (mantissa_.size() > 1 && mantissa_.back() == 0) {
        mantissa_.pop_back();
    }
}

void Decimal::carry_normalize(std::vector<uint16_t>* digits) {
    uint32_t carry = 0;
    for (size_t i = 0; i < digits->size() || carry > 0; ++i) {
        if (i >= digits->size()) {
            digits->push_back(0);
        }
        uint32_t sum = static_cast<uint32_t>((*digits)[i]) + carry;
        (*digits)[i] = static_cast<uint16_t>(sum % kWordBase);
        carry = sum / kWordBase;
    }
}

int Decimal::compare_abs(const Decimal& rhs) const {
    // 先对齐指数
    Decimal lhs_copy = *this;
    Decimal rhs_copy = rhs;
    align_exponents(&lhs_copy, &rhs_copy);

    // 从高位到低位比较
    if (lhs_copy.mantissa_.size() != rhs_copy.mantissa_.size()) {
        return lhs_copy.mantissa_.size() < rhs_copy.mantissa_.size() ? -1 : 1;
    }
    for (int i = static_cast<int>(lhs_copy.mantissa_.size()) - 1; i >= 0; --i) {
        if (lhs_copy.mantissa_[i] != rhs_copy.mantissa_[i]) {
            return lhs_copy.mantissa_[i] < rhs_copy.mantissa_[i] ? -1 : 1;
        }
    }
    return 0;
}

void Decimal::align_exponents(Decimal* lhs, Decimal* rhs) {
    if (lhs->exponent_ == rhs->exponent_) {
        return;
    }

    // 将指数较小的数（绝对值更大的数）进行缩放
    if (lhs->exponent_ < rhs->exponent_) {
        int32_t diff = rhs->exponent_ - lhs->exponent_;
        // lhs 需要乘以 10^diff
        for (int32_t i = 0; i < diff; ++i) {
            // 乘以 10
            uint32_t carry = 0;
            for (size_t j = 0; j < lhs->mantissa_.size(); ++j) {
                uint32_t prod = lhs->mantissa_[j] * 10 + carry;
                lhs->mantissa_[j] = static_cast<uint16_t>(prod % kWordBase);
                carry = prod / kWordBase;
            }
            if (carry > 0) {
                lhs->mantissa_.push_back(static_cast<uint16_t>(carry));
            }
        }
        lhs->exponent_ = rhs->exponent_;
    } else {
        int32_t diff = lhs->exponent_ - rhs->exponent_;
        for (int32_t i = 0; i < diff; ++i) {
            uint32_t carry = 0;
            for (size_t j = 0; j < rhs->mantissa_.size(); ++j) {
                uint32_t prod = rhs->mantissa_[j] * 10 + carry;
                rhs->mantissa_[j] = static_cast<uint16_t>(prod % kWordBase);
                carry = prod / kWordBase;
            }
            if (carry > 0) {
                rhs->mantissa_.push_back(static_cast<uint16_t>(carry));
            }
        }
        rhs->exponent_ = lhs->exponent_;
    }
}

void Decimal::shift_exponent(int32_t delta) {
    exponent_ += delta;
}

// ============================================================================
// 算术运算
// ============================================================================

Decimal Decimal::operator+(const Decimal& rhs) const {
    if (sign_ != rhs.sign_) {
        // 异号相加 = 减法
        Decimal rhs_neg = rhs;
        rhs_neg.sign_ = -rhs_neg.sign_;
        return *this - rhs_neg;
    }

    // 同号相加
    Decimal lhs_copy = *this;
    Decimal rhs_copy = rhs;
    align_exponents(&lhs_copy, &rhs_copy);

    // 对齐位数
    size_t max_size = std::max(lhs_copy.mantissa_.size(), rhs_copy.mantissa_.size());
    lhs_copy.mantissa_.resize(max_size, 0);
    rhs_copy.mantissa_.resize(max_size, 0);

    // 逐位相加
    std::vector<uint16_t> result;
    uint32_t carry = 0;
    for (size_t i = 0; i < max_size; ++i) {
        uint32_t sum = lhs_copy.mantissa_[i] + rhs_copy.mantissa_[i] + carry;
        result.push_back(static_cast<uint16_t>(sum % kWordBase));
        carry = sum / kWordBase;
    }
    if (carry > 0) {
        result.push_back(static_cast<uint16_t>(carry));
    }

    return Decimal(result, lhs_copy.exponent_, sign_);
}

Decimal Decimal::operator-(const Decimal& rhs) const {
    if (sign_ != rhs.sign_) {
        // 异号相减 = 加法
        Decimal rhs_neg = rhs;
        rhs_neg.sign_ = -rhs_neg.sign_;
        return *this + rhs_neg;
    }

    // 同号相减，比较绝对值
    int cmp = compare_abs(rhs);
    if (cmp == 0) {
        return Decimal();  // 结果为零
    }

    const Decimal* larger = (cmp > 0) ? this : &rhs;
    const Decimal* smaller = (cmp > 0) ? &rhs : this;
    int8_t result_sign = (cmp > 0) ? sign_ : -sign_;

    Decimal large_copy = *larger;
    Decimal small_copy = *smaller;
    align_exponents(&large_copy, &small_copy);

    // 对齐位数
    size_t max_size = std::max(large_copy.mantissa_.size(), small_copy.mantissa_.size());
    large_copy.mantissa_.resize(max_size, 0);
    small_copy.mantissa_.resize(max_size, 0);

    // 逐位相减
    std::vector<uint16_t> result;
    int32_t borrow = 0;
    for (size_t i = 0; i < max_size; ++i) {
        int32_t diff = static_cast<int32_t>(large_copy.mantissa_[i]) -
                       static_cast<int32_t>(small_copy.mantissa_[i]) - borrow;
        if (diff < 0) {
            diff += kWordBase;
            borrow = 1;
        } else {
            borrow = 0;
        }
        result.push_back(static_cast<uint16_t>(diff));
    }

    return Decimal(result, large_copy.exponent_, result_sign);
}

std::vector<uint16_t> Decimal::multiply_mantissa(
    const std::vector<uint16_t>& a,
    const std::vector<uint16_t>& b) {
    std::vector<uint32_t> temp(a.size() + b.size(), 0);

    for (size_t i = 0; i < a.size(); ++i) {
        uint32_t carry = 0;
        for (size_t j = 0; j < b.size() || carry > 0; ++j) {
            uint32_t prod = temp[i + j] + a[i] * (j < b.size() ? b[j] : 0) + carry;
            temp[i + j] = prod % kWordBase;
            carry = prod / kWordBase;
        }
    }

    // 转换回 uint16_t
    std::vector<uint16_t> result;
    for (uint32_t val : temp) {
        result.push_back(static_cast<uint16_t>(val));
    }

    return result;
}

Decimal Decimal::operator*(const Decimal& rhs) const {
    int8_t result_sign = sign_ * rhs.sign_;

    if (is_zero() || rhs.is_zero()) {
        return Decimal();
    }

    std::vector<uint16_t> result_mantissa = multiply_mantissa(mantissa_, rhs.mantissa_);
    int32_t result_exponent = exponent_ + rhs.exponent_;

    // 调整指数（因为两个尾数都包含了各自的基数）
    // 实际上指数已经正确反映了小数点位置

    Decimal result(result_mantissa, result_exponent, result_sign);

    // 限制精度
    int32_t current_scale = -result.exponent_;
    if (current_scale > max_precision_) {
        result = result.round(max_precision_);
    }

    return result;
}

std::vector<uint16_t> Decimal::divide_mantissa(
    const std::vector<uint16_t>& dividend,
    const std::vector<uint16_t>& divisor,
    std::vector<uint16_t>* remainder,
    size_t max_digits) {
    // 长除法实现
    std::vector<uint16_t> quotient;
    std::vector<uint16_t> rem = dividend;

    // 将除数和被除数转换为更大的基数以便比较
    auto compare = [](const std::vector<uint16_t>& a, const std::vector<uint16_t>& b) {
        size_t len_a = a.size();
        while (len_a > 1 && a[len_a - 1] == 0) --len_a;
        size_t len_b = b.size();
        while (len_b > 1 && b[len_b - 1] == 0) --len_b;

        if (len_a != len_b) {
            return len_a < len_b ? -1 : 1;
        }
        for (int i = static_cast<int>(len_a) - 1; i >= 0; --i) {
            if (a[i] != b[i]) {
                return a[i] < b[i] ? -1 : 1;
            }
        }
        return 0;
    };

    auto subtract = [](std::vector<uint16_t>& a, const std::vector<uint16_t>& b) {
        int32_t borrow = 0;
        for (size_t i = 0; i < b.size(); ++i) {
            int32_t diff = a[i] - b[i] - borrow;
            if (diff < 0) {
                diff += kWordBase;
                borrow = 1;
            } else {
                borrow = 0;
            }
            a[i] = static_cast<uint16_t>(diff);
        }
        for (size_t i = b.size(); i < a.size() && borrow > 0; ++i) {
            int32_t diff = a[i] - borrow;
            if (diff < 0) {
                diff += kWordBase;
                borrow = 1;
            } else {
                borrow = 0;
            }
            a[i] = static_cast<uint16_t>(diff);
        }
    };

    // 主除法循环
    for (size_t digit = 0; digit < max_digits; ++digit) {
        // 找到当前商的数字
        uint32_t low = 0, high = kWordBase;
        while (low < high) {
            uint32_t mid = (low + high + 1) / 2;
            // 计算 divisor * mid
            std::vector<uint16_t> prod;
            uint32_t carry = 0;
            for (size_t i = 0; i < divisor.size(); ++i) {
                uint32_t p = divisor[i] * mid + carry;
                prod.push_back(static_cast<uint16_t>(p % kWordBase));
                carry = p / kWordBase;
            }
            if (carry > 0) {
                prod.push_back(static_cast<uint16_t>(carry));
            }

            if (compare(prod, rem) <= 0) {
                low = mid;
            } else {
                high = mid - 1;
            }
        }

        quotient.push_back(static_cast<uint16_t>(low));

        // 计算余数
        if (low > 0) {
            std::vector<uint16_t> prod;
            uint32_t carry = 0;
            for (size_t i = 0; i < divisor.size(); ++i) {
                uint32_t p = divisor[i] * low + carry;
                prod.push_back(static_cast<uint16_t>(p % kWordBase));
                carry = p / kWordBase;
            }
            if (carry > 0) {
                prod.push_back(static_cast<uint16_t>(carry));
            }
            subtract(rem, prod);
        }

        // 余数乘以基数（相当于添加一个零）
        rem.insert(rem.begin(), 0);
    }

    if (remainder) {
        *remainder = rem;
    }

    return quotient;
}

Decimal Decimal::operator/(const Decimal& rhs) const {
    if (rhs.is_zero()) {
        throw std::runtime_error("Division by zero");
    }

    if (is_zero()) {
        return Decimal();
    }

    int8_t result_sign = sign_ * rhs.sign_;

    // 为了获得足够的精度，将被除数扩展
    std::vector<uint16_t> extended_dividend = mantissa_;
    // 添加足够的零以获得所需精度
    size_t extra_digits = (max_precision_ + kDigitsPerWord - 1) / kDigitsPerWord + 2;
    for (size_t i = 0; i < extra_digits; ++i) {
        extended_dividend.insert(extended_dividend.begin(), 0);
    }

    std::vector<uint16_t> quotient = divide_mantissa(
        extended_dividend, rhs.mantissa_, nullptr, mantissa_.size() + extra_digits);

    int32_t result_exponent = exponent_ - rhs.exponent_ - static_cast<int32_t>(extra_digits) * kDigitsPerWord;

    Decimal result(quotient, result_exponent, result_sign);

    // 限制精度
    int32_t current_scale = -result.exponent_;
    if (current_scale > max_precision_) {
        result = result.round(max_precision_);
    }

    return result;
}

Decimal& Decimal::operator+=(const Decimal& rhs) {
    *this = *this + rhs;
    return *this;
}

Decimal& Decimal::operator-=(const Decimal& rhs) {
    *this = *this - rhs;
    return *this;
}

Decimal& Decimal::operator*=(const Decimal& rhs) {
    *this = *this * rhs;
    return *this;
}

Decimal& Decimal::operator/=(const Decimal& rhs) {
    *this = *this / rhs;
    return *this;
}

Decimal Decimal::operator-() const {
    Decimal result = *this;
    if (!is_zero()) {
        result.sign_ = -sign_;
    }
    return result;
}

// ============================================================================
// 比较运算
// ============================================================================

bool Decimal::operator==(const Decimal& rhs) const {
    return sign_ == rhs.sign_ && exponent_ == rhs.exponent_ && mantissa_ == rhs.mantissa_;
}

bool Decimal::operator!=(const Decimal& rhs) const {
    return !(*this == rhs);
}

bool Decimal::operator<(const Decimal& rhs) const {
    if (sign_ != rhs.sign_) {
        return sign_ < rhs.sign_;
    }
    if (is_zero() && rhs.is_zero()) {
        return false;
    }
    int cmp = compare_abs(rhs);
    return sign_ > 0 ? cmp < 0 : cmp > 0;
}

bool Decimal::operator>(const Decimal& rhs) const {
    return rhs < *this;
}

bool Decimal::operator<=(const Decimal& rhs) const {
    return !(rhs < *this);
}

bool Decimal::operator>=(const Decimal& rhs) const {
    return !(*this < rhs);
}

// ============================================================================
// 类型转换
// ============================================================================

std::string Decimal::to_string() const {
    if (is_zero()) {
        return "0";
    }

    // 将尾数转换为字符串
    std::string digits;
    for (int i = static_cast<int>(mantissa_.size()) - 1; i >= 0; --i) {
        std::ostringstream oss;
        if (i < static_cast<int>(mantissa_.size()) - 1) {
            oss << std::setw(kDigitsPerWord) << std::setfill('0');
        }
        oss << mantissa_[i];
        digits += oss.str();
    }

    // 移除前导零
    size_t start = 0;
    while (start < digits.size() && digits[start] == '0') {
        ++start;
    }
    digits = digits.substr(start);

    if (digits.empty()) {
        return "0";
    }

    // 计算小数点位置
    int32_t total_digits = static_cast<int32_t>(digits.size());
    int32_t int_digits = total_digits + exponent_;

    std::string result;
    if (sign_ < 0) {
        result = "-";
    }

    if (int_digits > 0) {
        // 小数点在数字内部或之后
        if (int_digits >= total_digits) {
            // 没有小数部分，只有整数部分
            result += digits;
            for (int i = 0; i < int_digits - total_digits; ++i) {
                result += '0';
            }
        } else {
            // 有小数部分
            result += digits.substr(0, int_digits);
            result += '.';
            result += digits.substr(int_digits);
        }
    } else {
        // 小数点在数字之前，需要添加前导零
        result += "0.";
        for (int i = 0; i < -int_digits; ++i) {
            result += '0';
        }
        result += digits;
    }

    // 移除末尾的小数点后的零
    size_t dot_pos = result.find('.');
    if (dot_pos != std::string::npos) {
        size_t end = result.size();
        while (end > dot_pos + 1 && result[end - 1] == '0') {
            --end;
        }
        if (end == dot_pos + 1) {
            // 移除小数点
            result = result.substr(0, dot_pos);
        } else {
            result = result.substr(0, end);
        }
    }

    return result;
}

double Decimal::to_double() const {
    std::string str = to_string();
    try {
        return std::stod(str);
    } catch (...) {
        // 如果字符串太长，使用近似方法
        double result = 0.0;
        double factor = std::pow(10.0, exponent_);
        for (int i = static_cast<int>(mantissa_.size()) - 1; i >= 0; --i) {
            result = result * kWordBase + mantissa_[i];
        }
        return sign_ * result * factor;
    }
}

int64_t Decimal::to_int64() const {
    Decimal truncated = this->truncate(0);
    std::string str = truncated.to_string();
    try {
        return std::stoll(str);
    } catch (...) {
        throw std::overflow_error("Decimal value out of int64 range");
    }
}

bool Decimal::is_zero() const {
    return mantissa_.size() == 1 && mantissa_[0] == 0;
}

bool Decimal::is_integer() const {
    if (is_zero()) {
        return true;
    }
    // 检查指数是否 >= 0（即小数点在最右边或更右）
    if (exponent_ >= 0) {
        return true;
    }
    // 检查小数部分是否全为零
    std::string str = to_string();
    size_t dot_pos = str.find('.');
    if (dot_pos == std::string::npos) {
        return true;
    }
    for (size_t i = dot_pos + 1; i < str.size(); ++i) {
        if (str[i] != '0') {
            return false;
        }
    }
    return true;
}

bool Decimal::is_negative() const {
    return sign_ < 0 && !is_zero();
}

bool Decimal::is_positive() const {
    return sign_ > 0 && !is_zero();
}

// ============================================================================
// 精度控制
// ============================================================================

void Decimal::set_max_precision(int32_t precision) {
    if (precision < 0) {
        throw std::invalid_argument("Precision must be non-negative");
    }
    max_precision_ = precision;
}

int32_t Decimal::get_max_precision() {
    return max_precision_;
}

Decimal Decimal::round(int32_t decimal_places) const {
    if (decimal_places < 0) {
        throw std::invalid_argument("Decimal places must be non-negative");
    }

    int32_t current_scale = -exponent_;
    if (current_scale <= decimal_places) {
        return *this;  // 不需要舍入
    }

    // 需要舍入的位置
    int32_t drop_digits = current_scale - decimal_places;

    // 创建舍入用的 Decimal
    Decimal factor(std::vector<uint16_t>{1}, -drop_digits, 1);
    Decimal half_factor = factor * Decimal("0.5");

    // 加上半个舍入单位（用于四舍五入）
    Decimal adjusted = is_negative() ? *this - half_factor : *this + half_factor;

    // 截断
    return adjusted.truncate(decimal_places);
}

Decimal Decimal::truncate(int32_t decimal_places) const {
    if (decimal_places < 0) {
        throw std::invalid_argument("Decimal places must be non-negative");
    }

    std::string str = to_string();
    size_t dot_pos = str.find('.');

    if (dot_pos == std::string::npos) {
        return *this;  // 没有小数部分
    }

    size_t keep_digits = dot_pos + 1 + decimal_places;
    if (keep_digits >= str.size()) {
        return *this;  // 小数位数不足
    }

    std::string truncated = str.substr(0, keep_digits);
    // 移除末尾的小数点
    if (truncated.back() == '.') {
        truncated.pop_back();
    }

    return Decimal(truncated);
}

int32_t Decimal::scale() const {
    if (is_zero()) {
        return 0;
    }
    int32_t s = -exponent_;
    return s > 0 ? s : 0;
}

// ============================================================================
// 字面量支持
// ============================================================================

Decimal operator"" _dec(const char* str, size_t len) {
    return Decimal(std::string(str, len));
}

}  // namespace decimal

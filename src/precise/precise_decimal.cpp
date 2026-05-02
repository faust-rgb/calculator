// ============================================================================
// 精确小数实现
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

namespace {

constexpr int kPreciseDecimalDivisionDigits = 40;

std::string trim_leading_zeros(std::string digits) {
    std::size_t first_non_zero = digits.find_first_not_of('0');
    if (first_non_zero == std::string::npos) {
        return "0";
    }
    return digits.substr(first_non_zero);
}

int compare_unsigned_digits(const std::string& lhs, const std::string& rhs) {
    const std::string normalized_lhs = trim_leading_zeros(lhs);
    const std::string normalized_rhs = trim_leading_zeros(rhs);
    if (normalized_lhs.size() != normalized_rhs.size()) {
        return normalized_lhs.size() < normalized_rhs.size() ? -1 : 1;
    }
    if (normalized_lhs == normalized_rhs) {
        return 0;
    }
    return normalized_lhs < normalized_rhs ? -1 : 1;
}

std::string add_unsigned_digits(const std::string& lhs, const std::string& rhs) {
    std::string result;
    result.reserve(std::max(lhs.size(), rhs.size()) + 1);

    int carry = 0;
    int lhs_index = static_cast<int>(lhs.size()) - 1;
    int rhs_index = static_cast<int>(rhs.size()) - 1;
    while (lhs_index >= 0 || rhs_index >= 0 || carry != 0) {
        int sum = carry;
        if (lhs_index >= 0) {
            sum += lhs[static_cast<std::size_t>(lhs_index)] - '0';
            --lhs_index;
        }
        if (rhs_index >= 0) {
            sum += rhs[static_cast<std::size_t>(rhs_index)] - '0';
            --rhs_index;
        }
        result.push_back(static_cast<char>('0' + (sum % 10)));
        carry = sum / 10;
    }

    std::reverse(result.begin(), result.end());
    return trim_leading_zeros(result);
}

std::string subtract_unsigned_digits(const std::string& lhs, const std::string& rhs) {
    std::string result;
    result.reserve(lhs.size());

    int borrow = 0;
    int lhs_index = static_cast<int>(lhs.size()) - 1;
    int rhs_index = static_cast<int>(rhs.size()) - 1;
    while (lhs_index >= 0) {
        int diff = (lhs[static_cast<std::size_t>(lhs_index)] - '0') - borrow;
        if (rhs_index >= 0) {
            diff -= rhs[static_cast<std::size_t>(rhs_index)] - '0';
            --rhs_index;
        }
        if (diff < 0) {
            diff += 10;
            borrow = 1;
        } else {
            borrow = 0;
        }
        result.push_back(static_cast<char>('0' + diff));
        --lhs_index;
    }

    std::reverse(result.begin(), result.end());
    return trim_leading_zeros(result);
}

std::string multiply_unsigned_digits(const std::string& lhs, const std::string& rhs) {
    const std::string normalized_lhs = trim_leading_zeros(lhs);
    const std::string normalized_rhs = trim_leading_zeros(rhs);
    if (normalized_lhs == "0" || normalized_rhs == "0") {
        return "0";
    }

    std::vector<int> digits(normalized_lhs.size() + normalized_rhs.size(), 0);
    for (int i = static_cast<int>(normalized_lhs.size()) - 1; i >= 0; --i) {
        for (int j = static_cast<int>(normalized_rhs.size()) - 1; j >= 0; --j) {
            const int product =
                (normalized_lhs[static_cast<std::size_t>(i)] - '0') *
                (normalized_rhs[static_cast<std::size_t>(j)] - '0');
            const std::size_t pos_low =
                static_cast<std::size_t>(i + j + 1);
            const std::size_t pos_high =
                static_cast<std::size_t>(i + j);
            const int sum = digits[pos_low] + product;
            digits[pos_low] = sum % 10;
            digits[pos_high] += sum / 10;
        }
    }

    std::string result;
    result.reserve(digits.size());
    for (int digit : digits) {
        if (result.empty() && digit == 0) {
            continue;
        }
        result.push_back(static_cast<char>('0' + digit));
    }
    return result.empty() ? "0" : result;
}

std::string multiply_unsigned_digits_by_digit(const std::string& value, int digit) {
    if (digit < 0 || digit > 9) {
        throw std::runtime_error("internal error: invalid multiplication digit");
    }
    if (digit == 0 || value == "0") {
        return "0";
    }
    if (digit == 1) {
        return trim_leading_zeros(value);
    }

    std::string result;
    result.reserve(value.size() + 1);

    int carry = 0;
    for (int i = static_cast<int>(value.size()) - 1; i >= 0; --i) {
        const int product =
            (value[static_cast<std::size_t>(i)] - '0') * digit + carry;
        result.push_back(static_cast<char>('0' + (product % 10)));
        carry = product / 10;
    }
    while (carry != 0) {
        result.push_back(static_cast<char>('0' + (carry % 10)));
        carry /= 10;
    }

    std::reverse(result.begin(), result.end());
    return trim_leading_zeros(result);
}

void align_precise_scales(PreciseDecimal* lhs, PreciseDecimal* rhs) {
    if (lhs->scale == rhs->scale) {
        return;
    }
    if (lhs->scale < rhs->scale) {
        lhs->digits.append(static_cast<std::size_t>(rhs->scale - lhs->scale), '0');
        lhs->scale = rhs->scale;
        return;
    }
    rhs->digits.append(static_cast<std::size_t>(lhs->scale - rhs->scale), '0');
    rhs->scale = lhs->scale;
}

} // namespace

// ============================================================================
// PreciseDecimal 方法实现
// ============================================================================

void PreciseDecimal::normalize() {
    digits = trim_leading_zeros(digits);
    if (scale < 0) {
        digits.append(static_cast<std::size_t>(-scale), '0');
        scale = 0;
    }
    while (scale > 0 && digits.size() > 1 && digits.back() == '0') {
        digits.pop_back();
        --scale;
    }
    if (digits == "0") {
        negative = false;
        scale = 0;
    }
}

bool PreciseDecimal::is_zero() const {
    return digits == "0";
}

std::string PreciseDecimal::to_string() const {
    if (digits == "0") {
        return "0";
    }

    std::string result;
    if (negative) {
        result.push_back('-');
    }

    if (scale == 0) {
        result += digits;
        return result;
    }

    if (static_cast<int>(digits.size()) <= scale) {
        result += "0.";
        result.append(static_cast<std::size_t>(scale) - digits.size(), '0');
        result += digits;
        return result;
    }

    const std::size_t split =
        digits.size() - static_cast<std::size_t>(scale);
    result += digits.substr(0, split);
    result.push_back('.');
    result += digits.substr(split);
    return result;
}

double PreciseDecimal::to_double() const {
    return std::stod(to_string());
}

PreciseDecimal PreciseDecimal::from_digits(std::string raw_digits,
                                           int raw_scale,
                                           bool is_negative) {
    PreciseDecimal value;
    value.digits = std::move(raw_digits);
    value.scale = raw_scale;
    value.negative = is_negative;
    value.normalize();
    return value;
}

PreciseDecimal PreciseDecimal::from_integer_string(
    const std::string& integer_text,
    bool is_negative) {
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

PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs,
                                   const PreciseDecimal& rhs) {
    if (lhs.negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return subtract_precise_decimal(lhs, rhs_flipped);
    }

    PreciseDecimal lhs_copy = lhs;
    PreciseDecimal rhs_copy = rhs;
    align_precise_scales(&lhs_copy, &rhs_copy);
    return PreciseDecimal::from_digits(
        add_unsigned_digits(lhs_copy.digits, rhs_copy.digits),
        lhs_copy.scale,
        lhs_copy.negative);
}

PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs,
                                        const PreciseDecimal& rhs) {
    if (lhs.negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return add_precise_decimal(lhs, rhs_flipped);
    }

    PreciseDecimal lhs_copy = lhs;
    PreciseDecimal rhs_copy = rhs;
    align_precise_scales(&lhs_copy, &rhs_copy);

    const int comparison = compare_unsigned_digits(lhs_copy.digits, rhs_copy.digits);
    if (comparison == 0) {
        return {};
    }
    if (comparison > 0) {
        return PreciseDecimal::from_digits(
            subtract_unsigned_digits(lhs_copy.digits, rhs_copy.digits),
            lhs_copy.scale,
            lhs_copy.negative);
    }
    return PreciseDecimal::from_digits(
        subtract_unsigned_digits(rhs_copy.digits, lhs_copy.digits),
        lhs_copy.scale,
        !lhs_copy.negative);
}

PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs,
                                        const PreciseDecimal& rhs) {
    return PreciseDecimal::from_digits(
        multiply_unsigned_digits(lhs.digits, rhs.digits),
        lhs.scale + rhs.scale,
        lhs.negative != rhs.negative);
}

PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs,
                                      const PreciseDecimal& rhs) {
    if (rhs.is_zero()) {
        throw std::runtime_error("division by zero");
    }
    if (lhs.is_zero()) {
        return {};
    }

    std::string numerator = lhs.digits;
    numerator.append(static_cast<std::size_t>(rhs.scale), '0');
    const std::string denominator =
        rhs.digits + std::string(static_cast<std::size_t>(lhs.scale), '0');

    std::string integer_part;
    std::string fractional_part;
    std::string remainder = "0";
    bool emitted_integer_digit = false;

    for (char digit : numerator) {
        if (remainder == "0") {
            remainder.assign(1, digit);
        } else {
            remainder.push_back(digit);
        }
        remainder = trim_leading_zeros(remainder);

        int quotient_digit = 0;
        for (int candidate = 9; candidate >= 1; --candidate) {
            if (compare_unsigned_digits(
                    multiply_unsigned_digits_by_digit(denominator, candidate),
                    remainder) <= 0) {
                quotient_digit = candidate;
                break;
            }
        }
        if (quotient_digit != 0) {
            remainder = subtract_unsigned_digits(
                remainder,
                multiply_unsigned_digits_by_digit(denominator, quotient_digit));
        }
        if (quotient_digit != 0 || emitted_integer_digit) {
            integer_part.push_back(static_cast<char>('0' + quotient_digit));
            emitted_integer_digit = true;
        }
    }

    if (!emitted_integer_digit) {
        integer_part = "0";
    }

    int produced_fractional_digits = 0;
    while (remainder != "0" &&
           produced_fractional_digits < kPreciseDecimalDivisionDigits) {
        remainder.push_back('0');
        remainder = trim_leading_zeros(remainder);

        int quotient_digit = 0;
        for (int candidate = 9; candidate >= 1; --candidate) {
            if (compare_unsigned_digits(
                    multiply_unsigned_digits_by_digit(denominator, candidate),
                    remainder) <= 0) {
                quotient_digit = candidate;
                break;
            }
        }
        if (quotient_digit != 0) {
            remainder = subtract_unsigned_digits(
                remainder,
                multiply_unsigned_digits_by_digit(denominator, quotient_digit));
        }
        fractional_part.push_back(static_cast<char>('0' + quotient_digit));
        ++produced_fractional_digits;
    }

    return PreciseDecimal::from_digits(
        integer_part + fractional_part,
        static_cast<int>(fractional_part.size()),
        lhs.negative != rhs.negative);
}

int compare_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (lhs.is_zero() && rhs.is_zero()) return 0;
    if (lhs.negative != rhs.negative) {
        return lhs.negative ? -1 : 1;
    }

    PreciseDecimal l = lhs;
    PreciseDecimal r = rhs;
    align_precise_scales(&l, &r);

    int res = compare_unsigned_digits(l.digits, r.digits);
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

namespace {

constexpr double kDisplayZeroEps = 1e-16;
constexpr double kDisplayIntegerEps = 1e-9;

double normalize_display_decimal(double value) {
    if (mymath::abs(value) < kDisplayZeroEps) {
        return 0.0;
    }
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

} // namespace

std::string stored_value_precise_decimal_text(const StoredValue& value) {
    if (value.exact) {
        return rational_to_precise_decimal_text(value.rational);
    }
    if (value.has_precise_decimal_text) {
        return value.precise_decimal_text;
    }
    return format_decimal(normalize_display_decimal(value.decimal));
}

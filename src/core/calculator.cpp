#include "calculator.h"

#include "function_analysis.h"
#include "matrix.h"
#include "multivariable_integrator.h"
#include "mymath.h"
#include "ode_solver.h"
#include "polynomial.h"
#include "script_parser.h"
#include "symbolic_expression.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kDisplayZeroEps = mymath::kDoubleDenormMin;
constexpr double kDisplayIntegerEps = 1e-9;

class ExactModeUnsupported : public std::runtime_error {
public:
    explicit ExactModeUnsupported(const std::string& message)
        : std::runtime_error(message) {}
};

long long gcd_ll(long long a, long long b) {
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    while (b != 0) {
        const long long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

long long lcm_ll(long long a, long long b) {
    if (a == 0 || b == 0) {
        return 0;
    }
    const long long result = (a / gcd_ll(a, b)) * b;
    return result < 0 ? -result : result;
}

bool is_integer_double(double x, double eps = 1e-10) {
    return mymath::is_integer(x, eps);
}

long long round_to_long_long(double x) {
    return static_cast<long long>(x >= 0.0 ? x + 0.5 : x - 0.5);
}

long long trunc_to_long_long(double x) {
    return static_cast<long long>(x);
}

long long floor_to_long_long(double x) {
    long long truncated = static_cast<long long>(x);
    if (x < 0.0 && static_cast<double>(truncated) != x) {
        --truncated;
    }
    return truncated;
}

long long ceil_to_long_long(double x) {
    long long truncated = static_cast<long long>(x);
    if (x > 0.0 && static_cast<double>(truncated) != x) {
        ++truncated;
    }
    return truncated;
}

double normalize_display_decimal(double value) {
    if (mymath::is_near_zero(value, kDisplayZeroEps)) {
        return 0.0;
    }
    if (mymath::abs(value) > kDisplayIntegerEps &&
        is_integer_double(value, kDisplayIntegerEps)) {
        return static_cast<double>(round_to_long_long(value));
    }
    return value;
}

std::mt19937_64& global_rng() {
    static std::mt19937_64 engine(std::random_device{}());
    return engine;
}

bool lookup_builtin_constant(const std::string& name, double* value) {
    if (name == "pi") {
        *value = mymath::kPi;
        return true;
    }
    if (name == "e") {
        *value = mymath::kE;
        return true;
    }
    if (name == "c") {
        *value = 299792458.0;
        return true;
    }
    if (name == "G") {
        *value = 6.67430e-11;
        return true;
    }
    if (name == "h") {
        *value = 6.62607015e-34;
        return true;
    }
    if (name == "k") {
        *value = 1.380649e-23;
        return true;
    }
    if (name == "NA") {
        *value = 6.02214076e23;
        return true;
    }
    return false;
}

double degrees_to_radians(double value) {
    return value * mymath::kPi / 180.0;
}

double radians_to_degrees(double value) {
    return value * 180.0 / mymath::kPi;
}

double celsius_to_fahrenheit(double value) {
    return value * 9.0 / 5.0 + 32.0;
}

double fahrenheit_to_celsius(double value) {
    return (value - 32.0) * 5.0 / 9.0;
}

double normal_pdf(double x, double mean, double sigma) {
    if (sigma <= 0.0) {
        throw std::runtime_error("normal distribution sigma must be positive");
    }
    const double z = (x - mean) / sigma;
    return mymath::exp(-0.5 * z * z) /
           (sigma * mymath::sqrt(2.0 * mymath::kPi));
}

double normal_cdf(double x, double mean, double sigma) {
    if (sigma <= 0.0) {
        throw std::runtime_error("normal distribution sigma must be positive");
    }
    return 0.5 * (1.0 + mymath::erf((x - mean) /
                                    (sigma * mymath::sqrt(2.0))));
}

bool is_prime_ll(long long value) {
    if (value <= 1) {
        return false;
    }
    if (value <= 3) {
        return true;
    }
    if (value % 2 == 0 || value % 3 == 0) {
        return false;
    }
    for (long long i = 5; i * i <= value; i += 6) {
        if (value % i == 0 || value % (i + 2) == 0) {
            return false;
        }
    }
    return true;
}

long long next_prime_ll(long long value) {
    long long candidate = value <= 2 ? 2 : value + 1;
    if (candidate % 2 == 0 && candidate != 2) {
        ++candidate;
    }
    while (!is_prime_ll(candidate)) {
        candidate += (candidate == 2 ? 1 : 2);
    }
    return candidate;
}

double fibonacci_value(long long n) {
    if (n < 0) {
        throw std::runtime_error("fib only accepts non-negative integers");
    }
    if (n > 186) {
        throw std::runtime_error("fib is limited to n <= 186 to avoid overflow");
    }
    if (n == 0) {
        return 0.0;
    }
    long long a = 0;
    long long b = 1;
    for (long long i = 1; i < n; ++i) {
        const long long next = a + b;
        a = b;
        b = next;
    }
    return static_cast<double>(b);
}

constexpr unsigned kProgrammerBitWidth = 64;

std::uint64_t to_unsigned_bits(long long value) {
    return static_cast<std::uint64_t>(value);
}

long long from_unsigned_bits(std::uint64_t value) {
    return static_cast<long long>(value);
}

unsigned normalize_rotation_count(long long count) {
    if (count < 0) {
        throw std::runtime_error("rotate count cannot be negative");
    }
    return static_cast<unsigned>(count % static_cast<long long>(kProgrammerBitWidth));
}

std::uint64_t rotate_left_bits(std::uint64_t value, unsigned count) {
    if (count == 0) {
        return value;
    }
    return (value << count) | (value >> (kProgrammerBitWidth - count));
}

std::uint64_t rotate_right_bits(std::uint64_t value, unsigned count) {
    if (count == 0) {
        return value;
    }
    return (value >> count) | (value << (kProgrammerBitWidth - count));
}

int popcount_bits(std::uint64_t value) {
    int count = 0;
    while (value != 0) {
        value &= (value - 1);
        ++count;
    }
    return count;
}

int bit_length_bits(std::uint64_t value) {
    int length = 0;
    while (value != 0) {
        ++length;
        value >>= 1;
    }
    return length;
}

int trailing_zero_count_bits(std::uint64_t value) {
    if (value == 0) {
        return static_cast<int>(kProgrammerBitWidth);
    }
    int count = 0;
    while ((value & 1ULL) == 0ULL) {
        ++count;
        value >>= 1;
    }
    return count;
}

int leading_zero_count_bits(std::uint64_t value) {
    if (value == 0) {
        return static_cast<int>(kProgrammerBitWidth);
    }
    int count = 0;
    std::uint64_t mask = 1ULL << (kProgrammerBitWidth - 1);
    while ((value & mask) == 0ULL) {
        ++count;
        mask >>= 1;
    }
    return count;
}

int parity_bits(std::uint64_t value) {
    return popcount_bits(value) % 2;
}

std::uint64_t reverse_bits(std::uint64_t value) {
    std::uint64_t reversed = 0ULL;
    for (unsigned i = 0; i < kProgrammerBitWidth; ++i) {
        reversed = (reversed << 1) | (value & 1ULL);
        value >>= 1;
    }
    return reversed;
}

std::string factor_integer(long long value) {
    if (value == 0) {
        return "0";
    }
    if (value == 1) {
        return "1";
    }
    if (value == -1) {
        return "-1";
    }

    std::ostringstream out;
    bool first = true;
    if (value < 0) {
        out << "-1";
        first = false;
        value = -value;
    }

    for (long long p = 2; p * p <= value; ++p) {
        int exponent = 0;
        while (value % p == 0) {
            value /= p;
            ++exponent;
        }
        if (exponent == 0) {
            continue;
        }
        if (!first) {
            out << " * ";
        }
        first = false;
        out << p;
        if (exponent > 1) {
            out << "^" << exponent;
        }
    }

    if (value > 1) {
        if (!first) {
            out << " * ";
        }
        out << value;
    }

    return out.str();
}

int digit_value(char ch) {
    // 为 2/8/16 进制统一提供字符到数值的映射。
    if (ch >= '0' && ch <= '9') {
        return ch - '0';
    }
    if (ch >= 'a' && ch <= 'f') {
        return 10 + (ch - 'a');
    }
    if (ch >= 'A' && ch <= 'F') {
        return 10 + (ch - 'A');
    }
    return -1;
}

bool prefixed_base(char prefix, int* base) {
    // 支持 0b / 0o / 0x 三种前缀整数写法。
    if (prefix == 'b' || prefix == 'B') {
        *base = 2;
        return true;
    }
    if (prefix == 'o' || prefix == 'O') {
        *base = 8;
        return true;
    }
    if (prefix == 'x' || prefix == 'X') {
        *base = 16;
        return true;
    }
    return false;
}

long long parse_prefixed_integer_token(const std::string& token) {
    // 这里只处理带前缀的“整数字面量”，不处理小数或科学计数法。
    if (token.size() < 3 || token[0] != '0') {
        throw std::runtime_error("invalid prefixed integer literal");
    }

    int base = 10;
    if (!prefixed_base(token[1], &base)) {
        throw std::runtime_error("invalid prefixed integer literal");
    }

    long long value = 0;
    bool has_digit = false;
    for (std::size_t i = 2; i < token.size(); ++i) {
        const int digit = digit_value(token[i]);
        if (digit < 0 || digit >= base) {
            throw std::runtime_error("invalid digit in prefixed integer literal");
        }
        has_digit = true;
        value = value * static_cast<long long>(base) + static_cast<long long>(digit);
    }

    if (!has_digit) {
        throw std::runtime_error("prefixed integer literal requires digits");
    }

    return value;
}

struct Rational {
    long long numerator = 0;
    long long denominator = 1;

    Rational() = default;

    Rational(long long num, long long den) : numerator(num), denominator(den) {
        normalize();
    }

    void normalize() {
        if (denominator == 0) {
            throw std::runtime_error("division by zero");
        }
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }

        const long long divisor = gcd_ll(numerator, denominator);
        numerator /= divisor;
        denominator /= divisor;
    }

    bool is_integer() const {
        return denominator == 1;
    }

    std::string to_string() const {
        if (denominator == 1) {
            return std::to_string(numerator);
        }
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }
};

Rational operator+(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator * rhs.denominator + rhs.numerator * lhs.denominator,
                    lhs.denominator * rhs.denominator);
}

Rational operator-(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator * rhs.denominator - rhs.numerator * lhs.denominator,
                    lhs.denominator * rhs.denominator);
}

Rational operator*(const Rational& lhs, const Rational& rhs) {
    return Rational(lhs.numerator * rhs.numerator,
                    lhs.denominator * rhs.denominator);
}

Rational operator/(const Rational& lhs, const Rational& rhs) {
    if (rhs.numerator == 0) {
        throw std::runtime_error("division by zero");
    }
    return Rational(lhs.numerator * rhs.denominator,
                    lhs.denominator * rhs.numerator);
}

Rational pow_rational(Rational base, long long exponent) {
    if (exponent == 0) {
        return Rational(1, 1);
    }
    if (base.numerator == 0 && exponent < 0) {
        throw std::runtime_error("zero cannot be raised to a negative power");
    }

    bool negative = exponent < 0;
    unsigned long long power = negative
                                   ? static_cast<unsigned long long>(-exponent)
                                   : static_cast<unsigned long long>(exponent);
    Rational result(1, 1);
    while (power > 0) {
        if (power & 1ULL) {
            result = result * base;
        }
        base = base * base;
        power >>= 1ULL;
    }

    if (negative) {
        return Rational(result.denominator, result.numerator);
    }
    return result;
}

Rational abs_rational(Rational value) {
    if (value.numerator < 0) {
        value.numerator = -value.numerator;
    }
    return value;
}

std::string format_decimal(double value) {
    std::ostringstream out;
    out << std::setprecision(12) << value;
    return out.str();
}

bool try_make_simple_rational(double value,
                              int max_denominator,
                              Rational* rational) {
    if (rational == nullptr || !mymath::isfinite(value)) {
        return false;
    }

    long long numerator = 0;
    long long denominator = 1;
    if (!mymath::approximate_fraction(value,
                                      &numerator,
                                      &denominator,
                                      max_denominator,
                                      1e-10)) {
        return false;
    }

    *rational = Rational(numerator, denominator);
    return true;
}

std::string format_symbolic_number(double value) {
    value = mymath::is_near_zero(value, kDisplayZeroEps) ? 0.0 : value;
    if (mymath::abs(value) > kDisplayIntegerEps &&
        is_integer_double(value, kDisplayIntegerEps)) {
        return std::to_string(round_to_long_long(value));
    }

    Rational rational;
    if (try_make_simple_rational(value, 999, &rational)) {
        return rational.to_string();
    }

    return format_decimal(value);
}

double root_position_tolerance(double value) {
    return 1e-10 * std::max(1.0, mymath::abs(value));
}

double root_function_tolerance(double value) {
    return 1e-10 * std::max(1.0, mymath::abs(value));
}

double root_derivative_step(double value) {
    return 1e-6 * std::max(1.0, mymath::abs(value));
}

constexpr int kPreciseDecimalDivisionDigits = 40;

class PreciseDecimalUnsupported : public std::runtime_error {
public:
    explicit PreciseDecimalUnsupported(const std::string& message)
        : std::runtime_error(message) {}
};

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

struct PreciseDecimal {
    std::string digits = "0";
    int scale = 0;
    bool negative = false;

    void normalize() {
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

    bool is_zero() const {
        return digits == "0";
    }

    std::string to_string() const {
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

    double to_double() const {
        return std::stod(to_string());
    }

    static PreciseDecimal from_digits(std::string raw_digits,
                                      int raw_scale,
                                      bool is_negative) {
        PreciseDecimal value;
        value.digits = std::move(raw_digits);
        value.scale = raw_scale;
        value.negative = is_negative;
        value.normalize();
        return value;
    }

    static PreciseDecimal from_integer_string(const std::string& integer_text,
                                              bool is_negative) {
        return from_digits(integer_text, 0, is_negative);
    }

    static PreciseDecimal from_decimal_literal(const std::string& token) {
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
};

PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs,
                                        const PreciseDecimal& rhs);

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

double rational_to_double(const Rational& value) {
    return static_cast<double>(value.numerator) /
           static_cast<double>(value.denominator);
}

bool is_valid_variable_name(const std::string& name) {
    if (name.empty() ||
        !std::isalpha(static_cast<unsigned char>(name.front()))) {
        return false;
    }

    for (char ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            return false;
        }
    }

    return true;
}

std::string trim_copy(const std::string& text) {
    std::size_t start = 0;
    while (start < text.size() &&
           std::isspace(static_cast<unsigned char>(text[start]))) {
        ++start;
    }

    std::size_t end = text.size();
    while (end > start &&
           std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }

    return text.substr(start, end - start);
}

bool is_identifier_text(const std::string& text) {
    if (text.empty() ||
        !std::isalpha(static_cast<unsigned char>(text.front()))) {
        return false;
    }
    for (char ch : text) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            return false;
        }
    }
    return true;
}

bool is_string_literal(const std::string& text) {
    return text.size() >= 2 && text.front() == '"' && text.back() == '"';
}

std::string decode_escaped_string(const std::string& text) {
    std::string result;
    result.reserve(text.size());
    bool escaping = false;
    for (char ch : text) {
        if (!escaping) {
            if (ch == '\\') {
                escaping = true;
            } else {
                result.push_back(ch);
            }
            continue;
        }

        escaping = false;
        if (ch == 'n') {
            result.push_back('\n');
        } else if (ch == 't') {
            result.push_back('\t');
        } else {
            result.push_back(ch);
        }
    }
    if (escaping) {
        result.push_back('\\');
    }
    return result;
}

std::string parse_string_literal_value(const std::string& text) {
    if (!is_string_literal(text)) {
        throw std::runtime_error("expected string literal");
    }
    return decode_escaped_string(text.substr(1, text.size() - 2));
}

std::string encode_state_field(const std::string& text) {
    std::string encoded;
    encoded.reserve(text.size());
    for (char ch : text) {
        if (ch == '\\') {
            encoded += "\\\\";
        } else if (ch == '\n') {
            encoded += "\\n";
        } else if (ch == '\t') {
            encoded += "\\t";
        } else {
            encoded.push_back(ch);
        }
    }
    return encoded;
}

std::string decode_state_field(const std::string& text) {
    return decode_escaped_string(text);
}

std::vector<std::string> split_top_level_arguments(const std::string& text);
std::string expand_inline_function_commands(Calculator* calculator,
                                            const std::string& expression);

bool split_assignment(const std::string& expression,
                      std::string* lhs,
                      std::string* rhs) {
    // 只在最外层寻找 '='，避免把函数参数或括号内部内容误判为赋值。
    int paren_depth = 0;
    int bracket_depth = 0;
    bool in_string = false;
    bool escaping = false;
    for (std::size_t i = 0; i < expression.size(); ++i) {
        const char ch = expression[i];
        if (in_string) {
            if (escaping) {
                escaping = false;
            } else if (ch == '\\') {
                escaping = true;
            } else if (ch == '"') {
                in_string = false;
            }
            continue;
        }
        if (ch == '"') {
            in_string = true;
            continue;
        }
        if (ch == '(') {
            ++paren_depth;
        } else if (ch == '[') {
            ++bracket_depth;
        } else if (ch == ')') {
            --paren_depth;
        } else if (ch == ']') {
            --bracket_depth;
        } else if (ch == '=' && paren_depth == 0 && bracket_depth == 0) {
            const bool is_comparison =
                (i > 0 && (expression[i - 1] == '!' ||
                           expression[i - 1] == '<' ||
                           expression[i - 1] == '>' ||
                           expression[i - 1] == '=')) ||
                (i + 1 < expression.size() && expression[i + 1] == '=');
            if (is_comparison) {
                continue;
            }
            *lhs = trim_copy(expression.substr(0, i));
            *rhs = trim_copy(expression.substr(i + 1));
            return true;
        }
    }
    return false;
}

bool split_named_call(const std::string& expression,
                      const std::string& name,
                      std::string* inside) {
    // 这个工具函数用于识别“显示型函数”，例如 factor(...) 或 hex(...)。
    // 这类功能不会走普通数值显示，而是会返回格式化字符串。
    const std::string trimmed = trim_copy(expression);
    if (trimmed.size() < name.size() + 2 ||
        trimmed.compare(0, name.size(), name) != 0) {
        return false;
    }

    std::size_t pos = name.size();
    while (pos < trimmed.size() &&
           std::isspace(static_cast<unsigned char>(trimmed[pos]))) {
        ++pos;
    }

    if (pos >= trimmed.size() || trimmed[pos] != '(' || trimmed.back() != ')') {
        return false;
    }

    int depth = 0;
    for (std::size_t i = pos; i < trimmed.size(); ++i) {
        if (trimmed[i] == '(') {
            ++depth;
        } else if (trimmed[i] == ')') {
            --depth;
            if (depth == 0 && i != trimmed.size() - 1) {
                return false;
            }
        }
    }

    if (depth != 0) {
        return false;
    }

    *inside = trim_copy(trimmed.substr(pos + 1, trimmed.size() - pos - 2));
    return true;
}

bool split_named_call_with_arguments(const std::string& expression,
                                     const std::string& name,
                                     std::vector<std::string>* arguments) {
    std::string inside;
    if (!split_named_call(expression, name, &inside)) {
        return false;
    }
    *arguments = split_top_level_arguments(inside);
    return true;
}

std::vector<std::string> split_top_level_arguments(const std::string& text) {
    // 只在最外层逗号处分割参数，允许参数里再嵌套括号表达式。
    std::vector<std::string> arguments;
    int depth = 0;
    int bracket_depth = 0;
    bool in_string = false;
    bool escaping = false;
    std::size_t start = 0;

    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (in_string) {
            if (escaping) {
                escaping = false;
            } else if (ch == '\\') {
                escaping = true;
            } else if (ch == '"') {
                in_string = false;
            }
            continue;
        }
        if (ch == '"') {
            in_string = true;
        } else if (ch == '(') {
            ++depth;
        } else if (ch == '[') {
            ++bracket_depth;
        } else if (ch == ')') {
            --depth;
        } else if (ch == ']') {
            --bracket_depth;
        } else if (ch == ',' && depth == 0 && bracket_depth == 0) {
            arguments.push_back(trim_copy(text.substr(start, i - start)));
            start = i + 1;
        }
    }

    if (!text.empty()) {
        arguments.push_back(trim_copy(text.substr(start)));
    }

    return arguments;
}

bool is_inline_function_command_name(const std::string& name) {
    return name == "diff" ||
           name == "double_integral" ||
           name == "double_integral_cyl" ||
           name == "double_integral_polar" ||
           name == "integral" ||
           name == "limit" ||
           name == "ode" ||
           name == "ode_system" ||
           name == "ode_table" ||
           name == "ode_system_table" ||
           name == "taylor" ||
           name == "triple_integral" ||
           name == "triple_integral_cyl" ||
           name == "triple_integral_sph" ||
           name == "poly_add" ||
           name == "poly_sub" ||
           name == "poly_mul";
}

std::size_t find_matching_paren(const std::string& text, std::size_t open_pos) {
    if (open_pos >= text.size() || text[open_pos] != '(') {
        return std::string::npos;
    }

    int depth = 0;
    bool in_string = false;
    bool escaping = false;
    for (std::size_t i = open_pos; i < text.size(); ++i) {
        const char ch = text[i];
        if (in_string) {
            if (escaping) {
                escaping = false;
            } else if (ch == '\\') {
                escaping = true;
            } else if (ch == '"') {
                in_string = false;
            }
            continue;
        }
        if (ch == '"') {
            in_string = true;
            continue;
        }
        if (ch == '(') {
            ++depth;
        } else if (ch == ')') {
            --depth;
            if (depth == 0) {
                return i;
            }
        }
    }

    return std::string::npos;
}

std::string expand_inline_function_commands(Calculator* calculator,
                                            const std::string& expression) {
    std::string expanded;
    expanded.reserve(expression.size());

    for (std::size_t i = 0; i < expression.size();) {
        const char ch = expression[i];
        if (!std::isalpha(static_cast<unsigned char>(ch))) {
            expanded.push_back(ch);
            ++i;
            continue;
        }

        std::size_t name_end = i;
        while (name_end < expression.size()) {
            const char name_ch = expression[name_end];
            if (std::isalnum(static_cast<unsigned char>(name_ch)) || name_ch == '_') {
                ++name_end;
            } else {
                break;
            }
        }

        const std::string name = expression.substr(i, name_end - i);
        if (!is_inline_function_command_name(name)) {
            expanded.append(expression, i, name_end - i);
            i = name_end;
            continue;
        }

        std::size_t open_pos = name_end;
        while (open_pos < expression.size() &&
               std::isspace(static_cast<unsigned char>(expression[open_pos]))) {
            ++open_pos;
        }
        if (open_pos >= expression.size() || expression[open_pos] != '(') {
            expanded.append(expression, i, name_end - i);
            i = name_end;
            continue;
        }

        const std::size_t close_pos = find_matching_paren(expression, open_pos);
        if (close_pos == std::string::npos) {
            expanded.append(expression, i, name_end - i);
            i = name_end;
            continue;
        }

        const std::string inner =
            expand_inline_function_commands(calculator,
                                            expression.substr(open_pos + 1,
                                                              close_pos - open_pos - 1));
        const std::string rebuilt = name + "(" + inner + ")";
        std::string command_output;
        if (calculator->try_process_function_command(rebuilt, &command_output)) {
            expanded += "(" + command_output + ")";
        } else {
            expanded += rebuilt;
        }
        i = close_pos + 1;
    }

    return expanded;
}

struct StoredValue {
    // exact=true 时优先使用 rational 字段显示和参与 exact mode 计算；
    // 否则使用 decimal 字段。
    bool is_matrix = false;
    bool is_string = false;
    bool has_symbolic_text = false;
    bool has_precise_decimal_text = false;
    bool exact = false;
    Rational rational;
    double decimal = 0.0;
    std::string string_value;
    std::string symbolic_text;
    std::string precise_decimal_text;
    matrix::Matrix matrix;
};

struct CustomFunction {
    std::string parameter_name;
    std::string expression;
};

struct ScriptFunction {
    std::vector<std::string> parameter_names;
    std::shared_ptr<const script::BlockStatement> body;
};

using HasScriptFunctionCallback = std::function<bool(const std::string&)>;
using InvokeScriptFunctionDecimalCallback =
    std::function<double(const std::string&, const std::vector<double>&)>;

std::string rational_to_precise_decimal_text(const Rational& value) {
    PreciseDecimal numerator = PreciseDecimal::from_integer_string(
        std::to_string(value.numerator < 0 ? -value.numerator : value.numerator),
        value.numerator < 0);
    const PreciseDecimal denominator =
        PreciseDecimal::from_integer_string(std::to_string(value.denominator), false);
    return divide_precise_decimal(numerator, denominator).to_string();
}

std::string stored_value_precise_decimal_text(const StoredValue& value) {
    if (value.exact) {
        return rational_to_precise_decimal_text(value.rational);
    }
    if (value.has_precise_decimal_text) {
        return value.precise_decimal_text;
    }
    return format_decimal(normalize_display_decimal(value.decimal));
}

class PreciseDecimalParser {
public:
    PreciseDecimalParser(std::string source,
                         const std::map<std::string, StoredValue>* variables)
        : source_(std::move(source)),
          variables_(variables) {}

    PreciseDecimal parse() {
        PreciseDecimal value = parse_expression();
        skip_spaces();
        if (!is_at_end()) {
            throw PreciseDecimalUnsupported("unsupported token");
        }
        return value;
    }

private:
    PreciseDecimal parse_expression() {
        PreciseDecimal value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = add_precise_decimal(value, parse_term());
            } else if (match('-')) {
                value = subtract_precise_decimal(value, parse_term());
            } else {
                break;
            }
        }
        return value;
    }

    PreciseDecimal parse_term() {
        PreciseDecimal value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = multiply_precise_decimal(value, parse_unary());
            } else if (match('/')) {
                value = divide_precise_decimal(value, parse_unary());
            } else {
                break;
            }
        }
        return value;
    }

    PreciseDecimal parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            PreciseDecimal value = parse_unary();
            if (!value.is_zero()) {
                value.negative = !value.negative;
            }
            return value;
        }
        return parse_primary();
    }

    PreciseDecimal parse_primary() {
        skip_spaces();
        if (match('(')) {
            PreciseDecimal value = parse_expression();
            skip_spaces();
            if (!peek(')')) {
                throw PreciseDecimalUnsupported(
                    "only parenthesized four-operator expressions are supported");
            }
            expect(')');
            return value;
        }

        if (peek_is_alpha()) {
            const std::string name = parse_identifier();
            skip_spaces();
            if (peek('(')) {
                throw PreciseDecimalUnsupported("function calls are not supported");
            }
            return lookup_variable(name);
        }

        return parse_number();
    }

    PreciseDecimal parse_number() {
        skip_spaces();

        if (!is_at_end() &&
            source_[pos_] == '0' &&
            pos_ + 1 < source_.size()) {
            int base = 10;
            if (prefixed_base(source_[pos_ + 1], &base)) {
                const std::size_t start = pos_;
                pos_ += 2;
                while (!is_at_end()) {
                    const int digit = digit_value(source_[pos_]);
                    if (digit < 0 || digit >= base) {
                        break;
                    }
                    ++pos_;
                }
                const long long integer_value = parse_prefixed_integer_token(
                    source_.substr(start, pos_ - start));
                return PreciseDecimal::from_integer_string(
                    std::to_string(integer_value < 0 ? -integer_value : integer_value),
                    integer_value < 0);
            }
        }

        const std::size_t start = pos_;
        bool has_digit = false;
        bool seen_dot = false;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                has_digit = true;
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else {
                break;
            }
        }

        if (!is_at_end() && (source_[pos_] == 'e' || source_[pos_] == 'E')) {
            const std::size_t exponent_pos = pos_;
            ++pos_;
            if (!is_at_end() && (source_[pos_] == '+' || source_[pos_] == '-')) {
                ++pos_;
            }
            const std::size_t exponent_digits = pos_;
            while (!is_at_end() &&
                   std::isdigit(static_cast<unsigned char>(source_[pos_]))) {
                ++pos_;
            }
            if (exponent_digits == pos_) {
                pos_ = exponent_pos;
            }
        }

        if (!has_digit) {
            throw std::runtime_error("expected number");
        }

        return PreciseDecimal::from_decimal_literal(source_.substr(start, pos_ - start));
    }

    PreciseDecimal lookup_variable(const std::string& name) const {
        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            throw PreciseDecimalUnsupported("unknown variable for precise parsing");
        }
        if (it->second.is_matrix || it->second.is_string || it->second.has_symbolic_text) {
            throw PreciseDecimalUnsupported("unsupported variable type for precise parsing");
        }
        return PreciseDecimal::from_decimal_literal(
            stored_value_precise_decimal_text(it->second));
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
                ++pos_;
            } else {
                break;
            }
        }
        return source_.substr(start, pos_ - start);
    }

    bool peek_is_alpha() const {
        return !is_at_end() &&
               std::isalpha(static_cast<unsigned char>(source_[pos_]));
    }

    bool peek(char expected) const {
        return !is_at_end() && source_[pos_] == expected;
    }

    void skip_spaces() {
        while (!is_at_end() &&
               std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }

    bool match(char expected) {
        if (is_at_end() || source_[pos_] != expected) {
            return false;
        }
        ++pos_;
        return true;
    }

    void expect(char expected) {
        if (!match(expected)) {
            throw std::runtime_error(std::string("expected '") + expected + "'");
        }
    }

    bool is_at_end() const {
        return pos_ >= source_.size();
    }

    std::string source_;
    std::size_t pos_ = 0;
    const std::map<std::string, StoredValue>* variables_;
};

std::string format_stored_value(const StoredValue& value, bool symbolic_constants_mode) {
    if (value.is_matrix) {
        return value.matrix.to_string();
    }
    if (value.is_string) {
        std::ostringstream out;
        out << '"';
        for (char ch : value.string_value) {
            if (ch == '\\' || ch == '"') {
                out << '\\';
            } else if (ch == '\n') {
                out << "\\n";
                continue;
            } else if (ch == '\t') {
                out << "\\t";
                continue;
            }
            out << ch;
        }
        out << '"';
        return out.str();
    }
    if (symbolic_constants_mode && value.has_symbolic_text) {
        return value.symbolic_text;
    }
    if (value.exact) {
        return value.rational.to_string();
    }
    const double normalized_decimal = normalize_display_decimal(value.decimal);
    if (value.has_precise_decimal_text) {
        if (normalized_decimal != value.decimal) {
            return format_decimal(normalized_decimal);
        }
        return value.precise_decimal_text;
    }
    return format_decimal(normalized_decimal);
}

std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode) {
    if (value.is_string) {
        return value.string_value;
    }
    return format_stored_value(value, symbolic_constants_mode);
}

std::string format_symbolic_scalar(double value) {
    return format_symbolic_number(value);
}

double factorial_int(int n) {
    double result = 1.0;
    for (int i = 2; i <= n; ++i) {
        result *= static_cast<double>(i);
    }
    return result;
}

double factorial_value(long long n) {
    if (n < 0) {
        throw std::runtime_error("factorial only accepts non-negative integers");
    }
    if (n > 170) {
        throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    }
    double result = 1.0;
    for (long long i = 2; i <= n; ++i) {
        result *= static_cast<double>(i);
    }
    return result;
}

Rational factorial_rational(long long n) {
    if (n < 0) {
        throw std::runtime_error("factorial only accepts non-negative integers");
    }
    if (n > 170) {
        throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    }
    Rational result(1, 1);
    for (long long i = 2; i <= n; ++i) {
        result = result * Rational(i, 1);
    }
    return result;
}

double combination_value(long long n, long long r) {
    if (n < 0 || r < 0) {
        throw std::runtime_error("nCr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nCr requires r <= n");
    }
    const long long effective_r = r < (n - r) ? r : (n - r);
    double result = 1.0;
    for (long long i = 1; i <= effective_r; ++i) {
        result *= static_cast<double>(n - effective_r + i);
        result /= static_cast<double>(i);
    }
    return result;
}

Rational combination_rational(long long n, long long r) {
    if (n < 0 || r < 0) {
        throw std::runtime_error("nCr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nCr requires r <= n");
    }
    const long long effective_r = r < (n - r) ? r : (n - r);
    Rational result(1, 1);
    for (long long i = 1; i <= effective_r; ++i) {
        result = result * Rational(n - effective_r + i, i);
    }
    return result;
}

double permutation_value(long long n, long long r) {
    if (n < 0 || r < 0) {
        throw std::runtime_error("nPr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nPr requires r <= n");
    }
    double result = 1.0;
    for (long long i = 0; i < r; ++i) {
        result *= static_cast<double>(n - i);
    }
    return result;
}

Rational permutation_rational(long long n, long long r) {
    if (n < 0 || r < 0) {
        throw std::runtime_error("nPr only accepts non-negative integers");
    }
    if (r > n) {
        throw std::runtime_error("nPr requires r <= n");
    }
    Rational result(1, 1);
    for (long long i = 0; i < r; ++i) {
        result = result * Rational(n - i, 1);
    }
    return result;
}

std::string taylor_series_to_string(const std::vector<double>& coefficients,
                                    const std::string& variable_name,
                                    double center) {
    bool zero_center = mymath::is_near_zero(center, 1e-10);
    std::ostringstream out;
    bool first = true;

    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const double coefficient = coefficients[i];
        if (mymath::is_near_zero(coefficient, 1e-10)) {
            continue;
        }

        const bool negative = coefficient < 0.0;
        const double abs_value = negative ? -coefficient : coefficient;
        std::string term;

        if (i == 0) {
            term = format_symbolic_scalar(abs_value);
        } else {
            if (!mymath::is_near_zero(abs_value - 1.0, 1e-10)) {
                term += format_symbolic_scalar(abs_value) + " * ";
            }

            if (zero_center) {
                term += variable_name;
            } else {
                term += "(" + variable_name;
                term += center < 0.0 ? " + " : " - ";
                term += format_symbolic_scalar(mymath::abs(center)) + ")";
            }

            if (i > 1) {
                term += " ^ " + std::to_string(i);
            }
        }

        if (first) {
            out << (negative ? "-" : "") << term;
            first = false;
        } else {
            out << (negative ? " - " : " + ") << term;
        }
    }

    return first ? "0" : out.str();
}

std::string shifted_series_base(const std::string& variable_name, double center) {
    if (mymath::is_near_zero(center, 1e-10)) {
        return variable_name;
    }

    std::string base = "(" + variable_name;
    base += center < 0.0 ? " + " : " - ";
    base += format_symbolic_scalar(mymath::abs(center)) + ")";
    return base;
}

std::string generalized_series_to_string(const std::vector<double>& coefficients,
                                         const std::string& variable_name,
                                         double center,
                                         int denominator) {
    if (denominator <= 0) {
        throw std::runtime_error("series denominator must be positive");
    }

    const std::string base = shifted_series_base(variable_name, center);
    std::ostringstream out;
    bool first = true;

    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const double coefficient = coefficients[i];
        if (mymath::is_near_zero(coefficient, 1e-10)) {
            continue;
        }

        const bool negative = coefficient < 0.0;
        const double abs_value = negative ? -coefficient : coefficient;
        std::string term;

        if (i == 0) {
            term = format_symbolic_scalar(abs_value);
        } else {
            if (!mymath::is_near_zero(abs_value - 1.0, 1e-10)) {
                term += format_symbolic_scalar(abs_value) + " * ";
            }

            term += base;
            if (denominator == 1) {
                if (i > 1) {
                    term += " ^ " + std::to_string(i);
                }
            } else if (i % static_cast<std::size_t>(denominator) == 0) {
                const std::size_t exponent = i / static_cast<std::size_t>(denominator);
                if (exponent > 1) {
                    term += " ^ " + std::to_string(exponent);
                }
            } else {
                term += " ^ (" + std::to_string(i) + " / " +
                        std::to_string(denominator) + ")";
            }
        }

        if (first) {
            out << (negative ? "-" : "") << term;
            first = false;
        } else {
            out << (negative ? " - " : " + ") << term;
        }
    }

    return first ? "0" : out.str();
}

std::vector<double> solve_dense_linear_system(std::vector<std::vector<double>> matrix,
                                              std::vector<double> rhs,
                                              const std::string& context) {
    const std::size_t n = matrix.size();
    if (rhs.size() != n) {
        throw std::runtime_error(context + " linear system dimension mismatch");
    }
    for (const auto& row : matrix) {
        if (row.size() != n) {
            throw std::runtime_error(context + " linear system must be square");
        }
    }

    std::vector<std::vector<long double>> high_precision_matrix(
        n, std::vector<long double>(n, 0.0L));
    std::vector<long double> high_precision_rhs(n, 0.0L);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            high_precision_matrix[row][col] =
                static_cast<long double>(matrix[row][col]);
        }
        high_precision_rhs[row] = static_cast<long double>(rhs[row]);
    }

    for (std::size_t pivot = 0; pivot < n; ++pivot) {
        std::size_t best_row = pivot;
        long double best_value = mymath::abs_long_double(high_precision_matrix[pivot][pivot]);
        for (std::size_t row = pivot + 1; row < n; ++row) {
            const long double current = mymath::abs_long_double(high_precision_matrix[row][pivot]);
            if (current > best_value) {
                best_value = current;
                best_row = row;
            }
        }

        if (mymath::abs_long_double(best_value) <= 1e-12L) {
            throw std::runtime_error(context + " system is singular");
        }
        if (best_row != pivot) {
            std::swap(high_precision_matrix[pivot], high_precision_matrix[best_row]);
            std::swap(high_precision_rhs[pivot], high_precision_rhs[best_row]);
        }

        const long double pivot_value = high_precision_matrix[pivot][pivot];
        for (std::size_t col = pivot; col < n; ++col) {
            high_precision_matrix[pivot][col] /= pivot_value;
        }
        high_precision_rhs[pivot] /= pivot_value;

        for (std::size_t row = 0; row < n; ++row) {
            if (row == pivot) {
                continue;
            }
            const long double factor = high_precision_matrix[row][pivot];
            if (mymath::abs_long_double(factor) <= 1e-12L) {
                continue;
            }
            for (std::size_t col = pivot; col < n; ++col) {
                high_precision_matrix[row][col] -=
                    factor * high_precision_matrix[pivot][col];
            }
            high_precision_rhs[row] -= factor * high_precision_rhs[pivot];
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        rhs[i] = static_cast<double>(high_precision_rhs[i]);
    }
    return rhs;
}

bool is_reserved_function_name(const std::string& name) {
    static const std::vector<std::string> names = {
        "abs", "acos", "acosh", "acot", "acsc", "and", "asec", "asin",
        "asinh", "atan", "atanh", "base", "beta", "bin", "binom",
        "bessel", "bitlen", "c2f", "cbrt", "cdf_normal",
        "ceil", "cholesky", "clamp", "complex", "cond", "conj", "corr", "cos",
        "cos_deg", "cosh", "cot", "cov", "deg", "deg2rad", "diag",
        "diff", "double_integral",
        "double_integral_cyl", "double_integral_polar", "exp", "exp2", "extrema",
        "f2c", "factor", "factorial", "fahrenheit", "fib", "fixed_point",
        "floor", "gamma", "gcd", "get", "hadamard", "hessenberg",
        "hex", "identity", "imag", "integral", "inverse", "is_prime",
        "kelvin", "kron", "lagrange", "lcm", "least_squares",
        "linear_regression", "ln", "log", "log10", "log2", "lu_l", "lu_u", "mat", "max",
        "mean", "median", "min", "mod", "mode", "nCr", "nPr", "next_prime",
        "norm", "not", "null", "oct", "ode", "ode_table", "ode_system", "ode_system_table", "or", "outer",
        "parity", "pdf_normal", "percentile", "pinv", "polar", "poly_add", "poly_compose",
        "poly_deriv", "poly_div", "poly_eval", "poly_fit", "poly_gcd",
        "poly_integ", "poly_mul", "poly_sub", "polynomial_fit", "pow",
        "qr_q", "qr_r", "rad", "rad2deg", "rand", "randint", "randn",
        "rank", "rat", "real", "reshape", "resize", "rref", "rol", "root", "roots",
        "ror", "round",
        "schur", "sec", "secant", "set", "shl", "shr", "sign", "sin",
        "sin_deg", "sinh", "solve", "spline", "sqrt", "std", "sum",
        "svd", "svd_s", "svd_u", "svd_vt", "tan", "tanh", "taylor",
        "trace", "transpose", "triple_integral", "triple_integral_cyl",
        "triple_integral_sph", "trunc", "avg", "var", "vec", "xor", "zeta",
        "celsius", "delta", "heaviside", "impulse", "step",
        "quartile", "popcount", "ctz", "clz", "reverse_bits",
        "fourier", "ifourier", "inverse_fourier",
        "laplace", "ilaplace", "inverse_laplace",
        "ztrans", "iztrans", "z_transform", "inverse_z",
        "dft", "idft", "fft", "ifft", "conv", "convolve",
        "pade", "puiseux", "series_sum", "summation",
        "eig", "eigvals", "eigvecs",
        "lp_max", "lp_min", "ilp_max", "ilp_min",
        "milp_max", "milp_min", "bip_max", "bip_min", "binary_max", "binary_min"
    };

    for (const std::string& builtin : names) {
        if (builtin == name) {
            return true;
        }
    }
    return name == "pi" || name == "e";
}

bool split_function_definition(const std::string& expression,
                               std::string* function_name,
                               std::string* parameter_name,
                               std::string* body) {
    std::string lhs;
    std::string rhs;
    if (!split_assignment(expression, &lhs, &rhs)) {
        return false;
    }

    const std::size_t open = lhs.find('(');
    const std::size_t close = lhs.rfind(')');
    if (open == std::string::npos || close == std::string::npos || close <= open) {
        return false;
    }

    const std::string before = trim_copy(lhs.substr(0, open));
    const std::string inside = trim_copy(lhs.substr(open + 1, close - open - 1));
    const std::string after = trim_copy(lhs.substr(close + 1));
    if (!after.empty()) {
        return false;
    }
    if (!is_valid_variable_name(before) || !is_valid_variable_name(inside) ||
        rhs.empty()) {
        return false;
    }

    *function_name = before;
    *parameter_name = inside;
    *body = rhs;
    return true;
}

bool try_evaluate_matrix_expression(const std::string& expression,
                                    const std::map<std::string, StoredValue>* variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const HasScriptFunctionCallback& has_script_function,
                                    const InvokeScriptFunctionDecimalCallback& invoke_script_function,
                                    matrix::Value* value);

// DecimalParser 使用浮点数进行常规求值，支持全部数学函数。
class DecimalParser {
public:
    DecimalParser(std::string source,
                  const std::map<std::string, StoredValue>* variables,
                  const std::map<std::string, CustomFunction>* functions,
                  HasScriptFunctionCallback has_script_function = {},
                  InvokeScriptFunctionDecimalCallback invoke_script_function = {})
        : source_(std::move(source)),
          variables_(variables),
          functions_(functions),
          has_script_function_(std::move(has_script_function)),
          invoke_script_function_(std::move(invoke_script_function)) {}

    double parse() {
        double value = parse_comparison();
        skip_spaces();
        if (!is_at_end()) {
            throw std::runtime_error("unexpected token near: " + source_.substr(pos_, 1));
        }
        return value;
    }

private:
    double parse_comparison() {
        double value = parse_expression();
        while (true) {
            skip_spaces();
            if (match_string("==")) {
                value = mymath::is_near_zero(value - parse_expression(), 1e-10) ? 1.0 : 0.0;
            } else if (match_string("!=")) {
                value = mymath::is_near_zero(value - parse_expression(), 1e-10) ? 0.0 : 1.0;
            } else if (match_string("<=")) {
                value = value <= parse_expression() ? 1.0 : 0.0;
            } else if (match_string(">=")) {
                value = value >= parse_expression() ? 1.0 : 0.0;
            } else if (match('<')) {
                value = value < parse_expression() ? 1.0 : 0.0;
            } else if (match('>')) {
                value = value > parse_expression() ? 1.0 : 0.0;
            } else {
                break;
            }
        }
        return value;
    }

    double parse_expression() {
        double value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value += parse_term();
            } else if (match('-')) {
                value -= parse_term();
            } else {
                break;
            }
        }
        return value;
    }

    double parse_term() {
        double value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value *= parse_unary();
            } else if (match('/')) {
                const double divisor = parse_unary();
                if (mymath::is_near_zero(divisor)) {
                    throw std::runtime_error("division by zero");
                }
                value /= divisor;
            } else {
                break;
            }
        }
        return value;
    }

    double parse_power() {
        double value = parse_primary();
        skip_spaces();
        if (match('^')) {
            const double exponent = parse_unary();
            return mymath::pow(value, exponent);
        }
        return value;
    }

    double parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            return -parse_unary();
        }
        return parse_power();
    }

    double parse_primary() {
        skip_spaces();
        if (match('(')) {
            const double value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        if (peek_is_alpha()) {
            const std::string name = parse_identifier();

            skip_spaces();
            if (!peek('(')) {
                return lookup_variable(name);
            }

            skip_spaces();
            expect('(');
            const std::vector<double> arguments = parse_argument_list();
            expect(')');
            return apply_function(name, arguments);
        }

        return parse_number();
    }

    std::vector<double> parse_argument_list() {
        std::vector<double> arguments;
        skip_spaces();
        if (peek(')')) {
            return arguments;
        }

        while (true) {
            arguments.push_back(parse_expression());
            skip_spaces();
            if (!match(',')) {
                break;
            }
        }
        return arguments;
    }

    double parse_number() {
        skip_spaces();

        // 先尝试解析 0b/0o/0x 这类前缀整数。
        if (!is_at_end() &&
            source_[pos_] == '0' &&
            pos_ + 1 < source_.size()) {
            int base = 10;
            if (prefixed_base(source_[pos_ + 1], &base)) {
                const std::size_t start = pos_;
                pos_ += 2;
                while (!is_at_end()) {
                    const int digit = digit_value(source_[pos_]);
                    if (digit < 0 || digit >= base) {
                        break;
                    }
                    ++pos_;
                }
                return static_cast<double>(
                    parse_prefixed_integer_token(source_.substr(start, pos_ - start)));
            }
        }

        const std::size_t start = pos_;
        bool has_digit = false;
        bool seen_dot = false;

        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                has_digit = true;
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else {
                break;
            }
        }

        if (!is_at_end() && (source_[pos_] == 'e' || source_[pos_] == 'E')) {
            const std::size_t exponent_pos = pos_;
            ++pos_;
            if (!is_at_end() && (source_[pos_] == '+' || source_[pos_] == '-')) {
                ++pos_;
            }
            const std::size_t exponent_digits = pos_;
            while (!is_at_end() &&
                   std::isdigit(static_cast<unsigned char>(source_[pos_]))) {
                ++pos_;
            }
            if (exponent_digits == pos_) {
                pos_ = exponent_pos;
            }
        }

        if (!has_digit) {
            throw std::runtime_error("expected number");
        }

        return parse_decimal(source_.substr(start, pos_ - start));
    }

    static double parse_decimal(const std::string& token) {
        return std::stod(token);
    }

    double apply_function(const std::string& name, const std::vector<double>& arguments) {
        // 这里是“普通浮点路径”的函数分发中心。
        // 需要字符串结果的功能，例如 factor/bin/hex，不从这里返回。
        if (name == "pow") {
            return apply_pow(arguments);
        }
        if (name == "root") {
            return apply_root(arguments);
        }
        if (name == "and") {
            return apply_and(arguments);
        }
        if (name == "or") {
            return apply_or(arguments);
        }
        if (name == "xor") {
            return apply_xor(arguments);
        }
        if (name == "shl") {
            return apply_shl(arguments);
        }
        if (name == "shr") {
            return apply_shr(arguments);
        }
        if (name == "rol") {
            return apply_rol(arguments);
        }
        if (name == "ror") {
            return apply_ror(arguments);
        }
        if (name == "gcd") {
            return apply_gcd(arguments);
        }
        if (name == "lcm") {
            return apply_lcm(arguments);
        }
        if (name == "mod") {
            return apply_mod(arguments);
        }
        if (name == "min") {
            return apply_min(arguments);
        }
        if (name == "max") {
            return apply_max(arguments);
        }
        if (name == "clamp") {
            return apply_clamp(arguments);
        }
        if (name == "log") {
            return apply_log(arguments);
        }
        if (name == "sum") {
            return apply_sum(arguments);
        }
        if (name == "mean") {
            return apply_mean(arguments);
        }
        if (name == "avg") {
            return apply_avg(arguments);
        }
        if (name == "median") {
            return apply_median(arguments);
        }
        if (name == "mode") {
            return apply_mode(arguments);
        }
        if (name == "var") {
            return apply_variance(arguments);
        }
        if (name == "std") {
            return apply_stddev(arguments);
        }
        if (name == "percentile") {
            return apply_percentile(arguments);
        }
        if (name == "quartile") {
            return apply_quartile(arguments);
        }
        if (name == "factorial") {
            return apply_factorial(arguments);
        }
        if (name == "nCr") {
            return apply_ncr(arguments);
        }
        if (name == "binom") {
            return apply_ncr(arguments);
        }
        if (name == "nPr") {
            return apply_npr(arguments);
        }
        if (name == "fib") {
            return apply_fib(arguments);
        }
        if (name == "is_prime") {
            return apply_is_prime(arguments);
        }
        if (name == "next_prime") {
            return apply_next_prime(arguments);
        }
        if (name == "rand") {
            return apply_rand(arguments);
        }
        if (name == "randn") {
            return apply_randn(arguments);
        }
        if (name == "randint") {
            return apply_randint(arguments);
        }
        if (name == "beta") {
            return apply_beta(arguments);
        }
        if (name == "zeta") {
            return apply_zeta(arguments);
        }
        if (name == "bessel") {
            return apply_bessel(arguments);
        }
        if (name == "pdf_normal") {
            return apply_pdf_normal(arguments);
        }
        if (name == "cdf_normal") {
            return apply_cdf_normal(arguments);
        }
        if (name == "popcount") {
            return apply_popcount(arguments);
        }
        if (name == "bitlen") {
            return apply_bitlen(arguments);
        }
        if (name == "ctz") {
            return apply_ctz(arguments);
        }
        if (name == "clz") {
            return apply_clz(arguments);
        }
        if (name == "parity") {
            return apply_parity(arguments);
        }
        if (name == "reverse_bits") {
            return apply_reverse_bits(arguments);
        }

        const auto function_it = functions_->find(name);
        if (function_it != functions_->end()) {
            if (arguments.size() != 1) {
                throw std::runtime_error("custom function " + name +
                                         " expects exactly one argument");
            }

            std::map<std::string, StoredValue> scoped_variables = *variables_;
            StoredValue parameter_value;
            parameter_value.decimal = arguments[0];
            parameter_value.exact = false;
            scoped_variables[function_it->second.parameter_name] = parameter_value;

            DecimalParser nested_parser(function_it->second.expression,
                                        &scoped_variables,
                                        functions_,
                                        has_script_function_,
                                        invoke_script_function_);
            return nested_parser.parse();
        }

        if (has_script_function_ && has_script_function_(name)) {
            return invoke_script_function_(name, arguments);
        }

        if (arguments.size() != 1) {
            throw std::runtime_error("function " + name + " expects exactly one argument");
        }

        const double argument = arguments[0];
        if (name == "abs") {
            return mymath::abs(argument);
        }
        if (name == "not") {
            // 位运算函数只接受整数，浮点路径里也保持这个约束。
            if (!is_integer_double(argument)) {
                throw std::runtime_error("not only accepts integers");
            }
            return static_cast<double>(~round_to_long_long(argument));
        }
        if (name == "sign") {
            if (mymath::is_near_zero(argument)) {
                return 0.0;
            }
            return argument > 0.0 ? 1.0 : -1.0;
        }
        if (name == "floor") {
            return static_cast<double>(floor_to_long_long(argument));
        }
        if (name == "ceil") {
            return static_cast<double>(ceil_to_long_long(argument));
        }
        if (name == "round") {
            return static_cast<double>(round_to_long_long(argument));
        }
        if (name == "trunc") {
            return static_cast<double>(trunc_to_long_long(argument));
        }
        if (name == "cbrt") {
            return mymath::cbrt(argument);
        }
        if (name == "asinh") {
            return mymath::asinh(argument);
        }
        if (name == "acosh") {
            return mymath::acosh(argument);
        }
        if (name == "atanh") {
            return mymath::atanh(argument);
        }
        if (name == "sinh") {
            return mymath::sinh(argument);
        }
        if (name == "cosh") {
            return mymath::cosh(argument);
        }
        if (name == "tanh") {
            return mymath::tanh(argument);
        }
        if (name == "sec") {
            return mymath::sec(argument);
        }
        if (name == "csc") {
            return mymath::csc(argument);
        }
        if (name == "cot") {
            return mymath::cot(argument);
        }
        if (name == "sin") {
            return mymath::sin(argument);
        }
        if (name == "cos") {
            return mymath::cos(argument);
        }
        if (name == "tan") {
            return mymath::tan(argument);
        }
        if (name == "atan") {
            return mymath::atan(argument);
        }
        if (name == "asin") {
            return mymath::asin(argument);
        }
        if (name == "acos") {
            return mymath::acos(argument);
        }
        if (name == "asec") {
            return mymath::asec(argument);
        }
        if (name == "acsc") {
            return mymath::acsc(argument);
        }
        if (name == "acot") {
            return mymath::acot(argument);
        }
        if (name == "ln") {
            return mymath::ln(argument);
        }
        if (name == "log2") {
            return mymath::ln(argument) / mymath::ln(2.0);
        }
        if (name == "log10") {
            return mymath::log10(argument);
        }
        if (name == "exp") {
            return mymath::exp(argument);
        }
        if (name == "exp2") {
            return mymath::exp(argument * mymath::ln(2.0));
        }
        if (name == "gamma") {
            return mymath::gamma(argument);
        }
        if (name == "erf") {
            return mymath::erf(argument);
        }
        if (name == "erfc") {
            return mymath::erfc(argument);
        }
        if (name == "sqrt") {
            return mymath::sqrt(argument);
        }
        if (name == "step" || name == "u" || name == "heaviside") {
            return argument >= 0.0 ? 1.0 : 0.0;
        }
        if (name == "delta" || name == "impulse") {
            return mymath::is_near_zero(argument, 1e-10) ? 1.0 : 0.0;
        }
        if (name == "deg") {
            return radians_to_degrees(argument);
        }
        if (name == "rad") {
            return degrees_to_radians(argument);
        }
        if (name == "deg2rad") {
            return degrees_to_radians(argument);
        }
        if (name == "rad2deg") {
            return radians_to_degrees(argument);
        }
        if (name == "sin_deg") {
            return mymath::sin(degrees_to_radians(argument));
        }
        if (name == "cos_deg") {
            return mymath::cos(degrees_to_radians(argument));
        }
        if (name == "celsius") {
            return fahrenheit_to_celsius(argument);
        }
        if (name == "fahrenheit") {
            return celsius_to_fahrenheit(argument);
        }
        if (name == "kelvin") {
            return argument + 273.15;
        }
        if (name == "c2f") {
            return celsius_to_fahrenheit(argument);
        }
        if (name == "f2c") {
            return fahrenheit_to_celsius(argument);
        }

        throw std::runtime_error("unknown function: " + name);
    }

    static double apply_gcd(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("gcd expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("gcd only accepts integers");
        }
        return static_cast<double>(
            gcd_ll(round_to_long_long(arguments[0]), round_to_long_long(arguments[1])));
    }

    static double apply_pow(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("pow expects exactly two arguments");
        }
        return mymath::pow(arguments[0], arguments[1]);
    }

    static double apply_root(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("root expects exactly two arguments");
        }
        return mymath::root(arguments[0], arguments[1]);
    }

    static double apply_lcm(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("lcm expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("lcm only accepts integers");
        }
        return static_cast<double>(
            lcm_ll(round_to_long_long(arguments[0]), round_to_long_long(arguments[1])));
    }

    static double apply_mod(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("mod expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("mod only accepts integers");
        }

        const long long lhs = round_to_long_long(arguments[0]);
        const long long rhs = round_to_long_long(arguments[1]);
        if (rhs == 0) {
            throw std::runtime_error("mod divisor cannot be zero");
        }
        return static_cast<double>(lhs % rhs);
    }

    static double apply_min(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("min expects exactly two arguments");
        }
        return arguments[0] < arguments[1] ? arguments[0] : arguments[1];
    }

    static double apply_max(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("max expects exactly two arguments");
        }
        return arguments[0] > arguments[1] ? arguments[0] : arguments[1];
    }

    static double apply_clamp(const std::vector<double>& arguments) {
        if (arguments.size() != 3) {
            throw std::runtime_error("clamp expects exactly three arguments");
        }
        double lower = arguments[1];
        double upper = arguments[2];
        if (lower > upper) {
            std::swap(lower, upper);
        }
        if (arguments[0] < lower) {
            return lower;
        }
        if (arguments[0] > upper) {
            return upper;
        }
        return arguments[0];
    }

    static double apply_log(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("log expects exactly two arguments");
        }
        if (mymath::is_near_zero(arguments[1] - 1.0)) {
            throw std::runtime_error("log base cannot be 1");
        }
        return mymath::ln(arguments[0]) / mymath::ln(arguments[1]);
    }

    static double apply_sum(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("sum expects at least one argument");
        }
        double total = 0.0;
        for (double value : arguments) {
            total += value;
        }
        return total;
    }

    static double apply_avg(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("avg expects at least one argument");
        }
        return apply_sum(arguments) / static_cast<double>(arguments.size());
    }

    static double apply_mean(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("mean expects at least one argument");
        }
        return apply_sum(arguments) / static_cast<double>(arguments.size());
    }

    static double apply_median(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("median expects at least one argument");
        }

        std::vector<double> sorted = arguments;
        std::sort(sorted.begin(), sorted.end());
        const std::size_t middle = sorted.size() / 2;
        if (sorted.size() % 2 == 1) {
            return sorted[middle];
        }
        return (sorted[middle - 1] + sorted[middle]) / 2.0;
    }

    static double apply_mode(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("mode expects at least one argument");
        }
        std::vector<double> sorted = arguments;
        std::sort(sorted.begin(), sorted.end());
        double best_value = sorted.front();
        int best_count = 1;
        double current_value = sorted.front();
        int current_count = 1;
        for (std::size_t i = 1; i < sorted.size(); ++i) {
            if (mymath::is_near_zero(sorted[i] - current_value, 1e-10)) {
                ++current_count;
                continue;
            }
            if (current_count > best_count) {
                best_count = current_count;
                best_value = current_value;
            }
            current_value = sorted[i];
            current_count = 1;
        }
        if (current_count > best_count) {
            best_value = current_value;
        }
        return best_value;
    }

    static double apply_variance(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("var expects at least one argument");
        }
        const double mean = apply_mean(arguments);
        double sum = 0.0;
        for (double value : arguments) {
            const double delta = value - mean;
            sum += delta * delta;
        }
        return sum / static_cast<double>(arguments.size());
    }

    static double apply_stddev(const std::vector<double>& arguments) {
        return mymath::sqrt(apply_variance(arguments));
    }

    static double apply_percentile(const std::vector<double>& arguments) {
        if (arguments.size() < 2) {
            throw std::runtime_error("percentile expects p followed by at least one value");
        }
        const double p = arguments[0];
        if (p < 0.0 || p > 100.0) {
            throw std::runtime_error("percentile p must be in [0, 100]");
        }
        std::vector<double> values(arguments.begin() + 1, arguments.end());
        std::sort(values.begin(), values.end());
        if (values.size() == 1) {
            return values.front();
        }
        const double position =
            p * static_cast<double>(values.size() - 1) / 100.0;
        const long long lower_index = floor_to_long_long(position);
        const long long upper_index = ceil_to_long_long(position);
        if (lower_index == upper_index) {
            return values[static_cast<std::size_t>(lower_index)];
        }
        const double fraction = position - static_cast<double>(lower_index);
        const double lower = values[static_cast<std::size_t>(lower_index)];
        const double upper = values[static_cast<std::size_t>(upper_index)];
        return lower + (upper - lower) * fraction;
    }

    static double apply_quartile(const std::vector<double>& arguments) {
        if (arguments.size() < 2) {
            throw std::runtime_error("quartile expects q followed by at least one value");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("quartile q must be an integer");
        }
        const long long q = round_to_long_long(arguments[0]);
        if (q < 0 || q > 4) {
            throw std::runtime_error("quartile q must be between 0 and 4");
        }
        std::vector<double> percentile_arguments;
        percentile_arguments.reserve(arguments.size());
        percentile_arguments.push_back(static_cast<double>(q * 25));
        percentile_arguments.insert(percentile_arguments.end(),
                                    arguments.begin() + 1,
                                    arguments.end());
        return apply_percentile(percentile_arguments);
    }

    static double apply_factorial(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("factorial expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("factorial only accepts integers");
        }
        return factorial_value(round_to_long_long(arguments[0]));
    }

    static double apply_ncr(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("nCr expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("nCr only accepts integers");
        }
        return combination_value(round_to_long_long(arguments[0]),
                                 round_to_long_long(arguments[1]));
    }

    static double apply_npr(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("nPr expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("nPr only accepts integers");
        }
        return permutation_value(round_to_long_long(arguments[0]),
                                 round_to_long_long(arguments[1]));
    }

    static double apply_fib(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("fib expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("fib only accepts integers");
        }
        return fibonacci_value(round_to_long_long(arguments[0]));
    }

    static double apply_is_prime(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("is_prime expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("is_prime only accepts integers");
        }
        return is_prime_ll(round_to_long_long(arguments[0])) ? 1.0 : 0.0;
    }

    static double apply_next_prime(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("next_prime expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("next_prime only accepts integers");
        }
        return static_cast<double>(next_prime_ll(round_to_long_long(arguments[0])));
    }

    static double apply_rand(const std::vector<double>& arguments) {
        if (!arguments.empty()) {
            throw std::runtime_error("rand expects no arguments");
        }
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return distribution(global_rng());
    }

    static double apply_randn(const std::vector<double>& arguments) {
        if (!arguments.empty()) {
            throw std::runtime_error("randn expects no arguments");
        }
        std::normal_distribution<double> distribution(0.0, 1.0);
        return distribution(global_rng());
    }

    static double apply_randint(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("randint expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("randint only accepts integers");
        }
        long long left = round_to_long_long(arguments[0]);
        long long right = round_to_long_long(arguments[1]);
        if (left > right) {
            std::swap(left, right);
        }
        std::uniform_int_distribution<long long> distribution(left, right);
        return static_cast<double>(distribution(global_rng()));
    }

    static double apply_beta(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("beta expects exactly two arguments");
        }
        return mymath::beta(arguments[0], arguments[1]);
    }

    static double apply_zeta(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("zeta expects exactly one argument");
        }
        return mymath::zeta(arguments[0]);
    }

    static double apply_bessel(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("bessel expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("bessel order must be an integer");
        }
        return mymath::bessel_j(static_cast<int>(round_to_long_long(arguments[0])),
                                arguments[1]);
    }

    static double apply_pdf_normal(const std::vector<double>& arguments) {
        if (arguments.size() != 3) {
            throw std::runtime_error("pdf_normal expects exactly three arguments");
        }
        return normal_pdf(arguments[0], arguments[1], arguments[2]);
    }

    static double apply_cdf_normal(const std::vector<double>& arguments) {
        if (arguments.size() != 3) {
            throw std::runtime_error("cdf_normal expects exactly three arguments");
        }
        return normal_cdf(arguments[0], arguments[1], arguments[2]);
    }

    static double apply_and(const std::vector<double>& arguments) {
        // 位运算统一通过整数检查后再执行，避免隐式截断带来困惑。
        if (arguments.size() != 2) {
            throw std::runtime_error("and expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("and only accepts integers");
        }
        return static_cast<double>(
            round_to_long_long(arguments[0]) & round_to_long_long(arguments[1]));
    }

    static double apply_or(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("or expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("or only accepts integers");
        }
        return static_cast<double>(
            round_to_long_long(arguments[0]) | round_to_long_long(arguments[1]));
    }

    static double apply_xor(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("xor expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("xor only accepts integers");
        }
        return static_cast<double>(
            round_to_long_long(arguments[0]) ^ round_to_long_long(arguments[1]));
    }

    static double apply_shl(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("shl expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("shl only accepts integers");
        }
        const long long shift = round_to_long_long(arguments[1]);
        if (shift < 0) {
            throw std::runtime_error("shift count cannot be negative");
        }
        return static_cast<double>(round_to_long_long(arguments[0]) << shift);
    }

    static double apply_shr(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("shr expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("shr only accepts integers");
        }
        const long long shift = round_to_long_long(arguments[1]);
        if (shift < 0) {
            throw std::runtime_error("shift count cannot be negative");
        }
        return static_cast<double>(round_to_long_long(arguments[0]) >> shift);
    }

    static double apply_rol(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("rol expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("rol only accepts integers");
        }
        const std::uint64_t value = to_unsigned_bits(round_to_long_long(arguments[0]));
        const unsigned count = normalize_rotation_count(round_to_long_long(arguments[1]));
        return static_cast<double>(from_unsigned_bits(rotate_left_bits(value, count)));
    }

    static double apply_ror(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("ror expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("ror only accepts integers");
        }
        const std::uint64_t value = to_unsigned_bits(round_to_long_long(arguments[0]));
        const unsigned count = normalize_rotation_count(round_to_long_long(arguments[1]));
        return static_cast<double>(from_unsigned_bits(rotate_right_bits(value, count)));
    }

    static double apply_popcount(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("popcount expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("popcount only accepts integers");
        }
        return static_cast<double>(
            popcount_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_bitlen(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("bitlen expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("bitlen only accepts integers");
        }
        return static_cast<double>(
            bit_length_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_ctz(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("ctz expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("ctz only accepts integers");
        }
        return static_cast<double>(
            trailing_zero_count_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_clz(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("clz expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("clz only accepts integers");
        }
        return static_cast<double>(
            leading_zero_count_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_parity(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("parity expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("parity only accepts integers");
        }
        return static_cast<double>(
            parity_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_reverse_bits(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("reverse_bits expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("reverse_bits only accepts integers");
        }
        return static_cast<double>(
            from_unsigned_bits(reverse_bits(to_unsigned_bits(round_to_long_long(arguments[0])))));
    }

    double lookup_variable(const std::string& name) const {
        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            double constant_value = 0.0;
            if (lookup_builtin_constant(name, &constant_value)) {
                return constant_value;
            }
            throw std::runtime_error("unknown variable: " + name);
        }
        if (it->second.is_matrix) {
            throw std::runtime_error("matrix variable " + name + " cannot be used as a scalar");
        }
        if (it->second.is_string) {
            throw std::runtime_error("string variable " + name + " cannot be used as a number");
        }
        return it->second.exact ? rational_to_double(it->second.rational)
                                : it->second.decimal;
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
                ++pos_;
            } else {
                break;
            }
        }
        return source_.substr(start, pos_ - start);
    }

    bool peek_is_alpha() const {
        return !is_at_end() &&
               std::isalpha(static_cast<unsigned char>(source_[pos_]));
    }

    bool peek(char expected) const {
        return !is_at_end() && source_[pos_] == expected;
    }

    void skip_spaces() {
        while (!is_at_end() &&
               std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }

    bool match(char expected) {
        if (is_at_end() || source_[pos_] != expected) {
            return false;
        }
        ++pos_;
        return true;
    }

    bool match_string(const std::string& text) {
        if (source_.compare(pos_, text.size(), text) != 0) {
            return false;
        }
        pos_ += text.size();
        return true;
    }

    void expect(char expected) {
        if (!match(expected)) {
            throw std::runtime_error(std::string("expected '") + expected + "'");
        }
    }

    bool is_at_end() const {
        return pos_ >= source_.size();
    }

    std::string source_;
    std::size_t pos_ = 0;
    const std::map<std::string, StoredValue>* variables_;
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionDecimalCallback invoke_script_function_;
};

bool try_evaluate_matrix_expression(const std::string& expression,
                                    const std::map<std::string, StoredValue>* variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const HasScriptFunctionCallback& has_script_function,
                                    const InvokeScriptFunctionDecimalCallback& invoke_script_function,
                                    matrix::Value* value) {
    const matrix::ScalarEvaluator scalar_evaluator =
        [variables, functions, has_script_function, invoke_script_function](const std::string& text) {
            DecimalParser parser(text,
                                 variables,
                                 functions,
                                 has_script_function,
                                 invoke_script_function);
            const double scalar_value = parser.parse();
            return mymath::is_near_zero(scalar_value, 1e-10) ? 0.0 : scalar_value;
        };
    const matrix::MatrixLookup matrix_lookup =
        [variables](const std::string& name, matrix::Matrix* matrix_value) {
            const auto it = variables->find(name);
            if (it == variables->end() || !it->second.is_matrix) {
                return false;
            }
            *matrix_value = it->second.matrix;
            return true;
        };
    return matrix::try_evaluate_expression(expression,
                                           scalar_evaluator,
                                           matrix_lookup,
                                           value);
}

// ExactParser 只处理能够保持为有理数的表达式。
// 当遇到 sin、pi 或非整数指数这类无法精确表示为分数的情况时，
// 它会抛出 ExactModeUnsupported，调用方再回退到普通浮点模式。
class ExactParser {
public:
    ExactParser(std::string source,
                const std::map<std::string, StoredValue>* variables,
                const std::map<std::string, CustomFunction>* functions,
                HasScriptFunctionCallback has_script_function = {})
        : source_(std::move(source)),
          variables_(variables),
          functions_(functions),
          has_script_function_(std::move(has_script_function)) {}

    Rational parse() {
        Rational value = parse_comparison();
        skip_spaces();
        if (!is_at_end()) {
            throw std::runtime_error("unexpected token near: " + source_.substr(pos_, 1));
        }
        return value;
    }

private:
    Rational parse_comparison() {
        Rational value = parse_expression();
        while (true) {
            skip_spaces();
            if (match_string("==")) {
                const Rational rhs = parse_expression();
                value = Rational(value.numerator * rhs.denominator ==
                                         rhs.numerator * value.denominator
                                     ? 1
                                     : 0,
                                 1);
            } else if (match_string("!=")) {
                const Rational rhs = parse_expression();
                value = Rational(value.numerator * rhs.denominator !=
                                         rhs.numerator * value.denominator
                                     ? 1
                                     : 0,
                                 1);
            } else if (match_string("<=")) {
                const Rational rhs = parse_expression();
                value = Rational(rational_to_double(value) <= rational_to_double(rhs) ? 1 : 0, 1);
            } else if (match_string(">=")) {
                const Rational rhs = parse_expression();
                value = Rational(rational_to_double(value) >= rational_to_double(rhs) ? 1 : 0, 1);
            } else if (match('<')) {
                const Rational rhs = parse_expression();
                value = Rational(rational_to_double(value) < rational_to_double(rhs) ? 1 : 0, 1);
            } else if (match('>')) {
                const Rational rhs = parse_expression();
                value = Rational(rational_to_double(value) > rational_to_double(rhs) ? 1 : 0, 1);
            } else {
                break;
            }
        }
        return value;
    }

    Rational parse_expression() {
        Rational value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = value + parse_term();
            } else if (match('-')) {
                value = value - parse_term();
            } else {
                break;
            }
        }
        return value;
    }

    Rational parse_term() {
        Rational value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = value * parse_unary();
            } else if (match('/')) {
                value = value / parse_unary();
            } else {
                break;
            }
        }
        return value;
    }

    Rational parse_power() {
        Rational value = parse_primary();
        skip_spaces();
        if (match('^')) {
            const Rational exponent = parse_unary();
            if (!exponent.is_integer()) {
                throw ExactModeUnsupported("exact rational mode does not support non-integer exponents");
            }
            return pow_rational(value, exponent.numerator);
        }
        return value;
    }

    Rational parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            const Rational value = parse_unary();
            return Rational(-value.numerator, value.denominator);
        }
        return parse_power();
    }

    Rational parse_primary() {
        skip_spaces();
        if (match('(')) {
            const Rational value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        if (peek_is_alpha()) {
            const std::string name = parse_identifier();
            skip_spaces();
            if (!peek('(')) {
                return lookup_variable(name);
            }

            skip_spaces();
            expect('(');
            const std::vector<Rational> arguments = parse_argument_list();
            expect(')');
            return apply_function(name, arguments);
        }

        return parse_number();
    }

    std::vector<Rational> parse_argument_list() {
        std::vector<Rational> arguments;
        skip_spaces();
        if (peek(')')) {
            return arguments;
        }

        while (true) {
            arguments.push_back(parse_expression());
            skip_spaces();
            if (!match(',')) {
                break;
            }
        }
        return arguments;
    }

    Rational parse_number() {
        skip_spaces();

        // exact mode 也支持前缀整数，这样 0xFF 在分数模式里仍可作为整数参与计算。
        if (!is_at_end() &&
            source_[pos_] == '0' &&
            pos_ + 1 < source_.size()) {
            int base = 10;
            if (prefixed_base(source_[pos_ + 1], &base)) {
                const std::size_t start = pos_;
                pos_ += 2;
                while (!is_at_end()) {
                    const int digit = digit_value(source_[pos_]);
                    if (digit < 0 || digit >= base) {
                        break;
                    }
                    ++pos_;
                }
                return Rational(
                    parse_prefixed_integer_token(source_.substr(start, pos_ - start)), 1);
            }
        }

        const std::size_t start = pos_;
        bool has_digit = false;
        bool seen_dot = false;

        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                has_digit = true;
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else {
                break;
            }
        }

        if (!is_at_end() && (source_[pos_] == 'e' || source_[pos_] == 'E')) {
            const std::size_t exponent_pos = pos_;
            ++pos_;
            if (!is_at_end() && (source_[pos_] == '+' || source_[pos_] == '-')) {
                ++pos_;
            }
            const std::size_t exponent_digits = pos_;
            while (!is_at_end() &&
                   std::isdigit(static_cast<unsigned char>(source_[pos_]))) {
                ++pos_;
            }
            if (exponent_digits == pos_) {
                pos_ = exponent_pos;
            }
        }

        if (!has_digit) {
            throw std::runtime_error("expected number");
        }

        return parse_rational_literal(source_.substr(start, pos_ - start));
    }

    static Rational parse_rational_literal(const std::string& token) {
        std::string significand = token;
        long long exponent_adjust = 0;
        const std::size_t exponent_pos = token.find_first_of("eE");
        if (exponent_pos != std::string::npos) {
            significand = token.substr(0, exponent_pos);
            exponent_adjust = std::stoll(token.substr(exponent_pos + 1));
        }

        long long numerator = 0;
        long long denominator = 1;
        std::size_t idx = 0;

        while (idx < significand.size() && significand[idx] != '.') {
            numerator =
                numerator * 10 + static_cast<long long>(significand[idx] - '0');
            ++idx;
        }

        if (idx < significand.size() && significand[idx] == '.') {
            ++idx;
            while (idx < significand.size()) {
                numerator =
                    numerator * 10 + static_cast<long long>(significand[idx] - '0');
                denominator *= 10;
                ++idx;
            }
        }

        while (exponent_adjust > 0) {
            numerator *= 10;
            --exponent_adjust;
        }
        while (exponent_adjust < 0) {
            denominator *= 10;
            ++exponent_adjust;
        }

        return Rational(numerator, denominator);
    }

    Rational apply_function(const std::string& name, const std::vector<Rational>& arguments) {
        // exact mode 只实现“结果仍能保持为有理数”的函数。
        // 其他函数会抛出 ExactModeUnsupported，再由上层回退到浮点显示。
        if (functions_->find(name) != functions_->end()) {
            throw ExactModeUnsupported("custom function " + name +
                                       " is not supported exactly");
        }
        if (has_script_function_ && has_script_function_(name)) {
            throw ExactModeUnsupported("script function " + name +
                                       " is not supported exactly");
        }
        if (name == "pow") {
            if (arguments.size() != 2) {
                throw std::runtime_error("pow expects exactly two arguments");
            }
            if (!arguments[1].is_integer()) {
                throw ExactModeUnsupported("exact rational mode does not support non-integer exponents");
            }
            return pow_rational(arguments[0], arguments[1].numerator);
        }
        if (name == "abs") {
            if (arguments.size() != 1) {
                throw std::runtime_error("abs expects exactly one argument");
            }
            return abs_rational(arguments[0]);
        }
        if (name == "step" || name == "u" || name == "heaviside") {
            if (arguments.size() != 1) {
                throw std::runtime_error("step expects exactly one argument");
            }
            return Rational(arguments[0].numerator >= 0 ? 1 : 0, 1);
        }
        if (name == "delta" || name == "impulse") {
            if (arguments.size() != 1) {
                throw std::runtime_error("delta expects exactly one argument");
            }
            return Rational(arguments[0].numerator == 0 ? 1 : 0, 1);
        }
        if (name == "not") {
            if (arguments.size() != 1) {
                throw std::runtime_error("not expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw std::runtime_error("not only accepts integers");
            }
            return Rational(~arguments[0].numerator, 1);
        }
        if (name == "sign") {
            if (arguments.size() != 1) {
                throw std::runtime_error("sign expects exactly one argument");
            }
            if (arguments[0].numerator == 0) {
                return Rational(0, 1);
            }
            return Rational(arguments[0].numerator > 0 ? 1 : -1, 1);
        }
        if (name == "gcd") {
            if (arguments.size() != 2) {
                throw std::runtime_error("gcd expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("gcd only accepts integers");
            }
            return Rational(gcd_ll(arguments[0].numerator, arguments[1].numerator), 1);
        }
        if (name == "lcm") {
            if (arguments.size() != 2) {
                throw std::runtime_error("lcm expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("lcm only accepts integers");
            }
            return Rational(lcm_ll(arguments[0].numerator, arguments[1].numerator), 1);
        }
        if (name == "mod") {
            if (arguments.size() != 2) {
                throw std::runtime_error("mod expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("mod only accepts integers");
            }
            if (arguments[1].numerator == 0) {
                throw std::runtime_error("mod divisor cannot be zero");
            }
            return Rational(arguments[0].numerator % arguments[1].numerator, 1);
        }
        if (name == "rol") {
            if (arguments.size() != 2) {
                throw std::runtime_error("rol expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("rol only accepts integers");
            }
            const unsigned count = normalize_rotation_count(arguments[1].numerator);
            return Rational(
                from_unsigned_bits(rotate_left_bits(
                    to_unsigned_bits(arguments[0].numerator), count)),
                1);
        }
        if (name == "ror") {
            if (arguments.size() != 2) {
                throw std::runtime_error("ror expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("ror only accepts integers");
            }
            const unsigned count = normalize_rotation_count(arguments[1].numerator);
            return Rational(
                from_unsigned_bits(rotate_right_bits(
                    to_unsigned_bits(arguments[0].numerator), count)),
                1);
        }
        if (name == "floor") {
            if (arguments.size() != 1) {
                throw std::runtime_error("floor expects exactly one argument");
            }
            return Rational(floor_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "ceil") {
            if (arguments.size() != 1) {
                throw std::runtime_error("ceil expects exactly one argument");
            }
            return Rational(ceil_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "round") {
            if (arguments.size() != 1) {
                throw std::runtime_error("round expects exactly one argument");
            }
            return Rational(round_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "trunc") {
            if (arguments.size() != 1) {
                throw std::runtime_error("trunc expects exactly one argument");
            }
            return Rational(trunc_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "min") {
            if (arguments.size() != 2) {
                throw std::runtime_error("min expects exactly two arguments");
            }
            return rational_to_double(arguments[0]) < rational_to_double(arguments[1])
                       ? arguments[0]
                       : arguments[1];
        }
        if (name == "max") {
            if (arguments.size() != 2) {
                throw std::runtime_error("max expects exactly two arguments");
            }
            return rational_to_double(arguments[0]) > rational_to_double(arguments[1])
                       ? arguments[0]
                       : arguments[1];
        }
        if (name == "clamp") {
            if (arguments.size() != 3) {
                throw std::runtime_error("clamp expects exactly three arguments");
            }
            Rational lower = arguments[1];
            Rational upper = arguments[2];
            if (rational_to_double(lower) > rational_to_double(upper)) {
                std::swap(lower, upper);
            }
            if (rational_to_double(arguments[0]) < rational_to_double(lower)) {
                return lower;
            }
            if (rational_to_double(arguments[0]) > rational_to_double(upper)) {
                return upper;
            }
            return arguments[0];
        }
        if (name == "sum") {
            if (arguments.empty()) {
                throw std::runtime_error("sum expects at least one argument");
            }
            Rational total(0, 1);
            for (const Rational& value : arguments) {
                total = total + value;
            }
            return total;
        }
        if (name == "avg") {
            if (arguments.empty()) {
                throw std::runtime_error("avg expects at least one argument");
            }
            Rational total(0, 1);
            for (const Rational& value : arguments) {
                total = total + value;
            }
            return total / Rational(static_cast<long long>(arguments.size()), 1);
        }
        if (name == "mean") {
            if (arguments.empty()) {
                throw std::runtime_error("mean expects at least one argument");
            }
            Rational total(0, 1);
            for (const Rational& value : arguments) {
                total = total + value;
            }
            return total / Rational(static_cast<long long>(arguments.size()), 1);
        }
        if (name == "median") {
            if (arguments.empty()) {
                throw std::runtime_error("median expects at least one argument");
            }
            std::vector<Rational> sorted = arguments;
            std::sort(sorted.begin(), sorted.end(),
                      [](const Rational& lhs, const Rational& rhs) {
                          return rational_to_double(lhs) < rational_to_double(rhs);
                      });
            const std::size_t middle = sorted.size() / 2;
            if (sorted.size() % 2 == 1) {
                return sorted[middle];
            }
            return (sorted[middle - 1] + sorted[middle]) / Rational(2, 1);
        }
        if (name == "factorial") {
            if (arguments.size() != 1) {
                throw std::runtime_error("factorial expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw std::runtime_error("factorial only accepts integers");
            }
            return factorial_rational(arguments[0].numerator);
        }
        if (name == "nCr") {
            if (arguments.size() != 2) {
                throw std::runtime_error("nCr expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("nCr only accepts integers");
            }
            return combination_rational(arguments[0].numerator, arguments[1].numerator);
        }
        if (name == "binom") {
            if (arguments.size() != 2) {
                throw std::runtime_error("binom expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("binom only accepts integers");
            }
            return combination_rational(arguments[0].numerator, arguments[1].numerator);
        }
        if (name == "nPr") {
            if (arguments.size() != 2) {
                throw std::runtime_error("nPr expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("nPr only accepts integers");
            }
            return permutation_rational(arguments[0].numerator, arguments[1].numerator);
        }
        if (name == "popcount") {
            if (arguments.size() != 1) {
                throw std::runtime_error("popcount expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw std::runtime_error("popcount only accepts integers");
            }
            return Rational(popcount_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "bitlen") {
            if (arguments.size() != 1) {
                throw std::runtime_error("bitlen expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw std::runtime_error("bitlen only accepts integers");
            }
            return Rational(bit_length_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "ctz") {
            if (arguments.size() != 1) {
                throw std::runtime_error("ctz expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw std::runtime_error("ctz only accepts integers");
            }
            return Rational(trailing_zero_count_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "clz") {
            if (arguments.size() != 1) {
                throw std::runtime_error("clz expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw std::runtime_error("clz only accepts integers");
            }
            return Rational(leading_zero_count_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "parity") {
            if (arguments.size() != 1) {
                throw std::runtime_error("parity expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw std::runtime_error("parity only accepts integers");
            }
            return Rational(parity_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "reverse_bits") {
            if (arguments.size() != 1) {
                throw std::runtime_error("reverse_bits expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw std::runtime_error("reverse_bits only accepts integers");
            }
            return Rational(
                from_unsigned_bits(reverse_bits(to_unsigned_bits(arguments[0].numerator))),
                1);
        }
        if (name == "and") {
            if (arguments.size() != 2) {
                throw std::runtime_error("and expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("and only accepts integers");
            }
            return Rational(arguments[0].numerator & arguments[1].numerator, 1);
        }
        if (name == "or") {
            if (arguments.size() != 2) {
                throw std::runtime_error("or expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("or only accepts integers");
            }
            return Rational(arguments[0].numerator | arguments[1].numerator, 1);
        }
        if (name == "xor") {
            if (arguments.size() != 2) {
                throw std::runtime_error("xor expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("xor only accepts integers");
            }
            return Rational(arguments[0].numerator ^ arguments[1].numerator, 1);
        }
        if (name == "shl") {
            if (arguments.size() != 2) {
                throw std::runtime_error("shl expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("shl only accepts integers");
            }
            if (arguments[1].numerator < 0) {
                throw std::runtime_error("shift count cannot be negative");
            }
            return Rational(arguments[0].numerator << arguments[1].numerator, 1);
        }
        if (name == "shr") {
            if (arguments.size() != 2) {
                throw std::runtime_error("shr expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw std::runtime_error("shr only accepts integers");
            }
            if (arguments[1].numerator < 0) {
                throw std::runtime_error("shift count cannot be negative");
            }
            return Rational(arguments[0].numerator >> arguments[1].numerator, 1);
        }

        throw ExactModeUnsupported("function " + name + " is not supported exactly");
    }

    Rational lookup_variable(const std::string& name) const {
        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            double constant_value = 0.0;
            if (lookup_builtin_constant(name, &constant_value)) {
                throw ExactModeUnsupported("built-in constants are not rational");
            }
            throw std::runtime_error("unknown variable: " + name);
        }
        if (it->second.is_matrix) {
            throw ExactModeUnsupported("matrix variable " + name + " cannot be used exactly");
        }
        if (it->second.is_string) {
            throw ExactModeUnsupported("string variable " + name + " cannot be used exactly");
        }
        if (!it->second.exact) {
            throw ExactModeUnsupported("variable " + name + " is only stored approximately");
        }
        return it->second.rational;
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
                ++pos_;
            } else {
                break;
            }
        }
        return source_.substr(start, pos_ - start);
    }

    bool peek_is_alpha() const {
        return !is_at_end() &&
               std::isalpha(static_cast<unsigned char>(source_[pos_]));
    }

    bool peek(char expected) const {
        return !is_at_end() && source_[pos_] == expected;
    }

    void skip_spaces() {
        while (!is_at_end() &&
               std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }

    bool match(char expected) {
        if (is_at_end() || source_[pos_] != expected) {
            return false;
        }
        ++pos_;
        return true;
    }

    bool match_string(const std::string& text) {
        if (source_.compare(pos_, text.size(), text) != 0) {
            return false;
        }
        pos_ += text.size();
        return true;
    }

    void expect(char expected) {
        if (!match(expected)) {
            throw std::runtime_error(std::string("expected '") + expected + "'");
        }
    }

    bool is_at_end() const {
        return pos_ >= source_.size();
    }

    std::string source_;
    std::size_t pos_ = 0;
    const std::map<std::string, StoredValue>* variables_;
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
};

struct HexFormatOptions {
    bool prefix = false;
    bool uppercase = true;
};

bool convert_base_value(long long value,
                        int base,
                        const HexFormatOptions& hex_options,
                        std::string* output) {
    if (base < 2 || base > 16) {
        return false;
    }

    static const char upper_digits[] = "0123456789ABCDEF";
    static const char lower_digits[] = "0123456789abcdef";
    const char* digits = hex_options.uppercase ? upper_digits : lower_digits;
    if (value == 0) {
        *output = base == 16 && hex_options.prefix ? "0x0" : "0";
        return true;
    }

    bool negative = value < 0;
    unsigned long long current = negative
                                     ? static_cast<unsigned long long>(-(value + 1)) + 1ULL
                                     : static_cast<unsigned long long>(value);

    std::string reversed;
    while (current > 0) {
        reversed.push_back(digits[current % static_cast<unsigned long long>(base)]);
        current /= static_cast<unsigned long long>(base);
    }

    output->clear();
    if (negative) {
        output->push_back('-');
    }
    if (base == 16 && hex_options.prefix) {
        output->push_back('0');
        output->push_back('x');
    }
    for (std::size_t i = reversed.size(); i > 0; --i) {
        output->push_back(reversed[i - 1]);
    }
    return true;
}

bool try_base_conversion_expression(const std::string& expression,
                                    const std::map<std::string, StoredValue>* variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const HexFormatOptions& hex_options,
                                    std::string* output) {
    // 进制转换是“显示型功能”：
    // 先把参数当表达式求成整数，再格式化成目标进制字符串。
    std::string inside;
    std::string mode;

    if (split_named_call(expression, "bin", &inside)) {
        mode = "bin";
    } else if (split_named_call(expression, "oct", &inside)) {
        mode = "oct";
    } else if (split_named_call(expression, "hex", &inside)) {
        mode = "hex";
    } else if (split_named_call(expression, "base", &inside)) {
        mode = "base";
    } else {
        return false;
    }

    const std::vector<std::string> arguments = split_top_level_arguments(inside);
    int base = 10;

    if (mode == "bin" || mode == "oct" || mode == "hex") {
        if (arguments.size() != 1) {
            throw std::runtime_error(mode + " expects exactly one argument");
        }
        base = mode == "bin" ? 2 : (mode == "oct" ? 8 : 16);
    } else {
        if (arguments.size() != 2) {
            throw std::runtime_error("base expects exactly two arguments");
        }
        DecimalParser base_parser(arguments[1], variables, functions);
        const double base_value = base_parser.parse();
        if (!is_integer_double(base_value)) {
            throw std::runtime_error("base conversion requires an integer base");
        }
        base = static_cast<int>(round_to_long_long(base_value));
    }

    DecimalParser value_parser(arguments[0], variables, functions);
    const double value = value_parser.parse();
    if (!is_integer_double(value)) {
        throw std::runtime_error("base conversion only accepts integers");
    }

    if (!convert_base_value(round_to_long_long(value), base, hex_options, output)) {
        throw std::runtime_error("base must be in the range [2, 16]");
    }

    return true;
}

std::string scalar_value_expression_text(const StoredValue& value) {
    if (value.has_symbolic_text) {
        return value.symbolic_text;
    }
    if (value.exact) {
        return value.rational.to_string();
    }
    if (value.has_precise_decimal_text) {
        return value.precise_decimal_text;
    }
    return format_decimal(normalize_display_decimal(value.decimal));
}

std::string matrix_literal_expression(const matrix::Matrix& value) {
    std::ostringstream out;
    out << '[';
    for (std::size_t row = 0; row < value.rows; ++row) {
        if (row != 0) {
            out << "; ";
        }
        for (std::size_t col = 0; col < value.cols; ++col) {
            if (col != 0) {
                out << ", ";
            }
            out << format_decimal(normalize_display_decimal(value.at(row, col)));
        }
    }
    out << ']';
    return out.str();
}

bool is_supported_symbolic_unary_function(const std::string& name) {
    return name == "sin" || name == "cos" || name == "tan" ||
           name == "asin" || name == "acos" || name == "atan" ||
           name == "exp" || name == "ln" || name == "log10" ||
           name == "sqrt" || name == "abs" || name == "sign" ||
           name == "floor" || name == "ceil" || name == "cbrt" ||
           name == "step" || name == "delta";
}

class SymbolicRenderParser {
public:
    SymbolicRenderParser(std::string source,
                         const std::map<std::string, StoredValue>* variables,
                         const std::map<std::string, CustomFunction>* functions,
                         int depth = 0)
        : source_(std::move(source)),
          variables_(variables),
          functions_(functions),
          depth_(depth) {}

    bool parse(std::string* output, bool* used_symbolic_constant) {
        if (depth_ > 12) {
            return false;
        }
        try {
            used_symbolic_constant_ = false;
            const std::string text = parse_expression();
            skip_spaces();
            if (pos_ != source_.size()) {
                return false;
            }
            SymbolicExpression expression = SymbolicExpression::parse(text);
            *output = expression.to_string();
            *used_symbolic_constant = used_symbolic_constant_;
            return true;
        } catch (const std::exception&) {
            return false;
        }
    }

private:
    std::string parse_expression() {
        std::string value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = "(" + value + " + " + parse_term() + ")";
            } else if (match('-')) {
                value = "(" + value + " - " + parse_term() + ")";
            } else {
                return value;
            }
        }
    }

    std::string parse_term() {
        std::string value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = "(" + value + " * " + parse_unary() + ")";
            } else if (match('/')) {
                value = "(" + value + " / " + parse_unary() + ")";
            } else {
                return value;
            }
        }
    }

    std::string parse_power() {
        std::string value = parse_primary();
        skip_spaces();
        if (match('^')) {
            value = "(" + value + " ^ " + parse_unary() + ")";
        }
        return value;
    }

    std::string parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            return "(-" + parse_unary() + ")";
        }
        return parse_power();
    }

    std::string parse_primary() {
        skip_spaces();
        if (match('(')) {
            const std::string value = parse_expression();
            skip_spaces();
            expect(')');
            return "(" + value + ")";
        }
        if (peek_is_alpha()) {
            const std::string name = parse_identifier();
            skip_spaces();
            if (!peek('(')) {
                return render_identifier(name);
            }

            expect('(');
            const std::vector<std::string> arguments = parse_argument_list();
            expect(')');
            return render_function(name, arguments);
        }
        return parse_number_token();
    }

    std::vector<std::string> parse_argument_list() {
        std::vector<std::string> arguments;
        skip_spaces();
        if (peek(')')) {
            return arguments;
        }

        while (true) {
            arguments.push_back(parse_expression());
            skip_spaces();
            if (!match(',')) {
                break;
            }
        }
        return arguments;
    }

    std::string render_identifier(const std::string& name) {
        if (name == "pi" || name == "e") {
            used_symbolic_constant_ = true;
            return name;
        }

        double builtin_constant = 0.0;
        if (lookup_builtin_constant(name, &builtin_constant)) {
            return format_symbolic_scalar(builtin_constant);
        }

        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            return name;
        }
        if (it->second.is_matrix || it->second.is_string) {
            throw std::runtime_error("unsupported symbolic variable");
        }
        if (it->second.has_symbolic_text) {
            used_symbolic_constant_ = true;
        }
        return "(" + scalar_value_expression_text(it->second) + ")";
    }

    std::string render_function(const std::string& name,
                                const std::vector<std::string>& arguments) {
        if (name == "pow") {
            if (arguments.size() != 2) {
                throw std::runtime_error("pow expects two arguments");
            }
            return "((" + arguments[0] + ") ^ (" + arguments[1] + "))";
        }
        if (name == "root") {
            if (arguments.size() != 2) {
                throw std::runtime_error("root expects two arguments");
            }
            return "((" + arguments[0] + ") ^ (1 / (" + arguments[1] + ")))";
        }
        const auto function_it = functions_->find(name);
        if (function_it != functions_->end()) {
            if (arguments.size() != 1) {
                throw std::runtime_error("custom function expects one argument");
            }
            std::map<std::string, StoredValue> scoped_variables = *variables_;
            StoredValue parameter_value;
            parameter_value.has_symbolic_text = true;
            parameter_value.symbolic_text = arguments[0];
            scoped_variables[function_it->second.parameter_name] = parameter_value;

            SymbolicRenderParser nested(function_it->second.expression,
                                        &scoped_variables,
                                        functions_,
                                        depth_ + 1);
            std::string expanded;
            bool nested_symbolic = false;
            if (!nested.parse(&expanded, &nested_symbolic)) {
                throw std::runtime_error("unable to expand custom function symbolically");
            }
            used_symbolic_constant_ = used_symbolic_constant_ || nested_symbolic;
            return "(" + expanded + ")";
        }
        if (!is_supported_symbolic_unary_function(name) || arguments.size() != 1) {
            throw std::runtime_error("unsupported symbolic function");
        }
        return name + "(" + arguments[0] + ")";
    }

    std::string parse_number_token() {
        skip_spaces();

        if (!is_at_end() &&
            source_[pos_] == '0' &&
            pos_ + 1 < source_.size()) {
            int base = 10;
            if (prefixed_base(source_[pos_ + 1], &base)) {
                const std::size_t start = pos_;
                pos_ += 2;
                while (!is_at_end()) {
                    const int digit = digit_value(source_[pos_]);
                    if (digit < 0 || digit >= base) {
                        break;
                    }
                    ++pos_;
                }
                return format_decimal(static_cast<double>(
                    parse_prefixed_integer_token(source_.substr(start, pos_ - start))));
            }
        }

        const std::size_t start = pos_;
        bool has_digit = false;
        bool seen_dot = false;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                has_digit = true;
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else {
                break;
            }
        }
        if (!is_at_end() && (source_[pos_] == 'e' || source_[pos_] == 'E')) {
            const std::size_t exponent_pos = pos_;
            ++pos_;
            if (!is_at_end() && (source_[pos_] == '+' || source_[pos_] == '-')) {
                ++pos_;
            }
            const std::size_t exponent_digits = pos_;
            while (!is_at_end() &&
                   std::isdigit(static_cast<unsigned char>(source_[pos_]))) {
                ++pos_;
            }
            if (exponent_digits == pos_) {
                pos_ = exponent_pos;
            }
        }
        if (!has_digit) {
            throw std::runtime_error("expected number");
        }
        return format_decimal(std::stod(source_.substr(start, pos_ - start)));
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
                ++pos_;
            } else {
                break;
            }
        }
        return source_.substr(start, pos_ - start);
    }

    bool peek_is_alpha() const {
        return !is_at_end() &&
               std::isalpha(static_cast<unsigned char>(source_[pos_]));
    }

    bool peek(char expected) const {
        return !is_at_end() && source_[pos_] == expected;
    }

    void skip_spaces() {
        while (!is_at_end() &&
               std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }

    bool match(char expected) {
        if (is_at_end() || source_[pos_] != expected) {
            return false;
        }
        ++pos_;
        return true;
    }

    void expect(char expected) {
        if (!match(expected)) {
            throw std::runtime_error(std::string("expected '") + expected + "'");
        }
    }

    bool is_at_end() const {
        return pos_ >= source_.size();
    }

    std::string source_;
    const std::map<std::string, StoredValue>* variables_;
    const std::map<std::string, CustomFunction>* functions_;
    std::size_t pos_ = 0;
    int depth_ = 0;
    bool used_symbolic_constant_ = false;
};

bool try_symbolic_constant_expression(const std::string& expression,
                                      const std::map<std::string, StoredValue>* variables,
                                      const std::map<std::string, CustomFunction>* functions,
                                      std::string* output) {
    SymbolicRenderParser parser(expression, variables, functions);
    bool used_symbolic_constant = false;
    if (!parser.parse(output, &used_symbolic_constant)) {
        return false;
    }
    return used_symbolic_constant;
}

}  // namespace

struct Calculator::Impl {
    std::map<std::string, StoredValue> variables;
    std::map<std::string, CustomFunction> functions;
    std::map<std::string, ScriptFunction> script_functions;
    std::vector<std::map<std::string, StoredValue>> local_scopes;
    bool symbolic_constants_mode = false;
    bool hex_prefix_mode = false;
    bool hex_uppercase_mode = true;
};

std::map<std::string, StoredValue> visible_variables(const Calculator::Impl* impl) {
    std::map<std::string, StoredValue> merged = impl->variables;
    for (const auto& scope : impl->local_scopes) {
        for (const auto& [name, value] : scope) {
            merged[name] = value;
        }
    }
    return merged;
}

bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name) {
    return impl->script_functions.find(name) != impl->script_functions.end();
}

void assign_visible_variable(Calculator::Impl* impl,
                             const std::string& name,
                             const StoredValue& value) {
    for (auto it = impl->local_scopes.rbegin(); it != impl->local_scopes.rend(); ++it) {
        const auto existing = it->find(name);
        if (existing != it->end()) {
            (*it)[name] = value;
            return;
        }
    }

    const auto global_existing = impl->variables.find(name);
    if (global_existing != impl->variables.end()) {
        impl->variables[name] = value;
        return;
    }

    if (!impl->local_scopes.empty()) {
        impl->local_scopes.back()[name] = value;
        return;
    }

    impl->variables[name] = value;
}

struct ScriptSignal {
    enum class Kind {
        kNone,
        kReturn,
        kBreak,
        kContinue,
    };

    Kind kind = Kind::kNone;
    bool has_value = false;
    StoredValue value;

    static ScriptSignal make_return(const StoredValue& return_value) {
        ScriptSignal signal;
        signal.kind = Kind::kReturn;
        signal.has_value = true;
        signal.value = return_value;
        return signal;
    }

    static ScriptSignal make_break() {
        ScriptSignal signal;
        signal.kind = Kind::kBreak;
        return signal;
    }

    static ScriptSignal make_continue() {
        ScriptSignal signal;
        signal.kind = Kind::kContinue;
        return signal;
    }
};

StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode);
ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope);
ScriptSignal execute_script_block(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope);

script::StatementPtr clone_statement(const script::Statement& statement);
std::string render_script_statement(const script::Statement& statement, int indent);
std::string render_script_block(const script::BlockStatement& block, int indent);

std::string indent_text(int indent) {
    return std::string(static_cast<std::size_t>(indent) * 2, ' ');
}

std::unique_ptr<script::BlockStatement> clone_block_statement(const script::BlockStatement& block) {
    auto clone = std::make_unique<script::BlockStatement>();
    for (const auto& statement : block.statements) {
        clone->statements.push_back(clone_statement(*statement));
    }
    return clone;
}

std::string render_script_block(const script::BlockStatement& block, int indent) {
    std::ostringstream out;
    out << "{\n";
    for (const auto& statement : block.statements) {
        out << render_script_statement(*statement, indent + 1);
    }
    out << indent_text(indent) << "}";
    return out.str();
}

std::string render_script_statement(const script::Statement& statement, int indent) {
    const std::string prefix = indent_text(indent);
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return prefix + render_script_block(
                                static_cast<const script::BlockStatement&>(statement),
                                indent) + "\n";
        case script::Statement::Kind::kSimple:
            return prefix + static_cast<const script::SimpleStatement&>(statement).text + ";\n";
        case script::Statement::Kind::kIf: {
            const auto& source = static_cast<const script::IfStatement&>(statement);
            std::string rendered =
                prefix + "if (" + source.condition + ") " +
                render_script_statement(*source.then_branch, indent).substr(prefix.size());
            if (source.else_branch) {
                rendered.pop_back();
                rendered += prefix + "else " +
                            render_script_statement(*source.else_branch, indent).substr(prefix.size());
            }
            return rendered;
        }
        case script::Statement::Kind::kWhile: {
            const auto& source = static_cast<const script::WhileStatement&>(statement);
            return prefix + "while (" + source.condition + ") " +
                   render_script_statement(*source.body, indent).substr(prefix.size());
        }
        case script::Statement::Kind::kFor: {
            const auto& source = static_cast<const script::ForStatement&>(statement);
            return prefix + "for (" + source.initializer + "; " + source.condition + "; " +
                   source.step + ") " +
                   render_script_statement(*source.body, indent).substr(prefix.size());
        }
        case script::Statement::Kind::kFunction: {
            const auto& source = static_cast<const script::FunctionStatement&>(statement);
            std::ostringstream out;
            out << prefix << "fn " << source.name << "(";
            for (std::size_t i = 0; i < source.parameters.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << source.parameters[i];
            }
            out << ") " << render_script_block(*source.body, indent) << '\n';
            return out.str();
        }
        case script::Statement::Kind::kReturn: {
            const auto& source = static_cast<const script::ReturnStatement&>(statement);
            return prefix + "return" +
                   (source.has_expression ? " " + source.expression : "") +
                   ";\n";
        }
        case script::Statement::Kind::kBreak:
            return prefix + "break;\n";
        case script::Statement::Kind::kContinue:
            return prefix + "continue;\n";
    }

    throw std::runtime_error("unknown script statement kind");
}

script::StatementPtr clone_statement(const script::Statement& statement) {
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return clone_block_statement(static_cast<const script::BlockStatement&>(statement));
        case script::Statement::Kind::kSimple: {
            auto clone = std::make_unique<script::SimpleStatement>();
            clone->text = static_cast<const script::SimpleStatement&>(statement).text;
            return clone;
        }
        case script::Statement::Kind::kIf: {
            const auto& source = static_cast<const script::IfStatement&>(statement);
            auto clone = std::make_unique<script::IfStatement>();
            clone->condition = source.condition;
            clone->then_branch = clone_statement(*source.then_branch);
            if (source.else_branch) {
                clone->else_branch = clone_statement(*source.else_branch);
            }
            return clone;
        }
        case script::Statement::Kind::kWhile: {
            const auto& source = static_cast<const script::WhileStatement&>(statement);
            auto clone = std::make_unique<script::WhileStatement>();
            clone->condition = source.condition;
            clone->body = clone_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kFor: {
            const auto& source = static_cast<const script::ForStatement&>(statement);
            auto clone = std::make_unique<script::ForStatement>();
            clone->initializer = source.initializer;
            clone->condition = source.condition;
            clone->step = source.step;
            clone->body = clone_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kFunction: {
            const auto& source = static_cast<const script::FunctionStatement&>(statement);
            auto clone = std::make_unique<script::FunctionStatement>();
            clone->name = source.name;
            clone->parameters = source.parameters;
            clone->body = clone_block_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kReturn: {
            const auto& source = static_cast<const script::ReturnStatement&>(statement);
            auto clone = std::make_unique<script::ReturnStatement>();
            clone->has_expression = source.has_expression;
            clone->expression = source.expression;
            return clone;
        }
        case script::Statement::Kind::kBreak:
            return std::make_unique<script::BreakStatement>();
        case script::Statement::Kind::kContinue:
            return std::make_unique<script::ContinueStatement>();
    }

    throw std::runtime_error("unknown script statement kind");
}

bool truthy_value(const StoredValue& value) {
    if (value.is_matrix) {
        throw std::runtime_error("matrix values cannot be used as script conditions");
    }
    return !mymath::is_near_zero(value.exact
                                     ? rational_to_double(value.rational)
                                     : value.decimal,
                                 1e-10);
}

double invoke_script_function_decimal(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& name,
                                      const std::vector<double>& arguments) {
    auto it = impl->script_functions.find(name);
    if (it == impl->script_functions.end()) {
        throw std::runtime_error("unknown function: " + name);
    }

    const ScriptFunction& function = it->second;
    if (arguments.size() != function.parameter_names.size()) {
        throw std::runtime_error("script function " + name + " expects " +
                                 std::to_string(function.parameter_names.size()) +
                                 " arguments");
    }

    std::map<std::string, StoredValue> frame;
    for (std::size_t i = 0; i < arguments.size(); ++i) {
        StoredValue parameter_value;
        parameter_value.decimal = arguments[i];
        parameter_value.exact = false;
        frame[function.parameter_names[i]] = parameter_value;
    }

    impl->local_scopes.push_back(frame);
    std::string ignored_output;
    try {
        const ScriptSignal signal =
            execute_script_block(calculator, impl, *function.body, false, &ignored_output, false);
        impl->local_scopes.pop_back();

        if (signal.kind != ScriptSignal::Kind::kReturn || !signal.has_value) {
            throw std::runtime_error("script function " + name + " must return a value");
        }
        if (signal.value.is_matrix) {
            throw std::runtime_error("script function " + name +
                                     " cannot be used as a scalar expression");
        }
        if (signal.value.is_string) {
            throw std::runtime_error("script function " + name +
                                     " cannot be used as a numeric expression");
        }
        return signal.value.exact
                   ? rational_to_double(signal.value.rational)
                   : signal.value.decimal;
    } catch (...) {
        impl->local_scopes.pop_back();
        throw;
    }
}

StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode) {
    const std::string expanded_expression =
        expand_inline_function_commands(calculator, expression);
    const std::string trimmed = trim_copy(expanded_expression);
    const std::map<std::string, StoredValue> variables = visible_variables(impl);
    if (is_string_literal(trimmed)) {
        StoredValue stored;
        stored.is_string = true;
        stored.string_value = parse_string_literal_value(trimmed);
        return stored;
    }
    if (is_identifier_text(trimmed)) {
        const auto it = variables.find(trimmed);
        if (it != variables.end() && it->second.is_string) {
            return it->second;
        }
    }

    const HasScriptFunctionCallback has_script_function =
        [impl](const std::string& name) {
            return has_visible_script_function(impl, name);
        };
    const InvokeScriptFunctionDecimalCallback invoke_script_function =
        [calculator, impl](const std::string& name, const std::vector<double>& arguments) {
            return invoke_script_function_decimal(calculator, impl, name, arguments);
        };

    StoredValue stored;
    std::vector<std::string> rational_arguments;
    if (split_named_call_with_arguments(trimmed, "rat", &rational_arguments)) {
        if (rational_arguments.size() != 1 && rational_arguments.size() != 2) {
            throw std::runtime_error(
                "rat expects one argument or expression plus max_denominator");
        }

        const StoredValue value =
            evaluate_expression_value(calculator, impl, rational_arguments[0], false);
        if (value.is_matrix) {
            throw std::runtime_error("rat cannot approximate a matrix value");
        }
        if (value.is_string) {
            throw std::runtime_error("rat cannot approximate a string value");
        }

        long long max_denominator = 999;
        if (rational_arguments.size() == 2) {
            const StoredValue max_denominator_value =
                evaluate_expression_value(calculator, impl, rational_arguments[1], false);
            if (max_denominator_value.is_matrix || max_denominator_value.is_string) {
                throw std::runtime_error("rat max_denominator must be a positive integer");
            }
            const double scalar =
                max_denominator_value.exact
                    ? rational_to_double(max_denominator_value.rational)
                    : max_denominator_value.decimal;
            if (!is_integer_double(scalar) || scalar <= 0.0) {
                throw std::runtime_error("rat max_denominator must be a positive integer");
            }
            max_denominator = round_to_long_long(scalar);
        }

        if (value.exact && value.rational.denominator <= max_denominator) {
            return value;
        }

        const double decimal_value = value.exact
                                         ? rational_to_double(value.rational)
                                         : value.decimal;
        long long numerator = 0;
        long long denominator = 1;
        if (!mymath::best_rational_approximation(decimal_value,
                                                 &numerator,
                                                 &denominator,
                                                 max_denominator)) {
            throw std::runtime_error("rat could not compute a rational approximation");
        }

        stored.exact = true;
        stored.rational = Rational(numerator, denominator);
        stored.decimal = rational_to_double(stored.rational);
        return stored;
    }

    if (exact_mode) {
        try {
            ExactParser parser(trimmed, &variables, &impl->functions, has_script_function);
            stored.rational = parser.parse();
            stored.exact = true;
            stored.decimal = rational_to_double(stored.rational);
            return stored;
        } catch (const ExactModeUnsupported&) {
        }
    }

    matrix::Value matrix_value;
    if (try_evaluate_matrix_expression(trimmed,
                                       &variables,
                                       &impl->functions,
                                       has_script_function,
                                       invoke_script_function,
                                       &matrix_value)) {
        if (matrix_value.is_matrix) {
            stored.is_matrix = true;
            stored.matrix = matrix_value.matrix;
        } else {
            stored.decimal = matrix_value.scalar;
            stored.exact = false;
        }
        return stored;
    }

    if (!exact_mode) {
        try {
            PreciseDecimalParser parser(trimmed, &variables);
            const PreciseDecimal precise_value = parser.parse();
            stored.decimal = precise_value.to_double();
            stored.exact = false;
            stored.has_precise_decimal_text = true;
            stored.precise_decimal_text = precise_value.to_string();
            return stored;
        } catch (const PreciseDecimalUnsupported&) {
        }
    }

    DecimalParser parser(trimmed,
                         &variables,
                         &impl->functions,
                         has_script_function,
                         invoke_script_function);
    const double parsed_value = parser.parse();
    stored.decimal = parsed_value;
    stored.exact = false;
    if (impl->symbolic_constants_mode) {
        std::string symbolic_output;
        if (try_symbolic_constant_expression(trimmed,
                                             &variables,
                                             &impl->functions,
                                             &symbolic_output)) {
            stored.has_symbolic_text = true;
            stored.symbolic_text = symbolic_output;
        }
    }
    return stored;
}

std::string execute_simple_script_line(Calculator* calculator,
                                       Calculator::Impl* impl,
                                       const std::string& text,
                                       bool exact_mode) {
    std::vector<std::string> print_arguments;
    if (split_named_call_with_arguments(text, "print", &print_arguments)) {
        if (print_arguments.empty()) {
            throw std::runtime_error("print expects at least one argument");
        }
        std::ostringstream out;
        for (std::size_t i = 0; i < print_arguments.size(); ++i) {
            if (i != 0) {
                out << ' ';
            }
            out << format_print_value(
                evaluate_expression_value(calculator, impl, print_arguments[i], exact_mode),
                impl->symbolic_constants_mode);
        }
        return out.str();
    }

    std::string function_output;
    if (calculator->try_process_function_command(text, &function_output)) {
        return function_output;
    }

    std::string lhs;
    std::string rhs;
    if (split_assignment(text, &lhs, &rhs)) {
        if (!is_valid_variable_name(lhs)) {
            throw std::runtime_error("invalid variable name: " + lhs);
        }
        if (rhs.empty()) {
            throw std::runtime_error("assignment requires a value");
        }
        const StoredValue stored = evaluate_expression_value(calculator, impl, rhs, exact_mode);
        assign_visible_variable(impl, lhs, stored);
        return lhs + " = " + format_stored_value(stored, impl->symbolic_constants_mode);
    }

    return format_stored_value(evaluate_expression_value(calculator, impl, text, exact_mode),
                               impl->symbolic_constants_mode);
}

ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope) {
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return execute_script_block(calculator,
                                        impl,
                                        static_cast<const script::BlockStatement&>(statement),
                                        exact_mode,
                                        last_output,
                                        create_scope);
        case script::Statement::Kind::kSimple: {
            const auto& simple = static_cast<const script::SimpleStatement&>(statement);
            *last_output = execute_simple_script_line(calculator, impl, simple.text, exact_mode);
            return {};
        }
        case script::Statement::Kind::kIf: {
            const auto& if_statement = static_cast<const script::IfStatement&>(statement);
            if (truthy_value(evaluate_expression_value(calculator, impl, if_statement.condition, false))) {
                return execute_script_statement(calculator,
                                                impl,
                                                *if_statement.then_branch,
                                                exact_mode,
                                                last_output,
                                                true);
            }
            if (if_statement.else_branch) {
                return execute_script_statement(calculator,
                                                impl,
                                                *if_statement.else_branch,
                                                exact_mode,
                                                last_output,
                                                true);
            }
            return {};
        }
        case script::Statement::Kind::kWhile: {
            const auto& while_statement = static_cast<const script::WhileStatement&>(statement);
            while (truthy_value(evaluate_expression_value(calculator,
                                                          impl,
                                                          while_statement.condition,
                                                          false))) {
                const ScriptSignal signal =
                    execute_script_statement(calculator,
                                             impl,
                                             *while_statement.body,
                                             exact_mode,
                                             last_output,
                                             true);
                if (signal.kind == ScriptSignal::Kind::kReturn) {
                    return signal;
                }
                if (signal.kind == ScriptSignal::Kind::kBreak) {
                    break;
                }
                if (signal.kind == ScriptSignal::Kind::kContinue) {
                    continue;
                }
            }
            return {};
        }
        case script::Statement::Kind::kFor: {
            const auto& for_statement = static_cast<const script::ForStatement&>(statement);
            impl->local_scopes.push_back({});
            try {
                if (!for_statement.initializer.empty()) {
                    (void)execute_simple_script_line(calculator,
                                                     impl,
                                                     for_statement.initializer,
                                                     exact_mode);
                }
                while (for_statement.condition.empty() ||
                       truthy_value(evaluate_expression_value(calculator,
                                                              impl,
                                                              for_statement.condition,
                                                              false))) {
                    const ScriptSignal signal =
                        execute_script_statement(calculator,
                                                 impl,
                                                 *for_statement.body,
                                                 exact_mode,
                                                 last_output,
                                                 true);
                    if (signal.kind == ScriptSignal::Kind::kReturn) {
                        impl->local_scopes.pop_back();
                        return signal;
                    }
                    if (signal.kind == ScriptSignal::Kind::kBreak) {
                        break;
                    }
                    if (!for_statement.step.empty()) {
                        (void)execute_simple_script_line(calculator,
                                                         impl,
                                                         for_statement.step,
                                                         exact_mode);
                    }
                    if (signal.kind == ScriptSignal::Kind::kContinue) {
                        continue;
                    }
                }
                impl->local_scopes.pop_back();
                return {};
            } catch (...) {
                impl->local_scopes.pop_back();
                throw;
            }
        }
        case script::Statement::Kind::kFunction: {
            const auto& function_statement =
                static_cast<const script::FunctionStatement&>(statement);
            if (!is_valid_variable_name(function_statement.name)) {
                throw std::runtime_error("invalid function name: " + function_statement.name);
            }
            if (is_reserved_function_name(function_statement.name)) {
                throw std::runtime_error("function name is reserved: " +
                                         function_statement.name);
            }
            for (const std::string& parameter_name : function_statement.parameters) {
                if (!is_valid_variable_name(parameter_name)) {
                    throw std::runtime_error("invalid parameter name: " + parameter_name);
                }
            }
            ScriptFunction function;
            function.parameter_names = function_statement.parameters;
            function.body = std::shared_ptr<const script::BlockStatement>(
                clone_block_statement(*function_statement.body).release());
            impl->script_functions[function_statement.name] = function;
            *last_output = function_statement.name + "(...)";
            return {};
        }
        case script::Statement::Kind::kReturn: {
            const auto& return_statement =
                static_cast<const script::ReturnStatement&>(statement);
            if (!return_statement.has_expression) {
                ScriptSignal signal;
                signal.kind = ScriptSignal::Kind::kReturn;
                return signal;
            }
            return ScriptSignal::make_return(
                evaluate_expression_value(calculator, impl, return_statement.expression, exact_mode));
        }
        case script::Statement::Kind::kBreak:
            return ScriptSignal::make_break();
        case script::Statement::Kind::kContinue:
            return ScriptSignal::make_continue();
    }

    throw std::runtime_error("unknown script statement kind");
}

ScriptSignal execute_script_block(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope) {
    if (create_scope) {
        impl->local_scopes.push_back({});
    }

    try {
        for (const auto& statement : block.statements) {
            const ScriptSignal signal =
                execute_script_statement(calculator,
                                         impl,
                                         *statement,
                                         exact_mode,
                                         last_output,
                                         true);
            if (signal.kind != ScriptSignal::Kind::kNone) {
                if (create_scope) {
                    impl->local_scopes.pop_back();
                }
                return signal;
            }
        }
        if (create_scope) {
            impl->local_scopes.pop_back();
        }
        return {};
    } catch (...) {
        if (create_scope) {
            impl->local_scopes.pop_back();
        }
        throw;
    }
}

double Calculator::evaluate(const std::string& expression) {
    return normalize_result(evaluate_raw(expression));
}

double Calculator::evaluate_raw(const std::string& expression) {
    const StoredValue value = evaluate_expression_value(this, impl_.get(), expression, false);
    if (value.is_matrix) {
        throw std::runtime_error("matrix expression cannot be used as a scalar");
    }
    return value.decimal;
}

std::string Calculator::evaluate_for_display(const std::string& expression, bool exact_mode) {
    // 显示型功能优先于普通数值/分数显示，例如 hex(255) 应直接得到 "FF"。
    const std::map<std::string, StoredValue> variables = visible_variables(impl_.get());
    std::string converted;
    if (try_base_conversion_expression(expression,
                                       &variables,
                                       &impl_->functions,
                                       {impl_->hex_prefix_mode, impl_->hex_uppercase_mode},
                                       &converted)) {
        return converted;
    }

    if (impl_->symbolic_constants_mode) {
        std::string symbolic_output;
        if (try_symbolic_constant_expression(expression,
                                             &variables,
                                             &impl_->functions,
                                             &symbolic_output)) {
            return symbolic_output;
        }
    }

    return format_stored_value(
        evaluate_expression_value(this, impl_.get(), expression, exact_mode),
        impl_->symbolic_constants_mode);
}

std::string Calculator::process_line(const std::string& expression, bool exact_mode) {
    std::string lhs;
    std::string rhs;
    if (!split_assignment(expression, &lhs, &rhs)) {
        return evaluate_for_display(expression, exact_mode);
    }

    if (!is_valid_variable_name(lhs)) {
        throw std::runtime_error("invalid variable name: " + lhs);
    }
    if (rhs.empty()) {
        throw std::runtime_error("assignment requires a value");
    }

    const StoredValue stored = evaluate_expression_value(this, impl_.get(), rhs, exact_mode);
    assign_visible_variable(impl_.get(), lhs, stored);
    return lhs + " = " + format_stored_value(stored, impl_->symbolic_constants_mode);
}

std::string Calculator::execute_script(const std::string& source, bool exact_mode) {
    script::Program program = script::parse_program(source);
    std::string last_output;
    for (const auto& statement : program.statements) {
        const ScriptSignal signal =
            execute_script_statement(this, impl_.get(), *statement, exact_mode, &last_output, false);
        if (signal.kind == ScriptSignal::Kind::kReturn) {
            return signal.has_value ? format_stored_value(signal.value, impl_->symbolic_constants_mode)
                                    : (last_output.empty() ? "OK" : last_output);
        }
        if (signal.kind == ScriptSignal::Kind::kBreak ||
            signal.kind == ScriptSignal::Kind::kContinue) {
            throw std::runtime_error("break/continue can only be used inside loops");
        }
    }
    return last_output.empty() ? "OK" : last_output;
}

std::string Calculator::list_variables() const {
    if (impl_->variables.empty()) {
        return "No variables defined.";
    }

    std::ostringstream out;
    bool first = true;
    for (const auto& [name, value] : impl_->variables) {
        // std::map 保证变量按名字稳定排序，便于人读和测试断言。
        if (!first) {
            out << '\n';
        }
        first = false;
        out << name << " = " << format_stored_value(value, impl_->symbolic_constants_mode);
    }
    return out.str();
}

std::string Calculator::factor_expression(const std::string& expression) const {
    std::string inside;
    if (!split_named_call(expression, "factor", &inside)) {
        throw std::runtime_error("expected factor(expression)");
    }

    // 先允许 inside 是一个普通表达式或变量，再检查最终值是否为整数。
    DecimalParser parser(inside, &impl_->variables, &impl_->functions);
    const double value = normalize_result(parser.parse());
    if (!is_integer_double(value)) {
        throw std::runtime_error("factor only accepts integers");
    }

    return factor_integer(round_to_long_long(value));
}

std::string Calculator::base_conversion_expression(const std::string& expression) const {
    std::string converted;
    if (!try_base_conversion_expression(expression,
                                        &impl_->variables,
                                        &impl_->functions,
                                        {impl_->hex_prefix_mode, impl_->hex_uppercase_mode},
                                        &converted)) {
        throw std::runtime_error("expected bin(...), oct(...), hex(...), or base(value, base)");
    }
    return converted;
}

bool Calculator::try_process_function_command(const std::string& expression,
                                              std::string* output) {
    std::string function_name;
    std::string parameter_name;
    std::string body;
    if (split_function_definition(expression, &function_name, &parameter_name, &body)) {
        if (is_reserved_function_name(function_name)) {
            throw std::runtime_error("function name is reserved: " + function_name);
        }
        impl_->functions[function_name] = {parameter_name, body};
        *output = function_name + "(" + parameter_name + ") = " + body;
        return true;
    }

    const std::string trimmed = trim_copy(expression);
    if (trimmed == ":funcs") {
        if (impl_->functions.empty() && impl_->script_functions.empty()) {
            *output = "No custom functions defined.";
            return true;
        }

        std::ostringstream out;
        bool first = true;
        for (const auto& [name, function] : impl_->functions) {
            if (!first) {
                out << '\n';
            }
            first = false;
            out << name << "(" << function.parameter_name << ") = "
                << function.expression;
        }
        for (const auto& [name, function] : impl_->script_functions) {
            if (!first) {
                out << '\n';
            }
            first = false;
            out << name << "(";
            for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << function.parameter_names[i];
            }
            out << ") = { ... }";
        }
        *output = out.str();
        return true;
    }

    if (trimmed == ":clearfuncs") {
        impl_->functions.clear();
        impl_->script_functions.clear();
        *output = "Cleared all custom functions.";
        return true;
    }

    if (trimmed.rfind(":clearfunc ", 0) == 0) {
        const std::string name = trim_copy(trimmed.substr(11));
        const auto it = impl_->functions.find(name);
        if (it != impl_->functions.end()) {
            impl_->functions.erase(it);
            *output = "Cleared custom function: " + name;
            return true;
        }
        const auto script_it = impl_->script_functions.find(name);
        if (script_it == impl_->script_functions.end()) {
            throw std::runtime_error("unknown custom function: " + name);
        }
        impl_->script_functions.erase(script_it);
        *output = "Cleared custom function: " + name;
        return true;
    }

    auto build_symbolic_expression =
        [this](const std::string& name, std::string* variable_name) {
            const auto it = impl_->functions.find(name);
            if (it == impl_->functions.end()) {
                throw std::runtime_error("unknown custom function: " + name);
            }
            *variable_name = it->second.parameter_name;
            return SymbolicExpression::parse(it->second.expression);
        };

    std::function<void(const std::string&, std::string*, std::vector<double>*)>
        build_polynomial = [&](const std::string& argument,
                               std::string* variable_name,
                               std::vector<double>* coefficients) {
            const std::string trimmed_argument = trim_copy(argument);
            std::string nested_inside;
            if (split_named_call(trimmed_argument, "poly_add", &nested_inside) ||
                split_named_call(trimmed_argument, "poly_sub", &nested_inside) ||
                split_named_call(trimmed_argument, "poly_mul", &nested_inside) ||
                split_named_call(trimmed_argument, "poly_div", &nested_inside)) {
                const std::vector<std::string> nested_arguments =
                    split_top_level_arguments(nested_inside);
                if (nested_arguments.size() != 2) {
                    throw std::runtime_error(
                        "polynomial operations expect exactly two arguments");
                }

                std::string lhs_variable;
                std::string rhs_variable;
                std::vector<double> lhs_coefficients;
                std::vector<double> rhs_coefficients;
                build_polynomial(nested_arguments[0], &lhs_variable, &lhs_coefficients);
                build_polynomial(nested_arguments[1], &rhs_variable, &rhs_coefficients);

                *variable_name = lhs_variable;
                if (trimmed_argument.rfind("poly_add", 0) == 0) {
                    *coefficients = polynomial_add(lhs_coefficients, rhs_coefficients);
                    return;
                }
                if (trimmed_argument.rfind("poly_sub", 0) == 0) {
                    *coefficients = polynomial_subtract(lhs_coefficients, rhs_coefficients);
                    return;
                }
                if (trimmed_argument.rfind("poly_div", 0) == 0) {
                    const PolynomialDivisionResult division =
                        polynomial_divide(lhs_coefficients, rhs_coefficients);
                    bool zero_remainder = true;
                    for (double coefficient : division.remainder) {
                        if (!mymath::is_near_zero(coefficient, 1e-10)) {
                            zero_remainder = false;
                            break;
                        }
                    }
                    if (!zero_remainder) {
                        throw std::runtime_error(
                            "nested poly_div requires zero remainder");
                    }
                    *coefficients = division.quotient;
                    return;
                }

                *coefficients = polynomial_multiply(lhs_coefficients, rhs_coefficients);
                return;
            }

            // 命令层不直接处理原始字符串，而是统一走：
            // 自定义函数名 -> 符号表达式 -> 多项式系数
            //
            // 这样 poly_add/poly_mul/roots 等功能只依赖一个稳定的中间表示，
            // 不需要各自重复做“它是不是多项式”的判断。
            SymbolicExpression expression =
                build_symbolic_expression(trimmed_argument, variable_name);
            if (!expression.polynomial_coefficients(*variable_name, coefficients)) {
                throw std::runtime_error("custom function " + trimmed_argument +
                                         " is not a polynomial");
            }
        };

    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)>
        resolve_symbolic_expression =
            [&](const std::string& argument,
                bool require_single_variable,
                std::string* variable_name,
                SymbolicExpression* expression) {
                const std::string trimmed_argument = trim_copy(argument);
                std::string nested_inside;
                if (split_named_call(trimmed_argument, "diff", &nested_inside)) {
                    const std::vector<std::string> nested_arguments =
                        split_top_level_arguments(nested_inside);
                    if (nested_arguments.empty()) {
                        throw std::runtime_error(
                            "nested symbolic diff expects at least one argument");
                    }
                    SymbolicExpression nested_expression;
                    resolve_symbolic_expression(nested_arguments[0],
                                                nested_arguments.size() == 1,
                                                variable_name,
                                                &nested_expression);
                    if (nested_arguments.size() == 1) {
                        *expression =
                            nested_expression.derivative(*variable_name).simplify();
                    } else {
                        SymbolicExpression differentiated = nested_expression;
                        for (std::size_t i = 1; i < nested_arguments.size(); ++i) {
                            const std::string derivative_variable =
                                trim_copy(nested_arguments[i]);
                            if (!is_identifier_text(derivative_variable)) {
                                throw std::runtime_error(
                                    "nested symbolic diff variable arguments must be identifiers");
                            }
                            differentiated =
                                differentiated.derivative(derivative_variable).simplify();
                        }
                        *variable_name = trim_copy(nested_arguments[1]);
                        *expression = differentiated;
                    }
                    return;
                }
                if (split_named_call(trimmed_argument, "integral", &nested_inside)) {
                    const std::vector<std::string> nested_arguments =
                        split_top_level_arguments(nested_inside);
                    if (nested_arguments.size() != 1 && nested_arguments.size() != 2) {
                        throw std::runtime_error(
                            "nested symbolic integral expects expression and optional variable");
                    }
                    SymbolicExpression nested_expression;
                    resolve_symbolic_expression(nested_arguments[0],
                                                nested_arguments.size() == 1,
                                                variable_name,
                                                &nested_expression);
                    if (nested_arguments.size() == 2) {
                        const std::string integral_variable =
                            trim_copy(nested_arguments[1]);
                        if (!is_identifier_text(integral_variable)) {
                            throw std::runtime_error(
                                "nested symbolic integral variable must be an identifier");
                        }
                        *variable_name = integral_variable;
                    }
                    *expression =
                        nested_expression.integral(*variable_name).simplify();
                    return;
                }
                if (split_named_call(trimmed_argument, "poly_add", &nested_inside) ||
                    split_named_call(trimmed_argument, "poly_sub", &nested_inside) ||
                    split_named_call(trimmed_argument, "poly_mul", &nested_inside) ||
                    split_named_call(trimmed_argument, "poly_div", &nested_inside)) {
                    std::vector<double> coefficients;
                    build_polynomial(trimmed_argument, variable_name, &coefficients);
                    *expression = SymbolicExpression::parse(
                        polynomial_to_string(coefficients, *variable_name));
                    return;
                }
                const auto function_it = impl_->functions.find(trimmed_argument);
                if (function_it != impl_->functions.end()) {
                    *expression = build_symbolic_expression(trimmed_argument, variable_name);
                    return;
                }
                *expression = SymbolicExpression::parse(
                    expand_inline_function_commands(this, trimmed_argument));
                const std::vector<std::string> identifiers =
                    expression->identifier_variables();
                if (identifiers.size() == 1) {
                    *variable_name = identifiers[0];
                } else if (identifiers.empty()) {
                    *variable_name = "x";
                } else if (require_single_variable) {
                    throw std::runtime_error(
                        "symbolic expressions with multiple variables must use a custom function");
                } else {
                    *variable_name = identifiers.front();
                }
            };

    auto build_analysis = [&](const std::string& argument) {
        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(argument, true, &variable_name, &expression);
        FunctionAnalysis analysis(variable_name);
        analysis.define(expression.to_string());
        return analysis;
    };

    auto build_scoped_decimal_evaluator = [&](const std::string& argument) {
        const std::string scoped_expression =
            trim_copy(expand_inline_function_commands(this, argument));
        return [this, scoped_expression](
                   const std::vector<std::pair<std::string, double>>& assignments) {
            std::map<std::string, StoredValue> scoped_variables =
                visible_variables(impl_.get());
            for (const auto& [name, value] : assignments) {
                StoredValue stored;
                stored.decimal = normalize_display_decimal(value);
                stored.exact = false;
                scoped_variables[name] = stored;
            }

            const HasScriptFunctionCallback has_script_function =
                [this](const std::string& name) {
                    return has_visible_script_function(impl_.get(), name);
                };
            const InvokeScriptFunctionDecimalCallback invoke_script_function =
                [this](const std::string& name,
                       const std::vector<double>& arguments) {
                    return invoke_script_function_decimal(
                        this, impl_.get(), name, arguments);
                };

            DecimalParser parser(scoped_expression,
                                 &scoped_variables,
                                 &impl_->functions,
                                 has_script_function,
                                 invoke_script_function);
            return normalize_result(parser.parse());
        };
    };

    auto build_scoped_matrix_evaluator = [&](const std::string& argument) {
        const std::string scoped_expression =
            trim_copy(expand_inline_function_commands(this, argument));
        return [this, scoped_expression](
                   const std::vector<std::pair<std::string, StoredValue>>& assignments) {
            std::map<std::string, StoredValue> scoped_variables =
                visible_variables(impl_.get());
            for (const auto& [name, value] : assignments) {
                scoped_variables[name] = value;
            }

            const HasScriptFunctionCallback has_script_function =
                [this](const std::string& name) {
                    return has_visible_script_function(impl_.get(), name);
                };
            const InvokeScriptFunctionDecimalCallback invoke_script_function =
                [this](const std::string& name,
                       const std::vector<double>& arguments) {
                    return invoke_script_function_decimal(
                        this, impl_.get(), name, arguments);
                };

            matrix::Value value;
            if (!try_evaluate_matrix_expression(scoped_expression,
                                               &scoped_variables,
                                               &impl_->functions,
                                               has_script_function,
                                               invoke_script_function,
                                               &value) ||
                !value.is_matrix) {
                throw std::runtime_error("expected a matrix-valued expression");
            }
            return value.matrix;
        };
    };

    auto build_scoped_scalar_evaluator = [&](const std::string& argument) {
        const std::string scoped_expression =
            trim_copy(expand_inline_function_commands(this, argument));
        return [this, scoped_expression](
                   const std::vector<std::pair<std::string, StoredValue>>& assignments) {
            std::map<std::string, StoredValue> frame;
            for (const auto& [name, value] : assignments) {
                frame[name] = value;
            }

            impl_->local_scopes.push_back(frame);
            try {
                const StoredValue value =
                    evaluate_expression_value(this, impl_.get(), scoped_expression, false);
                impl_->local_scopes.pop_back();

                if (value.is_matrix || value.is_string) {
                    throw std::runtime_error("expected a scalar-valued expression");
                }
                return normalize_result(value.exact
                                            ? rational_to_double(value.rational)
                                            : value.decimal);
            } catch (...) {
                impl_->local_scopes.pop_back();
                throw;
            }
        };
    };

    auto parse_decimal_argument = [&](const std::string& argument) {
        DecimalParser parser(argument, &impl_->variables, &impl_->functions);
        return parser.parse();
    };

    auto parse_matrix_argument = [&](const std::string& argument,
                                     const std::string& context) {
        const StoredValue value =
            evaluate_expression_value(this, impl_.get(), argument, false);
        if (!value.is_matrix) {
            throw std::runtime_error(context + " expects a matrix or vector argument");
        }
        return value.matrix;
    };

    auto is_matrix_argument = [&](const std::string& argument) {
        const std::map<std::string, StoredValue> visible = visible_variables(impl_.get());
        const HasScriptFunctionCallback has_script_function =
            [this](const std::string& name) {
                return has_visible_script_function(impl_.get(), name);
            };
        const InvokeScriptFunctionDecimalCallback invoke_script_function =
            [this](const std::string& name, const std::vector<double>& arguments) {
                return invoke_script_function_decimal(this, impl_.get(), name, arguments);
            };
        matrix::Value value;
        return try_evaluate_matrix_expression(trim_copy(argument),
                                             &visible,
                                             &impl_->functions,
                                             has_script_function,
                                             invoke_script_function,
                                             &value) &&
               value.is_matrix;
    };

    auto matrix_to_vector_values = [&](const matrix::Matrix& value,
                                       const std::string& context) {
        if (!value.is_vector()) {
            throw std::runtime_error(context + " expects vector arguments");
        }
        const std::size_t size = value.rows == 1 ? value.cols : value.rows;
        std::vector<double> result(size, 0.0);
        for (std::size_t i = 0; i < size; ++i) {
            result[i] = value.rows == 1 ? value.at(0, i) : value.at(i, 0);
        }
        return result;
    };

    auto vector_to_column_matrix = [&](const std::vector<double>& values) {
        matrix::Matrix result(values.size(), 1, 0.0);
        for (std::size_t i = 0; i < values.size(); ++i) {
            result.at(i, 0) = normalize_result(values[i]);
        }
        return result;
    };

    auto make_scalar_stored = [&](double value) {
        StoredValue stored;
        stored.decimal = normalize_result(value);
        return stored;
    };

    auto append_parameter_assignments =
        [&](const StoredValue& parameter_value,
            std::vector<std::pair<std::string, StoredValue>>* assignments) {
            assignments->push_back({"p", parameter_value});
            if (!parameter_value.is_matrix || !parameter_value.matrix.is_vector()) {
                return;
            }

            const std::size_t size =
                parameter_value.matrix.rows == 1
                    ? parameter_value.matrix.cols
                    : parameter_value.matrix.rows;
            for (std::size_t i = 0; i < size; ++i) {
                const double component_value =
                    parameter_value.matrix.rows == 1
                        ? parameter_value.matrix.at(0, i)
                        : parameter_value.matrix.at(i, 0);
                assignments->push_back({"p" + std::to_string(i + 1),
                                        make_scalar_stored(component_value)});
            }
        };

    auto try_parse_positive_step_argument = [&](const std::string& argument,
                                                int* steps) {
        try {
            const double value = parse_decimal_argument(argument);
            if (!is_integer_double(value) || value <= 0.0) {
                return false;
            }
            *steps = static_cast<int>(round_to_long_long(value));
            return true;
        } catch (const std::exception&) {
            return false;
        }
    };

    auto dot_product = [](const std::vector<double>& lhs,
                          const std::vector<double>& rhs) {
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error("objective and solution dimension mismatch");
        }
        long double total = 0.0L;
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            total += static_cast<long double>(lhs[i]) *
                     static_cast<long double>(rhs[i]);
        }
        return static_cast<double>(total);
    };

    auto format_planning_result = [&](const std::vector<double>& solution,
                                      double objective) {
        return "x = " + matrix::Matrix::vector(solution).to_string() +
               "\nobjective = " + format_decimal(normalize_result(objective));
    };

    auto solve_linear_box_problem =
        [&](const std::vector<double>& objective,
            const matrix::Matrix& inequality_matrix,
            const std::vector<double>& inequality_rhs,
            const matrix::Matrix& equality_matrix,
            const std::vector<double>& equality_rhs,
            const std::vector<double>& lower_bounds,
            const std::vector<double>& upper_bounds,
            double planning_tolerance,
            std::vector<double>* solution,
            double* objective_value) {
            const std::size_t variable_count = objective.size();
            if (inequality_matrix.cols != variable_count ||
                inequality_rhs.size() != inequality_matrix.rows ||
                equality_matrix.cols != variable_count ||
                equality_rhs.size() != equality_matrix.rows ||
                lower_bounds.size() != variable_count ||
                upper_bounds.size() != variable_count) {
                throw std::runtime_error("planning dimension mismatch");
            }

            auto feasible_solution = [&](const std::vector<double>& x) {
                if (x.size() != variable_count) {
                    return false;
                }
                for (std::size_t col = 0; col < variable_count; ++col) {
                    if (x[col] < lower_bounds[col] - planning_tolerance ||
                        x[col] > upper_bounds[col] + planning_tolerance) {
                        return false;
                    }
                }
                for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                    long double total = 0.0L;
                    for (std::size_t col = 0; col < variable_count; ++col) {
                        total += static_cast<long double>(inequality_matrix.at(row, col)) *
                                 static_cast<long double>(x[col]);
                    }
                    if (total >
                        static_cast<long double>(inequality_rhs[row]) + planning_tolerance) {
                        return false;
                    }
                }
                for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                    long double total = 0.0L;
                    for (std::size_t col = 0; col < variable_count; ++col) {
                        total += static_cast<long double>(equality_matrix.at(row, col)) *
                                 static_cast<long double>(x[col]);
                    }
                    if (mymath::abs(static_cast<double>(total - equality_rhs[row])) >
                        planning_tolerance) {
                        return false;
                    }
                }
                return true;
            };

            if (variable_count == 0) {
                const std::vector<double> empty_solution;
                if (!feasible_solution(empty_solution)) {
                    return false;
                }
                *solution = empty_solution;
                *objective_value = 0.0;
                return true;
            }

            std::vector<std::vector<double>> all_rows;
            std::vector<double> all_rhs;
            all_rows.reserve(inequality_matrix.rows + 2 * equality_matrix.rows +
                             2 * variable_count);
            all_rhs.reserve(inequality_matrix.rows + 2 * equality_matrix.rows +
                            2 * variable_count);

            for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                std::vector<double> coefficients(variable_count, 0.0);
                for (std::size_t col = 0; col < variable_count; ++col) {
                    coefficients[col] = inequality_matrix.at(row, col);
                }
                all_rows.push_back(coefficients);
                all_rhs.push_back(inequality_rhs[row]);
            }
            for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                std::vector<double> positive_row(variable_count, 0.0);
                std::vector<double> negative_row(variable_count, 0.0);
                for (std::size_t col = 0; col < variable_count; ++col) {
                    positive_row[col] = equality_matrix.at(row, col);
                    negative_row[col] = -equality_matrix.at(row, col);
                }
                all_rows.push_back(positive_row);
                all_rhs.push_back(equality_rhs[row]);
                all_rows.push_back(negative_row);
                all_rhs.push_back(-equality_rhs[row]);
            }
            for (std::size_t col = 0; col < variable_count; ++col) {
                std::vector<double> upper_row(variable_count, 0.0);
                upper_row[col] = 1.0;
                all_rows.push_back(upper_row);
                all_rhs.push_back(upper_bounds[col]);

                std::vector<double> lower_row(variable_count, 0.0);
                lower_row[col] = -1.0;
                all_rows.push_back(lower_row);
                all_rhs.push_back(-lower_bounds[col]);
            }

            bool found = false;
            std::vector<double> best_solution(variable_count, 0.0);
            double best_value = 0.0;
            std::vector<std::size_t> selection(variable_count, 0);

            std::function<void(std::size_t, std::size_t)> enumerate_bases =
                [&](std::size_t start, std::size_t depth) {
                    if (depth == variable_count) {
                        matrix::Matrix basis(variable_count, variable_count, 0.0);
                        matrix::Matrix basis_rhs(variable_count, 1, 0.0);
                        for (std::size_t row = 0; row < variable_count; ++row) {
                            basis_rhs.at(row, 0) = all_rhs[selection[row]];
                            for (std::size_t col = 0; col < variable_count; ++col) {
                                basis.at(row, col) = all_rows[selection[row]][col];
                            }
                        }

                        try {
                            const matrix::Matrix solved = matrix::solve(basis, basis_rhs);
                            std::vector<double> candidate(variable_count, 0.0);
                            for (std::size_t i = 0; i < variable_count; ++i) {
                                const double value = solved.at(i, 0);
                                candidate[i] =
                                    mymath::abs(value) <= planning_tolerance ? 0.0 : value;
                            }
                            if (!feasible_solution(candidate)) {
                                return;
                            }

                            const double candidate_objective =
                                dot_product(objective, candidate);
                            if (!found ||
                                candidate_objective > best_value + planning_tolerance) {
                                found = true;
                                best_value = candidate_objective;
                                best_solution = candidate;
                            }
                        } catch (const std::exception&) {
                        }
                        return;
                    }

                    const std::size_t remaining = variable_count - depth;
                    for (std::size_t index = start;
                         index + remaining <= all_rows.size();
                         ++index) {
                        selection[depth] = index;
                        enumerate_bases(index + 1, depth + 1);
                    }
                };

            enumerate_bases(0, 0);
            if (!found) {
                return false;
            }

            *solution = best_solution;
            *objective_value = best_value;
            return true;
        };

    auto parse_subdivisions = [&](const std::vector<std::string>& arguments,
                                  std::size_t offset,
                                  const std::vector<int>& defaults) {
        std::vector<int> subdivisions = defaults;
        if (arguments.size() == offset) {
            return subdivisions;
        }
        if (arguments.size() != offset + defaults.size()) {
            throw std::runtime_error("unexpected subdivision argument count");
        }

        for (std::size_t i = 0; i < defaults.size(); ++i) {
            const double value = parse_decimal_argument(arguments[offset + i]);
            if (!is_integer_double(value) || value <= 0.0) {
                throw std::runtime_error(
                    "integration subdivision counts must be positive integers");
            }
            subdivisions[i] = static_cast<int>(round_to_long_long(value));
        }
        return subdivisions;
    };

    auto evaluate_symbolic_at = [&](const SymbolicExpression& expression,
                                    const std::string& variable_name,
                                    double point) {
        // Taylor 展开需要在某个点反复计算“当前导函数”的数值。
        // 为了复用现有 Calculator::evaluate，这里临时把变量写入会话变量表，
        // 求值完成后再恢复原状态，避免污染用户上下文。
        const auto existing = impl_->variables.find(variable_name);
        const bool had_existing = existing != impl_->variables.end();
        StoredValue backup;
        if (had_existing) {
            backup = existing->second;
        }

        StoredValue temporary;
        temporary.decimal = point;
        temporary.exact = false;
        impl_->variables[variable_name] = temporary;

        try {
            const double value = evaluate(expression.to_string());
            if (had_existing) {
                impl_->variables[variable_name] = backup;
            } else {
                impl_->variables.erase(variable_name);
            }
            return value;
        } catch (...) {
            if (had_existing) {
                impl_->variables[variable_name] = backup;
            } else {
                impl_->variables.erase(variable_name);
            }
            throw;
        }
    };

    auto build_taylor_coefficients = [&](const SymbolicExpression& expression,
                                         const std::string& variable_name,
                                         double center,
                                         int degree) {
        std::vector<double> coefficients;
        coefficients.reserve(static_cast<std::size_t>(degree + 1));
        SymbolicExpression current = expression;
        for (int order = 0; order <= degree; ++order) {
            const double derivative_value =
                evaluate_symbolic_at(current, variable_name, center);
            coefficients.push_back(derivative_value / factorial_int(order));
            if (order != degree) {
                current = current.derivative(variable_name).simplify();
            }
        }
        return coefficients;
    };

    auto simplify_symbolic_text = [&](const std::string& text) {
        return SymbolicExpression::parse(text).simplify().to_string();
    };

    auto symbolic_vector_to_string =
        [](const std::vector<SymbolicExpression>& values) {
            std::ostringstream out;
            out << "[";
            for (std::size_t i = 0; i < values.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << values[i].simplify().to_string();
            }
            out << "]";
            return out.str();
        };

    auto symbolic_matrix_to_string =
        [](const std::vector<std::vector<SymbolicExpression>>& values) {
            std::ostringstream out;
            out << "[";
            for (std::size_t row = 0; row < values.size(); ++row) {
                if (row != 0) {
                    out << ", ";
                }
                out << "[";
                for (std::size_t col = 0; col < values[row].size(); ++col) {
                    if (col != 0) {
                        out << ", ";
                    }
                    out << values[row][col].simplify().to_string();
                }
                out << "]";
            }
            out << "]";
            return out.str();
        };

    auto parse_symbolic_variable_arguments =
        [](const std::vector<std::string>& arguments,
           std::size_t start_index,
           const std::vector<std::string>& fallback_variables) {
            std::vector<std::string> variables;
            for (std::size_t i = start_index; i < arguments.size(); ++i) {
                const std::string variable = trim_copy(arguments[i]);
                if (!is_identifier_text(variable)) {
                    throw std::runtime_error(
                        "symbolic variable arguments must be identifiers");
                }
                variables.push_back(variable);
            }
            if (variables.empty()) {
                variables = fallback_variables;
            }
            if (variables.empty()) {
                variables.push_back("x");
            }
            return variables;
        };

    auto parse_symbolic_expression_list = [&](const std::string& argument) {
        std::string text = trim_copy(argument);
        if (text.size() < 2 || text.front() != '[' || text.back() != ']') {
            throw std::runtime_error(
                "jacobian expects its first argument to be a bracketed expression list");
        }
        text = trim_copy(text.substr(1, text.size() - 2));

        std::vector<std::string> expression_texts;
        int paren_depth = 0;
        int bracket_depth = 0;
        std::size_t start = 0;
        for (std::size_t i = 0; i < text.size(); ++i) {
            const char ch = text[i];
            if (ch == '(') {
                ++paren_depth;
            } else if (ch == ')') {
                --paren_depth;
            } else if (ch == '[') {
                ++bracket_depth;
            } else if (ch == ']') {
                --bracket_depth;
            } else if ((ch == ';' || ch == ',') &&
                       paren_depth == 0 &&
                       bracket_depth == 0) {
                expression_texts.push_back(trim_copy(text.substr(start, i - start)));
                start = i + 1;
            }
        }
        if (!text.empty()) {
            expression_texts.push_back(trim_copy(text.substr(start)));
        }
        if (expression_texts.empty()) {
            throw std::runtime_error("jacobian expression list cannot be empty");
        }

        std::vector<SymbolicExpression> expressions;
        expressions.reserve(expression_texts.size());
        for (const std::string& expression_text : expression_texts) {
            if (expression_text.empty()) {
                throw std::runtime_error("jacobian expression list cannot contain empty items");
            }
            expressions.push_back(SymbolicExpression::parse(
                expand_inline_function_commands(this, expression_text)));
        }
        return expressions;
    };

    std::string inside;
    if (split_named_call(trimmed, "poly_add", &inside) ||
        split_named_call(trimmed, "poly_sub", &inside) ||
        split_named_call(trimmed, "poly_mul", &inside) ||
        split_named_call(trimmed, "poly_div", &inside)) {
        // 这些命令只处理“已经定义好的多项式函数名”，
        // 例如 poly_mul(p, q)，而不是直接吃任意表达式。
        // 这样接口更清晰，也便于后续缓存或扩展更多多项式算法。
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 2) {
            throw std::runtime_error("polynomial operations expect exactly two arguments");
        }

        std::string lhs_variable;
        std::string rhs_variable;
        std::vector<double> lhs_coefficients;
        std::vector<double> rhs_coefficients;
        build_polynomial(arguments[0], &lhs_variable, &lhs_coefficients);
        build_polynomial(arguments[1], &rhs_variable, &rhs_coefficients);

        const std::string output_variable = lhs_variable;
        if (trimmed.rfind("poly_add", 0) == 0) {
            *output = polynomial_to_string(
                polynomial_add(lhs_coefficients, rhs_coefficients),
                output_variable);
            return true;
        }
        if (trimmed.rfind("poly_sub", 0) == 0) {
            *output = polynomial_to_string(
                polynomial_subtract(lhs_coefficients, rhs_coefficients),
                output_variable);
            return true;
        }
        if (trimmed.rfind("poly_mul", 0) == 0) {
            *output = polynomial_to_string(
                polynomial_multiply(lhs_coefficients, rhs_coefficients),
                output_variable);
            return true;
        }

        const PolynomialDivisionResult division =
            polynomial_divide(lhs_coefficients, rhs_coefficients);
        *output = "quotient: " +
                  polynomial_to_string(division.quotient, output_variable) +
                  ", remainder: " +
                  polynomial_to_string(division.remainder, output_variable);
        return true;
    }

    if (split_named_call(trimmed, "roots", &inside)) {
        // 求根目前限定为“多项式实根”。
        // 如果用户传进来的自定义函数不能化成多项式，会在 build_polynomial
        // 阶段直接报错，避免对非多项式做错误的根求解。
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 1) {
            throw std::runtime_error("roots expects exactly one argument");
        }

        std::string variable_name;
        std::vector<double> coefficients;
        build_polynomial(arguments[0], &variable_name, &coefficients);
        const std::vector<double> roots = polynomial_real_roots(coefficients);
        if (roots.empty()) {
            *output = "No real roots.";
            return true;
        }

        std::ostringstream out;
        for (std::size_t i = 0; i < roots.size(); ++i) {
            if (i != 0) {
                out << ", ";
            }
            const double root =
                is_integer_double(roots[i], 1e-6)
                    ? static_cast<double>(round_to_long_long(roots[i]))
                    : roots[i];
            out << format_symbolic_scalar(root);
        }
        *output = out.str();
        return true;
    }

    if (split_named_call(trimmed, "taylor", &inside)) {
        // Taylor 多项式按定义生成：
        //   f(a) + f'(a)(x-a) + f''(a)/2!(x-a)^2 + ...
        //
        // 这里采取一个务实实现：
        // - 导数本身走符号求导
        // - 每一阶导数在 a 处的值走现有数值求值
        //
        // 这样既能得到可读的展开式，也不需要再单独实现一套
        // “高阶导数闭式系数提取器”。
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3) {
            throw std::runtime_error("taylor expects exactly three arguments");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], true, &variable_name, &expression);
        DecimalParser center_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser degree_parser(arguments[2], &impl_->variables, &impl_->functions);
        const double center = center_parser.parse();
        const double degree_value = degree_parser.parse();
        if (!is_integer_double(degree_value) || degree_value < 0.0) {
            throw std::runtime_error("taylor degree must be a non-negative integer");
        }

        const int degree = static_cast<int>(round_to_long_long(degree_value));
        const std::vector<double> coefficients =
            build_taylor_coefficients(expression, variable_name, center, degree);
        *output = taylor_series_to_string(coefficients, variable_name, center);
        return true;
    }

    if (split_named_call(trimmed, "pade", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "pade expects expr, m, n or expr, center, m, n");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], true, &variable_name, &expression);

        const bool explicit_center = arguments.size() == 4;
        const double center = explicit_center
                                  ? parse_decimal_argument(arguments[1])
                                  : 0.0;
        const double numerator_degree_value = parse_decimal_argument(
            arguments[explicit_center ? 2 : 1]);
        const double denominator_degree_value = parse_decimal_argument(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(numerator_degree_value) ||
            numerator_degree_value < 0.0 ||
            !is_integer_double(denominator_degree_value) ||
            denominator_degree_value < 0.0) {
            throw std::runtime_error(
                "pade degrees must be non-negative integers");
        }

        const int numerator_degree =
            static_cast<int>(round_to_long_long(numerator_degree_value));
        const int denominator_degree =
            static_cast<int>(round_to_long_long(denominator_degree_value));
        if (numerator_degree == 0 && denominator_degree == 0) {
            throw std::runtime_error("pade requires at least one non-zero degree");
        }

        const std::vector<double> coefficients = build_taylor_coefficients(
            expression,
            variable_name,
            center,
            numerator_degree + denominator_degree);
        auto coefficient_at = [&](int index) {
            if (index < 0 ||
                index >= static_cast<int>(coefficients.size())) {
                return 0.0;
            }
            return coefficients[static_cast<std::size_t>(index)];
        };

        std::vector<double> denominator(denominator_degree + 1, 0.0);
        denominator[0] = 1.0;
        if (denominator_degree > 0) {
            std::vector<std::vector<double>> matrix(
                static_cast<std::size_t>(denominator_degree),
                std::vector<double>(static_cast<std::size_t>(denominator_degree), 0.0));
            std::vector<double> rhs(static_cast<std::size_t>(denominator_degree), 0.0);
            for (int row = 0; row < denominator_degree; ++row) {
                for (int col = 0; col < denominator_degree; ++col) {
                    matrix[static_cast<std::size_t>(row)]
                          [static_cast<std::size_t>(col)] =
                        coefficient_at(numerator_degree + row - col);
                }
                rhs[static_cast<std::size_t>(row)] =
                    -coefficient_at(numerator_degree + row + 1);
            }
            const std::vector<double> solved = solve_dense_linear_system(
                matrix, rhs, "pade");
            for (int i = 0; i < denominator_degree; ++i) {
                denominator[static_cast<std::size_t>(i + 1)] =
                    solved[static_cast<std::size_t>(i)];
            }
        }

        std::vector<double> numerator(numerator_degree + 1, 0.0);
        for (int i = 0; i <= numerator_degree; ++i) {
            double value = 0.0;
            for (int j = 0; j <= denominator_degree && j <= i; ++j) {
                value += denominator[static_cast<std::size_t>(j)] *
                         coefficient_at(i - j);
            }
            numerator[static_cast<std::size_t>(i)] = value;
        }

        const std::string base = shifted_series_base(variable_name, center);
        const std::string numerator_text =
            polynomial_to_string(numerator, base);
        const std::string denominator_text =
            polynomial_to_string(denominator, base);
        if (denominator_text == "1") {
            *output = simplify_symbolic_text(numerator_text);
        } else {
            *output = simplify_symbolic_text(
                "(" + numerator_text + ") / (" + denominator_text + ")");
        }
        return true;
    }

    if (split_named_call(trimmed, "puiseux", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "puiseux expects expr, degree, denominator or expr, center, degree, denominator");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], true, &variable_name, &expression);

        const bool explicit_center = arguments.size() == 4;
        const double center = explicit_center
                                  ? parse_decimal_argument(arguments[1])
                                  : 0.0;
        const double degree_value = parse_decimal_argument(
            arguments[explicit_center ? 2 : 1]);
        const double denominator_value = parse_decimal_argument(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(degree_value) || degree_value < 0.0) {
            throw std::runtime_error(
                "puiseux degree must be a non-negative integer");
        }
        if (!is_integer_double(denominator_value) || denominator_value <= 0.0) {
            throw std::runtime_error(
                "puiseux denominator must be a positive integer");
        }

        const int degree = static_cast<int>(round_to_long_long(degree_value));
        const int denominator =
            static_cast<int>(round_to_long_long(denominator_value));
        const std::string auxiliary_variable = "puiseux_t";
        const std::string replacement_text =
            mymath::is_near_zero(center, 1e-10)
                ? auxiliary_variable + " ^ " + std::to_string(denominator)
                : format_symbolic_scalar(center) + " + " +
                      auxiliary_variable + " ^ " +
                      std::to_string(denominator);
        const SymbolicExpression substituted = expression.substitute(
            variable_name,
            SymbolicExpression::parse(replacement_text));
        const std::vector<double> coefficients = build_taylor_coefficients(
            substituted, auxiliary_variable, 0.0, degree);
        *output = generalized_series_to_string(
            coefficients, variable_name, center, denominator);
        return true;
    }

    if (split_named_call(trimmed, "series_sum", &inside) ||
        split_named_call(trimmed, "summation", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 4) {
            throw std::runtime_error(
                "series_sum expects expr, index, lower, upper");
        }

        const std::string index_name = trim_copy(arguments[1]);
        if (!is_identifier_text(index_name)) {
            throw std::runtime_error("series_sum index must be an identifier");
        }

        SymbolicExpression summand =
            SymbolicExpression::parse(expand_inline_function_commands(this, arguments[0]));
        SymbolicExpression upper_expression;
        SymbolicExpression lower_expression;
        const std::string upper_text = trim_copy(arguments[3]);
        const bool upper_is_infinite =
            upper_text == "inf" || upper_text == "oo" ||
            upper_text == "infinity";
        if (!upper_is_infinite) {
            upper_expression = SymbolicExpression::parse(
                expand_inline_function_commands(this, upper_text));
        }
        lower_expression = SymbolicExpression::parse(
            expand_inline_function_commands(this, arguments[2]));

        auto make_polynomial_sum_primitive =
            [&](const std::vector<double>& coefficients) {
                if (coefficients.size() > 4) {
                    throw std::runtime_error(
                        "series_sum polynomial summands are currently supported up to degree 3");
                }

                std::vector<std::string> pieces;
                if (coefficients.size() >= 1 &&
                    !mymath::is_near_zero(coefficients[0], 1e-10)) {
                    pieces.push_back("(" + format_symbolic_scalar(coefficients[0]) +
                                     ") * (" + index_name + " + 1)");
                }
                if (coefficients.size() >= 2 &&
                    !mymath::is_near_zero(coefficients[1], 1e-10)) {
                    pieces.push_back("(" + format_symbolic_scalar(coefficients[1]) +
                                     ") * (" + index_name + " * (" + index_name +
                                     " + 1) / 2)");
                }
                if (coefficients.size() >= 3 &&
                    !mymath::is_near_zero(coefficients[2], 1e-10)) {
                    pieces.push_back("(" + format_symbolic_scalar(coefficients[2]) +
                                     ") * (" + index_name + " * (" + index_name +
                                     " + 1) * (2 * " + index_name + " + 1) / 6)");
                }
                if (coefficients.size() >= 4 &&
                    !mymath::is_near_zero(coefficients[3], 1e-10)) {
                    pieces.push_back("(" + format_symbolic_scalar(coefficients[3]) +
                                     ") * ((" + index_name + " * (" + index_name +
                                     " + 1) / 2) ^ 2)");
                }

                if (pieces.empty()) {
                    return SymbolicExpression::number(0.0);
                }
                std::ostringstream out;
                for (std::size_t i = 0; i < pieces.size(); ++i) {
                    if (i != 0) {
                        out << " + ";
                    }
                    out << pieces[i];
                }
                return SymbolicExpression::parse(out.str()).simplify();
            };

        auto finite_sum_from_primitive =
            [&](const SymbolicExpression& primitive) {
                const SymbolicExpression lower_minus_one =
                    SymbolicExpression::parse("(" + arguments[2] + ") - 1").simplify();
                return SymbolicExpression::parse(
                           "(" +
                           primitive.substitute(index_name, upper_expression)
                               .to_string() +
                           ") - (" +
                           primitive.substitute(index_name, lower_minus_one)
                               .to_string() +
                           ")")
                    .simplify()
                    .to_string();
            };

        std::vector<double> polynomial_coefficients;
        if (summand.polynomial_coefficients(index_name, &polynomial_coefficients)) {
            if (upper_is_infinite) {
                bool all_zero = true;
                for (double coefficient : polynomial_coefficients) {
                    if (!mymath::is_near_zero(coefficient, 1e-10)) {
                        all_zero = false;
                        break;
                    }
                }
                if (!all_zero) {
                    throw std::runtime_error(
                        "series_sum does not support infinite polynomial sums");
                }
                *output = "0";
                return true;
            }

            const SymbolicExpression primitive =
                make_polynomial_sum_primitive(polynomial_coefficients);
            *output = finite_sum_from_primitive(primitive);
            return true;
        }

        auto geometric_ratio = [&](double* coefficient, double* ratio) {
            const double s0 = evaluate_symbolic_at(summand, index_name, 0.0);
            const double s1 = evaluate_symbolic_at(summand, index_name, 1.0);
            const double s2 = evaluate_symbolic_at(summand, index_name, 2.0);
            const double s3 = evaluate_symbolic_at(summand, index_name, 3.0);
            if (mymath::is_near_zero(s0, 1e-10)) {
                return false;
            }
            const double candidate = s1 / s0;
            if (!mymath::is_near_zero(s2 - s1 * candidate, 1e-8) ||
                !mymath::is_near_zero(s3 - s2 * candidate, 1e-8)) {
                return false;
            }
            *coefficient = s0;
            *ratio = candidate;
            return true;
        };

        double geometric_coefficient = 0.0;
        double geometric_ratio_value = 0.0;
        if (!geometric_ratio(&geometric_coefficient, &geometric_ratio_value)) {
            throw std::runtime_error(
                "series_sum currently supports polynomial summands up to degree 3 and common geometric series");
        }

        const std::string coefficient_text =
            format_symbolic_scalar(geometric_coefficient);
        const std::string ratio_text =
            format_symbolic_scalar(geometric_ratio_value);

        if (upper_is_infinite) {
            if (mymath::abs(geometric_ratio_value) >= 1.0 - 1e-10) {
                throw std::runtime_error(
                    "series_sum infinite geometric series requires |r| < 1");
            }
            if (mymath::is_near_zero(geometric_ratio_value - 1.0, 1e-10)) {
                throw std::runtime_error(
                    "series_sum infinite geometric series diverges for r = 1");
            }
            *output = simplify_symbolic_text(
                "(" + coefficient_text + ") * (" + ratio_text + ") ^ (" +
                arguments[2] + ") / (1 - (" + ratio_text + "))");
            return true;
        }

        const std::string geometric_primitive_text =
            mymath::is_near_zero(geometric_ratio_value - 1.0, 1e-10)
                ? "(" + coefficient_text + ") * (" + index_name + " + 1)"
                : "(" + coefficient_text + ") * (1 - (" + ratio_text +
                      ") ^ (" + index_name + " + 1)) / (1 - (" +
                      ratio_text + "))";
        const SymbolicExpression primitive =
            SymbolicExpression::parse(geometric_primitive_text).simplify();
        *output = finite_sum_from_primitive(primitive);
        return true;
    }

    if (split_named_call(trimmed, "simplify", &inside)) {
        const std::string argument = trim_copy(inside);
        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(argument, false, &variable_name, &expression);
        *output = expression.simplify().to_string();
        return true;
    }

    std::string transform_inside;
    std::string transform_command;
    if (split_named_call(trimmed, "fourier", &transform_inside)) {
        transform_command = "fourier";
    } else if (split_named_call(trimmed, "ifourier", &transform_inside) ||
               split_named_call(trimmed, "inverse_fourier", &transform_inside)) {
        transform_command = "ifourier";
    } else if (split_named_call(trimmed, "laplace", &transform_inside)) {
        transform_command = "laplace";
    } else if (split_named_call(trimmed, "ilaplace", &transform_inside) ||
               split_named_call(trimmed, "inverse_laplace", &transform_inside)) {
        transform_command = "ilaplace";
    } else if (split_named_call(trimmed, "ztrans", &transform_inside) ||
               split_named_call(trimmed, "z_transform", &transform_inside)) {
        transform_command = "ztrans";
    } else if (split_named_call(trimmed, "iztrans", &transform_inside) ||
               split_named_call(trimmed, "inverse_z", &transform_inside)) {
        transform_command = "iztrans";
    }
    if (!transform_command.empty()) {
        const std::vector<std::string> arguments =
            split_top_level_arguments(transform_inside);
        if (arguments.size() != 1 && arguments.size() != 3) {
            throw std::runtime_error(
                transform_command +
                " expects either one symbolic expression or expression plus input/output variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], false, &variable_name, &expression);

        std::string input_variable;
        std::string output_variable;
        if (arguments.size() == 3) {
            input_variable = trim_copy(arguments[1]);
            output_variable = trim_copy(arguments[2]);
            if (!is_identifier_text(input_variable) || !is_identifier_text(output_variable)) {
                throw std::runtime_error(transform_command + " variable names must be identifiers");
            }
        } else if (transform_command == "fourier") {
            input_variable = variable_name.empty() ? "t" : variable_name;
            output_variable = "w";
        } else if (transform_command == "ifourier") {
            input_variable = variable_name.empty() ? "w" : variable_name;
            output_variable = "t";
        } else if (transform_command == "laplace") {
            input_variable = variable_name.empty() ? "t" : variable_name;
            output_variable = "s";
        } else if (transform_command == "ilaplace") {
            input_variable = variable_name.empty() ? "s" : variable_name;
            output_variable = "t";
        } else if (transform_command == "ztrans") {
            input_variable = variable_name.empty() ? "n" : variable_name;
            output_variable = "z";
        } else {
            input_variable = variable_name.empty() ? "z" : variable_name;
            output_variable = "n";
        }

        if (transform_command == "fourier") {
            *output = expression.fourier_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else if (transform_command == "ifourier") {
            *output = expression.inverse_fourier_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else if (transform_command == "laplace") {
            *output = expression.laplace_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else if (transform_command == "ilaplace") {
            *output = expression.inverse_laplace_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else if (transform_command == "ztrans") {
            *output = expression.z_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else {
            *output = expression.inverse_z_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        }
        return true;
    }

    if (split_named_call(trimmed, "gradient", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.empty()) {
            throw std::runtime_error(
                "gradient expects a symbolic expression and optional variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            parse_symbolic_variable_arguments(arguments,
                                              1,
                                              expression.identifier_variables());
        *output = symbolic_vector_to_string(expression.gradient(variables));
        return true;
    }

    if (split_named_call(trimmed, "hessian", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.empty()) {
            throw std::runtime_error(
                "hessian expects a symbolic expression and optional variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            parse_symbolic_variable_arguments(arguments,
                                              1,
                                              expression.identifier_variables());
        *output = symbolic_matrix_to_string(expression.hessian(variables));
        return true;
    }

    if (split_named_call(trimmed, "jacobian", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() < 2) {
            throw std::runtime_error(
                "jacobian expects a bracketed expression list and variable names");
        }

        const std::vector<SymbolicExpression> expressions =
            parse_symbolic_expression_list(arguments[0]);
        std::vector<std::string> fallback_variables;
        for (const SymbolicExpression& expression_item : expressions) {
            const std::vector<std::string> identifiers =
                expression_item.identifier_variables();
            fallback_variables.insert(fallback_variables.end(),
                                      identifiers.begin(),
                                      identifiers.end());
        }
        std::sort(fallback_variables.begin(), fallback_variables.end());
        fallback_variables.erase(std::unique(fallback_variables.begin(),
                                             fallback_variables.end()),
                                 fallback_variables.end());
        const std::vector<std::string> variables =
            parse_symbolic_variable_arguments(arguments, 1, fallback_variables);
        *output = symbolic_matrix_to_string(
            SymbolicExpression::jacobian(expressions, variables));
        return true;
    }

    if (split_named_call(trimmed, "diff", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.empty()) {
            throw std::runtime_error(
                "diff expects a symbolic expression and optional variable names");
        }
        bool symbolic_derivative = arguments.size() == 1;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            symbolic_derivative =
                symbolic_derivative || is_identifier_text(trim_copy(arguments[i]));
        }
        if (symbolic_derivative) {
            std::string variable_name;
            SymbolicExpression expression;
            resolve_symbolic_expression(arguments[0],
                                        arguments.size() == 1,
                                        &variable_name,
                                        &expression);
            SymbolicExpression differentiated = expression;
            if (arguments.size() == 1) {
                differentiated = differentiated.derivative(variable_name).simplify();
            } else {
                for (std::size_t i = 1; i < arguments.size(); ++i) {
                    const std::string derivative_variable = trim_copy(arguments[i]);
                    if (!is_identifier_text(derivative_variable)) {
                        throw std::runtime_error(
                            "diff variable arguments must be identifiers");
                    }
                    differentiated =
                        differentiated.derivative(derivative_variable).simplify();
                }
            }
            *output = differentiated.simplify().to_string();
            return true;
        }
        if (arguments.size() != 2) {
            throw std::runtime_error(
                "diff numeric evaluation expects exactly two arguments");
        }
        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        DecimalParser parser(arguments[1], &impl_->variables, &impl_->functions);
        *output = format_decimal(normalize_result(analysis.derivative(parser.parse())));
        return true;
    }

    if (split_named_call(trimmed, "limit", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 2 && arguments.size() != 3) {
            throw std::runtime_error(
                "limit expects 2 arguments for a two-sided limit or 3 with direction");
        }

        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        DecimalParser point_parser(arguments[1], &impl_->variables, &impl_->functions);
        int direction = 0;
        if (arguments.size() == 3) {
            DecimalParser direction_parser(arguments[2], &impl_->variables, &impl_->functions);
            const double direction_value = direction_parser.parse();
            if (!is_integer_double(direction_value)) {
                throw std::runtime_error("limit direction must be -1, 0, or 1");
            }
            direction = static_cast<int>(round_to_long_long(direction_value));
        }

        *output = format_decimal(normalize_result(
            analysis.limit(point_parser.parse(), direction)));
        return true;
    }

    if (split_named_call(trimmed, "integral", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 1 && arguments.size() != 2 &&
            arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "integral expects 1 argument for symbolic indefinite integral, "
                "2 arguments for indefinite value, "
                "3 for definite integral, or 4 for anchor and constant");
        }

        if (arguments.size() == 1 ||
            (arguments.size() == 2 && is_identifier_text(trim_copy(arguments[1])))) {
            std::string variable_name;
            SymbolicExpression expression;
            resolve_symbolic_expression(arguments[0],
                                        arguments.size() == 1,
                                        &variable_name,
                                        &expression);
            if (arguments.size() == 2) {
                variable_name = trim_copy(arguments[1]);
            }
            *output = expression.integral(variable_name).simplify().to_string() + " + C";
            return true;
        }
        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        if (arguments.size() == 2) {
            DecimalParser x_parser(arguments[1], &impl_->variables, &impl_->functions);
            *output = format_decimal(
                normalize_result(analysis.indefinite_integral_at(x_parser.parse())));
            return true;
        }
        if (arguments.size() == 3) {
            DecimalParser left_parser(arguments[1], &impl_->variables, &impl_->functions);
            DecimalParser right_parser(arguments[2], &impl_->variables, &impl_->functions);
            *output = format_decimal(normalize_result(
                analysis.definite_integral(left_parser.parse(), right_parser.parse())));
            return true;
        }

        DecimalParser x_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser anchor_parser(arguments[2], &impl_->variables, &impl_->functions);
        DecimalParser constant_parser(arguments[3], &impl_->variables, &impl_->functions);
        *output = format_decimal(normalize_result(
            analysis.indefinite_integral_at(x_parser.parse(),
                                            anchor_parser.parse(),
                                            constant_parser.parse())));
        return true;
    }

    std::string multivariable_inside;
    std::string multivariable_command;
    if (split_named_call(trimmed, "double_integral", &multivariable_inside)) {
        multivariable_command = "double_integral";
    } else if (split_named_call(trimmed, "double_integral_cyl", &multivariable_inside)) {
        multivariable_command = "double_integral_cyl";
    } else if (split_named_call(trimmed, "double_integral_polar", &multivariable_inside)) {
        multivariable_command = "double_integral_polar";
    } else if (split_named_call(trimmed, "triple_integral", &multivariable_inside)) {
        multivariable_command = "triple_integral";
    } else if (split_named_call(trimmed, "triple_integral_cyl", &multivariable_inside)) {
        multivariable_command = "triple_integral_cyl";
    } else if (split_named_call(trimmed, "triple_integral_sph", &multivariable_inside)) {
        multivariable_command = "triple_integral_sph";
    }
    if (!multivariable_command.empty()) {
        const std::vector<std::string> arguments =
            split_top_level_arguments(multivariable_inside);

        if (multivariable_command == "double_integral") {
            if (arguments.size() != 5 && arguments.size() != 7) {
                throw std::runtime_error(
                    "double_integral expects expr, x0, x1, y0, y1, and optional nx, ny");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            const std::vector<int> subdivisions =
                parse_subdivisions(arguments, 5, {32, 32});
            const MultivariableIntegrator integrator(
                [evaluate_expression](const std::vector<double>& point) {
                    return evaluate_expression(
                        {{"x", point[0]}, {"y", point[1]}});
                });
            *output = format_decimal(normalize_result(integrator.integrate(
                {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
                 {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])}},
                subdivisions)));
            return true;
        }

        if (multivariable_command == "double_integral_cyl" ||
            multivariable_command == "double_integral_polar") {
            if (arguments.size() != 5 && arguments.size() != 7) {
                throw std::runtime_error(
                    multivariable_command +
                    " expects expr, r0, r1, theta0, theta1, and optional nr, ntheta");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            const std::vector<int> subdivisions =
                parse_subdivisions(arguments, 5, {32, 32});
            const MultivariableIntegrator integrator(
                [evaluate_expression](const std::vector<double>& point) {
                    const double r = point[0];
                    const double theta = point[1];
                    const double x = r * mymath::cos(theta);
                    const double y = r * mymath::sin(theta);
                    return evaluate_expression(
                               {{"r", r},
                                {"theta", theta},
                                {"x", x},
                                {"y", y}}) *
                           r;
                });
            *output = format_decimal(normalize_result(integrator.integrate(
                {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
                 {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])}},
                subdivisions)));
            return true;
        }

        if (multivariable_command == "triple_integral") {
            if (arguments.size() != 7 && arguments.size() != 10) {
                throw std::runtime_error(
                    "triple_integral expects expr, x0, x1, y0, y1, z0, z1, and optional nx, ny, nz");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            const std::vector<int> subdivisions =
                parse_subdivisions(arguments, 7, {16, 16, 16});
            const MultivariableIntegrator integrator(
                [evaluate_expression](const std::vector<double>& point) {
                    return evaluate_expression(
                        {{"x", point[0]}, {"y", point[1]}, {"z", point[2]}});
                });
            *output = format_decimal(normalize_result(integrator.integrate(
                {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
                 {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])},
                 {parse_decimal_argument(arguments[5]), parse_decimal_argument(arguments[6])}},
                subdivisions)));
            return true;
        }

        if (multivariable_command == "triple_integral_cyl") {
            if (arguments.size() != 7 && arguments.size() != 10) {
                throw std::runtime_error(
                    "triple_integral_cyl expects expr, r0, r1, theta0, theta1, z0, z1, and optional nr, ntheta, nz");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            const std::vector<int> subdivisions =
                parse_subdivisions(arguments, 7, {16, 16, 16});
            const MultivariableIntegrator integrator(
                [evaluate_expression](const std::vector<double>& point) {
                    const double r = point[0];
                    const double theta = point[1];
                    const double z = point[2];
                    const double x = r * mymath::cos(theta);
                    const double y = r * mymath::sin(theta);
                    return evaluate_expression(
                               {{"r", r},
                                {"theta", theta},
                                {"z", z},
                                {"x", x},
                                {"y", y}}) *
                           r;
                });
            *output = format_decimal(normalize_result(integrator.integrate(
                {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
                 {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])},
                 {parse_decimal_argument(arguments[5]), parse_decimal_argument(arguments[6])}},
                subdivisions)));
            return true;
        }

        if (arguments.size() != 7 && arguments.size() != 10) {
            throw std::runtime_error(
                "triple_integral_sph expects expr, rho0, rho1, theta0, theta1, phi0, phi1, and optional nrho, ntheta, nphi");
        }
        const auto evaluate_expression =
            build_scoped_decimal_evaluator(arguments[0]);
        const std::vector<int> subdivisions =
            parse_subdivisions(arguments, 7, {16, 16, 16});
        const MultivariableIntegrator integrator(
            [evaluate_expression](const std::vector<double>& point) {
                const double rho = point[0];
                const double theta = point[1];
                const double phi = point[2];
                const double sin_phi = mymath::sin(phi);
                const double cos_phi = mymath::cos(phi);
                const double cos_theta = mymath::cos(theta);
                const double sin_theta = mymath::sin(theta);
                const double x = rho * sin_phi * cos_theta;
                const double y = rho * sin_phi * sin_theta;
                const double z = rho * cos_phi;
                return evaluate_expression(
                           {{"rho", rho},
                            {"theta", theta},
                            {"phi", phi},
                            {"x", x},
                            {"y", y},
                            {"z", z}}) *
                       rho * rho * sin_phi;
            });
        *output = format_decimal(normalize_result(integrator.integrate(
            {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
             {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])},
             {parse_decimal_argument(arguments[5]), parse_decimal_argument(arguments[6])}},
            subdivisions)));
        return true;
    }

    std::string ode_inside;
    std::string ode_command_name;
    if (split_named_call(trimmed, "ode", &ode_inside)) {
        ode_command_name = "ode";
    } else if (split_named_call(trimmed, "ode_table", &ode_inside)) {
        ode_command_name = "ode_table";
    }
    if (!ode_command_name.empty()) {
        const std::vector<std::string> arguments = split_top_level_arguments(ode_inside);
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                ode_command_name +
                " expects rhs, x0, y0, x1, optional steps, optional event, and optional params");
        }

        DecimalParser x0_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser y0_parser(arguments[2], &impl_->variables, &impl_->functions);
        DecimalParser x1_parser(arguments[3], &impl_->variables, &impl_->functions);
        int steps = ode_command_name == "ode" ? 100 : 10;

        std::size_t optional_index = 4;
        int parsed_steps = steps;
        if (optional_index < arguments.size() &&
            try_parse_positive_step_argument(arguments[optional_index], &parsed_steps)) {
            steps = parsed_steps;
            ++optional_index;
        }

        std::string event_expression;
        bool has_event = false;
        StoredValue parameter_value;
        bool has_parameter = false;
        if (optional_index < arguments.size()) {
            if (optional_index + 1 == arguments.size()) {
                if (is_matrix_argument(arguments[optional_index])) {
                    parameter_value = evaluate_expression_value(this,
                                                               impl_.get(),
                                                               arguments[optional_index],
                                                               false);
                    has_parameter = true;
                } else {
                    event_expression = arguments[optional_index];
                    has_event = true;
                }
                ++optional_index;
            } else {
                event_expression = arguments[optional_index];
                has_event = true;
                ++optional_index;
                parameter_value =
                    evaluate_expression_value(this, impl_.get(), arguments[optional_index], false);
                has_parameter = true;
                ++optional_index;
            }
        }

        if (optional_index != arguments.size()) {
            throw std::runtime_error(
                ode_command_name +
                " received too many optional arguments");
        }

        const auto evaluate_rhs = build_scoped_scalar_evaluator(arguments[0]);
        std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)> evaluate_event;
        if (has_event) {
            evaluate_event = build_scoped_scalar_evaluator(event_expression);
        }
        const ODESolver solver(
            [evaluate_rhs,
             has_parameter,
             parameter_value,
             &make_scalar_stored,
             &append_parameter_assignments](double x_value, double y_value) {
                std::vector<std::pair<std::string, StoredValue>> assignments;
                assignments.reserve(has_parameter ? 4 : 2);
                assignments.push_back({"x", make_scalar_stored(x_value)});
                assignments.push_back({"y", make_scalar_stored(y_value)});
                if (has_parameter) {
                    append_parameter_assignments(parameter_value, &assignments);
                }
                return evaluate_rhs(assignments);
            },
            has_event
                ? ODESolver::EventFunction(
                      [evaluate_event,
                       has_parameter,
                       parameter_value,
                       &make_scalar_stored,
                       &append_parameter_assignments](double x_value, double y_value) {
                          std::vector<std::pair<std::string, StoredValue>> assignments;
                          assignments.reserve(has_parameter ? 4 : 2);
                          assignments.push_back({"x", make_scalar_stored(x_value)});
                          assignments.push_back({"y", make_scalar_stored(y_value)});
                          if (has_parameter) {
                              append_parameter_assignments(parameter_value, &assignments);
                          }
                          return evaluate_event(assignments);
                      })
                : ODESolver::EventFunction());
        const double x0 = x0_parser.parse();
        const double y0 = y0_parser.parse();
        const double x1 = x1_parser.parse();

        if (ode_command_name == "ode") {
            *output = format_decimal(normalize_result(
                solver.solve(x0, y0, x1, steps)));
            return true;
        }

        const std::vector<ODEPoint> points =
            solver.solve_trajectory(x0, y0, x1, steps);
        matrix::Matrix table(points.size(), 2, 0.0);
        for (std::size_t i = 0; i < points.size(); ++i) {
            table.at(i, 0) = normalize_result(points[i].x);
            table.at(i, 1) = normalize_result(points[i].y);
        }
        *output = matrix_literal_expression(table);
        return true;
    }

    std::string ode_system_inside;
    std::string ode_system_command_name;
    if (split_named_call(trimmed, "ode_system", &ode_system_inside)) {
        ode_system_command_name = "ode_system";
    } else if (split_named_call(trimmed, "ode_system_table", &ode_system_inside)) {
        ode_system_command_name = "ode_system_table";
    }
    if (!ode_system_command_name.empty()) {
        const std::vector<std::string> arguments =
            split_top_level_arguments(ode_system_inside);
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                ode_system_command_name +
                " expects rhs_vector, x0, y0_vector, x1, optional steps, optional event, and optional params");
        }

        int steps = ode_system_command_name == "ode_system" ? 100 : 10;
        std::size_t optional_index = 4;
        int parsed_steps = steps;
        if (optional_index < arguments.size() &&
            try_parse_positive_step_argument(arguments[optional_index], &parsed_steps)) {
            steps = parsed_steps;
            ++optional_index;
        }

        std::string event_expression;
        bool has_event = false;
        StoredValue parameter_value;
        bool has_parameter = false;
        if (optional_index < arguments.size()) {
            if (optional_index + 1 == arguments.size()) {
                if (is_matrix_argument(arguments[optional_index])) {
                    parameter_value = evaluate_expression_value(this,
                                                               impl_.get(),
                                                               arguments[optional_index],
                                                               false);
                    has_parameter = true;
                } else {
                    event_expression = arguments[optional_index];
                    has_event = true;
                }
                ++optional_index;
            } else {
                event_expression = arguments[optional_index];
                has_event = true;
                ++optional_index;
                parameter_value =
                    evaluate_expression_value(this, impl_.get(), arguments[optional_index], false);
                has_parameter = true;
                ++optional_index;
            }
        }

        if (optional_index != arguments.size()) {
            throw std::runtime_error(
                ode_system_command_name +
                " received too many optional arguments");
        }

        const double x0 = parse_decimal_argument(arguments[1]);
        const double x1 = parse_decimal_argument(arguments[3]);
        const std::vector<double> initial_state =
            matrix_to_vector_values(parse_matrix_argument(arguments[2], ode_system_command_name),
                                    ode_system_command_name);
        const auto evaluate_rhs_matrix =
            build_scoped_matrix_evaluator(arguments[0]);
        std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)> evaluate_event;
        if (has_event) {
            evaluate_event = build_scoped_scalar_evaluator(event_expression);
        }

        const ODESystemSolver solver(
            [evaluate_rhs_matrix,
             has_parameter,
             parameter_value,
             &make_scalar_stored,
             &vector_to_column_matrix,
             &append_parameter_assignments](double x_value,
                                            const std::vector<double>& y_value) {
                std::vector<std::pair<std::string, StoredValue>> assignments;
                assignments.reserve(y_value.size() + (has_parameter ? 4 : 2));

                assignments.push_back({"x", make_scalar_stored(x_value)});

                StoredValue y_matrix_stored;
                y_matrix_stored.is_matrix = true;
                y_matrix_stored.matrix = vector_to_column_matrix(y_value);
                assignments.push_back({"y", y_matrix_stored});

                for (std::size_t i = 0; i < y_value.size(); ++i) {
                    assignments.push_back({"y" + std::to_string(i + 1),
                                           make_scalar_stored(y_value[i])});
                }
                if (has_parameter) {
                    append_parameter_assignments(parameter_value, &assignments);
                }

                const matrix::Matrix rhs_matrix = evaluate_rhs_matrix(assignments);
                if (!rhs_matrix.is_vector()) {
                    throw std::runtime_error("ODE system right-hand side must evaluate to a vector");
                }
                const std::size_t result_size =
                    rhs_matrix.rows == 1 ? rhs_matrix.cols : rhs_matrix.rows;
                if (result_size != y_value.size()) {
                    throw std::runtime_error("ODE system right-hand side dimension mismatch");
                }

                std::vector<double> result(result_size, 0.0);
                for (std::size_t i = 0; i < result_size; ++i) {
                    result[i] = rhs_matrix.rows == 1 ? rhs_matrix.at(0, i) : rhs_matrix.at(i, 0);
                }
                return result;
            },
            has_event
                ? ODESystemSolver::EventFunction(
                      [evaluate_event,
                       has_parameter,
                       parameter_value,
                       &make_scalar_stored,
                       &vector_to_column_matrix,
                       &append_parameter_assignments](double x_value,
                                                      const std::vector<double>& y_value) {
                          std::vector<std::pair<std::string, StoredValue>> assignments;
                          assignments.reserve(y_value.size() + (has_parameter ? 4 : 2));
                          assignments.push_back({"x", make_scalar_stored(x_value)});

                          StoredValue y_matrix_stored;
                          y_matrix_stored.is_matrix = true;
                          y_matrix_stored.matrix = vector_to_column_matrix(y_value);
                          assignments.push_back({"y", y_matrix_stored});

                          for (std::size_t i = 0; i < y_value.size(); ++i) {
                              assignments.push_back({"y" + std::to_string(i + 1),
                                                     make_scalar_stored(y_value[i])});
                          }
                          if (has_parameter) {
                              append_parameter_assignments(parameter_value, &assignments);
                          }
                          return evaluate_event(assignments);
                      })
                : ODESystemSolver::EventFunction());

        if (ode_system_command_name == "ode_system") {
            const std::vector<double> final_state =
                solver.solve(x0, initial_state, x1, steps);
            *output = matrix::Matrix::vector(final_state).to_string();
            return true;
        }

        const std::vector<ODESystemPoint> points =
            solver.solve_trajectory(x0, initial_state, x1, steps);
        matrix::Matrix table(points.size(), initial_state.size() + 1, 0.0);
        for (std::size_t row = 0; row < points.size(); ++row) {
            table.at(row, 0) = normalize_result(points[row].x);
            for (std::size_t col = 0; col < points[row].y.size(); ++col) {
                table.at(row, col + 1) = normalize_result(points[row].y[col]);
            }
        }
        *output = matrix_literal_expression(table);
        return true;
    }

    std::string planning_inside;
    std::string planning_command;
    if (split_named_call(trimmed, "lp_max", &planning_inside)) {
        planning_command = "lp_max";
    } else if (split_named_call(trimmed, "lp_min", &planning_inside)) {
        planning_command = "lp_min";
    } else if (split_named_call(trimmed, "ilp_max", &planning_inside)) {
        planning_command = "ilp_max";
    } else if (split_named_call(trimmed, "ilp_min", &planning_inside)) {
        planning_command = "ilp_min";
    } else if (split_named_call(trimmed, "milp_max", &planning_inside)) {
        planning_command = "milp_max";
    } else if (split_named_call(trimmed, "milp_min", &planning_inside)) {
        planning_command = "milp_min";
    } else if (split_named_call(trimmed, "bip_max", &planning_inside) ||
               split_named_call(trimmed, "binary_max", &planning_inside)) {
        planning_command = "bip_max";
    } else if (split_named_call(trimmed, "bip_min", &planning_inside) ||
               split_named_call(trimmed, "binary_min", &planning_inside)) {
        planning_command = "bip_min";
    }
    if (!planning_command.empty()) {
        const std::vector<std::string> arguments =
            split_top_level_arguments(planning_inside);
        const std::vector<double> objective = matrix_to_vector_values(
            parse_matrix_argument(arguments[0], planning_command), planning_command);
        const std::size_t variable_count = objective.size();
        const double planning_tolerance = 1e-8;
        const bool is_binary_program =
            planning_command == "bip_max" || planning_command == "bip_min";
        const bool is_mixed_integer =
            planning_command == "milp_max" || planning_command == "milp_min";
        const bool is_integer_program =
            planning_command == "ilp_max" || planning_command == "ilp_min";
        const bool maximize =
            planning_command == "lp_max" || planning_command == "ilp_max" ||
            planning_command == "milp_max" || planning_command == "bip_max";

        std::size_t argument_index = 1;
        const matrix::Matrix inequality_matrix =
            parse_matrix_argument(arguments[argument_index++], planning_command);
        const std::vector<double> inequality_rhs = matrix_to_vector_values(
            parse_matrix_argument(arguments[argument_index++], planning_command), planning_command);

        matrix::Matrix equality_matrix(0, variable_count, 0.0);
        std::vector<double> equality_rhs;
        std::vector<double> lower_bounds(variable_count, 0.0);
        std::vector<double> upper_bounds(variable_count, 0.0);
        std::vector<double> integrality(variable_count, 0.0);

        if (is_binary_program) {
            if (arguments.size() != 3 && arguments.size() != 5) {
                throw std::runtime_error(
                    planning_command +
                    " expects objective_vector, A, b, and optional Aeq, beq");
            }
            if (arguments.size() == 5) {
                equality_matrix =
                    parse_matrix_argument(arguments[argument_index++], planning_command);
                equality_rhs = matrix_to_vector_values(
                    parse_matrix_argument(arguments[argument_index++], planning_command),
                    planning_command);
            }
            std::fill(lower_bounds.begin(), lower_bounds.end(), 0.0);
            std::fill(upper_bounds.begin(), upper_bounds.end(), 1.0);
            std::fill(integrality.begin(), integrality.end(), 1.0);
        } else if (is_mixed_integer) {
            if (arguments.size() != 6 && arguments.size() != 8) {
                throw std::runtime_error(
                    planning_command +
                    " expects objective_vector, A, b, lower_bounds, upper_bounds, integrality, and optional Aeq, beq");
            }
            if (arguments.size() == 8) {
                equality_matrix =
                    parse_matrix_argument(arguments[argument_index++], planning_command);
                equality_rhs = matrix_to_vector_values(
                    parse_matrix_argument(arguments[argument_index++], planning_command),
                    planning_command);
            }
            lower_bounds = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
            upper_bounds = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
            integrality = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
        } else {
            if (arguments.size() != 5 && arguments.size() != 7) {
                throw std::runtime_error(
                    planning_command +
                    " expects objective_vector, A, b, lower_bounds, upper_bounds, and optional Aeq, beq");
            }
            if (arguments.size() == 7) {
                equality_matrix =
                    parse_matrix_argument(arguments[argument_index++], planning_command);
                equality_rhs = matrix_to_vector_values(
                    parse_matrix_argument(arguments[argument_index++], planning_command),
                    planning_command);
            }
            lower_bounds = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
            upper_bounds = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
            if (is_integer_program) {
                std::fill(integrality.begin(), integrality.end(), 1.0);
            }
        }

        if (argument_index != arguments.size()) {
            throw std::runtime_error(planning_command + " received an invalid argument count");
        }

        if (inequality_matrix.cols != variable_count ||
            inequality_rhs.size() != inequality_matrix.rows ||
            equality_matrix.cols != variable_count ||
            equality_rhs.size() != equality_matrix.rows ||
            lower_bounds.size() != variable_count ||
            upper_bounds.size() != variable_count ||
            integrality.size() != variable_count) {
            throw std::runtime_error(planning_command + " dimension mismatch");
        }

        for (std::size_t i = 0; i < variable_count; ++i) {
            if (lower_bounds[i] > upper_bounds[i]) {
                throw std::runtime_error(planning_command + " requires lower_bounds <= upper_bounds");
            }
            if (integrality[i] != 0.0 && !is_integer_double(lower_bounds[i])) {
                throw std::runtime_error(planning_command + " requires integer lower bounds for integer variables");
            }
            if (integrality[i] != 0.0 && !is_integer_double(upper_bounds[i])) {
                throw std::runtime_error(planning_command + " requires integer upper bounds for integer variables");
            }
        }

        std::vector<double> transformed_objective = objective;
        if (!maximize) {
            for (double& value : transformed_objective) {
                value = -value;
            }
        }

        std::vector<std::size_t> integer_indices;
        std::vector<std::size_t> continuous_indices;
        for (std::size_t i = 0; i < variable_count; ++i) {
            if (mymath::is_near_zero(integrality[i], planning_tolerance)) {
                continuous_indices.push_back(i);
            } else {
                integer_indices.push_back(i);
            }
        }

        if (integer_indices.empty()) {
            std::vector<double> best_solution;
            double best_value = 0.0;
            if (!solve_linear_box_problem(transformed_objective,
                                          inequality_matrix,
                                          inequality_rhs,
                                          equality_matrix,
                                          equality_rhs,
                                          lower_bounds,
                                          upper_bounds,
                                          planning_tolerance,
                                          &best_solution,
                                          &best_value)) {
                throw std::runtime_error(planning_command + " found no feasible bounded solution");
            }
            *output = format_planning_result(best_solution,
                                             dot_product(objective, best_solution));
            return true;
        }

        std::vector<long long> integer_lower(variable_count, 0);
        std::vector<long long> integer_upper(variable_count, 0);
        for (std::size_t index : integer_indices) {
            integer_lower[index] = round_to_long_long(lower_bounds[index]);
            integer_upper[index] = round_to_long_long(upper_bounds[index]);
        }

        bool found = false;
        double best_value = 0.0;
        std::vector<double> best_solution(variable_count, 0.0);
        std::vector<long long> current_integer_values(variable_count, 0);

        std::function<void(std::size_t, long double)> search_integer =
            [&](std::size_t depth, long double current_objective) {
                for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                    long double assigned_total = 0.0L;
                    for (std::size_t assigned_depth = 0; assigned_depth < depth; ++assigned_depth) {
                        const std::size_t col = integer_indices[assigned_depth];
                        assigned_total += static_cast<long double>(inequality_matrix.at(row, col)) *
                                          static_cast<long double>(current_integer_values[col]);
                    }

                    long double minimum_possible = assigned_total;
                    for (std::size_t remaining_depth = depth;
                         remaining_depth < integer_indices.size();
                         ++remaining_depth) {
                        const std::size_t col = integer_indices[remaining_depth];
                        const long double coefficient =
                            static_cast<long double>(inequality_matrix.at(row, col));
                        minimum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(integer_lower[col])
                                : coefficient * static_cast<long double>(integer_upper[col]);
                    }
                    for (std::size_t col : continuous_indices) {
                        const long double coefficient =
                            static_cast<long double>(inequality_matrix.at(row, col));
                        minimum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(lower_bounds[col])
                                : coefficient * static_cast<long double>(upper_bounds[col]);
                    }
                    if (minimum_possible >
                        static_cast<long double>(inequality_rhs[row]) + planning_tolerance) {
                        return;
                    }
                }

                for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                    long double assigned_total = 0.0L;
                    for (std::size_t assigned_depth = 0; assigned_depth < depth; ++assigned_depth) {
                        const std::size_t col = integer_indices[assigned_depth];
                        assigned_total += static_cast<long double>(equality_matrix.at(row, col)) *
                                          static_cast<long double>(current_integer_values[col]);
                    }

                    long double minimum_possible = assigned_total;
                    long double maximum_possible = assigned_total;
                    for (std::size_t remaining_depth = depth;
                         remaining_depth < integer_indices.size();
                         ++remaining_depth) {
                        const std::size_t col = integer_indices[remaining_depth];
                        const long double coefficient =
                            static_cast<long double>(equality_matrix.at(row, col));
                        minimum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(integer_lower[col])
                                : coefficient * static_cast<long double>(integer_upper[col]);
                        maximum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(integer_upper[col])
                                : coefficient * static_cast<long double>(integer_lower[col]);
                    }
                    for (std::size_t col : continuous_indices) {
                        const long double coefficient =
                            static_cast<long double>(equality_matrix.at(row, col));
                        minimum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(lower_bounds[col])
                                : coefficient * static_cast<long double>(upper_bounds[col]);
                        maximum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(upper_bounds[col])
                                : coefficient * static_cast<long double>(lower_bounds[col]);
                    }

                    const long double target = static_cast<long double>(equality_rhs[row]);
                    if (target < minimum_possible - planning_tolerance ||
                        target > maximum_possible + planning_tolerance) {
                        return;
                    }
                }

                long double optimistic_objective = current_objective;
                for (std::size_t remaining_depth = depth;
                     remaining_depth < integer_indices.size();
                     ++remaining_depth) {
                    const std::size_t col = integer_indices[remaining_depth];
                    const long double coefficient =
                        static_cast<long double>(transformed_objective[col]);
                    optimistic_objective +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(integer_upper[col])
                            : coefficient * static_cast<long double>(integer_lower[col]);
                }
                for (std::size_t col : continuous_indices) {
                    const long double coefficient =
                        static_cast<long double>(transformed_objective[col]);
                    optimistic_objective +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(upper_bounds[col])
                            : coefficient * static_cast<long double>(lower_bounds[col]);
                }
                if (found &&
                    optimistic_objective <=
                        static_cast<long double>(best_value) + planning_tolerance) {
                    return;
                }

                if (depth == integer_indices.size()) {
                    std::vector<double> candidate(variable_count, 0.0);
                    for (std::size_t col = 0; col < variable_count; ++col) {
                        candidate[col] = lower_bounds[col];
                    }
                    for (std::size_t col : integer_indices) {
                        candidate[col] = static_cast<double>(current_integer_values[col]);
                    }

                    if (continuous_indices.empty()) {
                        bool feasible = true;
                        for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                            long double total = 0.0L;
                            for (std::size_t col = 0; col < variable_count; ++col) {
                                total += static_cast<long double>(inequality_matrix.at(row, col)) *
                                         static_cast<long double>(candidate[col]);
                            }
                            if (total >
                                static_cast<long double>(inequality_rhs[row]) + planning_tolerance) {
                                feasible = false;
                                break;
                            }
                        }
                        if (!feasible) {
                            return;
                        }
                        for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                            long double total = 0.0L;
                            for (std::size_t col = 0; col < variable_count; ++col) {
                                total += static_cast<long double>(equality_matrix.at(row, col)) *
                                         static_cast<long double>(candidate[col]);
                            }
                            if (mymath::abs(static_cast<double>(total - equality_rhs[row])) >
                                planning_tolerance) {
                                feasible = false;
                                break;
                            }
                        }
                        if (!feasible) {
                            return;
                        }

                        const double objective_value = dot_product(transformed_objective, candidate);
                        if (!found || objective_value > best_value + planning_tolerance) {
                            found = true;
                            best_value = objective_value;
                            best_solution = candidate;
                        }
                        return;
                    }

                    matrix::Matrix reduced_inequality(inequality_matrix.rows,
                                                      continuous_indices.size(),
                                                      0.0);
                    std::vector<double> reduced_inequality_rhs(inequality_rhs.size(), 0.0);
                    for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                        long double rhs_adjustment = static_cast<long double>(inequality_rhs[row]);
                        for (std::size_t col : integer_indices) {
                            rhs_adjustment -=
                                static_cast<long double>(inequality_matrix.at(row, col)) *
                                static_cast<long double>(candidate[col]);
                        }
                        reduced_inequality_rhs[row] = static_cast<double>(rhs_adjustment);
                        for (std::size_t reduced_col = 0;
                             reduced_col < continuous_indices.size();
                             ++reduced_col) {
                            reduced_inequality.at(row, reduced_col) =
                                inequality_matrix.at(row, continuous_indices[reduced_col]);
                        }
                    }

                    matrix::Matrix reduced_equality(equality_matrix.rows,
                                                    continuous_indices.size(),
                                                    0.0);
                    std::vector<double> reduced_equality_rhs(equality_rhs.size(), 0.0);
                    for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                        long double rhs_adjustment = static_cast<long double>(equality_rhs[row]);
                        for (std::size_t col : integer_indices) {
                            rhs_adjustment -=
                                static_cast<long double>(equality_matrix.at(row, col)) *
                                static_cast<long double>(candidate[col]);
                        }
                        reduced_equality_rhs[row] = static_cast<double>(rhs_adjustment);
                        for (std::size_t reduced_col = 0;
                             reduced_col < continuous_indices.size();
                             ++reduced_col) {
                            reduced_equality.at(row, reduced_col) =
                                equality_matrix.at(row, continuous_indices[reduced_col]);
                        }
                    }

                    std::vector<double> reduced_objective(continuous_indices.size(), 0.0);
                    std::vector<double> reduced_lower(continuous_indices.size(), 0.0);
                    std::vector<double> reduced_upper(continuous_indices.size(), 0.0);
                    for (std::size_t reduced_col = 0;
                         reduced_col < continuous_indices.size();
                         ++reduced_col) {
                        const std::size_t original_col = continuous_indices[reduced_col];
                        reduced_objective[reduced_col] = transformed_objective[original_col];
                        reduced_lower[reduced_col] = lower_bounds[original_col];
                        reduced_upper[reduced_col] = upper_bounds[original_col];
                    }

                    std::vector<double> reduced_solution;
                    double reduced_objective_value = 0.0;
                    if (!solve_linear_box_problem(reduced_objective,
                                                  reduced_inequality,
                                                  reduced_inequality_rhs,
                                                  reduced_equality,
                                                  reduced_equality_rhs,
                                                  reduced_lower,
                                                  reduced_upper,
                                                  planning_tolerance,
                                                  &reduced_solution,
                                                  &reduced_objective_value)) {
                        return;
                    }

                    for (std::size_t reduced_col = 0;
                         reduced_col < continuous_indices.size();
                         ++reduced_col) {
                        candidate[continuous_indices[reduced_col]] = reduced_solution[reduced_col];
                    }

                    const double objective_value = dot_product(transformed_objective, candidate);
                    if (!found || objective_value > best_value + planning_tolerance) {
                        found = true;
                        best_value = objective_value;
                        best_solution = candidate;
                    }
                    return;
                }

                const std::size_t current_col = integer_indices[depth];
                const bool descending = transformed_objective[current_col] >= 0.0;
                if (descending) {
                    for (long long value = integer_upper[current_col];
                         value >= integer_lower[current_col];
                         --value) {
                        current_integer_values[current_col] = value;
                        search_integer(
                            depth + 1,
                            current_objective +
                                static_cast<long double>(transformed_objective[current_col]) *
                                    static_cast<long double>(value));
                    }
                } else {
                    for (long long value = integer_lower[current_col];
                         value <= integer_upper[current_col];
                         ++value) {
                        current_integer_values[current_col] = value;
                        search_integer(
                            depth + 1,
                            current_objective +
                                static_cast<long double>(transformed_objective[current_col]) *
                                    static_cast<long double>(value));
                    }
                }
            };

        search_integer(0, 0.0L);
        if (!found) {
            throw std::runtime_error(planning_command + " found no feasible mixed-integer solution");
        }

        *output = format_planning_result(best_solution,
                                         dot_product(objective, best_solution));
        return true;
    }

    std::string root_inside;
    std::string root_command_name;
    if (split_named_call(trimmed, "solve", &root_inside)) {
        root_command_name = "solve";
    } else if (split_named_call(trimmed, "bisect", &root_inside)) {
        root_command_name = "bisect";
    } else if (split_named_call(trimmed, "secant", &root_inside)) {
        root_command_name = "secant";
    } else if (split_named_call(trimmed, "fixed_point", &root_inside)) {
        root_command_name = "fixed_point";
    }
    if (!root_command_name.empty()) {
        const std::vector<std::string> arguments = split_top_level_arguments(root_inside);
        if (root_command_name == "solve") {
            if (arguments.size() == 2 &&
                !is_matrix_argument(arguments[0]) &&
                !is_matrix_argument(arguments[1])) {
                const auto evaluate_expression =
                    build_scoped_decimal_evaluator(arguments[0]);
                double x = parse_decimal_argument(arguments[1]);
                for (int iteration = 0; iteration < 64; ++iteration) {
                    const double fx = evaluate_expression({{"x", x}});
                    if (mymath::abs(fx) <= root_function_tolerance(fx)) {
                        *output = format_decimal(normalize_result(x));
                        return true;
                    }
                    const double h = root_derivative_step(x);
                    const long double derivative =
                        (static_cast<long double>(evaluate_expression({{"x", x + h}})) -
                         static_cast<long double>(evaluate_expression({{"x", x - h}}))) /
                        (2.0L * static_cast<long double>(h));
                    if (mymath::abs_long_double(derivative) <=
                        1e-12L * std::max(1.0L, mymath::abs_long_double(static_cast<long double>(fx)))) {
                        throw std::runtime_error("solve failed because the derivative vanished");
                    }
                    const double next = static_cast<double>(
                        static_cast<long double>(x) -
                        static_cast<long double>(fx) / derivative);
                    if (mymath::abs(next - x) <=
                        root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
                        *output = format_decimal(normalize_result(next));
                        return true;
                    }
                    x = next;
                }
                *output = format_decimal(normalize_result(x));
                return true;
            }
        } else if (root_command_name == "bisect") {
            if (arguments.size() != 3 || is_matrix_argument(arguments[0])) {
                throw std::runtime_error("bisect expects expression, a, b");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            double left = parse_decimal_argument(arguments[1]);
            double right = parse_decimal_argument(arguments[2]);
            if (left > right) {
                std::swap(left, right);
            }
            double left_value = evaluate_expression({{"x", left}});
            double right_value = evaluate_expression({{"x", right}});
            if (left_value * right_value > 0.0) {
                throw std::runtime_error("bisect requires f(a) and f(b) to have opposite signs");
            }
            for (int iteration = 0; iteration < 100; ++iteration) {
                const double mid = 0.5 * (left + right);
                const double mid_value = evaluate_expression({{"x", mid}});
                if (mymath::abs(mid_value) <= root_function_tolerance(mid_value) ||
                    mymath::abs(right - left) <=
                        root_position_tolerance(std::max(mymath::abs(left), mymath::abs(right)))) {
                    *output = format_decimal(normalize_result(mid));
                    return true;
                }
                if ((left_value < 0.0 && mid_value > 0.0) ||
                    (left_value > 0.0 && mid_value < 0.0)) {
                    right = mid;
                    right_value = mid_value;
                } else {
                    left = mid;
                    left_value = mid_value;
                }
            }
            *output = format_decimal(normalize_result(0.5 * (left + right)));
            return true;
        } else if (root_command_name == "secant") {
            if (arguments.size() != 3 || is_matrix_argument(arguments[0])) {
                throw std::runtime_error("secant expects expression, x0, x1");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            double x0 = parse_decimal_argument(arguments[1]);
            double x1 = parse_decimal_argument(arguments[2]);
            for (int iteration = 0; iteration < 64; ++iteration) {
                const double f0 = evaluate_expression({{"x", x0}});
                const double f1 = evaluate_expression({{"x", x1}});
                const long double denominator =
                    static_cast<long double>(f1) - static_cast<long double>(f0);
                if (mymath::abs_long_double(denominator) <=
                    1e-12L * std::max({1.0L,
                                       mymath::abs_long_double(static_cast<long double>(f0)),
                                       mymath::abs_long_double(static_cast<long double>(f1))})) {
                    throw std::runtime_error("secant failed because consecutive function values matched");
                }
                const double next = static_cast<double>(
                    static_cast<long double>(x1) -
                    static_cast<long double>(f1) *
                        (static_cast<long double>(x1) - static_cast<long double>(x0)) /
                        denominator);
                if (mymath::abs(next - x1) <=
                    root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x1)))) {
                    *output = format_decimal(normalize_result(next));
                    return true;
                }
                x0 = x1;
                x1 = next;
            }
            *output = format_decimal(normalize_result(x1));
            return true;
        } else if (root_command_name == "fixed_point") {
            if (arguments.size() != 2 || is_matrix_argument(arguments[0])) {
                throw std::runtime_error("fixed_point expects expression, x0");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            double x = parse_decimal_argument(arguments[1]);
            for (int iteration = 0; iteration < 128; ++iteration) {
                const double next = evaluate_expression({{"x", x}});
                if (mymath::abs(next - x) <=
                    root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
                    *output = format_decimal(normalize_result(next));
                    return true;
                }
                x = next;
            }
            *output = format_decimal(normalize_result(x));
            return true;
        }
    }

    std::string decomposition_inside;
    std::string decomposition_command;
    if (split_named_call(trimmed, "eig", &decomposition_inside)) {
        decomposition_command = "eig";
    } else if (split_named_call(trimmed, "svd", &decomposition_inside)) {
        decomposition_command = "svd";
    }
    if (!decomposition_command.empty()) {
        const std::vector<std::string> arguments = split_top_level_arguments(decomposition_inside);
        if (arguments.size() != 1 || !is_matrix_argument(arguments[0])) {
            throw std::runtime_error(decomposition_command + " expects exactly one matrix argument");
        }

        const matrix::Matrix matrix_value = evaluate_expression_value(this,
                                                                      impl_.get(),
                                                                      arguments[0],
                                                                      false).matrix;
        if (decomposition_command == "svd") {
            *output = "U: " + matrix::svd_u(matrix_value).to_string() +
                      "\nS: " + matrix::svd_s(matrix_value).to_string() +
                      "\nVt: " + matrix::svd_vt(matrix_value).to_string();
            return true;
        }

        try {
            *output = "values: " + matrix::eigenvalues(matrix_value).to_string() +
                      "\nvectors: " + matrix::eigenvectors(matrix_value).to_string();
            return true;
        } catch (const std::exception&) {
            if (matrix_value.rows == 2 && matrix_value.cols == 2) {
                const double trace = matrix_value.at(0, 0) + matrix_value.at(1, 1);
                const double det = matrix::determinant(matrix_value);
                const double discriminant = trace * trace - 4.0 * det;
                if (discriminant < 0.0) {
                    const double real = trace * 0.5;
                    const double imag = mymath::sqrt(-discriminant) * 0.5;
                    std::ostringstream out;
                    out << "values: [complex(" << format_decimal(real) << ", "
                        << format_decimal(imag) << "), complex("
                        << format_decimal(real) << ", " << format_decimal(-imag)
                        << ")]\nvectors: unavailable for complex eigenvalues";
                    *output = out.str();
                    return true;
                }
            }
            throw;
        }
    }

    if (split_named_call(trimmed, "extrema", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3) {
            throw std::runtime_error("extrema expects exactly three arguments");
        }

        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        DecimalParser left_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser right_parser(arguments[2], &impl_->variables, &impl_->functions);
        const std::vector<ExtremumPoint> points =
            analysis.solve_extrema(left_parser.parse(), right_parser.parse());

        if (points.empty()) {
            *output = "No extrema found in the given interval.";
            return true;
        }

        std::ostringstream out;
        auto display_value = [](double value) {
            if (is_integer_double(value, 1e-6)) {
                return format_decimal(static_cast<double>(round_to_long_long(value)));
            }
            return format_decimal(normalize_result(value));
        };
        for (std::size_t i = 0; i < points.size(); ++i) {
            if (i != 0) {
                out << '\n';
            }
            out << (points[i].is_maximum ? "max" : "min")
                << ": x = " << display_value(points[i].x)
                << ", f(x) = " << display_value(points[i].value);
        }
        *output = out.str();
        return true;
    }

    return false;
}

std::string Calculator::save_state(const std::string& path) const {
    for (const auto& [name, value] : impl_->variables) {
        (void)name;
        if (value.is_matrix) {
            throw std::runtime_error("save_state does not yet support matrix variables");
        }
    }

    const std::filesystem::path target_path(path);
    const std::filesystem::path temp_path =
        target_path.parent_path() /
        (target_path.filename().string() + ".tmp-save");
    std::ofstream out(temp_path);
    if (!out) {
        throw std::runtime_error("unable to open file for writing: " + path);
    }

    out << "STATE_V3\n";

    for (const auto& [name, value] : impl_->variables) {
        if (value.is_string) {
            out << "VAR\t" << encode_state_field(name)
                << "\tSTRING\t" << encode_state_field(value.string_value) << '\n';
        } else if (value.exact) {
            out << "VAR\t" << encode_state_field(name)
                << "\tEXACT\t" << value.rational.numerator
                << '\t' << value.rational.denominator
                << '\t' << std::setprecision(17) << value.decimal << '\n';
        } else {
            out << "VAR\t" << encode_state_field(name)
                << "\tDECIMAL\t" << std::setprecision(17) << value.decimal << '\n';
        }
        if (value.has_precise_decimal_text) {
            out << "PRECISE\t" << encode_state_field(name)
                << '\t' << encode_state_field(value.precise_decimal_text) << '\n';
        }
        if (value.has_symbolic_text) {
            out << "SYMBOLIC\t" << encode_state_field(name)
                << '\t' << encode_state_field(value.symbolic_text) << '\n';
        }
    }

    for (const auto& [name, function] : impl_->functions) {
        out << "EXPRFUNC\t"
            << encode_state_field(name + "(" + function.parameter_name + ") = " +
                                  function.expression)
            << '\n';
    }

    for (const auto& [name, function] : impl_->script_functions) {
        std::ostringstream source;
        source << "fn " << name << "(";
        for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
            if (i != 0) {
                source << ", ";
            }
            source << function.parameter_names[i];
        }
        source << ") " << render_script_block(*function.body, 0);
        out << "SCRIPT\t" << encode_state_field(source.str()) << '\n';
    }

    out.close();
    if (!out) {
        std::error_code remove_error;
        std::filesystem::remove(temp_path, remove_error);
        throw std::runtime_error("unable to finish writing state file: " + path);
    }

    std::error_code rename_error;
    std::filesystem::rename(temp_path, target_path, rename_error);
    if (rename_error) {
        std::error_code remove_existing_error;
        std::filesystem::remove(target_path, remove_existing_error);
        rename_error.clear();
        std::filesystem::rename(temp_path, target_path, rename_error);
    }
    if (rename_error) {
        std::error_code remove_error;
        std::filesystem::remove(temp_path, remove_error);
        throw std::runtime_error("unable to replace state file: " + path);
    }

    return "Saved variables to: " + path;
}

std::string Calculator::load_state(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("unable to open file for reading: " + path);
    }

    std::map<std::string, StoredValue> loaded;
    std::map<std::string, CustomFunction> loaded_functions;
    std::map<std::string, ScriptFunction> loaded_script_functions;
    std::string line;
    int state_version = 1;

    auto split_tab_fields = [](const std::string& row_text) {
        std::vector<std::string> parts;
        std::size_t start = 0;
        for (std::size_t i = 0; i <= row_text.size(); ++i) {
            if (i == row_text.size() || row_text[i] == '\t') {
                parts.push_back(row_text.substr(start, i - start));
                start = i + 1;
            }
        }
        return parts;
    };

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line == "STATE_V2") {
            state_version = 2;
            continue;
        }
        if (line == "STATE_V3") {
            state_version = 3;
            continue;
        }

        const std::vector<std::string> parts = split_tab_fields(line);
        if (state_version >= 2) {
            if (parts.empty()) {
                continue;
            }
            if (parts[0] == "VAR") {
                if (parts.size() < 4) {
                    throw std::runtime_error("invalid save file format");
                }

                StoredValue value;
                const std::string name = decode_state_field(parts[1]);
                if (parts[2] == "STRING") {
                    if (parts.size() != 4) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.is_string = true;
                    value.string_value = decode_state_field(parts[3]);
                } else if (parts[2] == "EXACT") {
                    if (parts.size() != 6) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.exact = true;
                    value.rational = Rational(std::stoll(parts[3]), std::stoll(parts[4]));
                    value.decimal = std::stod(parts[5]);
                } else if (parts[2] == "DECIMAL") {
                    if ((state_version == 2 && parts.size() != 4 && parts.size() != 5) ||
                        (state_version >= 3 && parts.size() != 4)) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.decimal = std::stod(parts[3]);
                    if (state_version == 2 && parts.size() == 5) {
                        value.has_precise_decimal_text = true;
                        value.precise_decimal_text = decode_state_field(parts[4]);
                    }
                } else {
                    throw std::runtime_error("invalid save file format");
                }
                loaded[name] = value;
                continue;
            }

            if (state_version >= 3 && parts[0] == "PRECISE") {
                if (parts.size() != 3) {
                    throw std::runtime_error("invalid save file format");
                }
                const std::string name = decode_state_field(parts[1]);
                auto it = loaded.find(name);
                if (it == loaded.end() || it->second.is_matrix || it->second.is_string) {
                    throw std::runtime_error("invalid save file format");
                }
                it->second.has_precise_decimal_text = true;
                it->second.precise_decimal_text = decode_state_field(parts[2]);
                continue;
            }

            if (state_version >= 3 && parts[0] == "SYMBOLIC") {
                if (parts.size() != 3) {
                    throw std::runtime_error("invalid save file format");
                }
                const std::string name = decode_state_field(parts[1]);
                auto it = loaded.find(name);
                if (it == loaded.end() || it->second.is_matrix || it->second.is_string) {
                    throw std::runtime_error("invalid save file format");
                }
                it->second.has_symbolic_text = true;
                it->second.symbolic_text = decode_state_field(parts[2]);
                continue;
            }

            if (parts[0] == "EXPRFUNC") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                std::string function_name;
                std::string parameter_name;
                std::string body;
                const std::string definition = decode_state_field(parts[1]);
                if (!split_function_definition(definition,
                                               &function_name,
                                               &parameter_name,
                                               &body)) {
                    throw std::runtime_error("invalid save file format");
                }
                loaded_functions[function_name] = {parameter_name, body};
                continue;
            }

            if (parts[0] == "SCRIPT") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                Calculator temp;
                temp.execute_script(decode_state_field(parts[1]), false);
                for (const auto& [name, function] : temp.impl_->script_functions) {
                    loaded_script_functions[name] = function;
                }
                continue;
            }

            throw std::runtime_error("invalid save file format");
        }

        if (parts.size() != 5) {
            throw std::runtime_error("invalid save file format");
        }

        StoredValue value;
        value.exact = std::stoi(parts[1]) != 0;
        value.rational = Rational(std::stoll(parts[2]), std::stoll(parts[3]));
        value.decimal = std::stod(parts[4]);
        loaded[parts[0]] = value;
    }

    impl_->variables = loaded;
    if (state_version >= 2) {
        impl_->functions = loaded_functions;
        impl_->script_functions = loaded_script_functions;
    }
    return "Loaded variables from: " + path;
}

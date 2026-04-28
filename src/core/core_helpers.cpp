#include "calculator_internal_types.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <random>
#include <sstream>
#include <utility>

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

bool is_integer_double(double x, double eps) {
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
    if (name == "inf" || name == "infinity" || name == "oo") {
        *value = mymath::infinity();
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

long long prev_prime_ll(long long value) {
    if (value <= 2) {
        throw std::runtime_error("prev_prime requires n > 2");
    }
    long long candidate = value - 1;
    if (candidate % 2 == 0) {
        --candidate;
    }
    while (candidate >= 2 && !is_prime_ll(candidate)) {
        candidate -= 2;
    }
    if (candidate < 2) {
        throw std::runtime_error("prev_prime requires n > 2");
    }
    return candidate;
}

long long euler_phi_ll(long long value) {
    if (value <= 0) {
        throw std::runtime_error("euler_phi only accepts positive integers");
    }
    long long n = value;
    long long result = value;
    for (long long p = 2; p * p <= n; ++p) {
        if (n % p != 0) {
            continue;
        }
        while (n % p == 0) {
            n /= p;
        }
        result -= result / p;
    }
    if (n > 1) {
        result -= result / n;
    }
    return result;
}

long long mobius_ll(long long value) {
    if (value <= 0) {
        throw std::runtime_error("mobius only accepts positive integers");
    }
    long long n = value;
    int prime_factor_count = 0;
    for (long long p = 2; p * p <= n; ++p) {
        if (n % p != 0) {
            continue;
        }
        n /= p;
        ++prime_factor_count;
        if (n % p == 0) {
            return 0;
        }
        while (n % p == 0) {
            n /= p;
        }
    }
    if (n > 1) {
        ++prime_factor_count;
    }
    return prime_factor_count % 2 == 0 ? 1 : -1;
}

long long prime_pi_ll(long long value) {
    if (value < 2) {
        return 0;
    }
    long long count = 0;
    for (long long n = 2; n <= value; ++n) {
        if (is_prime_ll(n)) {
            ++count;
        }
    }
    return count;
}

long long extended_gcd_ll(long long a, long long b, long long* x, long long* y) {
    long long old_r = a;
    long long r = b;
    long long old_s = 1;
    long long s = 0;
    long long old_t = 0;
    long long t = 1;

    while (r != 0) {
        const long long quotient = old_r / r;
        const long long next_r = old_r - quotient * r;
        old_r = r;
        r = next_r;

        const long long next_s = old_s - quotient * s;
        old_s = s;
        s = next_s;

        const long long next_t = old_t - quotient * t;
        old_t = t;
        t = next_t;
    }

    if (old_r < 0) {
        old_r = -old_r;
        old_s = -old_s;
        old_t = -old_t;
    }
    *x = old_s;
    *y = old_t;
    return old_r;
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

Rational::Rational(long long num, long long den)
    : numerator(num), denominator(den) {
    normalize();
}

void Rational::normalize() {
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

bool Rational::is_integer() const {
    return denominator == 1;
}

std::string Rational::to_string() const {
    if (denominator == 1) {
        return std::to_string(numerator);
    }
    return std::to_string(numerator) + "/" + std::to_string(denominator);
}

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
    return format_decimal(value, process_display_precision());
}

std::string format_decimal(double value, int precision) {
    precision = std::clamp(precision, kMinDisplayPrecision, kMaxDisplayPrecision);
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    return out.str();
}

namespace {

int& mutable_process_display_precision() {
    static int precision = kDefaultDisplayPrecision;
    return precision;
}

} // namespace

void set_process_display_precision(int precision) {
    mutable_process_display_precision() =
        std::clamp(precision, kMinDisplayPrecision, kMaxDisplayPrecision);
}

int process_display_precision() {
    return mutable_process_display_precision();
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

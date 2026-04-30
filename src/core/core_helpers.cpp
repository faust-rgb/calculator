#include "calculator_internal_types.h"
#include "../math/mymath.h"
#include "utils.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
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
        *value = mymath::kSpeedOfLight;
        return true;
    }
    if (name == "G") {
        *value = mymath::kGravitationalConstant;
        return true;
    }
    if (name == "h") {
        *value = mymath::kPlanckConstant;
        return true;
    }
    if (name == "k") {
        *value = mymath::kBoltzmannConstant;
        return true;
    }
    if (name == "NA") {
        *value = mymath::kAvogadroNumber;
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

// ============================================================================
// 有理数运算（实现在 precise/rational.cpp）
// ============================================================================

// Rational 方法实现已移至 precise/rational.cpp

std::string format_decimal(double value) {
    return format_decimal(value, process_display_precision());
}

std::string format_decimal(double value, int precision) {
    value = normalize_display_decimal(value);
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

// ============================================================================
// 精确小数运算（实现在 precise/precise_decimal.cpp）
// ============================================================================

// PreciseDecimal 实现已移至 precise/precise_decimal.cpp

// rational_to_double 已移至 precise/rational.cpp

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

bool is_identifier_text(const std::string& text) {
    if (text.empty() ||
        !std::isalpha(static_cast<unsigned char>(text.front()))) {
        return false;
    }
    for (char ch : text) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_' && ch != '\'') {
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
            *lhs = utils::trim_copy(expression.substr(0, i));
            *rhs = utils::trim_copy(expression.substr(i + 1));
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
    const std::string trimmed = utils::trim_copy(expression);
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

    *inside = utils::trim_copy(trimmed.substr(pos + 1, trimmed.size() - pos - 2));
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
            arguments.push_back(utils::trim_copy(text.substr(start, i - start)));
            start = i + 1;
        }
    }

    if (!text.empty()) {
        arguments.push_back(utils::trim_copy(text.substr(start)));
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

// rational_to_precise_decimal_text 和 stored_value_precise_decimal_text
// 已移至 precise/precise_decimal.cpp

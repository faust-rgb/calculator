#include "calculator_internal_types.h"

#include "conversion.h"
#include "matrix.h"
#include "mymath.h"
#include "symbolic_expression.h"

#include <algorithm>
#include <cctype>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

class ExactParserImpl {
public:
    ExactParserImpl(std::string source,
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

Rational parse_exact_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function) {
    ExactParserImpl parser(expression,
                           variables,
                           functions,
                           std::move(has_script_function));
    return parser.parse();
}

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
            out << format_decimal(normalize_display_decimal(numeric::to_double(value.at(row, col))));
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

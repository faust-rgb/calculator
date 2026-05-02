#include "exact_parser.h"
#include "calculator_internal_types.h"
#include "base_parser.h"
#include "mymath.h"
#include <algorithm>
#include <cctype>

class ExactParserImpl : public BaseParser {
public:
    ExactParserImpl(std::string_view source,
                    const VariableResolver& variables,
                    const std::map<std::string, CustomFunction>* functions,
                    HasScriptFunctionCallback has_script_function = {})
        : BaseParser(source),
          variables_(variables),
          functions_(functions),
          has_script_function_(std::move(has_script_function)) {}

    Rational parse() {
        Rational value = parse_comparison();
        skip_spaces();
        if (!is_at_end()) {
            throw_error_at_pos<SyntaxError>("unexpected token near: " + std::string(source_.substr(pos_, 1)));
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
                throw_error_at_pos<ExactModeUnsupported>("exact rational mode does not support non-integer exponents");
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
            const std::string name = std::string(parse_identifier());
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
                    parse_prefixed_integer_token(std::string(source_.substr(start, pos_ - start))), 1);
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
            throw_error_at_pos<std::runtime_error>("expected number");
        }

        return parse_rational_literal(std::string(source_.substr(start, pos_ - start)));
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
        if (functions_->find(name) != functions_->end()) {
            throw_error_at_pos<ExactModeUnsupported>("custom function " + name +
                                       " is not supported exactly");
        }
        if (has_script_function_ && has_script_function_(name)) {
            throw_error_at_pos<ExactModeUnsupported>("script function " + name +
                                       " is not supported exactly");
        }
        if (name == "pow") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("pow expects exactly two arguments");
            }
            if (!arguments[1].is_integer()) {
                throw_error_at_pos<ExactModeUnsupported>("exact rational mode does not support non-integer exponents");
            }
            return pow_rational(arguments[0], arguments[1].numerator);
        }
        if (name == "abs") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("abs expects exactly one argument");
            }
            return abs_rational(arguments[0]);
        }
        if (name == "step" || name == "u" || name == "heaviside") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("step expects exactly one argument");
            }
            return Rational(arguments[0].numerator >= 0 ? 1 : 0, 1);
        }
        if (name == "delta" || name == "impulse") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("delta expects exactly one argument");
            }
            return Rational(arguments[0].numerator == 0 ? 1 : 0, 1);
        }
        if (name == "not") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("not expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw_error_at_pos<std::runtime_error>("not only accepts integers");
            }
            return Rational(~arguments[0].numerator, 1);
        }
        if (name == "sign") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("sign expects exactly one argument");
            }
            if (arguments[0].numerator == 0) {
                return Rational(0, 1);
            }
            return Rational(arguments[0].numerator > 0 ? 1 : -1, 1);
        }
        if (name == "gcd") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("gcd expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("gcd only accepts integers");
            }
            return Rational(gcd_ll(arguments[0].numerator, arguments[1].numerator), 1);
        }
        if (name == "lcm") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("lcm expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("lcm only accepts integers");
            }
            return Rational(lcm_ll(arguments[0].numerator, arguments[1].numerator), 1);
        }
        if (name == "mod") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("mod expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("mod only accepts integers");
            }
            if (arguments[1].numerator == 0) {
                throw_error_at_pos<std::runtime_error>("mod divisor cannot be zero");
            }
            return Rational(arguments[0].numerator % arguments[1].numerator, 1);
        }
        if (name == "rol") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("rol expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("rol only accepts integers");
            }
            const unsigned count = normalize_rotation_count(arguments[1].numerator);
            return Rational(
                from_unsigned_bits(rotate_left_bits(
                    to_unsigned_bits(arguments[0].numerator), count)),
                1);
        }
        if (name == "ror") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("ror expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("ror only accepts integers");
            }
            const unsigned count = normalize_rotation_count(arguments[1].numerator);
            return Rational(
                from_unsigned_bits(rotate_right_bits(
                    to_unsigned_bits(arguments[0].numerator), count)),
                1);
        }
        if (name == "floor") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("floor expects exactly one argument");
            }
            return Rational(floor_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "ceil") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("ceil expects exactly one argument");
            }
            return Rational(ceil_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "round") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("round expects exactly one argument");
            }
            return Rational(round_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "trunc") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("trunc expects exactly one argument");
            }
            return Rational(trunc_to_long_long(rational_to_double(arguments[0])), 1);
        }
        if (name == "min") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("min expects exactly two arguments");
            }
            return rational_to_double(arguments[0]) < rational_to_double(arguments[1])
                       ? arguments[0]
                       : arguments[1];
        }
        if (name == "max") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("max expects exactly two arguments");
            }
            return rational_to_double(arguments[0]) > rational_to_double(arguments[1])
                       ? arguments[0]
                       : arguments[1];
        }
        if (name == "clamp") {
            if (arguments.size() != 3) {
                throw_error_at_pos<std::runtime_error>("clamp expects exactly three arguments");
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
                throw_error_at_pos<std::runtime_error>("sum expects at least one argument");
            }
            Rational total(0, 1);
            for (const Rational& value : arguments) {
                total = total + value;
            }
            return total;
        }
        if (name == "avg") {
            if (arguments.empty()) {
                throw_error_at_pos<std::runtime_error>("avg expects at least one argument");
            }
            Rational total(0, 1);
            for (const Rational& value : arguments) {
                total = total + value;
            }
            return total / Rational(static_cast<long long>(arguments.size()), 1);
        }
        if (name == "mean") {
            if (arguments.empty()) {
                throw_error_at_pos<std::runtime_error>("mean expects at least one argument");
            }
            Rational total(0, 1);
            for (const Rational& value : arguments) {
                total = total + value;
            }
            return total / Rational(static_cast<long long>(arguments.size()), 1);
        }
        if (name == "median") {
            if (arguments.empty()) {
                throw_error_at_pos<std::runtime_error>("median expects at least one argument");
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
                throw_error_at_pos<std::runtime_error>("factorial expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw_error_at_pos<std::runtime_error>("factorial only accepts integers");
            }
            return factorial_rational(arguments[0].numerator);
        }
        if (name == "nCr") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("nCr expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("nCr only accepts integers");
            }
            return combination_rational(arguments[0].numerator, arguments[1].numerator);
        }
        if (name == "binom") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("binom expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("binom only accepts integers");
            }
            return combination_rational(arguments[0].numerator, arguments[1].numerator);
        }
        if (name == "nPr") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("nPr expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("nPr only accepts integers");
            }
            return permutation_rational(arguments[0].numerator, arguments[1].numerator);
        }
        if (name == "popcount") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("popcount expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw_error_at_pos<std::runtime_error>("popcount only accepts integers");
            }
            return Rational(popcount_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "bitlen") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("bitlen expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw_error_at_pos<std::runtime_error>("bitlen only accepts integers");
            }
            return Rational(bit_length_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "ctz") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("ctz expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw_error_at_pos<std::runtime_error>("ctz only accepts integers");
            }
            return Rational(trailing_zero_count_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "clz") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("clz expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw_error_at_pos<std::runtime_error>("clz only accepts integers");
            }
            return Rational(leading_zero_count_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "parity") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("parity expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw_error_at_pos<std::runtime_error>("parity only accepts integers");
            }
            return Rational(parity_bits(to_unsigned_bits(arguments[0].numerator)), 1);
        }
        if (name == "reverse_bits") {
            if (arguments.size() != 1) {
                throw_error_at_pos<std::runtime_error>("reverse_bits expects exactly one argument");
            }
            if (!arguments[0].is_integer()) {
                throw_error_at_pos<std::runtime_error>("reverse_bits only accepts integers");
            }
            return Rational(
                from_unsigned_bits(reverse_bits(to_unsigned_bits(arguments[0].numerator))),
                1);
        }
        if (name == "and") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("and expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("and only accepts integers");
            }
            return Rational(arguments[0].numerator & arguments[1].numerator, 1);
        }
        if (name == "or") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("or expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("or only accepts integers");
            }
            return Rational(arguments[0].numerator | arguments[1].numerator, 1);
        }
        if (name == "xor") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("xor expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("xor only accepts integers");
            }
            return Rational(arguments[0].numerator ^ arguments[1].numerator, 1);
        }
        if (name == "shl") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("shl expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("shl only accepts integers");
            }
            if (arguments[1].numerator < 0) {
                throw_error_at_pos<std::runtime_error>("shift count cannot be negative");
            }
            return Rational(arguments[0].numerator << arguments[1].numerator, 1);
        }
        if (name == "shr") {
            if (arguments.size() != 2) {
                throw_error_at_pos<std::runtime_error>("shr expects exactly two arguments");
            }
            if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
                throw_error_at_pos<std::runtime_error>("shr only accepts integers");
            }
            if (arguments[1].numerator < 0) {
                throw_error_at_pos<std::runtime_error>("shift count cannot be negative");
            }
            return Rational(arguments[0].numerator >> arguments[1].numerator, 1);
        }

        throw_error_at_pos<ExactModeUnsupported>("function " + name + " is not supported exactly");
    }

    Rational lookup_variable(const std::string& name) const {
        const StoredValue* found = variables_.lookup(name);
        if (!found) {
            double constant_value = 0.0;
            if (lookup_builtin_constant(name, &constant_value)) {
                throw_error_at_pos<ExactModeUnsupported>("built-in constants are not rational");
            }
            throw_error_at_pos<std::runtime_error>("unknown variable: " + name);
        }
        if (found->is_matrix || found->is_complex) {
            throw_error_at_pos<ExactModeUnsupported>("matrix or complex variable " + name + " cannot be used exactly");
        }
        if (found->is_string) {
            throw_error_at_pos<ExactModeUnsupported>("string variable " + name + " cannot be used exactly");
        }
        if (!found->exact) {
            throw_error_at_pos<ExactModeUnsupported>("variable " + name + " is only stored approximately");
        }
        return found->rational;
    }

    VariableResolver variables_;
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
};

ExactParser::ExactParser(std::string source,
                         const VariableResolver& variables,
                         const std::map<std::string, CustomFunction>* functions,
                         HasScriptFunctionCallback has_script_function)
    : source_(std::move(source)),
      variables_(variables),
      functions_(functions),
      has_script_function_(std::move(has_script_function)) {}

Rational ExactParser::parse() {
    ExactParserImpl parser(source_, variables_, functions_, std::move(has_script_function_));
    return parser.parse();
}

Rational parse_exact_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function) {
    ExactParserImpl parser(expression,
                           variables,
                           functions,
                           std::move(has_script_function));
    return parser.parse();
}

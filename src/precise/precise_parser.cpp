// ============================================================================
// 精确小数表达式解析器
// ============================================================================

#include "types/precise_decimal.h"
#include "types/stored_value.h"
#include "core/base_parser.h"
#include "core/calculator_exceptions.h"
#include "math/mymath.h"

#include <algorithm>
#include <map>

namespace {

/**
 * @brief 使用牛顿法计算 sqrt(x)
 */
PreciseDecimal sqrt_precise_decimal_taylor(const PreciseDecimal& x, int iterations = 50) {
    if (x.is_zero()) return PreciseDecimal();
    if (x.negative) throw PreciseDecimalUnsupported("sqrt of negative number");

    double x_val = x.to_double();
    double y = x_val > 1.0 ? x_val : 1.0;
    for (int i = 0; i < iterations; ++i) {
        double next = (y + x_val / y) * 0.5;
        if (mymath::abs(next - y) < 1e-15 * mymath::abs(y)) break;
        y = next;
    }

    return PreciseDecimal::from_decimal_literal(std::to_string(y));
}

/**
 * @brief 使用泰勒级数计算 exp(x)
 */
PreciseDecimal exp_precise_decimal_taylor(const PreciseDecimal& x, int terms = 30) {
    double x_val = x.to_double();

    int n = 1;
    while (mymath::abs(x_val / n) > 1.0) ++n;

    double x_reduced = x_val / n;
    double result = 1.0;
    double term = 1.0;

    for (int i = 1; i <= terms; ++i) {
        term *= x_reduced / i;
        result += term;
        if (mymath::abs(term) < 1e-16 * mymath::abs(result)) break;
    }

    for (int i = 0; i < n; ++i) {
        result *= result;
    }

    return PreciseDecimal::from_decimal_literal(std::to_string(result));
}

/**
 * @brief 使用泰勒级数计算 ln(x)
 */
PreciseDecimal ln_precise_decimal_taylor(const PreciseDecimal& x, int terms = 50) {
    if (x.is_zero()) throw PreciseDecimalUnsupported("ln(0) is undefined");
    if (x.negative) throw PreciseDecimalUnsupported("ln of negative number");

    double x_val = x.to_double();

    double y = (x_val - 1.0) / (x_val + 1.0);
    double y_sq = y * y;
    double result = 0.0;
    double term = y;

    for (int i = 0; i < terms; ++i) {
        result += term / (2 * i + 1);
        term *= y_sq;
        if (mymath::abs(term) < 1e-16) break;
    }

    return PreciseDecimal::from_decimal_literal(std::to_string(2.0 * result));
}

/**
 * @brief 使用泰勒级数计算 sin(x)
 */
PreciseDecimal sin_precise_decimal_taylor(const PreciseDecimal& x, int terms = 30) {
    double x_val = x.to_double();

    const double pi = 3.14159265358979323846;
    while (x_val > pi) x_val -= 2 * pi;
    while (x_val < -pi) x_val += 2 * pi;

    double result = 0.0;
    double term = x_val;
    double x_sq = x_val * x_val;

    for (int i = 1; i <= terms; ++i) {
        result += term;
        term *= -x_sq / ((2 * i) * (2 * i + 1));
        if (mymath::abs(term) < 1e-16 * mymath::abs(result)) break;
    }

    return PreciseDecimal::from_decimal_literal(std::to_string(result));
}

/**
 * @brief 使用泰勒级数计算 cos(x)
 */
PreciseDecimal cos_precise_decimal_taylor(const PreciseDecimal& x, int terms = 30) {
    double x_val = x.to_double();

    const double pi = 3.14159265358979323846;
    while (x_val > pi) x_val -= 2 * pi;
    while (x_val < -pi) x_val += 2 * pi;

    double result = 0.0;
    double term = 1.0;
    double x_sq = x_val * x_val;

    for (int i = 0; i < terms; ++i) {
        result += term;
        term *= -x_sq / ((2 * i + 1) * (2 * i + 2));
        if (mymath::abs(term) < 1e-16 * mymath::abs(result)) break;
    }

    return PreciseDecimal::from_decimal_literal(std::to_string(result));
}

} // namespace

// ============================================================================
// 解析器实现
// ============================================================================

class PreciseDecimalParserImpl : public BaseParser {
public:
    PreciseDecimalParserImpl(std::string source,
                             const std::map<std::string, StoredValue>* variables)
        : BaseParser(std::move(source)),
          variables_(variables) {}

    PreciseDecimal parse() {
        PreciseDecimal value = parse_comparison();
        skip_spaces();
        if (!is_at_end()) {
            throw PreciseDecimalUnsupported("unsupported token");
        }
        return value;
    }

private:
    PreciseDecimal parse_comparison() {
        PreciseDecimal left = parse_expression();

        skip_spaces();
        bool result_bool = false;
        bool has_comparison = false;

        if (match_string("<=")) {
            PreciseDecimal right = parse_expression();
            result_bool = compare_precise_decimal(left, right) <= 0;
            has_comparison = true;
        } else if (match_string(">=")) {
            PreciseDecimal right = parse_expression();
            result_bool = compare_precise_decimal(left, right) >= 0;
            has_comparison = true;
        } else if (match_string("==")) {
            PreciseDecimal right = parse_expression();
            result_bool = compare_precise_decimal(left, right) == 0;
            has_comparison = true;
        } else if (match_string("!=")) {
            PreciseDecimal right = parse_expression();
            result_bool = compare_precise_decimal(left, right) != 0;
            has_comparison = true;
        } else if (match('<')) {
            PreciseDecimal right = parse_expression();
            result_bool = compare_precise_decimal(left, right) < 0;
            has_comparison = true;
        } else if (match('>')) {
            PreciseDecimal right = parse_expression();
            result_bool = compare_precise_decimal(left, right) > 0;
            has_comparison = true;
        }

        if (has_comparison) {
            return PreciseDecimal::from_integer_string(result_bool ? "1" : "0", false);
        }
        return left;
    }

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
        return parse_power();
    }

    PreciseDecimal parse_power() {
        PreciseDecimal value = parse_primary();
        skip_spaces();
        if (match('^')) {
            PreciseDecimal exponent = parse_unary();
            double exp_val = exponent.to_double();
            if (!mymath::is_integer(exp_val)) {
                throw PreciseDecimalUnsupported("only integer exponents are supported in precise mode");
            }
            return pow_precise_decimal(value, static_cast<long long>(exp_val));
        }
        return value;
    }

    PreciseDecimal parse_primary() {
        skip_spaces();
        if (match('(')) {
            PreciseDecimal value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        if (peek_is_alpha()) {
            const std::string name = parse_identifier();
            skip_spaces();
            if (peek() == '(') {
                expect('(');
                std::vector<PreciseDecimal> args;
                if (peek() != ')') {
                    args.push_back(parse_expression());
                    skip_spaces();
                    while (match(',')) {
                        args.push_back(parse_expression());
                        skip_spaces();
                    }
                }
                expect(')');
                return call_function(name, args);
            }
            return lookup_variable(name);
        }

        return parse_number();
    }

    PreciseDecimal call_function(const std::string& name, const std::vector<PreciseDecimal>& args) {
        if (name == "abs") {
            if (args.size() != 1) throw ArgumentError("abs expects 1 argument");
            PreciseDecimal res = args[0];
            res.negative = false;
            return res;
        }
        if (name == "min") {
            if (args.size() < 1) throw ArgumentError("min expects at least 1 argument");
            PreciseDecimal res = args[0];
            for (std::size_t i = 1; i < args.size(); ++i) {
                if (compare_precise_decimal(args[i], res) < 0) res = args[i];
            }
            return res;
        }
        if (name == "max") {
            if (args.size() < 1) throw ArgumentError("max expects at least 1 argument");
            PreciseDecimal res = args[0];
            for (std::size_t i = 1; i < args.size(); ++i) {
                if (compare_precise_decimal(args[i], res) > 0) res = args[i];
            }
            return res;
        }
        if (name == "sqrt") {
            if (args.size() != 1) throw ArgumentError("sqrt expects 1 argument");
            return sqrt_precise_decimal_taylor(args[0]);
        }
        if (name == "exp") {
            if (args.size() != 1) throw ArgumentError("exp expects 1 argument");
            return exp_precise_decimal_taylor(args[0]);
        }
        if (name == "ln" || name == "log") {
            if (args.size() != 1) throw ArgumentError("ln expects 1 argument");
            return ln_precise_decimal_taylor(args[0]);
        }
        if (name == "sin") {
            if (args.size() != 1) throw ArgumentError("sin expects 1 argument");
            return sin_precise_decimal_taylor(args[0]);
        }
        if (name == "cos") {
            if (args.size() != 1) throw ArgumentError("cos expects 1 argument");
            return cos_precise_decimal_taylor(args[0]);
        }
        if (name == "tan") {
            if (args.size() != 1) throw ArgumentError("tan expects 1 argument");
            PreciseDecimal c = cos_precise_decimal_taylor(args[0]);
            if (c.is_zero()) throw PreciseDecimalUnsupported("tan undefined at this point");
            PreciseDecimal s = sin_precise_decimal_taylor(args[0]);
            return divide_precise_decimal(s, c);
        }
        if (name == "floor") {
            if (args.size() != 1) throw ArgumentError("floor expects 1 argument");
            double val = args[0].to_double();
            val = mymath::floor(val);
            return PreciseDecimal::from_decimal_literal(std::to_string(val));
        }
        if (name == "ceil") {
            if (args.size() != 1) throw ArgumentError("ceil expects 1 argument");
            double val = args[0].to_double();
            val = mymath::ceil(val);
            return PreciseDecimal::from_decimal_literal(std::to_string(val));
        }
        if (name == "round") {
            if (args.size() != 1) throw ArgumentError("round expects 1 argument");
            double val = args[0].to_double();
            val = mymath::round(val);
            return PreciseDecimal::from_decimal_literal(std::to_string(val));
        }
        if (name == "sign" || name == "sgn") {
            if (args.size() != 1) throw ArgumentError("sign expects 1 argument");
            if (args[0].is_zero()) return PreciseDecimal::from_integer_string("0", false);
            return PreciseDecimal::from_integer_string("1", args[0].negative);
        }
        throw PreciseDecimalUnsupported("function '" + name + "' is not supported in precise mode");
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
            throw SyntaxError("expected number");
        }

        return PreciseDecimal::from_decimal_literal(source_.substr(start, pos_ - start));
    }

    PreciseDecimal lookup_variable(const std::string& name) const {
        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            throw UndefinedError("unknown variable: " + name);
        }
        if (it->second.is_matrix || it->second.is_complex ||
            it->second.is_string || it->second.has_symbolic_text) {
            throw MathError("unsupported variable type for precise parsing: " + name);
        }
        return PreciseDecimal::from_decimal_literal(
            stored_value_precise_decimal_text(it->second));
    }

    const std::map<std::string, StoredValue>* variables_;
};

// ============================================================================
// 公共接口
// ============================================================================

PreciseDecimal parse_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables) {
    PreciseDecimalParserImpl parser(expression, variables);
    return parser.parse();
}

// ============================================================================
// 精确小数表达式解析器
// ============================================================================

#include "precise_decimal.h"
#include "types/stored_value.h"
#include "parser/base_parser.h"
#include "core/calculator_exceptions.h"
#include "math/helpers/base_conversions.h"
#include "math/mymath.h"

#include <algorithm>
#include <map>
#include "symbolic/core/symbolic_expression_internal.h"

namespace {

class PreciseDecimalParserImpl : public BaseParser {
public:
    PreciseDecimalParserImpl(std::string_view source,
                            const std::map<std::string, StoredValue>* variables)
        : BaseParser(source), variables_(variables) {}

    PreciseDecimal parse() {
        skip_spaces();
        PreciseDecimal value = parse_expression();
        skip_spaces();
        if (!is_at_end()) {
            throw SyntaxError("unexpected characters at end of expression");
        }
        return value;
    }

private:
    PreciseDecimal parse_expression() {
        PreciseDecimal value = parse_term();
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

    PreciseDecimal parse_term() {
        PreciseDecimal value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value *= parse_unary();
            } else if (match('/')) {
                value /= parse_unary();
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
            return -parse_unary();
        }
        return parse_power();
    }

    PreciseDecimal parse_power() {
        PreciseDecimal value = parse_primary();
        skip_spaces();
        if (match('^')) {
            PreciseDecimal exponent = parse_unary();
            if (!exponent.negative && exponent.scale == 0) {
                const std::string exponent_text = exponent.to_string();
                constexpr std::size_t kMaxSafeLongLongDigits = 18;
                if (exponent_text.size() <= kMaxSafeLongLongDigits) {
                    return precise::pow(value, std::stoll(exponent_text));
                }
            }
            if (value <= PreciseDecimal(0LL)) {
                throw PreciseDecimalUnsupported(
                    "non-integer exponents for non-positive bases are not supported in precise mode");
            }
            return precise::exp(exponent * precise::ln(value));
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
            const std::string name = std::string(parse_identifier());
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
            if (name == "pi") return precise::pi();
            if (name == "e") return precise::e();
            return lookup_variable(name);
        }

        return parse_number();
    }

    PreciseDecimal call_function(const std::string& name, const std::vector<PreciseDecimal>& args) {
        if (name == "abs") {
            if (args.size() != 1) throw ArgumentError("abs expects 1 argument");
            return precise::abs(args[0]);
        }
        if (name == "min") {
            if (args.size() < 1) throw ArgumentError("min expects at least 1 argument");
            PreciseDecimal res = args[0];
            for (std::size_t i = 1; i < args.size(); ++i) {
                if (args[i] < res) res = args[i];
            }
            return res;
        }
        if (name == "max") {
            if (args.size() < 1) throw ArgumentError("max expects at least 1 argument");
            PreciseDecimal res = args[0];
            for (std::size_t i = 1; i < args.size(); ++i) {
                if (args[i] > res) res = args[i];
            }
            return res;
        }
        if (name == "sqrt") {
            if (args.size() != 1) throw ArgumentError("sqrt expects 1 argument");
            return precise::sqrt(args[0]);
        }
        if (name == "exp") {
            if (args.size() != 1) throw ArgumentError("exp expects 1 argument");
            return precise::exp(args[0]);
        }
        if (name == "ln" || name == "log") {
            if (args.size() != 1) throw ArgumentError("ln expects 1 argument");
            return precise::ln(args[0]);
        }
        if (name == "log10") {
            if (args.size() != 1) throw ArgumentError("log10 expects 1 argument");
            return precise::log10(args[0]);
        }
        if (name == "sin") {
            if (args.size() != 1) throw ArgumentError("sin expects 1 argument");
            return precise::sin(args[0]);
        }
        if (name == "cos") {
            if (args.size() != 1) throw ArgumentError("cos expects 1 argument");
            return precise::cos(args[0]);
        }
        if (name == "tan") {
            if (args.size() != 1) throw ArgumentError("tan expects 1 argument");
            return precise::tan(args[0]);
        }
        if (name == "atan") {
            if (args.size() != 1) throw ArgumentError("atan expects 1 argument");
            return precise::atan(args[0]);
        }
        if (name == "asin") {
            if (args.size() != 1) throw ArgumentError("asin expects 1 argument");
            return precise::asin(args[0]);
        }
        if (name == "acos") {
            if (args.size() != 1) throw ArgumentError("acos expects 1 argument");
            return precise::acos(args[0]);
        }
        if (name == "sinh") {
            if (args.size() != 1) throw ArgumentError("sinh expects 1 argument");
            return precise::sinh(args[0]);
        }
        if (name == "cosh") {
            if (args.size() != 1) throw ArgumentError("cosh expects 1 argument");
            return precise::cosh(args[0]);
        }
        if (name == "tanh") {
            if (args.size() != 1) throw ArgumentError("tanh expects 1 argument");
            return precise::tanh(args[0]);
        }
        if (name == "asinh") {
            if (args.size() != 1) throw ArgumentError("asinh expects 1 argument");
            return precise::asinh(args[0]);
        }
        if (name == "acosh") {
            if (args.size() != 1) throw ArgumentError("acosh expects 1 argument");
            return precise::acosh(args[0]);
        }
        if (name == "atanh") {
            if (args.size() != 1) throw ArgumentError("atanh expects 1 argument");
            return precise::atanh(args[0]);
        }
        if (name == "floor") {
            if (args.size() != 1) throw ArgumentError("floor expects 1 argument");
            return precise::floor(args[0]);
        }
        if (name == "ceil") {
            if (args.size() != 1) throw ArgumentError("ceil expects 1 argument");
            return precise::ceil(args[0]);
        }
        if (name == "round") {
            if (args.size() != 1) throw ArgumentError("round expects 1 argument");
            return precise::round(args[0]);
        }
        if (name == "sign" || name == "sgn") {
            if (args.size() != 1) throw ArgumentError("sign expects 1 argument");
            if (mymath::is_near_zero(args[0].to_double(), 1e-10)) return PreciseDecimal(0LL);
            return PreciseDecimal(args[0].negative ? -1LL : 1LL);
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
                    std::string(source_.substr(start, pos_ - start)));
                return PreciseDecimal(integer_value);
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

        return PreciseDecimal::from_decimal_literal(std::string(source_.substr(start, pos_ - start)));
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
        return PreciseDecimal(stored_value_precise_decimal_text(it->second));
    }

    const std::map<std::string, StoredValue>* variables_;
};

} // namespace

// ============================================================================
// 公共接口
// ============================================================================

PreciseDecimal parse_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables) {
    PreciseDecimalParserImpl parser(expression, variables);
    return parser.parse();
}

// ============================================================================
// 精确小数表达式解析器
// ============================================================================

#include "precise_decimal.h"
#include "core/value/stored_value.h"
#include "parser/infra/base_parser.h"
#include "core/common/calculator_exceptions.h"
#include "parser/ast/expression_ast.h"
#include "math/functions/conversion/base_conversions.h"
#include "math/mymath.h"

#include <algorithm>
#include <map>

namespace {

PreciseDecimal evaluate_precise_ast(const ExpressionAST* ast,
                                    const std::map<std::string, StoredValue>* variables);

PreciseDecimal parse_precise_literal(const std::string& text) {
    int base = 10;
    if (text.size() >= 2 && text[0] == '0') {
        if (text[1] == 'x' || text[1] == 'X') base = 16;
        else if (text[1] == 'b' || text[1] == 'B') base = 2;
        else if (text[1] == 'o' || text[1] == 'O') base = 8;
    }
    if (base == 10) return PreciseDecimal(text);
    PreciseDecimal value(0LL);
    for (std::size_t i = 2; i < text.size(); ++i) {
        const char ch = text[i];
        int digit = ch >= '0' && ch <= '9' ? ch - '0' :
                    ch >= 'a' && ch <= 'f' ? ch - 'a' + 10 : ch - 'A' + 10;
        if (digit < 0 || digit >= base) throw SyntaxError("invalid digit in base literal");
        value = value * PreciseDecimal(static_cast<long long>(base)) + PreciseDecimal(static_cast<long long>(digit));
    }
    if (text.size() == 2) throw SyntaxError("expected digits after base prefix");
    return value;
}

PreciseDecimal evaluate_precise_call(const ExpressionAST* ast,
                                     const std::map<std::string, StoredValue>* variables) {
    std::vector<PreciseDecimal> args;
    for (const auto& child : ast->children) args.push_back(evaluate_precise_ast(child.get(), variables));
    const auto one = [&]() {
        if (args.size() != 1) throw ArgumentError(ast->identifier + " expects 1 argument");
        return args[0];
    };
    if (ast->identifier == "abs") return precise::abs(one());
    if (ast->identifier == "sqrt") return precise::sqrt(one());
    if (ast->identifier == "exp") return precise::exp(one());
    if (ast->identifier == "ln" || ast->identifier == "log") return precise::ln(one());
    if (ast->identifier == "log10") return precise::log10(one());
    if (ast->identifier == "log2") return precise::log2(one());
    if (ast->identifier == "pow") {
        if (args.size() != 2) throw ArgumentError("pow expects 2 arguments");
        return precise::pow(args[0], args[1]);
    }
    if (ast->identifier == "sin") return precise::sin(one());
    if (ast->identifier == "cos") return precise::cos(one());
    if (ast->identifier == "tan") return precise::tan(one());
    if (ast->identifier == "atan") return precise::atan(one());
    if (ast->identifier == "asin") return precise::asin(one());
    if (ast->identifier == "acos") return precise::acos(one());
    if (ast->identifier == "sinh") return precise::sinh(one());
    if (ast->identifier == "cosh") return precise::cosh(one());
    if (ast->identifier == "tanh") return precise::tanh(one());
    if (ast->identifier == "asinh") return precise::asinh(one());
    if (ast->identifier == "acosh") return precise::acosh(one());
    if (ast->identifier == "atanh") return precise::atanh(one());
    if (ast->identifier == "floor") return precise::floor(one());
    if (ast->identifier == "ceil") return precise::ceil(one());
    if (ast->identifier == "round") return precise::round(one());
    if (ast->identifier == "sign" || ast->identifier == "sgn") {
        if (args.size() != 1) throw ArgumentError(ast->identifier + " expects 1 argument");
        return args[0].is_zero() ? PreciseDecimal(0LL) : PreciseDecimal(args[0].negative ? -1LL : 1LL);
    }
    if (ast->identifier == "min" || ast->identifier == "max") {
        if (args.empty()) throw ArgumentError(ast->identifier + " expects at least 1 argument");
        PreciseDecimal result = args[0];
        for (std::size_t i = 1; i < args.size(); ++i) {
            if ((ast->identifier == "min" && args[i] < result) ||
                (ast->identifier == "max" && args[i] > result)) result = args[i];
        }
        return result;
    }
    throw PreciseDecimalUnsupported("function '" + ast->identifier + "' is not supported in precise mode");
}

PreciseDecimal evaluate_precise_ast(const ExpressionAST* ast,
                                    const std::map<std::string, StoredValue>* variables) {
    if (!ast) throw SyntaxError("null expression AST");
    switch (ast->kind) {
        case ExprKind::kNumber: return parse_precise_literal(ast->string_value);
        case ExprKind::kVariable: {
            if (ast->identifier == "pi") return precise::pi();
            if (ast->identifier == "e") return precise::e();
            const auto it = variables->find(ast->identifier);
            if (it == variables->end()) throw UndefinedError("unknown variable: " + ast->identifier);
            if (it->second.is_matrix || it->second.is_complex || it->second.is_string || it->second.has_symbolic_text)
                throw PreciseDecimalUnsupported("unsupported variable type for precise parsing: " + ast->identifier);
            if (it->second.precise_decimal_value) return *it->second.precise_decimal_value;
            return PreciseDecimal(stored_value_precise_decimal_text(it->second));
        }
        case ExprKind::kImaginary:
            throw PreciseDecimalUnsupported("complex values are not supported in precise mode");
        case ExprKind::kUnaryOp: {
            const auto value = evaluate_precise_ast(ast->children.at(0).get(), variables);
            return ast->op_char == '-' ? -value : value;
        }
        case ExprKind::kBinaryOp: {
            const auto left = evaluate_precise_ast(ast->children.at(0).get(), variables);
            const auto right = evaluate_precise_ast(ast->children.at(1).get(), variables);
            if (ast->op_char == '+') return left + right;
            if (ast->op_char == '-') return left - right;
            if (ast->op_char == '*') return left * right;
            if (ast->op_char == '/') return left / right;
            if (ast->op_char == '%') throw PreciseDecimalUnsupported("modulo is not supported in precise mode");
            if (ast->op_char == '^') {
                if (!right.negative && right.scale == 0 && right.to_string().size() <= 18)
                    return precise::pow(left, std::stoll(right.to_string()));
                if (left <= PreciseDecimal(0LL)) throw PreciseDecimalUnsupported("non-integer exponent for non-positive base");
                return precise::exp(right * precise::ln(left));
            }
            throw PreciseDecimalUnsupported("unsupported precise operator");
        }
        case ExprKind::kPostfixOp: {
            const auto value = evaluate_precise_ast(ast->children.at(0).get(), variables);
            if (ast->postfix_op == "%") return value / PreciseDecimal(100LL);
            if (ast->postfix_op != "!" && ast->postfix_op != "!!") throw PreciseDecimalUnsupported("unsupported postfix operator");
            if (value.negative || value.scale != 0) throw PreciseDecimalUnsupported("factorial requires a non-negative integer");
            PreciseDecimal result(1LL);
            const long long n = std::stoll(value.to_string());
            const long long step = ast->postfix_op == "!!" ? 2 : 1;
            const long long start = ast->postfix_op == "!!" && n % 2 == 0 ? 2 : 1;
            for (long long i = n; i >= start; i -= step) result *= PreciseDecimal(i);
            return result;
        }
        case ExprKind::kComparison: {
            const auto left = evaluate_precise_ast(ast->children.at(0).get(), variables);
            const auto right = evaluate_precise_ast(ast->children.at(1).get(), variables);
            bool result = false;
            if (ast->comparison_op == "==") result = left == right;
            else if (ast->comparison_op == "!=") result = left != right;
            else if (ast->comparison_op == "<") result = left < right;
            else if (ast->comparison_op == ">") result = left > right;
            else if (ast->comparison_op == "<=") result = left <= right;
            else if (ast->comparison_op == ">=") result = left >= right;
            return PreciseDecimal(result ? 1LL : 0LL);
        }
        case ExprKind::kLogicalOp: {
            const auto left = evaluate_precise_ast(ast->children.at(0).get(), variables);
            if (ast->comparison_op == "&&" && left.is_zero()) return PreciseDecimal(0LL);
            if (ast->comparison_op == "||" && !left.is_zero()) return PreciseDecimal(1LL);
            const auto right = evaluate_precise_ast(ast->children.at(1).get(), variables);
            return PreciseDecimal((!right.is_zero()) == (ast->comparison_op == "||") ? 1LL : 0LL);
        }
        case ExprKind::kConditional:
            return !evaluate_precise_ast(ast->children.at(0).get(), variables).is_zero()
                ? evaluate_precise_ast(ast->children.at(1).get(), variables)
                : evaluate_precise_ast(ast->children.at(2).get(), variables);
        case ExprKind::kFunctionCall: return evaluate_precise_call(ast, variables);
        default: throw PreciseDecimalUnsupported("unsupported expression in precise mode");
    }
}

#if 0
// Retained only as migration documentation. Precise evaluation now consumes
// the canonical ExpressionAST below and no longer maintains a second grammar.
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
            if (args[0].is_zero()) return PreciseDecimal(0LL);
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
                pos_ += 2;
                PreciseDecimal val(0LL);
                PreciseDecimal base_dec(static_cast<long long>(base));
                bool has_any_digit = false;
                while (!is_at_end()) {
                    const int digit = digit_value(source_[pos_]);
                    if (digit < 0 || digit >= base) {
                        break;
                    }
                    val = val * base_dec + PreciseDecimal(static_cast<long long>(digit));
                    has_any_digit = true;
                    ++pos_;
                }
                if (!has_any_digit) {
                    throw SyntaxError("expected digits after base prefix");
                }
                return val;
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
        if (it->second.precise_decimal_value) {
            return *(it->second.precise_decimal_value);
        }
        return PreciseDecimal(stored_value_precise_decimal_text(it->second));
    }

    const std::map<std::string, StoredValue>* variables_;
};
#endif

} // namespace

// ============================================================================
// 公共接口
// ============================================================================

PreciseDecimal parse_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables) {
    const auto ast = compile_expression_ast_diagnostic(expression);
    if (!ast) throw SyntaxError("failed to compile precise expression");
    return evaluate_precise_ast(ast.get(), variables);
}

// ============================================================================
// 精确小数表达式解析器
// ============================================================================

#include "precise_decimal.h"
#include "core/value/stored_value.h"
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
            if (it->second.is_matrix() || it->second.is_complex() || it->second.is_string() || it->second.has_symbolic_text)
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
                std::string right_str = right.to_string();
                if (right.scale == 0 && right_str.size() <= (right.negative ? 19 : 18)) {
                    try {
                        long long exp_val = std::stoll(right_str);
                        return precise::pow(left, exp_val);
                    } catch (const std::exception&) {}
                }
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

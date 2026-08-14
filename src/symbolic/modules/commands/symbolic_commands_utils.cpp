/**
 * @file symbolic_commands_utils.cpp
 * @brief 符号命令工具函数实现
 *
 * 本文件实现了符号命令共用的工具函数：
 * - symbolic_vector_to_string: 符号向量格式化输出
 * - symbolic_matrix_to_string: 符号矩阵格式化输出
 * - 表达式列表解析
 *
 * 这些工具函数被各符号命令处理函数调用。
 */

#include "symbolic/modules/commands/symbolic_commands_internal.h"
#include "parser/grammars/unified_expression_parser.h"
#include "core/services/string_utils.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/calculus/risch/risch_algorithm.h"
#include "polynomial/polynomial.h"
#include "parser/grammars/command_parser.h"
#include <sstream>

namespace symbolic_commands {

namespace {
using namespace symbolic_expression_internal;
}

std::string symbolic_vector_to_string(const std::vector<SymbolicExpression>& values) {
    std::ostringstream out;
    out << "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) out << ", ";
        out << values[i].simplify().to_string();
    }
    out << "]";
    return out.str();
}

std::string symbolic_matrix_to_string(const std::vector<std::vector<SymbolicExpression>>& values) {
    std::ostringstream out;
    out << "[";
    for (std::size_t row = 0; row < values.size(); ++row) {
        if (row != 0) out << ", ";
        out << "[";
        for (std::size_t col = 0; col < values[row].size(); ++col) {
            if (col != 0) out << ", ";
            out << values[row][col].simplify().to_string();
        }
        out << "]";
    }
    out << "]";
    return out.str();
}

std::vector<std::string> split_top_level_semicolon_list(const std::string& text) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    int paren_depth = 0;
    int bracket_depth = 0;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char ch = text[i];
        if (ch == '(') ++paren_depth;
        else if (ch == ')') { if (paren_depth > 0) --paren_depth; }
        else if (ch == '[') ++bracket_depth;
        else if (ch == ']') { if (bracket_depth > 0) --bracket_depth; }
        else if (ch == ';' && paren_depth == 0 && bracket_depth == 0) {
            parts.push_back(trim_copy(text.substr(start, i - start)));
            start = i + 1;
        }
    }
    parts.push_back(trim_copy(text.substr(start)));
    return parts;
}

bool parse_symbolic_vector_literal(const SymbolicCommandContext& ctx,
                                   const std::string& text,
                                   std::vector<SymbolicExpression>* components) {
    const std::string trimmed = trim_copy(text);
    if (trimmed.size() < 2 || trimmed.front() != '[' || trimmed.back() != ']') return false;
    const std::string inner = trimmed.substr(1, trimmed.size() - 2);
    const std::vector<std::string> parts =
        inner.find(';') == std::string::npos ? split_top_level_arguments(inner) : split_top_level_semicolon_list(inner);
    components->clear();
    for (const std::string& part : parts) {
        if (trim_copy(part).empty()) continue;
        std::string var; SymbolicExpression expr;
        ctx.resolve_symbolic(part, false, &var, &expr);
        components->push_back(expr);
    }
    return true;
}

bool is_infinity_literal(const std::string& text) {
    std::string value = trim_copy(text);
    if (!value.empty() && (value.front() == '+' || value.front() == '-')) value = trim_copy(value.substr(1));
    return value == "inf" || value == "infinity" || value == "oo";
}

void resolve_symbolic_expression(const SymbolicResolverContext& ctx,
                                 const std::string& argument,
                                 bool require_single_variable,
                                 std::string* variable_name,
                                 SymbolicExpression* expression) {
    const std::string trimmed_argument = trim_copy(argument);
    CommandASTNode ast = parse_command(trimmed_argument);
    
    if (ast.kind == CommandKind::kFunctionCall) {
        const auto* call = ast.as_function_call();
        if (call->name == "diff") {
            SymbolicExpression nested;
            resolve_symbolic_expression(ctx, std::string(call->arguments[0].text), call->arguments.size() == 1, variable_name, &nested);
            if (call->arguments.size() == 1) *expression = nested.derivative(*variable_name).simplify();
            else {
                SymbolicExpression diffed = nested;
                for (size_t i = 1; i < call->arguments.size(); ++i) {
                    diffed = diffed.derivative(trim_copy(call->arguments[i].text)).simplify();
                }
                *variable_name = trim_copy(call->arguments[1].text);
                *expression = diffed;
            }
            return;
        }
        if (call->name == "integral") {
            SymbolicExpression nested;
            resolve_symbolic_expression(ctx, std::string(call->arguments[0].text), call->arguments.size() == 1, variable_name, &nested);
            if (call->arguments.size() == 2) *variable_name = trim_copy(call->arguments[1].text);
            auto res = RischAlgorithm::integrate_full(nested, *variable_name);
            if (!res.success || res.type != IntegralType::kElementary) throw std::runtime_error("Integration failed in nested integral");
            *expression = res.value.simplify();
            return;
        }
        if (call->name == "poly_add" || call->name == "poly_sub" || call->name == "poly_mul" || call->name == "poly_div") {
            polynomial_ops::PolynomialContext p_ctx;
            p_ctx.resolve_symbolic = [&](const std::string& n, std::string* v) {
                SymbolicExpression r; resolve_symbolic_expression(ctx, n, true, v, &r); return r;
            };
            const auto p = polynomial_ops::build_polynomial(p_ctx, trimmed_argument);
            *variable_name = p.variable_name;
            *expression = SymbolicExpression::parse(polynomial_to_string(p.coefficients, *variable_name));
            return;
        }
    }
    if (ctx.has_custom_function(trimmed_argument)) {
        *expression = ctx.resolve_custom_function(trimmed_argument, variable_name);
        return;
    }
    *expression = SymbolicExpression::parse(ctx.expand_inline(trimmed_argument));
    const auto ids = expression->identifier_variables();
    if (ids.size() == 1) *variable_name = ids[0];
    else if (ids.empty()) *variable_name = "x";
    else if (require_single_variable) throw std::runtime_error("Multiple variables in expression");
    else *variable_name = ids.front();
}

std::vector<std::string> parse_symbolic_variable_arguments(
    const std::vector<std::string>& arguments,
    std::size_t start_index,
    const std::vector<std::string>& fallback_variables) {
    std::vector<std::string> variables;
    for (std::size_t i = start_index; i < arguments.size(); ++i) {
        const std::string variable = trim_copy(arguments[i]);
        if (!is_identifier_text(variable)) throw std::runtime_error("Variable arguments must be identifiers");
        variables.push_back(variable);
    }
    if (variables.empty()) variables = fallback_variables;
    if (variables.empty()) variables.push_back("x");
    return variables;
}

std::vector<SymbolicExpression> parse_symbolic_expression_list(
    const std::string& argument,
    const std::function<std::string(const std::string&)>& expand_inline) {
    auto expand = expand_inline ? expand_inline : [](const std::string& s) { return s; };
    std::string text = trim_copy(argument);
    if (text.size() < 2 || text.front() != '[' || text.back() != ']') return {SymbolicExpression::parse(expand(text))};
    text = trim_copy(text.substr(1, text.size() - 2));
    std::vector<std::string> texts;
    int p_d = 0, b_d = 0; std::size_t s = 0;
    for (std::size_t i = 0; i < text.size(); ++i) {
        char ch = text[i];
        if (ch == '(') ++p_d; else if (ch == ')') --p_d;
        else if (ch == '[') ++b_d; else if (ch == ']') --b_d;
        else if ((ch == ';' || ch == ',') && p_d == 0 && b_d == 0) {
            texts.push_back(trim_copy(text.substr(s, i - s))); s = i + 1;
        }
    }
    if (!text.empty()) texts.push_back(trim_copy(text.substr(s)));
    if (texts.empty()) throw std::runtime_error("Empty expression list");
    std::vector<SymbolicExpression> exprs;
    for (const auto& t : texts) {
        if (t.empty()) throw std::runtime_error("Empty item in list");
        exprs.push_back(SymbolicExpression::parse(expand_inline(t)));
    }
    return exprs;
}

} // namespace symbolic_commands

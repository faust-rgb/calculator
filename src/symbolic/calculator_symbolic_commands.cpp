// ============================================================================
// 符号命令实现
// ============================================================================

#include "calculator_symbolic_commands.h"

#include "calculator_internal_types.h"
#include "polynomial.h"

#include <algorithm>
#include <sstream>

namespace symbolic_commands {

namespace {

std::string symbolic_vector_to_string(const std::vector<SymbolicExpression>& values) {
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
}

std::string symbolic_matrix_to_string(
    const std::vector<std::vector<SymbolicExpression>>& values) {
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
}

}  // namespace

bool is_symbolic_command(const std::string& command) {
    return command == "simplify" ||
           command == "gradient" ||
           command == "hessian" ||
           command == "jacobian" ||
           command == "diff" ||
           command == "integral";
}

void resolve_symbolic_expression(const SymbolicResolverContext& ctx,
                                 const std::string& argument,
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
        resolve_symbolic_expression(ctx,
                                    nested_arguments[0],
                                    nested_arguments.size() == 1,
                                    variable_name,
                                    &nested_expression);
        if (nested_arguments.size() == 1) {
            *expression = nested_expression.derivative(*variable_name).simplify();
        } else {
            SymbolicExpression differentiated = nested_expression;
            for (std::size_t i = 1; i < nested_arguments.size(); ++i) {
                const std::string derivative_variable = trim_copy(nested_arguments[i]);
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
        resolve_symbolic_expression(ctx,
                                    nested_arguments[0],
                                    nested_arguments.size() == 1,
                                    variable_name,
                                    &nested_expression);
        if (nested_arguments.size() == 2) {
            const std::string integral_variable = trim_copy(nested_arguments[1]);
            if (!is_identifier_text(integral_variable)) {
                throw std::runtime_error(
                    "nested symbolic integral variable must be an identifier");
            }
            *variable_name = integral_variable;
        }
        *expression = nested_expression.integral(*variable_name).simplify();
        return;
    }
    if (split_named_call(trimmed_argument, "poly_add", &nested_inside) ||
        split_named_call(trimmed_argument, "poly_sub", &nested_inside) ||
        split_named_call(trimmed_argument, "poly_mul", &nested_inside) ||
        split_named_call(trimmed_argument, "poly_div", &nested_inside)) {
        polynomial_ops::PolynomialContext polynomial_ctx;
        polynomial_ctx.resolve_symbolic =
            [&](const std::string& name, std::string* variable) {
                SymbolicExpression resolved;
                resolve_symbolic_expression(ctx, name, true, variable, &resolved);
                return resolved;
            };
        const polynomial_ops::PolynomialData polynomial =
            polynomial_ops::build_polynomial(polynomial_ctx, trimmed_argument);
        *variable_name = polynomial.variable_name;
        *expression = SymbolicExpression::parse(
            polynomial_to_string(polynomial.coefficients, *variable_name));
        return;
    }
    if (ctx.has_custom_function(trimmed_argument)) {
        *expression = ctx.resolve_custom_function(trimmed_argument, variable_name);
        return;
    }

    *expression = SymbolicExpression::parse(ctx.expand_inline(trimmed_argument));
    const std::vector<std::string> identifiers = expression->identifier_variables();
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
}

std::vector<std::string> parse_symbolic_variable_arguments(
    const std::vector<std::string>& arguments,
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
}

std::vector<SymbolicExpression> parse_symbolic_expression_list(
    const std::string& argument,
    const std::function<std::string(const std::string&)>& expand_inline) {
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
        expressions.push_back(SymbolicExpression::parse(expand_inline(expression_text)));
    }
    return expressions;
}

bool handle_symbolic_command(const SymbolicCommandContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "simplify") {
        const std::string argument = trim_copy(inside);
        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(argument, false, &variable_name, &expression);
        *output = expression.simplify().to_string();
        return true;
    }

    if (command == "gradient" || command == "hessian") {
        if (arguments.empty()) {
            throw std::runtime_error(
                command + " expects a symbolic expression and optional variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            ctx.parse_symbolic_variable_arguments(arguments,
                                                  1,
                                                  expression.identifier_variables());
        if (command == "gradient") {
            *output = symbolic_vector_to_string(expression.gradient(variables));
        } else {
            *output = symbolic_matrix_to_string(expression.hessian(variables));
        }
        return true;
    }

    if (command == "jacobian") {
        if (arguments.size() < 2) {
            throw std::runtime_error(
                "jacobian expects a bracketed expression list and variable names");
        }

        const std::vector<SymbolicExpression> expressions =
            ctx.parse_symbolic_expression_list(arguments[0]);
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
            ctx.parse_symbolic_variable_arguments(arguments, 1, fallback_variables);
        *output = symbolic_matrix_to_string(
            SymbolicExpression::jacobian(expressions, variables));
        return true;
    }

    if (command == "diff") {
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
            ctx.resolve_symbolic(arguments[0],
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
        const FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
        *output = format_decimal(ctx.normalize_result(
            analysis.derivative(ctx.parse_decimal(arguments[1]))));
        return true;
    }

    if (command == "integral") {
        if (arguments.size() != 1 && arguments.size() != 2 &&
            arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "integral expects 1 argument for symbolic indefinite integral, "
                "identifier arguments for chained symbolic integrals, "
                "2 arguments for indefinite value, "
                "3 for definite integral, or 4 for anchor and constant");
        }

        bool symbolic_integral = true;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            symbolic_integral =
                symbolic_integral && is_identifier_text(trim_copy(arguments[i]));
        }
        if (symbolic_integral) {
            std::string variable_name;
            SymbolicExpression expression;
            ctx.resolve_symbolic(arguments[0],
                                 arguments.size() == 1,
                                 &variable_name,
                                 &expression);
            if (arguments.size() == 2) {
                variable_name = trim_copy(arguments[1]);
            }
            SymbolicExpression integrated = expression;
            if (arguments.size() == 1) {
                integrated = integrated.integral(variable_name).simplify();
            } else {
                for (std::size_t i = 1; i < arguments.size(); ++i) {
                    integrated =
                        integrated.integral(trim_copy(arguments[i])).simplify();
                }
            }
            *output = integrated.simplify().to_string() + " + C";
            return true;
        }

        const FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
        if (arguments.size() == 2) {
            *output = format_decimal(ctx.normalize_result(
                analysis.indefinite_integral_at(ctx.parse_decimal(arguments[1]))));
            return true;
        }
        if (arguments.size() == 3) {
            *output = format_decimal(ctx.normalize_result(
                analysis.definite_integral(ctx.parse_decimal(arguments[1]),
                                           ctx.parse_decimal(arguments[2]))));
            return true;
        }

        *output = format_decimal(ctx.normalize_result(
            analysis.indefinite_integral_at(ctx.parse_decimal(arguments[1]),
                                            ctx.parse_decimal(arguments[2]),
                                            ctx.parse_decimal(arguments[3]))));
        return true;
    }

    return false;
}

}  // namespace symbolic_commands

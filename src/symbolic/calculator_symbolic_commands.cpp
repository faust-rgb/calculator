// ============================================================================
// 符号命令实现
// ============================================================================
//
// 本文件实现符号计算的各种命令处理器，包括：
//
// 1. 基本符号操作
//    - simplify: 表达式简化
//    - expand: 表达式展开
//
// 2. 微积分命令
//    - diff: 符号求导（支持高阶偏导）
//    - integral: 符号积分（不定积分）
//
// 3. 向量场分析
//    - gradient: 梯度计算
//    - hessian: Hessian 矩阵计算
//    - jacobian: Jacobian 矩阵计算
//    - divergence/div: 散度计算
//    - curl: 旋度计算
//    - laplacian: 拉普拉斯算子
//
// 4. 微分方程求解
//    - dsolve: 常微分方程符号解（支持线性 ODE）
//
// 命令通过 SymbolicCommandContext 接收解析上下文，
// 支持自定义函数展开和多变量表达式处理。
// ============================================================================

#include "symbolic/calculator_symbolic_commands.h"
#include "parser/unified_expression_parser.h"

#include "core/calculator_internal_types.h"
#include "analysis/calculator_integration.h"
#include "analysis/multivariable_integrator.h"
#include "analysis/multidim_integration.h"
#include "polynomial/polynomial.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/integration_engine.h"
#include "parser/command_parser.h"
#include "math/mymath.h"

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

bool is_infinity_literal(const std::string& text) {
    std::string value = trim_copy(text);
    if (!value.empty() && value.front() == '+') {
        value = trim_copy(value.substr(1));
    } else if (!value.empty() && value.front() == '-') {
        value = trim_copy(value.substr(1));
    }
    return value == "inf" || value == "infinity" || value == "oo";
}

}  // namespace

bool is_symbolic_command(const std::string& command) {
    return command == "simplify" ||
           command == "expand" ||
           command == "cse" ||
           command == "gradient" ||
           command == "numerical_gradient" ||
           command == "num_grad" ||
           command == "hessian" ||
           command == "jacobian" ||
           command == "divergence" ||
           command == "div" ||
           command == "curl" ||
           command == "curl_2d" ||
           command == "laplacian" ||
           command == "implicit_diff" ||
           command == "param_deriv" ||
           command == "directional" ||
           command == "line_integral" ||
           command == "line_integral_vector" ||
           command == "surface_integral" ||
           command == "greens_theorem" ||
           command == "stokes_theorem" ||
           command == "divergence_theorem" ||
           command == "integrate_region" ||
           command == "diff" ||
           command == "integral" ||
           command == "dsolve";
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
            if (call->arguments.empty()) {
                throw std::runtime_error(
                    "nested symbolic diff expects at least one argument");
            }
            SymbolicExpression nested_expression;
            resolve_symbolic_expression(ctx,
                                        std::string(call->arguments[0].text),
                                        call->arguments.size() == 1,
                                        variable_name,
                                        &nested_expression);
            if (call->arguments.size() == 1) {
                *expression = nested_expression.derivative(*variable_name).simplify();
            } else {
                SymbolicExpression differentiated = nested_expression;
                for (std::size_t i = 1; i < call->arguments.size(); ++i) {
                    const std::string derivative_variable = trim_copy(call->arguments[i].text);
                    if (!is_identifier_text(derivative_variable)) {
                        throw std::runtime_error(
                            "nested symbolic diff variable arguments must be identifiers");
                    }
                    differentiated =
                        differentiated.derivative(derivative_variable).simplify();
                }
                *variable_name = trim_copy(call->arguments[1].text);
                *expression = differentiated;
            }
            return;
        }
        if (call->name == "integral") {
            if (call->arguments.size() != 1 && call->arguments.size() != 2) {
                throw std::runtime_error(
                    "nested symbolic integral expects expression and optional variable");
            }
            SymbolicExpression nested_expression;
            resolve_symbolic_expression(ctx,
                                        std::string(call->arguments[0].text),
                                        call->arguments.size() == 1,
                                        variable_name,
                                        &nested_expression);
            if (call->arguments.size() == 2) {
                const std::string integral_variable = trim_copy(call->arguments[1].text);
                if (!is_identifier_text(integral_variable)) {
                    throw std::runtime_error(
                        "nested symbolic integral variable must be an identifier");
                }
                *variable_name = integral_variable;
            }
            *expression = nested_expression.integral(*variable_name).simplify();
            return;
        }
        if (call->name == "poly_add" || call->name == "poly_sub" ||
            call->name == "poly_mul" || call->name == "poly_div") {
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
        // If not a list, return as a single expression
        return {SymbolicExpression::parse(expand_inline(text))};
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

    if (command == "dsolve") {
        if (arguments.size() < 1 || arguments.size() > 3) {
            throw std::runtime_error("dsolve expects (rhs, [x, y])");
        }

        std::string x_var = "x";
        std::string y_var = "y";
        if (arguments.size() >= 2) x_var = trim_copy(arguments[1]);
        if (arguments.size() >= 3) y_var = trim_copy(arguments[2]);

        SymbolicExpression rhs = SymbolicExpression::parse(arguments[0]);

        auto contains = [](const SymbolicExpression& expr, const std::string& var) {
            auto vars = expr.identifier_variables();
            return std::find(vars.begin(), vars.end(), var) != vars.end();
        };

        // 1. Check if independent of y: y' = f(x)
        if (!contains(rhs, y_var)) {
            *output = y_var + " = " + rhs.integral(x_var).simplify().to_string() + " + C";
            return true;
        }

        // 2. Check if linear: y' = -P(x)y + Q(x) => y' + P(x)y = Q(x)
        SymbolicExpression neg_P = rhs.derivative(y_var).simplify();
        if (!contains(neg_P, y_var)) {
            SymbolicExpression Q = rhs.substitute(y_var, SymbolicExpression::number(0.0)).simplify();
            if (!contains(Q, y_var)) {
                // Linear ODE: y' + P(x)y = Q(x) where P = -neg_P
                SymbolicExpression P = (-neg_P).simplify();
                SymbolicExpression mu_exponent = P.integral(x_var).simplify();

                std::string mu_str = "exp(" + mu_exponent.to_string() + ")";
                std::string mu_inv_str = "exp(" + ((-mu_exponent).simplify()).to_string() + ")";

                SymbolicExpression mu = SymbolicExpression::parse(mu_str).simplify();
                SymbolicExpression mu_inv = SymbolicExpression::parse(mu_inv_str).simplify();

                SymbolicExpression integral_part = (mu * Q).simplify().integral(x_var).simplify();
                *output = y_var + " = " + (mu_inv * (integral_part)).simplify().to_string() + " + " + mu_inv.to_string() + " * C";
                return true;
            }
        }

        throw std::runtime_error("Only linear ODEs and y' = f(x) are currently supported by dsolve.");
    }
    if (command == "simplify" || command == "expand" || command == "cse") {
        const std::string argument = trim_copy(inside);
        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(argument, false, &variable_name, &expression);
        if (command == "expand") {
            *output = expression.expand().simplify().to_string();
        } else if (command == "cse") {
            auto subs = expression.common_subexpressions();
            if (subs.empty()) {
                *output = expression.simplify().to_string();
            } else {
                std::ostringstream out;
                SymbolicExpression current = expression;
                int t_index = 0;
                for (const auto& sub : subs) {
                    std::string t_name = "_t" + std::to_string(t_index++);
                    out << t_name << " = " << sub.first.simplify().to_string() << "\n";
                    current = current.substitute_expression(sub.first, SymbolicExpression::parse(t_name)).simplify();
                }
                out << "ans = " << current.simplify().to_string();
                *output = out.str();
            }
        } else {
            *output = expression.simplify().to_string();
        }
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

    if (command == "numerical_gradient" || command == "num_grad") {
        // numerical_gradient(expr, [vars], [point])
        // Example: numerical_gradient(x^2 + y^2, x, y, 3, 4) -> [6, 8]
        // Or: numerical_gradient(x^2, x, 2) -> [4]
        if (arguments.size() < 3) {
            throw std::runtime_error(
                "numerical_gradient expects (expression, variables..., point_values...)");
        }

        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &expression);

        // Determine how many variables we have
        // Variables are identifiers, point values are numbers
        std::vector<std::string> variables;
        std::vector<double> point_values;

        for (std::size_t i = 1; i < arguments.size(); ++i) {
            const std::string arg = trim_copy(arguments[i]);
            // Try to parse as number first
            try {
                double val = ctx.parse_decimal(arg);
                point_values.push_back(val);
            } catch (...) {
                // Not a number, must be a variable name
                if (!is_identifier_text(arg)) {
                    throw std::runtime_error(
                        "numerical_gradient: argument '" + arg + "' is neither a number nor a valid variable name");
                }
                variables.push_back(arg);
            }
        }

        // If no variables specified, use expression's identifiers
        if (variables.empty()) {
            variables = expression.identifier_variables();
        }
        if (variables.empty()) {
            variables.push_back("x");
        }

        if (variables.size() != point_values.size()) {
            throw std::runtime_error(
                "numerical_gradient: number of variables (" + std::to_string(variables.size()) +
                ") must match number of point values (" + std::to_string(point_values.size()) + ")");
        }

        // Compute gradient using forward-mode AutoDiff
        std::vector<double> gradient_values;
        for (std::size_t i = 0; i < variables.size(); ++i) {
            symbolic_expression_internal::DualEvaluationContext dual_ctx;
            dual_ctx.differentiation_variable = variables[i];
            dual_ctx.point_value = point_values[i];
            for (std::size_t j = 0; j < variables.size(); ++j) {
                if (j != i) {
                    dual_ctx.other_values[variables[j]] = point_values[j];
                }
            }

            mymath::dual<double> result;
            if (!symbolic_expression_internal::try_evaluate_dual_node(expression.node_, dual_ctx, &result)) {
                throw std::runtime_error(
                    "numerical_gradient: failed to evaluate expression with dual numbers");
            }
            gradient_values.push_back(result.derivative());
        }

        // Format output
        std::ostringstream out;
        out << "[";
        for (std::size_t i = 0; i < gradient_values.size(); ++i) {
            if (i != 0) {
                out << ", ";
            }
            out << gradient_values[i];
        }
        out << "]";
        *output = out.str();
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

    if (command == "divergence" || command == "div" || command == "curl" || command == "curl_2d") {
        if (arguments.size() < 2) {
            throw std::runtime_error(command + " expects a bracketed expression list and variable names");
        }
        const std::vector<SymbolicExpression> components =
            ctx.parse_symbolic_expression_list(arguments[0]);
        std::vector<std::string> fallback_variables;
        for (const SymbolicExpression& comp : components) {
            const std::vector<std::string> ids = comp.identifier_variables();
            fallback_variables.insert(fallback_variables.end(), ids.begin(), ids.end());
        }
        std::sort(fallback_variables.begin(), fallback_variables.end());
        fallback_variables.erase(std::unique(fallback_variables.begin(), fallback_variables.end()), fallback_variables.end());

        const std::vector<std::string> variables =
            ctx.parse_symbolic_variable_arguments(arguments, 1, fallback_variables);

        if (command == "curl") {
            *output = symbolic_vector_to_string(SymbolicExpression::curl(components, variables));
        } else if (command == "curl_2d") {
            *output = SymbolicExpression::curl_2d(components, variables).simplify().to_string();
        } else {
            *output = SymbolicExpression::divergence(components, variables).simplify().to_string();
        }
        return true;
    }

    if (command == "laplacian") {
        if (arguments.empty()) {
            throw std::runtime_error("laplacian expects an expression and optional variable names");
        }
        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            ctx.parse_symbolic_variable_arguments(arguments, 1, expression.identifier_variables());
        *output = expression.laplacian(variables).simplify().to_string();
        return true;
    }

    if (command == "implicit_diff") {
        if (arguments.size() != 3) {
            throw std::runtime_error("implicit_diff expects (F, y, x)");
        }
        const std::string y_var = trim_copy(arguments[1]);
        const std::string x_var = trim_copy(arguments[2]);
        if (!is_identifier_text(y_var) || !is_identifier_text(x_var)) {
            throw std::runtime_error("implicit_diff variables must be identifiers");
        }
        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &expression);
        const SymbolicExpression numerator =
            -expression.derivative(x_var).simplify();
        const SymbolicExpression denominator =
            expression.derivative(y_var).simplify();
        *output = (numerator / denominator).simplify().to_string();
        return true;
    }

    if (command == "param_deriv") {
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "param_deriv expects (x(t), y(t), t[, order])");
        }
        const std::string parameter = trim_copy(arguments[2]);
        if (!is_identifier_text(parameter)) {
            throw std::runtime_error("param_deriv parameter must be an identifier");
        }
        int order = 1;
        if (arguments.size() == 4) {
            const double order_value = ctx.parse_decimal(arguments[3]);
            if (!mymath::isfinite(order_value) || mymath::floor(order_value) != order_value ||
                order_value < 1.0) {
                throw std::runtime_error("param_deriv order must be a positive integer");
            }
            order = static_cast<int>(order_value);
        }
        SymbolicExpression x_expr = SymbolicExpression::parse(arguments[0]);
        SymbolicExpression y_expr = SymbolicExpression::parse(arguments[1]);
        SymbolicExpression dx_dt = x_expr.derivative(parameter).simplify();
        SymbolicExpression result =
            (y_expr.derivative(parameter).simplify() / dx_dt).simplify();
        for (int i = 2; i <= order; ++i) {
            result = (result.derivative(parameter).simplify() / dx_dt).simplify();
        }
        *output = result.simplify().to_string();
        return true;
    }

    if (command == "directional") {
        if (arguments.size() < 3 || ((arguments.size() - 1) % 2) != 0) {
            throw std::runtime_error(
                "directional expects (expr, vars..., direction_components...)");
        }
        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &expression);
        const std::size_t dimension = (arguments.size() - 1) / 2;
        std::vector<std::string> variables;
        variables.reserve(dimension);
        for (std::size_t i = 0; i < dimension; ++i) {
            const std::string variable = trim_copy(arguments[1 + i]);
            if (!is_identifier_text(variable)) {
                throw std::runtime_error(
                    "directional variable arguments must be identifiers");
            }
            variables.push_back(variable);
        }

        std::vector<SymbolicExpression> direction;
        direction.reserve(dimension);
        SymbolicExpression norm_squared = SymbolicExpression::number(0.0);
        for (std::size_t i = 0; i < dimension; ++i) {
            SymbolicExpression component =
                SymbolicExpression::parse(arguments[1 + dimension + i]).simplify();
            direction.push_back(component);
            norm_squared = (norm_squared + component * component).simplify();
        }
        const SymbolicExpression norm =
            SymbolicExpression::parse("sqrt(" + norm_squared.to_string() + ")").simplify();
        const std::vector<SymbolicExpression> gradient =
            expression.gradient(variables);
        SymbolicExpression result = SymbolicExpression::number(0.0);
        for (std::size_t i = 0; i < dimension; ++i) {
            result = (result + gradient[i] * direction[i] / norm).simplify();
        }
        *output = result.simplify().to_string();
        return true;
    }

    // ============================================================================
    // Line Integral Commands
    // ============================================================================

    if (command == "line_integral") {
        // line_integral(f_or_F, [x(t), y(t), ...], t, a, b, [coord_vars])
        // Optional coord_vars: list of coordinate variable names (default: x, y, z, ...)
        if (arguments.size() < 5) {
            throw std::runtime_error("line_integral expects (f_or_F, [x(t), y(t), ...], t, a, b, [coord_vars])");
        }
        const std::vector<SymbolicExpression> field_components = ctx.parse_symbolic_expression_list(arguments[0]);
        const std::vector<SymbolicExpression> curve_components = ctx.parse_symbolic_expression_list(arguments[1]);
        const std::string param = trim_copy(arguments[2]);
        if (!is_identifier_text(param)) throw std::runtime_error("line_integral parameter must be an identifier");
        const double a = ctx.parse_decimal(arguments[3]);
        const double b = ctx.parse_decimal(arguments[4]);
        const std::size_t dim = curve_components.size();

        // Parse optional coordinate variable names
        std::vector<std::string> coord_vars;
        if (arguments.size() > 5) {
            const std::vector<SymbolicExpression> coord_exprs = ctx.parse_symbolic_expression_list(arguments[5]);
            for (const auto& expr : coord_exprs) {
                std::string var_name = expr.to_string();
                // Clean up the variable name
                if (!var_name.empty() && is_identifier_text(var_name)) {
                    coord_vars.push_back(var_name);
                }
            }
        }
        // Fill with default names if not enough provided
        const std::vector<std::string> default_vars = {"x", "y", "z", "u", "v", "w", "a", "b", "c", "d"};
        while (coord_vars.size() < dim) {
            coord_vars.push_back(default_vars[coord_vars.size()]);
        }

        std::vector<SymbolicExpression> dr_dt(dim);
        for (std::size_t i = 0; i < dim; ++i) dr_dt[i] = curve_components[i].derivative(param).simplify();

        if (field_components.size() > 1) {
            if (field_components.size() != dim) throw std::runtime_error("line_integral: field and curve dimensions must match");
            SymbolicExpression integrand = SymbolicExpression::number(0.0);
            for (std::size_t i = 0; i < dim; ++i) {
                SymbolicExpression Fi = field_components[i];
                for (std::size_t j = 0; j < dim; ++j) Fi = Fi.substitute(coord_vars[j], curve_components[j]).simplify();
                integrand = (integrand + Fi * dr_dt[i]).simplify();
            }
            SymbolicExpression antideriv;
            bool symbolic_success = false;
            try {
                antideriv = integrand.integral(param).simplify();
                symbolic_success = antideriv.to_string().find("integral(") == std::string::npos;
            } catch (...) { symbolic_success = false; }

            if (symbolic_success) {
                SymbolicExpression res = (antideriv.substitute(param, SymbolicExpression::number(b)) - antideriv.substitute(param, SymbolicExpression::number(a))).simplify();
                double n; if (res.is_number(&n)) *output = format_decimal(ctx.normalize_result(n)); else *output = res.to_string();
            } else {
                const auto eval = ctx.build_scoped_evaluator(integrand.to_string());
                MultivariableIntegrator integrator([eval, param](const std::vector<double>& p) { return eval({{param, p[0]}}); });
                *output = format_decimal(ctx.normalize_result(integrator.integrate({[a, b](const std::vector<double>&) { return std::make_pair(a, b); }}, {256})));
            }
        } else {
            std::string var; SymbolicExpression f; ctx.resolve_symbolic(arguments[0], false, &var, &f);
            SymbolicExpression speed_sq = SymbolicExpression::number(0.0);
            for (std::size_t i = 0; i < dim; ++i) speed_sq = (speed_sq + dr_dt[i] * dr_dt[i]).simplify();
            SymbolicExpression speed = SymbolicExpression::parse("sqrt(" + speed_sq.to_string() + ")").simplify();
            SymbolicExpression f_r = f;
            for (std::size_t i = 0; i < dim; ++i) f_r = f_r.substitute(coord_vars[i], curve_components[i]).simplify();
            SymbolicExpression integrand = (f_r * speed).simplify();
            SymbolicExpression antideriv; bool symbolic_success = false;
            try {
                antideriv = integrand.integral(param).simplify();
                symbolic_success = antideriv.to_string().find("integral(") == std::string::npos;
            } catch (...) { symbolic_success = false; }
            if (symbolic_success) {
                SymbolicExpression res = (antideriv.substitute(param, SymbolicExpression::number(b)) - antideriv.substitute(param, SymbolicExpression::number(a))).simplify();
                double n; if (res.is_number(&n)) *output = format_decimal(ctx.normalize_result(n)); else *output = res.to_string();
            } else {
                const auto f_eval = ctx.build_scoped_evaluator(f.to_string());
                std::vector<std::function<double(const std::vector<std::pair<std::string, double>>&)>> c_evals, d_evals;
                for (const auto& c : curve_components) c_evals.push_back(ctx.build_scoped_evaluator(c.to_string()));
                for (const auto& d : dr_dt) d_evals.push_back(ctx.build_scoped_evaluator(d.to_string()));
                MultivariableIntegrator integrator([f_eval, c_evals, d_evals, coord_vars, param](const std::vector<double>& p) {
                    const std::vector<std::pair<std::string, double>> p_s = {{param, p[0]}};
                    std::vector<std::pair<std::string, double>> f_s = p_s;
                    for (std::size_t i = 0; i < c_evals.size(); ++i) f_s.push_back({coord_vars[i], c_evals[i](p_s)});
                    double s_v = 0.0; for (const auto& d : d_evals) { double v = d(p_s); s_v += v * v; }
                    return f_eval(f_s) * mymath::sqrt(s_v);
                });
                *output = format_decimal(ctx.normalize_result(integrator.integrate({[a, b](const std::vector<double>&) { return std::make_pair(a, b); }}, {256})));
            }
        }
        return true;
    }

    if (command == "line_integral_vector") {
        // line_integral_vector([F_x, F_y, ...], [x(t), y(t), ...], t, a, b)
        // Computes ∫_C F(r(t)) · r'(t) dt
        if (arguments.size() < 5) {
            throw std::runtime_error(
                "line_integral_vector expects ([F_x, F_y, ...], [x(t), y(t), ...], t, a, b)");
        }

        // Parse the vector field components [F_x, F_y, ...]
        const std::vector<SymbolicExpression> field_components =
            ctx.parse_symbolic_expression_list(arguments[0]);

        // Parse the parametric curve components [x(t), y(t), ...]
        const std::vector<SymbolicExpression> curve_components =
            ctx.parse_symbolic_expression_list(arguments[1]);

        if (field_components.size() != curve_components.size()) {
            throw std::runtime_error(
                "line_integral_vector: field and curve dimensions must match");
        }

        // Parse parameter name
        const std::string param = trim_copy(arguments[2]);
        if (!is_identifier_text(param)) {
            throw std::runtime_error("line_integral_vector parameter must be an identifier");
        }

        // Parse bounds
        const double a = ctx.parse_decimal(arguments[3]);
        const double b = ctx.parse_decimal(arguments[4]);

        const std::size_t dim = field_components.size();

        // Compute r'(t) components
        std::vector<SymbolicExpression> dr_dt(dim);
        for (std::size_t i = 0; i < dim; ++i) {
            dr_dt[i] = curve_components[i].derivative(param).simplify();
        }

        // Collect coordinate variable names (default: x, y, z, ...)
        std::vector<std::string> coord_vars;
        const char* default_vars[] = {"x", "y", "z", "u", "v", "w"};
        for (std::size_t i = 0; i < dim && i < 6; ++i) {
            coord_vars.push_back(default_vars[i]);
        }

        // Substitute r(t) into F and compute dot product F(r(t)) · r'(t)
        SymbolicExpression integrand = SymbolicExpression::number(0.0);
        for (std::size_t i = 0; i < dim; ++i) {
            SymbolicExpression F_i = field_components[i];
            for (std::size_t j = 0; j < dim; ++j) {
                F_i = F_i.substitute(coord_vars[j], curve_components[j]).simplify();
            }
            integrand = (integrand + F_i * dr_dt[i]).simplify();
        }

        // Try symbolic integration first, then fall back to numerical integration.
        SymbolicExpression antideriv;
        bool symbolic_success = false;
        try {
            antideriv = integrand.integral(param).simplify();
            const std::string antideriv_str = antideriv.to_string();
            symbolic_success = antideriv_str.find("integral(") == std::string::npos;
        } catch (const std::exception&) {
            symbolic_success = false;
        }

        if (symbolic_success) {
            SymbolicExpression result_at_b = antideriv.substitute(param, SymbolicExpression::number(b)).simplify();
            SymbolicExpression result_at_a = antideriv.substitute(param, SymbolicExpression::number(a)).simplify();
            SymbolicExpression result = (result_at_b - result_at_a).simplify();

            double numeric_result;
            if (result.is_number(&numeric_result)) {
                *output = format_decimal(ctx.normalize_result(numeric_result));
            } else {
                *output = result.to_string();
            }
        } else {
            // Fall back to numerical integration
            std::vector<std::function<double(const std::vector<std::pair<std::string, double>>&)>> field_evaluators;
            field_evaluators.reserve(dim);
            std::vector<std::function<double(const std::vector<std::pair<std::string, double>>&)>> curve_evaluators;
            curve_evaluators.reserve(dim);
            std::vector<std::function<double(const std::vector<std::pair<std::string, double>>&)>> derivative_evaluators;
            derivative_evaluators.reserve(dim);
            for (std::size_t i = 0; i < dim; ++i) {
                field_evaluators.push_back(ctx.build_scoped_evaluator(field_components[i].to_string()));
                curve_evaluators.push_back(ctx.build_scoped_evaluator(curve_components[i].to_string()));
                derivative_evaluators.push_back(ctx.build_scoped_evaluator(dr_dt[i].to_string()));
            }
            MultivariableIntegrator integrator([field_evaluators, curve_evaluators, derivative_evaluators, coord_vars, param](
                const std::vector<double>& point) {
                const std::vector<std::pair<std::string, double>> param_scope = {{param, point[0]}};
                std::vector<std::pair<std::string, double>> field_scope = param_scope;
                for (std::size_t i = 0; i < curve_evaluators.size(); ++i) {
                    field_scope.push_back({coord_vars[i], curve_evaluators[i](param_scope)});
                }
                double value = 0.0;
                for (std::size_t i = 0; i < field_evaluators.size(); ++i) {
                    value += field_evaluators[i](field_scope) *
                             derivative_evaluators[i](param_scope);
                }
                return value;
            });
            double numeric_result = integrator.integrate({
                [a, b](const std::vector<double>&) {
                    return std::make_pair(a, b);
                }
            }, {256});
            *output = format_decimal(ctx.normalize_result(numeric_result));
        }
        return true;
    }

    // ============================================================================
    // Surface Integral Commands
    // ============================================================================

    if (command == "surface_integral" && !arguments.empty() &&
        ctx.parse_symbolic_expression_list(arguments[0]).size() == 1) {
        // surface_integral(f, [x(u,v), y(u,v), z(u,v)], u, a, b, v, c, d)
        // Computes ∬_S f(r(u,v)) * |r_u × r_v| du dv
        if (arguments.size() < 8) {
            throw std::runtime_error(
                "surface_integral expects (f_or_F, [x(u,v), y(u,v), z(u,v)], u, a, b, v, c, d)");
        }

        // Parse the scalar field f
        std::string var_name;
        SymbolicExpression f_expr;
        ctx.resolve_symbolic(arguments[0], false, &var_name, &f_expr);

        // Parse the parametric surface components [x(u,v), y(u,v), z(u,v)]
        const std::vector<SymbolicExpression> components =
            ctx.parse_symbolic_expression_list(arguments[1]);

        if (components.size() != 3) {
            throw std::runtime_error(
                "surface_integral_scalar currently only supports 3D surfaces");
        }

        // Parse parameter names
        const std::string u_param = trim_copy(arguments[2]);
        const std::string v_param = trim_copy(arguments[5]);
        if (!is_identifier_text(u_param) || !is_identifier_text(v_param)) {
            throw std::runtime_error("surface_integral_scalar parameters must be identifiers");
        }

        // Parse bounds
        const double u_a = ctx.parse_decimal(arguments[3]);
        const double u_b = ctx.parse_decimal(arguments[4]);
        const double v_c = ctx.parse_decimal(arguments[6]);
        const double v_d = ctx.parse_decimal(arguments[7]);

        // Compute r_u and r_v (partial derivatives)
        std::vector<SymbolicExpression> r_u(3), r_v(3);
        for (int i = 0; i < 3; ++i) {
            r_u[i] = components[i].derivative(u_param).simplify();
            r_v[i] = components[i].derivative(v_param).simplify();
        }

        // Compute cross product r_u × r_v
        // (r_u × r_v)_x = r_u[1] * r_v[2] - r_u[2] * r_v[1]
        // (r_u × r_v)_y = r_u[2] * r_v[0] - r_u[0] * r_v[2]
        // (r_u × r_v)_z = r_u[0] * r_v[1] - r_u[1] * r_v[0]
        SymbolicExpression cross_x = (r_u[1] * r_v[2] - r_u[2] * r_v[1]).simplify();
        SymbolicExpression cross_y = (r_u[2] * r_v[0] - r_u[0] * r_v[2]).simplify();
        SymbolicExpression cross_z = (r_u[0] * r_v[1] - r_u[1] * r_v[0]).simplify();

        // |r_u × r_v| = sqrt(cross_x^2 + cross_y^2 + cross_z^2)
        SymbolicExpression cross_norm_squared = (
            cross_x * cross_x + cross_y * cross_y + cross_z * cross_z).simplify();
        SymbolicExpression cross_norm = SymbolicExpression::parse(
            "sqrt(" + cross_norm_squared.to_string() + ")").simplify();

        // Substitute r(u,v) into f: f(r(u,v))
        const std::vector<std::string> coord_vars = {"x", "y", "z"};
        SymbolicExpression f_of_r = f_expr;
        for (int i = 0; i < 3; ++i) {
            f_of_r = f_of_r.substitute(coord_vars[i], components[i]).simplify();
        }

        // Integrand: f(r(u,v)) * |r_u × r_v|
        SymbolicExpression integrand = (f_of_r * cross_norm).simplify();

        // Use the improved double_integral which supports adaptive integration
        integration_ops::IntegrationContext int_ctx;
        int_ctx.parse_decimal = [&](const std::string& s) { return ctx.parse_decimal(s); };
        int_ctx.build_scoped_evaluator = [&](const std::string& s) { return ctx.build_scoped_evaluator(s); };
        int_ctx.normalize_result = [&](double d) { return ctx.normalize_result(d); };

        double result = integration_ops::double_integral(
            int_ctx, integrand.to_string(),
            u_param, u_a, u_b,
            v_param, std::to_string(v_c), std::to_string(v_d),
            50, 50);

        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    if (command == "surface_integral") {
        // surface_integral([F_x, F_y, F_z], [x(u,v), y(u,v), z(u,v)], u, a, b, v, c(u), d(u), [orientation])
        // or surface_integral(f, [x(u,v), y(u,v), z(u,v)], u, a, b, v, c(u), d(u)) for scalar integral
        // c(u) and d(u) can be expressions depending on u, or constants
        // orientation: "outward" (default) or "inward" for flux integrals
        if (arguments.size() < 8) {
            throw std::runtime_error(
                "surface_integral expects (f_or_F, [x(u,v), y(u,v), z(u,v)], u, a, b, v, c, d, [orientation])");
        }

        const std::vector<SymbolicExpression> field_components =
            ctx.parse_symbolic_expression_list(arguments[0]);

        const std::vector<SymbolicExpression> surface_components =
            ctx.parse_symbolic_expression_list(arguments[1]);
        const std::size_t dim = surface_components.size();

        // Parse parameter names
        const std::string u_param = trim_copy(arguments[2]);
        const std::string v_param = trim_copy(arguments[5]);
        if (!is_identifier_text(u_param) || !is_identifier_text(v_param)) {
            throw std::runtime_error("surface_integral parameters must be identifiers");
        }

        // Parse bounds - v bounds can depend on u
        const double u_a = ctx.parse_decimal(arguments[3]);
        const double u_b = ctx.parse_decimal(arguments[4]);
        const std::string v_c_expr = arguments[6];  // Can be constant or expression in u
        const std::string v_d_expr = arguments[7];  // Can be constant or expression in u

        // Parse optional orientation
        std::string orientation = "outward";
        if (arguments.size() > 8) {
            orientation = trim_copy(arguments[8]);
        }

        // Compute r_u and r_v (partial derivatives)
        std::vector<SymbolicExpression> r_u(dim), r_v(dim);
        for (std::size_t i = 0; i < dim; ++i) {
            r_u[i] = surface_components[i].derivative(u_param).simplify();
            r_v[i] = surface_components[i].derivative(v_param).simplify();
        }

        // Build integrand based on field type
        SymbolicExpression integrand;

        if (field_components.size() > 1) {
            // Vector field: flux integral F · n dS = F · (r_u × r_v) du dv
            if (field_components.size() != dim) {
                throw std::runtime_error("surface_integral: field and surface dimensions must match");
            }
            if (dim != 3) {
                throw std::runtime_error("surface_integral flux currently only supports 3D surfaces");
            }

            // Compute cross product r_u × r_v
            SymbolicExpression cross_x = (r_u[1] * r_v[2] - r_u[2] * r_v[1]).simplify();
            SymbolicExpression cross_y = (r_u[2] * r_v[0] - r_u[0] * r_v[2]).simplify();
            SymbolicExpression cross_z = (r_u[0] * r_v[1] - r_u[1] * r_v[0]).simplify();

            // Apply orientation
            double sign = (orientation == "inward") ? -1.0 : 1.0;

            // Substitute r(u,v) into F
            const std::vector<std::string> coord_vars = {"x", "y", "z"};
            std::vector<SymbolicExpression> F_of_r(3);
            for (int i = 0; i < 3; ++i) {
                F_of_r[i] = field_components[i];
                for (int j = 0; j < 3; ++j) {
                    F_of_r[i] = F_of_r[i].substitute(coord_vars[j], surface_components[j]).simplify();
                }
            }

            // Integrand: F(r(u,v)) · (r_u × r_v)
            integrand = (F_of_r[0] * cross_x + F_of_r[1] * cross_y + F_of_r[2] * cross_z).simplify();
            if (sign < 0) {
                integrand = (SymbolicExpression::number(-1.0) * integrand).simplify();
            }
        } else {
            // Scalar field: surface integral f dS = f * |r_u × r_v| du dv
            std::string var;
            SymbolicExpression f_expr;
            ctx.resolve_symbolic(arguments[0], false, &var, &f_expr);

            if (dim == 3) {
                // Compute |r_u × r_v|
                SymbolicExpression cross_x = (r_u[1] * r_v[2] - r_u[2] * r_v[1]).simplify();
                SymbolicExpression cross_y = (r_u[2] * r_v[0] - r_u[0] * r_v[2]).simplify();
                SymbolicExpression cross_z = (r_u[0] * r_v[1] - r_u[1] * r_v[0]).simplify();
                SymbolicExpression cross_norm_sq = (cross_x * cross_x + cross_y * cross_y + cross_z * cross_z).simplify();
                SymbolicExpression cross_norm = SymbolicExpression::parse("sqrt(" + cross_norm_sq.to_string() + ")").simplify();

                const std::vector<std::string> coord_vars = {"x", "y", "z"};
                SymbolicExpression f_of_r = f_expr;
                for (int i = 0; i < 3; ++i) {
                    f_of_r = f_of_r.substitute(coord_vars[i], surface_components[i]).simplify();
                }
                integrand = (f_of_r * cross_norm).simplify();
            } else {
                // General dimension: compute sqrt(det(g)) where g is the metric tensor
                SymbolicExpression metric_det = SymbolicExpression::number(0.0);
                for (std::size_t i = 0; i < dim; ++i) {
                    for (std::size_t j = 0; j < dim; ++j) {
                        // g_ij = r_i · r_j
                        SymbolicExpression g_ij = SymbolicExpression::number(0.0);
                        for (std::size_t k = 0; k < dim; ++k) {
                            // For 2D surface, we compute the 2x2 metric tensor determinant
                        }
                    }
                }
                // For 2D surface in nD: area element = sqrt(|r_u|^2 * |r_v|^2 - (r_u·r_v)^2)
                SymbolicExpression r_u_sq = SymbolicExpression::number(0.0);
                SymbolicExpression r_v_sq = SymbolicExpression::number(0.0);
                SymbolicExpression r_u_dot_r_v = SymbolicExpression::number(0.0);
                for (std::size_t i = 0; i < dim; ++i) {
                    r_u_sq = (r_u_sq + r_u[i] * r_u[i]).simplify();
                    r_v_sq = (r_v_sq + r_v[i] * r_v[i]).simplify();
                    r_u_dot_r_v = (r_u_dot_r_v + r_u[i] * r_v[i]).simplify();
                }
                SymbolicExpression area_sq = (r_u_sq * r_v_sq - r_u_dot_r_v * r_u_dot_r_v).simplify();
                SymbolicExpression area_elem = SymbolicExpression::parse("sqrt(" + area_sq.to_string() + ")").simplify();

                // Build default coordinate variable names
                std::vector<std::string> coord_vars;
                const std::vector<std::string> default_vars = {"x", "y", "z", "u", "v", "w"};
                for (std::size_t i = 0; i < dim && i < default_vars.size(); ++i) {
                    coord_vars.push_back(default_vars[i]);
                }

                SymbolicExpression f_of_r = f_expr;
                for (std::size_t i = 0; i < dim && i < coord_vars.size(); ++i) {
                    f_of_r = f_of_r.substitute(coord_vars[i], surface_components[i]).simplify();
                }
                integrand = (f_of_r * area_elem).simplify();
            }
        }

        // Use double_integral with potentially non-constant v bounds
        integration_ops::IntegrationContext int_ctx;
        int_ctx.parse_decimal = [&](const std::string& s) { return ctx.parse_decimal(s); };
        int_ctx.build_scoped_evaluator = [&](const std::string& s) { return ctx.build_scoped_evaluator(s); };
        int_ctx.normalize_result = [&](double d) { return ctx.normalize_result(d); };

        double result = integration_ops::double_integral(
            int_ctx, integrand.to_string(),
            u_param, u_a, u_b,
            v_param, v_c_expr, v_d_expr,
            50, 50);

        *output = format_decimal(ctx.normalize_result(result));
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
        if (arguments.empty()) {
            throw std::runtime_error(
                "integral expects at least 1 argument for symbolic indefinite integral");
        }

        bool symbolic_integral = true;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            symbolic_integral =
                symbolic_integral &&
                is_identifier_text(trim_copy(arguments[i])) &&
                !is_infinity_literal(arguments[i]);
        }
        if (symbolic_integral) {
            std::string variable_name;
            SymbolicExpression expression;
            ctx.resolve_symbolic(arguments[0],
                                 arguments.size() == 1,
                                 &variable_name,
                                 &expression);
            if (arguments.size() >= 2) {
                variable_name = trim_copy(arguments[1]);
            }

            if (arguments.size() > 2) {
                SymbolicExpression integrated = expression;
                for (std::size_t i = 1; i < arguments.size(); ++i) {
                    const std::string current_variable = trim_copy(arguments[i]);
                    RischAlgorithm::IntegrationResult risch_result =
                        RischAlgorithm::integrate_full(integrated, current_variable);
                    if (!risch_result.success ||
                        risch_result.type != IntegralType::kElementary) {
                        throw std::runtime_error("Integration failed: Risch could not integrate with respect to " +
                                                 current_variable);
                    }
                    integrated = risch_result.value.simplify();
                }
                *output = integrated.simplify().to_string() + " + C";
                return true;
            }

            // 使用 IntegrationEngine 进行积分
            IntegrationEngine engine;
            IntegrationResult result = engine.integrate(expression, variable_name);

            if (result.success) {
                *output = result.value.simplify().to_string() + " + C";
                return true;
            }

            // 如果 IntegrationEngine 失败，回退到直接调用 integral()
            try {
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
            } catch (const std::runtime_error& e) {
                throw std::runtime_error("Integration failed: " + std::string(e.what()));
            }
        }

        if (arguments.size() != 2 && arguments.size() != 3 && arguments.size() != 4) {
             throw std::runtime_error(
                "numerical integral expects 2 arguments for indefinite value, "
                "3 for definite integral, or 4 for anchor and constant");
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

    // ============================================================================
    // Vector Field Theorem Commands
    // ============================================================================

    if (command == "greens_theorem") {
        // greens_theorem([Fx, Fy], [x(t), y(t)], t, a, b)
        // Green's theorem: ∮_C F·dr = ∬_D (∂Fy/∂x - ∂Fx/∂y) dA
        // Computes the line integral using the double integral of curl_2d
        if (arguments.size() < 5) {
            throw std::runtime_error(
                "greens_theorem expects ([Fx, Fy], [x(t), y(t)], t, a, b)");
        }

        // Parse the vector field [Fx, Fy]
        const std::vector<SymbolicExpression> field =
            ctx.parse_symbolic_expression_list(arguments[0]);
        if (field.size() != 2) {
            throw std::runtime_error("greens_theorem requires a 2D vector field");
        }

        // Parse the parametric curve [x(t), y(t)]
        const std::vector<SymbolicExpression> curve =
            ctx.parse_symbolic_expression_list(arguments[1]);
        if (curve.size() != 2) {
            throw std::runtime_error("greens_theorem requires a 2D curve");
        }

        // Parse parameter
        const std::string param = trim_copy(arguments[2]);
        if (!is_identifier_text(param)) {
            throw std::runtime_error("greens_theorem parameter must be an identifier");
        }

        // Parse bounds
        const double a = ctx.parse_decimal(arguments[3]);
        const double b = ctx.parse_decimal(arguments[4]);

        // Compute curl_2d = ∂Fy/∂x - ∂Fx/∂y
        const std::vector<std::string> vars = {"x", "y"};
        SymbolicExpression curl_2d_expr = SymbolicExpression::curl_2d(field, vars);

        // Compute the line integral directly: ∮ F·dr = ∫ (Fx*dx/dt + Fy*dy/dt) dt
        // This is the verification method
        SymbolicExpression dx_dt = curve[0].derivative(param).simplify();
        SymbolicExpression dy_dt = curve[1].derivative(param).simplify();

        // Substitute curve into field
        SymbolicExpression Fx_at_curve = field[0].substitute("x", curve[0]).substitute("y", curve[1]).simplify();
        SymbolicExpression Fy_at_curve = field[1].substitute("x", curve[0]).substitute("y", curve[1]).simplify();

        // Integrand: Fx * dx/dt + Fy * dy/dt
        SymbolicExpression integrand = (Fx_at_curve * dx_dt + Fy_at_curve * dy_dt).simplify();

        // Try symbolic integration
        SymbolicExpression antideriv;
        bool symbolic_success = false;
        try {
            antideriv = integrand.integral(param).simplify();
            const std::string antideriv_str = antideriv.to_string();
            symbolic_success = antideriv_str.find("integral(") == std::string::npos;
        } catch (const std::exception&) {
            symbolic_success = false;
        }

        if (symbolic_success) {
            SymbolicExpression result_at_b = antideriv.substitute(param, SymbolicExpression::number(b)).simplify();
            SymbolicExpression result_at_a = antideriv.substitute(param, SymbolicExpression::number(a)).simplify();
            SymbolicExpression result = (result_at_b - result_at_a).simplify();

            double numeric_result;
            if (result.is_number(&numeric_result)) {
                *output = format_decimal(ctx.normalize_result(numeric_result));
            } else {
                *output = result.to_string();
            }
        } else {
            // Fall back to numerical integration
            const auto evaluate_integrand = ctx.build_scoped_evaluator(integrand.to_string());
            MultivariableIntegrator integrator([&evaluate_integrand, param](const std::vector<double>& point) {
                return evaluate_integrand({{param, point[0]}});
            });

            std::vector<MultivariableIntegrator::BoundFunc> bounds;
            bounds.push_back([&a, &b](const std::vector<double>&) {
                return std::make_pair(a, b);
            });

            double result = integrator.integrate(bounds, {100});
            *output = format_decimal(ctx.normalize_result(result));
        }
        return true;
    }

    if (command == "stokes_theorem") {
        // stokes_theorem([Fx, Fy, Fz], [x(u,v), y(u,v), z(u,v)], u, a, b, v, c, d)
        // Stokes' theorem: ∮_C F·dr = ∬_S curl(F)·n dS
        // Computes the surface integral of curl(F) through the surface
        if (arguments.size() < 8) {
            throw std::runtime_error(
                "stokes_theorem expects ([Fx, Fy, Fz], [x(u,v), y(u,v), z(u,v)], u, a, b, v, c, d)");
        }

        // Parse the vector field
        const std::vector<SymbolicExpression> field =
            ctx.parse_symbolic_expression_list(arguments[0]);
        if (field.size() != 3) {
            throw std::runtime_error("stokes_theorem requires a 3D vector field");
        }

        // Parse the parametric surface
        const std::vector<SymbolicExpression> surface =
            ctx.parse_symbolic_expression_list(arguments[1]);
        if (surface.size() != 3) {
            throw std::runtime_error("stokes_theorem requires a 3D surface");
        }

        // Parse parameters
        const std::string u_param = trim_copy(arguments[2]);
        const std::string v_param = trim_copy(arguments[5]);
        if (!is_identifier_text(u_param) || !is_identifier_text(v_param)) {
            throw std::runtime_error("stokes_theorem parameters must be identifiers");
        }

        // Parse bounds
        const double u_a = ctx.parse_decimal(arguments[3]);
        const double u_b = ctx.parse_decimal(arguments[4]);
        const double v_c = ctx.parse_decimal(arguments[6]);
        const double v_d = ctx.parse_decimal(arguments[7]);

        // Compute curl(F)
        const std::vector<std::string> vars = {"x", "y", "z"};
        std::vector<SymbolicExpression> curl_F = SymbolicExpression::curl(field, vars);

        // Compute r_u and r_v
        std::vector<SymbolicExpression> r_u(3), r_v(3);
        for (int i = 0; i < 3; ++i) {
            r_u[i] = surface[i].derivative(u_param).simplify();
            r_v[i] = surface[i].derivative(v_param).simplify();
        }

        // Compute cross product r_u × r_v (this is n * dS)
        SymbolicExpression cross_x = (r_u[1] * r_v[2] - r_u[2] * r_v[1]).simplify();
        SymbolicExpression cross_y = (r_u[2] * r_v[0] - r_u[0] * r_v[2]).simplify();
        SymbolicExpression cross_z = (r_u[0] * r_v[1] - r_u[1] * r_v[0]).simplify();

        // Substitute surface into curl(F)
        std::vector<SymbolicExpression> curl_at_surface(3);
        for (int i = 0; i < 3; ++i) {
            curl_at_surface[i] = curl_F[i];
            for (int j = 0; j < 3; ++j) {
                curl_at_surface[i] = curl_at_surface[i].substitute(vars[j], surface[j]).simplify();
            }
        }

        // Integrand: curl(F) · (r_u × r_v)
        SymbolicExpression integrand = (
            curl_at_surface[0] * cross_x +
            curl_at_surface[1] * cross_y +
            curl_at_surface[2] * cross_z).simplify();

        // Use the improved double_integral (adaptive)
        integration_ops::IntegrationContext int_ctx;
        int_ctx.parse_decimal = [&](const std::string& s) { return ctx.parse_decimal(s); };
        int_ctx.build_scoped_evaluator = [&](const std::string& s) { return ctx.build_scoped_evaluator(s); };
        int_ctx.normalize_result = [&](double d) { return ctx.normalize_result(d); };

        double result = integration_ops::double_integral(
            int_ctx, integrand.to_string(),
            u_param, u_a, u_b,
            v_param, std::to_string(v_c), std::to_string(v_d),
            50, 50);

        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    if (command == "divergence_theorem") {
        // Supports two modes:
        // 1. Surface integral mode: ([Fx, Fy, Fz], [x,y,z], u, a, b, v, c, d)
        // 2. Volume integral mode: ([Fx, Fy, Fz], x, x0, x1, y, y0, y1, z, z0, z1)
        
        bool is_volume_mode = (arguments.size() == 10 && is_identifier_text(arguments[1]));

        if (is_volume_mode) {
            const std::vector<SymbolicExpression> field = ctx.parse_symbolic_expression_list(arguments[0]);
            if (field.size() != 3) throw std::runtime_error("divergence_theorem requires a 3D vector field");

            const std::string x_var = arguments[1];
            const double x0 = ctx.parse_decimal(arguments[2]);
            const double x1 = ctx.parse_decimal(arguments[3]);
            const std::string y_var = arguments[4];
            const std::string y0 = arguments[5];
            const std::string y1 = arguments[6];
            const std::string z_var = arguments[7];
            const std::string z0 = arguments[8];
            const std::string z1 = arguments[9];

            SymbolicExpression div_F = SymbolicExpression::divergence(field, {x_var, y_var, z_var});
            
            integration_ops::IntegrationContext int_ctx;
            int_ctx.parse_decimal = [&](const std::string& s) { return ctx.parse_decimal(s); };
            int_ctx.build_scoped_evaluator = [&](const std::string& s) { return ctx.build_scoped_evaluator(s); };
            int_ctx.normalize_result = [&](double d) { return ctx.normalize_result(d); };

            double result = integration_ops::triple_integral(
                int_ctx, div_F.to_string(),
                x_var, x0, x1, y_var, y0, y1, z_var, z0, z1,
                20, 20, 20);

            *output = format_decimal(ctx.normalize_result(result));
            return true;
        } else {
            if (arguments.size() < 8) {
                throw std::runtime_error(
                    "divergence_theorem expects surface or volume integral arguments");
            }

            const std::vector<SymbolicExpression> field = ctx.parse_symbolic_expression_list(arguments[0]);
            const std::vector<SymbolicExpression> surface = ctx.parse_symbolic_expression_list(arguments[1]);
            const std::string u_param = trim_copy(arguments[2]);
            const std::string v_param = trim_copy(arguments[5]);
            const double u_a = ctx.parse_decimal(arguments[3]);
            const double u_b = ctx.parse_decimal(arguments[4]);
            const double v_c = ctx.parse_decimal(arguments[6]);
            const double v_d = ctx.parse_decimal(arguments[7]);

            std::vector<SymbolicExpression> r_u(3), r_v(3);
            for (int i = 0; i < 3; ++i) {
                r_u[i] = surface[i].derivative(u_param).simplify();
                r_v[i] = surface[i].derivative(v_param).simplify();
            }
            SymbolicExpression nx = (r_u[1] * r_v[2] - r_u[2] * r_v[1]).simplify();
            SymbolicExpression ny = (r_u[2] * r_v[0] - r_u[0] * r_v[2]).simplify();
            SymbolicExpression nz = (r_u[0] * r_v[1] - r_u[1] * r_v[0]).simplify();

            std::vector<SymbolicExpression> F_at_surface(3);
            const std::vector<std::string> vars = {"x", "y", "z"};
            for (int i = 0; i < 3; ++i) {
                F_at_surface[i] = field[i];
                for (int j = 0; j < 3; ++j) {
                    F_at_surface[i] = F_at_surface[i].substitute(vars[j], surface[j]).simplify();
                }
            }

            SymbolicExpression integrand = (
                F_at_surface[0] * nx + F_at_surface[1] * ny + F_at_surface[2] * nz).simplify();

            integration_ops::IntegrationContext int_ctx;
            int_ctx.parse_decimal = [&](const std::string& s) { return ctx.parse_decimal(s); };
            int_ctx.build_scoped_evaluator = [&](const std::string& s) { return ctx.build_scoped_evaluator(s); };
            int_ctx.normalize_result = [&](double d) { return ctx.normalize_result(d); };

            double result = integration_ops::double_integral(
                int_ctx, integrand.to_string(),
                u_param, u_a, u_b,
                v_param, std::to_string(v_c), std::to_string(v_d),
                50, 50);

            *output = format_decimal(ctx.normalize_result(result));
            return true;
        }
    }

    if (command == "integrate_region") {
        // integrate_region(f, constraint, [x, x0, x1, y, y0, y1, ...], :samples N)
        // Integrates f over the implicit region defined by constraint(x) <= 0
        // Example: integrate_region(1, x^2 + y^2 - 1, x, -1, 1, y, -1, 1)
        //          computes the area of the unit disk
        if (arguments.size() < 4) {
            throw std::runtime_error(
                "integrate_region expects (f, constraint, x, x0, x1, y, y0, y1, ...)");
        }

        // Parse the integrand
        SymbolicExpression integrand_sym;
        std::string dummy_var;
        ctx.resolve_symbolic(arguments[0], false, &dummy_var, &integrand_sym);

        // Parse the constraint
        SymbolicExpression constraint_sym;
        ctx.resolve_symbolic(arguments[1], false, &dummy_var, &constraint_sym);

        // Parse bounds
        std::vector<std::string> var_names;
        std::vector<double> lower_bounds, upper_bounds;
        std::size_t i = 2;
        while (i + 2 < arguments.size()) {
            std::string var = trim_copy(arguments[i]);
            if (var.front() == ':') break;  // Option
            if (!is_identifier_text(var)) {
                throw std::runtime_error("integrate_region variable must be an identifier");
            }
            var_names.push_back(var);
            lower_bounds.push_back(ctx.parse_decimal(arguments[i + 1]));
            upper_bounds.push_back(ctx.parse_decimal(arguments[i + 2]));
            i += 3;
        }

        if (var_names.empty()) {
            throw std::runtime_error("integrate_region requires at least one variable");
        }

        // Parse options
        int num_samples = 50000;
        std::string method = "auto";
        while (i < arguments.size()) {
            std::string opt = trim_copy(arguments[i]);
            if (opt == ":samples" && i + 1 < arguments.size()) {
                num_samples = static_cast<int>(ctx.parse_decimal(arguments[++i]));
            } else if (opt == ":method" && i + 1 < arguments.size()) {
                method = trim_copy(arguments[++i]);
            }
            ++i;
        }

        // Build the integrand evaluator
        auto integrand_eval = ctx.build_scoped_evaluator(integrand_sym.to_string());
        auto constraint_eval = ctx.build_scoped_evaluator(constraint_sym.to_string());

        // Create the multidim function and constraint
        std::size_t dim = var_names.size();
        multidim::MultidimFunction f = [&](const std::vector<double>& point) {
            std::vector<std::pair<std::string, double>> scope;
            for (std::size_t d = 0; d < dim; ++d) {
                scope.push_back({var_names[d], point[d]});
            }
            return integrand_eval(scope);
        };

        multidim::RegionConstraint constraint = [&](const std::vector<double>& point) {
            std::vector<std::pair<std::string, double>> scope;
            for (std::size_t d = 0; d < dim; ++d) {
                scope.push_back({var_names[d], point[d]});
            }
            return constraint_eval(scope);  // <= 0 means inside region
        };

        // Build bounds
        std::vector<multidim::IntegrationBounds> bounds;
        for (std::size_t d = 0; d < dim; ++d) {
            bounds.push_back(multidim::IntegrationBounds(lower_bounds[d], upper_bounds[d]));
        }

        // Select method
        multidim::IntegrationOptions options;
        options.max_evaluations = num_samples;

        multidim::IntegrationResult result;
        if (method == "qmc" || method == "quasi_monte_carlo") {
            result = multidim::integrate_quasi_monte_carlo(f, bounds, constraint, num_samples);
        } else if (method == "mc" || method == "monte_carlo") {
            result = multidim::integrate_monte_carlo(f, bounds, constraint, num_samples);
        } else {
            // Auto: use quasi-monte carlo for low dimensions, mc for high
            result = multidim::integrate_implicit_region(f, bounds, constraint, options);
        }

        if (!result.converged) {
            *output = format_decimal(ctx.normalize_result(result.value)) + " (warning: may not have converged)";
        } else {
            *output = format_decimal(ctx.normalize_result(result.value));
        }
        return true;
    }

    return false;
}


std::string SymbolicModule::execute_args(const std::string& command,
                                        const std::vector<std::string>& args,
                                        const CoreServices& services) {
    SymbolicCommandContext ctx;
    ctx.resolve_symbolic = services.symbolic.resolve_symbolic;
    ctx.parse_symbolic_variable_arguments = services.parse_symbolic_vars;
    ctx.parse_symbolic_expression_list = services.symbolic.parse_symbolic_expr_list;
    ctx.build_analysis = services.symbolic.build_analysis;
    ctx.build_scoped_evaluator = services.evaluation.build_decimal_evaluator;
    ctx.parse_decimal = services.evaluation.parse_decimal;
    ctx.normalize_result = services.evaluation.normalize_result;

    std::string inside;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += args[i];
    }

    std::string output;
    if (handle_symbolic_command(ctx, command, inside, &output)) {
        return output;
    }
    throw std::runtime_error("Symbolic command failed: " + command);
}

std::string SymbolicModule::get_help_snippet(const std::string& topic) const {
    if (topic == "symbolic") {
        return "Symbolic Operations:\n"
               "  simplify(expr)         Simplify an algebraic expression\n"
               "  expand(expr)           Expand polynomial/algebraic expression\n"
               "  diff(expr, [var])      Symbolic derivative\n"
               "  integral(expr, [var])  Symbolic indefinite integral\n"
               "  dsolve(rhs, [x, y])    Solve simple linear ODEs";
    }
    return "";
}

}  // namespace symbolic_commands

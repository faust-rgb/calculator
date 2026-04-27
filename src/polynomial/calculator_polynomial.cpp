// ============================================================================
// 多项式操作命令实现
// ============================================================================

#include "calculator_polynomial.h"

#include "polynomial.h"
#include "mymath.h"

#include <sstream>

namespace polynomial_ops {

namespace {

// 递归构建多项式，支持嵌套调用
void build_polynomial_recursive(
    const PolynomialContext& ctx,
    const std::string& argument,
    std::string* variable_name,
    std::vector<double>* coefficients) {

    const std::string trimmed_argument = trim_copy(argument);
    std::string nested_inside;

    // 检查是否为嵌套的多项式操作
    if (split_named_call(trimmed_argument, "poly_add", &nested_inside) ||
        split_named_call(trimmed_argument, "poly_sub", &nested_inside) ||
        split_named_call(trimmed_argument, "poly_mul", &nested_inside) ||
        split_named_call(trimmed_argument, "poly_div", &nested_inside)) {

        const std::vector<std::string> nested_arguments =
            split_top_level_arguments(nested_inside);
        if (nested_arguments.size() != 2) {
            throw std::runtime_error(
                "polynomial operations expect exactly two arguments");
        }

        std::string lhs_variable;
        std::string rhs_variable;
        std::vector<double> lhs_coefficients;
        std::vector<double> rhs_coefficients;
        build_polynomial_recursive(ctx, nested_arguments[0],
                                   &lhs_variable, &lhs_coefficients);
        build_polynomial_recursive(ctx, nested_arguments[1],
                                   &rhs_variable, &rhs_coefficients);

        *variable_name = lhs_variable;
        if (trimmed_argument.rfind("poly_add", 0) == 0) {
            *coefficients = polynomial_add(lhs_coefficients, rhs_coefficients);
            return;
        }
        if (trimmed_argument.rfind("poly_sub", 0) == 0) {
            *coefficients = polynomial_subtract(lhs_coefficients, rhs_coefficients);
            return;
        }
        if (trimmed_argument.rfind("poly_div", 0) == 0) {
            const PolynomialDivisionResult division =
                polynomial_divide(lhs_coefficients, rhs_coefficients);
            bool zero_remainder = true;
            for (double coefficient : division.remainder) {
                if (!mymath::is_near_zero(coefficient, 1e-10)) {
                    zero_remainder = false;
                    break;
                }
            }
            if (!zero_remainder) {
                throw std::runtime_error(
                    "nested poly_div requires zero remainder");
            }
            *coefficients = division.quotient;
            return;
        }

        *coefficients = polynomial_multiply(lhs_coefficients, rhs_coefficients);
        return;
    }

    // 从符号表达式构建多项式
    SymbolicExpression expression =
        ctx.resolve_symbolic(trimmed_argument, variable_name);
    if (!expression.polynomial_coefficients(*variable_name, coefficients)) {
        throw std::runtime_error("custom function " + trimmed_argument +
                                 " is not a polynomial");
    }
}

}  // namespace

PolynomialData build_polynomial(const PolynomialContext& ctx,
                                const std::string& argument) {
    PolynomialData result;
    build_polynomial_recursive(ctx, argument,
                               &result.variable_name, &result.coefficients);
    return result;
}

std::string poly_add(const PolynomialData& lhs, const PolynomialData& rhs) {
    return polynomial_to_string(
        polynomial_add(lhs.coefficients, rhs.coefficients),
        lhs.variable_name);
}

std::string poly_sub(const PolynomialData& lhs, const PolynomialData& rhs) {
    return polynomial_to_string(
        polynomial_subtract(lhs.coefficients, rhs.coefficients),
        lhs.variable_name);
}

std::string poly_mul(const PolynomialData& lhs, const PolynomialData& rhs) {
    return polynomial_to_string(
        polynomial_multiply(lhs.coefficients, rhs.coefficients),
        lhs.variable_name);
}

std::string poly_div(const PolynomialData& lhs, const PolynomialData& rhs) {
    const PolynomialDivisionResult division =
        polynomial_divide(lhs.coefficients, rhs.coefficients);
    return "quotient: " +
           polynomial_to_string(division.quotient, lhs.variable_name) +
           ", remainder: " +
           polynomial_to_string(division.remainder, lhs.variable_name);
}

std::string roots(const PolynomialData& poly) {
    const std::vector<double> roots = polynomial_real_roots(poly.coefficients);
    if (roots.empty()) {
        return "No real roots.";
    }

    std::ostringstream out;
    for (std::size_t i = 0; i < roots.size(); ++i) {
        if (i != 0) {
            out << ", ";
        }
        const double root =
            is_integer_double(roots[i], 1e-6)
                ? static_cast<double>(round_to_long_long(roots[i]))
                : roots[i];
        out << format_symbolic_scalar(root);
    }
    return out.str();
}

bool is_polynomial_command(const std::string& command) {
    return command == "poly_add" ||
           command == "poly_sub" ||
           command == "poly_mul" ||
           command == "poly_div" ||
           command == "roots";
}

bool handle_polynomial_command(const PolynomialContext& ctx,
                               const std::string& command,
                               const std::string& inside,
                               std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "roots") {
        if (arguments.size() != 1) {
            throw std::runtime_error("roots expects exactly one argument");
        }
        PolynomialData poly = build_polynomial(ctx, arguments[0]);
        *output = roots(poly);
        return true;
    }

    // poly_add, poly_sub, poly_mul, poly_div
    if (arguments.size() != 2) {
        throw std::runtime_error("polynomial operations expect exactly two arguments");
    }

    PolynomialData lhs = build_polynomial(ctx, arguments[0]);
    PolynomialData rhs = build_polynomial(ctx, arguments[1]);

    if (command == "poly_add") {
        *output = poly_add(lhs, rhs);
    } else if (command == "poly_sub") {
        *output = poly_sub(lhs, rhs);
    } else if (command == "poly_mul") {
        *output = poly_mul(lhs, rhs);
    } else if (command == "poly_div") {
        *output = poly_div(lhs, rhs);
    } else {
        return false;
    }
    return true;
}

}  // namespace polynomial_ops

// ============================================================================
// 级数展开命令实现
// ============================================================================

#include "calculator_series.h"

#include "polynomial.h"
#include "mymath.h"

#include <sstream>

namespace series_ops {

namespace {

std::string simplify_symbolic_text(const std::string& text) {
    return SymbolicExpression::parse(text).simplify().to_string();
}

}  // namespace

std::string taylor(const SeriesContext& ctx,
                   const std::string& expr,
                   double center,
                   int degree) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

    const std::vector<double> coefficients =
        ctx.build_taylor_coefficients(expression, variable_name, center, degree);
    return taylor_series_to_string(coefficients, variable_name, center);
}

std::string pade(const SeriesContext& ctx,
                 const std::string& expr,
                 double center,
                 int numerator_degree,
                 int denominator_degree) {
    if (numerator_degree == 0 && denominator_degree == 0) {
        throw std::runtime_error("pade requires at least one non-zero degree");
    }

    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

    const std::vector<double> coefficients = ctx.build_taylor_coefficients(
        expression, variable_name, center, numerator_degree + denominator_degree);

    auto coefficient_at = [&](int index) {
        if (index < 0 || index >= static_cast<int>(coefficients.size())) {
            return 0.0;
        }
        return coefficients[static_cast<std::size_t>(index)];
    };

    std::vector<double> denominator(denominator_degree + 1, 0.0);
    denominator[0] = 1.0;
    if (denominator_degree > 0) {
        std::vector<std::vector<double>> matrix(
            static_cast<std::size_t>(denominator_degree),
            std::vector<double>(static_cast<std::size_t>(denominator_degree), 0.0));
        std::vector<double> rhs(static_cast<std::size_t>(denominator_degree), 0.0);
        for (int row = 0; row < denominator_degree; ++row) {
            for (int col = 0; col < denominator_degree; ++col) {
                matrix[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(col)] =
                    coefficient_at(numerator_degree + row - col);
            }
            rhs[static_cast<std::size_t>(row)] =
                -coefficient_at(numerator_degree + row + 1);
        }
        const std::vector<double> solved = solve_dense_linear_system(
            matrix, rhs, "pade");
        for (int i = 0; i < denominator_degree; ++i) {
            denominator[static_cast<std::size_t>(i + 1)] =
                solved[static_cast<std::size_t>(i)];
        }
    }

    std::vector<double> numerator(numerator_degree + 1, 0.0);
    for (int i = 0; i <= numerator_degree; ++i) {
        double value = 0.0;
        for (int j = 0; j <= denominator_degree && j <= i; ++j) {
            value += denominator[static_cast<std::size_t>(j)] *
                     coefficient_at(i - j);
        }
        numerator[static_cast<std::size_t>(i)] = value;
    }

    const std::string base = shifted_series_base(variable_name, center);
    const std::string numerator_text = polynomial_to_string(numerator, base);
    const std::string denominator_text = polynomial_to_string(denominator, base);
    if (denominator_text == "1") {
        return simplify_symbolic_text(numerator_text);
    } else {
        return simplify_symbolic_text(
            "(" + numerator_text + ") / (" + denominator_text + ")");
    }
}

std::string puiseux(const SeriesContext& ctx,
                    const std::string& expr,
                    double center,
                    int degree,
                    int denominator) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

    const std::string auxiliary_variable = "puiseux_t";
    const std::string replacement_text =
        mymath::is_near_zero(center, 1e-10)
            ? auxiliary_variable + " ^ " + std::to_string(denominator)
            : format_symbolic_scalar(center) + " + " +
                  auxiliary_variable + " ^ " +
                  std::to_string(denominator);
    const SymbolicExpression substituted = expression.substitute(
        variable_name, SymbolicExpression::parse(replacement_text));
    const std::vector<double> coefficients = ctx.build_taylor_coefficients(
        substituted, auxiliary_variable, 0.0, degree);
    return generalized_series_to_string(
        coefficients, variable_name, center, denominator);
}

std::string series_sum(const SeriesContext& ctx,
                       const std::string& expr,
                       const std::string& index_name,
                       const std::string& lower,
                       const std::string& upper) {
    // 这个函数的实现比较复杂，涉及符号求和和数值求和
    // 暂时保留在 calculator_commands.cpp 中，后续可以迁移
    throw std::runtime_error("series_sum not yet implemented in series_ops module");
}

bool is_series_command(const std::string& command) {
    return command == "taylor" ||
           command == "pade" ||
           command == "puiseux" ||
           command == "series_sum" ||
           command == "summation";
}

bool handle_series_command(const SeriesContext& ctx,
                           const std::string& command,
                           const std::string& inside,
                           std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "taylor") {
        if (arguments.size() != 3) {
            throw std::runtime_error("taylor expects exactly three arguments");
        }

        const double center = ctx.parse_decimal(arguments[1]);
        const double degree_value = ctx.parse_decimal(arguments[2]);
        if (!is_integer_double(degree_value) || degree_value < 0.0) {
            throw std::runtime_error("taylor degree must be a non-negative integer");
        }
        const int degree = static_cast<int>(round_to_long_long(degree_value));

        *output = taylor(ctx, arguments[0], center, degree);
        return true;
    }

    if (command == "pade") {
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "pade expects expr, m, n or expr, center, m, n");
        }

        const bool explicit_center = arguments.size() == 4;
        const double center = explicit_center
                                  ? ctx.parse_decimal(arguments[1])
                                  : 0.0;
        const double numerator_degree_value = ctx.parse_decimal(
            arguments[explicit_center ? 2 : 1]);
        const double denominator_degree_value = ctx.parse_decimal(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(numerator_degree_value) ||
            numerator_degree_value < 0.0 ||
            !is_integer_double(denominator_degree_value) ||
            denominator_degree_value < 0.0) {
            throw std::runtime_error(
                "pade degrees must be non-negative integers");
        }

        const int numerator_degree =
            static_cast<int>(round_to_long_long(numerator_degree_value));
        const int denominator_degree =
            static_cast<int>(round_to_long_long(denominator_degree_value));

        *output = pade(ctx, arguments[0], center, numerator_degree, denominator_degree);
        return true;
    }

    if (command == "puiseux") {
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "puiseux expects expr, degree, denominator or expr, center, degree, denominator");
        }

        const bool explicit_center = arguments.size() == 4;
        const double center = explicit_center
                                  ? ctx.parse_decimal(arguments[1])
                                  : 0.0;
        const double degree_value = ctx.parse_decimal(
            arguments[explicit_center ? 2 : 1]);
        const double denominator_value = ctx.parse_decimal(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(degree_value) || degree_value < 0.0) {
            throw std::runtime_error(
                "puiseux degree must be a non-negative integer");
        }
        if (!is_integer_double(denominator_value) || denominator_value <= 0.0) {
            throw std::runtime_error(
                "puiseux denominator must be a positive integer");
        }

        const int degree = static_cast<int>(round_to_long_long(degree_value));
        const int denom = static_cast<int>(round_to_long_long(denominator_value));

        *output = puiseux(ctx, arguments[0], center, degree, denom);
        return true;
    }

    // series_sum / summation 暂时返回 false，由 calculator_commands.cpp 处理
    return false;
}

}  // namespace series_ops

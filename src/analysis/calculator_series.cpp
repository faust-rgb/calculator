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

std::vector<double> build_taylor_coefficients(
    const SeriesContext& ctx,
    const SymbolicExpression& expression,
    const std::string& variable_name,
    double center,
    int degree) {
    struct TaylorDerivativeCacheEntry {
        SymbolicExpression derivative;
        double value = 0.0;
        bool has_value = false;
    };
    static thread_local std::map<std::string, TaylorDerivativeCacheEntry> derivative_cache;
    static constexpr std::size_t kMaxTaylorDerivativeCacheSize = 256;

    const std::string base_key =
        variable_name + "|" + format_symbolic_scalar(center) + "|" +
        expression.simplify().to_string();
    std::vector<double> coefficients;
    coefficients.reserve(static_cast<std::size_t>(degree + 1));
    SymbolicExpression current = expression;
    for (int order = 0; order <= degree; ++order) {
        const std::string order_key = base_key + "|" + std::to_string(order);
        auto found = derivative_cache.find(order_key);
        if (found == derivative_cache.end()) {
            if (derivative_cache.size() >= kMaxTaylorDerivativeCacheSize) {
                derivative_cache.clear();
            }
            TaylorDerivativeCacheEntry entry;
            entry.derivative = current.simplify();
            found = derivative_cache.emplace(order_key, entry).first;
        } else {
            current = found->second.derivative;
        }

        if (!found->second.has_value) {
            found->second.value =
                ctx.evaluate_at(found->second.derivative, variable_name, center);
            found->second.has_value = true;
        }
        const double derivative_value = found->second.value;
        coefficients.push_back(derivative_value / factorial_value(order));
        if (order != degree) {
            const std::string next_key = base_key + "|" + std::to_string(order + 1);
            auto next_found = derivative_cache.find(next_key);
            if (next_found != derivative_cache.end()) {
                current = next_found->second.derivative;
            } else {
                current = found->second.derivative.derivative(variable_name).simplify();
            }
        }
    }
    return coefficients;
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
        build_taylor_coefficients(ctx, expression, variable_name, center, degree);
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

    const std::vector<double> coefficients = build_taylor_coefficients(
        ctx, expression, variable_name, center, numerator_degree + denominator_degree);

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
    const std::vector<double> coefficients = build_taylor_coefficients(
        ctx, substituted, auxiliary_variable, 0.0, degree);
    return generalized_series_to_string(
        coefficients, variable_name, center, denominator);
}

std::string series_sum(const SeriesContext& ctx,
                       const std::string& expr,
                       const std::string& index_name,
                       const std::string& lower,
                       const std::string& upper) {
    SymbolicExpression summand = SymbolicExpression::parse(ctx.expand_inline(expr));
    SymbolicExpression upper_expression;
    const bool upper_is_infinite =
        upper == "inf" || upper == "oo" || upper == "infinity";
    if (!upper_is_infinite) {
        upper_expression = SymbolicExpression::parse(ctx.expand_inline(upper));
    }

    auto make_polynomial_sum_primitive =
        [&](const std::vector<double>& coefficients) {
            if (coefficients.size() > 4) {
                throw std::runtime_error(
                    "series_sum polynomial summands are currently supported up to degree 3");
            }

            std::vector<std::string> pieces;
            if (coefficients.size() >= 1 &&
                !mymath::is_near_zero(coefficients[0], 1e-10)) {
                pieces.push_back("(" + format_symbolic_scalar(coefficients[0]) +
                                 ") * (" + index_name + " + 1)");
            }
            if (coefficients.size() >= 2 &&
                !mymath::is_near_zero(coefficients[1], 1e-10)) {
                pieces.push_back("(" + format_symbolic_scalar(coefficients[1]) +
                                 ") * (" + index_name + " * (" + index_name +
                                 " + 1) / 2)");
            }
            if (coefficients.size() >= 3 &&
                !mymath::is_near_zero(coefficients[2], 1e-10)) {
                pieces.push_back("(" + format_symbolic_scalar(coefficients[2]) +
                                 ") * (" + index_name + " * (" + index_name +
                                 " + 1) * (2 * " + index_name + " + 1) / 6)");
            }
            if (coefficients.size() >= 4 &&
                !mymath::is_near_zero(coefficients[3], 1e-10)) {
                pieces.push_back("(" + format_symbolic_scalar(coefficients[3]) +
                                 ") * ((" + index_name + " * (" + index_name +
                                 " + 1) / 2) ^ 2)");
            }

            if (pieces.empty()) {
                return SymbolicExpression::number(0.0);
            }
            std::ostringstream out;
            for (std::size_t i = 0; i < pieces.size(); ++i) {
                if (i != 0) {
                    out << " + ";
                }
                out << pieces[i];
            }
            return SymbolicExpression::parse(out.str()).simplify();
        };

    auto finite_sum_from_primitive =
        [&](const SymbolicExpression& primitive) {
            const SymbolicExpression lower_minus_one =
                SymbolicExpression::parse("(" + lower + ") - 1").simplify();
            return SymbolicExpression::parse(
                       "(" +
                       primitive.substitute(index_name, upper_expression).to_string() +
                       ") - (" +
                       primitive.substitute(index_name, lower_minus_one).to_string() +
                       ")")
                .simplify()
                .to_string();
        };

    std::vector<double> polynomial_coefficients;
    if (summand.polynomial_coefficients(index_name, &polynomial_coefficients)) {
        if (upper_is_infinite) {
            bool all_zero = true;
            for (double coefficient : polynomial_coefficients) {
                if (!mymath::is_near_zero(coefficient, 1e-10)) {
                    all_zero = false;
                    break;
                }
            }
            if (!all_zero) {
                throw std::runtime_error(
                    "series_sum does not support infinite polynomial sums");
            }
            return "0";
        }

        const SymbolicExpression primitive =
            make_polynomial_sum_primitive(polynomial_coefficients);
        return finite_sum_from_primitive(primitive);
    }

    auto geometric_ratio = [&](double* coefficient, double* ratio) {
        const double s0 = ctx.evaluate_at(summand, index_name, 0.0);
        const double s1 = ctx.evaluate_at(summand, index_name, 1.0);
        const double s2 = ctx.evaluate_at(summand, index_name, 2.0);
        const double s3 = ctx.evaluate_at(summand, index_name, 3.0);
        if (mymath::is_near_zero(s0, 1e-10)) {
            return false;
        }
        const double candidate = s1 / s0;
        if (!mymath::is_near_zero(s2 - s1 * candidate, 1e-8) ||
            !mymath::is_near_zero(s3 - s2 * candidate, 1e-8)) {
            return false;
        }
        *coefficient = s0;
        *ratio = candidate;
        return true;
    };

    double geometric_coefficient = 0.0;
    double geometric_ratio_value = 0.0;
    if (!geometric_ratio(&geometric_coefficient, &geometric_ratio_value)) {
        throw std::runtime_error(
            "series_sum currently supports polynomial summands up to degree 3 and common geometric series");
    }

    const std::string coefficient_text =
        format_symbolic_scalar(geometric_coefficient);
    const std::string ratio_text = format_symbolic_scalar(geometric_ratio_value);

    if (upper_is_infinite) {
        if (mymath::abs(geometric_ratio_value) >= 1.0 - 1e-10) {
            throw std::runtime_error(
                "series_sum infinite geometric series requires |r| < 1");
        }
        if (mymath::is_near_zero(geometric_ratio_value - 1.0, 1e-10)) {
            throw std::runtime_error(
                "series_sum infinite geometric series diverges for r = 1");
        }
        return ctx.simplify_symbolic(
            "(" + coefficient_text + ") * (" + ratio_text + ") ^ (" +
            lower + ") / (1 - (" + ratio_text + "))");
    }

    const std::string geometric_primitive_text =
        mymath::is_near_zero(geometric_ratio_value - 1.0, 1e-10)
            ? "(" + coefficient_text + ") * (" + index_name + " + 1)"
            : "(" + coefficient_text + ") * (1 - (" + ratio_text +
                  ") ^ (" + index_name + " + 1)) / (1 - (" +
                  ratio_text + "))";
    const SymbolicExpression primitive =
        SymbolicExpression::parse(geometric_primitive_text).simplify();
    return finite_sum_from_primitive(primitive);
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

    if (command == "series_sum" || command == "summation") {
        if (arguments.size() != 4) {
            throw std::runtime_error(
                "series_sum expects expr, index, lower, upper");
        }

        const std::string index_name = trim_copy(arguments[1]);
        if (!is_identifier_text(index_name)) {
            throw std::runtime_error("series_sum index must be an identifier");
        }

        *output = series_sum(ctx,
                             arguments[0],
                             index_name,
                             arguments[2],
                             trim_copy(arguments[3]));
        return true;
    }

    return false;
}

}  // namespace series_ops

#include "analysis/series/puiseux_series.h"
#include "analysis/series/taylor_series.h"
#include "analysis/modules/series_module.h"
#include "core/scalar_type.h"
#include "core/format_utils.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "math/mymath.h"
#include "math/helpers/integer_helpers.h"
#include "statistics/probability.h"
#include "expression_utils.h"
#include "string_utils.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace series_ops {
namespace puiseux {

bool compute_newton_polygon(const std::vector<int>& x_powers,
                            const std::vector<int>& y_powers,
                            const std::vector<Scalar>& coefficients,
                            std::vector<std::tuple<Scalar, int, int, int, int>>* edges) {
    if (x_powers.empty() || y_powers.empty()) return false;

    std::vector<std::pair<int, int>> points;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        if (mymath::precise128::abs(coefficients[i]) >= Scalar(1e-15L)) {
            points.emplace_back(x_powers[i], y_powers[i]);
        }
    }

    if (points.empty()) return false;

    std::sort(points.begin(), points.end());

    int min_y = points[0].second;
    int min_x_at_min_y = points[0].first;
    for (const auto& p : points) {
        if (p.second < min_y || (p.second == min_y && p.first < min_x_at_min_y)) {
            min_y = p.second;
            min_x_at_min_y = p.first;
        }
    }

    int max_y = 0;
    for (const auto& p : points) {
        max_y = std::max(max_y, p.second);
    }

    if (max_y <= 0) return false;

    int min_x = mymath::kIntMax;
    for (const auto& p : points) {
        if (p.second == 0) {
            min_x = std::min(min_x, p.first);
        }
    }

    if (min_x == mymath::kIntMax) {
        Scalar min_ratio = mymath::kDoubleMax;
        for (const auto& p : points) {
            if (p.second > 0) {
                Scalar ratio = Scalar(p.first) / Scalar(p.second);
                min_ratio = std::min(min_ratio, ratio);
            }
        }
        if (min_ratio != mymath::kDoubleMax) {
            edges->emplace_back(min_ratio, 0, 0, 0, max_y);
            return true;
        }
        return false;
    }

    Scalar slope = Scalar(min_x) / Scalar(max_y);
    edges->emplace_back(slope, 0, min_x, 0, max_y);
    return true;
}

bool newton_puiseux_expand(const std::vector<Scalar>& poly_coeffs,
                           int degree,
                           std::vector<std::pair<Scalar, int>>* result) {
    int leading_power = -1;
    Scalar leading_coeff = Scalar(0);
    for (std::size_t i = 0; i < poly_coeffs.size(); ++i) {
        if (mymath::precise128::abs(poly_coeffs[i]) >= Scalar(1e-15L)) {
            leading_power = static_cast<int>(i);
            leading_coeff = poly_coeffs[i];
            break;
        }
    }

    if (leading_power < 0) {
        return false;
    }

    result->clear();
    result->reserve(degree + 1);

    Scalar sqrt_leading = mymath::precise128::sqrt(mymath::precise128::abs(leading_coeff));
    result->emplace_back(sqrt_leading, leading_power);

    if (leading_power + 1 >= static_cast<int>(poly_coeffs.size())) {
        return true;
    }

    std::vector<Scalar> normalized;
    normalized.reserve(poly_coeffs.size() - leading_power);
    for (std::size_t i = leading_power; i < poly_coeffs.size(); ++i) {
        normalized.push_back(poly_coeffs[i] / leading_coeff);
    }

    std::vector<Scalar> sqrt_coeffs;
    sqrt_coeffs.reserve(degree + 1);
    sqrt_coeffs.push_back(Scalar(1));

    for (int k = 1; k <= degree; ++k) {
        Scalar coeff = Scalar(1);
        for (int j = 1; j < k; ++j) {
            coeff = coeff * (Scalar(0.5L) - Scalar(j));
        }
        coeff = coeff / Scalar(prob::factorial(k));
        sqrt_coeffs.push_back(coeff);
    }

    std::vector<Scalar> t_series(degree + 1, Scalar(0));
    for (std::size_t i = 1; i < normalized.size() && i <= static_cast<std::size_t>(degree); ++i) {
        t_series[i] = normalized[i];
    }

    std::vector<Scalar> sqrt_result(degree + 1, Scalar(0));
    sqrt_result[0] = Scalar(1);

    std::vector<Scalar> t_power(degree + 1, Scalar(0));
    t_power[0] = Scalar(1);

    for (int k = 1; k <= degree; ++k) {
        std::vector<Scalar> new_t_power(degree + 1, Scalar(0));
        for (int i = 0; i <= degree; ++i) {
            if (mymath::precise128::abs(t_power[i]) < Scalar(1e-15L)) continue;
            for (int j = 0; i + j <= degree; ++j) {
                if (mymath::precise128::abs(t_series[j]) >= Scalar(1e-15L)) {
                    new_t_power[i + j] = new_t_power[i + j] + t_power[i] * t_series[j];
                }
            }
        }
        t_power = new_t_power;

        for (int i = 0; i <= degree; ++i) {
            sqrt_result[i] = sqrt_result[i] + sqrt_coeffs[k] * t_power[i];
        }
    }

    for (int i = 0; i <= degree && i < static_cast<int>(sqrt_result.size()); ++i) {
        if (mymath::precise128::abs(sqrt_result[i]) >= Scalar(1e-15L)) {
            result->emplace_back(sqrt_leading * sqrt_result[i], leading_power + i);
        }
    }

    return true;
}

std::string puiseux(const SeriesContext& ctx,
                    const std::string& expr,
                    Scalar center,
                    int degree,
                    int denominator) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

    bool is_sqrt_form = false;
    SymbolicExpression sqrt_arg;
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == "sqrt") {
        sqrt_arg = SymbolicExpression(expression.node_->left);
        is_sqrt_form = true;
    }

    if (is_sqrt_form && mymath::is_near_zero(center, 1e-10)) {
        SymbolicExpression arg_at_zero = sqrt_arg.substitute(variable_name, SymbolicExpression::number(0.0L)).simplify();
        Scalar val = 0.0L;
        if (arg_at_zero.is_number(&val) && mymath::is_near_zero(val, 1e-10)) {
            std::vector<Scalar> poly_coeffs;
            std::vector<Scalar> ld_poly_coeffs;
            if (sqrt_arg.polynomial_coefficients(variable_name, &ld_poly_coeffs)) {
                for (const auto& c : ld_poly_coeffs) {
                    poly_coeffs.push_back(Scalar(c));
                }
                int leading_power = -1;
                Scalar leading_coeff = Scalar(0);
                for (std::size_t i = 0; i < poly_coeffs.size(); ++i) {
                    if (mymath::precise128::abs(poly_coeffs[i]) >= Scalar(1e-10L)) {
                        leading_power = static_cast<int>(i);
                        leading_coeff = poly_coeffs[i];
                        break;
                    }
                }

                if (leading_power >= 1) {
                    std::vector<std::pair<Scalar, int>> puiseux_coeffs;
                    if (newton_puiseux_expand(poly_coeffs, degree, &puiseux_coeffs)) {
                        std::ostringstream result;
                        bool first = true;
                        for (const auto& [coeff, power] : puiseux_coeffs) {
                            if (mymath::precise128::abs(coeff) < Scalar(1e-15L)) continue;

                            Scalar coeff_ld = coeff.to_long_double();
                            if (!first) {
                                if (coeff_ld > 0) {
                                    result << " + ";
                                } else {
                                    result << " - ";
                                }
                            } else {
                                if (coeff_ld < 0) {
                                    result << "-";
                                }
                            }

                            Scalar abs_coeff = mymath::abs(coeff_ld);
                            int actual_power = power;
                            int actual_denom = 2;

                            if (actual_power % 2 == 0) {
                                actual_power /= 2;
                                actual_denom = 1;
                            }

                            const bool has_power = actual_power > 0;
                            const bool omit_unit_coeff =
                                has_power &&
                                mymath::is_near_zero(abs_coeff - 1.0L, 1e-9);

                            if (!omit_unit_coeff) {
                                if (!first && coeff_ld < 0) {
                                    result << format_symbolic_scalar(abs_coeff);
                                } else {
                                    result << format_symbolic_scalar(coeff_ld);
                                }
                            }

                            if (actual_denom == 1) {
                                if (actual_power == 1) {
                                    if (!omit_unit_coeff) {
                                        result << " * ";
                                    }
                                    result << variable_name;
                                } else if (actual_power > 0) {
                                    if (!omit_unit_coeff) {
                                        result << " * ";
                                    }
                                    result << variable_name << " ^ " << actual_power;
                                }
                            } else {
                                if (!omit_unit_coeff) {
                                    result << " * ";
                                }
                                result << variable_name << " ^ (" << actual_power
                                       << " / " << actual_denom << ")";
                            }

                            first = false;
                        }

                        std::string result_str = result.str();
                        if (result_str.empty()) {
                            return "0";
                        }
                        return result_str;
                    }

                    if (leading_power % 2 == 0) {
                        const std::string auxiliary_variable = "puiseux_t";
                        const std::string replacement_text =
                            auxiliary_variable + " ^ " + std::to_string(denominator);
                        const SymbolicExpression substituted = expression.substitute(
                            variable_name, SymbolicExpression::parse(replacement_text));
                        std::string substituted_text = substituted.to_string();
                        const std::string positive_aux_abs = "abs(" + auxiliary_variable + ")";
                        std::size_t abs_pos = 0;
                        while ((abs_pos = substituted_text.find(positive_aux_abs, abs_pos)) !=
                               std::string::npos) {
                            substituted_text.replace(abs_pos,
                                                    positive_aux_abs.size(),
                                                    auxiliary_variable);
                            abs_pos += auxiliary_variable.size();
                        }
                        const SymbolicExpression puiseux_expression =
                            SymbolicExpression::parse(substituted_text);
                        const std::vector<Scalar> coefficients = taylor::build_taylor_coefficients(
                            ctx, puiseux_expression, auxiliary_variable, 0.0L, degree);
                        std::vector<Scalar> ld_coefficients;
                        for (const auto& c : coefficients) {
                            ld_coefficients.push_back(c.to_long_double());
                        }
                        return generalized_series_to_string(
                            ld_coefficients, variable_name, center, denominator);
                    } else {
                        int effective_denominator = denominator;
                        if (denominator == 1) {
                            effective_denominator = 2;
                        }

                        const std::string auxiliary_variable = "puiseux_t";
                        const std::string replacement_text =
                            auxiliary_variable + " ^ " + std::to_string(effective_denominator);
                        const SymbolicExpression substituted = expression.substitute(
                            variable_name, SymbolicExpression::parse(replacement_text));
                        std::string substituted_text = substituted.to_string();
                        const std::string positive_aux_abs = "abs(" + auxiliary_variable + ")";
                        std::size_t abs_pos = 0;
                        while ((abs_pos = substituted_text.find(positive_aux_abs, abs_pos)) !=
                               std::string::npos) {
                            substituted_text.replace(abs_pos,
                                                    positive_aux_abs.size(),
                                                    auxiliary_variable);
                            abs_pos += auxiliary_variable.size();
                        }
                        const SymbolicExpression puiseux_expression =
                            SymbolicExpression::parse(substituted_text);
                        const std::vector<Scalar> raw_coeffs = taylor::build_taylor_coefficients(
                            ctx, puiseux_expression, auxiliary_variable, 0.0L, degree);

                        std::vector<Scalar> ld_raw_coeffs;
                        for (const auto& c : raw_coeffs) {
                            ld_raw_coeffs.push_back(c.to_long_double());
                        }
                        return generalized_series_to_string(
                            ld_raw_coeffs, variable_name, center, effective_denominator);
                    }
                }
            }
        }
    }

    const std::string auxiliary_variable = "puiseux_t";
    const std::string replacement_text =
        mymath::is_near_zero(center, 1e-10)
            ? auxiliary_variable + " ^ " + std::to_string(denominator)
            : format_symbolic_scalar(center) + " + " +
                  auxiliary_variable + " ^ " +
                  std::to_string(denominator);
    const SymbolicExpression substituted = expression.substitute(
        variable_name, SymbolicExpression::parse(replacement_text));
    std::string substituted_text = substituted.to_string();
    const std::string positive_aux_abs = "abs(" + auxiliary_variable + ")";
    std::size_t abs_pos = 0;
    while ((abs_pos = substituted_text.find(positive_aux_abs, abs_pos)) !=
           std::string::npos) {
        substituted_text.replace(abs_pos,
                                positive_aux_abs.size(),
                                auxiliary_variable);
        abs_pos += auxiliary_variable.size();
    }
    const SymbolicExpression puiseux_expression =
        SymbolicExpression::parse(substituted_text);
    const std::vector<Scalar> coefficients = taylor::build_taylor_coefficients(
        ctx, puiseux_expression, auxiliary_variable, 0.0L, degree);
    std::vector<Scalar> ld_coefficients;
    for (const auto& c : coefficients) {
        ld_coefficients.push_back(c.to_long_double());
    }
    return generalized_series_to_string(
        ld_coefficients, variable_name, center, denominator);
}

std::string general_series_auto(const SeriesContext& ctx,
                              const std::string& expression_str,
                              Scalar center,
                              int degree) {
    std::string variable_name;
    SymbolicExpression expr;
    ctx.resolve_symbolic(expression_str, true, &variable_name, &expr);

    try {
        SymbolicExpression val = expr.substitute(variable_name, SymbolicExpression::number(center)).simplify();
        Scalar n;
        if (val.is_number(&n) && !mymath::isinf(n)) {
            return taylor::taylor(ctx, expression_str, center, degree);
        }
    } catch (...) {}

    for (int k = 1; k <= 3; ++k) {
        SymbolicExpression factor =
            (SymbolicExpression::variable(variable_name) - SymbolicExpression::number(center)) ^
            SymbolicExpression::number(Scalar(k));
        SymbolicExpression shifted = (expr * factor).simplify();

        try {
            SymbolicExpression val = shifted.substitute(variable_name, SymbolicExpression::number(center)).simplify();
            Scalar n;
            if (val.is_number(&n) && !mymath::isinf(n)) {
                std::string t_str = taylor::taylor(ctx, shifted.to_string(), center, degree + k);
                return "(" + t_str + ") / (" + variable_name + (center == 0 ? "" : " - " + std::to_string(center.to_long_double())) + ")^" + std::to_string(k);
            }
        } catch (...) {}
    }

    return puiseux(ctx, expression_str, center, degree, 2);
}

}  // namespace puiseux
}  // namespace series_ops
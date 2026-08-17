// ============================================================================
// 极限求解器实现（直接基于 Scalar 类型）
// ============================================================================

#include "analysis/calculus/limit_solver.h"
#include "analysis/series/psa_engine.h"
#include "analysis/modules/series_module.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/public/symbolic_node_types.h"
#include "math/numeric/precision/tolerances.h"
#include "math/mymath.h"
#include <algorithm>
#include <vector>
#include <stdexcept>

namespace numeric {

bool is_infinite_probe(SymbolicLimitProbeKind kind) {
    return kind == SymbolicLimitProbeKind::kPositiveInfinity ||
           kind == SymbolicLimitProbeKind::kNegativeInfinity;
}

bool is_zero_probe(SymbolicLimitProbeKind kind, Scalar value) {
    return kind == SymbolicLimitProbeKind::kFinite &&
           mymath::is_near_zero(value, precision::default_relative_tolerance<Scalar>());
}

static bool is_effective_infinity_point(Scalar val) {
    return !mymath::isfinite(val);
}

Scalar handle_pole_limit(int shift, Scalar leading_coefficient, int direction) {
    if (direction == 0) {
        if (shift % 2 == 0) {
            return (leading_coefficient > Scalar(0)) ? mymath::infinity() : -mymath::infinity();
        } else {
            throw std::runtime_error("two-sided limit does not exist (pole with odd shift)");
        }
    } else if (direction == 1) {
        return (leading_coefficient > Scalar(0)) ? mymath::infinity() : -mymath::infinity();
    } else {
        bool flip_sign = (shift % 2 != 0);
        Scalar effective_c = flip_sign ? -leading_coefficient : leading_coefficient;
        return (effective_c > Scalar(0)) ? mymath::infinity() : -mymath::infinity();
    }
}

SymbolicLimitProbeKind probe_symbolic_value_at(
    SymbolicExpression expression,
    const std::string& variable_name,
    Scalar point,
    Scalar* finite_value) {
    try {
        if (is_effective_infinity_point(point)) {
            expression = expression.simplify();
            expression = expression.simplify();

            switch (expression.node_type()) {
                case NodeType::kNumber:
                case NodeType::kPi:
                case NodeType::kE:
                case NodeType::kInfinity: {
                    Scalar value = 0.0L;
                    if (expression.is_number(&value)) {
                        *finite_value = value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kVariable:
                    if (expression.node_text() == variable_name) {
                        return point > Scalar(0) ? SymbolicLimitProbeKind::kPositiveInfinity
                                                 : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                case NodeType::kNegate: {
                    Scalar child_value = Scalar(0);
                    const SymbolicLimitProbeKind child_kind =
                        probe_symbolic_value_at(expression.left_child(),
                                                variable_name,
                                                point,
                                                &child_value);
                    if (child_kind == SymbolicLimitProbeKind::kFinite) {
                        *finite_value = -child_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (child_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                        return SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    if (child_kind == SymbolicLimitProbeKind::kNegativeInfinity) {
                        return SymbolicLimitProbeKind::kPositiveInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kAdd:
                case NodeType::kSubtract: {
                    Scalar left_value = Scalar(0);
                    Scalar right_value = Scalar(0);
                    SymbolicLimitProbeKind left_kind =
                        probe_symbolic_value_at(expression.left_child(),
                                                variable_name,
                                                point,
                                                &left_value);
                    SymbolicLimitProbeKind right_kind =
                        probe_symbolic_value_at(expression.right_child(),
                                                variable_name,
                                                point,
                                                &right_value);
                    if (expression.node_type() == NodeType::kSubtract) {
                        if (right_kind == SymbolicLimitProbeKind::kFinite) {
                            right_value = -right_value;
                        } else if (right_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            right_kind = SymbolicLimitProbeKind::kNegativeInfinity;
                        } else if (right_kind == SymbolicLimitProbeKind::kNegativeInfinity) {
                            right_kind = SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                    }
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        right_kind == SymbolicLimitProbeKind::kFinite) {
                        *finite_value = left_value + right_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (left_kind == SymbolicLimitProbeKind::kFinite) return right_kind;
                    if (right_kind == SymbolicLimitProbeKind::kFinite) return left_kind;
                    if (left_kind == right_kind) return left_kind;
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kMultiply: {
                    Scalar left_value = Scalar(0);
                    Scalar right_value = Scalar(0);
                    const SymbolicLimitProbeKind left_kind =
                        probe_symbolic_value_at(expression.left_child(),
                                                variable_name,
                                                point,
                                                &left_value);
                    const SymbolicLimitProbeKind right_kind =
                        probe_symbolic_value_at(expression.right_child(),
                                                variable_name,
                                                point,
                                                &right_value);
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        right_kind == SymbolicLimitProbeKind::kFinite) {
                        *finite_value = left_value * right_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        mymath::is_near_zero(left_value, precision::default_relative_tolerance<Scalar>())) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    if (right_kind == SymbolicLimitProbeKind::kFinite &&
                        mymath::is_near_zero(right_value, precision::default_relative_tolerance<Scalar>())) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    auto sign_of = [](SymbolicLimitProbeKind kind, Scalar value) {
                        if (kind == SymbolicLimitProbeKind::kFinite) {
                            return value >= Scalar(0) ? 1 : -1;
                        }
                        return kind == SymbolicLimitProbeKind::kPositiveInfinity ? 1 : -1;
                    };
                    if ((left_kind == SymbolicLimitProbeKind::kFinite || is_infinite_probe(left_kind)) &&
                        (right_kind == SymbolicLimitProbeKind::kFinite || is_infinite_probe(right_kind))) {
                        const int sign = sign_of(left_kind, left_value) *
                                         sign_of(right_kind, right_value);
                        return sign >= 0 ? SymbolicLimitProbeKind::kPositiveInfinity
                                         : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kDivide: {
                    Scalar left_value = Scalar(0);
                    Scalar right_value = Scalar(0);
                    const SymbolicLimitProbeKind left_kind =
                        probe_symbolic_value_at(expression.left_child(),
                                                variable_name,
                                                point,
                                                &left_value);
                    const SymbolicLimitProbeKind right_kind =
                        probe_symbolic_value_at(expression.right_child(),
                                                variable_name,
                                                point,
                                                &right_value);
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        right_kind == SymbolicLimitProbeKind::kFinite &&
                        !mymath::is_near_zero(right_value, precision::default_relative_tolerance<Scalar>())) {
                        *finite_value = left_value / right_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        is_infinite_probe(right_kind)) {
                        *finite_value = Scalar(0);
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (is_infinite_probe(left_kind) &&
                        right_kind == SymbolicLimitProbeKind::kFinite &&
                        !mymath::is_near_zero(right_value, precision::default_relative_tolerance<Scalar>())) {
                        const bool positive =
                            (left_kind == SymbolicLimitProbeKind::kPositiveInfinity) ==
                            (right_value > Scalar(0));
                        return positive ? SymbolicLimitProbeKind::kPositiveInfinity
                                        : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kFunction: {
                    const std::string name = expression.node_text();
                    Scalar argument_value = Scalar(0);
                    const SymbolicLimitProbeKind argument_kind =
                        probe_symbolic_value_at(expression.left_child(),
                                                variable_name,
                                                point,
                                                &argument_value);
                    if (name == "ln") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite &&
                            argument_value > Scalar(0)) {
                            *finite_value = mymath::ln(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    if (name == "exp") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kNegativeInfinity) {
                            *finite_value = Scalar(0);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite) {
                            *finite_value = mymath::exp(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    if (name == "sin" || name == "cos") {
                        if (argument_kind == SymbolicLimitProbeKind::kFinite) {
                            *finite_value = (name == "sin") ? mymath::sin(argument_value)
                                                            : mymath::cos(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    if (name == "sqrt") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite &&
                            argument_value >= Scalar(0)) {
                            *finite_value = mymath::sqrt(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                default:
                    return SymbolicLimitProbeKind::kUnknown;
            }
        }

        SymbolicExpression sub_expr = expression.substitute(
            variable_name,
            SymbolicExpression::number(point)).simplify();
        Scalar value = 0.0L;
        if (!sub_expr.is_number(&value)) {
            return SymbolicLimitProbeKind::kUnknown;
        }
        if (mymath::isfinite(value)) {
            *finite_value = value;
            return SymbolicLimitProbeKind::kFinite;
        }
        return value > Scalar(0) ? SymbolicLimitProbeKind::kPositiveInfinity
                                 : SymbolicLimitProbeKind::kNegativeInfinity;
    } catch (...) {
        return SymbolicLimitProbeKind::kUnknown;
    }
}

bool try_symbolic_one_to_infinity_limit(
    const SymbolicExpression& base,
    const SymbolicExpression& exponent,
    const std::string& variable_name,
    Scalar point,
    Scalar* result) {
    Scalar exponent_value = Scalar(0);
    const SymbolicLimitProbeKind exponent_kind =
        probe_symbolic_value_at(exponent, variable_name, point, &exponent_value);
    if (!is_infinite_probe(exponent_kind)) {
        return false;
    }

    Scalar base_value = Scalar(0);
    const SymbolicLimitProbeKind base_kind =
        probe_symbolic_value_at(base, variable_name, point, &base_value);
    if (base_kind != SymbolicLimitProbeKind::kFinite ||
        !mymath::is_near_zero(base_value - Scalar(1), Scalar(1e-8L))) {
        return false;
    }

    const SymbolicExpression transformed =
        SymbolicExpression::function(
            "exp",
            ((base - SymbolicExpression::number(1.0L)) * exponent).simplify()).simplify();
    Scalar transformed_value = Scalar(0);
    const SymbolicLimitProbeKind transformed_kind =
        probe_symbolic_value_at(transformed, variable_name, point, &transformed_value);
    if (transformed_kind == SymbolicLimitProbeKind::kFinite) {
        *result = transformed_value;
        return true;
    }
    return false;
}

bool try_symbolic_lhopital_limit(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    Scalar point,
    int direction,
    Scalar* result,
    std::function<Scalar(const SymbolicExpression&, const std::string&, Scalar)> evaluate_at_override) {
    SymbolicExpression current = expression.simplify();
    if (current.node_type() != NodeType::kDivide) {
        return false;
    }

    SymbolicExpression numerator = current.left_child();
    SymbolicExpression denominator = current.right_child();
    Scalar numerator_value = Scalar(0);
    Scalar denominator_value = Scalar(0);
    const SymbolicLimitProbeKind numerator_kind =
        probe_symbolic_value_at(numerator, variable_name, point, &numerator_value);
    const SymbolicLimitProbeKind denominator_kind =
        probe_symbolic_value_at(denominator, variable_name, point, &denominator_value);

    const bool zero_over_zero =
        is_zero_probe(numerator_kind, numerator_value) &&
        is_zero_probe(denominator_kind, denominator_value);
    const bool infinity_over_infinity =
        is_infinite_probe(numerator_kind) &&
        is_infinite_probe(denominator_kind);

    if (!zero_over_zero && !infinity_over_infinity) {
        return false;
    }

    if (mymath::isfinite(point)) {
        series_ops::SeriesContext ctx;
        ctx.evaluate_at = evaluate_at_override ? evaluate_at_override :
            [](const SymbolicExpression& e, const std::string&, Scalar) {
                Scalar val = 0.0L;
                if (e.is_number(&val)) return val;
                return Scalar(0.0L);
            };

        struct LaurentInfo {
            int degree = 0;
            Scalar coefficient = Scalar(0);
            bool valid = false;
        };

        auto get_laurent_info = [&](const SymbolicExpression& e) -> LaurentInfo {
            LaurentInfo info;
            std::vector<Scalar> coeffs;
            try {
                if (series_ops::internal::evaluate_psa(e, variable_name, point, 4, coeffs, ctx)) {
                    for (int i = 0; i < static_cast<int>(coeffs.size()); ++i) {
                        if (!mymath::is_near_zero(coeffs[i], Scalar(1e-15L))) {
                            info.degree = i;
                            info.coefficient = coeffs[i];
                            info.valid = true;
                            return info;
                        }
                    }
                    info.degree = 100;
                    info.coefficient = Scalar(0);
                    info.valid = true;
                    return info;
                }
            } catch (const series_ops::internal::PoleException& ex) {
                info.degree = ex.shift;
                info.coefficient = ex.leading_coefficient;
                info.valid = true;
                return info;
            } catch (...) {
            }
            return info;
        };

        LaurentInfo infoN = get_laurent_info(numerator);
        LaurentInfo infoD = get_laurent_info(denominator);

        if (infoN.valid && infoD.valid) {
            int res_degree = infoN.degree - infoD.degree;
            Scalar res_coeff = infoN.coefficient / infoD.coefficient;

            if (res_degree > 0) {
                *result = Scalar(0);
                return true;
            } else if (res_degree == 0) {
                *result = res_coeff;
                return true;
            } else {
                try {
                    *result = handle_pole_limit(res_degree, res_coeff, direction);
                    return true;
                } catch (...) {
                    return false;
                }
            }
        }
    }

    static constexpr int kMaxLhopitalDepth = 5;
    SymbolicExpression iter_expr = current;
    for (int depth = 0; depth < kMaxLhopitalDepth; ++depth) {
        SymbolicExpression n = iter_expr.left_child();
        SymbolicExpression d = iter_expr.right_child();
        iter_expr = (n.derivative(variable_name).simplify() /
                     d.derivative(variable_name).simplify()).simplify();

        Scalar val = Scalar(0);
        const SymbolicLimitProbeKind kind =
            probe_symbolic_value_at(iter_expr, variable_name, point, &val);
        if (kind == SymbolicLimitProbeKind::kFinite) {
            *result = val;
            return true;
        }
        if (is_infinite_probe(kind)) {
            *result = (kind == SymbolicLimitProbeKind::kPositiveInfinity) ? mymath::infinity() : -mymath::infinity();
            return true;
        }
    }

    return false;
}

bool symbolic_limit_at_infinity(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    bool positive,
    Scalar* result) {
    if (expression.has_node() && expression.node_type() == NodeType::kPower) {
        Scalar transformed = Scalar(0);
        if (try_symbolic_one_to_infinity_limit(expression.left_child(),
                                               expression.right_child(),
                                               variable_name,
                                               positive ? mymath::infinity()
                                                        : -mymath::infinity(),
                                               &transformed)) {
            *result = transformed;
            return true;
        }
    }

    series_ops::SeriesContext ctx;
    ctx.evaluate_at = [](const SymbolicExpression& e, const std::string&, Scalar) {
        Scalar val = 0.0L;
        if (e.is_number(&val)) return val;
        return Scalar(0.0L);
    };

    SymbolicExpression t_var = SymbolicExpression::variable("t_limit_inf_tmp");
    SymbolicExpression inv_t;
    if (positive) {
        inv_t = SymbolicExpression::number(1.0L) / t_var;
    } else {
        inv_t = SymbolicExpression::number(-1.0L) / t_var;
    }
    SymbolicExpression substituted = expression.substitute(variable_name, inv_t).simplify();

    std::vector<Scalar> coeffs;
    try {
        if (series_ops::internal::evaluate_psa(substituted, "t_limit_inf_tmp", Scalar(0), 2, coeffs, ctx)) {
            if (!coeffs.empty() && mymath::isfinite(coeffs[0])) {
                *result = coeffs[0];
                return true;
            }
        }
    } catch (const series_ops::internal::PoleException& e) {
        Scalar lhopital_result = Scalar(0);
        if (try_symbolic_lhopital_limit(expression, variable_name,
                                         positive ? mymath::infinity()
                                                  : -mymath::infinity(),
                                         1,
                                         &lhopital_result)) {
            *result = lhopital_result;
            return true;
        }

        try {
            *result = handle_pole_limit(e.shift, e.leading_coefficient, 1);
            return true;
        } catch (...) {
        }
    }

    Scalar lhopital_result = Scalar(0);
    if (try_symbolic_lhopital_limit(expression, variable_name,
                                     positive ? mymath::infinity()
                                              : -mymath::infinity(),
                                     1,
                                     &lhopital_result)) {
        *result = lhopital_result;
        return true;
    }

    return false;
}

Scalar compute_numerical_limit(
    const std::function<Scalar(Scalar)>& evaluate,
    Scalar x,
    int direction) {
    auto compute_limit_at = [&](Scalar x_target, int side) -> Scalar {
        Scalar richardson[14][14] = {};
        Scalar best_value = Scalar(0);
        Scalar best_error = mymath::infinity();
        bool have_best = false;

        const Scalar base_h = is_effective_infinity_point(x_target)
                             ? Scalar(1e-2L)
                             : precision::sqrt_epsilon<Scalar>() * std::max(Scalar(1), mymath::abs(x_target));
        Scalar adaptive_h = base_h;
        int consecutive_bad = 0;
        constexpr int kMaxBadSamples = 3;

        Scalar prev_val = Scalar(0);
        bool have_prev = false;
        int oscillation_count = 0;
        Scalar total_amplitude = Scalar(0);

        for (int row = 0; row < 14; ++row) {
            const Scalar h = adaptive_h / mymath::pow(Scalar(2), Scalar(row + 4));
            Scalar sample_x;
            if (!is_effective_infinity_point(x_target)) {
                sample_x = x_target + Scalar(side) * h;
            } else {
                sample_x = (x_target > Scalar(0) ? Scalar(1) : Scalar(-1)) / h;
            }

            Scalar val = Scalar(0);
            try {
                val = evaluate(sample_x);
            } catch (...) {
                adaptive_h *= Scalar(0.5L);
                consecutive_bad++;
                if (consecutive_bad >= kMaxBadSamples) {
                    throw std::runtime_error("limit did not converge (sampling failures)");
                }
                continue;
            }

            if (!mymath::isfinite(val)) {
                if (have_prev && mymath::isfinite(prev_val)) {
                    if (prev_val > Scalar(1e10L)) return mymath::infinity();
                    else if (prev_val < Scalar(-1e10L)) return -mymath::infinity();
                }
                adaptive_h *= Scalar(0.5L);
                consecutive_bad++;
                if (consecutive_bad >= kMaxBadSamples) {
                    throw std::runtime_error("limit appears to be infinite (numerical evidence)");
                }
                continue;
            }

            if (have_prev) {
                const Scalar diff = mymath::abs(val - prev_val);
                total_amplitude += diff;
                if ((val > Scalar(0) && prev_val < Scalar(0)) ||
                    (val < Scalar(0) && prev_val > Scalar(0))) {
                    oscillation_count++;
                    if (oscillation_count >= 5) {
                        const Scalar avg_amp = total_amplitude / Scalar(row + 1);
                        if (avg_amp > Scalar(1e-2L)) {
                            throw std::runtime_error("limit does not exist (oscillation)");
                        }
                        adaptive_h *= Scalar(0.25L);
                        oscillation_count = 0;
                    }
                } else {
                    oscillation_count = std::max(0, oscillation_count - 1);
                }
            }
            prev_val = val;
            have_prev = true;

            if (have_best && row > 0) {
                const Scalar expected_change = best_error * Scalar(10) + Scalar(1e-10L);
                const Scalar actual_change = mymath::abs(val - best_value);
                if (actual_change > expected_change * Scalar(1e6L)) {
                    adaptive_h *= Scalar(0.5L);
                    row = -1;
                    consecutive_bad++;
                    if (consecutive_bad >= kMaxBadSamples) break;
                    continue;
                }
            }

            richardson[row][0] = val;
            for (int col = 1; col <= row; ++col) {
                Scalar p4 = mymath::pow(Scalar(2), Scalar(col));
                richardson[row][col] = (p4 * richardson[row][col - 1] - richardson[row - 1][col - 1]) / (p4 - Scalar(1));
            }

            if (row >= 1) {
                const Scalar current_error = mymath::abs(richardson[row][row] - richardson[row - 1][row - 1]);
                if (!have_best || current_error < best_error) {
                    best_value = richardson[row][row];
                    best_error = current_error;
                    have_best = true;
                }
                if (best_error < precision::default_relative_tolerance<Scalar>()) break;
            } else {
                best_value = richardson[0][0];
                have_best = true;
            }
        }

        if (!have_best) throw std::runtime_error("limit did not converge");
        return best_value;
    };

    if (direction == 0) {
        Scalar left = compute_limit_at(x, -1);
        Scalar right = compute_limit_at(x, 1);
        if (mymath::abs(left - right) > precision::default_relative_tolerance<Scalar>() &&
            mymath::isfinite(left) && mymath::isfinite(right)) {
            throw std::runtime_error("limit does not exist (left and right limits differ)");
        }
        return (left + right) * Scalar(0.5L);
    }
    return compute_limit_at(x, direction);
}

}  // namespace numeric

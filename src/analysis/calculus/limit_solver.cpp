// ============================================================================
// 极限求解器实现
// ============================================================================
//
// 本文件实现了符号极限的数值求解算法，包括：
// - 极限点探测（有限值、无穷大、振荡）
// - Puiseux 级数展开辅助
// - 极限值计算与验证
//
// 主要用于处理符号计算中无法直接代入的极限问题，
// 如 0/0 型、∞/∞ 型等不定式。
//
// 相关文件：
// - limit_solver.h: 极限求解器接口定义
// - psa_engine.h: Puiseux 级数展开引擎
// ============================================================================

#include "analysis/calculus/limit_solver.h"
#include "analysis/series/psa_engine.h"
#include "analysis/modules/series_module.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "math/mymath.h"
#include <algorithm>
#include <cmath>

namespace numeric {

template <typename T>
bool is_infinite_probe(SymbolicLimitProbeKind kind) {
    return kind == SymbolicLimitProbeKind::kPositiveInfinity ||
           kind == SymbolicLimitProbeKind::kNegativeInfinity;
}

template <typename T>
bool is_zero_probe(SymbolicLimitProbeKind kind, T value) {
    return kind == SymbolicLimitProbeKind::kFinite &&
           t_is_near_zero(value, precision::default_relative_tolerance<T>());
}

template <typename T>
bool t_is_effective_infinity_point(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return !std::isfinite(val);
    } else if constexpr (std::is_same_v<T, Scalar>) {
        return !mymath::isfinite(val);
    } else {
        return false;
    }
}

template <typename T>
T handle_pole_limit(int shift, T leading_coefficient, int direction) {
    if (direction == 0) {
        if (shift % 2 == 0) {
            return (leading_coefficient > T(static_cast<long long>(0))) ? t_infinity<T>() : -t_infinity<T>();
        } else {
            throw std::runtime_error("two-sided limit does not exist (pole with odd shift)");
        }
    } else if (direction == 1) {
        return (leading_coefficient > T(static_cast<long long>(0))) ? t_infinity<T>() : -t_infinity<T>();
    } else {
        bool flip_sign = (shift % 2 != 0);
        T effective_c = flip_sign ? -leading_coefficient : leading_coefficient;
        return (effective_c > T(static_cast<long long>(0))) ? t_infinity<T>() : -t_infinity<T>();
    }
}

template <typename T>
SymbolicLimitProbeKind probe_symbolic_value_at(
    SymbolicExpression expression,
    const std::string& variable_name,
    T point,
    T* finite_value) {
    try {
        if (t_is_effective_infinity_point(point)) {
            expression = expression.simplify();
            const std::shared_ptr<SymbolicExpression::Node>& node = expression.node_;

            switch (node->type) {
                case NodeType::kNumber:
                case NodeType::kPi:
                case NodeType::kE:
                case NodeType::kInfinity: {
                    Scalar value = 0.0L;
                    if (expression.is_number(&value)) {
                        *finite_value = T(value);
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kVariable:
                    if (node->text == variable_name) {
                        return point > T(static_cast<long long>(0)) ? SymbolicLimitProbeKind::kPositiveInfinity
                                           : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                case NodeType::kNegate: {
                    T child_value = T(static_cast<long long>(0));
                    const SymbolicLimitProbeKind child_kind =
                        probe_symbolic_value_at(SymbolicExpression(node->left),
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
                    T left_value = T(static_cast<long long>(0));
                    T right_value = T(static_cast<long long>(0));
                    SymbolicLimitProbeKind left_kind =
                        probe_symbolic_value_at(SymbolicExpression(node->left),
                                                variable_name,
                                                point,
                                                &left_value);
                    SymbolicLimitProbeKind right_kind =
                        probe_symbolic_value_at(SymbolicExpression(node->right),
                                                variable_name,
                                                point,
                                                &right_value);
                    if (node->type == NodeType::kSubtract) {
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
                    T left_value = T(static_cast<long long>(0));
                    T right_value = T(static_cast<long long>(0));
                    const SymbolicLimitProbeKind left_kind =
                        probe_symbolic_value_at(SymbolicExpression(node->left),
                                                variable_name,
                                                point,
                                                &left_value);
                    const SymbolicLimitProbeKind right_kind =
                        probe_symbolic_value_at(SymbolicExpression(node->right),
                                                variable_name,
                                                point,
                                                &right_value);
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        right_kind == SymbolicLimitProbeKind::kFinite) {
                        *finite_value = left_value * right_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        t_is_near_zero(left_value, precision::default_relative_tolerance<T>())) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    if (right_kind == SymbolicLimitProbeKind::kFinite &&
                        t_is_near_zero(right_value, precision::default_relative_tolerance<T>())) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    auto sign_of = [](SymbolicLimitProbeKind kind, T value) {
                        if (kind == SymbolicLimitProbeKind::kFinite) {
                            return value >= T(static_cast<long long>(0)) ? 1 : -1;
                        }
                        return kind == SymbolicLimitProbeKind::kPositiveInfinity ? 1 : -1;
                    };
                    if ((left_kind == SymbolicLimitProbeKind::kFinite || is_infinite_probe<T>(left_kind)) &&
                        (right_kind == SymbolicLimitProbeKind::kFinite || is_infinite_probe<T>(right_kind))) {
                        const int sign = sign_of(left_kind, left_value) *
                                         sign_of(right_kind, right_value);
                        return sign >= 0 ? SymbolicLimitProbeKind::kPositiveInfinity
                                         : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kDivide: {
                    T left_value = T(static_cast<long long>(0));
                    T right_value = T(static_cast<long long>(0));
                    const SymbolicLimitProbeKind left_kind =
                        probe_symbolic_value_at(SymbolicExpression(node->left),
                                                variable_name,
                                                point,
                                                &left_value);
                    const SymbolicLimitProbeKind right_kind =
                        probe_symbolic_value_at(SymbolicExpression(node->right),
                                                variable_name,
                                                point,
                                                &right_value);
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        right_kind == SymbolicLimitProbeKind::kFinite &&
                        !t_is_near_zero(right_value, precision::default_relative_tolerance<T>())) {
                        *finite_value = left_value / right_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        is_infinite_probe<T>(right_kind)) {
                        *finite_value = T(static_cast<long long>(0));
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (is_infinite_probe<T>(left_kind) &&
                        right_kind == SymbolicLimitProbeKind::kFinite &&
                        !t_is_near_zero(right_value, precision::default_relative_tolerance<T>())) {
                        const bool positive =
                            (left_kind == SymbolicLimitProbeKind::kPositiveInfinity) ==
                            (right_value > T(static_cast<long long>(0)));
                        return positive ? SymbolicLimitProbeKind::kPositiveInfinity
                                        : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kFunction: {
                    const std::string& name = node->text;
                    T argument_value = T(static_cast<long long>(0));
                    const SymbolicLimitProbeKind argument_kind =
                        probe_symbolic_value_at(SymbolicExpression(node->left),
                                                variable_name,
                                                point,
                                                &argument_value);
                    if (name == "ln") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite &&
                            argument_value > T(static_cast<long long>(0))) {
                            *finite_value = t_log(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    if (name == "exp") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kNegativeInfinity) {
                            *finite_value = T(static_cast<long long>(0));
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite) {
                            *finite_value = t_exp(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    if (name == "sin" || name == "cos") {
                        if (argument_kind == SymbolicLimitProbeKind::kFinite) {
                            *finite_value = (name == "sin") ? t_sin(argument_value)
                                                            : t_cos(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    if (name == "sqrt") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite &&
                            argument_value >= T(static_cast<long long>(0))) {
                            *finite_value = t_sqrt(argument_value);
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

        Scalar p_val;
        if constexpr (std::is_floating_point_v<T>) {
            p_val = Scalar(point);
        } else if constexpr (std::is_same_v<T, Scalar>) {
            p_val = point.to_long_double();
        } else {
            p_val = point.to_long_double();
        }

        SymbolicExpression sub_expr = expression.substitute(
            variable_name,
            SymbolicExpression::number(p_val)).simplify();
        Scalar value = 0.0L;
        if (!sub_expr.is_number(&value)) {
            return SymbolicLimitProbeKind::kUnknown;
        }
        if (mymath::isfinite(value)) {
            *finite_value = T(value);
            return SymbolicLimitProbeKind::kFinite;
        }
        return value > 0.0L ? SymbolicLimitProbeKind::kPositiveInfinity
                           : SymbolicLimitProbeKind::kNegativeInfinity;
    } catch (...) {
        return SymbolicLimitProbeKind::kUnknown;
    }
}

template <typename T>
bool try_symbolic_one_to_infinity_limit(
    const SymbolicExpression& base,
    const SymbolicExpression& exponent,
    const std::string& variable_name,
    T point,
    T* result) {
    T exponent_value = T(static_cast<long long>(0));
    const SymbolicLimitProbeKind exponent_kind =
        probe_symbolic_value_at(exponent, variable_name, point, &exponent_value);
    if (!is_infinite_probe<T>(exponent_kind)) {
        return false;
    }

    T base_value = T(static_cast<long long>(0));
    const SymbolicLimitProbeKind base_kind =
        probe_symbolic_value_at(base, variable_name, point, &base_value);
    if (base_kind != SymbolicLimitProbeKind::kFinite ||
        !t_is_near_zero(base_value - T(static_cast<long long>(1)),
                        T(1e-8))) {
        return false;
    }

    const SymbolicExpression transformed =
        symbolic_expression_internal::make_function(
            "exp",
            ((base - SymbolicExpression::number(1.0L)) * exponent).simplify()).simplify();
    T transformed_value = T(static_cast<long long>(0));
    const SymbolicLimitProbeKind transformed_kind =
        probe_symbolic_value_at(transformed, variable_name, point, &transformed_value);
    if (transformed_kind == SymbolicLimitProbeKind::kFinite) {
        *result = transformed_value;
        return true;
    }
    return false;
}

template <typename T>
bool try_symbolic_lhopital_limit(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    T point,
    int direction,
    T* result,
    std::function<Scalar(const SymbolicExpression&, const std::string&, Scalar)> evaluate_at_override) {
    SymbolicExpression current = expression.simplify();
    if (current.node_->type != NodeType::kDivide) {
        return false;
    }

    SymbolicExpression numerator(current.node_->left);
    SymbolicExpression denominator(current.node_->right);
    T numerator_value = T(static_cast<long long>(0));
    T denominator_value = T(static_cast<long long>(0));
    const SymbolicLimitProbeKind numerator_kind =
        probe_symbolic_value_at(numerator, variable_name, point, &numerator_value);
    const SymbolicLimitProbeKind denominator_kind =
        probe_symbolic_value_at(denominator, variable_name, point, &denominator_value);

    const bool zero_over_zero =
        is_zero_probe(numerator_kind, numerator_value) &&
        is_zero_probe(denominator_kind, denominator_value);
    const bool infinity_over_infinity =
        is_infinite_probe<T>(numerator_kind) &&
        is_infinite_probe<T>(denominator_kind);

    if (!zero_over_zero && !infinity_over_infinity) {
        return false;
    }

    if (t_isfinite(point)) {
        series_ops::SeriesContext ctx;
        ctx.evaluate_at = evaluate_at_override ? evaluate_at_override :
            [](const SymbolicExpression& e, const std::string&, Scalar) {
                Scalar val = 0.0L;
                if (e.is_number(&val)) return val;
                return Scalar(0.0L);
            };

        struct LaurentInfo {
            int degree = 0;
            T coefficient = T(static_cast<long long>(0));
            bool valid = false;
        };

        auto get_laurent_info = [&](const SymbolicExpression& e) -> LaurentInfo {
            LaurentInfo info;
            std::vector<Scalar> coeffs;
            try {
                Scalar p_val;
                if constexpr (std::is_floating_point_v<T>) p_val = Scalar(point);
                else if constexpr (std::is_same_v<T, Scalar>) p_val = point.to_long_double();
                else p_val = point.to_long_double();

                if (series_ops::internal::evaluate_psa(e, variable_name, p_val, 4, coeffs, ctx)) {
                    for (int i = 0; i < static_cast<int>(coeffs.size()); ++i) {
                        if (!mymath::is_near_zero(coeffs[i].to_long_double(), 1e-15)) {
                            info.degree = i;
                            info.coefficient = T(coeffs[i]);
                            info.valid = true;
                            return info;
                        }
                    }
                    info.degree = 100;
                    info.coefficient = T(static_cast<long long>(0));
                    info.valid = true;
                    return info;
                }
            } catch (const series_ops::internal::PoleException& ex) {
                info.degree = ex.shift;
                info.coefficient = T(ex.leading_coefficient);
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
            T res_coeff = infoN.coefficient / infoD.coefficient;

            if (res_degree > 0) {
                *result = T(static_cast<long long>(0));
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
        SymbolicExpression n(iter_expr.node_->left);
        SymbolicExpression d(iter_expr.node_->right);
        iter_expr = (n.derivative(variable_name).simplify() /
                     d.derivative(variable_name).simplify()).simplify();

        T val = T(static_cast<long long>(0));
        const SymbolicLimitProbeKind kind =
            probe_symbolic_value_at(iter_expr, variable_name, point, &val);
        if (kind == SymbolicLimitProbeKind::kFinite) {
            *result = val;
            return true;
        }
        if (is_infinite_probe<T>(kind)) {
            *result = (kind == SymbolicLimitProbeKind::kPositiveInfinity) ? t_infinity<T>() : -t_infinity<T>();
            return true;
        }
    }

    return false;
}

template <typename T>
bool symbolic_limit_at_infinity(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    bool positive,
    T* result) {
    if (expression.node_ && expression.node_->type == NodeType::kPower) {
        T transformed = T(static_cast<long long>(0));
        if (try_symbolic_one_to_infinity_limit(SymbolicExpression(expression.node_->left),
                                               SymbolicExpression(expression.node_->right),
                                               variable_name,
                                               positive ? t_infinity<T>()
                                                        : -t_infinity<T>(),
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
        if (series_ops::internal::evaluate_psa(substituted, "t_limit_inf_tmp", 0.0L, 2, coeffs, ctx)) {
            if (!coeffs.empty() && mymath::isfinite(coeffs[0])) {
                *result = T(coeffs[0]);
                return true;
            }
        }
    } catch (const series_ops::internal::PoleException& e) {
        T lhopital_result = T(static_cast<long long>(0));
        if (try_symbolic_lhopital_limit(expression, variable_name,
                                         positive ? t_infinity<T>()
                                                  : -t_infinity<T>(),
                                         1,
                                         &lhopital_result)) {
            *result = lhopital_result;
            return true;
        }

        try {
            *result = handle_pole_limit(e.shift, T(e.leading_coefficient), 1);
            return true;
        } catch (...) {
        }
    }

    T lhopital_result = T(static_cast<long long>(0));
    if (try_symbolic_lhopital_limit(expression, variable_name,
                                     positive ? t_infinity<T>()
                                              : -t_infinity<T>(),
                                     1,
                                     &lhopital_result)) {
        *result = lhopital_result;
        return true;
    }

    return false;
}

template <typename T>
T compute_numerical_limit(
    const std::function<T(T)>& evaluate,
    T x,
    int direction) {
    auto compute_limit_at = [&](T x_target, int side) -> T {
        T richardson[14][14] = {};
        //bool row_valid[14] = {};
        T best_value = T(static_cast<long long>(0));
        T best_error = t_infinity<T>();
        bool have_best = false;

        const T base_h = t_is_effective_infinity_point(x_target)
                             ? T(1e-2)
                             : precision::sqrt_epsilon<T>() * std::max(T(static_cast<long long>(1)), t_abs(x_target));
        T adaptive_h = base_h;
        int consecutive_bad = 0;
        constexpr int kMaxBadSamples = 3;

        T prev_val = T(static_cast<long long>(0));
        bool have_prev = false;
        int oscillation_count = 0;
        T total_amplitude = T(static_cast<long long>(0));

        for (int row = 0; row < 14; ++row) {
            const T h = adaptive_h / t_pow(T(static_cast<long long>(2)), T(static_cast<long long>(row + 4)));
            T sample_x;
            if (!t_is_effective_infinity_point(x_target)) {
                sample_x = x_target + T(static_cast<long long>(side)) * h;
            } else {
                sample_x = (x_target > T(static_cast<long long>(0)) ? T(static_cast<long long>(1)) : -T(static_cast<long long>(1))) / h;
            }

            T val = T(static_cast<long long>(0));
            try {
                val = evaluate(sample_x);
            } catch (...) {
                adaptive_h *= T(0.5L);
                consecutive_bad++;
                if (consecutive_bad >= kMaxBadSamples) {
                    throw std::runtime_error("limit did not converge (sampling failures)");
                }
                continue;
            }

            if (!t_isfinite(val)) {
                if (have_prev && t_isfinite(prev_val)) {
                    if (prev_val > T(1e10)) return t_infinity<T>();
                    else if (prev_val < T(-1e10)) return -t_infinity<T>();
                }
                adaptive_h *= T(0.5L);
                consecutive_bad++;
                if (consecutive_bad >= kMaxBadSamples) {
                    throw std::runtime_error("limit appears to be infinite (numerical evidence)");
                }
                continue;
            }

            if (have_prev) {
                const T diff = t_abs(val - prev_val);
                total_amplitude += diff;
                if ((val > T(static_cast<long long>(0)) && prev_val < T(static_cast<long long>(0))) ||
                    (val < T(static_cast<long long>(0)) && prev_val > T(static_cast<long long>(0)))) {
                    oscillation_count++;
                    if (oscillation_count >= 5) {
                        const T avg_amp = total_amplitude / T(static_cast<long long>(row + 1));
                        if (avg_amp > T(1e-2)) {
                            throw std::runtime_error("limit does not exist (oscillation)");
                        }
                        adaptive_h *= T(0.25);
                        oscillation_count = 0;
                    }
                } else {
                    oscillation_count = std::max(0, oscillation_count - 1);
                }
            }
            prev_val = val;
            have_prev = true;

            if (have_best && row > 0) {
                const T expected_change = best_error * T(static_cast<long long>(10)) + T(1e-10);
                const T actual_change = t_abs(val - best_value);
                if (actual_change > expected_change * T(1e6)) {
                    adaptive_h *= T(0.5L);
                    row = -1;
                    consecutive_bad++;
                    if (consecutive_bad >= kMaxBadSamples) break;
                    continue;
                }
            }

            richardson[row][0] = val;
            for (int col = 1; col <= row; ++col) {
                T p4 = t_pow(T(static_cast<long long>(2)), T(static_cast<long long>(col)));
                richardson[row][col] = (p4 * richardson[row][col - 1] - richardson[row - 1][col - 1]) / (p4 - T(static_cast<long long>(1)));
            }
            //row_valid[row] = true;

            if (row >= 1) {
                const T current_error = t_abs(richardson[row][row] - richardson[row - 1][row - 1]);
                if (!have_best || current_error < best_error) {
                    best_value = richardson[row][row];
                    best_error = current_error;
                    have_best = true;
                }
                if (best_error < precision::default_relative_tolerance<T>()) break;
            } else {
                best_value = richardson[0][0];
                have_best = true;
            }
        }

        if (!have_best) throw std::runtime_error("limit did not converge");
        return best_value;
    };

    if (direction == 0) {
        T left = compute_limit_at(x, -1);
        T right = compute_limit_at(x, 1);
        if (t_abs(left - right) > precision::default_relative_tolerance<T>() &&
            t_isfinite(left) && t_isfinite(right)) {
            throw std::runtime_error("limit does not exist (left and right limits differ)");
        }
        return (left + right) * T(0.5L);
    }
    return compute_limit_at(x, direction);
}

// 显式实例化
template SymbolicLimitProbeKind probe_symbolic_value_at<Scalar>(
    SymbolicExpression, const std::string&, Scalar, Scalar*);
template bool try_symbolic_lhopital_limit<Scalar>(
    const SymbolicExpression&, const std::string&, Scalar, int, Scalar*,
    std::function<Scalar(const SymbolicExpression&, const std::string&, Scalar)>);
template bool symbolic_limit_at_infinity<Scalar>(
    const SymbolicExpression&, const std::string&, bool, Scalar*);
template Scalar compute_numerical_limit<Scalar>(
    const std::function<Scalar(Scalar)>&, Scalar, int);

template SymbolicLimitProbeKind probe_symbolic_value_at<double>(
    SymbolicExpression, const std::string&, double, double*);
template bool try_symbolic_lhopital_limit<double>(
    const SymbolicExpression&, const std::string&, double, int, double*,
    std::function<Scalar(const SymbolicExpression&, const std::string&, Scalar)>);
template bool symbolic_limit_at_infinity<double>(
    const SymbolicExpression&, const std::string&, bool, double*);
template double compute_numerical_limit<double>(
    const std::function<double(double)>&, double, int);

}  // namespace numeric

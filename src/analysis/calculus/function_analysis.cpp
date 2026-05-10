/**
 * @file function_analysis.cpp
 * @brief 函数分析实现
 *
 * 实现数值微积分运算：
 * - 数值微分（自适应中心差分 + Richardson 外推）
 * - 极限计算（逐步逼近法）
 * - 数值积分（自适应 Gauss-Kronrod G7-K15）
 * - 极值点查找（导数变号检测 + 二分法）
 */

#include "analysis/calculus/function_analysis.h"
#include "analysis/base/precision_constants.h"

#include "core/api/calculator.h"
#include "app/scalar_type.h"
#include "math/mymath.h"
#include "analysis/modules/series_module.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "precise/precise_parser.h"
#include "statistics/probability.h"
#include "precise/precise_decimal.h"

#include <algorithm>
#include <cctype>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <type_traits>

namespace {

using Scalar = mymath::Scalar;

/** @brief 泛型绝对值函数 */

Scalar t_abs(const Scalar& val) { return mymath::abs(val); }

/** @brief 泛型平方根函数 */

Scalar t_sqrt(const Scalar& val) { return mymath::sqrt(val); }

/** @brief 泛型幂函数 */

Scalar t_pow(const Scalar& base, const Scalar& exponent) { return mymath::pow(base, exponent); }

/** @brief 泛型有限值检查 */

bool t_isfinite(const Scalar& val) { return mymath::isfinite(val); }

/** @brief 泛型正弦函数 */

Scalar t_sin(const Scalar& val) { return mymath::sin(val); }

/** @brief 泛型余弦函数 */

Scalar t_cos(const Scalar& val) { return mymath::cos(val); }

/** @brief 泛型正切函数 */

Scalar t_tan(const Scalar& val) { return mymath::tan(val); }

/** @brief 泛型自然对数函数 */

Scalar t_log(const Scalar& val) { return mymath::ln(val); }

/** @brief 泛型指数函数 */

Scalar t_exp(const Scalar& val) { return mymath::exp(val); }

/** @brief 泛型双曲正弦函数 */

Scalar t_sinh(const Scalar& val) { return mymath::sinh(val); }

/** @brief 泛型双曲余弦函数 */

Scalar t_cosh(const Scalar& val) { return mymath::cosh(val); }

/** @brief 泛型双曲正切函数 */

Scalar t_tanh(const Scalar& val) { return mymath::tanh(val); }

/** @brief 泛型圆周率 */

Scalar t_pi() { return mymath::pi(); }

/** @brief 泛型无穷大 */

Scalar t_infinity() { return Scalar("1e1000"); }


bool t_is_effective_infinity_point(const Scalar& val) {
    return !mymath::isfinite(val);
}


bool symbolic_limit_at_infinity(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                bool positive,
                                Scalar* result);

/** @brief 泛型整数判断 */

bool t_is_integer(const Scalar& val) { return mymath::is_integer(val.to_long_double()); }

/** @brief 泛型零判断 */

bool t_is_near_zero(const Scalar& val, const Scalar& eps) {
    return t_abs(val) <= eps;
}

/** @brief 泛型四舍五入到 long long */

long long t_llround(const Scalar& val) {
    return static_cast<long long>(std::llround(val.to_long_double()));
}

/** @brief 导数计算的基准步长 - 使用精度感知常量 */

Scalar kDerivativeBaseStep_v() {
    return precision::sqrt_epsilon<Scalar>();
}

/** @brief 极限计算的初始步长 */

Scalar kLimitInitialStep_v() {
    return Scalar(1e-1);
}

/** @brief 极限计算的收敛容差 - 使用精度感知常量 */

Scalar kLimitTolerance_v() {
    return precision::default_relative_tolerance<Scalar>();
}

/** @brief 根查找的收敛容差 - 使用精度感知常量 */

Scalar kRootTolerance_v() {
    return precision::newton_tolerance<Scalar>();
}

/** @brief 数值积分的精度要求 - 使用精度感知常量 */

Scalar kIntegralTolerance_v() {
    return precision::default_relative_tolerance<Scalar>();
}

/** @brief 自适应积分的最大递归深度 */
constexpr int kMaxIntegralDepth = 18;


std::string format_t(const Scalar& value) {
    std::ostringstream out;
    out << std::setprecision(20) << value.to_long_double();
    return out.str();
}


void compensated_add(Scalar value,
                     Scalar* sum,
                     Scalar* compensation) {
    const Scalar adjusted = value - *compensation;
    const Scalar next = *sum + adjusted;
    *compensation = (next - *sum) - adjusted;
    *sum = next;
}


Scalar compensated_pair_sum(Scalar lhs, Scalar rhs) {
    Scalar sum = Scalar(static_cast<long long>(0));
    Scalar compensation = Scalar(static_cast<long long>(0));
    compensated_add(lhs, &sum, &compensation);
    compensated_add(rhs, &sum, &compensation);
    return sum;
}


Scalar scale_aware_step(Scalar x) {
    const Scalar scale = std::max(Scalar(static_cast<long long>(1)), t_abs(x));
    return kDerivativeBaseStep_v() * scale;
}


Scalar central_difference_step_value(Scalar scale, Scalar factor) {
    Scalar base_step = std::max(precision::sqrt_epsilon<Scalar>(), Scalar(1e-6L));
    return std::max(base_step * scale, kDerivativeBaseStep_v() * scale * factor);
}


Scalar numeric_control_value(const char*, Scalar fallback) {
    return Scalar(fallback);
}


Scalar derivative_quarter_power_scale(const Scalar& value) {
    return t_pow(value, Scalar(0.25));
}


Scalar relative_tolerance(Scalar baseline, Scalar scale) {
    return baseline * std::max(Scalar(static_cast<long long>(1)), scale);
}


Scalar limit_step_scale(Scalar x) {
    return kLimitInitialStep_v() * std::max(Scalar(static_cast<long long>(1)), t_abs(x));
}


bool same_extremum_x(Scalar lhs, Scalar rhs) {
    const Scalar scale = std::max({Scalar(static_cast<long long>(1)), t_abs(lhs), t_abs(rhs)});
    return t_abs(lhs - rhs) <= numeric_control_value("1e-4", 1e-4) * scale;
}


Scalar require_finite_integral(Scalar value) {
    if (!t_isfinite(value)) {
        throw std::runtime_error("integral did not converge");
    }
    return value;
}


void reject_divergent_transformed_endpoint(
    const std::function<Scalar(Scalar)>& transformed,
    bool check_left,
    bool check_right) {
    const Scalar offsets[] = {
        numeric_control_value("1e-3", 1e-3),
        numeric_control_value("2e-3", 2e-3),
        numeric_control_value("5e-3", 5e-3),
        numeric_control_value("1e-4", 1e-4),
        numeric_control_value("5e-4", 5e-4),
        numeric_control_value("1e-5", 1e-5),
        numeric_control_value("1e-6", 1e-6),
    };
    auto check_at = [&](Scalar t) {
        Scalar value = transformed(t);
        if (!t_isfinite(value) || t_abs(value) > Scalar(static_cast<long long>(5000))) {
            throw std::runtime_error("integral did not converge (divergence detected at endpoint)");
        }
    };

    for (Scalar offset : offsets) {
        if (check_left) {
            check_at(offset);
        }
        if (check_right) {
            check_at(Scalar(static_cast<long long>(1)) - offset);
        }
    }
}


void reject_persistent_tail_oscillation(
    const std::function<Scalar(Scalar)>& function,
    Scalar start) {
    bool have_previous = false;
    Scalar previous = Scalar(0);
    int sign_changes = 0;
    int significant_samples = 0;
    Scalar max_abs = Scalar(0);

    for (int i = 0; i < 64; ++i) {
        const Scalar offset = Scalar(10) + Scalar(i) * (t_pi() / Scalar(2));
        const Scalar value = function(start + offset);
        if (!t_isfinite(value)) {
            throw std::runtime_error("integral did not converge (non-finite tail sample)");
        }

        const Scalar abs_value = t_abs(value);
        max_abs = std::max(max_abs, abs_value);
        if (abs_value > Scalar(1e-3L)) ++significant_samples;

        if (have_previous &&
            abs_value > Scalar(1e-3L) &&
            t_abs(previous) > Scalar(1e-3L) &&
            ((value > Scalar(0)) != (previous > Scalar(0)))) {
            ++sign_changes;
        }

        previous = value;
        have_previous = true;
    }
    if (sign_changes >= 12 &&
        significant_samples >= 32 &&
        max_abs > Scalar(1e-2L)) {
        throw std::runtime_error("integral did not converge (persistent tail oscillation detected)");
    }
}

// 自适应 Simpson 积分辅助函数（用于 callable，支持高精度）

Scalar simpson_rule_callable(const std::function<Scalar(Scalar)>& func, Scalar a, Scalar b) {
    const Scalar h = (b - a) / Scalar(static_cast<long long>(2));
    const Scalar fa = func(a);
    const Scalar fb = func(b);
    const Scalar fc = func((a + b) / Scalar(static_cast<long long>(2)));
    return h / Scalar(static_cast<long long>(3)) * (fa + Scalar(static_cast<long long>(4)) * fc + fb);
}


Scalar adaptive_simpson_callable_recursive(const std::function<Scalar(Scalar)>& func,
                                       Scalar a, Scalar b, Scalar whole, Scalar left, Scalar right, Scalar eps, int depth) {
    const Scalar c = (a + b) / Scalar(static_cast<long long>(2));
    const Scalar combined = left + right;
    const Scalar error = t_abs(combined - whole) / Scalar(static_cast<long long>(15));

    const Scalar scale = std::max(Scalar(static_cast<long long>(1)), t_abs(combined));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        return combined + (combined - whole) / Scalar(static_cast<long long>(15));
    }

    const Scalar d = (a + c) / Scalar(static_cast<long long>(2));
    const Scalar e = (c + b) / Scalar(static_cast<long long>(2));
    const Scalar left_left = simpson_rule_callable(func, a, d);
    const Scalar left_right = simpson_rule_callable(func, d, c);
    const Scalar right_left = simpson_rule_callable(func, c, e);
    const Scalar right_right = simpson_rule_callable(func, e, b);

    return adaptive_simpson_callable_recursive(func, a, c, left, left_left, left_right, eps / Scalar(static_cast<long long>(2)), depth - 1) +
           adaptive_simpson_callable_recursive(func, c, b, right, right_left, right_right, eps / Scalar(static_cast<long long>(2)), depth - 1);
}


Scalar adaptive_simpson_callable(const std::function<Scalar(Scalar)>& func, Scalar left, Scalar right, Scalar eps, int max_depth) {
    const Scalar c = (left + right) / Scalar(static_cast<long long>(2));
    const Scalar whole = simpson_rule_callable(func, left, right);
    const Scalar left_val = simpson_rule_callable(func, left, c);
    const Scalar right_val = simpson_rule_callable(func, c, right);
    return adaptive_simpson_callable_recursive(func, left, right, whole, left_val, right_val, eps, max_depth);
}


Scalar gauss_kronrod_15_callable(const std::function<Scalar(Scalar)>& function,
                                 Scalar left,
                                 Scalar right,
                                 Scalar* error_estimate) {
    static const Scalar kNodes[] = {
        Scalar(0.9914553711208126),
        Scalar(0.9491079123427585),
        Scalar(0.8648644233597691),
        Scalar(0.7415311855993945),
        Scalar(0.5860872354676911),
        Scalar(0.4058451513773972),
        Scalar(0.2077849550078985),
        Scalar(static_cast<long long>(0)),
    };
    static const Scalar kKronrodWeights[] = {
        Scalar(0.02293532201052922),
        Scalar(0.06309209262997855),
        Scalar(0.1047900103222502),
        Scalar(0.1406532597155259),
        Scalar(0.1690047266392679),
        Scalar(0.1903505780647854),
        Scalar(0.2044329400752989),
        Scalar(0.2094821410847278),
    };
    static const Scalar kGaussWeights[] = {
        Scalar(static_cast<long long>(0)),
        Scalar(0.1294849661688697),
        Scalar(static_cast<long long>(0)),
        Scalar(0.2797053914892767),
        Scalar(static_cast<long long>(0)),
        Scalar(0.3818300505051189),
        Scalar(static_cast<long long>(0)),
        Scalar(0.4179591836734694),
    };

    const Scalar center = (left + right) * Scalar(0.5L);
    const Scalar half_width = (right - left) * Scalar(0.5L);
    Scalar kronrod_sum = Scalar(static_cast<long long>(0));
    Scalar gauss_sum = Scalar(static_cast<long long>(0));
    Scalar kronrod_compensation = Scalar(static_cast<long long>(0));
    Scalar gauss_compensation = Scalar(static_cast<long long>(0));

    for (int i = 0; i < 8; ++i) {
        if (t_is_near_zero(kNodes[i], Scalar(static_cast<long long>(0)))) {
            const Scalar value = function(center);
            compensated_add(kKronrodWeights[i] * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(kGaussWeights[i] * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const Scalar offset = half_width * kNodes[i];
        const Scalar left_value = function(center - offset);
        const Scalar right_value = function(center + offset);
        const Scalar pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(kKronrodWeights[i] * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(kGaussWeights[i] * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const Scalar kronrod = half_width * kronrod_sum;
    const Scalar gauss = half_width * gauss_sum;
    *error_estimate = t_abs(kronrod - gauss);
    return kronrod;
}


Scalar adaptive_gauss_kronrod_callable_recursive(
    const std::function<Scalar(Scalar)>& function,
    Scalar left,
    Scalar right,
    Scalar eps,
    Scalar whole,
    Scalar error,
    int depth) {
    // 检查结果是否有效
    if (!t_isfinite(whole) || !t_isfinite(error)) {
        throw std::runtime_error("integral did not converge (non-finite value encountered)");
    }

    // 检查区间是否过小，避免数值问题
    const Scalar interval_width = t_abs(right - left);
    const Scalar interval_scale = std::max(t_abs(left), t_abs(right));
    const Scalar min_width = precision::min_step_size<Scalar>(interval_scale);
    if (interval_width < min_width) {
        return whole;
    }

    const Scalar scale = std::max(Scalar(static_cast<long long>(1)), t_abs(whole));
    const Scalar tol = relative_tolerance(eps, scale);
    if (error <= tol) {
        return whole;
    }
    if (depth <= 0) {
        if (error > tol * Scalar(1e4L)) { // 严重不收敛
            throw std::runtime_error("integral did not converge (max depth reached with large error)");
        }
        return whole;
    }

    const Scalar mid = (left + right) * Scalar(0.5L);
    Scalar left_error = Scalar(static_cast<long long>(0));
    Scalar right_error = Scalar(static_cast<long long>(0));
    const Scalar left_area =
        gauss_kronrod_15_callable(function, left, mid, &left_error);
    const Scalar right_area =
        gauss_kronrod_15_callable(function, mid, right, &right_error);

    // 检查子区间结果是否有效
    if (!t_isfinite(left_area) || !t_isfinite(right_area)) {
        throw std::runtime_error("integral did not converge (non-finite value in subinterval)");
    }

    const Scalar left_result =
        adaptive_gauss_kronrod_callable_recursive(function,
                                                  left,
                                                  mid,
                                                  eps * Scalar(0.5L),
                                                  left_area,
                                                  left_error,
                                                  depth - 1);
    const Scalar right_result =
        adaptive_gauss_kronrod_callable_recursive(function,
                                                  mid,
                                                  right,
                                                  eps * Scalar(0.5L),
                                                  right_area,
                                                  right_error,
                                                  depth - 1);
    return compensated_pair_sum(left_result, right_result);
}


Scalar adaptive_gauss_kronrod_callable(const std::function<Scalar(Scalar)>& function,
                                       Scalar left,
                                       Scalar right,
                                       Scalar eps,
                                       int depth) {
    Scalar error = Scalar(static_cast<long long>(0));
    const Scalar whole = gauss_kronrod_15_callable(function, left, right, &error);
    return require_finite_integral(
        adaptive_gauss_kronrod_callable_recursive(function,
                                                  left,
                                                  right,
                                                  eps,
                                                  whole,
                                                  error,
                                                  depth));
}

bool is_valid_analysis_variable_name(const std::string& name) {
    if (name.empty() ||
        !std::isalpha(static_cast<unsigned char>(name.front()))) {
        return false;
    }

    for (char ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            return false;
        }
    }

    return true;
}

enum class SymbolicLimitProbeKind {
    kFinite,
    kPositiveInfinity,
    kNegativeInfinity,
    kUnknown,
};


SymbolicLimitProbeKind probe_symbolic_value_at(
    SymbolicExpression expression,
    const std::string& variable_name,
    Scalar point,
    Scalar* finite_value);

bool is_infinite_probe(SymbolicLimitProbeKind kind) {
    return kind == SymbolicLimitProbeKind::kPositiveInfinity ||
           kind == SymbolicLimitProbeKind::kNegativeInfinity;
}


bool try_symbolic_one_to_infinity_limit(const SymbolicExpression& base,
                                        const SymbolicExpression& exponent,
                                        const std::string& variable_name,
                                        Scalar point,
                                        Scalar* result) {
    Scalar exponent_value = Scalar(static_cast<long long>(0));
    const SymbolicLimitProbeKind exponent_kind =
        probe_symbolic_value_at(exponent, variable_name, point, &exponent_value);
    if (!is_infinite_probe(exponent_kind)) {
        return false;
    }

    Scalar base_value = Scalar(static_cast<long long>(0));
    const SymbolicLimitProbeKind base_kind =
        probe_symbolic_value_at(base, variable_name, point, &base_value);
    if (base_kind != SymbolicLimitProbeKind::kFinite ||
        !t_is_near_zero(base_value - Scalar(static_cast<long long>(1)),
                        numeric_control_value("1e-8", 1e-8))) {
        return false;
    }

    const SymbolicExpression transformed =
        symbolic_expression_internal::make_function(
            "exp",
            ((base - SymbolicExpression::number(Scalar(1.0L))) * exponent).simplify()).simplify();
    const SymbolicExpression product =
        ((base - SymbolicExpression::number(Scalar(1.0L))) * exponent).simplify();
    Scalar product_limit = Scalar(static_cast<long long>(0));
    if (symbolic_limit_at_infinity(product, variable_name, point > Scalar(static_cast<long long>(0)), &product_limit)) {
        *result = t_exp(product_limit);
        return true;
    }

    Scalar transformed_value = Scalar(static_cast<long long>(0));
    const SymbolicLimitProbeKind transformed_kind =
        probe_symbolic_value_at(transformed, variable_name, point, &transformed_value);
    if (transformed_kind == SymbolicLimitProbeKind::kFinite) {
        *result = transformed_value;
        return true;
    }
    return false;
}


SymbolicLimitProbeKind probe_symbolic_value_at(
    SymbolicExpression expression,
    const std::string& variable_name,
    Scalar point,
    Scalar* finite_value) {
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
                        *finite_value = Scalar(value);
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kVariable:
                    if (node->text == variable_name) {
                        return point > Scalar(static_cast<long long>(0)) ? SymbolicLimitProbeKind::kPositiveInfinity
                                           : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                case NodeType::kNegate: {
                    Scalar child_value = Scalar(static_cast<long long>(0));
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
                    Scalar left_value = Scalar(static_cast<long long>(0));
                    Scalar right_value = Scalar(static_cast<long long>(0));
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
                    Scalar left_value = Scalar(static_cast<long long>(0));
                    Scalar right_value = Scalar(static_cast<long long>(0));
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
                        t_is_near_zero(left_value, kLimitTolerance_v())) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    if (right_kind == SymbolicLimitProbeKind::kFinite &&
                        t_is_near_zero(right_value, kLimitTolerance_v())) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    auto sign_of = [](SymbolicLimitProbeKind kind, Scalar value) {
                        if (kind == SymbolicLimitProbeKind::kFinite) {
                            return value >= Scalar(static_cast<long long>(0)) ? 1 : -1;
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
                    Scalar left_value = Scalar(static_cast<long long>(0));
                    Scalar right_value = Scalar(static_cast<long long>(0));
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
                        !t_is_near_zero(right_value, kLimitTolerance_v())) {
                        *finite_value = left_value / right_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        is_infinite_probe(right_kind)) {
                        *finite_value = Scalar(static_cast<long long>(0));
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (is_infinite_probe(left_kind) &&
                        right_kind == SymbolicLimitProbeKind::kFinite &&
                        !t_is_near_zero(right_value, kLimitTolerance_v())) {
                        const bool positive =
                            (left_kind == SymbolicLimitProbeKind::kPositiveInfinity) ==
                            (right_value > Scalar(static_cast<long long>(0)));
                        return positive ? SymbolicLimitProbeKind::kPositiveInfinity
                                        : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kPower: {
                    const SymbolicExpression base(node->left);
                    const SymbolicExpression exponent(node->right);
                    Scalar one_to_infinity_value = Scalar(static_cast<long long>(0));
                    if (try_symbolic_one_to_infinity_limit(base,
                                                           exponent,
                                                           variable_name,
                                                           point,
                                                           &one_to_infinity_value)) {
                        *finite_value = one_to_infinity_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }

                    Scalar exponent_number = 0.0L;
                    if (!exponent.is_number(&exponent_number)) {
                        Scalar base_value = Scalar(static_cast<long long>(0));
                        Scalar exponent_value = Scalar(static_cast<long long>(0));
                        const SymbolicLimitProbeKind base_kind =
                            probe_symbolic_value_at(base,
                                                    variable_name,
                                                    point,
                                                    &base_value);
                        const SymbolicLimitProbeKind exponent_kind =
                            probe_symbolic_value_at(exponent,
                                                    variable_name,
                                                    point,
                                                    &exponent_value);
                        if (base_kind == SymbolicLimitProbeKind::kFinite &&
                            exponent_kind == SymbolicLimitProbeKind::kFinite &&
                            base_value > Scalar(static_cast<long long>(0))) {
                            *finite_value = t_exp(exponent_value * t_log(base_value));
                            return t_isfinite(*finite_value)
                                       ? SymbolicLimitProbeKind::kFinite
                                       : SymbolicLimitProbeKind::kUnknown;
                        }
                        return SymbolicLimitProbeKind::kUnknown;
                    }

                    const Scalar exponent_value(exponent_number);
                    Scalar base_value = Scalar(static_cast<long long>(0));
                    const SymbolicLimitProbeKind base_kind =
                        probe_symbolic_value_at(base,
                                                variable_name,
                                                point,
                                                &base_value);
                    if (base_kind == SymbolicLimitProbeKind::kFinite) {
                        *finite_value = t_pow(base_value, exponent_value);
                        return t_isfinite(*finite_value)
                                   ? SymbolicLimitProbeKind::kFinite
                                   : SymbolicLimitProbeKind::kUnknown;
                    }
                    if (is_infinite_probe(base_kind)) {
                        if (exponent_number > 0.0L) {
                            if (base_kind == SymbolicLimitProbeKind::kNegativeInfinity &&
                                t_is_integer(exponent_value) &&
                                t_llround(exponent_value) % 2 != 0) {
                                return SymbolicLimitProbeKind::kNegativeInfinity;
                            }
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (exponent_number < 0.0L) {
                            *finite_value = Scalar(static_cast<long long>(0));
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kFunction: {
                    const std::string& name = node->text;
                    Scalar argument_value = Scalar(static_cast<long long>(0));
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
                            argument_value > Scalar(static_cast<long long>(0))) {
                            *finite_value = t_log(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    if (name == "exp") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kNegativeInfinity) {
                            *finite_value = Scalar(static_cast<long long>(0));
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
                    if (name == "tan") {
                        if (argument_kind == SymbolicLimitProbeKind::kFinite) {
                            *finite_value = t_tan(argument_value);
                            return t_isfinite(*finite_value)
                                       ? SymbolicLimitProbeKind::kFinite
                                       : SymbolicLimitProbeKind::kUnknown;
                        }
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    if (name == "sqrt") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite &&
                            argument_value >= Scalar(static_cast<long long>(0))) {
                            *finite_value = t_sqrt(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    if (name == "sinh") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kNegativeInfinity) {
                            return SymbolicLimitProbeKind::kNegativeInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite) {
                            *finite_value = t_sinh(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    if (name == "cosh") {
                        if (is_infinite_probe(argument_kind)) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite) {
                            *finite_value = t_cosh(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    if (name == "tanh") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            *finite_value = Scalar(static_cast<long long>(1));
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kNegativeInfinity) {
                            *finite_value = -Scalar(static_cast<long long>(1));
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite) {
                            *finite_value = t_tanh(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kVector:
                case NodeType::kTensor:
                case NodeType::kDifferentialOp:
                case NodeType::kRootOf:
                    return SymbolicLimitProbeKind::kUnknown;
            }
            return SymbolicLimitProbeKind::kUnknown;
        }

        Scalar p_val = point;

        SymbolicExpression sub_expr = expression.substitute(
            variable_name,
            SymbolicExpression::number(p_val)).simplify();
        Scalar value = 0.0L;
        if (!sub_expr.is_number(&value)) {
            return SymbolicLimitProbeKind::kUnknown;
        }
        if (mymath::isfinite(value)) {
            *finite_value = Scalar(value);
            return SymbolicLimitProbeKind::kFinite;
        }
        return value > 0.0L ? SymbolicLimitProbeKind::kPositiveInfinity
                           : SymbolicLimitProbeKind::kNegativeInfinity;
    } catch (...) {
        return SymbolicLimitProbeKind::kUnknown;
    }
}


bool is_zero_probe(SymbolicLimitProbeKind kind, Scalar value) {
    return kind == SymbolicLimitProbeKind::kFinite &&
           t_is_near_zero(value, kLimitTolerance_v());
}

/**
 * @brief 处理极点极限
 */

Scalar handle_pole_limit(int shift, Scalar leading_coefficient, int direction) {
    if (direction == 0) {
        // 双侧极限：只有当 shift 为偶数时才存在
        if (shift % 2 == 0) {
            throw std::runtime_error("limit diverges to infinity");
        } else {
            throw std::runtime_error("two-sided limit does not exist (pole with odd shift)");
        }
    } else if (direction == 1) {
        // 右极限：(x - x0) > 0，符号不变
        throw std::runtime_error("limit diverges to infinity");
    } else {
        // 左极限：(x - x0) < 0，奇数 shift 时符号翻转
        throw std::runtime_error("limit diverges to infinity");
    }
}


bool try_symbolic_lhopital_limit(const SymbolicExpression& expression,
                                 const std::string& variable_name,
                                 Scalar point,
                                 int direction,
                                 Scalar* result,
                                 std::function<Scalar(const SymbolicExpression&, const std::string&, Scalar)> evaluate_at_override = nullptr) {
    SymbolicExpression current = expression.simplify();
    if (current.node_->type != NodeType::kDivide) {
        return false;
    }

    SymbolicExpression numerator(current.node_->left);
    SymbolicExpression denominator(current.node_->right);
    Scalar numerator_value = Scalar(static_cast<long long>(0));
    Scalar denominator_value = Scalar(static_cast<long long>(0));
    const SymbolicLimitProbeKind numerator_kind =
        probe_symbolic_value_at(numerator,
                                variable_name,
                                point,
                                &numerator_value);
    const SymbolicLimitProbeKind denominator_kind =
        probe_symbolic_value_at(denominator,
                                variable_name,
                                point,
                                &denominator_value);

    const bool zero_over_zero =
        is_zero_probe(numerator_kind, numerator_value) &&
        is_zero_probe(denominator_kind, denominator_value);
    const bool infinity_over_infinity =
        is_infinite_probe(numerator_kind) &&
        is_infinite_probe(denominator_kind);

    if (!zero_over_zero && !infinity_over_infinity) {
        return false;
    }

    // 对于有限点，使用 PSA 提取 Laurent 信息
    if (t_isfinite(point)) {
        series_ops::SeriesContext ctx;
        // 使用传入的 evaluate_at 回调，或默认实现
        ctx.evaluate_at = evaluate_at_override ? evaluate_at_override :
            [](const SymbolicExpression& e, const std::string& /*v*/, Scalar /*p*/) {
                Scalar val = 0.0L;
                if (e.is_number(&val)) return val;
                return Scalar(0.0L);
            };

        struct LaurentInfo {
            int degree = 0;
            Scalar coefficient = Scalar(static_cast<long long>(0));
            bool valid = false;
        };

        auto get_laurent_info = [&](const SymbolicExpression& e) -> LaurentInfo {
            LaurentInfo info;
            std::vector<Scalar> coeffs;
            try {
                Scalar p_val = point.to_long_double();

                if (series_ops::internal::evaluate_psa(e, variable_name, p_val, 4, coeffs, ctx)) {
                    for (int i = 0; i < static_cast<int>(coeffs.size()); ++i) {
                        if (!mymath::is_near_zero(coeffs[i].to_long_double(), 1e-15)) {
                            info.degree = i;
                            info.coefficient = Scalar(coeffs[i]);
                            info.valid = true;
                            return info;
                        }
                    }
                    info.degree = 100;
                    info.coefficient = Scalar(static_cast<long long>(0));
                    info.valid = true;
                    return info;
                }
            } catch (const series_ops::internal::PoleException& ex) {
                info.degree = ex.shift;
                info.coefficient = Scalar(ex.leading_coefficient);
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
                *result = Scalar(static_cast<long long>(0));
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

    // 洛必达法则深度限制：增加到 5 层，对于复杂函数可以更多次求导
    // 对于 float128 精度，可以承受更多次数值求导的误差累积
    static constexpr int kMaxLhopitalDepth = 5;
    SymbolicExpression iter_expr = current;
    for (int depth = 0; depth < kMaxLhopitalDepth; ++depth) {
        SymbolicExpression n(iter_expr.node_->left);
        SymbolicExpression d(iter_expr.node_->right);
        iter_expr = (n.derivative(variable_name).simplify() /
                     d.derivative(variable_name).simplify()).simplify();

        Scalar val = Scalar(static_cast<long long>(0));
        const SymbolicLimitProbeKind kind =
            probe_symbolic_value_at(iter_expr, variable_name, point, &val);
        if (kind == SymbolicLimitProbeKind::kFinite) {
            *result = val;
            return true;
        }
        if (is_infinite_probe(kind)) {
            throw std::runtime_error("limit diverges to infinity");
        }
    }

    return false;
}


bool symbolic_limit_at_infinity(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                bool positive,
                                Scalar* result) {
    if (expression.node_ && expression.node_->type == NodeType::kPower) {
        Scalar transformed = Scalar(static_cast<long long>(0));
        if (try_symbolic_one_to_infinity_limit(SymbolicExpression(expression.node_->left),
                                               SymbolicExpression(expression.node_->right),
                                               variable_name,
                                               positive ? t_infinity()
                                                        : -t_infinity(),
                                               &transformed)) {
            *result = transformed;
            return true;
        }
    }

    series_ops::SeriesContext ctx;
    ctx.evaluate_at = [](const SymbolicExpression& e, const std::string& /*v*/, Scalar /*p*/) {
        Scalar val = 0.0L;
        if (e.is_number(&val)) return val;
        return Scalar(0.0L);
    };

    SymbolicExpression t_var = SymbolicExpression::variable("t_limit_inf_tmp");
    SymbolicExpression inv_t;
    if (positive) {
        inv_t = SymbolicExpression::number(Scalar(1.0L)) / t_var;
    } else {
        inv_t = SymbolicExpression::number(Scalar(-1.0L)) / t_var;
    }
    SymbolicExpression substituted = expression.substitute(variable_name, inv_t).simplify();

    std::vector<Scalar> coeffs;
    try {
        if (series_ops::internal::evaluate_psa(substituted, "t_limit_inf_tmp", Scalar(0.0L), 2, coeffs, ctx)) {
            if (!coeffs.empty() && mymath::isfinite(coeffs[0])) {
                *result = Scalar(coeffs[0]);
                return true;
            }
        }
    } catch (const series_ops::internal::PoleException& e) {
        Scalar lhopital_result = Scalar(static_cast<long long>(0));
        if (try_symbolic_lhopital_limit(expression, variable_name,
                                         positive ? t_infinity()
                                                  : -t_infinity(),
                                         1,
                                         &lhopital_result)) {
            *result = lhopital_result;
            return true;
        }

        try {
            Scalar large_x = positive ? 1e12 : -1e12;
            Scalar val = ctx.evaluate_at(expression, variable_name, large_x);
            if (mymath::isfinite(val) && mymath::abs(val) < 1e10) {
                return false;
            }
        } catch (...) {
        }

        try {
            *result = handle_pole_limit(e.shift, Scalar(e.leading_coefficient), 1);
            return true;
        } catch (...) {
        }
    }

    Scalar lhopital_result = Scalar(static_cast<long long>(0));
    if (try_symbolic_lhopital_limit(expression, variable_name,
                                     positive ? t_infinity()
                                              : -t_infinity(),
                                     1,
                                     &lhopital_result)) {
        *result = lhopital_result;
        return true;
    }

    return false;
}

}  // namespace


FunctionAnalysis::FunctionAnalysis(std::string variable_name)
    : variable_name_(std::move(variable_name)) {
    if (!is_valid_analysis_variable_name(variable_name_)) {
        throw std::runtime_error("invalid variable name for custom function");
    }
}


FunctionAnalysis::FunctionAnalysis(const FunctionAnalysis& other)
    : expression_(other.expression_),
      variable_name_(other.variable_name_),
      evaluator_(other.evaluator_),
      fallback_calculator_(other.fallback_calculator_),
      variable_lookup_(other.variable_lookup_) {}


FunctionAnalysis& FunctionAnalysis::operator=(const FunctionAnalysis& other) {
    if (this != &other) {
        expression_ = other.expression_;
        variable_name_ = other.variable_name_;
        evaluator_ = other.evaluator_;
        fallback_calculator_ = other.fallback_calculator_;
        variable_lookup_ = other.variable_lookup_;
        evaluation_cache_entries_.clear();
        evaluation_cache_index_.clear();
    }
    return *this;
}


FunctionAnalysis::FunctionAnalysis(FunctionAnalysis&& other) noexcept = default;


FunctionAnalysis& FunctionAnalysis::operator=(FunctionAnalysis&& other) noexcept = default;


FunctionAnalysis::~FunctionAnalysis() = default;


void FunctionAnalysis::define(const std::string& expression) {
    if (expression.empty()) {
        throw std::runtime_error("function expression cannot be empty");
    }

    expression_ = expression;
    evaluator_ = nullptr;
    fallback_calculator_.reset();
    evaluation_cache_entries_.clear();
    evaluation_cache_index_.clear();
}


void FunctionAnalysis::set_evaluator(std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> evaluator) {
    evaluator_ = std::move(evaluator);
    fallback_calculator_.reset();
    evaluation_cache_entries_.clear();
    evaluation_cache_index_.clear();
}


void FunctionAnalysis::set_variable_lookup(std::function<Scalar(const std::string&)> lookup) {
    variable_lookup_ = std::move(lookup);
}


Scalar FunctionAnalysis::evaluate(Scalar x) const {
    return evaluate_with_variable(x);
}


Scalar FunctionAnalysis::derivative(Scalar x) const {
    const Scalar scale = std::max(Scalar(1.0L), mymath::abs(x));
    const Scalar center = evaluate_with_variable(x);
    if (!mymath::isfinite(center)) {
        throw std::runtime_error("derivative is undefined at this point");
    }

    constexpr int kLayers = 4; // 恢复为 4 层，主要利用 float128 减少舍入误差

    // 计算曲率以确定步长
    const Scalar curvature_sample_step = scale * Scalar(1e-3L);
    const Scalar curvature_probe = evaluate_with_variable(x + curvature_sample_step) -
                                   Scalar(2.0L) * center +
                                   evaluate_with_variable(x - curvature_sample_step);
    const Scalar curvature_scale =
        std::max(Scalar(1.0L),
                 mymath::abs(curvature_probe) /
                     std::max(Scalar(1e-12L), mymath::abs(center)));

    Scalar base_step = central_difference_step_value(
        scale,
        Scalar(1.0L) / derivative_quarter_power_scale(curvature_scale));

    Scalar richardson[kLayers][kLayers] = {};
    bool row_valid[kLayers] = {};
    Scalar best_value = Scalar(0.0L);
    Scalar best_error = Scalar(1e300);  // Large number instead of infinity

    for (int row = 0; row < kLayers; ++row) {
        const Scalar step = base_step / mymath::pow(Scalar(2.0L), Scalar(static_cast<long long>(row)));
        const Scalar forward_x = x + step;
        const Scalar backward_x = x - step;
        const Scalar actual_step = (forward_x - backward_x) * Scalar(0.5L);
        if (actual_step <= Scalar(0.0L)) {
            continue;
        }
        const Scalar forward = evaluate_with_variable(forward_x);
        const Scalar backward = evaluate_with_variable(backward_x);
        if (!mymath::isfinite(forward) || !mymath::isfinite(backward)) {
            continue;
        }
        richardson[row][0] = (forward - backward) / (Scalar(2.0L) * actual_step);
        row_valid[row] = mymath::isfinite(richardson[row][0]);
        if (!row_valid[row]) {
            continue;
        }
        for (int col = 1; col <= row; ++col) {
            if (!row_valid[row - 1]) {
                row_valid[row] = false;
                break;
            }
            Scalar factor = mymath::pow(Scalar(4.0L), Scalar(static_cast<long long>(col)));
            richardson[row][col] =
                richardson[row][col - 1] +
                (richardson[row][col - 1] - richardson[row - 1][col - 1]) /
                    (factor - 1.0L);
            if (!mymath::isfinite(richardson[row][col])) {
                row_valid[row] = false;
                break;
            }
        }
        if (row > 0 && row_valid[row] && row_valid[row - 1]) {
            const Scalar candidate = richardson[row][row];
            const Scalar error_estimate = mymath::abs(candidate - richardson[row - 1][row - 1]);
            if (error_estimate < best_error && mymath::isfinite(candidate)) {
                best_error = error_estimate;
                best_value = candidate;

                Scalar tol = mymath::abs(best_value) * Scalar(1e-16L);
                if (tol < Scalar(1e-18L)) tol = Scalar(1e-18L);

                if (error_estimate < tol) break;
            }
        }
    }

    if (best_error < Scalar(1e300)) {
        const Scalar side_step = std::max(Scalar(1e-7L) * scale, base_step / Scalar(64.0L));
        const Scalar left_value = evaluate_with_variable(x - side_step);
        const Scalar right_value = evaluate_with_variable(x + side_step);
        if (!mymath::isfinite(left_value) || !mymath::isfinite(right_value)) {
            throw std::runtime_error("derivative is undefined at this point");
        }
        const Scalar left_slope = (center - left_value) / side_step;
        const Scalar right_slope = (right_value - center) / side_step;
        const Scalar slope_scale =
            std::max({Scalar(1.0L), mymath::abs(left_slope), mymath::abs(right_slope), mymath::abs(best_value)});
        if (mymath::abs(left_slope - right_slope) >
            std::max(Scalar(1e-4L), Scalar(1e-5L) * slope_scale)) {
            throw std::runtime_error("derivative does not exist at this point");
        }
        return best_value;
    }

    throw std::runtime_error("derivative did not converge");
}


Scalar FunctionAnalysis::limit(Scalar x, int direction) const {
    if (direction != -1 && direction != 0 && direction != 1) {
        throw std::runtime_error("limit direction must be -1, 0, or 1");
    }

    SymbolicExpression expr;
    try {
        expr = SymbolicExpression::parse(expression_);
    } catch (...) {
        return compute_numerical_limit(x, direction);
    }

    series_ops::SeriesContext ctx;
    ctx.evaluate_at = [this](const SymbolicExpression& e, const std::string& v, Scalar p) {
        // 如果是极限变量，返回极限点的值
        if (v == variable_name_) return p;

        // 尝试从表达式中提取数值
        Scalar val = Scalar(0.0L);
        if (e.is_number(&val)) return val;

        // 如果有外部变量查找函数，尝试使用它
        if (variable_lookup_) {
            // 尝试从表达式中提取变量名
            std::string var_name = e.to_string();
            // 简单变量名检查
            bool is_simple_var = true;
            for (char c : var_name) {
                if (!std::isalpha(c) && c != '_') {
                    is_simple_var = false;
                    break;
                }
            }
            if (is_simple_var && !var_name.empty()) {
                try {
                    Scalar lookup_val = variable_lookup_(var_name);
                    // 如果查找成功（返回非零或有效值），使用它
                    return lookup_val;
                } catch (...) {
                    // 查找失败，继续
                }
            }
        }

        // 无法确定值，返回 0（保守处理）
        return Scalar(0.0L);
    };

    if (t_is_effective_infinity_point(x)) {
        bool positive = x > Scalar(0.0L);
        Scalar inf_result = Scalar(0.0L);
        if (symbolic_limit_at_infinity(expr, variable_name_, positive, &inf_result)) {
            return inf_result;
        }
    } else if (t_isfinite(x)) {
        std::vector<Scalar> coeffs;
        try {
            Scalar p_val = x;
            if (series_ops::internal::evaluate_psa(expr, variable_name_, p_val, 2, coeffs, ctx)) {
                if (!coeffs.empty()) return coeffs[0];
            }
        } catch (const series_ops::internal::PoleException& e) {
            return handle_pole_limit(e.shift, Scalar(e.leading_coefficient), direction);
        }
    }

    Scalar lhopital_value = Scalar(0.0L);
    if (direction == 0 &&
        try_symbolic_lhopital_limit(expr,
                                    variable_name_,
                                    x,
                                    direction,
                                    &lhopital_value,
                                    ctx.evaluate_at)) {
        return lhopital_value;
    }

    return compute_numerical_limit(x, direction);
}


Scalar FunctionAnalysis::compute_numerical_limit(Scalar x, int direction) const {
    // For limits at infinity, use float128 precision for better accuracy
    // when computing expressions like (1 + 1/x)^x for large x
    if (t_is_effective_infinity_point(x)) {
        const bool positive = x > 0.0L;

        // Parse the expression and substitute x = 1/t
        SymbolicExpression expr;
        try {
            expr = SymbolicExpression::parse(expression_);
        } catch (...) {
            // Fall back to direct computation if parsing fails
            goto direct_computation;
        }

        // Create substitution: x -> 1/t or x -> -1/t
        SymbolicExpression t_var = SymbolicExpression::variable("t_limit_inf_subst");
        SymbolicExpression inv_t = positive
            ? SymbolicExpression::number(Scalar(1.0L)) / t_var
            : SymbolicExpression::number(Scalar(-1.0L)) / t_var;
        SymbolicExpression substituted = expr.substitute(variable_name_, inv_t).simplify();

        // Use float128 precision for the substituted expression
        FunctionAnalysis analysis_128("t_limit_inf_subst");
        analysis_128.define(substituted.to_string());

        // Compute limit as t -> 0 using float128 precision
        try {
            Scalar result_128 = analysis_128.limit(Scalar(0.0L), direction);
            if (mymath::isfinite(result_128)) {
                return result_128;
            }
        } catch (...) {
            // Fall through to direct computation
        }
    }

direct_computation:
    auto compute_limit_at = [this](Scalar x_target, int side) {
        Scalar richardson[14][14] = {};
        //bool row_valid[14] = {};
        Scalar best_value = Scalar(0.0L);
        Scalar best_error = Scalar(1e300);
        bool have_best = false;

        const Scalar base_h = t_is_effective_infinity_point(x_target)
                             ? Scalar(1e-2L)
                             : limit_step_scale(x_target);
        Scalar adaptive_h = base_h;
        int consecutive_bad = 0;
        constexpr int kMaxBadSamples = 3;

        Scalar prev_val = Scalar(0.0L);
        bool have_prev = false;
        int oscillation_count = 0;
        Scalar total_amplitude = Scalar(0.0L);

        for (int row = 0; row < 14; ++row) {
            const Scalar h = adaptive_h / mymath::pow(Scalar(2.0L), Scalar(static_cast<long long>(row + 4)));
            Scalar sample_x;
            if (!t_is_effective_infinity_point(x_target)) {
                sample_x = x_target + Scalar(static_cast<long long>(side)) * h;
            } else {
                sample_x = (x_target > Scalar(0.0L) ? Scalar(1.0L) : Scalar(-1.0L)) / h;
            }

            Scalar val = Scalar(0.0L);
            try {
                val = evaluate_with_variable(sample_x);
            } catch (...) {
                adaptive_h *= 0.5L;
                consecutive_bad++;
                if (consecutive_bad >= kMaxBadSamples) {
                    throw std::runtime_error("limit did not converge (sampling failures)");
                }
                continue;
            }

            if (!mymath::isfinite(val)) {
                if (have_prev && mymath::isfinite(prev_val)) {
                    if (prev_val > 1e10L) return Scalar(mymath::infinity());
                    else if (prev_val < -1e10L) return Scalar(-mymath::infinity());
                }
                adaptive_h *= 0.5L;
                consecutive_bad++;
                if (consecutive_bad >= kMaxBadSamples) {
                    throw std::runtime_error("limit appears to be infinite (numerical evidence)");
                }
                continue;
            }

            if (have_prev) {
                const Scalar diff = mymath::abs(val - prev_val);
                total_amplitude += diff;
                if ((val > 0.0L && prev_val < 0.0L) || (val < 0.0L && prev_val > 0.0L)) {
                    oscillation_count++;
                    if (oscillation_count >= 5) {
                        const Scalar avg_amp = total_amplitude / (row + 1);
                        if (avg_amp > 1e-2L) {
                            throw std::runtime_error("limit does not exist (oscillation)");
                        }
                        adaptive_h *= 0.25L;
                        oscillation_count = 0;
                    }
                } else {
                    oscillation_count = std::max(0, oscillation_count - 1);
                }
            }
            prev_val = val;
            have_prev = true;

            Scalar f_val(val);
            if (have_best && row > 0) {
                const Scalar expected_change =
                    best_error * 10.0L + 1e-10L;
                const Scalar actual_change = mymath::abs(f_val - best_value);
                if (actual_change > expected_change * 1e6L) {
                    adaptive_h *= 0.5L;
                    row = -1;
                    consecutive_bad++;
                    if (consecutive_bad >= kMaxBadSamples) break;
                    continue;
                }
            }

            richardson[row][0] = f_val;
            for (int col = 1; col <= row; ++col) {
                Scalar p4 = mymath::pow(Scalar(2.0L), Scalar(static_cast<long long>(col)));
                richardson[row][col] = (p4 * richardson[row][col - 1] - richardson[row - 1][col - 1]) / (p4 - 1.0L);
            }
            //row_valid[row] = true;

            if (row >= 1) {
                const Scalar current_error = mymath::abs(richardson[row][row] - richardson[row - 1][row - 1]);
                if (!have_best || current_error < best_error) {
                    best_value = richardson[row][row];
                    best_error = current_error;
                    have_best = true;
                }
                if (best_error < 1e-18L) break;
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
        if (mymath::abs(left - right) > 1e-7L && mymath::isfinite(left) && mymath::isfinite(right)) {
            throw std::runtime_error("limit does not exist (left and right limits differ)");
        }
        return (left + right) * 0.5L;
    }
    return compute_limit_at(x, direction);
}



Scalar FunctionAnalysis::definite_integral(Scalar lower_bound,
                                           Scalar upper_bound) const {
    const Scalar diff = lower_bound - upper_bound;
    if (t_isfinite(diff) && t_is_near_zero(diff, Scalar(1e-15L))) {
        return Scalar(static_cast<long long>(0));
    }
    if (lower_bound > upper_bound) {
        return -definite_integral(upper_bound, lower_bound);
    }
    const bool lower_is_infinite = !t_isfinite(lower_bound);
    const bool upper_is_infinite = !t_isfinite(upper_bound);
    if (lower_is_infinite || upper_is_infinite) {
        if (lower_is_infinite && upper_is_infinite) {
            if (lower_bound > Scalar(static_cast<long long>(0)) || upper_bound < Scalar(static_cast<long long>(0))) {
                throw std::runtime_error("invalid infinite integration bounds");
            }
            auto transformed = [this](Scalar t) {
                const Scalar angle = t_pi() * (t - Scalar(0.5L));
                const Scalar cos_angle = t_cos(angle);
                const Scalar x = t_tan(angle);
                return evaluate_with_variable(x) * t_pi() / (cos_angle * cos_angle);
            };
            reject_divergent_transformed_endpoint(std::function<Scalar(Scalar)>(transformed), true, true);
            return adaptive_gauss_kronrod_callable(std::function<Scalar(Scalar)>(transformed),
                                                   Scalar(static_cast<long long>(0)),
                                                   Scalar(static_cast<long long>(1)),
                                                   kIntegralTolerance_v(),
                                                   kMaxIntegralDepth);
        }

        if (lower_is_infinite) {
            if (lower_bound > Scalar(static_cast<long long>(0))) {
                throw std::runtime_error("invalid infinite integration bounds");
            }
            reject_persistent_tail_oscillation(
                [this](Scalar x) { return evaluate_with_variable(x); },
                -std::max(Scalar(static_cast<long long>(1)), t_abs(upper_bound)));
            auto transformed = [this, upper_bound](Scalar t) {
                const Scalar x = upper_bound - (Scalar(static_cast<long long>(1)) - t) / t;
                return evaluate_with_variable(x) / (t * t);
            };
            reject_divergent_transformed_endpoint(std::function<Scalar(Scalar)>(transformed), true, false);
            return adaptive_gauss_kronrod_callable(std::function<Scalar(Scalar)>(transformed),
                                                   Scalar(static_cast<long long>(0)),
                                                   Scalar(static_cast<long long>(1)),
                                                   kIntegralTolerance_v(),
                                                   kMaxIntegralDepth);
        }

        if (upper_bound < Scalar(static_cast<long long>(0))) {
            throw std::runtime_error("invalid infinite integration bounds");
        }
        reject_persistent_tail_oscillation(
            [this](Scalar x) { return evaluate_with_variable(x); },
            lower_bound);
        auto transformed = [this, lower_bound](Scalar t) {
            const Scalar one_minus_t = Scalar(static_cast<long long>(1)) - t;
            if (one_minus_t < Scalar(1e-18L)) return Scalar(0);
            const Scalar x = lower_bound + t / one_minus_t;
            return evaluate_with_variable(x) / (one_minus_t * one_minus_t);
        };
        // 增加对变换后函数的振荡检测
        Scalar prev_val_scan = transformed(Scalar(0.01L));
        int sign_changes = 0;
        for (int i = 2; i < 500; ++i) {
            Scalar val = transformed(Scalar(i) / Scalar(500.0L));
            if (t_isfinite(val) && t_isfinite(prev_val_scan)) {
                if ((val > Scalar(0)) != (prev_val_scan > Scalar(0))) sign_changes++;
            }
            prev_val_scan = val;
        }
        if (sign_changes > 20) {
            throw std::runtime_error("integral did not converge (excessive oscillations detected)");
        }
        reject_divergent_transformed_endpoint(std::function<Scalar(Scalar)>(transformed), false, true);
        return adaptive_gauss_kronrod_callable(std::function<Scalar(Scalar)>(transformed),
                                               Scalar(static_cast<long long>(0)),
                                               Scalar(static_cast<long long>(1)),
                                               kIntegralTolerance_v(),
                                               kMaxIntegralDepth);
    }
    const Scalar span = t_abs(upper_bound - lower_bound);
    const Scalar scaled_eps =
        relative_tolerance(kIntegralTolerance_v(), span + t_abs(lower_bound) + t_abs(upper_bound));
    bool left_singular = false;
    bool right_singular = false;
    try {
        left_singular = !t_isfinite(evaluate_with_variable(lower_bound));
    } catch (...) {
        left_singular = true;
    }
    try {
        right_singular = !t_isfinite(evaluate_with_variable(upper_bound));
    } catch (...) {
        right_singular = true;
    }

    if (left_singular || right_singular) {
        const Scalar width = upper_bound - lower_bound;
        auto transformed = [this, lower_bound, upper_bound, width,
                            left_singular, right_singular](Scalar t) {
            if (left_singular && right_singular) {
                const Scalar s = t * t * (Scalar(static_cast<long long>(3)) - Scalar(static_cast<long long>(2)) * t);
                const Scalar dx_dt = width * Scalar(static_cast<long long>(6)) * t * (Scalar(static_cast<long long>(1)) - t);
                return evaluate_with_variable(lower_bound + width * s) * dx_dt;
            }
            if (left_singular) {
                const Scalar x = lower_bound + width * t * t;
                return evaluate_with_variable(x) * Scalar(static_cast<long long>(2)) * width * t;
            }
            const Scalar one_minus_t = Scalar(static_cast<long long>(1)) - t;
            const Scalar x = upper_bound - width * one_minus_t * one_minus_t;
            return evaluate_with_variable(x) * Scalar(static_cast<long long>(2)) * width * one_minus_t;
        };
        reject_divergent_transformed_endpoint(std::function<Scalar(Scalar)>(transformed),
                                              left_singular,
                                              right_singular);
        return adaptive_gauss_kronrod_callable(std::function<Scalar(Scalar)>(transformed),
                                               Scalar(static_cast<long long>(0)),
                                               Scalar(static_cast<long long>(1)),
                                               scaled_eps,
                                               kMaxIntegralDepth);
    }

    Scalar prev_scan_val = evaluate_with_variable(lower_bound);
    for (int i = 1; i <= 40; ++i) {
        const Scalar x =
            lower_bound + (upper_bound - lower_bound) *
                              (Scalar((i)) / Scalar(40.0L));
        Scalar value = evaluate_with_variable(x);
        if (!t_isfinite(value)) {
            throw std::runtime_error("integral did not converge");
        }
        if (t_isfinite(prev_scan_val)) {
            Scalar prev_abs = t_abs(prev_scan_val);
            Scalar curr_abs = t_abs(value);
            bool sign_change = (value > Scalar(0)) != (prev_scan_val > Scalar(0));
            if (sign_change && (prev_abs > Scalar(1e6L) || curr_abs > Scalar(1e6L))) {
                throw std::runtime_error("potential internal discontinuity detected");
            }
        }
        prev_scan_val = value;
    }

    const int coarse_scan_points = 100;
    std::vector<std::pair<Scalar, Scalar>> suspicious_points;
    Scalar prev_x = lower_bound;
    prev_scan_val = evaluate_with_variable(lower_bound);
    for (int i = 1; i <= coarse_scan_points; ++i) {
        const Scalar x = lower_bound + (upper_bound - lower_bound) * Scalar((i)) / Scalar((coarse_scan_points));
        Scalar val;
        try {
            val = evaluate_with_variable(x);
        } catch (...) {
            suspicious_points.push_back({prev_x, x});
            prev_x = x;
            prev_scan_val = Scalar(static_cast<long long>(0));
            continue;
        }

        if (!t_isfinite(val)) {
            suspicious_points.push_back({prev_x, x});
        } else if (t_abs(val) > Scalar(1e9L)) {
            suspicious_points.push_back({prev_x, x});
        } else if (t_isfinite(prev_scan_val)) {
            Scalar jump = t_abs(val - prev_scan_val);
            Scalar avg = (t_abs(val) + t_abs(prev_scan_val)) / Scalar(2.0);
            
            // 检测可能的奇异点：值非常大且发生符号改变，或者跳跃巨大
            bool sign_change = (val > Scalar(0)) != (prev_scan_val > Scalar(0));
            if (sign_change && (t_abs(val) > Scalar(1e6L) || t_abs(prev_scan_val) > Scalar(1e6L))) {
                suspicious_points.push_back({prev_x, x});
            } else if (avg > Scalar(1e-10L) && jump > Scalar(50.0L) * avg) {
                suspicious_points.push_back({prev_x, x});
            }
        }
        prev_x = x;
        prev_scan_val = val;
    }

    for (const auto& susp : suspicious_points) {
        Scalar l = susp.first, r = susp.second;
        for (int iter = 0; iter < 30 && (r - l) > Scalar(1e-12L); ++iter) {
            Scalar mid = (l + r) * Scalar(0.5L);
            Scalar mid_val;
            try {
                mid_val = evaluate_with_variable(mid);
                if (!t_isfinite(mid_val) || t_abs(mid_val) > Scalar(1e10)) {
                    r = mid;
                } else {
                    Scalar left_mid = (l + mid) * Scalar(0.5L);
                    Scalar right_mid = (mid + r) * Scalar(0.5L);
                    Scalar left_val = evaluate_with_variable(left_mid);
                    Scalar right_val = evaluate_with_variable(right_mid);
                    if (!t_isfinite(left_val) || t_abs(left_val) > Scalar(1e10)) {
                        r = mid;
                    } else if (!t_isfinite(right_val) || t_abs(right_val) > Scalar(1e10)) {
                        l = mid;
                    } else {
                        break;
                    }
                }
            } catch (...) {
                r = mid;
            }
        }
        if ((r - l) < Scalar(1e-6)) {
            throw std::runtime_error("potential internal discontinuity detected near x = " +
                                     format_t((l + r) * Scalar(0.5L)));
        }
    }

    return require_finite_integral(
        adaptive_gauss_kronrod(lower_bound,
                               upper_bound,
                               scaled_eps,
                               kMaxIntegralDepth));
}


Scalar FunctionAnalysis::indefinite_integral_at(Scalar x,
                                                Scalar anchor,
                                                Scalar constant) const {
    return constant + definite_integral(anchor, x);
}


std::vector<ExtremumPoint> FunctionAnalysis::solve_extrema(Scalar left_bound,
                                                           Scalar right_bound,
                                                           int scan_segments) const {
    if (left_bound >= right_bound) {
        throw std::runtime_error("extrema search requires left_bound < right_bound");
    }
    if (scan_segments < 8) {
        throw std::runtime_error("scan_segments must be at least 8");
    }

    std::vector<ExtremumPoint> extrema;
    Scalar previous_x = left_bound;
    Scalar previous_derivative = derivative(previous_x);

    for (int i = 1; i <= scan_segments; ++i) {
        const Scalar current_x =
            left_bound +
            (right_bound - left_bound) * Scalar((i)) /
                Scalar((scan_segments));
        const Scalar current_derivative = derivative(current_x);

        if (t_is_near_zero(previous_derivative, Scalar(1e-5))) {
            const Scalar stationary_x = previous_x;
            const Scalar second = second_derivative(stationary_x);
            if (!t_is_near_zero(second, Scalar(1e-4))) {
                bool duplicate = false;
                for (const auto& point : extrema) {
                    if (same_extremum_x(point.x, stationary_x)) {
                        duplicate = true;
                        break;
                    }
                }
                if (!duplicate) {
                    extrema.push_back(
                        {stationary_x, evaluate_with_variable(stationary_x), second < Scalar(static_cast<long long>(0))});
                }
            }
        } else if ((previous_derivative < Scalar(static_cast<long long>(0)) && current_derivative > Scalar(static_cast<long long>(0))) ||
                   (previous_derivative > Scalar(static_cast<long long>(0)) && current_derivative < Scalar(static_cast<long long>(0)))) {
            const Scalar stationary_x =
                bisect_stationary_point(previous_x, current_x);
            const Scalar second = second_derivative(stationary_x);
            if (!t_is_near_zero(second, Scalar(1e-4))) {
                extrema.push_back(
                    {stationary_x, evaluate_with_variable(stationary_x), second < Scalar(static_cast<long long>(0))});
            }
        }

        previous_x = current_x;
        previous_derivative = current_derivative;
    }

    std::vector<ExtremumPoint> unique_extrema;
    for (const auto& point : extrema) {
        bool duplicate = false;
        for (const auto& kept : unique_extrema) {
            if (same_extremum_x(point.x, kept.x)) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            unique_extrema.push_back(point);
        }
    }

    return unique_extrema;
}


const std::string& FunctionAnalysis::expression() const {
    return expression_;
}


const std::string& FunctionAnalysis::variable_name() const {
    return variable_name_;
}


void FunctionAnalysis::ensure_evaluator_initialized() const {
    if (evaluator_) {
        return;
    }
    if (expression_.empty()) {
        throw std::runtime_error("function is not defined");
    }

    fallback_calculator_ = std::make_shared<Calculator>();
    const auto decimal_evaluator =
        fallback_calculator_->get_core_services().evaluation.build_decimal_evaluator(expression_);
    evaluator_ = [decimal_evaluator](
                     const std::vector<std::pair<std::string, Scalar>>& assignments) -> Scalar {
        std::vector<std::pair<std::string, Scalar>> decimal_assignments;
        decimal_assignments.reserve(assignments.size());
        for (const auto& [name, value] : assignments) {
            decimal_assignments.push_back({name, value.to_long_double()});
        }
        return Scalar(decimal_evaluator(decimal_assignments));
    };
}


Scalar FunctionAnalysis::evaluate_with_variable(Scalar x) const {
    if (expression_.empty()) {
        throw std::runtime_error("function is not defined");
    }

    ensure_evaluator_initialized();

    static constexpr std::size_t kMaxEvaluationCacheSize = 4096;
    const std::string cache_key =
        variable_name_ + "|" + expression_ + "|" + format_t(x);
    const auto found = evaluation_cache_index_.find(cache_key);
    if (found != evaluation_cache_index_.end()) {
        evaluation_cache_entries_.splice(evaluation_cache_entries_.begin(),
                                         evaluation_cache_entries_,
                                         found->second);
        return found->second->second;
    }

    Scalar value = Scalar(static_cast<long long>(0));
    value = evaluator_({{variable_name_, x}});
    
    evaluation_cache_entries_.push_front({cache_key, value});
    evaluation_cache_index_[cache_key] = evaluation_cache_entries_.begin();
    while (evaluation_cache_entries_.size() > kMaxEvaluationCacheSize) {
        evaluation_cache_index_.erase(evaluation_cache_entries_.back().first);
        evaluation_cache_entries_.pop_back();
    }
    return value;
}


Scalar FunctionAnalysis::second_derivative(Scalar x) const {
    const Scalar step = scale_aware_step(x);
    const Scalar center = evaluate_with_variable(x);
    const Scalar left_x = x - step;
    const Scalar right_x = x + step;
    const Scalar actual_step = (right_x - left_x) * Scalar(0.5L);
    const Scalar left = evaluate_with_variable(left_x);
    const Scalar right = evaluate_with_variable(right_x);
    const Scalar numerator = compensated_pair_sum(left - center, right - center);
    return numerator / (actual_step * actual_step);
}


Scalar FunctionAnalysis::bisect_stationary_point(Scalar left, Scalar right) const {
    Scalar left_derivative = derivative(left);

    for (int i = 0; i < 80; ++i) {
        const Scalar mid = (left + right) * Scalar(0.5L);
        const Scalar mid_derivative = derivative(mid);
        if (t_abs(mid_derivative) <= kRootTolerance_v() ||
            t_abs(right - left) <=
                relative_tolerance(kRootTolerance_v(),
                                   std::max(t_abs(left), t_abs(right)))) {
            return mid;
        }

        if ((left_derivative < Scalar(static_cast<long long>(0)) && mid_derivative > Scalar(static_cast<long long>(0))) ||
            (left_derivative > Scalar(static_cast<long long>(0)) && mid_derivative < Scalar(static_cast<long long>(0)))) {
            right = mid;
        } else {
            left = mid;
            left_derivative = mid_derivative;
        }
    }

    return (left + right) * Scalar(0.5L);
}

// 自适应 Simpson 积分辅助函数

Scalar FunctionAnalysis::simpson_rule(Scalar a, Scalar b) const {
    const Scalar h = (b - a) / Scalar(static_cast<long long>(2));
    const Scalar fa = evaluate_with_variable(a);
    const Scalar fb = evaluate_with_variable(b);
    const Scalar fc = evaluate_with_variable((a + b) / Scalar(static_cast<long long>(2)));
    return h / Scalar(static_cast<long long>(3)) * (fa + Scalar(static_cast<long long>(4)) * fc + fb);
}


Scalar FunctionAnalysis::adaptive_simpson_precise(Scalar a, Scalar b, Scalar eps, int max_depth) const {
    // 使用自适应 Simpson 方法进行高精度积分
    // 递归实现，当误差足够小或达到最大深度时停止
    int actual_max_depth = max_depth;

    const Scalar c = (a + b) / Scalar(static_cast<long long>(2));
    const Scalar whole = simpson_rule(a, b);
    const Scalar left = simpson_rule(a, c);
    const Scalar right = simpson_rule(c, b);

    return adaptive_simpson_recursive(a, b, whole, left, right, eps, actual_max_depth);
}


Scalar FunctionAnalysis::adaptive_simpson_recursive(Scalar a, Scalar b, Scalar whole, Scalar left, Scalar right, Scalar eps, int depth) const {
    const Scalar c = (a + b) / Scalar(static_cast<long long>(2));
    const Scalar combined = left + right;
    const Scalar error = t_abs(combined - whole) / Scalar(static_cast<long long>(15));

    const Scalar scale = std::max(Scalar(static_cast<long long>(1)), t_abs(combined));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        // Richardson 外推提高精度
        return combined + (combined - whole) / Scalar(static_cast<long long>(15));
    }

    // 继续细分
    const Scalar d = (a + c) / Scalar(static_cast<long long>(2));
    const Scalar e = (c + b) / Scalar(static_cast<long long>(2));
    const Scalar left_left = simpson_rule(a, d);
    const Scalar left_right = simpson_rule(d, c);
    const Scalar right_left = simpson_rule(c, e);
    const Scalar right_right = simpson_rule(e, b);

    return adaptive_simpson_recursive(a, c, left, left_left, left_right, eps / Scalar(static_cast<long long>(2)), depth - 1) +
           adaptive_simpson_recursive(c, b, right, right_left, right_right, eps / Scalar(static_cast<long long>(2)), depth - 1);
}


Scalar FunctionAnalysis::adaptive_gauss_kronrod(Scalar left,
                                                Scalar right,
                                                Scalar eps,
                                                int max_depth) const {
    Scalar error = Scalar(static_cast<long long>(0));
    const Scalar whole = gauss_kronrod_15(left, right, &error);
    return require_finite_integral(
        adaptive_gauss_kronrod_recursive(left,
                                         right,
                                         eps,
                                         whole,
                                         error,
                                         max_depth));
}


Scalar FunctionAnalysis::adaptive_gauss_kronrod_recursive(Scalar left,
                                                          Scalar right,
                                                          Scalar eps,
                                                          Scalar whole,
                                                          Scalar error,
                                                          int depth) const {
    // 检查结果是否有效
    if (!t_isfinite(whole) || !t_isfinite(error)) {
        throw std::runtime_error("integral did not converge (non-finite value encountered)");
    }

    // 检查区间是否过小，避免数值问题
    const Scalar interval_width = t_abs(right - left);
    const Scalar interval_scale = std::max(t_abs(left), t_abs(right));
    const Scalar min_width = precision::min_step_size<Scalar>(interval_scale);
    if (interval_width < min_width) {
        return whole;
    }

    const Scalar scale = std::max(Scalar(static_cast<long long>(1)), t_abs(whole));
    const Scalar tol = relative_tolerance(eps, scale);
    if (error <= tol) {
        return whole;
    }
    if (depth <= 0) {
        if (error > tol * Scalar(1e4L)) { // 严重不收敛
            throw std::runtime_error("integral did not converge (max depth reached with large error)");
        }
        return whole;
    }

    const Scalar mid = (left + right) * Scalar(0.5L);
    Scalar left_error = Scalar(static_cast<long long>(0));
    Scalar right_error = Scalar(static_cast<long long>(0));
    const Scalar left_area = gauss_kronrod_15(left, mid, &left_error);
    const Scalar right_area = gauss_kronrod_15(mid, right, &right_error);

    // 检查子区间结果是否有效
    if (!t_isfinite(left_area) || !t_isfinite(right_area)) {
        throw std::runtime_error("integral did not converge (non-finite value in subinterval)");
    }

    const Scalar left_result =
        adaptive_gauss_kronrod_recursive(left,
                                         mid,
                                         eps * Scalar(0.5L),
                                         left_area,
                                         left_error,
                                         depth - 1);
    const Scalar right_result =
        adaptive_gauss_kronrod_recursive(mid,
                                         right,
                                         eps * Scalar(0.5L),
                                         right_area,
                                         right_error,
                                         depth - 1);
    return compensated_pair_sum(left_result, right_result);
}


Scalar FunctionAnalysis::gauss_kronrod_15(Scalar left,
                                                       Scalar right,
                                                       Scalar* error_estimate) const {
    static const Scalar kNodes[] = {
        Scalar(0.0L),
        Scalar(0.20778495500789846760L),
        Scalar(0.40584515137739716691L),
        Scalar(0.58608723546769113029L),
        Scalar(0.74153118559939443986L),
        Scalar(0.86486442335976907279L),
        Scalar(0.94910791234275852453L),
        Scalar(0.99145537112081263921L)
    };

    static const Scalar kKronrodWeights[] = {
        Scalar(0.20948214108472782801L),
        Scalar(0.20443294007529889241L),
        Scalar(0.19035057806478540991L),
        Scalar(0.16900472663926790283L),
        Scalar(0.14065325971552591875L),
        Scalar(0.10479001032225018384L),
        Scalar(0.06309209262997855329L),
        Scalar(0.02293532201052922496L)
    };

    static const Scalar kGaussWeights[] = {
        Scalar(0.41795918367346938776L),
        Scalar(0.0L),
        Scalar(0.38183005050511894495L),
        Scalar(0.0L),
        Scalar(0.27970539148927666790L),
        Scalar(0.0L),
        Scalar(0.12948496616886969327L),
        Scalar(0.0L)
    };

    const Scalar center = (left + right) * Scalar(0.5L);
    const Scalar half_width = (right - left) * Scalar(0.5L);

    Scalar kronrod_sum = Scalar(0.0L);
    Scalar gauss_sum = Scalar(0.0L);

    for (int i = 0; i < 8; ++i) {
        if (mymath::is_near_zero(kNodes[i], Scalar(0.0L))) {
            Scalar value(evaluate_with_variable(center));
            kronrod_sum = kronrod_sum + kKronrodWeights[i] * value;
            gauss_sum = gauss_sum + kGaussWeights[i] * value;
            continue;
        }

        const Scalar offset = half_width * kNodes[i];
        Scalar left_val(evaluate_with_variable(center - offset));
        Scalar right_val(evaluate_with_variable(center + offset));
        Scalar pair_sum = left_val + right_val;

        kronrod_sum = kronrod_sum + kKronrodWeights[i] * pair_sum;
        if (kGaussWeights[i] != Scalar(0.0L)) {
            gauss_sum = gauss_sum + kGaussWeights[i] * pair_sum;
        }
    }

    Scalar kronrod = half_width * kronrod_sum;
    Scalar gauss = half_width * gauss_sum;

    *error_estimate = mymath::abs(kronrod - gauss);
    return kronrod;
}


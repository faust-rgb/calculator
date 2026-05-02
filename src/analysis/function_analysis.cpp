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

#include "analysis/function_analysis.h"

#include "core/calculator.h"
#include "math/mymath.h"
#include "analysis/calculator_series.h"
#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_expression_internal.h"
#include "statistics/probability.h"

#include <algorithm>
#include <cctype>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace {

/** @brief 数值微分基准步长 */
constexpr double kDerivativeBaseStep = 1e-4;

/** @brief 极限计算的初始步长 */
constexpr double kLimitInitialStep = 1e-1;

/** @brief 极限计算的收敛容差 */
constexpr double kLimitTolerance = 1e-10;

/** @brief 根查找的收敛容差 */
constexpr double kRootTolerance = 1e-7;

/** @brief 数值积分的精度要求 */
constexpr double kIntegralTolerance = 1e-8;

/** @brief 自适应积分的最大递归深度 */
constexpr int kMaxIntegralDepth = 18;

std::string format_double(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    std::string text = out.str();
    if (text.find_first_of("eE") == std::string::npos) {
        while (!text.empty() && text.back() == '0') {
            text.pop_back();
        }
        if (!text.empty() && text.back() == '.') {
            text.pop_back();
        }
    }
    if (text.empty() || text == "-0") {
        return "0";
    }
    return text;
}

long double to_long_double(double value) {
    return static_cast<long double>(value);
}

void compensated_add(long double value,
                     long double* sum,
                     long double* compensation) {
    const long double adjusted = value - *compensation;
    const long double next = *sum + adjusted;
    *compensation = (next - *sum) - adjusted;
    *sum = next;
}

long double compensated_pair_sum(long double lhs, long double rhs) {
    long double sum = 0.0L;
    long double compensation = 0.0L;
    compensated_add(lhs, &sum, &compensation);
    compensated_add(rhs, &sum, &compensation);
    return sum;
}

double scale_aware_step(double x) {
    const double scale = std::max(1.0, mymath::abs(x));
    return kDerivativeBaseStep * scale;
}

double central_difference_step_value(double scale, double factor) {
    return std::max(1e-7 * scale, kDerivativeBaseStep * scale * factor);
}

double relative_tolerance(double baseline, double scale) {
    return baseline * std::max(1.0, scale);
}

double limit_step_scale(double x) {
    return kLimitInitialStep * std::max(1.0, mymath::abs(x));
}

bool same_extremum_x(double lhs, double rhs) {
    return mymath::abs(lhs - rhs) <= 1e-5;
}

double require_finite_integral(double value) {
    if (!mymath::isfinite(value)) {
        throw std::runtime_error("integral did not converge");
    }
    return value;
}

void reject_divergent_transformed_endpoint(
    const std::function<double(double)>& transformed,
    bool check_left,
    bool check_right) {
    const double offsets[] = {1e-3, 1e-4, 1e-5};
    auto check_at = [&](double t) {
        double value = transformed(t);
        if (!mymath::isfinite(value) || mymath::abs(value) > 1e4) {
            throw std::runtime_error("integral did not converge");
        }
    };

    for (double offset : offsets) {
        if (check_left) {
            check_at(offset);
        }
        if (check_right) {
            check_at(1.0 - offset);
        }
    }
}

double gauss_kronrod_15_callable(const std::function<double(double)>& function,
                                 double left,
                                 double right,
                                 double* error_estimate) {
    static const double kNodes[] = {
        0.9914553711208126,
        0.9491079123427585,
        0.8648644233597691,
        0.7415311855993945,
        0.5860872354676911,
        0.4058451513773972,
        0.2077849550078985,
        0.0,
    };
    static const double kKronrodWeights[] = {
        0.02293532201052922,
        0.06309209262997855,
        0.1047900103222502,
        0.1406532597155259,
        0.1690047266392679,
        0.1903505780647854,
        0.2044329400752989,
        0.2094821410847278,
    };
    static const double kGaussWeights[] = {
        0.0,
        0.1294849661688697,
        0.0,
        0.2797053914892767,
        0.0,
        0.3818300505051189,
        0.0,
        0.4179591836734694,
    };

    const long double center =
        (to_long_double(left) + to_long_double(right)) * 0.5L;
    const long double half_width =
        (to_long_double(right) - to_long_double(left)) * 0.5L;
    long double kronrod_sum = 0.0L;
    long double gauss_sum = 0.0L;
    long double kronrod_compensation = 0.0L;
    long double gauss_compensation = 0.0L;

    for (int i = 0; i < 8; ++i) {
        if (mymath::is_near_zero(kNodes[i], 0.0)) {
            const long double value =
                to_long_double(function(static_cast<double>(center)));
            compensated_add(static_cast<long double>(kKronrodWeights[i]) * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(static_cast<long double>(kGaussWeights[i]) * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const long double offset = half_width * static_cast<long double>(kNodes[i]);
        const long double left_value =
            to_long_double(function(static_cast<double>(center - offset)));
        const long double right_value =
            to_long_double(function(static_cast<double>(center + offset)));
        const long double pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(static_cast<long double>(kKronrodWeights[i]) * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(static_cast<long double>(kGaussWeights[i]) * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const long double kronrod = half_width * kronrod_sum;
    const long double gauss = half_width * gauss_sum;
    *error_estimate = static_cast<double>(mymath::abs_long_double(kronrod - gauss));
    return static_cast<double>(kronrod);
}

double adaptive_gauss_kronrod_callable_recursive(
    const std::function<double(double)>& function,
    double left,
    double right,
    double eps,
    double whole,
    double error,
    int depth) {
    const double scale = std::max(1.0, mymath::abs(whole));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        return whole;
    }

    const double mid = (left + right) * 0.5;
    double left_error = 0.0;
    double right_error = 0.0;
    const double left_area =
        gauss_kronrod_15_callable(function, left, mid, &left_error);
    const double right_area =
        gauss_kronrod_15_callable(function, mid, right, &right_error);
    const long double left_result = to_long_double(
        adaptive_gauss_kronrod_callable_recursive(function,
                                                  left,
                                                  mid,
                                                  eps * 0.5,
                                                  left_area,
                                                  left_error,
                                                  depth - 1));
    const long double right_result = to_long_double(
        adaptive_gauss_kronrod_callable_recursive(function,
                                                  mid,
                                                  right,
                                                  eps * 0.5,
                                                  right_area,
                                                  right_error,
                                                  depth - 1));
    return static_cast<double>(compensated_pair_sum(left_result, right_result));
}

double adaptive_gauss_kronrod_callable(const std::function<double(double)>& function,
                                       double left,
                                       double right,
                                       double eps,
                                       int depth) {
    double error = 0.0;
    const double whole = gauss_kronrod_15_callable(function, left, right, &error);
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

bool is_infinite_probe(SymbolicLimitProbeKind kind) {
    return kind == SymbolicLimitProbeKind::kPositiveInfinity ||
           kind == SymbolicLimitProbeKind::kNegativeInfinity;
}

SymbolicLimitProbeKind probe_symbolic_value_at(
    SymbolicExpression expression,
    const std::string& variable_name,
    double point,
    double* finite_value) {
    try {
        if (!mymath::isfinite(point)) {
            expression = expression.simplify();
            const std::shared_ptr<SymbolicExpression::Node>& node = expression.node_;
            switch (node->type) {
                case NodeType::kNumber:
                case NodeType::kPi:
                case NodeType::kE:
                case NodeType::kInfinity: {
                    double value = 0.0;
                    if (expression.is_number(&value)) {
                        *finite_value = value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kVariable:
                    if (node->text == variable_name) {
                        return point > 0.0 ? SymbolicLimitProbeKind::kPositiveInfinity
                                           : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                case NodeType::kNegate: {
                    double child_value = 0.0;
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
                    double left_value = 0.0;
                    double right_value = 0.0;
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
                    double left_value = 0.0;
                    double right_value = 0.0;
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
                        mymath::is_near_zero(left_value, kLimitTolerance)) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    if (right_kind == SymbolicLimitProbeKind::kFinite &&
                        mymath::is_near_zero(right_value, kLimitTolerance)) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    auto sign_of = [](SymbolicLimitProbeKind kind, double value) {
                        if (kind == SymbolicLimitProbeKind::kFinite) {
                            return value >= 0.0 ? 1 : -1;
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
                    double left_value = 0.0;
                    double right_value = 0.0;
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
                        !mymath::is_near_zero(right_value, kLimitTolerance)) {
                        *finite_value = left_value / right_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        is_infinite_probe(right_kind)) {
                        *finite_value = 0.0;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (is_infinite_probe(left_kind) &&
                        right_kind == SymbolicLimitProbeKind::kFinite &&
                        !mymath::is_near_zero(right_value, kLimitTolerance)) {
                        const bool positive =
                            (left_kind == SymbolicLimitProbeKind::kPositiveInfinity) ==
                            (right_value > 0.0);
                        return positive ? SymbolicLimitProbeKind::kPositiveInfinity
                                        : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kPower: {
                    const SymbolicExpression base(node->left);
                    const SymbolicExpression exponent(node->right);
                    double exponent_value = 0.0;
                    if (!exponent.is_number(&exponent_value)) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    double base_value = 0.0;
                    const SymbolicLimitProbeKind base_kind =
                        probe_symbolic_value_at(base,
                                                variable_name,
                                                point,
                                                &base_value);
                    if (base_kind == SymbolicLimitProbeKind::kFinite) {
                        *finite_value = mymath::pow(base_value, exponent_value);
                        return mymath::isfinite(*finite_value)
                                   ? SymbolicLimitProbeKind::kFinite
                                   : SymbolicLimitProbeKind::kUnknown;
                    }
                    if (is_infinite_probe(base_kind)) {
                        if (exponent_value > 0.0) {
                            if (base_kind == SymbolicLimitProbeKind::kNegativeInfinity &&
                                mymath::is_integer(exponent_value) &&
                                static_cast<long long>(std::llround(exponent_value)) % 2 != 0) {
                                return SymbolicLimitProbeKind::kNegativeInfinity;
                            }
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (exponent_value < 0.0) {
                            *finite_value = 0.0;
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kFunction: {
                    const std::string& name = node->text;
                    double argument_value = 0.0;
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
                            argument_value > 0.0) {
                            *finite_value = mymath::log(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    if (name == "exp") {
                        if (argument_kind == SymbolicLimitProbeKind::kPositiveInfinity) {
                            return SymbolicLimitProbeKind::kPositiveInfinity;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kNegativeInfinity) {
                            *finite_value = 0.0;
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kFinite) {
                            *finite_value = mymath::exp(argument_value);
                            return SymbolicLimitProbeKind::kFinite;
                        }
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kVector:
                case NodeType::kTensor:
                case NodeType::kDifferentialOp:
                    return SymbolicLimitProbeKind::kUnknown;
            }
            return SymbolicLimitProbeKind::kUnknown;
        }

        expression = expression.substitute(
            variable_name,
            SymbolicExpression::number(point)).simplify();
        double value = 0.0;
        if (!expression.is_number(&value)) {
            return SymbolicLimitProbeKind::kUnknown;
        }
        if (mymath::isfinite(value)) {
            *finite_value = value;
            return SymbolicLimitProbeKind::kFinite;
        }
        return value > 0.0 ? SymbolicLimitProbeKind::kPositiveInfinity
                           : SymbolicLimitProbeKind::kNegativeInfinity;
    } catch (...) {
        return SymbolicLimitProbeKind::kUnknown;
    }
}

bool is_zero_probe(SymbolicLimitProbeKind kind, double value) {
    return kind == SymbolicLimitProbeKind::kFinite &&
           mymath::is_near_zero(value, kLimitTolerance);
}

bool try_symbolic_lhopital_limit(const SymbolicExpression& expression,
                                 const std::string& variable_name,
                                 double point,
                                 double* result) {
    SymbolicExpression current = expression.simplify();
    static constexpr int kMaxLhopitalDepth = 8;
    for (int depth = 0; depth < kMaxLhopitalDepth; ++depth) {
        if (current.node_->type != NodeType::kDivide) {
            return false;
        }

        SymbolicExpression numerator(current.node_->left);
        SymbolicExpression denominator(current.node_->right);
        double numerator_value = 0.0;
        double denominator_value = 0.0;
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

        current = (numerator.derivative(variable_name).simplify() /
                   denominator.derivative(variable_name).simplify()).simplify();

        double current_value = 0.0;
        const SymbolicLimitProbeKind current_kind =
            probe_symbolic_value_at(current,
                                    variable_name,
                                    point,
                                    &current_value);
        if (current_kind == SymbolicLimitProbeKind::kFinite) {
            *result = current_value;
            return true;
        }
    }
    return false;
}

/**
 * @brief 符号极限在无穷点
 *
 * 计算 limit(expr, var, ±inf) 的符号形式。
 * 使用代换 x = 1/t，将问题转化为 limit(expr(1/t), t, 0)。
 *
 * @param expression 符号表达式
 * @param variable_name 变量名
 * @param positive true 表示 +inf，false 表示 -inf
 * @param result 输出极限值
 * @return true 如果成功计算符号极限
 */
bool symbolic_limit_at_infinity(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                bool positive,
                                double* result) {
    series_ops::SeriesContext ctx;
    ctx.evaluate_at = [](const SymbolicExpression& e, const std::string& /*v*/, double /*p*/) {
        double val = 0.0;
        if (e.is_number(&val)) return val;
        return 0.0;
    };

    // 策略：对于 x -> inf，代换 x = 1/t，然后在 t -> 0 处展开
    // 注意：如果代换后的表达式在 t=0 处有极点，需要区分两种情况：
    // 1. 原表达式趋向无穷大（如 x^2 -> 1/t^2）
    // 2. 代换引入的极点（如 1/x -> 1/(1/t)），原极限是有限值

    SymbolicExpression t_var = SymbolicExpression::variable("t_limit_inf_tmp");
    SymbolicExpression inv_t;
    if (positive) {
        inv_t = SymbolicExpression::number(1.0) / t_var;
    } else {
        inv_t = SymbolicExpression::number(-1.0) / t_var;
    }
    SymbolicExpression substituted = expression.substitute(variable_name, inv_t).simplify();

    std::vector<double> coeffs;
    try {
        if (series_ops::internal::evaluate_psa(substituted, "t_limit_inf_tmp", 0.0, 2, coeffs, ctx)) {
            if (!coeffs.empty() && mymath::isfinite(coeffs[0])) {
                *result = coeffs[0];
                return true;
            }
        }
    } catch (const series_ops::internal::PoleException& e) {
        // 代换后的表达式在 t=0 处有极点
        // 检查原表达式的主导项来确定极限行为
        // 对于多项式类表达式，最高次项决定行为
        // 对于有理函数，比较分子分母的最高次

        // 尝试在非常大的点处展开，获取渐近行为
        double large_x = positive ? 1e6 : -1e6;
        try {
            std::vector<double> large_coeffs;
            if (series_ops::internal::evaluate_psa(expression, variable_name, large_x, 1, large_coeffs, ctx)) {
                if (!large_coeffs.empty()) {
                    // 使用线性近似估计极限趋势
                    double val_at_large = large_coeffs[0];
                    double slope = (large_coeffs.size() > 1) ? large_coeffs[1] : 0.0;

                    // 如果值已经很大且在增长，极限是无穷
                    if (mymath::abs(val_at_large) > 1e10) {
                        *result = (val_at_large > 0) ? mymath::infinity() : -mymath::infinity();
                        return true;
                    }

                    // 如果值很小且斜率也很小，极限可能是有限值
                    if (mymath::abs(val_at_large) < 1e-6 && mymath::abs(slope) < 1e-6) {
                        *result = val_at_large;
                        return true;
                    }
                }
            }
        } catch (...) {
            // 忽略，回退到其他方法
        }

        // 如果上述方法失败，返回 false 让数值方法处理
        return false;
    }

    // 尝试 L'Hopital 规则（对于无穷点）
    double lhopital_result = 0.0;
    if (try_symbolic_lhopital_limit(expression, variable_name,
                                     positive ? mymath::infinity()
                                              : -mymath::infinity(),
                                     &lhopital_result)) {
        *result = lhopital_result;
        return true;
    }

    return false;
}

/**
 * @brief 处理极点极限
 *
 * 根据 Laurent 级数的位移和前导系数判定无穷极限。
 *
 * @param shift Laurent 位移（负数表示极点阶数）
 * @param leading_coefficient 前导系数
 * @param direction 方向：-1 左极限，1 右极限，0 双侧极限
 * @return 极限值（+inf 或 -inf）
 * @throw std::runtime_error 当双侧极限不存在时
 */
double handle_pole_limit(int shift, double leading_coefficient, int direction) {
    if (direction == 0) {
        // 双侧极限：只有当 shift 为偶数时才存在
        if (shift % 2 == 0) {
            return (leading_coefficient > 0) ? mymath::infinity() : -mymath::infinity();
        } else {
            throw std::runtime_error("two-sided limit does not exist (pole with odd shift)");
        }
    } else if (direction == 1) {
        // 右极限：(x - x0) > 0，符号不变
        return (leading_coefficient > 0) ? mymath::infinity() : -mymath::infinity();
    } else {
        // 左极限：(x - x0) < 0，奇数 shift 时符号翻转
        bool flip_sign = (shift % 2 != 0);
        double effective_c = flip_sign ? -leading_coefficient : leading_coefficient;
        return (effective_c > 0) ? mymath::infinity() : -mymath::infinity();
    }
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
      evaluator_(other.evaluator_) {}

FunctionAnalysis& FunctionAnalysis::operator=(const FunctionAnalysis& other) {
    if (this != &other) {
        expression_ = other.expression_;
        variable_name_ = other.variable_name_;
        evaluator_ = other.evaluator_;
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
    evaluation_cache_entries_.clear();
    evaluation_cache_index_.clear();
}

void FunctionAnalysis::set_evaluator(std::function<double(const std::vector<std::pair<std::string, double>>&)> evaluator) {
    evaluator_ = std::move(evaluator);
}

double FunctionAnalysis::evaluate(double x) const {
    return evaluate_with_variable(x);
}

double FunctionAnalysis::derivative(double x) const {
    const double scale = std::max(1.0, mymath::abs(x));
    const double center = evaluate_with_variable(x);
    if (!mymath::isfinite(center)) {
        throw std::runtime_error("derivative is undefined at this point");
    }
    const double curvature_probe = evaluate_with_variable(x + scale * 1e-3) -
                                   2.0 * center +
                                   evaluate_with_variable(x - scale * 1e-3);
    const double curvature_scale =
        std::max(1.0, mymath::abs(curvature_probe) / std::max(1e-12, mymath::abs(center)));
    const double base_step =
        central_difference_step_value(scale, 1.0 / mymath::pow(curvature_scale, 0.25));

    long double richardson[4][4] = {};
    bool row_valid[4] = {};
    long double best_value = 0.0L;
    long double best_error = static_cast<long double>(mymath::infinity());
    for (int row = 0; row < 4; ++row) {
        const double step = base_step / mymath::pow(2.0, static_cast<double>(row));
        const long double forward_x = to_long_double(x) + to_long_double(step);
        const long double backward_x = to_long_double(x) - to_long_double(step);
        const long double actual_step =
            (forward_x - backward_x) * 0.5L;
        if (actual_step <= 0.0L) {
            continue;
        }
        const long double forward =
            to_long_double(evaluate_with_variable(static_cast<double>(forward_x)));
        const long double backward =
            to_long_double(evaluate_with_variable(static_cast<double>(backward_x)));
        if (!mymath::isfinite(static_cast<double>(forward)) ||
            !mymath::isfinite(static_cast<double>(backward))) {
            continue;
        }
        richardson[row][0] = (forward - backward) / (2.0L * actual_step);
        row_valid[row] = mymath::isfinite(static_cast<double>(richardson[row][0]));
        if (!row_valid[row]) {
            continue;
        }
        for (int col = 1; col <= row; ++col) {
            if (!row_valid[row - 1]) {
                row_valid[row] = false;
                break;
            }
            const long double factor =
                static_cast<long double>(mymath::pow(4.0, static_cast<double>(col)));
            richardson[row][col] =
                richardson[row][col - 1] +
                (richardson[row][col - 1] - richardson[row - 1][col - 1]) /
                    (factor - 1.0L);
            if (!mymath::isfinite(static_cast<double>(richardson[row][col]))) {
                row_valid[row] = false;
                break;
            }
        }
        if (row > 0 && row_valid[row] && row_valid[row - 1]) {
            const long double candidate = richardson[row][row];
            const long double error_estimate =
                mymath::abs_long_double(candidate - richardson[row - 1][row - 1]);
            if (error_estimate < best_error &&
                mymath::isfinite(static_cast<double>(candidate))) {
                best_error = error_estimate;
                best_value = candidate;
            }
        }
    }

    if (best_error < static_cast<long double>(mymath::infinity())) {
        const double side_step = std::max(1e-7 * scale, base_step / 64.0);
        const double left_value = evaluate_with_variable(x - side_step);
        const double right_value = evaluate_with_variable(x + side_step);
        if (!mymath::isfinite(left_value) || !mymath::isfinite(right_value)) {
            throw std::runtime_error("derivative is undefined at this point");
        }
        const double left_slope = (center - left_value) / side_step;
        const double right_slope = (right_value - center) / side_step;
        const double slope_scale =
            std::max({1.0, mymath::abs(left_slope), mymath::abs(right_slope),
                      mymath::abs(static_cast<double>(best_value))});
        if (mymath::abs(left_slope - right_slope) >
            std::max(1e-4, 1e-5 * slope_scale)) {
            throw std::runtime_error("derivative does not exist at this point");
        }
        return static_cast<double>(best_value);
    }
    for (int row = 3; row >= 0; --row) {
        if (row_valid[row]) {
            return static_cast<double>(richardson[row][row]);
        }
    }
    return static_cast<double>(richardson[3][3]);
}

double FunctionAnalysis::limit(double x, int direction) const {
    if (direction != -1 && direction != 0 && direction != 1) {
        throw std::runtime_error("limit direction must be -1, 0, or 1");
    }

    // --- 策略优先级：PSA > 洛必达 > 数值方法 ---
    // PSA 是主力方法，对于大多数情况能在微秒级完成
    // 洛必达作为备选，处理 PSA 无法处理的情况
    // 数值方法是最后手段

    SymbolicExpression expr;
    try {
        expr = SymbolicExpression::parse(expression_);
    } catch (...) {
        // 解析失败，直接使用数值方法
        return compute_numerical_limit(x, direction);
    }

    series_ops::SeriesContext ctx;
    ctx.evaluate_at = [this](const SymbolicExpression& e, const std::string& v, double p) {
        if (v == variable_name_) return p;
        double val = 0.0;
        if (e.is_number(&val)) return val;
        return 0.0;
    };

    // --- 第一优先级：PSA（幂级数展开）---
    if (mymath::isfinite(x)) {
        std::vector<double> coeffs;
        try {
            if (series_ops::internal::evaluate_psa(expr, variable_name_, x, 2, coeffs, ctx)) {
                if (!coeffs.empty()) return coeffs[0];
            }
        } catch (const series_ops::internal::PoleException& e) {
            // PSA 检测到极点，直接利用 Laurent 级数信息判定无穷极限
            return handle_pole_limit(e.shift, e.leading_coefficient, direction);
        }
    } else {
        // x 为无穷大：使用封装的 symbolic_limit_at_infinity
        bool positive = x > 0;
        double inf_result = 0.0;
        if (symbolic_limit_at_infinity(expr, variable_name_, positive, &inf_result)) {
            return inf_result;
        }
    }

    // --- 第二优先级：洛必达法则（仅对 0/0 或 inf/inf 型）---
    // 只有当 PSA 失败时才尝试洛必达
    double lhopital_value = 0.0;
    if (direction == 0 &&
        try_symbolic_lhopital_limit(expr,
                                    variable_name_,
                                    x,
                                    &lhopital_value)) {
        return lhopital_value;
    }

    // --- 第三优先级：数值方法（Richardson 外推）---
    return compute_numerical_limit(x, direction);
}

double FunctionAnalysis::compute_numerical_limit(double x, int direction) const {
    // 自适应步长的 Richardson 外推算法
    // 特点：
    // 1. 自适应步长控制：根据误差估计调整步长
    // 2. 步长拒绝机制：当检测到剧烈变化时拒绝当前步长
    // 3. 振荡检测：检测高频振荡函数并调整策略

    auto compute_limit_at = [this](double x_target, int side) {
        // side: -1 (左), 1 (右), 0 (中心或无穷大)
        long double richardson[14][14] = {};
        bool row_valid[14] = {};
        long double best_value = 0.0L;
        long double best_error = static_cast<long double>(mymath::infinity());
        bool have_best = false;

        // 自适应步长参数
        const double base_h = mymath::isfinite(x_target) ? limit_step_scale(x_target) : 1e-2;
        double adaptive_h = base_h;  // 当前自适应步长
        int consecutive_bad = 0;     // 连续不良采样计数
        constexpr int kMaxBadSamples = 3;  // 最大连续不良采样数

        // 用于振荡检测的历史值
        long double prev_val = 0.0L;
        bool have_prev = false;
        int oscillation_count = 0;

        for (int row = 0; row < 14; ++row) {
            // 使用自适应步长
            const double h = adaptive_h / mymath::pow(2.0, static_cast<double>(row + 4));
            double sample_x;
            if (mymath::isfinite(x_target)) {
                sample_x = x_target + static_cast<double>(side) * h;
            } else {
                // 无穷大：x = 1/h (或 -1/h)
                sample_x = (x_target > 0 ? 1.0 : -1.0) / h;
            }

            long double val = 0.0L;
            try {
                val = to_long_double(evaluate_with_variable(sample_x));
            } catch (...) {
                // 采样失败，调整步长并重试
                adaptive_h *= 0.5;
                consecutive_bad++;
                if (consecutive_bad >= kMaxBadSamples) {
                    throw std::runtime_error("limit did not converge (sampling failures)");
                }
                continue;
            }

            // 检测无穷大值（可能是极点）
            if (!mymath::isfinite(static_cast<double>(val))) {
                // 如果值趋向无穷，可能是无限极限
                if (have_prev && mymath::isfinite(static_cast<double>(prev_val))) {
                    // 从有限变为无穷，检查趋势
                    if (prev_val > 1e10) {
                        return mymath::infinity();
                    } else if (prev_val < -1e10) {
                        return -mymath::infinity();
                    }
                }
                adaptive_h *= 0.5;
                consecutive_bad++;
                if (consecutive_bad >= kMaxBadSamples) {
                    // 多次采样都得到无穷大，可能是极点
                    throw std::runtime_error("limit appears to be infinite (numerical evidence)");
                }
                continue;
            }

            // 振荡检测：检测符号频繁变化
            if (have_prev) {
                if ((val > 0 && prev_val < 0) || (val < 0 && prev_val > 0)) {
                    oscillation_count++;
                    if (oscillation_count >= 4) {
                        // 可能是高频振荡函数，使用更小的步长
                        adaptive_h *= 0.25;
                        oscillation_count = 0;
                    }
                } else {
                    oscillation_count = std::max(0, oscillation_count - 1);
                }
            }
            prev_val = val;
            have_prev = true;

            // 步长拒绝机制：检测剧烈变化
            if (have_best && row > 0) {
                const long double expected_change = best_error * 10.0L + 1e-10L;
                const long double actual_change = mymath::abs_long_double(val - best_value);
                if (actual_change > expected_change * 1e6L) {
                    // 变化过于剧烈，拒绝此步长
                    adaptive_h *= 0.5;
                    consecutive_bad++;
                    if (consecutive_bad >= kMaxBadSamples) {
                        break;  // 尝试用已有的 Richardson 表外推
                    }
                    continue;
                }
            }

            consecutive_bad = 0;  // 重置不良采样计数
            richardson[row][0] = val;
            row_valid[row] = true;

            for (int col = 1; col <= row; ++col) {
                if (!row_valid[row - 1]) {
                    row_valid[row] = false;
                    break;
                }
                // Richardson 外推公式
                const long double factor = static_cast<long double>(1LL << col);
                richardson[row][col] = richardson[row][col-1] +
                    (richardson[row][col-1] - richardson[row-1][col-1]) / (factor - 1.0L);
                if (!mymath::isfinite(static_cast<double>(richardson[row][col]))) {
                    row_valid[row] = false;
                    break;
                }
            }

            if (row > 0 && row_valid[row]) {
                long double err = mymath::abs_long_double(richardson[row][row] - richardson[row-1][row-1]);
                if (err < best_error) {
                    best_error = err;
                    best_value = richardson[row][row];
                    have_best = true;
                }

                // 自适应步长调整：误差下降慢时减小步长
                if (err > best_error * 0.9L && row > 3) {
                    adaptive_h *= 0.75;
                }

                // 如果误差已经非常小，提前退出
                if (err < 1e-15L) {
                    return static_cast<double>(richardson[row][row]);
                }
            } else if (row == 0 && row_valid[0]) {
                best_value = richardson[0][0];
            }
        }

        if (!have_best) {
            throw std::runtime_error("limit did not converge");
        }

        const long double scale =
            std::max(1.0L, mymath::abs_long_double(best_value));
        const long double acceptable_error =
            std::max(1e-8L, static_cast<long double>(kLimitTolerance) * 1000.0L * scale);
        if (best_error > acceptable_error) {
            throw std::runtime_error("limit did not converge");
        }
        return static_cast<double>(best_value);
    };

    if (!mymath::isfinite(x)) {
        return compute_limit_at(x, 0);
    }

    if (direction == -1) return compute_limit_at(x, -1);
    if (direction == 1) return compute_limit_at(x, 1);

    // 双侧极限
    const double left = compute_limit_at(x, -1);
    const double right = compute_limit_at(x, 1);

    // 如果两侧极限非常接近，或者是符号一致的无穷大
    if (mymath::abs(left - right) <= kLimitTolerance * 100.0 ||
        (!mymath::isfinite(left) && !mymath::isfinite(right) && ((left > 0) == (right > 0)))) {
        return (left + right) * 0.5;
    }

    throw std::runtime_error("two-sided limit does not exist");
}

double FunctionAnalysis::definite_integral(double lower_bound,
                                           double upper_bound) const {
    if (mymath::is_near_zero(lower_bound - upper_bound, 1e-15)) {
        return 0.0;
    }
    if (lower_bound > upper_bound) {
        return -definite_integral(upper_bound, lower_bound);
    }
    const bool lower_is_infinite = !mymath::isfinite(lower_bound);
    const bool upper_is_infinite = !mymath::isfinite(upper_bound);
    if (lower_is_infinite || upper_is_infinite) {
        if (lower_is_infinite && upper_is_infinite) {
            if (lower_bound > 0.0 || upper_bound < 0.0) {
                throw std::runtime_error("invalid infinite integration bounds");
            }
            auto transformed = [this](double t) {
                const double angle = mymath::kPi * (t - 0.5);
                const double cos_angle = mymath::cos(angle);
                const double x = mymath::tan(angle);
                return evaluate_with_variable(x) * mymath::kPi /
                       (cos_angle * cos_angle);
            };
            reject_divergent_transformed_endpoint(transformed, true, true);
            return adaptive_gauss_kronrod_callable(transformed,
                                                   0.0,
                                                   1.0,
                                                   kIntegralTolerance,
                                                   kMaxIntegralDepth);
        }

        if (lower_is_infinite) {
            if (lower_bound > 0.0) {
                throw std::runtime_error("invalid infinite integration bounds");
            }
            auto transformed = [this, upper_bound](double t) {
                const double x = upper_bound - (1.0 - t) / t;
                return evaluate_with_variable(x) / (t * t);
            };
            reject_divergent_transformed_endpoint(transformed, true, false);
            return adaptive_gauss_kronrod_callable(transformed,
                                                   0.0,
                                                   1.0,
                                                   kIntegralTolerance,
                                                   kMaxIntegralDepth);
        }

        if (upper_bound < 0.0) {
            throw std::runtime_error("invalid infinite integration bounds");
        }
        auto transformed = [this, lower_bound](double t) {
            const double one_minus_t = 1.0 - t;
            const double x = lower_bound + t / one_minus_t;
            return evaluate_with_variable(x) / (one_minus_t * one_minus_t);
        };
        reject_divergent_transformed_endpoint(transformed, false, true);
        return adaptive_gauss_kronrod_callable(transformed,
                                               0.0,
                                               1.0,
                                               kIntegralTolerance,
                                               kMaxIntegralDepth);
    }
    const double span = mymath::abs(upper_bound - lower_bound);
    const double scaled_eps =
        relative_tolerance(kIntegralTolerance, span + mymath::abs(lower_bound) + mymath::abs(upper_bound));
    bool left_singular = false;
    bool right_singular = false;
    try {
        left_singular = !mymath::isfinite(evaluate_with_variable(lower_bound));
    } catch (const std::exception&) {
        left_singular = true;
    }
    try {
        right_singular = !mymath::isfinite(evaluate_with_variable(upper_bound));
    } catch (const std::exception&) {
        right_singular = true;
    }

    if (left_singular || right_singular) {
        const double width = upper_bound - lower_bound;
        auto transformed = [this, lower_bound, upper_bound, width,
                            left_singular, right_singular](double t) {
            if (left_singular && right_singular) {
                const double s = t * t * (3.0 - 2.0 * t);
                const double dx_dt = width * 6.0 * t * (1.0 - t);
                return evaluate_with_variable(lower_bound + width * s) * dx_dt;
            }
            if (left_singular) {
                const double x = lower_bound + width * t * t;
                return evaluate_with_variable(x) * 2.0 * width * t;
            }
            const double one_minus_t = 1.0 - t;
            const double x = upper_bound - width * one_minus_t * one_minus_t;
            return evaluate_with_variable(x) * 2.0 * width * one_minus_t;
        };
        reject_divergent_transformed_endpoint(transformed,
                                              left_singular,
                                              right_singular);
        return adaptive_gauss_kronrod_callable(transformed,
                                               0.0,
                                               1.0,
                                               scaled_eps,
                                               kMaxIntegralDepth);
    }

    // Check for internal singularities by looking for large values or sign changes
    // that indicate a pole (like tan(x) near pi/2)
    double prev_scan_val = evaluate_with_variable(lower_bound);
    for (int i = 1; i <= 40; ++i) {
        const double x =
            lower_bound + (upper_bound - lower_bound) *
                              (static_cast<double>(i) / 40.0);
        double value = evaluate_with_variable(x);
        if (!mymath::isfinite(value)) {
            throw std::runtime_error("integral did not converge");
        }
        // Check for sign change through large values (characteristic of a pole)
        if (mymath::isfinite(prev_scan_val)) {
            double prev_abs = mymath::abs(prev_scan_val);
            double curr_abs = mymath::abs(value);
            // If values have opposite signs and at least one is large, likely a pole
            if ((prev_scan_val > 0) != (value > 0) &&
                (prev_abs > 100.0 || curr_abs > 100.0)) {
                throw std::runtime_error("potential internal discontinuity detected");
            }
        }
        prev_scan_val = value;
    }

    // Additional adaptive scan for internal discontinuities
    // Use a more thorough scan with adaptive refinement near suspicious points
    const int coarse_scan_points = 100;
    std::vector<std::pair<double, double>> suspicious_points;
    double prev_x = lower_bound;
    double prev_val = evaluate_with_variable(lower_bound);
    for (int i = 1; i <= coarse_scan_points; ++i) {
        const double x = lower_bound + (upper_bound - lower_bound) * static_cast<double>(i) / coarse_scan_points;
        double val;
        try {
            val = evaluate_with_variable(x);
        } catch (...) {
            suspicious_points.push_back({prev_x, x});
            prev_x = x;
            prev_val = 0.0;
            continue;
        }

        // Check for discontinuity indicators
        if (!mymath::isfinite(val)) {
            suspicious_points.push_back({prev_x, x});
        } else if (mymath::abs(val) > 1e6) {  // Lower threshold to catch poles
            suspicious_points.push_back({prev_x, x});
        } else if (mymath::isfinite(prev_val)) {
            // Check for large jump (potential discontinuity)
            double jump = mymath::abs(val - prev_val);
            double avg = (mymath::abs(val) + mymath::abs(prev_val)) / 2.0;
            // Lower threshold to catch sign changes through infinity (like tan near pi/2)
            if (avg > 1e-10 && jump > 10.0 * avg) {
                suspicious_points.push_back({prev_x, x});
            }
        }
        prev_x = x;
        prev_val = val;
    }

    // Refine suspicious regions with binary search
    for (const auto& [left, right] : suspicious_points) {
        // Binary search for the discontinuity point
        double l = left, r = right;
        for (int iter = 0; iter < 30 && (r - l) > 1e-12; ++iter) {
            double mid = (l + r) / 2.0;
            double mid_val;
            try {
                mid_val = evaluate_with_variable(mid);
                if (!mymath::isfinite(mid_val) || mymath::abs(mid_val) > 1e10) {
                    r = mid;
                } else {
                    // Check if left or right side has the issue
                    double left_mid = (l + mid) / 2.0;
                    double right_mid = (mid + r) / 2.0;
                    double left_val = evaluate_with_variable(left_mid);
                    double right_val = evaluate_with_variable(right_mid);
                    if (!mymath::isfinite(left_val) || mymath::abs(left_val) > 1e10) {
                        r = mid;
                    } else if (!mymath::isfinite(right_val) || mymath::abs(right_val) > 1e10) {
                        l = mid;
                    } else {
                        break;  // No clear discontinuity found
                    }
                }
            } catch (...) {
                r = mid;
            }
        }
        // If we found a narrow suspicious region, throw an error
        if ((r - l) < 1e-6) {
            throw std::runtime_error("potential internal discontinuity detected near x = " +
                                     std::to_string((l + r) / 2.0));
        }
    }

    return require_finite_integral(
        adaptive_gauss_kronrod(lower_bound,
                               upper_bound,
                               scaled_eps,
                               kMaxIntegralDepth));
}

double FunctionAnalysis::indefinite_integral_at(double x,
                                                double anchor,
                                                double constant) const {
    return constant + definite_integral(anchor, x);
}

std::vector<ExtremumPoint> FunctionAnalysis::solve_extrema(double left_bound,
                                                           double right_bound,
                                                           int scan_segments) const {
    if (left_bound >= right_bound) {
        throw std::runtime_error("extrema search requires left_bound < right_bound");
    }
    if (scan_segments < 8) {
        throw std::runtime_error("scan_segments must be at least 8");
    }

    std::vector<ExtremumPoint> extrema;
    double previous_x = left_bound;
    double previous_derivative = derivative(previous_x);

    for (int i = 1; i <= scan_segments; ++i) {
        const double current_x =
            left_bound +
            (right_bound - left_bound) * static_cast<double>(i) /
                static_cast<double>(scan_segments);
        const double current_derivative = derivative(current_x);

        if (mymath::is_near_zero(previous_derivative, 1e-5)) {
            const double stationary_x = previous_x;
            const double second = second_derivative(stationary_x);
            if (!mymath::is_near_zero(second, 1e-4)) {
                bool duplicate = false;
                for (const ExtremumPoint& point : extrema) {
                    if (same_extremum_x(point.x, stationary_x)) {
                        duplicate = true;
                        break;
                    }
                }
                if (!duplicate) {
                    extrema.push_back(
                        {stationary_x, evaluate_with_variable(stationary_x), second < 0.0});
                }
            }
        } else if ((previous_derivative < 0.0 && current_derivative > 0.0) ||
                   (previous_derivative > 0.0 && current_derivative < 0.0)) {
            const double stationary_x =
                bisect_stationary_point(previous_x, current_x);
            const double second = second_derivative(stationary_x);
            if (!mymath::is_near_zero(second, 1e-4)) {
                extrema.push_back(
                    {stationary_x, evaluate_with_variable(stationary_x), second < 0.0});
            }
        }

        previous_x = current_x;
        previous_derivative = current_derivative;
    }

    std::vector<ExtremumPoint> unique_extrema;
    for (const ExtremumPoint& point : extrema) {
        bool duplicate = false;
        for (const ExtremumPoint& kept : unique_extrema) {
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

double FunctionAnalysis::evaluate_with_variable(double x) const {
    if (expression_.empty()) {
        throw std::runtime_error("function is not defined");
    }

    static constexpr std::size_t kMaxEvaluationCacheSize = 4096;
    const std::string cache_key =
        variable_name_ + "|" + expression_ + "|" + format_double(x);
    const auto found = evaluation_cache_index_.find(cache_key);
    if (found != evaluation_cache_index_.end()) {
        evaluation_cache_entries_.splice(evaluation_cache_entries_.begin(),
                                         evaluation_cache_entries_,
                                         found->second);
        return found->second->second;
    }

    double value = 0.0;
    if (evaluator_) {
        value = evaluator_({{variable_name_, x}});
    } else {
        Calculator calculator;
        calculator.process_line(variable_name_ + " = " + format_double(x), false);
        value = calculator.evaluate_raw(expression_);
    }
    
    evaluation_cache_entries_.push_front({cache_key, value});
    evaluation_cache_index_[cache_key] = evaluation_cache_entries_.begin();
    while (evaluation_cache_entries_.size() > kMaxEvaluationCacheSize) {
        evaluation_cache_index_.erase(evaluation_cache_entries_.back().first);
        evaluation_cache_entries_.pop_back();
    }
    return value;
}

double FunctionAnalysis::second_derivative(double x) const {
    const double step = scale_aware_step(x);
    const long double center = to_long_double(evaluate_with_variable(x));
    const long double left_x = to_long_double(x) - to_long_double(step);
    const long double right_x = to_long_double(x) + to_long_double(step);
    const long double actual_step = (right_x - left_x) * 0.5L;
    const long double left =
        to_long_double(evaluate_with_variable(static_cast<double>(left_x)));
    const long double right =
        to_long_double(evaluate_with_variable(static_cast<double>(right_x)));
    const long double numerator =
        compensated_pair_sum(left - center, right - center);
    return static_cast<double>(
        numerator / (actual_step * actual_step));
}

double FunctionAnalysis::bisect_stationary_point(double left, double right) const {
    double left_derivative = derivative(left);

    for (int i = 0; i < 80; ++i) {
        const double mid = (left + right) * 0.5;
        const double mid_derivative = derivative(mid);
        if (mymath::abs(mid_derivative) <= kRootTolerance ||
            mymath::abs(right - left) <=
                relative_tolerance(kRootTolerance,
                                   std::max(mymath::abs(left), mymath::abs(right)))) {
            return mid;
        }

        if ((left_derivative < 0.0 && mid_derivative > 0.0) ||
            (left_derivative > 0.0 && mid_derivative < 0.0)) {
            right = mid;
        } else {
            left = mid;
            left_derivative = mid_derivative;
        }
    }

    return (left + right) * 0.5;
}

double FunctionAnalysis::adaptive_gauss_kronrod(double left,
                                                double right,
                                                double eps,
                                                int max_depth) const {
    double error = 0.0;
    const double whole = gauss_kronrod_15(left, right, &error);
    return require_finite_integral(
        adaptive_gauss_kronrod_recursive(left,
                                         right,
                                         eps,
                                         whole,
                                         error,
                                         max_depth));
}

double FunctionAnalysis::adaptive_gauss_kronrod_recursive(double left,
                                                          double right,
                                                          double eps,
                                                          double whole,
                                                          double error,
                                                          int depth) const {
    const double scale = std::max(1.0, mymath::abs(whole));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        return whole;
    }

    const double mid = (left + right) * 0.5;
    double left_error = 0.0;
    double right_error = 0.0;
    const double left_area = gauss_kronrod_15(left, mid, &left_error);
    const double right_area = gauss_kronrod_15(mid, right, &right_error);
    const long double left_result = to_long_double(
        adaptive_gauss_kronrod_recursive(left,
                                         mid,
                                         eps * 0.5,
                                         left_area,
                                         left_error,
                                         depth - 1));
    const long double right_result = to_long_double(
        adaptive_gauss_kronrod_recursive(mid,
                                         right,
                                         eps * 0.5,
                                         right_area,
                                         right_error,
                                         depth - 1));
    return static_cast<double>(compensated_pair_sum(left_result, right_result));
}

double FunctionAnalysis::gauss_kronrod_15(double left,
                                          double right,
                                          double* error_estimate) const {
    static const double kNodes[] = {
        0.9914553711208126,
        0.9491079123427585,
        0.8648644233597691,
        0.7415311855993945,
        0.5860872354676911,
        0.4058451513773972,
        0.2077849550078985,
        0.0,
    };
    static const double kKronrodWeights[] = {
        0.02293532201052922,
        0.06309209262997855,
        0.1047900103222502,
        0.1406532597155259,
        0.1690047266392679,
        0.1903505780647854,
        0.2044329400752989,
        0.2094821410847278,
    };
    static const double kGaussWeights[] = {
        0.0,
        0.1294849661688697,
        0.0,
        0.2797053914892767,
        0.0,
        0.3818300505051189,
        0.0,
        0.4179591836734694,
    };

    const long double center =
        (to_long_double(left) + to_long_double(right)) * 0.5L;
    const long double half_width =
        (to_long_double(right) - to_long_double(left)) * 0.5L;
    long double kronrod_sum = 0.0L;
    long double gauss_sum = 0.0L;
    long double kronrod_compensation = 0.0L;
    long double gauss_compensation = 0.0L;

    for (int i = 0; i < 8; ++i) {
        if (mymath::is_near_zero(kNodes[i], 0.0)) {
            const long double value =
                to_long_double(evaluate_with_variable(static_cast<double>(center)));
            compensated_add(static_cast<long double>(kKronrodWeights[i]) * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(static_cast<long double>(kGaussWeights[i]) * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const long double offset = half_width * static_cast<long double>(kNodes[i]);
        const long double left_value =
            to_long_double(evaluate_with_variable(static_cast<double>(center - offset)));
        const long double right_value =
            to_long_double(evaluate_with_variable(static_cast<double>(center + offset)));
        const long double pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(static_cast<long double>(kKronrodWeights[i]) * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(static_cast<long double>(kGaussWeights[i]) * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const long double kronrod = half_width * kronrod_sum;
    const long double gauss = half_width * gauss_sum;
    *error_estimate = static_cast<double>(mymath::abs_long_double(kronrod - gauss));
    return static_cast<double>(kronrod);
}

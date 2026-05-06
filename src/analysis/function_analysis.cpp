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

/** @brief 泛型绝对值函数 */
template <typename T>
T t_abs(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::abs(val);
    } else {
        return val < T(static_cast<long long>(0)) ? -val : val;
    }
}

/** @brief 泛型平方根函数 */
template <typename T>
T t_sqrt(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::sqrt(val);
    } else {
        throw std::runtime_error("t_sqrt not implemented for this type");
    }
}

/** @brief 泛型幂函数 */
template <typename T>
T t_pow(const T& base, const T& exponent) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::pow(base, exponent);
    } else {
        throw std::runtime_error("t_pow not implemented for this type");
    }
}

/** @brief 泛型有限值检查 */
template <typename T>
bool t_isfinite(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::isfinite(val);
    } else {
        return true;
    }
}

/** @brief 泛型正弦函数 */
template <typename T>
T t_sin(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::sin(val);
    } else {
        return T(mymath::sin(val.to_double()));
    }
}

/** @brief 泛型余弦函数 */
template <typename T>
T t_cos(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::cos(val);
    } else {
        return T(mymath::cos(val.to_double()));
    }
}

/** @brief 泛型正切函数 */
template <typename T>
T t_tan(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::tan(val);
    } else {
        return T(mymath::tan(val.to_double()));
    }
}

/** @brief 泛型自然对数函数 */
template <typename T>
T t_log(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::log(val);
    } else {
        return T(mymath::log(val.to_double()));
    }
}

/** @brief 泛型指数函数 */
template <typename T>
T t_exp(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::exp(val);
    } else {
        return T(mymath::exp(val.to_double()));
    }
}

/** @brief 泛型双曲正弦函数 */
template <typename T>
T t_sinh(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::sinh(val);
    } else {
        return T(mymath::sinh(val.to_double()));
    }
}

/** @brief 泛型双曲余弦函数 */
template <typename T>
T t_cosh(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::cosh(val);
    } else {
        return T(mymath::cosh(val.to_double()));
    }
}

/** @brief 泛型双曲正切函数 */
template <typename T>
T t_tanh(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::tanh(val);
    } else {
        return T(mymath::tanh(val.to_double()));
    }
}

/** @brief 泛型圆周率 */
template <typename T>
T t_pi() {
    if constexpr (std::is_same_v<T, long double>) {
        return 3.1415926535897932384626433832795029L;
    } else if constexpr (std::is_same_v<T, long double>) {
        return mymath::kPi;
    } else {
        return T(3.14159265358979323846);
    }
}

/** @brief 泛型无穷大 */
template <typename T>
T t_infinity() {
    if constexpr (std::is_floating_point_v<T>) {
        return std::numeric_limits<T>::infinity();
    } else {
        return T(1e300);
    }
}

template <typename T>
bool t_is_effective_infinity_point(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return !std::isfinite(val);
    } else {
        return false;
    }
}

template <typename T>
bool symbolic_limit_at_infinity(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                bool positive,
                                T* result);

/** @brief 泛型整数判断 */
template <typename T>
bool t_is_integer(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::floor(val) == val;
    } else {
        return mymath::is_integer(val.to_double());
    }
}

/** @brief 泛型零判断 */
template <typename T>
bool t_is_near_zero(const T& val, const T& eps) {
    return t_abs(val) <= eps;
}

/** @brief 泛型四舍五入到 long long */
template <typename T>
long long t_llround(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return static_cast<long long>(std::llround(val));
    } else {
        return static_cast<long long>(std::llround(val.to_double()));
    }
}

/** @brief 导数计算的基准步长 */
template <typename T>
T kDerivativeBaseStep_v() {
    if constexpr (std::is_same_v<T, long double>) {
        return 1e-5L;
    } else {
        return T(1e-4);
    }
}

/** @brief 极限计算的初始步长 */
template <typename T>
T kLimitInitialStep_v() {
    return T(1e-1);
}

/** @brief 极限计算的收敛容差 */
template <typename T>
T kLimitTolerance_v() {
    if constexpr (std::is_same_v<T, long double>) {
        return 1e-15L;
    } else {
        return T(1e-10L);
    }
}

/** @brief 根查找的收敛容差 */
template <typename T>
T kRootTolerance_v() {
    if constexpr (std::is_same_v<T, long double>) {
        return 1e-12L;
    } else {
        return T(1e-7L);
    }
}

/** @brief 数值积分的精度要求 */
template <typename T>
T kIntegralTolerance_v() {
    if constexpr (std::is_same_v<T, long double>) {
        return 1e-12L;
    } else {
        return T(1e-8L);
    }
}

/** @brief 自适应积分的最大递归深度 */
constexpr int kMaxIntegralDepth = 18;

template <typename T>
std::string format_t(const T& value) {
    if constexpr (std::is_floating_point_v<T>) {
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<T>::max_digits10) << value;
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
    } else {
        return value.to_string();
    }
}

template <typename T>
void compensated_add(T value,
                     T* sum,
                     T* compensation) {
    const T adjusted = value - *compensation;
    const T next = *sum + adjusted;
    *compensation = (next - *sum) - adjusted;
    *sum = next;
}

template <typename T>
T compensated_pair_sum(T lhs, T rhs) {
    T sum = T(static_cast<long long>(0));
    T compensation = T(static_cast<long long>(0));
    compensated_add(lhs, &sum, &compensation);
    compensated_add(rhs, &sum, &compensation);
    return sum;
}

template <typename T>
T scale_aware_step(T x) {
    const T scale = std::max(T(static_cast<long long>(1)), t_abs(x));
    return kDerivativeBaseStep_v<T>() * scale;
}

template <typename T>
T central_difference_step_value(T scale, T factor) {
    T base_step = T(1e-7L);
    return std::max(base_step * scale, kDerivativeBaseStep_v<T>() * scale * factor);
}

template <typename T>
T numeric_control_value(const char*, long double fallback) {
    return T(fallback);
}

template <typename T>
T derivative_quarter_power_scale(const T& value) {
    return t_pow(value, T(0.25));
}

template <typename T>
T relative_tolerance(T baseline, T scale) {
    return baseline * std::max(T(static_cast<long long>(1)), scale);
}

template <typename T>
T limit_step_scale(T x) {
    return kLimitInitialStep_v<T>() * std::max(T(static_cast<long long>(1)), t_abs(x));
}

template <typename T>
bool same_extremum_x(T lhs, T rhs) {
    return t_abs(lhs - rhs) <= T(1e-5);
}

template <typename T>
T require_finite_integral(T value) {
    if (!t_isfinite(value)) {
        throw std::runtime_error("integral did not converge");
    }
    return value;
}

template <typename T>
void reject_divergent_transformed_endpoint(
    const std::function<T(T)>& transformed,
    bool check_left,
    bool check_right) {
    const T offsets[] = {
        numeric_control_value<T>("1e-3", 1e-3),
        numeric_control_value<T>("1e-4", 1e-4),
        numeric_control_value<T>("1e-5", 1e-5),
    };
    auto check_at = [&](T t) {
        T value = transformed(t);
        if (!t_isfinite(value) || t_abs(value) > T(static_cast<long long>(10000))) {
            throw std::runtime_error("integral did not converge");
        }
    };

    for (T offset : offsets) {
        if (check_left) {
            check_at(offset);
        }
        if (check_right) {
            check_at(T(static_cast<long long>(1)) - offset);
        }
    }
}

// 自适应 Simpson 积分辅助函数（用于 callable，支持高精度）
template <typename T>
T simpson_rule_callable(const std::function<T(T)>& func, T a, T b) {
    const T h = (b - a) / T(static_cast<long long>(2));
    const T fa = func(a);
    const T fb = func(b);
    const T fc = func((a + b) / T(static_cast<long long>(2)));
    return h / T(static_cast<long long>(3)) * (fa + T(static_cast<long long>(4)) * fc + fb);
}

template <typename T>
T adaptive_simpson_callable_recursive(const std::function<T(T)>& func,
                                       T a, T b, T whole, T left, T right, T eps, int depth) {
    const T c = (a + b) / T(static_cast<long long>(2));
    const T combined = left + right;
    const T error = t_abs(combined - whole) / T(static_cast<long long>(15));

    const T scale = std::max(T(static_cast<long long>(1)), t_abs(combined));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        return combined + (combined - whole) / T(static_cast<long long>(15));
    }

    const T d = (a + c) / T(static_cast<long long>(2));
    const T e = (c + b) / T(static_cast<long long>(2));
    const T left_left = simpson_rule_callable(func, a, d);
    const T left_right = simpson_rule_callable(func, d, c);
    const T right_left = simpson_rule_callable(func, c, e);
    const T right_right = simpson_rule_callable(func, e, b);

    return adaptive_simpson_callable_recursive(func, a, c, left, left_left, left_right, eps / T(static_cast<long long>(2)), depth - 1) +
           adaptive_simpson_callable_recursive(func, c, b, right, right_left, right_right, eps / T(static_cast<long long>(2)), depth - 1);
}

template <typename T>
T adaptive_simpson_callable(const std::function<T(T)>& func, T left, T right, T eps, int max_depth) {
    const T c = (left + right) / T(static_cast<long long>(2));
    const T whole = simpson_rule_callable(func, left, right);
    const T left_val = simpson_rule_callable(func, left, c);
    const T right_val = simpson_rule_callable(func, c, right);
    return adaptive_simpson_callable_recursive(func, left, right, whole, left_val, right_val, eps, max_depth);
}

template <typename T>
T gauss_kronrod_15_callable(const std::function<T(T)>& function,
                                 T left,
                                 T right,
                                 T* error_estimate) {
    static const T kNodes[] = {
        T(0.9914553711208126),
        T(0.9491079123427585),
        T(0.8648644233597691),
        T(0.7415311855993945),
        T(0.5860872354676911),
        T(0.4058451513773972),
        T(0.2077849550078985),
        T(static_cast<long long>(0)),
    };
    static const T kKronrodWeights[] = {
        T(0.02293532201052922),
        T(0.06309209262997855),
        T(0.1047900103222502),
        T(0.1406532597155259),
        T(0.1690047266392679),
        T(0.1903505780647854),
        T(0.2044329400752989),
        T(0.2094821410847278),
    };
    static const T kGaussWeights[] = {
        T(static_cast<long long>(0)),
        T(0.1294849661688697),
        T(static_cast<long long>(0)),
        T(0.2797053914892767),
        T(static_cast<long long>(0)),
        T(0.3818300505051189),
        T(static_cast<long long>(0)),
        T(0.4179591836734694),
    };

    const T center = (left + right) * T(0.5L);
    const T half_width = (right - left) * T(0.5L);
    T kronrod_sum = T(static_cast<long long>(0));
    T gauss_sum = T(static_cast<long long>(0));
    T kronrod_compensation = T(static_cast<long long>(0));
    T gauss_compensation = T(static_cast<long long>(0));

    for (int i = 0; i < 8; ++i) {
        if (t_is_near_zero(kNodes[i], T(static_cast<long long>(0)))) {
            const T value = function(center);
            compensated_add(kKronrodWeights[i] * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(kGaussWeights[i] * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const T offset = half_width * kNodes[i];
        const T left_value = function(center - offset);
        const T right_value = function(center + offset);
        const T pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(kKronrodWeights[i] * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(kGaussWeights[i] * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const T kronrod = half_width * kronrod_sum;
    const T gauss = half_width * gauss_sum;
    *error_estimate = t_abs(kronrod - gauss);
    return kronrod;
}

template <typename T>
T adaptive_gauss_kronrod_callable_recursive(
    const std::function<T(T)>& function,
    T left,
    T right,
    T eps,
    T whole,
    T error,
    int depth) {
    const T scale = std::max(T(static_cast<long long>(1)), t_abs(whole));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        return whole;
    }

    const T mid = (left + right) * T(0.5L);
    T left_error = T(static_cast<long long>(0));
    T right_error = T(static_cast<long long>(0));
    const T left_area =
        gauss_kronrod_15_callable(function, left, mid, &left_error);
    const T right_area =
        gauss_kronrod_15_callable(function, mid, right, &right_error);
    const T left_result =
        adaptive_gauss_kronrod_callable_recursive(function,
                                                  left,
                                                  mid,
                                                  eps * T(0.5L),
                                                  left_area,
                                                  left_error,
                                                  depth - 1);
    const T right_result =
        adaptive_gauss_kronrod_callable_recursive(function,
                                                  mid,
                                                  right,
                                                  eps * T(0.5L),
                                                  right_area,
                                                  right_error,
                                                  depth - 1);
    return compensated_pair_sum(left_result, right_result);
}

template <typename T>
T adaptive_gauss_kronrod_callable(const std::function<T(T)>& function,
                                       T left,
                                       T right,
                                       T eps,
                                       int depth) {
    T error = T(static_cast<long long>(0));
    const T whole = gauss_kronrod_15_callable(function, left, right, &error);
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

template <typename T>
SymbolicLimitProbeKind probe_symbolic_value_at(
    SymbolicExpression expression,
    const std::string& variable_name,
    T point,
    T* finite_value);

bool is_infinite_probe(SymbolicLimitProbeKind kind) {
    return kind == SymbolicLimitProbeKind::kPositiveInfinity ||
           kind == SymbolicLimitProbeKind::kNegativeInfinity;
}

template <typename T>
bool try_symbolic_one_to_infinity_limit(const SymbolicExpression& base,
                                        const SymbolicExpression& exponent,
                                        const std::string& variable_name,
                                        T point,
                                        T* result) {
    T exponent_value = T(static_cast<long long>(0));
    const SymbolicLimitProbeKind exponent_kind =
        probe_symbolic_value_at(exponent, variable_name, point, &exponent_value);
    if (!is_infinite_probe(exponent_kind)) {
        return false;
    }

    T base_value = T(static_cast<long long>(0));
    const SymbolicLimitProbeKind base_kind =
        probe_symbolic_value_at(base, variable_name, point, &base_value);
    if (base_kind != SymbolicLimitProbeKind::kFinite ||
        !t_is_near_zero(base_value - T(static_cast<long long>(1)),
                        numeric_control_value<T>("1e-8", 1e-8))) {
        return false;
    }

    const SymbolicExpression transformed =
        symbolic_expression_internal::make_function(
            "exp",
            ((base - SymbolicExpression::number(1.0L)) * exponent).simplify()).simplify();
    const SymbolicExpression product =
        ((base - SymbolicExpression::number(1.0L)) * exponent).simplify();
    T product_limit = T(static_cast<long long>(0));
    if (symbolic_limit_at_infinity(product, variable_name, point > T(static_cast<long long>(0)), &product_limit)) {
        *result = t_exp(product_limit);
        return true;
    }

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
                    long double value = 0.0L;
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
                        t_is_near_zero(left_value, kLimitTolerance_v<T>())) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    if (right_kind == SymbolicLimitProbeKind::kFinite &&
                        t_is_near_zero(right_value, kLimitTolerance_v<T>())) {
                        return SymbolicLimitProbeKind::kUnknown;
                    }
                    auto sign_of = [](SymbolicLimitProbeKind kind, T value) {
                        if (kind == SymbolicLimitProbeKind::kFinite) {
                            return value >= T(static_cast<long long>(0)) ? 1 : -1;
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
                        !t_is_near_zero(right_value, kLimitTolerance_v<T>())) {
                        *finite_value = left_value / right_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (left_kind == SymbolicLimitProbeKind::kFinite &&
                        is_infinite_probe(right_kind)) {
                        *finite_value = T(static_cast<long long>(0));
                        return SymbolicLimitProbeKind::kFinite;
                    }
                    if (is_infinite_probe(left_kind) &&
                        right_kind == SymbolicLimitProbeKind::kFinite &&
                        !t_is_near_zero(right_value, kLimitTolerance_v<T>())) {
                        const bool positive =
                            (left_kind == SymbolicLimitProbeKind::kPositiveInfinity) ==
                            (right_value > T(static_cast<long long>(0)));
                        return positive ? SymbolicLimitProbeKind::kPositiveInfinity
                                        : SymbolicLimitProbeKind::kNegativeInfinity;
                    }
                    return SymbolicLimitProbeKind::kUnknown;
                }
                case NodeType::kPower: {
                    const SymbolicExpression base(node->left);
                    const SymbolicExpression exponent(node->right);
                    T one_to_infinity_value = T(static_cast<long long>(0));
                    if (try_symbolic_one_to_infinity_limit(base,
                                                           exponent,
                                                           variable_name,
                                                           point,
                                                           &one_to_infinity_value)) {
                        *finite_value = one_to_infinity_value;
                        return SymbolicLimitProbeKind::kFinite;
                    }

                    long double exponent_number = 0.0L;
                    if (!exponent.is_number(&exponent_number)) {
                        T base_value = T(static_cast<long long>(0));
                        T exponent_value = T(static_cast<long long>(0));
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
                            base_value > T(static_cast<long long>(0))) {
                            *finite_value = t_exp(exponent_value * t_log(base_value));
                            return t_isfinite(*finite_value)
                                       ? SymbolicLimitProbeKind::kFinite
                                       : SymbolicLimitProbeKind::kUnknown;
                        }
                        return SymbolicLimitProbeKind::kUnknown;
                    }

                    const T exponent_value(exponent_number);
                    T base_value = T(static_cast<long long>(0));
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
                            *finite_value = T(static_cast<long long>(0));
                            return SymbolicLimitProbeKind::kFinite;
                        }
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
                            argument_value >= T(static_cast<long long>(0))) {
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
                            *finite_value = T(static_cast<long long>(1));
                            return SymbolicLimitProbeKind::kFinite;
                        }
                        if (argument_kind == SymbolicLimitProbeKind::kNegativeInfinity) {
                            *finite_value = -T(static_cast<long long>(1));
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

        long double p_val;
        if constexpr (std::is_floating_point_v<T>) {
            p_val = static_cast<long double>(point);
        } else {
            p_val = point.to_double();
        }

        SymbolicExpression sub_expr = expression.substitute(
            variable_name,
            SymbolicExpression::number(p_val)).simplify();
        long double value = 0.0L;
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
bool is_zero_probe(SymbolicLimitProbeKind kind, T value) {
    return kind == SymbolicLimitProbeKind::kFinite &&
           t_is_near_zero(value, kLimitTolerance_v<T>());
}

/**
 * @brief 处理极点极限
 */
template <typename T>
T handle_pole_limit(int shift, T leading_coefficient, int direction) {
    if (direction == 0) {
        // 双侧极限：只有当 shift 为偶数时才存在
        if (shift % 2 == 0) {
            return (leading_coefficient > T(static_cast<long long>(0))) ? t_infinity<T>() : -t_infinity<T>();
        } else {
            throw std::runtime_error("two-sided limit does not exist (pole with odd shift)");
        }
    } else if (direction == 1) {
        // 右极限：(x - x0) > 0，符号不变
        return (leading_coefficient > T(static_cast<long long>(0))) ? t_infinity<T>() : -t_infinity<T>();
    } else {
        // 左极限：(x - x0) < 0，奇数 shift 时符号翻转
        bool flip_sign = (shift % 2 != 0);
        T effective_c = flip_sign ? -leading_coefficient : leading_coefficient;
        return (effective_c > T(static_cast<long long>(0))) ? t_infinity<T>() : -t_infinity<T>();
    }
}

template <typename T>
bool try_symbolic_lhopital_limit(const SymbolicExpression& expression,
                                 const std::string& variable_name,
                                 T point,
                                 int direction,
                                 T* result) {
    SymbolicExpression current = expression.simplify();
    if (current.node_->type != NodeType::kDivide) {
        return false;
    }

    SymbolicExpression numerator(current.node_->left);
    SymbolicExpression denominator(current.node_->right);
    T numerator_value = T(static_cast<long long>(0));
    T denominator_value = T(static_cast<long long>(0));
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
        ctx.evaluate_at = [](const SymbolicExpression& e, const std::string& /*v*/, long double /*p*/) {
            long double val = 0.0L;
            if (e.is_number(&val)) return val;
            return 0.0L;
        };

        struct LaurentInfo {
            int degree = 0;
            T coefficient = T(static_cast<long long>(0));
            bool valid = false;
        };

        auto get_laurent_info = [&](const SymbolicExpression& e) -> LaurentInfo {
            LaurentInfo info;
            std::vector<long double> coeffs;
            try {
                long double p_val;
                if constexpr (std::is_floating_point_v<T>) p_val = static_cast<long double>(point);
                else p_val = point.to_double();

                if (series_ops::internal::evaluate_psa(e, variable_name, p_val, 4, coeffs, ctx)) {
                    for (int i = 0; i < static_cast<int>(coeffs.size()); ++i) {
                        if (!mymath::is_near_zero(coeffs[i], 1e-15)) {
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

    static constexpr int kMaxLhopitalDepth = 3;
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
        if (is_infinite_probe(kind)) {
            *result = (kind == SymbolicLimitProbeKind::kPositiveInfinity) ? t_infinity<T>() : -t_infinity<T>();
            return true;
        }
    }

    return false;
}

template <typename T>
bool symbolic_limit_at_infinity(const SymbolicExpression& expression,
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
    ctx.evaluate_at = [](const SymbolicExpression& e, const std::string& /*v*/, long double /*p*/) {
        long double val = 0.0L;
        if (e.is_number(&val)) return val;
        return 0.0L;
    };

    SymbolicExpression t_var = SymbolicExpression::variable("t_limit_inf_tmp");
    SymbolicExpression inv_t;
    if (positive) {
        inv_t = SymbolicExpression::number(1.0L) / t_var;
    } else {
        inv_t = SymbolicExpression::number(-1.0L) / t_var;
    }
    SymbolicExpression substituted = expression.substitute(variable_name, inv_t).simplify();

    std::vector<long double> coeffs;
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
            long double large_x = positive ? 1e12 : -1e12;
            long double val = ctx.evaluate_at(expression, variable_name, large_x);
            if (mymath::isfinite(val) && mymath::abs(val) < 1e10) {
                return false;
            }
        } catch (...) {
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

}  // namespace

template <typename T>
TFunctionAnalysis<T>::TFunctionAnalysis(std::string variable_name)
    : variable_name_(std::move(variable_name)) {
    if (!is_valid_analysis_variable_name(variable_name_)) {
        throw std::runtime_error("invalid variable name for custom function");
    }
}

template <typename T>
TFunctionAnalysis<T>::TFunctionAnalysis(const TFunctionAnalysis& other)
    : expression_(other.expression_),
      variable_name_(other.variable_name_),
      evaluator_(other.evaluator_) {}

template <typename T>
TFunctionAnalysis<T>& TFunctionAnalysis<T>::operator=(const TFunctionAnalysis& other) {
    if (this != &other) {
        expression_ = other.expression_;
        variable_name_ = other.variable_name_;
        evaluator_ = other.evaluator_;
        evaluation_cache_entries_.clear();
        evaluation_cache_index_.clear();
    }
    return *this;
}

template <typename T>
TFunctionAnalysis<T>::TFunctionAnalysis(TFunctionAnalysis&& other) noexcept = default;

template <typename T>
TFunctionAnalysis<T>& TFunctionAnalysis<T>::operator=(TFunctionAnalysis&& other) noexcept = default;

template <typename T>
TFunctionAnalysis<T>::~TFunctionAnalysis() = default;

template <typename T>
void TFunctionAnalysis<T>::define(const std::string& expression) {
    if (expression.empty()) {
        throw std::runtime_error("function expression cannot be empty");
    }

    expression_ = expression;
    evaluation_cache_entries_.clear();
    evaluation_cache_index_.clear();
}

template <typename T>
void TFunctionAnalysis<T>::set_evaluator(std::function<T(const std::vector<std::pair<std::string, T>>&)> evaluator) {
    evaluator_ = std::move(evaluator);
}

template <typename T>
T TFunctionAnalysis<T>::evaluate(T x) const {
    return evaluate_with_variable(x);
}

template <typename T>
T TFunctionAnalysis<T>::derivative(T x) const {
    const T scale = std::max(T(static_cast<long long>(1)), t_abs(x));
    const T center = evaluate_with_variable(x);
    if (!t_isfinite(center)) {
        throw std::runtime_error("derivative is undefined at this point");
    }

    constexpr int kLayers = 4;

    // 计算曲率以确定步长
    const T curvature_sample_step = scale * T(1e-3);
    const T curvature_probe = evaluate_with_variable(x + curvature_sample_step) -
                                   T(static_cast<long long>(2)) * center +
                                   evaluate_with_variable(x - curvature_sample_step);
    const T curvature_scale =
        std::max(T(static_cast<long long>(1)),
                 t_abs(curvature_probe) /
                     std::max(T(1e-12L), t_abs(center)));

    T base_step = central_difference_step_value(
        scale,
        T(static_cast<long long>(1)) / derivative_quarter_power_scale(curvature_scale));

    std::vector<std::vector<T>> richardson(kLayers, std::vector<T>(kLayers, T(0)));
    std::vector<bool> row_valid(kLayers, false);
    T best_value = T(static_cast<long long>(0));
    T best_error = t_infinity<T>();

    for (int row = 0; row < kLayers; ++row) {
        const T step = base_step / t_pow(T(static_cast<long long>(2)), T(static_cast<long long>(row)));
        const T forward_x = x + step;
        const T backward_x = x - step;
        const T actual_step = (forward_x - backward_x) * T(0.5L);
        if (actual_step <= T(static_cast<long long>(0))) {
            continue;
        }
        const T forward = evaluate_with_variable(forward_x);
        const T backward = evaluate_with_variable(backward_x);
        if (!t_isfinite(forward) || !t_isfinite(backward)) {
            continue;
        }
        richardson[row][0] = (forward - backward) / (T(static_cast<long long>(2)) * actual_step);
        row_valid[row] = t_isfinite(richardson[row][0]);
        if (!row_valid[row]) {
            continue;
        }
        for (int col = 1; col <= row; ++col) {
            if (!row_valid[row - 1]) {
                row_valid[row] = false;
                break;
            }
            const T factor = t_pow(T(static_cast<long long>(4)), T(static_cast<long long>(col)));
            richardson[row][col] =
                richardson[row][col - 1] +
                (richardson[row][col - 1] - richardson[row - 1][col - 1]) /
                    (factor - T(static_cast<long long>(1)));
            if (!t_isfinite(richardson[row][col])) {
                row_valid[row] = false;
                break;
            }
        }
        if (row > 0 && row_valid[row] && row_valid[row - 1]) {
            const T candidate = richardson[row][row];
            const T error_estimate = t_abs(candidate - richardson[row - 1][row - 1]);
            if (error_estimate < best_error && t_isfinite(candidate)) {
                best_error = error_estimate;
                best_value = candidate;

                T tol = std::max(T(1e-12L), t_abs(best_value) * T(1e-14));

                if (error_estimate < tol) break;
            }
        }
    }

    if (best_error < t_infinity<T>()) {
        const T side_step = std::max(
            numeric_control_value<T>("1e-7", 1e-7) * scale,
            base_step / T(static_cast<long long>(64)));
        const T left_value = evaluate_with_variable(x - side_step);
        const T right_value = evaluate_with_variable(x + side_step);
        if (!t_isfinite(left_value) || !t_isfinite(right_value)) {
            throw std::runtime_error("derivative is undefined at this point");
        }
        const T left_slope = (center - left_value) / side_step;
        const T right_slope = (right_value - center) / side_step;
        const T slope_scale =
            std::max({T(static_cast<long long>(1)), t_abs(left_slope), t_abs(right_slope), t_abs(best_value)});
        if (t_abs(left_slope - right_slope) >
            std::max(numeric_control_value<T>("1e-4", 1e-4),
                     numeric_control_value<T>("1e-5", 1e-5) * slope_scale)) {
            throw std::runtime_error("derivative does not exist at this point");
        }
        return best_value;
    }
    for (int row = kLayers - 1; row >= 0; --row) {
        if (row_valid[row]) {
            return richardson[row][row];
        }
    }
    return richardson[0][0];
}

template <typename T>
T TFunctionAnalysis<T>::limit(T x, int direction) const {
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
    ctx.evaluate_at = [this](const SymbolicExpression& e, const std::string& v, long double p) {
        if (v == variable_name_) return p;
        long double val = 0.0L;
        if (e.is_number(&val)) return val;
        return 0.0L;
    };

    if (t_is_effective_infinity_point(x)) {
        bool positive = x > T(static_cast<long long>(0));
        T inf_result = T(static_cast<long long>(0));
        if (symbolic_limit_at_infinity(expr, variable_name_, positive, &inf_result)) {
            return inf_result;
        }
    } else if (t_isfinite(x)) {
        std::vector<long double> coeffs;
        try {
            long double p_val;
            if constexpr (std::is_floating_point_v<T>) p_val = static_cast<long double>(x);
            else p_val = x.to_double();

            if (series_ops::internal::evaluate_psa(expr, variable_name_, p_val, 2, coeffs, ctx)) {
                if (!coeffs.empty()) return T(coeffs[0]);
            }
        } catch (const series_ops::internal::PoleException& e) {
            return handle_pole_limit(e.shift, T(e.leading_coefficient), direction);
        }
    }

    T lhopital_value = T(static_cast<long long>(0));
    if (direction == 0 &&
        try_symbolic_lhopital_limit(expr,
                                    variable_name_,
                                    x,
                                    direction,
                                    &lhopital_value)) {
        return lhopital_value;
    }

    return compute_numerical_limit(x, direction);
}

template <typename T>
T TFunctionAnalysis<T>::compute_numerical_limit(T x, int direction) const {
    auto compute_limit_at = [this](T x_target, int side) {
        T richardson[14][14] = {};
        bool row_valid[14] = {};
        T best_value = T(static_cast<long long>(0));
        T best_error = t_infinity<T>();
        bool have_best = false;

        const T base_h = t_is_effective_infinity_point(x_target)
                             ? numeric_control_value<T>("1e-2", 1e-2)
                             : limit_step_scale(x_target);
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
                val = evaluate_with_variable(sample_x);
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
                    if (prev_val > numeric_control_value<T>("1e10", 1e10)) {
                        return t_infinity<T>();
                    } else if (prev_val < -numeric_control_value<T>("1e10", 1e10)) {
                        return -t_infinity<T>();
                    }
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
                if ((val > T(static_cast<long long>(0)) && prev_val < T(static_cast<long long>(0))) || (val < T(static_cast<long long>(0)) && prev_val > T(static_cast<long long>(0)))) {
                    oscillation_count++;
                    if (oscillation_count >= 5) {
                        const T avg_amp = total_amplitude / T(static_cast<long long>(row + 1));
                        if (avg_amp > numeric_control_value<T>("1e-2", 1e-2)) {
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
                const T expected_change =
                    best_error * T(static_cast<long long>(10)) +
                    numeric_control_value<T>("1e-10", 1e-10);
                const T actual_change = t_abs(val - best_value);
                if (actual_change > expected_change * T(1e6)) {
                    adaptive_h *= T(0.5L);
                    consecutive_bad++;
                    if (consecutive_bad >= kMaxBadSamples) {
                        break;
                    }
                    continue;
                }
            }

            consecutive_bad = 0;
            richardson[row][0] = val;
            row_valid[row] = true;

            for (int col = 1; col <= row; ++col) {
                if (!row_valid[row - 1]) {
                    row_valid[row] = false;
                    break;
                }
                const T factor = T(static_cast<long long>(1LL << col));
                richardson[row][col] = richardson[row][col-1] +
                    (richardson[row][col-1] - richardson[row-1][col-1]) / (factor - T(static_cast<long long>(1)));
                if (!t_isfinite(richardson[row][col])) {
                    row_valid[row] = false;
                    break;
                }
            }

            if (row > 0 && row_valid[row]) {
                T err = t_abs(richardson[row][row] - richardson[row-1][row-1]);
                if (err < best_error) {
                    best_error = err;
                    best_value = richardson[row][row];
                    have_best = true;
                }

                if (err > best_error * T(0.9) && row > 3) {
                    adaptive_h *= numeric_control_value<T>("0.75", 0.75);
                }

                if (err < numeric_control_value<T>("1e-15", 1e-15)) {
                    return richardson[row][row];
                }
            } else if (row == 0 && row_valid[0]) {
                best_value = richardson[0][0];
            }
        }

        if (!have_best) {
            throw std::runtime_error("limit did not converge");
        }

        const T scale = std::max(T(static_cast<long long>(1)), t_abs(best_value));
        const T acceptable_error =
            std::max(numeric_control_value<T>("1e-8", 1e-8),
                     kLimitTolerance_v<T>() * T(static_cast<long long>(1000)) * scale);
        if (best_error > acceptable_error) {
            throw std::runtime_error("limit did not converge");
        }
        return best_value;
    };

    if (t_is_effective_infinity_point(x)) {
        return compute_limit_at(x, 0);
    }

    if (direction == -1) return compute_limit_at(x, -1);
    if (direction == 1) return compute_limit_at(x, 1);

    T left = T(static_cast<long long>(0)), right = T(static_cast<long long>(0));
    bool left_ok = false, right_ok = false;
    std::string left_err, right_err;

    try {
        left = compute_limit_at(x, -1);
        left_ok = true;
    } catch (const std::exception& e) {
        left_err = e.what();
    }

    try {
        right = compute_limit_at(x, 1);
        right_ok = true;
    } catch (const std::exception& e) {
        right_err = e.what();
    }

    if (left_ok && right_ok) {
        if (t_abs(left - right) <= kLimitTolerance_v<T>() * T(static_cast<long long>(100)) ||
            (!t_isfinite(left) && !t_isfinite(right) && ((left > T(static_cast<long long>(0))) == (right > T(static_cast<long long>(0)))))) {
            return (left + right) * T(0.5L);
        }
        throw std::runtime_error("two-sided limit does not exist (left=" +
                                 format_t(left) + ", right=" + format_t(right) + ")");
    } else if (left_ok) {
        return left;
    } else if (right_ok) {
        return right;
    }

    throw std::runtime_error("limit did not converge on either side: " + left_err + " / " + right_err);
}

template <typename T>
T TFunctionAnalysis<T>::definite_integral(T lower_bound,
                                           T upper_bound) const {
    if (t_is_near_zero(lower_bound - upper_bound, T(1e-15L))) {
        return T(static_cast<long long>(0));
    }
    if (lower_bound > upper_bound) {
        return -definite_integral(upper_bound, lower_bound);
    }
    const bool lower_is_infinite = !t_isfinite(lower_bound);
    const bool upper_is_infinite = !t_isfinite(upper_bound);
    if (lower_is_infinite || upper_is_infinite) {
        if (lower_is_infinite && upper_is_infinite) {
            if (lower_bound > T(static_cast<long long>(0)) || upper_bound < T(static_cast<long long>(0))) {
                throw std::runtime_error("invalid infinite integration bounds");
            }
            auto transformed = [this](T t) {
                const T angle = t_pi<T>() * (t - T(0.5L));
                const T cos_angle = t_cos(angle);
                const T x = t_tan(angle);
                return evaluate_with_variable(x) * t_pi<T>() / (cos_angle * cos_angle);
            };
            reject_divergent_transformed_endpoint<T>(std::function<T(T)>(transformed), true, true);
            return adaptive_gauss_kronrod_callable<T>(std::function<T(T)>(transformed),
                                                   T(static_cast<long long>(0)),
                                                   T(static_cast<long long>(1)),
                                                   kIntegralTolerance_v<T>(),
                                                   kMaxIntegralDepth);
        }

        if (lower_is_infinite) {
            if (lower_bound > T(static_cast<long long>(0))) {
                throw std::runtime_error("invalid infinite integration bounds");
            }
            auto transformed = [this, upper_bound](T t) {
                const T x = upper_bound - (T(static_cast<long long>(1)) - t) / t;
                return evaluate_with_variable(x) / (t * t);
            };
            reject_divergent_transformed_endpoint<T>(std::function<T(T)>(transformed), true, false);
            return adaptive_gauss_kronrod_callable<T>(std::function<T(T)>(transformed),
                                                   T(static_cast<long long>(0)),
                                                   T(static_cast<long long>(1)),
                                                   kIntegralTolerance_v<T>(),
                                                   kMaxIntegralDepth);
        }

        if (upper_bound < T(static_cast<long long>(0))) {
            throw std::runtime_error("invalid infinite integration bounds");
        }
        auto transformed = [this, lower_bound](T t) {
            const T one_minus_t = T(static_cast<long long>(1)) - t;
            const T x = lower_bound + t / one_minus_t;
            return evaluate_with_variable(x) / (one_minus_t * one_minus_t);
        };
        reject_divergent_transformed_endpoint<T>(std::function<T(T)>(transformed), false, true);
        return adaptive_gauss_kronrod_callable<T>(std::function<T(T)>(transformed),
                                               T(static_cast<long long>(0)),
                                               T(static_cast<long long>(1)),
                                               kIntegralTolerance_v<T>(),
                                               kMaxIntegralDepth);
    }
    const T span = t_abs(upper_bound - lower_bound);
    const T scaled_eps =
        relative_tolerance(kIntegralTolerance_v<T>(), span + t_abs(lower_bound) + t_abs(upper_bound));
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
        const T width = upper_bound - lower_bound;
        auto transformed = [this, lower_bound, upper_bound, width,
                            left_singular, right_singular](T t) {
            if (left_singular && right_singular) {
                const T s = t * t * (T(static_cast<long long>(3)) - T(static_cast<long long>(2)) * t);
                const T dx_dt = width * T(static_cast<long long>(6)) * t * (T(static_cast<long long>(1)) - t);
                return evaluate_with_variable(lower_bound + width * s) * dx_dt;
            }
            if (left_singular) {
                const T x = lower_bound + width * t * t;
                return evaluate_with_variable(x) * T(static_cast<long long>(2)) * width * t;
            }
            const T one_minus_t = T(static_cast<long long>(1)) - t;
            const T x = upper_bound - width * one_minus_t * one_minus_t;
            return evaluate_with_variable(x) * T(static_cast<long long>(2)) * width * one_minus_t;
        };
        reject_divergent_transformed_endpoint<T>(std::function<T(T)>(transformed),
                                              left_singular,
                                              right_singular);
        return adaptive_gauss_kronrod_callable<T>(std::function<T(T)>(transformed),
                                               T(static_cast<long long>(0)),
                                               T(static_cast<long long>(1)),
                                               scaled_eps,
                                               kMaxIntegralDepth);
    }

    T prev_scan_val = evaluate_with_variable(lower_bound);
    for (int i = 1; i <= 40; ++i) {
        const T x =
            lower_bound + (upper_bound - lower_bound) *
                              (T(static_cast<long double>(i)) / T(40.0L));
        T value = evaluate_with_variable(x);
        if (!t_isfinite(value)) {
            throw std::runtime_error("integral did not converge");
        }
        if (t_isfinite(prev_scan_val)) {
            T prev_abs = t_abs(prev_scan_val);
            T curr_abs = t_abs(value);
            if ((prev_scan_val > T(static_cast<long long>(0))) != (value > T(static_cast<long long>(0))) &&
                (prev_abs > T(100.0L) || curr_abs > T(100.0L))) {
                throw std::runtime_error("potential internal discontinuity detected");
            }
        }
        prev_scan_val = value;
    }

    const int coarse_scan_points = 100;
    std::vector<std::pair<T, T>> suspicious_points;
    T prev_x = lower_bound;
    prev_scan_val = evaluate_with_variable(lower_bound);
    for (int i = 1; i <= coarse_scan_points; ++i) {
        const T x = lower_bound + (upper_bound - lower_bound) * T(static_cast<long double>(i)) / T(static_cast<long double>(coarse_scan_points));
        T val;
        try {
            val = evaluate_with_variable(x);
        } catch (...) {
            suspicious_points.push_back({prev_x, x});
            prev_x = x;
            prev_scan_val = T(static_cast<long long>(0));
            continue;
        }

        if (!t_isfinite(val)) {
            suspicious_points.push_back({prev_x, x});
        } else if (t_abs(val) > T(1e6)) {
            suspicious_points.push_back({prev_x, x});
        } else if (t_isfinite(prev_scan_val)) {
            T jump = t_abs(val - prev_scan_val);
            T avg = (t_abs(val) + t_abs(prev_scan_val)) / T(2.0);
            if (avg > T(1e-10L) && jump > T(10.0L) * avg) {
                suspicious_points.push_back({prev_x, x});
            }
        }
        prev_x = x;
        prev_scan_val = val;
    }

    for (const auto& susp : suspicious_points) {
        T l = susp.first, r = susp.second;
        for (int iter = 0; iter < 30 && (r - l) > T(1e-12L); ++iter) {
            T mid = (l + r) * T(0.5L);
            T mid_val;
            try {
                mid_val = evaluate_with_variable(mid);
                if (!t_isfinite(mid_val) || t_abs(mid_val) > T(1e10)) {
                    r = mid;
                } else {
                    T left_mid = (l + mid) * T(0.5L);
                    T right_mid = (mid + r) * T(0.5L);
                    T left_val = evaluate_with_variable(left_mid);
                    T right_val = evaluate_with_variable(right_mid);
                    if (!t_isfinite(left_val) || t_abs(left_val) > T(1e10)) {
                        r = mid;
                    } else if (!t_isfinite(right_val) || t_abs(right_val) > T(1e10)) {
                        l = mid;
                    } else {
                        break;
                    }
                }
            } catch (...) {
                r = mid;
            }
        }
        if ((r - l) < T(1e-6)) {
            throw std::runtime_error("potential internal discontinuity detected near x = " +
                                     format_t((l + r) * T(0.5L)));
        }
    }

    return require_finite_integral(
        adaptive_gauss_kronrod(lower_bound,
                               upper_bound,
                               scaled_eps,
                               kMaxIntegralDepth));
}

template <typename T>
T TFunctionAnalysis<T>::indefinite_integral_at(T x,
                                                T anchor,
                                                T constant) const {
    return constant + definite_integral(anchor, x);
}

template <typename T>
std::vector<TExtremumPoint<T>> TFunctionAnalysis<T>::solve_extrema(T left_bound,
                                                           T right_bound,
                                                           int scan_segments) const {
    if (left_bound >= right_bound) {
        throw std::runtime_error("extrema search requires left_bound < right_bound");
    }
    if (scan_segments < 8) {
        throw std::runtime_error("scan_segments must be at least 8");
    }

    std::vector<TExtremumPoint<T>> extrema;
    T previous_x = left_bound;
    T previous_derivative = derivative(previous_x);

    for (int i = 1; i <= scan_segments; ++i) {
        const T current_x =
            left_bound +
            (right_bound - left_bound) * T(static_cast<long double>(i)) /
                T(static_cast<long double>(scan_segments));
        const T current_derivative = derivative(current_x);

        if (t_is_near_zero(previous_derivative, T(1e-5))) {
            const T stationary_x = previous_x;
            const T second = second_derivative(stationary_x);
            if (!t_is_near_zero(second, T(1e-4))) {
                bool duplicate = false;
                for (const auto& point : extrema) {
                    if (same_extremum_x(point.x, stationary_x)) {
                        duplicate = true;
                        break;
                    }
                }
                if (!duplicate) {
                    extrema.push_back(
                        {stationary_x, evaluate_with_variable(stationary_x), second < T(static_cast<long long>(0))});
                }
            }
        } else if ((previous_derivative < T(static_cast<long long>(0)) && current_derivative > T(static_cast<long long>(0))) ||
                   (previous_derivative > T(static_cast<long long>(0)) && current_derivative < T(static_cast<long long>(0)))) {
            const T stationary_x =
                bisect_stationary_point(previous_x, current_x);
            const T second = second_derivative(stationary_x);
            if (!t_is_near_zero(second, T(1e-4))) {
                extrema.push_back(
                    {stationary_x, evaluate_with_variable(stationary_x), second < T(static_cast<long long>(0))});
            }
        }

        previous_x = current_x;
        previous_derivative = current_derivative;
    }

    std::vector<TExtremumPoint<T>> unique_extrema;
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

template <typename T>
const std::string& TFunctionAnalysis<T>::expression() const {
    return expression_;
}

template <typename T>
const std::string& TFunctionAnalysis<T>::variable_name() const {
    return variable_name_;
}

template <typename T>
T TFunctionAnalysis<T>::evaluate_with_variable(T x) const {
    if (expression_.empty()) {
        throw std::runtime_error("function is not defined");
    }

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

    T value = T(static_cast<long long>(0));
    if (evaluator_) {
        value = evaluator_({{variable_name_, x}});
    } else {
        Calculator calculator;
        calculator.process_line(variable_name_ + " = " + format_t(x), false);
        value = T(calculator.evaluate_raw(expression_));
    }
    
    evaluation_cache_entries_.push_front({cache_key, value});
    evaluation_cache_index_[cache_key] = evaluation_cache_entries_.begin();
    while (evaluation_cache_entries_.size() > kMaxEvaluationCacheSize) {
        evaluation_cache_index_.erase(evaluation_cache_entries_.back().first);
        evaluation_cache_entries_.pop_back();
    }
    return value;
}

template <typename T>
T TFunctionAnalysis<T>::second_derivative(T x) const {
    const T step = scale_aware_step(x);
    const T center = evaluate_with_variable(x);
    const T left_x = x - step;
    const T right_x = x + step;
    const T actual_step = (right_x - left_x) * T(0.5L);
    const T left = evaluate_with_variable(left_x);
    const T right = evaluate_with_variable(right_x);
    const T numerator = compensated_pair_sum(left - center, right - center);
    return numerator / (actual_step * actual_step);
}

template <typename T>
T TFunctionAnalysis<T>::bisect_stationary_point(T left, T right) const {
    T left_derivative = derivative(left);

    for (int i = 0; i < 80; ++i) {
        const T mid = (left + right) * T(0.5L);
        const T mid_derivative = derivative(mid);
        if (t_abs(mid_derivative) <= kRootTolerance_v<T>() ||
            t_abs(right - left) <=
                relative_tolerance(kRootTolerance_v<T>(),
                                   std::max(t_abs(left), t_abs(right)))) {
            return mid;
        }

        if ((left_derivative < T(static_cast<long long>(0)) && mid_derivative > T(static_cast<long long>(0))) ||
            (left_derivative > T(static_cast<long long>(0)) && mid_derivative < T(static_cast<long long>(0)))) {
            right = mid;
        } else {
            left = mid;
            left_derivative = mid_derivative;
        }
    }

    return (left + right) * T(0.5L);
}

// 自适应 Simpson 积分辅助函数
template <typename T>
T TFunctionAnalysis<T>::simpson_rule(T a, T b) const {
    const T h = (b - a) / T(static_cast<long long>(2));
    const T fa = evaluate_with_variable(a);
    const T fb = evaluate_with_variable(b);
    const T fc = evaluate_with_variable((a + b) / T(static_cast<long long>(2)));
    return h / T(static_cast<long long>(3)) * (fa + T(static_cast<long long>(4)) * fc + fb);
}

template <typename T>
T TFunctionAnalysis<T>::adaptive_simpson_precise(T a, T b, T eps, int max_depth) const {
    // 使用自适应 Simpson 方法进行高精度积分
    // 递归实现，当误差足够小或达到最大深度时停止
    int actual_max_depth = max_depth;

    const T c = (a + b) / T(static_cast<long long>(2));
    const T whole = simpson_rule(a, b);
    const T left = simpson_rule(a, c);
    const T right = simpson_rule(c, b);

    return adaptive_simpson_recursive(a, b, whole, left, right, eps, actual_max_depth);
}

template <typename T>
T TFunctionAnalysis<T>::adaptive_simpson_recursive(T a, T b, T whole, T left, T right, T eps, int depth) const {
    const T c = (a + b) / T(static_cast<long long>(2));
    const T combined = left + right;
    const T error = t_abs(combined - whole) / T(static_cast<long long>(15));

    const T scale = std::max(T(static_cast<long long>(1)), t_abs(combined));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        // Richardson 外推提高精度
        return combined + (combined - whole) / T(static_cast<long long>(15));
    }

    // 继续细分
    const T d = (a + c) / T(static_cast<long long>(2));
    const T e = (c + b) / T(static_cast<long long>(2));
    const T left_left = simpson_rule(a, d);
    const T left_right = simpson_rule(d, c);
    const T right_left = simpson_rule(c, e);
    const T right_right = simpson_rule(e, b);

    return adaptive_simpson_recursive(a, c, left, left_left, left_right, eps / T(static_cast<long long>(2)), depth - 1) +
           adaptive_simpson_recursive(c, b, right, right_left, right_right, eps / T(static_cast<long long>(2)), depth - 1);
}

template <typename T>
T TFunctionAnalysis<T>::adaptive_gauss_kronrod(T left,
                                                T right,
                                                T eps,
                                                int max_depth) const {
    T error = T(static_cast<long long>(0));
    const T whole = gauss_kronrod_15(left, right, &error);
    return require_finite_integral(
        adaptive_gauss_kronrod_recursive(left,
                                         right,
                                         eps,
                                         whole,
                                         error,
                                         max_depth));
}

template <typename T>
T TFunctionAnalysis<T>::adaptive_gauss_kronrod_recursive(T left,
                                                          T right,
                                                          T eps,
                                                          T whole,
                                                          T error,
                                                          int depth) const {
    const T scale = std::max(T(static_cast<long long>(1)), t_abs(whole));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        return whole;
    }

    const T mid = (left + right) * T(0.5L);
    T left_error = T(static_cast<long long>(0));
    T right_error = T(static_cast<long long>(0));
    const T left_area = gauss_kronrod_15(left, mid, &left_error);
    const T right_area = gauss_kronrod_15(mid, right, &right_error);
    const T left_result =
        adaptive_gauss_kronrod_recursive(left,
                                         mid,
                                         eps * T(0.5L),
                                         left_area,
                                         left_error,
                                         depth - 1);
    const T right_result =
        adaptive_gauss_kronrod_recursive(mid,
                                         right,
                                         eps * T(0.5L),
                                         right_area,
                                         right_error,
                                         depth - 1);
    return compensated_pair_sum(left_result, right_result);
}

template <typename T>
T TFunctionAnalysis<T>::gauss_kronrod_15(T left,
                                          T right,
                                          T* error_estimate) const {
    static const T kNodes[] = {
        T(0.9914553711208126),
        T(0.9491079123427585),
        T(0.8648644233597691),
        T(0.7415311855993945),
        T(0.5860872354676911),
        T(0.4058451513773972),
        T(0.2077849550078985),
        T(static_cast<long long>(0)),
    };
    static const T kKronrodWeights[] = {
        T(0.02293532201052922),
        T(0.06309209262997855),
        T(0.1047900103222502),
        T(0.1406532597155259),
        T(0.1690047266392679),
        T(0.1903505780647854),
        T(0.2044329400752989),
        T(0.2094821410847278),
    };
    static const T kGaussWeights[] = {
        T(static_cast<long long>(0)),
        T(0.1294849661688697),
        T(static_cast<long long>(0)),
        T(0.2797053914892767),
        T(static_cast<long long>(0)),
        T(0.3818300505051189),
        T(static_cast<long long>(0)),
        T(0.4179591836734694),
    };

    const T center = (left + right) * T(0.5L);
    const T half_width = (right - left) * T(0.5L);
    T kronrod_sum = T(static_cast<long long>(0));
    T gauss_sum = T(static_cast<long long>(0));
    T kronrod_compensation = T(static_cast<long long>(0));
    T gauss_compensation = T(static_cast<long long>(0));

    for (int i = 0; i < 8; ++i) {
        if (t_is_near_zero(kNodes[i], T(static_cast<long long>(0)))) {
            const T value = evaluate_with_variable(center);
            compensated_add(kKronrodWeights[i] * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(kGaussWeights[i] * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const T offset = half_width * kNodes[i];
        const T left_value = evaluate_with_variable(center - offset);
        const T right_value = evaluate_with_variable(center + offset);
        const T pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(kKronrodWeights[i] * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(kGaussWeights[i] * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const T kronrod = half_width * kronrod_sum;
    const T gauss = half_width * gauss_sum;
    *error_estimate = t_abs(kronrod - gauss);
    return kronrod;
}

// 显式模板实例化
template class TFunctionAnalysis<long double>;

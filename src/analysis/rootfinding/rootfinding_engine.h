// ============================================================================
// 求根算法引擎
// ============================================================================
//
// 本文件提供纯粹的数值求根算法实现，不依赖于计算器命令系统。
// 包括：
// - Newton 法（带回溯）
// - 二分法
// - 割线法
// - 不动点迭代
// - Brent 法
//
// 这些算法可以独立使用，便于单元测试。

#ifndef ROOTFINDING_ENGINE_H
#define ROOTFINDING_ENGINE_H

#include "app/scalar_type.h"
#include "precise/precise_decimal.h"
#include "math/core/basic_ops.h"
#include "math/core/roots_powers.h"
#include "analysis/base/precision_constants.h"

#include <functional>
#include <string>
#include <vector>
#include <type_traits>

namespace rootfinding_engine {

using Scalar = mymath::Scalar;

// ============================================================================
// 泛型辅助函数
// ============================================================================

/**
 * @brief 泛型绝对值函数
 */
template <typename T>
T t_abs(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::abs(val);
    } else {
        return mymath::abs(val);
    }
}

/**
 * @brief 泛型平方根函数
 */
template <typename T>
T t_sqrt(const T& val) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::sqrt(val);
    } else {
        return mymath::sqrt(val);
    }
}

template <typename T>
T t_max(const T& a, const T& b) {
    return a < b ? b : a;
}

// ============================================================================
// 内部类型映射
// ============================================================================

template <typename T>
struct InternalType {
    using type = T;
};

template <>
struct InternalType<Scalar> {
    using type = Scalar;
};

template <typename T>
using internal_t = typename InternalType<T>::type;

template <typename T>
internal_t<T> to_internal(T val) {
    return static_cast<internal_t<T>>(val);
}

template <typename T>
T from_internal(internal_t<T> val) {
    return static_cast<T>(val);
}

// 特化：对于 Scalar 类型，保持完整精度，不降级
template <>
inline Scalar from_internal<Scalar>(internal_t<Scalar> val) {
    return val;  // 直接返回，不调用 to_long_double()
}

// ============================================================================
// 容差计算函数
// ============================================================================

/**
 * @brief 计算函数值容差 - 使用精度感知常量
 */
template <typename T>
T root_function_tolerance(T fx) {
    return precision::newton_tolerance<T>() * t_max(T(static_cast<long long>(1)), t_abs(fx));
}

/**
 * @brief 计算位置容差 - 使用精度感知常量
 */
template <typename T>
T root_position_tolerance(T x) {
    return precision::default_absolute_tolerance<T>() * t_max(T(static_cast<long long>(1)), t_abs(x));
}

/**
 * @brief 计算数值导数的步长 - 使用精度感知常量
 */
template <typename T>
T root_derivative_step(T x) {
    return precision::optimal_derivative_step<T>(x);
}

/**
 * @brief 获取求根算法的最大迭代次数
 */
template <typename T>
constexpr int root_max_iterations() {
    if constexpr (std::is_same_v<T, Scalar>) {
        return 300;  // PreciseDecimal 需要更多迭代
    } else {
        return 100;  // double/long double
    }
}

// ============================================================================
// 求根算法
// ============================================================================

/**
 * @brief Newton 法求根
 *
 * 使用 Newton 法求解 f(x) = 0。
 * 如果未提供导数，则使用中心差分近似。
 * 包含回溯（backtracking）以保证收敛。
 *
 * @param evaluate 函数求值器
 * @param initial 初始值
 * @param normalize 结果归一化函数
 * @param evaluate_derivative 导数求值器（可选）
 * @param variable_name 变量名
 * @return 求得的根
 */
template <typename T>
T newton_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize,
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate_derivative = nullptr,
    const std::string& variable_name = "x");

/**
 * @brief 二分法求根
 *
 * 使用二分法求解 f(x) = 0。
 * 要求 f(a) 和 f(b) 异号。
 *
 * @param evaluate 函数求值器
 * @param left 左端点
 * @param right 右端点
 * @param normalize 结果归一化函数
 * @param variable_name 变量名
 * @return 求得的根
 */
template <typename T>
T bisection_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T left,
    T right,
    const std::function<T(T)>& normalize,
    const std::string& variable_name = "x");

/**
 * @brief 割线法求根
 *
 * 使用割线法求解 f(x) = 0。
 * 需要两个初始点 x0 和 x1。
 *
 * @param evaluate 函数求值器
 * @param x0 第一个初始点
 * @param x1 第二个初始点
 * @param normalize 结果归一化函数
 * @param variable_name 变量名
 * @return 求得的根
 */
template <typename T>
T secant_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T x0,
    T x1,
    const std::function<T(T)>& normalize,
    const std::string& variable_name = "x");

/**
 * @brief 不动点迭代
 *
 * 使用不动点迭代求解 x = f(x)。
 *
 * @param evaluate 函数求值器
 * @param initial 初始值
 * @param normalize 结果归一化函数
 * @param variable_name 变量名
 * @return 求得的不动点
 */
template <typename T>
T fixed_point_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T initial,
    const std::function<T(T)>& normalize,
    const std::string& variable_name = "x");

/**
 * @brief Brent 法求根
 *
 * 结合二分法、割线法和逆二次插值的高效求根方法。
 * 具有二分法的稳健性和割线法的快速收敛性。
 * 不需要导数信息，是求解一元方程的推荐方法。
 *
 * @param evaluate 函数求值器
 * @param left 左端点
 * @param right 右端点
 * @param normalize 结果归一化函数
 * @param variable_name 变量名
 * @return 求得的根
 */
template <typename T>
T brent_solve(
    const std::function<T(const std::vector<std::pair<std::string, T>>&)>& evaluate,
    T left,
    T right,
    const std::function<T(T)>& normalize,
    const std::string& variable_name = "x");

}  // namespace rootfinding_engine

#endif  // ROOTFINDING_ENGINE_H
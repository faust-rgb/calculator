// ============================================================================
// 精度感知常量工具
// ============================================================================
//
// 根据 Scalar 类型的精度动态计算最优常量，避免硬编码针对 double 优化的值
// ============================================================================

#ifndef PRECISION_CONSTANTS_H
#define PRECISION_CONSTANTS_H

#include <limits>
#include <cmath>
#include "app/scalar_type.h"
#include "app/default_precision.h"
#include "math/core/basic_ops.h"
#include "math/core/roots_powers.h"
#include "math/precise/precise_decimal.h"

namespace precision {

// ============================================================================
// 基础精度常量
// ============================================================================

/**
 * @brief 获取类型的机器 epsilon
 * 对于 double: ~2.2e-16
 * 对于 long double: ~1.1e-19 (80-bit) 或 ~1.9e-34 (128-bit)
 * 对于 float128_t: ~1.9e-34
 * 对于 PreciseDecimal: 基于 default_scale 动态计算
 */
template <typename T>
inline T epsilon() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        // float128_t 约有 113 位有效数字，epsilon ≈ 2^(-113)
        return T(1.0L) / T(10384593717069655257060992658440192.0L); // 2^113
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        // 基于 default_scale 动态计算，保留 2 位安全余量
        const int scale = app::get_default_scale();
        return PreciseDecimal("1e-" + std::to_string(std::max(scale - 2, 10)));
    } else {
        return std::numeric_limits<T>::epsilon();
    }
}

/**
 * @brief 获取 sqrt(epsilon) - 数值导数的最优步长基数
 * 对于 double: ~1.5e-8
 * 对于 float128_t: ~3e-17
 * 对于 PreciseDecimal: 基于 default_scale 动态计算
 */
template <typename T>
inline T sqrt_epsilon() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::sqrt(epsilon<T>());
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        // 基于 default_scale 动态计算步长基数
        const int scale = app::get_default_scale();
        const int sqrt_scale = std::max(scale / 2, 5);
        return PreciseDecimal("1e-" + std::to_string(sqrt_scale));
    } else {
        return std::sqrt(epsilon<T>());
    }
}

/**
 * @brief 获取 cbrt(epsilon) - 用于某些优化算法
 * 对于 double: ~6e-6
 * 对于 float128_t: ~1e-11
 * 对于 PreciseDecimal: 基于 default_scale 动态计算
 */
template <typename T>
inline T cbrt_epsilon() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        return mymath::cbrt(epsilon<T>());
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        // 基于 default_scale 动态计算
        const int scale = app::get_default_scale();
        const int cbrt_scale = std::max(scale / 3, 3);
        return PreciseDecimal("1e-" + std::to_string(cbrt_scale));
    } else {
        return std::cbrt(epsilon<T>());
    }
}

// ============================================================================
// 数值导数步长
// ============================================================================

/**
 * @brief 计算中心差分数值导数的最优步长
 *
 * 最优步长 h ≈ sqrt(epsilon) * max(1, |x|)
 * 这平衡了截断误差 (O(h^2)) 和舍入误差 (O(epsilon/h))
 *
 * @param x 求导点
 * @return 最优步长
 */
template <typename T>
inline T optimal_derivative_step(T x) {
    const T h_base = sqrt_epsilon<T>();
    const T abs_x = (x < T(0)) ? -x : x;
    return h_base * (abs_x > T(1) ? abs_x : T(1));
}

/**
 * @brief 计算前向差分数值导数的步长（精度稍低但计算量减半）
 */
template <typename T>
inline T forward_derivative_step(T x) {
    // 前向差分误差为 O(h)，需要更大的步长
    const T h_base = cbrt_epsilon<T>();
    const T abs_x = (x < T(0)) ? -x : x;
    return h_base * (abs_x > T(1) ? abs_x : T(1));
}

/**
 * @brief 计算数值 Jacobian 的步长
 */
template <typename T>
inline T jacobian_step(T y) {
    // Jacobian 计算使用中心差分
    return optimal_derivative_step(y);
}

// ============================================================================
// 收敛容差
// ============================================================================

/**
 * @brief 获取默认相对容差
 * 用于积分、ODE 求解等
 */
template <typename T>
inline T default_relative_tolerance() {
    // 相对容差设为 sqrt(epsilon) * 100，提供合理的精度和性能平衡
    return sqrt_epsilon<T>() * T(100);
}

/**
 * @brief 获取默认绝对容差
 */
template <typename T>
inline T default_absolute_tolerance() {
    // 绝对容差设为 epsilon * 10
    return epsilon<T>() * T(10);
}

/**
 * @brief 获取 Newton 迭代的收敛容差
 */
template <typename T>
inline T newton_tolerance() {
    if constexpr (std::is_same_v<T, PreciseDecimal>) {
        // 基于 default_scale 动态计算
        const int scale = app::get_default_scale();
        const int tol_scale = std::max(scale - 4, 8);
        return PreciseDecimal("1e-" + std::to_string(tol_scale));
    } else {
        // Newton 容差设为 epsilon^(1/3) * 10
        return cbrt_epsilon<T>() * T(10);
    }
}

/**
 * @brief 获取事件检测容差
 * 用于 ODE 事件检测
 */
template <typename T>
inline T event_tolerance() {
    // 事件检测容差设为 sqrt(epsilon) * 10
    return sqrt_epsilon<T>() * T(10);
}

/**
 * @brief 获取最小步长
 * @param segment 总区间长度
 */
template <typename T>
inline T min_step_size(T segment) {
    const T abs_segment = (segment < T(0)) ? -segment : segment;
    // 最小步长为 epsilon^(1/4) * segment
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        const T factor = mymath::sqrt(sqrt_epsilon<T>());
        return factor * abs_segment;
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const T factor = mymath::sqrt(sqrt_epsilon<T>());
        return factor * abs_segment;
    } else {
        return std::sqrt(sqrt_epsilon<T>()) * abs_segment;
    }
}

// ============================================================================
// 积分精度常量
// ============================================================================

/**
 * @brief 获取 Gauss-Kronrod 积分的阶数建议
 * 根据精度要求选择合适的积分规则
 */
template <typename T>
inline int recommended_gauss_kronrod_order() {
    if constexpr (std::is_same_v<T, PreciseDecimal>) {
        // 对于 PreciseDecimal，根据当前精度选择阶数
        const int scale = app::get_default_scale();
        if (scale >= 50) return 41;  // 高精度使用 G20-K41
        if (scale >= 30) return 31;  // 中等精度使用 G15-K31
        if (scale >= 20) return 21;  // 较低精度使用 G10-K21
        return 15;                    // 默认 G7-K15
    } else if constexpr (std::is_same_v<T, mymath::float128_t>) {
        // 对于 float128，使用更高阶的 G21-K41 或 G30-K61
        return 41;
    } else if constexpr (std::is_same_v<T, long double> && sizeof(long double) > 8) {
        // 对于 80-bit long double
        return 21;
    } else {
        // 对于 double，G7-K15 足够
        return 15;
    }
}

/**
 * @brief 获取多维积分的默认采样数
 */
template <typename T>
inline int default_integration_samples(int dimensions) {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        // 高精度需要更多采样
        return 100000 / dimensions;
    } else {
        return 10000 / dimensions;
    }
}

// ============================================================================
// 优化算法常量
// ============================================================================

/**
 * @brief 获取梯度下降的收敛阈值
 */
template <typename T>
inline T gradient_convergence_threshold() {
    return sqrt_epsilon<T>() * T(1000);
}

/**
 * @brief 获取 Hessian 矩阵正定性检测阈值
 */
template <typename T>
inline T positive_definite_threshold() {
    return sqrt_epsilon<T>() * T(100);
}

/**
 * @brief 获取线性搜索的最小步长
 */
template <typename T>
inline T line_search_min_step() {
    return epsilon<T>() * T(1000);
}

// ============================================================================
// 矩阵条件数阈值
// ============================================================================

/**
 * @brief 获取矩阵奇异值阈值（用于判断秩亏）
 */
template <typename T>
inline T singular_value_threshold() {
    return sqrt_epsilon<T>() * T(100);
}

/**
 * @brief 获取条件数警告阈值
 */
template <typename T>
inline T condition_number_warning_threshold() {
    return T(1) / sqrt_epsilon<T>();
}

// ============================================================================
// 级数求和收敛阈值
// ============================================================================

/**
 * @brief 获取级数求和的收敛阈值
 * 用于 Taylor 级数、连分式等迭代算法
 */
template <typename T>
inline T series_convergence_threshold() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        // float128 约 35 位精度，使用更严格的阈值
        return T(1e-40L);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const int scale = app::get_default_scale();
        return PreciseDecimal("1e-" + std::to_string(std::max(scale - 3, 12)));
    } else {
        return T(1e-18L);
    }
}

// ============================================================================
// 三角函数阈值
// ============================================================================

/**
 * @brief 获取三角函数奇异点检测阈值
 * 用于 tan, sec, csc, cot 等函数的奇异点检测
 */
template <typename T>
inline T trigonometric_singular_threshold() {
    if constexpr (std::is_same_v<T, mymath::float128_t>) {
        // float128 精度更高，使用更严格的阈值
        return T(1e-17L);
    } else if constexpr (std::is_same_v<T, PreciseDecimal>) {
        const int scale = app::get_default_scale();
        return PreciseDecimal("1e-" + std::to_string(std::max(scale / 2, 8)));
    } else {
        return T(1e-10L);
    }
}

// ============================================================================
// 矩阵算法阈值
// ============================================================================

/**
 * @brief 获取 Jacobi 特征值算法的安全阈值
 * 用于检测非对角元素是否足够小以跳过旋转
 */
template <typename T>
inline T jacobi_safety_threshold() {
    return sqrt_epsilon<T>() * T(1000);
}

/**
 * @brief 获取 Cholesky 分解的正定性检测阈值
 */
template <typename T>
inline T cholesky_positive_threshold() {
    return epsilon<T>() * T(100);
}

} // namespace precision

#endif // PRECISION_CONSTANTS_H

// ============================================================================
// 默认精度配置
// ============================================================================
//
// 统一管理计算器的默认精度设置，包括：
// - PreciseDecimal 的默认小数位数
// - 显示精度阈值
// - 数值计算容差
//
// 这些设置影响所有使用 Scalar 类型的数值计算。
// 所有阈值基于 default_scale 动态计算，适配不同精度类型。
// ============================================================================

#ifndef APP_DEFAULT_PRECISION_H
#define APP_DEFAULT_PRECISION_H

#include <cmath>
#include <algorithm>

namespace app {

// ============================================================================
// PreciseDecimal 默认精度
// ============================================================================

/**
 * @brief 默认小数精度（小数点后位数）
 *
 * 用于 PreciseDecimal 类型的默认精度设置。
 * 可以通过 PrecisionContext::set_default_scale() 在运行时修改。
 *
 * 典型值：
 * - 20: 标准高精度
 * - 40: 超高精度计算
 * - 10: 快速计算（较低精度）
 */
inline int& default_scale() {
    thread_local int scale = 40; // 默认小数位数为 40 位
    return scale;
}

/**
 * @brief 获取默认小数精度
 */
inline int get_default_scale() {
    return default_scale();
}

/**
 * @brief 设置默认小数精度
 */
inline void set_default_scale(int scale) {
    default_scale() = scale;
}

// ============================================================================
// 显示精度阈值
// ============================================================================

/**
 * @brief 默认十进制显示有效位数
 */
constexpr int kDefaultDisplayPrecision = 12;

/**
 * @brief 十进制显示有效位数范围
 */
constexpr int kMinDisplayPrecision = 1;
constexpr int kMaxDisplayPrecision = 50; // Increased from 17 to support high-precision display

/**
 * @brief 获取显示零的阈值
 *
 * 小于此值的数值在显示时会被视为零。
 * 用途：输出格式化时判断是否显示为 0
 */
inline long double display_zero_threshold() {
    return std::pow(10.0L, -std::max(get_default_scale() * 2, 35));
}

/**
 * @brief 获取显示整数的阈值
 *
 * 当数值与最近整数的差小于此值时，显示为整数。
 * 用途：输出格式化时判断是否显示为整数
 */
inline long double display_integer_threshold() {
    return std::pow(10.0L, -std::max(get_default_scale() / 2, 5));
}

// ============================================================================
// 基础数值计算容差
// ============================================================================

/**
 * @brief 获取机器 epsilon 级别的容差
 *
 * 用于需要最严格判断的场合，如判断两个数是否完全相等。
 * 用途：浮点数精确比较、判断是否为零的最严格标准
 */
inline long double epsilon_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 2, 10));
}

/**
 * @brief 获取数值计算的默认容差
 *
 * 用于一般数值计算中的零值判断。
 * 用途：多项式系数零值检测、数值比较（替代 1e-12）
 */
inline long double numeric_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 4, 10));
}

/**
 * @brief 获取整数判断容差
 *
 * 用于判断一个数是否为整数。
 * 用途：is_integer 判断、整数系数检测（替代 1e-10, 1e-9）
 */
inline long double integer_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() / 2, 8));
}

/**
 * @brief 获取宽松容差
 *
 * 用于不需要严格精度的场合，如初步筛选。
 * 用途：多项式除法余数判断、粗略零值检测（替代 1e-7, 1e-6）
 */
inline long double loose_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() / 3, 5));
}

// ============================================================================
// 数值分析容差
// ============================================================================

/**
 * @brief 获取默认相对容差
 *
 * 用于积分、ODE 求解等数值计算的相对误差控制。
 * 用途：自适应积分、ODE 求解器
 */
inline long double default_relative_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() / 2, 8));
}

/**
 * @brief 获取默认绝对容差
 *
 * 用于判断数值是否为零的绝对阈值。
 * 用途：收敛判断、零值检测
 */
inline long double default_absolute_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 2, 10));
}

/**
 * @brief 获取 Newton 迭代的收敛容差
 *
 * 用于 Newton-Raphson 等迭代方法的收敛判断。
 * 用途：求根算法、优化算法
 */
inline long double newton_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 4, 8));
}

/**
 * @brief 获取数值导数步长基数
 *
 * 最优步长 h ≈ sqrt(epsilon) * max(1, |x|)
 * 用途：数值微分、有限差分计算（替代 1e-8）
 */
inline long double derivative_step_base() {
    return std::pow(10.0L, -std::max(get_default_scale() / 2, 5));
}

/**
 * @brief 获取积分容差
 *
 * 用于数值积分的误差控制。
 * 用途：Simpson 积分、Gauss 积分（替代 1e-6）
 */
inline long double integration_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() / 2, 6));
}

/**
 * @brief 获取级数收敛容差
 *
 * 用于判断级数项是否足够小以停止求和。
 * 用途：Taylor 级数、Fourier 级数求和
 */
inline long double series_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 3, 10));
}

// ============================================================================
// 符号计算容差
// ============================================================================

/**
 * @brief 获取符号表达式格式化容差
 *
 * 用于符号表达式输出时的零值判断。
 * 用途：表达式简化、多项式格式化（替代 1e-12）
 */
inline long double format_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 4, 10));
}

/**
 * @brief 获取多项式系数零值容差
 *
 * 用于判断多项式系数是否为零。
 * 用途：多项式运算、GCD 计算（替代 1e-10）
 */
inline long double polynomial_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 4, 10));
}

/**
 * @brief 获取代数数判断容差
 *
 * 用于代数数的相等性判断。
 * 用途：代数数运算、域扩张（替代 1e-15）
 */
inline long double algebraic_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 5, 12));
}

/**
 * @brief 获取极限计算容差
 *
 * 用于极限计算中的数值判断。
 * 用途：极限求解、渐近分析
 */
inline long double limit_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() / 2, 8));
}

// ============================================================================
// 矩阵计算容差
// ============================================================================

/**
 * @brief 获取矩阵计算容差
 *
 * 用于矩阵运算中的零值判断。
 * 用途：矩阵求逆、LU 分解、秩计算（替代 1e-10）
 */
inline long double matrix_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 4, 10));
}

/**
 * @brief 获取矩阵主元相对容差
 *
 * 用于矩阵分解中的主元选择。
 * 用途：LU 分解、QR 分解（替代 1e-13）
 */
inline long double matrix_pivot_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 5, 12));
}

/**
 * @brief 获取特征值计算容差
 *
 * 用于特征值/特征向量计算的收敛判断。
 * 用途：特征值分解、SVD（替代 1e-15）
 */
inline long double eigenvalue_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 5, 12));
}

// ============================================================================
// 统计计算容差
// ============================================================================

/**
 * @brief 获取统计计算容差
 *
 * 用于统计分析中的零值判断。
 * 用途：方差计算、相关性分析（替代 1e-10）
 */
inline long double statistics_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 4, 10));
}

/**
 * @brief 获取概率计算容差
 *
 * 用于概率分布计算。
 * 用途：正态分布、t 分布等
 */
inline long double probability_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() / 2, 8));
}

// ============================================================================
// 特殊函数容差
// ============================================================================

/**
 * @brief 获取特殊函数计算容差
 *
 * 用于特殊函数（Gamma、Bessel 等）的数值计算。
 * 用途：Gamma 函数、Bessel 函数、误差函数
 */
inline long double special_function_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 5, 12));
}

/**
 * @brief 获取级数求和容差
 *
 * 用于级数求和时的项截断判断。
 * 用途：无穷级数求和、Taylor 展开（替代 1e-35）
 */
inline long double summation_tolerance() {
    return std::pow(10.0L, -std::max(get_default_scale() - 5, 15));
}

// ============================================================================
// 精度作用域守卫
// ============================================================================

/**
 * @class PrecisionGuard
 * @brief 临时修改精度的 RAII 守卫
 *
 * 在作用域内临时提高或降低精度，离开作用域时自动恢复。
 */
class PrecisionGuard {
public:
    explicit PrecisionGuard(int new_scale)
        : old_scale_(get_default_scale()) {
        set_default_scale(new_scale);
    }

    ~PrecisionGuard() {
        set_default_scale(old_scale_);
    }

    // 禁止拷贝
    PrecisionGuard(const PrecisionGuard&) = delete;
    PrecisionGuard& operator=(const PrecisionGuard&) = delete;

private:
    int old_scale_;
};

} // namespace app

#endif // APP_DEFAULT_PRECISION_H

// ============================================================================
// Runtime precision configuration
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

#ifndef MATH_RUNTIME_PRECISION_DEFAULT_PRECISION_H
#define MATH_RUNTIME_PRECISION_DEFAULT_PRECISION_H

#include <algorithm>
#include <cmath>

namespace math {
namespace config {

// ============================================================================
// PreciseDecimal 默认精度
// ============================================================================

/**
 * @brief 默认小数精度（小数点后位数）
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

constexpr int kDefaultDisplayPrecision = 12;
constexpr int kMinDisplayPrecision = 1;
constexpr int kMaxDisplayPrecision = 50;

// ============================================================================
// 精度作用域守卫
// ============================================================================

class PrecisionGuard {
public:
    explicit PrecisionGuard(int new_scale)
        : old_scale_(get_default_scale()) {
        set_default_scale(new_scale);
    }

    ~PrecisionGuard() {
        set_default_scale(old_scale_);
    }

    PrecisionGuard(const PrecisionGuard&) = delete;
    PrecisionGuard& operator=(const PrecisionGuard&) = delete;

private:
    int old_scale_;
};

} // namespace config
} // namespace math

namespace math {
namespace config {

// Compatibility wrappers for callers that still use app::* tolerance names.
inline long double display_zero_threshold() { return std::pow(10.0L, -std::max(get_default_scale() * 2, 35)); }
inline long double display_integer_threshold() { return std::pow(10.0L, -std::max(get_default_scale() / 2, 5)); }
inline long double epsilon_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 2, 10)); }
inline long double numeric_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 4, 10)); }
inline long double integer_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() / 2, 8)); }
inline long double loose_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() / 3, 5)); }
inline long double default_relative_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() / 2, 8)); }
inline long double default_absolute_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 2, 10)); }
inline long double newton_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 4, 8)); }
inline long double derivative_step_base() { return std::pow(10.0L, -std::max(get_default_scale() / 2, 5)); }
inline long double integration_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() / 2, 6)); }
inline long double series_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 3, 10)); }
inline long double format_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 4, 10)); }
inline long double polynomial_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 4, 10)); }
inline long double algebraic_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 5, 12)); }
inline long double limit_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() / 2, 8)); }
inline long double matrix_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 4, 10)); }
inline long double matrix_pivot_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 5, 12)); }
inline long double eigenvalue_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 5, 12)); }
inline long double statistics_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 4, 10)); }
inline long double probability_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() / 2, 8)); }
inline long double special_function_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 5, 12)); }
inline long double summation_tolerance() { return std::pow(10.0L, -std::max(get_default_scale() - 5, 15)); }

} // namespace config
} // namespace math

// 为了向后兼容，在 app 命名空间中暴露别名
namespace app {
    using namespace math::config;
}

#endif // MATH_RUNTIME_PRECISION_DEFAULT_PRECISION_H

// ============================================================================
// 表达式工具函数实现
// ============================================================================
//
// 提供级数格式化和数值容差计算功能。
// ============================================================================

#include "expression_utils.h"
#include "core/scalar_type.h"
#include "format_utils.h"
#include "math/mymath.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

// ============================================================================
// 数值容差
// ============================================================================

namespace {
using Scalar = mymath::Scalar;
} // namespace

/**
 * @brief 计算根位置容差
 * @param value 参考值
 * @return 基于参考值计算的容差，保证相对精度
 *
 * 使用自适应容差，对于大值放宽绝对误差要求。
 */
Scalar root_position_tolerance(Scalar value) {
    Scalar v(value);
    return (Scalar(1e-10) * mymath::precise128::fmax(Scalar(1), mymath::precise128::abs(v)));
}

/**
 * @brief 计算根函数容差
 * @param value 参考值
 * @return 基于参考值计算的函数值容差
 */
Scalar root_function_tolerance(Scalar value) {
    Scalar v(value);
    return (Scalar(1e-10) * mymath::precise128::fmax(Scalar(1), mymath::precise128::abs(v)));
}

/**
 * @brief 计算根导数步长
 * @param value 参考值
 * @return 数值微分使用的步长
 */
Scalar root_derivative_step(Scalar value) {
    Scalar v(value);
    return (Scalar(1e-6) * mymath::precise128::fmax(Scalar(1), mymath::precise128::abs(v)));
}

// ============================================================================
// 级数格式化
// ============================================================================

/**
 * @brief 生成移位级数的基表达式
 * @param variable_name 变量名
 * @param center 展开中心点
 * @return 格式化后的基表达式，如 "(x - 1)" 或 "x"
 *
 * 当中心点为零时，直接返回变量名；否则返回 "(变量名 - 中心点)" 形式。
 */
std::string shifted_series_base(const std::string& variable_name, Scalar center) {
    if (mymath::precise128::is_near_zero(center, 1e-12)) {
        return variable_name;
    }
    return "(" + variable_name + signed_center_text(center.to_long_double()) + ")";
}

/**
 * @brief 将广义级数转换为字符串
 * @param coefficients 系数数组
 * @param variable_name 变量名
 * @param center 展开中心点
 * @param denominator 幂次分母（用于 Puiseux 级数等）
 * @return 格式化后的级数字符串
 *
 * 支持任意分母的幂次，例如 denominator=2 时生成形如 "c0 + c1*(x-a)^(1/2) + ..." 的级数。
 * 自动跳过零系数项，并优化输出格式。
 */
std::string generalized_series_to_string(const std::vector<Scalar>& coefficients,
                                         const std::string& variable_name,
                                         Scalar center,
                                         int denominator) {
    if (denominator <= 0) {
        throw std::runtime_error("series denominator must be positive");
    }

    const std::string base = shifted_series_base(variable_name, center);
    std::vector<std::string> terms;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const Scalar coefficient = coefficients[i];
        if (mymath::is_near_zero(coefficient, 1e-12)) {
            continue;
        }
        const std::string factor =
            power_term(base, static_cast<int>(i), denominator);
        terms.push_back(format_term(coefficient.to_long_double(), factor));
    }

    if (terms.empty()) {
        return "0";
    }

    std::ostringstream out;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        if (i == 0) {
            out << terms[i];
        } else if (!terms[i].empty() && terms[i][0] == '-') {
            out << " - " << terms[i].substr(1);
        } else {
            out << " + " << terms[i];
        }
    }
    std::string result = out.str();
    if (result.rfind("1 * ", 0) == 0) {
        result.erase(0, 4);
    }
    return result;
}

/**
 * @brief 将泰勒级数转换为字符串
 * @param coefficients 系数数组
 * @param variable_name 变量名
 * @param center 展开中心点
 * @return 格式化后的泰勒级数字符串
 *
 * 泰勒级数的特化版本，幂次分母固定为 1。
 */
std::string taylor_series_to_string(const std::vector<Scalar>& coefficients,
                                    const std::string& variable_name,
                                    Scalar center) {
    return generalized_series_to_string(coefficients, variable_name, center, 1);
}

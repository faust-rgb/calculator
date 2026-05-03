// ============================================================================
// 表达式工具函数
// ============================================================================
//
// 提供表达式展开和级数处理的函数。
// 表达式拆分函数已移至 parser/unified_expression_parser.h
// 内联函数展开已移至 command/inline_expander.h
// ============================================================================

#ifndef CORE_EXPRESSION_UTILS_H
#define CORE_EXPRESSION_UTILS_H

#include <string>
#include <vector>

// ============================================================================
// 级数工具
// ============================================================================

/**
 * @brief 生成移位级数基
 */
std::string shifted_series_base(const std::string& variable_name, double center);

/**
 * @brief 将广义级数转换为字符串
 */
std::string generalized_series_to_string(const std::vector<double>& coefficients,
                                         const std::string& variable_name,
                                         double center,
                                         int denominator = 1);

/**
 * @brief 将泰勒级数转换为字符串
 */
std::string taylor_series_to_string(const std::vector<double>& coefficients,
                                    const std::string& variable_name,
                                    double center);

// ============================================================================
// 数值容差
// ============================================================================

/**
 * @brief 根据值计算根位置容差
 */
double root_position_tolerance(double value);

/**
 * @brief 根据值计算根函数容差
 */
double root_function_tolerance(double value);

/**
 * @brief 根据值计算根导数步长
 */
double root_derivative_step(double value);

#endif // CORE_EXPRESSION_UTILS_H

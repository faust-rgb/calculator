// ============================================================================
// 表达式工具函数
// ============================================================================
//
// 提供表达式拆分、分析和展开的函数。
// 从 utils.h 拆分以便更好地组织代码。
// ============================================================================

#ifndef CORE_EXPRESSION_UTILS_H
#define CORE_EXPRESSION_UTILS_H

#include <string>
#include <string_view>
#include <vector>

class Calculator;

// ============================================================================
// 表达式拆分
// ============================================================================

/**
 * @brief 拆分赋值表达式 (lhs = rhs)
 */
bool split_assignment(std::string_view expression,
                      std::string_view* lhs,
                      std::string_view* rhs);

/**
 * @brief 拆分命名调用 (name(args))
 */
bool split_named_call(std::string_view expression,
                      std::string* name,
                      std::string* arguments);

bool split_named_call(std::string_view expression,
                      std::string_view name,
                      std::string* arguments);

/**
 * @brief 拆分命名调用并返回参数列表
 */
bool split_named_call_with_arguments(std::string_view expression,
                                     std::string_view name,
                                     std::vector<std::string_view>* arguments);

// ============================================================================
// 内联函数展开
// ============================================================================

/**
 * @brief 检查是否是内联函数命令名
 */
bool is_inline_function_command_name(std::string_view name);

/**
 * @brief 展开内联函数命令
 */
std::string expand_inline_function_commands(Calculator* calculator,
                                            std::string_view expression);

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

// ============================================================================
// 有理数近似
// ============================================================================

/**
 * @brief 尝试将值表示为简单有理数
 */
bool try_make_simple_rational(double value,
                              long long* numerator,
                              long long* denominator);

#endif // CORE_EXPRESSION_UTILS_H

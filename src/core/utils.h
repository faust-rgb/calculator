#ifndef CALCULATOR_UTILS_H
#define CALCULATOR_UTILS_H

// ============================================================================
// 工具函数模块
// ============================================================================
//
// 本头文件提供统一的工具函数接口，并包含以下子模块：
// - string_utils.h: 字符串处理、标识符验证、字符串字面量
// - format_utils.h: 数值和存储值格式化
// - expression_utils.h: 表达式拆分、展开、级数处理
//
// 新代码可以直接包含 utils.h 获取所有功能，
// 或按需包含具体的子模块头文件。
// ============================================================================

#include "core/string_utils.h"
#include "core/format_utils.h"
#include "core/expression_utils.h"

#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <memory>

// 前向声明
struct StoredValue;
struct Rational;
class Calculator;
namespace matrix {
    class Matrix;
}

namespace utils {

/**
 * @brief 去除字符串首尾空白
 * @param text 输入字符串
 * @return 处理后的字符串视图
 */
std::string_view trim_view(std::string_view text);
std::string trim_copy(std::string_view text);

/**
 * @brief 检查字符串是否为合法的标识符
 * @param name 字符串
 * @return 是否合法
 */
bool is_valid_identifier(std::string_view name);

} // namespace utils

// ============================================================================
// 全局命名空间辅助函数（实现于 utils.cpp）
// ============================================================================

/** @brief 去除字符串首尾空白（全局版本） */
std::string_view trim_view(std::string_view text);
std::string trim_copy(std::string_view text);

/** @brief 检查是否为保留函数名 */
bool is_reserved_function_name(std::string_view name);

/** @brief 分割函数定义表达式 */
bool split_function_definition(std::string_view expression,
                               std::string_view* function_name,
                               std::string_view* parameter_name,
                               std::string_view* body);

// 进制转换
/** @brief 转换数值进制为字符串 */
bool convert_base_value(long long value,
                        int base,
                        const struct HexFormatOptions& hex_options,
                        std::string* output);

/** @brief 尝试在表达式中执行进制转换（如 bin(10)） */
bool try_base_conversion_expression(std::string_view expression,
                                    const class VariableResolver& variables,
                                    const std::map<std::string, struct CustomFunction>* functions,
                                    const struct HexFormatOptions& hex_options,
                                    std::string* output);

// 值格式化
/** @brief 规范化显示数值 */
double normalize_display_decimal(double value);

/** @brief 格式化小数 */
std::string format_decimal(double value);

/** @brief 按指定有效位数格式化小数 */
std::string format_decimal(double value, int precision);

/** @brief 格式化存储的值为字符串 */
std::string format_stored_value(const StoredValue& value, bool symbolic_constants_mode);

/** @brief 格式化显示值 */
std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode);

/** @brief 格式化符号数值 */
std::string format_symbolic_number(double value);

/** @brief 格式化符号标量 */
std::string format_symbolic_scalar(double value);

/** @brief 设置进程级显示有效位数 */
void set_process_display_precision(int precision);

/** @brief 查询进程级显示有效位数 */
int process_display_precision();

// 常量与有理数
/** @brief 查找内置常量 */
bool lookup_builtin_constant(const std::string& name, double* value);

/** @brief 尝试将 double 转换为简单的 Rational */
bool try_make_simple_rational(double value, int max_denominator, Rational* rational);

// 根求解容差
double root_position_tolerance(double value);
double root_function_tolerance(double value);
double root_derivative_step(double value);

// 字符串与标识符
bool is_valid_variable_name(std::string_view name);
bool is_identifier_text(std::string_view text);
bool is_string_literal(std::string_view text);
std::string parse_string_literal_value(std::string_view text);
std::string decode_escaped_string(std::string_view text);
std::string encode_state_field(const std::string& text);
std::string decode_state_field(const std::string& text);

// ============================================================================
// 级数格式化
// ============================================================================

std::string shifted_series_base(const std::string& variable_name, double center);

/** @brief 将广义级数系数转换为字符串 */
std::string generalized_series_to_string(const std::vector<double>& coefficients,
                                         const std::string& variable_name,
                                         double center,
                                         int denominator);

/** @brief 将泰勒级数系数转换为字符串 */
std::string taylor_series_to_string(const std::vector<double>& coefficients,
                                    const std::string& variable_name,
                                    double center);

#endif // CALCULATOR_UTILS_H

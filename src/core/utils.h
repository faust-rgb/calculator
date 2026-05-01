#ifndef CALCULATOR_UTILS_H
#define CALCULATOR_UTILS_H

#include <string>
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
 * @return 处理后的字符串
 */
std::string trim_copy(const std::string& text);

/**
 * @brief 检查字符串是否为合法的标识符
 * @param name 字符串
 * @return 是否合法
 */
bool is_valid_identifier(const std::string& name);

} // namespace utils

// ============================================================================
// 全局命名空间辅助函数（实现于 utils.cpp）
// ============================================================================

/** @brief 去除字符串首尾空白（全局版本） */
std::string trim_copy(const std::string& text);

/** @brief 检查是否为保留函数名 */
bool is_reserved_function_name(const std::string& name);

/** @brief 分割函数定义表达式 */
bool split_function_definition(const std::string& expression,
                               std::string* function_name,
                               std::string* parameter_name,
                               std::string* body);

// 进制转换
/** @brief 转换数值进制为字符串 */
bool convert_base_value(long long value,
                        int base,
                        const struct HexFormatOptions& hex_options,
                        std::string* output);

/** @brief 尝试在表达式中执行进制转换（如 bin(10)） */
bool try_base_conversion_expression(const std::string& expression,
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
bool is_valid_variable_name(const std::string& name);
bool is_identifier_text(const std::string& text);
bool is_string_literal(const std::string& text);
std::string parse_string_literal_value(const std::string& text);
std::string decode_escaped_string(const std::string& text);
std::string encode_state_field(const std::string& text);
std::string decode_state_field(const std::string& text);

// 表达式分割
bool split_assignment(const std::string& expression, std::string* lhs, std::string* rhs);
bool split_named_call(const std::string& expression, const std::string& name, std::string* inside);
bool split_named_call_with_arguments(const std::string& expression, const std::string& name, std::vector<std::string>* arguments);
std::vector<std::string> split_top_level_arguments(const std::string& text);

// 指令展开
bool is_inline_function_command_name(const std::string& name);
std::size_t find_matching_paren(const std::string& text, std::size_t open_pos);
std::string expand_inline_function_commands(Calculator* calculator, const std::string& expression);

// 级数格式化
/** @brief 获取级数基底文本 */
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

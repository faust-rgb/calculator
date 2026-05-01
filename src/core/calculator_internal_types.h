// ============================================================================
// 计算器内部类型定义
// ============================================================================
//
// 本文件定义 Calculator 类内部使用的所有数据结构和辅助函数。
// 这些类型不对外暴露，允许在不修改公共 API 的情况下更改实现。
//
// 主要组件：
// 1. StoredValue - 存储的值（标量、矩阵、字符串等）
// 2. CustomFunction/ScriptFunction - 用户定义的函数
// 3. Calculator::Impl - Pimpl 模式的实现类
// 4. ScriptSignal - 脚本执行控制流信号
// 5. 辅助函数声明
// ============================================================================

#ifndef CALCULATOR_INTERNAL_TYPES_H
#define CALCULATOR_INTERNAL_TYPES_H

#include "calculator.h"

#include "matrix.h"
#include "mymath.h"
#include "script_ast.h"
#include "types/rational.h"
#include "types/precise_decimal.h"
#include "types/stored_value.h"

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include "statistics/statistics.h"
#include "statistics/probability.h"

#include "calculator_exceptions.h"
#include "utils.h"

// ============================================================================
// 显示精度常量
// ============================================================================

/** @brief 判断数值是否为零的显示阈值 */
constexpr double kDisplayZeroEps = mymath::kDoubleDenormMin;

/** @brief 判断数值是否为整数的显示阈值 */
constexpr double kDisplayIntegerEps = 1e-9;

/** @brief 默认十进制显示有效位数 */
constexpr int kDefaultDisplayPrecision = 12;

/** @brief 十进制显示有效位数范围 */
constexpr int kMinDisplayPrecision = 1;
constexpr int kMaxDisplayPrecision = 17;

// ============================================================================
// 存储值类型（定义在 types/stored_value.h）
// ============================================================================

// StoredValue 结构体已移至 types/stored_value.h

// ============================================================================
// 函数类型
// ============================================================================

#include "variable_resolver.h"
#include "../parser/decimal_parser.h"
#include "../parser/exact_parser.h"
#include "../parser/symbolic_render_parser.h"

class CalculatorModule;

// ============================================================================
// Calculator 实现类
// ============================================================================

/**
 * @struct Calculator::Impl
 * @brief Calculator 的内部实现
 *
 * 存储计算器的所有状态：
 * - 变量表（全局和局部作用域）
 * - 函数表（简单函数和脚本函数）
 * - 显示选项
 */
struct Calculator::Impl {
    std::map<std::string, StoredValue> variables;          ///< 全局变量
    std::map<std::string, CustomFunction> functions;       ///< 简单函数
    std::map<std::string, ScriptFunction> script_functions; ///< 脚本函数
    std::vector<std::map<std::string, StoredValue>> local_scopes; ///< 局部作用域栈

    std::vector<std::shared_ptr<CalculatorModule>> registered_modules; ///< 已注册的数学模块
    std::vector<std::shared_ptr<CalculatorModule>> implicit_evaluation_modules; ///< 优化后的隐式求值模块列表

    std::map<std::string, std::function<double(const std::vector<double>&)>> scalar_functions; ///< 汇总的标量函数
    std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> matrix_functions; ///< 汇总的矩阵函数
    std::map<std::string, matrix::ValueFunction> value_functions; ///< 汇总的值多态函数

    // 汇总的模块元数据，用于补全和帮助
    std::vector<std::string> module_commands;
    std::vector<std::string> module_functions;
    std::map<std::string, std::vector<std::shared_ptr<CalculatorModule>>> help_topic_to_modules;
    std::map<std::string, std::shared_ptr<CalculatorModule>> command_to_module;

    bool symbolic_constants_mode = false;  ///< 符号常量模式（pi, e 保留符号形式）
    bool hex_prefix_mode = false;          ///< 十六进制输出前缀
    bool hex_uppercase_mode = true;        ///< 十六进制大写字母
    int display_precision = kDefaultDisplayPrecision; ///< 十进制显示有效位数
    int script_call_depth = 0;             ///< 脚本递归深度计数器
};

#include "../script/script_signal.h"

// ============================================================================
// 格式化选项
// ============================================================================

/**
 * @struct HexFormatOptions
 * @brief 十六进制格式化选项
 */
struct HexFormatOptions {
    bool prefix = false;    ///< 是否添加 0x 前缀
    bool uppercase = true;  ///< 是否使用大写字母
};

// ============================================================================
// 辅助函数声明
// ============================================================================

#include "math/helpers/integer_helpers.h"
#include "math/helpers/combinatorics.h"
#include "math/helpers/bitwise_helpers.h"
#include "math/helpers/unit_conversions.h"
#include "math/helpers/base_conversions.h"

// 显示格式化
void set_process_display_precision(int precision);    ///< 设置进程级显示有效位数
int process_display_precision();                      ///< 查询进程级显示有效位数

// 内置常量

// 有理数运算（实现在 precise/rational.cpp）
Rational pow_rational(Rational base, long long exponent); ///< 有理数幂
Rational abs_rational(Rational value);  ///< 有理数绝对值
double rational_to_double(const Rational& value); ///< 有理数转 double

// 精确小数运算（实现在 precise/precise_decimal.cpp）
PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal pow_precise_decimal(const PreciseDecimal& base, long long exponent);
int compare_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

// 字符串处理
bool is_valid_variable_name(const std::string& name);  ///< 检查变量名合法性
bool is_identifier_text(const std::string& text);      ///< 检查是否为标识符
bool is_string_literal(const std::string& text);       ///< 检查是否为字符串字面量
std::string parse_string_literal_value(const std::string& text); ///< 解析字符串字面量

// 状态持久化
std::string encode_state_field(const std::string& text); ///< 编码状态字段
std::string decode_state_field(const std::string& text); ///< 解码状态字段

// 表达式分割
bool split_assignment(const std::string& expression, std::string* lhs, std::string* rhs);
bool split_named_call(const std::string& expression, const std::string& name, std::string* inside);
bool split_named_call_with_arguments(const std::string& expression, const std::string& name, std::vector<std::string>* arguments);
std::vector<std::string> split_top_level_arguments(const std::string& text);

// 函数展开
std::string expand_inline_function_commands(Calculator* calculator, const std::string& expression);

// 精确小数处理（实现在 precise/precise_decimal.cpp 和 precise/precise_parser.cpp）
// stored_value_precise_decimal_text 和 parse_precise_decimal_expression
// 已移至 types/stored_value.h 和 precise/precise_parser.h

// 值格式化
void apply_calculator_display_precision(const Calculator::Impl* impl);

// 根求解容差

// 级数展开

// 线性方程组
std::vector<double> solve_dense_linear_system(std::vector<std::vector<double>> matrix, std::vector<double> rhs, const std::string& context);

// 函数定义解析

#include "../script/script_runtime.h"

// 模块注册
void register_standard_modules(Calculator* calculator);

#endif

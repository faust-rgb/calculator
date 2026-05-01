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

/**
 * @struct CustomFunction
 * @brief 单参数自定义函数
 *
 * 存储用户定义的简单函数，如 f(x) = x^2 + 1。
 * 参数名和表达式以字符串形式存储。
 */
struct CustomFunction {
    std::string parameter_name;  ///< 参数名
    std::string expression;      ///< 函数体表达式
};

/**
 * @struct ScriptFunction
 * @brief 脚本定义的函数
 *
 * 存储通过脚本语言定义的函数，支持多参数和复杂控制流。
 * 函数体为 AST 形式，可高效执行。
 */
struct ScriptFunction {
    std::vector<std::string> parameter_names;              ///< 参数名列表
    std::shared_ptr<const script::BlockStatement> body;    ///< 函数体 AST
};

// ============================================================================
// 回调类型定义
// ============================================================================

/** @brief 检查脚本函数是否存在的回调类型 */
using HasScriptFunctionCallback = std::function<bool(const std::string&)>;

/** @brief 调用脚本函数的回调类型 */
using InvokeScriptFunctionDecimalCallback =
    std::function<double(const std::string&, const std::vector<double>&)>;

/**
 * @class VariableResolver
 * @brief 变量查找解析器，避免 Map 复制
 *
 * 在脚本执行或表达式求值时，用于在局部作用域栈、全局变量表和内置常量中
 * 按优先级查找变量。通过持有引用或指针，避免了 visible_variables() 的 Map 复制开销。
 */
class VariableResolver {
public:
    VariableResolver() : global_vars_(nullptr), local_scopes_(nullptr), override_vars_(nullptr), parent_(nullptr), is_owned_(false) {}
    VariableResolver(const std::map<std::string, StoredValue>* global_vars,
                     const std::vector<std::map<std::string, StoredValue>>* local_scopes,
                     const std::map<std::string, StoredValue>* override_vars = nullptr,
                     const VariableResolver* parent = nullptr)
        : global_vars_(global_vars), local_scopes_(local_scopes), override_vars_(override_vars), parent_(parent), is_owned_(false) {}

    /** @brief 构造一个拥有所有权副本的解析器，用于 Lambda 捕获 */
    static VariableResolver make_owned(const VariableResolver& other);

    /** @brief 查找变量值，返回指针，未找到返回 nullptr */
    const StoredValue* lookup(const std::string& name) const;

    /** @brief 检查变量是否存在（包括内置常量） */
    bool contains(const std::string& name) const;

    /** @brief 获取所有可见变量的快照（用于兼容旧代码或显示） */
    std::map<std::string, StoredValue> snapshot() const;

private:
    const std::map<std::string, StoredValue>* global_vars_;
    const std::vector<std::map<std::string, StoredValue>>* local_scopes_;
    const std::map<std::string, StoredValue>* override_vars_;
    const VariableResolver* parent_;
    
    bool is_owned_;
    std::shared_ptr<std::map<std::string, StoredValue>> owned_global_vars_;
    std::shared_ptr<std::vector<std::map<std::string, StoredValue>>> owned_local_scopes_;
    std::shared_ptr<VariableResolver> owned_parent_;
};

// ============================================================================
// 数值表达式解析器
// ============================================================================

/**
 * @class DecimalParser
 * @brief 数值表达式解析器
 *
 * 解析并计算数值表达式，支持：
 * - 变量引用
 * - 函数调用（内置和自定义）
 * - 运算符优先级
 * - 括号分组
 */
class DecimalParser {
public:
    using ScalarFunction = std::function<double(const std::vector<double>&)>;

    DecimalParser(std::string source,
                  const VariableResolver& variables,
                  const std::map<std::string, CustomFunction>* functions,
                  const std::map<std::string, ScalarFunction>* scalar_functions = nullptr,
                  HasScriptFunctionCallback has_script_function = {},
                  InvokeScriptFunctionDecimalCallback invoke_script_function = {});

    double parse();

private:
    std::string source_;
    VariableResolver variables_;
    const std::map<std::string, CustomFunction>* functions_;
    const std::map<std::string, ScalarFunction>* scalar_functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionDecimalCallback invoke_script_function_;
};

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

    std::map<std::string, std::function<double(const std::vector<double>&)>> scalar_functions; ///< 汇总的标量函数
    std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> matrix_functions; ///< 汇总的矩阵函数
    std::map<std::string, matrix::ValueFunction> value_functions; ///< 汇总的值多态函数

    bool symbolic_constants_mode = false;  ///< 符号常量模式（pi, e 保留符号形式）
    bool hex_prefix_mode = false;          ///< 十六进制输出前缀
    bool hex_uppercase_mode = true;        ///< 十六进制大写字母
    int display_precision = kDefaultDisplayPrecision; ///< 十进制显示有效位数
};

// ============================================================================
// 脚本控制流信号
// ============================================================================

/**
 * @struct ScriptSignal
 * @brief 脚本执行的控制流信号
 *
 * 用于实现 return、break、continue 语句。
 * 语句执行返回信号，外层控制结构根据信号类型决定行为。
 */
struct ScriptSignal {
    enum class Kind {
        kNone,      ///< 无信号，正常执行
        kReturn,    ///< return 语句
        kBreak,     ///< break 语句
        kContinue,  ///< continue 语句
    };

    Kind kind = Kind::kNone;  ///< 信号类型
    bool has_value = false;   ///< 是否有返回值
    StoredValue value;        ///< 返回值（return 时使用）

    /** @brief 创建 return 信号 */
    static ScriptSignal make_return(const StoredValue& return_value);

    /** @brief 创建 break 信号 */
    static ScriptSignal make_break();

    /** @brief 创建 continue 信号 */
    static ScriptSignal make_continue();
};

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

// 整数运算
long long gcd_ll(long long a, long long b);           ///< 最大公约数
long long lcm_ll(long long a, long long b);           ///< 最小公倍数
bool is_integer_double(double x, double eps = 1e-10); ///< 检查 double 是否为整数
long long round_to_long_long(double x);               ///< 四舍五入
long long trunc_to_long_long(double x);               ///< 向零截断
long long floor_to_long_long(double x);               ///< 向下取整
long long ceil_to_long_long(double x);                ///< 向上取整

// 显示格式化
double normalize_display_decimal(double value);       ///< 规范化显示数值
std::string format_decimal(double value);             ///< 格式化小数
std::string format_decimal(double value, int precision); ///< 按指定有效位数格式化小数
std::string format_symbolic_number(double value);     ///< 格式化符号数值
void set_process_display_precision(int precision);    ///< 设置进程级显示有效位数
int process_display_precision();                      ///< 查询进程级显示有效位数

// 内置常量
bool lookup_builtin_constant(const std::string& name, double* value); ///< 查找内置常量

// 单位转换
double degrees_to_radians(double value);    ///< 角度转弧度
double radians_to_degrees(double value);    ///< 弧度转角度
double celsius_to_fahrenheit(double value); ///< 摄氏转华氏
double fahrenheit_to_celsius(double value); ///< 摄氏转华氏


// 数论函数
bool is_prime_ll(long long value);         ///< 素数判断
long long next_prime_ll(long long value);  ///< 下一个素数
long long prev_prime_ll(long long value);  ///< 上一个素数
long long euler_phi_ll(long long value);   ///< 欧拉函数
long long mobius_ll(long long value);      ///< 莫比乌斯函数
long long prime_pi_ll(long long value);    ///< 素数计数函数
long long extended_gcd_ll(long long a, long long b, long long* x, long long* y); ///< 扩展欧几里得

// 组合数学
double fibonacci_value(long long n);       ///< 斐波那契数
double factorial_value(long long n);       ///< 阶乘
Rational factorial_rational(long long n);  ///< 阶乘（有理数形式）
double combination_value(long long n, long long r);  ///< 组合数 C(n,r)
Rational combination_rational(long long n, long long r);
double permutation_value(long long n, long long r);  ///< 排列数 P(n,r)
Rational permutation_rational(long long n, long long r);

// 位运算
std::uint64_t to_unsigned_bits(long long value);      ///< 转无符号位表示
long long from_unsigned_bits(std::uint64_t value);    ///< 从无符号位表示转换
unsigned normalize_rotation_count(long long count);   ///< 规范化旋转次数
std::uint64_t rotate_left_bits(std::uint64_t value, unsigned count);  ///< 循环左移
std::uint64_t rotate_right_bits(std::uint64_t value, unsigned count); ///< 循环右移
int popcount_bits(std::uint64_t value);    ///< 1 的个数
int bit_length_bits(std::uint64_t value);  ///< 位长度
int trailing_zero_count_bits(std::uint64_t value); ///< 末尾零个数
int leading_zero_count_bits(std::uint64_t value);  ///< 前导零个数
int parity_bits(std::uint64_t value);      ///< 奇偶校验
std::uint64_t reverse_bits(std::uint64_t value);   ///< 位反转

// 因式分解
std::string factor_integer(long long value); ///< 整数因式分解

// 进制转换
int digit_value(char ch);                   ///< 字符转数字值
bool prefixed_base(char prefix, int* base); ///< 判断前缀对应的进制
long long parse_prefixed_integer_token(const std::string& token); ///< 解析带前缀整数

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
std::string trim_copy(const std::string& text);        ///< 去除首尾空白
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
std::string format_stored_value(const StoredValue& value, bool symbolic_constants_mode);
std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode);
std::string format_symbolic_scalar(double value);

// 根求解容差
double root_position_tolerance(double value);  ///< 位置容差
double root_function_tolerance(double value);   ///< 函数值容差
double root_derivative_step(double value);      ///< 导数步长

// 级数展开
std::string taylor_series_to_string(const std::vector<double>& coefficients, const std::string& variable_name, double center);
std::string shifted_series_base(const std::string& variable_name, double center);
std::string generalized_series_to_string(const std::vector<double>& coefficients, const std::string& variable_name, double center, int denominator);

// 线性方程组
std::vector<double> solve_dense_linear_system(std::vector<std::vector<double>> matrix, std::vector<double> rhs, const std::string& context);

// 函数定义解析
bool is_reserved_function_name(const std::string& name);
bool split_function_definition(const std::string& expression, std::string* function_name, std::string* parameter_name, std::string* body);

// 表达式求值
double parse_decimal_expression(const std::string& expression, const VariableResolver& variables, const std::map<std::string, CustomFunction>* functions, const std::map<std::string, DecimalParser::ScalarFunction>* scalar_functions = nullptr, HasScriptFunctionCallback has_script_function = {}, InvokeScriptFunctionDecimalCallback invoke_script_function = {});

bool try_evaluate_matrix_expression(const std::string& expression, const VariableResolver& variables, const std::map<std::string, CustomFunction>* functions, const std::map<std::string, DecimalParser::ScalarFunction>* scalar_functions, const std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>>* matrix_functions, const std::map<std::string, matrix::ValueFunction>* value_functions, const HasScriptFunctionCallback& has_script_function, const InvokeScriptFunctionDecimalCallback& invoke_script_function, matrix::Value* value);

Rational parse_exact_expression(const std::string& expression, const VariableResolver& variables, const std::map<std::string, CustomFunction>* functions, HasScriptFunctionCallback has_script_function = {});

bool try_base_conversion_expression(const std::string& expression, const VariableResolver& variables, const std::map<std::string, CustomFunction>* functions, const HexFormatOptions& hex_options, std::string* output);

std::string matrix_literal_expression(const matrix::Matrix& value);

bool try_symbolic_constant_expression(const std::string& expression, const VariableResolver& variables, const std::map<std::string, CustomFunction>* functions, std::string* output);

// 变量作用域
VariableResolver visible_variables(const Calculator::Impl* impl);
bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name);
void assign_visible_variable(Calculator::Impl* impl, const std::string& name, const StoredValue& value);

// 表达式求值（返回 StoredValue）- 新接口使用 ExpressionCache
struct ExpressionCache;
StoredValue evaluate_expression_value(Calculator* calculator, Calculator::Impl* impl, const std::string& expression, bool exact_mode, std::shared_ptr<ExpressionCache>* cache = nullptr);

// 脚本函数调用
double invoke_script_function_decimal(Calculator* calculator, Calculator::Impl* impl, const std::string& name, const std::vector<double>& arguments);

// 脚本执行
ScriptSignal execute_script_statement(Calculator* calculator, Calculator::Impl* impl, const script::Statement& statement, bool exact_mode, std::string* last_output, bool create_scope);
std::string render_script_block(const script::BlockStatement& block, int indent);

// 模块注册
void register_standard_modules(Calculator* calculator, Calculator::Impl* impl);

#endif

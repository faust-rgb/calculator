// ============================================================================
// 计算器内部类型定义
// ============================================================================
//
// 本文件定义 Calculator 类内部使用的所有数据结构和辅助函数。
// 这些类型不对外暴露，允许在不修改公共 API 的情况下更改实现。
//
// 主要组件：
// 1. Rational - 有理数表示，用于精确模式
// 2. PreciseDecimal - 精确小数表示，避免浮点误差
// 3. StoredValue - 存储的值（标量、矩阵、字符串等）
// 4. CustomFunction/ScriptFunction - 用户定义的函数
// 5. Calculator::Impl - Pimpl 模式的实现类
// 6. ScriptSignal - 脚本执行控制流信号
// 7. 辅助函数声明
// ============================================================================

#ifndef CALCULATOR_INTERNAL_TYPES_H
#define CALCULATOR_INTERNAL_TYPES_H

#include "calculator.h"

#include "matrix.h"
#include "mymath.h"
#include "script_ast.h"

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

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
// 异常类型
// ============================================================================

/**
 * @class ExactModeUnsupported
 * @brief 精确模式不支持异常
 *
 * 当表达式无法在精确模式下计算时抛出。
 */
class ExactModeUnsupported : public std::runtime_error {
public:
    explicit ExactModeUnsupported(const std::string& message)
        : std::runtime_error(message) {}
};

/**
 * @class PreciseDecimalUnsupported
 * @brief 精确小数不支持异常
 *
 * 当表达式无法用精确小数表示时抛出。
 */
class PreciseDecimalUnsupported : public std::runtime_error {
public:
    explicit PreciseDecimalUnsupported(const std::string& message)
        : std::runtime_error(message) {}
};

// ============================================================================
// 有理数类型
// ============================================================================

/**
 * @struct Rational
 * @brief 有理数表示
 *
 * 用于精确模式下的分数运算，避免浮点误差。
 * 自动规范化为最简分数（分母为正）。
 *
 * 例如：1/3 + 1/6 = 1/2（精确计算）
 */
struct Rational {
    long long numerator = 0;    ///< 分子
    long long denominator = 1;  ///< 分母（始终为正）

    Rational() = default;
    Rational(long long num, long long den);

    /** @brief 规范化为最简分数 */
    void normalize();

    /** @brief 检查是否为整数 */
    bool is_integer() const;

    /** @brief 转换为字符串，如 "1/2" 或 "3" */
    std::string to_string() const;
};

Rational operator+(const Rational& lhs, const Rational& rhs);
Rational operator-(const Rational& lhs, const Rational& rhs);
Rational operator*(const Rational& lhs, const Rational& rhs);
Rational operator/(const Rational& lhs, const Rational& rhs);

// ============================================================================
// 精确小数类型
// ============================================================================

/**
 * @struct PreciseDecimal
 * @brief 精确小数表示
 *
 * 使用字符串存储数字，避免浮点误差。
 * 适用于需要精确表示小数的场景，如货币计算。
 *
 * 内部表示：digits 存储有效数字，scale 表示小数点位置。
 * 例如：123.45 → digits="12345", scale=2
 */
struct PreciseDecimal {
    std::string digits = "0";  ///< 有效数字字符串
    int scale = 0;             ///< 小数点后的位数
    bool negative = false;     ///< 是否为负数

    /** @brief 规范化表示（去除前导零、末尾零） */
    void normalize();

    /** @brief 检查是否为零 */
    bool is_zero() const;

    /** @brief 转换为字符串 */
    std::string to_string() const;

    /** @brief 转换为 double（可能有精度损失） */
    double to_double() const;

    /** @brief 从原始数字构造 */
    static PreciseDecimal from_digits(std::string raw_digits,
                                      int raw_scale,
                                      bool is_negative);

    /** @brief 从整数字符串构造 */
    static PreciseDecimal from_integer_string(const std::string& integer_text,
                                              bool is_negative);

    /** @brief 从小数字面量构造 */
    static PreciseDecimal from_decimal_literal(const std::string& token);
};

// ============================================================================
// 存储值类型
// ============================================================================

/**
 * @struct StoredValue
 * @brief 计算器中存储的值
 *
 * 支持多种值类型：
 * - 标量（double 或有理数）
 * - 矩阵
 * - 字符串
 * - 符号表达式文本
 * - 精确小数文本
 *
 * 通过标志位区分值类型，实现统一的存储接口。
 */
struct StoredValue {
    bool is_matrix = false;              ///< 是否为矩阵
    bool is_complex = false;             ///< 是否为复数标量
    bool is_string = false;              ///< 是否为字符串
    bool has_symbolic_text = false;      ///< 是否有符号表达式文本
    bool has_precise_decimal_text = false; ///< 是否有精确小数文本
    bool exact = false;                  ///< 是否在精确模式下创建

    Rational rational;                   ///< 有理数值
    double decimal = 0.0;                ///< 浮点数值
    matrix::ComplexNumber complex;       ///< 复数值
    std::string string_value;            ///< 字符串值
    std::string symbolic_text;           ///< 符号表达式文本
    std::string precise_decimal_text;    ///< 精确小数文本
    matrix::Matrix matrix;               ///< 矩阵值
};

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
    DecimalParser(std::string source,
                  const std::map<std::string, StoredValue>* variables,
                  const std::map<std::string, CustomFunction>* functions,
                  HasScriptFunctionCallback has_script_function = {},
                  InvokeScriptFunctionDecimalCallback invoke_script_function = {});

    double parse();

private:
    std::string source_;
    const std::map<std::string, StoredValue>* variables_;
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionDecimalCallback invoke_script_function_;
};

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

// 随机数
std::mt19937_64& global_rng();  ///< 全局随机数生成器

// 内置常量
bool lookup_builtin_constant(const std::string& name, double* value); ///< 查找内置常量

// 单位转换
double degrees_to_radians(double value);    ///< 角度转弧度
double radians_to_degrees(double value);    ///< 弧度转角度
double celsius_to_fahrenheit(double value); ///< 摄氏转华氏
double fahrenheit_to_celsius(double value); ///< 华氏转摄氏

// 统计函数
double normal_pdf(double x, double mean, double sigma); ///< 正态分布 PDF
double normal_cdf(double x, double mean, double sigma); ///< 正态分布 CDF

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

// 有理数运算
Rational pow_rational(Rational base, long long exponent); ///< 有理数幂
Rational abs_rational(Rational value);  ///< 有理数绝对值
double rational_to_double(const Rational& value); ///< 有理数转 double

// 精确小数运算
PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);
PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs);

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

// 精确小数处理
std::string stored_value_precise_decimal_text(const StoredValue& value);
PreciseDecimal parse_precise_decimal_expression(const std::string& expression, const std::map<std::string, StoredValue>* variables);

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
double parse_decimal_expression(const std::string& expression, const std::map<std::string, StoredValue>* variables, const std::map<std::string, CustomFunction>* functions, HasScriptFunctionCallback has_script_function = {}, InvokeScriptFunctionDecimalCallback invoke_script_function = {});

bool try_evaluate_matrix_expression(const std::string& expression, const std::map<std::string, StoredValue>* variables, const std::map<std::string, CustomFunction>* functions, const HasScriptFunctionCallback& has_script_function, const InvokeScriptFunctionDecimalCallback& invoke_script_function, matrix::Value* value);

Rational parse_exact_expression(const std::string& expression, const std::map<std::string, StoredValue>* variables, const std::map<std::string, CustomFunction>* functions, HasScriptFunctionCallback has_script_function = {});

bool try_base_conversion_expression(const std::string& expression, const std::map<std::string, StoredValue>* variables, const std::map<std::string, CustomFunction>* functions, const HexFormatOptions& hex_options, std::string* output);

std::string matrix_literal_expression(const matrix::Matrix& value);

bool try_symbolic_constant_expression(const std::string& expression, const std::map<std::string, StoredValue>* variables, const std::map<std::string, CustomFunction>* functions, std::string* output);

// 变量作用域
std::map<std::string, StoredValue> visible_variables(const Calculator::Impl* impl);
bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name);
void assign_visible_variable(Calculator::Impl* impl, const std::string& name, const StoredValue& value);

// 表达式求值（返回 StoredValue）
StoredValue evaluate_expression_value(Calculator* calculator, Calculator::Impl* impl, const std::string& expression, bool exact_mode);

// 脚本函数调用
double invoke_script_function_decimal(Calculator* calculator, Calculator::Impl* impl, const std::string& name, const std::vector<double>& arguments);

// 脚本执行
ScriptSignal execute_script_statement(Calculator* calculator, Calculator::Impl* impl, const script::Statement& statement, bool exact_mode, std::string* last_output, bool create_scope);
std::string render_script_block(const script::BlockStatement& block, int indent);

#endif

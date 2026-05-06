// ============================================================================
// Calculator 实现类
// ============================================================================
//
// Calculator 类的内部实现，使用 Pimpl 模式隐藏细节。
// 存储计算器的所有状态：变量、函数、模块、显示选项等。
//
// 注意：此文件仅供 Calculator 内部使用，不应被外部模块直接引用。
// ============================================================================

#ifndef CORE_CALCULATOR_IMPL_H
#define CORE_CALCULATOR_IMPL_H

#include "core/calculator.h"
#include "core/scope.h"
#include "execution/command_registry.h"
#include "types/function.h"

#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// 前向声明
class CalculatorModule;
struct CoreServices;
struct CommandKey;

// ============================================================================
// 显示精度常量
// ============================================================================

/** @brief 判断数值是否为零的显示阈值 */
constexpr double kDisplayZeroEps = 1e-290;  // 近似 double 最小值

/** @brief 判断数值是否为整数的显示阈值 */
constexpr double kDisplayIntegerEps = 1e-9;

/** @brief 默认十进制显示有效位数 */
constexpr int kDefaultDisplayPrecision = 12;

/** @brief 十进制显示有效位数范围 */
constexpr int kMinDisplayPrecision = 1;
constexpr int kMaxDisplayPrecision = 17;

// ============================================================================
// 模块绑定
// ============================================================================

struct CommandBinding {
    std::shared_ptr<CalculatorModule> module;
    std::string dispatch_name;
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
 * - 模块注册
 * - 显示选项
 */
struct Calculator::Impl {
    // 变量存储
    std::map<std::string, StoredValue> variables;          ///< 全局变量
    FlatScopeStack flat_scopes;                            ///< 局部作用域栈

    // 函数存储
    std::map<std::string, CustomFunction> functions;       ///< 简单函数
    std::map<std::string, ScriptFunction> script_functions; ///< 脚本函数

    // 模块管理
    std::vector<std::shared_ptr<CalculatorModule>> registered_modules;
    std::vector<std::shared_ptr<CalculatorModule>> implicit_evaluation_modules;

    // 函数汇总（由模块注册）
    std::map<std::string, std::function<double(const std::vector<double>&)>> scalar_functions;
    std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> matrix_functions;
    std::map<std::string, matrix::ValueFunction> value_functions;
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> native_functions;

    // 模块元数据
    std::vector<std::string> module_commands;
    std::vector<std::string> module_functions;
    std::map<std::string, std::vector<std::shared_ptr<CalculatorModule>>> help_topic_to_modules;
    std::map<CommandKey, CommandBinding> command_to_module;

    // 命令注册表
    CommandRegistry command_registry;

    // 核心服务缓存
    std::unique_ptr<CoreServices> core_services;

    // 显示选项
    bool symbolic_constants_mode = false;  ///< 符号常量模式（pi, e 保留符号形式）
    bool hex_prefix_mode = false;          ///< 十六进制输出前缀
    bool hex_uppercase_mode = true;        ///< 十六进制大写字母
    int display_precision = kDefaultDisplayPrecision; ///< 十进制显示有效位数

    // 脚本执行状态
    int script_call_depth = 0;             ///< 脚本递归深度计数器
    std::vector<std::filesystem::path> script_file_stack; ///< 当前脚本文件栈
    std::set<std::filesystem::path> importing_script_files; ///< 正在导入的脚本，防循环
};

// ============================================================================
// 辅助函数声明
// ============================================================================

// 显示精度
void set_process_display_precision(int precision);
int process_display_precision();

// 字符串处理
bool is_valid_variable_name(std::string_view name);
bool is_identifier_text(std::string_view text);
bool is_string_literal(std::string_view text);
std::string parse_string_literal_value(std::string_view text);
bool is_reserved_user_function_name(const Calculator::Impl* impl, std::string_view name);

// 状态持久化
std::string encode_state_field(const std::string& text);
std::string decode_state_field(const std::string& text);

// 值格式化
void apply_calculator_display_precision(const Calculator::Impl* impl);

// 线性方程组
std::vector<double> solve_dense_linear_system(
    std::vector<std::vector<double>> matrix,
    std::vector<double> rhs,
    const std::string& context);

// 模块注册
void register_standard_modules(Calculator* calculator);

#endif // CORE_CALCULATOR_IMPL_H
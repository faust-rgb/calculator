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

#include "core/calculator.h"

#include "matrix/matrix.h"
#include "math/mymath.h"
#include "script/script_ast.h"
#include "precise/rational.h"
#include "precise/precise_decimal.h"
#include "types/stored_value.h"

#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>

#include "core/calculator_exceptions.h"
#include "module/calculator_module.h"
#include "command/command_registry.h"
#include "core/string_utils.h"
#include "core/format_utils.h"
#include "core/expression_utils.h"
#include "command/variable_resolver.h"
#include "parser/unified_expression_parser.h"
#include "parser/symbolic_render_parser.h"
#include "script/script_signal.h"
#include "command/builtin_constants.h"
#include "command/inline_expander.h"
#include "script/script_runtime.h"

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

class CalculatorModule;

struct CommandBinding {
    std::shared_ptr<CalculatorModule> module;
    std::string dispatch_name;
};

// ============================================================================
// Calculator 实现类
// ============================================================================

// ============================================================================
// 平坦作用域栈（高性能变量存储）
// ============================================================================

/**
 * @struct VariableSlot
 * @brief 变量槽位，用于平坦数组存储
 */
struct VariableSlot {
    std::string name;
    StoredValue value;
    int scope_level;  ///< 0 = 全局, 1+ = 局部作用域层级
};

/**
 * @struct FlatScopeStack
 * @brief 平坦作用域栈，用连续数组替代 std::vector<std::map>
 *
 * 优势：
 * - 连续内存布局，缓存友好
 * - 无 per-variable 堆分配
 * - O(n) 线性搜索，对于少量局部变量（通常 <10）比 O(log n) 的红黑树更快
 */
struct FlatScopeStack {
    std::vector<VariableSlot> slots;        ///< 所有变量连续存储
    std::vector<std::size_t> scope_starts;  ///< 每个作用域的起始索引
    int current_scope_level = 0;            ///< 当前作用域层级

    /// 进入新作用域
    void push_scope() {
        scope_starts.push_back(slots.size());
        ++current_scope_level;
    }

    /// 退出当前作用域
    void pop_scope() {
        if (scope_starts.empty()) return;
        std::size_t start = scope_starts.back();
        scope_starts.pop_back();
        slots.resize(start);
        --current_scope_level;
    }

    /// 从当前作用域向外搜索变量
    VariableSlot* find(const std::string& name) {
        // 从后向前搜索，优先匹配内层作用域
        for (auto it = slots.rbegin(); it != slots.rend(); ++it) {
            if (it->name == name) {
                return &(*it);
            }
        }
        return nullptr;
    }

    /// 从当前作用域向外搜索变量（const 版本）
    const VariableSlot* find(const std::string& name) const {
        for (auto it = slots.rbegin(); it != slots.rend(); ++it) {
            if (it->name == name) {
                return &(*it);
            }
        }
        return nullptr;
    }

    /// 仅在当前作用域查找
    VariableSlot* find_in_current_scope(const std::string& name) {
        if (scope_starts.empty()) return nullptr;
        std::size_t start = scope_starts.back();
        for (std::size_t i = start; i < slots.size(); ++i) {
            if (slots[i].name == name) {
                return &slots[i];
            }
        }
        return nullptr;
    }

    /// 设置变量值（在当前作用域）
    void set(const std::string& name, const StoredValue& value) {
        // 先在当前作用域查找是否已存在
        if (VariableSlot* existing = find_in_current_scope(name)) {
            existing->value = value;
            return;
        }
        // 不存在则新建
        slots.push_back({name, value, current_scope_level});
    }

    /// 设置或更新变量值（搜索所有作用域，找到则更新，否则在当前作用域创建）
    void set_or_create(const std::string& name, const StoredValue& value) {
        if (VariableSlot* existing = find(name)) {
            existing->value = value;
            return;
        }
        set(name, value);
    }

    /// 清空所有作用域
    void clear() {
        slots.clear();
        scope_starts.clear();
        current_scope_level = 0;
    }

    /// 获取当前作用域的变量数量
    std::size_t current_scope_size() const {
        if (scope_starts.empty()) return slots.size();
        return slots.size() - scope_starts.back();
    }

    /// 获取总变量数量
    std::size_t total_size() const {
        return slots.size();
    }

    /// 获取作用域深度
    int scope_depth() const {
        return static_cast<int>(scope_starts.size());
    }
};

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

    FlatScopeStack flat_scopes;  ///< 统一的平坦作用域栈（移除旧的 local_scopes）

    std::vector<std::shared_ptr<CalculatorModule>> registered_modules; ///< 已注册的数学模块
    std::vector<std::shared_ptr<CalculatorModule>> implicit_evaluation_modules; ///< 优化后的隐式求值模块列表

    std::map<std::string, std::function<double(const std::vector<double>&)>> scalar_functions; ///< 汇总的标量函数
    std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> matrix_functions; ///< 汇总的矩阵函数
    std::map<std::string, matrix::ValueFunction> value_functions; ///< 汇总的值多态函数
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> native_functions; ///< 原生函数

    // 汇总的模块元数据，用于补全和帮助
    std::vector<std::string> module_commands;
    std::vector<std::string> module_functions;
    std::map<std::string, std::vector<std::shared_ptr<CalculatorModule>>> help_topic_to_modules;
    std::map<CommandKey, CommandBinding> command_to_module;

    // 命令注册表（新架构）
    CommandRegistry command_registry;

    // 缓存的核心服务
    std::unique_ptr<CoreServices> core_services;

    bool symbolic_constants_mode = false;  ///< 符号常量模式（pi, e 保留符号形式）
    bool hex_prefix_mode = false;          ///< 十六进制输出前缀
    bool hex_uppercase_mode = true;        ///< 十六进制大写字母
    int display_precision = kDefaultDisplayPrecision; ///< 十进制显示有效位数
    int script_call_depth = 0;             ///< 脚本递归深度计数器
    std::vector<std::filesystem::path> script_file_stack; ///< 当前脚本文件栈，用于相对 import
    std::set<std::filesystem::path> importing_script_files; ///< 正在导入/执行的脚本，用于防循环
};

// ============================================================================
// 辅助函数声明（前向声明或关键辅助）
// ============================================================================

// 显示格式化
void set_process_display_precision(int precision);    ///< 设置进程级显示有效位数
int process_display_precision();                      ///< 查询进程级显示有效位数

// 内置常量

// 有理数和精确小数的前向声明已由 types/ 包含

// 字符串处理
bool is_valid_variable_name(std::string_view name);  ///< 检查变量名合法性
bool is_identifier_text(std::string_view text);      ///< 检查是否为标识符
bool is_string_literal(std::string_view text);       ///< 检查是否为字符串字面量
std::string parse_string_literal_value(std::string_view text); ///< 解析字符串字面量
bool is_reserved_user_function_name(const Calculator::Impl* impl, std::string_view name);

// 状态持久化
std::string encode_state_field(const std::string& text); ///< 编码状态字段
std::string decode_state_field(const std::string& text); ///< 解码状态字段

// 值格式化
void apply_calculator_display_precision(const Calculator::Impl* impl);

// 线性方程组
std::vector<double> solve_dense_linear_system(std::vector<std::vector<double>> matrix, std::vector<double> rhs, const std::string& context);

// 模块注册
void register_standard_modules(Calculator* calculator);

#endif

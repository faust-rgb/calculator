#ifndef CALCULATOR_MODULE_H
#define CALCULATOR_MODULE_H

#include "command_types.h"
#include "calculator_internal_types.h"
#include "symbolic_expression.h"
#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <map>
#include <array>

// Forward declarations
class FunctionAnalysis;

/**
 * @struct ServiceInterfaces
 * @brief 核心提供的分项服务接口
 */
struct IEvaluationService {
    std::function<double(const std::string&)> parse_decimal;
    std::function<StoredValue(const std::string&, bool)> evaluate_value;
    std::function<double(double)> normalize_result;

    // 作用域求值器构造
    std::function<std::function<double(const std::vector<std::pair<std::string, double>>&)>(const std::string&)> build_decimal_evaluator;
    std::function<std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_scalar_evaluator;
    std::function<std::function<matrix::Matrix(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_matrix_evaluator;
};

struct ISymbolicService {
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;
    std::function<std::string(const std::string&)> expand_inline;
    std::function<std::string(const std::string&)> simplify_symbolic;
    std::function<double(const SymbolicExpression&, const std::string&, double)> evaluate_symbolic_at;
    std::function<std::vector<SymbolicExpression>(const std::string&)> parse_symbolic_expr_list;
    std::function<FunctionAnalysis(const std::string&)> build_analysis;
};

struct IEnvironmentService {
    std::function<bool(const std::string&)> has_variable;
    std::function<bool(const std::string&)> has_function;
    std::function<std::string()> list_variables;
    std::function<std::string()> list_functions;
    std::function<std::string(const std::string&)> clear_variable;
    std::function<std::string(const std::string&)> clear_function;
    std::function<std::string()> clear_all_variables;
    std::function<std::string()> clear_all_functions;

    // 状态管理
    std::function<std::string(const std::string&)> save_state;
    std::function<std::string(const std::string&)> load_state;
    std::function<std::string(const std::string&)> export_variable;
    std::function<std::string(const std::string&, bool)> execute_script;

    // 模式与配置管理
    std::function<std::string(bool)> set_exact_mode;
    std::function<std::string(bool)> set_symbolic_mode;
    std::function<std::string(int)> set_precision;
    std::function<std::string(bool)> set_hex_prefix;
    std::function<std::string(bool)> set_hex_uppercase;
};

/**
 * @struct CoreServices
 * @brief 汇总后的核心服务，作为各模块访问核心功能的唯一入口
 */
struct CoreServices {
    IEvaluationService evaluation;
    ISymbolicService symbolic;
    IEnvironmentService env;

    // 参数解析与辅助（保留作为通用工具）
    std::function<std::vector<std::string>(const std::vector<std::string>&, std::size_t, const std::vector<std::string>&)> parse_symbolic_vars;
    std::function<bool(const std::string&)> is_matrix_argument;
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;
    std::function<std::string(const std::vector<std::string>&, bool)> render_plot;
    std::function<bool(double, double)> is_integer_double;
    std::function<long long(double)> round_to_long_long;
};

struct CommandSpec {
    CommandKey key;
    std::string dispatch_name;
};

/**
 * @struct CalculatorSettings
 * @brief 汇总计算器的全局配置状态
 */
struct CalculatorSettings {
    int display_precision;
    bool exact_mode;
    bool symbolic_constants_mode;
    bool hex_prefix_mode;
    bool hex_uppercase_mode;
};

/**
 * @class CalculatorModule
 * @brief 所有数学模块的基类，简化命令注册接口
 */
class CalculatorModule {
public:
    virtual ~CalculatorModule() = default;

    // 模块基本信息
    virtual std::string name() const = 0;

    /**
     * @brief 模块初始化钩子
     * @param services 核心提供的服务
     */
    virtual void initialize(const CoreServices& /*services*/) {}

    /**
     * @brief 暴露模块特定的服务接口（用于跨模块协作）
     * @param service_name 服务名
     * @return 服务指针，如果不支持则返回 nullptr
     */
    virtual void* query_service(const std::string& service_name) {
        (void)service_name;
        return nullptr;
    }

    /**
     * @brief 设置变化钩子（观察者模式）
     * @param settings 最新的全局配置
     */
    virtual void on_settings_changed(const CalculatorSettings& /*settings*/) {}

    // ============================================================================
    // 命令注册接口 - 简化版
    // ============================================================================

    /**
     * @brief 返回该模块支持的命令列表
     *
     * 命令格式：
     * - 元命令以冒号开头，如 ":vars", ":help"
     * - 函数命令不带冒号，如 "plot", "print"
     *
     * 基类的 get_command_specs() 会自动根据此列表生成 CommandSpec。
     */
    virtual std::vector<std::string> get_commands() const { return {}; }

    /**
     * @brief 自动生成命令规格（基类默认实现）
     *
     * 根据 get_commands() 自动派生：
     * - 冠号开头的命令 → MetaCommandKey
     * - 其他命令 → CallCommandKey
     */
    virtual std::vector<CommandSpec> get_command_specs() const {
        std::vector<CommandSpec> specs;
        for (const std::string& cmd : get_commands()) {
            bool is_meta = !cmd.empty() && cmd.front() == ':';
            std::string key_name = is_meta ? cmd.substr(1) : cmd;
            CommandKey key = is_meta ? meta_command_key(key_name)
                                     : call_command_key(key_name);
            specs.push_back({key, cmd});
        }
        return specs;
    }

    /**
     * @brief 执行具体命令（结构化版本，接收预拆分参数）
     *
     * 子模块只需实现此方法，基类提供默认的参数拼接逻辑。
     */
    virtual std::string execute_args(const std::string& command,
                                    const std::vector<std::string>& args,
                                    const CoreServices& services) {
        // 默认实现：将 args 拼接为 inside 字符串，调用旧版 execute
        std::string inside;
        for (std::size_t i = 0; i < args.size(); ++i) {
            if (i != 0) inside += ", ";
            inside += args[i];
        }
        return execute(command, inside, services);
    }

    /**
     * @brief 执行具体命令（string_view 版本，避免不必要的拷贝）
     *
     * 默认实现转换为 string 后调用 execute_args。
     * 子模块可重写此方法直接使用 string_view 以提升性能。
     */
    virtual std::string execute_args_view(std::string_view command,
                                          const std::vector<std::string_view>& args,
                                          const CoreServices& services) {
        std::string cmd(command);
        std::vector<std::string> string_args;
        string_args.reserve(args.size());
        for (auto arg : args) {
            string_args.emplace_back(arg);
        }
        return execute_args(cmd, string_args, services);
    }

    /**
     * @brief 执行具体命令（旧版接口，已废弃但保留兼容）
     */
    virtual std::string execute(const std::string& command,
                               const std::string& inside,
                               const CoreServices& services) {
        (void)command; (void)inside; (void)services;
        return "";
    }

    // ============================================================================
    // 隐式求值接口
    // ============================================================================

    /**
     * @brief 隐式求值触发特征（可选）
     * @return 返回该模块感兴趣的字符串特征（如 "0123456789."），
     *         核心据此过滤掉明显不相关的输入。
     */
    virtual std::string get_implicit_trigger_chars() const { return ""; }

    /**
     * @brief 是否支持隐式求值（用于核心层优化过滤）
     */
    virtual bool wants_implicit_evaluation() const { return false; }

    /**
     * @brief 获取缓存的触发字符查找表
     *
     * 首次调用时构建并缓存，后续直接返回缓存的表。
     * 用于优化热路径性能，避免重复构建。
     */
    const std::array<bool, 256>* get_cached_trigger_table() const {
        if (!trigger_table_cached_) {
            std::string triggers = get_implicit_trigger_chars();
            if (!triggers.empty()) {
                trigger_table_.fill(false);
                for (char c : triggers) {
                    trigger_table_[static_cast<unsigned char>(c)] = true;
                }
            }
            trigger_table_cached_ = true;
        }
        return trigger_table_cached_ && !get_implicit_trigger_chars().empty() ? &trigger_table_ : nullptr;
    }

    // 隐式求值钩子
    virtual bool try_evaluate_implicit(const std::string&,
                                      StoredValue*,
                                      const std::map<std::string, StoredValue>&) const { return false; }

    // ============================================================================
    // 函数注册接口
    // ============================================================================

    /**
     * @brief 获取该模块提供的各种函数映射
     */
    virtual std::map<std::string, std::function<double(const std::vector<double>&)>> get_scalar_functions() const { return {}; }
    virtual std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> get_matrix_functions() const { return {}; }

    using ValueFunction = matrix::ValueFunction;
    virtual std::map<std::string, ValueFunction> get_value_functions() const { return {}; }

    /**
     * @brief 自描述接口：返回该模块支持的内置函数列表，用于自动补全
     */
    virtual std::vector<std::string> get_functions() const { return {}; }

    // ============================================================================
    // 帮助接口
    // ============================================================================

    /**
     * @brief 自描述接口：返回特定主题的帮助文本
     * @param topic 主题名（如 "commands", "matrix" 等）
     * @return 帮助内容，如果不支持该主题则返回空字符串
     */
    virtual std::string get_help_snippet(const std::string& topic) const {
        (void)topic;
        return "";
    }

protected:
    // 缓存的触发字符查找表
    mutable std::array<bool, 256> trigger_table_{};
    mutable bool trigger_table_cached_ = false;
};

#endif
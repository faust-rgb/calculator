#ifndef CALCULATOR_MODULE_H
#define CALCULATOR_MODULE_H

#include "calculator_internal_types.h"
#include "symbolic_expression.h"
#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <map>

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
    std::function<bool(double, double)> is_integer_double;
    std::function<long long(double)> round_to_long_long;
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
 * @brief 所有数学模块的基类，增加了事件钩子和隐式求值元数据
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

    // 命令处理逻辑
    virtual bool can_handle(const std::string&) const { return false; }
    
    /**
     * @brief 执行具体命令（原始版本，接收内部字符串）
     */
    virtual std::string execute(const std::string& command, 
                               const std::string& inside, 
                               const CoreServices& services) { 
        (void)command; (void)inside; (void)services;
        return ""; 
    }

    /**
     * @brief 执行具体命令（结构化版本，接收预拆分参数）
     */
    virtual std::string execute_args(const std::string& command,
                                    const std::vector<std::string>& args,
                                    const CoreServices& services) {
        // 默认实现回退到原始 execute
        std::string inside;
        for (std::size_t i = 0; i < args.size(); ++i) {
            if (i != 0) inside += ", ";
            inside += args[i];
        }
        return execute(command, inside, services);
    }

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
                               
    // 隐式求值钩子
    virtual bool try_evaluate_implicit(const std::string&, 
                                      StoredValue*, 
                                      const std::map<std::string, StoredValue>&) const { return false; }

    /**
     * @brief 获取该模块提供的各种函数映射
     */
    virtual std::map<std::string, std::function<double(const std::vector<double>&)>> get_scalar_functions() const { return {}; }
    virtual std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> get_matrix_functions() const { return {}; }
    
    using ValueFunction = matrix::ValueFunction;
    virtual std::map<std::string, ValueFunction> get_value_functions() const { return {}; }

    /**
     * @brief 自描述接口：返回该模块支持的命令列表，用于自动补全
     */
    virtual std::vector<std::string> get_commands() const { return {}; }

    /**
     * @brief 自描述接口：返回该模块支持的内置函数列表，用于自动补全
     */
    virtual std::vector<std::string> get_functions() const { return {}; }

    /**
     * @brief 自描述接口：返回特定主题的帮助文本
     * @param topic 主题名（如 "commands", "matrix" 等）
     * @return 帮助内容，如果不支持该主题则返回空字符串
     */
    virtual std::string get_help_snippet(const std::string& topic) const { 
        (void)topic;
        return ""; 
    }
};

#endif

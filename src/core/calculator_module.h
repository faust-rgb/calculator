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
 * @struct CoreServices
 * @brief core 提供的基础设施服务，模块通过这些服务与核心交互
 */
struct CoreServices {
    // 基础求值服务
    std::function<double(const std::string&)> parse_decimal;
    std::function<StoredValue(const std::string&, bool)> evaluate_value;
    std::function<double(double)> normalize_result;
    
    // 符号计算服务
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;
    std::function<std::string(const std::string&)> expand_inline;
    std::function<std::string(const std::string&)> simplify_symbolic;
    std::function<double(const SymbolicExpression&, const std::string&, double)> evaluate_symbolic_at;
    
    // 作用域求值器构造服务
    std::function<std::function<double(const std::vector<std::pair<std::string, double>>&)>(const std::string&)> build_decimal_evaluator;
    std::function<std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_scalar_evaluator;
    std::function<std::function<matrix::Matrix(const std::vector<std::pair<std::string, StoredValue>>&)>(const std::string&)> build_matrix_evaluator;

    // 参数解析辅助
    std::function<std::vector<std::string>(const std::vector<std::string>&, std::size_t, const std::vector<std::string>&)> parse_symbolic_vars;
    std::function<std::vector<SymbolicExpression>(const std::string&)> parse_symbolic_expr_list;
    std::function<FunctionAnalysis(const std::string&)> build_analysis;
    
    // 矩阵与辅助
    std::function<bool(const std::string&)> is_matrix_argument;
    std::function<matrix::Matrix(const std::string&, const std::string&)> parse_matrix_argument;
    std::function<bool(double, double)> is_integer_double;
    std::function<long long(double)> round_to_long_long;

    // 环境查询与管理
    std::function<bool(const std::string&)> has_variable;
    std::function<bool(const std::string&)> has_function;
    std::function<std::string()> list_variables;
    std::function<std::string()> list_functions;
    std::function<std::string()> clear_all_variables;
    std::function<std::string(const std::string&)> clear_variable;
    std::function<std::string()> clear_all_functions;
    std::function<std::string(const std::string&)> clear_function;

    // 系统服务
    std::function<std::string(const std::string&)> save_state;
    std::function<std::string(const std::string&)> load_state;
    std::function<std::string(const std::string&)> export_variable;
    std::function<std::string()> get_history;
    std::function<std::string(const std::string&, bool)> execute_script;
};

/**
 * @class CalculatorModule
 * @brief 所有数学模块的基类
 */
class CalculatorModule {
public:
    virtual ~CalculatorModule() = default;
    
    // 模块名称
    virtual std::string name() const = 0;
    
    // 该模块支持的命令列表匹配逻辑
    virtual bool can_handle(const std::string&) const { return false; }
    
    // 执行具体命令
    virtual std::string execute(const std::string&, 
                               const std::string&, 
                               const CoreServices&) { return ""; }
                               
    // 隐式求值钩子（如高精度自动识别）
    virtual bool try_evaluate_implicit(const std::string&, 
                                      StoredValue*, 
                                      const std::map<std::string, StoredValue>&) const { return false; }

    /**
     * @brief 获取该模块提供的标量函数映射
     * @return 函数名 -> (参数列表 -> 结果)
     */
    virtual std::map<std::string, std::function<double(const std::vector<double>&)>> get_scalar_functions() const { return {}; }

    /**
     * @brief 获取该模块提供的矩阵函数映射
     * @return 函数名 -> (参数列表 -> 结果)
     */
    virtual std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>> get_matrix_functions() const { return {}; }

    /**
     * @brief 值多态函数类型：接受参数字符串列表和求值上下文，返回 Value
     *
     * 这种函数可以处理标量、复数、矩阵等多种输入类型，
     * 并根据输入类型返回相应的结果。
     *
     * @param arguments 参数字符串列表（未求值）
     * @param scalar_evaluator 标量求值器
     * @param matrix_lookup 矩阵查找函数
     * @param complex_lookup 复数查找函数
     * @param matrix_functions 矩阵函数表（用于递归调用）
     * @return 求值结果
     */
    using ValueFunction = matrix::ValueFunction;

    /**
     * @brief 获取该模块提供的值多态函数映射
     * @return 函数名 -> 值多态函数
     */
    virtual std::map<std::string, ValueFunction> get_value_functions() const { return {}; }
};

#endif

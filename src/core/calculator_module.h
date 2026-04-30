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

    // 环境查询
    std::function<bool(const std::string&)> has_variable;
    std::function<bool(const std::string&)> has_function;
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
};

#endif

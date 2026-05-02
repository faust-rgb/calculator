#ifndef CORE_SERVICE_INTERFACES_H
#define CORE_SERVICE_INTERFACES_H

#include "core/calculator_internal_types.h"
#include "symbolic/symbolic_expression.h"
#include <string>
#include <vector>
#include <functional>
#include <map>

// Forward declarations
class FunctionAnalysis;

/**
 * @struct IEvaluationService
 * @brief 核心提供的求值服务接口
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


#endif

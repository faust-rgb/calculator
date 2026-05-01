// ============================================================================
// 符号命令头文件
// ============================================================================
//
// 定义符号计算命令的接口和数据结构。主要包含：
//
// 1. SymbolicResolverContext - 符号表达式解析上下文
//    用于解析嵌套表达式和自定义函数
//
// 2. SymbolicCommandContext - 符号命令执行上下文
//    提供命令执行所需的回调函数接口
//
// 3. 命令处理函数
//    - is_symbolic_command: 判断是否为符号命令
//    - resolve_symbolic_expression: 解析符号表达式
//    - handle_symbolic_command: 分发处理具体命令
// ============================================================================

#ifndef CALCULATOR_SYMBOLIC_COMMANDS_H
#define CALCULATOR_SYMBOLIC_COMMANDS_H

#include "function_analysis.h"
#include "symbolic_expression.h"
#include "calculator_polynomial.h"

#include <functional>
#include <string>
#include <vector>

#include "../core/calculator_module.h"

namespace symbolic_commands {

/**
 * @class SymbolicModule
 * @brief 提供核心符号运算（微分、化简、极值等）功能的模块
 */
class SymbolicModule : public CalculatorModule {
public:
    std::string name() const override { return "Symbolic"; }

    std::vector<std::string> get_commands() const override {
        return {"simplify", "diff", "gradient", "jacobian", "hessian", "divergence", "div",
                "curl", "curl_2d", "laplacian", "implicit_diff", "param_deriv", "directional",
                "line_integral", "line_integral_vector", "surface_integral", "greens_theorem",
                "stokes_theorem", "divergence_theorem", "integrate_region", "integral", "dsolve"};
    }

    bool can_handle(const std::string& command) const override;

    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;

    std::string get_help_snippet(const std::string& topic) const override;
};

struct SymbolicResolverContext {
    std::function<SymbolicExpression(const std::string&, std::string*)>
        resolve_custom_function;
    std::function<bool(const std::string&)> has_custom_function;
    std::function<std::string(const std::string&)> expand_inline;
};

struct SymbolicCommandContext {
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)>
        resolve_symbolic;
    std::function<std::vector<std::string>(
        const std::vector<std::string>&,
        std::size_t,
        const std::vector<std::string>&)> parse_symbolic_variable_arguments;
    std::function<std::vector<SymbolicExpression>(const std::string&)>
        parse_symbolic_expression_list;
    std::function<FunctionAnalysis(const std::string&)> build_analysis;
    std::function<std::function<double(const std::vector<std::pair<std::string, double>>&)>(const std::string&)>
        build_scoped_evaluator;
    std::function<double(const std::string&)> parse_decimal;
    std::function<double(double)> normalize_result;
};

bool is_symbolic_command(const std::string& command);

void resolve_symbolic_expression(const SymbolicResolverContext& ctx,
                                 const std::string& argument,
                                 bool require_single_variable,
                                 std::string* variable_name,
                                 SymbolicExpression* expression);

std::vector<std::string> parse_symbolic_variable_arguments(
    const std::vector<std::string>& arguments,
    std::size_t start_index,
    const std::vector<std::string>& fallback_variables);

std::vector<SymbolicExpression> parse_symbolic_expression_list(
    const std::string& argument,
    const std::function<std::string(const std::string&)>& expand_inline);

bool handle_symbolic_command(const SymbolicCommandContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             std::string* output);

}  // namespace symbolic_commands

#endif  // CALCULATOR_SYMBOLIC_COMMANDS_H

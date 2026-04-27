// ============================================================================
// 符号命令
// ============================================================================

#ifndef CALCULATOR_SYMBOLIC_COMMANDS_H
#define CALCULATOR_SYMBOLIC_COMMANDS_H

#include "function_analysis.h"
#include "symbolic_expression.h"
#include "calculator_polynomial.h"

#include <functional>
#include <string>
#include <vector>

namespace symbolic_commands {

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

// ============================================================================
// 计算器命令辅助函数实现
// ============================================================================

#include "calculator_command_helpers.h"

#include <stdexcept>

namespace command_helpers {

SymbolicExpression build_symbolic_expression_from_function(
    const std::map<std::string, CustomFunction>& functions,
    const std::string& name,
    std::string* variable_name) {
    const auto it = functions.find(name);
    if (it == functions.end()) {
        throw std::runtime_error("unknown custom function: " + name);
    }
    *variable_name = it->second.parameter_name;
    return SymbolicExpression::parse(it->second.expression);
}

}  // namespace command_helpers

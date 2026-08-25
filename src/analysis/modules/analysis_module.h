/**
 * @file analysis_module.h
 * @brief 函数分析模块
 *
 * 本文件定义了函数分析模块：
 * - 极限计算：符号极限和数值极限
 * - 极值查找：局部极值和全局极值
 * - 函数性质分析：单调性、凸性等
 *
 * 支持符号分析和数值分析两种模式。
 */

#ifndef ANALYSIS_MODULE_H
#define ANALYSIS_MODULE_H

#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"
#include "types/scalar_type.h"

class ServiceLocator;

namespace analysis_cmds {

/**
 * @class AnalysisModule
 * @brief 提供高级函数 analysis 功能（极限、极值等）的模块
 */
class AnalysisModule : public CommandModuleBase {
public:
    std::string name() const override { return "Analysis"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args_view(std::string_view command,
                                  const std::vector<std::string_view>& args,
                                  ::ServiceLocator& locator) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct AnalysisContext {
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;
    std::function<std::vector<std::string>(const std::vector<std::string>&, std::size_t, const std::vector<std::string>&)> parse_symbolic_variable_arguments;
    std::function<Scalar(const std::string&)> parse_decimal;
    std::function<Scalar(Scalar)> normalize_result;
    std::function<FunctionAnalysis(const std::string&)> build_analysis;
};

bool is_analysis_command(const std::string& command);

bool handle_analysis_command(const AnalysisContext& ctx,
                             const std::string& command,
                             const std::vector<std::string>& arguments,
                             std::string* output);

}  // namespace analysis_cmds

#endif  // ANALYSIS_MODULE_H

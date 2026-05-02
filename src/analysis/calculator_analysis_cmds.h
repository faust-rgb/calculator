#ifndef CALCULATOR_ANALYSIS_CMDS_H
#define CALCULATOR_ANALYSIS_CMDS_H

#include <string>
#include <functional>
#include <vector>
#include "core/calculator_module.h"

namespace analysis_cmds {

/**
 * @class AnalysisModule
 * @brief 提供高级函数 analysis 功能（极限、极值等）的模块
 */
class AnalysisModule : public CalculatorModule {
public:
    std::string name() const override { return "Analysis"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct AnalysisContext {
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;
    std::function<std::vector<std::string>(const std::vector<std::string>&, std::size_t, const std::vector<std::string>&)> parse_symbolic_variable_arguments;
    std::function<double(const std::string&)> parse_decimal;
    std::function<double(double)> normalize_result;
    std::function<FunctionAnalysis(const std::string&)> build_analysis;
};

bool is_analysis_command(const std::string& command);

bool handle_analysis_command(const AnalysisContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             std::string* output);

}  // namespace analysis_cmds

#endif  // CALCULATOR_ANALYSIS_CMDS_H

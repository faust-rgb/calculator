#ifndef CALCULATOR_SERIES_H
#define CALCULATOR_SERIES_H

#include "calculator_internal_types.h"
#include "symbolic_expression.h"
#include <string>
#include <vector>
#include <map>
#include <functional>
#include "../core/calculator_module.h"

namespace series_ops {

/**
 * @class SeriesModule
 * @brief 提供级数展开和求和功能的模块
 */
class SeriesModule : public CalculatorModule {
public:
    std::string name() const override { return "Series"; }
    std::vector<std::string> get_commands() const override;
    bool can_handle(const std::string& command) const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct SeriesContext {
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;
    std::function<double(const std::string&)> parse_decimal;
    std::function<double(const SymbolicExpression&, const std::string&, double)> evaluate_at;
    std::function<std::string(const std::string&)> simplify_symbolic;
    std::function<std::string(const std::string&)> expand_inline;
};

bool is_series_command(const std::string& command);

bool handle_series_command(const SeriesContext& ctx,
                           const std::string& command,
                           const std::string& inside,
                           std::string* output);

namespace internal {
struct PoleException : public std::runtime_error {
    int shift;
    double leading_coefficient;
    PoleException(int s, double coeff) : std::runtime_error("Pole encountered"), shift(s), leading_coefficient(coeff) {}
};

bool evaluate_psa(const SymbolicExpression& expr, const std::string& var_name, double center, int degree, std::vector<double>& result, const SeriesContext& ctx);
}

}  // namespace series_ops

#endif  // CALCULATOR_SERIES_H

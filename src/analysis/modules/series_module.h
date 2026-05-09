#ifndef SERIES_MODULE_H
#define SERIES_MODULE_H

// ============================================================================
// 级数展开模块（重构版）
// ============================================================================
//
// 本模块提供级数展开和求和功能，核心算法已拆分到子模块：
// - series/psa_engine.h: 幂级数代数运算
// - series/taylor_series.h: Taylor 级数展开
// - series/pade_approximation.h: Pade 有理逼近
// - series/puiseux_series.h: Puiseux 级数
// - series/series_summation.h: 级数求和
//
// ============================================================================

#include "core/calculator_internal_types.h"
#include "core/scalar_type.h"
#include "symbolic/core/symbolic_expression.h"
#include "module/calculator_module.h"
#include "analysis/series/psa_engine.h"
#include <string>
#include <vector>
#include <map>
#include <functional>

namespace series_ops {

using Scalar = mymath::Scalar;

/**
 * @class SeriesModule
 * @brief 提供级数展开和求和功能的模块
 */
class SeriesModule : public CalculatorModule {
public:
    std::string name() const override { return "Series"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct SeriesContext {
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;
    std::function<Scalar(const std::string&)> parse_decimal;
    std::function<Scalar(const SymbolicExpression&, const std::string&, Scalar)> evaluate_at;
    std::function<std::string(const std::string&)> simplify_symbolic;
    std::function<std::string(const std::string&)> expand_inline;
};

bool is_series_command(const std::string& command);

bool handle_series_command(const SeriesContext& ctx,
                           const std::string& command,
                           const std::string& inside,
                           std::string* output);

// 内部命名空间中的 PoleException 和 evaluate_psa 已移至 series/psa_engine.h

}  // namespace series_ops

#endif  // SERIES_MODULE_H

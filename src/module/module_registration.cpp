// ============================================================================
// module_registration.cpp - 标准模块注册实现 (显式装配与确定性注册)
// ============================================================================

#include "module/calculator_module.h"
#include "core/api/calculator.h"

// 显式包含所有内建模块，消除静态库链接时未引用符号被剔除的风险
#include "analysis/modules/analysis_module.h"
#include "analysis/modules/integration_module.h"
#include "analysis/modules/ode_module.h"
#include "analysis/modules/optimization_module.h"
#include "analysis/modules/rootfinding_module.h"
#include "analysis/modules/series_module.h"
#include "dsp/dsp_module.h"
#include "io/io_module.h"
#include "math/modules/descriptive_stats_module.h"
#include "math/modules/integer_math_module.h"
#include "math/modules/precise_module.h"
#include "math/modules/probability_module.h"
#include "math/modules/standard_math_module.h"
#include "matrix/matrix_module.h"
#include "module/system_module.h"
#include "plot/plot_module.h"
#include "polynomial/calculator_polynomial.h"
#include "statistics/statistics_module.h"
#include "symbolic/modules/symbolic_calculus_module.h"
#include "symbolic/modules/symbolic_module.h"
#include "symbolic/modules/transform_module.h"
#include "time/time_module.h"

#include <set>

void register_standard_modules(Calculator* calculator) {
    if (!calculator) return;

    // 首先注册已通过 ModuleRegistry 登记的模块（支持第三方扩展）
    std::set<std::string> registered_names;
    for (const auto& factory : ModuleRegistry::instance().factories()) {
        auto mod = factory();
        if (mod) {
            registered_names.insert(mod->name());
            calculator->register_module(mod);
        }
    }

    // 显式装配内建模块（确定性兜底，防止静态库链接裁剪）
    auto register_if_missing = [&](const auto& create_fn, const std::string& name) {
        if (registered_names.find(name) == registered_names.end()) {
            calculator->register_module(create_fn());
            registered_names.insert(name);
        }
    };

    register_if_missing([]() { return std::make_shared<StandardMathModule>(); }, "standard_math");
    register_if_missing([]() { return std::make_shared<IntegerMathModule>(); }, "integer_math");
    register_if_missing([]() { return std::make_shared<PreciseModule>(); }, "precise");
    register_if_missing([]() { return std::make_shared<ProbabilityModule>(); }, "probability");
    register_if_missing([]() { return std::make_shared<DescriptiveStatsModule>(); }, "descriptive_stats");
    register_if_missing([]() { return std::make_shared<MatrixModule>(); }, "matrix");
    register_if_missing([]() { return std::make_shared<polynomial_ops::PolynomialModule>(); }, "polynomial");
    register_if_missing([]() { return std::make_shared<StatisticsModule>(); }, "statistics");
    register_if_missing([]() { return std::make_shared<DspModule>(); }, "dsp");
    register_if_missing([]() { return std::make_shared<symbolic_commands::SymbolicModule>(); }, "symbolic");
    register_if_missing([]() { return std::make_shared<symbolic_commands::SymbolicCalculusModule>(); }, "symbolic_calculus");
    register_if_missing([]() { return std::make_shared<transforms::TransformModule>(); }, "transform");
    register_if_missing([]() { return std::make_shared<analysis_cmds::AnalysisModule>(); }, "analysis");
    register_if_missing([]() { return std::make_shared<integration_ops::IntegrationModule>(); }, "integration");
    register_if_missing([]() { return std::make_shared<ode_ops::ODEModule>(); }, "ode");
    register_if_missing([]() { return std::make_shared<optimization::OptimizationModule>(); }, "optimization");
    register_if_missing([]() { return std::make_shared<rootfinding::RootfindingModule>(); }, "rootfinding");
    register_if_missing([]() { return std::make_shared<series_ops::SeriesModule>(); }, "series");
    register_if_missing([]() { return std::make_shared<PlotModule>(); }, "plot");
    register_if_missing([]() { return std::make_shared<TimeModule>(); }, "time");
    register_if_missing([]() { return std::make_shared<IoModule>(); }, "io");
    register_if_missing([]() { return std::make_shared<SystemModule>(); }, "system");
}

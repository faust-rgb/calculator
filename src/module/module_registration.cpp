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
#include "analysis/modules/pde_module.h"
#include "analysis/modules/rootfinding_module.h"
#include "analysis/modules/series_module.h"
#include "dsp/dsp_module.h"
#include "io/io_module.h"
#include "math/runtime/modules/integer_math_module.h"
#include "math/runtime/modules/precise_module.h"
#include "math/runtime/modules/standard_math_module.h"
#include "matrix/matrix_module.h"
#include "module/system_module.h"
#include "plot/plot_module.h"
#include "polynomial/calculator_polynomial.h"
#include "statistics/statistics_module.h"
#include "symbolic/modules/symbolic_calculus_module.h"
#include "symbolic/modules/symbolic_module.h"
#include "symbolic/modules/transform_module.h"
#include "time/time_module.h"

void register_standard_modules(Calculator* calculator) {
    if (!calculator) return;

    // Explicit assembly keeps startup order deterministic and avoids static
    // initialization and static-library dead-stripping issues.
    auto register_builtin = [&](const auto& create_fn) {
        calculator->register_module(create_fn());
    };

    register_builtin([]() { return std::make_shared<StandardMathModule>(); });
    register_builtin([]() { return std::make_shared<IntegerMathModule>(); });
    register_builtin([]() { return std::make_shared<PreciseModule>(); });
    register_builtin([]() { return std::make_shared<MatrixModule>(); });
    register_builtin([]() { return std::make_shared<polynomial_ops::PolynomialModule>(); });
    register_builtin([]() { return std::make_shared<StatisticsModule>(); });
    register_builtin([]() { return std::make_shared<DspModule>(); });
    register_builtin([]() { return std::make_shared<symbolic_commands::SymbolicModule>(); });
    register_builtin([]() { return std::make_shared<symbolic_commands::SymbolicCalculusModule>(); });
    register_builtin([]() { return std::make_shared<transforms::TransformModule>(); });
    register_builtin([]() { return std::make_shared<analysis_cmds::AnalysisModule>(); });
    register_builtin([]() { return std::make_shared<analysis::pde::PdeModule>(); });
    register_builtin([]() { return std::make_shared<integration_ops::IntegrationModule>(); });
    register_builtin([]() { return std::make_shared<ode_ops::ODEModule>(); });
    register_builtin([]() { return std::make_shared<optimization::OptimizationModule>(); });
    register_builtin([]() { return std::make_shared<rootfinding::RootfindingModule>(); });
    register_builtin([]() { return std::make_shared<series_ops::SeriesModule>(); });
    register_builtin([]() { return std::make_shared<PlotModule>(); });
    register_builtin([]() { return std::make_shared<TimeModule>(); });
    register_builtin([]() { return std::make_shared<IoModule>(); });
    register_builtin([]() { return std::make_shared<SystemModule>(); });
}

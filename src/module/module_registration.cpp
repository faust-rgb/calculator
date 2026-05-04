// ============================================================================
// 模块注册实现
// ============================================================================
//
// 本文件负责将所有标准模块注册到计算器实例。
// 模块按功能分组：
// 1. 基础数学与系统模块 - 标准数学函数、整数运算、精确计算、统计
// 2. 矩阵与 DSP 模块 - 矩阵运算、信号处理
// 3. 分析模块 - 级数展开、数值积分、求根、优化、ODE 求解
// 4. 符号计算与多项式 - 符号运算、多项式操作
// 5. 绘图模块 - 函数可视化
// ============================================================================

#include "module_registration.h"
#include "calculator.h"
#include "calculator_module.h"

// 基础数学与系统模块
#include "math/standard_math_module.h"
#include "math/integer_math_module.h"
#include "matrix/matrix_module.h"
#include "precise/precise_module.h"
#include "statistics/statistics_module.h"
#include "module/system_module.h"
#include "io/io_module.h"
#include "time/time_module.h"

// 分析模块
#include "analysis/calculator_series.h"
#include "analysis/calculator_integration.h"
#include "analysis/calculator_rootfinding.h"
#include "analysis/calculator_optimization.h"
#include "analysis/calculator_analysis_cmds.h"
#include "analysis/calculator_ode.h"

// 符号计算与多项式
#include "symbolic/calculator_symbolic_commands.h"
#include "symbolic/calculator_transforms.h"
#include "polynomial/calculator_polynomial.h"

// 信号处理
#include "dsp/dsp_module.h"
#include "plot/plot_module.h"

/// 注册所有标准模块到计算器实例
void register_standard_modules(Calculator* calculator) {
    // 注册基础数学模块
    calculator->register_module(std::make_shared<StandardMathModule>());
    calculator->register_module(std::make_shared<IntegerMathModule>());
    calculator->register_module(std::make_shared<PreciseModule>());
    calculator->register_module(std::make_shared<StatisticsModule>());
    
    // 注册矩阵与 DSP 模块
    calculator->register_module(std::make_shared<MatrixModule>());
    calculator->register_module(std::make_shared<DspModule>());
    calculator->register_module(std::make_shared<PlotModule>());

    // 注册核心功能模块
    calculator->register_module(std::make_shared<SystemModule>());
    calculator->register_module(std::make_shared<IoModule>());
    calculator->register_module(std::make_shared<TimeModule>());
    calculator->register_module(std::make_shared<polynomial_ops::PolynomialModule>());
    calculator->register_module(std::make_shared<series_ops::SeriesModule>());
    calculator->register_module(std::make_shared<transforms::TransformModule>());
    calculator->register_module(std::make_shared<integration_ops::IntegrationModule>());
    calculator->register_module(std::make_shared<rootfinding::RootfindingModule>());
    calculator->register_module(std::make_shared<optimization::OptimizationModule>());
    calculator->register_module(std::make_shared<symbolic_commands::SymbolicModule>());
    calculator->register_module(std::make_shared<analysis_cmds::AnalysisModule>());
    calculator->register_module(std::make_shared<ode_ops::ODEModule>());
}

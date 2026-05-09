// ============================================================================
// module_registration.cpp - 标准模块注册实现
// ============================================================================
//
// 本文件实现 register_standard_modules() 函数，负责将所有标准模块
// 注册到计算器实例。模块按功能分组注册：
//
// 1. 基础数学与系统模块 - 标准数学函数、整数运算、精确计算、统计
// 2. 矩阵与 DSP 模块 - 矩阵运算、信号处理
// 3. 分析模块 - 级数展开、数值积分、求根、优化、ODE 求解
// 4. 符号计算与多项式 - 符号运算、多项式操作
// 5. 绘图模块 - 函数可视化
// 6. 辅助模块 - I/O、时间
// ============================================================================

#include "module_registration.h"
#include "calculator.h"
#include "calculator_module.h"

// ==================== 基础数学与系统模块 ====================
#include "math/standard_math_module.h"    // 标准数学函数（sin, cos, exp 等）
#include "math/integer_math_module.h"     // 整数运算（阶乘、组合数等）
#include "matrix/matrix_module.h"         // 矩阵运算
#include "precise/precise_module.h"       // 高精度计算
#include "statistics/statistics_module.h" // 统计函数
#include "module/system_module.h"         // 系统命令（:vars, :clear 等）
#include "io/io_module.h"                 // 输入输出操作
#include "time/time_module.h"             // 时间相关功能

// ==================== 分析模块 ====================
#include "analysis/modules/series_module.h"       // 级数展开（泰勒、傅里叶）
#include "analysis/modules/integration_module.h"   // 数值积分
#include "analysis/modules/rootfinding_module.h"   // 方程求根
#include "analysis/modules/optimization_module.h"  // 数值优化
#include "analysis/modules/analysis_module.h" // 分析相关命令
#include "analysis/modules/ode_module.h"           // 常微分方程求解

// ==================== 符号计算与多项式 ====================
#include "symbolic/modules/symbolic_module.h"  // 符号运算命令
#include "symbolic/modules/transform_module.h" // 变换（傅里叶、拉普拉斯）
#include "polynomial/calculator_polynomial.h"      // 多项式操作

// ==================== 信号处理与绘图 ====================
#include "dsp/dsp_module.h"   // 数字信号处理
#include "plot/plot_module.h" // 函数绘图

/**
 * @brief 注册所有标准模块到计算器实例
 * @param calculator 指向计算器实例的指针
 *
 * 按顺序注册所有标准模块。注册顺序不影响功能，但建议保持一致的顺序
 * 以便于调试和维护。
 */
void register_standard_modules(Calculator* calculator) {
    // -------------------- 基础数学模块 --------------------
    calculator->register_module(std::make_shared<StandardMathModule>());  // 标准数学函数
    calculator->register_module(std::make_shared<IntegerMathModule>());   // 整数运算
    calculator->register_module(std::make_shared<PreciseModule>());       // 高精度计算
    calculator->register_module(std::make_shared<StatisticsModule>());    // 统计函数

    // -------------------- 矩阵与 DSP 模块 --------------------
    calculator->register_module(std::make_shared<MatrixModule>());  // 矩阵运算
    calculator->register_module(std::make_shared<DspModule>());     // 数字信号处理
    calculator->register_module(std::make_shared<PlotModule>());    // 绘图功能

    // -------------------- 核心功能模块 --------------------
    calculator->register_module(std::make_shared<SystemModule>());                               // 系统命令
    calculator->register_module(std::make_shared<IoModule>());                                   // 输入输出
    calculator->register_module(std::make_shared<TimeModule>());                                 // 时间功能
    calculator->register_module(std::make_shared<polynomial_ops::PolynomialModule>());           // 多项式操作
    calculator->register_module(std::make_shared<series_ops::SeriesModule>());                   // 级数展开
    calculator->register_module(std::make_shared<transforms::TransformModule>());                 // 变换
    calculator->register_module(std::make_shared<integration_ops::IntegrationModule>());         // 数值积分
    calculator->register_module(std::make_shared<rootfinding::RootfindingModule>());             // 方程求根
    calculator->register_module(std::make_shared<optimization::OptimizationModule>());           // 数值优化
    calculator->register_module(std::make_shared<symbolic_commands::SymbolicModule>());          // 符号运算
    calculator->register_module(std::make_shared<analysis_cmds::AnalysisModule>());              // 分析命令
    calculator->register_module(std::make_shared<ode_ops::ODEModule>());                         // ODE 求解
}

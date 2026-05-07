// ============================================================================
// module_registration.h - 标准模块注册接口
// ============================================================================
//
// 本头文件声明了标准模块的统一注册函数。
// 注册函数负责将所有内置模块（数学、矩阵、分析、绘图等）
// 注册到计算器实例中，使这些模块的功能可用于计算器。
//
// 使用方法：
//   Calculator calc;
//   register_standard_modules(&calc);
// ============================================================================

#ifndef MODULE_REGISTRATION_H
#define MODULE_REGISTRATION_H

class Calculator;

/**
 * @brief 注册所有标准模块到计算器实例
 * @param calculator 指向计算器实例的指针
 *
 * 此函数将所有内置模块注册到指定的计算器实例，包括：
 *
 * **基础数学模块：**
 * - StandardMathModule - 标准数学函数（sin, cos, log 等）
 * - IntegerMathModule - 整数运算（阶乘、组合数等）
 * - PreciseModule - 高精度计算
 * - StatisticsModule - 统计函数
 *
 * **矩阵与 DSP 模块：**
 * - MatrixModule - 矩阵运算
 * - DspModule - 数字信号处理
 *
 * **分析模块：**
 * - SeriesModule - 级数展开
 * - IntegrationModule - 数值积分
 * - RootfindingModule - 方程求根
 * - OptimizationModule - 数值优化
 * - ODEModule - 常微分方程求解
 *
 * **符号计算与多项式：**
 * - SymbolicModule - 符号运算
 * - TransformModule - 变换（傅里叶、拉普拉斯等）
 * - PolynomialModule - 多项式操作
 *
 * **其他模块：**
 * - PlotModule - 函数绘图
 * - SystemModule - 系统命令
 * - IoModule - 输入输出
 * - TimeModule - 时间相关功能
 */
void register_standard_modules(Calculator* calculator);

#endif

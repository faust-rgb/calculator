// ============================================================================
// module_registration.h - 标准模块注册接口
// ============================================================================
//
// 本头文件声明了标准模块的统一注册函数。
// 注册函数负责将所有内置模块（数学、矩阵、分析、绘图等）
// 注册到计算器实例中，使这些模块的功能可用于计算器。
//
// 注意：ModuleRegistry 类和 REGISTER_CALCULATOR_MODULE 宏已移至 calculator_module.h
//
// 使用方法：
//   Calculator calc;
//   register_standard_modules(&calc);
// ============================================================================

#ifndef MODULE_REGISTRATION_H
#define MODULE_REGISTRATION_H

#include "module/calculator_module.h"  // 引入 ModuleRegistry

class Calculator;

/**
 * @brief 注册标准模块到计算器实例的辅助函数
 * 遍历全局注册表，并将其中的所有模块注册到 Calculator 中。
 */
void register_standard_modules(Calculator* calculator);

#endif

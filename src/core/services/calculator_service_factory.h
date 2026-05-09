// ============================================================================
// 核心服务工厂
// ============================================================================
//
// 本模块提供 CoreServices 对象的构建逻辑，将服务构造与命令处理解耦。
// CoreServices 汇总了计算器核心提供的所有服务接口，供各模块使用。
// ============================================================================

#ifndef CALCULATOR_SERVICE_FACTORY_H
#define CALCULATOR_SERVICE_FACTORY_H

#include "module/calculator_module.h"
#include "core/calculator.h"
#include "core/calculator_internal_types.h"

namespace core {

/**
 * @brief 构建 CoreServices 对象的工厂函数
 *
 * @param calculator 计算器实例
 * @param impl 计算器实现对象
 * @return 配置完整的 CoreServices 对象
 *
 * 将服务构造逻辑与命令处理逻辑解耦，提高代码复用性并减少核心类的复杂性。
 */
CoreServices build_core_services(Calculator* calculator, Calculator::Impl* impl);

} // namespace core

#endif

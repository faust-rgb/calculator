#ifndef CALCULATOR_SERVICE_FACTORY_H
#define CALCULATOR_SERVICE_FACTORY_H

#include "calculator_module.h"
#include "calculator.h"
#include "calculator_internal_types.h"

namespace core {

/**
 * @brief 构建 CoreServices 对象的工厂函数
 * 
 * 将服务构造逻辑与命令处理逻辑解耦，提高代码复用性并减少核心类的复杂性。
 */
CoreServices build_core_services(Calculator* calculator, Calculator::Impl* impl);

} // namespace core

#endif

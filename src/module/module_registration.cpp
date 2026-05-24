// ============================================================================
// module_registration.cpp - 标准模块注册实现
// ============================================================================

#include "module_registration.h"
#include "core/api/calculator.h"

void register_standard_modules(Calculator* calculator) {
    for (const auto& factory : ModuleRegistry::instance().factories()) {
        calculator->register_module(factory());
    }
}

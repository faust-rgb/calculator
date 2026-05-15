// ============================================================================
// script_context.h - 脚本执行上下文
// ============================================================================

#ifndef EXECUTION_SCRIPT_CONTEXT_H
#define EXECUTION_SCRIPT_CONTEXT_H

#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"

/**
 * @struct ScriptContext
 * @brief 汇总脚本执行所需的全部服务
 */
struct ScriptContext {
    std::shared_ptr<IVariableManager> variables;
    std::shared_ptr<IFunctionManager> functions;
    std::shared_ptr<IConfigManager> config;
    std::shared_ptr<ICommandRegistry> commands;
    std::shared_ptr<IEvaluationEngine> engine;

    static ScriptContext from_locator(const ServiceLocator& locator) {
        return {
            locator.resolve<IVariableManager>(),
            locator.resolve<IFunctionManager>(),
            locator.resolve<IConfigManager>(),
            locator.resolve<ICommandRegistry>(),
            locator.resolve<IEvaluationEngine>()
        };
    }
};

#endif // EXECUTION_SCRIPT_CONTEXT_H

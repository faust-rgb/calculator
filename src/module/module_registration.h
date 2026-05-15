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

#include <vector>
#include <functional>
#include <memory>

class Calculator;
class CalculatorModule;

/**
 * @class ModuleRegistry
 * @brief 全局模块注册表，用于去中心化注册
 */
class ModuleRegistry {
public:
    using ModuleFactory = std::function<std::shared_ptr<CalculatorModule>()>;

    static ModuleRegistry& instance() {
        static ModuleRegistry registry;
        return registry;
    }

    void register_factory(ModuleFactory factory) {
        factories_.push_back(std::move(factory));
    }

    const std::vector<ModuleFactory>& factories() const {
        return factories_;
    }

private:
    std::vector<ModuleFactory> factories_;
};

/**
 * @brief 注册标准模块到计算器实例的辅助函数
 * 遍历全局注册表，并将其中的所有模块注册到 Calculator 中。
 */
void register_standard_modules(Calculator* calculator);

/**
 * @brief 模块注册宏，放在模块的 .cpp 文件末尾
 * 
 * 示例：
 * REGISTER_CALCULATOR_MODULE(MyCustomModule)
 */
#define REGISTER_CALCULATOR_MODULE(ModuleClass) \
    namespace { \
        struct ModuleRegistrar { \
            ModuleRegistrar() { \
                ModuleRegistry::instance().register_factory([]() { \
                    return std::make_shared<ModuleClass>(); \
                }); \
            } \
        }; \
        static ModuleRegistrar global_module_registrar; \
    }

#endif

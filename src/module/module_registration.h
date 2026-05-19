// ============================================================================
// module_registration.h - 标准模块注册接口
// ============================================================================
//
// 本头文件声明了标准模块的统一注册函数和注册宏。
// 注册函数负责将所有内置模块（数学、矩阵、分析、绘图等）
// 注册到计算器实例中，使这些模块的功能可用于计算器。
//
// 特性：
// - 支持模块依赖声明和拓扑排序初始化
// - 解决静态初始化顺序问题
// - 避免命名冲突
//
// 使用方法：
//   Calculator calc;
//   register_standard_modules(&calc);
//
// 在模块 .cpp 文件末尾：
//   REGISTER_CALCULATOR_MODULE(MyModule)
// ============================================================================

#ifndef MODULE_REGISTRATION_H
#define MODULE_REGISTRATION_H

#include <vector>
#include <functional>
#include <memory>
#include <string>
#include <map>
#include <set>

class Calculator;
class CalculatorModule;
struct ModuleMetadata;

/**
 * @struct ModuleRegistrationInfo
 * @brief 模块注册信息，包含工厂和元数据
 */
struct ModuleRegistrationInfo {
    std::function<std::shared_ptr<CalculatorModule>()> factory;
    std::string unique_id;  ///< 唯一标识符，避免命名冲突
    int registration_order; ///< 注册顺序，用于稳定排序
};

/**
 * @class ModuleRegistry
 * @brief 全局模块注册表，用于去中心化注册
 *
 * 支持模块依赖解析和拓扑排序初始化。
 */
class ModuleRegistry {
public:
    using ModuleFactory = std::function<std::shared_ptr<CalculatorModule>()>;

    static ModuleRegistry& instance() {
        static ModuleRegistry registry;
        return registry;
    }

    /**
     * @brief 注册模块工厂
     * @param unique_id 模块唯一标识符（避免命名冲突）
     * @param factory 模块工厂函数
     * @return 注册序号
     */
    int register_factory(const std::string& unique_id, ModuleFactory factory) {
        int order = next_order_++;
        registrations_.push_back({factory, unique_id, order});
        return order;
    }

    /**
     * @brief 获取所有注册信息
     */
    const std::vector<ModuleRegistrationInfo>& registrations() const {
        return registrations_;
    }

    /**
     * @brief 清除所有注册（用于测试）
     */
    void clear() {
        registrations_.clear();
        next_order_ = 0;
    }

private:
    ModuleRegistry() = default;
    std::vector<ModuleRegistrationInfo> registrations_;
    int next_order_ = 0;
};

/**
 * @brief 注册标准模块到计算器实例
 *
 * 此函数会：
 * 1. 收集所有已注册的模块
 * 2. 解析模块依赖关系
 * 3. 按拓扑顺序初始化模块
 *
 * @param calculator 计算器实例指针
 */
void register_standard_modules(Calculator* calculator);

/**
 * @brief 解析模块依赖并返回拓扑排序后的初始化顺序
 * @param modules 模块列表
 * @return 按依赖顺序排列的模块索引
 */
std::vector<size_t> resolve_initialization_order(
    const std::vector<std::shared_ptr<CalculatorModule>>& modules);

/**
 * @brief 模块注册宏，放在模块的 .cpp 文件末尾
 *
 * 使用 __COUNTER__ 生成唯一标识符，避免命名空间导致的标识符问题。
 *
 * 示例：
 * REGISTER_CALCULATOR_MODULE(MyCustomModule)
 *
 * 对于带命名空间的类：
 * REGISTER_CALCULATOR_MODULE(optimization::OptimizationModule)
 */
#define REGISTER_CALCULATOR_MODULE(ModuleClass) \
    REGISTER_CALCULATOR_MODULE_IMPL(ModuleClass, __COUNTER__)

#define REGISTER_CALCULATOR_MODULE_IMPL(ModuleClass, Counter) \
    namespace { \
        namespace module_registration_detail_##Counter { \
            struct Registrar { \
                Registrar() { \
                    ModuleRegistry::instance().register_factory(#ModuleClass, []() { \
                        return std::make_shared<ModuleClass>(); \
                    }); \
                } \
            }; \
            static Registrar instance; \
        } \
    }

#endif
// ============================================================================
// service_registry.h - 动态服务注册表
// ============================================================================
//
// 允许模块动态注册自己的服务回调，而不是在工厂中硬编码。
// 这解决了服务工厂的中心化耦合问题。
// ============================================================================

#ifndef CORE_SERVICE_REGISTRY_H
#define CORE_SERVICE_REGISTRY_H

#include "core/services/service_interfaces.h"
#include "core/api/calculator.h"
#include "core/api/calculator_impl.h"
#include <functional>
#include <map>
#include <string>
#include <vector>

/**
 * @class ServiceRegistry
 * @brief 动态服务注册表
 *
 * 模块可以通过此注册表注册自己的服务回调，而不是在
 * calculator_service_factory.cpp 中硬编码所有服务。
 */
class ServiceRegistry {
public:
    // 服务回调类型
    using ServiceCallback = std::function<void(CoreServices&)>;
    using ServiceBuilder = std::function<void(CoreServices&, Calculator*, Calculator::Impl*)>;

    /**
     * @brief 注册服务构建器
     * @param module_name 模块名称（用于日志）
     * @param builder 服务构建回调
     */
    void register_builder(const std::string& module_name, ServiceBuilder builder);

    /**
     * @brief 构建所有注册的服务
     * @param services 核心服务对象
     * @param calculator 计算器实例
     * @param impl 计算器实现
     */
    void build_all(CoreServices& services, Calculator* calculator, Calculator::Impl* impl);

    /**
     * @brief 获取所有已注册的模块名
     */
    std::vector<std::string> get_registered_modules() const;

private:
    std::vector<std::pair<std::string, ServiceBuilder>> builders_;
};

// 全局服务注册表单例
ServiceRegistry& global_service_registry();

/**
 * @brief 注册服务构建器的宏
 * @param module_name 模块名称
 * @param builder_func 构建函数
 */
#define REGISTER_SERVICE_BUILDER(module_name, builder_func) \
    namespace { \
        static bool _service_registered_##module_name = []() { \
            global_service_registry().register_builder(#module_name, builder_func); \
            return true; \
        }(); \
    }

#endif // CORE_SERVICE_REGISTRY_H
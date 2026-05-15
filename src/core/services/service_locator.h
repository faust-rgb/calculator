// ============================================================================
// service_locator.h - 服务定位器模式实现
// ============================================================================
//
// 提供依赖注入与服务发现的基础设施，用于解耦核心组件与功能模块。
// 各组件通过 ServiceLocator 注册和查询所需的接口，而不是直接引用具体实现。
// ============================================================================

#ifndef CORE_SERVICE_LOCATOR_H
#define CORE_SERVICE_LOCATOR_H

#include <any>
#include <functional>
#include <memory>
#include <string>
#include <typeindex>
#include <unordered_map>
#include <stdexcept>

/**
 * @class ServiceLocator
 * @brief 核心服务发现机制
 *
 * 支持通过类型索引或名称注册和解析服务。
 */
class ServiceLocator {
public:
    /**
     * @brief 注册服务实例
     * @tparam T 接口类型
     * @param instance 接口实现指针
     */
    template <typename T>
    void register_service(std::shared_ptr<T> instance) {
        services_[typeid(T)] = instance;
    }

    /**
     * @brief 获取服务实例
     * @tparam T 接口类型
     * @return 服务实例指针，如果未注册则抛出异常
     */
    template <typename T>
    std::shared_ptr<T> resolve() const {
        auto it = services_.find(typeid(T));
        if (it == services_.end()) {
            throw std::runtime_error(std::string("Service not registered: ") + typeid(T).name());
        }
        return std::any_cast<std::shared_ptr<T>>(it->second);
    }

    /**
     * @brief 检查服务是否已注册
     * @tparam T 接口类型
     * @return 如果已注册返回 true
     */
    template <typename T>
    bool has_service() const {
        return services_.find(typeid(T)) != services_.end();
    }

    /**
     * @brief 清除所有服务
     */
    void clear() {
        services_.clear();
    }

private:
    std::unordered_map<std::type_index, std::any> services_;
};

#endif // CORE_SERVICE_LOCATOR_H

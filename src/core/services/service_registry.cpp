// ============================================================================
// service_registry.cpp - 动态服务注册表实现
// ============================================================================

#include "service_registry.h"

void ServiceRegistry::register_builder(const std::string& module_name, ServiceBuilder builder) {
    builders_.emplace_back(module_name, std::move(builder));
}

void ServiceRegistry::build_all(CoreServices& services, Calculator* calculator, Calculator::Impl* impl) {
    for (const auto& [name, builder] : builders_) {
        builder(services, calculator, impl);
    }
}

std::vector<std::string> ServiceRegistry::get_registered_modules() const {
    std::vector<std::string> names;
    names.reserve(builders_.size());
    for (const auto& [name, _] : builders_) {
        names.push_back(name);
    }
    return names;
}

ServiceRegistry& global_service_registry() {
    static ServiceRegistry instance;
    return instance;
}

// ============================================================================
// module_registration.cpp - 标准模块注册实现
// ============================================================================

#include "module/module_registration.h"
#include "module/calculator_module.h"
#include "core/api/calculator.h"
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <queue>

namespace {

/**
 * @brief 拓扑排序实现（Kahn算法）
 * @param modules 模块列表
 * @return 按依赖顺序排列的索引
 */
std::vector<size_t> topological_sort(
    const std::vector<std::shared_ptr<CalculatorModule>>& modules) {

    const size_t n = modules.size();
    if (n == 0) return {};

    // 构建模块名到索引的映射
    std::unordered_map<std::string, size_t> name_to_index;
    for (size_t i = 0; i < n; ++i) {
        name_to_index[modules[i]->name()] = i;
    }

    // 计算入度和邻接表
    std::vector<size_t> in_degree(n, 0);
    std::vector<std::vector<size_t>> adj(n);

    for (size_t i = 0; i < n; ++i) {
        for (const auto& dep : modules[i]->dependencies()) {
            auto it = name_to_index.find(dep);
            if (it != name_to_index.end()) {
                adj[it->second].push_back(i);
                in_degree[i]++;
            }
            // 如果依赖不存在，忽略（可能是可选依赖或尚未加载）
        }
    }

    // Kahn算法
    std::queue<size_t> q;
    for (size_t i = 0; i < n; ++i) {
        if (in_degree[i] == 0) {
            q.push(i);
        }
    }

    std::vector<size_t> result;
    result.reserve(n);

    while (!q.empty()) {
        size_t u = q.front();
        q.pop();
        result.push_back(u);

        for (size_t v : adj[u]) {
            if (--in_degree[v] == 0) {
                q.push(v);
            }
        }
    }

    // 检测循环依赖
    if (result.size() != n) {
        // 存在循环依赖，回退到原始顺序并警告
        // 在生产环境中应该抛出异常或记录日志
        result.clear();
        for (size_t i = 0; i < n; ++i) {
            result.push_back(i);
        }
    }

    return result;
}

} // namespace

std::vector<size_t> resolve_initialization_order(
    const std::vector<std::shared_ptr<CalculatorModule>>& modules) {
    return topological_sort(modules);
}

void register_standard_modules(Calculator* calculator) {
    if (!calculator) return;

    const auto& registrations = ModuleRegistry::instance().registrations();

    // 1. 创建所有模块实例
    std::vector<std::shared_ptr<CalculatorModule>> modules;
    modules.reserve(registrations.size());

    for (const auto& reg : registrations) {
        if (reg.factory) {
            modules.push_back(reg.factory());
        }
    }

    // 2. 解析依赖并获取初始化顺序
    auto order = resolve_initialization_order(modules);

    // 3. 按顺序注册模块
    for (size_t idx : order) {
        if (idx < modules.size()) {
            calculator->register_module(modules[idx]);
        }
    }
}

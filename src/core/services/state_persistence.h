// ============================================================================
// state_persistence.h - 状态持久化服务
// ============================================================================
//
// 提供计算器状态（变量、函数）的保存和加载功能。
// 从 calculator_core.cpp 中提取，实现职责分离。
// ============================================================================

#ifndef CORE_STATE_PERSISTENCE_H
#define CORE_STATE_PERSISTENCE_H

#include "core_manager_interfaces.h"
#include "types/stored_value.h"
#include "types/function.h"
#include <string>
#include <map>
#include <memory>

// 前向声明
class Calculator;

/**
 * @class StatePersistenceService
 * @brief 状态持久化服务实现
 *
 * 使用管理器接口进行状态持久化，避免直接访问 Calculator::Impl。
 */
class StatePersistenceService : public IStatePersistence {
public:
    /**
     * @brief 构造函数
     * @param variables 变量管理器
     * @param functions 函数管理器
     */
    StatePersistenceService(
        std::shared_ptr<IVariableManager> variables,
        std::shared_ptr<IFunctionManager> functions);

    std::string save_state(const std::string& path) const override;
    std::string load_state(const std::string& path) override;

    /**
     * @brief 设置脚本执行器（用于加载脚本函数）
     * @param executor 脚本执行回调
     */
    void set_script_executor(std::function<void(const std::string&)> executor);

private:
    // 编码/解码辅助函数
    static std::string encode_field(const std::string& text);
    static std::string decode_field(const std::string& text);
    static std::vector<std::string> split_tab_fields(const std::string& row_text);

    std::shared_ptr<IVariableManager> variables_;
    std::shared_ptr<IFunctionManager> functions_;
    std::function<void(const std::string&)> script_executor_;
};

#endif // CORE_STATE_PERSISTENCE_H

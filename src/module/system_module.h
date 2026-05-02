// ============================================================================
// 系统模块 - 提供核心系统命令
// ============================================================================
//
// 本模块实现计算器的核心系统命令，包括：
// - :help - 获取帮助信息
// - :vars - 列出变量
// - :clear - 清除变量
// - :save/:load - 状态持久化
// - :export - 数据导出
// ============================================================================

#ifndef SYSTEM_MODULE_H
#define SYSTEM_MODULE_H

#include "calculator_module.h"

#include <string>
#include <vector>

/**
 * @class SystemModule
 * @brief 系统命令模块，提供核心系统级命令
 *
 * 继承自 CalculatorModule，实现系统管理相关命令。
 */
class SystemModule : public CalculatorModule {
public:
    /// 返回模块名称
    std::string name() const override;

    /// 返回支持的命令列表
    std::vector<std::string> get_commands() const override;

    /// 执行命令
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;

    /// 返回帮助文本
    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
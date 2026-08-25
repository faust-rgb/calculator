// ============================================================================
// system_module.h - 系统命令模块头文件
// ============================================================================
//
// 本头文件声明了 SystemModule 类，该模块提供计算器的核心系统命令。
// 系统命令用于管理计算器的状态、配置和持久化功能。
//
// 支持的命令包括：
// - :vars, :funcs - 列出变量和函数
// - :clear, :clearfunc, :clearfuncs - 清除变量和函数
// - :save, :load, :export - 状态持久化
// - :run - 执行脚本文件
// - :exact, :symbolic, :precision, :scale - 配置设置
// - :hexprefix, :hexcase - 十六进制格式设置
// - print - 打印值
// ============================================================================

#ifndef SYSTEM_MODULE_H
#define SYSTEM_MODULE_H

#include "calculator_module.h"

#include <string>
#include <vector>

class ServiceLocator;

/**
 * @class SystemModule
 * @brief 系统命令模块，提供核心系统级命令
 *
 * SystemModule 继承自 CalculatorModule，实现计算器的系统管理功能。
 * 这些功能包括：
 * - 变量和函数的查看与清除
 * - 状态的保存与加载
 * - 计算模式的配置（精确模式、符号模式等）
 * - 输出格式的设置（精度、十六进制格式等）
 *
 * 系统模块是计算器的核心组件，通常在初始化时最先注册。
 */
class SystemModule : public CommandModuleBase {
public:
    /**
     * @brief 返回模块名称
     * @return 模块名称字符串 "System"
     */
    std::string name() const override;

    /**
     * @brief 返回支持的命令列表
     * @return 命令名字符串向量
     *
     * 包括所有系统管理相关的命令。
     */
    std::vector<std::string> get_commands() const override;

    /**
     * @brief 执行系统命令
     * @param command 命令名（如 ":vars", ":clear"）
     * @param args 命令参数
     * @param locator 服务定位器
     * @return 命令执行结果字符串
     *
     * 根据命令名分发到相应的处理函数。
     */
    std::string execute_args_view(std::string_view command,
                                  const std::vector<std::string_view>& args,
                                  ServiceLocator& locator) override;
    std::vector<std::string> get_help_topics() const override;

    /**
     * @brief 返回指定主题的帮助文本
     * @param topic 帮助主题（如 "commands", "variables", "persistence"）
     * @return 帮助文本字符串
     *
     * 支持的主题：commands, variables, persistence, exact, examples
     */
    std::string get_help_snippet(const std::string& topic) const override;
};

#endif
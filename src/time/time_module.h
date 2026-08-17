/**
 * @file time_module.h
 * @brief 时间模块头文件
 *
 * 本文件定义了 TimeModule 类，提供计算器的时间相关功能。
 * 包括获取当前时间、格式化时间、解析时间字符串、计时器等功能。
 *
 * @author Calculator Team
 * @date 2024
 */

#ifndef TIME_MODULE_H
#define TIME_MODULE_H

#include "module/calculator_module.h"
#include "core/value/stored_value.h"

#include <chrono>
#include <map>

/**
 * @class TimeModule
 * @brief 提供时间相关功能的计算器模块
 *
 * 该模块实现了多种时间相关的函数，包括：
 * - 时间戳获取（now, clock）
 * - 时间格式化（time, strftime）
 * - 时间解析（strptime）
 * - 执行控制（sleep）
 * - 计时器管理（timer_start, timer_elapsed, timer_stop）
 *
 * 支持的函数：
 * - now() - 获取当前时间戳（秒，Unix 时间戳）
 * - time() - 获取当前时间字符串（默认格式 YYYY-MM-DD HH:MM:SS）
 * - strftime(format) - 按指定格式返回当前时间字符串
 * - strptime(time_string, format) - 解析时间字符串为时间戳
 * - clock() - 获取高精度计时器值（秒，用于性能测量）
 * - sleep(seconds) - 暂停执行指定秒数
 * - timer_start() - 启动一个计时器，返回计时器 ID
 * - timer_elapsed(id) - 获取计时器经过的时间（秒）
 * - timer_stop(id) - 停止计时器并返回经过的时间（秒）
 *
 * 时间格式说明（strftime/strptime）：
 * - %Y: 四位年份
 * - %m: 两位月份 (01-12)
 * - %d: 两位日期 (01-31)
 * - %H: 24小时制小时 (00-23)
 * - %M: 分钟 (00-59)
 * - %S: 秒 (00-59)
 * - %s: Unix 时间戳（秒）
 */
class TimeModule : public CalculatorModule {
public:
    /**
     * @brief 获取模块名称
     * @return 返回模块名称 "Time"
     */
    std::string name() const override { return "Time"; }

    /**
     * @brief 获取模块支持的命令列表
     * @return 返回所有支持的命令名称列表
     */
    std::vector<std::string> get_commands() const override {
        return {"now", "time", "strftime", "strptime", "clock",
                "sleep", "timer_start", "timer_elapsed", "timer_stop"};
    }

    /**
     * @brief 获取模块提供的原生函数映射
     * @return 返回函数名到函数实现的映射表
     */
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_functions_map() const override;

    /**
     * @brief 执行命令行参数形式的时间命令
     * @param command 命令名称
     * @param args 参数列表
     * @param locator 服务定位器
     * @return 返回命令执行结果的字符串表示
     */
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ServiceLocator& locator) override;

    /**
     * @brief 获取帮助信息片段
     * @param topic 帮助主题
     * @return 返回对应主题的帮助文本
     */
    std::string get_help_snippet(const std::string& topic) const override;

private:
    /** @brief 计时器映射表，存储计时器ID到启动时间的映射 */
    mutable std::map<int, std::chrono::steady_clock::time_point> timers_;
    /** @brief 下一个计时器ID，用于生成唯一的计时器标识 */
    mutable int next_timer_id_ = 1;
};

#endif // TIME_MODULE_H

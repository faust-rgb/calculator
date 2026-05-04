#ifndef TIME_MODULE_H
#define TIME_MODULE_H

#include "module/calculator_module.h"
#include "types/stored_value.h"

#include <chrono>
#include <map>

/**
 * @class TimeModule
 * @brief 提供时间相关功能的模块
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
    std::string name() const override { return "Time"; }

    std::vector<std::string> get_commands() const override {
        return {"now", "time", "strftime", "strptime", "clock",
                "sleep", "timer_start", "timer_elapsed", "timer_stop"};
    }

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_native_functions() const override;

    std::string get_help_snippet(const std::string& topic) const override;

private:
    mutable std::map<int, std::chrono::steady_clock::time_point> timers_;
    mutable int next_timer_id_ = 1;
};

#endif // TIME_MODULE_H

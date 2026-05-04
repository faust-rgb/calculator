#include "time_module.h"
#include "module/calculator_module.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

namespace {

double get_scalar(const StoredValue& val, const char* context) {
    if (val.is_matrix || val.is_complex || val.is_string || val.is_list || val.is_dict) {
        throw std::runtime_error(std::string(context) + " expects a scalar value");
    }
    return val.exact ? rational_to_double(val.rational) : val.decimal;
}

std::string format_time(const std::string& format, std::time_t timestamp) {
    std::tm* tm_info = std::localtime(&timestamp);
    if (!tm_info) {
        throw std::runtime_error("Failed to get local time");
    }

    // 处理自定义格式说明符
    std::string result;
    std::size_t i = 0;
    while (i < format.size()) {
        if (format[i] == '%' && i + 1 < format.size()) {
            char spec = format[i + 1];
            switch (spec) {
                case 'Y':
                    result += std::to_string(1900 + tm_info->tm_year);
                    break;
                case 'm':
                    result += (tm_info->tm_mon + 1 < 10 ? "0" : "") +
                              std::to_string(tm_info->tm_mon + 1);
                    break;
                case 'd':
                    result += (tm_info->tm_mday < 10 ? "0" : "") +
                              std::to_string(tm_info->tm_mday);
                    break;
                case 'H':
                    result += (tm_info->tm_hour < 10 ? "0" : "") +
                              std::to_string(tm_info->tm_hour);
                    break;
                case 'M':
                    result += (tm_info->tm_min < 10 ? "0" : "") +
                              std::to_string(tm_info->tm_min);
                    break;
                case 'S':
                    result += (tm_info->tm_sec < 10 ? "0" : "") +
                              std::to_string(tm_info->tm_sec);
                    break;
                case 's':
                    result += std::to_string(timestamp);
                    break;
                case '%':
                    result += '%';
                    break;
                default:
                    // 未知格式说明符，原样保留
                    result += format.substr(i, 2);
                    break;
            }
            i += 2;
        } else {
            result += format[i];
            ++i;
        }
    }
    return result;
}

std::time_t parse_time(const std::string& time_str, const std::string& format) {
    // 简单解析：仅支持常见格式
    std::tm tm_info = {};
    tm_info.tm_year = 70;  // 默认 1970
    tm_info.tm_mday = 1;   // 默认 1 日

    std::size_t ti = 0, fi = 0;
    while (ti < time_str.size() && fi < format.size()) {
        if (format[fi] == '%' && fi + 1 < format.size()) {
            char spec = format[fi + 1];
            switch (spec) {
                case 'Y': {
                    // 四位年份
                    if (ti + 4 <= time_str.size()) {
                        int year = std::stoi(time_str.substr(ti, 4));
                        tm_info.tm_year = year - 1900;
                        ti += 4;
                    }
                    break;
                }
                case 'm': {
                    // 两位月份
                    if (ti + 2 <= time_str.size()) {
                        int month = std::stoi(time_str.substr(ti, 2));
                        tm_info.tm_mon = month - 1;
                        ti += 2;
                    }
                    break;
                }
                case 'd': {
                    // 两位日期
                    if (ti + 2 <= time_str.size()) {
                        tm_info.tm_mday = std::stoi(time_str.substr(ti, 2));
                        ti += 2;
                    }
                    break;
                }
                case 'H': {
                    // 小时
                    if (ti + 2 <= time_str.size()) {
                        tm_info.tm_hour = std::stoi(time_str.substr(ti, 2));
                        ti += 2;
                    }
                    break;
                }
                case 'M': {
                    // 分钟
                    if (ti + 2 <= time_str.size()) {
                        tm_info.tm_min = std::stoi(time_str.substr(ti, 2));
                        ti += 2;
                    }
                    break;
                }
                case 'S': {
                    // 秒
                    if (ti + 2 <= time_str.size()) {
                        tm_info.tm_sec = std::stoi(time_str.substr(ti, 2));
                        ti += 2;
                    }
                    break;
                }
                case 's': {
                    // Unix 时间戳
                    std::size_t end = ti;
                    while (end < time_str.size() && std::isdigit(time_str[end])) ++end;
                    if (end > ti) {
                        return static_cast<std::time_t>(std::stoll(time_str.substr(ti, end - ti)));
                    }
                    break;
                }
                default:
                    break;
            }
            fi += 2;
        } else {
            // 字面字符匹配
            if (time_str[ti] == format[fi]) {
                ++ti;
            }
            ++fi;
        }
    }

    return std::mktime(&tm_info);
}

} // namespace

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> TimeModule::get_native_functions() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    // now() - 获取当前 Unix 时间戳（秒）
    funcs["now"] = [](const std::vector<StoredValue>& /*args*/) -> StoredValue {
        auto now = std::chrono::system_clock::now();
        auto timestamp = std::chrono::system_clock::to_time_t(now);
        StoredValue res;
        res.decimal = static_cast<double>(timestamp);
        res.exact = false;
        return res;
    };

    // time() - 获取当前时间字符串
    funcs["time"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        auto now = std::chrono::system_clock::now();
        auto timestamp = std::chrono::system_clock::to_time_t(now);

        std::string fmt = "%Y-%m-%d %H:%M:%S";
        if (!args.empty() && args[0].is_string) {
            fmt = args[0].string_value;
        }

        StoredValue res;
        res.is_string = true;
        res.string_value = format_time(fmt, timestamp);
        return res;
    };

    // strftime(format, [timestamp]) - 格式化时间
    funcs["strftime"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) {
            throw std::runtime_error("strftime expects at least 1 argument (format)");
        }
        if (!args[0].is_string) {
            throw std::runtime_error("strftime format must be a string");
        }

        std::string format = args[0].string_value;
        std::time_t timestamp;

        if (args.size() > 1) {
            double ts = get_scalar(args[1], "strftime timestamp");
            timestamp = static_cast<std::time_t>(ts);
        } else {
            auto now = std::chrono::system_clock::now();
            timestamp = std::chrono::system_clock::to_time_t(now);
        }

        StoredValue res;
        res.is_string = true;
        res.string_value = format_time(format, timestamp);
        return res;
    };

    // strptime(time_string, format) - 解析时间字符串为时间戳
    funcs["strptime"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2) {
            throw std::runtime_error("strptime expects 2 arguments (time_string, format)");
        }
        if (!args[0].is_string || !args[1].is_string) {
            throw std::runtime_error("strptime arguments must be strings");
        }

        std::time_t timestamp = parse_time(args[0].string_value, args[1].string_value);

        StoredValue res;
        res.decimal = static_cast<double>(timestamp);
        res.exact = false;
        return res;
    };

    // clock() - 高精度计时器（秒）
    funcs["clock"] = [](const std::vector<StoredValue>& /*args*/) -> StoredValue {
        auto now = std::chrono::steady_clock::now();
        auto duration = now.time_since_epoch();
        double seconds = std::chrono::duration<double>(duration).count();

        StoredValue res;
        res.decimal = seconds;
        res.exact = false;
        return res;
    };

    // sleep(seconds) - 暂停执行
    funcs["sleep"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) {
            throw std::runtime_error("sleep expects 1 argument (seconds)");
        }
        double seconds = get_scalar(args[0], "sleep duration");
        if (seconds < 0) {
            throw std::runtime_error("sleep duration must be non-negative");
        }

        std::this_thread::sleep_for(std::chrono::duration<double>(seconds));

        StoredValue res;
        res.decimal = 1.0;
        res.exact = false;
        return res;
    };

    // timer_start() - 启动计时器
    funcs["timer_start"] = [this](const std::vector<StoredValue>& /*args*/) -> StoredValue {
        int timer_id = next_timer_id_++;
        timers_[timer_id] = std::chrono::steady_clock::now();

        StoredValue res;
        res.decimal = timer_id;
        res.exact = false;
        return res;
    };

    // timer_elapsed(id) - 获取计时器经过时间
    funcs["timer_elapsed"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) {
            throw std::runtime_error("timer_elapsed expects 1 argument (timer_id)");
        }
        int timer_id = static_cast<int>(get_scalar(args[0], "timer_id"));

        auto it = timers_.find(timer_id);
        if (it == timers_.end()) {
            throw std::runtime_error("Invalid timer ID: " + std::to_string(timer_id));
        }

        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - it->second).count();

        StoredValue res;
        res.decimal = elapsed;
        res.exact = false;
        return res;
    };

    // timer_stop(id) - 停止计时器
    funcs["timer_stop"] = [this](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) {
            throw std::runtime_error("timer_stop expects 1 argument (timer_id)");
        }
        int timer_id = static_cast<int>(get_scalar(args[0], "timer_id"));

        auto it = timers_.find(timer_id);
        if (it == timers_.end()) {
            throw std::runtime_error("Invalid timer ID: " + std::to_string(timer_id));
        }

        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - it->second).count();
        timers_.erase(it);

        StoredValue res;
        res.decimal = elapsed;
        res.exact = false;
        return res;
    };

    return funcs;
}

std::string TimeModule::get_help_snippet(const std::string& topic) const {
    if (topic == "time" || topic == "timer" || topic == "now" || topic == "clock") {
        return R"(
Time Functions (time module):
  now()                - Current Unix timestamp (seconds since 1970-01-01)
  time()               - Current time as string "YYYY-MM-DD HH:MM:SS"
  time(format)         - Current time formatted (see strftime formats)
  strftime(fmt, [ts])  - Format timestamp (default: current time)
  strptime(str, fmt)   - Parse time string to timestamp
  clock()              - High-resolution timer value (seconds)
  sleep(seconds)       - Pause execution for given seconds
  timer_start()        - Start a timer, returns timer ID
  timer_elapsed(id)    - Get elapsed time for timer (seconds)
  timer_stop(id)       - Stop timer and return elapsed time

Time format specifiers:
  %Y - Year (4 digits)    %m - Month (01-12)    %d - Day (01-31)
  %H - Hour (00-23)       %M - Minute (00-59)   %S - Second (00-59)
  %s - Unix timestamp

Examples:
  now()                           # 1714833600
  time()                          # "2024-05-04 12:00:00"
  strftime("%Y-%m-%d")            # "2024-05-04"
  t = timer_start()               # Start timing
  # ... some computation ...
  elapsed = timer_stop(t)         # Get elapsed seconds
)";
    }
    return "";
}

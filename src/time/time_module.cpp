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

std::string format_time(const std::string& format, std::time_t timestamp, bool use_local = true) {
    std::tm* tm_info = use_local ? std::localtime(&timestamp) : std::gmtime(&timestamp);
    if (!tm_info) {
        throw std::runtime_error("Failed to get time info");
    }

    std::ostringstream oss;
    oss << std::put_time(tm_info, format.c_str());
    return oss.str();
}

std::time_t parse_time(const std::string& time_str, const std::string& format) {
    // 优先处理 Unix 时间戳
    if (format == "%s") {
        try {
            return static_cast<std::time_t>(std::stoll(time_str));
        } catch (...) {
            throw std::runtime_error("Failed to parse Unix timestamp: " + time_str);
        }
    }

    std::tm tm_info = {};
    std::istringstream iss(time_str);
    iss >> std::get_time(&tm_info, format.c_str());
    if (iss.fail()) {
        throw std::runtime_error("Failed to parse time string: '" + time_str + "' with format: '" + format + "'");
    }
    return std::mktime(&tm_info);
}

} // namespace

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> TimeModule::get_native_functions() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    // now() - 获取当前 Unix 时间戳（秒）
    funcs["now"] = [](const std::vector<StoredValue>& /*args*/) -> StoredValue {
        auto now = std::chrono::system_clock::now();
        auto duration = now.time_since_epoch();
        double seconds = std::chrono::duration<double>(duration).count();
        StoredValue res;
        res.decimal = seconds;
        res.exact = false;
        return res;
    };

    // time() - 获取当前本地时间字符串
    funcs["time"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        auto now = std::chrono::system_clock::now();
        auto timestamp = std::chrono::system_clock::to_time_t(now);

        std::string fmt = "%Y-%m-%d %H:%M:%S";
        if (!args.empty() && args[0].is_string) {
            fmt = args[0].string_value;
        }

        StoredValue res;
        res.is_string = true;
        res.string_value = format_time(fmt, timestamp, true);
        return res;
    };

    // utctime() - 获取当前 UTC 时间字符串
    funcs["utctime"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        auto now = std::chrono::system_clock::now();
        auto timestamp = std::chrono::system_clock::to_time_t(now);

        std::string fmt = "%Y-%m-%d %H:%M:%S";
        if (!args.empty() && args[0].is_string) {
            fmt = args[0].string_value;
        }

        StoredValue res;
        res.is_string = true;
        res.string_value = format_time(fmt, timestamp, false);
        return res;
    };

    // strftime(format, [timestamp]) - 格式化本地时间
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
        res.string_value = format_time(format, timestamp, true);
        return res;
    };

    // strftime_utc(format, [timestamp]) - 格式化 UTC 时间
    funcs["strftime_utc"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) {
            throw std::runtime_error("strftime_utc expects at least 1 argument (format)");
        }
        if (!args[0].is_string) {
            throw std::runtime_error("strftime_utc format must be a string");
        }

        std::string format = args[0].string_value;
        std::time_t timestamp;

        if (args.size() > 1) {
            double ts = get_scalar(args[1], "strftime_utc timestamp");
            timestamp = static_cast<std::time_t>(ts);
        } else {
            auto now = std::chrono::system_clock::now();
            timestamp = std::chrono::system_clock::to_time_t(now);
        }

        StoredValue res;
        res.is_string = true;
        res.string_value = format_time(format, timestamp, false);
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

/**
 * @file time_module.cpp
 * @brief 时间模块实现文件
 *
 * 本文件实现了 TimeModule 类，提供计算器的时间相关功能。
 * 包括获取当前时间、格式化时间、解析时间字符串、计时器等功能。
 *
 * @author Calculator Team
 * @date 2024
 */

#include "time_module.h"
#include "module/calculator_module.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

/**
 * @brief 匿名命名空间，包含内部辅助函数
 */
namespace {

/**
 * @brief 从 StoredValue 中提取标量值
 *
 * 检查值是否为标量类型（非矩阵、复数、字符串、列表或字典），
 * 并将其转换为 long double 类型的数值。
 *
 * @param val 要提取的存储值
 * @param context 错误消息中的上下文描述
 * @return 提取的标量值
 * @throws std::runtime_error 如果值不是标量类型
 */
long double get_scalar(const StoredValue& val, const char* context) {
    if (val.is_matrix || val.is_complex || val.is_string || val.is_list || val.is_dict) {
        throw std::runtime_error(std::string(context) + " expects a scalar value");
    }
    return val.exact ? rational_to_double(val.rational) : val.decimal;
}

/**
 * @brief 格式化时间戳为字符串
 *
 * 根据指定的格式字符串将时间戳转换为可读的时间字符串。
 * 支持本地时间和 UTC 时间。
 *
 * @param format 时间格式字符串（如 "%Y-%m-%d %H:%M:%S"）
 * @param timestamp Unix 时间戳
 * @param use_local 是否使用本地时间（true 为本地时间，false 为 UTC）
 * @return 格式化后的时间字符串
 * @throws std::runtime_error 如果获取时间信息失败
 */
std::string format_time(const std::string& format, std::time_t timestamp, bool use_local = true) {
    std::tm* tm_info = use_local ? std::localtime(&timestamp) : std::gmtime(&timestamp);
    if (!tm_info) {
        throw std::runtime_error("Failed to get time info");
    }

    std::ostringstream oss;
    oss << std::put_time(tm_info, format.c_str());
    return oss.str();
}

/**
 * @brief 解析时间字符串为时间戳
 *
 * 根据指定的格式字符串将时间字符串解析为 Unix 时间戳。
 * 特殊处理 "%s" 格式，直接解析为 Unix 时间戳。
 *
 * @param time_str 时间字符串
 * @param format 时间格式字符串
 * @return 解析后的 Unix 时间戳
 * @throws std::runtime_error 如果解析失败
 */
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

/**
 * @brief 获取模块提供的原生函数映射
 *
 * 实现并返回所有时间相关函数的映射表。
 * 每个函数都是一个 lambda 表达式，接收参数向量并返回 StoredValue。
 *
 * @return 函数名到函数实现的映射表
 */
std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> TimeModule::get_native_functions() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    // ==================== 时间戳函数 ====================

    /**
     * @brief now() - 获取当前 Unix 时间戳（秒）
     *
     * 返回自 1970-01-01 00:00:00 UTC 以来经过的秒数。
     * 使用高精度浮点数表示，可包含小数部分。
     *
     * @return 当前 Unix 时间戳（秒）
     */
    funcs["now"] = [](const std::vector<StoredValue>& /*args*/) -> StoredValue {
        auto now = std::chrono::system_clock::now();
        auto duration = now.time_since_epoch();
        long double seconds = std::chrono::duration<long double>(duration).count();
        StoredValue res;
        res.decimal = seconds;
        res.exact = false;
        return res;
    };

    // ==================== 时间格式化函数 ====================

    /**
     * @brief time() - 获取当前本地时间字符串
     *
     * 返回当前本地时间的字符串表示。
     * 可选参数指定时间格式，默认为 "YYYY-MM-DD HH:MM:SS"。
     *
     * @param args[0] 可选的时间格式字符串
     * @return 当前时间字符串
     */
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

    /**
     * @brief utctime() - 获取当前 UTC 时间字符串
     *
     * 返回当前 UTC 时间的字符串表示。
     * 可选参数指定时间格式，默认为 "YYYY-MM-DD HH:MM:SS"。
     *
     * @param args[0] 可选的时间格式字符串
     * @return 当前 UTC 时间字符串
     */
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

    /**
     * @brief strftime(format, [timestamp]) - 格式化本地时间
     *
     * 按照指定格式将时间戳转换为本地时间字符串。
     * 如果不提供时间戳，则使用当前时间。
     *
     * @param args[0] 时间格式字符串（必需）
     * @param args[1] 可选的 Unix 时间戳
     * @return 格式化后的时间字符串
     * @throws std::runtime_error 如果参数不足或类型错误
     */
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
            long double ts = get_scalar(args[1], "strftime timestamp");
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

    /**
     * @brief strftime_utc(format, [timestamp]) - 格式化 UTC 时间
     *
     * 按照指定格式将时间戳转换为 UTC 时间字符串。
     * 如果不提供时间戳，则使用当前时间。
     *
     * @param args[0] 时间格式字符串（必需）
     * @param args[1] 可选的 Unix 时间戳
     * @return 格式化后的 UTC 时间字符串
     * @throws std::runtime_error 如果参数不足或类型错误
     */
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
            long double ts = get_scalar(args[1], "strftime_utc timestamp");
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

    // ==================== 时间解析函数 ====================

    /**
     * @brief strptime(time_string, format) - 解析时间字符串为时间戳
     *
     * 根据指定格式解析时间字符串，返回 Unix 时间戳。
     *
     * @param args[0] 时间字符串（必需）
     * @param args[1] 时间格式字符串（必需）
     * @return 解析后的 Unix 时间戳（秒）
     * @throws std::runtime_error 如果参数不足、类型错误或解析失败
     */
    funcs["strptime"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2) {
            throw std::runtime_error("strptime expects 2 arguments (time_string, format)");
        }
        if (!args[0].is_string || !args[1].is_string) {
            throw std::runtime_error("strptime arguments must be strings");
        }

        std::time_t timestamp = parse_time(args[0].string_value, args[1].string_value);

        StoredValue res;
        res.decimal = static_cast<long double>(timestamp);
        res.exact = false;
        return res;
    };

    // ==================== 计时器函数 ====================

    /**
     * @brief clock() - 获取高精度计时器值（秒）
     *
     * 返回高精度计时器的当前值（秒）。
     * 主要用于性能测量，返回值本身无绝对意义，
     * 但两次调用的差值可表示经过的时间。
     *
     * @return 计时器值（秒）
     */
    funcs["clock"] = [](const std::vector<StoredValue>& /*args*/) -> StoredValue {
        auto now = std::chrono::steady_clock::now();
        auto duration = now.time_since_epoch();
        long double seconds = std::chrono::duration<long double>(duration).count();

        StoredValue res;
        res.decimal = seconds;
        res.exact = false;
        return res;
    };

    /**
     * @brief sleep(seconds) - 暂停执行
     *
     * 暂停当前线程执行指定的秒数。
     *
     * @param args[0] 暂停时间（秒，必需，必须非负）
     * @return 返回 1.0 表示执行完成
     * @throws std::runtime_error 如果参数不足、时间值为负
     */
    funcs["sleep"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) {
            throw std::runtime_error("sleep expects 1 argument (seconds)");
        }
        long double seconds = get_scalar(args[0], "sleep duration");
        if (seconds < 0) {
            throw std::runtime_error("sleep duration must be non-negative");
        }

        std::this_thread::sleep_for(std::chrono::duration<long double>(seconds));

        StoredValue res;
        res.decimal = 1.0L;
        res.exact = false;
        return res;
    };

    /**
     * @brief timer_start() - 启动计时器
     *
     * 启动一个新的计时器并返回其唯一标识符。
     * 可用于测量代码执行时间。
     *
     * @return 计时器 ID（整数）
     */
    funcs["timer_start"] = [this](const std::vector<StoredValue>& /*args*/) -> StoredValue {
        int timer_id = next_timer_id_++;
        timers_[timer_id] = std::chrono::steady_clock::now();

        StoredValue res;
        res.decimal = timer_id;
        res.exact = false;
        return res;
    };

    /**
     * @brief timer_elapsed(id) - 获取计时器经过时间
     *
     * 返回指定计时器自启动以来经过的时间（秒）。
     * 计时器继续运行，不会被停止。
     *
     * @param args[0] 计时器 ID（必需）
     * @return 经过的秒数
     * @throws std::runtime_error 如果参数不足或计时器 ID 无效
     */
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
        long double elapsed = std::chrono::duration<long double>(now - it->second).count();

        StoredValue res;
        res.decimal = elapsed;
        res.exact = false;
        return res;
    };

    /**
     * @brief timer_stop(id) - 停止计时器
     *
     * 停止指定的计时器并返回自启动以来经过的时间（秒）。
     * 计时器停止后将被删除，无法再次使用。
     *
     * @param args[0] 计时器 ID（必需）
     * @return 经过的秒数
     * @throws std::runtime_error 如果参数不足或计时器 ID 无效
     */
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
        long double elapsed = std::chrono::duration<long double>(now - it->second).count();
        timers_.erase(it);

        StoredValue res;
        res.decimal = elapsed;
        res.exact = false;
        return res;
    };

    return funcs;
}

/**
 * @brief 获取时间模块的帮助信息片段
 *
 * 根据请求的主题返回相应的帮助文本。
 * 目前支持 "time"、"timer"、"now"、"clock" 等主题。
 *
 * @param topic 帮助主题
 * @return 对应主题的帮助文本，如果主题不匹配则返回空字符串
 */
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

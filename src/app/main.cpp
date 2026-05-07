/**
 * @file main.cpp
 * @brief 计算器应用程序入口
 *
 * 提供交互式命令行界面，支持：
 * - 表达式求值
 * - 历史记录浏览（上下方向键）
 * - Tab 自动补全
 * - 脚本文件执行
 * - 交互式帮助系统
 */

#include "calculator.h"
#include "string_utils.h"

#include <algorithm>
#include <cctype>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <termios.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

/**
 * @class RawModeGuard
 * @brief 终端原始模式管理器（RAII）
 *
 * 在构造时将终端切换到原始模式（禁用行缓冲和回显），
 * 在析构时自动恢复原始设置。
 */
class RawModeGuard {
public:
    explicit RawModeGuard(int fd) : fd_(fd) {
        enabled_ = tcgetattr(fd_, &original_) == 0;
        if (!enabled_) {
            return;
        }

        termios raw = original_;
        raw.c_lflag &= static_cast<unsigned long>(~(ICANON | ECHO));
        raw.c_cc[VMIN] = 1;
        raw.c_cc[VTIME] = 0;
        enabled_ = tcsetattr(fd_, TCSAFLUSH, &raw) == 0;
    }

    ~RawModeGuard() {
        if (enabled_) {
            tcsetattr(fd_, TCSAFLUSH, &original_);
        }
    }

    bool enabled() const {
        return enabled_;
    }

private:
    int fd_;
    bool enabled_ = false;
    termios original_ {};
};

/**
 * @brief 重绘输入行
 *
 * 清除当前行并重新绘制提示符、输入内容和光标位置。
 * 用于在用户编辑输入时保持显示同步。
 *
 * @param prompt 命令行提示符字符串
 * @param line 当前输入行内容
 * @param cursor_pos 光标在输入行中的位置
 */
void redraw_input(const std::string& prompt, const std::string& line, std::size_t cursor_pos) {
    // \33[2K 清空当前行，\r 回到行首
    // 重画提示符和输入内容
    std::cout << "\33[2K\r" << prompt << line;

    // 将光标移回到正确的位置。使用 \r 后接 \33[C (向前移动)
    std::size_t offset = prompt.size() + cursor_pos;
    if (offset > 0) {
        std::cout << "\r\33[" << offset << "C";
    } else {
        std::cout << "\r";
    }
    std::cout << std::flush;
}

/**
 * @brief 格式化历史记录输出
 *
 * 将命令历史记录格式化为带编号的列表字符串。
 *
 * @param history 历史命令列表
 * @return 格式化后的历史记录字符串，每条命令带有编号
 */
std::string format_history(const std::vector<std::string>& history) {
    if (history.empty()) {
        return "No history.";
    }

    std::ostringstream out;
    for (std::size_t i = 0; i < history.size(); ++i) {
        if (i != 0) {
            out << '\n';
        }
        out << (i + 1) << ". " << history[i];
    }
    return out.str();
}

/**
 * @brief 获取帮助主题列表
 *
 * 返回计算器支持的所有帮助主题的静态列表。
 *
 * @return 帮助主题名称列表的常量引用
 */
const std::vector<std::string>& help_topics() {
    static const std::vector<std::string> topics = {
        "commands", "functions", "matrix", "symbolic", "analysis", "planning", "examples",
        "exact", "variables", "persistence", "programmer"
    };
    return topics;
}

/**
 * @brief 获取命令补全关键词列表
 *
 * 返回可用于Tab补全的核心命令名称列表。
 *
 * @return 命令名称列表的常量引用
 */
const std::vector<std::string>& command_completion_words() {
    static const std::vector<std::string> words = {
        ":help", ":exact", ":precision", ":symbolic",
        ":hexprefix", ":hexcase"
    };
    return words;
}

/**
 * @brief 提取当前正在输入的标记（token）
 *
 * 从输入行中提取光标所在位置的标记，用于Tab补全。
 * 标记由字母、数字、冒号和下划线组成。
 *
 * @param line 当前的输入行
 * @return 包含标记起始位置和标记内容的 pair
 */
std::pair<std::size_t, std::string> current_token(const std::string& line) {
    // 从行尾向前扫描，找出当前正在输入的 token。
    // 这样 Tab 补全只会替换最后一段命令/函数前缀。
    std::size_t start = line.size();
    while (start > 0) {
        const char ch = line[start - 1];
        if (std::isalnum(static_cast<unsigned char>(ch)) ||
            ch == ':' || ch == '_') {
            --start;
        } else {
            break;
        }
    }
    return {start, line.substr(start)};
}

/**
 * @brief 计算字符串列表的公共前缀
 *
 * 找出所有匹配字符串的最长公共前缀，用于Tab补全时的部分补全。
 *
 * @param matches 候选字符串列表
 * @return 最长公共前缀字符串，如果没有匹配则返回空字符串
 */
std::string common_prefix(const std::vector<std::string>& matches) {
    if (matches.empty()) {
        return "";
    }

    std::string prefix = matches.front();
    for (std::size_t i = 1; i < matches.size(); ++i) {
        std::size_t j = 0;
        while (j < prefix.size() &&
               j < matches[i].size() &&
               prefix[j] == matches[i][j]) {
            ++j;
        }
        prefix.resize(j);
        if (prefix.empty()) {
            break;
        }
    }
    return prefix;
}

/**
 * @brief 收集补全候选词
 *
 * 根据当前输入上下文收集所有可能用于Tab补全的候选词。
 * 根据输入内容不同，返回帮助主题、命令或函数名。
 *
 * @param calculator 计算器实例，用于获取模块函数和变量列表
 * @param line 当前的输入行，用于判断上下文
 * @param token 当前正在输入的标记
 * @return 候选补全词列表
 */
std::vector<std::string> gather_completion_words(const Calculator& calculator,
                                                 const std::string& line,
                                                 const std::string& token) {
    if (line.rfind(":help ", 0) == 0 && token.find(':') == std::string::npos) {
        return help_topics();
    }

    if (!token.empty() && token.front() == ':') {
        // 合并核心命令和模块提供的命令
        std::vector<std::string> words = command_completion_words();
        const std::vector<std::string> module_cmds = calculator.module_command_names();
        for (const auto& cmd : module_cmds) {
            words.push_back(":" + cmd);
        }
        std::sort(words.begin(), words.end());
        words.erase(std::unique(words.begin(), words.end()), words.end());
        return words;
    }

    // 合并核心函数、模块函数、变量和自定义函数
    std::vector<std::string> words = { "help", "exit", "quit" };

    const std::vector<std::string> module_funcs = calculator.module_function_names();
    for (const auto& name : module_funcs) {
        words.push_back(name + "(");
    }

    const std::vector<std::string> variable_names = calculator.variable_names();
    words.insert(words.end(), variable_names.begin(), variable_names.end());

    const std::vector<std::string> function_names = calculator.custom_function_names();
    for (const std::string& name : function_names) {
        words.push_back(name + "(");
    }

    std::sort(words.begin(), words.end());
    words.erase(std::unique(words.begin(), words.end()), words.end());
    return words;
}

/**
 * @brief 格式化补全候选列表
 *
 * 将多个补全候选词格式化为用户友好的显示字符串。
 *
 * @param matches 补全候选词列表
 * @return 格式化后的候选列表字符串
 */
std::string format_completion_candidates(const std::vector<std::string>& matches) {
    std::ostringstream out;
    out << "Candidates:\n";
    for (const std::string& match : matches) {
        out << "  " << match << '\n';
    }
    return out.str();
}

/**
 * @struct CompletionResult
 * @brief Tab补全结果结构体
 *
 * 存储Tab补全操作的结果信息，包括是否有匹配、
 * 是否已应用补全、是否存在歧义以及所有匹配项。
 */
struct CompletionResult {
    bool has_matches = false;    ///< 是否找到匹配项
    bool applied = false;        ///< 是否已应用补全到输入行
    bool ambiguous = false;      ///< 是否存在多个匹配项（歧义）
    std::vector<std::string> matches;  ///< 所有匹配的候选词列表
};

/**
 * @brief 应用Tab补全
 *
 * 根据当前输入内容执行Tab补全操作，查找并应用匹配的补全词。
 * 如果存在唯一匹配，则直接补全；如果存在多个匹配，则补全公共前缀。
 *
 * @param calculator 计算器实例
 * @param line 输入行的指针，补全结果会直接修改此字符串
 * @return 补全结果结构体，包含匹配状态和信息
 */
CompletionResult apply_completion(const Calculator& calculator, std::string* line) {
    CompletionResult result;

    const auto [start, token] = current_token(*line);
    if (token.empty()) {
        return result;
    }

    const std::vector<std::string> words =
        gather_completion_words(calculator, *line, token);
    for (const std::string& word : words) {
        if (word.rfind(token, 0) == 0) {
            result.matches.push_back(word);
        }
    }

    if (result.matches.empty()) {
        return result;
    }

    result.has_matches = true;
    result.ambiguous = result.matches.size() > 1;

    const std::string replacement =
        result.matches.size() == 1 ? result.matches.front() : common_prefix(result.matches);
    if (replacement.size() > token.size()) {
        line->replace(start, token.size(), replacement);
        result.applied = true;
    }
    return result;
}

/**
 * @brief 读取一行输入（支持历史记录和Tab补全）
 *
 * 交互式输入读取函数，支持：
 * - 方向键浏览历史记录
 * - Tab自动补全
 * - 光标移动和编辑
 * - Ctrl+A/E/K/D等快捷键
 *
 * @param calculator 计算器实例，用于Tab补全
 * @param prompt 命令行提示符
 * @param history 历史命令列表
 * @return 用户输入的命令行字符串
 */
std::string read_line_with_history(Calculator& calculator,
                                   const std::string& prompt,
                                   const std::vector<std::string>& history) {
    std::cout << prompt << std::flush;

    // 如果当前输入并不是交互式终端，就退回 getline，避免影响重定向输入。
    if (!isatty(STDIN_FILENO)) {
        std::string line;
        std::getline(std::cin, line);
        return line;
    }

    RawModeGuard raw_mode(STDIN_FILENO);
    if (!raw_mode.enabled()) {
        std::string line;
        std::getline(std::cin, line);
        return line;
    }

    std::string line;
    std::size_t cursor_pos = 0;
    std::size_t history_index = history.size();
    std::string last_tab_line;
    bool last_tab_was_ambiguous = false;

    while (true) {
        // 逐字符读取，这样才能自己处理方向键、Tab 和退格。
        char ch = '\0';
        if (read(STDIN_FILENO, &ch, 1) <= 0) {
            return line;
        }

        if (ch == '\n' || ch == '\r') {
            std::cout << '\n';
            return line;
        }

        if (ch == 1) { // Ctrl+A (Home)
            cursor_pos = 0;
            redraw_input(prompt, line, cursor_pos);
            continue;
        }

        if (ch == 5) { // Ctrl+E (End)
            cursor_pos = line.size();
            redraw_input(prompt, line, cursor_pos);
            continue;
        }

        if (ch == 11) { // Ctrl+K (Kill to end)
            if (cursor_pos < line.size()) {
                line.erase(cursor_pos);
                redraw_input(prompt, line, cursor_pos);
            }
            continue;
        }

        if (ch == 4) {
            // Ctrl+D: 在空行上视为 EOF，其他情况下删除光标处字符。
            if (line.empty()) {
                std::cout << '\n';
                return line;
            } else if (cursor_pos < line.size()) {
                line.erase(cursor_pos, 1);
                redraw_input(prompt, line, cursor_pos);
            }
            continue;
        }

        if (ch == 127 || ch == '\b') {
            if (cursor_pos > 0) {
                line.erase(cursor_pos - 1, 1);
                --cursor_pos;
                redraw_input(prompt, line, cursor_pos);
            }
            last_tab_was_ambiguous = false;
            continue;
        }

        if (ch == '\t') {
            const CompletionResult result = apply_completion(calculator, &line);
            if (result.applied) {
                cursor_pos = line.size();
                redraw_input(prompt, line, cursor_pos);
            }
            if (result.ambiguous && last_tab_was_ambiguous && last_tab_line == line) {
                std::cout << '\n' << format_completion_candidates(result.matches);
                redraw_input(prompt, line, cursor_pos);
                last_tab_was_ambiguous = false;
            } else {
                last_tab_was_ambiguous = result.ambiguous;
                last_tab_line = line;
            }
            continue;
        }

        last_tab_was_ambiguous = false;

        if (ch == '\33') {
            // 方向键通常会发出 ESC [ A / B / C / D 这样的序列。
            char seq[2] = {'\0', '\0'};
            if (read(STDIN_FILENO, &seq[0], 1) <= 0 ||
                read(STDIN_FILENO, &seq[1], 1) <= 0) {
                continue;
            }

            if (seq[0] == '[') {
                if (seq[1] == 'A') { // Up
                    if (!history.empty() && history_index > 0) {
                        --history_index;
                        line = history[history_index];
                        cursor_pos = line.size();
                        redraw_input(prompt, line, cursor_pos);
                    }
                } else if (seq[1] == 'B') { // Down
                    if (history_index < history.size()) {
                        ++history_index;
                        line = history_index == history.size() ? "" : history[history_index];
                        cursor_pos = line.size();
                        redraw_input(prompt, line, cursor_pos);
                    }
                } else if (seq[1] == 'C') { // Right
                    if (cursor_pos < line.size()) {
                        ++cursor_pos;
                        redraw_input(prompt, line, cursor_pos);
                    }
                } else if (seq[1] == 'D') { // Left
                    if (cursor_pos > 0) {
                        --cursor_pos;
                        redraw_input(prompt, line, cursor_pos);
                    }
                }
            }
            continue;
        }

        if (ch >= 32 && ch <= 126) {
            line.insert(cursor_pos, 1, ch);
            ++cursor_pos;
            redraw_input(prompt, line, cursor_pos);
        }
    }
}

/**
 * @brief 从输入流读取所有内容
 *
 * 将输入流中的所有内容读取到一个字符串中。
 *
 * @param in 输入流引用
 * @return 包含输入流所有内容的字符串
 */
std::string read_all(std::istream& in) {
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

/**
 * @brief 检查文件是否具有 .calc 扩展名
 *
 * 判断给定路径是否以 ".calc" 扩展名结尾。
 *
 * @param path 文件路径
 * @return 如果是 .calc 文件返回 true，否则返回 false
 */
bool has_calc_extension(const std::string& path) {
    constexpr std::string_view extension = ".calc";
    return path.size() >= extension.size() &&
           path.compare(path.size() - extension.size(), extension.size(), extension) == 0;
}

/**
 * @brief 运行脚本文件
 *
 * 执行指定路径的 .calc 脚本文件并输出结果。
 *
 * @param calculator 计算器实例
 * @param script_path 脚本文件路径
 * @param exact_mode 是否启用精确模式
 * @return 成功返回 0，失败返回 1
 */
int run_script_file(Calculator& calculator, const std::string& script_path, bool exact_mode) {
    if (!has_calc_extension(script_path)) {
        std::cerr << "Error: script file must use .calc extension\n";
        return 1;
    }

    try {
        const std::string output = calculator.execute_script_file(script_path, exact_mode);
        if (!output.empty()) {
            std::cout << output << '\n';
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}

/**
 * @brief 执行REPL单行命令
 *
 * 处理REPL中输入的单行命令，包括内置命令（help、exit等）
 * 和计算器表达式求值。
 *
 * @param calculator 计算器实例
 * @param raw_line 原始输入行
 * @param exact_mode 精确模式标志的指针，可能被修改
 * @param history 历史命令列表
 * @param should_exit 退出标志的指针，当用户输入exit时设为true
 * @return 命令执行结果字符串
 */
std::string execute_repl_line(Calculator& calculator,
                              const std::string& raw_line,
                              bool* exact_mode,
                              const std::vector<std::string>& history,
                              bool* should_exit) {
    const std::string line = utils::trim_copy(raw_line);
    if (line.empty()) {
        return "";
    }

    if (line == "exit" || line == "quit") {
        *should_exit = true;
        return "";
    }

    if (line == "help" || line == ":help") {
        return calculator.help_text();
    }
    if (line.rfind(":help ", 0) == 0) {
        return calculator.help_topic(utils::trim_copy(line.substr(6)));
    }

    // 处理特殊 REPL 本地指令
    if (line == ":history") {
        return format_history(history);
    }

    // 拦截 :exact 以便同步 REPL 的本地 exact_mode 状态
    if (line == ":exact on") { *exact_mode = true; }
    else if (line == ":exact off") { *exact_mode = false; }

    std::string output;
    if (calculator.try_process_function_command(line, &output, *exact_mode)) {
        return output;
    }

    return calculator.process_line(line, *exact_mode);
}

/**
 * @brief 运行非交互式输入
 *
 * 处理来自管道或重定向的非交互式输入，逐行执行命令。
 *
 * @param calculator 计算器实例
 * @param input 输入字符串
 * @return 成功返回 0，失败返回 1
 */
int run_non_interactive_input(Calculator& calculator, const std::string& input) {
    std::istringstream stream(input);
    std::string line;
    std::vector<std::string> history;
    bool exact_mode = false;

    while (std::getline(stream, line)) {
        const std::string trimmed = utils::trim_copy(line);
        if (trimmed.empty()) {
            continue;
        }

        bool should_exit = false;
        try {
            const std::string output =
                execute_repl_line(calculator, trimmed, &exact_mode, history, &should_exit);
            if (!should_exit) {
                history.push_back(trimmed);
            }
            if (!output.empty()) {
                std::cout << output << '\n';
            }
            if (should_exit) {
                break;
            }
        } catch (const std::exception& ex) {
            std::cerr << "Error: " << ex.what() << '\n';
            return 1;
        }
    }

    return 0;
}

}  // namespace

/**
 * @brief 程序入口函数
 *
 * 解析命令行参数，根据参数决定运行模式：
 * - 脚本文件执行模式
 * - 非交互式输入模式（管道/重定向）
 * - 交互式REPL模式
 *
 * @param argc 参数个数
 * @param argv 参数数组
 * @return 成功返回 0，失败返回非零值
 */
int main(int argc, char* argv[]) {
    Calculator calculator;
    std::vector<std::string> history;
    bool exact_mode = false;
    bool plain_mode = !isatty(STDIN_FILENO);
    std::string script_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--plain") {
            plain_mode = true;
        } else if (arg == "--version") {
            std::cout << "Calculator Version 2.0.0L\n";
            return 0;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [options] [script.calc]\n"
                      << "Options:\n"
                      << "  --plain      Minimal output (no prompts/headers)\n"
                      << "  --version    Show version info\n"
                      << "  --help, -h   Show this help\n";
            return 0;
        } else if (arg.size() > 5 && arg.substr(arg.size() - 5) == ".calc") {
            script_path = arg;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return 1;
        }
    }

    if (!script_path.empty()) {
        return run_script_file(calculator, script_path, exact_mode);
    }

    if (!isatty(STDIN_FILENO)) {
        const std::string input = read_all(std::cin);
        if (utils::trim_copy(input).empty()) {
            return 0;
        }
        try {
            return run_non_interactive_input(calculator, input);
        } catch (const std::exception& ex) {
            std::cerr << "Error: " << ex.what() << '\n';
            return 1;
        }
    }

    if (!plain_mode) {
        std::cout << "Command Line Calculator\n";
        std::cout << "Enter an expression, or type 'exit' to quit.\n";
        std::cout << "Type 'help' or ':help' to see available commands.\n";
        std::cout << "Use ':exact on' or ':exact off' to toggle exact fraction mode.\n";
        std::cout << "Use ':symbolic on' or ':symbolic off' to preserve pi/e in scalar results.\n";
        std::cout << "Use ':vars', ':clear name', or ':clear' to manage variables.\n";
        std::cout << "Use 'f(x) = ...', 'poly_*', 'diff(...)', 'taylor(...)', 'limit(...)', 'integral(...)', and 'extrema(...)' for custom functions.\n";
        std::cout << "Use ':save file', ':load file', or ':run file.calc' for files and scripts.\n";
    }

    while (true) {
        std::string line;
        if (plain_mode) {
            if (!std::getline(std::cin, line)) break;
        } else {
            line = read_line_with_history(calculator, "> ", history);
            if (!std::cin && line.empty()) break;
        }

        if (line.empty()) continue;

        try {
            bool should_exit = false;
            const std::string output =
                execute_repl_line(calculator, line, &exact_mode, history, &should_exit);
            if (!should_exit && !plain_mode) {
                history.push_back(line);
            }
            if (!output.empty()) {
                std::cout << output << '\n';
            }
            if (should_exit) {
                break;
            }
        } catch (const std::exception& ex) {
            if (plain_mode) {
                std::cerr << "Error: " << ex.what() << '\n';
            } else {
                std::cout << "Error: " << ex.what() << '\n';
            }
        }
    }

    return 0;
}

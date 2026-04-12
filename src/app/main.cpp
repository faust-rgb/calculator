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

#include <cctype>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
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

void redraw_input(const std::string& prompt, const std::string& line) {
    // \33[2K 清空当前行，\r 回到行首，然后重画提示符和输入内容。
    std::cout << "\33[2K\r" << prompt << line << std::flush;
}

std::string trim_copy(const std::string& text) {
    std::size_t start = 0;
    while (start < text.size() &&
           std::isspace(static_cast<unsigned char>(text[start]))) {
        ++start;
    }

    std::size_t end = text.size();
    while (end > start &&
           std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }

    return text.substr(start, end - start);
}

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

const std::vector<std::string>& completion_words() {
    // 这是一个固定的“轻量补全词典”。
    // 目前只覆盖高频命令和函数名，不尝试做上下文感知补全。
    static const std::vector<std::string> words = {
        "help", ":help", ":help commands", ":help functions", ":help matrix", ":help examples",
        ":exact", ":exact on", ":exact off",
        ":symbolic", ":symbolic on", ":symbolic off",
        ":vars", ":funcs", ":history", ":clear", ":clearfunc", ":clearfuncs", ":save", ":load",
        ":run",
        "sin(", "cos(", "tan(", "asin(", "acos(", "atan(",
        "exp(", "ln(", "log10(", "sqrt(", "cbrt(", "root(",
        "abs(", "sign(", "floor(", "ceil(",
        "min(", "max(", "gcd(", "lcm(", "mod(", "pow(", "factor(",
        "poly_add(", "poly_sub(", "poly_mul(", "poly_div(", "roots(",
        "diff(", "limit(", "integral(", "taylor(", "extrema(", "simplify(",
        "vec(", "mat(", "zeros(", "eye(", "identity(",
        "resize(", "append_row(", "append_col(", "transpose(", "inverse(",
        "dot(", "outer(", "null(", "least_squares(", "qr_q(", "qr_r(",
        "svd_u(", "svd_s(", "svd_vt(",
        "solve(", "get(", "set(",
        "norm(", "trace(", "det(", "rank(", "rref(", "eigvals(", "eigvecs(",
        "bin(", "oct(", "hex(", "base(",
        "and(", "or(", "xor(", "not(", "shl(", "shr(",
        "exit", "quit"
    };
    return words;
}

std::pair<std::size_t, std::string> current_token(const std::string& line) {
    // 从行尾向前扫描，找出当前正在输入的 token。
    // 这样 Tab 补全只会替换最后一段命令/函数前缀。
    std::size_t start = line.size();
    while (start > 0) {
        const char ch = line[start - 1];
        if (std::isalnum(static_cast<unsigned char>(ch)) ||
            ch == ':' || ch == '_' || ch == '(') {
            --start;
        } else {
            break;
        }
    }
    return {start, line.substr(start)};
}

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

bool apply_completion(std::string* line) {
    // 补全策略：
    // 1. 找到当前 token
    // 2. 用固定词典做前缀匹配
    // 3. 单候选则直接补全，多候选则补到它们的最长公共前缀
    const auto [start, token] = current_token(*line);
    if (token.empty()) {
        return false;
    }

    std::vector<std::string> matches;
    for (const std::string& word : completion_words()) {
        if (word.rfind(token, 0) == 0) {
            matches.push_back(word);
        }
    }

    if (matches.empty()) {
        return false;
    }

    const std::string replacement =
        matches.size() == 1 ? matches.front() : common_prefix(matches);
    if (replacement.size() <= token.size()) {
        return false;
    }

    line->replace(start, token.size(), replacement);
    return true;
}

std::string read_line_with_history(const std::string& prompt,
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
    std::size_t history_index = history.size();

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

        if (ch == 4) {
            // Ctrl+D: 在空行上视为 EOF，其他情况下忽略。
            if (line.empty()) {
                std::cout << '\n';
            }
            return line;
        }

        if (ch == 127 || ch == '\b') {
            if (!line.empty()) {
                line.pop_back();
                redraw_input(prompt, line);
            }
            continue;
        }

        if (ch == '\t') {
            // Tab 只负责就地补全，不自动执行输入。
            if (apply_completion(&line)) {
                redraw_input(prompt, line);
            }
            continue;
        }

        if (ch == '\33') {
            // 方向键通常会发出 ESC [ A / B / C / D 这样的序列。
            char seq[2] = {'\0', '\0'};
            if (read(STDIN_FILENO, &seq[0], 1) <= 0 ||
                read(STDIN_FILENO, &seq[1], 1) <= 0) {
                continue;
            }

            if (seq[0] == '[' && seq[1] == 'A') {
                // 上方向键：从最近一条历史开始向前浏览。
                if (!history.empty() && history_index > 0) {
                    --history_index;
                    line = history[history_index];
                    redraw_input(prompt, line);
                }
            } else if (seq[0] == '[' && seq[1] == 'B') {
                // 下方向键：向较新的历史移动，超过尾部则回到空输入。
                if (history_index < history.size()) {
                    ++history_index;
                    line = history_index == history.size() ? "" : history[history_index];
                    redraw_input(prompt, line);
                }
            }
            continue;
        }

        if (ch >= 32 && ch <= 126) {
            // 仅接受可打印 ASCII 字符，足够覆盖当前表达式输入需求。
            line.push_back(ch);
            std::cout << ch << std::flush;
        }
    }
}

std::string read_all(std::istream& in) {
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

bool looks_like_script(const std::string& input) {
    return input.find('{') != std::string::npos ||
           input.find("fn ") != std::string::npos ||
           input.find("if ") != std::string::npos ||
           input.find("while ") != std::string::npos ||
           input.find("for ") != std::string::npos;
}

bool line_requires_script_mode(const std::string& line) {
    const std::string trimmed = trim_copy(line);
    if (trimmed.empty()) {
        return false;
    }
    return trimmed.back() == ';' ||
           trimmed.find('{') != std::string::npos ||
           trimmed.find('}') != std::string::npos ||
           trimmed.rfind("fn ", 0) == 0 ||
           trimmed.rfind("if ", 0) == 0 ||
           trimmed.rfind("while ", 0) == 0 ||
           trimmed.rfind("for ", 0) == 0 ||
           trimmed == "break" ||
           trimmed == "continue" ||
           trimmed.rfind("return", 0) == 0 ||
           trimmed.rfind("print(", 0) == 0;
}

bool should_execute_as_script(const std::string& input) {
    if (looks_like_script(input)) {
        return true;
    }

    std::istringstream stream(input);
    std::string line;
    while (std::getline(stream, line)) {
        if (line_requires_script_mode(line)) {
            return true;
        }
    }
    return false;
}

std::string execute_repl_line(Calculator& calculator,
                              const std::string& raw_line,
                              bool* exact_mode,
                              const std::vector<std::string>& history,
                              bool* should_exit) {
    const std::string line = trim_copy(raw_line);
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
        return calculator.help_topic(line.substr(6));
    }
    if (line == ":exact on") {
        *exact_mode = true;
        return "Exact fraction mode: ON";
    }
    if (line == ":exact off") {
        *exact_mode = false;
        return "Exact fraction mode: OFF";
    }
    if (line == ":exact") {
        return std::string("Exact fraction mode: ") + (*exact_mode ? "ON" : "OFF");
    }
    if (line == ":symbolic on") {
        return calculator.set_symbolic_constants_mode(true);
    }
    if (line == ":symbolic off") {
        return calculator.set_symbolic_constants_mode(false);
    }
    if (line == ":symbolic") {
        return std::string("Symbolic constants mode: ") +
               (calculator.symbolic_constants_mode() ? "ON" : "OFF");
    }
    if (line == ":vars") {
        return calculator.list_variables();
    }
    if (line == ":funcs") {
        std::string output;
        return calculator.try_process_function_command(line, &output)
                   ? output
                   : "No custom functions defined.";
    }
    if (line == ":history") {
        return format_history(history);
    }
    if (line.rfind(":save ", 0) == 0) {
        return calculator.save_state(line.substr(6));
    }
    if (line.rfind(":load ", 0) == 0) {
        return calculator.load_state(line.substr(6));
    }
    if (line.rfind(":run ", 0) == 0) {
        std::ifstream in(line.substr(5));
        if (!in) {
            throw std::runtime_error("unable to open script file: " + line.substr(5));
        }
        return calculator.execute_script(read_all(in), *exact_mode);
    }
    if (line == ":clear") {
        return calculator.clear_all_variables();
    }
    if (line == ":clearfuncs") {
        std::string output;
        return calculator.try_process_function_command(line, &output)
                   ? output
                   : "Cleared all custom functions.";
    }
    if (line.rfind(":clearfunc ", 0) == 0) {
        std::string output;
        if (!calculator.try_process_function_command(line, &output)) {
            throw std::runtime_error("invalid custom function command");
        }
        return output;
    }
    if (line.rfind(":clear ", 0) == 0) {
        return calculator.clear_variable(line.substr(7));
    }

    std::string function_output;
    if (calculator.try_process_function_command(line, &function_output)) {
        return function_output;
    }
    if (line.rfind("factor", 0) == 0) {
        return calculator.factor_expression(line);
    }
    return calculator.process_line(line, *exact_mode);
}

int run_non_interactive_input(Calculator& calculator, const std::string& input) {
    if (should_execute_as_script(input)) {
        std::cout << calculator.execute_script(input, false) << '\n';
        return 0;
    }

    std::istringstream stream(input);
    std::string line;
    std::vector<std::string> history;
    bool exact_mode = false;

    while (std::getline(stream, line)) {
        const std::string trimmed = trim_copy(line);
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

int main() {
    // Calculator 封装了解析和求值逻辑，main 只负责与用户交互。
    Calculator calculator;
    std::vector<std::string> history;
    bool exact_mode = false;

    if (!isatty(STDIN_FILENO)) {
        const std::string input = read_all(std::cin);
        if (trim_copy(input).empty()) {
            return 0;
        }
        try {
            return run_non_interactive_input(calculator, input);
        } catch (const std::exception& ex) {
            std::cerr << "Error: " << ex.what() << '\n';
            return 1;
        }
    }

    std::cout << "Command Line Calculator\n";
    std::cout << "Enter an expression, or type 'exit' to quit.\n";
    std::cout << "Type 'help' or ':help' to see available commands.\n";
    std::cout << "Use ':exact on' or ':exact off' to toggle exact fraction mode.\n";
    std::cout << "Use ':symbolic on' or ':symbolic off' to preserve pi/e in scalar results.\n";
    std::cout << "Use ':vars', ':clear name', or ':clear' to manage variables.\n";
    std::cout << "Use 'f(x) = ...', 'poly_*', 'diff(...)', 'taylor(...)', 'limit(...)', 'integral(...)', and 'extrema(...)' for custom functions.\n";
    std::cout << "Use ':save file', ':load file', or ':run file' for files and scripts.\n";

    while (true) {
        // 自定义输入函数支持按上方向键回填历史命令。
        std::string line = read_line_with_history("> ", history);
        if (!std::cin && line.empty()) {
            break;
        }

        if (line.empty()) {
            // 空输入直接跳过，避免给出噪声错误。
            continue;
        }

        try {
            bool should_exit = false;
            const std::string output =
                execute_repl_line(calculator, line, &exact_mode, history, &should_exit);
            if (!should_exit) {
                history.push_back(line);
            }
            if (!output.empty()) {
                std::cout << output << '\n';
            }
            if (should_exit) {
                break;
            }
        } catch (const std::exception& ex) {
            // 解析错误、定义域错误和除零都会在这里统一展示。
            std::cout << "Error: " << ex.what() << '\n';
        }
    }

    return 0;
}

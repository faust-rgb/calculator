// ============================================================================
// repl_application.cpp - 交互式命令行应用程序实现
// ============================================================================

#ifndef APP_REPL_APPLICATION_CPP
#define APP_REPL_APPLICATION_CPP

#include "repl_application.h"
#include "core/services/string_utils.h"
#include <iostream>
#include <termios.h>
#include <unistd.h>
#include <algorithm>

namespace {

/**
 * @class RawModeGuard
 * @brief 终端原始模式管理器（RAII）
 */
class RawModeGuard {
public:
    explicit RawModeGuard(int fd) : fd_(fd) {
        enabled_ = tcgetattr(fd_, &original_) == 0;
        if (!enabled_) return;
        termios raw = original_;
        raw.c_lflag &= static_cast<unsigned long>(~(ICANON | ECHO));
        raw.c_cc[VMIN] = 1;
        raw.c_cc[VTIME] = 0;
        enabled_ = tcsetattr(fd_, TCSAFLUSH, &raw) == 0;
    }
    ~RawModeGuard() {
        if (enabled_) tcsetattr(fd_, TCSAFLUSH, &original_);
    }
    bool enabled() const { return enabled_; }
private:
    int fd_;
    bool enabled_ = false;
    termios original_ {};
};

void redraw_input(const std::string& prompt, const std::string& line, std::size_t cursor_pos) {
    std::cout << "\33[2K\r" << prompt << line;
    std::size_t offset = prompt.size() + cursor_pos;
    if (offset > 0) std::cout << "\r\33[" << offset << "C";
    else std::cout << "\r";
    std::cout << std::flush;
}

} // namespace

ReplApplication::ReplApplication(Calculator& calculator) : calculator_(calculator) {}

int ReplApplication::run(bool plain_mode) {
    if (!plain_mode) {
        std::cout << "Command Line Calculator (Clean Architecture Edition)\n";
        std::cout << "Enter an expression, or type 'exit' to quit.\n";
    }

    while (true) {
        std::string line;
        if (plain_mode) {
            if (!std::getline(std::cin, line)) break;
        } else {
            line = read_line("> ", false);
            if (!std::cin && line.empty()) break;
        }

        if (line.empty()) continue;

        try {
            bool should_exit = false;
            const std::string output = execute_line(line, &should_exit);
            if (!should_exit && !plain_mode) {
                history_.push_back(line);
            }
            if (!output.empty()) {
                std::cout << output << '\n';
            }
            if (should_exit) break;
        } catch (const std::exception& ex) {
            std::cout << "Error: " << ex.what() << '\n';
        }
    }
    return 0;
}

std::string ReplApplication::read_line(const std::string& prompt, bool plain_mode) {
    if (plain_mode || !isatty(STDIN_FILENO)) {
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

    std::cout << prompt << std::flush;
    std::string line;
    std::size_t cursor_pos = 0;
    std::size_t history_index = history_.size();

    while (true) {
        char ch = '\0';
        if (read(STDIN_FILENO, &ch, 1) <= 0) return line;
        if (ch == '\n' || ch == '\r') { std::cout << '\n'; return line; }
        
        if (ch == 127 || ch == '\b') {
            if (cursor_pos > 0) {
                line.erase(cursor_pos - 1, 1);
                --cursor_pos;
                redraw_input(prompt, line, cursor_pos);
            }
            continue;
        }

        if (ch == '\33') {
            char seq[2];
            if (read(STDIN_FILENO, &seq[0], 1) <= 0 || read(STDIN_FILENO, &seq[1], 1) <= 0) continue;
            if (seq[0] == '[') {
                if (seq[1] == 'A') { // Up
                    if (!history_.empty() && history_index > 0) {
                        --history_index;
                        line = history_[history_index];
                        cursor_pos = line.size();
                        redraw_input(prompt, line, cursor_pos);
                    }
                } else if (seq[1] == 'B') { // Down
                    if (history_index < history_.size()) {
                        ++history_index;
                        line = history_index == history_.size() ? "" : history_[history_index];
                        cursor_pos = line.size();
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

std::string ReplApplication::execute_line(const std::string& raw_line, bool* should_exit) {
    const std::string line = utils::trim_copy(raw_line);
    if (line.empty()) return "";
    if (line == "exit" || line == "quit") { *should_exit = true; return ""; }

    // 同步精确模式状态
    if (line == ":exact on") exact_mode_ = true;
    else if (line == ":exact off") exact_mode_ = false;

    return calculator_.process_line(line, exact_mode_);
}

#endif // APP_REPL_APPLICATION_CPP

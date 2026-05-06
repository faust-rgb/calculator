#include "parser/parser_utils.h"
#include "core/string_utils.h"

namespace parser_utils {

// ============================================================================
// 内部辅助：字符串状态追踪器
// ============================================================================

/**
 * @class StringScanState
 * @brief 追踪字符串扫描状态，统一处理转义和引号
 *
 * 使用方法：
 * - 每次迭代调用 update(ch)
 * - 检查 in_string() 判断是否在字符串中
 * - 检查 is_top_level() 判断是否在顶层（不在任何嵌套结构中）
 */
class StringScanState {
public:
    StringScanState() : in_string_(false), escaping_(false), depth_(0) {}

    /// 更新状态
    void update(char ch) {
        if (in_string_) {
            if (escaping_) {
                escaping_ = false;
            } else if (ch == '\\') {
                escaping_ = true;
            } else if (ch == '"') {
                in_string_ = false;
            }
            return;
        }

        if (ch == '"') {
            in_string_ = true;
            return;
        }

        if (ch == '(' || ch == '[' || ch == '{') depth_++;
        else if (ch == ')' || ch == ']' || ch == '}') depth_--;
    }

    /// 是否在字符串中
    bool in_string() const { return in_string_; }

    /// 是否在顶层（不在任何嵌套结构中）
    bool is_top_level() const { return depth_ == 0 && !in_string_; }

    /// 当前嵌套深度
    int depth() const { return depth_; }

    /// 增加深度（用于特定字符检测）
    void inc_depth() { depth_++; }

    /// 减少深度（用于特定字符检测）
    void dec_depth() { depth_--; }

private:
    bool in_string_;
    bool escaping_;
    int depth_;
};

// ============================================================================
// 公共 API 实现
// ============================================================================

BalancedScanResult scan_balanced(std::string_view text) {
    BalancedScanResult res;
    StringScanState state;

    for (std::size_t i = 0; i < text.size(); ++i) {
        char ch = text[i];

        // 先处理字符串状态
        if (state.in_string()) {
            state.update(ch);
            continue;
        }

        state.update(ch);

        // 检查不匹配的关闭括号
        if (ch == ')') {
            if (res.paren_depth == 0) {
                res.balanced = false;
                res.first_mismatch_pos = i;
                return res;
            }
            res.paren_depth--;
        } else if (ch == ']') {
            if (res.bracket_depth == 0) {
                res.balanced = false;
                res.first_mismatch_pos = i;
                return res;
            }
            res.bracket_depth--;
        } else if (ch == '}') {
            if (res.brace_depth == 0) {
                res.balanced = false;
                res.first_mismatch_pos = i;
                return res;
            }
            res.brace_depth--;
        } else if (ch == '(') {
            res.paren_depth++;
        } else if (ch == '[') {
            res.bracket_depth++;
        } else if (ch == '{') {
            res.brace_depth++;
        }
    }

    if (res.paren_depth != 0 || res.bracket_depth != 0 || res.brace_depth != 0 || state.in_string()) {
        res.balanced = false;
    }
    return res;
}

bool is_wrapped_by(std::string_view text, char open, char close) {
    std::string_view trimmed = utils::trim_view(text);
    if (trimmed.size() < 2 || trimmed.front() != open || trimmed.back() != close) return false;

    // Check if the wrapping is actually balanced
    StringScanState state;
    int match_depth = 0;

    for (std::size_t i = 0; i < trimmed.size() - 1; ++i) {
        char ch = trimmed[i];
        state.update(ch);

        if (!state.in_string()) {
            if (ch == open) match_depth++;
            else if (ch == close) {
                if (--match_depth == 0) return false; // Closed too early
            }
        }
    }
    return match_depth == 1;
}

std::size_t find_top_level(std::string_view text, char target) {
    StringScanState state;
    for (std::size_t i = 0; i < text.size(); ++i) {
        char ch = text[i];
        state.update(ch);
        if (state.is_top_level() && ch == target) return i;
    }
    return std::string_view::npos;
}

std::vector<std::string> split_top_level(std::string_view text, char delimiter) {
    std::vector<std::string> result;
    std::size_t last = 0;
    StringScanState state;

    for (std::size_t i = 0; i < text.size(); ++i) {
        char ch = text[i];
        state.update(ch);
        if (state.is_top_level() && ch == delimiter) {
            result.push_back(std::string(text.substr(last, i - last)));
            last = i + 1;
        }
    }
    result.push_back(std::string(text.substr(last)));
    return result;
}

bool contains_script_syntax(std::string_view text) {
    StringScanState state;
    for (char ch : text) {
        state.update(ch);
        if (!state.in_string() && (ch == '{' || ch == '}')) return true;
    }
    return false;
}

} // namespace parser_utils
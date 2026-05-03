#ifndef CALCULATOR_BASE_PARSER_H
#define CALCULATOR_BASE_PARSER_H

#include <string>
#include <string_view>
#include <vector>
#include <cctype>
#include <stdexcept>
#include <sstream>
#include "calculator_exceptions.h"

/**
 * @class BaseParser
 * @brief 基础解析器类，提供通用的词法分析工具
 */
class BaseParser {
public:
    /**
     * @struct NestingTracker
     * @brief 追踪嵌套深度（(), [], {}）
     */
    struct NestingTracker {
        int paren_depth = 0;
        int bracket_depth = 0;
        int brace_depth = 0;

        bool is_balanced() const { return paren_depth == 0 && bracket_depth == 0 && brace_depth == 0; }
        bool is_top_level() const { return is_balanced(); }

        void update(char ch) {
            if (ch == '(') paren_depth++;
            else if (ch == ')') { if (paren_depth > 0) paren_depth--; }
            else if (ch == '[') bracket_depth++;
            else if (ch == ']') { if (bracket_depth > 0) bracket_depth--; }
            else if (ch == '{') brace_depth++;
            else if (ch == '}') { if (brace_depth > 0) brace_depth--; }
        }
    };

protected:
    std::string_view source_;
    std::size_t pos_ = 0;
    bool skip_comments_ = false;  ///< 是否跳过注释

    explicit BaseParser(std::string_view source, bool skip_comments = false)
        : source_(source), skip_comments_(skip_comments) {}

    /// 判断字符是否为空白
    static constexpr bool is_space_char(char ch) {
        return ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r';
    }

    /// 跳过空白字符
    void skip_spaces() {
        while (pos_ < source_.size() && is_space_char(source_[pos_])) {
            ++pos_;
        }
    }

    /// 跳过空白和注释
    void skip_ignorable() {
        while (pos_ < source_.size()) {
            if (is_space_char(source_[pos_])) {
                ++pos_;
                continue;
            }
            if (skip_comments_ && source_[pos_] == '#') {
                while (pos_ < source_.size() && source_[pos_] != '\n') {
                    ++pos_;
                }
                continue;
            }
            if (skip_comments_ && pos_ + 1 < source_.size() &&
                source_[pos_] == '/' && source_[pos_ + 1] == '/') {
                pos_ += 2;
                while (pos_ < source_.size() && source_[pos_] != '\n') {
                    ++pos_;
                }
                continue;
            }
            break;
        }
    }

    bool is_at_end() const {
        return pos_ >= source_.size();
    }

    char peek() const {
        if (is_at_end()) return '\0';
        return source_[pos_];
    }

    char peek_next() const {
        if (pos_ + 1 >= source_.size()) return '\0';
        return source_[pos_ + 1];
    }

    bool match(char expected) {
        if (peek() == expected) {
            ++pos_;
            return true;
        }
        return false;
    }

    bool match_string(std::string_view expected) {
        if (source_.substr(pos_, expected.size()) == expected) {
            pos_ += expected.size();
            return true;
        }
        return false;
    }

    void expect(char expected) {
        if (!match(expected)) {
            throw_error_at_pos(std::string("expected '") + expected + "'");
        }
    }

    /// 解析标识符
    std::string_view parse_identifier() {
        std::size_t start = pos_;
        while (!is_at_end() && (std::isalnum(static_cast<unsigned char>(source_[pos_])) ||
                                source_[pos_] == '_' || source_[pos_] == '\'')) {
            ++pos_;
        }
        return source_.substr(start, pos_ - start);
    }

    /// 解析字符串字面量（处理转义）
    std::string parse_string_literal(std::string* raw_text = nullptr) {
        std::size_t start = pos_;
        if (!match('"')) throw_error_at_pos("expected '\"'");

        std::string content;
        bool escaping = false;
        while (!is_at_end()) {
            const char ch = source_[pos_++];
            if (escaping) {
                switch (ch) {
                    case 'n': content.push_back('\n'); break;
                    case 't': content.push_back('\t'); break;
                    case '\\': content.push_back('\\'); break;
                    case '"': content.push_back('"'); break;
                    default: content.push_back(ch); break;
                }
                escaping = false;
            } else if (ch == '\\') {
                escaping = true;
            } else if (ch == '"') {
                if (raw_text) *raw_text = std::string(source_.substr(start, pos_ - start));
                return content;
            } else {
                content.push_back(ch);
            }
        }
        throw_error_at_pos("unterminated string literal");
    }

    /// 解析通用数字记号（支持进制前缀和虚数后缀）
    std::string_view parse_number_token_view() {
        const std::size_t start = pos_;

        // 进制前缀 (0x, 0b, 0o)
        if (peek() == '0' && pos_ + 1 < source_.size()) {
            char base_char = std::tolower(static_cast<unsigned char>(source_[pos_ + 1]));
            if (base_char == 'x' || base_char == 'b' || base_char == 'o') {
                pos_ += 2;
                while (!is_at_end() && std::isxdigit(static_cast<unsigned char>(peek()))) {
                    pos_++;
                }
                goto end_num;
            }
        }

        // 整数/小数/指数
        if (peek() == '.') pos_++;
        while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) pos_++;
        if (!is_at_end() && peek() == '.') {
            pos_++;
            while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) pos_++;
        }
        if (!is_at_end() && (peek() == 'e' || peek() == 'E')) {
            pos_++;
            if (!is_at_end() && (peek() == '+' || peek() == '-')) pos_++;
            while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) pos_++;
        }
        
        // 虚数后缀
        if (!is_at_end() && peek() == 'i') pos_++;

    end_num:
        return source_.substr(start, pos_ - start);
    }

    template <typename ExceptionType = SyntaxError>
    [[noreturn]] void throw_error_at_pos(const std::string& message, std::size_t error_pos = std::string::npos) const {
        if (error_pos == std::string::npos) {
            error_pos = pos_ < source_.size() ? pos_ : source_.size();
        }
        std::ostringstream oss;
        oss << message << " at position " << error_pos << "\n";
        oss << "  " << source_ << "\n";
        oss << "  " << std::string(error_pos, ' ') << "^";
        throw ExceptionType(oss.str());
    }

    bool peek_is_digit() const {
        return std::isdigit(static_cast<unsigned char>(peek()));
    }

    bool peek_is_alpha() const {
        return std::isalpha(static_cast<unsigned char>(peek()));
    }

    bool peek_is_identifier_start() const {
        return !is_at_end() &&
               (std::isalpha(static_cast<unsigned char>(source_[pos_])) ||
                source_[pos_] == '_');
    }

};

#endif

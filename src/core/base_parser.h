#ifndef CALCULATOR_BASE_PARSER_H
#define CALCULATOR_BASE_PARSER_H

#include <string>
#include <string_view>
#include <vector>
#include <cctype>
#include <stdexcept>
#include "calculator_exceptions.h"
#include "calculator_internal_types.h"

/**
 * @class BaseParser
 * @brief 基础解析器类，提供通用的词法分析工具
 *
 * 所有解析器（DecimalParser、ScriptParser、MatrixExpressionParser 等）
 * 都应继承此类，以确保词法分析的一致性。
 */
class BaseParser {
protected:
    std::string_view source_;
    std::size_t pos_ = 0;
    bool skip_comments_ = false;  ///< 是否跳过注释

    explicit BaseParser(std::string_view source, bool skip_comments = false)
        : source_(source), skip_comments_(skip_comments) {}

    /// 判断字符是否为空白（内联比较，避免 std::isspace 的 locale 开销）
    static constexpr bool is_space_char(char ch) {
        return ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r';
    }

    /// 跳过空白字符
    void skip_spaces() {
        while (pos_ < source_.size() && is_space_char(source_[pos_])) {
            ++pos_;
        }
    }

    /// 跳过空白和注释（用于脚本解析）
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

    /// 检查下一个字符是否匹配（不消费）
    bool peek(char expected) const {
        return peek() == expected;
    }

    /// 尝试匹配并消费指定字符
    bool match(char expected) {
        if (peek() == expected) {
            ++pos_;
            return true;
        }
        return false;
    }

    /// 尝试匹配并消费指定字符串
    bool match_string(const std::string& expected) {
        if (source_.compare(pos_, expected.size(), expected) == 0) {
            pos_ += expected.size();
            return true;
        }
        return false;
    }

    /// 匹配并消费指定字符，失败则抛出异常
    void expect(char expected) {
        if (!match(expected)) {
            std::string msg = "expected '";
            msg += expected;
            msg += "' but found '";
            msg += (is_at_end() ? "EOF" : std::string(1, peek()));
            msg += "'";
            throw SyntaxError(msg);
        }
    }

    /// 检查当前字符是否为字母或下划线
    bool peek_is_alpha() const {
        return std::isalpha(static_cast<unsigned char>(peek())) || peek() == '_';
    }

    /// 检查当前字符是否为数字
    bool peek_is_digit() const {
        return std::isdigit(static_cast<unsigned char>(peek()));
    }

    /// 检查当前字符是否为标识符起始字符
    bool peek_is_identifier_start() const {
        return !is_at_end() &&
               (std::isalpha(static_cast<unsigned char>(source_[pos_])) ||
                source_[pos_] == '_');
    }

    /// 解析标识符（字母或下划线开头，后跟字母、数字或下划线）
    std::string_view parse_identifier() {
        std::size_t start = pos_;
        while (!is_at_end() && (std::isalnum(static_cast<unsigned char>(source_[pos_])) ||
                                source_[pos_] == '_' || source_[pos_] == '\'')) {
            ++pos_;
        }
        return source_.substr(start, pos_ - start);
    }

    /// 尝试匹配关键字（确保后面不是标识符字符）
    bool match_keyword(const std::string& keyword) {
        if (skip_comments_) skip_ignorable();
        if (source_.compare(pos_, keyword.size(), keyword) != 0) {
            return false;
        }
        const std::size_t after = pos_ + keyword.size();
        if (after < source_.size()) {
            const char next = source_[after];
            if (std::isalnum(static_cast<unsigned char>(next)) || next == '_') {
                return false;
            }
        }
        pos_ = after;
        return true;
    }

    /// 跳过空白后尝试匹配字符
    bool match_with_skip(char expected) {
        skip_spaces();
        return match(expected);
    }

    /// 跳过空白后尝试匹配字符（支持注释跳过）
    bool match_ignorable(char expected) {
        if (skip_comments_) skip_ignorable();
        else skip_spaces();
        return match(expected);
    }

    /// 跳过空白后期望字符
    void expect_with_skip(char expected) {
        skip_spaces();
        expect(expected);
    }

    /// 跳过空白后期望字符（支持注释跳过）
    void expect_ignorable(char expected) {
        if (skip_comments_) skip_ignorable();
        else skip_spaces();
        expect(expected);
    }
};

#endif

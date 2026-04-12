/**
 * @file script_parser.cpp
 * @brief 脚本语言解析器实现
 *
 * 使用递归下降算法解析自定义脚本语言。
 * 支持：代码块、函数定义、条件语句、循环、控制流等。
 */

#include "script_parser.h"

#include <cctype>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace script {
namespace {

/**
 * @brief 去除字符串首尾空白字符
 * @param text 输入字符串
 * @return 去除空白后的字符串
 */
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

/**
 * @class Parser
 * @brief 递归下降解析器
 *
 * 将源代码解析为抽象语法树（AST）。
 * 处理词法分析和语法分析两个阶段。
 */
class Parser {
public:
    explicit Parser(std::string source) : source_(std::move(source)) {}

    Program parse() {
        Program program;
        skip_ignorable();
        while (!is_at_end()) {
            program.statements.push_back(parse_statement());
            skip_ignorable();
        }
        return program;
    }

private:
    StatementPtr parse_statement() {
        skip_ignorable();
        if (match('{')) {
            return parse_block_after_open();
        }
        if (match_keyword("fn")) {
            return parse_function();
        }
        if (match_keyword("if")) {
            return parse_if();
        }
        if (match_keyword("while")) {
            return parse_while();
        }
        if (match_keyword("for")) {
            return parse_for();
        }
        if (match_keyword("return")) {
            return parse_return();
        }
        if (match_keyword("break")) {
            auto statement = std::make_unique<BreakStatement>();
            expect(';');
            return statement;
        }
        if (match_keyword("continue")) {
            auto statement = std::make_unique<ContinueStatement>();
            expect(';');
            return statement;
        }
        return parse_simple_statement();
    }

    StatementPtr parse_block_after_open() {
        auto block = std::make_unique<BlockStatement>();
        skip_ignorable();
        while (!is_at_end() && !peek('}')) {
            block->statements.push_back(parse_statement());
            skip_ignorable();
        }
        expect('}');
        return block;
    }

    StatementPtr parse_function() {
        skip_ignorable();
        const std::string name = parse_identifier();
        skip_ignorable();
        const std::string parameter_text = parse_group('(', ')');
        skip_ignorable();
        expect('{');

        auto statement = std::make_unique<FunctionStatement>();
        statement->name = name;
        statement->parameters = split_top_level(parameter_text, ',');

        for (std::string& parameter : statement->parameters) {
            parameter = trim_copy(parameter);
            if (parameter.empty()) {
                throw std::runtime_error("function parameters must be identifiers");
            }
        }

        StatementPtr body = parse_block_after_open();
        statement->body.reset(static_cast<BlockStatement*>(body.release()));
        return statement;
    }

    StatementPtr parse_if() {
        auto statement = std::make_unique<IfStatement>();
        skip_ignorable();
        statement->condition = trim_copy(parse_group('(', ')'));
        skip_ignorable();
        statement->then_branch = parse_statement();
        skip_ignorable();
        if (match_keyword("else")) {
            skip_ignorable();
            statement->else_branch = parse_statement();
        }
        return statement;
    }

    StatementPtr parse_while() {
        auto statement = std::make_unique<WhileStatement>();
        skip_ignorable();
        statement->condition = trim_copy(parse_group('(', ')'));
        skip_ignorable();
        statement->body = parse_statement();
        return statement;
    }

    StatementPtr parse_for() {
        auto statement = std::make_unique<ForStatement>();
        skip_ignorable();
        const std::string header = parse_group('(', ')');
        const std::vector<std::string> parts = split_top_level(header, ';');
        if (parts.size() != 3) {
            throw std::runtime_error("for expects initializer; condition; step");
        }
        statement->initializer = trim_copy(parts[0]);
        statement->condition = trim_copy(parts[1]);
        statement->step = trim_copy(parts[2]);
        skip_ignorable();
        statement->body = parse_statement();
        return statement;
    }

    StatementPtr parse_return() {
        auto statement = std::make_unique<ReturnStatement>();
        const std::string expression = trim_copy(parse_until_semicolon());
        if (!expression.empty()) {
            statement->has_expression = true;
            statement->expression = expression;
        }
        return statement;
    }

    StatementPtr parse_simple_statement() {
        auto statement = std::make_unique<SimpleStatement>();
        statement->text = trim_copy(parse_until_semicolon());
        if (statement->text.empty()) {
            throw std::runtime_error("empty statement");
        }
        return statement;
    }

    std::string parse_group(char open, char close) {
        expect(open);
        const std::size_t start = pos_;
        int depth = 1;
        bool in_string = false;
        bool escaping = false;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (in_string) {
                if (escaping) {
                    escaping = false;
                } else if (ch == '\\') {
                    escaping = true;
                } else if (ch == '"') {
                    in_string = false;
                }
            } else if (ch == '"') {
                in_string = true;
            } else if (ch == open) {
                ++depth;
            } else if (ch == close) {
                --depth;
                if (depth == 0) {
                    const std::string content = source_.substr(start, pos_ - start);
                    ++pos_;
                    return content;
                }
            }
            ++pos_;
        }
        throw std::runtime_error("unterminated grouped expression");
    }

    std::string parse_until_semicolon() {
        const std::size_t start = pos_;
        int paren_depth = 0;
        int bracket_depth = 0;
        bool in_string = false;
        bool escaping = false;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (in_string) {
                if (escaping) {
                    escaping = false;
                } else if (ch == '\\') {
                    escaping = true;
                } else if (ch == '"') {
                    in_string = false;
                }
            } else if (ch == '"') {
                in_string = true;
            } else if (ch == '(') {
                ++paren_depth;
            } else if (ch == ')') {
                --paren_depth;
            } else if (ch == '[') {
                ++bracket_depth;
            } else if (ch == ']') {
                --bracket_depth;
            } else if (ch == ';' && paren_depth == 0 && bracket_depth == 0) {
                const std::string text = source_.substr(start, pos_ - start);
                ++pos_;
                return text;
            }
            ++pos_;
        }
        throw std::runtime_error("expected ';'");
    }

    std::string parse_identifier() {
        if (is_at_end() ||
            !std::isalpha(static_cast<unsigned char>(source_[pos_]))) {
            throw std::runtime_error("expected identifier");
        }

        const std::size_t start = pos_;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
                ++pos_;
            } else {
                break;
            }
        }
        return source_.substr(start, pos_ - start);
    }

    std::vector<std::string> split_top_level(const std::string& text, char delimiter) const {
        std::vector<std::string> parts;
        int paren_depth = 0;
        int bracket_depth = 0;
        bool in_string = false;
        bool escaping = false;
        std::size_t start = 0;

        for (std::size_t i = 0; i < text.size(); ++i) {
            const char ch = text[i];
            if (in_string) {
                if (escaping) {
                    escaping = false;
                } else if (ch == '\\') {
                    escaping = true;
                } else if (ch == '"') {
                    in_string = false;
                }
                continue;
            }
            if (ch == '"') {
                in_string = true;
            } else if (ch == '(') {
                ++paren_depth;
            } else if (ch == ')') {
                --paren_depth;
            } else if (ch == '[') {
                ++bracket_depth;
            } else if (ch == ']') {
                --bracket_depth;
            } else if (ch == delimiter && paren_depth == 0 && bracket_depth == 0) {
                parts.push_back(text.substr(start, i - start));
                start = i + 1;
            }
        }

        if (!text.empty() || delimiter == ',') {
            parts.push_back(text.substr(start));
        }
        if (text.empty() && delimiter == ',') {
            parts.clear();
        }
        return parts;
    }

    void skip_ignorable() {
        while (!is_at_end()) {
            if (std::isspace(static_cast<unsigned char>(source_[pos_]))) {
                ++pos_;
                continue;
            }
            if (source_[pos_] == '#') {
                while (!is_at_end() && source_[pos_] != '\n') {
                    ++pos_;
                }
                continue;
            }
            if (source_[pos_] == '/' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '/') {
                pos_ += 2;
                while (!is_at_end() && source_[pos_] != '\n') {
                    ++pos_;
                }
                continue;
            }
            break;
        }
    }

    bool match_keyword(const std::string& keyword) {
        skip_ignorable();
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

    bool match(char ch) {
        skip_ignorable();
        if (is_at_end() || source_[pos_] != ch) {
            return false;
        }
        ++pos_;
        return true;
    }

    bool peek(char ch) const {
        return !is_at_end() && source_[pos_] == ch;
    }

    void expect(char ch) {
        skip_ignorable();
        if (is_at_end() || source_[pos_] != ch) {
            throw std::runtime_error(std::string("expected '") + ch + "'");
        }
        ++pos_;
    }

    bool is_at_end() const {
        return pos_ >= source_.size();
    }

    std::string source_;
    std::size_t pos_ = 0;
};

}  // namespace

Program parse_program(const std::string& source) {
    Parser parser(source);
    return parser.parse();
}

}  // namespace script

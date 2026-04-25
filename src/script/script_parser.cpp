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
#include <sstream>
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

std::vector<std::string> split_top_level_text(const std::string& text, char delimiter) {
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

std::string strip_hash_comment(const std::string& line) {
    bool in_string = false;
    bool escaping = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char ch = line[i];
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
        } else if (ch == '#') {
            return line.substr(0, i);
        }
    }
    return line;
}

int indentation_width(const std::string& line) {
    int width = 0;
    for (const char ch : line) {
        if (ch == ' ') {
            ++width;
        } else if (ch == '\t') {
            width += 4;
        } else {
            break;
        }
    }
    return width;
}

bool starts_with_keyword(const std::string& text, const std::string& keyword) {
    if (text.compare(0, keyword.size(), keyword) != 0) {
        return false;
    }
    if (text.size() == keyword.size()) {
        return true;
    }
    const char next = text[keyword.size()];
    return !std::isalnum(static_cast<unsigned char>(next)) && next != '_';
}

bool has_legacy_block_markers(const std::string& source) {
    bool in_string = false;
    bool escaping = false;
    for (const char ch : source) {
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
        } else if (ch == '{' || ch == '}') {
            return true;
        }
    }
    return false;
}

std::string parenthesize_condition(const std::string& condition) {
    return "(" + trim_copy(condition) + ")";
}

std::string translate_range_for(const std::string& header) {
    const std::size_t in_pos = header.find(" in ");
    if (in_pos == std::string::npos) {
        throw std::runtime_error("for expects 'name in range(...)' or a legacy for header");
    }

    const std::string variable = trim_copy(header.substr(0, in_pos));
    const std::string range_call = trim_copy(header.substr(in_pos + 4));
    if (variable.empty() || range_call.rfind("range(", 0) != 0 || range_call.back() != ')') {
        throw std::runtime_error("for expects 'name in range(...)'");
    }

    std::vector<std::string> parts =
        split_top_level_text(range_call.substr(6, range_call.size() - 7), ',');
    for (std::string& part : parts) {
        part = trim_copy(part);
    }
    if (parts.empty() || parts.size() > 3) {
        throw std::runtime_error("range expects stop, start/stop, or start/stop/step");
    }

    std::string start = "0";
    std::string stop = parts[0];
    std::string step = "1";
    if (parts.size() >= 2) {
        start = parts[0];
        stop = parts[1];
    }
    if (parts.size() == 3) {
        step = parts[2];
    }
    if (stop.empty() || start.empty() || step.empty()) {
        throw std::runtime_error("range arguments cannot be empty");
    }

    const std::string comparison = step.front() == '-' ? " > " : " < ";
    return variable + " = " + start + "; " + variable + comparison + stop + "; " +
           variable + " = " + variable + " + " + step;
}

std::string normalize_python_like_script(const std::string& source) {
    if (has_legacy_block_markers(source)) {
        return source;
    }

    std::istringstream input(source);
    std::ostringstream output;
    std::vector<int> indent_stack = {0};
    bool pending_block = false;
    std::string line;

    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        const std::string without_comment = strip_hash_comment(line);
        const std::string stripped = trim_copy(without_comment);
        if (stripped.empty()) {
            continue;
        }

        const int indent = indentation_width(without_comment);
        if (indent > indent_stack.back()) {
            if (!pending_block) {
                throw std::runtime_error("unexpected indentation");
            }
            indent_stack.push_back(indent);
            pending_block = false;
        } else {
            if (pending_block) {
                throw std::runtime_error("expected an indented block");
            }
            while (indent < indent_stack.back()) {
                output << "}\n";
                indent_stack.pop_back();
            }
            if (indent != indent_stack.back()) {
                throw std::runtime_error("inconsistent indentation");
            }
        }

        if (stripped == "pass") {
            continue;
        }

        if (!stripped.empty() && stripped.back() == ':') {
            const std::string header = trim_copy(stripped.substr(0, stripped.size() - 1));
            if (starts_with_keyword(header, "def")) {
                output << "fn " << trim_copy(header.substr(3)) << " {\n";
            } else if (starts_with_keyword(header, "fn")) {
                output << "fn " << trim_copy(header.substr(2)) << " {\n";
            } else if (starts_with_keyword(header, "if")) {
                output << "if " << parenthesize_condition(header.substr(2)) << " {\n";
            } else if (starts_with_keyword(header, "elif")) {
                output << "else if " << parenthesize_condition(header.substr(4)) << " {\n";
            } else if (header == "else") {
                output << "else {\n";
            } else if (starts_with_keyword(header, "while")) {
                output << "while " << parenthesize_condition(header.substr(5)) << " {\n";
            } else if (starts_with_keyword(header, "for")) {
                const std::string for_header = trim_copy(header.substr(3));
                if (!for_header.empty() && for_header.front() == '(' && for_header.back() == ')') {
                    output << "for " << for_header << " {\n";
                } else {
                    output << "for (" << translate_range_for(for_header) << ") {\n";
                }
            } else {
                throw std::runtime_error("unknown block header: " + header);
            }
            pending_block = true;
            continue;
        }

        if (starts_with_keyword(stripped, "return")) {
            output << "return " << trim_copy(stripped.substr(6)) << ";\n";
        } else if (stripped == "break" || stripped == "continue") {
            output << stripped << ";\n";
        } else {
            output << stripped;
            if (stripped.back() != ';') {
                output << ';';
            }
            output << '\n';
        }
    }

    if (pending_block) {
        throw std::runtime_error("expected an indented block");
    }
    while (indent_stack.size() > 1) {
        output << "}\n";
        indent_stack.pop_back();
    }

    return output.str();
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
    explicit Parser(std::string source)
        : source_(normalize_python_like_script(std::move(source))) {}

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
        return split_top_level_text(text, delimiter);
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

/**
 * @file script_parser.cpp
 * @brief 脚本语言解析器实现
 *
 * 使用递归下降算法解析自定义脚本语言。
 * 支持：代码块、函数定义、条件语句、循环、控制流等。
 */

#include "script_parser.h"
#include "parser/base_parser.h"
#include "command/expression_compiler.h"
#include "parser/lazy_token_stream.h"
#include "string_utils.h"

#include <cctype>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace script {
namespace {

using utils::trim_copy;

/**
 * @brief 创建表达式缓存（在解析阶段预编译）
 * @param expression 表达式字符串
 * @return 预编译的表达式缓存
 */
std::shared_ptr<ExpressionCache> create_expression_cache(const std::string& expression) {
    auto cache = std::make_shared<ExpressionCache>();
    cache->expanded = expression;  // 展开会在运行时进行
    cache->hint = analyze_expression_hint(expression);
    cache->features = analyze_expression_features(expression);
    return cache;
}

/**
 * @brief 基于 Token 流的参数拆分（替代手动字符扫描）
 * @param text 输入文本
 * @param delimiter 分隔符
 * @return 拆分后的参数列表
 *
 * 使用 LazyTokenStream 进行正确的词法分析，
 * 自动处理字符串、注释和嵌套结构。
 */
std::vector<std::string> split_top_level_by_tokens(const std::string& text, char delimiter) {
    std::vector<std::string> parts;
    LazyTokenStream tokens(text);

    std::size_t start_pos = 0;
    int paren_depth = 0;
    int bracket_depth = 0;
    int brace_depth = 0;

    while (!tokens.is_at_end()) {
        ::Token tok = tokens.advance();

        // 跟踪嵌套深度
        switch (tok.kind) {
            case TokenKind::kLParen:
                paren_depth++;
                break;
            case TokenKind::kRParen:
                if (paren_depth > 0) paren_depth--;
                break;
            case TokenKind::kLBracket:
                bracket_depth++;
                break;
            case TokenKind::kRBracket:
                if (bracket_depth > 0) bracket_depth--;
                break;
            case TokenKind::kLBrace:
                brace_depth++;
                break;
            case TokenKind::kRBrace:
                if (brace_depth > 0) brace_depth--;
                break;
            case TokenKind::kComma:
                // 只在顶层拆分
                if (paren_depth == 0 && bracket_depth == 0 && brace_depth == 0 && delimiter == ',') {
                    std::string part = text.substr(start_pos, tok.position - start_pos);
                    parts.push_back(trim_copy(part));
                    start_pos = tok.position + 1;
                }
                break;
            case TokenKind::kSemicolon:
                // 只在顶层拆分
                if (paren_depth == 0 && bracket_depth == 0 && brace_depth == 0 && delimiter == ';') {
                    std::string part = text.substr(start_pos, tok.position - start_pos);
                    parts.push_back(trim_copy(part));
                    start_pos = tok.position + 1;
                }
                break;
            default:
                break;
        }
    }

    // 添加最后一个部分
    if (start_pos < text.size()) {
        std::string part = text.substr(start_pos);
        parts.push_back(trim_copy(part));
    } else if (!parts.empty() && delimiter == ',') {
        // 尾部逗号情况
        parts.push_back("");
    }

    return parts;
}

/**
 * @struct Token
 * @brief 脚本标记
 */
struct Token {
    enum class Kind {
        kIdentifier,
        kKeyword,
        kSymbol,
        kNumber,
        kString,
        kNewline,
        kIndent,
        kDedent,
        kEof,
    };

    Kind kind;
    std::string_view text;
    int line = 1;
};

/**
 * @class ScriptLexer
 * @brief 脚本词法分析器，支持缩进处理和行号追踪
 */
class ScriptLexer : public BaseParser {
public:
    ScriptLexer(std::string_view source)
        : BaseParser(source), current_line_{1}, indent_stack_{0} {}


    std::vector<Token> tokenize() {
        std::vector<Token> tokens;
        while (!is_at_end()) {
            // 处理行首缩进
            if (pos_ == 0 || (pos_ > 0 && source_[pos_ - 1] == '\n')) {
                handle_indentation(tokens);
            }

            skip_spaces_except_newline();
            if (is_at_end()) break;

            const char ch = peek();
            if (ch == '\n') {
                tokens.push_back({Token::Kind::kNewline, "\n", current_line_});
                ++pos_;
                ++current_line_;
                continue;
            }

            if (ch == '#' || (ch == '/' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '/')) {
                skip_comment();
                continue;
            }

            if (ch == '"') {
                tokens.push_back(parse_string_literal_token());
                continue;
            }

            if (std::isdigit(static_cast<unsigned char>(ch)) || 
                (ch == '.' && pos_ + 1 < source_.size() && std::isdigit(static_cast<unsigned char>(source_[pos_ + 1])))) {
                tokens.push_back(parse_number_literal_token());
                continue;
            }

            if (peek_is_identifier_start()) {
                const std::string_view id = parse_identifier();
                if (is_keyword(id)) {
                    tokens.push_back({Token::Kind::kKeyword, id, current_line_});
                } else {
                    tokens.push_back({Token::Kind::kIdentifier, id, current_line_});
                }
            } else if (std::ispunct(static_cast<unsigned char>(ch))) {
                // 处理多字符符号
                int start_line = current_line_;
                const std::size_t start = pos_;
                if (ch == '=' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '=') {
                    pos_ += 2;
                    tokens.push_back({Token::Kind::kSymbol, source_.substr(start, 2), start_line});
                } else if (ch == '!' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '=') {
                    pos_ += 2;
                    tokens.push_back({Token::Kind::kSymbol, source_.substr(start, 2), start_line});
                } else if (ch == '<' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '=') {
                    pos_ += 2;
                    tokens.push_back({Token::Kind::kSymbol, source_.substr(start, 2), start_line});
                } else if (ch == '>' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '=') {
                    pos_ += 2;
                    tokens.push_back({Token::Kind::kSymbol, source_.substr(start, 2), start_line});
                } else {
                    ++pos_;
                    tokens.push_back({Token::Kind::kSymbol, source_.substr(start, 1), start_line});
                }
            } else {
                const std::size_t start = pos_;
                ++pos_;
                tokens.push_back({Token::Kind::kSymbol, source_.substr(start, 1), current_line_});
            }
        }

        // 补齐末尾的 Dedent
        while (indent_stack_.size() > 1) {
            tokens.push_back({Token::Kind::kDedent, "", current_line_});
            indent_stack_.pop_back();
        }
        tokens.push_back({Token::Kind::kEof, "", current_line_});
        return tokens;
    }

private:
    void skip_spaces_except_newline() {
        while (pos_ < source_.size() && std::isspace(static_cast<unsigned char>(source_[pos_])) && source_[pos_] != '\n') {
            ++pos_;
        }
    }

    void skip_comment() {
        while (pos_ < source_.size() && source_[pos_] != '\n') {
            ++pos_;
        }
    }

    Token parse_string_literal_token() {
        const std::size_t start = pos_;
        int start_line = current_line_;
        expect('"');
        bool escaping = false;
        while (!is_at_end()) {
            const char ch = source_[pos_++];
            if (ch == '\n') ++current_line_;
            if (escaping) {
                escaping = false;
            } else if (ch == '\\') {
                escaping = true;
            } else if (ch == '"') {
                return {Token::Kind::kString, source_.substr(start, pos_ - start), start_line};
            }
        }
        throw std::runtime_error("Line " + std::to_string(start_line) + ": Unterminated string literal");
    }

    Token parse_number_literal_token() {
        const std::size_t start = pos_;
        int start_line = current_line_;
        bool seen_dot = false;
        if (peek() == '.') {
            seen_dot = true;
            ++pos_;
        }
        while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) {
            ++pos_;
        }
        if (!seen_dot && !is_at_end() && peek() == '.') {
            seen_dot = true;
            ++pos_;
            while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) {
                ++pos_;
            }
        }
        // 指数部分
        if (!is_at_end() && (peek() == 'e' || peek() == 'E')) {
            ++pos_;
            if (!is_at_end() && (peek() == '+' || peek() == '-')) ++pos_;
            while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) {
                ++pos_;
            }
        }
        return {Token::Kind::kNumber, source_.substr(start, pos_ - start), start_line};
    }

    void handle_indentation(std::vector<Token>& tokens) {
        int current_indent = 0;
        while (pos_ < source_.size() && (source_[pos_] == ' ' || source_[pos_] == '\t')) {
            if (source_[pos_] == ' ') current_indent += 1;
            else current_indent += 4;
            ++pos_;
        }

        if (is_at_end() || peek() == '\n' || peek() == '#' || (peek() == '/' && pos_ + 1 < source_.size() && source_[pos_ + 1] == '/')) {
            return; // 忽略空行或只有注释的行
        }

        if (current_indent > indent_stack_.back()) {
            indent_stack_.push_back(current_indent);
            tokens.push_back({Token::Kind::kIndent, "", current_line_});
        } else {
            while (current_indent < indent_stack_.back()) {
                tokens.push_back({Token::Kind::kDedent, "", current_line_});
                indent_stack_.pop_back();
            }
            if (current_indent != indent_stack_.back()) {
                throw std::runtime_error("Line " + std::to_string(current_line_) + ": Inconsistent indentation");
            }
        }
    }

    bool is_keyword(std::string_view id) {
        static const std::vector<std::string_view> keywords = {
            "fn", "def", "if", "elif", "else", "while", "for", "return", "break", "continue", "pass", "in"
        };
        for (const auto& k : keywords) if (k == id) return true;
        return false;
    }

    int current_line_;
    std::vector<int> indent_stack_;
};

/**
 * @brief 辅助函数：拆分顶层文本（使用 Token 流实现）
 *
 * 保留旧函数签名以兼容，内部使用新的 Token 流实现。
 */
std::vector<std::string> split_top_level_text(const std::string& text, char delimiter) {
    return split_top_level_by_tokens(text, delimiter);
}

/**
 * @brief 辅助函数：转换 Python 风格的 for 循环为 C 风格（使用纯 Token 流）
 */
std::string translate_range_for(const std::string& header) {
    // 使用 Token 流查找 "in" 关键字
    LazyTokenStream tokens(header);
    std::string variable;
    std::string range_args_text;
    bool found_in = false;

    // 第一遍：找到变量名和 "in" 关键字
    std::string last_identifier;
    while (!tokens.is_at_end()) {
        ::Token tok = tokens.advance();

        if (tok.kind == TokenKind::kIdentifier) {
            if (tok.text == "in" && !last_identifier.empty()) {
                // 找到了 "in"，last_identifier 是循环变量
                variable = last_identifier;
                found_in = true;
                break;
            }
            last_identifier = std::string(tok.text);
        } else if (tok.kind != TokenKind::kEnd) {
            // 重置（非标识符 token）
            if (tok.kind != TokenKind::kLParen && tok.kind != TokenKind::kRParen) {
                // 保留 last_identifier 以支持 "for i in range(...)" 格式
            }
        }
    }

    if (!found_in) {
        return header; // 可能是已经转换过的或者是 C 风格的 (init; cond; step)
    }

    // 收集 "in" 之后的所有 token 直到末尾
    std::vector<::Token> range_tokens;
    while (!tokens.is_at_end()) {
        ::Token tok = tokens.advance();
        if (tok.kind == TokenKind::kEnd) break;
        range_tokens.push_back(tok);
    }

    // 验证格式：应该是 "range(...)"
    if (range_tokens.size() < 2 ||
        range_tokens[0].kind != TokenKind::kIdentifier ||
        range_tokens[0].text != "range" ||
        range_tokens[1].kind != TokenKind::kLParen) {
        throw std::runtime_error("for expects 'name in range(...)'");
    }

    // 找到右括号
    if (range_tokens.back().kind != TokenKind::kRParen) {
        throw std::runtime_error("for expects 'name in range(...)'");
    }

    // 提取 range 参数（跳过 "range(" 和最后的 ")"）
    std::vector<std::string> parts;
    std::string current_arg;
    int paren_depth = 0;

    for (std::size_t i = 2; i < range_tokens.size() - 1; ++i) {
        const ::Token& tok = range_tokens[i];

        if (tok.kind == TokenKind::kLParen) {
            paren_depth++;
            current_arg += "(";
        } else if (tok.kind == TokenKind::kRParen) {
            paren_depth--;
            current_arg += ")";
        } else if (tok.kind == TokenKind::kComma && paren_depth == 0) {
            parts.push_back(trim_copy(current_arg));
            current_arg.clear();
        } else {
            current_arg += std::string(tok.text);
        }
    }

    // 添加最后一个参数
    if (!current_arg.empty() || !parts.empty()) {
        parts.push_back(trim_copy(current_arg));
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

    const std::string comparison = (!step.empty() && step[0] == '-') ? " > " : " < ";
    return variable + " = " + start + "; " + variable + comparison + stop + "; " +
           variable + " = " + variable + " + " + step;
}

std::string strip_enclosing_parentheses(const std::string& text) {
    const std::string trimmed = trim_copy(text);
    if (trimmed.size() < 2 || trimmed.front() != '(' || trimmed.back() != ')') {
        return trimmed;
    }

    int depth = 0;
    bool in_string = false;
    bool escaping = false;
    for (std::size_t i = 0; i < trimmed.size(); ++i) {
        const char ch = trimmed[i];
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
            ++depth;
        } else if (ch == ')') {
            --depth;
            if (depth == 0 && i != trimmed.size() - 1) {
                return trimmed;
            }
        }
    }

    return depth == 0 ? trim_copy(trimmed.substr(1, trimmed.size() - 2)) : trimmed;
}

/**
 * @class ScriptParserImpl
 * @brief 递归下降解析器
 */
class ScriptParserImpl {
public:
    explicit ScriptParserImpl(std::vector<Token> tokens)
        : tokens_(std::move(tokens)), pos_(0), recursion_depth_(0) {}

    Program parse() {
        Program program;
        skip_newlines();
        while (!is_at_end()) {
            program.statements.push_back(parse_statement());
            skip_newlines();
        }
        return program;
    }

private:
    std::vector<Token> tokens_;
    std::size_t pos_;
    int recursion_depth_;
    static constexpr int kMaxRecursionDepth = 256;

    void check_recursion_depth() {
        if (recursion_depth_ > kMaxRecursionDepth) {
            throw std::runtime_error("Line " + std::to_string(peek().line) + ": Script too deeply nested (recursion limit exceeded)");
        }
    }

    bool is_at_end() const {
        return pos_ >= tokens_.size() || tokens_[pos_].kind == Token::Kind::kEof;
    }

    const Token& peek() const {
        return tokens_[pos_];
    }

    Token advance() {
        if (!is_at_end()) ++pos_;
        return tokens_[pos_ - 1];
    }

    bool match(Token::Kind kind) {
        if (peek().kind == kind) { advance(); return true; }
        return false;
    }

    bool match_keyword(const std::string& text) {
        if (peek().kind == Token::Kind::kKeyword && peek().text == text) {
            advance();
            return true;
        }
        return false;
    }

    bool match_symbol(const std::string& text) {
        if (peek().kind == Token::Kind::kSymbol && peek().text == text) {
            advance();
            return true;
        }
        return false;
    }

    void expect(Token::Kind kind, const std::string& msg = "Unexpected token") {
        if (!match(kind)) {
            throw std::runtime_error("Line " + std::to_string(peek().line) + ": " + msg);
        }
    }

    void expect_symbol(const std::string& text) {
        if (!match_symbol(text)) {
            throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected symbol '" + text + "'");
        }
    }

    void skip_newlines() {
        while (!is_at_end() && peek().kind == Token::Kind::kNewline) {
            advance();
        }
    }

    void skip_brace_block_separators() {
        while (!is_at_end() &&
               (peek().kind == Token::Kind::kNewline ||
                peek().kind == Token::Kind::kIndent ||
                peek().kind == Token::Kind::kDedent)) {
            advance();
        }
    }

    StatementPtr parse_statement() {
        skip_newlines();
        recursion_depth_++;
        check_recursion_depth();
        
        StatementPtr result;
        if (match_symbol("{")) {
            result = parse_block_after_open();
        } else if (match_keyword("fn") || match_keyword("def")) {
            result = parse_function();
        } else if (match_keyword("if")) {
            result = parse_if();
        } else if (match_keyword("while")) {
            result = parse_while();
        } else if (match_keyword("for")) {
            result = parse_for();
        } else if (match_keyword("return")) {
            result = parse_return();
        } else if (match_keyword("break")) {
            auto statement = std::make_unique<BreakStatement>();
            if (!is_at_end() && peek().kind == Token::Kind::kSymbol && peek().text == ";") advance();
            result = std::move(statement);
        } else if (match_keyword("continue")) {
            auto statement = std::make_unique<ContinueStatement>();
            if (!is_at_end() && peek().kind == Token::Kind::kSymbol && peek().text == ";") advance();
            result = std::move(statement);
        } else if (match_keyword("pass")) {
            advance();
            result = std::make_unique<SimpleStatement>(); 
        } else {
            result = parse_simple_statement();
        }
        
        recursion_depth_--;
        return result;
    }

    StatementPtr parse_block_after_open() {
        auto block = std::make_unique<BlockStatement>();
        skip_brace_block_separators();
        while (!is_at_end() && (peek().kind != Token::Kind::kSymbol || peek().text != "}")) {
            block->statements.push_back(parse_statement());
            skip_brace_block_separators();
        }
        expect_symbol("}");
        return block;
    }

    StatementPtr parse_indented_block() {
        if (match_symbol("{")) return parse_block_after_open();
        
        expect(Token::Kind::kIndent, "Expected indented block");
        auto block = std::make_unique<BlockStatement>();
        while (!is_at_end() && peek().kind != Token::Kind::kDedent) {
            block->statements.push_back(parse_statement());
            skip_newlines();
        }
        expect(Token::Kind::kDedent, "Expected dedent at end of block");
        return block;
    }

    StatementPtr parse_function() {
        const int line = peek().line;
        if (is_at_end() || peek().kind != Token::Kind::kIdentifier) {
            throw std::runtime_error("Line " + std::to_string(line) + ": Expected function name");
        }
        const std::string name = std::string(advance().text);
        expect_symbol("(");
        std::vector<std::string> params;
        if (!match_symbol(")")) {
            while (true) {
                if (peek().kind != Token::Kind::kIdentifier) {
                    throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected parameter name");
                }
                params.push_back(std::string(advance().text));
                if (!match_symbol(",")) break;
            }
            expect_symbol(")");
        }
        
        if (match_symbol(":")) {
            skip_newlines();
            return finish_function(name, params, parse_indented_block());
        } else {
            expect_symbol("{");
            return finish_function(name, params, parse_block_after_open());
        }
    }

    StatementPtr finish_function(const std::string& name, const std::vector<std::string>& params, StatementPtr body) {
        auto statement = std::make_unique<FunctionStatement>();
        statement->name = name;
        statement->parameters = params;
        statement->body.reset(static_cast<BlockStatement*>(body.release()));
        return statement;
    }

    StatementPtr parse_if() {
        auto statement = std::make_unique<IfStatement>();
        statement->condition = parse_expression_string();
        statement->cache = create_expression_cache(statement->condition);

        if (match_symbol(":")) {
            skip_newlines();
            statement->then_branch = parse_indented_block();
        } else {
            statement->then_branch = parse_statement();
        }

        skip_newlines();
        if (match_keyword("else")) {
            if (match_symbol(":")) {
                skip_newlines();
                statement->else_branch = parse_indented_block();
            } else {
                statement->else_branch = parse_statement();
            }
        } else if (match_keyword("elif")) {
            statement->else_branch = parse_if();
        }
        return statement;
    }

    StatementPtr parse_while() {
        auto statement = std::make_unique<WhileStatement>();
        statement->condition = parse_expression_string();
        statement->cache = create_expression_cache(statement->condition);
        if (match_symbol(":")) {
            skip_newlines();
            statement->body = parse_indented_block();
        } else {
            statement->body = parse_statement();
        }
        return statement;
    }

    StatementPtr parse_for() {
        std::string header = strip_enclosing_parentheses(parse_expression_string());

        // 检测是否为 Python 风格的 for 循环
        const std::size_t in_pos = header.find(" in ");
        if (in_pos != std::string::npos) {
            // Python 风格：for i in range(...)
            return parse_for_range(header, in_pos);
        }

        // C 风格：for (init; cond; step)
        return parse_for_c_style(header);
    }

    StatementPtr parse_for_range(const std::string& header, std::size_t in_pos) {
        auto statement = std::make_unique<ForRangeStatement>();

        const std::string variable = trim_copy(header.substr(0, in_pos));
        const std::string range_call = trim_copy(header.substr(in_pos + 4));

        if (variable.empty()) {
            throw std::runtime_error("for expects variable name before 'in'");
        }
        if (range_call.rfind("range(", 0) != 0 || range_call.back() != ')') {
            throw std::runtime_error("for expects 'name in range(...)'");
        }

        statement->variable = variable;

        // 解析 range 参数
        std::vector<std::string> parts =
            split_top_level_text(range_call.substr(6, range_call.size() - 7), ',');
        for (std::string& part : parts) {
            part = trim_copy(part);
        }
        if (parts.empty() || parts.size() > 3) {
            throw std::runtime_error("range expects stop, start/stop, or start/stop/step");
        }

        // 设置 start, stop, step
        statement->start_expr = "0";
        statement->stop_expr = parts[0];
        statement->step_expr = "1";

        if (parts.size() >= 2) {
            statement->start_expr = parts[0];
            statement->stop_expr = parts[1];
        }
        if (parts.size() == 3) {
            statement->step_expr = parts[2];
        }

        // 检测步长是否为负数
        if (!statement->step_expr.empty() && statement->step_expr[0] == '-') {
            statement->step_is_negative = true;
        }

        // 创建缓存
        statement->start_cache = create_expression_cache(statement->start_expr);
        statement->stop_cache = create_expression_cache(statement->stop_expr);
        statement->step_cache = create_expression_cache(statement->step_expr);

        // 解析循环体
        if (match_symbol(":")) {
            skip_newlines();
            statement->body = parse_indented_block();
        } else {
            statement->body = parse_statement();
        }
        return statement;
    }

    StatementPtr parse_for_c_style(const std::string& header) {
        auto statement = std::make_unique<ForStatement>();
        std::string translated = translate_range_for(header);

        std::vector<std::string> parts = split_top_level_text(translated, ';');
        if (parts.size() == 3) {
            statement->initializer = trim_copy(parts[0]);
            statement->condition = trim_copy(parts[1]);
            statement->step = trim_copy(parts[2]);
            if (!statement->initializer.empty()) statement->init_cache = create_expression_cache(statement->initializer);
            if (!statement->condition.empty()) statement->cond_cache = create_expression_cache(statement->condition);
            if (!statement->step.empty()) statement->step_cache = create_expression_cache(statement->step);
        } else {
            statement->condition = translated;
            if (!statement->condition.empty()) statement->cond_cache = create_expression_cache(statement->condition);
        }

        if (match_symbol(":")) {
            skip_newlines();
            statement->body = parse_indented_block();
        } else {
            statement->body = parse_statement();
        }
        return statement;
    }

    StatementPtr parse_return() {
        auto statement = std::make_unique<ReturnStatement>();
        std::string expr = parse_expression_string();
        if (!expr.empty()) {
            statement->has_expression = true;
            statement->expression = expr;
            statement->cache = create_expression_cache(expr);
        }
        return statement;
    }

    StatementPtr parse_simple_statement() {
        auto statement = std::make_unique<SimpleStatement>();
        const int start_line = peek().line;
        statement->text = parse_expression_string();
        if (statement->text.empty()) {
             // 可能是空行或分号
             return statement;
        }
        try {
            statement->cache = create_expression_cache(statement->text);
        } catch (const std::exception& e) {
            throw std::runtime_error("Line " + std::to_string(start_line) + ": Invalid expression: " + e.what());
        }
        return statement;
    }

    std::string parse_expression_string() {
        int depth = 0;
        int start_line = peek().line;
        if (is_at_end()) return "";

        const char* start_ptr = peek().text.data();
        const char* end_ptr = start_ptr;

        while (!is_at_end()) {
            const auto& t = peek();
            if (depth == 0 && (t.kind == Token::Kind::kNewline || 
                               (t.kind == Token::Kind::kSymbol &&
                                (t.text == ":" || t.text == ";" ||
                                 t.text == "{" || t.text == "}")))) {
                break;
            }
            
            if (t.kind == Token::Kind::kSymbol) {
                if (t.text == "(" || t.text == "[" || t.text == "{") ++depth;
                else if (t.text == ")" || t.text == "]" || t.text == "}") --depth;
            }
            
            end_ptr = t.text.data() + t.text.size();
            advance();
        }

        if (depth != 0) {
            throw std::runtime_error("Line " + std::to_string(start_line) + ": Unbalanced parentheses/brackets");
        }

        std::string result;
        if (start_ptr && end_ptr >= start_ptr) {
            result = std::string(start_ptr, end_ptr - start_ptr);
        }

        if (!is_at_end() && peek().kind == Token::Kind::kSymbol && peek().text == ";") advance();
        return trim_copy(result);
    }
};

}  // namespace

Program parse_program(const std::string& source) {
    ScriptLexer lexer(source);
    ScriptParserImpl parser(lexer.tokenize());
    return parser.parse();
}

}  // namespace script

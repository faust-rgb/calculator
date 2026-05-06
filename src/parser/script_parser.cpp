/**
 * @file script_parser.cpp
 * @brief 脚本语言解析器实现
 */

#include "script_parser.h"
#include "parser/base_parser.h"
#include "parser/expression_compiler.h"
#include "parser/lazy_token_stream.h"
#include "parser/parser_utils.h"
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
 * @brief 检查表达式是否为索引赋值（如 x[i] = value 或 m[r, c] = value）
 */
bool looks_like_index_assignment(std::string_view text) {
    std::size_t equal_pos = parser_utils::find_top_level(text, '=');
    if (equal_pos == std::string_view::npos) return false;
    
    // Check for compound assignments or comparison operators
    if (equal_pos + 1 < text.size() && text[equal_pos + 1] == '=') return false; // ==
    if (equal_pos > 0 && (text[equal_pos - 1] == '!' || text[equal_pos - 1] == '<' || text[equal_pos - 1] == '>')) return false; // !=, <=, >=

    const std::string lhs = trim_copy(std::string(text.substr(0, equal_pos)));
    if (lhs.empty() || lhs.back() != ']' || lhs.find('[') == std::string::npos) return false;

    const std::size_t first_bracket = lhs.find('[');
    if (first_bracket == 0) return false;
    const char first_char = lhs[0];
    if (!std::isalpha(static_cast<unsigned char>(first_char)) && first_char != '_') return false;

    return true;
}

struct Token {
    enum class Kind {
        kIdentifier, kKeyword, kSymbol, kNumber, kString, kNewline, kIndent, kDedent, kEof,
    };
    Kind kind;
    std::string_view text;
    int line = 1;
};

class ScriptLexer {
public:
    explicit ScriptLexer(std::string_view source)
        : source_(source),
          tokens_(source, LazyTokenStream::WhitespaceMode::kScript),
          current_line_{1},
          line_start_pos_{0},
          at_line_start_{true},
          indent_stack_{0} {}

    std::vector<Token> tokenize() {
        std::vector<Token> out;
        while (tokens_.peek().kind != TokenKind::kEnd) {
            if (at_line_start_) {
                handle_indentation(out);
            }

            const ::Token raw = tokens_.advance();
            if (raw.kind == TokenKind::kEnd) {
                break;
            }
            if (raw.kind == TokenKind::kNewline) {
                out.push_back({Token::Kind::kNewline, raw.text, current_line_});
                ++current_line_;
                line_start_pos_ = raw.position + raw.text.size();
                at_line_start_ = true;
                continue;
            }

            out.push_back(convert_token(raw));
            at_line_start_ = false;
        }
        while (indent_stack_.size() > 1) { out.push_back({Token::Kind::kDedent, "", current_line_}); indent_stack_.pop_back(); }
        out.push_back({Token::Kind::kEof, "", current_line_});
        return out;
    }
private:
    Token convert_token(const ::Token& raw) const {
        switch (raw.kind) {
            case TokenKind::kIdentifier:
                if (raw.text == "_") return {Token::Kind::kSymbol, raw.text, current_line_};
                return {is_keyword(raw.text) ? Token::Kind::kKeyword : Token::Kind::kIdentifier,
                        raw.text,
                        current_line_};
            case TokenKind::kNumber:
                return {Token::Kind::kNumber, raw.text, current_line_};
            case TokenKind::kString:
                return {Token::Kind::kString, raw.text, current_line_};
            case TokenKind::kNewline:
                return {Token::Kind::kNewline, raw.text, current_line_};
            case TokenKind::kEnd:
                return {Token::Kind::kEof, raw.text, current_line_};
            default:
                return {Token::Kind::kSymbol, raw.text, current_line_};
        }
    }

    void handle_indentation(std::vector<Token>& out) {
        const ::Token& next = tokens_.peek();
        if (next.kind == TokenKind::kEnd || next.kind == TokenKind::kNewline) {
            return;
        }

        int cur = 0;
        const std::size_t indent_end = next.position;
        for (std::size_t i = line_start_pos_; i < indent_end && i < source_.size(); ++i) {
            if (source_[i] == ' ') cur += 1;
            else if (source_[i] == '\t') cur += 4;
        }

        if (cur > indent_stack_.back()) {
            indent_stack_.push_back(cur);
            out.push_back({Token::Kind::kIndent, "", current_line_});
        } else {
            while (cur < indent_stack_.back()) {
                out.push_back({Token::Kind::kDedent, "", current_line_});
                indent_stack_.pop_back();
            }
            if (cur != indent_stack_.back()) {
                throw std::runtime_error("Line " + std::to_string(current_line_) + ": Indent error");
            }
        }
    }

    static bool is_keyword(std::string_view id) {
        static const std::vector<std::string_view> kw = { "fn", "def", "if", "elif", "else", "while", "for", "return", "break", "continue", "pass", "import", "in", "match", "case" };
        for (const auto& k : kw) {
            if (k == id) return true;
        }
        return false;
    }

    std::string_view source_;
    LazyTokenStream tokens_;
    int current_line_;
    std::size_t line_start_pos_;
    bool at_line_start_;
    std::vector<int> indent_stack_;
};

class ScriptParserImpl {
public:
    explicit ScriptParserImpl(std::vector<Token> tokens,
                              CommandParser::IsCommandCallback is_command)
        : tokens_(std::move(tokens)),
          pos_(0),
          depth_(0),
          is_command_(std::move(is_command)) {}
    Program parse() { Program p; skip_newlines(); while (!is_at_end()) { p.statements.push_back(parse_statement()); skip_newlines(); } return p; }
private:
    std::vector<Token> tokens_; std::size_t pos_; int depth_; CommandParser::IsCommandCallback is_command_;
    bool is_at_end() const { return pos_ >= tokens_.size() || tokens_[pos_].kind == Token::Kind::kEof; }
    const Token& peek() const { return tokens_[pos_]; }
    Token advance() { if (!is_at_end()) ++pos_; return tokens_[pos_ - 1]; }
    bool match_keyword(const std::string& t) { if (peek().kind == Token::Kind::kKeyword && peek().text == t) { advance(); return true; } return false; }
    bool match_symbol(const std::string& t) { if (peek().kind == Token::Kind::kSymbol && peek().text == t) { advance(); return true; } return false; }
    void expect_symbol(const std::string& t) { if (!match_symbol(t)) throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected '" + t + "'"); }
    void skip_newlines() { while (!is_at_end() && peek().kind == Token::Kind::kNewline) advance(); }
    void skip_block_sep() { while (!is_at_end() && (peek().kind == Token::Kind::kNewline || peek().kind == Token::Kind::kIndent || peek().kind == Token::Kind::kDedent)) advance(); }
    CommandASTNode parse_script_command(const std::string& text) { return parse_command(text, is_command_); }

    StatementPtr parse_statement() {
        skip_newlines(); depth_++; if (depth_ > 256) throw std::runtime_error("Line " + std::to_string(peek().line) + ": Nested too deep");
        StatementPtr res;
        if (match_symbol("{")) res = parse_block_after_open();
        else if (match_keyword("fn") || match_keyword("def")) res = parse_function();
        else if (match_keyword("if")) res = parse_if();
        else if (match_keyword("while")) res = parse_while();
        else if (match_keyword("for")) res = parse_for();
        else if (match_keyword("match")) res = parse_match();
        else if (match_keyword("return")) res = parse_return();
        else if (match_keyword("break")) { auto s = std::make_unique<BreakStatement>(); s->line = peek().line; match_symbol(";"); res = std::move(s); }
        else if (match_keyword("continue")) { auto s = std::make_unique<ContinueStatement>(); s->line = peek().line; match_symbol(";"); res = std::move(s); }
        else if (match_keyword("import")) res = parse_import();
        else if (match_keyword("pass")) { auto s = std::make_unique<SimpleStatement>(); s->line = peek().line; res = std::move(s); }
        else res = parse_simple_statement();
        depth_--; return res;
    }

    StatementPtr parse_block_after_open() {
        auto b = std::make_unique<BlockStatement>(); skip_block_sep();
        while (!is_at_end() && (peek().kind != Token::Kind::kSymbol || peek().text != "}")) { b->statements.push_back(parse_statement()); skip_block_sep(); }
        expect_symbol("}"); return b;
    }

    StatementPtr parse_indented_block() {
        if (match_symbol("{")) return parse_block_after_open();
        if (peek().kind != Token::Kind::kIndent) throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected indent");
        advance(); auto b = std::make_unique<BlockStatement>();
        while (!is_at_end() && peek().kind != Token::Kind::kDedent) { b->statements.push_back(parse_statement()); skip_newlines(); }
        if (peek().kind != Token::Kind::kDedent) throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected dedent");
        advance(); return b;
    }

    StatementPtr parse_simple_statement() {
        auto s = std::make_unique<SimpleStatement>(); s->line = peek().line;
        s->text = parse_expr_str();
        if (!s->text.empty() && !looks_like_index_assignment(s->text)) {
            s->command_ast = parse_script_command(s->text);
        }
        return s;
    }

    StatementPtr parse_function() {
        int l = peek().line; if (peek().kind != Token::Kind::kIdentifier) throw std::runtime_error("Line " + std::to_string(l) + ": Expected name");
        std::string n = std::string(advance().text); expect_symbol("("); std::vector<std::string> p;
        if (!match_symbol(")")) {
            while (true) {
                if (peek().kind != Token::Kind::kIdentifier) throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected param");
                p.push_back(std::string(advance().text)); if (!match_symbol(",")) break;
            }
            expect_symbol(")");
        }
        std::shared_ptr<const BlockStatement> body;
        if (match_symbol(":")) { skip_newlines(); body = std::shared_ptr<const BlockStatement>(static_cast<BlockStatement*>(parse_indented_block().release())); }
        else { expect_symbol("{"); body = std::shared_ptr<const BlockStatement>(static_cast<BlockStatement*>(parse_block_after_open().release())); }
        auto s = std::make_unique<FunctionStatement>(); s->line = l; s->name = n; s->parameters = p; s->body = std::move(body);
        return s;
    }

    StatementPtr parse_if() {
        auto s = std::make_unique<IfStatement>(); s->line = peek().line; s->condition_ast = parse_script_command(parse_expr_str());
        if (match_symbol(":")) { skip_newlines(); s->then_branch = parse_indented_block(); } else s->then_branch = parse_statement();
        skip_newlines();
        if (match_keyword("else")) { if (match_symbol(":")) { skip_newlines(); s->else_branch = parse_indented_block(); } else s->else_branch = parse_statement(); }
        else if (match_keyword("elif")) s->else_branch = parse_if();
        return s;
    }

    StatementPtr parse_match() {
        int l = peek().line;
        auto s = std::make_unique<MatchStatement>();
        s->line = l;

        // 解析匹配主体表达式
        std::string subject_expr = parse_expr_str();
        s->subject_ast = parse_script_command(subject_expr);

        // 期望冒号或缩进块开始
        if (!match_symbol(":")) {
            throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected ':' after match expression");
        }
        skip_newlines();

        // 解析 case 分支列表
        // case 分支必须缩进
        if (peek().kind != Token::Kind::kIndent) {
            throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected indent after 'match:'");
        }
        advance();  // 消费 kIndent

        while (!is_at_end() && peek().kind != Token::Kind::kDedent) {
            skip_newlines();

            // 检查是否是 case 关键字
            if (!match_keyword("case")) {
                // 可能是嵌套的语句或到达 dedent
                if (peek().kind == Token::Kind::kDedent) break;
                throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected 'case' in match block");
            }

            CaseClause clause;

            // 解析模式
            // 特殊情况：case _ 表示默认分支
            if (peek().kind == Token::Kind::kSymbol && peek().text == "_") {
                advance();  // 消费 _
                clause.is_default = true;
            } else {
                // 解析模式表达式
                std::string pattern_expr = parse_case_pattern_expr_str();
                if (pattern_expr.empty()) {
                    throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected case pattern");
                }
                clause.pattern_ast = parse_script_command(pattern_expr);
            }

            // 检查是否有守卫条件 (if condition)
            if (match_keyword("if")) {
                clause.is_guarded = true;
                std::string guard_expr = parse_expr_str();
                if (guard_expr.empty()) {
                    throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected guard expression");
                }
                clause.guard_ast = parse_script_command(guard_expr);
            }

            // 期望冒号
            if (!match_symbol(":")) {
                throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected ':' after case pattern");
            }
            skip_newlines();

            // 解析 case 体
            // case 体可以是单行语句或缩进块
            if (peek().kind == Token::Kind::kIndent) {
                clause.body = parse_indented_block();
            } else {
                clause.body = parse_statement();
            }

            s->cases.push_back(std::move(clause));
            skip_newlines();
        }

        if (peek().kind != Token::Kind::kDedent) {
            throw std::runtime_error("Line " + std::to_string(peek().line) + ": Expected dedent after match cases");
        }
        advance();  // 消费 kDedent

        return s;
    }

    StatementPtr parse_while() {
        auto s = std::make_unique<WhileStatement>(); s->line = peek().line; s->condition_ast = parse_script_command(parse_expr_str());
        if (match_symbol(":")) { skip_newlines(); s->body = parse_indented_block(); } else s->body = parse_statement();
        return s;
    }

    struct ExprFragmentOptions {
        bool stop_at_case_guard = false;
        bool stop_at_comma = false;
        bool stop_at_rparen = false;
        bool consume_semicolon = true;
    };

    StatementPtr parse_for() {
        int l = peek().line;

        // 检查是否是 for-in-range 循环：for variable in range(...)
        // 使用 Token 流识别 'in' 关键字，而不是字符串查找
        std::size_t saved_pos = pos_;
        std::string variable;

        // 尝试解析 for-in-range 格式
        if (peek().kind == Token::Kind::kIdentifier) {
            variable = std::string(advance().text);

            // 检查 'in' 关键字
            if (match_keyword("in")) {
                // 检查是否是 range(...) 调用
                if (peek().kind == Token::Kind::kIdentifier && peek().text == "range") {
                    advance();  // 消费 'range'
                    expect_symbol("(");

                    std::vector<std::string> range_args;
                    if (!match_symbol(")")) {
                        while (true) {
                            range_args.push_back(parse_expr_fragment({ false, true, true, false }));
                            if (match_symbol(",")) continue;
                            if (match_symbol(")")) break;
                            throw std::runtime_error("Line " + std::to_string(l) + ": Unmatched '(' in range()");
                        }
                    }

                    auto s = std::make_unique<ForRangeStatement>();
                    s->line = l;
                    s->variable = variable;

                    // 根据 range 参数数量设置 start/stop/step
                    if (range_args.size() == 1) {
                        s->start_ast = parse_script_command("0");
                        s->stop_ast = parse_script_command(range_args[0]);
                        s->step_ast = parse_script_command("1");
                    } else if (range_args.size() == 2) {
                        s->start_ast = parse_script_command(range_args[0]);
                        s->stop_ast = parse_script_command(range_args[1]);
                        s->step_ast = parse_script_command("1");
                    } else if (range_args.size() == 3) {
                        s->start_ast = parse_script_command(range_args[0]);
                        s->stop_ast = parse_script_command(range_args[1]);
                        s->step_ast = parse_script_command(range_args[2]);
                    } else {
                        throw std::runtime_error("Line " + std::to_string(l) + ": range() takes 1-3 arguments");
                    }

                    if (match_symbol(":")) { skip_newlines(); s->body = parse_indented_block(); }
                    else s->body = parse_statement();
                    return s;
                } else {
                    // general for-in loop
                    auto s = std::make_unique<ForInStatement>();
                    s->line = l;
                    s->variable = variable;
                    s->iterable_ast = parse_script_command(parse_expr_str());
                    if (match_symbol(":")) { skip_newlines(); s->body = parse_indented_block(); }
                    else s->body = parse_statement();
                    return s;
                }
            }
        }

        // 不是 for-in-range，回退并解析 C 风格 for 循环
        pos_ = saved_pos;
        std::string h = strip_paren(parse_expr_str());

        auto s = std::make_unique<ForStatement>();
        s->line = l;
        std::vector<std::string> pts;
        std::size_t st = 0;
        int d = 0;
        for (std::size_t i = 0; i < h.size(); ++i) {
            if (h[i] == '(') d++;
            else if (h[i] == ')') d--;
            else if (h[i] == ';' && d == 0) {
                pts.push_back(trim_copy(h.substr(st, i - st)));
                st = i + 1;
            }
        }
        pts.push_back(trim_copy(h.substr(st)));

        if (pts.size() == 3) {
            s->init_ast = parse_script_command(pts[0]);
            s->cond_ast = parse_script_command(pts[1]);
            s->step_ast = parse_script_command(pts[2]);
        } else {
            s->cond_ast = parse_script_command(h);
        }

        if (match_symbol(":")) { skip_newlines(); s->body = parse_indented_block(); }
        else s->body = parse_statement();
        return s;
    }

    StatementPtr parse_return() {
        auto s = std::make_unique<ReturnStatement>(); s->line = peek().line; std::string e = parse_expr_str();
        if (!e.empty()) { s->has_expression = true; s->expr_ast = parse_script_command(e); }
        return s;
    }

    StatementPtr parse_import() {
        auto s = std::make_unique<ImportStatement>();
        s->line = peek().line;
        s->path_text = parse_expr_str();
        if (s->path_text.empty()) {
            throw std::runtime_error("Line " + std::to_string(s->line) + ": Expected import path");
        }
        s->path_ast = parse_script_command(s->path_text);
        return s;
    }

    std::string strip_paren(const std::string& t) {
        std::string s = trim_copy(t); if (s.size() >= 2 && s.front() == '(' && s.back() == ')') {
            int d = 0; for (std::size_t i = 0; i < s.size(); ++i) { if (s[i] == '(') d++; else if (s[i] == ')') { d--; if (d == 0 && i != s.size() - 1) return s; } }
            return trim_copy(s.substr(1, s.size() - 2));
        }
        return s;
    }

    std::string parse_expr_str() {
        return parse_expr_fragment({ false, false, false, true });
    }

    std::string parse_case_pattern_expr_str() {
        return parse_expr_fragment({ true, false, false, false });
    }

    std::string parse_expr_fragment(const ExprFragmentOptions& options) {
        int d = 0; if (is_at_end()) return ""; const char* s_ptr = peek().text.data(); const char* e_ptr = s_ptr;
        std::string_view previous_token;
        while (!is_at_end()) {
            const auto& t = peek();
            if (d == 0) {
                if (options.stop_at_case_guard && t.kind == Token::Kind::kKeyword && t.text == "if") break;
                const bool brace_starts_literal =
                    t.kind == Token::Kind::kSymbol && t.text == "{" &&
                    (previous_token.empty() || previous_token == "=" || previous_token == "(" ||
                     previous_token == "[" || previous_token == "," || previous_token == ":");
                if (t.kind == Token::Kind::kNewline ||
                    (t.kind == Token::Kind::kSymbol &&
                     (t.text == ":" || t.text == ";" || (t.text == "{" && !brace_starts_literal) || t.text == "}" ||
                      (options.stop_at_comma && t.text == ",") ||
                      (options.stop_at_rparen && t.text == ")")))) {
                    break;
                }
            }
            if (t.kind == Token::Kind::kSymbol) {
                if (t.text == "(" || t.text == "[" || t.text == "{") ++d;
                else if (t.text == ")" || t.text == "]" || t.text == "}") --d;
            }
            previous_token = t.text;
            e_ptr = t.text.data() + t.text.size(); advance();
        }
        std::string res; if (s_ptr && e_ptr >= s_ptr) res = std::string(s_ptr, e_ptr - s_ptr);
        if (options.consume_semicolon && !is_at_end() && peek().kind == Token::Kind::kSymbol && peek().text == ";") advance();
        return trim_copy(res);
    }
};

} // namespace

Program parse_program(const std::string& source, CommandParser::IsCommandCallback is_command) {
    ScriptLexer lexer(source);
    ScriptParserImpl parser(lexer.tokenize(), std::move(is_command));
    return parser.parse();
}

} // namespace script

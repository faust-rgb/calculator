// ============================================================================
// 命令解析器实现
// ============================================================================

#include "command_parser.h"
#include "utils.h"
#include <algorithm>
#include <cctype>

// ============================================================================
// CommandASTNode 实现
// ============================================================================

CommandASTNode CommandASTNode::make_empty() {
    CommandASTNode node;
    node.kind = CommandKind::kEmpty;
    node.data = std::monostate{};
    return node;
}

CommandASTNode CommandASTNode::make_meta_command(std::string_view cmd,
                                                  const std::vector<std::string_view>& args) {
    CommandASTNode node;
    node.kind = CommandKind::kMetaCommand;
    node.data = MetaCommandInfo{cmd, args};
    return node;
}

CommandASTNode CommandASTNode::make_function_definition(std::string_view name,
                                                        std::string_view param,
                                                        std::string_view body) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionDefinition;
    node.data = FunctionDefinitionInfo{name, param, body};
    return node;
}

CommandASTNode CommandASTNode::make_function_call(std::string_view name,
                                                  const std::vector<std::string_view>& args) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionCall;
    node.data = FunctionCallInfo{name, args};
    return node;
}

CommandASTNode CommandASTNode::make_assignment(std::string_view var,
                                               std::string_view expr) {
    CommandASTNode node;
    node.kind = CommandKind::kAssignment;
    node.data = AssignmentInfo{var, expr};
    return node;
}

CommandASTNode CommandASTNode::make_expression(std::string_view expr) {
    CommandASTNode node;
    node.kind = CommandKind::kExpression;
    node.data = expr;
    return node;
}

CommandASTNode CommandASTNode::make_string_literal(const std::string& value) {
    CommandASTNode node;
    node.kind = CommandKind::kStringLiteral;
    node.data = value;
    return node;
}

// ============================================================================
// CommandParser 实现 - 词法分析
// ============================================================================

CommandParser::CommandParser(std::string_view source)
    : source_(source) {}

void CommandParser::skip_spaces() {
    while (pos_ < source_.size() && is_space_char(source_[pos_])) {
        ++pos_;
    }
}

char CommandParser::peek() const {
    if (pos_ >= source_.size()) return '\0';
    return source_[pos_];
}

char CommandParser::peek_next() const {
    if (pos_ + 1 >= source_.size()) return '\0';
    return source_[pos_ + 1];
}

char CommandParser::advance() {
    return source_[pos_++];
}

bool CommandParser::is_at_end() const {
    return pos_ >= source_.size();
}

CommandToken CommandParser::next_token() {
    skip_spaces();
    if (is_at_end()) {
        CommandToken tok;
        tok.type = CommandTokenType::kEnd;
        tok.position = pos_;
        return tok;
    }

    const char ch = peek();
    const std::size_t start_pos = pos_;

    // 字符串字面量
    if (ch == '"') {
        return parse_string_token();
    }

    // 数字字面量
    if (std::isdigit(static_cast<unsigned char>(ch)) ||
        (ch == '.' && std::isdigit(static_cast<unsigned char>(peek_next())))) {
        return parse_number_token();
    }

    // 标识符
    if (std::isalpha(static_cast<unsigned char>(ch)) || ch == '_') {
        return parse_identifier_token();
    }

    // 单字符 Token 和运算符
    switch (ch) {
        case '(': {
            advance();
            CommandToken tok;
            tok.type = CommandTokenType::kLParen;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ')': {
            advance();
            CommandToken tok;
            tok.type = CommandTokenType::kRParen;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '[': {
            advance();
            CommandToken tok;
            tok.type = CommandTokenType::kLBracket;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ']': {
            advance();
            CommandToken tok;
            tok.type = CommandTokenType::kRBracket;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ',': {
            advance();
            CommandToken tok;
            tok.type = CommandTokenType::kComma;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ':': {
            advance();
            CommandToken tok;
            tok.type = CommandTokenType::kColon;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '=': {
            // 检查是否是 ==
            if (peek_next() == '=') {
                return parse_operator_token();
            }
            advance();
            CommandToken tok;
            tok.type = CommandTokenType::kEqual;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '+': case '-': case '*': case '/': case '^': case '%':
        case '<': case '>': case '!': case '&': case '|':
            return parse_operator_token();
        default:
            // 未知字符，作为运算符处理
            return parse_operator_token();
    }
}

CommandToken CommandParser::parse_identifier_token() {
    const std::size_t start = pos_;
    while (!is_at_end() && (std::isalnum(static_cast<unsigned char>(peek())) ||
                            peek() == '_' || peek() == '\'')) {
        advance();
    }
    CommandToken tok;
    tok.type = CommandTokenType::kIdentifier;
    tok.text = source_.substr(start, pos_ - start);
    tok.position = start;
    return tok;
}

CommandToken CommandParser::parse_number_token() {
    const std::size_t start = pos_;

    // 整数部分
    while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) {
        advance();
    }

    // 小数部分
    if (!is_at_end() && peek() == '.') {
        advance();
        while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) {
            advance();
        }
    }

    // 指数部分
    if (!is_at_end() && (peek() == 'e' || peek() == 'E')) {
        advance();
        if (!is_at_end() && (peek() == '+' || peek() == '-')) {
            advance();
        }
        while (!is_at_end() && std::isdigit(static_cast<unsigned char>(peek()))) {
            advance();
        }
    }

    CommandToken tok;
    tok.type = CommandTokenType::kNumber;
    tok.text = source_.substr(start, pos_ - start);
    tok.position = start;

    // 解析数值
    try {
        tok.number_value = std::stod(std::string(tok.text));
    } catch (...) {
        tok.number_value = 0.0;
    }

    return tok;
}

CommandToken CommandParser::parse_string_token() {
    const std::size_t start = pos_;
    advance(); // 消费开头的 "

    std::string content;
    bool escaping = false;

    while (!is_at_end()) {
        const char ch = advance();

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
            // 字符串结束
            CommandToken tok;
            tok.type = CommandTokenType::kString;
            tok.text = source_.substr(start, pos_ - start);  // 包含引号
            tok.string_value = std::move(content);
            tok.position = start;
            return tok;
        } else {
            content.push_back(ch);
        }
    }

    // 未闭合的字符串
    throw SyntaxError("unterminated string literal starting at position " + std::to_string(start));
}

CommandToken CommandParser::parse_operator_token() {
    const std::size_t start = pos_;
    const char ch = advance();

    // 检查双字符运算符
    if (!is_at_end()) {
        const char next = peek();
        if ((ch == '=' && next == '=') ||
            (ch == '!' && next == '=') ||
            (ch == '<' && next == '=') ||
            (ch == '>' && next == '=') ||
            (ch == '&' && next == '&') ||
            (ch == '|' && next == '|')) {
            advance();
        }
    }

    CommandToken tok;
    tok.type = CommandTokenType::kOperator;
    tok.text = source_.substr(start, pos_ - start);
    tok.position = start;
    return tok;
}

// ============================================================================
// CommandParser 实现 - 语法分析
// ============================================================================

const CommandToken& CommandParser::peek_token() {
    if (!tokens_scanned_) {
        // 单次扫描生成所有 Token
        while (true) {
            CommandToken tok = next_token();
            tokens_.push_back(tok);
            if (tok.type == CommandTokenType::kEnd) {
                break;
            }
        }
        tokens_scanned_ = true;
    }

    if (token_pos_ >= tokens_.size()) {
        static CommandToken end_token;
        end_token.type = CommandTokenType::kEnd;
        end_token.position = source_.size();
        return end_token;
    }
    return tokens_[token_pos_];
}

CommandToken CommandParser::advance_token() {
    const CommandToken& tok = peek_token();
    if (token_pos_ < tokens_.size()) {
        ++token_pos_;
    }
    return tok;
}

bool CommandParser::check_token(CommandTokenType type) const {
    if (token_pos_ >= tokens_.size()) return false;
    return tokens_[token_pos_].type == type;
}

bool CommandParser::match_token(CommandTokenType type) {
    if (check_token(type)) {
        advance_token();
        return true;
    }
    return false;
}

CommandToken CommandParser::expect_token(CommandTokenType type, const char* message) {
    if (!check_token(type)) {
        throw_syntax_error(message);
    }
    return advance_token();
}

void CommandParser::throw_syntax_error(const std::string& message) {
    const std::size_t pos = token_pos_ < tokens_.size() ? tokens_[token_pos_].position : source_.size();
    throw SyntaxError(message + " at position " + std::to_string(pos));
}

CommandASTNode CommandParser::parse() {
    // 初始化 Token 流
    peek_token();

    skip_spaces();
    if (is_at_end() && tokens_.empty()) {
        return CommandASTNode::make_empty();
    }

    return parse_command();
}

CommandASTNode CommandParser::parse_command() {
    const CommandToken& tok = peek_token();

    // 1. 空输入
    if (tok.type == CommandTokenType::kEnd) {
        return CommandASTNode::make_empty();
    }

    // 2. 元命令（以冒号开头）
    if (tok.type == CommandTokenType::kColon) {
        return parse_meta_command();
    }

    // 3. 字符串字面量 - 需要检查后面是否还有内容
    if (tok.type == CommandTokenType::kString) {
        CommandToken str_tok = advance_token();
        skip_spaces();

        // 检查后面是否还有内容（如 "abc" + "def"）
        if (peek_token().type != CommandTokenType::kEnd) {
            // 后面还有内容，作为表达式处理
            return parse_expression();
        }

        return CommandASTNode::make_string_literal(str_tok.string_value);
    }

    // 4. 标识符开头：可能是函数定义、赋值、函数调用或表达式
    if (tok.type == CommandTokenType::kIdentifier) {
        return parse_definition_or_assignment(tok);
    }

    // 5. 其他情况作为表达式
    return parse_expression();
}

CommandASTNode CommandParser::parse_meta_command() {
    expect_token(CommandTokenType::kColon, "expected ':' for meta command");

    CommandToken cmd_tok = expect_token(CommandTokenType::kIdentifier,
                                         "expected command name after ':'");

    std::vector<std::string_view> arguments;

    // 收集所有剩余参数
    while (peek_token().type != CommandTokenType::kEnd) {
        const CommandToken& arg_tok = peek_token();

        if (arg_tok.type == CommandTokenType::kIdentifier ||
            arg_tok.type == CommandTokenType::kNumber ||
            arg_tok.type == CommandTokenType::kString) {
            arguments.push_back(arg_tok.text);
            advance_token();
        } else if (arg_tok.type == CommandTokenType::kComma) {
            advance_token(); // 跳过逗号
        } else {
            // 其他 Token 作为参数文本
            arguments.push_back(arg_tok.text);
            advance_token();
        }
    }

    return CommandASTNode::make_meta_command(cmd_tok.text, arguments);
}

CommandASTNode CommandParser::parse_definition_or_assignment(const CommandToken& id_tok) {
    // 消费标识符
    advance_token();

    const CommandToken& next = peek_token();

    // 情况 1: f(...) 可能是函数定义或函数调用
    if (next.type == CommandTokenType::kLParen) {
        // 先保存位置，尝试解析为函数定义
        std::size_t saved_token_pos = token_pos_;
        advance_token(); // 消费 (

        // 检查是否是函数定义格式: f(x) = ...
        // 要求: 单个标识符参数，紧接着 ) 和 =
        if (check_token(CommandTokenType::kIdentifier)) {
            CommandToken param_tok = advance_token();

            if (check_token(CommandTokenType::kRParen)) {
                advance_token(); // 消费 )

                if (check_token(CommandTokenType::kEqual)) {
                    advance_token(); // 消费 =

                    // 这是函数定义
                    std::string_view body;
                    if (peek_token().type != CommandTokenType::kEnd) {
                        body = source_.substr(peek_token().position);
                        while (peek_token().type != CommandTokenType::kEnd) {
                            advance_token();
                        }
                    }

                    return CommandASTNode::make_function_definition(id_tok.text, param_tok.text, body);
                }
            }
        }

        // 不是函数定义，回退并作为函数调用处理
        token_pos_ = saved_token_pos;
        return parse_function_call(id_tok);
    }

    // 情况 2: x = ... 赋值
    if (next.type == CommandTokenType::kEqual) {
        advance_token(); // 消费 =

        // 收集表达式
        std::string_view expr;
        if (peek_token().type != CommandTokenType::kEnd) {
            expr = source_.substr(peek_token().position);
            // 消费所有剩余 Token
            while (peek_token().type != CommandTokenType::kEnd) {
                advance_token();
            }
        }

        return CommandASTNode::make_assignment(id_tok.text, expr);
    }

    // 情况 3: 其他情况作为表达式（单独的标识符或表达式的一部分）
    std::string_view expr = id_tok.text;
    if (peek_token().type != CommandTokenType::kEnd) {
        // 后面还有内容，收集完整表达式
        expr = source_.substr(id_tok.position);
        while (peek_token().type != CommandTokenType::kEnd) {
            advance_token();
        }
    }
    return CommandASTNode::make_expression(expr);
}

CommandASTNode CommandParser::parse_function_call(const CommandToken& id_tok) {
    // parse_definition_or_assignment has already consumed the identifier and
    // rewound to the opening parenthesis after ruling out a function definition.
    expect_token(CommandTokenType::kLParen, "expected '(' in function call");

    std::vector<std::string_view> arguments;

    // 解析参数列表
    if (!check_token(CommandTokenType::kRParen)) {
        // 使用位置信息来提取参数
        std::size_t arg_start = peek_token().position;
        int paren_depth = 1;
        int bracket_depth = 0;

        while (peek_token().type != CommandTokenType::kEnd) {
            const CommandToken& tok = peek_token();

            if (tok.type == CommandTokenType::kLParen) {
                ++paren_depth;
            } else if (tok.type == CommandTokenType::kRParen) {
                --paren_depth;
                if (paren_depth == 0 && bracket_depth == 0) {
                    // 保存最后一个参数
                    if (tok.position > arg_start) {
                        arguments.push_back(source_.substr(arg_start, tok.position - arg_start));
                    }
                    advance_token(); // 消费 )
                    break;
                }
            } else if (tok.type == CommandTokenType::kLBracket) {
                ++bracket_depth;
            } else if (tok.type == CommandTokenType::kRBracket) {
                --bracket_depth;
            } else if (tok.type == CommandTokenType::kComma && paren_depth == 1 && bracket_depth == 0) {
                // 顶级逗号，分割参数
                if (tok.position > arg_start) {
                    arguments.push_back(source_.substr(arg_start, tok.position - arg_start));
                }
                advance_token();
                arg_start = peek_token().position;
                continue;
            }

            advance_token();
        }

        // 检查括号是否匹配
        if (paren_depth > 0) {
            throw_syntax_error("unmatched '(' in function call");
        }
    } else {
        advance_token(); // 消费空的 )
    }

    // 检查后面是否还有内容
    if (peek_token().type != CommandTokenType::kEnd) {
        // 后面还有内容，整个输入是以函数调用开头的表达式。
        const std::string_view expr = source_.substr(id_tok.position);
        while (peek_token().type != CommandTokenType::kEnd) {
            advance_token();
        }
        return CommandASTNode::make_expression(expr);
    }

    return CommandASTNode::make_function_call(id_tok.text, arguments);
}

CommandASTNode CommandParser::parse_expression() {
    // 收集所有剩余 Token 作为表达式
    if (peek_token().type == CommandTokenType::kEnd) {
        return CommandASTNode::make_expression("");
    }

    std::size_t start_pos = peek_token().position;

    // 消费所有剩余 Token
    while (peek_token().type != CommandTokenType::kEnd) {
        advance_token();
    }

    std::string_view expr = source_.substr(start_pos);
    return CommandASTNode::make_expression(expr);
}

// ============================================================================
// 便捷函数实现
// ============================================================================

CommandASTNode parse_command(std::string_view source) {
    CommandParser parser(source);
    return parser.parse();
}

std::vector<std::string_view> split_command_arguments(std::string_view text) {
    return split_top_level_arguments_view(text);
}

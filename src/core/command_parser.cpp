// ============================================================================
// 命令解析器实现
// ============================================================================

#include "command_parser.h"
#include "utils.h"
#include <algorithm>
#include <cctype>
#include <sstream>

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
                                                        const std::vector<std::string_view>& params,
                                                        std::string_view body) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionDefinition;
    node.data = FunctionDefinitionInfo{name, params, body};
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
    : BaseParser(source) {}

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
    if (peek_is_digit() || (ch == '.' && std::isdigit(static_cast<unsigned char>(peek_next())))) {
        return parse_number_token();
    }

    // 标识符
    if (peek_is_identifier_start()) {
        return parse_identifier_token();
    }

    // 单字符 Token 和运算符
    switch (ch) {
        case '(': {
            pos_++;
            CommandToken tok;
            tok.type = CommandTokenType::kLParen;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ')': {
            pos_++;
            CommandToken tok;
            tok.type = CommandTokenType::kRParen;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case '[': {
            pos_++;
            CommandToken tok;
            tok.type = CommandTokenType::kLBracket;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ']': {
            pos_++;
            CommandToken tok;
            tok.type = CommandTokenType::kRBracket;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ',': {
            pos_++;
            CommandToken tok;
            tok.type = CommandTokenType::kComma;
            tok.text = source_.substr(start_pos, 1);
            tok.position = start_pos;
            return tok;
        }
        case ':': {
            pos_++;
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
            pos_++;
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
    parse_identifier(); // 使用 BaseParser 的工具
    CommandToken tok;
    tok.type = CommandTokenType::kIdentifier;
    tok.text = source_.substr(start, pos_ - start);
    tok.position = start;
    return tok;
}

CommandToken CommandParser::parse_number_token() {
    const std::size_t start = pos_;

    // 处理进制前缀 (0x, 0b, 0o)
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

    // 整数部分
    while (!is_at_end() && peek_is_digit()) {
        pos_++;
    }

    // 小数部分
    if (!is_at_end() && peek() == '.') {
        pos_++;
        while (!is_at_end() && peek_is_digit()) {
            pos_++;
        }
    }

    // 指数部分
    if (!is_at_end() && (peek() == 'e' || peek() == 'E')) {
        pos_++;
        if (!is_at_end() && (peek() == '+' || peek() == '-')) {
            pos_++;
        }
        while (!is_at_end() && peek_is_digit()) {
            pos_++;
        }
    }
    
    // 虚数后缀
    if (!is_at_end() && peek() == 'i') {
        pos_++;
    }

end_num:
    CommandToken tok;
    tok.type = CommandTokenType::kNumber;
    tok.text = source_.substr(start, pos_ - start);
    tok.position = start;

    return tok;
}

CommandToken CommandParser::parse_string_token() {
    const std::size_t start = pos_;
    pos_++; // 消费开头的 "

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
    throw_syntax_error("unterminated string literal");
}

CommandToken CommandParser::parse_operator_token() {
    const std::size_t start = pos_;
    const char ch = source_[pos_++];

    // 检查双字符运算符
    if (!is_at_end()) {
        const char next = peek();
        if ((ch == '=' && next == '=') ||
            (ch == '!' && next == '=') ||
            (ch == '<' && next == '=') ||
            (ch == '>' && next == '=') ||
            (ch == '&' && next == '&') ||
            (ch == '|' && next == '|')) {
            pos_++;
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
    if (token_pos_ >= tokens_.size()) return type == CommandTokenType::kEnd;
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
    const std::size_t err_pos = token_pos_ < tokens_.size() ? tokens_[token_pos_].position : source_.size();
    
    std::ostringstream oss;
    oss << message << " at position " << err_pos << "\n";
    oss << "  " << source_ << "\n";
    oss << "  " << std::string(err_pos, ' ') << "^";
    
    throw SyntaxError(oss.str());
}

CommandASTNode CommandParser::parse() {
    // 触发词法分析
    peek_token();

    if (tokens_.empty() || tokens_[0].type == CommandTokenType::kEnd) {
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
        
        // 检查后面是否还有内容（如 "abc" + "def"）
        if (peek_token().type != CommandTokenType::kEnd) {
            // 后面还有内容，回退到字符串开始，作为表达式处理
            token_pos_--; 
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
    if (peek_token().type != CommandTokenType::kEnd) {
        arguments = parse_argument_list_by_tokens(false);
    }

    return CommandASTNode::make_meta_command(cmd_tok.text, arguments);
}

CommandASTNode CommandParser::parse_definition_or_assignment(const CommandToken& id_tok) {
    // 消费标识符
    advance_token();

    const CommandToken& next = peek_token();

    // 情况 1: f(...) 可能是函数定义或函数调用
    if (next.type == CommandTokenType::kLParen) {
        // 保存当前位置，以便回退
        std::size_t saved_token_pos = token_pos_;
        advance_token(); // 消费 (

        // 尝试解析参数名列表 (x, y, z)
        std::vector<std::string_view> params;
        bool is_valid_def_params = true;
        
        if (!check_token(CommandTokenType::kRParen)) {
            while (true) {
                if (check_token(CommandTokenType::kIdentifier)) {
                    params.push_back(advance_token().text);
                } else {
                    is_valid_def_params = false;
                    break;
                }
                
                if (match_token(CommandTokenType::kComma)) {
                    continue;
                } else if (check_token(CommandTokenType::kRParen)) {
                    break;
                } else {
                    is_valid_def_params = false;
                    break;
                }
            }
        }
        
        if (is_valid_def_params && match_token(CommandTokenType::kRParen)) {
            if (match_token(CommandTokenType::kEqual)) {
                // 这是函数定义
                std::string_view body;
                if (peek_token().type != CommandTokenType::kEnd) {
                    body = trim_view(source_.substr(peek_token().position));
                    while (peek_token().type != CommandTokenType::kEnd) {
                        advance_token();
                    }
                }
                return CommandASTNode::make_function_definition(id_tok.text, params, body);
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
            expr = trim_view(source_.substr(peek_token().position));
            while (peek_token().type != CommandTokenType::kEnd) {
                advance_token();
            }
        }
        return CommandASTNode::make_assignment(id_tok.text, expr);
    }

    // 情况 3: 其他情况作为表达式
    token_pos_--; // 回退标识符
    return parse_expression();
}

CommandASTNode CommandParser::parse_function_call(const CommandToken& id_tok) {
    expect_token(CommandTokenType::kLParen, "expected '(' in function call");

    std::vector<std::string_view> arguments;

    // 解析参数列表
    if (!check_token(CommandTokenType::kRParen)) {
        arguments = parse_argument_list_by_tokens(true);
        if (!match_token(CommandTokenType::kRParen)) {
            throw_syntax_error("unmatched '(' in function call");
        }
    } else {
        advance_token(); // 消费空的 )
    }

    // 检查后面是否还有内容
    if (peek_token().type != CommandTokenType::kEnd) {
        // 后面还有内容，整个输入是以函数调用开头的表达式。
        token_pos_ = 0; // 完全回退，交给 parse_expression
        return parse_expression();
    }

    return CommandASTNode::make_function_call(id_tok.text, arguments);
}

std::vector<std::string_view> CommandParser::parse_argument_list_by_tokens(bool stop_at_rparen) {
    std::vector<std::string_view> arguments;
    if (peek_token().type == CommandTokenType::kEnd) {
        return arguments;
    }

    std::size_t start_pos = peek_token().position;
    int paren_depth = 0;
    int bracket_depth = 0;

    while (peek_token().type != CommandTokenType::kEnd) {
        const CommandToken& tok = peek_token();

        if (tok.type == CommandTokenType::kLParen) paren_depth++;
        else if (tok.type == CommandTokenType::kRParen) {
            if (paren_depth == 0 && stop_at_rparen) {
                // Reached the end of the argument list for a function call
                break;
            }
            if (paren_depth > 0) paren_depth--;
        }
        else if (tok.type == CommandTokenType::kLBracket) bracket_depth++;
        else if (tok.type == CommandTokenType::kRBracket) {
            if (bracket_depth > 0) bracket_depth--;
        }
        else if (tok.type == CommandTokenType::kComma && paren_depth == 0 && bracket_depth == 0) {
            // Found a top-level comma
            std::size_t end_pos = tok.position;
            arguments.push_back(trim_view(source_.substr(start_pos, end_pos - start_pos)));
            advance_token(); // consume comma
            if (peek_token().type != CommandTokenType::kEnd) {
                start_pos = peek_token().position;
            } else {
                start_pos = source_.size();
            }
            continue;
        }

        advance_token();
    }

    // push the last argument
    std::size_t end_pos = peek_token().type == CommandTokenType::kEnd ? source_.size() : peek_token().position;
    if (start_pos < end_pos) {
        arguments.push_back(trim_view(source_.substr(start_pos, end_pos - start_pos)));
    } else if (!arguments.empty()) {
        // trailing comma case like `f(a, )`
        arguments.push_back("");
    }

    return arguments;
}

CommandASTNode CommandParser::parse_expression() {
    if (peek_token().type == CommandTokenType::kEnd) {
        return CommandASTNode::make_expression("");
    }

    std::size_t start_pos = tokens_[token_pos_].position;

    // 消费所有剩余 Token
    while (peek_token().type != CommandTokenType::kEnd) {
        advance_token();
    }

    std::string_view expr = trim_view(source_.substr(start_pos));

    // 基本语法预检：检测明显的语法错误
    if (has_obvious_syntax_error(expr)) {
        // 返回表达式节点，但错误会在求值阶段报告
        // 这里只是记录警告，不阻止解析
    }

    return CommandASTNode::make_expression(expr);
}

// ============================================================================
// 表达式预检验证
// ============================================================================

bool CommandParser::has_obvious_syntax_error(std::string_view expr) const {
    if (expr.empty()) return false;

    // 检查括号匹配
    if (!check_paren_balance(expr)) {
        return true;
    }

    // 检查连续运算符（如 "++", "--" 不是合法的运算符组合）
    // 注意："--" 可能是合法的负负，所以不检查

    // 检查表达式是否以非法字符结尾
    char last = expr.back();
    if (last == '+' || last == '-' || last == '*' || last == '/' ||
        last == '^' || last == '=') {
        return true;
    }

    return false;
}

bool CommandParser::check_paren_balance(std::string_view expr) const {
    int paren_depth = 0;
    int bracket_depth = 0;
    bool in_string = false;
    bool escaping = false;

    for (std::size_t i = 0; i < expr.size(); ++i) {
        const char ch = expr[i];

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
            continue;
        }

        if (ch == '(') ++paren_depth;
        else if (ch == ')') {
            if (paren_depth == 0) return false;  // 多余的右括号
            --paren_depth;
        }
        else if (ch == '[') ++bracket_depth;
        else if (ch == ']') {
            if (bracket_depth == 0) return false;  // 多余的右方括号
            --bracket_depth;
        }
    }

    return paren_depth == 0 && bracket_depth == 0 && !in_string;
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

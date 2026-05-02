// ============================================================================
// 命令解析器实现
// ============================================================================

#include "parser/command_parser.h"
#include "parser/syntax_validator.h"
#include "core/utils.h"
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
    FunctionDefinitionInfo info;
    info.name = name;
    info.parameters = params;
    info.body = ExpressionInfo(body);
    node.data = std::move(info);
    return node;
}

CommandASTNode CommandASTNode::make_function_call(std::string_view name,
                                                  const std::vector<std::string_view>& args) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionCall;
    FunctionCallInfo info;
    info.name = name;
    for (const auto& arg : args) {
        info.arguments.push_back(ExpressionInfo(arg));
    }
    node.data = std::move(info);
    return node;
}

CommandASTNode CommandASTNode::make_assignment(std::string_view var,
                                               std::string_view expr) {
    CommandASTNode node;
    node.kind = CommandKind::kAssignment;
    AssignmentInfo info;
    info.variable = var;
    info.expression = ExpressionInfo(expr);
    node.data = std::move(info);
    return node;
}

CommandASTNode CommandASTNode::make_expression(std::string_view expr) {
    CommandASTNode node;
    node.kind = CommandKind::kExpression;
    node.data = ExpressionInfo(expr);
    return node;
}

CommandASTNode CommandASTNode::make_string_literal(const std::string& value) {
    CommandASTNode node;
    node.kind = CommandKind::kStringLiteral;
    node.data = value;
    return node;
}

// ============================================================================
// CommandParser 实现
// ============================================================================

CommandParser::CommandParser(std::string_view source)
    : tokens_(source) {}

// ============================================================================
// Token 访问
// ============================================================================

const Token& CommandParser::peek_token() {
    return tokens_.peek();
}

const Token& CommandParser::peek_token(std::size_t offset) {
    return tokens_.peek(offset);
}

Token CommandParser::advance_token() {
    return tokens_.advance();
}

bool CommandParser::check_token(TokenKind kind) const {
    // 需要通过 tokens_ 的 peek 方法检查
    // 由于 tokens_ 不是 const 方法，这里用一种间接方式
    return const_cast<CommandParser*>(this)->tokens_.check(kind);
}

bool CommandParser::match_token(TokenKind kind) {
    return tokens_.match(kind);
}

Token CommandParser::expect_token(TokenKind kind, const char* message) {
    return tokens_.expect(kind, message);
}

// ============================================================================
// 回溯支持
// ============================================================================

LazyTokenStream::Checkpoint CommandParser::save_checkpoint() {
    return tokens_.save_checkpoint();
}

void CommandParser::restore_checkpoint(const LazyTokenStream::Checkpoint& cp) {
    tokens_.restore_checkpoint(cp);
}

std::size_t CommandParser::current_position() const {
    return const_cast<CommandParser*>(this)->peek_token().position;
}

// ============================================================================
// 语法分析
// ============================================================================

void CommandParser::throw_syntax_error(const std::string& message) {
    Token tok = peek_token();
    std::ostringstream oss;
    oss << message << " at position " << tok.position << "\n";
    oss << "  " << tokens_.source() << "\n";
    oss << "  " << std::string(tok.position, ' ') << "^";

    throw SyntaxError(oss.str());
}

CommandASTNode CommandParser::parse() {
    const Token& tok = peek_token();

    if (tok.kind == TokenKind::kEnd) {
        return CommandASTNode::make_empty();
    }

    return parse_command();
}

CommandASTNode CommandParser::parse_command() {
    const Token& tok = peek_token();

    // 1. 空输入
    if (tok.kind == TokenKind::kEnd) {
        return CommandASTNode::make_empty();
    }

    // 2. 元命令（以冒号开头）
    if (tok.kind == TokenKind::kColon) {
        return parse_meta_command();
    }

    // 3. 字符串字面量 - 需要检查后面是否还有内容
    if (tok.kind == TokenKind::kString) {
        Token str_tok = advance_token();

        // 检查后面是否还有内容（如 "abc" + "def"）
        if (peek_token().kind != TokenKind::kEnd) {
            // 后面还有内容，回退到字符串开始，作为表达式处理
            restore_checkpoint(save_checkpoint());  // 重置到开始
            return parse_expression();
        }

        return CommandASTNode::make_string_literal(str_tok.string_value);
    }

    // 4. 标识符开头：可能是函数定义、赋值、函数调用或表达式
    if (tok.kind == TokenKind::kIdentifier) {
        return parse_definition_or_assignment(tok);
    }

    // 5. 其他情况作为表达式
    return parse_expression();
}

CommandASTNode CommandParser::parse_meta_command() {
    expect_token(TokenKind::kColon, "expected ':' for meta command");

    Token cmd_tok = expect_token(TokenKind::kIdentifier,
                                 "expected command name after ':'");

    std::vector<std::string_view> arguments;
    if (peek_token().kind != TokenKind::kEnd) {
        arguments = parse_argument_list_by_tokens(false);
    }

    return CommandASTNode::make_meta_command(cmd_tok.text, arguments);
}

CommandASTNode CommandParser::parse_definition_or_assignment(Token id_tok) {
    // 保存检查点，用于回退
    LazyTokenStream::Checkpoint checkpoint = save_checkpoint();

    // 消费标识符
    advance_token();

    const Token& next = peek_token();

    // 情况 1: f(...) 可能是函数定义或函数调用
    if (next.kind == TokenKind::kLParen) {
        // 保存第二个检查点（在 '(' 之前）
        LazyTokenStream::Checkpoint paren_checkpoint = save_checkpoint();

        advance_token(); // 消费 (

        // 尝试解析参数名列表 (x, y, z)
        std::vector<std::string_view> params;
        bool is_valid_def_params = true;

        if (!check_token(TokenKind::kRParen)) {
            while (true) {
                if (check_token(TokenKind::kIdentifier)) {
                    params.push_back(advance_token().text);
                } else {
                    is_valid_def_params = false;
                    break;
                }

                if (match_token(TokenKind::kComma)) {
                    continue;
                } else if (check_token(TokenKind::kRParen)) {
                    break;
                } else {
                    is_valid_def_params = false;
                    break;
                }
            }
        }

        if (is_valid_def_params && match_token(TokenKind::kRParen)) {
            if (match_token(TokenKind::kEqual)) {
                // 这是函数定义
                std::string_view body;
                if (peek_token().kind != TokenKind::kEnd) {
                    body = trim_view(tokens_.source_from(peek_token().position));
                    while (peek_token().kind != TokenKind::kEnd) {
                        advance_token();
                    }
                }
                return CommandASTNode::make_function_definition(id_tok.text, params, body);
            }
        }

        // 不是函数定义，回退到 '(' 之前并作为函数调用处理
        restore_checkpoint(paren_checkpoint);
        return parse_function_call(id_tok);
    }

    // 情况 2: x = ... 赋值
    if (next.kind == TokenKind::kEqual) {
        advance_token(); // 消费 =

        // 收集表达式
        std::string_view expr;
        if (peek_token().kind != TokenKind::kEnd) {
            expr = trim_view(tokens_.source_from(peek_token().position));
            while (peek_token().kind != TokenKind::kEnd) {
                advance_token();
            }
        }
        return CommandASTNode::make_assignment(id_tok.text, expr);
    }

    // 情况 3: 其他情况作为表达式
    restore_checkpoint(checkpoint);
    return parse_expression();
}

CommandASTNode CommandParser::parse_function_call(Token id_tok) {
    expect_token(TokenKind::kLParen, "expected '(' in function call");

    std::vector<std::string_view> arguments;

    // 解析参数列表
    if (!check_token(TokenKind::kRParen)) {
        arguments = parse_argument_list_by_tokens(true);
        if (!match_token(TokenKind::kRParen)) {
            throw_syntax_error("unmatched '(' in function call");
        }
    } else {
        advance_token(); // 消费空的 )
    }

    // 检查后面是否还有内容
    if (peek_token().kind != TokenKind::kEnd) {
        // 后面还有内容，整个输入是以函数调用开头的表达式。
        // 完全回退，交给 parse_expression
        tokens_.reset();
        return parse_expression();
    }

    return CommandASTNode::make_function_call(id_tok.text, arguments);
}

std::vector<std::string_view> CommandParser::parse_argument_list_by_tokens(bool stop_at_rparen) {
    std::vector<std::string_view> arguments;
    if (peek_token().kind == TokenKind::kEnd) {
        return arguments;
    }

    std::size_t start_pos = peek_token().position;
    int paren_depth = 0;
    int bracket_depth = 0;

    while (peek_token().kind != TokenKind::kEnd) {
        const Token& tok = peek_token();

        if (tok.kind == TokenKind::kLParen) paren_depth++;
        else if (tok.kind == TokenKind::kRParen) {
            if (paren_depth == 0 && stop_at_rparen) {
                // Reached the end of the argument list for a function call
                break;
            }
            if (paren_depth > 0) paren_depth--;
        }
        else if (tok.kind == TokenKind::kLBracket) bracket_depth++;
        else if (tok.kind == TokenKind::kRBracket) {
            if (bracket_depth > 0) bracket_depth--;
        }
        else if (tok.kind == TokenKind::kComma && paren_depth == 0 && bracket_depth == 0) {
            // Found a top-level comma
            std::size_t end_pos = tok.position;
            arguments.push_back(trim_view(tokens_.source().substr(start_pos, end_pos - start_pos)));
            advance_token(); // consume comma
            if (peek_token().kind != TokenKind::kEnd) {
                start_pos = peek_token().position;
            } else {
                start_pos = tokens_.source().size();
            }
            continue;
        }

        advance_token();
    }

    // push the last argument
    std::size_t end_pos = peek_token().kind == TokenKind::kEnd ? tokens_.source().size() : peek_token().position;
    if (start_pos < end_pos) {
        arguments.push_back(trim_view(tokens_.source().substr(start_pos, end_pos - start_pos)));
    } else if (!arguments.empty()) {
        // trailing comma case like `f(a, )`
        arguments.push_back("");
    }

    return arguments;
}

CommandASTNode CommandParser::parse_expression() {
    if (peek_token().kind == TokenKind::kEnd) {
        return CommandASTNode::make_expression("");
    }

    std::size_t start_pos = peek_token().position;

    // 消费所有剩余 Token
    while (peek_token().kind != TokenKind::kEnd) {
        advance_token();
    }

    std::string_view expr = trim_view(tokens_.source_from(start_pos));

    // 增强的语法预检：使用 SyntaxValidator
    static SyntaxValidator validator;
    if (validator.has_errors(expr)) {
        std::string error_msg = validator.get_first_error(expr);
        // 记录错误信息，但不阻止解析
        // 错误将在求值阶段报告
    }

    // 创建带预编译缓存的 ExpressionInfo
    ExpressionInfo info(expr);
    compile_expression_info(info);

    return CommandASTNode(CommandKind::kExpression, std::move(info));
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

// ============================================================================
// 表达式预编译实现
// ============================================================================

void compile_expression_info(ExpressionInfo& info) {
    if (info.text.empty() || info.cache) {
        return;  // 已缓存或空表达式
    }

    info.cache = std::make_shared<ExpressionCache>();
    info.cache->expanded = std::string(info.text);
    info.cache->hint = analyze_expression_hint(info.cache->expanded);
    info.cache->features = analyze_expression_features(info.cache->expanded);

    // 尝试编译 AST（仅对标量表达式）
    if (info.cache->hint == ExpressionHint::kScalar) {
        info.cache->compiled_ast = compile_expression_ast(info.cache->expanded);
        if (info.cache->compiled_ast) {
            info.cache->is_compiled = true;
        }
    }
}

void compile_command_expressions(CommandASTNode& node) {
    switch (node.kind) {
        case CommandKind::kEmpty:
        case CommandKind::kMetaCommand:
        case CommandKind::kStringLiteral:
            // 无表达式需要编译
            break;

        case CommandKind::kFunctionDefinition: {
            auto* info = const_cast<FunctionDefinitionInfo*>(node.as_function_definition());
            if (info) {
                compile_expression_info(info->body);
            }
            break;
        }

        case CommandKind::kFunctionCall: {
            auto* info = const_cast<FunctionCallInfo*>(node.as_function_call());
            if (info) {
                for (auto& arg : info->arguments) {
                    compile_expression_info(arg);
                }
            }
            break;
        }

        case CommandKind::kAssignment: {
            auto* info = const_cast<AssignmentInfo*>(node.as_assignment());
            if (info) {
                compile_expression_info(info->expression);
            }
            break;
        }

        case CommandKind::kExpression: {
            auto* info = const_cast<ExpressionInfo*>(node.as_expression());
            if (info) {
                compile_expression_info(*info);
            }
            break;
        }
    }
}

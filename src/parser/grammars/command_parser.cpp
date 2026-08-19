// ============================================================================
// 命令解析器实现
// ============================================================================

#include "parser/grammars/command_parser.h"
#include "parser/infra/syntax_validator.h"
#include "core/services/string_utils.h"
#include "parser/infra/parser_utils.h"
#include <algorithm>
#include <cctype>
#include <sstream>

// ============================================================================
// 内部函数
// ============================================================================

namespace {

void compile_expression_info(ExpressionInfo& info) {
    if (info.text.empty()) {
        return;  // 空表达式
    }

    // 获取或创建缓存
    ExpressionCache* cache = info.get_or_create_cache();

    // 如果已经分析过，直接返回
    if (cache->analyzed) {
        return;
    }

    // 获取有效文本（展开后的或原始的）
    std::string_view expr = cache->effective_text();

    // 分析表达式
    cache->hint = analyze_expression_hint(std::string(expr));
    cache->features = analyze_expression_features(std::string(expr));
    cache->analyzed = true;

    // 所有表达式域统一预编译；矩阵/复数不再被排除在 AST 缓存之外。
    if (!parser_utils::contains_script_syntax(expr)) {
        cache->compiled_ast = compile_expression_ast(std::string(expr));
        if (cache->compiled_ast) {
            cache->is_compiled = true;
        }
    }
}

void validate_expression_text(std::string_view expr, LazyTokenStream& tokens, std::size_t source_offset) {
    static SyntaxValidator validator;
    std::vector<SyntaxErrorInfo> syntax_errors = validator.validate(expr);
    for (const auto& err : syntax_errors) {
        if (err.severity == Severity::kError) {
            const std::size_t error_pos = source_offset + err.position;
            std::ostringstream oss;
            oss << err.message << " at position " << error_pos << "\n";
            oss << "  " << tokens.source() << "\n";
            oss << "  " << std::string(error_pos, ' ') << "^";
            throw SyntaxError(oss.str());
        }
    }
}

void attach_source_owner(CommandASTNode& node,
                         const std::shared_ptr<const std::string>& source_owner) {
    node.source_owner = source_owner;
    if (node.kind == CommandKind::kSequence) {
        auto* nodes = &std::get<std::vector<CommandASTNode>>(node.data);
        for (auto& child : *nodes) {
            attach_source_owner(child, source_owner);
        }
    }
    // Note: For kFunctionCall, kFunctionDefinition, kAssignment, kIndexAssignment, kExpression,
    // the string_views in the data field point into source_owner, which is now held by the node.
    // No further recursion needed as these are leaf nodes in terms of CommandASTNode structure.
}

} // namespace

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
    compile_expression_info(info.body);
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
        ExpressionInfo expression(arg);
        compile_expression_info(expression);
        info.arguments.push_back(std::move(expression));
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
    compile_expression_info(info.expression);
    node.data = std::move(info);
    return node;
}

CommandASTNode CommandASTNode::make_index_assignment(
    std::string_view var,
    const std::vector<std::string_view>& indices,
    std::string_view value) {
    CommandASTNode node;
    node.kind = CommandKind::kIndexAssignment;
    IndexAssignmentInfo info;
    info.variable = var;
    for (const auto& idx : indices) {
        info.indices.push_back(ExpressionInfo(idx));
    }
    info.value = ExpressionInfo(value);
    compile_expression_info(info.value);
    node.data = std::move(info);
    return node;
}

CommandASTNode CommandASTNode::make_function_call_with_named(
    std::string_view name,
    const std::vector<std::string_view>& positional_args,
    const std::vector<std::pair<std::string_view, std::string_view>>& named_args) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionCall;
    FunctionCallInfo info;
    info.name = name;
    for (const auto& arg : positional_args) {
        ExpressionInfo expression(arg);
        compile_expression_info(expression);
        info.arguments.push_back(std::move(expression));
    }
    for (const auto& [name, value] : named_args) {
        ExpressionInfo expression(value);
        compile_expression_info(expression);
        info.named_args.push_back(NamedArgument{name, std::move(expression)});
    }
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

CommandASTNode CommandASTNode::make_sequence(std::vector<CommandASTNode> nodes) {
    CommandASTNode node;
    node.kind = CommandKind::kSequence;
    node.data = std::move(nodes);
    return node;
}

// ============================================================================
// CommandParser 实现
// ============================================================================

CommandParser::CommandParser(std::string_view source, IsCommandCallback is_command)
    : tokens_(source), is_command_(std::move(is_command)) {}

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
    if (peek_token().kind == TokenKind::kEnd) {
        return CommandASTNode::make_empty();
    }

    std::vector<CommandASTNode> nodes;
    while (peek_token().kind != TokenKind::kEnd) {
        // 跳过可能的分号
        while (match_token(TokenKind::kSemicolon)) {
            if (peek_token().kind == TokenKind::kEnd) break;
        }
        if (peek_token().kind == TokenKind::kEnd) break;

        nodes.push_back(parse_command(false));

        // 必须以分号分隔或到达末尾
        if (peek_token().kind != TokenKind::kEnd && !match_token(TokenKind::kSemicolon)) {
            throw_syntax_error("expected ';' or end of statement");
        }
    }

    if (nodes.empty()) return CommandASTNode::make_empty();
    if (nodes.size() == 1) return std::move(nodes[0]);
    return CommandASTNode::make_sequence(std::move(nodes));
}

CommandASTNode CommandParser::parse_command(bool single_statement) {
    LazyTokenStream::Checkpoint statement_checkpoint = save_checkpoint();
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

        // 检查后面是否还有内容
        const TokenKind next_kind = peek_token().kind;
        if (next_kind == TokenKind::kEnd || next_kind == TokenKind::kSemicolon) {
            return CommandASTNode::make_string_literal(str_tok.string_value);
        }

        // 后面还有内容，回退并作为表达式处理
        restore_checkpoint(statement_checkpoint);
        return parse_expression(single_statement);
    }

    // 4. 标识符开头：可能是函数定义、赋值、函数调用或表达式
    if (tok.kind == TokenKind::kIdentifier) {
        return parse_definition_or_assignment(tok, single_statement);
    }

    // 5. 其他情况作为表达式
    return parse_expression(single_statement);
}

CommandASTNode CommandParser::parse_meta_command() {
    expect_token(TokenKind::kColon, "expected ':' for meta command");

    Token cmd_tok = expect_token(TokenKind::kIdentifier,
                                 "expected command name after ':'");

    std::vector<std::string_view> arguments;
    const TokenKind next = peek_token().kind;
    if (next != TokenKind::kEnd && next != TokenKind::kSemicolon) {
        arguments = parse_argument_list_by_tokens(false);
    }

    return CommandASTNode::make_meta_command(cmd_tok.text, arguments);
}

CommandASTNode CommandParser::parse_definition_or_assignment(Token id_tok, bool single_statement) {
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
                std::string_view body = collect_statement_expression();
                if (body.empty()) {
                    throw_syntax_error("expected expression after '='");
                }
                validate_expression_text(body, tokens_, body.data() - tokens_.source().data());
                return CommandASTNode::make_function_definition(id_tok.text, params, body);
            }
        }

        // 不是函数定义，回退到 '(' 之前
        restore_checkpoint(paren_checkpoint);

        // 有命令注册表时只把注册命令识别为命令调用；没有回调时保留通用解析能力，
        // 供 base/rat/poly_* 等 helper 直接解析完整 name(args)。
        if (!is_command_ || is_command_(id_tok.text)) {
            return parse_function_call(id_tok, single_statement, checkpoint);
        }

        // 否则回退，作为普通表达式解析
        restore_checkpoint(checkpoint);
        return parse_expression(single_statement);
    }

    // 情况 2: x[...] = value 索引赋值
    if (next.kind == TokenKind::kLBracket) {
        return parse_index_assignment(id_tok, single_statement, checkpoint);
    }

    // 情况 3: x = ... 赋值
    if (next.kind == TokenKind::kEqual) {
        advance_token(); // 消费 =

        std::string_view expr = collect_statement_expression();
        if (expr.empty()) {
            throw_syntax_error("expected expression after '='");
        }
        validate_expression_text(expr, tokens_, expr.data() - tokens_.source().data());
        return CommandASTNode::make_assignment(id_tok.text, expr);
    }

    // 情况 4: 其他情况作为表达式
    restore_checkpoint(checkpoint);
    return parse_expression(single_statement);
}

CommandASTNode CommandParser::parse_index_assignment(
    Token id_tok,
    bool single_statement,
    const LazyTokenStream::Checkpoint& expression_checkpoint) {

    expect_token(TokenKind::kLBracket, "expected '[' in index assignment");

    // 解析索引列表
    std::vector<std::string_view> indices;
    if (!check_token(TokenKind::kRBracket)) {
        indices = parse_argument_list_by_tokens(true);
        if (!match_token(TokenKind::kRBracket)) {
            throw_syntax_error("unmatched '[' in index assignment");
        }
    } else {
        advance_token(); // 消费空的 ]
        // An empty subscript is never a valid assignment target.  Empty
        // list/dictionary literals are handled by the expression parser.
        if (check_token(TokenKind::kEqual)) {
            throw_syntax_error("index assignment requires an index");
        }
    }

    // 期望等号
    if (!match_token(TokenKind::kEqual)) {
        // 不是索引赋值，回退作为表达式解析
        restore_checkpoint(expression_checkpoint);
        return parse_expression(single_statement);
    }

    // 解析赋值值
    std::string_view value = collect_statement_expression();
    if (value.empty()) {
        throw_syntax_error("expected expression after '='");
    }
    validate_expression_text(value, tokens_, value.data() - tokens_.source().data());

    return CommandASTNode::make_index_assignment(id_tok.text, indices, value);
}

CommandASTNode CommandParser::parse_function_call(
    Token id_tok,
    bool single_statement,
    const LazyTokenStream::Checkpoint& expression_checkpoint) {
    expect_token(TokenKind::kLParen, "expected '(' in function call");

    std::vector<std::string_view> arguments;

    // 解析参数列表
    std::vector<std::pair<std::string_view, std::string_view>> named_arguments;
    if (!check_token(TokenKind::kRParen)) {
        arguments = parse_argument_list_by_tokens(true);
        if (!match_token(TokenKind::kRParen)) {
            throw_syntax_error("unmatched '(' in function call");
        }
    } else {
        advance_token(); // 消费空的 )
    }

    // Split `name=value` only when the equals sign is at the argument's
    // top level and the left side is a valid identifier.
    std::vector<std::string_view> positional_arguments;
    for (std::string_view argument : arguments) {
        const std::size_t equal = parser_utils::find_top_level(argument, '=');
        if (equal == std::string_view::npos) {
            positional_arguments.push_back(argument);
            continue;
        }
        const std::string_view name = trim_view(argument.substr(0, equal));
        const std::string_view value = trim_view(argument.substr(equal + 1));
        const bool valid_name = !name.empty() &&
            (std::isalpha(static_cast<unsigned char>(name.front())) || name.front() == '_') &&
            std::all_of(name.begin() + 1, name.end(), [](char ch) {
                return std::isalnum(static_cast<unsigned char>(ch)) || ch == '_';
            });
        if (!valid_name) {
            positional_arguments.push_back(argument);
            continue;
        }
        if (value.empty()) {
            throw_syntax_error("named argument requires a value");
        }
        named_arguments.emplace_back(name, value);
    }

    // 检查后面是否还有内容
    const TokenKind next = peek_token().kind;
    if (next != TokenKind::kEnd && next != TokenKind::kSemicolon) {
        // 后面还有内容，当前语句是以函数调用开头的普通表达式。
        // 回到当前语句起点，而不是重置整个输入；否则在分号序列中会
        // 反复跳回第一条语句。
        restore_checkpoint(expression_checkpoint);
        return parse_expression(single_statement);
    }

    if (!named_arguments.empty()) {
        return CommandASTNode::make_function_call_with_named(
            id_tok.text, positional_arguments, named_arguments);
    }
    return CommandASTNode::make_function_call(id_tok.text, positional_arguments);
}

std::vector<std::string_view> CommandParser::parse_argument_list_by_tokens(bool stop_at_rparen) {
    std::vector<std::string_view> arguments;
    if (peek_token().kind == TokenKind::kEnd) {
        return arguments;
    }

    std::size_t start_pos = peek_token().position;
    int paren_depth = 0;
    int bracket_depth = 0;
    int brace_depth = 0;

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
            if (bracket_depth == 0) {
                // Reached the end of the index list for index assignment
                break;
            }
            if (bracket_depth > 0) bracket_depth--;
        }
        else if (tok.kind == TokenKind::kLBrace) brace_depth++;
        else if (tok.kind == TokenKind::kRBrace) {
            if (brace_depth > 0) brace_depth--;
        }
        else if (tok.kind == TokenKind::kComma && paren_depth == 0 && bracket_depth == 0 && brace_depth == 0) {
            // Found a top-level comma
            std::size_t end_pos = tok.position;
            const std::string_view argument =
                trim_view(tokens_.source().substr(start_pos, end_pos - start_pos));
            if (argument.empty()) {
                throw_syntax_error("empty argument");
            }
            arguments.push_back(argument);
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
        const std::string_view argument =
            trim_view(tokens_.source().substr(start_pos, end_pos - start_pos));
        if (argument.empty()) {
            throw_syntax_error("empty argument");
        }
        arguments.push_back(argument);
    } else if (!arguments.empty()) {
        throw_syntax_error("empty argument after ','");
    }

    return arguments;
}

CommandASTNode CommandParser::parse_expression(bool single_statement) {
    (void)single_statement;
    if (peek_token().kind == TokenKind::kEnd) {
        return CommandASTNode::make_expression("");
    }

    std::string_view expr = collect_statement_expression();

    validate_expression_text(expr, tokens_, expr.data() - tokens_.source().data());

    // 创建带预编译缓存的 ExpressionInfo
    ExpressionInfo info(expr);
    compile_expression_info(info);

    return CommandASTNode(CommandKind::kExpression, std::move(info));
}

std::string_view CommandParser::collect_statement_expression() {
    std::size_t start_pos = peek_token().position;
    int paren_depth = 0;
    int bracket_depth = 0;
    int brace_depth = 0;

    while (peek_token().kind != TokenKind::kEnd) {
        const Token& tok = peek_token();
        if (tok.kind == TokenKind::kSemicolon && paren_depth == 0 && bracket_depth == 0 && brace_depth == 0) {
            break; // 遇到顶层分号，停止当前表达式解析
        }
        if (tok.kind == TokenKind::kLParen) paren_depth++;
        else if (tok.kind == TokenKind::kRParen) { if (paren_depth > 0) paren_depth--; }
        else if (tok.kind == TokenKind::kLBracket) bracket_depth++;
        else if (tok.kind == TokenKind::kRBracket) { if (bracket_depth > 0) bracket_depth--; }
        else if (tok.kind == TokenKind::kLBrace) brace_depth++;
        else if (tok.kind == TokenKind::kRBrace) { if (brace_depth > 0) brace_depth--; }
        
        advance_token();
    }

    std::size_t end_pos = (peek_token().kind == TokenKind::kEnd) ? tokens_.source().size() : peek_token().position;
    return trim_view(tokens_.source().substr(start_pos, end_pos - start_pos));
}

// ============================================================================
// 便捷函数实现
// ============================================================================

CommandASTNode parse_command(std::string_view source, CommandParser::IsCommandCallback is_command) {
    auto source_owner = std::make_shared<const std::string>(source);
    CommandParser parser(*source_owner, std::move(is_command));
    CommandASTNode node = parser.parse();
    attach_source_owner(node, source_owner);
    return node;
}

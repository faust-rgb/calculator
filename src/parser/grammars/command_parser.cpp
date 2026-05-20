// ============================================================================
// 命令解析器实现
// ============================================================================

#include "parser/grammars/command_parser.h"
#include "parser/infra/syntax_validator.h"
#include "core/services/string_utils.h"
#include "core/services/core_manager_interfaces.h"
#include "parser/infra/parser_utils.h"
#include <algorithm>
#include <cctype>
#include <sstream>

// ============================================================================
// 内部函数
// ============================================================================

namespace {

void compile_expression_info(ExpressionInfo& info) {
    if (info.text.empty()) return;
    ExpressionCache* cache = info.get_or_create_cache();
    if (cache->analyzed) return;
    std::string_view expr = cache->effective_text();
    cache->hint = analyze_expression_hint(std::string(expr));
    cache->features = analyze_expression_features(std::string(expr));
    cache->analyzed = true;
    if (cache->hint == ExpressionHint::kScalar && !parser_utils::contains_script_syntax(expr)) {
        cache->compiled_ast = compile_expression_ast(std::string(expr));
        if (cache->compiled_ast) cache->is_compiled = true;
    }
}

void validate_expression_text(std::string_view expr, LazyTokenStream& tokens, std::size_t source_offset) {
    SyntaxValidator validator;
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
        for (auto& child : *nodes) attach_source_owner(child, source_owner);
    } else if (node.kind == CommandKind::kMetaCommand) {
        auto* meta = &std::get<MetaCommandInfo>(node.data);
        for (auto& arg : meta->arguments) attach_source_owner(*arg, source_owner);
    } else if (node.kind == CommandKind::kFunctionCall) {
        auto* call = &std::get<FunctionCallInfo>(node.data);
        for (auto& arg : call->arguments) attach_source_owner(*arg, source_owner);
        for (auto& n_arg : call->named_args) attach_source_owner(*n_arg.value, source_owner);
    } else if (node.kind == CommandKind::kFunctionDefinition) {
        auto* def = &std::get<FunctionDefinitionInfo>(node.data);
        if (def->body_expr) attach_source_owner(*def->body_expr, source_owner);
    } else if (node.kind == CommandKind::kAssignment) {
        auto* assign = &std::get<AssignmentInfo>(node.data);
        if (assign->value_expr) attach_source_owner(*assign->value_expr, source_owner);
    } else if (node.kind == CommandKind::kIndexAssignment) {
        auto* idx_assign = &std::get<IndexAssignmentInfo>(node.data);
        for (auto& idx : idx_assign->index_exprs) attach_source_owner(*idx, source_owner);
        if (idx_assign->value_expr) attach_source_owner(*idx_assign->value_expr, source_owner);
    }
}

std::vector<std::unique_ptr<CommandASTNode>> to_unique_ptrs(std::vector<CommandASTNode> nodes) {
    std::vector<std::unique_ptr<CommandASTNode>> result;
    for (auto& node : nodes) {
        result.push_back(std::make_unique<CommandASTNode>(std::move(node)));
    }
    return result;
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
                                                  std::vector<CommandASTNode> args) {
    CommandASTNode node;
    node.kind = CommandKind::kMetaCommand;
    node.data = MetaCommandInfo{cmd, to_unique_ptrs(std::move(args))};
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
                                                  std::vector<CommandASTNode> args) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionCall;
    FunctionCallInfo info;
    info.name = name;
    info.arguments = to_unique_ptrs(std::move(args));
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
    for (const auto& idx : indices) info.indices.emplace_back(idx);
    info.value = ExpressionInfo(value);
    compile_expression_info(info.value);
    node.data = std::move(info);
    return node;
}

CommandASTNode CommandASTNode::make_function_call_with_named(
    std::string_view name,
    std::vector<CommandASTNode> positional_args,
    std::vector<std::pair<std::string_view, CommandASTNode>> named_args) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionCall;
    FunctionCallInfo info;
    info.name = name;
    info.arguments = to_unique_ptrs(std::move(positional_args));
    for (auto& pair : named_args) {
        info.named_args.push_back(NamedArgument{pair.first, std::make_unique<CommandASTNode>(std::move(pair.second))});
    }
    node.data = std::move(info);
    return node;
}

CommandASTNode CommandASTNode::make_function_definition_ast(std::string_view name,
                                                            const std::vector<std::string_view>& params,
                                                            CommandASTNode body_expr) {
    CommandASTNode node;
    node.kind = CommandKind::kFunctionDefinition;
    FunctionDefinitionInfo info;
    info.name = name;
    info.parameters = params;
    // 保留文本视图用于显示
    if (body_expr.kind == CommandKind::kExpression && body_expr.as_expression()) {
        info.body = ExpressionInfo(body_expr.as_expression()->text);
    }
    info.body_expr = std::make_unique<CommandASTNode>(std::move(body_expr));
    node.data = std::move(info);
    return node;
}

CommandASTNode CommandASTNode::make_assignment_ast(std::string_view var,
                                                   CommandASTNode value_expr) {
    CommandASTNode node;
    node.kind = CommandKind::kAssignment;
    AssignmentInfo info;
    info.variable = var;
    // 保留文本视图用于显示
    if (value_expr.kind == CommandKind::kExpression && value_expr.as_expression()) {
        info.expression = ExpressionInfo(value_expr.as_expression()->text);
    }
    info.value_expr = std::make_unique<CommandASTNode>(std::move(value_expr));
    node.data = std::move(info);
    return node;
}

CommandASTNode CommandASTNode::make_index_assignment_ast(
    std::string_view var,
    std::vector<CommandASTNode> index_exprs,
    CommandASTNode value_expr) {
    CommandASTNode node;
    node.kind = CommandKind::kIndexAssignment;
    IndexAssignmentInfo info;
    info.variable = var;
    // 保留文本视图用于显示
    for (const auto& idx_expr : index_exprs) {
        if (idx_expr.kind == CommandKind::kExpression && idx_expr.as_expression()) {
            info.indices.emplace_back(idx_expr.as_expression()->text);
        }
    }
    if (value_expr.kind == CommandKind::kExpression && value_expr.as_expression()) {
        info.value = ExpressionInfo(value_expr.as_expression()->text);
    }
    for (auto& idx : index_exprs) {
        info.index_exprs.push_back(std::make_unique<CommandASTNode>(std::move(idx)));
    }
    info.value_expr = std::make_unique<CommandASTNode>(std::move(value_expr));
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

// CommandConfig 实现
CommandParser::CommandConfig CommandParser::CommandConfig::from_registry(const ICommandRegistry& registry) {
    CommandConfig config;
    config.exact_commands = registry.get_all_commands();
    return config;
}

CommandParser::CommandParser(std::string_view source, const CommandConfig& config)
    : tokens_(source), config_(config) {}

CommandParser::CommandParser(std::string_view source, IsCommandCallback is_command)
    : tokens_(source), is_command_(std::move(is_command)) {}

const Token& CommandParser::peek_token() { return tokens_.peek(); }
const Token& CommandParser::peek_token(std::size_t offset) { return tokens_.peek(offset); }
Token CommandParser::advance_token() { return tokens_.advance(); }
bool CommandParser::check_token(TokenKind kind) const { return const_cast<CommandParser*>(this)->tokens_.check(kind); }
bool CommandParser::match_token(TokenKind kind) { return tokens_.match(kind); }
Token CommandParser::expect_token(TokenKind kind, const char* message) { return tokens_.expect(kind, message); }
LazyTokenStream::Checkpoint CommandParser::save_checkpoint() { return tokens_.save_checkpoint(); }
void CommandParser::restore_checkpoint(const LazyTokenStream::Checkpoint& cp) { tokens_.restore_checkpoint(cp); }
std::size_t CommandParser::current_position() const { return const_cast<CommandParser*>(this)->peek_token().position; }

void CommandParser::throw_syntax_error(const std::string& message) {
    Token tok = peek_token();
    std::ostringstream oss;
    oss << message << " at position " << tok.position << "\n" << "  " << tokens_.source() << "\n" << "  " << std::string(tok.position, ' ') << "^";
    throw SyntaxError(oss.str());
}

CommandASTNode CommandParser::parse() {
    if (peek_token().kind == TokenKind::kEnd) return CommandASTNode::make_empty();
    std::vector<CommandASTNode> nodes;
    while (peek_token().kind != TokenKind::kEnd) {
        while (match_token(TokenKind::kSemicolon)) { if (peek_token().kind == TokenKind::kEnd) break; }
        if (peek_token().kind == TokenKind::kEnd) break;
        nodes.push_back(parse_command(false));
    }
    if (nodes.empty()) return CommandASTNode::make_empty();
    if (nodes.size() == 1) return std::move(nodes[0]);
    return CommandASTNode::make_sequence(std::move(nodes));
}

CommandASTNode CommandParser::parse_command(bool single_statement) {
    LazyTokenStream::Checkpoint statement_checkpoint = save_checkpoint();
    const Token& tok = peek_token();
    if (tok.kind == TokenKind::kEnd) return CommandASTNode::make_empty();
    if (tok.kind == TokenKind::kColon) return parse_meta_command();
    if (tok.kind == TokenKind::kString) {
        Token str_tok = advance_token();
        const TokenKind next_kind = peek_token().kind;
        if (next_kind == TokenKind::kEnd || next_kind == TokenKind::kSemicolon) return CommandASTNode::make_string_literal(str_tok.string_value);
        restore_checkpoint(statement_checkpoint);
        return parse_expression(single_statement);
    }
    if (tok.kind == TokenKind::kIdentifier) return parse_definition_or_assignment(tok, single_statement);
    return parse_expression(single_statement);
}

CommandASTNode CommandParser::parse_meta_command() {
    expect_token(TokenKind::kColon, "expected ':' for meta command");
    Token cmd_tok = expect_token(TokenKind::kIdentifier, "expected command name after ':'");
    std::vector<CommandASTNode> arguments;
    if (peek_token().kind != TokenKind::kEnd && peek_token().kind != TokenKind::kSemicolon) arguments = parse_argument_list_by_tokens(false);
    return CommandASTNode::make_meta_command(cmd_tok.text, std::move(arguments));
}

CommandASTNode CommandParser::parse_definition_or_assignment(Token id_tok, bool single_statement) {
    LazyTokenStream::Checkpoint checkpoint = save_checkpoint();
    advance_token();
    const Token& next = peek_token();
    if (next.kind == TokenKind::kLParen) {
        LazyTokenStream::Checkpoint paren_checkpoint = save_checkpoint();
        advance_token();
        std::vector<std::string_view> params;
        bool is_valid_def_params = true;
        if (!check_token(TokenKind::kRParen)) {
            while (true) {
                if (check_token(TokenKind::kIdentifier)) params.push_back(advance_token().text);
                else { is_valid_def_params = false; break; }
                if (match_token(TokenKind::kComma)) continue;
                else if (check_token(TokenKind::kRParen)) break;
                else { is_valid_def_params = false; break; }
            }
        }
        if (is_valid_def_params && match_token(TokenKind::kRParen)) {
            if (match_token(TokenKind::kEqual)) {
                std::string_view body = collect_statement_expression();
                validate_expression_text(body, tokens_, body.data() - tokens_.source().data());
                // 使用新的带完整 AST 的工厂方法
                CommandASTNode body_expr = CommandASTNode::make_expression(body);
                compile_expression_info(*const_cast<ExpressionInfo*>(body_expr.as_expression()));
                return CommandASTNode::make_function_definition_ast(id_tok.text, params, std::move(body_expr));
            }
        }
        restore_checkpoint(paren_checkpoint);
        
        // 尝试解析为函数调用（如果是顶层唯一语句）
        // 不再强依赖 is_command 回调或配置表，让执行层决定是命令还是函数
        return parse_function_call(id_tok, single_statement, checkpoint);
    }
    if (next.kind == TokenKind::kLBracket) return parse_index_assignment(id_tok, single_statement, checkpoint);
    if (next.kind == TokenKind::kEqual) {
        advance_token();
        std::string_view expr = collect_statement_expression();
        validate_expression_text(expr, tokens_, expr.data() - tokens_.source().data());
        // 使用新的带完整 AST 的工厂方法
        CommandASTNode value_expr = CommandASTNode::make_expression(expr);
        compile_expression_info(*const_cast<ExpressionInfo*>(value_expr.as_expression()));
        return CommandASTNode::make_assignment_ast(id_tok.text, std::move(value_expr));
    }
    restore_checkpoint(checkpoint);
    return parse_expression(single_statement);
}

CommandASTNode CommandParser::parse_index_assignment(Token id_tok, bool single_statement, const LazyTokenStream::Checkpoint& expression_checkpoint) {
    expect_token(TokenKind::kLBracket, "expected '[' in index assignment");
    std::vector<std::string_view> indices;
    std::vector<CommandASTNode> index_exprs;
    if (!check_token(TokenKind::kRBracket)) {
        std::size_t start_pos = peek_token().position;
        int depth = 0;
        while (peek_token().kind != TokenKind::kEnd) {
            if (peek_token().kind == TokenKind::kLBracket) depth++;
            else if (peek_token().kind == TokenKind::kRBracket) { if (depth == 0) break; depth--; }
            else if (peek_token().kind == TokenKind::kComma && depth == 0) {
                std::string_view idx_text = trim_view(tokens_.source().substr(start_pos, peek_token().position - start_pos));
                indices.push_back(idx_text);
                // 创建完整的 AST 节点
                CommandASTNode idx_expr = CommandASTNode::make_expression(idx_text);
                compile_expression_info(*const_cast<ExpressionInfo*>(idx_expr.as_expression()));
                index_exprs.push_back(std::move(idx_expr));
                advance_token();
                start_pos = peek_token().position;
                continue;
            }
            advance_token();
        }
        std::string_view idx_text = trim_view(tokens_.source().substr(start_pos, peek_token().position - start_pos));
        indices.push_back(idx_text);
        // 创建完整的 AST 节点
        CommandASTNode idx_expr = CommandASTNode::make_expression(idx_text);
        compile_expression_info(*const_cast<ExpressionInfo*>(idx_expr.as_expression()));
        index_exprs.push_back(std::move(idx_expr));
        if (!match_token(TokenKind::kRBracket)) throw_syntax_error("unmatched '[' in index assignment");
    } else advance_token();
    if (!match_token(TokenKind::kEqual)) { restore_checkpoint(expression_checkpoint); return parse_expression(single_statement); }
    std::string_view value = collect_statement_expression();
    validate_expression_text(value, tokens_, value.data() - tokens_.source().data());
    // 使用新的带完整 AST 的工厂方法
    CommandASTNode value_expr = CommandASTNode::make_expression(value);
    compile_expression_info(*const_cast<ExpressionInfo*>(value_expr.as_expression()));
    return CommandASTNode::make_index_assignment_ast(id_tok.text, std::move(index_exprs), std::move(value_expr));
}

CommandASTNode CommandParser::parse_function_call(Token id_tok, bool single_statement, const LazyTokenStream::Checkpoint& expression_checkpoint) {
    expect_token(TokenKind::kLParen, "expected '(' in function call");
    std::vector<CommandASTNode> arguments;
    if (!check_token(TokenKind::kRParen)) {
        arguments = parse_argument_list_by_tokens(true);
        if (!match_token(TokenKind::kRParen)) throw_syntax_error("unmatched '(' in function call");
    } else advance_token();
    const TokenKind next = peek_token().kind;
    if (next != TokenKind::kEnd && next != TokenKind::kSemicolon) { restore_checkpoint(expression_checkpoint); return parse_expression(single_statement); }
    return CommandASTNode::make_function_call(id_tok.text, std::move(arguments));
}

std::vector<CommandASTNode> CommandParser::parse_argument_list_by_tokens(bool stop_at_rparen) {
    std::vector<CommandASTNode> arguments;
    if (peek_token().kind == TokenKind::kEnd) return arguments;
    while (peek_token().kind != TokenKind::kEnd) {
        if (stop_at_rparen && peek_token().kind == TokenKind::kRParen) break;
        arguments.push_back(parse_expression(false));
        if (match_token(TokenKind::kComma)) {
            if (stop_at_rparen && peek_token().kind == TokenKind::kRParen) { arguments.push_back(CommandASTNode::make_empty()); break; }
            continue;
        } else break;
    }
    return arguments;
}

CommandASTNode CommandParser::parse_expression(bool single_statement) {
    (void)single_statement;
    if (peek_token().kind == TokenKind::kEnd) return CommandASTNode::make_expression("");
    std::string_view expr = collect_statement_expression();
    validate_expression_text(expr, tokens_, expr.data() - tokens_.source().data());
    ExpressionInfo info(expr);
    compile_expression_info(info);
    return CommandASTNode(CommandKind::kExpression, std::move(info));
}

std::string_view CommandParser::collect_statement_expression() {
    std::size_t start_pos = peek_token().position;
    int paren_depth = 0, bracket_depth = 0, brace_depth = 0;
    while (peek_token().kind != TokenKind::kEnd) {
        const Token& tok = peek_token();
        if (tok.kind == TokenKind::kSemicolon && paren_depth == 0 && bracket_depth == 0 && brace_depth == 0) break;
        if (tok.kind == TokenKind::kLParen) paren_depth++;
        else if (tok.kind == TokenKind::kRParen) {
            if (paren_depth == 0) break;  // 在顶层右括号处停止
            paren_depth--;
        }
        else if (tok.kind == TokenKind::kLBracket) bracket_depth++;
        else if (tok.kind == TokenKind::kRBracket) {
            if (bracket_depth == 0) break;  // 在顶层右方括号处停止
            bracket_depth--;
        }
        else if (tok.kind == TokenKind::kLBrace) brace_depth++;
        else if (tok.kind == TokenKind::kRBrace) {
            if (brace_depth == 0) break;  // 在顶层右花括号处停止
            brace_depth--;
        }
        else if (tok.kind == TokenKind::kComma && paren_depth == 0 && bracket_depth == 0 && brace_depth == 0) break;  // 在顶层逗号处停止
        advance_token();
    }
    return trim_view(tokens_.source().substr(start_pos, peek_token().position - start_pos));
}

CommandASTNode parse_command(std::string_view source, const CommandParser::CommandConfig& config) {
    auto source_owner = std::make_shared<const std::string>(source);
    CommandParser parser(*source_owner, config);
    CommandASTNode node = parser.parse();
    attach_source_owner(node, source_owner);
    return node;
}

CommandASTNode parse_command(std::string_view source, CommandParser::IsCommandCallback is_command) {
    auto source_owner = std::make_shared<const std::string>(source);
    CommandParser parser(*source_owner, std::move(is_command));
    CommandASTNode node = parser.parse();
    attach_source_owner(node, source_owner);
    return node;
}

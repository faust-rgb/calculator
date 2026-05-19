// ============================================================================
// 命令解析器 - 统一的命令词法分析与语法树构建
// ============================================================================
//
// 本模块实现统一的命令解析框架：
//
// 1. 命令 AST 节点 (CommandASTNode)
//    - 分类命令类型：元命令、函数定义、函数调用、表达式
//    - 存储解析后的结构化信息，避免重复字符串操作
//    - 支持表达式预编译（ExpressionInfo）
//
// 2. 命令解析器 (CommandParser)
//    - 基于 LazyTokenStream 的惰性词法分析
//    - 统一的语法错误处理
//    - 检查点机制用于高效回溯
//
// 设计目标：
// - 按需生成 Token，避免全量词法分析开销
// - O(N) 单次扫描解析
// - 统一的异常处理策略
// - 精确的位置信息用于错误报告
// - 支持进制前缀和虚数后缀
// ============================================================================

#ifndef COMMAND_PARSER_H
#define COMMAND_PARSER_H

#include "core/common/calculator_exceptions.h"
#include "parser/lexer/lazy_token_stream.h"
#include "parser/lexer/token_types.h"
#include "parser/ast/expression_compiler.h"
#include "parser/ast/expression_ast.h"
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <variant>

// 前向声明
class ICommandRegistry;

// ============================================================================
// 命令 AST 节点类型
// ============================================================================

/**
 * @enum CommandKind
 * @brief 命令 AST 节点的类型分类
 */
enum class CommandKind {
    kEmpty,              ///< 空输入
    kMetaCommand,        ///< 元命令，如 :help, :vars
    kFunctionDefinition, ///< 函数定义，如 f(x,y) = x^2+y^2
    kFunctionCall,       ///< 函数调用，如 sin(pi/4)
    kAssignment,         ///< 变量赋值，如 x = 5
    kIndexAssignment,    ///< 索引赋值，如 x[i] = value 或 m[r, c] = value
    kExpression,         ///< 纯表达式，如 2 + 3
    kStringLiteral,      ///< 字符串字面量，如 "hello"
    kSequence            ///< 语句序列，如 x=1; y=2; x+y
};

/**
 * @struct ExpressionInfo
 * @brief 表达式信息（支持预编译）
 *
 * 注意：此结构正在逐步被废弃。新设计使用 CommandASTNode 子节点
 * 来表示表达式，实现完全解析。保留此结构仅用于向后兼容。
 */
struct ExpressionInfo {
    std::string_view text;                    ///< 原始文本视图
    mutable std::shared_ptr<ExpressionCache> cache;  ///< 预编译缓存

    ExpressionInfo() = default;
    explicit ExpressionInfo(std::string_view t) : text(t) {}

    ExpressionCache* get_or_create_cache() const {
        if (!cache) cache = std::make_shared<ExpressionCache>(text);
        return cache.get();
    }

    bool has_cache() const { return cache != nullptr; }
};

// 前向声明
class CommandASTNode;

/**
 * @struct FunctionDefinitionInfo
 * @brief 函数定义信息
 *
 * 新设计：body_expr 存储完全解析的表达式 AST 节点，
 * body 保留文本视图用于显示和向后兼容。
 */
struct FunctionDefinitionInfo {
    std::string_view name;
    std::vector<std::string_view> parameters;
    ExpressionInfo body;           ///< 向后兼容：原始文本
    std::unique_ptr<CommandASTNode> body_expr;  ///< 新增：完全解析的表达式 AST
};

/**
 * @struct AssignmentInfo
 * @brief 变量赋值信息
 *
 * 新设计：value_expr 存储完全解析的表达式 AST 节点，
 * expression 保留文本视图用于显示和向后兼容。
 */
struct AssignmentInfo {
    std::string_view variable;
    ExpressionInfo expression;     ///< 向后兼容：原始文本
    std::unique_ptr<CommandASTNode> value_expr;  ///< 新增：完全解析的表达式 AST
};

/**
 * @struct IndexAssignmentInfo
 * @brief 索引赋值信息
 *
 * 新设计：index_exprs 和 value_expr 存储完全解析的 AST 节点，
 * indices 和 value 保留文本视图用于显示和向后兼容。
 */
struct IndexAssignmentInfo {
    std::string_view variable;
    std::vector<ExpressionInfo> indices;  ///< 向后兼容：原始文本
    ExpressionInfo value;                 ///< 向后兼容：原始文本
    std::vector<std::unique_ptr<CommandASTNode>> index_exprs;  ///< 新增：完全解析的索引 AST
    std::unique_ptr<CommandASTNode> value_expr;                ///< 新增：完全解析的值 AST
};

/**
 * @struct MetaCommandInfo
 * @brief 元命令信息
 */
struct MetaCommandInfo {
    std::string_view command;
    std::vector<std::unique_ptr<CommandASTNode>> arguments;
};

/**
 * @struct NamedArgument
 * @brief 命名参数信息
 */
struct NamedArgument {
    std::string_view name;
    std::unique_ptr<CommandASTNode> value;
};

/**
 * @struct FunctionCallInfo
 * @brief 函数调用信息
 */
struct FunctionCallInfo {
    std::string_view name;
    std::vector<std::unique_ptr<CommandASTNode>> arguments;
    std::vector<NamedArgument> named_args;
};

/**
 * @class CommandASTNode
 * @brief 命令抽象语法树节点
 */
class CommandASTNode {
public:
    CommandKind kind;
    std::shared_ptr<const std::string> source_owner;

    std::variant<
        std::monostate,              // kEmpty
        MetaCommandInfo,             // kMetaCommand
        FunctionDefinitionInfo,      // kFunctionDefinition
        FunctionCallInfo,            // kFunctionCall
        AssignmentInfo,              // kAssignment
        IndexAssignmentInfo,         // kIndexAssignment
        ExpressionInfo,              // kExpression
        std::string,                 // kStringLiteral
        std::vector<CommandASTNode>  // kSequence
    > data;

    CommandASTNode() : kind(CommandKind::kEmpty), data(std::monostate{}) {}

    // 移动支持
    CommandASTNode(CommandASTNode&&) = default;
    CommandASTNode& operator=(CommandASTNode&&) = default;

    // 禁用拷贝
    CommandASTNode(const CommandASTNode&) = delete;
    CommandASTNode& operator=(const CommandASTNode&) = delete;

    CommandASTNode(CommandKind k, ExpressionInfo info)
        : kind(k), data(std::move(info)) {}

    static CommandASTNode make_empty();
    static CommandASTNode make_meta_command(std::string_view cmd,
                                            std::vector<CommandASTNode> args);
    static CommandASTNode make_function_definition(std::string_view name,
                                                   const std::vector<std::string_view>& params,
                                                   std::string_view body);
    /// 新增：带完整 AST 的函数定义
    static CommandASTNode make_function_definition_ast(std::string_view name,
                                                       const std::vector<std::string_view>& params,
                                                       CommandASTNode body_expr);
    static CommandASTNode make_function_call(std::string_view name,
                                             std::vector<CommandASTNode> args);
    static CommandASTNode make_function_call_with_named(
        std::string_view name,
        std::vector<CommandASTNode> positional_args,
        std::vector<std::pair<std::string_view, CommandASTNode>> named_args);
    static CommandASTNode make_assignment(std::string_view var,
                                          std::string_view expr);
    /// 新增：带完整 AST 的赋值
    static CommandASTNode make_assignment_ast(std::string_view var,
                                              CommandASTNode value_expr);
    static CommandASTNode make_index_assignment(
        std::string_view var,
        const std::vector<std::string_view>& indices,
        std::string_view value);
    /// 新增：带完整 AST 的索引赋值
    static CommandASTNode make_index_assignment_ast(
        std::string_view var,
        std::vector<CommandASTNode> index_exprs,
        CommandASTNode value_expr);
    static CommandASTNode make_expression(std::string_view expr);
    static CommandASTNode make_string_literal(const std::string& value);
    static CommandASTNode make_sequence(std::vector<CommandASTNode> nodes);

    const MetaCommandInfo* as_meta_command() const {
        return kind == CommandKind::kMetaCommand ? &std::get<MetaCommandInfo>(data) : nullptr;
    }
    const FunctionDefinitionInfo* as_function_definition() const {
        return kind == CommandKind::kFunctionDefinition ? &std::get<FunctionDefinitionInfo>(data) : nullptr;
    }
    const FunctionCallInfo* as_function_call() const {
        return kind == CommandKind::kFunctionCall ? &std::get<FunctionCallInfo>(data) : nullptr;
    }
    const AssignmentInfo* as_assignment() const {
        return kind == CommandKind::kAssignment ? &std::get<AssignmentInfo>(data) : nullptr;
    }
    const IndexAssignmentInfo* as_index_assignment() const {
        return kind == CommandKind::kIndexAssignment ? &std::get<IndexAssignmentInfo>(data) : nullptr;
    }
    const ExpressionInfo* as_expression() const {
        return kind == CommandKind::kExpression ? &std::get<ExpressionInfo>(data) : nullptr;
    }
    const std::string* as_string_literal() const {
        return kind == CommandKind::kStringLiteral ? &std::get<std::string>(data) : nullptr;
    }
    const std::vector<CommandASTNode>* as_sequence() const {
        return kind == CommandKind::kSequence ? &std::get<std::vector<CommandASTNode>>(data) : nullptr;
    }
};

/**
 * @class CommandParser
 * @brief 统一的命令解析器
 *
 * 新设计：解析器不再依赖外部回调来判断命令。
 * 使用内置的命令配置表来驱动解析决策，实现解析器与注册表的解耦。
 */
class CommandParser {
public:
    /// @deprecated 使用 CommandConfig 替代
    using IsCommandCallback = std::function<bool(std::string_view)>;

    /**
     * @struct CommandConfig
     * @brief 命令配置表，用于驱动解析决策
     *
     * 解耦设计：解析器使用此配置表来判断标识符是否可能是命令，
     * 而不直接依赖 CommandRegistry。
     */
    struct CommandConfig {
        std::vector<std::string> exact_commands;    ///< 精确匹配的命令名
        std::vector<std::string> prefix_commands;   ///< 前缀匹配的命令名

        /// 检查标识符是否可能是命令
        bool could_be_command(std::string_view name) const {
            for (const auto& cmd : exact_commands) {
                if (name == cmd) return true;
            }
            for (const auto& prefix : prefix_commands) {
                if (name.size() >= prefix.size() &&
                    name.substr(0, prefix.size()) == prefix) return true;
            }
            return false;
        }

        /// 从 CommandRegistry 构建配置
        static CommandConfig from_registry(const ICommandRegistry& registry);
    };

    /// 使用配置表构造（推荐）
    explicit CommandParser(std::string_view source,
                           const CommandConfig& config = {});

    /// 使用回调构造（向后兼容，已废弃）
    explicit CommandParser(std::string_view source, IsCommandCallback is_command);

    CommandASTNode parse();
    std::size_t current_position() const;
    std::string_view source() const { return tokens_.source(); }

private:
    const Token& peek_token();
    const Token& peek_token(std::size_t offset);
    Token advance_token();
    bool check_token(TokenKind kind) const;
    bool match_token(TokenKind kind);
    Token expect_token(TokenKind kind, const char* message);

    LazyTokenStream::Checkpoint save_checkpoint();
    void restore_checkpoint(const LazyTokenStream::Checkpoint& cp);

    CommandASTNode parse_command(bool single_statement = false);
    CommandASTNode parse_meta_command();
    CommandASTNode parse_definition_or_assignment(Token id_token, bool single_statement);
    CommandASTNode parse_function_call(Token id_token,
                                       bool single_statement,
                                       const LazyTokenStream::Checkpoint& expression_checkpoint);
    CommandASTNode parse_index_assignment(Token id_token,
                                          bool single_statement,
                                          const LazyTokenStream::Checkpoint& expression_checkpoint);

    std::vector<CommandASTNode> parse_argument_list_by_tokens(bool stop_at_rparen);
    CommandASTNode parse_expression(bool single_statement);
    std::string_view collect_statement_expression();
    [[noreturn]] void throw_syntax_error(const std::string& message);

    LazyTokenStream tokens_;
    CommandConfig config_;           ///< 新增：命令配置表
    IsCommandCallback is_command_;   ///< 向后兼容：回调函数（已废弃）
};

/// 使用配置表解析命令（推荐）
CommandASTNode parse_command(std::string_view source,
                             const CommandParser::CommandConfig& config);

/// 使用回调解析命令（向后兼容）
CommandASTNode parse_command(std::string_view source,
                             CommandParser::IsCommandCallback is_command);

#endif // COMMAND_PARSER_H

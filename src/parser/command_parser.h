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

#include "core/calculator_exceptions.h"
#include "parser/lazy_token_stream.h"
#include "parser/token_types.h"
#include "parser/expression_compiler.h"
#include "parser/expression_ast.h"
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <variant>

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
    kExpression,         ///< 纯表达式，如 2 + 3
    kStringLiteral,      ///< 字符串字面量，如 "hello"
    kSequence            ///< 语句序列，如 x=1; y=2; x+y
};

/**
 * @struct ExpressionInfo
 * @brief 表达式信息（支持预编译）
 *
 * 优化设计：
 * - text 是原始文本视图，不拥有内存
 * - cache 是延迟分配的，只有需要时才创建
 * - 避免重复存储文本
 */
struct ExpressionInfo {
    std::string_view text;                    ///< 原始文本视图
    mutable std::shared_ptr<ExpressionCache> cache;  ///< 预编译缓存（延迟分配）

    ExpressionInfo() = default;
    explicit ExpressionInfo(std::string_view t) : text(t) {}

    // 获取缓存（延迟创建）
    ExpressionCache* get_or_create_cache() const {
        if (!cache) {
            cache = std::make_shared<ExpressionCache>(text);
        }
        return cache.get();
    }

    // 检查是否有缓存
    bool has_cache() const { return cache != nullptr; }
};

/**
 * @struct FunctionDefinitionInfo
 * @brief 函数定义信息
 */
struct FunctionDefinitionInfo {
    std::string_view name;                      ///< 函数名
    std::vector<std::string_view> parameters;   ///< 参数名列表
    ExpressionInfo body;                        ///< 函数体表达式（支持预编译）
};

/**
 * @struct AssignmentInfo
 * @brief 变量赋值信息
 */
struct AssignmentInfo {
    std::string_view variable;       ///< 变量名
    ExpressionInfo expression;       ///< 赋值表达式（支持预编译）
};

/**
 * @struct MetaCommandInfo
 * @brief 元命令信息
 */
struct MetaCommandInfo {
    std::string_view command;                    ///< 命令名（不含冒号）
    std::vector<std::string_view> arguments;     ///< 参数列表
};

/**
 * @struct FunctionCallInfo
 * @brief 函数调用信息
 */
struct FunctionCallInfo {
    std::string_view name;                    ///< 函数名
    std::vector<ExpressionInfo> arguments;    ///< 参数列表（支持预编译）
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
        ExpressionInfo,              // kExpression (改为 ExpressionInfo 支持预编译)
        std::string,                 // kStringLiteral
        std::vector<CommandASTNode>  // kSequence
    > data;

    CommandASTNode() : kind(CommandKind::kEmpty), data(std::monostate{}) {}

    // 直接接受 ExpressionInfo 的构造函数（支持预编译）
    CommandASTNode(CommandKind k, ExpressionInfo info)
        : kind(k), data(std::move(info)) {}

    static CommandASTNode make_empty();
    static CommandASTNode make_meta_command(std::string_view cmd,
                                            const std::vector<std::string_view>& args);
    static CommandASTNode make_function_definition(std::string_view name,
                                                   const std::vector<std::string_view>& params,
                                                   std::string_view body);
    static CommandASTNode make_function_call(std::string_view name,
                                             const std::vector<std::string_view>& args);
    static CommandASTNode make_assignment(std::string_view var,
                                          std::string_view expr);
    static CommandASTNode make_expression(std::string_view expr);
    static CommandASTNode make_string_literal(const std::string& value);
    static CommandASTNode make_sequence(std::vector<CommandASTNode> nodes);

    // 类型安全访问器
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

// ============================================================================
// 命令解析器
// ============================================================================

/**
 * @class CommandParser
 * @brief 统一的命令解析器
 *
 * 使用 LazyTokenStream 实现按需词法分析。
 * 支持检查点机制用于高效回溯。
 *
 * 特性：
 * - 支持进制前缀（0x, 0b, 0o）
 * - 支持虚数后缀（i）
 * - 保留数字原始文本
 * - 基于 Token 的参数解析
 * - 检查点回溯机制
 */
class CommandParser {
public:
    /// 回调类型：判断一个标识符是否为已注册的命令
    using IsCommandCallback = std::function<bool(std::string_view)>;

    explicit CommandParser(std::string_view source, IsCommandCallback is_command = nullptr);

    /**
     * @brief 解析命令，构建 AST
     * @return 命令 AST 节点
     */
    CommandASTNode parse();

    /**
     * @brief 获取当前解析位置（用于错误报告）
     */
    std::size_t current_position() const;

    /**
     * @brief 获取源字符串（用于错误报告）
     */
    std::string_view source() const { return tokens_.source(); }

private:
    // ========================================================================
    // Token 访问
    // ========================================================================

    /// 查看当前 Token（不消费）
    const Token& peek_token();

    /// 查看指定偏移量的 Token
    const Token& peek_token(std::size_t offset);

    /// 消费并返回当前 Token
    Token advance_token();

    /// 检查当前 Token 类型
    bool check_token(TokenKind kind) const;

    /// 匹配并消费指定类型的 Token
    bool match_token(TokenKind kind);

    /// 期望指定类型的 Token，失败则抛出异常
    Token expect_token(TokenKind kind, const char* message);

    // ========================================================================
    // 回溯支持
    // ========================================================================

    /// 保存当前解析状态
    LazyTokenStream::Checkpoint save_checkpoint();

    /// 恢复到指定解析状态
    void restore_checkpoint(const LazyTokenStream::Checkpoint& cp);

    // ========================================================================
    // 语法分析
    // ========================================================================

    /// 解析命令
    CommandASTNode parse_command(bool single_statement = false);

    /// 解析元命令
    CommandASTNode parse_meta_command();

    /// 解析函数定义或赋值
    CommandASTNode parse_definition_or_assignment(Token id_token, bool single_statement);

    /// 解析函数调用
    CommandASTNode parse_function_call(Token id_token,
                                       bool single_statement,
                                       const LazyTokenStream::Checkpoint& expression_checkpoint);

    /// 基于 Token 流解析参数列表（不使用字符串切分）
    std::vector<std::string_view> parse_argument_list_by_tokens(bool stop_at_rparen);

    /// 解析表达式（兜底）
    CommandASTNode parse_expression(bool single_statement);

    /// 收集当前语句尾部的表达式片段，停在顶层分号或输入结束
    std::string_view collect_statement_expression();

    /// 抛出带有上下文信息的语法错误
    [[noreturn]] void throw_syntax_error(const std::string& message);

    // ========================================================================
    // 成员变量
    // ========================================================================

    LazyTokenStream tokens_;
    IsCommandCallback is_command_;
};

// ============================================================================
// 便捷函数
// ============================================================================

/**
 * @brief 解析命令字符串，返回 AST
 */
CommandASTNode parse_command(std::string_view source,
                             CommandParser::IsCommandCallback is_command = nullptr);

#endif // COMMAND_PARSER_H

// ============================================================================
// 命令解析器 - 统一的命令词法分析与语法树构建
// ============================================================================
//
// 本模块实现统一的命令解析框架：
//
// 1. 命令 Token 类型 (CommandToken)
//    - 标识符、字符串字面量、运算符、括号等
//    - 单次扫描生成 Token 流，避免重复遍历
//
// 2. 命令 AST 节点 (CommandASTNode)
//    - 分类命令类型：元命令、函数定义、函数调用、表达式
//    - 存储解析后的结构化信息，避免重复字符串操作
//
// 3. 命令解析器 (CommandParser)
//    - 基于 Token 流的递归下降解析
//    - 统一的语法错误处理
//
// 设计目标：
// - O(N) 单次扫描解析
// - 统一的异常处理策略
// - 精确的位置信息用于错误报告
// ============================================================================

#ifndef COMMAND_PARSER_H
#define COMMAND_PARSER_H

#include "calculator_exceptions.h"
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <variant>

// ============================================================================
// 命令 Token 类型
// ============================================================================

/**
 * @enum CommandTokenType
 * @brief 命令 Token 类型
 */
enum class CommandTokenType {
    kEnd,           ///< 输入结束
    kIdentifier,    ///< 标识符（变量名、函数名）
    kString,        ///< 字符串字面量（已解码转义）
    kNumber,        ///< 数字字面量
    kLParen,        ///< 左圆括号 (
    kRParen,        ///< 右圆括号 )
    kLBracket,      ///< 左方括号 [
    kRBracket,      ///< 右方括号 ]
    kComma,         ///< 逗号 ,
    kEqual,         ///< 等号 =（赋值或定义）
    kColon,         ///< 冒号 :
    kOperator,      ///< 运算符 + - * / ^ 等
};

/**
 * @struct CommandToken
 * @brief 命令 Token 结构
 */
struct CommandToken {
    CommandTokenType type;
    std::string_view text;      ///< 原始文本视图
    std::string string_value;   ///< 字符串值（仅 kString 类型使用）
    double number_value = 0.0;  ///< 数值（仅 kNumber 类型使用）
    std::size_t position = 0;   ///< 在源字符串中的起始位置

    CommandToken() : type(CommandTokenType::kEnd) {}
};

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
    kFunctionDefinition, ///< 函数定义，如 f(x) = x^2
    kFunctionCall,       ///< 函数调用，如 sin(pi/4)
    kAssignment,         ///< 变量赋值，如 x = 5
    kExpression,         ///< 纯表达式，如 2 + 3
    kStringLiteral       ///< 字符串字面量，如 "hello"
};

/**
 * @struct FunctionDefinitionInfo
 * @brief 函数定义信息
 */
struct FunctionDefinitionInfo {
    std::string_view name;           ///< 函数名
    std::string_view parameter;      ///< 参数名
    std::string_view body;           ///< 函数体表达式
};

/**
 * @struct FunctionCallInfo
 * @brief 函数调用信息
 */
struct FunctionCallInfo {
    std::string_view name;                    ///< 函数名
    std::vector<std::string_view> arguments;  ///< 参数列表（原始字符串形式）
};

/**
 * @struct AssignmentInfo
 * @brief 变量赋值信息
 */
struct AssignmentInfo {
    std::string_view variable;       ///< 变量名
    std::string_view expression;     ///< 赋值表达式
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
 * @class CommandASTNode
 * @brief 命令抽象语法树节点
 */
class CommandASTNode {
public:
    CommandKind kind;

    std::variant<
        std::monostate,              // kEmpty
        MetaCommandInfo,             // kMetaCommand
        FunctionDefinitionInfo,      // kFunctionDefinition
        FunctionCallInfo,            // kFunctionCall
        AssignmentInfo,              // kAssignment
        std::string_view,            // kExpression
        std::string                  // kStringLiteral
    > data;

    CommandASTNode() : kind(CommandKind::kEmpty), data(std::monostate{}) {}

    static CommandASTNode make_empty();
    static CommandASTNode make_meta_command(std::string_view cmd,
                                            const std::vector<std::string_view>& args);
    static CommandASTNode make_function_definition(std::string_view name,
                                                   std::string_view param,
                                                   std::string_view body);
    static CommandASTNode make_function_call(std::string_view name,
                                             const std::vector<std::string_view>& args);
    static CommandASTNode make_assignment(std::string_view var,
                                          std::string_view expr);
    static CommandASTNode make_expression(std::string_view expr);
    static CommandASTNode make_string_literal(const std::string& value);

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

    const std::string_view* as_expression() const {
        return kind == CommandKind::kExpression ? &std::get<std::string_view>(data) : nullptr;
    }

    const std::string* as_string_literal() const {
        return kind == CommandKind::kStringLiteral ? &std::get<std::string>(data) : nullptr;
    }
};

// ============================================================================
// 命令解析器
// ============================================================================

/**
 * @class CommandParser
 * @brief 统一的命令解析器
 *
 * 基于 Token 流实现 O(N) 单次扫描解析：
 * 1. 词法分析生成 Token 流
 * 2. 递归下降解析构建 AST
 * 3. 统一的语法错误处理
 */
class CommandParser {
public:
    explicit CommandParser(std::string_view source);

    /**
     * @brief 解析命令，构建 AST
     * @return 命令 AST 节点
     */
    CommandASTNode parse();

    /**
     * @brief 获取当前解析位置（用于错误报告）
     */
    std::size_t current_position() const { return pos_; }

    /**
     * @brief 获取源字符串（用于错误报告）
     */
    std::string_view source() const { return source_; }

private:
    // ========================================================================
    // 词法分析
    // ========================================================================

    /// 判断字符是否为空白
    static constexpr bool is_space_char(char ch) {
        return ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r';
    }

    /// 跳过空白字符
    void skip_spaces();

    /// 查看当前字符
    char peek() const;

    /// 查看下一个字符（不消费）
    char peek_next() const;

    /// 消费当前字符
    char advance();

    /// 是否到达末尾
    bool is_at_end() const;

    /// 解析下一个 Token
    CommandToken next_token();

    /// 解析标识符
    CommandToken parse_identifier_token();

    /// 解析数字
    CommandToken parse_number_token();

    /// 解析字符串字面量
    CommandToken parse_string_token();

    /// 解析运算符
    CommandToken parse_operator_token();

    // ========================================================================
    // 语法分析
    // ========================================================================

    /// 查看当前 Token（不消费）
    const CommandToken& peek_token();

    /// 消费并返回当前 Token
    CommandToken advance_token();

    /// 检查当前 Token 类型
    bool check_token(CommandTokenType type) const;

    /// 匹配并消费指定类型的 Token
    bool match_token(CommandTokenType type);

    /// 期望指定类型的 Token，失败则抛出异常
    CommandToken expect_token(CommandTokenType type, const char* message);

    /// 解析命令
    CommandASTNode parse_command();

    /// 解析元命令
    CommandASTNode parse_meta_command();

    /// 解析函数定义或赋值
    CommandASTNode parse_definition_or_assignment(const CommandToken& id_token);

    /// 解析函数调用
    CommandASTNode parse_function_call(const CommandToken& id_token);

    /// 解析参数列表
    std::vector<std::string_view> parse_argument_list();

    /// 解析表达式（兜底）
    CommandASTNode parse_expression();

    /// 抛出语法错误
    [[noreturn]] void throw_syntax_error(const std::string& message);

    // ========================================================================
    // 成员变量
    // ========================================================================

    std::string_view source_;
    std::size_t pos_ = 0;
    std::vector<CommandToken> tokens_;
    std::size_t token_pos_ = 0;
    bool tokens_scanned_ = false;
};

// ============================================================================
// 便捷函数
// ============================================================================

/**
 * @brief 解析命令字符串，返回 AST
 */
CommandASTNode parse_command(std::string_view source);

/**
 * @brief 拆分顶层参数（正确处理嵌套和字符串）
 */
std::vector<std::string_view> split_command_arguments(std::string_view text);

#endif // COMMAND_PARSER_H

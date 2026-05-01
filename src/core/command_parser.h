// ============================================================================
// 命令解析器 - 统一的命令词法分析与语法树构建
// ============================================================================
//
// 本模块实现统一的命令解析框架：
//
// 1. 命令词法分析器 (CommandTokenizer)
//    - 基于 BaseParser 衍生，复用成熟的词法分析工具
//    - 正确处理嵌套括号、字符串字面量、转义字符
//    - 支持元命令 (:cmd)、函数调用 (func(...))、表达式等
//
// 2. 命令 AST 节点 (CommandASTNode)
//    - 分类命令类型：元命令、函数定义、函数调用、表达式
//    - 存储解析后的结构化信息，避免重复字符串操作
//
// 3. 命令解析器 (CommandParser)
//    - 递归下降解析，构建 AST
//    - 统一处理命令边界判定，消除硬编码字符串切分
//
// 设计目标：
// - 废弃直接的字符串截取（find('('), substr 等）
// - 通过 AST 节点判定命令类型
// - 提高解析性能和正确性
// ============================================================================

#ifndef COMMAND_PARSER_H
#define COMMAND_PARSER_H

#include "base_parser.h"
#include <memory>
#include <string>
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
    std::string name;           ///< 函数名
    std::string parameter;      ///< 参数名
    std::string body;           ///< 函数体表达式
};

/**
 * @struct FunctionCallInfo
 * @brief 函数调用信息
 */
struct FunctionCallInfo {
    std::string name;                    ///< 函数名
    std::vector<std::string> arguments;  ///< 参数列表（原始字符串形式）
};

/**
 * @struct AssignmentInfo
 * @brief 变量赋值信息
 */
struct AssignmentInfo {
    std::string variable;       ///< 变量名
    std::string expression;     ///< 赋值表达式
};

/**
 * @struct MetaCommandInfo
 * @brief 元命令信息
 */
struct MetaCommandInfo {
    std::string command;                    ///< 命令名（不含冒号）
    std::vector<std::string> arguments;     ///< 参数列表
};

/**
 * @class CommandASTNode
 * @brief 命令抽象语法树节点
 *
 * 使用 std::variant 存储不同类型命令的具体信息，
 * 避免运行时类型检查和 void* 转换。
 */
class CommandASTNode {
public:
    CommandKind kind;

    // 使用 variant 存储具体信息
    std::variant<
        std::monostate,              // kEmpty
        MetaCommandInfo,             // kMetaCommand
        FunctionDefinitionInfo,      // kFunctionDefinition
        FunctionCallInfo,            // kFunctionCall
        AssignmentInfo,              // kAssignment
        std::string                  // kExpression, kStringLiteral
    > data;

    CommandASTNode() : kind(CommandKind::kEmpty), data(std::monostate{}) {}

    // 便捷构造函数
    static CommandASTNode make_empty();
    static CommandASTNode make_meta_command(const std::string& cmd,
                                            const std::vector<std::string>& args);
    static CommandASTNode make_function_definition(const std::string& name,
                                                   const std::string& param,
                                                   const std::string& body);
    static CommandASTNode make_function_call(const std::string& name,
                                             const std::vector<std::string>& args);
    static CommandASTNode make_assignment(const std::string& var,
                                          const std::string& expr);
    static CommandASTNode make_expression(const std::string& expr);
    static CommandASTNode make_string_literal(const std::string& value);

    // 便捷访问器
    const MetaCommandInfo* as_meta_command() const;
    const FunctionDefinitionInfo* as_function_definition() const;
    const FunctionCallInfo* as_function_call() const;
    const AssignmentInfo* as_assignment() const;
    const std::string* as_expression() const;
    const std::string* as_string_literal() const;
};

// ============================================================================
// 命令解析器
// ============================================================================

/**
 * @class CommandParser
 * @brief 统一的命令解析器
 *
 * 基于 BaseParser 实现，提供：
 * - 正确的嵌套括号处理
 * - 字符串字面量识别
 * - AST 构建
 *
 * 解析优先级：
 * 1. 空输入 → kEmpty
 * 2. 元命令 (:cmd ...) → kMetaCommand
 * 3. 函数定义 (f(x) = ...) → kFunctionDefinition
 * 4. 变量赋值 (x = ...) → kAssignment
 * 5. 函数调用 (func(...)) → kFunctionCall
 * 6. 纯表达式 → kExpression
 */
class CommandParser : public BaseParser {
public:
    explicit CommandParser(std::string source);

    /**
     * @brief 解析命令，构建 AST
     * @return 命令 AST 节点
     */
    CommandASTNode parse();

private:
    /**
     * @brief 解析函数调用或表达式
     * @return AST 节点
     */
    CommandASTNode parse_function_call_or_expression();

    /**
     * @brief 解析参数列表
     * @return 参数字符串列表
     */
    std::vector<std::string> parse_argument_list();

    /**
     * @brief 解析字符串字面量
     * @return 字符串内容（不含引号）
     */
    std::string parse_string_literal_content();
};

// ============================================================================
// 便捷函数
// ============================================================================

/**
 * @brief 解析命令字符串，返回 AST
 * @param source 命令字符串
 * @return 命令 AST 节点
 */
CommandASTNode parse_command(const std::string& source);

/**
 * @brief 拆分顶层参数（正确处理嵌套和字符串）
 * @param text 参数列表文本
 * @return 参数字符串列表
 */
std::vector<std::string> split_command_arguments(const std::string& text);

#endif // COMMAND_PARSER_H

#ifndef SCRIPT_AST_H
#define SCRIPT_AST_H

#include <memory>
#include <string>
#include <vector>

namespace script {

/**
 * @file script_ast.h
 * @brief 脚本语言的抽象语法树（AST）定义
 *
 * 定义了脚本语言的所有语句类型，形成树形结构。
 * 使用访问者模式的基础结构（通过基类虚析构）。
 */

/**
 * @struct Statement
 * @brief 所有语句类型的抽象基类
 */
struct Statement {
    /** @brief 语句类型枚举 */
    enum class Kind {
        kBlock,      ///< 代码块 { ... }
        kSimple,     ///< 简单语句（表达式）
        kIf,         ///< 条件语句 if-else
        kWhile,      ///< while 循环
        kFor,        ///< for 循环
        kFunction,   ///< 函数定义
        kReturn,     ///< 返回语句
        kBreak,      ///< 跳出循环
        kContinue,   ///< 继续下一次循环
    };

    /** @brief 构造函数，指定语句类型 */
    explicit Statement(Kind kind_value) : kind(kind_value) {}

    /** @brief 虚析构函数，支持多态删除 */
    virtual ~Statement() = default;

    Kind kind;  ///< 语句类型标识
};

/** @brief 语句的智能指针类型 */
using StatementPtr = std::unique_ptr<Statement>;

/**
 * @struct BlockStatement
 * @brief 代码块语句，包含多条子语句
 */
struct BlockStatement : Statement {
    BlockStatement() : Statement(Kind::kBlock) {}

    std::vector<StatementPtr> statements;  ///< 代码块内的语句列表
};

/**
 * @struct SimpleStatement
 * @brief 简单语句，通常是表达式
 */
struct SimpleStatement : Statement {
    SimpleStatement() : Statement(Kind::kSimple) {}

    std::string text;  ///< 语句文本（表达式字符串）
};

/**
 * @struct IfStatement
 * @brief 条件语句
 */
struct IfStatement : Statement {
    IfStatement() : Statement(Kind::kIf) {}

    std::string condition;    ///< 条件表达式
    StatementPtr then_branch; ///< then 分支
    StatementPtr else_branch; ///< else 分支（可选）
};

/**
 * @struct WhileStatement
 * @brief while 循环语句
 */
struct WhileStatement : Statement {
    WhileStatement() : Statement(Kind::kWhile) {}

    std::string condition;  ///< 循环条件
    StatementPtr body;      ///< 循环体
};

/**
 * @struct ForStatement
 * @brief for 循环语句
 */
struct ForStatement : Statement {
    ForStatement() : Statement(Kind::kFor) {}

    std::string initializer;  ///< 初始化表达式
    std::string condition;    ///< 循环条件
    std::string step;         ///< 步进表达式
    StatementPtr body;        ///< 循环体
};

/**
 * @struct FunctionStatement
 * @brief 函数定义语句
 */
struct FunctionStatement : Statement {
    FunctionStatement() : Statement(Kind::kFunction) {}

    std::string name;                  ///< 函数名
    std::vector<std::string> parameters;  ///< 参数名列表
    std::unique_ptr<BlockStatement> body; ///< 函数体
};

/**
 * @struct ReturnStatement
 * @brief 返回语句
 */
struct ReturnStatement : Statement {
    ReturnStatement() : Statement(Kind::kReturn) {}

    bool has_expression = false;  ///< 是否有返回值
    std::string expression;       ///< 返回表达式（可选）
};

/**
 * @struct BreakStatement
 * @brief break 语句，跳出循环
 */
struct BreakStatement : Statement {
    BreakStatement() : Statement(Kind::kBreak) {}
};

/**
 * @struct ContinueStatement
 * @brief continue 语句，跳过当前迭代
 */
struct ContinueStatement : Statement {
    ContinueStatement() : Statement(Kind::kContinue) {}
};

/**
 * @struct Program
 * @brief 完整的程序，包含顶层语句列表
 */
struct Program {
    std::vector<StatementPtr> statements;  ///< 顶层语句列表
};

}  // namespace script

#endif

#ifndef SCRIPT_AST_H
#define SCRIPT_AST_H

#include "parser/command_parser.h"

#include <memory>
#include <string>
#include <vector>

struct ExpressionCache;

namespace script {

struct Statement;
struct BlockStatement;

using StatementPtr = std::unique_ptr<Statement>;
using SharedBlockPtr = std::shared_ptr<const BlockStatement>;

/**
 * @struct Statement
 * @brief 所有语句类型的抽象基类
 */
struct Statement {
    enum class Kind {
        kBlock,
        kSimple,
        kIf,
        kWhile,
        kFor,
        kForRange,
        kForIn,
        kFunction,
        kReturn,
        kBreak,
        kContinue,
        kImport,
        kMatch,
    };

    explicit Statement(Kind kind_value) : kind(kind_value), line(0) {}
    virtual ~Statement() = default;

    Kind kind;
    int line;
};

/**
 * @struct BlockStatement
 * @brief 语句块
 */
struct BlockStatement : Statement {
    BlockStatement() : Statement(Kind::kBlock) {}

    std::vector<StatementPtr> statements;
};

/**
 * @struct SimpleStatement
 * @brief 简单语句，集成 CommandASTNode
 */
struct SimpleStatement : Statement {
    SimpleStatement() : Statement(Kind::kSimple) {}

    CommandASTNode command_ast;
    std::string text;
};

/**
 * @struct IfStatement
 * @brief 条件语句
 */
struct IfStatement : Statement {
    IfStatement() : Statement(Kind::kIf) {}

    CommandASTNode condition_ast;
    StatementPtr then_branch;
    StatementPtr else_branch;
};

/**
 * @struct WhileStatement
 * @brief while 循环
 */
struct WhileStatement : Statement {
    WhileStatement() : Statement(Kind::kWhile) {}

    CommandASTNode condition_ast;
    StatementPtr body;
};

/**
 * @struct ForStatement
 * @brief C 风格 for 循环
 */
struct ForStatement : Statement {
    ForStatement() : Statement(Kind::kFor) {}

    CommandASTNode init_ast;
    CommandASTNode cond_ast;
    CommandASTNode step_ast;
    StatementPtr body;
};

/**
 * @struct ForRangeStatement
 * @brief range 风格 for 循环
 */
struct ForRangeStatement : Statement {
    ForRangeStatement() : Statement(Kind::kForRange) {}

    std::string variable;
    CommandASTNode start_ast;
    CommandASTNode stop_ast;
    CommandASTNode step_ast;
    bool step_is_negative = false;
    StatementPtr body;
};

/**
 * @struct ForInStatement
 * @brief iterable 风格 for 循环
 */
struct ForInStatement : Statement {
    ForInStatement() : Statement(Kind::kForIn) {}

    std::string variable;
    CommandASTNode iterable_ast;
    StatementPtr body;
};

/**
 * @struct FunctionStatement
 * @brief 函数定义语句
 */
struct FunctionStatement : Statement {
    FunctionStatement() : Statement(Kind::kFunction) {}

    std::string name;
    std::vector<std::string> parameters;
    SharedBlockPtr body;
};

/**
 * @struct ReturnStatement
 * @brief 返回语句
 */
struct ReturnStatement : Statement {
    ReturnStatement() : Statement(Kind::kReturn) {}

    bool has_expression = false;
    CommandASTNode expr_ast;
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
 * @struct ImportStatement
 * @brief import 语句，执行另一个脚本文件
 */
struct ImportStatement : Statement {
    ImportStatement() : Statement(Kind::kImport) {}

    CommandASTNode path_ast;
    std::string path_text;
};

/**
 * @struct CaseClause
 * @brief match 语句的单个 case 分支
 */
struct CaseClause {
    CommandASTNode pattern_ast;      ///< 匹配模式（可以是值或条件）
    StatementPtr body;              ///< case 分支的执行体
    bool is_default = false;        ///< 是否是默认分支 (case _)
    bool is_guarded = false;        ///< 是否有守卫条件 (case x if x > 0)
    CommandASTNode guard_ast;       ///< 守卫条件 AST
};

/**
 * @struct MatchStatement
 * @brief match-case 语句，类似 Python 的模式匹配
 */
struct MatchStatement : Statement {
    MatchStatement() : Statement(Kind::kMatch) {}

    CommandASTNode subject_ast;     ///< 要匹配的表达式
    std::vector<CaseClause> cases;  ///< case 分支列表
};

/**
 * @struct Program
 * @brief 完整的程序，包含顶层语句列表
 */
struct Program {
    std::vector<StatementPtr> statements;
};

inline SharedBlockPtr make_shared_block(std::unique_ptr<BlockStatement> block) {
    return SharedBlockPtr(block.release());
}

}  // namespace script

#endif

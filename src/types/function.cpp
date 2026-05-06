// ============================================================================
// 用户函数类型实现
// ============================================================================

#include "function.h"
#include "parser/expression_ast.h"

std::shared_ptr<ExpressionAST> CustomFunction::get_or_compile_ast() const {
    if (ast_compiled) {
        return cached_ast;
    }

    // 编译表达式
    try {
        auto ast = compile_expression_ast(expression);
        if (ast) {
            cached_ast = std::move(ast);
        }
    } catch (...) {
        // 编译失败，保持 cached_ast 为空
    }

    ast_compiled = true;
    return cached_ast;
}

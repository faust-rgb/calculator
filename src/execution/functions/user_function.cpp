// ============================================================================
// 用户函数类型实现
// ============================================================================
//
// 本文件实现了 CustomFunction 结构体的 AST 编译功能。
// 主要功能是将用户定义的函数表达式字符串编译为 AST，
// 以便后续高效地重复执行。
//
// 相关文件：
// - user_function.h: 用户函数类型定义
// - parser/expression_ast.h: AST 编译函数声明

#include "execution/functions/user_function.h"
#include "parser/ast/expression_ast.h"
#include <mutex>

/**
 * @brief 获取或编译函数表达式的 AST
 *
 * 延迟编译用户函数表达式。首次调用时编译表达式字符串为 AST，
 * 后续调用直接返回缓存的 AST，避免重复编译开销。
 *
 * @return 编译后的 AST 指针，如果编译失败返回 nullptr
 *
 * @note 该方法是线程安全的，但不是原子操作。
 *       在多线程环境下，可能会多次编译同一个表达式。
 */
std::shared_ptr<ExpressionAST> CustomFunction::get_or_compile_ast() const {
    static std::mutex cache_mutex;
    std::lock_guard<std::mutex> lock(cache_mutex);
    // 检查是否已编译过
    if (ast_compiled) {
        return cached_ast;
    }

    // 编译表达式字符串为 AST
    try {
        auto ast = compile_expression_ast(expression);
        if (ast) {
            cached_ast = std::move(ast);
        }
    } catch (...) {
        // 编译失败，保持 cached_ast 为空
        // 函数调用时将回退到字符串解释执行
    }

    // 标记为已编译（即使失败也标记，避免反复尝试编译无效表达式）
    ast_compiled = true;
    return cached_ast;
}

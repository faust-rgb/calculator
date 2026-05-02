// ============================================================================
// 统一表达式解析器
// ============================================================================
//
// 整合所有表达式解析功能：
// 1. 特征分析（使用 UnifiedParserFactory）
// 2. 十进制解析
// 3. 精确有理数解析
// 4. 矩阵表达式解析
//
// 设计目标：
// - 单一入口点进行表达式分析
// - 避免重复的特征检测
// - 统一的错误处理
// - 支持缓存机制
//

#ifndef PARSER_UNIFIED_EXPRESSION_PARSER_H
#define PARSER_UNIFIED_EXPRESSION_PARSER_H

#include "parser/unified_parser_factory.h"
#include "precise/rational.h"
#include "types/stored_value.h"
#include "matrix.h"
#include <string>
#include <map>
#include <functional>
#include <memory>

// 前向声明
class VariableResolver;
struct CustomFunction;
struct ExpressionAST;
struct ExpressionCache;

// ============================================================================
// 统一表达式解析器
// ============================================================================

/**
 * @class UnifiedExpressionParser
 * @brief 统一的表达式解析和求值入口
 *
 * 使用方法：
 * 1. 创建解析器实例，传入必要的上下文
 * 2. 调用 analyze() 分析表达式特征
 * 3. 调用 evaluate() 或 evaluate_stored() 求值
 */
class UnifiedExpressionParser {
public:
    using ScalarFunction = std::function<double(const std::vector<double>&)>;
    using MatrixFunction = std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>;
    using HasScriptFunctionCallback = std::function<bool(const std::string&)>;
    using InvokeScriptFunctionCallback = std::function<double(const std::string&, const std::vector<double>&)>;

    /**
     * @brief 构造解析器
     */
    UnifiedExpressionParser(
        const VariableResolver& variables,
        const std::map<std::string, CustomFunction>* functions = nullptr,
        const std::map<std::string, ScalarFunction>* scalar_functions = nullptr,
        const std::map<std::string, MatrixFunction>* matrix_functions = nullptr,
        const std::map<std::string, matrix::ValueFunction>* value_functions = nullptr,
        HasScriptFunctionCallback has_script_function = {},
        InvokeScriptFunctionCallback invoke_script_function = {});

    // ========================================================================
    // 特征分析
    // ========================================================================

    /**
     * @brief 分析表达式特征
     */
    UnifiedParserFactory::AnalysisResult analyze(const std::string& expression);

    /**
     * @brief 判断表达式是否可以编译为 AST
     */
    bool can_compile_to_ast(const std::string& expression);

    /**
     * @brief 获取表达式类型提示
     */
    ExpressionHint get_hint(const std::string& expression);

    // ========================================================================
    // 求值
    // ========================================================================

    /**
     * @brief 求值表达式（十进制）
     */
    double evaluate(const std::string& expression);

    /**
     * @brief 求值表达式（精确有理数）
     */
    Rational evaluate_exact(const std::string& expression);

    /**
     * @brief 尝试求值可能包含矩阵的表达式
     */
    bool try_evaluate_value(const std::string& expression, matrix::Value* value);

    /**
     * @brief 求值表达式并返回 StoredValue（完整功能）
     * @param expression 表达式字符串
     * @param exact_mode 是否精确模式
     * @param symbolic_mode 是否符号模式
     * @param cache 表达式缓存（可选）
     * @return StoredValue 结果
     */
    StoredValue evaluate_stored(const std::string& expression,
                                bool exact_mode = false,
                                bool symbolic_mode = false,
                                ExpressionCache* cache = nullptr);

    /**
     * @brief 求值已编译的 AST
     */
    double evaluate_ast(const ExpressionAST* ast);

    /**
     * @brief 求值已编译的 AST（精确模式）
     */
    Rational evaluate_ast_exact(const ExpressionAST* ast);

    // ========================================================================
    // 编译
    // ========================================================================

    /**
     * @brief 编译表达式为 AST
     */
    std::unique_ptr<ExpressionAST> compile(const std::string& expression);

private:
    const VariableResolver& variables_;
    const std::map<std::string, CustomFunction>* functions_;
    const std::map<std::string, ScalarFunction>* scalar_functions_;
    const std::map<std::string, MatrixFunction>* matrix_functions_;
    const std::map<std::string, matrix::ValueFunction>* value_functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionCallback invoke_script_function_;

    UnifiedParserFactory factory_;
};

// ============================================================================
// 便捷函数
// ============================================================================

/**
 * @brief 快速求值十进制表达式
 */
double evaluate_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions = nullptr,
    const std::map<std::string, std::function<double(const std::vector<double>&)>>* scalar_functions = nullptr);

/**
 * @brief 快速求值精确有理数表达式
 */
Rational evaluate_expression_exact(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions = nullptr);

#endif // PARSER_UNIFIED_EXPRESSION_PARSER_H
// ============================================================================
// 编译表达式 AST
// ============================================================================
//
// 提供表达式的 AST 表示，用于：
// 1. 一次编译，多次求值（解决循环中重复解析问题）
// 2. 编译期变量绑定优化
// 3. 直接执行，跳过字符串解析开销
// ============================================================================

#ifndef EXPRESSION_AST_H
#define EXPRESSION_AST_H

#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <map>

// 前向声明
class VariableResolver;
struct CustomFunction;
struct StoredValue;

// 回调类型定义
using HasScriptFunctionCallback = std::function<bool(const std::string&)>;
using InvokeScriptFunctionDecimalCallback = std::function<double(const std::string&, const std::vector<double>&)>;

/**
 * @enum ExprKind
 * @brief 表达式 AST 节点类型
 */
enum class ExprKind {
    kNumber,       ///< 数字字面量
    kVariable,     ///< 变量引用
    kBinaryOp,     ///< 二元运算 (+, -, *, /, ^)
    kUnaryOp,      ///< 一元运算 (-, +)
    kFunctionCall, ///< 函数调用
    kComparison,   ///< 比较运算 (==, !=, <, >, <=, >=)
    kConditional,  ///< 三元条件 (condition ? then : else)
};

/**
 * @struct ExpressionAST
 * @brief 编译后的表达式 AST 节点
 */
struct ExpressionAST {
    ExprKind kind;

    // 数值（kNumber）
    double number_value = 0.0;
    std::string string_value; // 用于精确解析数字的原始文本

    // 标识符/函数名（kVariable, kFunctionCall）
    std::string identifier;

    // 源字符串位置，用于求值时精准报错
    std::size_t position = 0;

    // 子节点（kBinaryOp, kUnaryOp, kFunctionCall, kComparison, kConditional）
    std::vector<std::unique_ptr<ExpressionAST>> children;

    // 运算符字符（kBinaryOp, kUnaryOp）
    char op_char = '\0';

    // 比较运算符字符串（kComparison）
    std::string comparison_op;

    // 函数调用标志
    bool is_script_function = false;

    // 变量绑定优化（编译期绑定到槽位索引）
    int variable_slot_index = -1;  ///< -1 = 未绑定, >= 0 = 槽位索引
    bool is_slot_bound = false;    ///< 是否已绑定到槽位

    // 传统变量绑定
    int variable_scope_level = -1;  ///< -1 = 未解析, 0 = 全局, 1+ = 局部作用域
    bool is_builtin_constant = false;

    ExpressionAST() = default;
    explicit ExpressionAST(ExprKind k) : kind(k) {}
};

/**
 * @brief 判断表达式是否可以被编译为 AST
 * @param expression 表达式字符串
 * @return 如果可以编译返回 true
 */
bool can_compile_to_ast(const std::string& expression);

/**
 * @brief 编译表达式为 AST
 * @param expression 表达式字符串
 * @return 编译后的 AST，失败返回 nullptr
 */
std::unique_ptr<ExpressionAST> compile_expression_ast(const std::string& expression);

/**
 * @brief 绑定 AST 中的变量到槽位索引（编译期优化）
 * @param ast AST 节点
 * @param variables 变量解析器
 * @return 是否成功绑定
 */
bool bind_variable_slots(ExpressionAST* ast, const VariableResolver& variables);

/**
 * @brief 求值编译后的 AST
 */
double evaluate_compiled_ast(
    const ExpressionAST* ast,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, std::function<double(const std::vector<double>&)>>* scalar_functions,
    const HasScriptFunctionCallback& has_script_function,
    const InvokeScriptFunctionDecimalCallback& invoke_script_function);

#endif // EXPRESSION_AST_H

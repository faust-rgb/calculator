// ============================================================================
// 编译表达式 AST
// ============================================================================
//
// 提供表达式的 AST 表示，用于：
// 1. 一次编译，多次求值（解决循环中重复解析问题）
// 2. 编译期变量绑定优化
// 3. 直接执行，跳过字符串解析开销
// 4. 表达式特征分析和缓存
//
// 本文件整合了表达式编译、缓存和分析功能，避免功能分散。
// ============================================================================

#ifndef EXPRESSION_AST_H
#define EXPRESSION_AST_H

#include "parser/token_types.h"

#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <functional>
#include <map>
#include <cstdint>

// 前向声明
class VariableResolver;
struct CustomFunction;
struct StoredValue;

// 回调类型定义（来自 types/function.h）
using HasScriptFunctionCallback = std::function<bool(const std::string&)>;
using InvokeScriptFunctionDecimalCallback = std::function<double(const std::string&, const std::vector<double>&)>;

/**
 * @enum ExprKind
 * @brief 表达式 AST 节点类型
 */
enum class ExprKind {
    kNumber,       ///< 数字字面量
    kString,       ///< 字符串字面量
    kVariable,     ///< 变量引用
    kBinaryOp,     ///< 二元运算 (+, -, *, /, ^)
    kUnaryOp,      ///< 一元运算 (-, +)
    kFunctionCall, ///< 函数调用
    kComparison,   ///< 比较运算 (==, !=, <, >, <=, >=)
    kLogicalOp,    ///< 逻辑运算 (&&, ||)
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

// ============================================================================
// 表达式缓存
// ============================================================================

/**
 * @struct ExpressionCache
 * @brief 表达式缓存结构，用于存储预编译结果
 *
 * 优化设计：
 * - 原始文本通过 string_view 引用，不复制
 * - 只有在需要展开时才分配 expanded 字符串
 * - hint 和 features 合并到编译阶段一起计算
 */
struct ExpressionCache {
    // 原始文本存储。缓存可能从临时 std::string 创建，必须拥有文本。
    std::string owned_original;

    // 原始文本视图
    std::string_view original_text;

    // 展开后的表达式（仅当需要展开时才分配，如 inline 命令）
    std::string expanded;

    // 是否使用了展开后的文本
    bool has_expanded = false;

    // 类型提示和特征（延迟计算）
    ExpressionHint hint = ExpressionHint::kUnknown;
    ExpressionFeature features = ExpressionFeature::kNone;
    bool analyzed = false;

    // 编译后的 AST（用于快速求值）
    std::unique_ptr<ExpressionAST> compiled_ast;
    bool is_compiled = false;

    ExpressionCache() = default;

    explicit ExpressionCache(std::string_view text)
        : owned_original(text), original_text(owned_original), has_expanded(false),
          hint(ExpressionHint::kUnknown),
          features(ExpressionFeature::kNone),
          analyzed(false), is_compiled(false) {}

    // 获取当前有效的文本（展开后的或原始的）
    std::string_view effective_text() const {
        return has_expanded ? std::string_view(expanded) : original_text;
    }

    // 设置展开后的文本
    void set_expanded(std::string text) {
        expanded = std::move(text);
        has_expanded = true;
        // 展开后需要重新分析
        analyzed = false;
        is_compiled = false;
    }
};

// ============================================================================
// 表达式分析
// ============================================================================

/**
 * @brief 分析表达式类型（不完整编译）
 * @param expression 表达式字符串
 * @return 类型提示
 */
ExpressionHint analyze_expression_hint(const std::string& expression);

/**
 * @brief 分析表达式特征
 * @param expression 表达式字符串
 * @return 特征位掩码
 */
ExpressionFeature analyze_expression_features(const std::string& expression);

#endif // EXPRESSION_AST_H

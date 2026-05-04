#ifndef SYMBOLIC_EXPRESSION_INTERNAL_H
#define SYMBOLIC_EXPRESSION_INTERNAL_H

#include "symbolic/symbolic_expression.h"
#include "math/mymath_dual.h"

#include <initializer_list>
#include <memory>
#include <string>
#include <unordered_map>
#include <map>
#include <set>
#include <vector>

// 前向声明
class SymbolicPolynomial;

// ============================================================================
// 节点类型枚举
// ============================================================================

/**
 * @enum NodeType
 * @brief 表达式树节点类型
 *
 * 每种类型对应一种表达式构造：
 * - kNumber: 数值常量，如 3.14, 2.0
 * - kVariable: 变量，如 x, y, pi, e
 * - kAdd/kSubtract/kMultiply/kDivide/kPower: 二元运算
 * - kNegate: 一元取负
 * - kFunction: 函数调用，如 sin(x), exp(x)
 * - kVector: 向量表达式，children 存储各分量
 * - kTensor: 张量（矩阵）表达式，children 存储各行
 * - kDifferentialOp: 微分算子，text 存储算子名（grad, div, curl, laplacian）
 */
enum class NodeType {
    kNumber,         ///< 数值常量节点
    kVariable,       ///< 变量节点
    kPi,             ///< 精确常数 pi
    kE,              ///< 精确常数 e
    kInfinity,       ///< 无穷大节点 (+inf 或 -inf)
    kAdd,            ///< 加法节点: left + right
    kSubtract,       ///< 减法节点: left - right
    kMultiply,       ///< 乘法节点: left * right
    kDivide,         ///< 除法节点: left / right
    kPower,          ///< 幂运算节点: left ^ right
    kNegate,         ///< 取负节点: -left
    kFunction,       ///< 函数调用节点: text(left)，如 sin(x)
    kVector,         ///< 向量节点: children 存储各分量
    kTensor,         ///< 张量节点: children 存储各行（每行为 kVector）
    kDifferentialOp, ///< 微分算子节点: text 为算子名，left 为操作数
    kRootOf,         ///< 代数数节点: RootOf(poly, var, root_index)
                     ///< - children[0]: 最小多项式（多项式表达式）
                     ///< - text: 变量名
                     ///< - number_value: 根编号（第几个实根）
};

// ============================================================================
// 表达式节点结构
// ============================================================================

/**
 * @struct SymbolicExpression::Node
 * @brief 表达式树的内部节点表示
 *
 * 采用多态联合体设计，通过 type 字段区分节点类型：
 * - kNumber: 使用 number_value 存储数值
 * - kVariable: 使用 text 存储变量名
 * - kFunction: 使用 text 存储函数名，left 存储参数
 * - 二元运算: 使用 left/right 存储操作数
 * - kNegate: 使用 left 存储被取负的表达式
 * - kVector: 使用 children 存储各分量
 * - kTensor: 使用 children 存储各行，shape 存储 [rows, cols]
 * - kDifferentialOp: 使用 text 存储算子名，left 存储操作数
 *
 * 节点通过 intern_node() 实现驻留（interning），相同结构的节点共享，
 * 以减少内存分配并加速结构比较。
 */
struct SymbolicExpression::Node {
    NodeType type = NodeType::kNumber;  ///< 节点类型
    double number_value = 0.0;          ///< 数值（仅 kNumber 类型使用）
    std::string text;                   ///< 变量名、函数名或微分算子名
    std::shared_ptr<Node> left;         ///< 左子节点/唯一子节点
    std::shared_ptr<Node> right;        ///< 右子节点（二元运算使用）
    std::vector<std::shared_ptr<Node>> children;  ///< 子节点列表（kVector, kTensor）
    std::vector<std::size_t> shape;     ///< 形状信息（kTensor: [rows, cols]）
    mutable std::string structural_key_cache;  ///< 结构键缓存，用于快速比较

    /** @brief 默认构造函数，创建零节点 */
    Node() = default;

    /**
     * @brief 从数值构造节点
     * @param value 数值
     */
    explicit Node(double value) : type(NodeType::kNumber), number_value(value) {}
};

// ============================================================================
// 边界/端点解析协议
// ============================================================================

/**
 * @enum BoundKind
 * @brief 边界值类型
 */
enum class BoundKind {
    kFinite,   ///< 有限数值
    kPosInf,   ///< 正无穷大 (+inf)
    kNegInf,   ///< 负无穷大 (-inf)
};

/**
 * @struct BoundArgument
 * @brief 边界参数解析结果
 *
 * 统一表示有限值和无穷大，用于 integral、limit 等命令的边界解析。
 */
struct BoundArgument {
    BoundKind kind = BoundKind::kFinite;
    double value = 0.0;  ///< 仅当 kind == kFinite 时有效

    /** @brief 是否为有限值 */
    bool is_finite() const { return kind == BoundKind::kFinite; }

    /** @brief 是否为无穷大（正或负） */
    bool is_infinite() const { return kind != BoundKind::kFinite; }

    /** @brief 是否为正无穷 */
    bool is_pos_inf() const { return kind == BoundKind::kPosInf; }

    /** @brief 是否为负无穷 */
    bool is_neg_inf() const { return kind == BoundKind::kNegInf; }

    /** @brief 获取数值（无穷大返回 ±inf） */
    double to_double() const;

    /** @brief 创建有限边界 */
    static BoundArgument finite(double v);

    /** @brief 创建正无穷边界 */
    static BoundArgument pos_inf();

    /** @brief 创建负无穷边界 */
    static BoundArgument neg_inf();
};

/**
 * @brief 解析边界参数字符串
 * @param text 边界参数文本（如 "0", "inf", "+inf", "-inf", "infinity", "oo"）
 * @return 解析后的 BoundArgument
 *
 * 支持的格式：
 * - 有限数值：直接解析为 double
 * - 正无穷：inf, +inf, infinity, +infinity, oo, +oo
 * - 负无穷：-inf, -infinity, -oo
 */
BoundArgument parse_bound_argument(const std::string& text);

/**
 * @brief 检查字符串是否为无穷大字面量
 * @param text 待检查文本
 * @return true 如果是 inf/infinity/oo（不含符号）
 */
bool is_infinity_literal(const std::string& text);

// ============================================================================
// 内部实现命名空间
// ============================================================================

namespace symbolic_expression_internal {

// ============================================================================
// 常量定义
// ============================================================================

/**
 * @brief 格式化精度阈值
 *
 * 用于判断两个浮点数是否"足够接近"，应用于：
 * - 数值比较（如判断是否为零、是否为整数）
 * - 简化规则中的等式检验
 * - 多项式系数的零值裁剪
 */
constexpr double kFormatEps = 1e-12;

// ============================================================================
// 节点构造函数
// ============================================================================

/**
 * @brief 创建数值节点（带驻留）
 * @param value 数值
 * @return 驻留后的节点指针
 */
std::shared_ptr<SymbolicExpression::Node> make_number(double value);

/**
 * @brief 创建变量节点（带驻留）
 * @param name 变量名
 * @return 驻留后的节点指针
 */
std::shared_ptr<SymbolicExpression::Node> make_variable(const std::string& name);

/**
 * @brief 创建无穷大节点（带驻留）
 * @param positive true 表示 +inf，false 表示 -inf
 * @return 驻留后的节点指针
 */
std::shared_ptr<SymbolicExpression::Node> make_infinity(bool positive = true);

/**
 * @brief 创建一元节点（带驻留）
 * @param type 节点类型（kNegate 或 kFunction）
 * @param operand 操作数节点
 * @param text 函数名（仅 kFunction 使用）
 * @return 驻留后的节点指针
 */
std::shared_ptr<SymbolicExpression::Node> make_unary(
    NodeType type,
    std::shared_ptr<SymbolicExpression::Node> operand,
    const std::string& text = "");

/**
 * @brief 创建二元节点（带驻留）
 * @param type 节点类型（kAdd/kSubtract/kMultiply/kDivide/kPower）
 * @param left 左操作数
 * @param right 右操作数
 * @return 驻留后的节点指针
 */
std::shared_ptr<SymbolicExpression::Node> make_binary(
    NodeType type,
    std::shared_ptr<SymbolicExpression::Node> left,
    std::shared_ptr<SymbolicExpression::Node> right);

// ============================================================================
// 表达式构造辅助函数
// ============================================================================

/** @brief 创建取负表达式 -expression */
SymbolicExpression make_negate(SymbolicExpression expression);

/** @brief 创建加法表达式 lhs + rhs */
SymbolicExpression make_add(SymbolicExpression lhs, SymbolicExpression rhs);

/** @brief 创建减法表达式 lhs - rhs */
SymbolicExpression make_subtract(SymbolicExpression lhs, SymbolicExpression rhs);

/** @brief 创建乘法表达式 lhs * rhs */
SymbolicExpression make_multiply(SymbolicExpression lhs, SymbolicExpression rhs);

/** @brief 创建除法表达式 lhs / rhs */
SymbolicExpression make_divide(SymbolicExpression lhs, SymbolicExpression rhs);

/** @brief 创建幂运算表达式 lhs ^ rhs */
SymbolicExpression make_power(SymbolicExpression lhs, SymbolicExpression rhs);

/**
 * @brief 创建函数调用表达式
 * @param name 函数名（如 sin, cos, exp, ln）
 * @param argument 函数参数
 * @return 函数调用表达式
 */
SymbolicExpression make_function(const std::string& name, SymbolicExpression argument);

/**
 * @brief 创建 RootOf 表达式（代数数）
 * @param polynomial 最小多项式（作为表达式）
 * @param var_name 变量名
 * @param root_index 根编号（第几个实根，从 0 开始）
 * @return RootOf 表达式，表示多项式的指定根
 *
 * 例如: make_rootof(x^3 - 2, "x", 0) 表示 x^3 - 2 = 0 的第一个实根
 */
SymbolicExpression make_rootof(
    const SymbolicExpression& polynomial,
    const std::string& var_name,
    int root_index);

/**
 * @brief 从多项式创建 RootOf 表达式
 * @param polynomial 符号多项式
 * @param root_index 根编号
 * @return RootOf 表达式
 */
SymbolicExpression make_rootof_from_polynomial(
    const SymbolicPolynomial& polynomial,
    int root_index);

// ============================================================================
// 字符串转换与结构键
// ============================================================================

/**
 * @brief 将表达式节点转换为字符串
 * @param node 表达式节点
 * @param parent_precedence 父节点优先级，用于决定是否添加括号
 * @return 表达式字符串，如 "x ^ 2 + 2 * x + 1"
 */
std::string to_string_impl(const std::shared_ptr<SymbolicExpression::Node>& node,
                           int parent_precedence);

/**
 * @brief 计算节点的结构键
 * @param node 表达式节点
 * @return 唯一标识表达式结构的字符串
 *
 * 结构键用于：
 * - 节点驻留时的查重
 * - 简化缓存和导数缓存的键
 * - 表达式结构相等性比较
 *
 * 例如 x + 1 和 1 + x 有不同的结构键，但数学上等价。
 */
std::string node_structural_key(const std::shared_ptr<SymbolicExpression::Node>& node);

// ============================================================================
// 表达式谓词函数
// ============================================================================

/** @brief 检查表达式是否为指定变量 */
bool expr_is_variable(const SymbolicExpression& expression, const std::string& name);

/** @brief 检查表达式是否为零 */
bool expr_is_zero(const SymbolicExpression& expression);

/** @brief 检查表达式是否为一 */
bool expr_is_one(const SymbolicExpression& expression);

/** @brief 检查表达式是否为负一 */
bool expr_is_minus_one(const SymbolicExpression& expression);

/** @brief 检查表达式是否为数值，可选输出该数值 */
bool expr_is_number(const SymbolicExpression& expression, double* value = nullptr);

/**
 * @brief 检查表达式是否为无穷大
 * @param expression 表达式
 * @param positive 输出参数，true 表示 +inf，false 表示 -inf（可选）
 * @return true 如果表达式是无穷大节点或数值无穷大
 */
bool expr_is_infinity(const SymbolicExpression& expression, bool* positive = nullptr);

/**
 * @brief 尝试数值求值
 * @param node 表达式节点
 * @param value 输出参数，求值结果
 * @return true 如果表达式可完全求值为数值
 *
 * 递归求值表达式树，遇到变量则失败。
 */
bool try_evaluate_numeric_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                               double* value);

/**
 * @brief Context for dual number evaluation (forward-mode AutoDiff)
 */
struct DualEvaluationContext {
    std::string differentiation_variable;  ///< Variable to differentiate w.r.t.
    double point_value;                     ///< Value of the variable at the point
    std::unordered_map<std::string, double> other_values;  ///< Values of other variables
};

/**
 * @brief Evaluate expression node with dual numbers for automatic differentiation.
 *
 * For the differentiation variable, returns dual(point_value, 1.0).
 * For other variables, returns dual(their_value, 0.0).
 * All operations propagate derivatives using dual arithmetic.
 */
bool try_evaluate_dual_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                            const DualEvaluationContext& ctx,
                            mymath::dual<double>* result);

// ============================================================================
// Expression size monitoring
// ============================================================================

/**
 * @brief Count nodes in expression tree
 * @param node Expression root
 * @return Total number of nodes (including leaves)
 */
std::size_t count_nodes(const std::shared_ptr<SymbolicExpression::Node>& node);

// ============================================================================
// 多项式操作
// ============================================================================

/** @brief 裁剪多项式系数向量末尾的零 */
void trim_polynomial_coefficients(std::vector<double>* coefficients);

/** @brief 多项式加法（系数向量形式） */
std::vector<double> polynomial_add_impl(const std::vector<double>& lhs,
                                        const std::vector<double>& rhs);

/** @brief 多项式减法（系数向量形式） */
std::vector<double> polynomial_subtract_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs);

/** @brief 多项式乘法（系数向量形式），朴素 O(n²) 算法 */
std::vector<double> polynomial_multiply_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs);

/** @brief 检查表达式是否为关于指定变量的多项式（支持符号系数） */
bool is_symbolic_polynomial(const SymbolicExpression& expression,
                            const std::string& variable_name);

/**
 * @brief 从简化后的表达式提取符号多项式系数
 * @param expression 已简化的表达式
 * @param variable_name 多项式变量
 * @param coefficients 输出系数向量（低次到高次，系数本身为表达式）
 * @return true 如果表达式是关于该变量的多项式
 */
bool symbolic_polynomial_coefficients_from_simplified(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    std::vector<SymbolicExpression>* coefficients);

bool symbolic_laurent_coefficients(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    std::map<int, SymbolicExpression>* coefficients);

/** @brief 裁剪符号多项式系数向量末尾的零 */
void trim_symbolic_polynomial_coefficients(std::vector<SymbolicExpression>* coefficients);

/**
 * @brief 从简化后的表达式提取多项式系数
 * @param expression 已简化的表达式
 * @param variable_name 多项式变量
 * @param coefficients 输出系数向量（低次到高次）
 * @return true 如果表达式是关于该变量的多项式
 */
bool polynomial_coefficients_from_simplified(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    std::vector<double>* coefficients);

/**
 * @brief 从系数向量构建多项式表达式
 * @param coefficients 系数向量（低次到高次）
 * @param variable_name 变量名
 * @return 多项式表达式
 */
SymbolicExpression build_polynomial_expression_from_coefficients(
    const std::vector<double>& coefficients,
    const std::string& variable_name);

// ============================================================================
// 变量分析
// ============================================================================

/**
 * @brief 收集表达式中的所有标识符变量
 * @param expression 表达式
 * @param names 输出变量名列表（可能重复）
 *
 * 排除保留符号 pi, e, i。
 */
void collect_identifier_variables(const SymbolicExpression& expression,
                                  std::vector<std::string>* names);

/**
 * @brief 获取表达式中唯一的标识符变量
 * @return 变量名（如果只有一个），否则返回空串
 */
std::string unique_identifier_variable(const SymbolicExpression& expression);

/**
 * @brief 检查字符串是否为合法的标识符变量名
 * @param name 待检查的字符串
 * @return true 如果是合法标识符（字母开头，仅含字母数字下划线）
 */
bool is_identifier_variable_name(const std::string& name);

// ============================================================================
// 简化与替换
// ============================================================================

/**
 * @brief 简化的内部实现（带多轮迭代）
 * @param expression 待简化表达式
 * @return 简化后的表达式
 *
 * 最多进行 16 轮简化迭代，直到表达式结构不再变化。
 */
SymbolicExpression simplify_impl(const SymbolicExpression& expression);

/**
 * @brief Simplify with a node count budget
 * @param expression Expression to simplify
 * @param max_nodes Maximum allowed node count
 * @return Simplified expression, or partial result if budget exceeded
 */
SymbolicExpression simplify_with_budget_impl(const SymbolicExpression& expression,
                                             std::size_t max_nodes);

/** @brief 强制完全展开的内部实现 */
SymbolicExpression expand_impl(const SymbolicExpression& expression);

/**
 * @brief 单轮简化
 * @param expression 待简化表达式
 * @return 单轮简化后的表达式
 *
 * 应用所有简化规则一次，不保证达到最优形式。
 */
SymbolicExpression simplify_once(const SymbolicExpression& expression);

/**
 * @brief 变量替换的内部实现
 * @param expression 表达式
 * @param variable_name 被替换的变量
 * @param replacement 替换表达式
 * @return 替换后的表达式
 */
SymbolicExpression substitute_impl(const SymbolicExpression& expression,
                                   const std::string& variable_name,
                                   const SymbolicExpression& replacement);

// ============================================================================
// 表达式匹配与分解
// ============================================================================

/** @brief 检查两个表达式是否结构相同（基于结构键） */
bool expressions_match(const SymbolicExpression& lhs, const SymbolicExpression& rhs);

/**
 * @brief 检查表达式是否已知为正
 * @return true 对于正数、exp、sqrt、abs、pi、e 等
 */
bool is_known_positive_expression(const SymbolicExpression& expression);

/**
 * @brief 分解常数乘积表达式
 * @param expression 表达式
 * @param variable_name 变量名
 * @param constant 输出常数因子
 * @param remainder 输出剩余表达式
 * @return true 如果 expression = constant * remainder
 */
bool decompose_constant_times_expression(const SymbolicExpression& expression,
                                         const std::string& variable_name,
                                         double* constant,
                                         SymbolicExpression* remainder);

bool decompose_linear(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      double* a,
                      double* b);

/** @brief 分解符号线性表达式 a*x + b */
bool symbolic_decompose_linear(const SymbolicExpression& expression,
                               const std::string& variable_name,
                               SymbolicExpression* a,
                               SymbolicExpression* b);

/**
 * @brief 分解幂因子
 * @param expression 表达式
 * @param base 输出底数
 * @param exponent 输出指数
 * @return true（总是成功，非幂表达式返回 exponent=1）
 */
bool decompose_power_factor(const SymbolicExpression& expression,
                            SymbolicExpression* base,
                            double* exponent);

// ============================================================================
// 乘积项收集与重建
// ============================================================================

/**
 * @brief 收集乘积项
 * @param expression 表达式
 * @param numeric_factor 输出数值因子（累乘）
 * @param symbolic_factors 输出符号因子列表（追加）
 *
 * 例如 2 * x * 3 * y 收集为 numeric_factor=6, symbolic_factors=[x, y]
 */
void collect_multiplicative_terms(const SymbolicExpression& expression,
                                  double* numeric_factor,
                                  std::vector<SymbolicExpression>* symbolic_factors);

/**
 * @brief 收集除法因子（展开幂为重复因子）
 * @param expression 表达式
 * @param numeric_factor 输出数值因子
 * @param symbolic_factors 输出符号因子列表
 *
 * 例如 x^3 收集为 [x, x, x]，便于后续约分。
 */
void collect_division_factors(const SymbolicExpression& expression,
                              double* numeric_factor,
                              std::vector<SymbolicExpression>* symbolic_factors);

/**
 * @brief 收集加法项
 * @param expression 表达式
 * @param terms 输出项列表
 */
void collect_additive_expressions(const SymbolicExpression& expression,
                                  std::vector<SymbolicExpression>* terms);

/**
 * @brief 从因子列表重建乘积表达式
 * @param numeric_factor 数值因子
 * @param factors 符号因子列表
 * @return 乘积表达式
 */
SymbolicExpression rebuild_product_expression(
    double numeric_factor,
    const std::vector<SymbolicExpression>& factors);

/**
 * @brief 构建排序后的和式
 * @param terms 项列表
 * @return 按规范顺序排列的和式
 */
SymbolicExpression make_sorted_sum(std::vector<SymbolicExpression> terms);

/**
 * @brief 构建排序后的积式（合并相同底数的幂）
 * @param numeric_factor 数值因子
 * @param factors 因子列表
 * @return 规范化的乘积表达式
 */
SymbolicExpression make_sorted_product(double numeric_factor,
                                       std::vector<SymbolicExpression> factors);

/**
 * @brief 重建幂表达式（处理负指数）
 * @param base 底数
 * @param exponent 指数
 * @return base^exponent 或 1/base^(-exponent)
 */
SymbolicExpression rebuild_power_difference(const SymbolicExpression& base,
                                            double exponent);

// ============================================================================
// 多项式简化规则
// ============================================================================

/** @brief 检查表达式是否为单变量多项式 */
bool is_single_variable_polynomial(const SymbolicExpression& expression);

/**
 * @brief 尝试将多项式规范化
 * @param expression 表达式
 * @return 规范化后的多项式（如果成功），否则原样返回
 */
SymbolicExpression maybe_canonicalize_polynomial(const SymbolicExpression& expression);

/** @brief 检查表达式是否为多项式，输出多项式形式 */
bool polynomial_expression(const SymbolicExpression& expression,
                           const std::string& variable_name,
                           SymbolicExpression* polynomial);

// ============================================================================
// 简化规则辅助函数
// ============================================================================

/**
 * @brief 尝试合并相似项
 * @param left 左项
 * @param right 右项
 * @param right_sign 右项符号（+1 或 -1）
 * @param combined 输出合并结果
 * @return true 如果成功合并
 *
 * 例如 2*x + 3*x → 5*x
 */
bool try_combine_like_terms(const SymbolicExpression& left,
                            const SymbolicExpression& right,
                            double right_sign,
                            SymbolicExpression* combined);

/**
 * @brief 尝试提取公因子
 * @param left 左项
 * @param right 右项
 * @param right_sign 右项符号
 * @param combined 输出提取结果
 * @return true 如果成功提取
 *
 * 例如 2*x + 2*y → 2*(x + y)
 */
bool try_factor_common_terms(const SymbolicExpression& left,
                             const SymbolicExpression& right,
                             double right_sign,
                             SymbolicExpression* combined);

/**
 * @brief 尝试多项式商约分（整除情况）
 * @param left 分子
 * @param right 分母
 * @param reduced 输出约分结果
 * @return true 如果可整除
 *
 * 例如 (x^2 - 1) / (x - 1) → x + 1
 */
bool try_reduce_polynomial_quotient(const SymbolicExpression& left,
                                    const SymbolicExpression& right,
                                    SymbolicExpression* reduced);

/**
 * @brief 尝试多项式 GCD 约分
 * @param left 分子
 * @param right 分母
 * @param reduced 输出约分结果
 * @return true 如果存在非常数 GCD
 *
 * 例如 (x^3 - x) / (x^2 - 2x + 1) → (x^2 + x) / (x - 1)
 */
bool try_reduce_polynomial_gcd_quotient(const SymbolicExpression& left,
                                        const SymbolicExpression& right,
                                        SymbolicExpression* reduced);

/**
 * @brief 尝试规范因子商约分
 * @param numerator 分子
 * @param denominator 分母
 * @param quotient 输出结果
 * @return true 如果进行了约分
 */
bool try_canonical_factor_quotient(const SymbolicExpression& numerator,
                                   const SymbolicExpression& denominator,
                                   SymbolicExpression* quotient);

// ============================================================================
// 三角函数辅助
// ============================================================================

/**
 * @brief 检查是否为函数的平方
 * @param expression 表达式
 * @param function_name 函数名（如 sin, cos）
 * @param argument_key 输出参数的结构键
 * @return true 如果 expression = function(argument)^2
 */
bool is_squared_function(const SymbolicExpression& expression,
                         const std::string& function_name,
                         std::string* argument_key);

/**
 * @brief 分解符号的数值倍数
 * @param expression 表达式
 * @param symbol_name 符号名（如 pi）
 * @param coefficient 输出系数
 * @return true 如果 expression = coefficient * symbol
 */
bool decompose_numeric_multiple_of_symbol(const SymbolicExpression& expression,
                                          const std::string& symbol_name,
                                          double* coefficient);

/**
 * @brief 检查数值是否匹配候选列表之一
 * @param value 待检查值
 * @param candidates 候选值列表
 * @return true 如果 value 近似等于某个候选值
 */
bool numeric_matches_any(double value,
                         const std::initializer_list<double>& candidates);

/** @brief 创建 sqrt(3) 符号表达式 */
SymbolicExpression sqrt3_symbol();

/** @brief 创建 1/2 符号表达式 */
SymbolicExpression half_symbol();

// ============================================================================
// 积分辅助
// ============================================================================

/**
 * @brief 计算阶乘（浮点结果）
 * @param exponent 非负整数
 * @return exponent! 作为 double
 */
double factorial_double(int exponent);

/**
 * @brief 积分多项式乘函数
 * @param polynomial 多项式部分
 * @param function_name 函数名（sin/cos/exp）
 * @param argument 函数参数
 * @param variable_name 积分变量
 * @param integrated 输出积分结果
 * @return true 如果成功
 *
 * 使用递推公式处理 x^n * sin(ax+b) 等形式。
 */
bool integrate_polynomial_times_function(const SymbolicExpression& polynomial,
                                         const std::string& function_name,
                                         const SymbolicExpression& argument,
                                         const std::string& variable_name,
                                         SymbolicExpression* integrated);

/**
 * @brief 符号逆二次式积分（常数分子）
 *
 * 计算 A * ∫ dx / (ax^2 + bx + c)^n
 * 使用符号公式，返回 atan 或 atanh/log 形式。
 */
SymbolicExpression integrate_symbolic_inverse_quadratic(
    const SymbolicExpression& a,
    const SymbolicExpression& b,
    const SymbolicExpression& c,
    const SymbolicExpression& coeff,
    int power,
    const std::string& variable_name);

/**
 * @brief 符号逆二次式积分（线性分子）
 *
 * 计算 B * ∫ x dx / (ax^2 + bx + c)^n
 */
SymbolicExpression integrate_symbolic_inverse_quadratic_linear(
    const SymbolicExpression& a,
    const SymbolicExpression& b,
    const SymbolicExpression& c,
    const SymbolicExpression& coeff,
    int power,
    const std::string& variable_name);

/**
 * @brief 符号部分分式积分
 *
 * 使用系数恒等式方法求解部分分式系数。
 */
bool integrate_symbolic_partial_fractions(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated);

// ============================================================================
// 特殊函数构造
// ============================================================================

/**
 * @brief 创建单位阶跃函数 step(variable - location)
 * @param variable_name 变量名
 * @param location 跳跃位置
 * @return step 表达式
 */
SymbolicExpression make_step_expression(const std::string& variable_name,
                                        double location);

/**
 * @brief 创建 Dirac delta 函数 delta(variable - shift)
 * @param variable_name 变量名
 * @param shift 平移量
 * @return delta 表达式
 */
SymbolicExpression make_delta_expression(const std::string& variable_name,
                                         double shift);

// ============================================================================
// 积分变换实现
// ============================================================================

/** @brief Fourier 变换内部实现 */
SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable);

/** @brief 逆 Fourier 变换内部实现 */
SymbolicExpression inverse_fourier_transform_impl(
    const SymbolicExpression& expression,
    const std::string& frequency_variable,
    const std::string& time_variable);

/** @brief Laplace 变换内部实现 */
SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable);

/** @brief 逆 Laplace 变换内部实现 */
SymbolicExpression inverse_laplace_transform_impl(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& time_variable);

/** @brief z 变换内部实现 */
SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable);

/** @brief 逆 z 变换内部实现 */
SymbolicExpression inverse_z_transform_impl(
    const SymbolicExpression& expression,
    const std::string& transform_variable,
    const std::string& index_variable);

}  // namespace symbolic_expression_internal

#endif

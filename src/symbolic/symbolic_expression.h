#ifndef SYMBOLIC_EXPRESSION_H
#define SYMBOLIC_EXPRESSION_H

#include <memory>
#include <string>
#include <vector>

/**
 * @file symbolic_expression.h
 * @brief 符号表达式库
 *
 * 提供符号数学表达式的表示、解析和操作功能。
 * 支持符号微分、积分、简化和多项式分析。
 *
 * 使用共享指针实现表达式树的共享和惰性求值。
 */

/**
 * @class SymbolicExpression
 * @brief 符号表达式类
 *
 * 表示一个数学表达式的符号形式，可以包含：
 * - 数字常量
 * - 变量
 * - 二元运算（+、-、*、/、^）
 * - 函数调用（sin、cos、exp、ln 等）
 */
class SymbolicExpression {
public:
    /** @brief 表达式节点的前向声明（Pimpl 模式） */
    struct Node;

    /** @brief 默认构造函数，创建空表达式 */
    SymbolicExpression();

    /**
     * @brief 从字符串解析表达式
     * @param text 表达式字符串，如 "x^2 + 2*x + 1"
     * @return 解析后的符号表达式
     * @throw std::runtime_error 当解析失败时抛出
     */
    static SymbolicExpression parse(const std::string& text);

    /**
     * @brief 创建数字常量表达式
     * @param value 数值
     * @return 常量表达式
     */
    static SymbolicExpression number(double value);

    /**
     * @brief 设置符号表达式字符串输出中的十进制显示有效位数
     * @param precision 有效位数，范围 1..17
     */
    static void set_display_precision(int precision);

    /**
     * @brief 创建变量表达式
     * @param name 变量名
     * @return 变量表达式
     */
    static SymbolicExpression variable(const std::string& name);

    /**
     * @brief 转换为字符串表示
     * @return 表达式字符串，如 "x ^ 2 + 2 * x + 1"
     */
    std::string to_string() const;

    /**
     * @brief 计算符号导数
     * @param variable_name 求导变量名
     * @return 导数表达式
     *
     * 使用符号微分规则（链式法则、乘积法则等）。
     */
    SymbolicExpression derivative(const std::string& variable_name) const;

    /**
     * @brief 计算多个变量上的梯度
     * @param variable_names 求导变量名列表
     * @return 按变量顺序排列的一阶偏导表达式
     */
    std::vector<SymbolicExpression> gradient(
        const std::vector<std::string>& variable_names) const;

    /**
     * @brief 计算多个变量上的 Hessian 矩阵
     * @param variable_names 求导变量名列表
     * @return Hessian 矩阵表达式
     */
    std::vector<std::vector<SymbolicExpression>> hessian(
        const std::vector<std::string>& variable_names) const;

    /**
     * @brief 计算表达式列表对变量列表的 Jacobian 矩阵
     * @param expressions 表达式列表
     * @param variable_names 求导变量名列表
     * @return Jacobian 矩阵表达式
     */
    static std::vector<std::vector<SymbolicExpression>> jacobian(
        const std::vector<SymbolicExpression>& expressions,
        const std::vector<std::string>& variable_names);

    /**
     * @brief 计算符号积分
     * @param variable_name 积分变量名
     * @return 积分表达式（不含常数 C）
     *
     * 支持基本积分规则和简单替换。
     */
    SymbolicExpression integral(const std::string& variable_name) const;

    /**
     * @brief 计算符号 Fourier 变换
     * @param time_variable 时域变量
     * @param frequency_variable 频域变量
     * @return 变换后的表达式
     */
    SymbolicExpression fourier_transform(const std::string& time_variable,
                                         const std::string& frequency_variable) const;

    /**
     * @brief 计算符号逆 Fourier 变换
     * @param frequency_variable 频域变量
     * @param time_variable 时域变量
     * @return 逆变换后的表达式
     */
    SymbolicExpression inverse_fourier_transform(
        const std::string& frequency_variable,
        const std::string& time_variable) const;

    /**
     * @brief 计算符号 Laplace 变换
     * @param time_variable 时域变量
     * @param transform_variable 复频域变量
     * @return 变换后的表达式
     */
    SymbolicExpression laplace_transform(const std::string& time_variable,
                                         const std::string& transform_variable) const;

    /**
     * @brief 计算符号逆 Laplace 变换
     * @param transform_variable 复频域变量
     * @param time_variable 时域变量
     * @return 逆变换后的表达式
     */
    SymbolicExpression inverse_laplace_transform(
        const std::string& transform_variable,
        const std::string& time_variable) const;

    /**
     * @brief 计算符号 z 变换
     * @param index_variable 序列索引变量
     * @param transform_variable z 域变量
     * @return 变换后的表达式
     */
    SymbolicExpression z_transform(const std::string& index_variable,
                                   const std::string& transform_variable) const;

    /**
     * @brief 计算符号逆 z 变换
     * @param transform_variable z 域变量
     * @param index_variable 序列索引变量
     * @return 逆变换后的表达式
     */
    SymbolicExpression inverse_z_transform(const std::string& transform_variable,
                                           const std::string& index_variable) const;

    /**
     * @brief 用另一个表达式替换指定变量
     * @param variable_name 被替换的变量名
     * @param replacement 替换表达式
     * @return 代换后的表达式
     */
    SymbolicExpression substitute(const std::string& variable_name,
                                  const SymbolicExpression& replacement) const;

    /**
     * @brief 简化表达式
     * @return 简化后的表达式
     *
     * 执行代数简化：
     * - 常数折叠（如 2 + 3 → 5）
     * - 恒等式消去（如 x + 0 → x, x * 1 → x）
     * - 幂运算简化
     */
    SymbolicExpression simplify() const;

    /**
     * @brief 强制完全展开表达式（乘法分配律与多项式幂次展开）
     * @return 展开后的表达式
     */
    SymbolicExpression expand() const;

    /**
     * @brief 检查表达式是否不依赖于指定变量
     * @param variable_name 变量名
     * @return true 如果表达式中不包含该变量
     */
    bool is_constant(const std::string& variable_name) const;

    /**
     * @brief 检查表达式是否为纯数字
     * @param value 可选的输出参数，用于获取数值
     * @return true 如果表达式是数字常量
     */
    bool is_number(double* value = nullptr) const;

    /**
     * @brief 检查表达式是否为指定变量
     * @param variable_name 变量名
     * @return true 如果表达式就是该变量
     */
    bool is_variable_named(const std::string& variable_name) const;

    /**
     * @brief 提取多项式系数
     * @param variable_name 多项式变量名
     * @param coefficients 输出参数，系数向量（低次到高次）
     * @return true 如果表达式是关于该变量的多项式
     *
     * 例如：x^2 + 2x + 1 的系数为 [1, 2, 1]
     */
    bool polynomial_coefficients(const std::string& variable_name,
                                 std::vector<double>* coefficients) const;

    /**
     * @brief 获取表达式中出现的普通标识符变量名（不含 pi/e）
     * @return 去重后的变量名列表
     */
    std::vector<std::string> identifier_variables() const;

    /**
     * @brief 提取公共子表达式 (CSE)
     * @return 包含重复出现的非平凡子表达式及其出现次数的列表
     */
    std::vector<std::pair<SymbolicExpression, int>> common_subexpressions() const;

    /**
     * @brief 从现有节点构造表达式（内部使用）
     * @param node 表达式节点
     */
    explicit SymbolicExpression(std::shared_ptr<Node> node);

    std::shared_ptr<Node> node_;  ///< 表达式树根节点
};

// ============================================================================
// 运算符重载
// ============================================================================

SymbolicExpression operator+(const SymbolicExpression& lhs, const SymbolicExpression& rhs);
SymbolicExpression operator-(const SymbolicExpression& lhs, const SymbolicExpression& rhs);
SymbolicExpression operator*(const SymbolicExpression& lhs, const SymbolicExpression& rhs);
SymbolicExpression operator/(const SymbolicExpression& lhs, const SymbolicExpression& rhs);
SymbolicExpression operator^(const SymbolicExpression& lhs, const SymbolicExpression& rhs);
SymbolicExpression operator-(const SymbolicExpression& expr);

#endif

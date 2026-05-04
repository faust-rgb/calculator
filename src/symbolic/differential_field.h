#ifndef DIFFERENTIAL_FIELD_H
#define DIFFERENTIAL_FIELD_H

#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_polynomial.h"
#include "symbolic/risch_types.h"
#include <string>
#include <vector>
#include <set>
#include <functional>

/**
 * @file differential_field.h
 * @brief 微分域的形式化抽象
 *
 * 提供微分域 K、常数域 C = {k ∈ K : D(k) = 0}、微分塔 K_0 ⊂ K_1 ⊂ ... ⊂ K_n
 * 的形式化表示和运算。
 *
 * 在 Risch 算法中，我们需要：
 * 1. 判断元素是否在常数域 C 中 (D(k) = 0)
 * 2. 判断元素在哪一层域 K_i 中
 * 3. 计算导数 D: K → K
 * 4. 判断元素是否在基域 K_0 中
 */

/**
 * @brief 精确有理数类 - 用于代数数的隔离区间边界
 *
 * 避免浮点精度问题，用于 Sturm 序列的实根隔离
 */
class ExactRational {
public:
    int64_t numerator;
    int64_t denominator;

    ExactRational() : numerator(0), denominator(1) {}

    ExactRational(int64_t n) : numerator(n), denominator(1) {}

    ExactRational(int64_t n, int64_t d) : numerator(n), denominator(d) {
        normalize();
    }

    // 从 double 构造 (近似)
    static ExactRational from_double(double value, int64_t max_den = 1000000);

    void normalize() {
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
        // GCD 约简
        int64_t g = gcd(std::abs(numerator), std::abs(denominator));
        if (g > 1) {
            numerator /= g;
            denominator /= g;
        }
    }

    ExactRational add(const ExactRational& other) const {
        return ExactRational(
            numerator * other.denominator + other.numerator * denominator,
            denominator * other.denominator
        );
    }

    ExactRational subtract(const ExactRational& other) const {
        return ExactRational(
            numerator * other.denominator - other.numerator * denominator,
            denominator * other.denominator
        );
    }

    ExactRational multiply(const ExactRational& other) const {
        return ExactRational(numerator * other.numerator, denominator * other.denominator);
    }

    ExactRational divide(const ExactRational& other) const {
        return ExactRational(numerator * other.denominator, denominator * other.numerator);
    }

    ExactRational negate() const {
        return ExactRational(-numerator, denominator);
    }

    int compare(const ExactRational& other) const {
        int64_t left = numerator * other.denominator;
        int64_t right = other.numerator * denominator;
        if (left < right) return -1;
        if (left > right) return 1;
        return 0;
    }

    bool operator<(const ExactRational& other) const { return compare(other) < 0; }
    bool operator>(const ExactRational& other) const { return compare(other) > 0; }
    bool operator<=(const ExactRational& other) const { return compare(other) <= 0; }
    bool operator>=(const ExactRational& other) const { return compare(other) >= 0; }
    bool operator==(const ExactRational& other) const { return compare(other) == 0; }
    bool operator!=(const ExactRational& other) const { return compare(other) != 0; }

    double to_double() const {
        return static_cast<double>(numerator) / static_cast<double>(denominator);
    }

    bool is_zero() const { return numerator == 0; }
    bool is_positive() const { return numerator > 0; }
    bool is_negative() const { return numerator < 0; }

    std::string to_string() const {
        if (denominator == 1) {
            return std::to_string(numerator);
        }
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }

private:
    static int64_t gcd(int64_t a, int64_t b) {
        while (b != 0) {
            int64_t t = b;
            b = a % b;
            a = t;
        }
        return a;
    }
};

// Note: Use ExactRational directly to avoid conflict with precise/rational.h

/**
 * @brief 微分域类
 *
 * 表示微分域 K = K_0(t_1)(t_2)...(t_n) 其中每个 t_i 是对数、指数或代数扩展
 */
class DifferentialField {
public:
    // 基域变量名 (积分变量)
    std::string base_variable;

    // 扩展塔 K_0 ⊂ K_1 ⊂ ... ⊂ K_n
    // K_0 = C(x) 是有理函数域
    // K_i = K_{i-1}(t_i) 是第 i 层扩展
    std::vector<DifferentialExtension> tower;

    // 已知常数表 (这些元素的导数为零)
    std::set<std::string> known_constants;

    // 默认构造函数
    DifferentialField() = default;

    // 从表达式构建微分域
    static DifferentialField from_expression(
        const SymbolicExpression& expr,
        const std::string& x_var);

    // ==================== 核心操作 ====================

    /**
     * @brief 计算导数 D(expr)
     *
     * 对于 f ∈ K，计算 D(f) = df/dx
     * 使用链式法则: D(f) = ∂f/∂x + Σ (∂f/∂t_i) * D(t_i)
     */
    SymbolicExpression derive(const SymbolicExpression& expr) const;

    /**
     * @brief 判断元素是否在常数域 C 中
     *
     * C = {k ∈ K : D(k) = 0}
     *
     * 判断方法:
     * 1. 语法检查: 不含积分变量和塔变量
     * 2. 导数检查: D(expr) = 0
     * 3. 已知常数表: pi, e, sin(1) 等
     */
    bool is_constant(const SymbolicExpression& expr) const;

    /**
     * @brief 判断元素在哪一层域中
     *
     * 返回最小的 i 使得 expr ∈ K_i
     * 如果 expr ∈ K_0，返回 0
     * 如果 expr 不在任何 K_i 中，返回 -1
     */
    int field_level(const SymbolicExpression& expr) const;

    /**
     * @brief 判断元素是否在基域 K_0 中
     *
     * K_0 = C(x) 是有理函数域
     */
    bool is_in_base_field(const SymbolicExpression& expr) const;

    /**
     * @brief 判断元素是否在第 i 层域 K_i 中
     */
    bool is_in_field(const SymbolicExpression& expr, int level) const;

    /**
     * @brief 获取第 i 层扩展的导数 t_i'
     */
    SymbolicExpression get_extension_derivative(int level) const;

    /**
     * @brief 获取第 i 层扩展的变量名 t_i
     */
    std::string get_extension_variable(int level) const;

    /**
     * @brief 获取塔的高度
     */
    int tower_height() const { return static_cast<int>(tower.size()); }

    /**
     * @brief 检查表达式是否包含塔中的任何变量
     */
    bool contains_tower_variable(const SymbolicExpression& expr) const;

    /**
     * @brief 获取表达式依赖的所有塔变量
     */
    std::set<std::string> get_tower_dependencies(const SymbolicExpression& expr) const;

    // ==================== 辅助操作 ====================

    /**
     * @brief 添加已知常数
     */
    void add_known_constant(const std::string& name);

    /**
     * @brief 初始化标准常数表 (pi, e, 等)
     */
    void initialize_standard_constants();

    /**
     * @brief 将表达式中的塔变量替换回原始形式
     */
    SymbolicExpression substitute_back(const SymbolicExpression& expr) const;

private:
    // 内部导数计算
    SymbolicExpression derive_impl(const SymbolicExpression& expr,
                                    const std::set<std::string>& active_vars) const;

    // 检查是否是已知常数
    bool is_known_constant(const SymbolicExpression& expr) const;

    // 检查表达式是否只包含指定层及以下的变量
    bool only_uses_vars_up_to_level(const SymbolicExpression& expr, int level) const;
};

/**
 * @brief 带域关联的多项式
 *
 * 表示多项式环 K_i[t] 中的元素，其中系数在 K_i 中
 */
class PolynomialOverField {
public:
    SymbolicPolynomial poly;
    const DifferentialField* field;
    int field_level;  // 系数在哪层域中

    PolynomialOverField() : field(nullptr), field_level(0) {}

    PolynomialOverField(const SymbolicPolynomial& p,
                        const DifferentialField& f,
                        int level = 0)
        : poly(p), field(&f), field_level(level) {}

    // 算术运算
    PolynomialOverField add(const PolynomialOverField& other) const;
    PolynomialOverField subtract(const PolynomialOverField& other) const;
    PolynomialOverField multiply(const PolynomialOverField& other) const;

    // 多项式除法
    bool divide(const PolynomialOverField& divisor,
                PolynomialOverField* quotient,
                PolynomialOverField* remainder) const;

    /**
     * @brief 计算微分 D(P)
     *
     * D(P) = dP/dx + Σ (dP/dt_i) * D(t_i)
     *
     * 这是微分域上多项式的正确导数计算
     */
    PolynomialOverField differential() const;

    /**
     * @brief 计算形式导数 dP/dt
     *
     * 这只是对多项式变量的偏导数，不考虑微分结构
     */
    PolynomialOverField formal_derivative() const;

    // 属性
    int degree() const { return poly.degree(); }
    bool is_zero() const { return poly.is_zero(); }
    SymbolicExpression leading_coefficient() const { return poly.leading_coefficient(); }
    SymbolicExpression coefficient(int i) const { return poly.coefficient(i); }

    // 转换
    SymbolicExpression to_expression() const;
    static PolynomialOverField from_expression(const SymbolicExpression& expr,
                                                const std::string& var,
                                                const DifferentialField& field,
                                                int level = 0);
};

/**
 * @brief 微分塔构建器
 *
 * 从表达式构建微分塔，处理代数独立性
 */
class DifferentialTowerBuilder {
public:
    /**
     * @brief 构建微分塔
     *
     * @param expr 要积分的表达式
     * @param x_var 积分变量
     * @return 构建好的微分域
     */
    static DifferentialField build(const SymbolicExpression& expr,
                                   const std::string& x_var);

private:
    // 收集表达式中的所有扩展
    static void collect_extensions(
        const SymbolicExpression& expr,
        const std::string& x_var,
        std::vector<std::pair<SymbolicExpression, DifferentialExtension::Kind>>& extensions);

    // 检查代数独立性
    static IndependenceCheck check_independence(
        const SymbolicExpression& arg,
        DifferentialExtension::Kind kind,
        const DifferentialField& current_field);

    // 拓扑排序
    static void topological_sort(std::vector<DifferentialExtension>& tower);
};

#endif // DIFFERENTIAL_FIELD_H

#ifndef SYMBOLIC_SUM_H
#define SYMBOLIC_SUM_H

#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_expression_internal.h"

#include <string>
#include <vector>
#include <optional>

namespace symbolic_sum {

// ============================================================================
// 求和类型枚举
// ============================================================================

enum class TermType {
    kPolynomial,      // 多项式项
    kGeometric,       // 几何级数项 a*r^n
    kHypergeometric,  // 超几何项
    kTelescoping,     // 伸缩级数项 f(n+1) - f(n)
    kExponential,     // 指数相关
    kTrigonometric,   // 三角相关
    kRational,        // 有理函数
    kUnknown
};

enum class SumType {
    kFinite,          // 有限求和
    kInfiniteUpper,   // 上限无穷
    kInfiniteBoth,    // 双边无穷
    kDivergent        // 发散
};

// ============================================================================
// 求和结果结构
// ============================================================================

struct SumResult {
    SymbolicExpression value;
    bool is_convergent = true;    // 是否收敛（无穷级数）
    bool is_elementary = true;    // 是否为初等函数
    SumType sum_type = SumType::kFinite;
    TermType term_type = TermType::kUnknown;
    std::string method_used;

    static SumResult elementary(const SymbolicExpression& val, TermType tt, const std::string& method) {
        SumResult r;
        r.value = val;
        r.term_type = tt;
        r.method_used = method;
        return r;
    }

    static SumResult divergent(const std::string& reason) {
        SumResult r;
        r.is_convergent = false;
        r.sum_type = SumType::kDivergent;
        r.method_used = "divergent: " + reason;
        return r;
    }

    static SumResult special_function(const std::string& name, const SymbolicExpression& arg) {
        SumResult r;
        r.value = symbolic_expression_internal::make_function(name, arg);
        r.is_elementary = false;
        r.method_used = "special_function: " + name;
        return r;
    }

    static SumResult unknown() {
        SumResult r;
        r.is_convergent = false;
        r.method_used = "unknown";
        return r;
    }
};

// ============================================================================
// 符号求和引擎
// ============================================================================

class SymbolicSumEngine {
public:
    /**
     * @brief 计算符号求和
     * @param term 求和项
     * @param var 求和变量
     * @param lower 下限
     * @param upper 上限
     * @return 求和结果
     */
    SumResult compute_sum(
        const SymbolicExpression& term,
        const std::string& var,
        const BoundArgument& lower,
        const BoundArgument& upper);

    /**
     * @brief 计算符号求和（字符串接口）
     */
    static std::optional<SymbolicExpression> sum(
        const std::string& term_str,
        const std::string& var,
        const BoundArgument& lower,
        const BoundArgument& upper);

private:
    // 分类求和项类型
    TermType classify_term(
        const SymbolicExpression& term,
        const std::string& var);

    // 多项式求和 (Faulhaber 公式)
    SumResult sum_polynomial(
        const SymbolicExpression& term,
        const std::string& var,
        int lower,
        int upper);

    // 几何级数求和
    SumResult sum_geometric(
        const SymbolicExpression& term,
        const std::string& var,
        const BoundArgument& lower,
        const BoundArgument& upper);

    // 有理函数求和
    SumResult sum_rational(
        const SymbolicExpression& term,
        const std::string& var,
        int lower,
        int upper);

    // Gosper 算法（超几何级数）
    bool gosper_algorithm(
        const SymbolicExpression& term,
        const std::string& var,
        int lower,
        int upper,
        SumResult* result);

    // 检测伸缩级数
    bool detect_telescoping(
        const SymbolicExpression& term,
        const std::string& var,
        SymbolicExpression* f);

    // 已知无穷级数
    SumResult known_infinite_series(
        const SymbolicExpression& term,
        const std::string& var);

    // Faulhaber 公式实现
    SymbolicExpression faulhaber_sum(int m, const SymbolicExpression& n);

    // 获取 Bernoulli 数
    std::pair<long long, long long> get_bernoulli(int n);

    // 检查是否为几何级数项
    bool is_geometric_term(
        const SymbolicExpression& term,
        const std::string& var,
        SymbolicExpression* a,
        SymbolicExpression* r);

    // 提取多项式幂次
    int get_polynomial_degree(
        const SymbolicExpression& term,
        const std::string& var);

    // 计算有限数值求和（回退）
    SymbolicExpression numerical_sum(
        const SymbolicExpression& term,
        const std::string& var,
        int lower,
        int upper);
};

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 解析求和参数
 */
bool parse_sum_arguments(
    const std::vector<std::string>& args,
    SymbolicExpression* term,
    std::string* var,
    BoundArgument* lower,
    BoundArgument* upper);

/**
 * @brief 格式化求和结果
 */
std::string format_sum_result(const SumResult& result);

}  // namespace symbolic_sum

#endif  // SYMBOLIC_SUM_H

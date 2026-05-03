#ifndef RISCH_TYPES_H
#define RISCH_TYPES_H

#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_polynomial.h"
#include <string>
#include <set>
#include <vector>

// 最大递归深度限制
#define RISCH_MAX_RECURSION_DEPTH 100

/**
 * @file risch_types.h
 * @brief Risch积分算法的类型定义
 *
 * 这个头文件定义了Risch积分算法所需的所有数据类型，
 * 包括特殊函数类型、微分扩展、积分结果等。
 */

// 特殊函数类型
enum class SpecialFunction {
    kEi,      // 指数积分 Ei(x)
    kErf,     // 误差函数 erf(x)
    kSi,      // 正弦积分 Si(x)
    kCi,      // 余弦积分 Ci(x)
    kLi,      // 对数积分 li(x)
    kGamma,   // Gamma 函数
    kPolyLog  // 多对数函数
};

// 复数根表示，用于 Rothstein-Trager 算法
struct ComplexRoot {
    SymbolicExpression real_part;      // 实部
    SymbolicExpression imag_part;      // 虚部
    bool is_complex;                   // true 如果有非零虚部
    bool is_conjugate_pair;            // true 如果表示共轭对 a+bi 和 a-bi

    static ComplexRoot real(const SymbolicExpression& r) {
        return {r, SymbolicExpression::number(0.0), false, false};
    }
    static ComplexRoot complex(const SymbolicExpression& re, const SymbolicExpression& im, bool conjugate = true) {
        return {re, im, true, conjugate};
    }
};

// 微分扩展类型
struct DifferentialExtension {
    enum class Kind { kLogarithmic, kExponential, kAlgebraic, kTrigonometric, kNone };
    Kind kind;
    SymbolicExpression argument;
    SymbolicExpression derivation;
    std::string t_name;
    int dependency_depth;
    std::set<std::string> dependencies;
    std::string original_function_name;
};

// 积分类型
enum class IntegralType {
    kElementary,
    kNonElementary,
    kSpecialFunction,
    kUnknown
};

// 积分结果
struct RischIntegrationResult {
    bool success;
    IntegralType type;
    SymbolicExpression value;
    std::string message;
    SpecialFunction special_func;

    static RischIntegrationResult elementary(const SymbolicExpression& value) {
        return {true, IntegralType::kElementary, value, "", SpecialFunction::kEi};
    }
    static RischIntegrationResult non_elementary(const std::string& msg = "") {
        return {false, IntegralType::kNonElementary, SymbolicExpression::number(0.0), msg, SpecialFunction::kEi};
    }
    static RischIntegrationResult special_function(SpecialFunction func, const SymbolicExpression& arg) {
        return {true, IntegralType::kSpecialFunction, arg, "", func};
    }
    static RischIntegrationResult unknown(const std::string& msg = "") {
        return {false, IntegralType::kUnknown, SymbolicExpression::number(0.0), msg, SpecialFunction::kEi};
    }
};

// RDE 解的类型
struct RDESolution {
    bool has_polynomial_part;
    SymbolicPolynomial polynomial_part;
    bool has_logarithmic_part;
    SymbolicExpression logarithmic_part;
    bool has_special_part;
    SpecialFunction special_type;
    SymbolicExpression special_arg;

    bool is_valid() const { return has_polynomial_part || has_logarithmic_part || has_special_part; }
};

// 积分缓存键
struct CacheKey {
    std::string expression_str;
    std::string variable;

    bool operator==(const CacheKey& other) const {
        return expression_str == other.expression_str && variable == other.variable;
    }
};

// 缓存键哈希函数
struct CacheKeyHash {
    std::size_t operator()(const CacheKey& k) const {
        return std::hash<std::string>()(k.expression_str) ^ std::hash<std::string>()(k.variable);
    }
};

#endif // RISCH_TYPES_H

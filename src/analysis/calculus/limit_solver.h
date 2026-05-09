#ifndef ANALYSIS_NUMERIC_LIMIT_SOLVER_H
#define ANALYSIS_NUMERIC_LIMIT_SOLVER_H

#include "core/common/scalar_type.h"
#include "symbolic/core/symbolic_expression.h"
#include "analysis/calculus/numerical_calculus.h"
#include <functional>
#include <string>

namespace numeric {

// ============================================================================
// 符号极限探测
// ============================================================================

enum class SymbolicLimitProbeKind {
    kFinite,
    kPositiveInfinity,
    kNegativeInfinity,
    kUnknown,
};

template <typename T>
SymbolicLimitProbeKind probe_symbolic_value_at(
    SymbolicExpression expression,
    const std::string& variable_name,
    T point,
    T* finite_value);

template <typename T>
bool is_infinite_probe(SymbolicLimitProbeKind kind);

template <typename T>
bool is_zero_probe(SymbolicLimitProbeKind kind, T value);

// ============================================================================
// 洛必达法则
// ============================================================================

template <typename T>
bool try_symbolic_lhopital_limit(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    T point,
    int direction,
    T* result,
    std::function<Scalar(const SymbolicExpression&, const std::string&, Scalar)> evaluate_at_override = nullptr);

template <typename T>
bool try_symbolic_one_to_infinity_limit(
    const SymbolicExpression& base,
    const SymbolicExpression& exponent,
    const std::string& variable_name,
    T point,
    T* result);

// ============================================================================
// 极点极限处理
// ============================================================================

template <typename T>
T handle_pole_limit(int shift, T leading_coefficient, int direction);

// ============================================================================
// 无穷远处极限
// ============================================================================

template <typename T>
bool symbolic_limit_at_infinity(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    bool positive,
    T* result);

// ============================================================================
// 数值极限计算
// ============================================================================

template <typename T>
T compute_numerical_limit(
    const std::function<T(T)>& evaluate,
    T x,
    int direction);

}  // namespace numeric

#endif  // ANALYSIS_NUMERIC_LIMIT_SOLVER_H
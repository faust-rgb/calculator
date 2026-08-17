/**
 * @file limit_solver.h
 * @brief 极限求解器（直接基于 Scalar 类型）
 *
 * 本文件定义了极限求解算法：
 * - 数值极限：逐步逼近法计算极限
 * - 符号极限探测：识别无穷极限和振荡极限
 * - L'Hôpital 法则：处理不定式极限
 */

#ifndef ANALYSIS_NUMERIC_LIMIT_SOLVER_H
#define ANALYSIS_NUMERIC_LIMIT_SOLVER_H

#include "types/scalar_type.h"
#include "symbolic/core/symbolic_expression.h"
#include <functional>
#include <string>

namespace numeric {

using Scalar = mymath::Scalar;

// ============================================================================
// 符号极限探测
// ============================================================================

enum class SymbolicLimitProbeKind {
    kFinite,
    kPositiveInfinity,
    kNegativeInfinity,
    kUnknown,
};

SymbolicLimitProbeKind probe_symbolic_value_at(
    SymbolicExpression expression,
    const std::string& variable_name,
    Scalar point,
    Scalar* finite_value);

bool is_infinite_probe(SymbolicLimitProbeKind kind);

bool is_zero_probe(SymbolicLimitProbeKind kind, Scalar value);

// ============================================================================
// 洛必达法则与特殊极限
// ============================================================================

bool try_symbolic_lhopital_limit(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    Scalar point,
    int direction,
    Scalar* result,
    std::function<Scalar(const SymbolicExpression&, const std::string&, Scalar)> evaluate_at_override = nullptr);

bool try_symbolic_one_to_infinity_limit(
    const SymbolicExpression& base,
    const SymbolicExpression& exponent,
    const std::string& variable_name,
    Scalar point,
    Scalar* result);

// ============================================================================
// 极点与无穷远处极限
// ============================================================================

Scalar handle_pole_limit(int shift, Scalar leading_coefficient, int direction);

bool symbolic_limit_at_infinity(
    const SymbolicExpression& expression,
    const std::string& variable_name,
    bool positive,
    Scalar* result);

// ============================================================================
// 数值极限计算
// ============================================================================

Scalar compute_numerical_limit(
    const std::function<Scalar(Scalar)>& evaluate,
    Scalar x,
    int direction);

}  // namespace numeric

#endif  // ANALYSIS_NUMERIC_LIMIT_SOLVER_H
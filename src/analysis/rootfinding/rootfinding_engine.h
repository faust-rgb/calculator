// ============================================================================
// 求根算法引擎 (基于 Scalar 类型)
// ============================================================================
//
// 本文件提供纯粹的数值求根算法实现：
// - Newton 法（带回溯）
// - 二分法
// - 割线法
// - 不动点迭代
// - Brent 法
// ============================================================================

#ifndef ROOTFINDING_ENGINE_H
#define ROOTFINDING_ENGINE_H

#include "types/scalar_type.h"
#include "math/functions/scalar/basic_ops.h"
#include "math/functions/scalar/basic_ops.h"
#include "math/numeric/precision/tolerances.h"
#include "math/mymath.h"

#include <functional>
#include <string>
#include <vector>

namespace rootfinding_engine {

using Scalar = mymath::Scalar;

// ============================================================================
// 容差计算函数
// ============================================================================

inline Scalar root_function_tolerance(Scalar fx) {
    return precision::newton_tolerance<Scalar>() * std::max(Scalar(1), mymath::abs(fx));
}

inline Scalar root_position_tolerance(Scalar x) {
    return precision::default_absolute_tolerance<Scalar>() * std::max(Scalar(1), mymath::abs(x));
}

inline Scalar root_derivative_step(Scalar x) {
    return precision::optimal_derivative_step<Scalar>(x);
}

constexpr int root_max_iterations() {
    return 300;
}

// ============================================================================
// 求根算法接口
// ============================================================================

Scalar newton_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar initial,
    const std::function<Scalar(Scalar)>& normalize,
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate_derivative = nullptr,
    const std::string& variable_name = "x");

Scalar bisection_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar left,
    Scalar right,
    const std::function<Scalar(Scalar)>& normalize,
    const std::string& variable_name = "x");

Scalar secant_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar x0,
    Scalar x1,
    const std::function<Scalar(Scalar)>& normalize,
    const std::string& variable_name = "x");

Scalar fixed_point_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar initial,
    const std::function<Scalar(Scalar)>& normalize,
    const std::string& variable_name = "x");

Scalar brent_solve(
    const std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>& evaluate,
    Scalar left,
    Scalar right,
    const std::function<Scalar(Scalar)>& normalize,
    const std::string& variable_name = "x");

}  // namespace rootfinding_engine

#endif  // ROOTFINDING_ENGINE_H

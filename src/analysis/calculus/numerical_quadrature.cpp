// ============================================================================
// numerical_quadrature.cpp - 自适应数值积分求解器实现
// ============================================================================

#include "numerical_quadrature.h"
#include "analysis/calculus/numerical_calculus.h"

namespace analysis {

Scalar adaptive_simpson(
    const std::function<Scalar(Scalar)>& func,
    Scalar left,
    Scalar right,
    Scalar eps,
    int max_depth) {
    return numeric::adaptive_simpson_callable(func, left, right, eps, max_depth);
}

Scalar adaptive_gauss_kronrod(
    const std::function<Scalar(Scalar)>& function,
    Scalar left,
    Scalar right,
    Scalar eps,
    int max_depth) {
    return numeric::adaptive_gauss_kronrod_callable(function, left, right, eps, max_depth);
}

Scalar tanh_sinh_quadrature(
    const std::function<Scalar(Scalar)>& func,
    Scalar left,
    Scalar right,
    Scalar eps,
    int max_levels) {
    return numeric::tanh_sinh_quadrature_callable(func, left, right, eps, max_levels);
}

} // namespace analysis

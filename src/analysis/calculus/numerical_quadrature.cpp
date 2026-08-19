// ============================================================================
// numerical_quadrature.cpp - 自适应数值积分求解器实现
// ============================================================================

#include "numerical_quadrature.h"
#include "math/mymath.h"
#include <algorithm>
#include <stdexcept>

namespace analysis {

namespace {


inline void compensated_add(Scalar value, Scalar* sum, Scalar* compensation) {
    const Scalar y = value - *compensation;
    const Scalar t = *sum + y;
    *compensation = (t - *sum) - y;
    *sum = t;
}

inline Scalar compensated_pair_sum(Scalar lhs, Scalar rhs) {
    const Scalar sum = lhs + rhs;
    const Scalar compensation = (sum - lhs) - rhs;
    return sum - compensation;
}

inline Scalar relative_tolerance(Scalar baseline, Scalar scale) {
    const Scalar resolved_scale = std::max(Scalar(1.0L), scale);
    return baseline * resolved_scale;
}

Scalar simpson_rule(const std::function<Scalar(Scalar)>& func, Scalar a, Scalar b) {
    const Scalar h = (b - a) * Scalar(0.5L);
    const Scalar fa = func(a);
    const Scalar fb = func(b);
    const Scalar fc = func((a + b) * Scalar(0.5L));
    return h / Scalar(3.0L) * (fa + Scalar(4.0L) * fc + fb);
}

Scalar adaptive_simpson_recursive(
    const std::function<Scalar(Scalar)>& func,
    Scalar a, Scalar b, Scalar whole, Scalar left, Scalar right, Scalar eps, int depth) {
    const Scalar c = (a + b) * Scalar(0.5L);
    const Scalar combined = left + right;
    const Scalar error = mymath::abs(combined - whole) / Scalar(15.0L);

    const Scalar scale = std::max(Scalar(1.0L), mymath::abs(combined));
    if (depth <= 0 || error <= relative_tolerance(eps, scale)) {
        return combined + (combined - whole) / Scalar(15.0L);
    }

    const Scalar d = (a + c) * Scalar(0.5L);
    const Scalar e = (c + b) * Scalar(0.5L);
    const Scalar left_left = simpson_rule(func, a, d);
    const Scalar left_right = simpson_rule(func, d, c);
    const Scalar right_left = simpson_rule(func, c, e);
    const Scalar right_right = simpson_rule(func, e, b);

    return adaptive_simpson_recursive(func, a, c, left, left_left, left_right, eps * Scalar(0.5L), depth - 1) +
           adaptive_simpson_recursive(func, c, b, right, right_left, right_right, eps * Scalar(0.5L), depth - 1);
}

Scalar gauss_kronrod_15(
    const std::function<Scalar(Scalar)>& function,
    Scalar left,
    Scalar right,
    Scalar* error_estimate) {
    static const Scalar kNodes[] = {
        Scalar(0.9914553711208126),
        Scalar(0.9491079123427585),
        Scalar(0.8648644233597691),
        Scalar(0.7415311855993945),
        Scalar(0.5860872354676911),
        Scalar(0.4058451513773972),
        Scalar(0.2077849550078985),
        Scalar(0.0L),
    };
    static const Scalar kKronrodWeights[] = {
        Scalar(0.02293532201052922),
        Scalar(0.06309209262997855),
        Scalar(0.1047900103222502),
        Scalar(0.1406532597155259),
        Scalar(0.1690047266392679),
        Scalar(0.1903505780647854),
        Scalar(0.2044329400752989),
        Scalar(0.2094821410847278),
    };
    static const Scalar kGaussWeights[] = {
        Scalar(0.0L),
        Scalar(0.1294849661688697),
        Scalar(0.0L),
        Scalar(0.2797053914892767),
        Scalar(0.0L),
        Scalar(0.3818300505051189),
        Scalar(0.0L),
        Scalar(0.4179591836734694),
    };

    const Scalar center = (left + right) * Scalar(0.5L);
    const Scalar half_width = (right - left) * Scalar(0.5L);
    Scalar kronrod_sum = Scalar(0.0L);
    Scalar gauss_sum = Scalar(0.0L);
    Scalar kronrod_compensation = Scalar(0.0L);
    Scalar gauss_compensation = Scalar(0.0L);

    for (int i = 0; i < 8; ++i) {
        if (mymath::is_near_zero(kNodes[i], Scalar(0.0L))) {
            const Scalar value = function(center);
            compensated_add(kKronrodWeights[i] * value,
                            &kronrod_sum,
                            &kronrod_compensation);
            compensated_add(kGaussWeights[i] * value,
                            &gauss_sum,
                            &gauss_compensation);
            continue;
        }

        const Scalar offset = half_width * kNodes[i];
        const Scalar left_value = function(center - offset);
        const Scalar right_value = function(center + offset);
        const Scalar pair_sum = compensated_pair_sum(left_value, right_value);
        compensated_add(kKronrodWeights[i] * pair_sum,
                        &kronrod_sum,
                        &kronrod_compensation);
        compensated_add(kGaussWeights[i] * pair_sum,
                        &gauss_sum,
                        &gauss_compensation);
    }

    const Scalar kronrod = half_width * kronrod_sum;
    const Scalar gauss = half_width * gauss_sum;
    if (!mymath::isfinite(kronrod) || !mymath::isfinite(gauss)) {
        throw std::runtime_error("integral did not converge (non-finite value encountered)");
    }
    *error_estimate = mymath::abs(kronrod - gauss);
    return kronrod;
}

Scalar adaptive_gauss_kronrod_recursive(
    const std::function<Scalar(Scalar)>& function,
    Scalar left,
    Scalar right,
    Scalar eps,
    Scalar whole,
    Scalar error,
    int depth) {
    (void)whole;
    const Scalar center = (left + right) * Scalar(0.5L);
    Scalar left_error = Scalar(0.0L);
    Scalar right_error = Scalar(0.0L);
    const Scalar left_result = gauss_kronrod_15(function, left, center, &left_error);
    const Scalar right_result = gauss_kronrod_15(function, center, right, &right_error);
    const Scalar sub_error = left_error + right_error;
    const Scalar sub_result = left_result + right_result;
    if (!mymath::isfinite(sub_result) || !mymath::isfinite(sub_error)) {
        throw std::runtime_error("integral did not converge (non-finite value encountered)");
    }

    const Scalar scale = std::max(Scalar(1.0L), mymath::abs(sub_result));
    if (depth <= 0 ||
        sub_error <= relative_tolerance(eps, scale) ||
        error <= relative_tolerance(eps, scale)) {
        return sub_result;
    }

    return adaptive_gauss_kronrod_recursive(
               function, left, center, eps * Scalar(0.5L), left_result, left_error, depth - 1) +
           adaptive_gauss_kronrod_recursive(
               function, center, right, eps * Scalar(0.5L), right_result, right_error, depth - 1);
}

} // namespace

Scalar adaptive_simpson(
    const std::function<Scalar(Scalar)>& func,
    Scalar left,
    Scalar right,
    Scalar eps,
    int max_depth) {
    const Scalar c = (left + right) * Scalar(0.5L);
    const Scalar whole = simpson_rule(func, left, right);
    const Scalar left_val = simpson_rule(func, left, c);
    const Scalar right_val = simpson_rule(func, c, right);
    return adaptive_simpson_recursive(func, left, right, whole, left_val, right_val, eps, max_depth);
}

Scalar adaptive_gauss_kronrod(
    const std::function<Scalar(Scalar)>& function,
    Scalar left,
    Scalar right,
    Scalar eps,
    int max_depth) {
    Scalar error = Scalar(0.0L);
    const Scalar whole = gauss_kronrod_15(function, left, right, &error);
    return adaptive_gauss_kronrod_recursive(function, left, right, eps, whole, error, max_depth);
}

} // namespace analysis

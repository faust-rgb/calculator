/**
 * @file hyperbolic.cpp
 * @brief Implementation of hyperbolic functions
 */

#include "hyperbolic.h"
#include "math/core/floating_point.h"
#include <stdexcept>

namespace mymath {

// ============================================================================
// Hyperbolic Functions
// ============================================================================

long double sinh(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (isinf(x)) {
        return x;
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::sinh(Scalar(x));
    return static_cast<long double>(result);
}

long double cosh(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (isinf(x)) {
        return infinity();
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::cosh(Scalar(x));
    return static_cast<long double>(result);
}

long double tanh(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (isinf(x)) {
        return x > 0.0L ? 1.0L : -1.0L;
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::tanh(Scalar(x));
    return static_cast<long double>(result);
}

// ============================================================================
// Inverse Hyperbolic Functions
// ============================================================================

long double asinh(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (isinf(x)) {
        return x;
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::asinh(Scalar(x));
    return static_cast<long double>(result);
}

long double acosh(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (isinf(x)) {
        return x;
    }
    if (x < 1.0L) {
        throw std::domain_error("acosh is only defined for x >= 1");
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::acosh(Scalar(x));
    return static_cast<long double>(result);
}

long double atanh(long double x) {
    if (isnan(x)) {
        return x;
    }
    if (x <= -1.0L || x >= 1.0L) {
        throw std::domain_error("atanh is only defined for values in (-1, 1)");
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::atanh(Scalar(x));
    return static_cast<long double>(result);
}

}  // namespace mymath
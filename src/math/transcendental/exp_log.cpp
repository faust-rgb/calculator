/**
 * @file exp_log.cpp
 * @brief Implementation of exponential and logarithmic functions
 */

#include "exp_log.h"
#include "math/core/constants.h"
#include "math/core/floating_point.h"
#include <stdexcept>

namespace mymath {

// ============================================================================
// Exponential Functions
// ============================================================================

long double exp(long double x) {
    if (x >= kLnDoubleMax) {
        return infinity();
    }
    if (x <= kLnDoubleDenormMin) {
        return 0.0L;
    }

    // Use precise for high precision calculation
    Scalar result = precise::exp(Scalar(x));
    return static_cast<long double>(result);
}

// ============================================================================
// Logarithmic Functions
// ============================================================================

long double ln(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("ln is only defined for positive numbers");
    }

    // Use precise for high precision calculation
    Scalar result = precise::ln(Scalar(x));
    return static_cast<long double>(result);
}

long double log(long double x) {
    return ln(x);
}

long double log1p(long double x) {
    // Use ln(1+x) for high precision calculation
    Scalar result = precise::ln(Scalar(1.0L) + Scalar(x));
    return static_cast<long double>(result);
}

long double log10(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("log10 is only defined for positive numbers");
    }

    // Use precise for high precision calculation
    Scalar result = precise::log10(Scalar(x));
    return static_cast<long double>(result);
}

long double log2(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("log2 is only defined for positive numbers");
    }

    // Use ln(x)/ln(2) for high precision calculation
    Scalar result = precise::ln(Scalar(x)) / precise::ln(Scalar(2.0L));
    return static_cast<long double>(result);
}

}  // namespace mymath
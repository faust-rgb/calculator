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

    // Use scalar version for high precision calculation
    Scalar result = mymath::exp(Scalar(x));
    return static_cast<long double>(result);
}

// ============================================================================
// Logarithmic Functions
// ============================================================================

long double ln(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("ln is only defined for positive numbers");
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::ln(Scalar(x));
    return static_cast<long double>(result);
}

long double log(long double x) {
    return ln(x);
}

long double log1p(long double x) {
    // Use scalar version for high precision calculation
    Scalar result = mymath::log1p(Scalar(x));
    return static_cast<long double>(result);
}

long double log10(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("log10 is only defined for positive numbers");
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::log10(Scalar(x));
    return static_cast<long double>(result);
}

long double log2(long double x) {
    if (x <= 0.0L) {
        throw std::domain_error("log2 is only defined for positive numbers");
    }

    // Use scalar version for high precision calculation
    Scalar result = mymath::log2(Scalar(x));
    return static_cast<long double>(result);
}

}  // namespace mymath
/**
 * @file mymath.h
 * @brief Aggregated header for the modular math library
 *
 * This file includes all math library components for backward compatibility.
 * New code should include specific headers from the appropriate subdirectories:
 * - math/types/ for numeric types (float128, complex, dual)
 * - math/core/ for constants and basic operations
 * - math/transcendental/ for trig, hyperbolic, exp/log functions
 * - math/special/ for gamma, beta, bessel, error functions
 */

#ifndef MYMATH_H
#define MYMATH_H

// Types - numeric type system
#include "app/scalar_type.h"
#include "math/types/float128.h"
#include "math/types/complex.h"
#include "math/types/dual.h"

// Core - constants and basic operations
#include "math/core/constants.h"
#include "math/core/floating_point.h"
#include "math/core/basic_ops.h"
#include "math/core/roots_powers.h"
#include "math/core/scalar_traits.h"

// Transcendental - trig, hyperbolic, exp/log
#include "math/transcendental/transcendental.h"

// Special - gamma, beta, bessel, error functions
#include "math/special/special_functions.h"

#endif // MYMATH_H

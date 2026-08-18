/**
 * @file mymath.h
 * @brief Aggregated header for the modular math library
 *
 * This file includes all math library components for backward compatibility.
 * New code should include specific headers from the appropriate subdirectories:
 * - math/numeric/types/ for numeric types (float128, complex, dual)
 * - math/numeric/scalar/ for Scalar backend dispatch
 * - math/functions/elementary/ for trig, hyperbolic, exp/log functions
 * - math/functions/special/ for gamma, beta, bessel, error functions
 */

#ifndef MYMATH_H
#define MYMATH_H

// Types - numeric type system
#include "types/scalar_type.h"
#include "math/numeric/float128/float128.h"
#include "math/numeric/types/complex.h"
#include "math/numeric/types/dual.h"

// Core - constants and basic operations
#include "math/numeric/constants/numeric.h"
#include "math/numeric/precision/predicates.h"
#include "math/numeric/scalar/dispatch.h"
#include "math/numeric/precision/tolerances.h"
#include "math/runtime/precision/default_precision.h"

// Public operations by numeric domain
#include "math/functions/scalar/basic_ops.h"
#include "math/functions/scalar/basic_ops.h"
#include "math/functions/integer/integer_helpers.h"
#include "math/functions/conversion/rational_approximation.h"

// Transcendental - trig, hyperbolic, exp/log

// Special - gamma, beta, bessel, error functions
#include "math/functions/special/special_functions.h"

#endif // MYMATH_H

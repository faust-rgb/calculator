#ifndef CORE_SCALAR_TYPE_H
#define CORE_SCALAR_TYPE_H

#include "math/types/float128.h"
#include "precise/precise_decimal.h"

namespace mymath {

/**
 * @brief Global scalar type definition
 *
 * Uses PreciseDecimal by default for higher precision.
 * Modify this to switch to other data types.
 */
//using Scalar = PreciseDecimal;
using Scalar = float128_t;

} // namespace mymath

// Export Scalar to global namespace for convenience
using mymath::Scalar;

#endif // CORE_SCALAR_TYPE_H
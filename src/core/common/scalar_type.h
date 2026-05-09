#ifndef CORE_SCALAR_TYPE_H
#define CORE_SCALAR_TYPE_H

#include "math/types/float128.h"

namespace mymath {

/**
 * @brief Global scalar type definition
 *
 * Uses float128_t by default for higher precision.
 * Modify this to switch to other data types.
 */
using Scalar = float128_t;

} // namespace mymath

// Export Scalar to global namespace for convenience
using mymath::Scalar;

#endif // CORE_SCALAR_TYPE_H
#ifndef NUMERIC_CONVERSION_H
#define NUMERIC_CONVERSION_H

#include "number.h"
#include "precision_context.h"

namespace numeric {

// Convert Number to double for interoperability
// Returns best double approximation
double to_double(const Number& value);

// Check if a Number can be exactly represented as a double
// Returns true if the value fits in double range without overflow
bool can_convert_to_double(const Number& value);

// Convert double to Number (exact representation if possible)
Number from_double(double value);

// Get a long double approximation (higher precision than double)
long double to_long_double(const Number& value);

}  // namespace numeric

#endif

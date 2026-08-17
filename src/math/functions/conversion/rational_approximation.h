#pragma once

#include "types/scalar_type.h"

namespace mymath {
bool approximate_fraction(Scalar value, long long* numerator, long long* denominator,
                          int max_denominator = 999, Scalar eps = Scalar(1e-30L));
bool best_rational_approximation(long double value, long long* numerator,
                                 long long* denominator, long long max_denominator = 999);
}

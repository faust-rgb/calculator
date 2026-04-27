#ifndef NUMERIC_FUNCTIONS_H
#define NUMERIC_FUNCTIONS_H

#include "number.h"

namespace numeric {

Number abs(const Number& value);
Number negate(const Number& value);
bool is_zero(const Number& value);

Number sqrt(const Number& value);
Number cbrt(const Number& value);
Number root(const Number& value, const Number& degree);
Number pow(const Number& base, const Number& exponent);

Number pi();
Number exp(const Number& value);
Number ln(const Number& value);
Number sin(const Number& value);
Number cos(const Number& value);
Number tan(const Number& value);
Number asin(const Number& value);
Number acos(const Number& value);
Number atan(const Number& value);
Number sinh(const Number& value);
Number cosh(const Number& value);
Number tanh(const Number& value);
Number asinh(const Number& value);
Number acosh(const Number& value);
Number atanh(const Number& value);
Number sec(const Number& value);
Number csc(const Number& value);
Number cot(const Number& value);
Number asec(const Number& value);
Number acsc(const Number& value);
Number acot(const Number& value);
Number erf(const Number& value);
Number erfc(const Number& value);
Number log10(const Number& value);
Number gamma(const Number& value);
Number beta(const Number& a, const Number& b);
Number zeta(const Number& value);
Number bessel_j(int order, const Number& x);

// Utility functions for migration from mymath
bool is_near_zero(const Number& value);
bool is_near_zero(double value, double eps = 1e-10);
bool is_integer_value(const Number& value);
bool is_integer_value(double value, double eps = 1e-10);
bool approximate_fraction(const Number& value,
                          BigInt* numerator,
                          BigInt* denominator,
                          long long max_denominator);
bool approximate_fraction(double value,
                          long long* numerator,
                          long long* denominator,
                          long long max_denominator,
                          double eps = 1e-10);

}  // namespace numeric

#endif

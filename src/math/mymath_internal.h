#ifndef MYMATH_INTERNAL_H
#define MYMATH_INTERNAL_H

namespace mymath {
namespace internal {

double log_gamma_positive(double x);
double finite_or_infinity_from_log(double log_value);

}  // namespace internal
}  // namespace mymath

#endif

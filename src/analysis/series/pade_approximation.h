#ifndef ANALYSIS_SERIES_PADE_APPROXIMATION_H
#define ANALYSIS_SERIES_PADE_APPROXIMATION_H

#include "core/common/scalar_type.h"
#include <string>
#include <vector>
#include <functional>

namespace series_ops {

struct SeriesContext;

namespace pade {

std::string pade(const SeriesContext& ctx,
                 const std::string& expr,
                 Scalar center,
                 int numerator_degree,
                 int denominator_degree);

std::string pade_from_coeffs(const std::vector<Scalar>& coefficients,
                              int numerator_degree,
                              int denominator_degree);

bool solve_pade_denominator(std::function<Scalar(int)> c,
                             int numerator_degree,
                             int denominator_degree,
                             std::vector<Scalar>& q);

bool solve_tohplitz_stable(std::function<Scalar(int)> c, int n, std::vector<Scalar>& q);

std::string format_pade_result(const std::vector<Scalar>& numerator,
                                const std::vector<Scalar>& denominator);

std::string format_simple_pade(const std::vector<Scalar>& numerator,
                                const std::vector<Scalar>& denominator);

}  // namespace pade
}  // namespace series_ops

#endif  // ANALYSIS_SERIES_PADE_APPROXIMATION_H
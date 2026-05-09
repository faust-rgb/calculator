#ifndef ANALYSIS_SERIES_PUISEUX_SERIES_H
#define ANALYSIS_SERIES_PUISEUX_SERIES_H

#include "core/scalar_type.h"
#include <string>
#include <vector>
#include <tuple>

namespace series_ops {

struct SeriesContext;

namespace puiseux {

bool compute_newton_polygon(const std::vector<int>& x_powers,
                            const std::vector<int>& y_powers,
                            const std::vector<Scalar>& coefficients,
                            std::vector<std::tuple<Scalar, int, int, int, int>>* edges);

bool newton_puiseux_expand(const std::vector<Scalar>& poly_coeffs,
                           int degree,
                           std::vector<std::pair<Scalar, int>>* result);

std::string puiseux(const SeriesContext& ctx,
                    const std::string& expr,
                    Scalar center,
                    int degree,
                    int denominator);

std::string general_series_auto(const SeriesContext& ctx,
                              const std::string& expression_str,
                              Scalar center,
                              int degree);

}  // namespace puiseux
}  // namespace series_ops

#endif  // ANALYSIS_SERIES_PUISEUX_SERIES_H
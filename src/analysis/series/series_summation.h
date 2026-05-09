#ifndef ANALYSIS_SERIES_SERIES_SUMMATION_H
#define ANALYSIS_SERIES_SERIES_SUMMATION_H

#include "core/common/scalar_type.h"
#include "symbolic/core/symbolic_expression.h"
#include <string>
#include <vector>
#include <functional>

namespace series_ops {

struct SeriesContext;

namespace summation {

std::string series_sum(const SeriesContext& ctx,
                       const std::string& expr,
                       const std::string& index_name,
                       const std::string& lower,
                       const std::string& upper);

bool detect_geometric_ratio_symbolic(const SymbolicExpression& summand,
                                      const std::string& index_name,
                                      Scalar* coefficient,
                                      Scalar* ratio);

bool detect_arith_geo_series(const SymbolicExpression& summand,
                              const std::string& index_name,
                              std::vector<Scalar>* poly_coeffs,
                              Scalar* ratio);

std::string compute_power_times_geo_sum(int m, const std::string& r,
                                         const std::string& index_name);

std::string arith_geo_sum(const std::vector<Scalar>& poly_coeffs,
                           Scalar ratio,
                           const std::string& index_name,
                           const SymbolicExpression& upper_expr,
                           const std::string& lower);

bool detect_telescoping(const SymbolicExpression& summand,
                        const std::string& index_name,
                        SymbolicExpression* g_function);

bool wynn_epsilon_acceleration(const std::vector<Scalar>& partial_sums,
                                int max_iterations,
                                Scalar tolerance,
                                Scalar* result);

bool euler_transform(const std::vector<Scalar>& terms,
                     Scalar tolerance,
                     Scalar* result);

bool levin_transform(const std::vector<Scalar>& partial_sums,
                     const std::vector<Scalar>& terms,
                     Scalar tolerance,
                     Scalar* result);

bool numerical_infinite_sum(const SeriesContext& ctx,
                            const SymbolicExpression& summand,
                            const std::string& index_name,
                            const std::string& lower,
                            int max_terms,
                            Scalar* result);

bool try_infinite_series_value(const SymbolicExpression& summand,
                                const std::string& index_name,
                                const std::string& lower,
                                std::string* result);

}  // namespace summation
}  // namespace series_ops

#endif  // ANALYSIS_SERIES_SERIES_SUMMATION_H
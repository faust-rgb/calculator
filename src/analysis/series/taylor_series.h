#ifndef ANALYSIS_SERIES_TAYLOR_SERIES_H
#define ANALYSIS_SERIES_TAYLOR_SERIES_H

#include "core/common/scalar_type.h"
#include "symbolic/core/symbolic_expression.h"
#include "analysis/series/psa_engine.h"
#include <string>
#include <vector>
#include <functional>

namespace series_ops {

struct SeriesContext;

namespace taylor {

std::vector<Scalar> build_taylor_coefficients(
    const SeriesContext& ctx,
    const SymbolicExpression& expression,
    const std::string& variable_name,
    Scalar center,
    int degree);

std::string taylor(const SeriesContext& ctx,
                   const std::string& expr,
                   Scalar center,
                   int degree);

bool compute_taylor_coefficients_ad(const SeriesContext& ctx,
                                    const SymbolicExpression& expression,
                                    const std::string& variable_name,
                                    Scalar center,
                                    int degree,
                                    std::vector<Scalar>& result);

namespace predefined_series {

bool try_predefined_taylor(const SymbolicExpression& expr,
                            const std::string& var_name,
                            Scalar center,
                            int degree,
                            std::vector<Scalar>& result);

}  // namespace predefined_series

}  // namespace taylor
}  // namespace series_ops

#endif  // ANALYSIS_SERIES_TAYLOR_SERIES_H

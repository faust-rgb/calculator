/**
 * @file series_summation.h
 * @brief 级数求和
 *
 * 本文件定义了符号级数求和算法：
 * - 几何级数求和：检测等比数列并计算闭式解
 * - 调和级数：识别常见的调和级数形式
 * - 符号求和：使用符号方法计算级数和
 */

#ifndef ANALYSIS_SERIES_SERIES_SUMMATION_H
#define ANALYSIS_SERIES_SERIES_SUMMATION_H

#include "app/scalar_type.h"
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
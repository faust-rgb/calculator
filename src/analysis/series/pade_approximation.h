/**
 * @file pade_approximation.h
 * @brief Padé 逼近
 *
 * 本文件定义了 Padé 逼近算法：
 * - 有理逼近：用有理函数逼近幂级数
 * - 分子分母计算：求解 Padé 系数
 * - 收敛加速：比 Taylor 级数收敛更快
 */

#ifndef ANALYSIS_SERIES_PADE_APPROXIMATION_H
#define ANALYSIS_SERIES_PADE_APPROXIMATION_H

#include "types/scalar_type.h"
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
#ifndef ANALYSIS_SERIES_PSA_ENGINE_H
#define ANALYSIS_SERIES_PSA_ENGINE_H

#include "core/common/scalar_type.h"
#include "symbolic/core/symbolic_expression.h"
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>

namespace series_ops {

using Scalar = mymath::Scalar;

struct SeriesContext;

namespace internal {

struct PoleException : public std::runtime_error {
    int shift;
    Scalar leading_coefficient;
    PoleException(int s, Scalar coeff);
};

struct LaurentSeries {
    std::vector<Scalar> coefficients;
    int low_power;
    int high_power;

    LaurentSeries();
    LaurentSeries(int low, int high, const std::vector<Scalar>& coeffs);
    bool empty() const;
    Scalar get_coefficient(int power) const;
    bool to_taylor(std::vector<Scalar>& result, int degree) const;
};

std::vector<Scalar> ps_add(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree);
std::vector<Scalar> ps_sub(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree);
std::vector<Scalar> ps_mul(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree);
std::vector<Scalar> ps_div_with_laurent(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree, int* laurent_shift);
std::vector<Scalar> ps_div(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree);
std::vector<Scalar> ps_exp(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_ln(const std::vector<Scalar>& a, int degree);
void ps_sincos(const std::vector<Scalar>& a, int degree, std::vector<Scalar>& sin_res, std::vector<Scalar>& cos_res);
std::vector<Scalar> ps_sin(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_cos(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_pow_const(const std::vector<Scalar>& a, Scalar n, int degree);
std::vector<Scalar> ps_scale(const std::vector<Scalar>& a, Scalar scale, int degree);
std::vector<Scalar> ps_derivative(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_integral(const std::vector<Scalar>& a, Scalar constant_term, int degree);
std::vector<Scalar> ps_tan(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_asin(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_acos(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_atan(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_sinh(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_cosh(const std::vector<Scalar>& a, int degree);
std::vector<Scalar> ps_tanh(const std::vector<Scalar>& a, int degree);

bool evaluate_psa(const SymbolicExpression& expr, const std::string& var_name, Scalar center, int degree, std::vector<Scalar>& result, const SeriesContext& ctx);

}  // namespace internal

}  // namespace series_ops

#endif  // ANALYSIS_SERIES_PSA_ENGINE_H

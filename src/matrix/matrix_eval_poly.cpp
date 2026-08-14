// ============================================================================
// matrix_eval_poly.cpp - 矩阵/向量多项式与插值函数求值实现
// ============================================================================

#include "matrix_eval_poly.h"
#include "matrix_internal.h"
#include "polynomial/polynomial.h"
#include "mymath.h"
#include <stdexcept>

namespace matrix {
namespace internal {

bool try_evaluate_poly_function(
    const std::string& name,
    const std::vector<std::string>& arguments,
    const ScalarEvaluator& scalar_eval,
    const std::function<Matrix(const std::string&, const std::string&)>& require_matrix,
    Value* result) {
    if (!result) return false;

    if (name == "lagrange") {
        if (arguments.size() != 3) {
            throw std::runtime_error("lagrange expects x samples, y samples, and xi");
        }
        *result = Value::from_scalar(lagrange_interpolate(
            as_vector_values(require_matrix(arguments[0], "lagrange"), "lagrange"),
            as_vector_values(require_matrix(arguments[1], "lagrange"), "lagrange"),
            scalar_eval(arguments[2])));
        return true;
    }

    if (name == "spline") {
        if (arguments.size() != 3) {
            throw std::runtime_error("spline expects x samples, y samples, and xi");
        }
        *result = Value::from_scalar(spline_interpolate(
            as_vector_values(require_matrix(arguments[0], "spline"), "spline"),
            as_vector_values(require_matrix(arguments[1], "spline"), "spline"),
            scalar_eval(arguments[2])));
        return true;
    }

    if (name == "linear_regression") {
        if (arguments.size() != 2) {
            throw std::runtime_error("linear_regression expects exactly two vector arguments");
        }
        const auto fit = linear_regression_fit(
            as_vector_values(require_matrix(arguments[0], "linear_regression"),
                             "linear_regression"),
            as_vector_values(require_matrix(arguments[1], "linear_regression"),
                             "linear_regression"));
        *result = Value::from_matrix(Matrix::vector({fit.first, fit.second}));
        return true;
    }

    if (name == "poly_fit" || name == "polynomial_fit") {
        if (arguments.size() != 3) {
            throw std::runtime_error(name + " expects x samples, y samples, and degree");
        }
        const Scalar degree_value = scalar_eval(arguments[2]);
        if (!mymath::is_integer(degree_value) || degree_value < 0.0L) {
            throw std::runtime_error(name + " degree must be a non-negative integer");
        }
        *result = Value::from_matrix(Matrix::vector(polynomial_fit(
            as_vector_values(require_matrix(arguments[0], name), name),
            as_vector_values(require_matrix(arguments[1], name), name),
            static_cast<int>(degree_value + 0.5))));
        return true;
    }

    if (name == "poly_eval") {
        if (arguments.size() != 2) {
            throw std::runtime_error("poly_eval expects coefficient vector and x");
        }
        *result = Value::from_scalar(polynomial_evaluate(
            as_vector_values(require_matrix(arguments[0], "poly_eval"), "poly_eval"),
            scalar_eval(arguments[1])));
        return true;
    }

    if (name == "poly_deriv") {
        if (arguments.size() != 1) {
            throw std::runtime_error("poly_deriv expects exactly one coefficient vector");
        }
        *result = Value::from_matrix(Matrix::vector(polynomial_derivative(
            as_vector_values(require_matrix(arguments[0], "poly_deriv"), "poly_deriv"))));
        return true;
    }

    if (name == "poly_integ") {
        if (arguments.size() != 1) {
            throw std::runtime_error("poly_integ expects exactly one coefficient vector");
        }
        *result = Value::from_matrix(Matrix::vector(polynomial_integral(
            as_vector_values(require_matrix(arguments[0], "poly_integ"), "poly_integ"))));
        return true;
    }

    if (name == "poly_compose") {
        if (arguments.size() != 2) {
            throw std::runtime_error("poly_compose expects exactly two coefficient vectors");
        }
        *result = Value::from_matrix(Matrix::vector(polynomial_compose(
            as_vector_values(require_matrix(arguments[0], "poly_compose"), "poly_compose"),
            as_vector_values(require_matrix(arguments[1], "poly_compose"), "poly_compose"))));
        return true;
    }

    if (name == "poly_gcd") {
        if (arguments.size() != 2) {
            throw std::runtime_error("poly_gcd expects exactly two coefficient vectors");
        }
        *result = Value::from_matrix(Matrix::vector(polynomial_gcd(
            as_vector_values(require_matrix(arguments[0], "poly_gcd"), "poly_gcd"),
            as_vector_values(require_matrix(arguments[1], "poly_gcd"), "poly_gcd"))));
        return true;
    }

    return false;
}

} // namespace internal
} // namespace matrix

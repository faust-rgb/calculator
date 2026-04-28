/**
 * @file polynomial.cpp
 * @brief 多项式运算实现
 */

#include "polynomial.h"

#include "matrix.h"
#include "mymath.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace {

/** @brief 多项式运算的数值精度阈值 */
constexpr double kPolynomialEps = 1e-10;

/** @brief 根隔离过程的精度阈值 */
constexpr double kRootIsolationEps = 1e-8;

void trim_trailing_zeros(std::vector<double>* coefficients) {
    while (coefficients->size() > 1 &&
           mymath::is_near_zero(coefficients->back(), kPolynomialEps)) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(0.0);
    }
}

std::string format_coefficient(double value) {
    if (mymath::is_integer(value, 1e-10)) {
        long long rounded =
            static_cast<long long>(value >= 0.0 ? value + 0.5 : value - 0.5);
        return std::to_string(rounded);
    }

    long long numerator = 0;
    long long denominator = 1;
    if (mymath::approximate_fraction(value,
                                     &numerator,
                                     &denominator,
                                     999,
                                     1e-10)) {
        if (value < 0.0) {
            numerator = -numerator;
        }
        if (denominator == 1) {
            return std::to_string(numerator);
        }
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }

    std::string text = std::to_string(value);
    while (!text.empty() && text.back() == '0') {
        text.pop_back();
    }
    if (!text.empty() && text.back() == '.') {
        text.pop_back();
    }
    return text;
}

double polynomial_evaluate_impl(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    for (std::size_t i = coefficients.size(); i > 0; --i) {
        result = result * x + coefficients[i - 1];
    }
    return result;
}

std::vector<double> polynomial_derivative_impl(const std::vector<double>& coefficients) {
    if (coefficients.size() <= 1) {
        return {0.0};
    }
    std::vector<double> derivative(coefficients.size() - 1, 0.0);
    for (std::size_t i = 1; i < coefficients.size(); ++i) {
        derivative[i - 1] = coefficients[i] * static_cast<double>(i);
    }
    trim_trailing_zeros(&derivative);
    return derivative;
}

double polynomial_root_bound(const std::vector<double>& coefficients) {
    const double leading = coefficients.back();
    double bound = 0.0;
    for (std::size_t i = 0; i + 1 < coefficients.size(); ++i) {
        const double ratio = mymath::abs(coefficients[i] / leading);
        if (ratio > bound) bound = ratio;
    }
    return 1.0 + bound;
}

double bisect_root(const std::vector<double>& coefficients, double left, double right) {
    double left_value = polynomial_evaluate_impl(coefficients, left);
    for (int i = 0; i < 100; ++i) {
        const double mid = (left + right) * 0.5;
        const double mid_value = polynomial_evaluate_impl(coefficients, mid);
        if (mymath::is_near_zero(mid_value, kRootIsolationEps) ||
            mymath::abs(right - left) <= kRootIsolationEps) {
            return mid;
        }
        if ((left_value < 0.0) == (mid_value < 0.0)) {
            left = mid;
            left_value = mid_value;
        } else {
            right = mid;
        }
    }
    return (left + right) * 0.5;
}

void add_unique_root(std::vector<double>* roots, double candidate) {
    for (double existing : *roots) {
        if (mymath::abs(existing - candidate) <= 1e-6) return;
    }
    roots->push_back(candidate);
}

}  // namespace

double polynomial_evaluate(const std::vector<double>& coefficients, double x) {
    return polynomial_evaluate_impl(coefficients, x);
}

std::vector<double> polynomial_derivative(const std::vector<double>& coefficients) {
    return polynomial_derivative_impl(coefficients);
}

std::vector<double> polynomial_add(const std::vector<double>& lhs,
                                   const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) result[i] += lhs[i];
    for (std::size_t i = 0; i < rhs.size(); ++i) result[i] += rhs[i];
    trim_trailing_zeros(&result);
    return result;
}

std::vector<double> polynomial_subtract(const std::vector<double>& lhs,
                                        const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) result[i] += lhs[i];
    for (std::size_t i = 0; i < rhs.size(); ++i) result[i] -= rhs[i];
    trim_trailing_zeros(&result);
    return result;
}

std::vector<double> polynomial_multiply(const std::vector<double>& lhs,
                                        const std::vector<double>& rhs) {
    std::vector<double> result(lhs.size() + rhs.size() - 1, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result[i + j] += lhs[i] * rhs[j];
        }
    }
    trim_trailing_zeros(&result);
    return result;
}

PolynomialDivisionResult polynomial_divide(const std::vector<double>& dividend,
                                           const std::vector<double>& divisor) {
    std::vector<double> normalized_dividend = dividend;
    std::vector<double> normalized_divisor = divisor;
    trim_trailing_zeros(&normalized_dividend);
    trim_trailing_zeros(&normalized_divisor);

    if (normalized_divisor.size() == 1 &&
        mymath::is_near_zero(normalized_divisor[0], kPolynomialEps)) {
        throw std::runtime_error("polynomial divisor cannot be zero");
    }

    if (normalized_dividend.size() < normalized_divisor.size()) {
        return {{0.0}, normalized_dividend};
    }

    std::vector<double> quotient(
        normalized_dividend.size() - normalized_divisor.size() + 1, 0.0);
    std::vector<double> remainder = normalized_dividend;

    while (remainder.size() >= normalized_divisor.size() &&
           !(remainder.size() == 1 && mymath::is_near_zero(remainder[0], kPolynomialEps))) {
        const std::size_t degree_diff = remainder.size() - normalized_divisor.size();
        const double factor = remainder.back() / normalized_divisor.back();
        quotient[degree_diff] = factor;
        for (std::size_t i = 0; i < normalized_divisor.size(); ++i) {
            remainder[degree_diff + i] -= factor * normalized_divisor[i];
        }
        trim_trailing_zeros(&remainder);
    }

    trim_trailing_zeros(&quotient);
    trim_trailing_zeros(&remainder);
    return {quotient, remainder};
}

std::vector<double> polynomial_real_roots(const std::vector<double>& coefficients) {
    std::vector<double> normalized = coefficients;
    trim_trailing_zeros(&normalized);

    if (normalized.size() <= 1) {
        throw std::runtime_error("constant polynomial does not have isolated roots");
    }

    if (normalized.size() == 2) {
        double r = -normalized[0] / normalized[1];
        long long num, den;
        if (mymath::approximate_fraction(r, &num, &den, 100, 1e-9)) {
            if (r < 0.0) {
                num = -num;
            }
            r = static_cast<double>(num) / static_cast<double>(den);
        }
        return {r};
    }

    const std::vector<double> derivative = polynomial_derivative(normalized);
    std::vector<double> critical_points = polynomial_real_roots(derivative);
    std::sort(critical_points.begin(), critical_points.end());

    const double bound = polynomial_root_bound(normalized);
    std::vector<double> points;
    points.push_back(-bound);
    for (double point : critical_points) points.push_back(point);
    points.push_back(bound);

    std::vector<double> roots;
    for (double point : critical_points) {
        if (mymath::is_near_zero(polynomial_evaluate(normalized, point), 1e-6)) {
            add_unique_root(&roots, point);
        }
    }

    for (std::size_t i = 1; i < points.size(); ++i) {
        const double left = points[i - 1];
        const double right = points[i];
        const double left_value = polynomial_evaluate(normalized, left);
        const double right_value = polynomial_evaluate(normalized, right);

        if (mymath::is_near_zero(left_value, 1e-6)) {
            add_unique_root(&roots, left);
            continue;
        }
        if (mymath::is_near_zero(right_value, 1e-6)) {
            add_unique_root(&roots, right);
            continue;
        }

        if ((left_value < 0.0) != (right_value < 0.0)) {
            double r = bisect_root(normalized, left, right);
            long long num, den;
            if (mymath::approximate_fraction(r, &num, &den, 100, 1e-9)) {
                if (r < 0.0) {
                    num = -num;
                }
                r = static_cast<double>(num) / static_cast<double>(den);
            }
            add_unique_root(&roots, r);
        }
    }

    std::sort(roots.begin(), roots.end());
    return roots;
}

std::vector<double> polynomial_integral(const std::vector<double>& coefficients) {
    std::vector<double> integral(coefficients.size() + 1, 0.0);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        integral[i + 1] = coefficients[i] / static_cast<double>(i + 1);
    }
    trim_trailing_zeros(&integral);
    return integral;
}

std::vector<double> polynomial_compose(const std::vector<double>& outer,
                                       const std::vector<double>& inner) {
    std::vector<double> result = {0.0};
    for (std::size_t i = outer.size(); i > 0; --i) {
        result = polynomial_multiply(result, inner);
        if (result.empty()) result.push_back(0.0);
        result[0] += outer[i - 1];
        trim_trailing_zeros(&result);
    }
    return result;
}

std::vector<double> polynomial_gcd(const std::vector<double>& lhs,
                                   const std::vector<double>& rhs) {
    std::vector<double> a = lhs;
    std::vector<double> b = rhs;
    trim_trailing_zeros(&a);
    trim_trailing_zeros(&b);
    
    // 动态容差：根据输入多项式系数的最大绝对值决定
    auto get_max_coeff = [](const std::vector<double>& p) {
        double m = 0.0;
        for (double c : p) m = std::max(m, mymath::abs(c));
        return m;
    };
    const double scale = std::max(get_max_coeff(a), get_max_coeff(b));
    const double dynamic_eps = std::max(kPolynomialEps, scale * 1e-12);

    while (!(b.size() == 1 && mymath::is_near_zero(b[0], dynamic_eps))) {
        const PolynomialDivisionResult division = polynomial_divide(a, b);
        a = b;
        b = division.remainder;
        trim_trailing_zeros(&a);
        trim_trailing_zeros(&b);
    }
    
    if (a.empty() || (a.size() == 1 && mymath::is_near_zero(a[0], dynamic_eps))) {
        return {0.0};
    }

    const double leading = a.back();
    if (!mymath::is_near_zero(leading, dynamic_eps)) {
        for (double& coefficient : a) coefficient /= leading;
    }
    trim_trailing_zeros(&a);
    return a;
}

std::vector<double> polynomial_fit(const std::vector<double>& x_samples,
                                   const std::vector<double>& y_samples,
                                   int degree) {
    if (degree < 0) throw std::runtime_error("polynomial degree must be non-negative");
    if (x_samples.size() != y_samples.size() || x_samples.empty()) {
        throw std::runtime_error("polynomial_fit requires non-empty sample vectors of the same length");
    }
    const std::size_t m_vars = static_cast<std::size_t>(degree + 1);
    if (x_samples.size() < m_vars) throw std::runtime_error("polynomial_fit requires at least degree + 1 samples");

    const std::size_t n = x_samples.size();
    long double x_sum = 0.0L;
    for (double x : x_samples) x_sum += static_cast<long double>(x);
    const double center = static_cast<double>(x_sum / static_cast<long double>(n));
    double scale = 0.0;
    for (double x : x_samples) {
        const double d = mymath::abs(x - center);
        if (d > scale) scale = d;
    }
    if (scale < 1e-9) scale = 1.0;

    matrix::Matrix A(n, m_vars);
    for (std::size_t i = 0; i < n; ++i) {
        const double sx = (x_samples[i] - center) / scale;
        double p = 1.0;
        for (std::size_t j = 0; j < m_vars; ++j) {
            A.at(i, j) = p;
            p *= sx;
        }
    }
    matrix::Matrix b = matrix::Matrix::vector(y_samples);
    try {
        matrix::Matrix solution = matrix::least_squares(A, b);
        std::vector<double> scaled_coeffs(m_vars);
        for (std::size_t i = 0; i < m_vars; ++i) scaled_coeffs[i] = solution.at(i, 0);
        const std::vector<double> linear_map = {-center / scale, 1.0 / scale};
        std::vector<double> coefficients = polynomial_compose(scaled_coeffs, linear_map);
        trim_trailing_zeros(&coefficients);
        return coefficients;
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("polynomial_fit failed: ") + e.what());
    }
}

std::string polynomial_to_string(const std::vector<double>& coefficients,
                                 const std::string& variable_name) {
    std::vector<double> normalized = coefficients;
    trim_trailing_zeros(&normalized);
    if (normalized.size() == 1 && mymath::is_near_zero(normalized[0], kPolynomialEps)) return "0";
    std::string result;
    bool first = true;
    for (std::size_t index = normalized.size(); index > 0; --index) {
        const std::size_t degree = index - 1;
        const double coefficient = normalized[degree];
        if (mymath::is_near_zero(coefficient, kPolynomialEps)) continue;
        const bool negative = coefficient < 0.0;
        const double abs_value = negative ? -coefficient : coefficient;
        std::string term;
        if (degree == 0) term = format_coefficient(abs_value);
        else {
            if (!mymath::is_near_zero(abs_value - 1.0, kPolynomialEps)) term += format_coefficient(abs_value) + " * ";
            term += variable_name;
            if (degree > 1) term += " ^ " + std::to_string(degree);
        }
        if (first) {
            result += negative ? "-" + term : term;
            first = false;
        } else result += negative ? " - " + term : " + " + term;
    }
    return result.empty() ? "0" : result;
}

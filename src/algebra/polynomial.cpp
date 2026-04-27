/**
 * @file polynomial.cpp
 * @brief 多项式运算实现（高精度版本）
 *
 * 实现多项式的基本运算（加减乘除）和实根查找。
 * 使用系数向量表示多项式，索引 i 对应 x^i 的系数。
 * 使用 numeric::Number 实现任意精度计算。
 */

#include "polynomial.h"

#include "functions.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace {

using numeric::Number;
using numeric::BigInt;
using numeric::BigDecimal;
using numeric::PrecisionContext;

PrecisionContext default_polynomial_context() {
    PrecisionContext ctx;
    ctx.digits = 50;
    ctx.max_iterations = 1000;
    return ctx;
}

bool is_number_zero(const Number& value) {
    return numeric::is_zero(value);
}

bool is_number_one(const Number& value) {
    if (value.is_integer()) {
        const BigInt& i = value.as_integer();
        return !i.is_zero() && (i - BigInt(1)).is_zero();
    }
    if (value.is_rational()) {
        const auto& r = value.as_rational();
        return r.numerator() == BigInt(1) && r.denominator() == BigInt(1);
    }
    if (value.is_decimal()) {
        const auto& d = value.as_decimal();
        BigInt remainder = d.coefficient().abs() % BigDecimal::pow10(d.scale());
        return remainder.is_zero() && d.coefficient() == BigInt(1);
    }
    return false;
}

Number number_abs(const Number& value) {
    return numeric::abs(value);
}

Number number_sqrt(const Number& value) {
    return numeric::sqrt(value);
}

void trim_trailing_zeros(std::vector<Number>* coefficients) {
    while (coefficients->size() > 1 && is_number_zero(coefficients->back())) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(Number(0));
    }
}

std::string format_coefficient(const Number& value) {
    if (value.is_integer()) {
        return value.as_integer().to_string();
    }
    if (value.is_rational()) {
        const auto& r = value.as_rational();
        if (r.denominator() == BigInt(1)) {
            return r.numerator().to_string();
        }
        return r.numerator().to_string() + "/" + r.denominator().to_string();
    }
    return value.to_string();
}

Number polynomial_evaluate_impl(const std::vector<Number>& coefficients, const Number& x) {
    Number result(0);
    for (std::size_t i = coefficients.size(); i > 0; --i) {
        result = result * x + coefficients[i - 1];
    }
    return result;
}

std::vector<Number> polynomial_derivative_impl(const std::vector<Number>& coefficients) {
    if (coefficients.size() <= 1) {
        return {Number(0)};
    }

    std::vector<Number> derivative(coefficients.size() - 1);
    for (std::size_t i = 1; i < coefficients.size(); ++i) {
        derivative[i - 1] = coefficients[i] * Number(static_cast<long long>(i));
    }
    trim_trailing_zeros(&derivative);
    return derivative;
}

Number polynomial_root_bound(const std::vector<Number>& coefficients) {
    const Number leading = coefficients.back();
    Number bound(0);
    for (std::size_t i = 0; i + 1 < coefficients.size(); ++i) {
        const Number ratio = number_abs(coefficients[i] / leading);
        if (ratio.compare(bound) > 0) {
            bound = ratio;
        }
    }
    return Number(1) + bound;
}

Number bisect_root(const std::vector<Number>& coefficients, const Number& left, const Number& right) {
    Number left_value = polynomial_evaluate_impl(coefficients, left);
    const Number eps(BigDecimal::from_string("1e-40"));
    const int max_iterations = 200;

    for (int i = 0; i < max_iterations; ++i) {
        const Number mid = (left + right) / Number(2);
        const Number mid_value = polynomial_evaluate_impl(coefficients, mid);
        if (is_number_zero(mid_value) || (right - mid).compare(eps) <= 0) {
            return mid;
        }

        const bool left_negative = left_value.compare(Number(0)) < 0;
        const bool mid_negative = mid_value.compare(Number(0)) < 0;
        if (left_negative != mid_negative) {
            return bisect_root(coefficients, left, mid);
        }
        return bisect_root(coefficients, mid, right);
    }

    return (left + right) / Number(2);
}

void add_unique_root(std::vector<Number>* roots, const Number& candidate) {
    const Number eps(BigDecimal::from_string("1e-30"));
    for (const Number& existing : *roots) {
        if (number_abs(existing - candidate).compare(eps) <= 0) {
            return;
        }
    }
    roots->push_back(candidate);
}

}  // namespace

Number polynomial_evaluate(const std::vector<Number>& coefficients, const Number& x) {
    return polynomial_evaluate_impl(coefficients, x);
}

std::vector<Number> polynomial_derivative(const std::vector<Number>& coefficients) {
    return polynomial_derivative_impl(coefficients);
}

std::vector<Number> polynomial_add(const std::vector<Number>& lhs,
                                   const std::vector<Number>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<Number> result(size, Number(0));
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] = result[i] + lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] = result[i] + rhs[i];
    }
    trim_trailing_zeros(&result);
    return result;
}

std::vector<Number> polynomial_subtract(const std::vector<Number>& lhs,
                                        const std::vector<Number>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<Number> result(size, Number(0));
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] = result[i] + lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] = result[i] - rhs[i];
    }
    trim_trailing_zeros(&result);
    return result;
}

std::vector<Number> polynomial_multiply(const std::vector<Number>& lhs,
                                        const std::vector<Number>& rhs) {
    std::vector<Number> result(lhs.size() + rhs.size() - 1, Number(0));
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result[i + j] = result[i + j] + lhs[i] * rhs[j];
        }
    }
    trim_trailing_zeros(&result);
    return result;
}

PolynomialDivisionResult polynomial_divide(const std::vector<Number>& dividend,
                                           const std::vector<Number>& divisor) {
    std::vector<Number> normalized_dividend = dividend;
    std::vector<Number> normalized_divisor = divisor;
    trim_trailing_zeros(&normalized_dividend);
    trim_trailing_zeros(&normalized_divisor);

    if (normalized_divisor.size() == 1 && is_number_zero(normalized_divisor[0])) {
        throw std::runtime_error("polynomial divisor cannot be zero");
    }

    if (normalized_dividend.size() < normalized_divisor.size()) {
        return {{Number(0)}, normalized_dividend};
    }

    std::vector<Number> quotient(
        normalized_dividend.size() - normalized_divisor.size() + 1, Number(0));
    std::vector<Number> remainder = normalized_dividend;

    while (remainder.size() >= normalized_divisor.size() &&
           !(remainder.size() == 1 && is_number_zero(remainder[0]))) {
        const std::size_t degree_diff =
            remainder.size() - normalized_divisor.size();
        const Number factor =
            remainder.back() / normalized_divisor.back();
        quotient[degree_diff] = factor;

        for (std::size_t i = 0; i < normalized_divisor.size(); ++i) {
            remainder[degree_diff + i] = remainder[degree_diff + i] - factor * normalized_divisor[i];
        }
        trim_trailing_zeros(&remainder);
        if (remainder.size() < normalized_divisor.size()) {
            break;
        }
    }

    trim_trailing_zeros(&quotient);
    trim_trailing_zeros(&remainder);
    return {quotient, remainder};
}

std::vector<Number> polynomial_real_roots(const std::vector<Number>& coefficients) {
    std::vector<Number> normalized = coefficients;
    trim_trailing_zeros(&normalized);

    if (normalized.size() <= 1) {
        throw std::runtime_error("constant polynomial does not have isolated roots");
    }

    if (normalized.size() == 2) {
        return {-normalized[0] / normalized[1]};
    }

    const std::vector<Number> derivative = polynomial_derivative(normalized);
    std::vector<Number> critical_points = polynomial_real_roots(derivative);
    std::sort(critical_points.begin(), critical_points.end(),
              [](const Number& a, const Number& b) { return a.compare(b) < 0; });

    const Number bound = polynomial_root_bound(normalized);
    std::vector<Number> points;
    points.push_back(-bound);
    for (const Number& point : critical_points) {
        points.push_back(point);
    }
    points.push_back(bound);

    std::vector<Number> roots;
    const Number eps(BigDecimal::from_string("1e-25"));
    for (const Number& point : critical_points) {
        if (is_number_zero(polynomial_evaluate(normalized, point))) {
            add_unique_root(&roots, point);
        }
    }

    for (std::size_t i = 1; i < points.size(); ++i) {
        const Number left = points[i - 1];
        const Number right = points[i];
        const Number left_value = polynomial_evaluate(normalized, left);
        const Number right_value = polynomial_evaluate(normalized, right);

        if (is_number_zero(left_value)) {
            add_unique_root(&roots, left);
            continue;
        }
        if (is_number_zero(right_value)) {
            add_unique_root(&roots, right);
            continue;
        }

        const bool left_negative = left_value.compare(Number(0)) < 0;
        const bool right_negative = right_value.compare(Number(0)) < 0;
        if (left_negative != right_negative) {
            add_unique_root(&roots, bisect_root(normalized, left, right));
        }
    }

    std::sort(roots.begin(), roots.end(),
              [](const Number& a, const Number& b) { return a.compare(b) < 0; });
    return roots;
}

std::vector<Number> polynomial_integral(const std::vector<Number>& coefficients) {
    std::vector<Number> integral(coefficients.size() + 1, Number(0));
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        integral[i + 1] = coefficients[i] / Number(static_cast<long long>(i + 1));
    }
    trim_trailing_zeros(&integral);
    return integral;
}

std::vector<Number> polynomial_compose(const std::vector<Number>& outer,
                                       const std::vector<Number>& inner) {
    std::vector<Number> result = {Number(0)};
    for (std::size_t i = outer.size(); i > 0; --i) {
        result = polynomial_multiply(result, inner);
        if (result.empty()) {
            result.push_back(Number(0));
        }
        result[0] = result[0] + outer[i - 1];
        trim_trailing_zeros(&result);
    }
    return result;
}

std::vector<Number> polynomial_gcd(const std::vector<Number>& lhs,
                                   const std::vector<Number>& rhs) {
    std::vector<Number> a = lhs;
    std::vector<Number> b = rhs;
    trim_trailing_zeros(&a);
    trim_trailing_zeros(&b);

    while (!(b.size() == 1 && is_number_zero(b[0]))) {
        const PolynomialDivisionResult division = polynomial_divide(a, b);
        a = b;
        b = division.remainder;
        trim_trailing_zeros(&a);
        trim_trailing_zeros(&b);
    }

    if (a.empty()) {
        return {Number(0)};
    }
    const Number leading = a.back();
    if (!is_number_zero(leading)) {
        for (Number& coefficient : a) {
            coefficient = coefficient / leading;
        }
    }
    trim_trailing_zeros(&a);
    return a;
}

std::vector<Number> polynomial_fit(const std::vector<Number>& x_samples,
                                   const std::vector<Number>& y_samples,
                                   int degree) {
    if (degree < 0) {
        throw std::runtime_error("polynomial degree must be non-negative");
    }
    if (x_samples.size() != y_samples.size() || x_samples.empty()) {
        throw std::runtime_error("polynomial_fit requires non-empty sample vectors of the same length");
    }
    if (x_samples.size() < static_cast<std::size_t>(degree + 1)) {
        throw std::runtime_error("polynomial_fit requires at least degree + 1 samples");
    }

    // For polynomial fitting, we use a simplified approach with Number
    // This is a basic least squares implementation
    const std::size_t n = x_samples.size();
    const std::size_t m = static_cast<std::size_t>(degree + 1);

    // Build Vandermonde matrix and solve using normal equations
    // A^T * A * c = A^T * y
    std::vector<std::vector<Number>> ata(m, std::vector<Number>(m, Number(0)));
    std::vector<Number> aty(m, Number(0));

    for (std::size_t row = 0; row < n; ++row) {
        std::vector<Number> powers(m, Number(1));
        for (std::size_t j = 1; j < m; ++j) {
            powers[j] = powers[j - 1] * x_samples[row];
        }

        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = 0; j < m; ++j) {
                ata[i][j] = ata[i][j] + powers[i] * powers[j];
            }
            aty[i] = aty[i] + powers[i] * y_samples[row];
        }
    }

    // Gaussian elimination with partial pivoting
    std::vector<Number> coefficients(m, Number(0));
    std::vector<std::vector<Number>> aug = ata;
    for (std::size_t i = 0; i < m; ++i) {
        aug[i].push_back(aty[i]);
    }

    for (std::size_t col = 0; col < m; ++col) {
        // Find pivot
        std::size_t max_row = col;
        for (std::size_t row = col + 1; row < m; ++row) {
            if (number_abs(aug[row][col]).compare(number_abs(aug[max_row][col])) > 0) {
                max_row = row;
            }
        }
        std::swap(aug[col], aug[max_row]);

        if (is_number_zero(aug[col][col])) {
            throw std::runtime_error("polynomial_fit design matrix is rank deficient");
        }

        // Eliminate
        for (std::size_t row = col + 1; row < m; ++row) {
            const Number factor = aug[row][col] / aug[col][col];
            for (std::size_t j = col; j <= m; ++j) {
                aug[row][j] = aug[row][j] - factor * aug[col][j];
            }
        }
    }

    // Back substitution
    for (std::size_t i = m; i-- > 0;) {
        coefficients[i] = aug[i][m];
        for (std::size_t j = i + 1; j < m; ++j) {
            coefficients[i] = coefficients[i] - aug[i][j] * coefficients[j];
        }
        coefficients[i] = coefficients[i] / aug[i][i];
    }

    trim_trailing_zeros(&coefficients);
    return coefficients;
}

std::string polynomial_to_string(const std::vector<Number>& coefficients,
                                 const std::string& variable_name) {
    std::vector<Number> normalized = coefficients;
    trim_trailing_zeros(&normalized);
    if (normalized.size() == 1 && is_number_zero(normalized[0])) {
        return "0";
    }

    std::string result;
    bool first = true;
    for (std::size_t index = normalized.size(); index > 0; --index) {
        const std::size_t degree = index - 1;
        const Number coefficient = normalized[degree];
        if (is_number_zero(coefficient)) {
            continue;
        }

        const bool negative = coefficient.compare(Number(0)) < 0;
        const Number abs_value = number_abs(coefficient);

        std::string term;
        if (degree == 0) {
            term = format_coefficient(abs_value);
        } else {
            if (!is_number_one(abs_value)) {
                term += format_coefficient(abs_value) + " * ";
            }
            term += variable_name;
            if (degree > 1) {
                term += " ^ " + std::to_string(degree);
            }
        }

        if (first) {
            result += negative ? "-" + term : term;
            first = false;
        } else {
            result += negative ? " - " + term : " + " + term;
        }
    }

    return result.empty() ? "0" : result;
}

// ============================================================================
// 兼容性接口：double 版本
// ============================================================================

namespace {

std::vector<Number> doubles_to_numbers(const std::vector<double>& values) {
    std::vector<Number> result;
    result.reserve(values.size());
    for (double v : values) {
        result.push_back(Number(BigDecimal::from_string(std::to_string(v))));
    }
    return result;
}

std::vector<double> numbers_to_doubles(const std::vector<Number>& values) {
    std::vector<double> result;
    result.reserve(values.size());
    for (const Number& v : values) {
        result.push_back(static_cast<double>(std::stod(v.to_string())));
    }
    return result;
}

}  // namespace

std::vector<double> polynomial_add_double(const std::vector<double>& lhs,
                                          const std::vector<double>& rhs) {
    return numbers_to_doubles(polynomial_add(doubles_to_numbers(lhs), doubles_to_numbers(rhs)));
}

std::vector<double> polynomial_subtract_double(const std::vector<double>& lhs,
                                               const std::vector<double>& rhs) {
    return numbers_to_doubles(polynomial_subtract(doubles_to_numbers(lhs), doubles_to_numbers(rhs)));
}

std::vector<double> polynomial_multiply_double(const std::vector<double>& lhs,
                                               const std::vector<double>& rhs) {
    return numbers_to_doubles(polynomial_multiply(doubles_to_numbers(lhs), doubles_to_numbers(rhs)));
}

double polynomial_evaluate_double(const std::vector<double>& coefficients, double x) {
    Number result = polynomial_evaluate(doubles_to_numbers(coefficients),
                                        Number(BigDecimal::from_string(std::to_string(x))));
    return static_cast<double>(std::stod(result.to_string()));
}

std::vector<double> polynomial_derivative_double(const std::vector<double>& coefficients) {
    return numbers_to_doubles(polynomial_derivative(doubles_to_numbers(coefficients)));
}

std::vector<double> polynomial_real_roots_double(const std::vector<double>& coefficients) {
    return numbers_to_doubles(polynomial_real_roots(doubles_to_numbers(coefficients)));
}

PolynomialDivisionResultDouble polynomial_divide_double(const std::vector<double>& dividend,
                                                        const std::vector<double>& divisor) {
    PolynomialDivisionResult result = polynomial_divide(doubles_to_numbers(dividend),
                                                        doubles_to_numbers(divisor));
    PolynomialDivisionResultDouble double_result;
    double_result.quotient = numbers_to_doubles(result.quotient);
    double_result.remainder = numbers_to_doubles(result.remainder);
    return double_result;
}

std::string polynomial_to_string_double(const std::vector<double>& coefficients,
                                        const std::string& variable_name) {
    return polynomial_to_string(doubles_to_numbers(coefficients), variable_name);
}

std::vector<double> polynomial_integral_double(const std::vector<double>& coefficients) {
    return numbers_to_doubles(polynomial_integral(doubles_to_numbers(coefficients)));
}

std::vector<double> polynomial_fit_double(const std::vector<double>& x_samples,
                                          const std::vector<double>& y_samples,
                                          int degree) {
    return numbers_to_doubles(polynomial_fit(doubles_to_numbers(x_samples),
                                             doubles_to_numbers(y_samples),
                                             degree));
}

std::vector<double> polynomial_compose_double(const std::vector<double>& outer,
                                              const std::vector<double>& inner) {
    return numbers_to_doubles(polynomial_compose(doubles_to_numbers(outer),
                                                 doubles_to_numbers(inner)));
}

std::vector<double> polynomial_gcd_double(const std::vector<double>& lhs,
                                          const std::vector<double>& rhs) {
    return numbers_to_doubles(polynomial_gcd(doubles_to_numbers(lhs),
                                             doubles_to_numbers(rhs)));
}

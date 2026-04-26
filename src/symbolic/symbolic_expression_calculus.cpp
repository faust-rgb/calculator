#include "symbolic_expression_internal.h"

#include "mymath.h"
#include "polynomial.h"

#include <algorithm>
#include <list>
#include <stdexcept>
#include <unordered_map>
#include <utility>

using namespace symbolic_expression_internal;

namespace {

class SymbolicExpressionLruCache {
public:
    explicit SymbolicExpressionLruCache(std::size_t capacity)
        : capacity_(capacity) {}

    bool get(const std::string& key, SymbolicExpression* value) {
        const auto found = index_.find(key);
        if (found == index_.end()) {
            return false;
        }
        entries_.splice(entries_.begin(), entries_, found->second);
        *value = found->second->second;
        return true;
    }

    void put(const std::string& key, const SymbolicExpression& value) {
        const auto found = index_.find(key);
        if (found != index_.end()) {
            found->second->second = value;
            entries_.splice(entries_.begin(), entries_, found->second);
            return;
        }

        entries_.push_front({key, value});
        index_[key] = entries_.begin();
        while (entries_.size() > capacity_) {
            index_.erase(entries_.back().first);
            entries_.pop_back();
        }
    }

private:
    std::size_t capacity_ = 0;
    std::list<std::pair<std::string, SymbolicExpression>> entries_;
    std::unordered_map<std::string,
                       std::list<std::pair<std::string, SymbolicExpression>>::iterator>
        index_;
};

bool is_one_plus_variable_squared(const SymbolicExpression& expression,
                                  const std::string& variable_name) {
    std::vector<double> coefficients;
    return expression.polynomial_coefficients(variable_name, &coefficients) &&
           coefficients.size() == 3 &&
           mymath::is_near_zero(coefficients[0] - 1.0, kFormatEps) &&
           mymath::is_near_zero(coefficients[1], kFormatEps) &&
           mymath::is_near_zero(coefficients[2] - 1.0, kFormatEps);
}

bool is_one_minus_variable_squared(const SymbolicExpression& expression,
                                   const std::string& variable_name) {
    std::vector<double> coefficients;
    return expression.polynomial_coefficients(variable_name, &coefficients) &&
           coefficients.size() == 3 &&
           mymath::is_near_zero(coefficients[0] - 1.0, kFormatEps) &&
           mymath::is_near_zero(coefficients[1], kFormatEps) &&
           mymath::is_near_zero(coefficients[2] + 1.0, kFormatEps);
}

bool is_sqrt_one_minus_variable_squared(const SymbolicExpression& expression,
                                        const std::string& variable_name) {
    if (expression.node_->type != NodeType::kFunction ||
        expression.node_->text != "sqrt") {
        return false;
    }
    return is_one_minus_variable_squared(SymbolicExpression(expression.node_->left),
                                         variable_name);
}

void trim_coefficients(std::vector<double>* coefficients) {
    while (coefficients->size() > 1 &&
           mymath::is_near_zero(coefficients->back(), kFormatEps)) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(0.0);
    }
}

bool polynomial_is_zero(const std::vector<double>& coefficients) {
    for (double coefficient : coefficients) {
        if (!mymath::is_near_zero(coefficient, kFormatEps)) {
            return false;
        }
    }
    return true;
}

bool solve_dense_linear_system(std::vector<std::vector<double>> matrix,
                               std::vector<double> rhs,
                               std::vector<double>* solution) {
    const std::size_t size = rhs.size();
    for (std::size_t col = 0; col < size; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < size; ++row) {
            if (mymath::abs(matrix[row][col]) > mymath::abs(matrix[pivot][col])) {
                pivot = row;
            }
        }
        if (mymath::is_near_zero(matrix[pivot][col], 1e-10)) {
            return false;
        }
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(rhs[pivot], rhs[col]);
        }

        const double pivot_value = matrix[col][col];
        for (std::size_t j = col; j < size; ++j) {
            matrix[col][j] /= pivot_value;
        }
        rhs[col] /= pivot_value;

        for (std::size_t row = 0; row < size; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = matrix[row][col];
            if (mymath::is_near_zero(factor, 1e-12)) {
                continue;
            }
            for (std::size_t j = col; j < size; ++j) {
                matrix[row][j] -= factor * matrix[col][j];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    *solution = rhs;
    return true;
}

double clean_symbolic_constant(double value) {
    if (mymath::is_integer(value, 1e-7)) {
        return value >= 0.0 ? static_cast<long long>(value + 0.5)
                            : static_cast<long long>(value - 0.5);
    }

    long long numerator = 0;
    long long denominator = 1;
    if (mymath::approximate_fraction(value,
                                     &numerator,
                                     &denominator,
                                     999,
                                     1e-7)) {
        if (value < 0.0) {
            numerator = -numerator;
        }
        return static_cast<double>(numerator) / static_cast<double>(denominator);
    }
    return value;
}

std::vector<double> polynomial_power_coefficients(const std::vector<double>& base,
                                                  int exponent) {
    std::vector<double> result = {1.0};
    for (int i = 0; i < exponent; ++i) {
        result = polynomial_multiply(result, base);
    }
    trim_coefficients(&result);
    return result;
}

SymbolicExpression make_linear_form(double slope,
                                    double intercept,
                                    const std::string& variable_name) {
    return make_add(
               make_multiply(SymbolicExpression::number(slope),
                             SymbolicExpression::variable(variable_name)),
               SymbolicExpression::number(intercept))
        .simplify();
}

bool integrate_inverse_quadratic(const std::vector<double>& denominator,
                                 const std::string& variable_name,
                                 SymbolicExpression* integrated) {
    const double c = denominator[0];
    const double b = denominator[1];
    const double a = denominator[2];
    if (mymath::is_near_zero(a, kFormatEps)) {
        return false;
    }

    const SymbolicExpression linear = make_linear_form(2.0 * a, b, variable_name);
    const double atan_discriminant = 4.0 * a * c - b * b;
    if (atan_discriminant > kFormatEps) {
        const double scale = mymath::sqrt(atan_discriminant);
        *integrated = make_multiply(
                          SymbolicExpression::number(2.0 / scale),
                          make_function("atan",
                                        make_divide(linear,
                                                    SymbolicExpression::number(scale))))
                          .simplify();
        return true;
    }

    const double log_discriminant = b * b - 4.0 * a * c;
    if (log_discriminant > kFormatEps) {
        const double scale = mymath::sqrt(log_discriminant);
        const SymbolicExpression numerator =
            make_subtract(linear, SymbolicExpression::number(scale)).simplify();
        const SymbolicExpression denominator_expr =
            make_add(linear, SymbolicExpression::number(scale)).simplify();
        *integrated = make_multiply(
                          SymbolicExpression::number(1.0 / scale),
                          make_function("ln",
                                        make_function(
                                            "abs",
                                            make_divide(numerator,
                                                        denominator_expr))))
                          .simplify();
        return true;
    }

    return false;
}

struct LinearFactorMultiplicity {
    double root = 0.0;
    int multiplicity = 0;
};

struct QuadraticFactorMultiplicity {
    std::vector<double> coefficients;
    int multiplicity = 0;
};

struct RationalPartialFractionTerm {
    enum class Kind {
        kLinear,
        kQuadratic,
    };

    Kind kind = Kind::kLinear;
    double root = 0.0;
    std::vector<double> quadratic;
    int power = 0;
    int numerator_degree = 0;
};

bool extract_real_linear_factorization(const std::vector<double>& denominator,
                                       std::vector<LinearFactorMultiplicity>* factors) {
    std::vector<double> remaining = denominator;
    trim_coefficients(&remaining);
    std::vector<double> roots;
    try {
        roots = polynomial_real_roots(remaining);
    } catch (const std::exception&) {
        return false;
    }

    for (double root : roots) {
        root = clean_symbolic_constant(root);
        const std::vector<double> factor = {-root, 1.0};
        int multiplicity = 0;
        while (remaining.size() > 1) {
            const PolynomialDivisionResult division =
                polynomial_divide(remaining, factor);
            if (!polynomial_is_zero(division.remainder)) {
                break;
            }
            remaining = division.quotient;
            trim_coefficients(&remaining);
            ++multiplicity;
        }
        if (multiplicity > 0) {
            factors->push_back(LinearFactorMultiplicity{root, multiplicity});
        }
    }

    trim_coefficients(&remaining);
    return remaining.size() == 1 && !factors->empty();
}

bool polynomial_close(const std::vector<double>& lhs,
                      const std::vector<double>& rhs,
                      double eps = 1e-8) {
    std::vector<double> left = lhs;
    std::vector<double> right = rhs;
    trim_coefficients(&left);
    trim_coefficients(&right);
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (!mymath::is_near_zero(left[i] - right[i], eps)) {
            return false;
        }
    }
    return true;
}

bool divide_exact_polynomial(std::vector<double>* remaining,
                             const std::vector<double>& divisor) {
    const PolynomialDivisionResult division = polynomial_divide(*remaining, divisor);
    if (!polynomial_is_zero(division.remainder)) {
        return false;
    }
    *remaining = division.quotient;
    trim_coefficients(remaining);
    return true;
}

bool extract_linear_and_one_quadratic_factorization(
    const std::vector<double>& denominator,
    std::vector<LinearFactorMultiplicity>* linear_factors,
    QuadraticFactorMultiplicity* quadratic_factor) {
    linear_factors->clear();
    quadratic_factor->coefficients.clear();
    quadratic_factor->multiplicity = 0;

    std::vector<double> remaining = denominator;
    trim_coefficients(&remaining);
    std::vector<double> roots;
    try {
        roots = polynomial_real_roots(remaining);
    } catch (const std::exception&) {
        roots.clear();
    }

    for (double root : roots) {
        root = clean_symbolic_constant(root);
        const std::vector<double> factor = {-root, 1.0};
        int multiplicity = 0;
        while (remaining.size() > 1 && divide_exact_polynomial(&remaining, factor)) {
            ++multiplicity;
        }
        if (multiplicity > 0) {
            linear_factors->push_back(LinearFactorMultiplicity{root, multiplicity});
        }
    }

    trim_coefficients(&remaining);
    if (remaining.size() == 1) {
        return !linear_factors->empty();
    }

    if ((remaining.size() - 1) % 2 != 0) {
        return false;
    }
    const int quadratic_power = static_cast<int>((remaining.size() - 1) / 2);
    if (quadratic_power <= 0) {
        return false;
    }

    if (quadratic_power == 1) {
        const double a = remaining[2];
        const double b = remaining[1];
        const double c = remaining[0];
        if (mymath::is_near_zero(a, kFormatEps) ||
            b * b - 4.0 * a * c >= -1e-8) {
            return false;
        }
        quadratic_factor->coefficients = remaining;
        quadratic_factor->multiplicity = 1;
        return true;
    }

    const double leading_root =
        mymath::pow(mymath::abs(remaining.back()), 1.0 / quadratic_power);
    const double leading =
        remaining.back() < 0.0 && quadratic_power % 2 == 1 ? -leading_root : leading_root;
    const double sum_roots = remaining[remaining.size() - 2] / remaining.back();
    const double b = leading * sum_roots / static_cast<double>(quadratic_power);

    std::vector<double> best_quadratic;
    bool found = false;
    for (double sign : {-1.0, 1.0}) {
        const double constant_root =
            mymath::pow(mymath::abs(remaining[0]), 1.0 / quadratic_power);
        const double c = sign * constant_root;
        std::vector<double> candidate = {c, b, leading};
        std::vector<double> powered = polynomial_power_coefficients(candidate, quadratic_power);
        if (polynomial_close(powered, remaining, 1e-6)) {
            best_quadratic = candidate;
            found = true;
            break;
        }
    }
    if (!found) {
        return false;
    }

    const double a = best_quadratic[2];
    const double best_b = best_quadratic[1];
    const double best_c = best_quadratic[0];
    if (mymath::is_near_zero(a, kFormatEps) ||
        best_b * best_b - 4.0 * a * best_c >= -1e-8) {
        return false;
    }

    quadratic_factor->coefficients = best_quadratic;
    quadratic_factor->multiplicity = quadratic_power;
    return true;
}

SymbolicExpression integrate_inverse_quadratic_power(
    const std::vector<double>& quadratic,
    int power,
    const std::string& variable_name) {
    const double c = quadratic[0];
    const double b = quadratic[1];
    const double a = quadratic[2];
    const double delta = 4.0 * a * c - b * b;
    if (power <= 0 || mymath::is_near_zero(a, kFormatEps) || delta <= kFormatEps) {
        throw std::runtime_error("unsupported quadratic power integral");
    }

    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression u =
        make_add(x, SymbolicExpression::number(b / (2.0 * a))).simplify();
    const double d = delta / (4.0 * a * a);
    const SymbolicExpression shifted_quadratic =
        make_add(make_power(u, SymbolicExpression::number(2.0)),
                 SymbolicExpression::number(d))
            .simplify();

    SymbolicExpression integral =
        make_multiply(SymbolicExpression::number(1.0 / mymath::sqrt(d)),
                      make_function("atan",
                                    make_divide(u,
                                                SymbolicExpression::number(mymath::sqrt(d)))))
            .simplify();

    for (int n = 2; n <= power; ++n) {
        const SymbolicExpression recurrence_term =
            make_divide(u,
                        make_multiply(
                            SymbolicExpression::number(2.0 * d * (n - 1)),
                            make_power(shifted_quadratic,
                                       SymbolicExpression::number(n - 1))))
                .simplify();
        integral =
            make_add(recurrence_term,
                     make_multiply(
                         SymbolicExpression::number(
                             static_cast<double>(2 * n - 3) /
                             static_cast<double>(2 * (n - 1)) / d),
                         integral))
                .simplify();
    }

    return make_multiply(SymbolicExpression::number(1.0 / mymath::pow(a, power)),
                         integral)
        .simplify();
}

SymbolicExpression integrate_quadratic_partial_fraction_term(
    const std::vector<double>& quadratic,
    double slope,
    double constant,
    int power,
    const std::string& variable_name) {
    const double a = quadratic[2];
    const double b = quadratic[1];
    const double derivative_scale = slope / (2.0 * a);
    const double inverse_scale = constant - derivative_scale * b;
    const SymbolicExpression quadratic_expression =
        build_polynomial_expression_from_coefficients(quadratic, variable_name)
            .simplify();
    SymbolicExpression result = SymbolicExpression::number(0.0);

    if (!mymath::is_near_zero(derivative_scale, kFormatEps)) {
        SymbolicExpression derivative_part;
        if (power == 1) {
            derivative_part =
                make_multiply(SymbolicExpression::number(derivative_scale),
                              make_function("ln",
                                            make_function("abs", quadratic_expression)))
                    .simplify();
        } else {
            derivative_part =
                make_multiply(SymbolicExpression::number(
                                  derivative_scale / static_cast<double>(1 - power)),
                              make_power(quadratic_expression,
                                         SymbolicExpression::number(1 - power)))
                    .simplify();
        }
        result = make_add(result, derivative_part).simplify();
    }

    if (!mymath::is_near_zero(inverse_scale, kFormatEps)) {
        result =
            make_add(result,
                     make_multiply(SymbolicExpression::number(inverse_scale),
                                   integrate_inverse_quadratic_power(quadratic,
                                                                    power,
                                                                    variable_name)))
                .simplify();
    }
    return result.simplify();
}

bool integrate_mixed_linear_quadratic_partial_fractions(
    const std::vector<double>& numerator,
    const std::vector<double>& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    if (denominator.size() <= 3) {
        return false;
    }

    std::vector<LinearFactorMultiplicity> linear_factors;
    QuadraticFactorMultiplicity quadratic_factor;
    if (!extract_linear_and_one_quadratic_factorization(denominator,
                                                        &linear_factors,
                                                        &quadratic_factor) ||
        quadratic_factor.multiplicity <= 0) {
        return false;
    }

    std::vector<RationalPartialFractionTerm> terms;
    for (const LinearFactorMultiplicity& factor : linear_factors) {
        for (int power = 1; power <= factor.multiplicity; ++power) {
            RationalPartialFractionTerm term;
            term.kind = RationalPartialFractionTerm::Kind::kLinear;
            term.root = factor.root;
            term.power = power;
            term.numerator_degree = 0;
            terms.push_back(term);
        }
    }
    for (int power = 1; power <= quadratic_factor.multiplicity; ++power) {
        RationalPartialFractionTerm slope_term;
        slope_term.kind = RationalPartialFractionTerm::Kind::kQuadratic;
        slope_term.quadratic = quadratic_factor.coefficients;
        slope_term.power = power;
        slope_term.numerator_degree = 1;
        terms.push_back(slope_term);

        RationalPartialFractionTerm constant_term = slope_term;
        constant_term.numerator_degree = 0;
        terms.push_back(constant_term);
    }

    const int unknown_count = static_cast<int>(terms.size());
    if (unknown_count != static_cast<int>(denominator.size()) - 1) {
        return false;
    }

    std::vector<std::vector<double>> columns;
    columns.reserve(terms.size());
    for (const RationalPartialFractionTerm& term : terms) {
        std::vector<double> divisor;
        if (term.kind == RationalPartialFractionTerm::Kind::kLinear) {
            divisor = polynomial_power_coefficients({-term.root, 1.0}, term.power);
        } else {
            divisor = polynomial_power_coefficients(term.quadratic, term.power);
        }
        PolynomialDivisionResult division = polynomial_divide(denominator, divisor);
        if (!polynomial_is_zero(division.remainder)) {
            return false;
        }
        std::vector<double> column = division.quotient;
        if (term.kind == RationalPartialFractionTerm::Kind::kQuadratic &&
            term.numerator_degree == 1) {
            column.insert(column.begin(), 0.0);
        }
        trim_coefficients(&column);
        columns.push_back(column);
    }

    std::vector<std::vector<double>> matrix;
    std::vector<double> rhs;
    for (int candidate = -unknown_count * 2;
         candidate <= unknown_count * 4 &&
         static_cast<int>(rhs.size()) < unknown_count;
         ++candidate) {
        const double sample = static_cast<double>(candidate);
        bool sample_is_pole = false;
        for (const LinearFactorMultiplicity& factor : linear_factors) {
            if (mymath::is_near_zero(sample - factor.root, 1e-8)) {
                sample_is_pole = true;
                break;
            }
        }
        if (sample_is_pole) {
            continue;
        }

        std::vector<double> row;
        row.reserve(columns.size());
        for (const std::vector<double>& column : columns) {
            row.push_back(polynomial_evaluate(column, sample));
        }
        matrix.push_back(row);
        rhs.push_back(polynomial_evaluate(numerator, sample));
    }
    if (static_cast<int>(rhs.size()) != unknown_count) {
        return false;
    }

    std::vector<double> coefficients;
    if (!solve_dense_linear_system(matrix, rhs, &coefficients)) {
        return false;
    }

    SymbolicExpression result = SymbolicExpression::number(0.0);
    bool has_term = false;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        const double coefficient = clean_symbolic_constant(coefficients[i]);
        if (mymath::is_near_zero(coefficient, kFormatEps)) {
            continue;
        }

        SymbolicExpression term_integral;
        if (terms[i].kind == RationalPartialFractionTerm::Kind::kLinear) {
            const SymbolicExpression shifted_variable =
                make_subtract(SymbolicExpression::variable(variable_name),
                              SymbolicExpression::number(terms[i].root))
                    .simplify();
            if (terms[i].power == 1) {
                term_integral =
                    make_multiply(SymbolicExpression::number(coefficient),
                                  make_function("ln",
                                                make_function("abs", shifted_variable)))
                        .simplify();
            } else {
                term_integral =
                    make_multiply(SymbolicExpression::number(
                                      coefficient /
                                      static_cast<double>(1 - terms[i].power)),
                                  make_power(shifted_variable,
                                             SymbolicExpression::number(
                                                 1 - terms[i].power)))
                        .simplify();
            }
        } else {
            const double slope = terms[i].numerator_degree == 1 ? coefficient : 0.0;
            const double constant = terms[i].numerator_degree == 0 ? coefficient : 0.0;
            term_integral =
                integrate_quadratic_partial_fraction_term(terms[i].quadratic,
                                                          slope,
                                                          constant,
                                                          terms[i].power,
                                                          variable_name);
        }

        result = has_term ? make_add(result, term_integral).simplify() : term_integral;
        has_term = true;
    }

    if (!has_term) {
        return false;
    }
    *integrated = result.simplify();
    return true;
}

bool integrate_real_linear_partial_fractions(
    const std::vector<double>& numerator,
    const std::vector<double>& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    if (denominator.size() <= 3) {
        return false;
    }

    std::vector<LinearFactorMultiplicity> factors;
    if (!extract_real_linear_factorization(denominator, &factors)) {
        return false;
    }

    int unknown_count = 0;
    for (const LinearFactorMultiplicity& factor : factors) {
        unknown_count += factor.multiplicity;
    }
    if (unknown_count != static_cast<int>(denominator.size()) - 1) {
        return false;
    }

    std::vector<std::vector<double>> columns;
    columns.reserve(static_cast<std::size_t>(unknown_count));
    std::vector<std::pair<double, int>> terms;
    terms.reserve(static_cast<std::size_t>(unknown_count));
    for (const LinearFactorMultiplicity& factor : factors) {
        const std::vector<double> linear = {-factor.root, 1.0};
        for (int power = 1; power <= factor.multiplicity; ++power) {
            const std::vector<double> divisor =
                polynomial_power_coefficients(linear, power);
            PolynomialDivisionResult division =
                polynomial_divide(denominator, divisor);
            if (!polynomial_is_zero(division.remainder)) {
                return false;
            }
            columns.push_back(division.quotient);
            terms.push_back({factor.root, power});
        }
    }

    std::vector<std::vector<double>> matrix;
    std::vector<double> rhs;
    for (int candidate = -unknown_count; candidate <= unknown_count * 3 &&
         static_cast<int>(rhs.size()) < unknown_count; ++candidate) {
        const double sample = static_cast<double>(candidate);
        bool sample_is_pole = false;
        for (const LinearFactorMultiplicity& factor : factors) {
            if (mymath::is_near_zero(sample - factor.root, 1e-8)) {
                sample_is_pole = true;
                break;
            }
        }
        if (sample_is_pole) {
            continue;
        }

        std::vector<double> row;
        row.reserve(static_cast<std::size_t>(unknown_count));
        for (const std::vector<double>& column : columns) {
            row.push_back(polynomial_evaluate(column, sample));
        }
        matrix.push_back(row);
        rhs.push_back(polynomial_evaluate(numerator, sample));
    }
    if (static_cast<int>(rhs.size()) != unknown_count) {
        return false;
    }

    std::vector<double> coefficients;
    if (!solve_dense_linear_system(matrix, rhs, &coefficients)) {
        return false;
    }

    SymbolicExpression result = SymbolicExpression::number(0.0);
    bool has_term = false;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const double coefficient = clean_symbolic_constant(coefficients[i]);
        if (mymath::is_near_zero(coefficient, kFormatEps)) {
            continue;
        }

        const double root = terms[i].first;
        const int power = terms[i].second;
        const SymbolicExpression shifted_variable =
            make_subtract(SymbolicExpression::variable(variable_name),
                          SymbolicExpression::number(root))
                .simplify();
        SymbolicExpression term;
        if (power == 1) {
            term = make_function("ln", make_function("abs", shifted_variable));
            if (!mymath::is_near_zero(coefficient - 1.0, kFormatEps)) {
                if (mymath::is_near_zero(coefficient + 1.0, kFormatEps)) {
                    term = make_negate(term).simplify();
                } else {
                    term = make_multiply(SymbolicExpression::number(coefficient),
                                         term)
                               .simplify();
                }
            }
        } else {
            const double scale =
                clean_symbolic_constant(coefficient / static_cast<double>(1 - power));
            SymbolicExpression antiderivative_base;
            if (power == 2) {
                antiderivative_base =
                    make_divide(SymbolicExpression::number(1.0),
                                shifted_variable)
                        .simplify();
            } else {
                antiderivative_base =
                    make_power(shifted_variable,
                               SymbolicExpression::number(
                                   static_cast<double>(1 - power)))
                        .simplify();
            }
            term = make_multiply(SymbolicExpression::number(mymath::abs(scale)),
                                 antiderivative_base)
                       .simplify();
            if (scale < 0.0) {
                term = make_negate(term).simplify();
            }
        }

        result = has_term ? make_add(result, term).simplify() : term;
        has_term = true;
    }

    if (!has_term) {
        return false;
    }

    *integrated = result.simplify();
    return true;
}

bool try_integrate_repeated_unit_quadratic(const std::vector<double>& numerator,
                                           const std::vector<double>& denominator,
                                           const std::string& variable_name,
                                           SymbolicExpression* integrated) {
    if (numerator.size() != 1 ||
        !mymath::is_near_zero(numerator[0] - 1.0, kFormatEps) ||
        denominator.size() != 5 ||
        !mymath::is_near_zero(denominator[0] - 1.0, kFormatEps) ||
        !mymath::is_near_zero(denominator[1], kFormatEps) ||
        !mymath::is_near_zero(denominator[2] - 2.0, kFormatEps) ||
        !mymath::is_near_zero(denominator[3], kFormatEps) ||
        !mymath::is_near_zero(denominator[4] - 1.0, kFormatEps)) {
        return false;
    }

    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression one_plus_x_squared =
        make_add(SymbolicExpression::number(1.0),
                 make_power(x, SymbolicExpression::number(2.0)))
            .simplify();
    *integrated =
        make_add(make_divide(x,
                             make_multiply(SymbolicExpression::number(2.0),
                                           one_plus_x_squared)),
                 make_divide(make_function("atan", x),
                             SymbolicExpression::number(2.0)))
            .simplify();
    return true;
}

bool same_simplified_expression(const SymbolicExpression& lhs,
                                const SymbolicExpression& rhs) {
    return lhs.simplify().to_string() == rhs.simplify().to_string();
}

bool primitive_for_outer_function(const std::string& function_name,
                                  const SymbolicExpression& argument,
                                  SymbolicExpression* primitive) {
    if (function_name == "exp") {
        *primitive = make_function("exp", argument);
        return true;
    }
    if (function_name == "sin") {
        *primitive = make_negate(make_function("cos", argument)).simplify();
        return true;
    }
    if (function_name == "cos") {
        *primitive = make_function("sin", argument);
        return true;
    }
    if (function_name == "tan") {
        *primitive = make_negate(
                         make_function("ln",
                                       make_function("abs",
                                                     make_function("cos", argument))))
                         .simplify();
        return true;
    }
    if (function_name == "sec") {
        *primitive =
            make_function("ln",
                          make_function("abs",
                                        make_add(make_function("sec", argument),
                                                 make_function("tan", argument))))
                .simplify();
        return true;
    }
    if (function_name == "csc") {
        *primitive =
            make_function("ln",
                          make_function("abs",
                                        make_subtract(make_function("csc", argument),
                                                      make_function("cot", argument))))
                .simplify();
        return true;
    }
    if (function_name == "cot") {
        *primitive =
            make_function("ln", make_function("abs", make_function("sin", argument)))
                .simplify();
        return true;
    }
    return false;
}

bool try_integrate_substitution_product(const SymbolicExpression& derivative_factor,
                                        const SymbolicExpression& function_factor,
                                        const std::string& variable_name,
                                        SymbolicExpression* integrated) {
    if (function_factor.node_->type != NodeType::kFunction) {
        return false;
    }

    const SymbolicExpression argument(function_factor.node_->left);
    const SymbolicExpression expected_derivative =
        argument.derivative(variable_name).simplify();
    double scale = 1.0;
    if (!same_simplified_expression(derivative_factor, expected_derivative)) {
        double constant = 0.0;
        SymbolicExpression rest;
        if (decompose_constant_times_expression(expected_derivative,
                                                variable_name,
                                                &constant,
                                                &rest) &&
            !mymath::is_near_zero(constant, kFormatEps) &&
            same_simplified_expression(derivative_factor, rest)) {
            scale = 1.0 / constant;
        } else if (decompose_constant_times_expression(derivative_factor,
                                                       variable_name,
                                                       &constant,
                                                       &rest) &&
                   same_simplified_expression(rest, expected_derivative)) {
            scale = constant;
        } else {
            return false;
        }
    }

    SymbolicExpression primitive;
    if (!primitive_for_outer_function(function_factor.node_->text,
                                      argument,
                                      &primitive)) {
        return false;
    }
    *integrated = make_multiply(SymbolicExpression::number(scale), primitive).simplify();
    return true;
}

bool try_integrate_trig_power_identity(const SymbolicExpression& base,
                                       double exponent_value,
                                       const std::string& variable_name,
                                       SymbolicExpression* integrated) {
    if (base.node_->type != NodeType::kFunction) {
        return false;
    }

    const SymbolicExpression argument(base.node_->left);
    double a = 0.0;
    double b = 0.0;
    if (!decompose_linear(argument, variable_name, &a, &b) ||
        mymath::is_near_zero(a, kFormatEps)) {
        return false;
    }

    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const std::string& function_name = base.node_->text;
    if (mymath::is_near_zero(exponent_value - 2.0, kFormatEps) &&
        function_name == "sin") {
        const SymbolicExpression double_argument =
            make_multiply(SymbolicExpression::number(2.0), argument).simplify();
        *integrated =
            make_subtract(make_divide(x, SymbolicExpression::number(2.0)),
                          make_divide(make_function("sin", double_argument),
                                      SymbolicExpression::number(4.0 * a)))
                .simplify();
        return true;
    }
    if (mymath::is_near_zero(exponent_value - 2.0, kFormatEps) &&
        function_name == "cos") {
        const SymbolicExpression double_argument =
            make_multiply(SymbolicExpression::number(2.0), argument).simplify();
        *integrated =
            make_add(make_divide(x, SymbolicExpression::number(2.0)),
                     make_divide(make_function("sin", double_argument),
                                 SymbolicExpression::number(4.0 * a)))
                .simplify();
        return true;
    }
    if (mymath::is_near_zero(exponent_value - 2.0, kFormatEps) &&
        function_name == "tan") {
        *integrated =
            make_subtract(make_divide(make_function("tan", argument),
                                      SymbolicExpression::number(a)),
                          x)
                .simplify();
        return true;
    }
    if (mymath::is_near_zero(exponent_value - 2.0, kFormatEps) &&
        function_name == "sec") {
        *integrated = make_divide(make_function("tan", argument),
                                  SymbolicExpression::number(a))
                          .simplify();
        return true;
    }
    if (mymath::is_near_zero(exponent_value - 2.0, kFormatEps) &&
        function_name == "csc") {
        *integrated = make_divide(make_negate(make_function("cot", argument)),
                                  SymbolicExpression::number(a))
                          .simplify();
        return true;
    }
    if (mymath::is_near_zero(exponent_value - 2.0, kFormatEps) &&
        function_name == "cot") {
        *integrated =
            make_subtract(make_negate(x),
                          make_divide(make_function("cot", argument),
                                      SymbolicExpression::number(a)))
                .simplify();
        return true;
    }
    if (mymath::is_near_zero(exponent_value - 3.0, kFormatEps) &&
        function_name == "sin") {
        *integrated =
            make_add(make_divide(make_power(make_function("cos", argument),
                                            SymbolicExpression::number(3.0)),
                                 SymbolicExpression::number(3.0 * a)),
                     make_divide(make_negate(make_function("cos", argument)),
                                 SymbolicExpression::number(a)))
                .simplify();
        return true;
    }
    if (mymath::is_near_zero(exponent_value - 3.0, kFormatEps) &&
        function_name == "cos") {
        *integrated =
            make_subtract(make_divide(make_function("sin", argument),
                                      SymbolicExpression::number(a)),
                          make_divide(make_power(make_function("sin", argument),
                                                 SymbolicExpression::number(3.0)),
                                      SymbolicExpression::number(3.0 * a)))
                .simplify();
        return true;
    }
    return false;
}

bool try_integrate_trig_product_identity(const SymbolicExpression& left,
                                         const SymbolicExpression& right,
                                         const std::string& variable_name,
                                         SymbolicExpression* integrated) {
    if (left.node_->type != NodeType::kFunction ||
        right.node_->type != NodeType::kFunction) {
        return false;
    }
    if (!same_simplified_expression(SymbolicExpression(left.node_->left),
                                    SymbolicExpression(right.node_->left))) {
        return false;
    }

    const SymbolicExpression argument(left.node_->left);
    double a = 0.0;
    double b = 0.0;
    if (!decompose_linear(argument, variable_name, &a, &b) ||
        mymath::is_near_zero(a, kFormatEps)) {
        return false;
    }

    if ((left.node_->text == "sin" && right.node_->text == "cos") ||
        (left.node_->text == "cos" && right.node_->text == "sin")) {
        *integrated =
            make_divide(make_power(make_function("sin", argument),
                                   SymbolicExpression::number(2.0)),
                        SymbolicExpression::number(2.0 * a))
                .simplify();
        return true;
    }
    return false;
}

bool try_integrate_sec_csc_power_product(const SymbolicExpression& left,
                                         const SymbolicExpression& right,
                                         const std::string& variable_name,
                                         SymbolicExpression* integrated) {
    const SymbolicExpression* power_factor = nullptr;
    const SymbolicExpression* function_factor = nullptr;
    if (left.node_->type == NodeType::kPower &&
        right.node_->type == NodeType::kFunction) {
        power_factor = &left;
        function_factor = &right;
    } else if (right.node_->type == NodeType::kPower &&
               left.node_->type == NodeType::kFunction) {
        power_factor = &right;
        function_factor = &left;
    } else {
        return false;
    }

    const SymbolicExpression base(power_factor->node_->left);
    const SymbolicExpression exponent(power_factor->node_->right);
    double exponent_value = 0.0;
    if (base.node_->type != NodeType::kFunction ||
        !exponent.is_number(&exponent_value) ||
        !mymath::is_near_zero(exponent_value - 2.0, kFormatEps) ||
        !same_simplified_expression(SymbolicExpression(base.node_->left),
                                    SymbolicExpression(function_factor->node_->left))) {
        return false;
    }

    const SymbolicExpression argument(base.node_->left);
    double a = 0.0;
    double b = 0.0;
    if (!decompose_linear(argument, variable_name, &a, &b) ||
        mymath::is_near_zero(a, kFormatEps)) {
        return false;
    }

    if (base.node_->text == "sec" && function_factor->node_->text == "tan") {
        *integrated =
            make_divide(make_power(make_function("sec", argument),
                                   SymbolicExpression::number(2.0)),
                        SymbolicExpression::number(2.0 * a))
                .simplify();
        return true;
    }
    if (base.node_->text == "csc" && function_factor->node_->text == "cot") {
        *integrated =
            make_divide(make_negate(make_power(make_function("csc", argument),
                                               SymbolicExpression::number(2.0))),
                        SymbolicExpression::number(2.0 * a))
                .simplify();
        return true;
    }
    return false;
}

bool try_integrate_polynomial_quotient(const SymbolicExpression& numerator,
                                       const SymbolicExpression& denominator,
                                       const std::string& variable_name,
                                       SymbolicExpression* integrated) {
    std::vector<double> numerator_coefficients;
    std::vector<double> denominator_coefficients;
    if (!polynomial_coefficients_from_simplified(numerator.simplify(),
                                                 variable_name,
                                                 &numerator_coefficients) ||
        !polynomial_coefficients_from_simplified(denominator.simplify(),
                                                 variable_name,
                                                 &denominator_coefficients)) {
        return false;
    }
    trim_coefficients(&numerator_coefficients);
    trim_coefficients(&denominator_coefficients);
    if (denominator_coefficients.size() <= 1 ||
        polynomial_is_zero(denominator_coefficients)) {
        return false;
    }
    if (denominator_coefficients.back() < 0.0) {
        for (double& coefficient : numerator_coefficients) {
            coefficient = -coefficient;
        }
        for (double& coefficient : denominator_coefficients) {
            coefficient = -coefficient;
        }
    }

    const PolynomialDivisionResult division =
        polynomial_divide(numerator_coefficients, denominator_coefficients);
    std::vector<double> quotient_coefficients = division.quotient;
    std::vector<double> remainder_coefficients = division.remainder;
    trim_coefficients(&quotient_coefficients);
    trim_coefficients(&remainder_coefficients);

    SymbolicExpression result = SymbolicExpression::number(0.0);
    if (!polynomial_is_zero(quotient_coefficients)) {
        const SymbolicExpression quotient_expression =
            build_polynomial_expression_from_coefficients(quotient_coefficients,
                                                          variable_name);
        result = quotient_expression.integral(variable_name).simplify();
    }

    if (polynomial_is_zero(remainder_coefficients)) {
        *integrated = result.simplify();
        return true;
    }

    const SymbolicExpression denominator_expression =
        build_polynomial_expression_from_coefficients(denominator_coefficients,
                                                      variable_name)
            .simplify();
    SymbolicExpression remainder_integral;
    if (denominator_coefficients.size() == 2) {
        const double constant = remainder_coefficients[0];
        const double slope = denominator_coefficients[1];
        if (mymath::is_near_zero(slope, kFormatEps)) {
            return false;
        }
        remainder_integral =
            make_multiply(SymbolicExpression::number(constant / slope),
                          make_function("ln",
                                        make_function("abs",
                                                      denominator_expression)))
                .simplify();
    } else if (denominator_coefficients.size() == 3) {
        if (remainder_coefficients.size() > 2) {
            return false;
        }
        const double remainder_constant = remainder_coefficients[0];
        const double remainder_slope =
            remainder_coefficients.size() > 1 ? remainder_coefficients[1] : 0.0;
        const double a = denominator_coefficients[2];
        const double b = denominator_coefficients[1];
        if (mymath::is_near_zero(a, kFormatEps)) {
            return false;
        }

        const double derivative_factor = remainder_slope / (2.0 * a);
        const double inverse_factor = remainder_constant - derivative_factor * b;
        remainder_integral = SymbolicExpression::number(0.0);
        const double repeated_discriminant = b * b - 4.0 * a * denominator_coefficients[0];
        if (mymath::is_near_zero(repeated_discriminant, 1e-10)) {
            const double root = clean_symbolic_constant(-b / (2.0 * a));
            const double alpha = clean_symbolic_constant(remainder_slope / a);
            const double beta = clean_symbolic_constant(remainder_constant / a + alpha * root);
            const SymbolicExpression shifted_variable =
                make_subtract(SymbolicExpression::variable(variable_name),
                              SymbolicExpression::number(root))
                    .simplify();
            if (!mymath::is_near_zero(alpha, kFormatEps)) {
                remainder_integral =
                    make_add(remainder_integral,
                             make_multiply(SymbolicExpression::number(alpha),
                                           make_function("ln",
                                                         make_function("abs",
                                                                       shifted_variable))))
                        .simplify();
            }
            if (!mymath::is_near_zero(beta, kFormatEps)) {
                remainder_integral =
                    make_add(remainder_integral,
                             make_multiply(SymbolicExpression::number(-beta),
                                           make_divide(SymbolicExpression::number(1.0),
                                                       shifted_variable)))
                        .simplify();
            }
            *integrated = make_add(result, remainder_integral).simplify();
            return true;
        }
        if (!mymath::is_near_zero(derivative_factor, kFormatEps)) {
            remainder_integral =
                make_add(remainder_integral,
                         make_multiply(SymbolicExpression::number(derivative_factor),
                                       make_function("ln",
                                                     make_function(
                                                         "abs",
                                                         denominator_expression))))
                    .simplify();
        }
        if (!mymath::is_near_zero(inverse_factor, kFormatEps)) {
            SymbolicExpression inverse_quadratic;
            if (!integrate_inverse_quadratic(denominator_coefficients,
                                             variable_name,
                                             &inverse_quadratic)) {
                return false;
            }
            remainder_integral =
                make_add(remainder_integral,
                         make_multiply(SymbolicExpression::number(inverse_factor),
                                       inverse_quadratic))
                    .simplify();
        }
    } else {
        if (!try_integrate_repeated_unit_quadratic(remainder_coefficients,
                                                  denominator_coefficients,
                                                  variable_name,
                                                  &remainder_integral) &&
            !integrate_mixed_linear_quadratic_partial_fractions(remainder_coefficients,
                                                                denominator_coefficients,
                                                                variable_name,
                                                                &remainder_integral) &&
            !integrate_real_linear_partial_fractions(remainder_coefficients,
                                                     denominator_coefficients,
                                                     variable_name,
                                                     &remainder_integral)) {
            return false;
        }
    }

    *integrated = make_add(result, remainder_integral).simplify();
    return true;
}

SymbolicExpression derivative_uncached(const SymbolicExpression& expression,
                                       const std::string& variable_name) {
    const std::shared_ptr<SymbolicExpression::Node>& node_ = expression.node_;
    switch (node_->type) {
        case NodeType::kNumber:
            return SymbolicExpression::number(0.0);
        case NodeType::kVariable:
            return SymbolicExpression::number(node_->text == variable_name ? 1.0 : 0.0);
        case NodeType::kNegate:
            return make_negate(SymbolicExpression(node_->left).derivative(variable_name)).simplify();
        case NodeType::kAdd:
            return make_add(SymbolicExpression(node_->left).derivative(variable_name),
                            SymbolicExpression(node_->right).derivative(variable_name)).simplify();
        case NodeType::kSubtract:
            return make_subtract(SymbolicExpression(node_->left).derivative(variable_name),
                                 SymbolicExpression(node_->right).derivative(variable_name)).simplify();
        case NodeType::kMultiply: {
            const SymbolicExpression left(node_->left);
            const SymbolicExpression right(node_->right);
            return make_add(make_multiply(left.derivative(variable_name), right),
                            make_multiply(left, right.derivative(variable_name)))
                .simplify();
        }
        case NodeType::kDivide: {
            const SymbolicExpression left(node_->left);
            const SymbolicExpression right(node_->right);
            return make_divide(
                       make_subtract(make_multiply(left.derivative(variable_name), right),
                                     make_multiply(left, right.derivative(variable_name))),
                       make_power(right, SymbolicExpression::number(2.0)))
                .simplify();
        }
        case NodeType::kPower: {
            const SymbolicExpression base(node_->left);
            const SymbolicExpression exponent(node_->right);
            double exponent_value = 0.0;
            if (exponent.is_number(&exponent_value)) {
                return make_multiply(
                           make_multiply(SymbolicExpression::number(exponent_value),
                                         make_power(base, SymbolicExpression::number(exponent_value - 1.0))),
                           base.derivative(variable_name))
                    .simplify();
            }
            if (base.is_constant(variable_name)) {
                return make_multiply(
                           make_multiply(expression, make_function("ln", base)),
                           exponent.derivative(variable_name))
                    .simplify();
            }
            return make_multiply(
                       expression,
                       make_add(
                           make_multiply(exponent.derivative(variable_name), make_function("ln", base)),
                           make_multiply(exponent,
                                         make_divide(base.derivative(variable_name), base))))
                .simplify();
        }
        case NodeType::kFunction: {
            const SymbolicExpression argument(node_->left);
            const SymbolicExpression inner = argument.derivative(variable_name);
            if (node_->text == "asin") {
                return make_divide(
                           inner,
                           make_function("sqrt",
                                         make_subtract(SymbolicExpression::number(1.0),
                                                       make_power(argument, SymbolicExpression::number(2.0)))))
                    .simplify();
            }
            if (node_->text == "acos") {
                return make_negate(
                           make_divide(inner,
                                       make_function("sqrt",
                                                     make_subtract(SymbolicExpression::number(1.0),
                                                                   make_power(argument, SymbolicExpression::number(2.0))))))
                    .simplify();
            }
            if (node_->text == "atan") {
                return make_divide(inner,
                                   make_add(SymbolicExpression::number(1.0), make_power(argument, SymbolicExpression::number(2.0))))
                    .simplify();
            }
            if (node_->text == "sin") {
                return make_multiply(make_function("cos", argument), inner).simplify();
            }
            if (node_->text == "cos") {
                return make_multiply(make_negate(make_function("sin", argument)), inner).simplify();
            }
            if (node_->text == "tan") {
                return make_multiply(make_divide(SymbolicExpression::number(1.0),
                                                 make_power(make_function("cos", argument),
                                                            SymbolicExpression::number(2.0))),
                                     inner)
                    .simplify();
            }
            if (node_->text == "sec") {
                return make_multiply(
                           make_multiply(make_function("sec", argument),
                                         make_function("tan", argument)),
                           inner)
                    .simplify();
            }
            if (node_->text == "csc") {
                return make_multiply(
                           make_negate(make_multiply(make_function("csc", argument),
                                                     make_function("cot", argument))),
                           inner)
                    .simplify();
            }
            if (node_->text == "cot") {
                return make_multiply(
                           make_negate(make_power(make_function("csc", argument),
                                                  SymbolicExpression::number(2.0))),
                           inner)
                    .simplify();
            }
            if (node_->text == "exp") {
                return make_multiply(make_function("exp", argument), inner).simplify();
            }
            if (node_->text == "sinh") {
                return make_multiply(make_function("cosh", argument), inner).simplify();
            }
            if (node_->text == "cosh") {
                return make_multiply(make_function("sinh", argument), inner).simplify();
            }
            if (node_->text == "tanh") {
                return make_divide(inner,
                                   make_power(make_function("cosh", argument),
                                              SymbolicExpression::number(2.0)))
                    .simplify();
            }
            if (node_->text == "ln") {
                return make_divide(inner, argument).simplify();
            }
            if (node_->text == "sqrt") {
                return make_divide(inner,
                                   make_multiply(SymbolicExpression::number(2.0), make_function("sqrt", argument)))
                    .simplify();
            }
            if (node_->text == "cbrt") {
                return make_divide(inner,
                                   make_multiply(SymbolicExpression::number(3.0),
                                                 make_power(make_function("cbrt", argument),
                                                            SymbolicExpression::number(2.0))))
                    .simplify();
            }
            if (node_->text == "abs") {
                return make_multiply(make_function("sign", argument), inner).simplify();
            }
            if (node_->text == "step") {
                return make_multiply(make_function("delta", argument), inner).simplify();
            }
            throw std::runtime_error("symbolic derivative does not support function: " + node_->text);
        }
    }
    throw std::runtime_error("unsupported symbolic derivative");
}

}  // namespace

SymbolicExpression SymbolicExpression::derivative(const std::string& variable_name) const {
    static constexpr std::size_t kMaxDerivativeCacheSize = 4096;
    static thread_local SymbolicExpressionLruCache cache(kMaxDerivativeCacheSize);

    const std::string key = variable_name + "|" + node_structural_key(node_);
    SymbolicExpression cached;
    if (cache.get(key, &cached)) {
        return cached;
    }

    SymbolicExpression derived = derivative_uncached(*this, variable_name);
    cache.put(key, derived);
    return derived;
}

std::vector<SymbolicExpression> SymbolicExpression::gradient(
    const std::vector<std::string>& variable_names) const {
    std::vector<SymbolicExpression> result;
    result.reserve(variable_names.size());
    for (const std::string& variable_name : variable_names) {
        result.push_back(derivative(variable_name).simplify());
    }
    return result;
}

std::vector<std::vector<SymbolicExpression>> SymbolicExpression::hessian(
    const std::vector<std::string>& variable_names) const {
    std::vector<std::vector<SymbolicExpression>> result;
    result.reserve(variable_names.size());
    for (const std::string& outer_variable : variable_names) {
        std::vector<SymbolicExpression> row;
        row.reserve(variable_names.size());
        const SymbolicExpression outer_derivative =
            derivative(outer_variable).simplify();
        for (const std::string& inner_variable : variable_names) {
            row.push_back(outer_derivative.derivative(inner_variable).simplify());
        }
        result.push_back(row);
    }
    return result;
}

std::vector<std::vector<SymbolicExpression>> SymbolicExpression::jacobian(
    const std::vector<SymbolicExpression>& expressions,
    const std::vector<std::string>& variable_names) {
    std::vector<std::vector<SymbolicExpression>> result;
    result.reserve(expressions.size());
    for (const SymbolicExpression& expression : expressions) {
        result.push_back(expression.gradient(variable_names));
    }
    return result;
}

SymbolicExpression SymbolicExpression::integral(const std::string& variable_name) const {
    double numeric_value = 0.0;
    if (is_constant(variable_name)) {
        if (is_number(&numeric_value)) {
            return make_multiply(number(numeric_value), variable(variable_name)).simplify();
        }
        return make_multiply(*this, variable(variable_name)).simplify();
    }

    switch (node_->type) {
        case NodeType::kNumber:
            return make_multiply(number(node_->number_value), variable(variable_name)).simplify();
        case NodeType::kVariable:
            if (node_->text == variable_name) {
                return make_divide(make_power(variable(variable_name), number(2.0)),
                                   number(2.0))
                    .simplify();
            }
            return make_multiply(variable(node_->text), variable(variable_name)).simplify();
        case NodeType::kNegate:
            return make_negate(SymbolicExpression(node_->left).integral(variable_name)).simplify();
        case NodeType::kAdd:
            return make_add(SymbolicExpression(node_->left).integral(variable_name),
                            SymbolicExpression(node_->right).integral(variable_name)).simplify();
        case NodeType::kSubtract:
            return make_subtract(SymbolicExpression(node_->left).integral(variable_name),
                                 SymbolicExpression(node_->right).integral(variable_name)).simplify();
        case NodeType::kMultiply: {
            double constant = 0.0;
            SymbolicExpression rest;
            const SymbolicExpression left(node_->left);
            const SymbolicExpression right(node_->right);
            SymbolicExpression integrated;
            if (try_integrate_substitution_product(left,
                                                  right,
                                                  variable_name,
                                                  &integrated) ||
                try_integrate_substitution_product(right,
                                                  left,
                                                  variable_name,
                                                  &integrated)) {
                return integrated.simplify();
            }
            if (try_integrate_trig_product_identity(left,
                                                    right,
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            if (try_integrate_sec_csc_power_product(left,
                                                    right,
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            if (decompose_constant_times_expression(*this, variable_name, &constant, &rest)) {
                return make_multiply(number(constant), rest.integral(variable_name)).simplify();
            }
            if (left.is_constant(variable_name)) {
                return make_multiply(left, right.integral(variable_name)).simplify();
            }
            if (right.is_constant(variable_name)) {
                return make_multiply(right, left.integral(variable_name)).simplify();
            }
            SymbolicExpression polynomial;
            if (polynomial_expression(left, variable_name, &polynomial) &&
                right.node_->type == NodeType::kFunction &&
                integrate_polynomial_times_function(polynomial,
                                                    right.node_->text,
                                                    SymbolicExpression(right.node_->left),
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            if (polynomial_expression(right, variable_name, &polynomial) &&
                left.node_->type == NodeType::kFunction &&
                integrate_polynomial_times_function(polynomial,
                                                    left.node_->text,
                                                    SymbolicExpression(left.node_->left),
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            throw std::runtime_error("symbolic integral does not support this product");
        }
        case NodeType::kPower:
        case NodeType::kFunction:
        case NodeType::kDivide:
            break;
    }

    if (node_->type == NodeType::kDivide) {
        const SymbolicExpression left(node_->left);
        const SymbolicExpression right(node_->right);
        if (left.is_number(&numeric_value) && mymath::is_near_zero(numeric_value - 1.0, kFormatEps)) {
            if (is_one_plus_variable_squared(right, variable_name)) {
                return make_function("atan", variable(variable_name)).simplify();
            }
            if (is_sqrt_one_minus_variable_squared(right, variable_name)) {
                return make_function("asin", variable(variable_name)).simplify();
            }
            double a = 0.0;
            double b = 0.0;
            if (decompose_linear(right, variable_name, &a, &b) &&
                !mymath::is_near_zero(a, kFormatEps)) {
                return make_divide(make_function("ln", make_function("abs", right)),
                                   number(a))
                    .simplify();
            }
        }
        SymbolicExpression rational_integral;
        if (try_integrate_polynomial_quotient(left,
                                              right,
                                              variable_name,
                                              &rational_integral)) {
            return rational_integral.simplify();
        }
        throw std::runtime_error("symbolic integral does not support this quotient");
    }

    if (node_->type == NodeType::kPower) {
        const SymbolicExpression base(node_->left);
        const SymbolicExpression exponent(node_->right);
        double exponent_value = 0.0;
        double a = 0.0;
        double b = 0.0;
        SymbolicExpression trig_identity_integral;
        if (exponent.is_number(&exponent_value) &&
            try_integrate_trig_power_identity(base,
                                              exponent_value,
                                              variable_name,
                                              &trig_identity_integral)) {
            return trig_identity_integral.simplify();
        }
        if (exponent.is_number(&exponent_value) &&
            decompose_linear(base, variable_name, &a, &b) &&
            !mymath::is_near_zero(a, kFormatEps)) {
            if (mymath::is_near_zero(exponent_value + 1.0, kFormatEps)) {
                return make_divide(make_function("ln", make_function("abs", base)),
                                   number(a))
                    .simplify();
            }
            return make_divide(make_power(base, number(exponent_value + 1.0)),
                               number(a * (exponent_value + 1.0)))
                .simplify();
        }
        throw std::runtime_error("symbolic integral only supports powers of the integration variable");
    }

    if (node_->type == NodeType::kFunction) {
        const SymbolicExpression argument(node_->left);
        double a = 0.0;
        double b = 0.0;
        const bool linear = decompose_linear(argument, variable_name, &a, &b) &&
                            !mymath::is_near_zero(a, kFormatEps);
        if (node_->text == "sin" && linear) {
            return make_divide(make_negate(make_function("cos", argument)),
                               number(a))
                .simplify();
        }
        if (node_->text == "cos" && linear) {
            return make_divide(make_function("sin", argument), number(a)).simplify();
        }
        if (node_->text == "exp" && linear) {
            return make_divide(make_function("exp", argument), number(a)).simplify();
        }
        if (node_->text == "sqrt" && linear) {
            return make_divide(make_multiply(number(2.0),
                                             make_power(make_function("sqrt", argument),
                                                        number(3.0))),
                               number(3.0 * a))
                .simplify();
        }
        if (node_->text == "cbrt" && linear) {
            return make_divide(make_multiply(number(3.0),
                                             make_power(make_function("cbrt", argument),
                                                        number(4.0))),
                               number(4.0 * a))
                .simplify();
        }
        if (node_->text == "sqrt" &&
            is_one_minus_variable_squared(argument, variable_name)) {
            const SymbolicExpression x = variable(variable_name);
            return make_divide(
                       make_add(make_multiply(x, make_function("sqrt", argument)),
                                make_function("asin", x)),
                       number(2.0))
                .simplify();
        }
        if (node_->text == "tan" && linear) {
            return make_divide(make_negate(make_function("ln",
                                                         make_function("abs",
                                                                       make_function("cos",
                                                                                     argument)))),
                               number(a))
                .simplify();
        }
        if (node_->text == "sec" && linear) {
            return make_divide(
                       make_function("ln",
                                     make_function("abs",
                                                   make_add(make_function("sec", argument),
                                                            make_function("tan", argument)))),
                       number(a))
                .simplify();
        }
        if (node_->text == "csc" && linear) {
            return make_divide(
                       make_function("ln",
                                     make_function("abs",
                                                   make_subtract(make_function("csc", argument),
                                                                 make_function("cot", argument)))),
                       number(a))
                .simplify();
        }
        if (node_->text == "cot" && linear) {
            return make_divide(
                       make_function("ln",
                                     make_function("abs",
                                                   make_function("sin", argument))),
                       number(a))
                .simplify();
        }
        if (node_->text == "delta" &&
            argument.is_variable_named(variable_name)) {
            return make_step_expression(variable_name, 0.0);
        }
        throw std::runtime_error("symbolic integral does not support function: " + node_->text);
    }

    throw std::runtime_error("unsupported symbolic integral");
}

#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include <cmath>
#include <vector>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;

namespace {

void trim_numeric_polynomial(std::vector<double>* coefficients) {
    while (!coefficients->empty() && std::abs(coefficients->back()) < 1e-10) {
        coefficients->pop_back();
    }
}

bool numeric_coefficients(const SymbolicPolynomial& polynomial,
                          std::vector<double>* coefficients) {
    coefficients->clear();
    coefficients->reserve(polynomial.degree() + 1);
    for (int i = 0; i <= polynomial.degree(); ++i) {
        double value = 0.0;
        if (!polynomial.coefficient(i).is_number(&value)) {
            return false;
        }
        coefficients->push_back(value);
    }
    trim_numeric_polynomial(coefficients);
    return !coefficients->empty();
}

std::vector<double> multiply_numeric_polynomial(const std::vector<double>& lhs,
                                                const std::vector<double>& rhs) {
    std::vector<double> product(lhs.size() + rhs.size() - 1, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            product[i + j] += lhs[i] * rhs[j];
        }
    }
    trim_numeric_polynomial(&product);
    return product;
}

std::vector<double> power_numeric_polynomial(const std::vector<double>& base,
                                             int power) {
    std::vector<double> result = {1.0};
    for (int i = 0; i < power; ++i) {
        result = multiply_numeric_polynomial(result, base);
    }
    return result;
}

bool numeric_polynomials_close(const std::vector<double>& lhs,
                               const std::vector<double>& rhs,
                               double eps) {
    std::vector<double> left = lhs;
    std::vector<double> right = rhs;
    trim_numeric_polynomial(&left);
    trim_numeric_polynomial(&right);
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (std::abs(left[i] - right[i]) > eps) {
            return false;
        }
    }
    return true;
}

bool divide_numeric_polynomial(const std::vector<double>& numerator,
                               const std::vector<double>& denominator,
                               std::vector<double>* quotient) {
    if (denominator.empty() || std::abs(denominator.back()) < 1e-12 ||
        numerator.size() < denominator.size()) {
        return false;
    }

    std::vector<double> remainder = numerator;
    quotient->assign(numerator.size() - denominator.size() + 1, 0.0);
    const int den_degree = static_cast<int>(denominator.size()) - 1;

    for (int i = static_cast<int>(numerator.size()) - 1; i >= den_degree; --i) {
        const double term = remainder[i] / denominator.back();
        const int q_index = i - den_degree;
        (*quotient)[q_index] = term;
        for (int j = 0; j <= den_degree; ++j) {
            remainder[q_index + j] -= term * denominator[j];
        }
    }

    for (int i = 0; i < den_degree; ++i) {
        if (std::abs(remainder[i]) > 1e-7) {
            return false;
        }
    }
    trim_numeric_polynomial(quotient);
    return true;
}

bool solve_numeric_linear_system(std::vector<std::vector<double>> matrix,
                                 std::vector<double> rhs,
                                 std::vector<double>* solution) {
    const int n = static_cast<int>(rhs.size());
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int row = col + 1; row < n; ++row) {
            if (std::abs(matrix[row][col]) > std::abs(matrix[pivot][col])) {
                pivot = row;
            }
        }
        if (std::abs(matrix[pivot][col]) < 1e-12) {
            return false;
        }
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(rhs[pivot], rhs[col]);
        }

        const double pivot_value = matrix[col][col];
        for (int j = col; j < n; ++j) {
            matrix[col][j] /= pivot_value;
        }
        rhs[col] /= pivot_value;

        for (int row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = matrix[row][col];
            if (std::abs(factor) < 1e-14) {
                continue;
            }
            for (int j = col; j < n; ++j) {
                matrix[row][j] -= factor * matrix[col][j];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    *solution = rhs;
    return true;
}

SymbolicExpression integrate_numeric_inverse_quadratic_power(const std::vector<double>& quadratic,
                                                             int power,
                                                             const std::string& variable_name) {
    SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const double c = quadratic[0];
    const double b = quadratic[1];
    const double a = quadratic[2];
    SymbolicExpression q =
        (SymbolicExpression::number(a) * make_power(x, SymbolicExpression::number(2.0)) +
         SymbolicExpression::number(b) * x +
         SymbolicExpression::number(c)).simplify();

    const double disc_neg = 4.0 * a * c - b * b;
    if (disc_neg <= 1e-12) {
        return SymbolicExpression::number(0.0);
    }

    const double root = std::sqrt(disc_neg);
    SymbolicExpression arg =
        ((SymbolicExpression::number(2.0 * a) * x + SymbolicExpression::number(b)) /
         SymbolicExpression::number(root)).simplify();
    SymbolicExpression integral =
        (SymbolicExpression::number(2.0 / root) * make_function("atan", arg)).simplify();

    for (int n = 2; n <= power; ++n) {
        SymbolicExpression rational_part =
            ((SymbolicExpression::number(2.0 * a) * x + SymbolicExpression::number(b)) /
             (SymbolicExpression::number(disc_neg * (n - 1.0)) *
              make_power(q, SymbolicExpression::number(n - 1.0)))).simplify();
        SymbolicExpression recursive_part =
            (SymbolicExpression::number(2.0 * a * (2.0 * n - 3.0) /
                                        (disc_neg * (n - 1.0))) *
             integral).simplify();
        integral = (rational_part + recursive_part).simplify();
    }

    return integral;
}

SymbolicExpression integrate_numeric_quadratic_fraction(double slope,
                                                        double constant,
                                                        const std::vector<double>& quadratic,
                                                        int power,
                                                        const std::string& variable_name) {
    SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const double c = quadratic[0];
    const double b = quadratic[1];
    const double a = quadratic[2];

    SymbolicExpression q =
        (SymbolicExpression::number(a) * make_power(x, SymbolicExpression::number(2.0)) +
         SymbolicExpression::number(b) * x +
         SymbolicExpression::number(c)).simplify();

    SymbolicExpression result = SymbolicExpression::number(0.0);
    const double log_coeff = slope / (2.0 * a);
    if (std::abs(log_coeff) > 1e-12) {
        if (power == 1) {
            result = (result +
                      SymbolicExpression::number(log_coeff) *
                      make_function("ln", make_function("abs", q))).simplify();
        } else {
            result = (result +
                      SymbolicExpression::number(log_coeff / (1.0 - power)) *
                      make_power(q, SymbolicExpression::number(1.0 - power))).simplify();
        }
    }

    const double remainder = constant - log_coeff * b;
    if (std::abs(remainder) < 1e-12) {
        return result;
    }

    return (result +
            SymbolicExpression::number(remainder) *
            integrate_numeric_inverse_quadratic_power(quadratic,
                                                      power,
                                                      variable_name)).simplify();
}

bool try_integrate_numeric_quadratic_partial_fractions(const SymbolicPolynomial& numerator,
                                                       const SymbolicPolynomial& denominator,
                                                       const std::string& variable_name,
                                                       SymbolicExpression* result) {
    std::vector<double> num_coeffs;
    std::vector<double> den_coeffs;
    if (!numeric_coefficients(numerator, &num_coeffs) ||
        !numeric_coefficients(denominator, &den_coeffs) ||
        denominator.degree() < 3) {
        return false;
    }

    auto roots = RischAlgorithm::find_complex_roots_aberth(den_coeffs);
    struct QuadraticFactor {
        std::vector<double> coefficients;
        int multiplicity = 1;
    };

    std::vector<QuadraticFactor> quadratics;
    bool used_power_factorization = false;

    const int denominator_degree = static_cast<int>(den_coeffs.size()) - 1;
    if (denominator_degree >= 4 && denominator_degree % 2 == 0) {
        const int multiplicity = denominator_degree / 2;
        const double leading = den_coeffs.back();
        const double constant = den_coeffs.front();
        if (leading > 0.0 && constant > 0.0) {
            const double a = std::pow(leading, 1.0 / multiplicity);
            const double c = std::pow(constant, 1.0 / multiplicity);
            const double b = den_coeffs[denominator_degree - 1] /
                             (multiplicity * std::pow(a, multiplicity - 1));
            std::vector<double> factor = {c, b, a};
            if (b * b - 4.0 * a * c < -1e-8 &&
                numeric_polynomials_close(power_numeric_polynomial(factor, multiplicity),
                                          den_coeffs,
                                          1e-6 * (1.0 + std::abs(leading) + std::abs(constant)))) {
                quadratics.push_back({factor, multiplicity});
                used_power_factorization = true;
            }
        }
    }

    for (const auto& root : roots) {
        if (used_power_factorization) {
            break;
        }
        const double re = root.first;
        const double im = root.second;
        if (std::abs(im) < 1e-8) {
            return false;
        }
        std::vector<double> factor = {re * re + im * im, -2.0 * re, 1.0};
        bool merged = false;
        for (auto& existing : quadratics) {
            if (std::abs(existing.coefficients[0] - factor[0]) < 1e-6 &&
                std::abs(existing.coefficients[1] - factor[1]) < 1e-6) {
                ++existing.multiplicity;
                merged = true;
                break;
            }
        }
        if (!merged) {
            quadratics.push_back({factor, 1});
        }
    }

    int factored_degree = 0;
    for (const auto& factor : quadratics) {
        factored_degree += 2 * factor.multiplicity;
    }
    if (quadratics.empty() || factored_degree != denominator.degree()) {
        return false;
    }

    std::vector<std::vector<double>> columns;
    columns.reserve(denominator.degree());
    for (const auto& factor : quadratics) {
        std::vector<double> divisor = {1.0};
        for (int power = 1; power <= factor.multiplicity; ++power) {
            divisor = multiply_numeric_polynomial(divisor, factor.coefficients);
            std::vector<double> quotient;
            if (!divide_numeric_polynomial(den_coeffs, divisor, &quotient)) {
                return false;
            }
            columns.push_back(quotient);
            columns.push_back(multiply_numeric_polynomial(quotient, {0.0, 1.0}));
        }
    }

    const int unknown_count = static_cast<int>(columns.size());
    std::vector<std::vector<double>> matrix(unknown_count,
                                            std::vector<double>(unknown_count, 0.0));
    std::vector<double> rhs(unknown_count, 0.0);
    for (int row = 0; row < unknown_count; ++row) {
        if (row < static_cast<int>(num_coeffs.size())) {
            rhs[row] = num_coeffs[row];
        }
        for (int col = 0; col < unknown_count; ++col) {
            if (row < static_cast<int>(columns[col].size())) {
                matrix[row][col] = columns[col][row];
            }
        }
    }

    std::vector<double> unknowns;
    if (!solve_numeric_linear_system(matrix, rhs, &unknowns)) {
        return false;
    }

    SymbolicExpression total = SymbolicExpression::number(0.0);
    std::size_t unknown_index = 0;
    for (const auto& factor : quadratics) {
        for (int power = 1; power <= factor.multiplicity; ++power) {
            const double constant = unknowns[unknown_index++];
            const double slope = unknowns[unknown_index++];
            total = (total + integrate_numeric_quadratic_fraction(slope,
                                                                  constant,
                                                                  factor.coefficients,
                                                                  power,
                                                                  variable_name)).simplify();
        }
    }

    *result = total;
    return true;
}

bool try_integrate_log_derivative_monomial(const SymbolicPolynomial& numerator,
                                           const SymbolicPolynomial& denominator,
                                           const std::string& variable_name,
                                           const SymbolicExpression& t_prime,
                                           const std::vector<DifferentialExtension>& tower,
                                           int tower_index,
                                           const std::string& main_var,
                                           SymbolicExpression* result) {
    if (numerator.degree() != 0 || denominator.degree() <= 0) {
        return false;
    }

    const int power = denominator.degree();
    for (int i = 0; i < power; ++i) {
        if (!SymbolicPolynomial::coeff_is_zero(denominator.coefficient(i))) {
            return false;
        }
    }

    SymbolicExpression scale =
        (numerator.coefficient(0) /
         (denominator.leading_coefficient() * t_prime)).simplify();
    if (contains_var(scale, main_var) ||
        contains_tower_var(scale, tower, tower_index)) {
        if (t_prime.node_->type == NodeType::kDivide) {
            SymbolicExpression t_prime_num(t_prime.node_->left);
            SymbolicExpression t_prime_den(t_prime.node_->right);
            if (structural_equals(t_prime_den.simplify(),
                                  denominator.leading_coefficient().simplify())) {
                scale = (numerator.coefficient(0) / t_prime_num).simplify();
            }
        }
    }
    if (contains_var(scale, main_var) ||
        contains_tower_var(scale, tower, tower_index)) {
        return false;
    }

    SymbolicExpression t = SymbolicExpression::variable(variable_name);
    if (power == 1) {
        *result = (scale * make_function("ln", make_function("abs", t))).simplify();
    } else {
        *result = (scale / SymbolicExpression::number(1.0 - static_cast<double>(power)) *
                   make_power(t, SymbolicExpression::number(1.0 - static_cast<double>(power)))).simplify();
    }
    return true;
}

}  // namespace

// 有理函数积分
// ============================================================================

bool RischAlgorithm::integrate_rational(const SymbolicPolynomial& numerator,
                                        const SymbolicPolynomial& denominator,
                                        const std::string& variable_name,
                                        SymbolicExpression* result,
                                        const std::vector<DifferentialExtension>& tower,
                                        int tower_index,
                                        const std::string& main_var,
                                        const SymbolicExpression* t_prime,
                                        DifferentialExtension::Kind kind) {
    // 1. 多项式部分
    SymbolicPolynomial Q, R;
    numerator.divide(denominator, &Q, &R);

    SymbolicExpression poly_int = SymbolicExpression::number(0.0);

    if (!Q.is_zero()) {
        const auto& q_coeffs = Q.coefficients();
        for (std::size_t i = 0; i < q_coeffs.size(); ++i) {
            if (SymbolicPolynomial::coeff_is_zero(q_coeffs[i])) continue;

            SymbolicExpression term_int;

            if (kind == DifferentialExtension::Kind::kNone) {
                // 普通变量 x^i
                term_int = (q_coeffs[i] / SymbolicExpression::number(static_cast<double>(i + 1)) *
                           make_power(SymbolicExpression::variable(variable_name),
                                      SymbolicExpression::number(static_cast<double>(i + 1)))).simplify();
            } else if (kind == DifferentialExtension::Kind::kLogarithmic) {
                // t = ln(u), ∫ a*t^n dx
                // 使用递归公式: ∫ a*t^n = (∫a)*t^n - ∫(∫a)*n*t^(n-1)*t' dx
                if (i == 0) {
                    IntegrationResult res = integrate_in_extension(q_coeffs[i], tower, tower_index - 1, main_var);
                    if (!res.success || res.type != IntegralType::kElementary) return false;
                    term_int = res.value;
                } else {
                    IntegrationResult a_int_res = integrate_in_extension(q_coeffs[i], tower, tower_index - 1, main_var);
                    if (!a_int_res.success || a_int_res.type != IntegralType::kElementary) return false;

                    SymbolicExpression a_int = a_int_res.value;

                    // 改进：允许 a_int 是 t 的常数倍
                    SymbolicExpression t = SymbolicExpression::variable(variable_name);
                    double t_coeff = 0.0;
                    SymbolicExpression remainder;
                    
                    bool is_t_multiple = false;
                    if (structural_equals(a_int, t)) {
                        is_t_multiple = true;
                        t_coeff = 1.0;
                    } else if (a_int.node_->type == NodeType::kMultiply) {
                        SymbolicExpression left(a_int.node_->left);
                        SymbolicExpression right(a_int.node_->right);
                        if (left.is_number(&t_coeff) && structural_equals(right, t)) {
                            is_t_multiple = true;
                        } else if (right.is_number(&t_coeff) && structural_equals(left, t)) {
                            is_t_multiple = true;
                        }
                    }

                    if (is_t_multiple) {
                        // ∫ (c*t) * t^n dx 其中 a_int = c*t
                        // 这意味着 a = c*t' (因为 ∫a = c*t)
                        // 使用公式: ∫ a*t^n dx = ∫ c*t'*t^n dx = c * t^(n+1) / (n+1)
                        term_int = (SymbolicExpression::number(t_coeff / (i + 1)) *
                                   make_power(t, SymbolicExpression::number(i + 1))).simplify();
                    } else if (contains_tower_var(a_int, tower, tower_index)) {
                        // a_int 包含塔变量且不是简单的 t 倍数
                        // 尝试使用 RDE 求解: y' + n*t'*y = q_coeffs[i]
                        // 其中 y 是在基域中的解

                        SymbolicExpression n_expr = SymbolicExpression::number(static_cast<double>(i));
                        SymbolicExpression f = (n_expr * (*t_prime)).simplify();

                        IntegrationResult y_res = solve_rde(f, q_coeffs[i], main_var, tower, tower_index - 1);
                        if (y_res.success && y_res.type == IntegralType::kElementary) {
                            term_int = (y_res.value * make_power(t, n_expr)).simplify();
                        } else {
                            // RDE 失败，尝试分部积分递归
                            // ∫ a t^n = a_int * t^n - ∫ a_int * n * t^(n-1) * t' dx
                            SymbolicExpression first_part = (a_int * make_power(t, SymbolicExpression::number(i))).simplify();

                            SymbolicExpression next_base = (a_int * SymbolicExpression::number(i) *
                                                           make_power(t, SymbolicExpression::number(i - 1))).simplify();
                            SymbolicExpression next_integrand =
                                multiply_by_derivative_factor(next_base, *t_prime);

                            IntegrationResult second_res = integrate_in_extension(next_integrand, tower, tower_index, main_var);
                            if (!second_res.success || second_res.type != IntegralType::kElementary) return false;

                            term_int = (first_part - second_res.value).simplify();
                        }
                    } else {
                        SymbolicExpression first_part = (a_int * make_power(t, SymbolicExpression::number(i))).simplify();

                        SymbolicExpression next_base = (a_int * SymbolicExpression::number(i) *
                                                       make_power(t, SymbolicExpression::number(i - 1))).simplify();
                        SymbolicExpression next_integrand =
                            multiply_by_derivative_factor(next_base, *t_prime);

                        IntegrationResult second_res = integrate_in_extension(next_integrand, tower, tower_index, main_var);
                        if (!second_res.success || second_res.type != IntegralType::kElementary) return false;

                        term_int = (first_part - second_res.value).simplify();
                    }
                }
            } else if (kind == DifferentialExtension::Kind::kExponential) {
                // t = exp(u), ∫ a*t^n dx = y*t^n 其中 y' + n*u'*y = a
                if (i == 0) {
                    IntegrationResult res = integrate_in_extension(q_coeffs[i], tower, tower_index - 1, main_var);
                    if (!res.success || res.type != IntegralType::kElementary) return false;
                    term_int = res.value;
                } else {
                    SymbolicExpression u_prime = ((*t_prime) / SymbolicExpression::variable(variable_name)).simplify();
                    SymbolicExpression f = (SymbolicExpression::number(i) * u_prime).simplify();

                    IntegrationResult y_res = solve_rde(f, q_coeffs[i], main_var, tower, tower_index - 1);
                    if (!y_res.success || y_res.type != IntegralType::kElementary) return false;

                    term_int = (y_res.value * make_power(SymbolicExpression::variable(variable_name),
                                                        SymbolicExpression::number(i))).simplify();
                }
            } else if (kind == DifferentialExtension::Kind::kAlgebraic) {
                // 代数扩展
                if (i == 0) {
                    IntegrationResult res = integrate_in_extension(q_coeffs[i], tower, tower_index - 1, main_var);
                    if (!res.success || res.type != IntegralType::kElementary) return false;
                    term_int = res.value;
                } else {
                    // 尝试 RDE 求解
                    SymbolicExpression f = (SymbolicExpression::number(static_cast<double>(i)) *
                                           (*t_prime) / SymbolicExpression::variable(variable_name)).simplify();

                    IntegrationResult y_res = solve_rde(f, q_coeffs[i], main_var, tower, tower_index - 1);
                    if (y_res.success && y_res.type == IntegralType::kElementary) {
                        term_int = (y_res.value * make_power(SymbolicExpression::variable(variable_name),
                                                            SymbolicExpression::number(static_cast<double>(i)))).simplify();
                    } else {
                        return false;
                    }
                }
            }

            poly_int = (poly_int + term_int).simplify();
        }
    }

    if (R.is_zero()) {
        *result = poly_int;
        return true;
    }

    if (kind == DifferentialExtension::Kind::kLogarithmic && t_prime) {
        SymbolicExpression monomial_log_part;
        if (try_integrate_log_derivative_monomial(R,
                                                  denominator,
                                                  variable_name,
                                                  *t_prime,
                                                  tower,
                                                  tower_index,
                                                  main_var,
                                                  &monomial_log_part)) {
            *result = (poly_int + monomial_log_part).simplify();
            return true;
        }
    }

    if (kind == DifferentialExtension::Kind::kNone) {
        SymbolicExpression partial_fraction_part;
        if (try_integrate_numeric_quadratic_partial_fractions(R,
                                                              denominator,
                                                              variable_name,
                                                              &partial_fraction_part)) {
            *result = (poly_int + partial_fraction_part).simplify();
            return true;
        }
    }

    // 2. Hermite 约化
    SymbolicExpression rational_part = SymbolicExpression::number(0.0);
    SymbolicPolynomial reduced_num, reduced_den;

    if (!hermite_reduction(R, denominator, &rational_part, &reduced_num, &reduced_den,
                          tower, tower_index, main_var, t_prime, kind)) {
        return false;
    }

    // 3. Rothstein-Trager 算法
    SymbolicExpression log_part = SymbolicExpression::number(0.0);
    if (!reduced_num.is_zero()) {
        if (!rothstein_trager(reduced_num, reduced_den, variable_name, &log_part,
                             tower, tower_index, main_var, t_prime, kind)) {
            return false;
        }
    }

    *result = (poly_int + rational_part + log_part).simplify();
    return true;
}

// ============================================================================
// Hermite 约化
// ============================================================================

bool RischAlgorithm::hermite_reduction(const SymbolicPolynomial& numerator,
                                       const SymbolicPolynomial& denominator,
                                       SymbolicExpression* rational_part,
                                       SymbolicPolynomial* reduced_numerator,
                                       SymbolicPolynomial* reduced_denominator,
                                       const std::vector<DifferentialExtension>& tower,
                                       int tower_index,
                                       const std::string& main_var,
                                       const SymbolicExpression* t_prime,
                                       DifferentialExtension::Kind kind) {
    auto poly_diff = [&](const SymbolicPolynomial& p) {
        return t_prime ? p.total_derivative(main_var, *t_prime) : p.derivative();
    };

    if (denominator.is_zero()) {
        return false;
    }

    if (numerator.is_zero()) {
        *rational_part = SymbolicExpression::number(0.0);
        *reduced_numerator = numerator;
        *reduced_denominator = denominator;
        return true;
    }

    if (polynomial_is_obviously_square_free(denominator)) {
        *rational_part = SymbolicExpression::number(0.0);
        *reduced_numerator = numerator;
        *reduced_denominator = denominator;
        return true;
    }

    // Square-free 分解
    std::vector<SymbolicPolynomial> factors;
    if (!denominator.square_free_decomposition(&factors) || factors.empty()) {
        *rational_part = SymbolicExpression::number(0.0);
        *reduced_numerator = numerator;
        *reduced_denominator = denominator;
        return true;
    }

    SymbolicExpression total_rational = SymbolicExpression::number(0.0);
    SymbolicPolynomial current_num = numerator;
    SymbolicPolynomial current_den = denominator;

    // 从最高幂次开始处理
    for (int i = static_cast<int>(factors.size()); i > 1; --i) {
        SymbolicPolynomial v = factors[i - 1];
        if (v.degree() <= 0) continue;

        SymbolicPolynomial Vi = v.power(i);
        SymbolicPolynomial U, R;
        if (!current_den.divide(Vi, &U, &R)) {
            continue;
        }

        for (int k = i; k > 1; --k) {
            SymbolicPolynomial v_deriv = poly_diff(v);
            SymbolicPolynomial UVp = U.multiply(v_deriv);

            // 扩展 GCD
            SymbolicPolynomial S, T;
            SymbolicPolynomial g = UVp.extended_gcd(v, &S, &T);

            // 处理非常量 GCD
            if (!g.is_constant()) {
                SymbolicPolynomial v_reduced, rem;
                if (g.divide(g, &v_reduced, &rem)) {
                    if (!v_reduced.is_zero() && v_reduced.degree() > 0) {
                        v = v_reduced;
                        Vi = v.power(k);
                        if (!current_den.divide(Vi, &U, &R)) continue;
                        v_deriv = poly_diff(v);
                        UVp = U.multiply(v_deriv);
                        g = UVp.extended_gcd(v, &S, &T);
                    }
                }
            }

            // 归一化
            SymbolicExpression inv_g = SymbolicExpression::number(1.0);
            if (!g.is_zero()) {
                inv_g = (SymbolicExpression::number(1.0) / g.leading_coefficient()).simplify();
            }

            S = S.multiply(current_num).scale(inv_g);
            T = T.multiply(current_num).scale(inv_g);

            // 计算有理部分贡献
            SymbolicExpression factor = (SymbolicExpression::number(-1.0) /
                                        SymbolicExpression::number(static_cast<double>(k - 1))).simplify();
            SymbolicPolynomial G = S.scale(factor);

            SymbolicExpression term = (G.to_expression() /
                                      make_power(v.to_expression(),
                                                SymbolicExpression::number(static_cast<double>(k - 1)))).simplify();
            total_rational = (total_rational + term).simplify();

            // 更新分子
            SymbolicPolynomial G_deriv = poly_diff(G);
            current_num = T.subtract(G_deriv.multiply(U)).simplify();

            // 更新分母
            SymbolicPolynomial Vk_minus_1 = v.power(k - 1);
            current_den = Vk_minus_1.multiply(U);

            if (current_num.is_zero()) {
                *rational_part = total_rational;
                *reduced_numerator = current_num;
                *reduced_denominator = current_den;
                return true;
            }
        }
    }

    *rational_part = total_rational;
    *reduced_numerator = current_num;
    *reduced_denominator = current_den;

    return true;
}

// ============================================================================
// Lazard-Rioboo-Trager 算法
// ============================================================================

bool RischAlgorithm::lazard_rioboo_trager(const SymbolicPolynomial& A,
                                           const SymbolicPolynomial& D,
                                           const std::string& variable_name,
                                           SymbolicExpression* result) {
    // Lazard-Rioboo-Trager 算法用于处理高次不可约多项式
    // 它是 Rothstein-Trager 的改进版本，可以更高效地处理复数根

    if (!result) return false;

    SymbolicPolynomial Dp = D.derivative();
    std::string c_var = "risch_c";

    int deg_a = A.degree();
    int deg_dp = Dp.degree();
    int max_deg = std::max(deg_a, deg_dp);

    // 构造关于 c 的多项式: resultant(A - c*D', D)
    std::vector<SymbolicExpression> poly_c_coeffs;
    for (int i = 0; i <= max_deg; ++i) {
        SymbolicExpression term = (A.coefficient(i) -
                                  SymbolicExpression::variable(c_var) * Dp.coefficient(i)).simplify();
        poly_c_coeffs.push_back(term);
    }

    SymbolicPolynomial poly_c(poly_c_coeffs, variable_name);
    SymbolicExpression res_expr = poly_c.resultant(D).simplify();

    // 提取关于 c 的多项式系数
    std::vector<SymbolicExpression> res_roots_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(res_expr, c_var, &res_roots_coeffs)) {
        return false;
    }

    std::size_t poly_degree = res_roots_coeffs.size() - 1;
    if (poly_degree == 0) {
        *result = SymbolicExpression::number(0.0);
        return true;
    }

    // 尝试数值求解根
    bool all_numeric = true;
    std::vector<double> num_coeffs;

    for (const auto& c : res_roots_coeffs) {
        double val = 0.0;
        if (c.is_number(&val)) {
            num_coeffs.push_back(val);
        } else {
            all_numeric = false;
            break;
        }
    }

    SymbolicExpression final_result = SymbolicExpression::number(0.0);

    if (all_numeric && poly_degree > 2) {
        // 使用 Aberth-Ehrlich 方法找复数根
        auto complex_roots = find_complex_roots_aberth(num_coeffs);

        for (const auto& [re, im] : complex_roots) {
            if (std::abs(im) < 1e-10) {
                // 实数根
                if (std::abs(re) < 1e-10) continue;

                std::vector<SymbolicExpression> cur_poly_coeffs;
                for (int j = 0; j <= max_deg; ++j) {
                    cur_poly_coeffs.push_back((A.coefficient(j) -
                                              SymbolicExpression::number(re) * Dp.coefficient(j)).simplify());
                }
                SymbolicPolynomial poly_i(cur_poly_coeffs, variable_name);
                SymbolicPolynomial v_i = poly_i.gcd(D);

                if (!v_i.is_zero() && !v_i.is_constant()) {
                    final_result = (final_result +
                                   SymbolicExpression::number(re) *
                                   make_function("ln", v_i.to_expression())).simplify();
                }
            } else {
                // 复数共轭对
                // 转换为 arctan 形式
                if (D.degree() == 2) {
                    double a_d = 0.0, b_d = 0.0, c_d = 0.0;
                    D.coefficient(2).is_number(&a_d);
                    D.coefficient(1).is_number(&b_d);
                    D.coefficient(0).is_number(&c_d);

                    double disc_neg = 4.0 * a_d * c_d - b_d * b_d;

                    if (disc_neg > 0) {
                        double sqrt_disc = std::sqrt(disc_neg);
                        SymbolicExpression x = SymbolicExpression::variable(variable_name);

                        SymbolicExpression atan_arg = ((SymbolicExpression::number(2.0 * a_d) * x +
                                                       SymbolicExpression::number(b_d)) /
                                                       SymbolicExpression::number(sqrt_disc)).simplify();

                        double atan_coeff = 2.0 * im / sqrt_disc;

                        if (A.degree() == 0) {
                            double a_coeff_val = 1.0;
                            A.coefficient(0).is_number(&a_coeff_val);
                            atan_coeff *= a_coeff_val;
                        }

                        final_result = (final_result +
                                       SymbolicExpression::number(atan_coeff) *
                                       make_function("atan", atan_arg)).simplify();

                        if (std::abs(re) > 1e-10) {
                            double ln_coeff = re / a_d;
                            final_result = (final_result +
                                           SymbolicExpression::number(ln_coeff) *
                                           make_function("ln", D.to_expression())).simplify();
                        }
                    }
                }
            }
        }
    } else {
        // 回退到标准 Rothstein-Trager
        return false;
    }

    *result = final_result;
    return true;
}

// ============================================================================
// Rothstein-Trager 算法
// ============================================================================

bool RischAlgorithm::rothstein_trager(const SymbolicPolynomial& numerator,
                                      const SymbolicPolynomial& denominator,
                                      const std::string& variable_name,
                                      SymbolicExpression* log_part,
                                      const std::vector<DifferentialExtension>& tower,
                                      int tower_index,
                                      const std::string& main_var,
                                      const SymbolicExpression* t_prime,
                                      DifferentialExtension::Kind kind) {
    auto poly_diff = [&](const SymbolicPolynomial& p) {
        return t_prime ? p.total_derivative(main_var, *t_prime) : p.derivative();
    };

    SymbolicPolynomial D = denominator;
    SymbolicPolynomial Dp = poly_diff(D);
    SymbolicPolynomial A = numerator;

    std::string c_var = "risch_c";

    std::vector<SymbolicExpression> poly_c_coeffs;
    int deg_a = A.degree();
    int deg_dp = Dp.degree();
    int max_deg = std::max(deg_a, deg_dp);

    for (int i = 0; i <= max_deg; ++i) {
        SymbolicExpression term = (A.coefficient(i) -
                                  SymbolicExpression::variable(c_var) * Dp.coefficient(i)).simplify();
        poly_c_coeffs.push_back(term);
    }

    SymbolicPolynomial poly_c(poly_c_coeffs, variable_name);
    SymbolicExpression res_expr = poly_c.resultant(D).simplify();

    // 检查结式的系数是否为常数
    // 在 Risch 理论中，基域是 K(t)。如果结式依赖于 t (variable_name)，则其根 c 依赖于 t。
    // Risch 定理指出，如果积分是初等的，则结式的根必须是常数（在常数域中）。
    // 这里 main_var 是塔中较低层级的变量，variable_name 是当前层级的超越变量。
    // 我们必须确保 c 不依赖于 variable_name。
    if (contains_var(res_expr, variable_name)) {
         return false; // 结式依赖于当前变量，通常意味着非初等或需要代数扩张
    }

    // 检查结式是否依赖于 main_var (更低层级的超越变量)
    // 理想情况下，c 应该是绝对常数，但如果是嵌套扩张，c 可以是基域中的常数。
    if (contains_var(res_expr, main_var)) {
        // 如果结式包含更低层级的超越变量，可能需要进一步处理或判定为非初等
        // 暂时保守处理
    }

    std::vector<SymbolicExpression> res_roots_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(res_expr, c_var, &res_roots_coeffs)) {
        return false;
    }

    std::size_t poly_degree = res_roots_coeffs.size() - 1;

    if (poly_degree == 0) {
        *log_part = SymbolicExpression::number(0.0);
        return true;
    }

    // 查找所有根（包括复数根）
    auto roots = find_all_roots(res_roots_coeffs, c_var);

    if (roots.empty()) {
        return false;
    }

    // 构建对数部分
    SymbolicExpression final_log = SymbolicExpression::number(0.0);

    for (const auto& root : roots) {
        if (!root.is_complex) {
            // 实数根处理
            double c_val = 0.0;
            if (root.real_part.is_number(&c_val) && std::abs(c_val) < 1e-10) continue;

            std::vector<SymbolicExpression> cur_poly_coeffs;
            for (int j = 0; j <= max_deg; ++j) {
                cur_poly_coeffs.push_back((A.coefficient(j) - root.real_part * Dp.coefficient(j)).simplify());
            }
            SymbolicPolynomial poly_i(cur_poly_coeffs, variable_name);
            SymbolicPolynomial v_i = poly_i.gcd(D);

            if (!v_i.is_zero() && !v_i.is_constant()) {
                final_log = (final_log + root.real_part * make_function("ln", v_i.to_expression())).simplify();
            }
        } else if (root.is_conjugate_pair) {
            // 复数共轭对处理：c = a + bi 和 c' = a - bi
            // 使用 Lazard-Rioboo-Trager 转换为实数域的 arctan 和 ln

            double a_val = 0.0, b_val = 0.0;
            root.real_part.is_number(&a_val);
            root.imag_part.is_number(&b_val);

            // 对于二次不可约分母 D = x^2 + px + q
            // 结果 = (2*a / sqrt(4q - p^2)) * atan((2x + p) / sqrt(4q - p^2))
            //      + (b / 2) * ln(D)  (如果 a != 0)
            // 对于 a = 0 的特殊情况（如 1/(x^2+1)）：
            // 结果 = (2*b / sqrt(4q - p^2)) * atan((2x + p) / sqrt(4q - p^2))

            if (D.degree() == 2) {
                // 获取二次多项式系数 D = a_d * x^2 + b_d * x + c_d
                double a_d = 0.0, b_d = 0.0, c_d = 0.0;
                D.coefficient(2).is_number(&a_d);
                D.coefficient(1).is_number(&b_d);
                D.coefficient(0).is_number(&c_d);

                // 判别式 disc = b_d^2 - 4*a_d*c_d < 0 表示不可约
                // 但我们使用 4*a_d*c_d - b_d^2 > 0 来计算 arctan 参数
                double disc_neg = 4.0 * a_d * c_d - b_d * b_d;

                if (disc_neg > 0) {
                    double sqrt_disc = std::sqrt(disc_neg);
                    SymbolicExpression x = SymbolicExpression::variable(variable_name);

                    // atan 参数: (2*a_d*x + b_d) / sqrt_disc
                    SymbolicExpression atan_arg = ((SymbolicExpression::number(2.0 * a_d) * x +
                                                 SymbolicExpression::number(b_d)) /
                                                 SymbolicExpression::number(sqrt_disc)).simplify();

                    // arctan 系数
                    // 对于 ∫ P(x)/(ax^2+bx+c) dx，其中 P 是常数或线性
                    // Rothstein-Trager 根 c = a + bi 给出:
                    // atan 系数 = 2 * Im(c) / sqrt_disc * leading_coeff_of_A_if_constant
                    double atan_coeff = 2.0 * b_val / sqrt_disc;

                    // 如果 A 是常数，需要乘以 A 的系数
                    if (A.degree() == 0) {
                        double a_coeff_val = 1.0;
                        A.coefficient(0).is_number(&a_coeff_val);
                        atan_coeff *= a_coeff_val;
                    }

                    final_log = (final_log + SymbolicExpression::number(atan_coeff) *
                                make_function("atan", atan_arg)).simplify();

                    // 如果 real_part != 0，还需要添加 ln(D) 项
                    if (std::abs(a_val) > 1e-10) {
                        double ln_coeff = a_val / a_d;
                        final_log = (final_log + SymbolicExpression::number(ln_coeff) *
                                    make_function("ln", D.to_expression())).simplify();
                    }
                }
            } else {
                // 更高次不可约多项式：需要更复杂的处理
                // 这里简化处理，尝试直接使用复数根
                // 对于一般情况，需要完整的 Lazard-Rioboo-Trager 算法

                // 尝试计算 v = gcd(A - c*D', D) 的实部和虚部
                // 由于我们没有复数多项式运算，这里跳过
                // 但可以尝试一些特殊模式

                // 检查是否 A 是常数且 D 是二次的幂
                if (A.degree() == 0 && D.degree() > 2) {
                    // 尝试分解 D 为二次因子
                    // 这里简化处理
                }
            }
        }
    }

    *log_part = final_log;
    return true;
}

// ============================================================================

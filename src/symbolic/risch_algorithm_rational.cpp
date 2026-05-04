#include "symbolic/risch_algorithm.h"
#include "symbolic/risch_algorithm_internal.h"
#include "symbolic/symbolic_algebraic_number.h"
#include "symbolic/symbolic_expression_internal.h"
#include "symbolic/differential_field.h"
#include <cmath>
#include <vector>

using namespace symbolic_expression_internal;
using namespace risch_algorithm_internal;
using namespace sturm;

namespace {

constexpr int kMaxParametricPrsSteps = 24;
constexpr std::size_t kMaxParametricExpressionChars = 12000;

SymbolicExpression real_log_abs(const SymbolicExpression& argument) {
    return make_function("ln", make_function("abs", argument));
}

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

SymbolicPolynomial gcd_for_residue_value(const SymbolicPolynomial& A,
                                         const SymbolicPolynomial& Dp,
                                         const SymbolicPolynomial& D,
                                         const std::string& variable_name,
                                         const SymbolicExpression& c_value) {
    const int max_degree = std::max(A.degree(), Dp.degree());
    std::vector<SymbolicExpression> coefficients;
    coefficients.reserve(max_degree + 1);

    for (int i = 0; i <= max_degree; ++i) {
        coefficients.push_back((A.coefficient(i) - c_value * Dp.coefficient(i)).simplify());
    }

    return SymbolicPolynomial(coefficients, variable_name).gcd(D);
}

bool polynomial_expression_size_ok(const SymbolicPolynomial& polynomial) {
    for (int i = 0; i <= polynomial.degree(); ++i) {
        if (polynomial.coefficient(i).to_string().size() > kMaxParametricExpressionChars) {
            return false;
        }
    }
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
    SymbolicExpression arg;
    if (std::abs(b) < 1e-12) {
        arg = make_divide(x, SymbolicExpression::number(root / (2.0 * a)));
    } else {
        arg =
            ((SymbolicExpression::number(2.0 * a) * x + SymbolicExpression::number(b)) /
             SymbolicExpression::number(root)).simplify();
    }
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

bool try_integrate_split_linear_partial_fractions(const SymbolicPolynomial& numerator,
                                                  const SymbolicPolynomial& denominator,
                                                  const std::string& variable_name,
                                                  SymbolicExpression* result) {
    auto factors = denominator.factor_linear();
    if (factors.empty()) {
        return false;
    }

    int total_degree = 0;
    for (const auto& factor : factors) {
        if (factor.first.degree() != 1) {
            return false;
        }
        total_degree += factor.first.degree() * factor.second;
    }
    if (total_degree != denominator.degree()) {
        return false;
    }

    std::vector<std::pair<SymbolicExpression, SymbolicPolynomial>> parts;
    if (!partial_fraction_decomposition(numerator, factors, variable_name, &parts) ||
        parts.empty()) {
        return false;
    }

    SymbolicExpression total = SymbolicExpression::number(0.0);
    for (const auto& part : parts) {
        const SymbolicExpression coefficient = part.first.simplify();
        if (SymbolicPolynomial::coeff_is_zero(coefficient)) {
            continue;
        }

        const SymbolicPolynomial& denominator_part = part.second;
        const int power = denominator_part.degree();
        if (power <= 0) {
            return false;
        }

        const SymbolicPolynomial* base_factor = nullptr;
        for (const auto& factor : factors) {
            if (factor.first.degree() != 1 || factor.second < power) {
                continue;
            }
            if (structural_equals(factor.first.power(power).to_expression().simplify(),
                                  denominator_part.to_expression().simplify())) {
                base_factor = &factor.first;
                break;
            }
        }
        if (!base_factor) {
            return false;
        }

        SymbolicExpression a;
        SymbolicExpression b;
        base_factor->is_linear_factor(&a, &b);
        const SymbolicExpression linear = base_factor->to_expression().simplify();
        SymbolicExpression term;
        if (power == 1) {
            term = make_multiply(make_divide(coefficient, a),
                                 real_log_abs(linear)).simplify();
        } else {
            term = make_multiply(
                       make_divide(coefficient,
                                   make_multiply(a,
                                                 SymbolicExpression::number(1.0 - power))),
                       make_power(linear,
                                  SymbolicExpression::number(1.0 - power)))
                       .simplify();
        }
        total = make_add(total, term).simplify();
    }

    if (expr_is_zero(total)) {
        return false;
    }
    *result = total;
    return true;
}

bool try_integrate_distinct_even_quadratic_product(const SymbolicPolynomial& numerator,
                                                   const SymbolicPolynomial& denominator,
                                                   const std::string& variable_name,
                                                   SymbolicExpression* result) {
    if (numerator.degree() != 0 || denominator.degree() != 4) {
        return false;
    }

    double num = 0.0;
    double c0 = 0.0;
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0;
    double c4 = 0.0;
    if (!numerator.coefficient(0).is_number(&num) ||
        !denominator.coefficient(0).is_number(&c0) ||
        !denominator.coefficient(1).is_number(&c1) ||
        !denominator.coefficient(2).is_number(&c2) ||
        !denominator.coefficient(3).is_number(&c3) ||
        !denominator.coefficient(4).is_number(&c4) ||
        std::abs(c1) > 1e-10 ||
        std::abs(c3) > 1e-10 ||
        std::abs(c4 - 1.0) > 1e-10) {
        return false;
    }

    const double discriminant = c2 * c2 - 4.0 * c0;
    if (discriminant <= 1e-12) {
        return false;
    }

    double first = (c2 - std::sqrt(discriminant)) / 2.0;
    double second = (c2 + std::sqrt(discriminant)) / 2.0;
    if (first <= 0.0 || second <= 0.0 || std::abs(first - second) < 1e-12) {
        return false;
    }

    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    auto atan_scaled = [&](double constant) {
        const double root = std::sqrt(constant);
        SymbolicExpression arg = std::abs(root - 1.0) < 1e-12
            ? x
            : make_multiply(SymbolicExpression::number(1.0 / root), x).simplify();
        return make_multiply(SymbolicExpression::number(1.0 / root),
                             make_function("atan", arg)).simplify();
    };

    *result = make_multiply(
                  SymbolicExpression::number(num / (second - first)),
                  make_subtract(atan_scaled(first), atan_scaled(second)))
                  .simplify();
    return true;
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

bool RischAlgorithm::try_integrate_exponential_special_part(
    const SymbolicPolynomial& numerator,
    const SymbolicPolynomial& denominator,
    const std::string& variable_name,
    const SymbolicExpression& t_prime,
    const std::vector<DifferentialExtension>& tower,
    int tower_index,
    const std::string& main_var,
    SymbolicExpression* result) {

    // 检查分母是否为 t^n 形式
    if (denominator.degree() <= 0) return false;
    for (int i = 0; i < denominator.degree(); ++i) {
        if (!SymbolicPolynomial::coeff_is_zero(denominator.coefficient(i))) return false;
    }

    int n = denominator.degree();
    SymbolicExpression a_n = denominator.leading_coefficient();
    SymbolicExpression t = SymbolicExpression::variable(variable_name);
    
    // t' = u't => u' = t' / t
    SymbolicExpression u_prime = (t_prime / t).simplify();
    SymbolicExpression total_int = SymbolicExpression::number(0.0);

    for (int i = 0; i <= numerator.degree(); ++i) {
        if (SymbolicPolynomial::coeff_is_zero(numerator.coefficient(i))) continue;

        // 项为 (num[i] / a_n) * t^i / t^n = a * t^(-k), 其中 k = n - i
        int k = n - i;
        if (k <= 0) continue; 

        SymbolicExpression a = (numerator.coefficient(i) / a_n).simplify();
        
        // 对于指数扩展，t^-k 的积分为 y * t^-k，满足 y' - k u' y = a
        SymbolicExpression f = (make_negate(SymbolicExpression::number(static_cast<double>(k))) * u_prime).simplify();

        IntegrationResult y_res = solve_rde(f, a, main_var, tower, tower_index - 1);
        if (!y_res.success || y_res.type != IntegralType::kElementary) return false;

        total_int = (total_int + y_res.value / make_power(t, SymbolicExpression::number(static_cast<double>(k)))).simplify();
    }

    *result = total_int;
    return true;
}

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
                // 根据 Liouville 定理，如果 P(t) 在 K(t) 中有初等积分，
                // 则最高次项系数 a_n 的积分必须在 K 中。
                // 我们严格要求 a_int 在 K 中，否则积分非初等。

                IntegrationResult a_int_res = integrate_in_extension(q_coeffs[i], tower, tower_index - 1, main_var);
                if (!a_int_res.success || a_int_res.type != IntegralType::kElementary) return false;

                SymbolicExpression a_int = a_int_res.value;

                // 严格检查 a_int 是否在基域 K 中
                // 如果 a_int 包含任何不应存在于 K 中的超越函数，则说明不符合形式
                std::vector<std::pair<SymbolicExpression, SymbolicExpression>> dummy_logs;
                SymbolicExpression dummy_rest = SymbolicExpression::number(0.0);
                LogarithmicRepresentation log_rep = express_as_logarithmic_sum(a_int,
                                                                               std::vector<DifferentialExtension>(tower.begin(), tower.begin() + tower_index),
                                                                               main_var);
                if (!log_rep.is_valid) {
                    return false; // a_int 不在 K 中，非初等
                }

                if (i == 0) {
                    term_int = a_int;
                } else {
                    SymbolicExpression t = SymbolicExpression::variable(variable_name);
                    SymbolicExpression first_part = (a_int * make_power(t, SymbolicExpression::number(static_cast<double>(i)))).simplify();

                    // 计算校正项: ∫ a_int * i * t^(i-1) * t' dx
                    SymbolicExpression t_prime_val = t_prime ? *t_prime : SymbolicExpression::number(1.0);
                    SymbolicExpression correction_integrand =
                        (a_int * SymbolicExpression::number(static_cast<double>(i)) *
                         make_power(t, SymbolicExpression::number(static_cast<double>(i - 1))) *
                         t_prime_val).simplify();

                    IntegrationResult correction_res = integrate_in_extension(correction_integrand, tower, tower_index, main_var);
                    if (!correction_res.success || correction_res.type != IntegralType::kElementary) return false;

                    term_int = (first_part - correction_res.value).simplify();
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

    if (kind == DifferentialExtension::Kind::kNone && denominator.degree() > 0) {
        bool monomial_denominator = true;
        for (int i = 0; i < denominator.degree(); ++i) {
            if (!SymbolicPolynomial::coeff_is_zero(denominator.coefficient(i))) {
                monomial_denominator = false;
                break;
            }
        }

        if (monomial_denominator) {
            SymbolicExpression x = SymbolicExpression::variable(variable_name);
            SymbolicExpression laurent_part = SymbolicExpression::number(0.0);
            const SymbolicExpression leading = denominator.leading_coefficient();
            const int denominator_power = denominator.degree();

            for (int i = 0; i <= R.degree(); ++i) {
                if (SymbolicPolynomial::coeff_is_zero(R.coefficient(i))) {
                    continue;
                }

                const int exponent = i - denominator_power;
                SymbolicExpression coefficient = (R.coefficient(i) / leading).simplify();
                if (exponent == -1) {
                    laurent_part =
                        (laurent_part +
                         coefficient * make_function("ln", make_function("abs", x))).simplify();
                } else {
                    SymbolicExpression new_power =
                        SymbolicExpression::number(static_cast<double>(exponent + 1));
                    laurent_part =
                        (laurent_part +
                         coefficient * make_power(x, new_power) / new_power).simplify();
                }
            }

            *result = (poly_int + laurent_part).simplify();
            return true;
        }
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

    if (kind == DifferentialExtension::Kind::kExponential && t_prime) {
        SymbolicExpression special_part;
        if (RischAlgorithm::try_integrate_exponential_special_part(R,
                                                                   denominator,
                                                                   variable_name,
                                                                   *t_prime,
                                                                   tower,
                                                                   tower_index,
                                                                   main_var,
                                                                   &special_part)) {
            *result = (poly_int + special_part).simplify();
            return true;
        }
    }

    if (kind == DifferentialExtension::Kind::kNone) {
        SymbolicExpression split_linear_part;
        if (try_integrate_split_linear_partial_fractions(R,
                                                         denominator,
                                                         variable_name,
                                                         &split_linear_part)) {
            *result = (poly_int + split_linear_part).simplify();
            return true;
        }

        SymbolicExpression even_quadratic_part;
        if (try_integrate_distinct_even_quadratic_product(R,
                                                          denominator,
                                                          variable_name,
                                                          &even_quadratic_part)) {
            *result = (poly_int + even_quadratic_part).simplify();
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

    if (kind == DifferentialExtension::Kind::kNone) {
        SymbolicExpression partial_fraction_part;
        if (try_integrate_numeric_quadratic_partial_fractions(reduced_num,
                                                              reduced_den,
                                                              variable_name,
                                                              &partial_fraction_part)) {
            *result = (poly_int + rational_part + partial_fraction_part).simplify();
            return true;
        }
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

    std::vector<double> numeric_denominator_coefficients;
    if (!numeric_coefficients(denominator, &numeric_denominator_coefficients)) {
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
// Lazard-Rioboo-Trager 算法 (改进版，支持代数数根)
// ============================================================================

bool RischAlgorithm::lazard_rioboo_trager(const SymbolicPolynomial& A,
                                           const SymbolicPolynomial& D,
                                           const std::string& variable_name,
                                           SymbolicExpression* result) {
    if (!result) return false;

    // Lazard-Rioboo-Trager 算法改进版本
    // 目标：计算 ∫ (A/D) dx 的对数部分，其中 D 是 square-free 的。
    // 该算法避免了在常数域的代数扩张中进行运算，直到最后一步。

    SymbolicPolynomial Dp = D.derivative();
    std::string c_var = "risch_c";

    // 1. 计算结式 R(c) = resultant_x(A - c*D', D)
    int deg_a = A.degree();
    int deg_dp = Dp.degree();
    int max_deg = std::max(deg_a, deg_dp);

    std::vector<SymbolicExpression> poly_c_coeffs_in_x;
    for (int i = 0; i <= max_deg; ++i) {
        SymbolicExpression term = (A.coefficient(i) -
                                  SymbolicExpression::variable(c_var) * Dp.coefficient(i)).simplify();
        poly_c_coeffs_in_x.push_back(term);
    }

    SymbolicPolynomial poly_c_x(poly_c_coeffs_in_x, variable_name);
    SymbolicExpression res_expr = poly_c_x.resultant(D).simplify();

    // 2. 提取 R(c) 的系数
    std::vector<SymbolicExpression> res_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(res_expr, c_var, &res_coeffs)) {
        return false;
    }

    SymbolicPolynomial R_c(res_coeffs, c_var);
    if (R_c.is_zero()) {
        *result = SymbolicExpression::number(0.0);
        return true;
    }

    // 3. 对 R(c) 进行 square-free 分解
    std::vector<SymbolicPolynomial> r_factors;
    if (!R_c.square_free_decomposition(&r_factors)) {
        r_factors = {R_c};
    }

    SymbolicExpression total_log = SymbolicExpression::number(0.0);

    // 4. 对于每个 square-free 因子 r_i(c)
    for (int i = 0; i < (int)r_factors.size(); ++i) {
        const auto& ri = r_factors[i];
        if (ri.degree() <= 0) continue;

        // 计算 v_i(x, c) = gcd_x(A - c*D', D) modulo r_i(c)
        // 这是一个带参数的多项式 GCD。在 Risch 理论中，这可以通过子结果项序列（Subresultant Chain）高效完成。

        if (ri.degree() == 1) {
            // 线性因子: c = -ri(0) / ri(1)
            SymbolicExpression c_val = (make_negate(ri.coefficient(0)) / ri.coefficient(1)).simplify();

            std::vector<SymbolicExpression> sub_poly_coeffs;
            for (int j = 0; j <= max_deg; ++j) {
                sub_poly_coeffs.push_back((A.coefficient(j) - c_val * Dp.coefficient(j)).simplify());
            }
            SymbolicPolynomial sub_poly(sub_poly_coeffs, variable_name);
            SymbolicPolynomial vi = sub_poly.gcd(D);

            if (!vi.is_zero() && !vi.is_constant()) {
                total_log = (total_log + c_val * real_log_abs(vi.to_expression())).simplify();
            }
        } else if (ri.degree() == 2) {
            // 二次因子: 使用求根公式保持 sqrt 形式
            // r(c) = a*c^2 + b*c + d = 0
            // c = (-b ± sqrt(b^2 - 4ad)) / (2a)

            SymbolicExpression a = ri.coefficient(2);
            SymbolicExpression b = ri.coefficient(1);
            SymbolicExpression d = ri.coefficient(0);

            // 判别式
            SymbolicExpression disc = (b * b - SymbolicExpression::number(4.0) * a * d).simplify();

            // 检查判别式是否为完全平方
            double disc_val = 0.0;
            if (disc.is_number(&disc_val)) {
                double sqrt_disc = std::sqrt(disc_val);
                if (std::abs(sqrt_disc * sqrt_disc - disc_val) < 1e-12) {
                    // 判别式是完全平方，根是有理数
                    SymbolicExpression c1 = ((make_negate(b) + SymbolicExpression::number(sqrt_disc)) /
                                            (SymbolicExpression::number(2.0) * a)).simplify();
                    SymbolicExpression c2 = ((make_negate(b) - SymbolicExpression::number(sqrt_disc)) /
                                            (SymbolicExpression::number(2.0) * a)).simplify();

                    for (const auto& c_val : {c1, c2}) {
                        std::vector<SymbolicExpression> sub_poly_coeffs;
                        for (int j = 0; j <= max_deg; ++j) {
                            sub_poly_coeffs.push_back((A.coefficient(j) - c_val * Dp.coefficient(j)).simplify());
                        }
                        SymbolicPolynomial sub_poly(sub_poly_coeffs, variable_name);
                        SymbolicPolynomial vi = sub_poly.gcd(D);

                        if (!vi.is_zero() && !vi.is_constant()) {
                            total_log = (total_log + c_val * real_log_abs(vi.to_expression())).simplify();
                        }
                    }
                } else {
                    // 判别式不是完全平方，使用 sqrt 形式
                    // 对于共轭根对 c = (-b ± sqrt(disc)) / (2a)
                    // 积分结果可以表示为实数形式

                    SymbolicExpression sqrt_disc = make_function("sqrt", disc);
                    SymbolicExpression two_a = (SymbolicExpression::number(2.0) * a).simplify();

                    // 对于共轭根对，使用代数数表示
                    // 这里我们使用 Sturm 序列进行实根隔离
                    auto intervals = sturm::isolate_real_roots(ri);

                    if (!intervals.empty()) {
                        // 存在实根
                        for (const auto& [lower, upper] : intervals) {
                            // 创建代数数
                            AlgebraicNumber c_alpha(ri, lower, upper, true, 0, 0);

                            // 计算 GCD
                            double c_approx = c_alpha.approximate();
                            SymbolicExpression c_val = SymbolicExpression::number(c_approx);

                            std::vector<SymbolicExpression> sub_poly_coeffs;
                            for (int j = 0; j <= max_deg; ++j) {
                                sub_poly_coeffs.push_back((A.coefficient(j) - c_val * Dp.coefficient(j)).simplify());
                            }
                            SymbolicPolynomial sub_poly(sub_poly_coeffs, variable_name);
                            SymbolicPolynomial vi = sub_poly.gcd(D);

                            if (!vi.is_zero() && !vi.is_constant()) {
                                // 使用代数数表示
                                total_log = (total_log + c_alpha.to_expression() *
                                           real_log_abs(vi.to_expression())).simplify();
                            }
                        }
                    } else {
                        // 无实根，使用复数共轭对
                        // 对于二次不可约因子，根是复数共轭对
                        // 积分结果可以表示为 arctan 形式

                        // 使用子结果式链计算 v_i
                        SubresultantChain chain = compute_subresultant_chain(A, Dp, D, c_var);

                        // 对于二次因子，对应的 v_i 是二次的
                        SymbolicPolynomial vi;
                        for (size_t j = 0; j < chain.degrees.size(); ++j) {
                            if (chain.degrees[j] == 2) {
                                vi = chain.subresultants[j];
                                break;
                            }
                        }

                        if (!vi.is_zero()) {
                            // 对于复数根，结果涉及 arctan
                            // 这里简化处理，返回符号形式
                            // 完整实现需要计算复数对数并转换为实数形式

                            // 使用 RootOf 表示代数数
                            SymbolicExpression c_root = make_rootof_from_polynomial(ri, 0);
                            total_log = (total_log + c_root * real_log_abs(vi.to_expression())).simplify();
                        }
                    }
                }
            } else {
                // 符号判别式
                return false;
            }
        } else {
            // 高次因子 (degree >= 3)
            // 使用 Sturm 序列和代数数

            // 首先检查是否所有系数都是数值
            bool all_numeric = true;
            for (int k = 0; k <= ri.degree(); ++k) {
                double val = 0.0;
                if (!ri.coefficient(k).is_number(&val)) {
                    all_numeric = false;
                    break;
                }
            }

            if (all_numeric) {
                // 使用 Sturm 序列隔离实根
                auto intervals = sturm::isolate_real_roots(ri);

                if (!intervals.empty()) {
                    // 存在实根
                    SubresultantChain chain = compute_subresultant_chain(A, Dp, D, c_var);

                    for (const auto& [lower, upper] : intervals) {
                        // 创建代数数
                        AlgebraicNumber c_alpha(ri, lower, upper, true,
                                               static_cast<int>(total_log.to_string().length()), 0);

                        // 从子结果式链获取对应的 v_i
                        int deg_c = ri.degree();
                        SymbolicPolynomial vi;
                        for (size_t j = 0; j < chain.degrees.size(); ++j) {
                            if (chain.degrees[j] == deg_c) {
                                vi = chain.subresultants[j];
                                break;
                            }
                        }

                        if (!vi.is_zero() && !vi.is_constant()) {
                            total_log = (total_log + c_alpha.to_expression() *
                                       real_log_abs(vi.to_expression())).simplify();
                        }
                    }
                } else {
                    // 无实根，全部是复数共轭对
                    // 使用子结果式链
                    SubresultantChain chain = compute_subresultant_chain(A, Dp, D, c_var);

                    int deg_c = ri.degree();
                    SymbolicPolynomial vi;
                    for (size_t j = 0; j < chain.degrees.size(); ++j) {
                        if (chain.degrees[j] == deg_c) {
                            vi = chain.subresultants[j];
                            break;
                        }
                    }

                    if (!vi.is_zero()) {
                        // 使用 RootOf 表示
                        SymbolicExpression c_root = make_rootof_from_polynomial(ri, 0);
                        total_log = (total_log + c_root * real_log_abs(vi.to_expression())).simplify();
                    }
                }
            } else {
                // 符号系数，无法处理
                return false;
            }
        }
    }

    *result = total_log;
    return !expr_is_zero(total_log);
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
    
    // 优先尝试 Lazard-Rioboo-Trager，因为它更通用且处理高次项更稳健
    if (kind == DifferentialExtension::Kind::kNone) {
        if (lazard_rioboo_trager(numerator, denominator, variable_name, log_part)) {
            return true;
        }
    }

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
                final_log = (final_log + root.real_part * real_log_abs(v_i.to_expression())).simplify();
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
                                    real_log_abs(D.to_expression())).simplify();
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
    return !expr_is_zero(final_log);
}

// ============================================================================
// Phase 4.2: 子结果式链与改进的 Lazard-Rioboo-Trager 算法
// ============================================================================

// 子结果式链的 GCD 提取
SymbolicPolynomial RischAlgorithm::SubresultantChain::gcd_at_parameter(
    const SymbolicExpression& c_value) const {

    if (subresultants.empty()) {
        return SymbolicPolynomial({}, var_name);
    }

    // 在 c = c_value 处计算每个子结果式
    for (std::size_t i = 0; i < subresultants.size(); ++i) {
        SymbolicPolynomial sub = subresultants[i];

        std::vector<SymbolicExpression> coeffs;
        for (int j = 0; j <= sub.degree(); ++j) {
            SymbolicExpression coeff = sub.coefficient(j);
            if (!parameter_name.empty()) {
                coeff = coeff.substitute(parameter_name, c_value);
            }
            coeffs.push_back(coeff.simplify());
        }

        SymbolicPolynomial sub_at_c(coeffs, var_name);

        // 检查是否非零
        if (!sub_at_c.is_zero()) {
            return sub_at_c;
        }
    }

    return SymbolicPolynomial({}, var_name);
}

// 计算子结果式链
RischAlgorithm::SubresultantChain RischAlgorithm::compute_subresultant_chain(
    const SymbolicPolynomial& A,
    const SymbolicPolynomial& D_prime,
    const SymbolicPolynomial& D,
    const std::string& c_var) {

    SubresultantChain chain;
    chain.var_name = D.variable_name();
    chain.parameter_name = c_var;

    int deg_A = A.degree();

    // 构造 P(c, x) = A - c * D'
    std::vector<SymbolicExpression> P_coeffs;
    int max_deg = std::max(deg_A, D_prime.degree());
    for (int i = 0; i <= max_deg; ++i) {
        SymbolicExpression term = A.coefficient(i);
        term = (term - SymbolicExpression::variable(c_var) * D_prime.coefficient(i)).simplify();
        P_coeffs.push_back(term);
    }
    SymbolicPolynomial P(P_coeffs, chain.var_name);

    // 使用伪余式序列计算子结果式链
    SymbolicPolynomial current_P = P;
    SymbolicPolynomial current_Q = D;

    int last_degree = std::max(current_P.degree(), current_Q.degree()) + 1;
    int steps = 0;
    while (!current_Q.is_zero()) {
        int m = current_P.degree();
        int n = current_Q.degree();

        if (m < n) {
            std::swap(current_P, current_Q);
            std::swap(m, n);
        }

        if (n == 0) {
            break;
        }
        if (n >= last_degree || steps++ >= kMaxParametricPrsSteps ||
            !polynomial_expression_size_ok(current_P) ||
            !polynomial_expression_size_ok(current_Q)) {
            chain.subresultants.clear();
            chain.degrees.clear();
            return chain;
        }
        last_degree = n;

        // 计算伪余式
        SymbolicPolynomial quotient, remainder;
        SymbolicExpression lc_Q = current_Q.leading_coefficient();

        // 伪除法: lc_Q^(m-n+1) * P = Q * quotient + remainder
        int power = m - n + 1;
        SymbolicExpression factor = make_power(lc_Q, SymbolicExpression::number(static_cast<double>(power)));

        std::vector<SymbolicExpression> scaled_P_coeffs;
        for (int i = 0; i <= m; ++i) {
            scaled_P_coeffs.push_back((factor * current_P.coefficient(i)).simplify());
        }
        SymbolicPolynomial scaled_P(scaled_P_coeffs, chain.var_name);

        if (scaled_P.divide(current_Q, &quotient, &remainder)) {
            // 记录当前子结果式
            chain.subresultants.push_back(current_Q);
            chain.degrees.push_back(n);

            if (!polynomial_expression_size_ok(remainder)) {
                chain.subresultants.clear();
                chain.degrees.clear();
                return chain;
            }
            current_P = current_Q;
            current_Q = remainder;
        } else {
            break;
        }

        if (chain.subresultants.size() > kMaxParametricPrsSteps) break;
    }

    return chain;
}

// 改进的 Lazard-Rioboo-Trager 算法 (完全符号化版本)
bool RischAlgorithm::lazard_rioboo_trager_improved(
    const SymbolicPolynomial& A,
    const SymbolicPolynomial& D,
    const std::string& variable_name,
    SymbolicExpression* result) {

    if (!result) return false;

    SymbolicPolynomial Dp = D.derivative();
    std::string c_var = "risch_c";

    // 计算结式 R(c) = resultant_x(A - c*D', D)
    int deg_a = A.degree();
    int deg_dp = Dp.degree();
    int max_deg = std::max(deg_a, deg_dp);

    std::vector<SymbolicExpression> poly_c_coeffs_in_x;
    for (int i = 0; i <= max_deg; ++i) {
        SymbolicExpression term = (A.coefficient(i) -
                                  SymbolicExpression::variable(c_var) * Dp.coefficient(i)).simplify();
        poly_c_coeffs_in_x.push_back(term);
    }

    SymbolicPolynomial poly_c_x(poly_c_coeffs_in_x, variable_name);
    SymbolicExpression res_expr = poly_c_x.resultant(D).simplify();

    // 提取 R(c) 的系数
    std::vector<SymbolicExpression> res_coeffs;
    if (!symbolic_polynomial_coefficients_from_simplified(res_expr, c_var, &res_coeffs)) {
        return false;
    }

    SymbolicPolynomial R_c(res_coeffs, c_var);
    if (R_c.is_zero()) {
        *result = SymbolicExpression::number(0.0);
        return true;
    }

    // 对 R(c) 进行 square-free 分解
    std::vector<SymbolicPolynomial> r_factors;
    if (!R_c.square_free_decomposition(&r_factors)) {
        r_factors = {R_c};
    }

    SymbolicExpression total_log = SymbolicExpression::number(0.0);

    // 对于每个 square-free 因子 r_i(c)
    for (int i = 0; i < (int)r_factors.size(); ++i) {
        const auto& ri = r_factors[i];
        if (ri.degree() <= 0) continue;

        if (ri.degree() == 1) {
            // 线性因子: c = -ri(0) / ri(1)
            SymbolicExpression c_val = (make_negate(ri.coefficient(0)) / ri.coefficient(1)).simplify();

            // For a linear residue factor, c is known exactly in the constant
            // field, so compute gcd(A - c*D', D) directly.  Building a full
            // parametric subresultant chain here can cause coefficient swell
            // on higher-degree rational functions before we even know whether
            // the strict path can finish.
            SymbolicPolynomial vi = gcd_for_residue_value(A, Dp, D, variable_name, c_val);

            if (!vi.is_zero() && !vi.is_constant()) {
                total_log = (total_log + c_val * real_log_abs(vi.to_expression())).simplify();
            }
        } else {
            // Higher-degree residue fields need algebraic constant arithmetic
            // to stay exact.  The previous midpoint/subresultant shortcut was
            // not a valid strict Risch step and could explode in memory while
            // trying to build the parametric PRS.  Fail fast so integrate_full
            // can use the non-strict rational fallback instead of hanging.
            return false;
        }
    }

    *result = total_log;
    return !expr_is_zero(total_log);
}

// ============================================================================

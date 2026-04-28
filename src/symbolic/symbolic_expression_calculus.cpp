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

bool is_general_quadratic(const SymbolicExpression& expression,
                          const std::string& variable_name,
                          SymbolicExpression* a,
                          SymbolicExpression* b,
                          SymbolicExpression* c) {
    std::vector<SymbolicExpression> coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(expression, variable_name, &coeffs) &&
        coeffs.size() <= 3) {
        *c = (coeffs.size() > 0) ? coeffs[0] : SymbolicExpression::number(0.0);
        *b = (coeffs.size() > 1) ? coeffs[1] : SymbolicExpression::number(0.0);
        *a = (coeffs.size() > 2) ? coeffs[2] : SymbolicExpression::number(0.0);
        return true;
    }
    return false;
}

bool is_pure_quadratic(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       SymbolicExpression* constant_term,
                       SymbolicExpression* x2_coeff) {
    SymbolicExpression a, b, c;
    if (is_general_quadratic(expression, variable_name, &a, &b, &c) &&
        expr_is_zero(b) && !expr_is_zero(a)) {
        *constant_term = c;
        *x2_coeff = a;
        return true;
    }
    return false;
}

[[maybe_unused]] bool is_sqrt_one_minus_variable_squared(const SymbolicExpression& expression,
                                        const std::string& variable_name) {
    if (expression.node_->type != NodeType::kFunction ||
        expression.node_->text != "sqrt") {
        return false;
    }
    const SymbolicExpression inner(expression.node_->left);
    SymbolicExpression a, b, c;
    return is_general_quadratic(inner, variable_name, &a, &b, &c) &&
           expr_is_one(c) && expr_is_zero(b) && expr_is_minus_one(a);
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

SymbolicExpression clean_symbolic_constant(double value) {
    if (mymath::is_integer(value, 1e-7)) {
        return SymbolicExpression::number(value >= 0.0 ? static_cast<long long>(value + 0.5)
                                                      : static_cast<long long>(value - 0.5));
    }

    // Try common radicals: sqrt(2), sqrt(3), and their reciprocal/half multiples
    const double sqrt2 = 1.4142135623730951;
    const double sqrt3 = 1.7320508075688772;
    const std::vector<std::pair<double, std::string>> candidates = {
        {sqrt2, "sqrt(2)"}, {sqrt3, "sqrt(3)"},
        {1.0/sqrt2, "1/sqrt(2)"}, {1.0/sqrt3, "1/sqrt(3)"},
        {sqrt2/2.0, "sqrt(2)/2"}, {sqrt3/2.0, "sqrt(3)/2"},
        {2.0/sqrt3, "2/sqrt(3)"}
    };

    for (const auto& cand : candidates) {
        double scale = value / cand.first;
        if (mymath::is_integer(scale, 1e-6)) {
            int int_scale = static_cast<int>(std::round(scale));
            if (int_scale == 1) return SymbolicExpression::parse(cand.second);
            if (int_scale == -1) return make_negate(SymbolicExpression::parse(cand.second));
            return make_multiply(SymbolicExpression::number(int_scale), SymbolicExpression::parse(cand.second)).simplify();
        }
        // Try common fractions of radicals
        for (int den = 2; den <= 6; ++den) {
            for (int num = 1; num <= 10; ++num) {
                double f = static_cast<double>(num) / den;
                if (mymath::abs(value - f * cand.first) < 1e-7) {
                    return make_multiply(make_divide(SymbolicExpression::number(num), SymbolicExpression::number(den)),
                                         SymbolicExpression::parse(cand.second)).simplify();
                }
                if (mymath::abs(value + f * cand.first) < 1e-7) {
                    return make_negate(make_multiply(make_divide(SymbolicExpression::number(num), SymbolicExpression::number(den)),
                                                     SymbolicExpression::parse(cand.second))).simplify();
                }
            }
        }
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
        return make_divide(SymbolicExpression::number(numerator), SymbolicExpression::number(denominator)).simplify();
    }
    return SymbolicExpression::number(value);
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

/*
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
*/

struct LinearFactorMultiplicity {
    double root = 0.0;
    int multiplicity = 0;
};

struct QuadraticFactorMultiplicity {
    std::vector<double> coefficients;
    int multiplicity = 0;
};

bool polynomial_close(const std::vector<double>& lhs,
                      const std::vector<double>& rhs,
                      double eps = 1e-8);

bool divide_exact_polynomial(std::vector<double>* remaining,
                             const std::vector<double>& divisor);

SymbolicExpression integrate_quadratic_partial_fraction_term(
    const std::vector<double>& quadratic,
    double slope,
    double constant,
    int power,
    const std::string& variable_name);

bool extract_general_rational_factorization(
    const std::vector<double>& denominator,
    std::vector<LinearFactorMultiplicity>* linear_factors,
    std::vector<QuadraticFactorMultiplicity>* quadratic_factors) {
    linear_factors->clear();
    quadratic_factors->clear();

    std::vector<double> remaining = denominator;
    trim_coefficients(&remaining);
    
    // 1. 提取线性因子
    std::vector<double> roots;
    try {
        roots = polynomial_real_roots(remaining);
    } catch (...) {
        roots.clear();
    }

    for (double root : roots) {
        SymbolicExpression cleaned_root_expr = clean_symbolic_constant(root);
        double cleaned_root = root;
        expr_is_number(cleaned_root_expr, &cleaned_root);
        const std::vector<double> factor = {-cleaned_root, 1.0};
        int multiplicity = 0;
        while (remaining.size() > 1 && divide_exact_polynomial(&remaining, factor)) {
            ++multiplicity;
        }
        if (multiplicity > 0) {
            linear_factors->push_back(LinearFactorMultiplicity{root, multiplicity});
        }
    }

    trim_coefficients(&remaining);
    if (remaining.size() == 1) return true;

    // 2. 提取二次因子
    // 策略：对于剩余多项式，如果它不含实根，则它由不可约二次因子组成（代数基本定理的推论）
    // 目前我们支持检测形如 (ax^2 + bx + c)^n 的重复因子，或者不同的二次因子。
    // 由于缺乏复数求根器，我们尝试对一些常见的二次因子进行模糊匹配或启发式搜索。
    
    // 简易策略：对于 (x^2 + k)^n 这种常见情况
    if (remaining.size() >= 3 && (remaining.size() - 1) % 2 == 0) {
        // 尝试检测重复的二次因子
        const int degree = static_cast<int>(remaining.size() - 1);
        for (int m = degree / 2; m >= 1; --m) {
            if (degree % m == 0) {
                const int q_deg = 2;
                if (m * q_deg == degree) {
                    // 假设形如 (a x^2 + b x + c)^m
                    double a = mymath::pow(mymath::abs(remaining.back()), 1.0 / m);
                    if (remaining.back() < 0 && m % 2 != 0) a = -a;
                    double c = mymath::pow(mymath::abs(remaining[0]), 1.0 / m);
                    // 尝试不同的 c 符号 (对于不可约，c/a 必须 > (b/2a)^2)
                    for (double c_sign : {1.0, -1.0}) {
                        double trial_c = c * c_sign;
                        // 估算 b：通过一次项或最高次项次一项
                        double b = 0.0;
                        if (degree > 2) {
                            // (ax^2 + bx + c)^m = (ax^2)^m + m(ax^2)^{m-1}(bx) + ...
                            // coeff[deg-1] = m * a^{m-1} * b
                            b = remaining[degree - 1] / (m * mymath::pow(a, m - 1));
                        } else {
                            b = remaining[1] / m; // 实际上 m=1 时就是 remaining[1]
                        }
                        
                        std::vector<double> q = {trial_c, b, a};
                        if (b * b - 4.0 * a * trial_c < -1e-7) {
                            std::vector<double> powered = polynomial_power_coefficients(q, m);
                            // 由于 LC 可能不对，我们需要标准化比较
                            double scale = remaining.back() / powered.back();
                            for (double& val : powered) val *= scale;
                            
                            if (polynomial_close(powered, remaining, 1e-4)) {
                                quadratic_factors->push_back({q, m});
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }

    // 如果无法提取完整的二次分解，返回 false 以便回退到数值积分
    return remaining.size() == 1;
}

bool polynomial_close(const std::vector<double>& lhs,
                      const std::vector<double>& rhs,
                      double eps) {
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

bool integrate_general_partial_fractions(
    const std::vector<double>& numerator,
    const std::vector<double>& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    if (denominator.size() < 2) {
        return false;
    }

    std::vector<LinearFactorMultiplicity> linear_factors;
    std::vector<QuadraticFactorMultiplicity> quadratic_factors;
    if (!extract_general_rational_factorization(denominator,
                                                &linear_factors,
                                                &quadratic_factors)) {
        // 如果无法完全分解，至少尝试处理已经提取出的部分
        if (linear_factors.empty() && quadratic_factors.empty()) {
            return false;
        }
    }

    std::vector<RationalPartialFractionTerm> terms;
    // ... (保持项生成逻辑不变)
    for (const auto& factor : linear_factors) {
        for (int p = 1; p <= factor.multiplicity; ++p) {
            RationalPartialFractionTerm term;
            term.kind = RationalPartialFractionTerm::Kind::kLinear;
            term.root = factor.root;
            term.power = p;
            terms.push_back(term);
        }
    }
    for (const auto& factor : quadratic_factors) {
        for (int p = 1; p <= factor.multiplicity; ++p) {
            RationalPartialFractionTerm slope_term;
            slope_term.kind = RationalPartialFractionTerm::Kind::kQuadratic;
            slope_term.quadratic = factor.coefficients;
            slope_term.power = p;
            slope_term.numerator_degree = 1;
            terms.push_back(slope_term);

            RationalPartialFractionTerm constant_term = slope_term;
            constant_term.numerator_degree = 0;
            terms.push_back(constant_term);
        }
    }

    const int unknown_count = static_cast<int>(terms.size());
    if (unknown_count == 0) return false;

    // 构建系数矩阵。我们通过在多个点采样并代入恒等式来构建线性方程组。
    // P(x)/Q(x) = sum(A_i / (x-r_i)^p_i) + sum((B_j x + C_j) / (q_j(x))^p_j)
    // P(x) = sum(A_i * Q(x)/(x-r_i)^p_i) + sum((B_j x + C_j) * Q(x)/(q_j(x))^p_j)
    
    std::vector<std::vector<double>> columns;
    for (const auto& term : terms) {
        std::vector<double> divisor;
        if (term.kind == RationalPartialFractionTerm::Kind::kLinear) {
            divisor = polynomial_power_coefficients({-term.root, 1.0}, term.power);
        } else {
            divisor = polynomial_power_coefficients(term.quadratic, term.power);
        }
        PolynomialDivisionResult div = polynomial_divide(denominator, divisor);
        std::vector<double> col = div.quotient;
        if (term.kind == RationalPartialFractionTerm::Kind::kQuadratic && term.numerator_degree == 1) {
            col.insert(col.begin(), 0.0);
        }
        columns.push_back(col);
    }

    std::vector<std::vector<double>> matrix;
    std::vector<double> rhs;
    for (int c = -unknown_count * 2; static_cast<int>(rhs.size()) < unknown_count && c < 1000; ++c) {
        double x = static_cast<double>(c);
        bool is_pole = false;
        for (const auto& f : linear_factors) {
            if (mymath::is_near_zero(x - f.root, 1e-8)) { is_pole = true; break; }
        }
        if (is_pole) continue;

        std::vector<double> row;
        for (const auto& col : columns) {
            double val = polynomial_evaluate(col, x);
            SymbolicExpression cleaned = clean_symbolic_constant(val);
            double cleaned_val = val;
            expr_is_number(cleaned, &cleaned_val);
            row.push_back(cleaned_val);
        }
        matrix.push_back(row);
        double rhs_val = polynomial_evaluate(numerator, x);
        SymbolicExpression cleaned_rhs = clean_symbolic_constant(rhs_val);
        double cleaned_rhs_val = rhs_val;
        expr_is_number(cleaned_rhs, &cleaned_rhs_val);
        rhs.push_back(cleaned_rhs_val);
    }

    std::vector<double> coeffs;
    if (!solve_dense_linear_system(matrix, rhs, &coeffs)) return false;

    SymbolicExpression result = SymbolicExpression::number(0.0);
    std::vector<SymbolicExpression> term_exprs;
    for (size_t i = 0; i < terms.size(); ++i) {
        SymbolicExpression val = clean_symbolic_constant(coeffs[i]);
        if (expr_is_zero(val)) continue;

        SymbolicExpression term_int;
        if (terms[i].kind == RationalPartialFractionTerm::Kind::kLinear) {
            const SymbolicExpression v = make_subtract(SymbolicExpression::variable(variable_name),
                                                       SymbolicExpression::number(terms[i].root)).simplify();
            if (terms[i].power == 1) {
                term_int = make_multiply(val,
                                         make_function("ln", make_function("abs", v))).simplify();
            } else {
                double p_val = static_cast<double>(terms[i].power);
                term_int = make_multiply(make_divide(val, SymbolicExpression::number(1.0 - p_val)),
                                         make_power(v, SymbolicExpression::number(1.0 - p_val))).simplify();
            }
        } else {
            double slope_val = (terms[i].numerator_degree == 1) ? coeffs[i] : 0.0;
            double const_val = (terms[i].numerator_degree == 0) ? coeffs[i] : 0.0;
            term_int = integrate_quadratic_partial_fraction_term(terms[i].quadratic, slope_val, const_val, terms[i].power, variable_name);
        }
        term_exprs.push_back(term_int);
    }

    if (term_exprs.empty()) return false;
    
    result = term_exprs[0];
    for (size_t i = 1; i < term_exprs.size(); ++i) {
        result = make_add(result, term_exprs[i]);
    }
    *integrated = result.simplify();
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
        SymbolicExpression clean_derivative_scale = clean_symbolic_constant(derivative_scale);
        SymbolicExpression derivative_part;
        if (power == 1) {
            derivative_part =
                make_multiply(clean_derivative_scale,
                              make_function("ln",
                                            make_function("abs", quadratic_expression)))
                    .simplify();
        } else {
            double p_val = static_cast<double>(power);
            derivative_part =
                make_multiply(make_divide(clean_derivative_scale, SymbolicExpression::number(1.0 - p_val)),
                              make_power(quadratic_expression,
                                         SymbolicExpression::number(1.0 - p_val)))
                    .simplify();
        }
        result = make_add(result, derivative_part).simplify();
    }

    if (!mymath::is_near_zero(inverse_scale, kFormatEps)) {
        result =
            make_add(result,
                     make_multiply(clean_symbolic_constant(inverse_scale),
                                   integrate_inverse_quadratic_power(quadratic,
                                                                    power,
                                                                    variable_name)))
                .simplify();
    }
    return result.simplify();
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
        *primitive = make_function("ln", make_function("abs", make_function("sin", argument))).simplify();
        return true;
    }
    if (function_name == "ln") {
        *primitive = make_subtract(make_multiply(argument, make_function("ln", argument)), argument).simplify();
        return true;
    }
    if (function_name == "asin") {
        *primitive = make_add(make_multiply(argument, make_function("asin", argument)),
                              make_function("sqrt", make_subtract(SymbolicExpression::number(1.0),
                                                                make_power(argument, SymbolicExpression::number(2.0)))))
                        .simplify();
        return true;
    }
    if (function_name == "acos") {
        *primitive = make_subtract(make_multiply(argument, make_function("acos", argument)),
                                   make_function("sqrt", make_subtract(SymbolicExpression::number(1.0),
                                                                     make_power(argument, SymbolicExpression::number(2.0)))))
                        .simplify();
        return true;
    }
    if (function_name == "atan") {
        *primitive = make_subtract(make_multiply(argument, make_function("atan", argument)),
                                   make_multiply(SymbolicExpression::number(0.5),
                                               make_function("ln", make_add(SymbolicExpression::number(1.0),
                                                                          make_power(argument, SymbolicExpression::number(2.0))))))
                        .simplify();
        return true;
    }
    if (function_name == "sinh") {
        *primitive = make_function("cosh", argument);
        return true;
    }
    if (function_name == "cosh") {
        *primitive = make_function("sinh", argument);
        return true;
    }
    if (function_name == "tanh") {
        *primitive = make_function("ln", make_function("cosh", argument)).simplify();
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

bool try_integrate_by_parts(const SymbolicExpression& left,
                            const SymbolicExpression& right,
                            const std::string& variable_name,
                            SymbolicExpression* integrated) {
    struct IbpState {
        std::string original_key;
        SymbolicExpression uv_sum = SymbolicExpression::number(0.0);
    };
    static thread_local std::vector<IbpState> ibp_stack;

    const SymbolicExpression original_expr = make_multiply(left, right).simplify();
    const std::string original_key = node_structural_key(original_expr.node_);

    // 检查是否已经在栈中（循环检测）
    for (const auto& state : ibp_stack) {
        if (state.original_key == original_key) {
            // 发现循环！我们不能直接积分，而是抛出一个特定异常或返回 false
            // 实际上，这里需要通过代数手段解决。
            // 简单起见，如果是在第二层递归发现循环，我们可以标记它。
            return false; 
        }
    }

    // LIATE 优先级
    auto get_priority = [&](const SymbolicExpression& e) {
        if (e.node_->type == NodeType::kFunction) {
            if (e.node_->text == "ln") return 1;
            if (e.node_->text == "asin" || e.node_->text == "atan") return 2;
            if (e.node_->text == "sin" || e.node_->text == "cos" || e.node_->text == "tan") return 4;
            if (e.node_->text == "exp") return 5;
        }
        if (is_symbolic_polynomial(e, variable_name)) return 3;
        return 6;
    };

    SymbolicExpression u = left;
    SymbolicExpression dv = right;
    if (get_priority(right) < get_priority(left)) {
        u = right;
        dv = left;
    }

    SymbolicExpression v;
    try {
        v = dv.integral(variable_name);
    } catch (...) {
        u = (u.node_ == left.node_) ? right : left;
        dv = (dv.node_ == right.node_) ? left : right;
        try {
            v = dv.integral(variable_name);
        } catch (...) {
            return false;
        }
    }

    const SymbolicExpression du = u.derivative(variable_name).simplify();
    const SymbolicExpression v_du = make_multiply(v, du).simplify();

    // 检查 v_du 是否是原始表达式的常数倍 (I = uv - k*I => I = uv/(1+k))
    double k = 0.0;
    SymbolicExpression remainder;
    if (decompose_constant_times_expression(v_du, variable_name, &k, &remainder) &&
        node_structural_key(remainder.node_) == original_key) {
        if (!mymath::is_near_zero(1.0 + k, kFormatEps)) {
            *integrated = make_divide(make_multiply(u, v),
                                      SymbolicExpression::number(1.0 + k)).simplify();
            return true;
        }
    }

    try {
        *integrated = make_subtract(make_multiply(u, v),
                                    v_du.integral(variable_name))
                          .simplify();
        return true;
    } catch (...) {
    }

    // 如果没有直接发现循环，尝试递归积分 v_du
    if (ibp_stack.size() >= 2) return false;

    ibp_stack.push_back({original_key, make_multiply(u, v)});
    try {
        SymbolicExpression second_integral;
        // 在这里，我们需要处理二阶循环： I = uv - (u2v2 - kI)
        // 这需要更复杂的逻辑。暂时先处理一阶循环。
        // 为了支持 exp(x)*cos(x)，我们需要两步 IBP。
        
        // 我们尝试手动展开一轮
        // I = u*v - integral(v*du)
        // 令 I2 = integral(v*du)
        SymbolicExpression u2, dv2;
        // 简单假设 v_du 也是乘积
        if (v_du.node_->type == NodeType::kMultiply) {
            u2 = SymbolicExpression(v_du.node_->left);
            dv2 = SymbolicExpression(v_du.node_->right);
            
            SymbolicExpression v2;
            try {
                v2 = dv2.integral(variable_name);
                SymbolicExpression du2 = u2.derivative(variable_name).simplify();
                SymbolicExpression v2_du2 = make_multiply(v2, du2).simplify();
                
                // 检查 v2_du2 是否回到 I
                double k2 = 0.0;
                SymbolicExpression remainder2;
                if (decompose_constant_times_expression(v2_du2, variable_name, &k2, &remainder2) &&
                    node_structural_key(remainder2.node_) == original_key) {
                    // I = uv - (u2v2 - k2*I) => I = (uv - u2v2) / (1 - k2)
                    if (!mymath::is_near_zero(1.0 - k2, kFormatEps)) {
                        *integrated = make_divide(make_subtract(make_multiply(u, v),
                                                                make_multiply(u2, v2)),
                                                  SymbolicExpression::number(1.0 - k2)).simplify();
                        ibp_stack.pop_back();
                        return true;
                    }
                }
            } catch (...) {}
        }

        ibp_stack.pop_back();
        return false;
    } catch (...) {
        ibp_stack.pop_back();
        return false;
    }
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
    if (left.node_->type == NodeType::kFunction &&
        right.node_->type == NodeType::kFunction &&
        same_simplified_expression(SymbolicExpression(left.node_->left),
                                   SymbolicExpression(right.node_->left))) {
        const SymbolicExpression argument(left.node_->left);
        double a = 0.0;
        double b = 0.0;
        if (decompose_linear(argument, variable_name, &a, &b) &&
            !mymath::is_near_zero(a, kFormatEps)) {
            if ((left.node_->text == "sec" && right.node_->text == "tan") ||
                (left.node_->text == "tan" && right.node_->text == "sec")) {
                *integrated =
                    make_divide(make_function("sec", argument),
                                SymbolicExpression::number(a))
                        .simplify();
                return true;
            }
            if ((left.node_->text == "csc" && right.node_->text == "cot") ||
                (left.node_->text == "cot" && right.node_->text == "csc")) {
                *integrated =
                    make_divide(make_negate(make_function("csc", argument)),
                                SymbolicExpression::number(a))
                        .simplify();
                return true;
            }
        }
    }

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
    std::vector<double> num_coeffs;
    std::vector<double> den_coeffs;
    if (!polynomial_coefficients_from_simplified(numerator.simplify(),
                                                 variable_name,
                                                 &num_coeffs) ||
        !polynomial_coefficients_from_simplified(denominator.simplify(),
                                                 variable_name,
                                                 &den_coeffs)) {
        return false;
    }
    trim_coefficients(&num_coeffs);
    trim_coefficients(&den_coeffs);
    if (den_coeffs.size() <= 1 || polynomial_is_zero(den_coeffs)) return false;

    if (den_coeffs.back() < 0.0) {
        for (double& c : num_coeffs) c = -c;
        for (double& c : den_coeffs) c = -c;
    }

    SymbolicExpression poly_part = SymbolicExpression::number(0.0);
    bool has_poly = false;
    if (num_coeffs.size() >= den_coeffs.size()) {
        const PolynomialDivisionResult div = polynomial_divide(num_coeffs, den_coeffs);
        if (!polynomial_is_zero(div.quotient)) {
            poly_part = build_polynomial_expression_from_coefficients(div.quotient, variable_name)
                            .integral(variable_name).simplify();
            has_poly = true;
        }
        num_coeffs = div.remainder;
        trim_coefficients(&num_coeffs);
    }

    if (polynomial_is_zero(num_coeffs)) {
        *integrated = poly_part;
        return has_poly;
    }

    SymbolicExpression special;
    if (try_integrate_repeated_unit_quadratic(num_coeffs, den_coeffs, variable_name, &special)) {
        *integrated = has_poly ? make_add(poly_part, special).simplify() : special;
        return true;
    }

    SymbolicExpression partial;
    if (integrate_general_partial_fractions(num_coeffs, den_coeffs, variable_name, &partial)) {
        *integrated = has_poly ? make_add(poly_part, partial).simplify() : partial;
        return true;
    }

    // 回退：处理简单的线性分母
    if (den_coeffs.size() == 2) {
        const double constant = num_coeffs[0];
        const double slope = den_coeffs[1];
        if (!mymath::is_near_zero(slope, kFormatEps)) {
            const SymbolicExpression den_expr = 
                build_polynomial_expression_from_coefficients(den_coeffs, variable_name);
            const SymbolicExpression linear_int = 
                make_multiply(SymbolicExpression::number(constant / slope),
                              make_function("ln", make_function("abs", den_expr))).simplify();
            *integrated = has_poly ? make_add(poly_part, linear_int).simplify() : linear_int;
            return true;
        }
    }

    return false;
}


bool try_integrate_weierstrass_substitution(const SymbolicExpression& expression,
                                             const std::string& variable_name,
                                             SymbolicExpression* integrated) {
    // 检查是否为 1 / (trig expression)
    if (expression.node_->type != NodeType::kDivide) {
        return false;
    }

    const SymbolicExpression numerator = SymbolicExpression(expression.node_->left);
    const SymbolicExpression denominator = SymbolicExpression(expression.node_->right);

    if (!numerator.is_constant(variable_name)) {
        return false;
    }

    // 尝试识别 A + B*sin(ax+b) + C*cos(ax+b)
    // 这里我们先做一个简化版的：检查分母是否只包含 sin/cos 的线性组合
    std::vector<SymbolicExpression> terms;
    collect_additive_expressions(denominator, &terms);

    double A = 0.0;
    SymbolicExpression B_expr = SymbolicExpression::number(0.0);
    SymbolicExpression C_expr = SymbolicExpression::number(0.0);
    SymbolicExpression argument;
    bool found_argument = false;

    for (const auto& term : terms) {
        double val = 0.0;
        if (term.is_number(&val)) {
            A += val;
        } else if (term.node_->type == NodeType::kFunction && term.node_->text == "sin") {
            if (!found_argument) {
                argument = SymbolicExpression(term.node_->left);
                found_argument = true;
            } else if (!expressions_match(argument, SymbolicExpression(term.node_->left))) {
                return false;
            }
            B_expr = make_add(B_expr, SymbolicExpression::number(1.0)).simplify();
        } else if (term.node_->type == NodeType::kFunction && term.node_->text == "cos") {
            if (!found_argument) {
                argument = SymbolicExpression(term.node_->left);
                found_argument = true;
            } else if (!expressions_match(argument, SymbolicExpression(term.node_->left))) {
                return false;
            }
            C_expr = make_add(C_expr, SymbolicExpression::number(1.0)).simplify();
        } else if (term.node_->type == NodeType::kMultiply) {
            // 处理 B * sin(x) 或 C * cos(x)
            double coeff = 0.0;
            SymbolicExpression rest;
            if (decompose_constant_times_expression(term, variable_name, &coeff, &rest)) {
                if (rest.node_->type == NodeType::kFunction && rest.node_->text == "sin") {
                    if (!found_argument) {
                        argument = SymbolicExpression(rest.node_->left);
                        found_argument = true;
                    } else if (!expressions_match(argument, SymbolicExpression(rest.node_->left))) {
                        return false;
                    }
                    B_expr = make_add(B_expr, SymbolicExpression::number(coeff)).simplify();
                } else if (rest.node_->type == NodeType::kFunction && rest.node_->text == "cos") {
                    if (!found_argument) {
                        argument = SymbolicExpression(rest.node_->left);
                        found_argument = true;
                    } else if (!expressions_match(argument, SymbolicExpression(rest.node_->left))) {
                        return false;
                    }
                    C_expr = make_add(C_expr, SymbolicExpression::number(coeff)).simplify();
                } else {
                    return false;
                }
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    if (!found_argument) return false;

    double a = 0.0, b = 0.0;
    if (!decompose_linear(argument, variable_name, &a, &b) || !mymath::is_near_zero(b, kFormatEps)) {
        // 目前仅支持 sin(ax), cos(ax)
        return false;
    }

    // Weierstrass 置换: t = tan(ax/2)
    // sin(ax) = 2t / (1+t^2)
    // cos(ax) = (1-t^2) / (1+t^2)
    // dx = 2 / (a * (1+t^2)) dt
    
    // 积分变为: integral( (2/a) / (A*(1+t^2) + B*(2t) + C*(1-t^2)), t )
    // 分母多项式: (A-C)t^2 + 2Bt + (A+C)
    double b_val = 0.0, c_val = 0.0;
    if (!B_expr.is_number(&b_val) || !C_expr.is_number(&c_val)) return false;

    std::vector<double> poly_numerator = { 2.0 / a };
    std::vector<double> poly_denominator = { A + c_val, 2.0 * b_val, A - c_val };

    SymbolicExpression t_variable = SymbolicExpression::variable("_t");
    SymbolicExpression t_numerator = build_polynomial_expression_from_coefficients(poly_numerator, "_t");
    SymbolicExpression t_denominator = build_polynomial_expression_from_coefficients(poly_denominator, "_t");

    SymbolicExpression t_integrated;
    if (try_integrate_polynomial_quotient(t_numerator, t_denominator, "_t", &t_integrated)) {
        // 替换 t 为 tan(ax/2)
        SymbolicExpression substitution = make_function("tan", make_divide(argument, SymbolicExpression::number(2.0)));
        *integrated = substitute_impl(t_integrated, "_t", substitution).simplify();
        return true;
    }

    return false;
}

SymbolicExpression derivative_uncached(const SymbolicExpression& expression,
                                       const std::string& variable_name) {
    const std::shared_ptr<SymbolicExpression::Node>& node_ = expression.node_;
    switch (node_->type) {
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
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
            if (node_->text == "erf") {
                const SymbolicExpression pi_val = SymbolicExpression::variable("pi");
                const SymbolicExpression exp_part = make_function("exp", make_negate(make_power(argument, SymbolicExpression::number(2.0))));
                const SymbolicExpression factor = make_divide(SymbolicExpression::number(2.0), make_function("sqrt", pi_val));
                return make_multiply(make_multiply(factor, exp_part), inner).simplify();
            }
            if (node_->text == "erfc") {
                const SymbolicExpression pi_val = SymbolicExpression::variable("pi");
                const SymbolicExpression exp_part = make_function("exp", make_negate(make_power(argument, SymbolicExpression::number(2.0))));
                const SymbolicExpression factor = make_divide(SymbolicExpression::number(-2.0), make_function("sqrt", pi_val));
                return make_multiply(make_multiply(factor, exp_part), inner).simplify();
            }
            if (node_->text == "erfi") {
                const SymbolicExpression pi_val = SymbolicExpression::variable("pi");
                const SymbolicExpression exp_part = make_function("exp", make_power(argument, SymbolicExpression::number(2.0)));
                const SymbolicExpression factor = make_divide(SymbolicExpression::number(2.0), make_function("sqrt", pi_val));
                return make_multiply(make_multiply(factor, exp_part), inner).simplify();
            }
            if (node_->text == "Ei") {
                return make_multiply(make_divide(make_function("exp", argument), argument), inner).simplify();
            }
            if (node_->text == "Si") {
                return make_multiply(make_divide(make_function("sin", argument), argument), inner).simplify();
            }
            if (node_->text == "Ci") {
                return make_multiply(make_divide(make_function("cos", argument), argument), inner).simplify();
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
    const std::size_t n = variable_names.size();
    std::vector<SymbolicExpression> first_derivatives;
    first_derivatives.reserve(n);
    for (const std::string& variable_name : variable_names) {
        first_derivatives.push_back(derivative(variable_name).simplify());
    }

    std::vector<std::vector<SymbolicExpression>> result(
        n, std::vector<SymbolicExpression>(n, SymbolicExpression::number(0.0)));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = row; col < n; ++col) {
            const SymbolicExpression value =
                first_derivatives[row].derivative(variable_names[col]).simplify();
            result[row][col] = value;
            if (row != col) {
                result[col][row] = value;
            }
        }
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

SymbolicExpression SymbolicExpression::divergence(
    const std::vector<SymbolicExpression>& components,
    const std::vector<std::string>& variable_names) {
    if (components.size() != variable_names.size()) {
        throw std::runtime_error("divergence requires vector field components to match variable names count");
    }
    SymbolicExpression result = number(0.0);
    for (std::size_t i = 0; i < components.size(); ++i) {
        result = make_add(result, components[i].derivative(variable_names[i])).simplify();
    }
    return result;
}

std::vector<SymbolicExpression> SymbolicExpression::curl(
    const std::vector<SymbolicExpression>& components,
    const std::vector<std::string>& variable_names) {
    if (components.size() != 3 || variable_names.size() != 3) {
        throw std::runtime_error("curl is currently only supported for 3D vector fields");
    }
    
    // components: [Fx, Fy, Fz], variables: [x, y, z]
    // curl = [dFz/dy - dFy/dz, dFx/dz - dFz/dx, dFy/dx - dFx/dy]
    return {
        make_subtract(components[2].derivative(variable_names[1]), components[1].derivative(variable_names[2])).simplify(),
        make_subtract(components[0].derivative(variable_names[2]), components[2].derivative(variable_names[0])).simplify(),
        make_subtract(components[1].derivative(variable_names[0]), components[0].derivative(variable_names[1])).simplify()
    };
}

SymbolicExpression SymbolicExpression::laplacian(const std::vector<std::string>& variable_names) const {
    SymbolicExpression result = number(0.0);
    for (const std::string& variable_name : variable_names) {
        result = make_add(result, derivative(variable_name).derivative(variable_name)).simplify();
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
        case NodeType::kPi:
        case NodeType::kE:
            return make_multiply(*this, variable(variable_name)).simplify();
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
            if (try_integrate_by_parts(left,
                                       right,
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
        
        if (left.is_constant(variable_name)) {
            SymbolicExpression c_term, x2_coeff;
            // Case 1: 1 / (c + a*x^2) -> atan
            if (is_pure_quadratic(right, variable_name, &c_term, &x2_coeff)) {
                double a_val, c_val;
                if (x2_coeff.is_number(&a_val) && c_term.is_number(&c_val)) {
                    if (a_val * c_val > 0) {
                        const double factor = 1.0 / mymath::sqrt(a_val * c_val);
                        const double internal = mymath::sqrt(a_val / c_val);
                        return make_multiply(left,
                            make_multiply(number(factor),
                                make_function("atan", make_multiply(number(internal), variable(variable_name)))))
                            .simplify();
                    } else if (a_val * c_val < 0) {
                        // 1 / (x^2 - 1) -> partial fractions (ln)
                        // This is handled by try_integrate_polynomial_quotient below
                    }
                } else if (expr_is_one(x2_coeff) && expr_is_one(c_term)) {
                    return make_multiply(left, make_function("atan", variable(variable_name))).simplify();
                }
            }
            
            // Case 2: 1 / sqrt(c - a*x^2) -> asin
            if (right.node_->type == NodeType::kFunction && right.node_->text == "sqrt") {
                const SymbolicExpression inner(right.node_->left);
                if (is_pure_quadratic(inner, variable_name, &c_term, &x2_coeff)) {
                    double a_val, c_val;
                    if (x2_coeff.is_number(&a_val) && c_term.is_number(&c_val)) {
                        if (c_val > 0 && a_val < 0) {
                            const double abs_a = -a_val;
                            const double factor = 1.0 / mymath::sqrt(abs_a);
                            const double internal = mymath::sqrt(abs_a / c_val);
                            return make_multiply(left,
                                make_multiply(number(factor),
                                    make_function("asin", make_multiply(number(internal), variable(variable_name)))))
                                .simplify();
                        }
                    } else if (expr_is_one(c_term) && expr_is_minus_one(x2_coeff)) {
                        return make_multiply(left, make_function("asin", variable(variable_name))).simplify();
                    }
                }
            }

            // Case 3: 1 / (ax + b) -> (1/a) * ln|ax + b|
            SymbolicExpression a_expr, b_expr;
            if (symbolic_decompose_linear(right, variable_name, &a_expr, &b_expr) && !expr_is_zero(a_expr)) {
                return make_divide(make_multiply(left, make_function("ln", make_function("abs", right))),
                                   a_expr)
                    .simplify();
            }
        }

        if (right.node_->type == NodeType::kFunction && right.node_->text == "sqrt") {
            const SymbolicExpression inner(right.node_->left);
            const SymbolicExpression expected_derivative =
                inner.derivative(variable_name).simplify();
            double scale = 1.0;
            bool matched = same_simplified_expression(left, expected_derivative);
            if (!matched) {
                double constant = 0.0;
                SymbolicExpression rest;
                if (decompose_constant_times_expression(expected_derivative,
                                                        variable_name,
                                                        &constant,
                                                        &rest) &&
                    !mymath::is_near_zero(constant, kFormatEps) &&
                    same_simplified_expression(left, rest)) {
                    scale = 1.0 / constant;
                    matched = true;
                } else if (decompose_constant_times_expression(left,
                                                               variable_name,
                                                               &constant,
                                                               &rest) &&
                           same_simplified_expression(rest, expected_derivative)) {
                    scale = constant;
                    matched = true;
                }
            }
            if (matched) {
                return make_multiply(SymbolicExpression::number(2.0 * scale),
                                     make_function("sqrt", inner))
                    .simplify();
            }

            SymbolicExpression c_term;
            SymbolicExpression x2_coeff;
            if (is_pure_quadratic(inner, variable_name, &c_term, &x2_coeff)) {
                double a_value = 0.0;
                if (x2_coeff.is_number(&a_value) &&
                    !mymath::is_near_zero(a_value, kFormatEps) &&
                    left.is_variable_named(variable_name)) {
                    return make_divide(make_function("sqrt", inner),
                                       SymbolicExpression::number(a_value))
                        .simplify();
                }
            }
        }
        
        if (left.node_->type == NodeType::kFunction && left.node_->text == "exp" &&
            right.is_variable_named(variable_name)) {
            const SymbolicExpression argument(left.node_->left);
            if (argument.is_variable_named(variable_name)) {
                return make_function("Ei", argument);
            }
        }
        if (left.node_->type == NodeType::kFunction && left.node_->text == "sin" &&
            right.is_variable_named(variable_name)) {
            const SymbolicExpression argument(left.node_->left);
            if (argument.is_variable_named(variable_name)) {
                return make_function("Si", argument);
            }
        }
        if (left.node_->type == NodeType::kFunction && left.node_->text == "cos" &&
            right.is_variable_named(variable_name)) {
            const SymbolicExpression argument(left.node_->left);
            if (argument.is_variable_named(variable_name)) {
                return make_function("Ci", argument);
            }
        }

        SymbolicExpression rational_integral;
        if (try_integrate_polynomial_quotient(left,
                                              right,
                                              variable_name,
                                              &rational_integral)) {
            return rational_integral.simplify();
        }
        if (try_integrate_weierstrass_substitution(*this,
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
        SymbolicExpression trig_identity_integral;
        if (exponent.is_number(&exponent_value) &&
            try_integrate_trig_power_identity(base,
                                              exponent_value,
                                              variable_name,
                                              &trig_identity_integral)) {
            return trig_identity_integral.simplify();
        }
        
        SymbolicExpression a_expr, b_expr;
        if (exponent.is_constant(variable_name) &&
            symbolic_decompose_linear(base, variable_name, &a_expr, &b_expr) &&
            !expr_is_zero(a_expr)) {
            if (expr_is_minus_one(exponent)) {
                return make_divide(make_function("ln", make_function("abs", base)),
                                   a_expr)
                    .simplify();
            }
            const SymbolicExpression new_exponent = make_add(exponent, number(1.0)).simplify();
            return make_divide(make_power(base, new_exponent),
                               make_multiply(a_expr, new_exponent))
                .simplify();
        }
        throw std::runtime_error("symbolic integral only supports powers of linear terms or certain trig identities");
    }

    if (node_->type == NodeType::kFunction) {
        const SymbolicExpression argument(node_->left);
        double a = 0.0;
        double b = 0.0;
        const bool linear = decompose_linear(argument, variable_name, &a, &b) &&
                            !mymath::is_near_zero(a, kFormatEps);
        
        SymbolicExpression primitive;
        if (linear && primitive_for_outer_function(node_->text, argument, &primitive)) {
            return make_divide(primitive, number(a)).simplify();
        }

        SymbolicExpression symbolic_a;
        SymbolicExpression symbolic_b;
        if (!linear &&
            symbolic_decompose_linear(argument, variable_name, &symbolic_a, &symbolic_b) &&
            !expr_is_zero(symbolic_a) &&
            primitive_for_outer_function(node_->text, argument, &primitive)) {
            return make_divide(primitive, symbolic_a).simplify();
        }

        if (node_->text == "sin" && linear) {
            return make_divide(make_negate(make_function("cos", argument)),
                               number(a))
                .simplify();
        }
        if (node_->text == "cos" && linear) {
            return make_divide(make_function("sin", argument), number(a)).simplify();
        }
        if (node_->text == "exp") {
            if (linear) {
                return make_divide(make_function("exp", argument), number(a)).simplify();
            }
            SymbolicExpression c_term, x2_coeff;
            if (is_pure_quadratic(argument, variable_name, &c_term, &x2_coeff)) {
                double a_val;
                if (x2_coeff.is_number(&a_val)) {
                    const SymbolicExpression pi_val = variable("pi");
                    if (a_val < 0) {
                        const double pos_a = -a_val;
                        const SymbolicExpression factor = make_divide(
                            make_multiply(make_function("exp", c_term), make_function("sqrt", pi_val)),
                            make_multiply(number(2.0), make_function("sqrt", number(pos_a)))
                        );
                        return make_multiply(
                            factor,
                            make_function("erf", make_multiply(make_function("sqrt", number(pos_a)), variable(variable_name)))
                        ).simplify();
                    } else if (a_val > 0) {
                        const SymbolicExpression factor = make_divide(
                            make_multiply(make_function("exp", c_term), make_function("sqrt", pi_val)),
                            make_multiply(number(2.0), make_function("sqrt", number(a_val)))
                        );
                        return make_multiply(
                            factor,
                            make_function("erfi", make_multiply(make_function("sqrt", number(a_val)), variable(variable_name)))
                        ).simplify();
                    }
                }
            }
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
        if (node_->text == "sqrt") {
            SymbolicExpression a_quad, b_quad, c_quad;
            if (is_general_quadratic(argument, variable_name, &a_quad, &b_quad, &c_quad) &&
                expr_is_one(c_quad) && expr_is_zero(b_quad) && expr_is_minus_one(a_quad)) {
                const SymbolicExpression x = variable(variable_name);
                return make_divide(
                           make_add(make_multiply(x, make_function("sqrt", argument)),
                                    make_function("asin", x)),
                           number(2.0))
                    .simplify();
            }
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

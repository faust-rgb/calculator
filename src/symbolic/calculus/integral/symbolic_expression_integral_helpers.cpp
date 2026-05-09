// ============================================================================
// 符号积分辅助函数
// ============================================================================
//
// 本文件包含符号积分所需的辅助函数：
// - 部分分式分解
// - 换元积分法辅助
// - 分部积分法辅助
// - 三角函数积分辅助
// - Weierstrass 置换辅助
// ============================================================================

#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/algebra/polynomial/symbolic_polynomial.h"
#include "symbolic/calculus/integral/symbolic_expression_integral_helpers.h"
#include "symbolic/calculus/integral/symbolic_expression_integral_internal.h"

#include "math/mymath.h"
#include "polynomial/polynomial.h"

#include <algorithm>
#include <stdexcept>

namespace symbolic_expression_internal {

namespace {

// ============================================================================// ============================================================================

/**
 * @struct LinearFactorMultiplicity
 * @brief 线性因子及其重数
 */
struct LinearFactorMultiplicity {
    Scalar root = 0.0L;
    int multiplicity = 0;
};

/**
 * @struct QuadraticFactorMultiplicity
 * @brief 二次因子及其重数
 */
struct QuadraticFactorMultiplicity {
    std::vector<Scalar> coefficients;
    int multiplicity = 0;
};
} // namespace

// ============================================================================
// 多项式因子提取
// ============================================================================

/**
 * @brief 检查两个多项式系数向量是否近似相等
 */
bool polynomial_close(const std::vector<Scalar>& lhs,
                      const std::vector<Scalar>& rhs,
                      Scalar eps = 1e-8);

/**
 * @brief 多项式整除检测
 */
bool divide_exact_polynomial(std::vector<Scalar>* remaining,
                             const std::vector<Scalar>& divisor);

/**
 * @brief 容差多项式整除检测
 */
bool divide_near_polynomial(std::vector<Scalar>* remaining,
                            const std::vector<Scalar>& divisor,
                            Scalar eps = 1e-7);

/**
 * @brief 积分二次式部分分式项
 */
SymbolicExpression integrate_quadratic_partial_fraction_term(
    const std::vector<Scalar>& quadratic,
    Scalar slope,
    Scalar constant,
    int power,
    const std::string& variable_name);

/**
 * @brief 提取有理式的一般因子分解
 *
 * 将分母多项式分解为线性因子和不可约二次因子的乘积。
 * 用于部分分式积分。
 */
bool extract_general_rational_factorization(
    const std::vector<Scalar>& denominator,
    std::vector<LinearFactorMultiplicity>* linear_factors,
    std::vector<QuadraticFactorMultiplicity>* quadratic_factors) {
    linear_factors->clear();
    quadratic_factors->clear();

    std::vector<Scalar> remaining = denominator;
    trim_coefficients(&remaining);
    
    // 1. 提取线性因子
    std::vector<Scalar> roots;
    try {
        roots = polynomial_real_roots(remaining);
    } catch (...) {
        roots.clear();
    }

    for (Scalar root : roots) {
        SymbolicExpression cleaned_root_expr = clean_symbolic_constant(root);
        Scalar cleaned_root = root;
        expr_is_number(cleaned_root_expr, &cleaned_root);
        const std::vector<Scalar> factor = {-cleaned_root, 1.0L};
        int multiplicity = 0;
        while (remaining.size() > 1 && divide_exact_polynomial(&remaining, factor)) {
            ++multiplicity;
        }
        if (multiplicity > 0) {
            linear_factors->push_back(LinearFactorMultiplicity{cleaned_root, multiplicity});
        }
    }

    trim_coefficients(&remaining);
    if (remaining.size() == 1) return true;

    // 2. 使用复根共轭配对提取一般实不可约二次因子。
    // 对实系数多项式，每一对非实共轭根 a±bi 对应
    // x^2 - 2a*x + (a^2+b^2)。这条路径覆盖多个不同二次因子
    // 及其重复因子，避免只识别 (x^2+k)^n 这类特殊形状。
    try {
        std::vector<mymath::complex<Scalar>> complex_roots =
            polynomial_complex_roots(remaining);
        std::vector<bool> used(complex_roots.size(), false);
        std::sort(complex_roots.begin(), complex_roots.end(),
                  [](const auto& lhs, const auto& rhs) {
                      if (mymath::abs(lhs.real() - rhs.real()) > 1e-8) {
                          return lhs.real() < rhs.real();
                      }
                      return lhs.imag() < rhs.imag();
                  });

        std::vector<std::vector<Scalar>> extracted_quadratics;
        for (std::size_t i = 0; i < complex_roots.size(); ++i) {
            if (used[i] || mymath::abs(complex_roots[i].imag()) <= 1e-7) {
                continue;
            }
            std::size_t pair_index = complex_roots.size();
            const mymath::complex<Scalar> conjugate_root =
                mymath::conj(complex_roots[i]);
            for (std::size_t j = i + 1; j < complex_roots.size(); ++j) {
                if (used[j]) {
                    continue;
                }
                if (mymath::abs(complex_roots[j].real() - conjugate_root.real()) <= 1e-6 &&
                    mymath::abs(complex_roots[j].imag() - conjugate_root.imag()) <= 1e-6) {
                    pair_index = j;
                    break;
                }
            }
            if (pair_index == complex_roots.size()) {
                continue;
            }

            const Scalar real = (complex_roots[i].real() +
                                 complex_roots[pair_index].real()) * 0.5;
            const Scalar imag_abs =
                (mymath::abs(complex_roots[i].imag()) +
                 mymath::abs(complex_roots[pair_index].imag())) * 0.5;
            if (imag_abs <= 1e-7) {
                continue;
            }
            std::vector<Scalar> factor = {
                real * real + imag_abs * imag_abs,
                -2.0 * real,
                1.0L,
            };
            for (Scalar& coefficient : factor) {
                SymbolicExpression cleaned = clean_symbolic_constant(coefficient);
                Scalar cleaned_value = coefficient;
                if (expr_is_number(cleaned, &cleaned_value)) {
                    coefficient = cleaned_value;
                }
            }
            if (factor[1] * factor[1] - 4.0 * factor[2] * factor[0] >= -1e-6) {
                continue;
            }
            extracted_quadratics.push_back(factor);
            used[i] = true;
            used[pair_index] = true;
        }

        for (const auto& factor : extracted_quadratics) {
            int multiplicity = 0;
            while (remaining.size() > 1 &&
                   divide_near_polynomial(&remaining, factor, 1e-5)) {
                ++multiplicity;
            }
            if (multiplicity > 0) {
                bool merged = false;
                for (auto& existing : *quadratic_factors) {
                    if (polynomial_close(existing.coefficients, factor, 1e-5)) {
                        existing.multiplicity += multiplicity;
                        merged = true;
                        break;
                    }
                }
                if (!merged) {
                    quadratic_factors->push_back({factor, multiplicity});
                }
            }
        }

        trim_coefficients(&remaining);
        if (remaining.size() == 1) return true;
    } catch (...) {
        // Fall through to the older shape-specific quadratic detector.
    }

    // 3. 提取二次因子
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
                    Scalar a = mymath::pow(mymath::abs(remaining.back()), 1.0L / m);
                    if (remaining.back() < 0 && m % 2 != 0) a = -a;
                    Scalar c = mymath::pow(mymath::abs(remaining[0]), 1.0L / m);
                    // 尝试不同的 c 符号 (对于不可约，c/a 必须 > (b/2a)^2)
                    for (Scalar c_sign : {1.0L, -1.0L}) {
                        Scalar trial_c = c * c_sign;
                        // 估算 b：通过一次项或最高次项次一项
                        Scalar b = 0.0L;
                        if (degree > 2) {
                            // (ax^2 + bx + c)^m = (ax^2)^m + m(ax^2)^{m-1}(bx) + ...
                            // coeff[deg-1] = m * a^{m-1} * b
                            b = remaining[degree - 1] / (m * mymath::pow(a, m - 1));
                        } else {
                            b = remaining[1] / m; // 实际上 m=1 时就是 remaining[1]
                        }
                        
                        std::vector<Scalar> q = {trial_c, b, a};
                        if (b * b - 4.0 * a * trial_c < -1e-7) {
                            std::vector<Scalar> powered = polynomial_power_coefficients(q, m);
                            // 由于 LC 可能不对，我们需要标准化比较
                            Scalar scale = remaining.back() / powered.back();
                            for (Scalar& val : powered) val *= scale;
                            
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

/**
 * @brief 检查两个多项式系数向量是否近似相等
 */
bool polynomial_close(const std::vector<Scalar>& lhs,
                      const std::vector<Scalar>& rhs,
                      Scalar eps) {
    std::vector<Scalar> left = lhs;
    std::vector<Scalar> right = rhs;
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

/**
 * @brief 多项式整除检测
 *
 * 检查 remaining 是否能被 divisor 整除，
 * 如果能，返回商并更新 remaining。
 */
bool divide_exact_polynomial(std::vector<Scalar>* remaining,
                             const std::vector<Scalar>& divisor) {
    const PolynomialDivisionResult division = polynomial_divide(*remaining, divisor);
    if (!polynomial_is_zero(division.remainder)) {
        return false;
    }
    *remaining = division.quotient;
    trim_coefficients(remaining);
    return true;
}

bool divide_near_polynomial(std::vector<Scalar>* remaining,
                            const std::vector<Scalar>& divisor,
                            Scalar eps) {
    PolynomialDivisionResult division = polynomial_divide(*remaining, divisor);
    trim_coefficients(&division.remainder);
    for (Scalar coefficient : division.remainder) {
        if (!mymath::is_near_zero(coefficient, eps)) {
            return false;
        }
    }
    *remaining = division.quotient;
    trim_coefficients(remaining);
    return true;
}

/**
 * @struct RationalPartialFractionTerm
 * @brief 部分分式分解项
 *
 * 表示部分分式分解中的单个项：
 * - kLinear: A / (x-r)^p 形式
 * - kQuadratic: (Bx+C) / q(x)^p 形式
 */
struct RationalPartialFractionTerm {
    enum class Kind {
        kLinear,
        kQuadratic,
    };

    Kind kind = Kind::kLinear;
    Scalar root = 0.0L;
    std::vector<Scalar> quadratic;
    int power = 0;
    int numerator_degree = 0;
};

// ============================================================================
// 部分分式积分
// ============================================================================

/**
 * @brief 积分一般有理式（部分分式分解法）
 *
 * 将有理式 P(x)/Q(x) 分解为部分分式后逐项积分：
 * 1. 提取分母因子（线性因子和不可约二次因子）
 * 2. 计算部分分式系数（求解线性方程组）
 * 3. 逐项积分
 *
 * @param numerator 分子多项式系数
 * @param denominator 分母多项式系数
 * @param variable_name 积分变量
 * @param integrated 输出积分结果
 * @return true 如果成功积分
 */
bool integrate_general_partial_fractions(
    const std::vector<Scalar>& numerator,
    const std::vector<Scalar>& denominator,
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
    
    std::vector<std::vector<Scalar>> columns;
    for (const auto& term : terms) {
        std::vector<Scalar> divisor;
        if (term.kind == RationalPartialFractionTerm::Kind::kLinear) {
            divisor = polynomial_power_coefficients({-term.root, 1.0L}, term.power);
        } else {
            divisor = polynomial_power_coefficients(term.quadratic, term.power);
        }
        PolynomialDivisionResult div = polynomial_divide(denominator, divisor);
        std::vector<Scalar> col = div.quotient;
        if (term.kind == RationalPartialFractionTerm::Kind::kQuadratic && term.numerator_degree == 1) {
            col.insert(col.begin(), 0.0L);
        }
        columns.push_back(col);
    }

    std::vector<std::vector<Scalar>> matrix;
    std::vector<Scalar> rhs;
    for (int c = -unknown_count * 2; static_cast<int>(rhs.size()) < unknown_count && c < 1000; ++c) {
        Scalar x = (c);
        bool is_pole = false;
        for (const auto& f : linear_factors) {
            if (mymath::is_near_zero(x - f.root, 1e-8)) { is_pole = true; break; }
        }
        if (is_pole) continue;

        std::vector<Scalar> row;
        for (const auto& col : columns) {
            Scalar val = polynomial_evaluate(col, x);
            SymbolicExpression cleaned = clean_symbolic_constant(val);
            Scalar cleaned_val = val;
            expr_is_number(cleaned, &cleaned_val);
            row.push_back(cleaned_val);
        }
        matrix.push_back(row);
        Scalar rhs_val = polynomial_evaluate(numerator, x);
        SymbolicExpression cleaned_rhs = clean_symbolic_constant(rhs_val);
        Scalar cleaned_rhs_val = rhs_val;
        expr_is_number(cleaned_rhs, &cleaned_rhs_val);
        rhs.push_back(cleaned_rhs_val);
    }

    std::vector<Scalar> coeffs;
    if (!solve_dense_linear_system(matrix, rhs, &coeffs)) return false;

    SymbolicExpression result = SymbolicExpression::number(0.0L);
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
                Scalar p_val = (terms[i].power);
                term_int = make_multiply(make_divide(val, SymbolicExpression::number(1.0L - p_val)),
                                         make_power(v, SymbolicExpression::number(1.0L - p_val))).simplify();
            }
        } else {
            Scalar slope_val = (terms[i].numerator_degree == 1) ? coeffs[i] : 0.0L;
            Scalar const_val = (terms[i].numerator_degree == 0) ? coeffs[i] : 0.0L;
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

/**
 * @brief 积分逆二次式幂次
 *
 * 计算 1 / (ax^2 + bx + c)^n 的积分，其中判别式小于零（不可约）。
 * 使用递推公式：
 * I_n = 2ax / [(4ac-b^2)(n-1)(ax^2+bx+c)^(n-1)] + (4n-6)/[(2n-2)(4ac-b^2)] * I_{n-1}
 */
SymbolicExpression integrate_inverse_quadratic_power(
    const std::vector<Scalar>& quadratic,
    int power,
    const std::string& variable_name) {
    const Scalar c = quadratic[0];
    const Scalar b = quadratic[1];
    const Scalar a = quadratic[2];
    const Scalar delta = 4.0 * a * c - b * b;
    if (power <= 0 || mymath::is_near_zero(a, kFormatEps) || delta <= kFormatEps) {
        throw std::runtime_error("unsupported quadratic power integral");
    }

    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression u =
        make_add(x, SymbolicExpression::number(b / (2.0 * a))).simplify();
    const Scalar d = delta / (4.0 * a * a);
    const SymbolicExpression shifted_quadratic =
        make_add(make_power(u, SymbolicExpression::number(2.0)),
                 SymbolicExpression::number(d))
            .simplify();

    SymbolicExpression integral =
        make_multiply(SymbolicExpression::number(1.0L / mymath::sqrt(d)),
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
                             (2 * n - 3) /
                             (2 * (n - 1)) / d),
                         integral))
                .simplify();
    }

    return make_multiply(SymbolicExpression::number(1.0L / mymath::pow(a, power)),
                         integral)
        .simplify();
}

/**
 * @brief 积分二次式部分分式项
 *
 * 处理 (slope*x + constant) / (ax^2 + bx + c)^power 形式的积分。
 * 分解为导数部分和逆二次式部分分别处理。
 */
SymbolicExpression integrate_quadratic_partial_fraction_term(
    const std::vector<Scalar>& quadratic,
    Scalar slope,
    Scalar constant,
    int power,
    const std::string& variable_name) {
    const Scalar a = quadratic[2];
    const Scalar b = quadratic[1];
    const Scalar derivative_scale = slope / (2.0 * a);
    const Scalar inverse_scale = constant - derivative_scale * b;
    const SymbolicExpression quadratic_expression =
        build_polynomial_expression_from_coefficients(quadratic, variable_name)
            .simplify();
    SymbolicExpression result = SymbolicExpression::number(0.0L);

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
            Scalar p_val = (power);
            derivative_part =
                make_multiply(make_divide(clean_derivative_scale, SymbolicExpression::number(1.0L - p_val)),
                              make_power(quadratic_expression,
                                         SymbolicExpression::number(1.0L - p_val)))
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

// ============================================================================
// 特殊有理式积分
// ============================================================================

/**
 * @brief 尝试积分 1/(1+x^2)^2 形式
 *
 * 特殊处理 (1+x^2)^n 分母的有理式。
 */
bool try_integrate_repeated_unit_quadratic(const std::vector<Scalar>& numerator,
                                           const std::vector<Scalar>& denominator,
                                           const std::string& variable_name,
                                           SymbolicExpression* integrated) {
    if (numerator.size() != 1 ||
        !mymath::is_near_zero(numerator[0] - 1.0L, kFormatEps) ||
        denominator.size() != 5 ||
        !mymath::is_near_zero(denominator[0] - 1.0L, kFormatEps) ||
        !mymath::is_near_zero(denominator[1], kFormatEps) ||
        !mymath::is_near_zero(denominator[2] - 2.0, kFormatEps) ||
        !mymath::is_near_zero(denominator[3], kFormatEps) ||
        !mymath::is_near_zero(denominator[4] - 1.0L, kFormatEps)) {
        return false;
    }

    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression one_plus_x_squared =
        make_add(SymbolicExpression::number(1.0L),
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

// ============================================================================
// 积分辅助函数
// ============================================================================

/**
 * @brief 检查两个表达式的简化形式是否相同
 */
bool same_simplified_expression(const SymbolicExpression& lhs,
                                const SymbolicExpression& rhs) {
    return lhs.simplify().to_string() == rhs.simplify().to_string();
}

/**
 * @brief 获取函数的原函数
 *
 * 返回常见函数的不定积分：
 * - exp(x) → exp(x)
 * - sin(x) → -cos(x)
 * - cos(x) → sin(x)
 * - tan(x) → -ln|cos(x)|
 * - ln(x) → x*ln(x) - x
 * - asin(x), acos(x), atan(x) 等反三角函数
 * - sinh(x), cosh(x), tanh(x) 双曲函数
 */
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
                              make_function("sqrt", make_subtract(SymbolicExpression::number(1.0L),
                                                                make_power(argument, SymbolicExpression::number(2.0)))))
                        .simplify();
        return true;
    }
    if (function_name == "acos") {
        *primitive = make_subtract(make_multiply(argument, make_function("acos", argument)),
                                   make_function("sqrt", make_subtract(SymbolicExpression::number(1.0L),
                                                                     make_power(argument, SymbolicExpression::number(2.0)))))
                        .simplify();
        return true;
    }
    if (function_name == "atan") {
        *primitive = make_subtract(make_multiply(argument, make_function("atan", argument)),
                                   make_multiply(SymbolicExpression::number(0.5),
                                               make_function("ln", make_add(SymbolicExpression::number(1.0L),
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

/**
 * @brief 尝试使用换元积分法处理乘积
 *
 * 检测 f(x) * g'(x) 形式的乘积，其中 g' 是某函数的导数。
 * 例如：x * exp(x^2) → 0.5 * exp(x^2)
 */
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
    Scalar scale = 1.0L;
    if (!same_simplified_expression(derivative_factor, expected_derivative)) {
        Scalar constant = 0.0L;
        SymbolicExpression rest;
        if (decompose_constant_times_expression(expected_derivative,
                                                variable_name,
                                                &constant,
                                                &rest) &&
            !mymath::is_near_zero(constant, kFormatEps) &&
            same_simplified_expression(derivative_factor, rest)) {
            scale = 1.0L / constant;
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

/**
 * @brief 尝试使用三角恒等式积分幂次形式
 *
 * 处理三角函数的幂次积分：
 * - sin^2(x) → x/2 - sin(2x)/4
 * - cos^2(x) → x/2 + sin(2x)/4
 * - tan^2(x) → tan(x)/a - x
 * - sin^3(x), cos^3(x) 等高次幂
 */
bool try_integrate_trig_power_identity(const SymbolicExpression& base,
                                       Scalar exponent_value,
                                       const std::string& variable_name,
                                       SymbolicExpression* integrated) {
    if (base.node_->type != NodeType::kFunction) {
        return false;
    }

    const SymbolicExpression argument(base.node_->left);
    Scalar a = 0.0L;
    Scalar b = 0.0L;
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

/**
 * @brief 尝试分部积分法
 *
 * 使用 LIATE 规则选择 u 和 dv：
 * L: 对数函数 (ln)
 * I: 反三角函数 (asin, atan)
 * A: 代数函数 (多项式)
 * T: 三角函数 (sin, cos)
 * E: 指数函数 (exp)
 *
 * 支持循环积分的代数求解（如 exp(x)*sin(x)）。
 */
bool try_integrate_by_parts(const SymbolicExpression& left,
                            const SymbolicExpression& right,
                            const std::string& variable_name,
                            SymbolicExpression* integrated) {
    struct IbpState {
        std::string original_key;
        SymbolicExpression uv_sum = SymbolicExpression::number(0.0L);
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
    Scalar k = 0.0L;
    SymbolicExpression remainder;
    if (decompose_constant_times_expression(v_du, variable_name, &k, &remainder) &&
        node_structural_key(remainder.node_) == original_key) {
        if (!mymath::is_near_zero(1.0L + k, kFormatEps)) {
            *integrated = make_divide(make_multiply(u, v),
                                      SymbolicExpression::number(1.0L + k)).simplify();
            return true;
        }
    }

    if (ibp_stack.size() >= 4) return false;
    ibp_stack.push_back({original_key, make_multiply(u, v)});
    try {
        *integrated = make_subtract(make_multiply(u, v),
                                    v_du.integral(variable_name))
                          .simplify();
        ibp_stack.pop_back();
        return true;
    } catch (...) {
        ibp_stack.pop_back();
    }

    // 如果没有直接发现循环，尝试递归积分 v_du
    if (ibp_stack.size() >= 4) return false;

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
                Scalar k2 = 0.0L;
                SymbolicExpression remainder2;
                if (decompose_constant_times_expression(v2_du2, variable_name, &k2, &remainder2) &&
                    node_structural_key(remainder2.node_) == original_key) {
                    // I = uv - (u2v2 - k2*I) => I = (uv - u2v2) / (1 - k2)
                    if (!mymath::is_near_zero(1.0L - k2, kFormatEps)) {
                        *integrated = make_divide(make_subtract(make_multiply(u, v),
                                                                make_multiply(u2, v2)),
                                                  SymbolicExpression::number(1.0L - k2)).simplify();
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

/**
 * @brief 尝试使用三角乘积恒等式积分
 *
 * 处理 sin * cos 乘积：
 * - sin(x) * cos(x) → sin^2(x) / 2
 */
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
    Scalar a = 0.0L;
    Scalar b = 0.0L;
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

/**
 * @brief 尝试积分 sec/csc 幂次与 tan/cot 的乘积
 *
 * 处理形如 sec^2(x) * tan(x) 的积分：
 * - sec(x) * tan(x) → sec(x)
 * - csc(x) * cot(x) → -csc(x)
 * - sec^2(x) * tan(x) → sec^2(x) / 2
 */
bool try_integrate_sec_csc_power_product(const SymbolicExpression& left,
                                         const SymbolicExpression& right,
                                         const std::string& variable_name,
                                         SymbolicExpression* integrated) {
    if (left.node_->type == NodeType::kFunction &&
        right.node_->type == NodeType::kFunction &&
        same_simplified_expression(SymbolicExpression(left.node_->left),
                                   SymbolicExpression(right.node_->left))) {
        const SymbolicExpression argument(left.node_->left);
        Scalar a = 0.0L;
        Scalar b = 0.0L;
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
    Scalar exponent_value = 0.0L;
    if (base.node_->type != NodeType::kFunction ||
        !exponent.is_number(&exponent_value) ||
        !mymath::is_near_zero(exponent_value - 2.0, kFormatEps) ||
        !same_simplified_expression(SymbolicExpression(base.node_->left),
                                    SymbolicExpression(function_factor->node_->left))) {
        return false;
    }

    const SymbolicExpression argument(base.node_->left);
    Scalar a = 0.0L;
    Scalar b = 0.0L;
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

// ============================================================================
// 有理式积分
// ============================================================================

/**
 * @brief 尝试积分多项式商
 *
 * 处理多项式有理式 P(x)/Q(x) 的积分：
 * 1. 如果分子次数 ≥ 分母次数，进行多项式除法
 * 2. 对真分式进行部分分式分解
 * 3. 特殊处理简单分母
 */
bool try_integrate_polynomial_quotient(const SymbolicExpression& numerator,
                                        const SymbolicExpression& denominator,
                                        const std::string& variable_name,
                                        SymbolicExpression* integrated) {
    std::vector<Scalar> num_coeffs;
    std::vector<Scalar> den_coeffs;
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

    if (den_coeffs.back() < 0.0L) {
        for (Scalar& c : num_coeffs) c = -c;
        for (Scalar& c : den_coeffs) c = -c;
    }

    SymbolicExpression poly_part = SymbolicExpression::number(0.0L);
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
        const Scalar constant = num_coeffs[0];
        const Scalar slope = den_coeffs[1];
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

bool is_symbolic_pure_monic_quadratic(const SymbolicExpression& expression,
                                      const std::string& variable_name,
                                      SymbolicExpression* constant_term) {
    std::vector<SymbolicExpression> coefficients;
    if (!symbolic_polynomial_coefficients_from_simplified(expression.simplify(),
                                                          variable_name,
                                                          &coefficients) ||
        coefficients.size() != 3 ||
        !expr_is_zero(coefficients[1]) ||
        !expr_is_one(coefficients[2])) {
        return false;
    }
    *constant_term = coefficients[0].simplify();
    return !constant_term->is_constant(variable_name) ? false : true;
}

bool collect_two_symbolic_linear_denominator_factors(
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* first_slope,
    SymbolicExpression* first_intercept,
    SymbolicExpression* second_slope,
    SymbolicExpression* second_intercept) {
    Scalar numeric_factor = 1.0L;
    std::vector<SymbolicExpression> factors;
    collect_multiplicative_terms(denominator.simplify(),
                                 &numeric_factor,
                                 &factors);
    if (!mymath::is_near_zero(numeric_factor - 1.0L, kFormatEps) ||
        factors.size() != 2) {
        return false;
    }
    return symbolic_decompose_linear(factors[0],
                                     variable_name,
                                     first_slope,
                                     first_intercept) &&
           !expr_is_zero(*first_slope) &&
           symbolic_decompose_linear(factors[1],
                                     variable_name,
                                     second_slope,
                                     second_intercept) &&
           !expr_is_zero(*second_slope);
}

bool collect_symbolic_linear_and_pure_quadratic_factors(
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* slope,
    SymbolicExpression* intercept,
    SymbolicExpression* quadratic_constant) {
    Scalar numeric_factor = 1.0L;
    std::vector<SymbolicExpression> factors;
    collect_multiplicative_terms(denominator.simplify(),
                                 &numeric_factor,
                                 &factors);
    if (!mymath::is_near_zero(numeric_factor - 1.0L, kFormatEps) ||
        factors.size() != 2) {
        return false;
    }

    auto try_order = [&](const SymbolicExpression& linear,
                         const SymbolicExpression& quadratic) {
        return symbolic_decompose_linear(linear,
                                         variable_name,
                                         slope,
                                         intercept) &&
               !expr_is_zero(*slope) &&
               is_symbolic_pure_monic_quadratic(quadratic,
                                                variable_name,
                                                quadratic_constant);
    };

    return try_order(factors[0], factors[1]) ||
           try_order(factors[1], factors[0]);
}

bool collect_symbolic_two_linear_and_pure_quadratic_factors(
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* first_slope,
    SymbolicExpression* first_intercept,
    SymbolicExpression* second_slope,
    SymbolicExpression* second_intercept,
    SymbolicExpression* quadratic_constant) {
    Scalar numeric_factor = 1.0L;
    std::vector<SymbolicExpression> factors;
    collect_multiplicative_terms(denominator.simplify(),
                                 &numeric_factor,
                                 &factors);
    if (!mymath::is_near_zero(numeric_factor - 1.0L, kFormatEps) ||
        factors.size() != 3) {
        return false;
    }

    for (std::size_t quadratic_index = 0; quadratic_index < factors.size(); ++quadratic_index) {
        if (!is_symbolic_pure_monic_quadratic(factors[quadratic_index],
                                              variable_name,
                                              quadratic_constant)) {
            continue;
        }
        std::vector<std::size_t> linear_indices;
        for (std::size_t i = 0; i < factors.size(); ++i) {
            if (i != quadratic_index) {
                linear_indices.push_back(i);
            }
        }
        if (linear_indices.size() != 2) {
            return false;
        }
        if (symbolic_decompose_linear(factors[linear_indices[0]],
                                      variable_name,
                                      first_slope,
                                      first_intercept) &&
            !expr_is_zero(*first_slope) &&
            symbolic_decompose_linear(factors[linear_indices[1]],
                                      variable_name,
                                      second_slope,
                                      second_intercept) &&
            !expr_is_zero(*second_slope)) {
            return true;
        }
    }
    return false;
}

bool collect_symbolic_repeated_linear_and_linear_factors(
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* repeated_slope,
    SymbolicExpression* repeated_intercept,
    SymbolicExpression* simple_slope,
    SymbolicExpression* simple_intercept) {
    Scalar numeric_factor = 1.0L;
    std::vector<SymbolicExpression> factors;
    collect_multiplicative_terms(denominator.simplify(),
                                 &numeric_factor,
                                 &factors);
    if (!mymath::is_near_zero(numeric_factor - 1.0L, kFormatEps) ||
        factors.size() != 2) {
        return false;
    }

    auto try_order = [&](const SymbolicExpression& repeated_factor,
                         const SymbolicExpression& simple_factor) {
        Scalar exponent_value = 0.0L;
        if (repeated_factor.node_->type != NodeType::kPower ||
            !SymbolicExpression(repeated_factor.node_->right).is_number(&exponent_value) ||
            !mymath::is_near_zero(exponent_value - 2.0, kFormatEps)) {
            return false;
        }
        return symbolic_decompose_linear(SymbolicExpression(repeated_factor.node_->left),
                                         variable_name,
                                         repeated_slope,
                                         repeated_intercept) &&
               !expr_is_zero(*repeated_slope) &&
               symbolic_decompose_linear(simple_factor,
                                         variable_name,
                                         simple_slope,
                                         simple_intercept) &&
               !expr_is_zero(*simple_slope);
    };

    return try_order(factors[0], factors[1]) ||
           try_order(factors[1], factors[0]);
}

bool collect_symbolic_repeated_linear_and_pure_quadratic_factors(
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* repeated_slope,
    SymbolicExpression* repeated_intercept,
    SymbolicExpression* quadratic_constant) {
    Scalar numeric_factor = 1.0L;
    std::vector<SymbolicExpression> factors;
    collect_multiplicative_terms(denominator.simplify(),
                                 &numeric_factor,
                                 &factors);
    if (!mymath::is_near_zero(numeric_factor - 1.0L, kFormatEps) ||
        factors.size() != 2) {
        return false;
    }

    auto try_order = [&](const SymbolicExpression& repeated_factor,
                         const SymbolicExpression& quadratic_factor) {
        Scalar exponent_value = 0.0L;
        if (repeated_factor.node_->type != NodeType::kPower ||
            !SymbolicExpression(repeated_factor.node_->right).is_number(&exponent_value) ||
            !mymath::is_near_zero(exponent_value - 2.0, kFormatEps)) {
            return false;
        }
        return symbolic_decompose_linear(SymbolicExpression(repeated_factor.node_->left),
                                         variable_name,
                                         repeated_slope,
                                         repeated_intercept) &&
               !expr_is_zero(*repeated_slope) &&
               is_symbolic_pure_monic_quadratic(quadratic_factor,
                                                variable_name,
                                                quadratic_constant);
    };

    return try_order(factors[0], factors[1]) ||
           try_order(factors[1], factors[0]);
}

bool collect_symbolic_repeated_pure_quadratic_factor(
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* constant_term) {
    const SymbolicExpression simplified = denominator.simplify();
    Scalar exponent_value = 0.0L;
    return simplified.node_->type == NodeType::kPower &&
           SymbolicExpression(simplified.node_->right).is_number(&exponent_value) &&
           mymath::is_near_zero(exponent_value - 2.0, kFormatEps) &&
           is_symbolic_pure_monic_quadratic(SymbolicExpression(simplified.node_->left),
                                            variable_name,
                                            constant_term);
}

bool symbolic_numerator_coefficients_up_to(
    const SymbolicExpression& numerator,
    const std::string& variable_name,
    std::size_t max_degree,
    std::vector<SymbolicExpression>* coefficients) {
    if (!symbolic_polynomial_coefficients_from_simplified(numerator.simplify(),
                                                          variable_name,
                                                          coefficients)) {
        return false;
    }
    trim_symbolic_polynomial_coefficients(coefficients);
    if (coefficients->size() > max_degree + 1) {
        return false;
    }
    while (coefficients->size() < max_degree + 1) {
        coefficients->push_back(SymbolicExpression::number(0.0L));
    }
    return true;
}

SymbolicExpression integrate_inverse_symbolic_pure_quadratic(
    const SymbolicExpression& constant_term,
    const std::string& variable_name);

SymbolicExpression symbolic_linear_value_at_root(
    const SymbolicExpression& slope,
    const SymbolicExpression& intercept,
    const SymbolicExpression& root) {
    return make_add(make_multiply(slope, root), intercept).simplify();
}

SymbolicExpression symbolic_quadratic_value_at(
    const SymbolicExpression& constant_term,
    const SymbolicExpression& root) {
    return make_add(make_power(root, SymbolicExpression::number(2.0)),
                    constant_term)
        .simplify();
}

bool try_integrate_symbolic_two_linear_factors(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    SymbolicExpression a1, b1, a2, b2;
    if (!collect_two_symbolic_linear_denominator_factors(denominator,
                                                         variable_name,
                                                         &a1,
                                                         &b1,
                                                         &a2,
                                                         &b2)) {
        return false;
    }

    std::vector<SymbolicExpression> numerator_coefficients;
    if (!symbolic_numerator_coefficients_up_to(numerator,
                                               variable_name,
                                               1,
                                               &numerator_coefficients)) {
        return false;
    }
    const SymbolicExpression q = numerator_coefficients[0].simplify();
    const SymbolicExpression p = numerator_coefficients[1].simplify();
    const SymbolicExpression determinant =
        make_subtract(make_multiply(a2, b1),
                      make_multiply(a1, b2))
            .simplify();
    if (expr_is_zero(determinant)) {
        return false;
    }

    const SymbolicExpression first_coefficient =
        make_divide(
            make_subtract(make_multiply(p, b1),
                          make_multiply(a1, q)),
            determinant)
            .simplify();
    const SymbolicExpression second_coefficient =
        make_divide(
            make_subtract(make_multiply(a2, q),
                          make_multiply(p, b2)),
            determinant)
            .simplify();
    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression first_linear =
        make_add(make_multiply(a1, x), b1).simplify();
    const SymbolicExpression second_linear =
        make_add(make_multiply(a2, x), b2).simplify();

    *integrated =
        make_add(
            make_multiply(make_divide(first_coefficient, a1),
                          make_function("ln", make_function("abs", first_linear))),
            make_multiply(make_divide(second_coefficient, a2),
                          make_function("ln", make_function("abs", second_linear))))
            .simplify();
    return true;
}

bool try_integrate_symbolic_repeated_linear_and_linear(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    SymbolicExpression a1, b1, a2, b2;
    if (!collect_symbolic_repeated_linear_and_linear_factors(denominator,
                                                             variable_name,
                                                             &a1,
                                                             &b1,
                                                             &a2,
                                                             &b2)) {
        return false;
    }

    std::vector<SymbolicExpression> numerator_coefficients;
    if (!symbolic_numerator_coefficients_up_to(numerator,
                                               variable_name,
                                               2,
                                               &numerator_coefficients)) {
        return false;
    }

    const SymbolicExpression repeated_root =
        make_divide(make_negate(b1), a1).simplify();
    const SymbolicExpression simple_root =
        make_divide(make_negate(b2), a2).simplify();
    const SymbolicExpression repeated_linear =
        make_add(make_multiply(a1, SymbolicExpression::variable(variable_name)), b1).simplify();
    const SymbolicExpression simple_linear =
        make_add(make_multiply(a2, SymbolicExpression::variable(variable_name)), b2).simplify();
    const SymbolicExpression simple_at_repeated =
        symbolic_linear_value_at_root(a2, b2, repeated_root);
    const SymbolicExpression repeated_at_simple =
        symbolic_linear_value_at_root(a1, b1, simple_root);
    if (expr_is_zero(simple_at_repeated) || expr_is_zero(repeated_at_simple)) {
        return false;
    }

    const SymbolicExpression numerator_at_repeated =
        numerator.substitute(variable_name, repeated_root).simplify();
    const SymbolicExpression numerator_at_simple =
        numerator.substitute(variable_name, simple_root).simplify();
    const SymbolicExpression numerator_prime_at_repeated =
        numerator.derivative(variable_name)
            .simplify()
            .substitute(variable_name, repeated_root)
            .simplify();

    const SymbolicExpression B =
        make_divide(numerator_at_repeated, simple_at_repeated).simplify();
    const SymbolicExpression C =
        make_divide(numerator_at_simple,
                    make_power(repeated_at_simple, SymbolicExpression::number(2.0)))
            .simplify();
    const SymbolicExpression A =
        make_divide(
            make_subtract(numerator_prime_at_repeated,
                          make_multiply(B, a2)),
            make_multiply(a1, simple_at_repeated))
            .simplify();

    *integrated =
        make_add(
            make_add(
                make_multiply(make_divide(A, a1),
                              make_function("ln", make_function("abs", repeated_linear))),
                make_multiply(make_negate(make_divide(B, a1)),
                              make_divide(SymbolicExpression::number(1.0L), repeated_linear))),
            make_multiply(make_divide(C, a2),
                          make_function("ln", make_function("abs", simple_linear))))
            .simplify();
    return true;
}

bool try_integrate_symbolic_two_linear_times_pure_quadratic(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    SymbolicExpression a1, b1, a2, b2, c;
    if (!collect_symbolic_two_linear_and_pure_quadratic_factors(denominator,
                                                                variable_name,
                                                                &a1,
                                                                &b1,
                                                                &a2,
                                                                &b2,
                                                                &c)) {
        return false;
    }

    std::vector<SymbolicExpression> numerator_coefficients;
    if (!symbolic_numerator_coefficients_up_to(numerator,
                                               variable_name,
                                               3,
                                               &numerator_coefficients)) {
        return false;
    }

    const SymbolicExpression first_root =
        make_divide(make_negate(b1), a1).simplify();
    const SymbolicExpression second_root =
        make_divide(make_negate(b2), a2).simplify();
    const SymbolicExpression first_linear =
        make_add(make_multiply(a1, SymbolicExpression::variable(variable_name)), b1).simplify();
    const SymbolicExpression second_linear =
        make_add(make_multiply(a2, SymbolicExpression::variable(variable_name)), b2).simplify();
    const SymbolicExpression simple_at_first =
        symbolic_linear_value_at_root(a2, b2, first_root);
    const SymbolicExpression simple_at_second =
        symbolic_linear_value_at_root(a1, b1, second_root);
    const SymbolicExpression quadratic_at_first =
        symbolic_quadratic_value_at(c, first_root);
    const SymbolicExpression quadratic_at_second =
        symbolic_quadratic_value_at(c, second_root);
    if (expr_is_zero(simple_at_first) || expr_is_zero(simple_at_second) ||
        expr_is_zero(quadratic_at_first) || expr_is_zero(quadratic_at_second)) {
        return false;
    }

    const SymbolicExpression A =
        make_divide(numerator.substitute(variable_name, first_root).simplify(),
                    make_multiply(simple_at_first, quadratic_at_first))
            .simplify();
    const SymbolicExpression B =
        make_divide(numerator.substitute(variable_name, second_root).simplify(),
                    make_multiply(simple_at_second, quadratic_at_second))
            .simplify();

    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression quadratic =
        make_add(make_power(x, SymbolicExpression::number(2.0)), c).simplify();
    const SymbolicExpression remaining =
        make_subtract(
            make_subtract(numerator,
                          make_multiply(A,
                                        make_multiply(second_linear, quadratic))),
            make_multiply(B,
                          make_multiply(first_linear, quadratic)))
            .simplify();
    std::vector<SymbolicExpression> remaining_coefficients;
    if (!symbolic_numerator_coefficients_up_to(remaining,
                                               variable_name,
                                               3,
                                               &remaining_coefficients)) {
        return false;
    }

    const SymbolicExpression linear_pair_x =
        make_multiply(a1, a2).simplify();
    if (expr_is_zero(linear_pair_x)) {
        return false;
    }
    const SymbolicExpression linear_pair_mid =
        make_add(make_multiply(a1, b2),
                 make_multiply(a2, b1))
            .simplify();
    const SymbolicExpression C =
        make_divide(remaining_coefficients[3], linear_pair_x).simplify();
    const SymbolicExpression D =
        make_divide(make_subtract(remaining_coefficients[2],
                                  make_multiply(C, linear_pair_mid)),
                    linear_pair_x)
            .simplify();

    *integrated =
        make_add(
            make_add(
                make_add(
                    make_multiply(make_divide(A, a1),
                                  make_function("ln", make_function("abs", first_linear))),
                    make_multiply(make_divide(B, a2),
                                  make_function("ln", make_function("abs", second_linear)))),
                make_multiply(make_divide(C, SymbolicExpression::number(2.0)),
                              make_function("ln", make_function("abs", quadratic)))),
            make_multiply(D, integrate_inverse_symbolic_pure_quadratic(c, variable_name)))
            .simplify();
    return true;
}

bool try_integrate_symbolic_linear_times_pure_quadratic(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    SymbolicExpression a, b, c;
    if (!collect_symbolic_linear_and_pure_quadratic_factors(denominator,
                                                            variable_name,
                                                            &a,
                                                            &b,
                                                            &c)) {
        return false;
    }

    std::vector<SymbolicExpression> numerator_coefficients;
    if (!symbolic_numerator_coefficients_up_to(numerator,
                                               variable_name,
                                               2,
                                               &numerator_coefficients)) {
        return false;
    }
    const SymbolicExpression n0 = numerator_coefficients[0].simplify();
    const SymbolicExpression n1 = numerator_coefficients[1].simplify();
    const SymbolicExpression n2 = numerator_coefficients[2].simplify();
    const SymbolicExpression denominator_scale =
        make_add(make_multiply(make_multiply(a, a), c),
                 make_multiply(b, b))
            .simplify();
    if (expr_is_zero(denominator_scale)) {
        return false;
    }

    const SymbolicExpression B =
        make_divide(
            make_subtract(
                make_add(make_multiply(make_multiply(a, c), n2),
                         make_multiply(b, n1)),
                make_multiply(a, n0)),
            denominator_scale)
            .simplify();
    const SymbolicExpression C =
        make_divide(make_subtract(n1, make_multiply(b, B)), a).simplify();
    const SymbolicExpression A =
        make_subtract(n2, make_multiply(a, B)).simplify();

    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression linear =
        make_add(make_multiply(a, x), b).simplify();
    const SymbolicExpression quadratic =
        make_add(make_power(x, SymbolicExpression::number(2.0)), c).simplify();
    const SymbolicExpression inverse_quadratic_integral =
        integrate_inverse_symbolic_pure_quadratic(c, variable_name);

    *integrated =
        make_add(
            make_add(
                make_multiply(make_divide(A, a),
                              make_function("ln", make_function("abs", linear))),
                make_multiply(make_divide(B, SymbolicExpression::number(2.0)),
                              make_function("ln", make_function("abs", quadratic)))),
            make_multiply(C, inverse_quadratic_integral))
            .simplify();
    return true;
}

bool try_integrate_symbolic_repeated_linear_times_pure_quadratic(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    SymbolicExpression a, b, c;
    if (!collect_symbolic_repeated_linear_and_pure_quadratic_factors(denominator,
                                                                     variable_name,
                                                                     &a,
                                                                     &b,
                                                                     &c)) {
        return false;
    }

    std::vector<SymbolicExpression> numerator_coefficients;
    if (!symbolic_numerator_coefficients_up_to(numerator,
                                               variable_name,
                                               3,
                                               &numerator_coefficients)) {
        return false;
    }

    const SymbolicExpression repeated_root =
        make_divide(make_negate(b), a).simplify();
    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression linear =
        make_add(make_multiply(a, x), b).simplify();
    const SymbolicExpression quadratic =
        make_add(make_power(x, SymbolicExpression::number(2.0)), c).simplify();
    const SymbolicExpression quadratic_at_root =
        symbolic_quadratic_value_at(c, repeated_root);
    if (expr_is_zero(quadratic_at_root)) {
        return false;
    }

    const SymbolicExpression numerator_at_root =
        numerator.substitute(variable_name, repeated_root).simplify();
    const SymbolicExpression numerator_prime_at_root =
        numerator.derivative(variable_name)
            .simplify()
            .substitute(variable_name, repeated_root)
            .simplify();
    const SymbolicExpression B =
        make_divide(numerator_at_root, quadratic_at_root).simplify();
    const SymbolicExpression quadratic_prime_at_root =
        make_multiply(SymbolicExpression::number(2.0), repeated_root).simplify();
    const SymbolicExpression A =
        make_divide(make_subtract(numerator_prime_at_root,
                                  make_multiply(B, quadratic_prime_at_root)),
                    make_multiply(a, quadratic_at_root))
            .simplify();

    const SymbolicExpression remaining =
        make_subtract(
            make_subtract(numerator,
                          make_multiply(A,
                                        make_multiply(linear, quadratic))),
            make_multiply(B, quadratic))
            .simplify();
    std::vector<SymbolicExpression> remaining_coefficients;
    if (!symbolic_numerator_coefficients_up_to(remaining,
                                               variable_name,
                                               3,
                                               &remaining_coefficients)) {
        return false;
    }

    const SymbolicExpression repeated_square_x2 =
        make_multiply(a, a).simplify();
    if (expr_is_zero(repeated_square_x2)) {
        return false;
    }
    const SymbolicExpression repeated_square_x =
        make_multiply(SymbolicExpression::number(2.0),
                      make_multiply(a, b))
            .simplify();
    const SymbolicExpression C =
        make_divide(remaining_coefficients[3], repeated_square_x2).simplify();
    const SymbolicExpression D =
        make_divide(make_subtract(remaining_coefficients[2],
                                  make_multiply(C, repeated_square_x)),
                    repeated_square_x2)
            .simplify();

    *integrated =
        make_add(
            make_add(
                make_add(
                    make_multiply(make_divide(A, a),
                                  make_function("ln", make_function("abs", linear))),
                    make_multiply(make_negate(make_divide(B, a)),
                                  make_divide(SymbolicExpression::number(1.0L), linear))),
                make_multiply(make_divide(C, SymbolicExpression::number(2.0)),
                              make_function("ln", make_function("abs", quadratic)))),
            make_multiply(D, integrate_inverse_symbolic_pure_quadratic(c, variable_name)))
            .simplify();
    return true;
}

bool try_integrate_symbolic_repeated_pure_quadratic(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    SymbolicExpression constant_term;
    if (!collect_symbolic_repeated_pure_quadratic_factor(denominator,
                                                         variable_name,
                                                         &constant_term)) {
        return false;
    }

    std::vector<SymbolicExpression> numerator_coefficients;
    if (!symbolic_numerator_coefficients_up_to(numerator,
                                               variable_name,
                                               1,
                                               &numerator_coefficients)) {
        return false;
    }
    const SymbolicExpression c0 = numerator_coefficients[0].simplify();
    const SymbolicExpression c1 = numerator_coefficients[1].simplify();
    const SymbolicExpression x = SymbolicExpression::variable(variable_name);
    const SymbolicExpression quadratic =
        make_add(make_power(x, SymbolicExpression::number(2.0)), constant_term).simplify();
    const SymbolicExpression inverse_quadratic_integral =
        integrate_inverse_symbolic_pure_quadratic(constant_term, variable_name);

    if (!expr_is_zero(c0)) {
        const SymbolicExpression first_term =
            make_multiply(
                c0,
                make_divide(x,
                            make_multiply(
                                make_multiply(SymbolicExpression::number(2.0), constant_term),
                                quadratic)))
                .simplify();
        const SymbolicExpression second_term =
            make_multiply(
                c0,
                make_divide(inverse_quadratic_integral,
                            make_multiply(SymbolicExpression::number(2.0), constant_term)))
                .simplify();
        *integrated = make_add(first_term, second_term).simplify();
        if (!expr_is_zero(c1)) {
            *integrated =
                make_add(
                    *integrated,
                    make_multiply(
                        c1,
                        make_negate(
                            make_divide(SymbolicExpression::number(1.0L),
                                        make_multiply(SymbolicExpression::number(2.0),
                                                      quadratic)))))
                    .simplify();
        }
        return true;
    }

    if (!expr_is_zero(c1)) {
        *integrated =
            make_multiply(
                c1,
                make_negate(
                    make_divide(SymbolicExpression::number(1.0L),
                                make_multiply(SymbolicExpression::number(2.0),
                                              quadratic))))
                .simplify();
        return true;
    }

    return false;
}

bool collect_two_pure_quadratic_denominator_factors(
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* first_constant,
    SymbolicExpression* second_constant) {
    Scalar numeric_factor = 1.0L;
    std::vector<SymbolicExpression> factors;
    collect_multiplicative_terms(denominator.simplify(),
                                 &numeric_factor,
                                 &factors);
    if (!mymath::is_near_zero(numeric_factor - 1.0L, kFormatEps) ||
        factors.size() != 2) {
        return false;
    }
    return is_symbolic_pure_monic_quadratic(factors[0],
                                            variable_name,
                                            first_constant) &&
           is_symbolic_pure_monic_quadratic(factors[1],
                                            variable_name,
                                            second_constant);
}

SymbolicExpression integrate_inverse_symbolic_pure_quadratic(
    const SymbolicExpression& constant_term,
    const std::string& variable_name) {
    const SymbolicExpression root =
        make_function("sqrt", constant_term).simplify();
    return make_divide(
               make_function(
                   "atan",
                   make_divide(SymbolicExpression::variable(variable_name), root)),
               root)
        .simplify();
}

bool try_integrate_symbolic_two_pure_quadratics(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& variable_name,
    SymbolicExpression* integrated) {
    SymbolicExpression first_constant;
    SymbolicExpression second_constant;
    if (!collect_two_pure_quadratic_denominator_factors(denominator,
                                                        variable_name,
                                                        &first_constant,
                                                        &second_constant)) {
        return false;
    }

    const SymbolicExpression difference =
        make_subtract(second_constant, first_constant).simplify();
    if (expr_is_zero(difference)) {
        return false;
    }

    std::vector<SymbolicExpression> numerator_coefficients;
    if (!symbolic_polynomial_coefficients_from_simplified(numerator.simplify(),
                                                          variable_name,
                                                          &numerator_coefficients)) {
        return false;
    }
    trim_symbolic_polynomial_coefficients(&numerator_coefficients);

    if (numerator_coefficients.size() == 1) {
        const SymbolicExpression coefficient = numerator_coefficients[0].simplify();
        const SymbolicExpression first_integral =
            integrate_inverse_symbolic_pure_quadratic(first_constant,
                                                      variable_name);
        const SymbolicExpression second_integral =
            integrate_inverse_symbolic_pure_quadratic(second_constant,
                                                      variable_name);
        *integrated =
            make_multiply(
                coefficient,
                make_divide(make_subtract(first_integral, second_integral),
                            difference))
                .simplify();
        return true;
    }

    if (numerator_coefficients.size() == 2 &&
        expr_is_zero(numerator_coefficients[0])) {
        const SymbolicExpression coefficient = numerator_coefficients[1].simplify();
        const SymbolicExpression x = SymbolicExpression::variable(variable_name);
        const SymbolicExpression first_quadratic =
            make_add(make_power(x, SymbolicExpression::number(2.0)),
                     first_constant)
                .simplify();
        const SymbolicExpression second_quadratic =
            make_add(make_power(x, SymbolicExpression::number(2.0)),
                     second_constant)
                .simplify();
        *integrated =
            make_multiply(
                coefficient,
                make_divide(
                    make_subtract(
                        make_function("ln", make_function("abs", first_quadratic)),
                        make_function("ln", make_function("abs", second_quadratic))),
                    make_multiply(SymbolicExpression::number(2.0), difference)))
                .simplify();
        return true;
    }

    return false;
}

// ============================================================================
// Weierstrass 置换
// ============================================================================

/**
 * @brief 尝试 Weierstrass 置换积分
 *
 * 对于含三角函数的分式积分，使用万能代换 t = tan(x/2)：
 * - sin(x) = 2t / (1+t^2)
 * - cos(x) = (1-t^2) / (1+t^2)
 * - dx = 2 / (1+t^2) dt
 *
 * 将三角积分转化为有理式积分。
 */
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

    Scalar A = 0.0L;
    SymbolicExpression B_expr = SymbolicExpression::number(0.0L);
    SymbolicExpression C_expr = SymbolicExpression::number(0.0L);
    SymbolicExpression argument;
    bool found_argument = false;

    for (const auto& term : terms) {
        Scalar val = 0.0L;
        if (term.is_number(&val)) {
            A += val;
        } else if (term.node_->type == NodeType::kFunction && term.node_->text == "sin") {
            if (!found_argument) {
                argument = SymbolicExpression(term.node_->left);
                found_argument = true;
            } else if (!expressions_match(argument, SymbolicExpression(term.node_->left))) {
                return false;
            }
            B_expr = make_add(B_expr, SymbolicExpression::number(1.0L)).simplify();
        } else if (term.node_->type == NodeType::kFunction && term.node_->text == "cos") {
            if (!found_argument) {
                argument = SymbolicExpression(term.node_->left);
                found_argument = true;
            } else if (!expressions_match(argument, SymbolicExpression(term.node_->left))) {
                return false;
            }
            C_expr = make_add(C_expr, SymbolicExpression::number(1.0L)).simplify();
        } else if (term.node_->type == NodeType::kMultiply) {
            // 处理 B * sin(x) 或 C * cos(x)
            Scalar coeff = 0.0L;
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

    Scalar a = 0.0L, b = 0.0L;
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
    Scalar b_val = 0.0L, c_val = 0.0L;
    if (!B_expr.is_number(&b_val) || !C_expr.is_number(&c_val)) return false;

    std::vector<Scalar> poly_numerator = { 2.0 / a };
    std::vector<Scalar> poly_denominator = { A + c_val, 2.0 * b_val, A - c_val };

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

}  // namespace symbolic_expression_internal

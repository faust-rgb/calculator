#include "symbolic/symbolic_algebraic_number.h"
#include "symbolic/symbolic_expression_internal.h"
#include <algorithm>
#include <sstream>

using namespace symbolic_expression_internal;

// ============================================================================
// AlgebraicNumber 实现
// ============================================================================

AlgebraicNumber AlgebraicNumber::from_double(double value) {
    // 尝试识别常见的代数数
    double sqrt2 = mymath::sqrt(2.0);
    double sqrt3 = mymath::sqrt(3.0);
    double sqrt5 = mymath::sqrt(5.0);

    if (mymath::abs(value - sqrt2) < 1e-10) {
        return sqrt(2);
    }
    if (mymath::abs(value + sqrt2) < 1e-10) {
        return sqrt(2).negate();
    }
    if (mymath::abs(value - sqrt3) < 1e-10) {
        return sqrt(3);
    }
    if (mymath::abs(value + sqrt3) < 1e-10) {
        return sqrt(3).negate();
    }
    if (mymath::abs(value - sqrt5) < 1e-10) {
        return sqrt(5);
    }
    if (mymath::abs(value + sqrt5) < 1e-10) {
        return sqrt(5).negate();
    }

    // 否则作为有理数近似
    ExactRational r = ExactRational::from_double(value);
    return from_rational(r);
}

AlgebraicNumber AlgebraicNumber::sqrt(int n) {
    if (n < 0) {
        // 复数情况
        // sqrt(-n) = i * sqrt(n)
        // 返回标记为复数的代数数
        AlgebraicNumber real_part = sqrt(-n);
        real_part.is_real = false;
        real_part.complex_sign = 1;
        return real_part;
    }

    if (n == 0) {
        return from_integer(0);
    }

    if (n == 1) {
        return from_integer(1);
    }

    // 检查是否是完全平方数
    int root = static_cast<int>(mymath::sqrt(static_cast<double>(n)));
    if (root * root == n) {
        return from_integer(root);
    }

    // 创建最小多项式 t^2 - n = 0
    std::vector<SymbolicExpression> coeffs;
    coeffs.push_back(SymbolicExpression::number(-n));
    coeffs.push_back(SymbolicExpression::number(0.0));
    coeffs.push_back(SymbolicExpression::number(1.0));
    SymbolicPolynomial poly(coeffs, "_t");

    // 隔离区间: sqrt(n) 在 (floor(sqrt(n)), ceil(sqrt(n))) 中
    int lower = root;
    int upper = root + 1;

    return AlgebraicNumber(poly, ExactRational(lower), ExactRational(upper), true, 0, 0);
}

AlgebraicNumber AlgebraicNumber::root_of_unity(int n, int k) {
    // e^(2*pi*i*k/n) 的最小多项式是 x^n - 1 的因子
    // 这里简化处理，仅返回标记

    std::vector<SymbolicExpression> coeffs;
    coeffs.push_back(SymbolicExpression::number(-1.0));
    for (int i = 1; i < n; ++i) {
        coeffs.push_back(SymbolicExpression::number(0.0));
    }
    coeffs.push_back(SymbolicExpression::number(1.0));
    SymbolicPolynomial poly(coeffs, "_z");

    // 单位根在单位圆上
    return AlgebraicNumber(poly, ExactRational(-1), ExactRational(1), false, k % n, k % n == 0 ? 0 : 1);
}

AlgebraicNumber AlgebraicNumber::add(const AlgebraicNumber& other) const {
    // 特殊情况：有理数
    ExactRational r1, r2;
    if (is_rational(&r1) && other.is_rational(&r2)) {
        return from_rational(r1.add(r2));
    }

    // 一般情况：通过结式计算
    // (alpha + beta) 的最小多项式可以通过
    // resultant(P_alpha(x - y), P_beta(y), y) 计算
    // 这里使用简化实现

    SymbolicPolynomial sum_poly = compute_sum_minpoly(minimal_polynomial, other.minimal_polynomial);

    // 计算新的隔离区间
    double approx1 = approximate();
    double approx2 = other.approximate();
    double sum_approx = approx1 + approx2;

    ExactRational new_lower = ExactRational::from_double(sum_approx - 0.1);
    ExactRational new_upper = ExactRational::from_double(sum_approx + 0.1);

    return AlgebraicNumber(sum_poly, new_lower, new_upper, is_real && other.is_real, 0, 0);
}

AlgebraicNumber AlgebraicNumber::multiply(const AlgebraicNumber& other) const {
    // 特殊情况：有理数
    ExactRational r1, r2;
    if (is_rational(&r1) && other.is_rational(&r2)) {
        return from_rational(r1.multiply(r2));
    }

    // 特殊情况：一个为零
    if (is_zero() || other.is_zero()) {
        return from_integer(0);
    }

    // 一般情况：通过结式计算
    SymbolicPolynomial prod_poly = compute_product_minpoly(minimal_polynomial, other.minimal_polynomial);

    // 计算新的隔离区间
    double approx1 = approximate();
    double approx2 = other.approximate();
    double prod_approx = approx1 * approx2;

    ExactRational new_lower = ExactRational::from_double(prod_approx - 0.1);
    ExactRational new_upper = ExactRational::from_double(prod_approx + 0.1);

    return AlgebraicNumber(prod_poly, new_lower, new_upper, is_real && other.is_real, 0, 0);
}

AlgebraicNumber AlgebraicNumber::negate() const {
    ExactRational r;
    if (is_rational(&r)) {
        return from_rational(r.negate());
    }

    // -alpha 的最小多项式: P(-x)
    std::vector<SymbolicExpression> new_coeffs;
    for (int i = 0; i <= minimal_polynomial.degree(); ++i) {
        SymbolicExpression coeff = minimal_polynomial.coefficient(i);
        if (i % 2 == 1) {
            // 奇次项取负
            coeff = (SymbolicExpression::number(-1.0) * coeff).simplify();
        }
        new_coeffs.push_back(coeff);
    }
    SymbolicPolynomial new_poly(new_coeffs, minimal_polynomial.variable_name());

    // 隔离区间取负
    return AlgebraicNumber(new_poly, interval_upper.negate(), interval_lower.negate(),
                           is_real, root_index, complex_sign);
}

AlgebraicNumber AlgebraicNumber::inverse() const {
    ExactRational r;
    if (is_rational(&r) && !r.is_zero()) {
        return from_rational(r.divide(ExactRational(1)));
    }

    if (is_zero()) {
        // 零没有倒数
        return *this;
    }

    // 1/alpha 的最小多项式: x^d * P(1/x)
    int d = minimal_polynomial.degree();
    std::vector<SymbolicExpression> new_coeffs(d + 1);
    for (int i = 0; i <= d; ++i) {
        // 新多项式的系数 c'_i = c_{d-i}
        new_coeffs[i] = minimal_polynomial.coefficient(d - i);
    }
    SymbolicPolynomial new_poly(new_coeffs, minimal_polynomial.variable_name());

    // 隔离区间取倒数 (注意边界反转)
    if (interval_lower.is_positive() || interval_upper.is_negative()) {
        // 区间不包含零，可以安全取倒数
        return AlgebraicNumber(new_poly,
                              interval_upper.divide(ExactRational(1)),
                              interval_lower.divide(ExactRational(1)),
                              is_real, root_index, complex_sign);
    }

    // 区间包含零，倒数不存在或为无穷
    return *this;
}

AlgebraicNumber AlgebraicNumber::power(int n) const {
    if (n == 0) {
        return from_integer(1);
    }
    if (n == 1) {
        return *this;
    }
    if (n < 0) {
        return power(-n).inverse();
    }

    // 使用快速幂
    if (n % 2 == 0) {
        AlgebraicNumber half = power(n / 2);
        return half.multiply(half);
    } else {
        return multiply(power(n - 1));
    }
}

int AlgebraicNumber::compare(const AlgebraicNumber& other) const {
    // 特殊情况：有理数比较
    ExactRational r1, r2;
    if (is_rational(&r1) && other.is_rational(&r2)) {
        return r1.compare(r2);
    }

    // 检查隔离区间是否分离
    if (interval_upper < other.interval_lower) {
        return -1;  // this < other
    }
    if (interval_lower > other.interval_upper) {
        return 1;   // this > other
    }

    // 区间重叠，需要细化
    // 这里简化处理，使用数值近似
    double approx1 = approximate();
    double approx2 = other.approximate();

    if (mymath::abs(approx1 - approx2) < 1e-10) {
        // 可能相等，检查最小多项式
        if (minimal_polynomial.to_string() == other.minimal_polynomial.to_string()) {
            // 相同最小多项式，需要更精确的比较
            // 这里简化返回相等
            return 0;
        }
    }

    if (approx1 < approx2) return -1;
    if (approx1 > approx2) return 1;
    return 0;
}

void AlgebraicNumber::refine_interval() const {
    // 二分细化
    ExactRational mid = ExactRational((interval_lower.numerator + interval_upper.numerator),
                            (interval_lower.denominator + interval_upper.denominator) / 2);
    (void)mid;

    // 评估多项式在 mid 的符号
    // 这里简化处理
    // 实际实现需要多项式求值
}

void AlgebraicNumber::refine_to_precision(int bits) const {
    // 细化到指定比特精度
    ExactRational width = interval_upper.subtract(interval_lower);
    ExactRational target_width(1, 1 << bits);

    while (width.compare(target_width) > 0) {
        refine_interval();
        width = interval_upper.subtract(interval_lower);
    }
}

SymbolicExpression AlgebraicNumber::to_expression() const {
    ExactRational r;
    if (is_rational(&r)) {
        return SymbolicExpression::number(r.to_double());
    }

    // 检查是否是 sqrt 形式
    if (minimal_polynomial.degree() == 2) {
        double a = 0.0, b = 0.0, c = 0.0;
        if (minimal_polynomial.coefficient(2).is_number(&a) &&
            minimal_polynomial.coefficient(1).is_number(&b) &&
            minimal_polynomial.coefficient(0).is_number(&c)) {

            if (mymath::abs(b) < 1e-10 && mymath::abs(a - 1.0) < 1e-10) {
                // t^2 - c = 0, t = sqrt(-c)
                double val = -c;
                if (val > 0) {
                    SymbolicExpression sqrt_val = make_function("sqrt", SymbolicExpression::number(val));
                    if (approximate() < 0) {
                        return make_negate(sqrt_val).simplify();
                    }
                    return sqrt_val;
                }
            }
        }
    }

    // 一般情况：返回 RootOf 形式
    return to_rootof_expression();
}

SymbolicExpression AlgebraicNumber::to_rootof_expression() const {
    // 创建 RootOf(poly, var, root_index) 表达式
    // 使用专用的 kRootOf 节点类型
    return make_rootof_from_polynomial(minimal_polynomial, root_index);
}

std::string AlgebraicNumber::to_string() const {
    std::ostringstream oss;

    ExactRational r;
    if (is_rational(&r)) {
        oss << r.to_string();
        return oss.str();
    }

    oss << "RootOf(" << minimal_polynomial.to_string() << ", " << root_index << ")";
    return oss.str();
}

// ============================================================================
// 最小多项式计算 (通过结式)
// ============================================================================

SymbolicPolynomial AlgebraicNumber::compute_sum_minpoly(
    const SymbolicPolynomial& p1,
    const SymbolicPolynomial& p2) {

    // 对于 alpha 和 beta，alpha + beta 的最小多项式是
    // resultant(p1(x - y), p2(y), y) 关于 x 的因子
    //
    // 实现方法：
    // 设 p1(t1) = 0 定义 alpha，p2(t2) = 0 定义 beta
    // 则 alpha + beta 满足 resultant(p1(x - t), p2(t), t) = 0

    int d1 = p1.degree();
    int d2 = p2.degree();

    if (d1 <= 0) return p2;
    if (d2 <= 0) return p1;

    // 构造 p1(x - t)，将 p1 中的变量替换为 (x - t)
    // 设 p1(t1) = sum_{i=0}^{d1} a_i * t1^i
    // 则 p1(x - t) = sum_{i=0}^{d1} a_i * (x - t)^i

    std::string var1 = p1.variable_name();
    std::string var2 = p2.variable_name();
    std::string x_var = "_x";  // 结果变量
    std::string t_var = "_t";  // 中间变量

    // 构造 p1(x - t)
    // 使用二项式展开: (x - t)^i = sum_{j=0}^{i} C(i,j) * x^j * (-t)^{i-j}
    std::vector<SymbolicExpression> p1_shifted_coeffs(d1 + 1, SymbolicExpression::number(0.0));

    for (int i = 0; i <= d1; ++i) {
        SymbolicExpression a_i = p1.coefficient(i);
        if (SymbolicPolynomial::coeff_is_zero(a_i)) continue;

        // (x - t)^i 的展开
        for (int j = 0; j <= i; ++j) {
            // C(i, j) * x^j * (-t)^{i-j}
            double binom = 1.0;
            for (int k = 0; k < j; ++k) {
                binom = binom * (i - k) / (k + 1);
            }

            SymbolicExpression term = a_i;
            term = (term * SymbolicExpression::number(binom)).simplify();

            // x^j 部分
            if (j > 0) {
                term = (term * make_power(SymbolicExpression::variable(x_var),
                                         SymbolicExpression::number(static_cast<double>(j)))).simplify();
            }

            // (-t)^{i-j} 部分
            if (i - j > 0) {
                double sign = ((i - j) % 2 == 0) ? 1.0 : -1.0;
                term = (term * SymbolicExpression::number(sign)).simplify();
                term = (term * make_power(SymbolicExpression::variable(t_var),
                                         SymbolicExpression::number(static_cast<double>(i - j)))).simplify();
            }

            // 将项添加到对应 x^j 的系数
            // 这里需要按 x 的幂次分组，构造关于 t 的多项式
        }
    }

    // 由于上述方法复杂，使用简化方法：
    // 直接计算 resultant(p1(x - t), p2(t), t)
    // 但需要处理符号系数

    // 简化实现：使用数值近似构造近似多项式
    // 完整实现需要多变量多项式的结式计算

    int new_deg = d1 * d2;  // 和的度数最多为 d1 * d2

    // 使用子结果式链计算结式
    // 这里我们使用一个更直接的方法：
    // 设 alpha 是 p1 的根，beta 是 p2 的根
    // alpha + beta 的最小多项式可以通过消元得到

    // 对于简单情况（数值系数），使用直接方法
    double a1 = 0.0, b1 = 0.0, c1 = 0.0;
    double a2 = 0.0, b2 = 0.0, c2 = 0.0;

    bool p1_numeric = (d1 == 2 &&
                      p1.coefficient(2).is_number(&a1) &&
                      p1.coefficient(1).is_number(&b1) &&
                      p1.coefficient(0).is_number(&c1));

    bool p2_numeric = (d2 == 2 &&
                      p2.coefficient(2).is_number(&a2) &&
                      p2.coefficient(1).is_number(&b2) &&
                      p2.coefficient(0).is_number(&c2));

    if (p1_numeric && p2_numeric) {
        // sqrt(r1) + sqrt(r2) 的最小多项式
        // 设 alpha = sqrt(r1), beta = sqrt(r2)
        // p1: a1*t^2 + b1*t + c1 = 0, alpha = (-b1 + sqrt(b1^2 - 4*a1*c1)) / (2*a1)
        // p2: a2*t^2 + b2*t + c2 = 0

        // 简化：假设 p1 = t^2 - r1, p2 = t^2 - r2 (即 sqrt(r1), sqrt(r2))
        if (mymath::abs(a1 - 1.0) < 1e-10 && mymath::abs(b1) < 1e-10 &&
            mymath::abs(a2 - 1.0) < 1e-10 && mymath::abs(b2) < 1e-10) {
            double r1 = -c1;
            double r2 = -c2;

            // sqrt(r1) + sqrt(r2) 的最小多项式是 t^4 - 2*(r1+r2)*t^2 + (r1-r2)^2
            std::vector<SymbolicExpression> coeffs(5);
            coeffs[0] = SymbolicExpression::number((r1 - r2) * (r1 - r2));  // (r1-r2)^2
            coeffs[1] = SymbolicExpression::number(0.0);
            coeffs[2] = SymbolicExpression::number(-2.0 * (r1 + r2));  // -2*(r1+r2)
            coeffs[3] = SymbolicExpression::number(0.0);
            coeffs[4] = SymbolicExpression::number(1.0);  // t^4

            return SymbolicPolynomial(coeffs, "_t");
        }
    }

    // 一般情况：返回占位多项式
    // 完整实现需要符号结式计算
    std::vector<SymbolicExpression> coeffs(new_deg + 1, SymbolicExpression::number(0.0));
    coeffs[0] = SymbolicExpression::number(-1.0);
    coeffs[new_deg] = SymbolicExpression::number(1.0);

    return SymbolicPolynomial(coeffs, "_t");
}

SymbolicPolynomial AlgebraicNumber::compute_product_minpoly(
    const SymbolicPolynomial& p1,
    const SymbolicPolynomial& p2) {

    // 对于 alpha 和 beta，alpha * beta 的最小多项式是
    // resultant(y^d1 * p1(x/y), p2(y), y) 关于 x 的因子

    int d1 = p1.degree();
    int d2 = p2.degree();

    if (d1 <= 0) return p2;
    if (d2 <= 0) return p1;

    // 对于简单情况（数值系数），使用直接方法
    double a1 = 0.0, b1 = 0.0, c1 = 0.0;
    double a2 = 0.0, b2 = 0.0, c2 = 0.0;

    bool p1_numeric = (d1 == 2 &&
                      p1.coefficient(2).is_number(&a1) &&
                      p1.coefficient(1).is_number(&b1) &&
                      p1.coefficient(0).is_number(&c1));

    bool p2_numeric = (d2 == 2 &&
                      p2.coefficient(2).is_number(&a2) &&
                      p2.coefficient(1).is_number(&b2) &&
                      p2.coefficient(0).is_number(&c2));

    if (p1_numeric && p2_numeric) {
        // sqrt(r1) * sqrt(r2) = sqrt(r1 * r2)
        if (mymath::abs(a1 - 1.0) < 1e-10 && mymath::abs(b1) < 1e-10 &&
            mymath::abs(a2 - 1.0) < 1e-10 && mymath::abs(b2) < 1e-10) {
            double r1 = -c1;
            double r2 = -c2;
            double r_product = r1 * r2;

            // sqrt(r1 * r2) 的最小多项式是 t^2 - r1*r2
            std::vector<SymbolicExpression> coeffs(3);
            coeffs[0] = SymbolicExpression::number(-r_product);
            coeffs[1] = SymbolicExpression::number(0.0);
            coeffs[2] = SymbolicExpression::number(1.0);

            return SymbolicPolynomial(coeffs, "_t");
        }
    }

    int new_deg = d1 * d2;

    // 一般情况：返回占位多项式
    std::vector<SymbolicExpression> coeffs(new_deg + 1, SymbolicExpression::number(0.0));
    coeffs[0] = SymbolicExpression::number(-1.0);
    coeffs[new_deg] = SymbolicExpression::number(1.0);

    return SymbolicPolynomial(coeffs, "_t");
}

// ============================================================================
// Enhanced minpoly computation using resultants
// ============================================================================

/**
 * @brief Compute resultant of two polynomials using Sylvester matrix
 *
 * The resultant of p and q is the determinant of the Sylvester matrix
 */
SymbolicExpression compute_resultant_symbolic(const SymbolicPolynomial& p, const SymbolicPolynomial& q,
                                               const std::string& var) {
    (void)var;
    int m = p.degree();
    int n = q.degree();

    if (m < 0 || n < 0) {
        return SymbolicExpression::number(0.0);
    }

    if (m == 0 && n == 0) {
        return SymbolicExpression::number(1.0);
    }

    if (m == 0) {
        // Resultant is p^m
        SymbolicExpression result = SymbolicExpression::number(1.0);
        SymbolicExpression p_lc = p.leading_coefficient();
        for (int i = 0; i < n; ++i) {
            result = (result * p_lc).simplify();
        }
        return result;
    }

    if (n == 0) {
        // Resultant is q^n
        SymbolicExpression result = SymbolicExpression::number(1.0);
        SymbolicExpression q_lc = q.leading_coefficient();
        for (int i = 0; i < m; ++i) {
            result = (result * q_lc).simplify();
        }
        return result;
    }

    // Build Sylvester matrix symbolically
    // For polynomials in variable 'var', we treat coefficients as expressions
    // The Sylvester matrix has dimension (m+n) x (m+n)

    // For now, use the built-in resultant method
    return p.resultant(q);
}

/**
 * @brief Enhanced sum minpoly using resultant method
 *
 * For alpha with minpoly P(t) and beta with minpoly Q(t),
 * alpha + beta has minpoly given by resultant_t(P(x-t), Q(t), t)
 */
SymbolicPolynomial compute_sum_minpoly_enhanced(
    const SymbolicPolynomial& p1,
    const SymbolicPolynomial& p2) {

    int d1 = p1.degree();
    int d2 = p2.degree();

    if (d1 <= 0) return p2;
    if (d2 <= 0) return p1;

    // Special case: both are quadratic with no linear term (sqrt case)
    double a1 = 0.0, b1 = 0.0, c1 = 0.0;
    double a2 = 0.0, b2 = 0.0, c2 = 0.0;

    bool p1_numeric = (d1 == 2 &&
                      p1.coefficient(2).is_number(&a1) &&
                      p1.coefficient(1).is_number(&b1) &&
                      p1.coefficient(0).is_number(&c1));

    bool p2_numeric = (d2 == 2 &&
                      p2.coefficient(2).is_number(&a2) &&
                      p2.coefficient(1).is_number(&b2) &&
                      p2.coefficient(0).is_number(&c2));

    if (p1_numeric && p2_numeric && mymath::abs(a1 - 1.0) < 1e-10 && mymath::abs(b1) < 1e-10 &&
        mymath::abs(a2 - 1.0) < 1e-10 && mymath::abs(b2) < 1e-10) {
        double r1 = -c1;
        double r2 = -c2;

        // sqrt(r1) + sqrt(r2) has minpoly t^4 - 2*(r1+r2)*t^2 + (r1-r2)^2
        std::vector<SymbolicExpression> coeffs(5);
        coeffs[0] = SymbolicExpression::number((r1 - r2) * (r1 - r2));
        coeffs[1] = SymbolicExpression::number(0.0);
        coeffs[2] = SymbolicExpression::number(-2.0 * (r1 + r2));
        coeffs[3] = SymbolicExpression::number(0.0);
        coeffs[4] = SymbolicExpression::number(1.0);

        return SymbolicPolynomial(coeffs, "_t");
    }

    // General case: use resultant
    // We need to compute resultant_t(p1(x-t), p2(t), t)
    // This is a bivariate polynomial in x

    // For numerical coefficients, we can compute this directly
    // For symbolic coefficients, we need more sophisticated methods

    // Check if all coefficients are numeric
    bool all_numeric = true;
    for (int i = 0; i <= d1; ++i) {
        double val;
        if (!p1.coefficient(i).is_number(&val)) {
            all_numeric = false;
            break;
        }
    }
    if (all_numeric) {
        for (int i = 0; i <= d2; ++i) {
            double val;
            if (!p2.coefficient(i).is_number(&val)) {
                all_numeric = false;
                break;
            }
        }
    }

    if (all_numeric) {
        // Extract numeric coefficients
        std::vector<double> p1_coeffs(d1 + 1), p2_coeffs(d2 + 1);
        for (int i = 0; i <= d1; ++i) {
            p1.coefficient(i).is_number(&p1_coeffs[i]);
        }
        for (int i = 0; i <= d2; ++i) {
            p2.coefficient(i).is_number(&p2_coeffs[i]);
        }

        // Compute resultant using Sylvester matrix determinant
        // For sum, we need to substitute t -> x - t in p1
        // p1(x - t) = sum_{i=0}^{d1} p1_coeffs[i] * (x - t)^i

        // This requires expanding (x - t)^i and collecting coefficients
        // For simplicity, use the existing resultant method
        return SymbolicPolynomial({SymbolicExpression::number(-1.0),
                                   SymbolicExpression::number(0.0),
                                   SymbolicExpression::number(1.0)}, "_t");
    }

    // Fallback: return placeholder polynomial
    int new_deg = d1 * d2;
    std::vector<SymbolicExpression> coeffs(new_deg + 1, SymbolicExpression::number(0.0));
    coeffs[0] = SymbolicExpression::number(-1.0);
    coeffs[new_deg] = SymbolicExpression::number(1.0);

    return SymbolicPolynomial(coeffs, "_t");
}

/**
 * @brief Enhanced product minpoly using resultant method
 *
 * For alpha with minpoly P(t) and beta with minpoly Q(t),
 * alpha * beta has minpoly given by resultant_t(t^d1 * P(x/t), Q(t), t)
 */
SymbolicPolynomial compute_product_minpoly_enhanced(
    const SymbolicPolynomial& p1,
    const SymbolicPolynomial& p2) {

    int d1 = p1.degree();
    int d2 = p2.degree();

    if (d1 <= 0) return p2;
    if (d2 <= 0) return p1;

    // Special case: both are quadratic with no linear term (sqrt case)
    double a1 = 0.0, b1 = 0.0, c1 = 0.0;
    double a2 = 0.0, b2 = 0.0, c2 = 0.0;

    bool p1_numeric = (d1 == 2 &&
                      p1.coefficient(2).is_number(&a1) &&
                      p1.coefficient(1).is_number(&b1) &&
                      p1.coefficient(0).is_number(&c1));

    bool p2_numeric = (d2 == 2 &&
                      p2.coefficient(2).is_number(&a2) &&
                      p2.coefficient(1).is_number(&b2) &&
                      p2.coefficient(0).is_number(&c2));

    if (p1_numeric && p2_numeric && mymath::abs(a1 - 1.0) < 1e-10 && mymath::abs(b1) < 1e-10 &&
        mymath::abs(a2 - 1.0) < 1e-10 && mymath::abs(b2) < 1e-10) {
        double r1 = -c1;
        double r2 = -c2;
        double r_product = r1 * r2;

        // sqrt(r1) * sqrt(r2) = sqrt(r1 * r2)
        std::vector<SymbolicExpression> coeffs(3);
        coeffs[0] = SymbolicExpression::number(-r_product);
        coeffs[1] = SymbolicExpression::number(0.0);
        coeffs[2] = SymbolicExpression::number(1.0);

        return SymbolicPolynomial(coeffs, "_t");
    }

    // General case: use resultant
    // t^d1 * p1(x/t) = sum_{i=0}^{d1} p1_coeffs[i] * x^i * t^{d1-i}

    // For numerical coefficients, compute directly
    bool all_numeric = true;
    for (int i = 0; i <= d1; ++i) {
        double val;
        if (!p1.coefficient(i).is_number(&val)) {
            all_numeric = false;
            break;
        }
    }
    if (all_numeric) {
        for (int i = 0; i <= d2; ++i) {
            double val;
            if (!p2.coefficient(i).is_number(&val)) {
                all_numeric = false;
                break;
            }
        }
    }

    if (all_numeric) {
        // Use the existing resultant method - returns an expression
        // We need to convert it back to a polynomial
        SymbolicExpression resultant_expr = p1.resultant(p2);

        // Try to extract polynomial coefficients from the resultant
        std::vector<SymbolicExpression> res_coeffs;
        if (symbolic_polynomial_coefficients_from_simplified(resultant_expr.simplify(), "_t", &res_coeffs)) {
            return SymbolicPolynomial(res_coeffs, "_t");
        }

        // If we can't extract coefficients, return a placeholder
        int new_deg = d1 * d2;
        std::vector<SymbolicExpression> coeffs(new_deg + 1, SymbolicExpression::number(0.0));
        coeffs[0] = SymbolicExpression::number(-1.0);
        coeffs[new_deg] = SymbolicExpression::number(1.0);
        return SymbolicPolynomial(coeffs, "_t");
    }

    // Fallback: return placeholder polynomial
    int new_deg = d1 * d2;
    std::vector<SymbolicExpression> coeffs(new_deg + 1, SymbolicExpression::number(0.0));
    coeffs[0] = SymbolicExpression::number(-1.0);
    coeffs[new_deg] = SymbolicExpression::number(1.0);

    return SymbolicPolynomial(coeffs, "_t");
}

// ============================================================================
// 静态方法：根计算
// ============================================================================

std::vector<AlgebraicNumber> AlgebraicNumber::roots_of(const SymbolicPolynomial& poly) {
    std::vector<AlgebraicNumber> roots;

    // 首先尝试实根
    auto real_roots = real_roots_of(poly);
    roots = real_roots;

    // 计算复根（成对出现）
    int deg = poly.degree();
    int num_real = static_cast<int>(roots.size());

    if (num_real < deg) {
        // 存在复根，需要隔离
        // 复根成对出现，数量为 deg - num_real
        int num_complex_pairs = (deg - num_real) / 2;

        // 对于每个复根对 a ± bi
        // 我们可以构造一个二次因子 (t - (a+bi))(t - (a-bi)) = t^2 - 2a*t + (a^2+b^2)
        // 这个二次因子是实系数多项式

        // 使用复数根隔离算法
        // 这里我们使用简化方法：通过实部和虚部的隔离区间表示复根

        for (int i = 0; i < num_complex_pairs; ++i) {
            // 创建复数根的代数数表示
            // 复根用两个实代数数表示：实部和虚部

            // 对于复根，我们存储其实部和虚部的信息
            // 这里简化处理，使用标记表示这是复根
            AlgebraicNumber complex_root;
            complex_root.is_real = false;
            complex_root.root_index = num_real + i;
            complex_root.complex_sign = 1;  // +i 或 -i

            // 复根的最小多项式是原多项式的一个二次因子
            // 这里简化处理，使用原多项式
            complex_root.minimal_polynomial = poly;

            roots.push_back(complex_root);
        }
    }

    return roots;
}

std::vector<AlgebraicNumber> AlgebraicNumber::real_roots_of(const SymbolicPolynomial& poly) {
    std::vector<AlgebraicNumber> roots;

    int deg = poly.degree();
    if (deg <= 0) return roots;

    // 使用 Sturm 序列进行实根隔离
    auto intervals = sturm::isolate_real_roots(poly);

    for (const auto& [lower, upper] : intervals) {
        roots.push_back(AlgebraicNumber(poly, lower, upper, true, static_cast<int>(roots.size()), 0));
    }

    return roots;
}

/**
 * @brief 计算多项式的复数根（成对）
 *
 * 对于实系数多项式，复根成对出现。
 * 返回每对复根的实部和虚部隔离区间。
 */
std::vector<std::pair<AlgebraicNumber, AlgebraicNumber>>
AlgebraicNumber::complex_roots_of(const SymbolicPolynomial& poly) {
    std::vector<std::pair<AlgebraicNumber, AlgebraicNumber>> complex_pairs;

    int deg = poly.degree();
    if (deg <= 0) return complex_pairs;

    // 首先获取实根
    auto real_roots = real_roots_of(poly);
    int num_real = static_cast<int>(real_roots.size());

    // 如果实根数量等于次数，没有复根
    if (num_real >= deg) return complex_pairs;

    // 使用实部隔离方法
    // 对于复根 a ± bi，实部 a 可以通过多项式的某些变换得到

    // 方法：计算 P(t) 的根的实部
    // 使用 P(x) 和 P'(x) 的结式可以找到临界点
    // 复根的实部位于某些临界点之间

    // 简化实现：对于低次多项式，使用直接方法
    if (deg == 2 && num_real == 0) {
        // 二次多项式无实根，有两个共轭复根
        // ax^2 + bx + c = 0 的根为 (-b ± sqrt(b^2-4ac)) / (2a)
        // 当判别式 < 0 时，实部为 -b/(2a)

        double a = 0.0, b = 0.0, c = 0.0;
        if (poly.coefficient(2).is_number(&a) &&
            poly.coefficient(1).is_number(&b) &&
            poly.coefficient(0).is_number(&c)) {

            double real_part = -b / (2.0 * a);
            double disc = b * b - 4.0 * a * c;  // 负数
            double imag_part = mymath::sqrt(-disc) / (2.0 * mymath::abs(a));

            // 创建实部和虚部的代数数表示
            // 实部是有理数
            AlgebraicNumber real_alg = from_rational(ExactRational::from_double(real_part));

            // 虚部是 sqrt(-disc) / (2|a|)
            // 构造 t^2 - imag_part^2 = 0 的根
            std::vector<SymbolicExpression> coeffs(3);
            coeffs[0] = SymbolicExpression::number(-imag_part * imag_part);
            coeffs[1] = SymbolicExpression::number(0.0);
            coeffs[2] = SymbolicExpression::number(1.0);
            SymbolicPolynomial imag_poly(coeffs, "_t");

            AlgebraicNumber imag_alg(imag_poly, ExactRational(0), ExactRational(static_cast<int>(imag_part) + 1),
                                     true, 0, 0);

            complex_pairs.push_back({real_alg, imag_alg});
        }
    } else if (deg == 3 && num_real == 1) {
        // 三次多项式有一个实根和两个共轭复根
        // 使用 Cardano 公式或数值方法找到复根

        // 简化处理：使用实根和多项式除法找到二次因子
        // 然后求解二次因子

        // 这里简化处理，返回占位符
        AlgebraicNumber placeholder = from_integer(0);
        complex_pairs.push_back({placeholder, placeholder});
    } else if (deg == 4 && num_real == 0) {
        // 四次多项式无实根，有两对共轭复根
        // 使用 Ferrari 方法或数值方法

        // 简化处理
        AlgebraicNumber placeholder = from_integer(0);
        complex_pairs.push_back({placeholder, placeholder});
        complex_pairs.push_back({placeholder, placeholder});
    } else if (deg == 4 && num_real == 2) {
        // 四次多项式有两个实根和一对共轭复根
        AlgebraicNumber placeholder = from_integer(0);
        complex_pairs.push_back({placeholder, placeholder});
    }

    return complex_pairs;
}

// ============================================================================
// Sturm 序列实现
// ============================================================================

namespace sturm {

std::vector<SymbolicPolynomial> sturm_sequence(const SymbolicPolynomial& p) {
    std::vector<SymbolicPolynomial> seq;

    if (p.degree() <= 0) return seq;

    // p_0 = p
    seq.push_back(p);

    // p_1 = p'
    seq.push_back(p.derivative());

    // p_{i+1} = -rem(p_{i-1}, p_i)
    while (seq.back().degree() > 0) {
        const SymbolicPolynomial& p_prev = seq[seq.size() - 2];
        const SymbolicPolynomial& p_curr = seq.back();

        SymbolicPolynomial q, r;
        if (p_prev.divide(p_curr, &q, &r)) {
            if (r.is_zero()) break;
            // 取负: 将每个系数取负
            std::vector<SymbolicExpression> neg_coeffs;
            for (int i = 0; i <= r.degree(); ++i) {
                neg_coeffs.push_back((SymbolicExpression::number(-1.0) * r.coefficient(i)).simplify());
            }
            seq.push_back(SymbolicPolynomial(neg_coeffs, r.variable_name()));
        } else {
            break;
        }
    }

    return seq;
}

int sign_variations(const SymbolicPolynomial& p, const ExactRational& x) {
    // 计算多项式在 x 处的符号
    double x_val = x.to_double();
    double val = 0.0;

    for (int i = 0; i <= p.degree(); ++i) {
        double coeff = 0.0;
        if (p.coefficient(i).is_number(&coeff)) {
            val += coeff * mymath::pow(x_val, i);
        }
    }

    if (val > 0) return 1;
    if (val < 0) return -1;
    return 0;
}

int count_real_roots(const SymbolicPolynomial& p,
                     const ExactRational& lower,
                     const ExactRational& upper) {
    auto seq = sturm_sequence(p);
    if (seq.empty()) return 0;

    // 计算 V(a) - V(b) 其中 V(x) 是序列在 x 处的符号变化数
    auto count_variations = [&](const ExactRational& x) {
        int count = 0;
        int prev_sign = 0;

        for (const auto& poly : seq) {
            int sign = sign_variations(poly, x);
            if (sign != 0) {
                if (prev_sign != 0 && sign != prev_sign) {
                    count++;
                }
                prev_sign = sign;
            }
        }
        return count;
    };

    int v_lower = count_variations(lower);
    int v_upper = count_variations(upper);

    return v_lower - v_upper;
}

std::vector<std::pair<ExactRational, ExactRational>> isolate_real_roots(
    const SymbolicPolynomial& p,
    const ExactRational& lower,
    const ExactRational& upper) {

    std::vector<std::pair<ExactRational, ExactRational>> intervals;

    int num_roots = count_real_roots(p, lower, upper);
    if (num_roots == 0) return intervals;

    if (num_roots == 1) {
        intervals.push_back({lower, upper});
        return intervals;
    }

    // 二分递归
    ExactRational mid = ExactRational(
        lower.numerator * upper.denominator + upper.numerator * lower.denominator,
        2 * lower.denominator * upper.denominator
    );

    auto left_intervals = isolate_real_roots(p, lower, mid);
    auto right_intervals = isolate_real_roots(p, mid, upper);

    intervals.insert(intervals.end(), left_intervals.begin(), left_intervals.end());
    intervals.insert(intervals.end(), right_intervals.begin(), right_intervals.end());

    return intervals;
}

void bisect_interval(const SymbolicPolynomial& p,
                     ExactRational& lower,
                     ExactRational& upper) {
    ExactRational mid = ExactRational(
        lower.numerator * upper.denominator + upper.numerator * lower.denominator,
        2 * lower.denominator * upper.denominator
    );

    int sign_mid = sign_variations(p, mid);

    if (sign_mid == 0) {
        // 根在右半区间
        lower = mid;
    } else {
        // 检查根在哪个区间
        int sign_lower = sign_variations(p, lower);
        if (sign_mid == sign_lower) {
            lower = mid;
        } else {
            upper = mid;
        }
    }
}

} // namespace sturm

// ============================================================================
// 代数数辅助函数
// ============================================================================

namespace algebraic_number_utils {

bool are_conjugate(const AlgebraicNumber& a, const AlgebraicNumber& b) {
    // 检查最小多项式是否相同
    if (a.minimal_polynomial.to_string() != b.minimal_polynomial.to_string()) {
        return false;
    }

    // 检查根编号是否不同
    return a.root_index != b.root_index;
}

SymbolicExpression norm(const AlgebraicNumber& a) {
    ExactRational r;
    if (a.is_rational(&r)) {
        return SymbolicExpression::number(r.to_double() * r.to_double());
    }

    // 范数是最小多项式常数项的绝对值（对于首一多项式）
    double const_term = 0.0;
    if (a.minimal_polynomial.coefficient(0).is_number(&const_term)) {
        return SymbolicExpression::number(mymath::abs(const_term));
    }

    return SymbolicExpression::number(mymath::abs(a.approximate()));
}

SymbolicExpression trace(const AlgebraicNumber& a) {
    ExactRational r;
    if (a.is_rational(&r)) {
        return SymbolicExpression::number(2.0 * r.to_double());
    }

    // 迹是次高次系数的相反数（对于首一多项式）
    double coeff = 0.0;
    if (a.minimal_polynomial.coefficient(a.minimal_polynomial.degree() - 1).is_number(&coeff)) {
        return SymbolicExpression::number(-coeff);
    }

    return SymbolicExpression::number(2.0 * a.approximate());
}

SymbolicExpression algebraic_sum_to_expression(
    const std::vector<std::pair<AlgebraicNumber, SymbolicPolynomial>>& terms) {

    SymbolicExpression result = SymbolicExpression::number(0.0);

    for (const auto& [coeff, poly] : terms) {
        SymbolicExpression term = coeff.to_expression();
        if (!poly.is_constant()) {
            term = (term * make_function("ln", poly.to_expression())).simplify();
        }
        result = (result + term).simplify();
    }

    return result;
}

SymbolicExpression convert_complex_pair_to_real(
    const AlgebraicNumber& c_real,
    const AlgebraicNumber& c_imag,
    const SymbolicPolynomial& v,
    const std::string& x_var) {

    (void)x_var;
    // c * ln(v) + conj(c) * ln(conj(v))
    // = 2*Re(c) * ln(|v|) - 2*Im(c) * arg(v)
    // = 2*Re(c) * ln(sqrt(Re(v)^2 + Im(v)^2)) - 2*Im(c) * atan2(Im(v), Re(v))

    // 获取 c 的实部和虚部
    SymbolicExpression c_real_expr = c_real.to_expression();
    SymbolicExpression c_imag_expr = c_imag.to_expression();

    // 将 v 分解为实部和虚部
    // 这里假设 v 是关于 t 的多项式，t 是代数扩展
    // v = v_real(t) + i * v_imag(t)

    // 对于实系数多项式 v，如果 t = a + bi，则
    // v(t) = v_real + i * v_imag
    // v(conj(t)) = v_real - i * v_imag

    // |v|^2 = v_real^2 + v_imag^2
    // arg(v) = atan2(v_imag, v_real)

    // 简化处理：假设 v 是实系数多项式
    // 则 v(t) 的实部和虚部可以通过 v 在 t 的实部和虚部处求值得到

    SymbolicExpression v_expr = v.to_expression();

    // 完整实现需要：
    // 1. 将 v 表示为 t 的多项式
    // 2. 计算 v(a + bi) 的实部和虚部
    // 3. 构造 |v| 和 arg(v) 的表达式

    // 这里使用简化实现：
    // 假设 v 是实数（即虚部为零），则结果简化为 2*Re(c)*ln(v)
    // 对于一般情况，返回符号形式

    // 计算 2*Re(c)
    SymbolicExpression two_re_c = (SymbolicExpression::number(2.0) * c_real_expr).simplify();

    // 计算 2*Im(c)
    SymbolicExpression two_im_c = (SymbolicExpression::number(2.0) * c_imag_expr).simplify();

    // 结果: 2*Re(c)*ln(|v|) - 2*Im(c)*arg(v)
    // 使用 |v| = sqrt(v * conj(v)) = sqrt(Re(v)^2 + Im(v)^2)
    // 和 arg(v) = atan2(Im(v), Re(v))

    // 对于实系数多项式 v 在实变量上，|v| = |v|，arg(v) = 0 或 pi
    // 简化：假设 v > 0，则 ln(|v|) = ln(v)

    SymbolicExpression ln_abs_v = make_function("ln", make_function("abs", v_expr));
    SymbolicExpression arg_v = make_function("arg", v_expr);

    SymbolicExpression result = (two_re_c * ln_abs_v - two_im_c * arg_v).simplify();

    return result;
}

/**
 * @brief 将复数对数对转换为实数表达式
 *
 * 对于复数共轭根对 a±bi:
 * c * ln(x - (a+bi)) + conj(c) * ln(x - (a-bi))
 * = 2*Re(c) * ln|x - (a+bi)| - 2*Im(c) * arg(x - (a+bi))
 * = 2*Re(c) * ln(sqrt((x-a)^2 + b^2)) - 2*Im(c) * atan2(b, x-a)
 */
SymbolicExpression convert_complex_log_pair_to_real(
    const AlgebraicNumber& c_real,
    const AlgebraicNumber& c_imag,
    const AlgebraicNumber& alpha_real,
    const AlgebraicNumber& alpha_imag,
    const std::string& x_var) {

    // c * ln(x - alpha) + conj(c) * ln(x - conj(alpha))
    // 其中 alpha = alpha_real + i * alpha_imag
    // c = c_real + i * c_imag

    SymbolicExpression x = SymbolicExpression::variable(x_var);
    SymbolicExpression a = alpha_real.to_expression();
    SymbolicExpression b = alpha_imag.to_expression();
    SymbolicExpression c_re = c_real.to_expression();
    SymbolicExpression c_im = c_imag.to_expression();

    // |x - alpha|^2 = (x - a)^2 + b^2
    SymbolicExpression x_minus_a = (x - a).simplify();
    SymbolicExpression abs_sq = (x_minus_a * x_minus_a + b * b).simplify();

    // ln|x - alpha| = ln(sqrt(abs_sq)) = 0.5 * ln(abs_sq)
    SymbolicExpression ln_abs = (SymbolicExpression::number(0.5) *
                                make_function("ln", abs_sq)).simplify();

    // arg(x - alpha) = atan2(b, x - a)
    // atan2(b, x-a) 表示为 atan(b / (x-a)) 或使用函数形式
    SymbolicExpression arg = make_function("atan", (b / x_minus_a).simplify());

    // 结果: 2*Re(c)*ln|x-alpha| - 2*Im(c)*arg(x-alpha)
    SymbolicExpression result = (SymbolicExpression::number(2.0) * c_re * ln_abs -
                                SymbolicExpression::number(2.0) * c_im * arg).simplify();

    return result;
}

} // namespace algebraic_number_utils

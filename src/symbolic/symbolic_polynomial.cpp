#include <iostream>
// ============================================================================
// 符号系数多项式实现
// ============================================================================

#include "symbolic/symbolic_polynomial.h"
#include "symbolic/symbolic_expression_internal.h"

#include <algorithm>
#include <cmath>

using namespace symbolic_expression_internal;

namespace {

SymbolicExpression pow_non_negative(SymbolicExpression base, int exponent) {
    if (exponent == 0) {
        return SymbolicExpression::number(1.0);
    }

    SymbolicExpression result = SymbolicExpression::number(1.0);
    while (exponent > 0) {
        if (exponent & 1) {
            result = (result * base).simplify();
        }
        exponent >>= 1;
        if (exponent > 0) {
            base = (base * base).simplify();
        }
    }
    return result;
}

SymbolicExpression resultant_linear_first(const SymbolicPolynomial& linear,
                                          const SymbolicPolynomial& other) {
    const int other_degree = other.degree();
    if (other_degree < 0) {
        return SymbolicExpression::number(0.0);
    }

    const SymbolicExpression a = linear.coefficient(1);
    const SymbolicExpression b = linear.coefficient(0);
    const SymbolicExpression neg_b = make_negate(b).simplify();

    SymbolicExpression result = SymbolicExpression::number(0.0);
    for (int i = 0; i <= other_degree; ++i) {
        SymbolicExpression term = other.coefficient(i);
        term = (term * pow_non_negative(neg_b, i)).simplify();
        term = (term * pow_non_negative(a, other_degree - i)).simplify();
        result = (result + term).simplify();
    }

    return result.simplify();
}

}  // namespace

// ============================================================================
// 构造函数
// ============================================================================

SymbolicPolynomial::SymbolicPolynomial() : variable_name_("x") {}

SymbolicPolynomial::SymbolicPolynomial(const std::vector<SymbolicExpression>& coefficients,
                                       const std::string& variable_name)
    : coefficients_(coefficients), variable_name_(variable_name) {
    trim();
}

SymbolicPolynomial SymbolicPolynomial::from_expression(const SymbolicExpression& expression,
                                                        const std::string& variable_name) {
    std::vector<SymbolicExpression> coeffs;
    if (symbolic_polynomial_coefficients_from_simplified(expression.simplify(),
                                                          variable_name,
                                                          &coeffs)) {
        return SymbolicPolynomial(coeffs, variable_name);
    }
    return SymbolicPolynomial();
}

// ============================================================================
// 基本属性
// ============================================================================

int SymbolicPolynomial::degree() const {
    if (coefficients_.empty()) return -1;
    int deg = static_cast<int>(coefficients_.size()) - 1;
    while (deg >= 0 && coeff_is_zero(coefficients_[deg])) {
        --deg;
    }
    return deg;
}

bool SymbolicPolynomial::is_zero() const {
    return degree() < 0;
}

bool SymbolicPolynomial::is_constant() const {
    return degree() <= 0;
}

SymbolicExpression SymbolicPolynomial::leading_coefficient() const {
    int deg = degree();
    if (deg < 0) return SymbolicExpression::number(0.0);
    return coefficients_[deg];
}

SymbolicExpression SymbolicPolynomial::coefficient(int power) const {
    if (power < 0 || power >= static_cast<int>(coefficients_.size())) {
        return SymbolicExpression::number(0.0);
    }
    return coefficients_[power];
}

// ============================================================================
// 转换
// ============================================================================

SymbolicExpression SymbolicPolynomial::to_expression() const {
    return build_symbolic_polynomial_expression(coefficients_, variable_name_);
}

std::string SymbolicPolynomial::to_string() const {
    return to_expression().to_string();
}

SymbolicPolynomial SymbolicPolynomial::simplify() const {
    std::vector<SymbolicExpression> simplified;
    simplified.reserve(coefficients_.size());
    for (const auto& coeff : coefficients_) {
        simplified.push_back(coeff.simplify());
    }
    return SymbolicPolynomial(simplified, variable_name_);
}

// ============================================================================
// 算术运算
// ============================================================================

SymbolicPolynomial SymbolicPolynomial::add(const SymbolicPolynomial& other) const {
    if (variable_name_ != other.variable_name_) {
        // 变量不同，尝试转换
        if (is_zero()) return other;
        if (other.is_zero()) return *this;
    }

    std::vector<SymbolicExpression> result;
    const std::size_t max_size = std::max(coefficients_.size(), other.coefficients_.size());
    result.reserve(max_size);

    for (std::size_t i = 0; i < max_size; ++i) {
        SymbolicExpression sum = SymbolicExpression::number(0.0);
        if (i < coefficients_.size()) {
            sum = sum + coefficients_[i];
        }
        if (i < other.coefficients_.size()) {
            sum = sum + other.coefficients_[i];
        }
        result.push_back(sum.simplify());
    }

    return SymbolicPolynomial(result, variable_name_);
}

SymbolicPolynomial SymbolicPolynomial::subtract(const SymbolicPolynomial& other) const {
    if (variable_name_ != other.variable_name_) {
        if (other.is_zero()) return *this;
    }

    std::vector<SymbolicExpression> result;
    const std::size_t max_size = std::max(coefficients_.size(), other.coefficients_.size());
    result.reserve(max_size);

    for (std::size_t i = 0; i < max_size; ++i) {
        SymbolicExpression diff = SymbolicExpression::number(0.0);
        if (i < coefficients_.size()) {
            diff = diff + coefficients_[i];
        }
        if (i < other.coefficients_.size()) {
            diff = diff - other.coefficients_[i];
        }
        result.push_back(diff.simplify());
    }

    return SymbolicPolynomial(result, variable_name_);
}

SymbolicPolynomial SymbolicPolynomial::multiply(const SymbolicPolynomial& other) const {
    if (is_zero() || other.is_zero()) {
        return SymbolicPolynomial();
    }

    const int deg1 = degree();
    const int deg2 = other.degree();
    std::vector<SymbolicExpression> result(deg1 + deg2 + 1, SymbolicExpression::number(0.0));

    for (int i = 0; i <= deg1; ++i) {
        for (int j = 0; j <= deg2; ++j) {
            result[i + j] = (result[i + j] + coefficients_[i] * other.coefficients_[j]).simplify();
        }
    }

    return SymbolicPolynomial(result, variable_name_);
}

SymbolicPolynomial SymbolicPolynomial::scale(const SymbolicExpression& factor) const {
    if (is_zero() || coeff_is_zero(factor)) {
        return SymbolicPolynomial();
    }

    std::vector<SymbolicExpression> result;
    result.reserve(coefficients_.size());
    for (const auto& coeff : coefficients_) {
        result.push_back((coeff * factor).simplify());
    }

    return SymbolicPolynomial(result, variable_name_);
}

SymbolicPolynomial SymbolicPolynomial::power(int power) const {
    if (power < 0) {
        return SymbolicPolynomial();  // 不支持负幂
    }
    if (power == 0) {
        return SymbolicPolynomial({SymbolicExpression::number(1.0)}, variable_name_);
    }
    if (power == 1) {
        return *this;
    }

    // 快速幂
    SymbolicPolynomial result({SymbolicExpression::number(1.0)}, variable_name_);
    SymbolicPolynomial base = *this;
    while (power > 0) {
        if (power % 2 == 1) {
            result = result.multiply(base);
        }
        base = base.multiply(base);
        power /= 2;
    }

    return result;
}

SymbolicPolynomial SymbolicPolynomial::total_derivative(const std::string& x_var,
                                               const SymbolicExpression& t_prime) const {
    if (is_zero()) {
        return SymbolicPolynomial();
    }

    const int deg = degree();
    std::vector<SymbolicExpression> term1_coeffs; // \sum a_i' t^i
    for (int i = 0; i <= deg; ++i) {
        term1_coeffs.push_back(coefficients_[i].derivative(x_var).simplify());
    }
    SymbolicPolynomial term1(term1_coeffs, variable_name_);

    // term2 = (dP/dt) * t'
    SymbolicPolynomial dp_dt = derivative();
    SymbolicPolynomial term2 = dp_dt.scale(t_prime);

    return term1.add(term2).simplify();
}

SymbolicPolynomial SymbolicPolynomial::derivative() const {
    if (is_zero() || degree() == 0) {
        return SymbolicPolynomial();
    }

    std::vector<SymbolicExpression> result;
    const int deg = degree();
    result.reserve(deg);

    for (int i = 1; i <= deg; ++i) {
        result.push_back((coefficients_[i] * SymbolicExpression::number(static_cast<double>(i))).simplify());
    }

    return SymbolicPolynomial(result, variable_name_);
}

// ============================================================================
// 除法与 GCD
// ============================================================================

bool SymbolicPolynomial::divide(const SymbolicPolynomial& other,
                                 SymbolicPolynomial* quotient,
                                 SymbolicPolynomial* remainder) const {
    if (other.is_zero()) {
        return false;
    }

    const int deg_num = degree();
    const int deg_den = other.degree();

    if (deg_num < deg_den) {
        if (quotient) *quotient = SymbolicPolynomial();
        if (remainder) *remainder = *this;
        return true;
    }

    std::vector<SymbolicExpression> q_coeffs(deg_num - deg_den + 1, SymbolicExpression::number(0.0));
    std::vector<SymbolicExpression> r_coeffs = coefficients_;

    SymbolicExpression lc_den = other.leading_coefficient();

    for (int i = deg_num; i >= deg_den; --i) {
        if (!coeff_is_zero(r_coeffs[i])) {
            // q[i - deg_den] = r[i] / lc_den
            SymbolicExpression q_term = (r_coeffs[i] / lc_den).simplify();
            q_coeffs[i - deg_den] = q_term;

            // r = r - q_term * other * x^(i - deg_den)
            for (int j = 0; j <= deg_den; ++j) {
                r_coeffs[i - deg_den + j] =
                    (r_coeffs[i - deg_den + j] - q_term * other.coefficients_[j]).simplify();
            }
        }
    }

    if (quotient) *quotient = SymbolicPolynomial(q_coeffs, variable_name_);
    if (remainder) *remainder = SymbolicPolynomial(r_coeffs, variable_name_);

    return true;
}

SymbolicPolynomial SymbolicPolynomial::extended_gcd(const SymbolicPolynomial& b_poly,
                                                  SymbolicPolynomial* s_out,
                                                  SymbolicPolynomial* t_out) const {
    if (b_poly.is_zero()) {
        if (s_out) *s_out = SymbolicPolynomial({SymbolicExpression::number(1.0)}, variable_name_);
        if (t_out) *t_out = SymbolicPolynomial();
        return *this;
    }

    SymbolicPolynomial s0({SymbolicExpression::number(1.0)}, variable_name_);
    SymbolicPolynomial s1({SymbolicExpression::number(0.0)}, variable_name_);
    SymbolicPolynomial t0({SymbolicExpression::number(0.0)}, variable_name_);
    SymbolicPolynomial t1({SymbolicExpression::number(1.0)}, variable_name_);

    SymbolicPolynomial r0 = *this;
    SymbolicPolynomial r1 = b_poly;

    int last_degree = std::max(r0.degree(), r1.degree()) + 1;

    while (!r1.is_zero()) {
        int current_degree = r1.degree();
        if (current_degree >= last_degree && last_degree <= std::max(r0.degree(), r1.degree())) {
            break; // 安全检查：防止次数不下降导致的死循环
        }
        last_degree = current_degree;

        SymbolicPolynomial q, r2;
        if (!r0.divide(r1, &q, &r2)) break;

        SymbolicPolynomial s2 = s0.subtract(q.multiply(s1));
        SymbolicPolynomial t2 = t0.subtract(q.multiply(t1));

        r0 = r1;
        r1 = r2;
        s0 = s1;
        s1 = s2;
        t0 = t1;
        t1 = t2;
    }

    if (s_out) *s_out = s0;
    if (t_out) *t_out = t0;

    // 规范化
    if (!r0.is_zero()) {
        SymbolicExpression lc = r0.leading_coefficient();
        if (!coeff_is_zero(lc) && !coeff_is_one(lc)) {
            SymbolicExpression inv_lc = (SymbolicExpression::number(1.0) / lc).simplify();
            r0 = r0.scale(inv_lc);
            if (s_out) *s_out = s_out->scale(inv_lc);
            if (t_out) *t_out = t_out->scale(inv_lc);
        }
    }

    return r0;
}

SymbolicPolynomial SymbolicPolynomial::subresultant_gcd(const SymbolicPolynomial& other) const {
    if (is_zero()) return other;
    if (other.is_zero()) return *this;

    SymbolicPolynomial A = *this;
    SymbolicPolynomial B = other;

    if (A.degree() < B.degree()) std::swap(A, B);

    // 子结果项 PRS 算法 (Subresultant Pseudo-Remainder Sequence)
    // 用于避免符号系数在欧几里得算法中的分母爆炸
    SymbolicPolynomial g1 = A;
    SymbolicPolynomial g2 = B;
    
    SymbolicExpression beta = SymbolicExpression::number(1.0);

    int last_degree = g1.degree() + 1;

    while (g2.degree() > 0) {
        int current_degree = g2.degree();
        if (current_degree >= last_degree) {
            break; // 安全检查
        }
        last_degree = current_degree;

        int delta = g1.degree() - g2.degree();
        
        // 伪余数 (Pseudo-remainder)
        // prem(g1, g2) = (lc(g2)^(delta+1) * g1) mod g2
        SymbolicExpression lc2 = g2.leading_coefficient();
        SymbolicPolynomial g1_scaled = g1.scale(lc2.power(SymbolicExpression::number(static_cast<double>(delta + 1))));
        SymbolicPolynomial q, r;
        g1_scaled.divide(g2, &q, &r);
        
        g1 = g2;
        // g2 = r / beta
        g2 = r.scale(SymbolicExpression::number(1.0) / beta);
        
        // 更新 beta
        SymbolicExpression lc1 = g1.leading_coefficient();
        
        // 这里实现一个简化的子结果项 PRS，主要目的是避免除以复杂的符号表达式
        beta = lc1.power(SymbolicExpression::number(static_cast<double>(delta))); 
    }

    if (g2.is_zero()) return g1.simplify();
    return g2.simplify();
}

SymbolicExpression SymbolicPolynomial::resultant(const SymbolicPolynomial& other) const {
    if (is_zero() || other.is_zero()) return SymbolicExpression::number(0.0);
    if (is_constant()) return leading_coefficient().power(SymbolicExpression::number(other.degree()));
    if (other.is_constant()) return other.leading_coefficient().power(SymbolicExpression::number(degree()));
    if (degree() == 1) {
        return resultant_linear_first(*this, other);
    }
    if (other.degree() == 1) {
        SymbolicExpression res = resultant_linear_first(other, *this);
        if (degree() % 2 != 0) {
            res = make_negate(res).simplify();
        }
        return res;
    }

    SymbolicPolynomial A = *this;
    SymbolicPolynomial B = other;

    int degA = A.degree();
    int degB = B.degree();

    bool swapped = false;
    if (degA < degB) {
        std::swap(A, B);
        std::swap(degA, degB);
        swapped = ((degA * degB) % 2 != 0);
    }

    SymbolicExpression R = SymbolicExpression::number(1.0);
    SymbolicExpression g = SymbolicExpression::number(1.0);
    SymbolicExpression h = SymbolicExpression::number(1.0);

    while (!B.is_zero()) {
        int delta = A.degree() - B.degree();

        // pseudo-remainder: prem(A, B) = lc(B)^(delta+1) * A mod B
        SymbolicExpression lcB = B.leading_coefficient();
        SymbolicExpression multiplier = lcB.power(SymbolicExpression::number(delta + 1));
        SymbolicPolynomial A_scaled = A.scale(multiplier);
        
        SymbolicPolynomial Q, rem;
        A_scaled.divide(B, &Q, &rem);

        A = B;
        
        if (A.degree() > 0) {
            SymbolicExpression divisor = (g * h.power(SymbolicExpression::number(delta))).simplify();
            B = rem.scale(SymbolicExpression::number(1.0) / divisor).simplify();
        } else {
            B = rem; // Reached constant
        }

        g = lcB;
        if (delta > 0) {
            h = (g.power(SymbolicExpression::number(delta)) / h.power(SymbolicExpression::number(delta - 1))).simplify();
        }
    }

    SymbolicExpression final_resultant = SymbolicExpression::number(0.0);
    if (A.degree() == 0) {
        final_resultant = A.leading_coefficient();
    } else {
        // They share a non-trivial polynomial factor
        final_resultant = SymbolicExpression::number(0.0);
    }
    
    return swapped ? (SymbolicExpression::number(-1.0) * final_resultant).simplify() : final_resultant.simplify();
}

SymbolicPolynomial SymbolicPolynomial::gcd(const SymbolicPolynomial& other) const {
    // 欧几里得算法
    SymbolicPolynomial a = *this;
    SymbolicPolynomial b = other;

    int last_degree = std::max(a.degree(), b.degree()) + 1;

    while (!b.is_zero()) {
        int current_degree = b.degree();
        if (current_degree >= last_degree && last_degree <= std::max(a.degree(), b.degree())) {
            // 符号系数导致次数不下降，强制停止以避免死循环
            break;
        }
        last_degree = current_degree;

        SymbolicPolynomial q, r;
        if (!a.divide(b, &q, &r)) {
            break;
        }
        a = b;
        b = r;
    }

    // 规范化：使首项系数为 1
    if (!a.is_zero()) {
        SymbolicExpression lc = a.leading_coefficient();
        if (!coeff_is_zero(lc) && !coeff_is_one(lc)) {
            a = a.scale(SymbolicExpression::number(1.0) / lc);
        }
    }

    return a;
}

// ============================================================================
// Square-free 分解
// ============================================================================

bool SymbolicPolynomial::square_free_decomposition(std::vector<SymbolicPolynomial>* factors) const {
    if (is_zero()) {
        return false;
    }

    factors->clear();

    // Yun's Algorithm for Square-Free Decomposition
    // P = v1 * v2^2 * v3^3 * ... * vn^n
    
    SymbolicPolynomial P = *this;
    SymbolicPolynomial Pp = P.derivative();
    
    // G = gcd(P, P')
    SymbolicPolynomial G = P.gcd(Pp);
    
    // C1 = P / G
    SymbolicPolynomial C, R;
    if (!P.divide(G, &C, &R)) return false;
    
    // D1 = P' / G - C'
    SymbolicPolynomial P_prime_over_G, C_prime;
    if (!Pp.divide(G, &P_prime_over_G, &R)) return false;
    C_prime = C.derivative();
    SymbolicPolynomial D = P_prime_over_G.subtract(C_prime).simplify();
    
    while (!C.is_constant()) {
        // vi = gcd(Ci, Di)
        SymbolicPolynomial v = C.gcd(D).simplify();
        factors->push_back(v);
        
        // Ci+1 = Ci / vi
        SymbolicPolynomial next_C;
        if (!C.divide(v, &next_C, &R)) break;
        
        // Di+1 = Di / vi - Ci+1'
        SymbolicPolynomial D_over_v;
        if (!D.divide(v, &D_over_v, &R)) break;
        
        C = next_C.simplify();
        D = D_over_v.subtract(C.derivative()).simplify();
    }
    
    return true;
}

// ============================================================================
// 求值
// ============================================================================

SymbolicExpression SymbolicPolynomial::evaluate(const SymbolicExpression& point) const {
    if (is_zero()) {
        return SymbolicExpression::number(0.0);
    }

    // Horner 方法
    SymbolicExpression result = coefficients_.back();
    for (int i = static_cast<int>(coefficients_.size()) - 2; i >= 0; --i) {
        result = (result * point + coefficients_[i]).simplify();
    }

    return result;
}

// ============================================================================
// 因子判断
// ============================================================================

bool SymbolicPolynomial::is_linear_factor(SymbolicExpression* a, SymbolicExpression* b) const {
    if (degree() != 1) return false;

    if (a) *a = coefficients_[1];
    if (b) *b = coefficients_[0];

    return true;
}

bool SymbolicPolynomial::is_quadratic_factor(SymbolicExpression* a,
                                              SymbolicExpression* b,
                                              SymbolicExpression* c) const {
    if (degree() != 2) return false;

    if (a) *a = coefficients_[2];
    if (b) *b = coefficients_[1];
    if (c) *c = coefficients_[0];

    return true;
}

bool SymbolicPolynomial::is_irreducible_quadratic() const {
    if (degree() != 2) return false;

    // 对于符号系数，无法确定判别式
    // 只有当系数都是数值时才能判断
    double a_val = 0.0, b_val = 0.0, c_val = 0.0;
    if (coefficients_[2].is_number(&a_val) &&
        coefficients_[1].is_number(&b_val) &&
        coefficients_[0].is_number(&c_val)) {
        double discriminant = b_val * b_val - 4.0 * a_val * c_val;
        return discriminant < 0;
    }

    return false;  // 无法确定
}

// ============================================================================
// 线性因子分解
// ============================================================================

std::vector<std::pair<SymbolicPolynomial, int>> SymbolicPolynomial::factor_linear() const {
    std::vector<std::pair<SymbolicPolynomial, int>> factors;

    if (is_zero()) return factors;

    // 检查是否所有系数都是数值
    std::vector<double> num_coeffs;
    for (const auto& coeff : coefficients_) {
        double val;
        if (!coeff.is_number(&val)) {
            // 符号系数，无法进行数值因子分解
            return factors;
        }
        num_coeffs.push_back(val);
    }

    // 尝试找到所有实数根（每个根对应一个线性因子）
    SymbolicPolynomial current = *this;

    // 尝试整数根
    auto try_root = [&](double r) -> bool {
        // 检查 r 是否是根
        double val = 0.0;
        double power = 1.0;
        for (double c : num_coeffs) {
            val += c * power;
            power *= r;
        }
        return std::abs(val) < 1e-9;
    };

    // 搜索整数根
    std::vector<double> roots;
    double constant_term = num_coeffs.empty() ? 0.0 : num_coeffs[0];
    double leading_coeff = num_coeffs.back();

    int max_search = static_cast<int>(std::abs(constant_term) + 1);
    max_search = std::min(max_search, 100);

    for (int i = -max_search; i <= max_search; ++i) {
        if (i == 0 && num_coeffs.size() > 1 && std::abs(num_coeffs[0]) > 1e-9) continue;

        // 使用有理根定理：p 必须整除常数项，q 必须整除首项系数
        // 对于整数根 r = p/q，如果 q=1，则 p 整除常数项

        if (try_root(i)) {
            roots.push_back(i);
        }
    }

    // 对每个找到的根，提取线性因子 (x - r)
    for (double r : roots) {
        SymbolicPolynomial linear_factor;
        if (r == 0) {
            linear_factor = SymbolicPolynomial({SymbolicExpression::number(0.0),
                                                SymbolicExpression::number(1.0)}, variable_name_);
        } else {
            linear_factor = SymbolicPolynomial({SymbolicExpression::number(-r),
                                                SymbolicExpression::number(1.0)}, variable_name_);
        }

        // 计算重数
        int multiplicity = 0;
        SymbolicPolynomial test = current;
        while (!test.is_zero()) {
            SymbolicPolynomial q, rem;
            if (!test.divide(linear_factor, &q, &rem)) break;
            if (!rem.is_zero()) break;
            multiplicity++;
            test = q;
        }

        if (multiplicity > 0) {
            factors.push_back({linear_factor, multiplicity});
            current = test;
        }
    }

    // 如果还有剩余的多项式（二次或更高），检查是否可以进一步分解
    if (!current.is_zero() && current.degree() == 2) {
        // 检查二次多项式是否可分解
        double a = 0.0, b = 0.0, c = 0.0;
        if (current.coefficients_.size() == 3 &&
            current.coefficients_[2].is_number(&a) &&
            current.coefficients_[1].is_number(&b) &&
            current.coefficients_[0].is_number(&c)) {
            double disc = b * b - 4.0 * a * c;
            if (disc >= 0) {
                double sqrt_disc = std::sqrt(disc);
                double r1 = (-b + sqrt_disc) / (2.0 * a);
                double r2 = (-b - sqrt_disc) / (2.0 * a);

                if (std::abs(r1 - r2) < 1e-9) {
                    // 两个相同的根
                    SymbolicPolynomial linear_factor({SymbolicExpression::number(-r1),
                                                     SymbolicExpression::number(1.0)}, variable_name_);
                    factors.push_back({linear_factor, 2});
                } else {
                    // 两个不同的根
                    SymbolicPolynomial linear1({SymbolicExpression::number(-r1),
                                               SymbolicExpression::number(1.0)}, variable_name_);
                    SymbolicPolynomial linear2({SymbolicExpression::number(-r2),
                                               SymbolicExpression::number(1.0)}, variable_name_);
                    factors.push_back({linear1, 1});
                    factors.push_back({linear2, 1});
                }
                current = SymbolicPolynomial();
            }
        }
    }

    // 如果还有剩余部分且不是常数，将其作为不可约因子添加
    if (!current.is_zero() && !current.is_constant()) {
        factors.push_back({current, 1});
    }

    return factors;
}

// ============================================================================
// 私有方法
// ============================================================================

void SymbolicPolynomial::trim() {
    while (!coefficients_.empty() && coeff_is_zero(coefficients_.back())) {
        coefficients_.pop_back();
    }
}

bool SymbolicPolynomial::coeff_is_zero(const SymbolicExpression& coeff) {
    return expr_is_zero(coeff.simplify());
}

bool SymbolicPolynomial::coeff_is_one(const SymbolicExpression& coeff) {
    return expr_is_one(coeff.simplify());
}

bool SymbolicPolynomial::coeff_equals(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return expressions_match(lhs.simplify(), rhs.simplify());
}

// ============================================================================
// 运算符重载
// ============================================================================

SymbolicPolynomial operator+(const SymbolicPolynomial& lhs, const SymbolicPolynomial& rhs) {
    return lhs.add(rhs);
}

SymbolicPolynomial operator-(const SymbolicPolynomial& lhs, const SymbolicPolynomial& rhs) {
    return lhs.subtract(rhs);
}

SymbolicPolynomial operator*(const SymbolicPolynomial& lhs, const SymbolicPolynomial& rhs) {
    return lhs.multiply(rhs);
}

SymbolicPolynomial operator*(const SymbolicPolynomial& poly, const SymbolicExpression& expr) {
    return poly.scale(expr);
}

// ============================================================================
// 辅助函数
// ============================================================================

SymbolicExpression build_symbolic_polynomial_expression(
    const std::vector<SymbolicExpression>& coefficients,
    const std::string& variable_name) {

    if (coefficients.empty()) {
        return SymbolicExpression::number(0.0);
    }

    SymbolicExpression result = SymbolicExpression::number(0.0);
    SymbolicExpression x = SymbolicExpression::variable(variable_name);

    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        if (!SymbolicPolynomial::coeff_is_zero(coefficients[i])) {
            if (i == 0) {
                result = (result + coefficients[i]).simplify();
            } else if (i == 1) {
                result = (result + coefficients[i] * x).simplify();
            } else {
                result = (result + coefficients[i] *
                          make_power(x, SymbolicExpression::number(static_cast<double>(i)))).simplify();
            }
        }
    }

    return result;
}

bool solve_coefficient_identity(
    const std::vector<SymbolicExpression>& identity_coeffs,
    const std::vector<std::vector<SymbolicExpression>>& term_coeffs,
    std::vector<SymbolicExpression>* unknowns) {

    const std::size_t num_unknowns = term_coeffs.size();
    if (num_unknowns == 0) return false;

    // 找到最大次数
    std::size_t max_degree = identity_coeffs.size();
    for (const auto& term : term_coeffs) {
        max_degree = std::max(max_degree, term.size());
    }

    // 构建系数矩阵
    // 每行对应一个幂次，每列对应一个未知数
    std::vector<std::vector<SymbolicExpression>> matrix(max_degree, std::vector<SymbolicExpression>(num_unknowns, SymbolicExpression::number(0.0)));
    std::vector<SymbolicExpression> rhs(max_degree, SymbolicExpression::number(0.0));

    for (std::size_t i = 0; i < identity_coeffs.size(); ++i) {
        rhs[i] = identity_coeffs[i];
    }

    for (std::size_t j = 0; j < num_unknowns; ++j) {
        for (std::size_t i = 0; i < term_coeffs[j].size(); ++i) {
            matrix[i][j] = term_coeffs[j][i];
        }
    }

    // 高斯消元法求解符号线性方程组
    // 注意：这是一个简化版本，对于符号系数可能不完全正确
    unknowns->assign(num_unknowns, SymbolicExpression::number(0.0));

    // 对于简单情况（对角占优），直接求解
    // 这里使用一个简化的方法：假设矩阵是方阵且可解

    if (num_unknowns == 1) {
        // 单变量情况
        for (std::size_t i = 0; i < max_degree; ++i) {
            if (!SymbolicPolynomial::coeff_is_zero(matrix[i][0])) {
                (*unknowns)[0] = (rhs[i] / matrix[i][0]).simplify();
                return true;
            }
        }
        return false;
    }

    // 多变量情况：使用高斯消元
    // 这里实现一个简化版本
    std::vector<std::vector<SymbolicExpression>> aug_matrix(max_degree, std::vector<SymbolicExpression>(num_unknowns + 1));
    for (std::size_t i = 0; i < max_degree; ++i) {
        for (std::size_t j = 0; j < num_unknowns; ++j) {
            aug_matrix[i][j] = matrix[i][j];
        }
        aug_matrix[i][num_unknowns] = rhs[i];
    }

    // 前向消元
    std::size_t rank = 0;
    for (std::size_t col = 0; col < num_unknowns && rank < max_degree; ++col) {
        // 找主元
        std::size_t pivot = rank;
        for (std::size_t row = rank + 1; row < max_degree; ++row) {
            if (!SymbolicPolynomial::coeff_is_zero(aug_matrix[row][col]) &&
                SymbolicPolynomial::coeff_is_zero(aug_matrix[pivot][col])) {
                pivot = row;
            }
        }

        if (SymbolicPolynomial::coeff_is_zero(aug_matrix[pivot][col])) {
            continue;
        }

        // 交换行
        if (pivot != rank) {
            std::swap(aug_matrix[pivot], aug_matrix[rank]);
        }

        // 消元
        SymbolicExpression pivot_val = aug_matrix[rank][col];
        for (std::size_t row = rank + 1; row < max_degree; ++row) {
            if (!SymbolicPolynomial::coeff_is_zero(aug_matrix[row][col])) {
                SymbolicExpression factor = (aug_matrix[row][col] / pivot_val).simplify();
                for (std::size_t j = col; j <= num_unknowns; ++j) {
                    aug_matrix[row][j] = (aug_matrix[row][j] - factor * aug_matrix[rank][j]).simplify();
                }
            }
        }

        ++rank;
    }

    // 回代
    for (int col = static_cast<int>(num_unknowns) - 1; col >= 0; --col) {
        if (static_cast<std::size_t>(col) >= rank) continue;

        // 找到对应的行
        std::size_t row = col;
        for (std::size_t r = 0; r < rank; ++r) {
            if (!SymbolicPolynomial::coeff_is_zero(aug_matrix[r][col])) {
                // 检查这是否是第一个非零元素
                bool is_first = true;
                for (std::size_t c = 0; c < static_cast<std::size_t>(col); ++c) {
                    if (!SymbolicPolynomial::coeff_is_zero(aug_matrix[r][c])) {
                        is_first = false;
                        break;
                    }
                }
                if (is_first) {
                    row = r;
                    break;
                }
            }
        }

        if (SymbolicPolynomial::coeff_is_zero(aug_matrix[row][col])) {
            continue;
        }

        SymbolicExpression sum = aug_matrix[row][num_unknowns];
        for (std::size_t j = static_cast<std::size_t>(col) + 1; j < num_unknowns; ++j) {
            sum = (sum - aug_matrix[row][j] * (*unknowns)[j]).simplify();
        }
        (*unknowns)[col] = (sum / aug_matrix[row][col]).simplify();
    }

    return true;
}

// ============================================================================
// 部分分式分解
// ============================================================================

bool partial_fraction_decomposition(
    const SymbolicPolynomial& numerator,
    const std::vector<std::pair<SymbolicPolynomial, int>>& denominator_factors,
    const std::string& variable_name,
    std::vector<std::pair<SymbolicExpression, SymbolicPolynomial>>* partial_fractions) {

    partial_fractions->clear();

    if (denominator_factors.empty()) {
        return false;
    }

    // 计算分母多项式
    SymbolicPolynomial denominator({SymbolicExpression::number(1.0)}, variable_name);
    for (const auto& [factor, power] : denominator_factors) {
        SymbolicPolynomial factor_power = factor.power(power);
        denominator = denominator.multiply(factor_power);
    }

    // 检查分子次数是否小于分母次数
    // 如果不是，先进行多项式除法
    SymbolicPolynomial proper_num = numerator;
    if (numerator.degree() >= denominator.degree()) {
        SymbolicPolynomial quotient, remainder;
        if (!numerator.divide(denominator, &quotient, &remainder)) {
            return false;
        }
        proper_num = remainder;
        // 商的部分可以作为多项式项添加
        // 这里我们只处理真分式部分
    }

    // 构建部分分式的未知数
    // 对于每个因子 (x - r)^k，我们有 k 个项：A1/(x-r), A2/(x-r)^2, ..., Ak/(x-r)^k
    // 对于二次因子 (ax^2 + bx + c)^k，我们有 k 个项：(B1x + C1)/(ax^2+bx+c), ...

    std::vector<SymbolicPolynomial> term_denominators;
    std::vector<std::vector<SymbolicExpression>> term_coefficients;

    for (const auto& [factor, power] : denominator_factors) {
        for (int p = 1; p <= power; ++p) {
            SymbolicPolynomial term_denom = factor.power(p);
            term_denominators.push_back(term_denom);

            // 确定分子的形式
            if (factor.degree() == 1) {
                // 线性因子：分子是常数
                term_coefficients.push_back({SymbolicExpression::number(1.0)});
            } else if (factor.degree() == 2) {
                // 二次因子：分子是线性式 Bx + C
                term_coefficients.push_back({
                    SymbolicExpression::variable(variable_name),
                    SymbolicExpression::number(1.0)
                });
            } else {
                // 更高次因子：分子是 (degree-1) 次多项式
                std::vector<SymbolicExpression> coeffs;
                for (int d = 0; d < factor.degree(); ++d) {
                    if (d == 0) {
                        coeffs.push_back(SymbolicExpression::number(1.0));
                    } else {
                        coeffs.push_back(SymbolicExpression::variable(variable_name).power(
                            SymbolicExpression::number(static_cast<double>(d))));
                    }
                }
                term_coefficients.push_back(coeffs);
            }
        }
    }

    // 计算未知数总数
    int num_unknowns = 0;
    for (const auto& coeffs : term_coefficients) {
        num_unknowns += static_cast<int>(coeffs.size());
    }

    if (num_unknowns == 0) {
        return false;
    }

    // 构建恒等式：proper_num = Σ (unknown_i * term_i)
    // 其中 term_i = coeff_j / denominator_j

    // 计算每个未知数对应的项的系数
    std::vector<std::vector<SymbolicExpression>> identity_terms;
    SymbolicPolynomial product_denom({SymbolicExpression::number(1.0)}, variable_name);

    // 计算公分母（所有项的分母的乘积）
    for (const auto& term_denom : term_denominators) {
        product_denom = product_denom.multiply(term_denom);
    }

    // 对每个未知数，计算其对应的项在通分后的分子
    int unknown_idx = 0;
    for (std::size_t i = 0; i < term_denominators.size(); ++i) {
        const auto& term_denom = term_denominators[i];
        const auto& coeffs = term_coefficients[i];

        // 计算补充分母（公分母除以当前项的分母）
        SymbolicPolynomial complement_denom;
        SymbolicPolynomial remainder;
        if (!product_denom.divide(term_denom, &complement_denom, &remainder)) {
            continue;
        }

        for (std::size_t j = 0; j < coeffs.size(); ++j) {
            // 项 = coeff_j * complement_denom
            SymbolicPolynomial term_poly;
            if (coeffs[j].is_number(nullptr)) {
                term_poly = complement_denom.scale(coeffs[j]);
            } else {
                // 符号系数，需要乘法
                std::vector<SymbolicExpression> term_coeffs;
                for (const auto& c : complement_denom.coefficients()) {
                    term_coeffs.push_back((c * coeffs[j]).simplify());
                }
                term_poly = SymbolicPolynomial(term_coeffs, variable_name);
            }

            // 将多项式系数添加到恒等式矩阵
            std::vector<SymbolicExpression> poly_coeffs(product_denom.degree() + 1, SymbolicExpression::number(0.0));
            for (int d = 0; d <= term_poly.degree(); ++d) {
                if (d < static_cast<int>(poly_coeffs.size())) {
                    poly_coeffs[d] = term_poly.coefficient(d);
                }
            }
            identity_terms.push_back(poly_coeffs);
            unknown_idx++;
        }
    }

    // 构建右边（proper_num 的系数）
    std::vector<SymbolicExpression> rhs(product_denom.degree() + 1, SymbolicExpression::number(0.0));
    for (int d = 0; d <= proper_num.degree(); ++d) {
        if (d < static_cast<int>(rhs.size())) {
            rhs[d] = proper_num.coefficient(d);
        }
    }

    // 求解线性方程组
    std::vector<SymbolicExpression> unknowns;
    if (!solve_coefficient_identity(rhs, identity_terms, &unknowns)) {
        return false;
    }

    // 构建部分分式结果
    unknown_idx = 0;
    for (std::size_t i = 0; i < term_denominators.size(); ++i) {
        const auto& coeffs = term_coefficients[i];
        const auto& term_denom = term_denominators[i];

        // 构建分子
        SymbolicExpression numer = SymbolicExpression::number(0.0);
        for (std::size_t j = 0; j < coeffs.size() && unknown_idx < unknowns.size(); ++j) {
            numer = (numer + unknowns[unknown_idx] * coeffs[j]).simplify();
            unknown_idx++;
        }

        partial_fractions->push_back({numer, term_denom});
    }

    return true;
}

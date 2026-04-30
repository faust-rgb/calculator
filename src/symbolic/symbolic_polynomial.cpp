// ============================================================================
// 符号系数多项式实现
// ============================================================================

#include "symbolic_polynomial.h"
#include "symbolic_expression_internal.h"

#include <algorithm>

using namespace symbolic_expression_internal;

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

SymbolicPolynomial SymbolicPolynomial::gcd(const SymbolicPolynomial& other) const {
    // 欧几里得算法
    // 注意：对于符号系数，这可能不精确

    SymbolicPolynomial a = *this;
    SymbolicPolynomial b = other;

    while (!b.is_zero()) {
        SymbolicPolynomial q, r;
        if (!a.divide(b, &q, &r)) {
            break;
        }
        a = b;
        b = r;

        // 防止无限循环
        if (b.degree() >= a.degree() && !b.is_zero()) {
            // 可能是符号系数导致的无法整除
            break;
        }
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

    // Yun 算法
    SymbolicPolynomial f = *this;
    SymbolicPolynomial fp = f.derivative();

    SymbolicPolynomial g;
    if (!f.gcd(fp).is_zero()) {
        g = f.gcd(fp);
    } else {
        g = SymbolicPolynomial({SymbolicExpression::number(1.0)}, variable_name_);
    }

    // f = g * h, 其中 h 是 square-free 的
    SymbolicPolynomial h, r;
    if (!f.divide(g, &h, &r) || !r.is_zero()) {
        // 无法整除，可能不是完全 square-free
        factors->push_back(f);
        return true;
    }

    // 递归分解 g
    if (!g.is_constant() && g.degree() > 0) {
        std::vector<SymbolicPolynomial> sub_factors;
        g.square_free_decomposition(&sub_factors);
        factors->insert(factors->end(), sub_factors.begin(), sub_factors.end());
    }

    if (!h.is_constant() && h.degree() > 0) {
        factors->push_back(h);
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

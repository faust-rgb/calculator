/**
 * @file polynomial.cpp
 * @brief 多项式运算实现
 */

#include "polynomial.h"

#include "types/scalar_type.h"
#include "core/services/format_utils.h"
#include "matrix.h"
#include "mymath.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace {

using Scalar = mymath::Scalar;

// ============================================================================
// 内部辅助函数
// ============================================================================

/** @brief 多项式运算的数值精度阈值 */
const Scalar kPolynomialEpsLD = Scalar(1e-10L);

/** @brief 根隔离过程的精度阈值 */
const Scalar kRootIsolationEpsLD = Scalar(1e-8L);

/**
 * @brief 移除系数向量末尾的零系数
 * @param coefficients 系数向量（会被原地修改）
 *
 * 多项式运算后可能产生尾随零，此函数用于规范化结果。
 * 注意：零多项式保留为 [0.0L]。
 */
void trim_trailing_zeros(std::vector<Scalar>* coefficients) {
    while (coefficients->size() > 1 &&
           mymath::is_near_zero(coefficients->back(), kPolynomialEpsLD)) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(0.0L);
    }
}



std::vector<Scalar> interpolate_polynomial(const std::vector<Scalar>& x_samples,
                                           const std::vector<Scalar>& y_samples) {
    std::vector<Scalar> coefficients(1, Scalar(0.0L));
    for (std::size_t i = 0; i < x_samples.size(); ++i) {
        std::vector<Scalar> basis(1, Scalar(1.0L));
        Scalar denominator = Scalar(1.0L);
        for (std::size_t j = 0; j < x_samples.size(); ++j) {
            if (i == j) continue;
            const Scalar delta = x_samples[i] - x_samples[j];
            if (mymath::is_near_zero(delta, kPolynomialEpsLD)) {
                throw std::runtime_error("polynomial_fit requires distinct x samples");
            }
            basis = polynomial_multiply(basis, {-x_samples[j], Scalar(1.0L)});
            denominator *= delta;
        }
        const Scalar scale = y_samples[i] / denominator;
        for (Scalar& coefficient : basis) coefficient *= scale;
        coefficients = polynomial_add(coefficients, basis);
    }
    trim_trailing_zeros(&coefficients);
    return coefficients;
}



/**
 * @brief 格式化单个系数为字符串
 * @param value 系数值
 * @return 格式化后的字符串
 *
 * 处理整数、分数和小数三种情况：
 * - 整数：直接显示整数形式
 * - 可近似为分数：显示分数形式
 * - 其他：显示小数形式（去除尾部多余的零）
 */
std::string format_coefficient(Scalar value) {
    const Scalar kEps = Scalar(1e-10L);
    if (mymath::is_integer(value, kEps)) {
        long long rounded =
            static_cast<long long>(mymath::round(value));
        return std::to_string(rounded);
    }

    long long numerator = 0;
    long long denominator = 1;
    if (mymath::approximate_fraction(value,
                                     &numerator,
                                     &denominator,
                                     999,
                                     kEps)) {
        if (denominator == 1) {
            return std::to_string(numerator);
        }
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }

    std::string text = format_decimal(value);
    while (!text.empty() && text.back() == '0') {
        text.pop_back();
    }
    if (!text.empty() && text.back() == '.') {
        text.pop_back();
    }
    return text;
}

/**
 * @brief 使用 Horner 法计算多项式值（Scalar版本）
 * @param coefficients 多项式系数（低次到高次，Scalar）
 * @param x 求值点
 * @return p(x) 的值
 */
Scalar polynomial_evaluate_scalar(const std::vector<Scalar>& coefficients, Scalar x) {
    Scalar result = Scalar(0.0L);
    for (std::size_t i = coefficients.size(); i > 0; --i) {
        result = result * x + coefficients[i - 1];
    }
    return result;
}

/**
 * @brief 使用 Horner 法计算多项式值（内部实现）
 * @param coefficients 多项式系数（低次到高次）
 * @param x 求值点
 * @return p(x) 的值
 */
Scalar polynomial_evaluate_impl(const std::vector<Scalar>& coefficients, Scalar x) {
    return polynomial_evaluate_scalar(coefficients, x);
}

/**
 * @brief 计算多项式导数（Scalar版本）
 * @param coefficients 原多项式系数
 * @return 导数系数
 */
std::vector<Scalar> polynomial_derivative_scalar(const std::vector<Scalar>& coefficients) {
    if (coefficients.size() <= 1) {
        return {Scalar(0.0L)};
    }
    std::vector<Scalar> derivative(coefficients.size() - 1, Scalar(0.0L));
    for (std::size_t i = 1; i < coefficients.size(); ++i) {
        derivative[i - 1] = coefficients[i] * Scalar(static_cast<long long>(i));
    }
    return derivative;
}

std::vector<Scalar> polynomial_derivative_impl(const std::vector<Scalar>& coefficients) {
    if (coefficients.size() <= 1) {
        return {Scalar(0.0L)};
    }
    std::vector<Scalar> derivative(coefficients.size() - 1, Scalar(0.0L));
    for (std::size_t i = 1; i < coefficients.size(); ++i) {
        derivative[i - 1] = coefficients[i] * Scalar(static_cast<long long>(i));
    }
    trim_trailing_zeros(&derivative);
    return derivative;
}

/**
 * @brief 计算多项式实根的柯西界
 * @param coefficients 多项式系数（Scalar）
 * @return 根的上界值
 *
 * 使用柯西界公式：所有实根的绝对值不超过 1 + max|a_i/a_n|
 * 其中 a_n 是首项系数。
 */
Scalar polynomial_root_bound(const std::vector<Scalar>& coefficients) {
    const Scalar leading = coefficients.back();
    Scalar bound = Scalar(0.0L);
    for (std::size_t i = 0; i + 1 < coefficients.size(); ++i) {
        const Scalar ratio = mymath::abs(coefficients[i] / leading);
        if (ratio > bound) bound = ratio;
    }
    return Scalar(1.0L) + bound;
}

/**
 * @brief 使用二分法在区间内求根
 * @param coefficients 多项式系数（Scalar）
 * @param left 区间左端点
 * @param right 区间右端点
 * @return 找到的根（近似值）
 *
 * 前提条件：多项式在 left 和 right 处的值符号相反。
 * 进行最多 100 次迭代，直到区间宽度或函数值足够小。
 */
Scalar bisect_root(const std::vector<Scalar>& coefficients, Scalar left, Scalar right) {
    Scalar left_value = polynomial_evaluate_scalar(coefficients, left);
    for (int i = 0; i < 100; ++i) {
        const Scalar mid = (left + right) * Scalar(0.5L);
        const Scalar mid_value = polynomial_evaluate_scalar(coefficients, mid);
        if (mymath::abs(mid_value) < Scalar(kRootIsolationEpsLD) ||
            mymath::abs(right - left) <= Scalar(kRootIsolationEpsLD)) {
            return mid;
        }
        if ((left_value < Scalar(0.0L)) == (mid_value < Scalar(0.0L))) {
            left = mid;
            left_value = mid_value;
        } else {
            right = mid;
        }
    }
    return (left + right) * Scalar(0.5L);
}

/**
 * @brief 向根列表中添加唯一的根
 * @param roots 现有根列表
 * @param candidate 候选根
 *
 * 如果候选根与已有根的差值在 1e-8 以内，则不重复添加。
 */
void add_unique_root(std::vector<Scalar>* roots, Scalar candidate) {
    for (Scalar existing : *roots) {
        if (mymath::abs(existing - candidate) <= Scalar(1e-8L)) return;
    }
    roots->push_back(candidate);
}

}  // namespace

// ============================================================================
// 公共接口实现
// ============================================================================

Scalar polynomial_evaluate(const std::vector<Scalar>& coefficients, Scalar x) {
    return polynomial_evaluate_impl(coefficients, x);
}

std::vector<Scalar> polynomial_derivative(const std::vector<Scalar>& coefficients) {
    return polynomial_derivative_impl(coefficients);
}

/**
 * @brief 多项式加法
 * @param lhs 左操作数系数
 * @param rhs 右操作数系数
 * @return 和的系数
 *
 * 对应系数相加，结果长度取两多项式长度的较大值。
 */
std::vector<Scalar> polynomial_add(const std::vector<Scalar>& lhs,
                                   const std::vector<Scalar>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<Scalar> result(size, Scalar(0.0L));
    for (std::size_t i = 0; i < lhs.size(); ++i) result[i] += lhs[i];
    for (std::size_t i = 0; i < rhs.size(); ++i) result[i] += rhs[i];
    trim_trailing_zeros(&result);
    return result;
}

/**
 * @brief 多项式减法
 * @param lhs 左操作数系数（被减数）
 * @param rhs 右操作数系数（减数）
 * @return 差的系数
 *
 * 对应系数相减，结果长度取两多项式长度的较大值。
 */
std::vector<Scalar> polynomial_subtract(const std::vector<Scalar>& lhs,
                                        const std::vector<Scalar>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<Scalar> result(size, Scalar(0.0L));
    for (std::size_t i = 0; i < lhs.size(); ++i) result[i] += lhs[i];
    for (std::size_t i = 0; i < rhs.size(); ++i) result[i] -= rhs[i];
    trim_trailing_zeros(&result);
    return result;
}

/**
 * @brief 多项式乘法
 * @param lhs 左操作数系数
 * @param rhs 右操作数系数
 * @return 积的系数
 *
 * 使用直接卷积算法，时间复杂度 O(n*m)。
 * 结果的次数为两多项式次数之和。
 */
std::vector<Scalar> polynomial_multiply(const std::vector<Scalar>& lhs,
                                        const std::vector<Scalar>& rhs) {
    std::vector<Scalar> result(lhs.size() + rhs.size() - 1, Scalar(0.0L));
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result[i + j] += lhs[i] * rhs[j];
        }
    }
    trim_trailing_zeros(&result);
    return result;
}

/**
 * @brief 多项式除法（长除法）
 * @param dividend 被除数系数
 * @param divisor 除数系数
 * @return 包含商和余的除法结果
 * @throw std::runtime_error 当除数为零多项式时抛出
 *
 * 使用标准的多项式长除法算法，迭代地消去最高次项。
 */
PolynomialDivisionResult polynomial_divide(const std::vector<Scalar>& dividend,
                                           const std::vector<Scalar>& divisor) {
    std::vector<Scalar> normalized_dividend = dividend;
    std::vector<Scalar> normalized_divisor = divisor;
    trim_trailing_zeros(&normalized_dividend);
    trim_trailing_zeros(&normalized_divisor);

    if (normalized_divisor.size() == 1 &&
        mymath::is_near_zero(normalized_divisor[0], kPolynomialEpsLD)) {
        throw std::runtime_error("polynomial divisor cannot be zero");
    }

    if (normalized_dividend.size() < normalized_divisor.size()) {
        return {{0.0L}, normalized_dividend};
    }

    std::vector<Scalar> quotient(
        normalized_dividend.size() - normalized_divisor.size() + 1, Scalar(0.0L));
    std::vector<Scalar> remainder = normalized_dividend;

    while (remainder.size() >= normalized_divisor.size()) {
        trim_trailing_zeros(&remainder);
        if (remainder.size() < normalized_divisor.size()) break;
        if (remainder.size() == 1 && mymath::abs(remainder[0]) < Scalar(kPolynomialEpsLD)) break;

        const std::size_t degree_diff = remainder.size() - normalized_divisor.size();
        const Scalar factor = remainder.back() / normalized_divisor.back();
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

/**
 * @brief 计算多项式的所有实根（内部Scalar版本）
 * @param coefficients 多项式系数（Scalar）
 * @return 按升序排列的实根列表
 */
std::vector<Scalar> polynomial_real_roots_scalar(const std::vector<Scalar>& coefficients) {
    std::vector<Scalar> normalized = coefficients;
    trim_trailing_zeros(&normalized);

    if (normalized.size() <= 1) {
        throw std::runtime_error("constant polynomial does not have isolated roots");
    }

    if (normalized.size() == 2) {
        Scalar r = -normalized[0] / normalized[1];
        return {r};
    }

    const std::vector<Scalar> derivative = polynomial_derivative_scalar(normalized);
    std::vector<Scalar> critical_points = polynomial_real_roots_scalar(derivative);
    std::sort(critical_points.begin(), critical_points.end(),
              [](const Scalar& a, const Scalar& b) { return a < b; });

    const Scalar bound = polynomial_root_bound(normalized);
    std::vector<Scalar> points;
    points.push_back(-bound);
    for (const Scalar& point : critical_points) points.push_back(point);
    points.push_back(bound);

    std::vector<Scalar> roots;
    for (const Scalar& point : critical_points) {
        if (mymath::abs(polynomial_evaluate_scalar(normalized, point)) < Scalar(kRootIsolationEpsLD)) {
            add_unique_root(&roots, point);
        }
    }

    for (std::size_t i = 1; i < points.size(); ++i) {
        const Scalar left = points[i - 1];
        const Scalar right = points[i];
        const Scalar left_value = polynomial_evaluate_scalar(normalized, left);
        const Scalar right_value = polynomial_evaluate_scalar(normalized, right);

        if (mymath::abs(left_value) < Scalar(kRootIsolationEpsLD)) {
            add_unique_root(&roots, left);
            continue;
        }
        if (mymath::abs(right_value) < Scalar(kRootIsolationEpsLD)) {
            add_unique_root(&roots, right);
            continue;
        }

        if ((left_value < Scalar(0.0L)) != (right_value < Scalar(0.0L))) {
            Scalar r = bisect_root(normalized, left, right);
            add_unique_root(&roots, r);
        }
    }

    std::sort(roots.begin(), roots.end(),
              [](const Scalar& a, const Scalar& b) { return a < b; });
    return roots;
}

/**
 * @brief 计算多项式的所有实根
 * @param coefficients 多项式系数
 * @return 按升序排列的实根列表
 * @throw std::runtime_error 当多项式为常数时抛出
 */
std::vector<Scalar> polynomial_real_roots(const std::vector<Scalar>& coefficients) {
    return polynomial_real_roots_scalar(coefficients);
}

/**
 * @brief 计算多项式的全部复根
 * @param coefficients 多项式系数
 * @return 按实部、虚部排序的复根列表
 * @throw std::runtime_error 当多项式为常数时抛出
 */
/**
 * @brief 计算多项式最大公因式（鲁棒单项式化欧几里得算法）
 * @param lhs 左多项式系数
 * @param rhs 右多项式系数
 * @return 单位首项化后的最大公因式系数
 */
std::vector<Scalar> polynomial_gcd(const std::vector<Scalar>& lhs,
                                   const std::vector<Scalar>& rhs) {
    std::vector<Scalar> a = lhs;
    std::vector<Scalar> b = rhs;
    trim_trailing_zeros(&a);
    trim_trailing_zeros(&b);

    if (a.empty() || (a.size() == 1 && mymath::is_near_zero(a[0], kPolynomialEpsLD))) {
        if (b.empty() || (b.size() == 1 && mymath::is_near_zero(b[0], kPolynomialEpsLD))) {
            return {Scalar(0.0L)};
        }
        const Scalar lead = b.back();
        for (Scalar& c : b) c /= lead;
        trim_trailing_zeros(&b);
        return b;
    }
    if (b.empty() || (b.size() == 1 && mymath::is_near_zero(b[0], kPolynomialEpsLD))) {
        const Scalar lead = a.back();
        for (Scalar& c : a) c /= lead;
        trim_trailing_zeros(&a);
        return a;
    }

    auto get_max_coeff = [](const std::vector<Scalar>& p) {
        Scalar m = Scalar(0.0L);
        for (const Scalar& c : p) m = std::max(m, mymath::abs(c));
        return m;
    };
    const Scalar scale = std::max(get_max_coeff(a), get_max_coeff(b));
    const Scalar dynamic_eps = std::max(kPolynomialEpsLD, scale * Scalar(1e-12L));

    while (true) {
        trim_trailing_zeros(&b);
        if (b.empty() || (b.size() == 1 && mymath::is_near_zero(b[0], dynamic_eps))) {
            break;
        }

        // 首一化除数以防止系数膨胀与数值漂移
        const Scalar lead_b = b.back();
        if (!mymath::is_near_zero(lead_b, dynamic_eps)) {
            for (Scalar& c : b) c /= lead_b;
        }

        const PolynomialDivisionResult division = polynomial_divide(a, b);
        a = b;
        b = division.remainder;
        trim_trailing_zeros(&a);
        trim_trailing_zeros(&b);

        // 检查余数是否全接近 0
        bool all_zero = true;
        for (const Scalar& c : b) {
            if (!mymath::is_near_zero(c, dynamic_eps)) {
                all_zero = false;
                break;
            }
        }
        if (all_zero) break;
    }

    if (a.empty() || (a.size() == 1 && mymath::is_near_zero(a[0], dynamic_eps))) {
        return {Scalar(1.0L)};
    }

    const Scalar leading = a.back();
    if (!mymath::is_near_zero(leading, dynamic_eps)) {
        for (Scalar& coefficient : a) coefficient /= leading;
    }
    trim_trailing_zeros(&a);
    return a;
}

/**
 * @brief 计算多项式的全部复根（Square-free + Aberth-Ehrlich + Newton 精炼）
 * @param coefficients 多项式系数
 * @return 按实部、虚部排序的复根列表
 * @throw std::runtime_error 当多项式为常数时抛出
 */
std::vector<mymath::complex<Scalar>> polynomial_complex_roots(
    const std::vector<Scalar>& coefficients) {
    std::vector<Scalar> normalized = coefficients;
    trim_trailing_zeros(&normalized);

    if (normalized.size() <= 1) {
        throw std::runtime_error("constant polynomial does not have isolated roots");
    }

    const std::size_t degree = normalized.size() - 1;
    if (degree == 1) {
        return {mymath::complex<Scalar>(-normalized[0] / normalized[1], Scalar(0.0L))};
    }

    // 无平方因子分解（Square-free Factorization）预处理以消除重根停滞
    if (degree >= 2) {
        const std::vector<Scalar> deriv = polynomial_derivative(normalized);
        const std::vector<Scalar> gcd_poly = polynomial_gcd(normalized, deriv);
        if (gcd_poly.size() > 1 && gcd_poly.size() - 1 < degree) {
            const PolynomialDivisionResult sff = polynomial_divide(normalized, gcd_poly);
            std::vector<Scalar> p0 = sff.quotient;
            trim_trailing_zeros(&p0);
            if (p0.size() > 1 && p0.size() - 1 < degree) {
                std::vector<mymath::complex<Scalar>> roots_p0 = polynomial_complex_roots(p0);
                std::vector<mymath::complex<Scalar>> roots_gcd = polynomial_complex_roots(gcd_poly);
                roots_p0.insert(roots_p0.end(), roots_gcd.begin(), roots_gcd.end());
                std::sort(roots_p0.begin(), roots_p0.end(), [](const auto& lhs, const auto& rhs) {
                    if (mymath::abs(lhs.real() - rhs.real()) > Scalar(1e-8L)) {
                        return lhs.real() < rhs.real();
                    }
                    return lhs.imag() < rhs.imag();
                });
                return roots_p0;
            }
        }
    }

    const Scalar leading_128(normalized.back());
    const Scalar bound_128(polynomial_root_bound(normalized));
    const Scalar radius_128 = (bound_128 > Scalar(1.0L)) ? bound_128 : Scalar(1.0L);

    std::vector<mymath::complex<Scalar>> roots_128;
    roots_128.reserve(degree);
    for (std::size_t k = 0; k < degree; ++k) {
        const Scalar angle_128 =
            Scalar(2.0L) * mymath::constants::pi<Scalar>() *
            (Scalar(static_cast<long long>(k)) + Scalar(0.25L)) /
            Scalar(static_cast<long long>(degree));
        roots_128.emplace_back(radius_128 * mymath::cos(angle_128),
                               radius_128 * mymath::sin(angle_128));
    }

    auto evaluate_complex_128 = [&](const mymath::complex<Scalar>& x) {
        mymath::complex<Scalar> result(Scalar(0.0L), Scalar(0.0L));
        for (std::size_t i = normalized.size(); i > 0; --i) {
            result = result * x;
            result = result + mymath::complex<Scalar>(normalized[i - 1] / leading_128, Scalar(0.0L));
        }
        return result;
    };

    auto evaluate_deriv_complex_128 = [&](const mymath::complex<Scalar>& x) {
        mymath::complex<Scalar> result(Scalar(0.0L), Scalar(0.0L));
        for (std::size_t i = normalized.size() - 1; i > 0; --i) {
            result = result * x;
            result = result + mymath::complex<Scalar>(normalized[i] * Scalar(static_cast<long long>(i)) / leading_128, Scalar(0.0L));
        }
        return result;
    };

    for (int iteration = 0; iteration < 2000; ++iteration) {
        Scalar max_delta_128(0.0L);
        for (std::size_t i = 0; i < roots_128.size(); ++i) {
            mymath::complex<Scalar> denominator(Scalar(1.0L), Scalar(0.0L));
            for (std::size_t j = 0; j < roots_128.size(); ++j) {
                if (i == j) continue;
                denominator = denominator * (roots_128[i] - roots_128[j]);
            }
            if (mymath::abs(denominator) <= Scalar(1e-24L)) {
                denominator = mymath::complex<Scalar>(Scalar(1e-12L), Scalar(1e-12L));
            }
            const mymath::complex<Scalar> delta = evaluate_complex_128(roots_128[i]) / denominator;
            roots_128[i] = roots_128[i] - delta;
            const Scalar abs_delta = mymath::abs(delta);
            max_delta_128 = (abs_delta > max_delta_128) ? abs_delta : max_delta_128;
        }
        if (max_delta_128 <= Scalar(1e-12L)) {
            break;
        }
    }

    // Newton-Raphson 根精炼
    for (std::size_t i = 0; i < roots_128.size(); ++i) {
        for (int iter = 0; iter < 5; ++iter) {
            const mymath::complex<Scalar> f_val = evaluate_complex_128(roots_128[i]);
            const mymath::complex<Scalar> df_val = evaluate_deriv_complex_128(roots_128[i]);
            if (mymath::abs(df_val) > Scalar(1e-14L)) {
                const mymath::complex<Scalar> step = f_val / df_val;
                roots_128[i] = roots_128[i] - step;
                if (mymath::abs(step) < Scalar(1e-14L)) break;
            }
        }
    }

    std::vector<mymath::complex<Scalar>> roots;
    roots.reserve(roots_128.size());
    for (const auto& root_128 : roots_128) {
        Scalar real = root_128.real();
        Scalar imag = root_128.imag();
        if (mymath::is_near_zero(real, 1e-9)) real = Scalar(0.0L);
        if (mymath::is_near_zero(imag, 1e-9)) imag = Scalar(0.0L);
        if (mymath::is_integer(real, 1e-9)) {
            real = mymath::round(real);
        }
        if (mymath::is_integer(imag, 1e-9)) {
            imag = mymath::round(imag);
        }
        roots.emplace_back(real, imag);
    }

    std::sort(roots.begin(), roots.end(), [](const auto& lhs, const auto& rhs) {
        if (mymath::abs(lhs.real() - rhs.real()) > Scalar(1e-8L)) {
            return lhs.real() < rhs.real();
        }
        return lhs.imag() < rhs.imag();
    });
    return roots;
}

/**
 * @brief 计算多项式积分系数
 * @param coefficients 原多项式系数
 * @return 不定积分系数（积分常数为 0）
 *
 * 使用公式：integral(a_i * x^i) = a_i * x^(i+1) / (i+1)
 */
std::vector<Scalar> polynomial_integral(const std::vector<Scalar>& coefficients) {
    std::vector<Scalar> integral(coefficients.size() + 1, Scalar(0.0L));
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        integral[i + 1] = coefficients[i] / Scalar(static_cast<long long>(i + 1));
    }
    trim_trailing_zeros(&integral);
    return integral;
}

/**
 * @brief 计算多项式复合 p(q(x))
 * @param outer 外层多项式 p 的系数
 * @param inner 内层多项式 q 的系数
 * @return 复合后的多项式系数
 *
 * 使用 Horner 法的推广形式：
 * p(q(x)) = (...((a_n * q + a_{n-1}) * q + a_{n-2}) * q + ...) + a_0
 */
std::vector<Scalar> polynomial_compose(const std::vector<Scalar>& outer,
                                       const std::vector<Scalar>& inner) {
    std::vector<Scalar> result = {Scalar(0.0L)};
    for (std::size_t i = outer.size(); i > 0; --i) {
        result = polynomial_multiply(result, inner);
        if (result.empty()) result.push_back(Scalar(0.0L));
        result[0] += outer[i - 1];
        trim_trailing_zeros(&result);
    }
    return result;
}

/**
 * @brief 使用最小二乘法进行多项式拟合
 * @param x_samples x 坐标样本点
 * @param y_samples y 坐标样本点
 * @param degree 拟合多项式的次数
 * @return 拟合多项式的系数（低次到高次）
 * @throw std::runtime_error 当参数无效时抛出
 *
 * 算法步骤：
 * 1. 对 x 坐标进行中心化和缩放，改善数值稳定性
 * 2. 构建范德蒙德矩阵
 * 3. 使用最小二乘法求解
 * 4. 将结果转换回原始坐标系
 */
std::vector<Scalar> polynomial_fit(const std::vector<Scalar>& x_samples,
                                   const std::vector<Scalar>& y_samples,
                                   int degree) {
    if (degree < 0) throw std::runtime_error("polynomial degree must be non-negative");
    if (x_samples.size() != y_samples.size() || x_samples.empty()) {
        throw std::runtime_error("polynomial_fit requires non-empty sample vectors of the same length");
    }
    const std::size_t m_vars = static_cast<std::size_t>(degree + 1);
    if (x_samples.size() < m_vars) throw std::runtime_error("polynomial_fit requires at least degree + 1 samples");

    const std::size_t n = x_samples.size();
    Scalar x_sum = Scalar(0.0L);
    for (Scalar x : x_samples) x_sum += x;
    const Scalar center = x_sum / Scalar(static_cast<long long>(n));
    Scalar scale = Scalar(0.0L);
    for (Scalar x : x_samples) {
        const Scalar d = mymath::abs(x - center);
        if (d > scale) scale = d;
    }
    if (scale < Scalar(1e-9L)) scale = Scalar(1.0L);

    matrix::Matrix A(n, m_vars);
    std::vector<Scalar> scaled_x(n);
    for (std::size_t i = 0; i < n; ++i) {
        const Scalar sx = (x_samples[i] - center) / scale;
        scaled_x[i] = sx;
        Scalar p = Scalar(1.0L);
        for (std::size_t j = 0; j < m_vars; ++j) {
            A.at(i, j) = p;
            p *= sx;
        }
    }
    std::vector<Scalar> y_samples_scalar(y_samples.size());
    for (std::size_t i = 0; i < y_samples.size(); ++i) {
        y_samples_scalar[i] = y_samples[i];
    }
    matrix::Matrix b = matrix::Matrix::vector(y_samples_scalar);
    try {
        std::vector<Scalar> scaled_coeffs(m_vars);
        if (n == m_vars) {
            scaled_coeffs = interpolate_polynomial(scaled_x, y_samples_scalar);
        } else {
            matrix::Matrix solution = matrix::least_squares(A, b);
            for (std::size_t i = 0; i < m_vars; ++i) scaled_coeffs[i] = solution.at(i, 0);
        }
        const std::vector<Scalar> linear_map = {-center / scale, 1.0L / scale};
        std::vector<Scalar> coefficients = polynomial_compose(scaled_coeffs, linear_map);
        trim_trailing_zeros(&coefficients);
        return coefficients;
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("polynomial_fit failed: ") + e.what());
    }
}

/**
 * @brief 将多项式转换为可读字符串
 * @param coefficients 多项式系数（低次到高次）
 * @param variable_name 变量名
 * @return 格式化后的多项式字符串
 *
 * 自动处理：
 * - 系数为 1 或 -1 时的简化（如 "x" 而非 "1 * x"）
 * - 整数系数的格式化
 * - 符号连接（使用 "+" 和 "-"）
 */
std::string polynomial_to_string(const std::vector<Scalar>& coefficients,
                                 const std::string& variable_name) {
    std::vector<Scalar> normalized = coefficients;
    trim_trailing_zeros(&normalized);
    if (normalized.size() == 1 && mymath::is_near_zero(normalized[0], kPolynomialEpsLD)) return "0";
    std::string result;
    bool first = true;
    for (std::size_t index = normalized.size(); index > 0; --index) {
        const std::size_t degree = index - 1;
        const Scalar coefficient = normalized[degree];
        if (mymath::is_near_zero(coefficient, kPolynomialEpsLD)) continue;
        const bool negative = coefficient < 0.0L;
        const Scalar abs_value = negative ? -coefficient : coefficient;
        std::string term;
        if (degree == 0) term = format_coefficient(abs_value);
        else {
            if (!mymath::is_near_zero(abs_value - 1.0L, kPolynomialEpsLD)) term += format_coefficient(abs_value) + " * ";
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

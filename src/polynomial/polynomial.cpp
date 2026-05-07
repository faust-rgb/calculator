/**
 * @file polynomial.cpp
 * @brief 多项式运算实现
 */

#include "polynomial.h"

#include "matrix.h"
#include "mymath.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace {

// ============================================================================
// 内部辅助函数
// ============================================================================

/** @brief 多项式运算的数值精度阈值 */
constexpr long double kPolynomialEps = 1e-10;

/** @brief 根隔离过程的精度阈值 */
constexpr long double kRootIsolationEps = 1e-8;

/**
 * @brief 移除系数向量末尾的零系数
 * @param coefficients 系数向量（会被原地修改）
 *
 * 多项式运算后可能产生尾随零，此函数用于规范化结果。
 * 注意：零多项式保留为 [0.0L]。
 */
void trim_trailing_zeros(std::vector<long double>* coefficients) {
    while (coefficients->size() > 1 &&
           mymath::is_near_zero(coefficients->back(), kPolynomialEps)) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(0.0L);
    }
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
std::string format_coefficient(long double value) {
    if (mymath::is_integer(value, 1e-10)) {
        long long rounded =
            static_cast<long long>(value >= 0.0L ? value + 0.5 : value - 0.5);
        return std::to_string(rounded);
    }

    long long numerator = 0;
    long long denominator = 1;
    if (mymath::approximate_fraction(value,
                                     &numerator,
                                     &denominator,
                                     999,
                                     1e-10)) {
        if (denominator == 1) {
            return std::to_string(numerator);
        }
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }

    std::string text = std::to_string(value);
    while (!text.empty() && text.back() == '0') {
        text.pop_back();
    }
    if (!text.empty() && text.back() == '.') {
        text.pop_back();
    }
    return text;
}

/**
 * @brief 使用 Horner 法计算多项式值（内部实现）
 * @param coefficients 多项式系数（低次到高次）
 * @param x 求值点
 * @return p(x) 的值
 *
 * Horner 法将 p(x) = a_0 + a_1*x + ... + a_n*x^n 重写为：
 * p(x) = (...((a_n * x + a_{n-1}) * x + a_{n-2}) * x + ...) + a_0
 * 这种形式减少了乘法次数，提高计算效率。
 */
long double polynomial_evaluate_impl(const std::vector<long double>& coefficients, long double x) {
    long double result = 0.0L;
    for (std::size_t i = coefficients.size(); i > 0; --i) {
        result = result * x + coefficients[i - 1];
    }
    return result;
}

/**
 * @brief 计算多项式导数（内部实现）
 * @param coefficients 原多项式系数
 * @return 导数系数
 *
 * 使用公式：d/dx(a_i * x^i) = i * a_i * x^(i-1)
 * 常数项的导数为零。
 */
std::vector<long double> polynomial_derivative_impl(const std::vector<long double>& coefficients) {
    if (coefficients.size() <= 1) {
        return {0.0L};
    }
    std::vector<long double> derivative(coefficients.size() - 1, 0.0L);
    for (std::size_t i = 1; i < coefficients.size(); ++i) {
        derivative[i - 1] = coefficients[i] * static_cast<long double>(i);
    }
    trim_trailing_zeros(&derivative);
    return derivative;
}

/**
 * @brief 计算多项式实根的柯西界
 * @param coefficients 多项式系数
 * @return 根的上界值
 *
 * 使用柯西界公式：所有实根的绝对值不超过 1 + max|a_i/a_n|
 * 其中 a_n 是首项系数。
 */
long double polynomial_root_bound(const std::vector<long double>& coefficients) {
    const long double leading = coefficients.back();
    long double bound = 0.0L;
    for (std::size_t i = 0; i + 1 < coefficients.size(); ++i) {
        const long double ratio = mymath::abs(coefficients[i] / leading);
        if (ratio > bound) bound = ratio;
    }
    return 1.0L + bound;
}

/**
 * @brief 使用二分法在区间内求根
 * @param coefficients 多项式系数
 * @param left 区间左端点
 * @param right 区间右端点
 * @return 找到的根（近似值）
 *
 * 前提条件：多项式在 left 和 right 处的值符号相反。
 * 进行最多 100 次迭代，直到区间宽度或函数值足够小。
 */
long double bisect_root(const std::vector<long double>& coefficients, long double left, long double right) {
    long double left_value = polynomial_evaluate_impl(coefficients, left);
    for (int i = 0; i < 100; ++i) {
        const long double mid = (left + right) * 0.5;
        const long double mid_value = polynomial_evaluate_impl(coefficients, mid);
        if (mymath::is_near_zero(mid_value, kRootIsolationEps) ||
            mymath::abs(right - left) <= kRootIsolationEps) {
            return mid;
        }
        if ((left_value < 0.0L) == (mid_value < 0.0L)) {
            left = mid;
            left_value = mid_value;
        } else {
            right = mid;
        }
    }
    return (left + right) * 0.5;
}

/**
 * @brief 向根列表中添加唯一的根
 * @param roots 现有根列表
 * @param candidate 候选根
 *
 * 如果候选根与已有根的差值在 1e-6 以内，则不重复添加。
 */
void add_unique_root(std::vector<long double>* roots, long double candidate) {
    for (long double existing : *roots) {
        if (mymath::abs(existing - candidate) <= 1e-6) return;
    }
    roots->push_back(candidate);
}

}  // namespace

// ============================================================================
// 公共接口实现
// ============================================================================

long double polynomial_evaluate(const std::vector<long double>& coefficients, long double x) {
    return polynomial_evaluate_impl(coefficients, x);
}

std::vector<long double> polynomial_derivative(const std::vector<long double>& coefficients) {
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
std::vector<long double> polynomial_add(const std::vector<long double>& lhs,
                                   const std::vector<long double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<long double> result(size, 0.0L);
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
std::vector<long double> polynomial_subtract(const std::vector<long double>& lhs,
                                        const std::vector<long double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<long double> result(size, 0.0L);
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
 * 使用 float128 进行中间计算以减少累积舍入误差。
 */
std::vector<long double> polynomial_multiply(const std::vector<long double>& lhs,
                                        const std::vector<long double>& rhs) {
    std::vector<mymath::float128_t> result_128(lhs.size() + rhs.size() - 1, mymath::float128_t(0.0L));
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result_128[i + j] = result_128[i + j] + mymath::float128_t(lhs[i]) * mymath::float128_t(rhs[j]);
        }
    }
    std::vector<long double> result(result_128.size());
    for (std::size_t i = 0; i < result_128.size(); ++i) {
        result[i] = result_128[i].to_long_double();
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
 * 使用 float128 进行中间计算以减少累积舍入误差。
 */
PolynomialDivisionResult polynomial_divide(const std::vector<long double>& dividend,
                                           const std::vector<long double>& divisor) {
    std::vector<long double> normalized_dividend = dividend;
    std::vector<long double> normalized_divisor = divisor;
    trim_trailing_zeros(&normalized_dividend);
    trim_trailing_zeros(&normalized_divisor);

    if (normalized_divisor.size() == 1 &&
        mymath::is_near_zero(normalized_divisor[0], kPolynomialEps)) {
        throw std::runtime_error("polynomial divisor cannot be zero");
    }

    if (normalized_dividend.size() < normalized_divisor.size()) {
        return {{0.0L}, normalized_dividend};
    }

    std::vector<mymath::float128_t> quotient_128(
        normalized_dividend.size() - normalized_divisor.size() + 1, mymath::float128_t(0.0L));
    std::vector<mymath::float128_t> remainder_128(normalized_dividend.size());
    for (std::size_t i = 0; i < normalized_dividend.size(); ++i) {
        remainder_128[i] = mymath::float128_t(normalized_dividend[i]);
    }
    std::vector<mymath::float128_t> divisor_128(normalized_divisor.size());
    for (std::size_t i = 0; i < normalized_divisor.size(); ++i) {
        divisor_128[i] = mymath::float128_t(normalized_divisor[i]);
    }

    auto trim_128 = [](std::vector<mymath::float128_t>* coeffs) {
        while (coeffs->size() > 1 &&
               mymath::precise128::abs(coeffs->back()) < mymath::float128_t(kPolynomialEps)) {
            coeffs->pop_back();
        }
        if (coeffs->empty()) {
            coeffs->push_back(mymath::float128_t(0.0L));
        }
    };

    while (remainder_128.size() >= normalized_divisor.size()) {
        trim_128(&remainder_128);
        if (remainder_128.size() < normalized_divisor.size()) break;
        if (remainder_128.size() == 1 && mymath::precise128::abs(remainder_128[0]) < mymath::float128_t(kPolynomialEps)) break;

        const std::size_t degree_diff = remainder_128.size() - normalized_divisor.size();
        const mymath::float128_t factor = remainder_128.back() / divisor_128.back();
        quotient_128[degree_diff] = factor;
        for (std::size_t i = 0; i < divisor_128.size(); ++i) {
            remainder_128[degree_diff + i] = remainder_128[degree_diff + i] - factor * divisor_128[i];
        }
        trim_128(&remainder_128);
    }

    trim_128(&quotient_128);
    trim_128(&remainder_128);

    std::vector<long double> quotient(quotient_128.size());
    for (std::size_t i = 0; i < quotient_128.size(); ++i) {
        quotient[i] = quotient_128[i].to_long_double();
    }
    std::vector<long double> remainder(remainder_128.size());
    for (std::size_t i = 0; i < remainder_128.size(); ++i) {
        remainder[i] = remainder_128[i].to_long_double();
    }

    trim_trailing_zeros(&quotient);
    trim_trailing_zeros(&remainder);
    return {quotient, remainder};
}

/**
 * @brief 计算多项式的所有实根
 * @param coefficients 多项式系数
 * @return 按升序排列的实根列表
 * @throw std::runtime_error 当多项式为常数时抛出
 *
 * 算法步骤：
 * 1. 计算导数的实根（临界点）
 * 2. 使用柯西根界定出搜索区间
 * 3. 用临界点将区间分段，每段内多项式单调
 * 4. 在符号变化的区间使用二分法求根
 */
std::vector<long double> polynomial_real_roots(const std::vector<long double>& coefficients) {
    std::vector<long double> normalized = coefficients;
    trim_trailing_zeros(&normalized);

    if (normalized.size() <= 1) {
        throw std::runtime_error("constant polynomial does not have isolated roots");
    }

    if (normalized.size() == 2) {
        long double r = -normalized[0] / normalized[1];
        long long num, den;
        if (mymath::approximate_fraction(r, &num, &den, 100, 1e-9)) {
            r = static_cast<long double>(num) / static_cast<long double>(den);
        }
        return {r};
    }

    const std::vector<long double> derivative = polynomial_derivative(normalized);
    std::vector<long double> critical_points = polynomial_real_roots(derivative);
    std::sort(critical_points.begin(), critical_points.end());

    const long double bound = polynomial_root_bound(normalized);
    std::vector<long double> points;
    points.push_back(-bound);
    for (long double point : critical_points) points.push_back(point);
    points.push_back(bound);

    std::vector<long double> roots;
    for (long double point : critical_points) {
        if (mymath::is_near_zero(polynomial_evaluate(normalized, point), 1e-6)) {
            add_unique_root(&roots, point);
        }
    }

    for (std::size_t i = 1; i < points.size(); ++i) {
        const long double left = points[i - 1];
        const long double right = points[i];
        const long double left_value = polynomial_evaluate(normalized, left);
        const long double right_value = polynomial_evaluate(normalized, right);

        if (mymath::is_near_zero(left_value, 1e-6)) {
            add_unique_root(&roots, left);
            continue;
        }
        if (mymath::is_near_zero(right_value, 1e-6)) {
            add_unique_root(&roots, right);
            continue;
        }

        if ((left_value < 0.0L) != (right_value < 0.0L)) {
            long double r = bisect_root(normalized, left, right);
            long long num, den;
            if (mymath::approximate_fraction(r, &num, &den, 100, 1e-9)) {
                r = static_cast<long double>(num) / static_cast<long double>(den);
            }
            add_unique_root(&roots, r);
        }
    }

    std::sort(roots.begin(), roots.end());
    return roots;
}

/**
 * @brief 计算多项式的全部复根
 * @param coefficients 多项式系数
 * @return 按实部、虚部排序的复根列表
 * @throw std::runtime_error 当多项式为常数时抛出
 *
 * 使用 Aberth-Ehrlich 方法：
 * 1. 在半径为根界圆上均匀分布初始猜测
 * 2. 迭代更新每个根的估计值
 * 3. 收敛后清理接近整数的实部和虚部
 * 使用 float128 进行中间计算以减少累积舍入误差。
 */
std::vector<mymath::complex<long double>> polynomial_complex_roots(
    const std::vector<long double>& coefficients) {
    std::vector<long double> normalized = coefficients;
    trim_trailing_zeros(&normalized);

    if (normalized.size() <= 1) {
        throw std::runtime_error("constant polynomial does not have isolated roots");
    }

    const std::size_t degree = normalized.size() - 1;
    if (degree == 1) {
        return {mymath::complex<long double>(-normalized[0] / normalized[1], 0.0L)};
    }

    const mymath::float128_t leading_128(normalized.back());
    const mymath::float128_t bound_128(polynomial_root_bound(normalized));
    const mymath::float128_t radius_128 = (bound_128 > mymath::float128_t(1.0L)) ? bound_128 : mymath::float128_t(1.0L);

    std::vector<mymath::complex<mymath::float128_t>> roots_128;
    roots_128.reserve(degree);
    for (std::size_t k = 0; k < degree; ++k) {
        const mymath::float128_t angle_128 =
            mymath::float128_t(2.0L) * mymath::precise128::pi() *
            (mymath::float128_t(static_cast<long double>(k)) + mymath::float128_t(0.25L)) /
            mymath::float128_t(static_cast<long double>(degree));
        roots_128.emplace_back(radius_128 * mymath::precise128::cos(angle_128),
                               radius_128 * mymath::precise128::sin(angle_128));
    }

    auto evaluate_complex_128 = [&](const mymath::complex<mymath::float128_t>& x) {
        mymath::complex<mymath::float128_t> result(mymath::float128_t(0.0L), mymath::float128_t(0.0L));
        for (std::size_t i = normalized.size(); i > 0; --i) {
            result = result * x;
            result = result + mymath::complex<mymath::float128_t>(mymath::float128_t(normalized[i - 1]) / leading_128, mymath::float128_t(0.0L));
        }
        return result;
    };

    for (int iteration = 0; iteration < 2000; ++iteration) {
        mymath::float128_t max_delta_128(0.0L);
        for (std::size_t i = 0; i < roots_128.size(); ++i) {
            mymath::complex<mymath::float128_t> denominator(mymath::float128_t(1.0L), mymath::float128_t(0.0L));
            for (std::size_t j = 0; j < roots_128.size(); ++j) {
                if (i == j) continue;
                denominator = denominator * (roots_128[i] - roots_128[j]);
            }
            if (mymath::abs(denominator) <= mymath::float128_t(1e-24L)) {
                denominator = mymath::complex<mymath::float128_t>(mymath::float128_t(1e-12L), mymath::float128_t(1e-12L));
            }
            const mymath::complex<mymath::float128_t> delta = evaluate_complex_128(roots_128[i]) / denominator;
            roots_128[i] = roots_128[i] - delta;
            const mymath::float128_t abs_delta = mymath::abs(delta);
            max_delta_128 = (abs_delta > max_delta_128) ? abs_delta : max_delta_128;
        }
        if (max_delta_128 <= mymath::float128_t(1e-12L)) {
            break;
        }
    }

    std::vector<mymath::complex<long double>> roots;
    roots.reserve(roots_128.size());
    for (const auto& root_128 : roots_128) {
        long double real = root_128.real().to_long_double();
        long double imag = root_128.imag().to_long_double();
        if (mymath::is_near_zero(real, 1e-9)) real = 0.0L;
        if (mymath::is_near_zero(imag, 1e-9)) imag = 0.0L;
        if (mymath::is_integer(real)) {
            real = mymath::round(real);
        }
        if (mymath::is_integer(imag)) {
            imag = mymath::round(imag);
        }
        roots.emplace_back(real, imag);
    }

    std::sort(roots.begin(), roots.end(), [](const auto& lhs, const auto& rhs) {
        if (mymath::abs(lhs.real() - rhs.real()) > 1e-8) {
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
std::vector<long double> polynomial_integral(const std::vector<long double>& coefficients) {
    std::vector<long double> integral(coefficients.size() + 1, 0.0L);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        integral[i + 1] = coefficients[i] / static_cast<long double>(i + 1);
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
std::vector<long double> polynomial_compose(const std::vector<long double>& outer,
                                       const std::vector<long double>& inner) {
    std::vector<long double> result = {0.0L};
    for (std::size_t i = outer.size(); i > 0; --i) {
        result = polynomial_multiply(result, inner);
        if (result.empty()) result.push_back(0.0L);
        result[0] += outer[i - 1];
        trim_trailing_zeros(&result);
    }
    return result;
}

/**
 * @brief 计算多项式最大公因式（欧几里得算法）
 * @param lhs 左多项式系数
 * @param rhs 右多项式系数
 * @return 单位首项化后的最大公因式系数
 *
 * 使用多项式版本的欧几里得算法：
 * gcd(a, b) = gcd(b, a mod b)
 *
 * 采用动态容差以处理不同量级的系数。
 */
std::vector<long double> polynomial_gcd(const std::vector<long double>& lhs,
                                   const std::vector<long double>& rhs) {
    std::vector<long double> a = lhs;
    std::vector<long double> b = rhs;
    trim_trailing_zeros(&a);
    trim_trailing_zeros(&b);
    
    // 动态容差：根据输入多项式系数的最大绝对值决定
    auto get_max_coeff = [](const std::vector<long double>& p) {
        long double m = 0.0L;
        for (long double c : p) m = std::max(m, mymath::abs(c));
        return m;
    };
    const long double scale = std::max(get_max_coeff(a), get_max_coeff(b));
    const long double dynamic_eps = std::max(kPolynomialEps, scale * 1e-12);

    while (!(b.size() == 1 && mymath::is_near_zero(b[0], dynamic_eps))) {
        const PolynomialDivisionResult division = polynomial_divide(a, b);
        a = b;
        b = division.remainder;
        trim_trailing_zeros(&a);
        trim_trailing_zeros(&b);
    }
    
    if (a.empty() || (a.size() == 1 && mymath::is_near_zero(a[0], dynamic_eps))) {
        return {0.0L};
    }

    const long double leading = a.back();
    if (!mymath::is_near_zero(leading, dynamic_eps)) {
        for (long double& coefficient : a) coefficient /= leading;
    }
    trim_trailing_zeros(&a);
    return a;
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
std::vector<long double> polynomial_fit(const std::vector<long double>& x_samples,
                                   const std::vector<long double>& y_samples,
                                   int degree) {
    if (degree < 0) throw std::runtime_error("polynomial degree must be non-negative");
    if (x_samples.size() != y_samples.size() || x_samples.empty()) {
        throw std::runtime_error("polynomial_fit requires non-empty sample vectors of the same length");
    }
    const std::size_t m_vars = static_cast<std::size_t>(degree + 1);
    if (x_samples.size() < m_vars) throw std::runtime_error("polynomial_fit requires at least degree + 1 samples");

    const std::size_t n = x_samples.size();
    long double x_sum = 0.0L;
    for (long double x : x_samples) x_sum += static_cast<long double>(x);
    const long double center = static_cast<long double>(x_sum / static_cast<long double>(n));
    long double scale = 0.0L;
    for (long double x : x_samples) {
        const long double d = mymath::abs(x - center);
        if (d > scale) scale = d;
    }
    if (scale < 1e-9) scale = 1.0L;

    matrix::Matrix A(n, m_vars);
    for (std::size_t i = 0; i < n; ++i) {
        const long double sx = (x_samples[i] - center) / scale;
        long double p = 1.0L;
        for (std::size_t j = 0; j < m_vars; ++j) {
            A.at(i, j) = p;
            p *= sx;
        }
    }
    matrix::Matrix b = matrix::Matrix::vector(y_samples);
    try {
        matrix::Matrix solution = matrix::least_squares(A, b);
        std::vector<long double> scaled_coeffs(m_vars);
        for (std::size_t i = 0; i < m_vars; ++i) scaled_coeffs[i] = solution.at(i, 0);
        const std::vector<long double> linear_map = {-center / scale, 1.0L / scale};
        std::vector<long double> coefficients = polynomial_compose(scaled_coeffs, linear_map);
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
std::string polynomial_to_string(const std::vector<long double>& coefficients,
                                 const std::string& variable_name) {
    std::vector<long double> normalized = coefficients;
    trim_trailing_zeros(&normalized);
    if (normalized.size() == 1 && mymath::is_near_zero(normalized[0], kPolynomialEps)) return "0";
    std::string result;
    bool first = true;
    for (std::size_t index = normalized.size(); index > 0; --index) {
        const std::size_t degree = index - 1;
        const long double coefficient = normalized[degree];
        if (mymath::is_near_zero(coefficient, kPolynomialEps)) continue;
        const bool negative = coefficient < 0.0L;
        const long double abs_value = negative ? -coefficient : coefficient;
        std::string term;
        if (degree == 0) term = format_coefficient(abs_value);
        else {
            if (!mymath::is_near_zero(abs_value - 1.0L, kPolynomialEps)) term += format_coefficient(abs_value) + " * ";
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

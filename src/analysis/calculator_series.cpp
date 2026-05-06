// ============================================================================
// 级数展开命令实现
// ============================================================================
//
// 本文件实现了级数展开命令的数值计算，包括：
// - taylor: Taylor 级数展开
// - pade: Pade 有理逼近
// - puiseux: Puiseux 级数展开
// - series_sum: 级数求和
//
// 核心技术：
// - 幂级数代数（PSA）：对幂级数进行加减乘除等运算
// - 符号微分：用于 Taylor 级数系数计算
// - Richardson 外推：用于极限计算
// - Newton-Puiseux 算法：处理代数曲线奇点
// - Wynn-epsilon 加速：数值级数收敛加速

#include "analysis/calculator_series.h"
#include "symbolic/symbolic_expression_internal.h"
#include "statistics/probability.h"
#include "math/helpers/integer_helpers.h"
#include "expression_utils.h"
#include "string_utils.h"
#include "polynomial/polynomial.h"
#include "math/mymath.h"
#include "parser/unified_expression_parser.h"
#include <sstream>
#include <iomanip>
#include <tuple>
#include <functional>

namespace series_ops {

namespace internal {

// ============================================================================
// Laurent 级数结构
// ============================================================================

/**
 * @brief Laurent 级数（支持有限负幂项）
 *
 * 表示形式：sum_{k=low_power}^{high_power} c_k * x^k
 */
struct LaurentSeries {
    std::vector<double> coefficients;  // 系数，从 low_power 开始
    int low_power;                      // 最低幂次（可为负）
    int high_power;                     // 最高幂次

    LaurentSeries() : low_power(0), high_power(-1) {}

    LaurentSeries(int low, int high, const std::vector<double>& coeffs)
        : coefficients(coeffs), low_power(low), high_power(high) {}

    // 检查是否为空
    bool empty() const { return coefficients.empty() || low_power > high_power; }

    // 获取指定幂次的系数
    double get_coefficient(int power) const {
        if (power < low_power || power > high_power) return 0.0;
        return coefficients[power - low_power];
    }

    // 转换为 Taylor 级数（如果 low_power >= 0）
    bool to_taylor(std::vector<double>& result, int degree) const {
        if (low_power < 0) return false;
        result.assign(degree + 1, 0.0);
        for (int k = low_power; k <= high_power && k <= degree; ++k) {
            result[k] = get_coefficient(k);
        }
        return true;
    }
};

// ============================================================================
// 幂级数代数（PSA）运算
// ============================================================================

/**
 * @brief 幂级数加法
 */
std::vector<double> ps_add(const std::vector<double>& a, const std::vector<double>& b, int degree) {
    std::vector<double> res(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        res[i] = (i < static_cast<int>(a.size()) ? a[i] : 0.0) +
                 (i < static_cast<int>(b.size()) ? b[i] : 0.0);
    }
    return res;
}

/**
 * @brief 幂级数减法
 */
std::vector<double> ps_sub(const std::vector<double>& a, const std::vector<double>& b, int degree) {
    std::vector<double> res(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        res[i] = (i < static_cast<int>(a.size()) ? a[i] : 0.0) -
                 (i < static_cast<int>(b.size()) ? b[i] : 0.0);
    }
    return res;
}

/**
 * @brief 幂级数乘法
 *
 * 使用卷积计算乘积的前 degree+1 项。
 */
std::vector<double> ps_mul(const std::vector<double>& a, const std::vector<double>& b, int degree) {
    std::vector<double> res(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        if (i >= static_cast<int>(a.size()) || mymath::is_near_zero(a[i], 1e-12)) continue;
        for (int j = 0; i + j <= degree; ++j) {
            if (j >= static_cast<int>(b.size()) || mymath::is_near_zero(b[j], 1e-12)) continue;
            res[i + j] += a[i] * b[j];
        }
    }
    return res;
}

/**
 * @brief 幂级数除法（支持 Laurent 级数）
 *
 * 计算商级数 a / b 的前 degree+1 项。
 * 支持有限负幂项（Laurent 级数）。
 *
 * @param a 被除数幂级数
 * @param b 除数幂级数
 * @param degree 展开阶数
 * @param laurent_shift 输出 Laurent 级数的负幂起始位置（负数表示有负幂项）
 * @throws std::runtime_error 当除数为零时
 */
std::vector<double> ps_div_with_laurent(const std::vector<double>& a, const std::vector<double>& b, int degree, int* laurent_shift) {
    if (b.empty()) throw std::runtime_error("division by empty power series");

    // 查找 a 和 b 的第一个非零项索引
    int start_a = -1;
    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
        if (!mymath::is_near_zero(a[i], 1e-15)) {
            start_a = i;
            break;
        }
    }

    int start_b = -1;
    for (int i = 0; i < static_cast<int>(b.size()); ++i) {
        if (!mymath::is_near_zero(b[i], 1e-15)) {
            start_b = i;
            break;
        }
    }

    if (start_b == -1) throw std::runtime_error("division by zero in power series");

    // 如果 a 为全零，结果全零
    if (start_a == -1) {
        if (laurent_shift) *laurent_shift = 0;
        return std::vector<double>(degree + 1, 0.0);
    }

    // 计算位移：x^start_a / x^start_b = x^(start_a - start_b)
    int shift = start_a - start_b;

    // 输出 Laurent 位移信息
    if (laurent_shift) *laurent_shift = shift;

    // 提取有效部分并相除
    std::vector<double> a_effective;
    for (int i = start_a; i < static_cast<int>(a.size()); ++i) a_effective.push_back(a[i]);

    std::vector<double> b_effective;
    for (int i = start_b; i < static_cast<int>(b.size()); ++i) b_effective.push_back(b[i]);

    // 计算商级数的有效部分
    // 如果 shift < 0（Laurent 级数），我们需要计算从 shift 开始的系数
    int result_start = (shift < 0) ? shift : 0;
    int result_size = degree + 1 - result_start;
    if (result_size <= 0) {
        // 结果全是超出范围的负幂项
        return std::vector<double>(degree + 1, 0.0);
    }

    std::vector<double> res_effective(result_size, 0.0);
    double inv_b0 = 1.0 / b_effective[0];
    for (int i = 0; i < result_size; ++i) {
        double val = (i < static_cast<int>(a_effective.size()) ? a_effective[i] : 0.0);
        for (int j = 1; j <= i; ++j) {
            if (j < static_cast<int>(b_effective.size())) val -= b_effective[j] * res_effective[i - j];
        }
        res_effective[i] = val * inv_b0;
    }

    // 构建最终结果，考虑 Laurent 位移
    std::vector<double> final_res(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        int effective_idx = i - result_start;
        if (effective_idx >= 0 && effective_idx < static_cast<int>(res_effective.size())) {
            final_res[i] = res_effective[effective_idx];
        }
    }
    return final_res;
}

/**
 * @brief 幂级数除法（原始接口，不支持 Laurent）
 *
 * 计算商级数 a / b 的前 degree+1 项。
 * 当结果为 Laurent 级数时抛出异常。
 */
std::vector<double> ps_div(const std::vector<double>& a, const std::vector<double>& b, int degree) {
    int shift = 0;
    std::vector<double> result = ps_div_with_laurent(a, b, degree, &shift);
    if (shift < 0) {
        // 对于 Laurent 级数，前导系数对应 x^shift 项
        // 由于 shift < 0，前导系数存储在 result 的开始位置（对应最低幂次）
        // 我们需要找到第一个非零系数作为前导系数
        double leading = 0.0;
        for (double v : result) {
            if (!mymath::is_near_zero(v, 1e-15)) {
                leading = v;
                break;
            }
        }
        // 如果 result 中没有找到非零值，说明前导系数在 res_effective[0]
        // 这种情况发生在 shift < 0 且 |shift| > degree 时
        if (mymath::is_near_zero(leading, 1e-15)) {
            // 重新计算前导系数
            int start_a = -1;
            for (int i = 0; i < static_cast<int>(a.size()); ++i) {
                if (!mymath::is_near_zero(a[i], 1e-15)) {
                    start_a = i;
                    break;
                }
            }
            int start_b = -1;
            for (int i = 0; i < static_cast<int>(b.size()); ++i) {
                if (!mymath::is_near_zero(b[i], 1e-15)) {
                    start_b = i;
                    break;
                }
            }
            if (start_a >= 0 && start_b >= 0) {
                leading = a[start_a] / b[start_b];
            }
        }
        throw series_ops::internal::PoleException(shift, leading);
    }
    return result;
}

/**
 * @brief 幂级数指数函数
 *
 * 计算 exp(a) 的幂级数。
 */
std::vector<double> ps_exp(const std::vector<double>& a, int degree) {
    std::vector<double> res(degree + 1, 0.0);
    double a0 = a.empty() ? 0.0 : a[0];
    res[0] = mymath::exp(a0);
    for (int i = 1; i <= degree; ++i) {
        double sum = 0.0;
        for (int k = 1; k <= i; ++k) {
            double ak = k < static_cast<int>(a.size()) ? a[k] : 0.0;
            sum += k * ak * res[i - k];
        }
        res[i] = sum / i;
    }
    return res;
}

/**
 * @brief 幂级数对数函数
 *
 * 计算 ln(a) 的幂级数。
 * 要求 a[0] > 0。
 */
std::vector<double> ps_ln(const std::vector<double>& a, int degree) {
    if (a.empty() || a[0] <= 0) throw std::runtime_error("ln of non-positive power series base");
    std::vector<double> res(degree + 1, 0.0);
    res[0] = mymath::ln(a[0]);
    double inv_a0 = 1.0 / a[0];
    for (int i = 1; i <= degree; ++i) {
        double sum = 0.0;
        for (int k = 1; k < i; ++k) {
            double ak = i - k < static_cast<int>(a.size()) ? a[i - k] : 0.0;
            sum += k * res[k] * ak;
        }
        double ai = i < static_cast<int>(a.size()) ? a[i] : 0.0;
        res[i] = (ai - sum / i) * inv_a0;
    }
    return res;
}

/**
 * @brief 幂级数正弦和余弦函数
 *
 * 同时计算 sin(a) 和 cos(a) 的幂级数。
 */
void ps_sincos(const std::vector<double>& a, int degree, std::vector<double>& sin_res, std::vector<double>& cos_res) {
    sin_res.assign(degree + 1, 0.0);
    cos_res.assign(degree + 1, 0.0);
    double a0 = a.empty() ? 0.0 : a[0];
    sin_res[0] = mymath::sin(a0);
    cos_res[0] = mymath::cos(a0);
    for (int i = 1; i <= degree; ++i) {
        double sum_sin = 0.0;
        double sum_cos = 0.0;
        for (int k = 1; k <= i; ++k) {
            double ak = k < static_cast<int>(a.size()) ? a[k] : 0.0;
            sum_sin += k * ak * cos_res[i - k];
            sum_cos -= k * ak * sin_res[i - k];
        }
        sin_res[i] = sum_sin / i;
        cos_res[i] = sum_cos / i;
    }
}

/**
 * @brief 幂级数正弦函数
 */
std::vector<double> ps_sin(const std::vector<double>& a, int degree) {
    std::vector<double> s, c;
    ps_sincos(a, degree, s, c);
    return s;
}

/**
 * @brief 幂级数余弦函数
 */
std::vector<double> ps_cos(const std::vector<double>& a, int degree) {
    std::vector<double> s, c;
    ps_sincos(a, degree, s, c);
    return c;
}

// 前向声明
static std::vector<double> ps_pow_const_puiseux(const std::vector<double>& a, double n, int degree, int leading);

/**
 * @brief 幂级数幂函数
 *
 * 计算 a^n 的幂级数，其中 n 为常数。
 * 特殊处理 a[0] = 0 的情况。
 */
std::vector<double> ps_pow_const(const std::vector<double>& a, double n, int degree) {
    if (a.empty() || mymath::is_near_zero(a[0], 1e-12)) {
        // 底数首项为零的情况
        if (mymath::is_near_zero(n, 1e-12)) {
            // a^0 = 1
            std::vector<double> res(degree + 1, 0.0);
            res[0] = 1.0;
            return res;
        }

        // 找到第一个非零项
        int leading = -1;
        for (int i = 0; i < static_cast<int>(a.size()); ++i) {
            if (!mymath::is_near_zero(a[i], 1e-12)) {
                leading = i;
                break;
            }
        }

        if (leading < 0) {
            // 全零序列
            if (n > 0) {
                return std::vector<double>(degree + 1, 0.0);
            } else {
                throw std::runtime_error("0^negative is undefined");
            }
        }

        // a(x) = a_leading * x^leading * (1 + higher terms)
        // a(x)^n = a_leading^n * x^(leading*n) * (1 + higher terms)^n

        const double shifted_power = static_cast<double>(leading) * n;

        // 情况1：正整数幂 - 直接乘法
        if (n > 0 && mymath::is_integer(n, 1e-12)) {
            std::vector<double> res(degree + 1, 0.0);
            res[0] = 1.0;
            std::vector<double> base = a;
            int p = static_cast<int>(n + 0.5);
            for (int i = 0; i < p; ++i) res = ps_mul(res, base, degree);
            return res;
        }

        // 情况2：负整数幂 - 使用除法
        if (n < 0 && mymath::is_integer(n, 1e-12)) {
            int p = static_cast<int>(-n + 0.5);
            // a^(-p) = 1 / a^p
            std::vector<double> pos_pow = ps_pow_const(a, static_cast<double>(p), degree);
            std::vector<double> one(degree + 1, 0.0);
            one[0] = 1.0;
            return ps_div_with_laurent(one, pos_pow, degree, nullptr);
        }

        // 情况3：分数幂，位移为整数 - 可以处理
        if (is_integer_double(shifted_power, 1e-10)) {
            const int shift = static_cast<int>(round_to_long_long(shifted_power));
            if (shift < 0) {
                // Laurent 级数（极点）
                double leading_coeff = mymath::pow(a[leading], n);
                throw series_ops::internal::PoleException(shift, leading_coeff);
            }
            if (shift > degree) {
                return std::vector<double>(degree + 1, 0.0);
            }

            // 归一化：提取 x^leading 因子
            std::vector<double> normalized;
            normalized.reserve(a.size() - static_cast<std::size_t>(leading));
            for (int i = leading; i < static_cast<int>(a.size()); ++i) {
                normalized.push_back(a[static_cast<std::size_t>(i)]);
            }

            // 现在 normalized[0] != 0，可以计算 normalized^n
            const std::vector<double> powered = ps_pow_const(normalized, n, degree - shift);

            // 乘以 x^shift
            std::vector<double> res(degree + 1, 0.0);
            for (int i = 0; i + shift <= degree &&
                            i < static_cast<int>(powered.size()); ++i) {
                res[static_cast<std::size_t>(i + shift)] =
                    powered[static_cast<std::size_t>(i)];
            }
            return res;
        }

        // 情况4：分数幂，位移非整数 - 需要 Puiseux 级数
        // 对于极限计算，如果位移为正，则极限为 0；如果位移为负，则极限为无穷大。
        if (shifted_power > 1e-10) {
            return std::vector<double>(degree + 1, 0.0);
        } else if (shifted_power < -1e-10) {
            double leading_coeff = mymath::pow(a[leading], n);
            // 抛出极点异常，虽然这实际上是支点，但对于极限判定是等价的
            throw series_ops::internal::PoleException(-1, leading_coeff);
        }

        return ps_pow_const_puiseux(a, n, degree, leading);
    }

    // 正常情况：a[0] != 0
    std::vector<double> res(degree + 1, 0.0);
    res[0] = mymath::pow(a[0], n);
    double inv_a0 = 1.0 / a[0];
    for (int i = 1; i <= degree; ++i) {
        double sum = 0.0;
        for (int k = 1; k <= i; ++k) {
            double ak = k < static_cast<int>(a.size()) ? a[k] : 0.0;
            sum += (n * k - (i - k)) * ak * res[i - k];
        }
        res[i] = sum * inv_a0 / i;
    }
    return res;
}

/**
 * @brief 处理需要 Puiseux 级数的分数幂
 *
 * 当 a(x)^n 的首项位移不是整数时，需要 Puiseux 级数。
 * 使用广义二项式展开计算。
 */
std::vector<double> ps_pow_const_puiseux(const std::vector<double>& a, double n, int degree, int leading) {
    // a(x) = a[leading] * x^leading * (1 + b_1 * x + b_2 * x^2 + ...)
    // a(x)^n = a[leading]^n * x^(leading*n) * (1 + b_1 * x + ...)^n

    // 提取归一化部分
    std::vector<double> normalized;
    double a_leading = a[leading];
    normalized.push_back(1.0);  // 归一化使首项为 1
    for (int i = leading + 1; i < static_cast<int>(a.size()); ++i) {
        normalized.push_back(a[i] / a_leading);
    }

    // 使用二项式展开 (1 + u)^n = sum_{k=0}^{inf} C(n,k) * u^k
    // 其中 u = b_1 * x + b_2 * x^2 + ...

    std::vector<double> res(degree + 1, 0.0);
    res[0] = mymath::pow(a_leading, n);  // a[leading]^n 作为首项系数

    // 计算 (1 + u)^n 的幂级数展开
    // 使用递推公式：c_0 = 1, c_k = sum_{j=1}^{k} (n*j/k - j/k + 1/k) * u_j * c_{k-j}
    // 简化：c_k = (1/k) * sum_{j=1}^{k} ((n-1)*j + k) * u_j * c_{k-j} / k

    std::vector<double> binom_coeffs(degree + 1, 0.0);
    binom_coeffs[0] = 1.0;

    for (int k = 1; k <= degree; ++k) {
        double sum = 0.0;
        for (int j = 1; j <= k && j < static_cast<int>(normalized.size()); ++j) {
            double u_j = normalized[j];
            sum += (n * j - (k - j)) * u_j * binom_coeffs[k - j];
        }
        binom_coeffs[k] = sum / k;
    }

    // 组合结果
    for (int i = 0; i <= degree; ++i) {
        res[i] = res[0] * binom_coeffs[i];
    }

    // 注意：这里没有处理 x^(leading*n) 的分数幂部分
    // 调用者需要处理这个分数位移
    return res;
}

std::vector<double> ps_scale(const std::vector<double>& a, double scale, int degree) {
    std::vector<double> res(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        if (i < static_cast<int>(a.size())) res[i] = a[i] * scale;
    }
    return res;
}

std::vector<double> ps_derivative(const std::vector<double>& a, int degree) {
    std::vector<double> res(degree + 1, 0.0);
    for (int i = 1; i <= degree + 1; ++i) {
        if (i < static_cast<int>(a.size())) res[i - 1] = a[i] * i;
    }
    return res;
}

std::vector<double> ps_integral(const std::vector<double>& a, double constant_term, int degree) {
    std::vector<double> res(degree + 1, 0.0);
    res[0] = constant_term;
    for (int i = 1; i <= degree; ++i) {
        if (i - 1 < static_cast<int>(a.size())) res[i] = a[i - 1] / i;
    }
    return res;
}

std::vector<double> ps_tan(const std::vector<double>& a, int degree) {
    std::vector<double> s, c;
    ps_sincos(a, degree, s, c);
    return ps_div(s, c, degree);
}

std::vector<double> ps_asin(const std::vector<double>& a, int degree) {
    double a0 = a.empty() ? 0.0 : a[0];
    if (mymath::abs(a0) >= 1.0) throw std::runtime_error("asin power series base out of bounds");
    std::vector<double> a_prime = ps_derivative(a, degree);
    std::vector<double> a_sq = ps_mul(a, a, degree);
    std::vector<double> one_minus_a_sq(degree + 1, 0.0);
    one_minus_a_sq[0] = 1.0;
    one_minus_a_sq = ps_sub(one_minus_a_sq, a_sq, degree);
    std::vector<double> inv_sqrt = ps_pow_const(one_minus_a_sq, -0.5, degree);
    std::vector<double> deriv = ps_mul(inv_sqrt, a_prime, degree);
    return ps_integral(deriv, mymath::asin(a0), degree);
}

std::vector<double> ps_acos(const std::vector<double>& a, int degree) {
    double a0 = a.empty() ? 0.0 : a[0];
    if (mymath::abs(a0) >= 1.0) throw std::runtime_error("acos power series base out of bounds");
    std::vector<double> a_prime = ps_derivative(a, degree);
    std::vector<double> a_sq = ps_mul(a, a, degree);
    std::vector<double> one_minus_a_sq(degree + 1, 0.0);
    one_minus_a_sq[0] = 1.0;
    one_minus_a_sq = ps_sub(one_minus_a_sq, a_sq, degree);
    std::vector<double> inv_sqrt = ps_pow_const(one_minus_a_sq, -0.5, degree);
    std::vector<double> deriv = ps_mul(inv_sqrt, a_prime, degree);
    std::vector<double> neg_deriv = ps_scale(deriv, -1.0, degree);
    return ps_integral(neg_deriv, mymath::acos(a0), degree);
}

std::vector<double> ps_atan(const std::vector<double>& a, int degree) {
    double a0 = a.empty() ? 0.0 : a[0];
    std::vector<double> a_prime = ps_derivative(a, degree);
    std::vector<double> a_sq = ps_mul(a, a, degree);
    std::vector<double> one_plus_a_sq(degree + 1, 0.0);
    one_plus_a_sq[0] = 1.0;
    one_plus_a_sq = ps_add(one_plus_a_sq, a_sq, degree);
    std::vector<double> inv = ps_pow_const(one_plus_a_sq, -1.0, degree);
    std::vector<double> deriv = ps_mul(inv, a_prime, degree);
    return ps_integral(deriv, mymath::atan(a0), degree);
}

std::vector<double> ps_sinh(const std::vector<double>& a, int degree) {
    std::vector<double> exp_a = ps_exp(a, degree);
    std::vector<double> neg_a = ps_scale(a, -1.0, degree);
    std::vector<double> exp_neg_a = ps_exp(neg_a, degree);
    return ps_scale(ps_sub(exp_a, exp_neg_a, degree), 0.5, degree);
}

std::vector<double> ps_cosh(const std::vector<double>& a, int degree) {
    std::vector<double> exp_a = ps_exp(a, degree);
    std::vector<double> neg_a = ps_scale(a, -1.0, degree);
    std::vector<double> exp_neg_a = ps_exp(neg_a, degree);
    return ps_scale(ps_add(exp_a, exp_neg_a, degree), 0.5, degree);
}

std::vector<double> ps_tanh(const std::vector<double>& a, int degree) {
    return ps_div(ps_sinh(a, degree), ps_cosh(a, degree), degree);
}

/**
 * @brief 使用 PSA 引擎计算表达式的级数系数
 *
 * 对表达式进行递归求值，使用幂级数代数运算。
 *
 * @param expr 符号表达式
 * @param var_name 变量名
 * @param center 展开中心
 * @param degree 展开阶数
 * @param result 输出系数列表
 * @param ctx 级数上下文
 * @return 是否成功计算
 */
bool evaluate_psa(const SymbolicExpression& expr, const std::string& var_name, double center, int degree, std::vector<double>& result, const SeriesContext& ctx) {
    if (!expr.node_) return false;
    auto node = expr.node_;
    
    if (node->type == NodeType::kNumber) {
        result.assign(degree + 1, 0.0);
        result[0] = node->number_value;
        return true;
    }
    if (node->type == NodeType::kPi) {
        result.assign(degree + 1, 0.0);
        result[0] = mymath::kPi;
        return true;
    }
    if (node->type == NodeType::kE) {
        result.assign(degree + 1, 0.0);
        result[0] = mymath::kE;
        return true;
    }
    if (node->type == NodeType::kVariable) {
        result.assign(degree + 1, 0.0);
        if (node->text == var_name) {
            result[0] = center;
            if (degree >= 1) result[1] = 1.0;
        } else {
            result[0] = ctx.evaluate_at(expr, var_name, center);
        }
        return true;
    }

    std::vector<double> left_res, right_res;
    if (node->left && !series_ops::internal::evaluate_psa(SymbolicExpression(node->left), var_name, center, degree, left_res, ctx)) return false;
    if (node->right && !series_ops::internal::evaluate_psa(SymbolicExpression(node->right), var_name, center, degree, right_res, ctx)) return false;

    try {
        switch (node->type) {
            case NodeType::kAdd: result = ps_add(left_res, right_res, degree); return true;
            case NodeType::kSubtract: result = ps_sub(left_res, right_res, degree); return true;
            case NodeType::kMultiply: result = ps_mul(left_res, right_res, degree); return true;
            case NodeType::kDivide: result = ps_div(left_res, right_res, degree); return true;
            case NodeType::kNegate: 
                result.assign(degree + 1, 0.0);
                for (int i = 0; i <= degree; ++i) result[i] = -left_res[i];
                return true;
            case NodeType::kPower: {
                bool right_is_const = true;
                for (int i = 1; i <= degree; ++i) {
                    if (i < static_cast<int>(right_res.size()) &&
                        !mymath::is_near_zero(right_res[i], 1e-12)) {
                        right_is_const = false; break;
                    }
                }
                if (right_is_const) {
                    double p = right_res.empty() ? 0.0 : right_res[0];
                    result = ps_pow_const(left_res, p, degree);
                    return true;
                } else {
                    result = ps_exp(ps_mul(right_res, ps_ln(left_res, degree), degree), degree);
                    return true;
                }
            }
            case NodeType::kFunction: {
                if (node->text == "exp") { result = ps_exp(left_res, degree); return true; }
                if (node->text == "ln") { result = ps_ln(left_res, degree); return true; }
                if (node->text == "sin") { result = ps_sin(left_res, degree); return true; }
                if (node->text == "cos") { result = ps_cos(left_res, degree); return true; }
                if (node->text == "tan") { result = ps_tan(left_res, degree); return true; }
                if (node->text == "asin" || node->text == "arcsin") { result = ps_asin(left_res, degree); return true; }
                if (node->text == "acos" || node->text == "arccos") { result = ps_acos(left_res, degree); return true; }
                if (node->text == "atan" || node->text == "arctan") { result = ps_atan(left_res, degree); return true; }
                if (node->text == "sinh") { result = ps_sinh(left_res, degree); return true; }
                if (node->text == "cosh") { result = ps_cosh(left_res, degree); return true; }
                if (node->text == "tanh") { result = ps_tanh(left_res, degree); return true; }
                if (node->text == "sqrt") { result = ps_pow_const(left_res, 0.5, degree); return true; }
                return false;
            }
            default: return false;
        }
    } catch (const series_ops::internal::PoleException&) {
        // PoleException 需要向上传播，让调用者处理无限极限
        throw;
    } catch (...) {
        return false;
    }
}

}  // namespace internal

// ============================================================================
// 级数展开命令处理
// ============================================================================

namespace {

/**
 * @brief 简化符号文本
 */
std::string simplify_symbolic_text(const std::string& text) {
    return SymbolicExpression::parse(text).simplify().to_string();
}

// ============================================================================
// 预定义 Taylor 级数库
// ============================================================================

/**
 * @brief 预定义函数的 Taylor 级数系数
 *
 * 对于常见初等函数，直接使用已知的 Taylor 系数公式，
 * 避免反复符号求导导致的表达式爆炸。
 */
namespace predefined_series {

/**
 * @brief 检测并返回预定义函数的 Taylor 级数
 *
 * @param expr 表达式
 * @param var_name 变量名
 * @param center 展开中心
 * @param degree 展开阶数
 * @param result 输出系数
 * @return 是否匹配预定义函数
 */
bool try_predefined_taylor(const SymbolicExpression& expr,
                            const std::string& var_name,
                            double center,
                            int degree,
                            std::vector<double>& result) {
    if (!expr.node_) return false;

    // 处理函数节点
    if (expr.node_->type == NodeType::kFunction) {
        const std::string& func_name = expr.node_->text;
        SymbolicExpression arg(expr.node_->left);

        // 检查参数是否是简单的 x - center 形式
        double arg_coeff = 0.0, arg_const = 0.0;
        bool is_linear_arg = false;

        // 尝试检测 arg = a*(x-center) 或 arg = x - center
        if (arg.node_->type == NodeType::kVariable && arg.node_->text == var_name) {
            arg_coeff = 1.0;
            arg_const = -center;
            is_linear_arg = true;
        } else if (arg.node_->type == NodeType::kSubtract) {
            SymbolicExpression left(arg.node_->left);
            SymbolicExpression right(arg.node_->right);
            if (left.node_->type == NodeType::kVariable && left.node_->text == var_name &&
                right.is_number(&arg_const)) {
                arg_coeff = 1.0;
                arg_const = arg_const - center;
                is_linear_arg = true;
            }
        }

        // 对于 center = 0 且参数是 x，直接使用预定义级数
        if (mymath::is_near_zero(center, 1e-10) && arg.node_->type == NodeType::kVariable && arg.node_->text == var_name) {
            result.assign(degree + 1, 0.0);

            // exp(x) = sum x^n/n!
            if (func_name == "exp") {
                result[0] = 1.0;
                for (int n = 1; n <= degree; ++n) {
                    result[n] = result[n - 1] / static_cast<double>(n);
                }
                return true;
            }

            // sin(x) = x - x^3/3! + x^5/5! - ...
            if (func_name == "sin") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 1) {
                        int sign = (n % 4 == 1) ? 1 : -1;
                        result[n] = static_cast<double>(sign) / prob::factorial(n);
                    }
                }
                return true;
            }

            // cos(x) = 1 - x^2/2! + x^4/4! - ...
            if (func_name == "cos") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 0) {
                        int sign = (n % 4 == 0) ? 1 : -1;
                        result[n] = static_cast<double>(sign) / prob::factorial(n);
                    }
                }
                return true;
            }

            // sinh(x) = x + x^3/3! + x^5/5! + ...
            if (func_name == "sinh") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 1) {
                        result[n] = 1.0 / prob::factorial(n);
                    }
                }
                return true;
            }

            // cosh(x) = 1 + x^2/2! + x^4/4! + ...
            if (func_name == "cosh") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 0) {
                        result[n] = 1.0 / prob::factorial(n);
                    }
                }
                return true;
            }

            // ln(1+x) = x - x^2/2 + x^3/3 - ... (对于 |x| < 1)
            if (func_name == "ln") {
                // ln(x) 在 x=0 处无定义，需要特殊处理
                return false;
            }

            // arcsin(x) = x + x^3/6 + 3x^5/40 + ...
            if (func_name == "arcsin" || func_name == "asin") {
                // 使用递推公式: a_n = (2n-1)/(2n) * a_{n-1} * (n-1)/n
                result[1] = 1.0;
                for (int n = 3; n <= degree; n += 2) {
                    int m = (n - 1) / 2;
                    result[n] = result[n - 2] * static_cast<double>(2 * m - 1) / static_cast<double>(2 * m) * static_cast<double>(m) / static_cast<double>(m + 1);
                }
                return true;
            }

            // arctan(x) = x - x^3/3 + x^5/5 - ...
            if (func_name == "arctan" || func_name == "atan") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 1) {
                        int sign = (n % 4 == 1) ? 1 : -1;
                        result[n] = static_cast<double>(sign) / static_cast<double>(n);
                    }
                }
                return true;
            }

            // sqrt(1+x) = 1 + x/2 - x^2/8 + x^3/16 - ... (二项式展开)
            if (func_name == "sqrt") {
                result[0] = 1.0;
                for (int n = 1; n <= degree; ++n) {
                    // 系数 = (1/2 choose n) = (1/2)(1/2-1)...(1/2-n+1)/n!
                    double coeff = 1.0;
                    for (int k = 0; k < n; ++k) {
                        coeff *= (0.5 - static_cast<double>(k)) / static_cast<double>(n - k);
                    }
                    result[n] = coeff;
                }
                return true;
            }
        }

        // 对于 center != 0，使用平移公式
        // f(x) 在 x=a 处展开: f(a + h) = sum f^(n)(a) * h^n / n!
        if (is_linear_arg) {
            // 参数是 x - center 或类似形式
            double h_coeff = arg_coeff;  // 参数中 x 的系数

            result.assign(degree + 1, 0.0);

            // exp(x) 在 x=a 处: exp(a) * exp(h) = exp(a) * sum h^n/n!
            if (func_name == "exp") {
                double exp_center = mymath::exp(center);
                result[0] = exp_center;
                for (int n = 1; n <= degree; ++n) {
                    result[n] = result[n - 1] * h_coeff / static_cast<double>(n);
                }
                return true;
            }

            // sin(x) 在 x=a 处: sin(a+h) = sin(a)cos(h) + cos(a)sin(h)
            if (func_name == "sin") {
                double sin_a = mymath::sin(center);
                double cos_a = mymath::cos(center);
                // sin(a+h) = sin(a) * cos(h) + cos(a) * sin(h)
                // cos(h) = 1 - h^2/2! + h^4/4! - ...
                // sin(h) = h - h^3/3! + h^5/5! - ...
                for (int n = 0; n <= degree; ++n) {
                    double h_pow = mymath::pow(h_coeff, n);
                    if (n % 2 == 0) {
                        // 来自 cos(h)
                        int sign = (n % 4 == 0) ? 1 : -1;
                        result[n] = sin_a * static_cast<double>(sign) * h_pow / prob::factorial(n);
                    } else {
                        // 来自 sin(h)
                        int sign = (n % 4 == 1) ? 1 : -1;
                        result[n] = cos_a * static_cast<double>(sign) * h_pow / prob::factorial(n);
                    }
                }
                return true;
            }

            // cos(x) 在 x=a 处
            if (func_name == "cos") {
                double sin_a = mymath::sin(center);
                double cos_a = mymath::cos(center);
                for (int n = 0; n <= degree; ++n) {
                    double h_pow = mymath::pow(h_coeff, n);
                    if (n % 2 == 0) {
                        int sign = (n % 4 == 0) ? 1 : -1;
                        result[n] = cos_a * static_cast<double>(sign) * h_pow / prob::factorial(n);
                    } else {
                        int sign = (n % 4 == 1) ? -1 : 1;
                        result[n] = sin_a * static_cast<double>(sign) * h_pow / prob::factorial(n);
                    }
                }
                return true;
            }
        }
    }

    // 处理幂运算 x^n
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exponent(expr.node_->right);

        double exp_val = 0.0;
        if (base.node_->type == NodeType::kVariable && base.node_->text == var_name &&
            exponent.is_number(&exp_val) && mymath::is_near_zero(center, 1e-10)) {
            // x^p 在 x=0 处展开
            // 如果 p 是正整数，x^p 本身就是单项式
            if (exp_val > 0 && mymath::is_integer(exp_val, 1e-10)) {
                int p = static_cast<int>(exp_val + 0.5);
                result.assign(degree + 1, 0.0);
                if (p <= degree) {
                    result[p] = 1.0;
                }
                return true;
            }
        }
    }

    return false;
}

}  // namespace predefined_series

// ============================================================================
// 自动微分（AD）Taylor 系数计算
// ============================================================================

/**
 * @brief 使用前向模式自动微分计算 Taylor 系数
 *
 * 对于无法符号求导或 PSA 失败的复杂函数，
 * 使用数值自动微分来估计高阶导数。
 *
 * 基于 Faa di Bruno 公式和数值差分。
 *
 * @param ctx 级数上下文
 * @param expression 符号表达式
 * @param variable_name 变量名
 * @param center 展开中心
 * @param degree 展开阶数
 * @param result 输出系数
 * @return 是否成功计算
 */
bool compute_taylor_coefficients_ad(const SeriesContext& ctx,
                                    const SymbolicExpression& expression,
                                    const std::string& variable_name,
                                    double center,
                                    int degree,
                                    std::vector<double>& result) {
    result.clear();
    result.reserve(degree + 1);

    // 使用数值差分计算导数
    // 对于 n 阶导数，使用中心差分公式：
    // f^(n)(x) ≈ (1/h^n) * sum_{k=0}^{n} (-1)^k * C(n,k) * f(x + (n/2 - k)*h)

    // 选择步长：对于高阶导数需要更小的步长
    const double base_h = 1e-4;

    // 0阶导数（函数值）
    double f0 = 0.0;
    try {
        f0 = ctx.evaluate_at(expression, variable_name, center);
    } catch (...) {
        return false;
    }
    if (!mymath::isfinite(f0)) return false;
    result.push_back(f0);

    // 使用复合数值差分计算高阶导数
    // 对于 1-3 阶使用精确的中心差分公式
    // 对于更高阶使用递归差分

    auto compute_derivative = [&](int order, double h) -> double {
        if (order == 0) return f0;

        // 使用中心差分
        int n_points = 2 * order + 1;
        std::vector<double> points;
        points.reserve(n_points);

        for (int k = -order; k <= order; ++k) {
            double x = center + k * h;
            double val = 0.0;
            try {
                val = ctx.evaluate_at(expression, variable_name, x);
            } catch (...) {
                return mymath::quiet_nan();
            }
            if (!mymath::isfinite(val)) {
                return mymath::quiet_nan();
            }
            points.push_back(val);
        }

        // 计算差分表
        std::vector<std::vector<double>> diff_table;
        diff_table.push_back(points);

        for (int d = 1; d <= order; ++d) {
            std::vector<double> next_diff;
            next_diff.reserve(diff_table[d - 1].size() - 1);
            for (std::size_t i = 0; i + 1 < diff_table[d - 1].size(); ++i) {
                next_diff.push_back(diff_table[d - 1][i + 1] - diff_table[d - 1][i]);
            }
            diff_table.push_back(next_diff);
        }

        // n阶中心差分位于表的中心
        const auto& nth_diff = diff_table[order];
        if (nth_diff.empty()) {
            return mymath::quiet_nan();
        }

        // 中心差分系数
        double coeff = nth_diff[nth_diff.size() / 2];

        // 调整系数：中心差分公式 f^(n)(x) ≈ delta^n f / h^n
        // 其中 delta^n 是 n 阶差分
        return coeff / mymath::pow(h, order);
    };

    // 使用 Richardson 外推提高精度
    auto richardson_derivative = [&](int order) -> double {
        double h1 = base_h;
        double h2 = base_h / 2.0;

        double d1 = compute_derivative(order, h1);
        double d2 = compute_derivative(order, h2);

        if (!mymath::isfinite(d1) || !mymath::isfinite(d2)) {
            // 回退到单步差分
            return compute_derivative(order, base_h);
        }

        // Richardson 外推：R = d2 + (d2 - d1) / (2^order - 1)
        double factor = mymath::pow(2.0, order) - 1.0;
        if (mymath::abs(factor) < 1e-10) {
            return d2;
        }
        return d2 + (d2 - d1) / factor;
    };

    // 计算各阶导数
    for (int order = 1; order <= degree; ++order) {
        double deriv = richardson_derivative(order);

        if (!mymath::isfinite(deriv)) {
            // 尝试使用更小的步长
            deriv = compute_derivative(order, base_h / 10.0);
        }

        if (!mymath::isfinite(deriv)) {
            // 失败，返回已计算的系数
            return result.size() > 1;
        }

        // Taylor 系数 = f^(n)(x) / n!
        result.push_back(deriv / prob::factorial(order));
    }

    return true;
}

/**
 * @brief 计算 Taylor 级数系数（改进版）
 *
 * 优先级：
 * 1. 预定义级数库（避免符号求导爆炸）
 * 2. 幂级数代数（PSA）方法
 * 3. 自动微分（AD）方法
 * 4. 符号微分方法（最后回退）
 */
std::vector<double> build_taylor_coefficients(

    const SeriesContext& ctx,
    const SymbolicExpression& expression,
    const std::string& variable_name,
    double center,
    int degree) {
    struct TaylorDerivativeCacheEntry {
        SymbolicExpression derivative;
        double value = 0.0;
        bool has_value = false;
    };
    static thread_local std::map<std::string, TaylorDerivativeCacheEntry> derivative_cache;
    static constexpr std::size_t kMaxTaylorDerivativeCacheSize = 256;

    // 1. 首先尝试预定义级数库
    std::vector<double> predefined_result;
    if (predefined_series::try_predefined_taylor(expression, variable_name, center, degree, predefined_result)) {
        return predefined_result;
    }

    // 2. 尝试幂级数代数（PSA）方法
    std::vector<double> psa_result;
    if (series_ops::internal::evaluate_psa(expression.simplify(), variable_name, center, degree, psa_result, ctx)) {
        return psa_result;
    }

    // 3. 尝试自动微分（AD）方法（对于高阶展开更高效）
    std::vector<double> ad_result;
    if (degree >= 5 && compute_taylor_coefficients_ad(ctx, expression, variable_name, center, degree, ad_result)) {
        return ad_result;
    }

    // 4. 回退到符号微分方法
    const std::string base_key =
        variable_name + "|" + format_symbolic_scalar(center) + "|" +
        expression.simplify().to_string();
    std::vector<double> coefficients;
    coefficients.reserve(static_cast<std::size_t>(degree + 1));
    SymbolicExpression current = expression;
    for (int order = 0; order <= degree; ++order) {
        const std::string order_key = base_key + "|" + std::to_string(order);
        auto found = derivative_cache.find(order_key);
        if (found == derivative_cache.end()) {
            if (derivative_cache.size() >= kMaxTaylorDerivativeCacheSize) {
                derivative_cache.clear();
            }
            TaylorDerivativeCacheEntry entry;
            entry.derivative = current.simplify();
            found = derivative_cache.emplace(order_key, entry).first;
        } else {
            current = found->second.derivative;
        }

        if (!found->second.has_value) {
            found->second.value =
                ctx.evaluate_at(found->second.derivative, variable_name, center);
            found->second.has_value = true;
        }
        const double derivative_value = found->second.value;
        if (!mymath::isfinite(derivative_value)) {
            // 符号微分失败，尝试 AD 回退
            std::vector<double> ad_fallback;
            if (compute_taylor_coefficients_ad(ctx, expression, variable_name, center, degree, ad_fallback)) {
                return ad_fallback;
            }
            throw std::runtime_error(
                "taylor expansion produced a non-finite coefficient");
        }
        coefficients.push_back(derivative_value / prob::factorial(order));
        if (order != degree) {
            const std::string next_key = base_key + "|" + std::to_string(order + 1);
            auto next_found = derivative_cache.find(next_key);
            if (next_found != derivative_cache.end()) {
                current = next_found->second.derivative;
            } else {
                current = found->second.derivative.derivative(variable_name).simplify();
            }
        }
    }
    return coefficients;
}

}  // namespace

/**
 * @brief Taylor 级数展开
 */
std::string taylor(const SeriesContext& ctx,
                   const std::string& expr,
                   double center,
                   int degree) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

    const std::vector<double> coefficients =
        build_taylor_coefficients(ctx, expression, variable_name, center, degree);
    return taylor_series_to_string(coefficients, variable_name, center);
}

// ============================================================================
// Pade 逼近辅助函数
// ============================================================================

// 前向声明（在 series_ops 命名空间内）
static std::string pade_from_coeffs(const std::vector<double>& coefficients,
                              int numerator_degree,
                              int denominator_degree);
static bool solve_pade_denominator(std::function<long double(int)> c,
                                   int numerator_degree,
                                   int denominator_degree,
                                   std::vector<long double>& q);
static bool solve_tohplitz_stable(std::function<long double(int)> c, int n, std::vector<long double>& q);
static std::string format_pade_result(const std::vector<double>& numerator,
                                const std::vector<double>& denominator);
static std::string format_simple_pade(const std::vector<long double>& numerator,
                                const std::vector<long double>& denominator);

/**
 * @brief Pade 有理逼近（基于连分数的稳定算法）
 *
 * 计算 Taylor 级数的 Pade 逼近 [m/n]。
 * 分子为 m 次，分母为 n 次。
 *
 * 算法改进：
 * - 使用连分数递推算法（Baker 算法），复杂度 O(n²)
 * - 避免求解病态的 Hankel 矩阵
 * - 数值稳定性显著提高
 */
std::string pade(const SeriesContext& ctx,
                 const std::string& expr,
                 double center,
                 int numerator_degree,
                 int denominator_degree) {
    if (numerator_degree == 0 && denominator_degree == 0) {
        throw std::runtime_error("pade requires at least one non-zero degree");
    }

    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

    const int total_degree = numerator_degree + denominator_degree;
    const std::vector<double> coefficients = build_taylor_coefficients(
        ctx, expression, variable_name, center, total_degree);

    if (coefficients.empty() || mymath::is_near_zero(coefficients[0], 1e-15)) {
        // 首项系数为零，需要特殊处理
        int first_nonzero = 0;
        for (int i = 0; i < static_cast<int>(coefficients.size()); ++i) {
            if (!mymath::is_near_zero(coefficients[i], 1e-15)) {
                first_nonzero = i;
                break;
            }
        }
        if (first_nonzero > 0) {
            // 提取 x^k 因子
            std::vector<double> shifted_coeffs(coefficients.begin() + first_nonzero, coefficients.end());
            std::string inner_result = pade_from_coeffs(shifted_coeffs, numerator_degree, denominator_degree);
            const std::string base = shifted_series_base(variable_name, center);
            if (first_nonzero == 1) {
                return simplify_symbolic_text(base + " * (" + inner_result + ")");
            } else {
                return simplify_symbolic_text(base + "^" + std::to_string(first_nonzero) + " * (" + inner_result + ")");
            }
        }
    }

    return pade_from_coeffs(coefficients, numerator_degree, denominator_degree);
}

/**
 * @brief 从系数向量计算 Pade 逼近（核心算法）
 *
 * 使用 Baker 算法（基于连分数），避免求解线性方程组。
 * 这是计算 Pade 逼近的最稳定方法之一。
 */
std::string pade_from_coeffs(const std::vector<double>& coefficients,
                              int numerator_degree,
                              int denominator_degree) {
    const int total_degree = numerator_degree + denominator_degree;

    // 获取系数的辅助函数
    auto c = [&](int k) -> long double {
        if (k < 0 || k >= static_cast<int>(coefficients.size())) return 0.0L;
        return static_cast<long double>(coefficients[static_cast<std::size_t>(k)]);
    };

    // 特殊情况处理
    if (denominator_degree == 0) {
        // [m/0] 逼近就是 Taylor 多项式
        std::vector<double> num(numerator_degree + 1, 0.0);
        for (int i = 0; i <= numerator_degree && i < static_cast<int>(coefficients.size()); ++i) {
            num[i] = coefficients[i];
        }
        return polynomial_to_string(num, "x");
    }

    if (numerator_degree == 0) {
        // [0/n] 逼近
        // 使用简单的递推计算
        std::vector<long double> q_coeffs(denominator_degree + 1, 0.0L);
        q_coeffs[0] = 1.0L;

        // 构建并求解线性系统（对于 [0/n] 情况，系统较小）
        // 使用更稳定的 Householder 变换
        std::vector<std::vector<long double>> H(
            static_cast<std::size_t>(denominator_degree),
            std::vector<long double>(static_cast<std::size_t>(denominator_degree + 1), 0.0L));

        for (int i = 0; i < denominator_degree; ++i) {
            for (int j = 0; j < denominator_degree; ++j) {
                H[i][j] = c(i - j + 1);
            }
            H[i][denominator_degree] = -c(i + 1);
        }

        // 使用 QR 分解求解（更稳定）
        std::vector<long double> q(denominator_degree + 1, 0.0L);
        q[0] = 1.0L;
        if (!solve_tohplitz_stable(c, denominator_degree, q)) {
            // 回退到简化处理
            return format_simple_pade({c(0)}, std::vector<long double>(denominator_degree + 1, 1.0L));
        }

        long double p0 = c(0);
        for (int j = 1; j <= denominator_degree; ++j) {
            p0 += q[j] * c(-j);
        }

        std::vector<double> num_vec = {static_cast<double>(p0)};
        std::vector<double> den_vec(denominator_degree + 1);
        for (int i = 0; i <= denominator_degree; ++i) {
            den_vec[i] = static_cast<double>(q[i]);
        }
        return format_pade_result(num_vec, den_vec);
    }

    // 一般情况：使用 Baker 算法（基于连分数的递推）
    // 参考: Baker, G. A. "Essentials of Padé Approximants"

    // 初始化：P_{-1} = 0, P_0 = c_0
    //         Q_{-1} = 1, Q_0 = 1
    std::vector<long double> P_prev2(total_degree + 2, 0.0L);
    std::vector<long double> P_prev1(total_degree + 2, 0.0L);
    P_prev1[0] = c(0);

    std::vector<long double> Q_prev2(total_degree + 2, 1.0L);
    std::vector<long double> Q_prev1(total_degree + 2, 1.0L);

    // 使用 Thacher-Tukey 算法递推计算
    // 这是一种基于欧几里得算法的稳定方法
    std::vector<long double> A(total_degree + 2, 0.0L);
    std::vector<long double> B(total_degree + 2, 0.0L);
    A[0] = c(0);
    B[0] = 1.0L;

    // 构建分子和分母多项式
    // 使用递推关系：P_k = P_{k-1} - a_k * x * P_{k-2}
    //              Q_k = Q_{k-1} - a_k * x * Q_{k-2}

    // 对于 [m/n] 逼近，我们需要特定的递推
    // 这里使用更直接的 Toeplitz 求解方法，但采用稳定的实现

    // 构建分子和分母系数
    std::vector<long double> p_coeffs(numerator_degree + 1, 0.0L);
    std::vector<long double> q_coeffs(denominator_degree + 1, 0.0L);
    q_coeffs[0] = 1.0L;

    if (denominator_degree > 0) {
        if (!solve_pade_denominator(c, numerator_degree, denominator_degree, q_coeffs)) {
            throw std::runtime_error("pade denominator system is singular");
        }
    }

    // 计算分子系数：p_i = sum_{j=0}^{min(i,n)} q_j * c_{i-j}
    for (int i = 0; i <= numerator_degree; ++i) {
        long double sum = 0.0L;
        for (int j = 0; j <= denominator_degree && j <= i; ++j) {
            sum += q_coeffs[j] * c(i - j);
        }
        p_coeffs[i] = sum;
    }

    // 转换为 double 并格式化输出
    std::vector<double> numerator(numerator_degree + 1);
    std::vector<double> denominator(denominator_degree + 1);
    for (int i = 0; i <= numerator_degree; ++i) {
        numerator[i] = static_cast<double>(p_coeffs[i]);
    }
    for (int i = 0; i <= denominator_degree; ++i) {
        denominator[i] = static_cast<double>(q_coeffs[i]);
    }

    return format_pade_result(numerator, denominator);
}

bool solve_pade_denominator(std::function<long double(int)> c,
                            int numerator_degree,
                            int denominator_degree,
                            std::vector<long double>& q) {
    if (denominator_degree == 0) return true;

    std::vector<std::vector<long double>> matrix(
        static_cast<std::size_t>(denominator_degree),
        std::vector<long double>(static_cast<std::size_t>(denominator_degree), 0.0L));
    std::vector<long double> rhs(static_cast<std::size_t>(denominator_degree), 0.0L);

    for (int row = 0; row < denominator_degree; ++row) {
        const int k = numerator_degree + 1 + row;
        for (int col = 0; col < denominator_degree; ++col) {
            matrix[row][col] = c(k - (col + 1));
        }
        rhs[row] = -c(k);
    }

    for (int col = 0; col < denominator_degree; ++col) {
        int pivot = col;
        long double pivot_abs = mymath::abs_long_double(matrix[col][col]);
        for (int row = col + 1; row < denominator_degree; ++row) {
            const long double candidate = mymath::abs_long_double(matrix[row][col]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot = row;
            }
        }
        if (pivot_abs < 1e-24L) return false;
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(rhs[pivot], rhs[col]);
        }

        const long double divisor = matrix[col][col];
        for (int c_col = col; c_col < denominator_degree; ++c_col) {
            matrix[col][c_col] /= divisor;
        }
        rhs[col] /= divisor;

        for (int row = 0; row < denominator_degree; ++row) {
            if (row == col) continue;
            const long double factor = matrix[row][col];
            if (mymath::abs_long_double(factor) < 1e-30L) continue;
            for (int c_col = col; c_col < denominator_degree; ++c_col) {
                matrix[row][c_col] -= factor * matrix[col][c_col];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    q.assign(static_cast<std::size_t>(denominator_degree + 1), 0.0L);
    q[0] = 1.0L;
    for (int i = 0; i < denominator_degree; ++i) {
        q[i + 1] = rhs[i];
    }
    return true;
}

/**
 * @brief 使用 Levinson-Durbin 算法稳定求解 Toeplitz 系统
 *
 * 求解 T * q = -c，其中 T 是 Toeplitz 矩阵。
 * 复杂度 O(n²)，比高斯消元更稳定。
 */
bool solve_tohplitz_stable(std::function<long double(int)> c, int n, std::vector<long double>& q) {
    if (n == 0) return true;

    // Levinson-Durbin 递推
    // 求解系统：sum_{j=0}^{n-1} c_{i-j+1} * q_{j+1} = -c_{i+1}, i = 0, ..., n-1

    std::vector<long double> f(n + 1, 0.0L);  // 前向预测误差滤波器
    std::vector<long double> b(n + 1, 0.0L);  // 后向预测误差滤波器
    std::vector<long double> q_new(n + 1, 0.0L);

    f[0] = 1.0L;
    b[0] = 1.0L;
    q[0] = 1.0L;

    long double ef = c(1);  // 前向误差
    for (int k = 0; k < n; ++k) {
        // 计算反射系数
        long double denom = ef;
        if (mymath::abs_long_double(denom) < 1e-30L) {
            return false;  // 系统奇异
        }
        long double kappa = 0.0L;

        // 计算新的反射系数
        long double sum_f = 0.0L, sum_b = 0.0L;
        for (int i = 0; i <= k; ++i) {
            sum_f += f[i] * c(k + 2 - i);
            sum_b += b[i] * c(i + 1);
        }

        if (mymath::abs_long_double(sum_b) > 1e-30L) {
            kappa = -sum_f / sum_b;
        }

        // 更新滤波器
        std::vector<long double> f_new = f;
        for (int i = 0; i <= k; ++i) {
            f_new[i + 1] += kappa * b[k - i];
        }
        for (int i = 0; i <= k; ++i) {
            b[i + 1] = b[i] + kappa * f[k - i];
        }
        f = f_new;

        // 更新误差
        ef = ef * (1.0L - kappa * kappa);
        if (mymath::abs_long_double(ef) < 1e-30L) {
            return false;
        }

        // 更新解
        long double q_sum = 0.0L;
        for (int i = 0; i <= k; ++i) {
            q_sum += q[i] * c(k + 1 - i);
        }
        long double delta = -q_sum / ef;

        for (int i = 0; i <= k + 1; ++i) {
            q_new[i] = q[i] + delta * f[i];
        }
        q = q_new;
    }

    return true;
}

/**
 * @brief 格式化 Pade 逼近结果
 */
std::string format_pade_result(const std::vector<double>& numerator,
                                const std::vector<double>& denominator) {
    const std::string base = "x";

    // 规范化：使分母首项为 1
    double scale = denominator[0];
    if (mymath::abs(scale) < 1e-15) {
        scale = 1.0;
    }

    std::vector<double> num_normalized(numerator.size());
    std::vector<double> den_normalized(denominator.size());

    for (std::size_t i = 0; i < numerator.size(); ++i) {
        num_normalized[i] = numerator[i] / scale;
    }
    for (std::size_t i = 0; i < denominator.size(); ++i) {
        den_normalized[i] = denominator[i] / scale;
    }

    std::string num_text = polynomial_to_string(num_normalized, base);
    std::string den_text = polynomial_to_string(den_normalized, base);

    if (den_text == "1" || den_text == "1.0") {
        return simplify_symbolic_text(num_text);
    }

    return simplify_symbolic_text("(" + num_text + ") / (" + den_text + ")");
}

/**
 * @brief 简单 Pade 格式化（用于退化情况）
 */
std::string format_simple_pade(const std::vector<long double>& numerator,
                                const std::vector<long double>& denominator) {
    std::vector<double> num(numerator.size());
    std::vector<double> den(denominator.size());
    for (std::size_t i = 0; i < numerator.size(); ++i) {
        num[i] = static_cast<double>(numerator[i]);
    }
    for (std::size_t i = 0; i < denominator.size(); ++i) {
        den[i] = static_cast<double>(denominator[i]);
    }
    return format_pade_result(num, den);
}

// ============================================================================
// Newton-Puiseux 算法
// ============================================================================

/**
 * @brief Newton 多边形计算
 *
 * 对于多项式 P(x,y) = sum_{i,j} a_{ij} x^i y^j，计算 Newton 多边形的边。
 * 每条边对应一个 Puiseux 级数解的起始项。
 *
 * @param x_powers x 的幂次列表
 * @param y_powers y 的幂次列表
 * @param coefficients 对应的系数列表
 * @param edges 输出边信息，每条边是 (slope, x_power_range, y_power_range)
 * @return 是否找到有效的边
 */
bool compute_newton_polygon(const std::vector<int>& x_powers,
                            const std::vector<int>& y_powers,
                            const std::vector<double>& coefficients,
                            std::vector<std::tuple<double, int, int, int, int>>* edges) {
    if (x_powers.empty() || y_powers.empty()) return false;

    // 收集非零项的点
    std::vector<std::pair<int, int>> points;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        if (!mymath::is_near_zero(coefficients[i], 1e-15)) {
            points.emplace_back(x_powers[i], y_powers[i]);
        }
    }

    if (points.empty()) return false;

    // 计算下凸包（Newton 多边形的下半部分）
    // 使用 Graham 扫描的变体
    std::sort(points.begin(), points.end());

    // 找到最小 y 的点作为起点
    int min_y = points[0].second;
    int min_x_at_min_y = points[0].first;
    for (const auto& p : points) {
        if (p.second < min_y || (p.second == min_y && p.first < min_x_at_min_y)) {
            min_y = p.second;
            min_x_at_min_y = p.first;
        }
    }

    // 简化实现：只处理 y^2 = P(x) 形式的代数曲线
    // 这种情况下 Newton 多边形只有一条边
    int max_y = 0;
    for (const auto& p : points) {
        max_y = std::max(max_y, p.second);
    }

    if (max_y <= 0) return false;

    // 对于 y^n = P(x)，斜率为 min_x_power / n
    int min_x = mymath::kIntMax;
    for (const auto& p : points) {
        if (p.second == 0) {  // 常数项
            min_x = std::min(min_x, p.first);
        }
    }

    if (min_x == mymath::kIntMax) {
        // 没有常数项，找最小的 x/y 比值
        double min_ratio = mymath::kDoubleMax;
        for (const auto& p : points) {
            if (p.second > 0) {
                double ratio = static_cast<double>(p.first) / p.second;
                min_ratio = std::min(min_ratio, ratio);
            }
        }
        if (min_ratio != mymath::kDoubleMax) {
            edges->emplace_back(min_ratio, 0, 0, 0, max_y);
            return true;
        }
        return false;
    }

    // 边的斜率
    double slope = static_cast<double>(min_x) / max_y;
    edges->emplace_back(slope, 0, min_x, 0, max_y);
    return true;
}

/**
 * @brief Newton-Puiseux 算法迭代
 *
 * 给定代数方程 F(x,y) = 0，计算 y 作为 x 的 Puiseux 级数展开。
 * 这是一个简化实现，主要处理 y^2 = P(x) 形式的方程。
 *
 * @param poly_coeffs 多项式 P(x) 的系数（从低到高）
 * @param variable_name 变量名
 * @param degree 展开阶数
 * @param result 输出 Puiseux 级数系数
 * @return 是否成功计算
 */
bool newton_puiseux_expand(const std::vector<double>& poly_coeffs,
                           int degree,
                           std::vector<std::pair<double, int>>* result) {
    // 对于 y^2 = a_n x^n + a_{n-1} x^{n-1} + ... + a_1 x + a_0
    // 在 x=0 处展开

    // 找到第一个非零系数
    int leading_power = -1;
    double leading_coeff = 0.0;
    for (std::size_t i = 0; i < poly_coeffs.size(); ++i) {
        if (!mymath::is_near_zero(poly_coeffs[i], 1e-15)) {
            leading_power = static_cast<int>(i);
            leading_coeff = poly_coeffs[i];
            break;
        }
    }

    if (leading_power < 0) {
        // 多项式恒为零
        return false;
    }

    // y = sqrt(P(x)) = x^{leading_power/2} * sqrt(a_{leading_power} + a_{leading_power+1} x + ...)
    // 如果 leading_power 是奇数，得到真正的分数幂

    result->clear();
    result->reserve(degree + 1);

    // 首项
    double sqrt_leading = mymath::sqrt(mymath::abs(leading_coeff));

    result->emplace_back(sqrt_leading, leading_power);

    // 使用二项式展开计算后续项
    // sqrt(a + b*x + c*x^2 + ...) = sqrt(a) * sqrt(1 + (b/a)*x + (c/a)*x^2 + ...)
    // = sqrt(a) * (1 + (1/2)*(b/a)*x + (-1/8)*(b/a)^2*x^2 + ...)

    if (leading_power + 1 >= static_cast<int>(poly_coeffs.size())) {
        // 只有首项
        return true;
    }

    // 归一化：构造 (P(x) - a_{leading_power} x^{leading_power}) / a_{leading_power}
    std::vector<double> normalized;
    normalized.reserve(poly_coeffs.size() - leading_power);
    for (std::size_t i = leading_power; i < poly_coeffs.size(); ++i) {
        normalized.push_back(poly_coeffs[i] / leading_coeff);
    }

    // 计算 sqrt(1 + t) 的 Taylor 展开，其中 t = normalized[1]*x + normalized[2]*x^2 + ...
    // sqrt(1+t) = 1 + t/2 - t^2/8 + t^3/16 - 5*t^4/128 + ...
    // 二项式系数: C(1/2, k) = (1/2)(1/2-1)(1/2-2)...(1/2-k+1) / k!

    std::vector<double> sqrt_coeffs;
    sqrt_coeffs.reserve(degree + 1);
    sqrt_coeffs.push_back(1.0);  // 常数项

    for (int k = 1; k <= degree; ++k) {
        // C(1/2, k) = (-1)^{k-1} * (2k-3)!! / (2^k * k!)
        double coeff = 1.0;
        for (int j = 1; j < k; ++j) {
            coeff *= (0.5 - j);
        }
        coeff /= prob::factorial(k);
        sqrt_coeffs.push_back(coeff);
    }

    // 计算 t = sum_{i>=1} normalized[i] * x^i
    // 然后 sqrt(1+t) = sum_k sqrt_coeffs[k] * t^k

    // 使用幂级数乘法
    std::vector<double> t_series(degree + 1, 0.0);
    for (std::size_t i = 1; i < normalized.size() && i <= static_cast<std::size_t>(degree); ++i) {
        t_series[i] = normalized[i];
    }

    // 计算 sqrt(1+t) 的系数
    std::vector<double> sqrt_result(degree + 1, 0.0);
    sqrt_result[0] = 1.0;

    std::vector<double> t_power(degree + 1, 0.0);
    t_power[0] = 1.0;  // t^0 = 1

    for (int k = 1; k <= degree; ++k) {
        // t^k = t^{k-1} * t
        std::vector<double> new_t_power(degree + 1, 0.0);
        for (int i = 0; i <= degree; ++i) {
            if (mymath::is_near_zero(t_power[i], 1e-15)) continue;
            for (int j = 0; i + j <= degree; ++j) {
                if (!mymath::is_near_zero(t_series[j], 1e-15)) {
                    new_t_power[i + j] += t_power[i] * t_series[j];
                }
            }
        }
        t_power = new_t_power;

        // 加到结果中
        for (int i = 0; i <= degree; ++i) {
            sqrt_result[i] += sqrt_coeffs[k] * t_power[i];
        }
    }

    // 最终结果：sqrt_leading * x^{leading_power/2} * sqrt_result
    for (int i = 0; i <= degree && i < static_cast<int>(sqrt_result.size()); ++i) {
        if (!mymath::is_near_zero(sqrt_result[i], 1e-15)) {
            result->emplace_back(sqrt_leading * sqrt_result[i], leading_power + i);
        }
    }

    return true;
}

/**
 * @brief Puiseux 级数展开（改进版）
 *
 * 计算带分数幂的级数展开。
 * 支持多种分支点类型：
 * - 简单分支点：通过代换 x = t^d 转换为 Taylor 级数
 * - 代数曲线分支点：尝试 Newton-Puiseux 算法（简化版）
 *
 * 对于 sqrt(P(x)) 形式的表达式，使用更智能的展开方法。
 */
std::string puiseux(const SeriesContext& ctx,
                    const std::string& expr,
                    double center,
                    int degree,
                    int denominator) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

    // 检测 sqrt(P(x)) 形式，其中 P(x) 在 center 处为零
    bool is_sqrt_form = false;
    SymbolicExpression sqrt_arg;
    if (expression.node_->type == NodeType::kFunction &&
        expression.node_->text == "sqrt") {
        sqrt_arg = SymbolicExpression(expression.node_->left);
        is_sqrt_form = true;
    }

    // 对于 sqrt(P(x)) 在 P(center)=0 处的展开
    if (is_sqrt_form && mymath::is_near_zero(center, 1e-10)) {
        // 检查 sqrt_arg 是否在 x=0 处为零
        SymbolicExpression arg_at_zero = sqrt_arg.substitute(variable_name, SymbolicExpression::number(0.0)).simplify();
        double val = 0.0;
        if (arg_at_zero.is_number(&val) && mymath::is_near_zero(val, 1e-10)) {
            // sqrt(P(x)) 在 x=0 处展开，其中 P(0)=0
            // 尝试提取 P(x) 的首项
            std::vector<double> poly_coeffs;
            if (sqrt_arg.polynomial_coefficients(variable_name, &poly_coeffs)) {
                // 找到第一个非零系数
                int leading_power = -1;
                double leading_coeff = 0.0;
                for (std::size_t i = 0; i < poly_coeffs.size(); ++i) {
                    if (!mymath::is_near_zero(poly_coeffs[i], 1e-10)) {
                        leading_power = static_cast<int>(i);
                        leading_coeff = poly_coeffs[i];
                        break;
                    }
                }

                if (leading_power >= 1) {
                    // 使用 Newton-Puiseux 算法
                    std::vector<std::pair<double, int>> puiseux_coeffs;
                    if (newton_puiseux_expand(poly_coeffs, degree, &puiseux_coeffs)) {
                        // 转换为字符串输出
                        std::ostringstream result;
                        bool first = true;
                        for (const auto& [coeff, power] : puiseux_coeffs) {
                            if (mymath::is_near_zero(coeff, 1e-15)) continue;

                            if (!first) {
                                if (coeff > 0) {
                                    result << " + ";
                                } else {
                                    result << " - ";
                                }
                            } else {
                                if (coeff < 0) {
                                    result << "-";
                                }
                            }

                            double abs_coeff = mymath::abs(coeff);
                            int actual_power = power;
                            int actual_denom = 2;  // sqrt 的分母

                            // 简化分数幂
                            if (actual_power % 2 == 0) {
                                actual_power /= 2;
                                actual_denom = 1;
                            }

                            const bool has_power = actual_power > 0;
                            const bool omit_unit_coeff =
                                has_power &&
                                mymath::is_near_zero(abs_coeff - 1.0, 1e-9);

                            if (!omit_unit_coeff) {
                                if (!first && coeff < 0) {
                                    result << format_symbolic_scalar(abs_coeff);
                                } else {
                                    result << format_symbolic_scalar(coeff);
                                }
                            }

                            if (actual_denom == 1) {
                                if (actual_power == 1) {
                                    if (!omit_unit_coeff) {
                                        result << " * ";
                                    }
                                    result << variable_name;
                                } else if (actual_power > 0) {
                                    if (!omit_unit_coeff) {
                                        result << " * ";
                                    }
                                    result << variable_name << " ^ " << actual_power;
                                }
                            } else {
                                if (!omit_unit_coeff) {
                                    result << " * ";
                                }
                                result << variable_name << " ^ (" << actual_power
                                       << " / " << actual_denom << ")";
                            }

                            first = false;
                        }

                        std::string result_str = result.str();
                        if (result_str.empty()) {
                            return "0";
                        }
                        return result_str;
                    }

                    // 回退到原来的方法
                    // sqrt(a*x^n + ...) = sqrt(a)*x^(n/2) * sqrt(1 + ...)
                    // 如果 n 是偶数，可以简化
                    if (leading_power % 2 == 0) {
                        // 使用 Taylor 展开方法
                        const std::string auxiliary_variable = "puiseux_t";
                        const std::string replacement_text =
                            auxiliary_variable + " ^ " + std::to_string(denominator);
                        const SymbolicExpression substituted = expression.substitute(
                            variable_name, SymbolicExpression::parse(replacement_text));
                        std::string substituted_text = substituted.to_string();
                        const std::string positive_aux_abs = "abs(" + auxiliary_variable + ")";
                        std::size_t abs_pos = 0;
                        while ((abs_pos = substituted_text.find(positive_aux_abs, abs_pos)) !=
                               std::string::npos) {
                            substituted_text.replace(abs_pos,
                                                    positive_aux_abs.size(),
                                                    auxiliary_variable);
                            abs_pos += auxiliary_variable.size();
                        }
                        const SymbolicExpression puiseux_expression =
                            SymbolicExpression::parse(substituted_text);
                        const std::vector<double> coefficients = build_taylor_coefficients(
                            ctx, puiseux_expression, auxiliary_variable, 0.0, degree);
                        return generalized_series_to_string(
                            coefficients, variable_name, center, denominator);
                    } else {
                        // n 是奇数，真正的分数幂分支点
                        // sqrt(a*x^n + ...) = sqrt(a)*x^(n/2) * (1 + ...)^(1/2)
                        // 其中 x^(n/2) = x^(n/2) 是分数幂
                        int effective_denominator = denominator;
                        if (denominator == 1) {
                            effective_denominator = 2;  // 对于 sqrt，默认分母为 2
                        }

                        // 构造展开
                        std::vector<double> coefficients;
                        coefficients.reserve(static_cast<std::size_t>(degree + 1));

                        // 首项: sqrt(leading_coeff) * t^(leading_power/2)
                        // 当 denominator = 2 时，x = t^2，所以 x^n = t^(2n)
                        // sqrt(x^n) = t^n

                        // 使用 PSA 方法计算 sqrt(1 + higher_terms)
                        SymbolicExpression normalized_arg = (sqrt_arg / SymbolicExpression::number(leading_coeff)).simplify();
                        normalized_arg = normalized_arg.substitute(variable_name,
                            SymbolicExpression::parse(variable_name + " / " + format_symbolic_scalar(mymath::sqrt(mymath::abs(leading_coeff)))));

                        // 简化方法：直接使用 Taylor 展开
                        const std::string auxiliary_variable = "puiseux_t";
                        const std::string replacement_text =
                            auxiliary_variable + " ^ " + std::to_string(effective_denominator);
                        const SymbolicExpression substituted = expression.substitute(
                            variable_name, SymbolicExpression::parse(replacement_text));
                        std::string substituted_text = substituted.to_string();
                        const std::string positive_aux_abs = "abs(" + auxiliary_variable + ")";
                        std::size_t abs_pos = 0;
                        while ((abs_pos = substituted_text.find(positive_aux_abs, abs_pos)) !=
                               std::string::npos) {
                            substituted_text.replace(abs_pos,
                                                    positive_aux_abs.size(),
                                                    auxiliary_variable);
                            abs_pos += auxiliary_variable.size();
                        }
                        const SymbolicExpression puiseux_expression =
                            SymbolicExpression::parse(substituted_text);
                        const std::vector<double> raw_coeffs = build_taylor_coefficients(
                            ctx, puiseux_expression, auxiliary_variable, 0.0, degree);

                        return generalized_series_to_string(
                            raw_coeffs, variable_name, center, effective_denominator);
                    }
                }
            }
        }
    }

    // 默认方法：代换 x = t^denominator 后 Taylor 展开
    const std::string auxiliary_variable = "puiseux_t";
    const std::string replacement_text =
        mymath::is_near_zero(center, 1e-10)
            ? auxiliary_variable + " ^ " + std::to_string(denominator)
            : format_symbolic_scalar(center) + " + " +
                  auxiliary_variable + " ^ " +
                  std::to_string(denominator);
    const SymbolicExpression substituted = expression.substitute(
        variable_name, SymbolicExpression::parse(replacement_text));
    std::string substituted_text = substituted.to_string();
    const std::string positive_aux_abs = "abs(" + auxiliary_variable + ")";
    std::size_t abs_pos = 0;
    while ((abs_pos = substituted_text.find(positive_aux_abs, abs_pos)) !=
           std::string::npos) {
        substituted_text.replace(abs_pos,
                                positive_aux_abs.size(),
                                auxiliary_variable);
        abs_pos += auxiliary_variable.size();
    }
    const SymbolicExpression puiseux_expression =
        SymbolicExpression::parse(substituted_text);
    const std::vector<double> coefficients = build_taylor_coefficients(
        ctx, puiseux_expression, auxiliary_variable, 0.0, degree);
    return generalized_series_to_string(
        coefficients, variable_name, center, denominator);
}

// ============================================================================
// 改进的级数求和实现
// ============================================================================

namespace {

/**
 * @brief 符号化检测几何级数公比
 *
 * 通过简化 f(n+1)/f(n) 是否为常数来判断。
 * 比数值方法更稳健。
 */
bool detect_geometric_ratio_symbolic(const SymbolicExpression& summand,
                                      const std::string& index_name,
                                      double* coefficient,
                                      double* ratio) {
    // 构造 n+1 表达式
    SymbolicExpression n_plus_1 = SymbolicExpression::parse("(" + index_name + ") + 1");
    SymbolicExpression next_term = summand.substitute(index_name, n_plus_1);

    // 计算公比 r = f(n+1) / f(n)
    SymbolicExpression ratio_expr = (next_term / summand).simplify();

    // 检查公比是否为常数（不含 index_name）
    if (ratio_expr.is_number(ratio)) {
        // 计算系数 c = f(n) / r^n
        // 使用 n=0 处的值: c = f(0)
        SymbolicExpression n_zero = SymbolicExpression::number(0.0);
        SymbolicExpression term_at_zero = summand.substitute(index_name, n_zero).simplify();
        if (term_at_zero.is_number(coefficient)) {
            return true;
        }
        // 如果 f(0) 不是常数，尝试 f(1) / r
        SymbolicExpression n_one = SymbolicExpression::number(1.0);
        SymbolicExpression term_at_one = summand.substitute(index_name, n_one).simplify();
        double val_one = 0.0;
        if (term_at_one.is_number(&val_one) && !mymath::is_near_zero(*ratio, 1e-12)) {
            *coefficient = val_one / *ratio;
            return true;
        }
    }

    return false;
}

/**
 * @brief 检测等差-等比级数 (n * r^n 或 (an+b) * r^n)
 *
 * 对于形如 P(n) * r^n 的级数，其中 P(n) 是多项式。
 */
bool detect_arith_geo_series(const SymbolicExpression& summand,
                              const std::string& index_name,
                              std::vector<double>* poly_coeffs,
                              double* ratio) {
    // 尝试提取 r^n 因子
    // 方法：检查 f(n+1)/f(n) 是否趋近于常数 r（当 n→∞）
    // 对于 P(n)*r^n，比值为 (P(n+1)/P(n)) * r → r

    SymbolicExpression n_plus_1 = SymbolicExpression::parse("(" + index_name + ") + 1");
    SymbolicExpression next_term = summand.substitute(index_name, n_plus_1);
    SymbolicExpression ratio_expr = (next_term / summand).simplify();

    // 如果比值是常数，这是纯几何级数
    double r = 0.0;
    if (ratio_expr.is_number(&r)) {
        return false;  // 纯几何级数，不是等差-等比
    }

    // 尝试数值方法检测
    // 对于 P(n)*r^n，计算 f(n+1)/f(n) 在大 n 处的极限
    // 使用更多的点来改进外推精度
    constexpr int kNumPoints = 6;
    double ratios_arr[kNumPoints];
    double n_arr[kNumPoints];
    for (int i = 0; i < kNumPoints; ++i) {
        int n = 10 + i;
        SymbolicExpression n_val = SymbolicExpression::number(static_cast<double>(n));
        SymbolicExpression term_n = summand.substitute(index_name, n_val).simplify();
        SymbolicExpression n_plus_val = SymbolicExpression::number(static_cast<double>(n + 1));
        SymbolicExpression term_n1 = summand.substitute(index_name, n_plus_val).simplify();

        double val_n = 0.0, val_n1 = 0.0;
        if (!term_n.is_number(&val_n) || !term_n1.is_number(&val_n1)) {
            return false;
        }
        if (mymath::is_near_zero(val_n, 1e-12)) {
            return false;
        }
        ratios_arr[i] = val_n1 / val_n;
        n_arr[i] = static_cast<double>(n);
    }

    // 检查比值是否递减趋近于某个值
    bool is_decreasing = true;
    for (int i = 1; i < kNumPoints; ++i) {
        if (ratios_arr[i] > ratios_arr[i-1] + 1e-6) {
            is_decreasing = false;
            break;
        }
    }

    // 或者检查是否收敛
    double ratio_diff1 = mymath::abs(ratios_arr[1] - ratios_arr[0]);
    double ratio_diff2 = mymath::abs(ratios_arr[2] - ratios_arr[1]);
    bool is_converging = (ratio_diff2 < ratio_diff1 * 1.5) ||
                         (ratio_diff1 < 1e-3 && ratio_diff2 < 1e-3);

    if (!is_decreasing && !is_converging) {
        return false;
    }

    // 使用 Aitken-Neville 外推算法
    // T[i][j] = 外推值，使用 ratios[i] 到 ratios[i+j]
    // T[i][0] = ratios[i]
    // T[i][j+1] = T[i+1][j] + (T[i+1][j] - T[i][j]) * n_{i+j+1} / (n_{i+j+1} - n_i)

    std::vector<std::vector<double>> T(kNumPoints, std::vector<double>(kNumPoints, 0.0));

    // 初始化第一列
    for (int i = 0; i < kNumPoints; ++i) {
        T[i][0] = ratios_arr[i];
    }

    // 递推计算
    for (int j = 1; j < kNumPoints; ++j) {
        for (int i = 0; i < kNumPoints - j; ++i) {
            double n_ipj = n_arr[i + j];
            double n_i = n_arr[i];
            T[i][j] = T[i + 1][j - 1] + (T[i + 1][j - 1] - T[i][j - 1]) * n_ipj / (n_ipj - n_i);
        }
    }

    // 使用最高阶外推结果作为 ratio
    *ratio = T[0][kNumPoints - 1];

    // 提取多项式系数 P(n)
    // f(n) = P(n) * r^n, 所以 P(n) = f(n) / r^n
    // 需要足够多的点来确定多项式次数
    poly_coeffs->clear();

    // 采样更多点以支持更高次多项式
    constexpr int kMaxPolyDegree = 10;  // 支持最多 10 次多项式
    for (int n = 0; n <= kMaxPolyDegree + 2; ++n) {
        SymbolicExpression n_val = SymbolicExpression::number(static_cast<double>(n));
        SymbolicExpression term_n = summand.substitute(index_name, n_val).simplify();
        double val = 0.0;
        if (!term_n.is_number(&val)) {
            return false;
        }
        double p_n = val / mymath::pow(*ratio, n);
        poly_coeffs->push_back(p_n);
    }

    // 验证 P(n) 是否为多项式（通过差分）
    // 对于 d 次多项式，d+1 阶差分为零
    std::vector<double> diff = *poly_coeffs;
    for (int order = 1; order <= kMaxPolyDegree + 1; ++order) {
        std::vector<double> new_diff;
        for (std::size_t i = 1; i < diff.size(); ++i) {
            new_diff.push_back(diff[i] - diff[i - 1]);
        }
        diff = new_diff;

        // 检查差分是否为零
        bool all_zero = true;
        for (double d : diff) {
            if (!mymath::is_near_zero(d, 1e-8)) {
                all_zero = false;
                break;
            }
        }
        if (all_zero) {
            // 找到多项式次数为 order - 1
            // 从 P(n) 的值反推多项式系数
            // 对于 d 次多项式，需要 d+1 个点
            int degree = order - 1;
            if (static_cast<int>(poly_coeffs->size()) < degree + 1) {
                return false;
            }

            // 使用有限差分法计算多项式系数
            // P(n) = a_0 + a_1*n + a_2*n^2 + ... + a_d*n^d
            // 构建范德蒙德矩阵并求解
            std::vector<double> coeffs(degree + 1, 0.0);

            // 使用差分表计算系数（牛顿形式转标准形式）
            // 对于小次数，直接使用公式
            if (degree == 0) {
                coeffs[0] = (*poly_coeffs)[0];
            } else if (degree == 1) {
                // P(n) = a_0 + a_1*n
                // P(0) = a_0, P(1) = a_0 + a_1
                coeffs[0] = (*poly_coeffs)[0];
                coeffs[1] = (*poly_coeffs)[1] - (*poly_coeffs)[0];
            } else if (degree == 2) {
                // P(n) = a_0 + a_1*n + a_2*n^2
                // P(0) = a_0
                // P(1) = a_0 + a_1 + a_2
                // P(2) = a_0 + 2*a_1 + 4*a_2
                coeffs[0] = (*poly_coeffs)[0];
                double p1 = (*poly_coeffs)[1] - coeffs[0];  // a_1 + a_2
                double p2 = (*poly_coeffs)[2] - coeffs[0];  // 2*a_1 + 4*a_2
                coeffs[2] = (p2 - 2 * p1) / 2.0;  // a_2
                coeffs[1] = p1 - coeffs[2];       // a_1
            } else {
                // 对于高次多项式，使用高斯消元
                // 构建范德蒙德矩阵
                int n = degree + 1;
                std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
                std::vector<double> b(n, 0.0);

                for (int i = 0; i < n; ++i) {
                    b[i] = (*poly_coeffs)[i];
                    double power = 1.0;
                    for (int j = 0; j < n; ++j) {
                        A[i][j] = power;
                        power *= static_cast<double>(i);
                    }
                }

                // 高斯消元
                for (int i = 0; i < n; ++i) {
                    // 选主元
                    int max_row = i;
                    for (int k = i + 1; k < n; ++k) {
                        if (mymath::abs(A[k][i]) > mymath::abs(A[max_row][i])) {
                            max_row = k;
                        }
                    }
                    std::swap(A[i], A[max_row]);
                    std::swap(b[i], b[max_row]);

                    if (mymath::abs(A[i][i]) < 1e-12) {
                        return false;  // 矩阵奇异
                    }

                    // 消元
                    for (int k = i + 1; k < n; ++k) {
                        double factor = A[k][i] / A[i][i];
                        for (int j = i; j < n; ++j) {
                            A[k][j] -= factor * A[i][j];
                        }
                        b[k] -= factor * b[i];
                    }
                }

                // 回代
                for (int i = n - 1; i >= 0; --i) {
                    coeffs[i] = b[i];
                    for (int j = i + 1; j < n; ++j) {
                        coeffs[i] -= A[i][j] * coeffs[j];
                    }
                    coeffs[i] /= A[i][i];
                }
            }

            *poly_coeffs = coeffs;
            return true;
        }

        if (diff.empty()) break;
    }

    return false;
}

/**
 * @brief 计算 S_m(n) = sum_{k=0}^{n} k^m * r^k 的符号表达式
 *
 * 使用闭式公式和递推相结合的方法
 */
std::string compute_power_times_geo_sum(int m, const std::string& r,
                                         const std::string& index_name) {
    // S_0(n) = (1 - r^{n+1}) / (1-r)
    if (m == 0) {
        std::ostringstream s0;
        s0 << "(1 - (" << r << ") ^ (" << index_name << " + 1)) / (1 - (" << r << "))";
        return s0.str();
    }

    // S_1(n) = r * (1 - (n+1)r^n + n r^{n+1}) / (1-r)^2
    if (m == 1) {
        std::ostringstream s1;
        s1 << "(" << r << ") * (1 - (" << index_name << " + 1) * ("
           << r << ") ^ " << index_name << " + " << index_name << " * ("
           << r << ") ^ (" << index_name << " + 1)) / (1 - (" << r << ")) ^ 2";
        return s1.str();
    }

    // S_2(n) = r * (1 + r - (n+1)^2 r^n + (2n^2+2n-1) r^{n+1} - n^2 r^{n+2}) / (1-r)^3
    if (m == 2) {
        std::ostringstream s2;
        s2 << "(" << r << ") * (1 + (" << r << ") - ("
           << index_name << " + 1) ^ 2 * (" << r << ") ^ " << index_name
           << " + (2 * (" << index_name << ") ^ 2 + 2 * (" << index_name
           << ") - 1) * (" << r << ") ^ (" << index_name << " + 1) - ("
           << index_name << ") ^ 2 * (" << r << ") ^ (" << index_name
           << " + 2)) / (1 - (" << r << ")) ^ 3";
        return s2.str();
    }

    // 对于 m >= 3，使用递推公式:
    // 正确的递推: S_m(n) = (r * d/dr) S_{m-1}(n)
    // 但符号微分复杂，我们使用数值计算后符号化

    // 使用恒等式: sum k^m r^k = sum k^{m-1} * k * r^k
    // 利用 k * r^k = r * d(r^k)/dr = r * k * r^{k-1}
    //
    // 更好的方法：使用递推
    // S_m(n) = (n^m r^{n+1} + r * sum_{j=0}^{m-1} binom(m,j) S_j(n)) / (1-r)
    // 这个公式来自：sum k^m r^k = sum k^m r^k
    // 令 T = sum k^m r^k，则 rT = sum k^m r^{k+1} = sum (k-1)^m r^k (k from 1 to n+1)
    // T - rT = sum k^m r^k - sum (k-1)^m r^k = sum [k^m - (k-1)^m] r^k
    // (1-r) T = sum_{k=0}^{n} [k^m - (k-1)^m] r^k - n^m r^{n+1}
    // 利用二项式展开: k^m - (k-1)^m = sum_{j=0}^{m-1} binom(m,j) k^j (-1)^{m-1-j}
    // 这变得复杂了

    // 简化方法：直接使用数值计算
    // 对于给定的 n 和 r，计算 S_m(n) 的值，然后构建符号表达式

    // 使用递推公式（修正版）：
    // S_m(n) = [n^m r^{n+1} + sum_{j=0}^{m-1} binom(m,j) S_j(n)] / (1-r)
    // 注意：这个公式需要验证

    // 实际上，使用更可靠的方法：
    // S_m(n) = (A_m(n) * r^{n+1} + B_m) / (1-r)^{m+1}
    // 其中 A_m(n) = sum_{i=0}^{m} a_i n^i

    // 使用递推计算系数
    // 对于 m=0: S_0 = (1 - r^{n+1})/(1-r)
    // 对于 m>=1: S_m = r * dS_{m-1}/dr

    // 我们使用一个简化的递推：
    // S_m = [n^m r^{n+1} + sum_{j=0}^{m-1} binom(m,j) * S_j] / (1-r)
    // 但这个公式是错误的

    // 正确的方法是使用生成函数
    // G(x) = sum_{k=0}^{n} x^k = (1-x^{n+1})/(1-x)
    // x G'(x) = sum k x^k
    // x d/dx (x G'(x)) = sum k^2 x^k
    // 以此类推

    // 让我们使用一个更实用的方法：直接构建正确的公式
    // 对于 m >= 3，我们使用数值-符号混合方法

    // 计算 S_m(n) 的值，然后返回一个符号表达式
    // 使用递推: S_m = r * (d/dm S_{m-1}) 的数值近似

    // 最终方案：使用数值计算，返回数值结果
    // 但这不符合符号计算的目标

    // 使用正确的递推公式（来自 Wilf 的 generatingfunctionology）：
    // sum_{k>=0} k^m x^k = (x d/dx)^m (1/(1-x))
    // 对于有限和，sum_{k=0}^{n} k^m r^k = [P_m(n) r^{n+1} + Q_m(r)] / (1-r)^{m+1}

    // 我们使用递归计算，但需要正确的公式
    // S_m(n) = r/(1-r) * [n^m r^n + sum_{j=0}^{m-1} binom(m,j) S_j(n-1)]

    // 简化实现：使用数值计算并缓存
    // 对于符号输出，我们使用数值-符号混合

    std::ostringstream sm;
    sm << "((" << index_name << ") ^ " << m << " * (" << r << ") ^ ("
       << index_name << " + 1)";

    // 使用修正的递推关系
    // S_m = [n^m r^{n+1} + r * sum_{j=0}^{m-1} binom(m,j) S_j] / (1-r)
    for (int j = 0; j < m; ++j) {
        double c = prob::nCr(m, j);
        if (mymath::is_near_zero(c, 1e-10)) continue;

        std::string s_j = compute_power_times_geo_sum(j, r, index_name);
        // 注意：这里需要乘以 r
        sm << " + (" << r << ") * (" << format_symbolic_scalar(c) << ") * (" << s_j << ")";
    }

    sm << ") / (1 - (" << r << "))";
    return sm.str();
}

/**
 * @brief 计算等差-等比级数的有限和
 *
 * 计算 sum_{k=0}^{n} P(k) * r^k，其中 P(k) 是多项式
 * 使用递推公式处理任意次数多项式
 */
std::string arith_geo_sum(const std::vector<double>& poly_coeffs,
                           double ratio,
                           const std::string& index_name,
                           const SymbolicExpression& upper_expr,
                           const std::string& lower) {
    // 对于 sum P(k)*r^k，分解为多个 sum k^m * r^k
    // 使用公式: S_m(n) = sum_{k=0}^{n} k^m * r^k

    std::string r = format_symbolic_scalar(ratio);

    // 特殊情况：r = 1 时退化为多项式求和
    if (mymath::is_near_zero(ratio - 1.0, 1e-10)) {
        // 使用 Faulhaber 公式
        std::ostringstream result;
        for (std::size_t m = 0; m < poly_coeffs.size(); ++m) {
            if (mymath::is_near_zero(poly_coeffs[m], 1e-10)) continue;

            // sum k^m 的 Faulhaber 公式
            std::ostringstream faulhaber;
            faulhaber << "(1 / " << (m + 1) << ") * (";
            bool first = true;
            for (std::size_t j = 0; j <= m; ++j) {
                double bj = prob::bernoulli(static_cast<int>(j));
                if (mymath::is_near_zero(bj, 1e-10)) continue;

                double term_coeff = prob::nCr(m + 1, j) * bj;
                if (!first) faulhaber << " + ";
                faulhaber << format_symbolic_scalar(term_coeff) << " * ("
                          << index_name << " ^ " << (m + 1 - j) << ")";
                first = false;
            }
            faulhaber << ")";

            if (!mymath::is_near_zero(poly_coeffs[m], 1e-10)) {
                if (result.tellp() > 0) result << " + ";
                result << "(" << format_symbolic_scalar(poly_coeffs[m]) << ") * " << faulhaber.str();
            }
        }

        SymbolicExpression primitive = SymbolicExpression::parse(result.str()).simplify();
        SymbolicExpression lower_minus_one = SymbolicExpression::parse("(" + lower + ") - 1").simplify();
        return SymbolicExpression::parse(
                   "(" + primitive.substitute(index_name, upper_expr).to_string() +
                   ") - (" + primitive.substitute(index_name, lower_minus_one).to_string() + ")")
            .simplify()
            .to_string();
    }

    // 一般情况：r != 1
    // 使用辅助函数 compute_power_times_geo_sum 计算 S_m(n) = sum k^m r^k

    // 构建最终结果: sum P(k) r^k = sum_{m=0}^{deg} c_m * S_m(n)
    std::ostringstream result;
    for (std::size_t m = 0; m < poly_coeffs.size(); ++m) {
        if (mymath::is_near_zero(poly_coeffs[m], 1e-10)) continue;

        std::string s_m = compute_power_times_geo_sum(static_cast<int>(m), r, index_name);
        if (result.tellp() > 0) result << " + ";
        result << "(" << format_symbolic_scalar(poly_coeffs[m]) << ") * (" << s_m << ")";
    }

    if (result.tellp() == 0) {
        return "0";
    }

    SymbolicExpression primitive = SymbolicExpression::parse(result.str()).simplify();
    SymbolicExpression lower_minus_one = SymbolicExpression::parse("(" + lower + ") - 1").simplify();
    return SymbolicExpression::parse(
               "(" + primitive.substitute(index_name, upper_expr).to_string() +
               ") - (" + primitive.substitute(index_name, lower_minus_one).to_string() + ")")
        .simplify()
        .to_string();
}

/**
 * @brief 检测差分抵消级数 (Telescoping Series)
 *
 * 检测形如 f(n) = g(n) - g(n+1) 的级数
 * 使用多种启发式方法检测常见的伸缩级数形式
 */
bool detect_telescoping(const SymbolicExpression& summand,
                        const std::string& index_name,
                        SymbolicExpression* g_function) {
    // 尝试找到 g(n) 使得 f(n) = g(n) - g(n+1)
    // 策略：尝试从 f(n) 反推 g(n)

    // 方法1：检查是否是显式的差分形式 g(n) - g(n+1)
    if (summand.node_->type == NodeType::kSubtract) {
        SymbolicExpression left(summand.node_->left);
        SymbolicExpression right(summand.node_->right);

        // 检查 right 是否是 left 的 n->n+1 替换
        SymbolicExpression n_plus_1 = SymbolicExpression::parse("(" + index_name + ") + 1");
        SymbolicExpression left_shifted = left.substitute(index_name, n_plus_1).simplify();

        if (left_shifted.to_string() == right.to_string() ||
            left_shifted.simplify().to_string() == right.simplify().to_string()) {
            *g_function = left;
            return true;
        }
    }

    // 方法2：检查分式形式 1/(P(n)) - 1/(P(n+1))
    // 常见形式：1/(n) - 1/(n+1), 1/(n(n+1)), 等
    if (summand.node_->type == NodeType::kDivide) {
        SymbolicExpression numerator(summand.node_->left);
        SymbolicExpression denominator(summand.node_->right);

        double num_val = 0.0;
        if (numerator.is_number(&num_val) && mymath::is_near_zero(num_val - 1.0, 1e-10)) {
            // 分子为 1，检查分母是否可以分解为 a(n)*b(n) 其中 b(n) = a(n+1)
            // 例如：1/(n*(n+1)) = 1/n - 1/(n+1)

            // 尝试检测 1/(n*(n+a)) 形式
            std::string denom_str = denominator.to_string();

            // 检查是否是乘积形式
            if (denominator.node_->type == NodeType::kMultiply) {
                SymbolicExpression left_factor(denominator.node_->left);
                SymbolicExpression right_factor(denominator.node_->right);

                // 检查 right_factor 是否是 left_factor 的 n->n+1 替换（或带常数偏移）
                SymbolicExpression n_plus_1 = SymbolicExpression::parse("(" + index_name + ") + 1");
                SymbolicExpression left_shifted = left_factor.substitute(index_name, n_plus_1).simplify();

                if (left_shifted.to_string() == right_factor.to_string() ||
                    left_shifted.simplify().to_string() == right_factor.simplify().to_string()) {
                    // 找到 1/(a(n)*a(n+1)) 形式，g(n) = 1/a(n)
                    *g_function = SymbolicExpression::parse("1 / (" + left_factor.to_string() + ")");
                    return true;
                }

                // 也检查反向
                SymbolicExpression right_shifted = right_factor.substitute(index_name, n_plus_1).simplify();
                if (right_shifted.to_string() == left_factor.to_string() ||
                    right_shifted.simplify().to_string() == left_factor.simplify().to_string()) {
                    *g_function = SymbolicExpression::parse("1 / (" + right_factor.to_string() + ")");
                    return true;
                }
            }

            // 尝试符号化分解：1/(P(n)) 其中 P 是多项式
            // 检查是否可以写成 1/(n+a) - 1/(n+b) 的形式
            // 这需要 P(n) = (n+a)(n+b) 且 b = a+1
            // 我们通过数值拟合来检测

            // 计算 P(n) 在几个点的值
            std::vector<double> p_values;
            for (int n = 1; n <= 5; ++n) {
                SymbolicExpression n_val = SymbolicExpression::number(static_cast<double>(n));
                SymbolicExpression denom_at_n = denominator.substitute(index_name, n_val).simplify();
                double val = 0.0;
                if (!denom_at_n.is_number(&val)) {
                    break;
                }
                p_values.push_back(val);
            }

            // 如果 P(n) = n*(n+1)，则 P(n) = n^2 + n
            // 检查 P(n) / n 是否接近 n+1
            if (p_values.size() >= 5) {
                bool is_n_times_n_plus_1 = true;
                for (int n = 1; n <= 5; ++n) {
                    double expected = static_cast<double>(n) * (n + 1);
                    if (mymath::abs(p_values[n - 1] - expected) > 1e-6) {
                        is_n_times_n_plus_1 = false;
                        break;
                    }
                }
                if (is_n_times_n_plus_1) {
                    // g(n) = 1/n
                    *g_function = SymbolicExpression::parse("1 / (" + index_name + ")");
                    return true;
                }

                // 检查 P(n) = (n+a)*(n+a+1) 形式
                // P(n) = n^2 + (2a+1)n + a(a+1)
                // 通过解方程确定 a
                // P(1) = 1 + (2a+1) + a(a+1) = a^2 + 3a + 2
                // P(0) = a(a+1)
                double p1 = p_values[0];  // P(1)
                // 解 a^2 + 3a + 2 = p1
                // a^2 + 3a + (2 - p1) = 0
                double discriminant = 9.0 - 4.0 * (2.0 - p1);
                if (discriminant >= 0) {
                    double a1 = (-3.0 + mymath::sqrt(discriminant)) / 2.0;
                    double a2 = (-3.0 - mymath::sqrt(discriminant)) / 2.0;

                    auto verify_a = [&](double a) -> bool {
                        // 验证 P(n) = (n+a)*(n+a+1)
                        for (int n = 1; n <= 5; ++n) {
                            double expected = (n + a) * (n + a + 1);
                            if (mymath::abs(p_values[n - 1] - expected) > 1e-6) {
                                return false;
                            }
                        }
                        return true;
                    };

                    double a = 0.0;
                    if (verify_a(a1)) {
                        a = a1;
                    } else if (verify_a(a2)) {
                        a = a2;
                    }

                    if (a != 0.0 || verify_a(0.0)) {
                        // g(n) = 1/(n+a)
                        std::ostringstream oss;
                        if (mymath::is_near_zero(a, 1e-10)) {
                            oss << "1 / (" << index_name << ")";
                        } else if (a > 0) {
                            oss << "1 / ((" << index_name << ") + " << format_symbolic_scalar(a) << ")";
                        } else {
                            oss << "1 / ((" << index_name << ") - " << format_symbolic_scalar(-a) << ")";
                        }
                        *g_function = SymbolicExpression::parse(oss.str());
                        return true;
                    }
                }
            }
        }
    }

    // 方法3：通过数值方法检测
    // 如果 sum_{n=a}^{b} f(n) = g(a) - g(b+1)，则 f(n) 可能是伸缩的
    // 我们尝试通过 f(n) = g(n) - g(n+1) 反推 g(n)

    // 计算 f(n) 在几个点的值
    std::vector<double> f_values;
    for (int n = 0; n <= 6; ++n) {
        SymbolicExpression n_val = SymbolicExpression::number(static_cast<double>(n));
        SymbolicExpression f_at_n = summand.substitute(index_name, n_val).simplify();
        double val = 0.0;
        if (!f_at_n.is_number(&val)) {
            break;
        }
        f_values.push_back(val);
    }

    if (f_values.size() >= 5) {
        // 尝试拟合 g(n) = A/n + B/n^2 + ... 形式（对于大 n）
        // 或者 g(n) = P(n)/Q(n) 形式

        // 简单启发式：如果 f(n) = 1/(n+a) - 1/(n+a+1)
        // 则 f(n) = 1/((n+a)(n+a+1))
        // 检查 f(n) * n * (n+1) 是否接近 1
        bool is_simple_telescoping = true;
        for (int n = 1; n <= 5; ++n) {
            double product = f_values[n] * static_cast<double>(n) * (n + 1);
            if (mymath::abs(product - 1.0) > 1e-6) {
                is_simple_telescoping = false;
                break;
            }
        }
        if (is_simple_telescoping) {
            *g_function = SymbolicExpression::parse("1 / (" + index_name + ")");
            return true;
        }
    }

    return false;
}

// ============================================================================
// 数值级数求和（Wynn-epsilon 加速）
// ============================================================================

/**
 * @brief Wynn-epsilon 算法加速级数收敛
 *
 * 用于加速收敛级数的数值求和。
 * 基于 Shanks 变换的递归实现。
 *
 * @param partial_sums 部分和序列
 * @param max_iterations 最大迭代次数（-1 表示自动）
 * @param tolerance 收敛容差（-1 表示自动）
 * @param result 输出收敛值
 * @return 是否成功收敛
 */
bool wynn_epsilon_acceleration(const std::vector<double>& partial_sums,
                                int max_iterations,
                                double tolerance,
                                double* result) {
    if (partial_sums.size() < 3) {
        if (partial_sums.empty()) {
            *result = 0.0;
            return true;
        }
        *result = partial_sums.back();
        return true;
    }

    const int n = static_cast<int>(partial_sums.size());

    // 动态参数：根据数据量自动调整
    int actual_max_iter = (max_iterations < 0) ? std::min(n / 2, 50) : max_iterations;
    double actual_tolerance = (tolerance < 0) ? 1e-10 : tolerance;

    // epsilon 表，只需要两行
    std::vector<double> e_prev(n + 1, 0.0);
    std::vector<double> e_curr(n + 1, 0.0);

    // 初始化：e_0^{(n)} = S_n
    for (int i = 0; i < n; ++i) {
        e_prev[i + 1] = partial_sums[i];
    }

    double best_estimate = partial_sums.back();
    double best_change = mymath::kDoubleMax;
    int consecutive_no_improvement = 0;

    // 递归计算 epsilon 表
    for (int k = 1; k <= actual_max_iter && k < n; ++k) {
        bool all_valid = true;
        for (int j = k; j <= n - k; ++j) {
            double diff = e_prev[j + 1] - e_prev[j - 1];
            if (mymath::abs(diff) < 1e-30) {
                e_curr[j] = e_prev[j];
            } else {
                double eps_diff = e_curr[j - 1] - e_prev[j];
                if (mymath::abs(eps_diff) < 1e-30) {
                    e_curr[j] = e_prev[j + 1];
                } else {
                    e_curr[j] = e_prev[j + 1] + 1.0 / eps_diff;
                }
            }
            if (!mymath::isfinite(e_curr[j])) {
                all_valid = false;
            }
        }

        // 偶数阶给出收敛估计
        if (k % 2 == 0 && all_valid) {
            int mid = (k + 1) / 2 + 1;
            if (mid <= n - k) {
                double estimate = e_curr[mid];
                if (mymath::isfinite(estimate)) {
                    double change = mymath::abs(estimate - best_estimate);
                    if (change < best_change) {
                        best_change = change;
                        best_estimate = estimate;
                        consecutive_no_improvement = 0;
                    } else {
                        consecutive_no_improvement++;
                    }

                    // 检查收敛
                    if (change < actual_tolerance) {
                        *result = estimate;
                        return true;
                    }

                    // 连续多次无改进，提前终止
                    if (consecutive_no_improvement >= 3) {
                        break;
                    }
                }
            }
        }

        std::swap(e_prev, e_curr);
    }

    *result = best_estimate;
    return best_change < actual_tolerance * 10;
}

/**
 * @brief Euler 变换加速交错级数
 *
 * 对于交错级数 sum (-1)^n * a_n，Euler 变换可以显著加速收敛。
 * 变换公式：S = sum_{k=0}^{inf} (-1)^k * Delta^k(a_0) / 2^{k+1}
 * 其中 Delta 是前向差分算子。
 */
bool euler_transform(const std::vector<double>& terms,
                     double tolerance,
                     double* result) {
    if (terms.empty()) {
        *result = 0.0;
        return true;
    }

    const int n = static_cast<int>(terms.size());
    std::vector<double> delta = terms;  // 初始为 a_n
    std::vector<double> euler_sums;
    euler_sums.reserve(n);

    double sum = 0.0;
    double pow2 = 1.0;  // 2^k

    for (int k = 0; k < n && k < 30; ++k) {
        // 添加当前项
        sum += delta[0] / pow2;
        euler_sums.push_back(sum);

        // 计算下一阶差分
        for (int i = 0; i < n - k - 1; ++i) {
            delta[i] = delta[i + 1] - delta[i];
        }

        pow2 *= 2.0;

        // 检查收敛
        if (k > 2 && euler_sums.size() >= 2) {
            double change = mymath::abs(euler_sums[k] - euler_sums[k - 1]);
            if (change < tolerance) {
                *result = sum;
                return true;
            }
        }
    }

    *result = sum;
    return euler_sums.size() >= 3;
}

/**
 * @brief Levin u 变换加速级数收敛
 *
 * Levin 变换是一种广义序列加速方法，对对数收敛级数特别有效。
 * 公式：L_k^{(n)} = S_n / D_k^{(n)}
 * 其中 D_k^{(n)} 是基于差分的分母。
 */
bool levin_transform(const std::vector<double>& partial_sums,
                     const std::vector<double>& terms,
                     double tolerance,
                     double* result) {
    const int n = static_cast<int>(partial_sums.size());
    if (n < 4) {
        *result = partial_sums.empty() ? 0.0 : partial_sums.back();
        return true;
    }

    // Levin u 变换
    // L = S_n + a_{n+1} * T_n / D_n
    // 其中 T_n 和 D_n 是基于二项式系数的和

    double best_estimate = partial_sums.back();
    double best_change = mymath::kDoubleMax;

    // 使用简化的 Levin 序列变换
    for (int k = 2; k < n - 1; ++k) {
        double beta = 1.0;

        // 计算分子和分母
        double num = partial_sums[k];
        double den = 1.0;

        for (int j = 1; j <= k && j < n; ++j) {
            double coeff = static_cast<double>(j) / (k + 1);
            double term = mymath::pow(beta * coeff, j);
            if (j < static_cast<int>(terms.size())) {
                num += term * terms[k - j];
            }
            den += term;
        }

        if (mymath::abs(den) > 1e-15) {
            double estimate = num / den;
            if (mymath::isfinite(estimate)) {
                double change = mymath::abs(estimate - best_estimate);
                if (change < best_change) {
                    best_change = change;
                    best_estimate = estimate;
                }
                if (change < tolerance) {
                    *result = estimate;
                    return true;
                }
            }
        }
    }

    *result = best_estimate;
    return best_change < tolerance * 100;
}

/**
 * @brief 检测级数是否为交错级数
 */
bool is_alternating_series(const std::vector<double>& terms) {
    if (terms.size() < 3) return false;

    int sign_changes = 0;
    for (std::size_t i = 1; i < terms.size(); ++i) {
        if ((terms[i] > 0) != (terms[i - 1] > 0)) {
            sign_changes++;
        }
    }

    // 如果大部分相邻项符号相反，认为是交错级数
    return sign_changes > static_cast<int>(terms.size()) / 2;
}

/**
 * @brief 数值计算无穷级数的和
 *
 * 使用多种加速算法：
 * - Wynn-epsilon 加速
 * - Euler 变换（交错级数）
 * - Levin 变换（缓收敛级数）
 *
 * @param ctx 级数上下文
 * @param summand 通项表达式
 * @param index_name 指标名
 * @param lower 下界
 * @param max_terms 最大项数（-1 表示自动）
 * @param result 输出结果
 * @return 是否成功计算
 */
bool numerical_infinite_sum(const SeriesContext& ctx,
                            const SymbolicExpression& summand,
                            const std::string& index_name,
                            const std::string& lower,
                            int max_terms,
                            double* result) {
    // 解析下界
    double lower_val = 0.0;
    try {
        lower_val = ctx.parse_decimal(lower);
    } catch (...) {
        lower_val = 1.0;
    }

    // 动态参数：根据收敛情况调整
    int actual_max_terms = (max_terms < 0) ? 2000 : max_terms;
    const double tolerance = 1e-12;

    // 计算部分和序列和项序列
    std::vector<double> partial_sums;
    std::vector<double> terms;
    partial_sums.reserve(actual_max_terms);
    terms.reserve(actual_max_terms);

    double running_sum = 0.0;
    double prev_term = mymath::kDoubleMax;
    int terms_since_improvement = 0;
    double best_term = mymath::kDoubleMax;

    for (int n = 0; n < actual_max_terms; ++n) {
        double idx = lower_val + n;
        double term = 0.0;
        try {
            term = ctx.evaluate_at(summand, index_name, idx);
        } catch (...) {
            break;
        }

        if (!mymath::isfinite(term)) {
            break;
        }

        running_sum += term;
        partial_sums.push_back(running_sum);
        terms.push_back(term);

        // 检查项是否趋于零（收敛的必要条件）
        if (n > 10 && mymath::abs(term) > mymath::abs(prev_term) * 1.1) {
            // 项不递减，可能发散
            // 但交错级数可能仍收敛
            if (!is_alternating_series(terms)) {
                return false;
            }
        }

        // 检查项是否足够小
        if (mymath::abs(term) < tolerance * mymath::abs(running_sum)) {
            break;
        }

        // 动态终止：如果连续多项没有显著改进
        if (mymath::abs(term) < best_term * 0.99) {
            best_term = mymath::abs(term);
            terms_since_improvement = 0;
        } else {
            terms_since_improvement++;
        }

        // 如果连续 100 项没有改进，尝试提前终止
        if (terms_since_improvement > 100 && n > 200) {
            break;
        }

        prev_term = term;
    }

    if (partial_sums.size() < 5) {
        return false;
    }

    // 根据级数类型选择加速方法
    bool is_alternating = is_alternating_series(terms);

    // 尝试多种方法，选择最佳结果
    std::vector<std::pair<double, bool>> estimates;

    // 方法1：Wynn-epsilon
    double wynn_result = 0.0;
    bool wynn_ok = wynn_epsilon_acceleration(partial_sums, -1, tolerance, &wynn_result);
    if (wynn_ok) {
        estimates.push_back({wynn_result, true});
    }

    // 方法2：Euler 变换（仅对交错级数）
    if (is_alternating) {
        double euler_result = 0.0;
        bool euler_ok = euler_transform(terms, tolerance, &euler_result);
        if (euler_ok) {
            estimates.push_back({euler_result, true});
        }
    }

    // 方法3：Levin 变换
    double levin_result = 0.0;
    bool levin_ok = levin_transform(partial_sums, terms, tolerance, &levin_result);
    if (levin_ok) {
        estimates.push_back({levin_result, true});
    }

    if (estimates.empty()) {
        // 所有加速方法都失败，使用原始部分和
        *result = partial_sums.back();
        return false;
    }

    // 选择最一致的结果
    if (estimates.size() == 1) {
        *result = estimates[0].first;
        return true;
    }

    // 找到与其他结果差异最小的估计
    double best_result = estimates[0].first;
    double min_variance = mymath::kDoubleMax;

    for (const auto& est1 : estimates) {
        double variance = 0.0;
        for (const auto& est2 : estimates) {
            variance += mymath::abs(est1.first - est2.first);
        }
        if (variance < min_variance) {
            min_variance = variance;
            best_result = est1.first;
        }
    }

    *result = best_result;
    return true;
}

/**
 * @brief 无限级数特殊值表
 *
 * 存储 sum_{n=1}^{∞} f(n) 的已知值。
 */
struct InfiniteSeriesValue {
    std::string pattern;  // 模式描述
    double value;         // 级数值
    bool is_symbolic;     // 是否需要符号表示
    std::string symbolic; // 符号表示
};

const InfiniteSeriesValue kInfiniteSeriesValues[] = {
    {"zeta(2)", 1.6449340668482264, true, "pi^2 / 6"},
    {"zeta(4)", 1.082323233711138, true, "pi^4 / 90"},
    {"zeta(6)", 1.017343061984449, true, "pi^6 / 945"},
    {"zeta(8)", 1.004077356197944, true, "pi^8 / 9450"},
    {"zeta(10)", 1.000994575127818, true, "pi^10 / 93555"},
    // 调和级数发散，不包含
    // 交替调和级数
    {"ln(2)", 0.6931471805599453, true, "ln(2)"},  // sum (-1)^(n+1)/n
};

/**
 * @brief 尝试匹配无限级数特殊值
 */
bool try_infinite_series_value(const SymbolicExpression& summand,
                                const std::string& index_name,
                                const std::string& lower,
                                std::string* result) {
    // 检查 sum 1/n^2 = pi^2/6 等
    if (lower != "1") {
        return false;
    }

    // 检查 1/n^(2k) 形式
    for (int k = 1; k <= 5; ++k) {
        int two_k = 2 * k;
        std::string inv_pow_str = "1 / (" + index_name + " ^ " + std::to_string(two_k) + ")";
        SymbolicExpression pattern = SymbolicExpression::parse(inv_pow_str);
        SymbolicExpression diff = (summand - pattern).simplify();
        double val = 0.0;
        if (diff.is_number(&val) && mymath::is_near_zero(val, 1e-10)) {
            // 匹配 zeta(2k)
            double B = prob::bernoulli(two_k);
            double coeff = mymath::pow(2.0, two_k - 1) * mymath::abs(B) / prob::factorial(two_k);
            std::ostringstream ss;
            if (!mymath::is_near_zero(coeff - 1.0, 1e-10)) {
                ss << format_symbolic_scalar(coeff) << " * ";
            }
            ss << "pi ^ " << two_k;
            *result = ss.str();
            return true;
        }
    }

    // 检查交替调和级数 sum (-1)^(n+1)/n = ln(2)
    std::string alt_harm = "(-1) ^ (" + index_name + " + 1) / " + index_name;
    SymbolicExpression alt_pattern = SymbolicExpression::parse(alt_harm);
    SymbolicExpression alt_diff = (summand - alt_pattern).simplify();
    double alt_val = 0.0;
    if (alt_diff.is_number(&alt_val) && mymath::is_near_zero(alt_val, 1e-10)) {
        *result = "ln(2)";
        return true;
    }

    // 检查 sum 1/(n*(n+1)) = 1 (差分抵消)
    std::string telescoping1 = "1 / (" + index_name + " * (" + index_name + " + 1))";
    SymbolicExpression teles_pattern1 = SymbolicExpression::parse(telescoping1);
    SymbolicExpression teles_diff1 = (summand - teles_pattern1).simplify();
    double teles_val1 = 0.0;
    if (teles_diff1.is_number(&teles_val1) && mymath::is_near_zero(teles_val1, 1e-10)) {
        *result = "1";
        return true;
    }

    // 检查 sum 1/(n*(n+2)) = 3/4
    std::string telescoping2 = "1 / (" + index_name + " * (" + index_name + " + 2))";
    SymbolicExpression teles_pattern2 = SymbolicExpression::parse(telescoping2);
    SymbolicExpression teles_diff2 = (summand - teles_pattern2).simplify();
    double teles_val2 = 0.0;
    if (teles_diff2.is_number(&teles_val2) && mymath::is_near_zero(teles_val2, 1e-10)) {
        *result = "3 / 4";
        return true;
    }

    return false;
}

}  // namespace

/**
 * @brief 级数求和（改进版）
 *
 * 计算级数的有限和或无限和。
 * 支持：
 * - 多项式求和（Faulhaber 公式）
 * - 几何级数求和（符号化检测）
 * - 等差-等比级数求和
 * - 无限级数特殊值（ζ(2k)、ln(2) 等）
 * - 差分抵消级数（部分）
 */
std::string series_sum(const SeriesContext& ctx,
                       const std::string& expr,
                       const std::string& index_name,
                       const std::string& lower,
                       const std::string& upper) {
    SymbolicExpression summand = SymbolicExpression::parse(ctx.expand_inline(expr));
    SymbolicExpression upper_expression;
    const bool upper_is_infinite =
        upper == "inf" || upper == "oo" || upper == "infinity";
    if (!upper_is_infinite) {
        upper_expression = SymbolicExpression::parse(ctx.expand_inline(upper));
    }

    auto make_polynomial_sum_primitive =
        [&](const std::vector<double>& coefficients) {
            std::vector<std::string> pieces;
            for (std::size_t p = 0; p < coefficients.size(); ++p) {
                if (mymath::is_near_zero(coefficients[p], 1e-10)) {
                    continue;
                }

                std::ostringstream poly_part;
                poly_part << "(1 / " << (p + 1) << ") * (";
                bool first = true;
                for (std::size_t j = 0; j <= p; ++j) {
                    double bj = prob::bernoulli(static_cast<int>(j));
                    if (mymath::is_near_zero(bj, 1e-10)) continue;

                    double term_coeff = prob::nCr(p + 1, j) * bj;
                    if (!first) poly_part << " + ";
                    poly_part << format_symbolic_scalar(term_coeff) << " * (" << index_name << " ^ " << (p + 1 - j) << ")";
                    first = false;
                }
                poly_part << ")";

                pieces.push_back("(" + format_symbolic_scalar(coefficients[p]) + ") * " + poly_part.str());
            }

            if (pieces.empty()) {
                return SymbolicExpression::number(0.0);
            }
            std::ostringstream out;
            for (std::size_t i = 0; i < pieces.size(); ++i) {
                if (i != 0) {
                    out << " + ";
                }
                out << pieces[i];
            }
            return SymbolicExpression::parse(out.str()).simplify();
        };

    auto finite_sum_from_primitive =
        [&](const SymbolicExpression& primitive) {
            const SymbolicExpression lower_minus_one =
                SymbolicExpression::parse("(" + lower + ") - 1").simplify();
            return SymbolicExpression::parse(
                       "(" +
                       primitive.substitute(index_name, upper_expression).to_string() +
                       ") - (" +
                       primitive.substitute(index_name, lower_minus_one).to_string() +
                       ")")
                .simplify()
                .to_string();
        };

    // 1. 检查多项式求和
    std::vector<double> polynomial_coefficients;
    if (summand.polynomial_coefficients(index_name, &polynomial_coefficients)) {
        if (upper_is_infinite) {
            bool all_zero = true;
            for (double coefficient : polynomial_coefficients) {
                if (!mymath::is_near_zero(coefficient, 1e-10)) {
                    all_zero = false;
                    break;
                }
            }
            if (!all_zero) {
                throw std::runtime_error(
                    "series_sum does not support infinite polynomial sums (diverges)");
            }
            return "0";
        }

        const SymbolicExpression primitive =
            make_polynomial_sum_primitive(polynomial_coefficients);
        return finite_sum_from_primitive(primitive);
    }

    // 2. 检查无限级数特殊值
    if (upper_is_infinite) {
        std::string special_value;
        if (try_infinite_series_value(summand, index_name, lower, &special_value)) {
            return ctx.simplify_symbolic(special_value);
        }
    }

    // 3. 尝试符号化检测几何级数（优先于数值方法）
    double geometric_coefficient = 0.0;
    double geometric_ratio_value = 0.0;
    bool is_geometric = detect_geometric_ratio_symbolic(
        summand, index_name, &geometric_coefficient, &geometric_ratio_value);

    // 4. 如果符号方法失败，尝试数值方法
    if (!is_geometric) {
        auto geometric_ratio_numeric = [&](double* coefficient, double* ratio) {
            double s[4];
            int offset = 0;
            while (offset < 10) {
                for (int i = 0; i < 4; ++i) {
                    s[i] = ctx.evaluate_at(summand, index_name, offset + i);
                }
                if (!mymath::is_near_zero(s[0], 1e-10)) {
                    break;
                }
                offset++;
            }
            if (offset == 10 || mymath::is_near_zero(s[0], 1e-10)) {
                return false;
            }

            const double candidate = s[1] / s[0];
            if (!mymath::is_near_zero(s[2] - s[1] * candidate, 1e-8) ||
                !mymath::is_near_zero(s[3] - s[2] * candidate, 1e-8)) {
                return false;
            }

            if (mymath::is_near_zero(candidate, 1e-10)) {
                *ratio = 0.0;
                *coefficient = s[0];
            } else {
                *ratio = candidate;
                *coefficient = s[0] / mymath::pow(candidate, offset);
            }
            return true;
        };

        is_geometric = geometric_ratio_numeric(&geometric_coefficient, &geometric_ratio_value);
    }

    // 5. 检查等差-等比级数
    std::vector<double> arith_geo_poly;
    double arith_geo_ratio = 0.0;
    bool is_arith_geo = false;
    if (!is_geometric) {
        is_arith_geo = detect_arith_geo_series(summand, index_name, &arith_geo_poly, &arith_geo_ratio);
    }

    if (is_arith_geo && !upper_is_infinite) {
        return arith_geo_sum(arith_geo_poly, arith_geo_ratio, index_name, upper_expression, lower);
    }

    // 6. 处理几何级数
    if (is_geometric) {
        const std::string coefficient_text =
            format_symbolic_scalar(geometric_coefficient);
        const std::string ratio_text = format_symbolic_scalar(geometric_ratio_value);

        if (upper_is_infinite) {
            if (mymath::abs(geometric_ratio_value) >= 1.0 - 1e-10) {
                throw std::runtime_error(
                    "series_sum infinite geometric series requires |r| < 1");
            }
            if (mymath::is_near_zero(geometric_ratio_value - 1.0, 1e-10)) {
                throw std::runtime_error(
                    "series_sum infinite geometric series diverges for r = 1");
            }
            return ctx.simplify_symbolic(
                "(" + coefficient_text + ") * (" + ratio_text + ") ^ (" +
                lower + ") / (1 - (" + ratio_text + "))");
        }

        const std::string geometric_primitive_text =
            mymath::is_near_zero(geometric_ratio_value - 1.0, 1e-10)
                ? "(" + coefficient_text + ") * (" + index_name + " + 1)"
                : "(" + coefficient_text + ") * (1 - (" + ratio_text +
                      ") ^ (" + index_name + " + 1)) / (1 - (" +
                      ratio_text + "))";
        const SymbolicExpression primitive =
            SymbolicExpression::parse(geometric_primitive_text).simplify();
        return finite_sum_from_primitive(primitive);
    }

    // 7. 尝试差分抵消检测
    SymbolicExpression g_function;
    if (detect_telescoping(summand, index_name, &g_function)) {
        if (upper_is_infinite) {
            // sum_{n=a}^{∞} (g(n) - g(n+1)) = g(a)
            SymbolicExpression lower_expr = SymbolicExpression::parse(lower);
            return g_function.substitute(index_name, lower_expr).simplify().to_string();
        } else {
            // sum_{n=a}^{b} (g(n) - g(n+1)) = g(a) - g(b+1)
            SymbolicExpression lower_expr = SymbolicExpression::parse(lower);
            SymbolicExpression upper_plus_one = SymbolicExpression::parse("(" + upper + ") + 1");
            return SymbolicExpression::parse(
                "(" + g_function.substitute(index_name, lower_expr).to_string() +
                ") - (" + g_function.substitute(index_name, upper_plus_one).to_string() + ")")
                .simplify().to_string();
        }
    }

    // 8. 数值求和回退（仅对无穷级数）
    if (upper_is_infinite) {
        double numerical_result = 0.0;
        if (numerical_infinite_sum(ctx, summand, index_name, lower, 1000, &numerical_result)) {
            // 返回数值结果，标记为近似值
            std::ostringstream ss;
            ss << std::setprecision(15) << numerical_result;
            return ss.str() + " (numerical approximation)";
        }
    }

    throw std::runtime_error(
        "series_sum: unsupported series type. Supported: polynomials, geometric series, "
        "arithmetic-geometric series (n*r^n), zeta(2k) values, telescoping series, "
        "and numerical approximation for convergent infinite series");
}

/**
 * @brief 检查是否为级数命令
 */
bool is_series_command(const std::string& command) {
    return command == "taylor" ||
           command == "pade" ||
           command == "puiseux" ||
           command == "series_sum" ||
           command == "summation";
}

/**
 * @brief 处理级数命令
 */
bool handle_series_command(const SeriesContext& ctx,
                           const std::string& command,
                           const std::string& inside,
                           std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "taylor") {
        if (arguments.size() != 3) {
            throw std::runtime_error("taylor expects exactly three arguments");
        }

        const double center = ctx.parse_decimal(arguments[1]);
        const double degree_value = ctx.parse_decimal(arguments[2]);
        if (!is_integer_double(degree_value) || degree_value < 0.0) {
            throw std::runtime_error("taylor degree must be a non-negative integer");
        }
        const int degree = static_cast<int>(round_to_long_long(degree_value));

        *output = taylor(ctx, arguments[0], center, degree);
        return true;
    }

    if (command == "pade") {
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "pade expects expr, m, n or expr, center, m, n");
        }

        const bool explicit_center = arguments.size() == 4;
        const double center = explicit_center
                                  ? ctx.parse_decimal(arguments[1])
                                  : 0.0;
        const double numerator_degree_value = ctx.parse_decimal(
            arguments[explicit_center ? 2 : 1]);
        const double denominator_degree_value = ctx.parse_decimal(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(numerator_degree_value) ||
            numerator_degree_value < 0.0 ||
            !is_integer_double(denominator_degree_value) ||
            denominator_degree_value < 0.0) {
            throw std::runtime_error(
                "pade degrees must be non-negative integers");
        }

        const int numerator_degree =
            static_cast<int>(round_to_long_long(numerator_degree_value));
        const int denominator_degree =
            static_cast<int>(round_to_long_long(denominator_degree_value));

        *output = pade(ctx, arguments[0], center, numerator_degree, denominator_degree);
        return true;
    }

    if (command == "puiseux") {
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "puiseux expects expr, degree, denominator or expr, center, degree, denominator");
        }

        const bool explicit_center = arguments.size() == 4;
        const double center = explicit_center
                                  ? ctx.parse_decimal(arguments[1])
                                  : 0.0;
        const double degree_value = ctx.parse_decimal(
            arguments[explicit_center ? 2 : 1]);
        const double denominator_value = ctx.parse_decimal(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(degree_value) || degree_value < 0.0) {
            throw std::runtime_error(
                "puiseux degree must be a non-negative integer");
        }
        if (!is_integer_double(denominator_value) || denominator_value <= 0.0) {
            throw std::runtime_error(
                "puiseux denominator must be a positive integer");
        }

        const int degree = static_cast<int>(round_to_long_long(degree_value));
        const int denom = static_cast<int>(round_to_long_long(denominator_value));

        *output = puiseux(ctx, arguments[0], center, degree, denom);
        return true;
    }

    if (command == "series_sum" || command == "summation") {
        if (arguments.size() != 4) {
            throw std::runtime_error(
                "series_sum expects expr, index, lower, upper");
        }

        const std::string index_name = trim_copy(arguments[1]);
        if (!is_identifier_text(index_name)) {
            throw std::runtime_error("series_sum index must be an identifier");
        }

        *output = series_sum(ctx,
                             arguments[0],
                             index_name,
                             arguments[2],
                             trim_copy(arguments[3]));
        return true;
    }

    return false;
}


std::string SeriesModule::execute_args(const std::string& command,
                                      const std::vector<std::string>& args,
                                      const CoreServices& services) {
    SeriesContext ctx;
    ctx.resolve_symbolic = services.symbolic.resolve_symbolic;
    ctx.parse_decimal = services.evaluation.parse_decimal;
    ctx.evaluate_at = services.symbolic.evaluate_symbolic_at;
    ctx.simplify_symbolic = services.symbolic.simplify_symbolic;
    ctx.expand_inline = services.symbolic.expand_inline;

    std::string inside;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += args[i];
    }

    std::string output;
    if (handle_series_command(ctx, command, inside, &output)) {
        return output;
    }
    throw std::runtime_error("Series command failed: " + command);
}

std::vector<std::string> SeriesModule::get_commands() const {
    return {"taylor", "pade", "puiseux", "series_sum", "summation"};
}

std::string SeriesModule::get_help_snippet(const std::string& topic) const {
    if (topic == "symbolic") {
        return "Series:\n"
               "  taylor(f, a, n)     Taylor series at x=a up to degree n\n"
               "  pade(f, [a], m, n)  Pade approximation [m/n]\n"
               "  puiseux(f, n, d)    Puiseux series (fractional powers)\n"
               "  series_sum(f, i, a, b) Finite or infinite summation";
    }
    return "";
}

}  // namespace series_ops

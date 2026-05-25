/**
 * @file psa_engine.cpp
 * @brief 幂级数分析（PSA）引擎实现
 *
 * 本文件实现了幂级数分析的核心算法：
 * - Laurent 级数：支持负幂项的级数表示
 * - 级数运算：加、减、乘、除、复合
 * - 极点检测：识别级数的极点位置和阶数
 * - 奇点分析：分析函数的奇点类型
 *
 * PSA 引擎是级数展开和奇点分析的基础设施。
 */

#include "analysis/series/psa_engine.h"
#include "analysis/modules/series_module.h"
#include "math/base/precision_constants.h"
#include "symbolic/core/symbolic_expression.h"
#include "math/mymath.h"
#include "math/helpers/integer_helpers.h"
#include <algorithm>

namespace series_ops {
namespace internal {

// 获取 PSA 引擎的数值容差
inline Scalar psa_tolerance() {
    return precision::epsilon<Scalar>();
}

// 获取 PSA 引擎的系数截断容差（稍宽松）
inline Scalar psa_coefficient_threshold() {
    return precision::sqrt_epsilon<Scalar>() * Scalar(10);
}

PoleException::PoleException(int s, Scalar coeff)
    : std::runtime_error("Pole encountered"), shift(s), leading_coefficient(coeff) {}

LaurentSeries::LaurentSeries() : low_power(0), high_power(-1) {}

LaurentSeries::LaurentSeries(int low, int high, const std::vector<Scalar>& coeffs)
    : coefficients(coeffs), low_power(low), high_power(high) {}

bool LaurentSeries::empty() const { return coefficients.empty() || low_power > high_power; }

Scalar LaurentSeries::get_coefficient(int power) const {
    if (power < low_power || power > high_power) return Scalar(0);
    return coefficients[power - low_power];
}

bool LaurentSeries::to_taylor(std::vector<Scalar>& result, int degree) const {
    if (low_power < 0) return false;
    result.assign(degree + 1, Scalar(0));
    for (int k = low_power; k <= high_power && k <= degree; ++k) {
        result[k] = get_coefficient(k);
    }
    return true;
}

std::vector<Scalar> ps_add(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree) {
    std::vector<Scalar> res(degree + 1, Scalar(0));
    for (int i = 0; i <= degree; ++i) {
        res[i] = (i < static_cast<int>(a.size()) ? a[i] : Scalar(0)) +
                 (i < static_cast<int>(b.size()) ? b[i] : Scalar(0));
    }
    return res;
}

std::vector<Scalar> ps_sub(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree) {
    std::vector<Scalar> res(degree + 1, Scalar(0));
    for (int i = 0; i <= degree; ++i) {
        res[i] = (i < static_cast<int>(a.size()) ? a[i] : Scalar(0)) -
                 (i < static_cast<int>(b.size()) ? b[i] : Scalar(0));
    }
    return res;
}

std::vector<Scalar> ps_mul(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree) {
    std::vector<Scalar> res(degree + 1, Scalar(0));
    for (int i = 0; i <= degree; ++i) {
        if (i >= static_cast<int>(a.size()) || mymath::abs(a[i]) < psa_tolerance()) continue;
        for (int j = 0; i + j <= degree; ++j) {
            if (j >= static_cast<int>(b.size()) || mymath::abs(b[j]) < psa_tolerance()) continue;
            res[i + j] = res[i + j] + a[i] * b[j];
        }
    }
    return res;
}

std::vector<Scalar> ps_div_with_laurent(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree, int* laurent_shift) {
    if (b.empty()) throw std::runtime_error("division by empty power series");

    int start_a = -1;
    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
        if (mymath::abs(a[i]) >= psa_tolerance() * Scalar(1e-3)) {
            start_a = i;
            break;
        }
    }

    int start_b = -1;
    for (int i = 0; i < static_cast<int>(b.size()); ++i) {
        if (mymath::abs(b[i]) >= psa_tolerance() * Scalar(1e-3)) {
            start_b = i;
            break;
        }
    }

    if (start_b == -1) throw std::runtime_error("division by zero in power series");

    if (start_a == -1) {
        if (laurent_shift) *laurent_shift = 0;
        return std::vector<Scalar>(degree + 1, Scalar(0));
    }

    int shift = start_a - start_b;
    if (laurent_shift) *laurent_shift = shift;

    std::vector<Scalar> a_effective;
    for (int i = start_a; i < static_cast<int>(a.size()); ++i) a_effective.push_back(a[i]);

    std::vector<Scalar> b_effective;
    for (int i = start_b; i < static_cast<int>(b.size()); ++i) b_effective.push_back(b[i]);

    int result_start = (shift < 0) ? shift : 0;
    int result_size = degree + 1 - result_start;
    if (result_size <= 0) {
        return std::vector<Scalar>(degree + 1, Scalar(0));
    }

    std::vector<Scalar> res_effective(result_size, Scalar(0));
    Scalar inv_b0 = Scalar(1) / b_effective[0];
    for (int i = 0; i < result_size; ++i) {
        Scalar val = (i < static_cast<int>(a_effective.size()) ? a_effective[i] : Scalar(0));
        for (int j = 1; j <= i; ++j) {
            if (j < static_cast<int>(b_effective.size()))
                val = val - b_effective[j] * res_effective[i - j];
        }
        res_effective[i] = val * inv_b0;
    }

    std::vector<Scalar> final_res(degree + 1, Scalar(0));
    for (int i = 0; i <= degree; ++i) {
        int effective_idx = i - result_start;
        if (effective_idx >= 0 && effective_idx < static_cast<int>(res_effective.size())) {
            final_res[i] = res_effective[effective_idx];
        }
    }
    return final_res;
}

std::vector<Scalar> ps_div(const std::vector<Scalar>& a, const std::vector<Scalar>& b, int degree) {
    int shift = 0;
    std::vector<Scalar> result = ps_div_with_laurent(a, b, degree, &shift);
    if (shift < 0) {
        Scalar leading = Scalar(0);
        for (const Scalar& v : result) {
            if (mymath::abs(v) >= psa_tolerance()) {
                leading = v;
                break;
            }
        }
        if (mymath::abs(leading) < psa_tolerance()) {
            int start_a = -1;
            for (int i = 0; i < static_cast<int>(a.size()); ++i) {
                if (mymath::abs(a[i]) >= psa_tolerance()) {
                    start_a = i;
                    break;
                }
            }
            int start_b = -1;
            for (int i = 0; i < static_cast<int>(b.size()); ++i) {
                if (mymath::abs(b[i]) >= psa_tolerance()) {
                    start_b = i;
                    break;
                }
            }
            if (start_a >= 0 && start_b >= 0) {
                leading = a[start_a] / b[start_b];
            }
        }
        throw PoleException(shift, leading.to_long_double());
    }
    return result;
}

std::vector<Scalar> ps_exp(const std::vector<Scalar>& a, int degree) {
    std::vector<Scalar> res(degree + 1, Scalar(0));
    Scalar a0 = a.empty() ? Scalar(0) : a[0];
    res[0] = mymath::exp(a0);
    for (int i = 1; i <= degree; ++i) {
        Scalar sum = Scalar(0);
        for (int k = 1; k <= i; ++k) {
            Scalar ak = k < static_cast<int>(a.size()) ? a[k] : Scalar(0);
            sum = sum + Scalar(k) * ak * res[i - k];
        }
        res[i] = sum / Scalar(i);
    }
    return res;
}

std::vector<Scalar> ps_ln(const std::vector<Scalar>& a, int degree) {
    if (a.empty() || a[0] <= Scalar(0)) throw std::runtime_error("ln of non-positive power series base");
    std::vector<Scalar> res(degree + 1, Scalar(0));
    res[0] = mymath::ln(a[0]);
    Scalar inv_a0 = Scalar(1) / a[0];
    for (int i = 1; i <= degree; ++i) {
        Scalar sum = Scalar(0);
        for (int k = 1; k < i; ++k) {
            Scalar ak = i - k < static_cast<int>(a.size()) ? a[i - k] : Scalar(0);
            sum = sum + Scalar(k) * res[k] * ak;
        }
        Scalar ai = i < static_cast<int>(a.size()) ? a[i] : Scalar(0);
        res[i] = (ai - sum / Scalar(i)) * inv_a0;
    }
    return res;
}

void ps_sincos(const std::vector<Scalar>& a, int degree, std::vector<Scalar>& sin_res, std::vector<Scalar>& cos_res) {
    sin_res.assign(degree + 1, Scalar(0));
    cos_res.assign(degree + 1, Scalar(0));
    Scalar a0 = a.empty() ? Scalar(0) : a[0];
    sin_res[0] = mymath::sin(a0);
    cos_res[0] = mymath::cos(a0);
    for (int i = 1; i <= degree; ++i) {
        Scalar sum_sin = Scalar(0);
        Scalar sum_cos = Scalar(0);
        for (int k = 1; k <= i; ++k) {
            Scalar ak = k < static_cast<int>(a.size()) ? a[k] : Scalar(0);
            sum_sin = sum_sin + Scalar(k) * ak * cos_res[i - k];
            sum_cos = sum_cos - Scalar(k) * ak * sin_res[i - k];
        }
        sin_res[i] = sum_sin / Scalar(i);
        cos_res[i] = sum_cos / Scalar(i);
    }
}

std::vector<Scalar> ps_sin(const std::vector<Scalar>& a, int degree) {
    std::vector<Scalar> s, c;
    ps_sincos(a, degree, s, c);
    return s;
}

std::vector<Scalar> ps_cos(const std::vector<Scalar>& a, int degree) {
    std::vector<Scalar> s, c;
    ps_sincos(a, degree, s, c);
    return c;
}

std::vector<Scalar> ps_pow_const(const std::vector<Scalar>& a, Scalar n, int degree) {
    if (a.empty() || mymath::abs(a[0]) < psa_tolerance()) {
        if (mymath::abs(Scalar(n)) < psa_tolerance()) {
            std::vector<Scalar> res(degree + 1, Scalar(0));
            res[0] = Scalar(1);
            return res;
        }

        int leading = -1;
        for (int i = 0; i < static_cast<int>(a.size()); ++i) {
            if (mymath::abs(a[i]) >= psa_tolerance() * Scalar(1e-3)) {
                leading = i;
                break;
            }
        }

        if (leading < 0) {
            if (n > 0) {
                return std::vector<Scalar>(degree + 1, Scalar(0));
            } else {
                throw std::runtime_error("0^negative is undefined");
            }
        }

        const Scalar shifted_power = Scalar(leading) * n;

        if (n > 0 && mymath::is_integer(n, psa_tolerance())) {
            std::vector<Scalar> res(degree + 1, Scalar(0));
            res[0] = Scalar(1);
            std::vector<Scalar> base = a;
            int p = static_cast<int>(n + 0.5);
            for (int i = 0; i < p; ++i) res = ps_mul(res, base, degree);
            return res;
        }

        if (n < 0 && mymath::is_integer(n, psa_tolerance())) {
            int p = static_cast<int>(-n + 0.5);
            std::vector<Scalar> pos_pow = ps_pow_const(a, Scalar(p), degree);
            std::vector<Scalar> one(degree + 1, Scalar(0));
            one[0] = Scalar(1);
            return ps_div_with_laurent(one, pos_pow, degree, nullptr);
        }

        if (is_integer_double(shifted_power, psa_coefficient_threshold())) {
            const int shift = static_cast<int>(round_to_long_long(shifted_power));
            if (shift < 0) {
                Scalar leading_coeff = mymath::pow(a[leading], Scalar(n));
                throw PoleException(shift, leading_coeff.to_long_double());
            }
            if (shift > degree) {
                return std::vector<Scalar>(degree + 1, Scalar(0));
            }

            std::vector<Scalar> normalized;
            normalized.reserve(a.size() - static_cast<std::size_t>(leading));
            for (int i = leading; i < static_cast<int>(a.size()); ++i) {
                normalized.push_back(a[static_cast<std::size_t>(i)]);
            }

            const std::vector<Scalar> powered = ps_pow_const(normalized, n, degree - shift);

            std::vector<Scalar> res(degree + 1, Scalar(0));
            for (int i = 0; i + shift <= degree && i < static_cast<int>(powered.size()); ++i) {
                res[static_cast<std::size_t>(i + shift)] = powered[static_cast<std::size_t>(i)];
            }
            return res;
        }

        if (shifted_power > psa_coefficient_threshold()) {
            return std::vector<Scalar>(degree + 1, Scalar(0));
        } else if (shifted_power < -psa_coefficient_threshold()) {
            Scalar leading_coeff = mymath::pow(a[leading], Scalar(n));
            throw PoleException(-1, leading_coeff.to_long_double());
        }
    }

    std::vector<Scalar> res(degree + 1, Scalar(0));
    Scalar n_scalar(n);
    res[0] = mymath::pow(a[0], n_scalar);
    Scalar inv_a0 = Scalar(1) / a[0];
    for (int i = 1; i <= degree; ++i) {
        Scalar sum = Scalar(0);
        for (int k = 1; k <= i; ++k) {
            Scalar ak = k < static_cast<int>(a.size()) ? a[k] : Scalar(0);
            sum = sum + (n_scalar * Scalar(k) - Scalar(i - k)) * ak * res[i - k];
        }
        res[i] = sum * inv_a0 / Scalar(i);
    }
    return res;
}

std::vector<Scalar> ps_scale(const std::vector<Scalar>& a, Scalar scale, int degree) {
    std::vector<Scalar> res(degree + 1, Scalar(0));
    for (int i = 0; i <= degree && i < static_cast<int>(a.size()); ++i) {
        res[i] = a[i] * scale;
    }
    return res;
}

std::vector<Scalar> ps_derivative(const std::vector<Scalar>& a, int degree) {
    std::vector<Scalar> res(degree + 1, Scalar(0));
    for (int i = 0; i < degree && i + 1 < static_cast<int>(a.size()); ++i) {
        res[i] = a[i + 1] * Scalar(i + 1);
    }
    return res;
}

std::vector<Scalar> ps_integral(const std::vector<Scalar>& a, Scalar constant_term, int degree) {
    std::vector<Scalar> res(degree + 1, Scalar(0));
    res[0] = constant_term;
    for (int i = 1; i <= degree && i < static_cast<int>(a.size()) + 1; ++i) {
        res[i] = (i - 1 < static_cast<int>(a.size()) ? a[i - 1] : Scalar(0)) / Scalar(i);
    }
    return res;
}

std::vector<Scalar> ps_tan(const std::vector<Scalar>& a, int degree) {
    std::vector<Scalar> s, c;
    ps_sincos(a, degree, s, c);
    return ps_div(s, c, degree);
}

std::vector<Scalar> ps_asin(const std::vector<Scalar>& a, int degree) {
    if (a.empty()) return std::vector<Scalar>(degree + 1, Scalar(0));
    std::vector<Scalar> res(degree + 1, Scalar(0));
    res[0] = mymath::asin(a[0]);

    std::vector<Scalar> one_minus_x2(degree + 1, Scalar(0));
    one_minus_x2[0] = Scalar(1);
    if (degree >= 2) one_minus_x2[2] = -Scalar(1);

    std::vector<Scalar> sqrt_denom = ps_pow_const(one_minus_x2, Scalar(-0.5L), degree);
    std::vector<Scalar> deriv = ps_derivative(a, degree);
    return ps_integral(ps_mul(deriv, sqrt_denom, degree), res[0], degree);
}

std::vector<Scalar> ps_acos(const std::vector<Scalar>& a, int degree) {
    if (a.empty()) return std::vector<Scalar>(degree + 1, Scalar(0));
    std::vector<Scalar> res(degree + 1, Scalar(0));
    res[0] = mymath::acos(a[0]);

    std::vector<Scalar> one_minus_x2(degree + 1, Scalar(0));
    one_minus_x2[0] = Scalar(1);
    if (degree >= 2) one_minus_x2[2] = -Scalar(1);

    std::vector<Scalar> sqrt_denom = ps_pow_const(one_minus_x2, Scalar(-0.5L), degree);
    std::vector<Scalar> deriv = ps_derivative(a, degree);
    return ps_integral(ps_mul(ps_scale(deriv, Scalar(-1), degree), sqrt_denom, degree), res[0], degree);
}

std::vector<Scalar> ps_atan(const std::vector<Scalar>& a, int degree) {
    if (a.empty()) return std::vector<Scalar>(degree + 1, Scalar(0));
    std::vector<Scalar> res(degree + 1, Scalar(0));
    res[0] = mymath::atan(a[0]);

    std::vector<Scalar> one_plus_x2(degree + 1, Scalar(0));
    one_plus_x2[0] = Scalar(1);
    if (degree >= 2) one_plus_x2[2] = Scalar(1);

    std::vector<Scalar> denom = ps_pow_const(one_plus_x2, Scalar(-1), degree);
    std::vector<Scalar> deriv = ps_derivative(a, degree);
    return ps_integral(ps_mul(deriv, denom, degree), res[0], degree);
}

std::vector<Scalar> ps_sinh(const std::vector<Scalar>& a, int degree) {
    std::vector<Scalar> res(degree + 1, Scalar(0));
    Scalar a0 = a.empty() ? Scalar(0) : a[0];
    res[0] = Scalar(mymath::sinh(a0.to_long_double()));
    for (int i = 1; i <= degree; ++i) {
        Scalar sum = Scalar(0);
        for (int k = 1; k <= i; ++k) {
            Scalar ak = k < static_cast<int>(a.size()) ? a[k] : Scalar(0);
            sum = sum + Scalar(k) * ak * res[i - k];
        }
        res[i] = sum / Scalar(i);
    }
    return res;
}

std::vector<Scalar> ps_cosh(const std::vector<Scalar>& a, int degree) {
    std::vector<Scalar> res(degree + 1, Scalar(0));
    Scalar a0 = a.empty() ? Scalar(0) : a[0];
    res[0] = Scalar(mymath::cosh(a0.to_long_double()));
    for (int i = 1; i <= degree; ++i) {
        Scalar sum = Scalar(0);
        for (int k = 1; k <= i; ++k) {
            Scalar ak = k < static_cast<int>(a.size()) ? a[k] : Scalar(0);
            sum = sum + Scalar(k) * ak * res[i - k];
        }
        res[i] = sum / Scalar(i);
    }
    return res;
}

std::vector<Scalar> ps_tanh(const std::vector<Scalar>& a, int degree) {
    std::vector<Scalar> sh = ps_sinh(a, degree);
    std::vector<Scalar> ch = ps_cosh(a, degree);
    return ps_div(sh, ch, degree);
}

bool evaluate_psa(const SymbolicExpression& expr, const std::string& var_name, Scalar center, int degree, std::vector<Scalar>& result, const SeriesContext& ctx) {
    if (!expr.has_node()) return false;

    if (expr.node_type() == NodeType::kNumber) {
        result.assign(degree + 1, Scalar(0));
        result[0] = Scalar(expr.node_numeric_value());
        return true;
    }
    if (expr.node_type() == NodeType::kPi) {
        result.assign(degree + 1, Scalar(0));
        result[0] = Scalar(mymath::kPi);
        return true;
    }
    if (expr.node_type() == NodeType::kE) {
        result.assign(degree + 1, Scalar(0));
        result[0] = Scalar(mymath::kE);
        return true;
    }
    if (expr.node_type() == NodeType::kVariable) {
        result.assign(degree + 1, Scalar(0));
        if (expr.node_text() == var_name) {
            result[0] = Scalar(center);
            if (degree >= 1) result[1] = Scalar(1);
        } else {
            result[0] = Scalar(ctx.evaluate_at(expr, var_name, center));
        }
        return true;
    }

    std::vector<Scalar> left_res, right_res;
    SymbolicExpression left_expr = expr.left_child();
    SymbolicExpression right_expr = expr.right_child();
    if (left_expr.has_node() && !evaluate_psa(left_expr, var_name, center, degree, left_res, ctx)) return false;
    if (right_expr.has_node() && !evaluate_psa(right_expr, var_name, center, degree, right_res, ctx)) return false;

    try {
        switch (expr.node_type()) {
            case NodeType::kAdd: result = ps_add(left_res, right_res, degree); return true;
            case NodeType::kSubtract: result = ps_sub(left_res, right_res, degree); return true;
            case NodeType::kMultiply: result = ps_mul(left_res, right_res, degree); return true;
            case NodeType::kDivide: result = ps_div(left_res, right_res, degree); return true;
            case NodeType::kNegate:
                result.assign(degree + 1, Scalar(0));
                for (int i = 0; i <= degree; ++i) result[i] = -left_res[i];
                return true;
            case NodeType::kPower: {
                bool right_is_const = true;
                for (int i = 1; i <= degree; ++i) {
                    if (i < static_cast<int>(right_res.size()) &&
                        mymath::abs(right_res[i]) >= psa_tolerance()) {
                        right_is_const = false; break;
                    }
                }
                if (right_is_const) {
                    Scalar p = right_res.empty() ? 0.0L : right_res[0].to_long_double();
                    result = ps_pow_const(left_res, p, degree);
                    return true;
                } else {
                    result = ps_exp(ps_mul(right_res, ps_ln(left_res, degree), degree), degree);
                    return true;
                }
            }
            case NodeType::kFunction: {
                const std::string func_name = expr.node_text();
                if (func_name == "exp") { result = ps_exp(left_res, degree); return true; }
                if (func_name == "ln") { result = ps_ln(left_res, degree); return true; }
                if (func_name == "sin") { result = ps_sin(left_res, degree); return true; }
                if (func_name == "cos") { result = ps_cos(left_res, degree); return true; }
                if (func_name == "tan") { result = ps_tan(left_res, degree); return true; }
                if (func_name == "asin" || func_name == "arcsin") { result = ps_asin(left_res, degree); return true; }
                if (func_name == "acos" || func_name == "arccos") { result = ps_acos(left_res, degree); return true; }
                if (func_name == "atan" || func_name == "arctan") { result = ps_atan(left_res, degree); return true; }
                if (func_name == "sinh") { result = ps_sinh(left_res, degree); return true; }
                if (func_name == "cosh") { result = ps_cosh(left_res, degree); return true; }
                if (func_name == "tanh") { result = ps_tanh(left_res, degree); return true; }
                if (func_name == "sqrt") { result = ps_pow_const(left_res, 0.5, degree); return true; }
                return false;
            }
            default: return false;
        }
    } catch (const PoleException&) {
        throw;
    } catch (...) {
        return false;
    }
}

}  // namespace internal
}  // namespace series_ops

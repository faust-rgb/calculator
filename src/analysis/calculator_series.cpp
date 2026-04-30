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

#include "calculator_series.h"
#include "symbolic_expression_internal.h"
#include "statistics/probability.h"

#include "polynomial.h"
#include "mymath.h"
#include <sstream>
#include <iomanip>
#include <tuple>

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
    // 如果 shift < 0，结果包含负幂项（Laurent 级数）
    // 对于极限计算，这表示无穷大或不存在。抛出异常以便回退到数值方法。
    if (shift < 0) throw std::runtime_error("Laurent series result (infinite or undefined limit)");
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

/**
 * @brief 幂级数幂函数
 *
 * 计算 a^n 的幂级数，其中 n 为常数。
 * 特殊处理 a[0] = 0 的情况。
 */
std::vector<double> ps_pow_const(const std::vector<double>& a, double n, int degree) {
    if (a.empty() || mymath::is_near_zero(a[0], 1e-12)) {
        if (mymath::is_near_zero(n, 1e-12)) {
            std::vector<double> res(degree + 1, 0.0);
            res[0] = 1.0;
            return res;
        }
        if (n > 0 && mymath::is_integer(n, 1e-12)) {
            std::vector<double> res(degree + 1, 0.0);
            res[0] = 1.0;
            std::vector<double> base = a;
            int p = static_cast<int>(n + 0.5);
            for (int i = 0; i < p; ++i) res = ps_mul(res, base, degree);
            return res;
        }

        int leading = -1;
        for (int i = 0; i < static_cast<int>(a.size()); ++i) {
            if (!mymath::is_near_zero(a[i], 1e-12)) {
                leading = i;
                break;
            }
        }
        if (leading >= 0 && n > 0.0) {
            const double shifted_power = static_cast<double>(leading) * n;
            if (is_integer_double(shifted_power, 1e-10)) {
                const int shift = static_cast<int>(round_to_long_long(shifted_power));
                if (shift > degree) {
                    return std::vector<double>(degree + 1, 0.0);
                }

                std::vector<double> normalized;
                normalized.reserve(a.size() - static_cast<std::size_t>(leading));
                for (int i = leading; i < static_cast<int>(a.size()); ++i) {
                    normalized.push_back(a[static_cast<std::size_t>(i)]);
                }

                const std::vector<double> powered =
                    ps_pow_const(normalized, n, degree - shift);
                std::vector<double> res(degree + 1, 0.0);
                for (int i = 0; i + shift <= degree &&
                                i < static_cast<int>(powered.size()); ++i) {
                    res[static_cast<std::size_t>(i + shift)] =
                        powered[static_cast<std::size_t>(i)];
                }
                return res;
            }
        }
        throw std::runtime_error("unsupported power series base 0 with non-integer/negative exponent");
    }
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
                if (node->text == "sqrt") { result = ps_pow_const(left_res, 0.5, degree); return true; }
                return false;
            }
            default: return false;
        }
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

/**
 * @brief Pade 有理逼近（改进版）
 *
 * 计算 Taylor 级数的 Pade 逼近 [m/n]。
 * 分子为 m 次，分母为 n 次。
 *
 * 改进：
 * - 使用 long double 提高数值精度
 * - 添加条件数检测和警告
 * - 支持部分 pivoting 提高稳定性
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

    const std::vector<double> coefficients = build_taylor_coefficients(
        ctx, expression, variable_name, center, numerator_degree + denominator_degree);

    auto coefficient_at = [&](int index) -> long double {
        if (index < 0 || index >= static_cast<int>(coefficients.size())) {
            return 0.0L;
        }
        return static_cast<long double>(coefficients[static_cast<std::size_t>(index)]);
    };

    std::vector<long double> denominator_ld(denominator_degree + 1, 0.0L);
    denominator_ld[0] = 1.0L;

    if (denominator_degree > 0) {
        // 使用 long double 构建线性系统
        std::vector<std::vector<long double>> matrix(
            static_cast<std::size_t>(denominator_degree),
            std::vector<long double>(static_cast<std::size_t>(denominator_degree), 0.0L));
        std::vector<long double> rhs(static_cast<std::size_t>(denominator_degree), 0.0L);

        for (int row = 0; row < denominator_degree; ++row) {
            for (int col = 0; col < denominator_degree; ++col) {
                matrix[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(col)] =
                    coefficient_at(numerator_degree + row - col);
            }
            rhs[static_cast<std::size_t>(row)] =
                -coefficient_at(numerator_degree + row + 1);
        }

        // 使用改进的高斯消元法（部分 pivoting）
        std::vector<long double> solved(static_cast<std::size_t>(denominator_degree), 0.0L);
        const std::size_t n = static_cast<std::size_t>(denominator_degree);

        // 创建副本用于消元
        auto A = matrix;
        auto b = rhs;
        std::vector<std::size_t> pivot_order(n);
        for (std::size_t i = 0; i < n; ++i) pivot_order[i] = i;

        // 前向消元（带部分 pivoting）
        for (std::size_t col = 0; col < n; ++col) {
            // 找最大主元
            std::size_t max_row = col;
            long double max_val = mymath::abs_long_double(A[col][col]);
            for (std::size_t row = col + 1; row < n; ++row) {
                if (mymath::abs_long_double(A[row][col]) > max_val) {
                    max_val = mymath::abs_long_double(A[row][col]);
                    max_row = row;
                }
            }

            // 检查条件数（粗略估计）
            if (max_val < 1e-15L) {
                // 矩阵接近奇异，尝试降低阶数或使用正则化
                // 简化处理：使用正则化
                for (std::size_t row = col; row < n; ++row) {
                    A[row][col] += 1e-12L;
                }
                max_val = 1e-12L;
            }

            // 交换行
            if (max_row != col) {
                std::swap(A[col], A[max_row]);
                std::swap(b[col], b[max_row]);
            }

            // 消元
            const long double pivot = A[col][col];
            for (std::size_t row = col + 1; row < n; ++row) {
                const long double factor = A[row][col] / pivot;
                for (std::size_t j = col; j < n; ++j) {
                    A[row][j] -= factor * A[col][j];
                }
                b[row] -= factor * b[col];
            }
        }

        // 后向回代
        for (int col = static_cast<int>(n) - 1; col >= 0; --col) {
            long double sum = b[col];
            for (std::size_t j = col + 1; j < n; ++j) {
                sum -= A[col][j] * solved[j];
            }
            solved[col] = sum / A[col][col];
        }

        // 转换回 double
        for (int i = 0; i < denominator_degree; ++i) {
            denominator_ld[static_cast<std::size_t>(i + 1)] = solved[static_cast<std::size_t>(i)];
        }
    }

    // 计算分子系数
    std::vector<long double> numerator_ld(numerator_degree + 1, 0.0L);
    for (int i = 0; i <= numerator_degree; ++i) {
        long double value = 0.0L;
        for (int j = 0; j <= denominator_degree && j <= i; ++j) {
            value += denominator_ld[static_cast<std::size_t>(j)] *
                     coefficient_at(i - j);
        }
        numerator_ld[static_cast<std::size_t>(i)] = value;
    }

    // 转换为 double 并输出
    std::vector<double> denominator(denominator_degree + 1);
    std::vector<double> numerator(numerator_degree + 1);
    for (std::size_t i = 0; i <= static_cast<std::size_t>(denominator_degree); ++i) {
        denominator[i] = static_cast<double>(denominator_ld[i]);
    }
    for (std::size_t i = 0; i <= static_cast<std::size_t>(numerator_degree); ++i) {
        numerator[i] = static_cast<double>(numerator_ld[i]);
    }

    const std::string base = shifted_series_base(variable_name, center);
    const std::string numerator_text = polynomial_to_string(numerator, base);
    const std::string denominator_text = polynomial_to_string(denominator, base);
    if (denominator_text == "1") {
        return simplify_symbolic_text(numerator_text);
    } else {
        return simplify_symbolic_text(
            "(" + numerator_text + ") / (" + denominator_text + ")");
    }
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
    double ratios[3];
    for (int n = 10; n <= 12; ++n) {
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
        ratios[n - 10] = val_n1 / val_n;
    }

    // 检查比值是否稳定
    if (mymath::abs(ratios[2] - ratios[1]) > 1e-6 || mymath::abs(ratios[1] - ratios[0]) > 1e-6) {
        return false;
    }

    *ratio = ratios[2];

    // 提取多项式系数 P(n)
    // f(n) = P(n) * r^n, 所以 P(n) = f(n) / r^n
    poly_coeffs->clear();
    for (int n = 0; n <= 3; ++n) {
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
    // 对于多项式，高阶差分最终为零
    std::vector<double> diff = *poly_coeffs;
    for (int order = 0; order < 4; ++order) {
        bool all_zero = true;
        for (double d : diff) {
            if (!mymath::is_near_zero(d, 1e-8)) {
                all_zero = false;
                break;
            }
        }
        if (all_zero && order > 0) {
            return true;  // 是多项式
        }
        std::vector<double> new_diff;
        for (std::size_t i = 1; i < diff.size(); ++i) {
            new_diff.push_back(diff[i] - diff[i - 1]);
        }
        diff = new_diff;
        if (diff.empty()) break;
    }

    return false;
}

/**
 * @brief 计算等差-等比级数的有限和
 *
 * 计算 sum_{k=0}^{n} k * r^k 或更一般的 P(k) * r^k
 */
std::string arith_geo_sum(const std::vector<double>& poly_coeffs,
                           double ratio,
                           const std::string& index_name,
                           const SymbolicExpression& upper_expr,
                           const std::string& lower) {
    // 对于 sum k*r^k，公式为 r(1 - (n+1)r^n + nr^{n+1}) / (1-r)^2
    // 对于 sum P(k)*r^k，可以分解为多个 sum k^m * r^k

    // 简化实现：只处理 P(k) = k 的情况
    if (poly_coeffs.size() >= 2 &&
        mymath::is_near_zero(poly_coeffs[0], 1e-8) &&
        mymath::is_near_zero(poly_coeffs[1] - 1.0, 1e-8)) {
        // sum k * r^k
        std::string r = format_symbolic_scalar(ratio);
        std::string lower_val = lower;

        // 公式: r * (1 - (n+1)*r^n + n*r^{n+1}) / (1-r)^2
        std::ostringstream formula;
        formula << "(" << r << " * (1 - (" << index_name << " + 1) * (" << r
                << ") ^ " << index_name << " + " << index_name << " * (" << r
                << ") ^ (" << index_name << " + 1))) / (1 - (" << r << ")) ^ 2";

        SymbolicExpression primitive = SymbolicExpression::parse(formula.str()).simplify();

        SymbolicExpression lower_minus_one = SymbolicExpression::parse("(" + lower + ") - 1").simplify();
        return SymbolicExpression::parse(
                   "(" + primitive.substitute(index_name, upper_expr).to_string() +
                   ") - (" + primitive.substitute(index_name, lower_minus_one).to_string() + ")")
            .simplify()
            .to_string();
    }

    throw std::runtime_error("arithmetic-geometric series with higher degree polynomials not yet supported");
}

/**
 * @brief 检测差分抵消级数 (Telescoping Series)
 *
 * 检测形如 f(n) = g(n) - g(n+1) 的级数
 */
bool detect_telescoping(const SymbolicExpression& summand,
                        const std::string& /* index_name */,
                        SymbolicExpression* /* g_function */) {
    // 尝试找到 g(n) 使得 f(n) = g(n) - g(n+1)
    // 这是一个困难的问题，我们使用启发式方法

    // 检查是否是分式形式，分子可能抵消
    if (summand.node_->type != NodeType::kDivide) {
        return false;
    }

    // 尝试 g(n) = 1/(P(n)) 形式
    SymbolicExpression numerator(summand.node_->left);
    SymbolicExpression denominator(summand.node_->right);

    // 如果分子是常数，检查分母是否是 P(n) 形式
    double num_val = 0.0;
    if (numerator.is_number(&num_val) && mymath::is_near_zero(num_val - 1.0, 1e-10)) {
        // 检查是否是 1/(n+a) - 1/(n+b) 形式
        // 这需要更复杂的分析，暂时跳过
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
 * @param max_iterations 最大迭代次数
 * @param tolerance 收敛容差
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
    // epsilon 表，只需要两行
    std::vector<double> e_prev(n + 1, 0.0);
    std::vector<double> e_curr(n + 1, 0.0);

    // 初始化：e_0^{(n)} = S_n
    for (int i = 0; i < n; ++i) {
        e_prev[i + 1] = partial_sums[i];
    }

    double best_estimate = partial_sums.back();
    double best_change = mymath::kDoubleMax;

    // 递归计算 epsilon 表
    for (int k = 1; k <= max_iterations && k < n; ++k) {
        bool all_valid = true;
        for (int j = k; j <= n - k; ++j) {
            double diff = e_prev[j + 1] - e_prev[j - 1];
            if (mymath::abs(diff) < 1e-30) {
                // 避免除零
                e_curr[j] = e_prev[j];
            } else {
                // epsilon 递推公式
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
                    }
                    // 检查收敛
                    if (change < tolerance) {
                        *result = estimate;
                        return true;
                    }
                }
            }
        }

        std::swap(e_prev, e_curr);
    }

    *result = best_estimate;
    return best_change < tolerance * 10;  // 放宽容差要求
}

/**
 * @brief 数值计算无穷级数的和
 *
 * 使用部分和序列 + Wynn-epsilon 加速来估计收敛值。
 *
 * @param ctx 级数上下文
 * @param summand 通项表达式
 * @param index_name 指标名
 * @param lower 下界
 * @param max_terms 最大项数
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
        lower_val = 1.0;  // 默认从 1 开始
    }

    // 计算部分和序列
    std::vector<double> partial_sums;
    partial_sums.reserve(max_terms);
    double running_sum = 0.0;
    double prev_term = mymath::kDoubleMax;

    for (int n = 0; n < max_terms; ++n) {
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

        // 检查项是否趋于零（收敛的必要条件）
        if (n > 10 && mymath::abs(term) > mymath::abs(prev_term) * 1.1) {
            // 项不递减，可能发散
            return false;
        }

        // 检查项是否足够小
        if (mymath::abs(term) < 1e-15 * mymath::abs(running_sum)) {
            break;
        }

        prev_term = term;
    }

    if (partial_sums.size() < 5) {
        return false;
    }

    // 使用 Wynn-epsilon 加速
    return wynn_epsilon_acceleration(partial_sums, 20, 1e-10, result);
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

}  // namespace series_ops

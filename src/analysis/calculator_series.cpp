// ============================================================================
// 级数展开命令实现
// ============================================================================

#include "calculator_series.h"
#include "symbolic_expression_internal.h"
#include "statistics/probability.h"

#include "polynomial.h"
#include "mymath.h"

#include <sstream>

namespace series_ops {

namespace internal {

std::vector<double> ps_add(const std::vector<double>& a, const std::vector<double>& b, int degree) {
    std::vector<double> res(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        res[i] = (i < static_cast<int>(a.size()) ? a[i] : 0.0) +
                 (i < static_cast<int>(b.size()) ? b[i] : 0.0);
    }
    return res;
}

std::vector<double> ps_sub(const std::vector<double>& a, const std::vector<double>& b, int degree) {
    std::vector<double> res(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        res[i] = (i < static_cast<int>(a.size()) ? a[i] : 0.0) -
                 (i < static_cast<int>(b.size()) ? b[i] : 0.0);
    }
    return res;
}

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

std::vector<double> ps_div(const std::vector<double>& a, const std::vector<double>& b, int degree) {
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
    if (start_a == -1) return std::vector<double>(degree + 1, 0.0);

    // 计算位移：x^start_a / x^start_b = x^(start_a - start_b)
    int shift = start_a - start_b;
    
    // 如果 shift < 0，结果包含负幂项（Laurent 级数），目前 Taylor PSA 不支持。
    // 但对于极限计算，这表示无穷大。这里抛出异常以便回退到数值方法。
    if (shift < 0) throw std::runtime_error("Laurent series result (infinite limit)");

    // 提取有效部分并相除
    std::vector<double> a_effective;
    for (int i = start_a; i < static_cast<int>(a.size()); ++i) a_effective.push_back(a[i]);
    
    std::vector<double> b_effective;
    for (int i = start_b; i < static_cast<int>(b.size()); ++i) b_effective.push_back(b[i]);

    std::vector<double> res_effective(degree + 1, 0.0);
    double inv_b0 = 1.0 / b_effective[0];
    for (int i = 0; i <= degree; ++i) {
        double val = (i < static_cast<int>(a_effective.size()) ? a_effective[i] : 0.0);
        for (int j = 1; j <= i; ++j) {
            if (j < static_cast<int>(b_effective.size())) val -= b_effective[j] * res_effective[i - j];
        }
        res_effective[i] = val * inv_b0;
    }

    // 处理 shift > 0 的情况：在前面补零
    std::vector<double> final_res(degree + 1, 0.0);
    for (int i = 0; i <= degree; ++i) {
        if (i >= shift && i - shift < static_cast<int>(res_effective.size())) {
            final_res[i] = res_effective[i - shift];
        }
    }
    return final_res;
}

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

std::vector<double> ps_sin(const std::vector<double>& a, int degree) {
    std::vector<double> s, c;
    ps_sincos(a, degree, s, c);
    return s;
}

std::vector<double> ps_cos(const std::vector<double>& a, int degree) {
    std::vector<double> s, c;
    ps_sincos(a, degree, s, c);
    return c;
}

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

} // namespace internal

namespace {

std::string simplify_symbolic_text(const std::string& text) {
    return SymbolicExpression::parse(text).simplify().to_string();
}

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

    std::vector<double> psa_result;
    if (series_ops::internal::evaluate_psa(expression.simplify(), variable_name, center, degree, psa_result, ctx)) {
        return psa_result;
    }

    // Fallback to symbolic derivative
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

    auto coefficient_at = [&](int index) {
        if (index < 0 || index >= static_cast<int>(coefficients.size())) {
            return 0.0;
        }
        return coefficients[static_cast<std::size_t>(index)];
    };

    std::vector<double> denominator(denominator_degree + 1, 0.0);
    denominator[0] = 1.0;
    if (denominator_degree > 0) {
        std::vector<std::vector<double>> matrix(
            static_cast<std::size_t>(denominator_degree),
            std::vector<double>(static_cast<std::size_t>(denominator_degree), 0.0));
        std::vector<double> rhs(static_cast<std::size_t>(denominator_degree), 0.0);
        for (int row = 0; row < denominator_degree; ++row) {
            for (int col = 0; col < denominator_degree; ++col) {
                matrix[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(col)] =
                    coefficient_at(numerator_degree + row - col);
            }
            rhs[static_cast<std::size_t>(row)] =
                -coefficient_at(numerator_degree + row + 1);
        }
        const std::vector<double> solved = solve_dense_linear_system(
            matrix, rhs, "pade");
        for (int i = 0; i < denominator_degree; ++i) {
            denominator[static_cast<std::size_t>(i + 1)] =
                solved[static_cast<std::size_t>(i)];
        }
    }

    std::vector<double> numerator(numerator_degree + 1, 0.0);
    for (int i = 0; i <= numerator_degree; ++i) {
        double value = 0.0;
        for (int j = 0; j <= denominator_degree && j <= i; ++j) {
            value += denominator[static_cast<std::size_t>(j)] *
                     coefficient_at(i - j);
        }
        numerator[static_cast<std::size_t>(i)] = value;
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

std::string puiseux(const SeriesContext& ctx,
                    const std::string& expr,
                    double center,
                    int degree,
                    int denominator) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

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
            auto get_bernoulli = [](int m) -> double {
                static std::vector<double> B = {1.0, 0.5};
                while (B.size() <= static_cast<std::size_t>(m)) {
                    int n = B.size();
                    double sum = 0.0;
                    for (int k = 0; k < n; ++k) {
                        sum += prob::nCr(n + 1, k) * B[k];
                    }
                    B.push_back((n + 1.0 - sum) / prob::nCr(n + 1, n));
                }
                return B[m];
            };

            std::vector<std::string> pieces;
            for (std::size_t p = 0; p < coefficients.size(); ++p) {
                if (mymath::is_near_zero(coefficients[p], 1e-10)) {
                    continue;
                }
                
                std::ostringstream poly_part;
                poly_part << "(1 / " << (p + 1) << ") * (";
                bool first = true;
                for (std::size_t j = 0; j <= p; ++j) {
                    double bj = get_bernoulli(j);
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
                    "series_sum does not support infinite polynomial sums");
            }
            return "0";
        }

        const SymbolicExpression primitive =
            make_polynomial_sum_primitive(polynomial_coefficients);
        return finite_sum_from_primitive(primitive);
    }

    if (upper_is_infinite && lower == "1") {
        for (int k = 1; k <= 5; ++k) {
            int two_k = 2 * k;
            std::string inv_pow_str = "1 / (" + index_name + " ^ " + std::to_string(two_k) + ")";
            SymbolicExpression diff = SymbolicExpression::parse("(" + summand.to_string() + ") - (" + inv_pow_str + ")").simplify();
            double val = 0.0;
            if (diff.is_number(&val) && mymath::is_near_zero(val, 1e-10)) {
                // zeta(2k) = (-1)^(k+1) (2pi)^(2k) B_{2k} / (2(2k)!)
                double B = prob::bernoulli(two_k);
                double coeff = mymath::pow(2.0, two_k - 1) * mymath::abs(B) / prob::factorial(two_k);
                std::ostringstream ss;
                if (!mymath::is_near_zero(coeff - 1.0, 1e-10)) {
                    ss << format_symbolic_scalar(coeff) << " * ";
                }
                ss << "pi ^ " << two_k;
                return ctx.simplify_symbolic(ss.str());
            }
        }
    }

    auto geometric_ratio = [&](double* coefficient, double* ratio) {
        SymbolicExpression n_plus_1 = SymbolicExpression::parse("(" + index_name + ") + 1");
        SymbolicExpression next_term = summand.substitute(index_name, n_plus_1);
        SymbolicExpression ratio_expr = SymbolicExpression::parse("(" + next_term.to_string() + ") / (" + summand.to_string() + ")").simplify();
        double symbolic_ratio = 0.0;
        if (ratio_expr.is_number(&symbolic_ratio)) {
            *ratio = symbolic_ratio;
            SymbolicExpression coeff_expr = SymbolicExpression::parse("(" + summand.to_string() + ") / ((" + format_symbolic_scalar(symbolic_ratio) + ") ^ (" + index_name + "))").simplify();
            double symbolic_coeff = 0.0;
            if (coeff_expr.is_number(&symbolic_coeff)) {
                *coefficient = symbolic_coeff;
                return true;
            }
        }

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

    double geometric_coefficient = 0.0;
    double geometric_ratio_value = 0.0;
    if (!geometric_ratio(&geometric_coefficient, &geometric_ratio_value)) {
        throw std::runtime_error(
            "series_sum currently supports polynomial summands up to degree 3 and common geometric series");
    }

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

bool is_series_command(const std::string& command) {
    return command == "taylor" ||
           command == "pade" ||
           command == "puiseux" ||
           command == "series_sum" ||
           command == "summation";
}

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

// ============================================================================
// Taylor 级数展开实现
// ============================================================================
//
// 本文件实现了 Taylor 级数展开算法，包括：
// - 预定义函数的 Taylor 展开（sin, cos, exp, ln 等）
// - 符号表达式的 Taylor 展开
// - 级数系数计算与格式化
//
// Taylor 级数是函数在某点附近的幂级数展开：
// f(x) = Σ f^(n)(a) / n! * (x-a)^n
//
// 相关文件：
// - taylor_series.h: Taylor 级数接口定义
// - psa_engine.cpp: Puiseux 级数引擎（更通用的级数展开）
// ============================================================================

#include "analysis/series/taylor_series.h"
#include "analysis/modules/series_module.h"
#include "core/scalar_type.h"
#include "core/format_utils.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "statistics/probability.h"
#include "math/mymath.h"
#include "math/helpers/integer_helpers.h"
#include "expression_utils.h"
#include "string_utils.h"
#include <sstream>
#include <iomanip>
#include <map>

namespace series_ops {
namespace taylor {

using internal::PoleException;

namespace predefined_series {

bool try_predefined_taylor(const SymbolicExpression& expr,
                            const std::string& var_name,
                            Scalar center,
                            int degree,
                            std::vector<Scalar>& result) {
    if (!expr.node_) return false;

    if (expr.node_->type == NodeType::kFunction) {
        const std::string& func_name = expr.node_->text;
        SymbolicExpression arg(expr.node_->left);

        Scalar arg_coeff = 0.0L, arg_const = 0.0L;
        bool is_linear_arg = false;

        if (arg.node_->type == NodeType::kVariable && arg.node_->text == var_name) {
            arg_coeff = 1.0L;
            arg_const = -center;
            is_linear_arg = true;
        } else if (arg.node_->type == NodeType::kSubtract) {
            SymbolicExpression left(arg.node_->left);
            SymbolicExpression right(arg.node_->right);
            if (left.node_->type == NodeType::kVariable && left.node_->text == var_name &&
                right.is_number(&arg_const)) {
                arg_coeff = 1.0L;
                arg_const = arg_const - center;
                is_linear_arg = true;
            }
        }

        if (mymath::is_near_zero(center, 1e-10) && arg.node_->type == NodeType::kVariable && arg.node_->text == var_name) {
            result.assign(degree + 1, Scalar(0));

            if (func_name == "exp") {
                result[0] = Scalar(1);
                for (int n = 1; n <= degree; ++n) {
                    result[n] = result[n - 1] / Scalar(n);
                }
                return true;
            }

            if (func_name == "sin") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 1) {
                        int sign = (n % 4 == 1) ? 1 : -1;
                        result[n] = Scalar(sign) / Scalar(prob::factorial(n));
                    }
                }
                return true;
            }

            if (func_name == "cos") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 0) {
                        int sign = (n % 4 == 0) ? 1 : -1;
                        result[n] = Scalar(sign) / Scalar(prob::factorial(n));
                    }
                }
                return true;
            }

            if (func_name == "sinh") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 1) {
                        result[n] = Scalar(1) / Scalar(prob::factorial(n));
                    }
                }
                return true;
            }

            if (func_name == "cosh") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 0) {
                        result[n] = Scalar(1) / Scalar(prob::factorial(n));
                    }
                }
                return true;
            }

            if (func_name == "ln") {
                return false;
            }

            if (func_name == "arcsin" || func_name == "asin") {
                result[1] = Scalar(1);
                for (int n = 3; n <= degree; n += 2) {
                    int m = (n - 1) / 2;
                    result[n] = result[n - 2] * Scalar(2 * m - 1) / Scalar(2 * m) * Scalar(m) / Scalar(m + 1);
                }
                return true;
            }

            if (func_name == "arctan" || func_name == "atan") {
                for (int n = 0; n <= degree; ++n) {
                    if (n % 2 == 1) {
                        int sign = (n % 4 == 1) ? 1 : -1;
                        result[n] = Scalar(sign) / Scalar(n);
                    }
                }
                return true;
            }

            if (func_name == "sqrt") {
                result[0] = Scalar(1);
                for (int n = 1; n <= degree; ++n) {
                    Scalar coeff = Scalar(1);
                    for (int k = 0; k < n; ++k) {
                        coeff = coeff * (Scalar(0.5L) - Scalar(k)) / Scalar(n - k);
                    }
                    result[n] = coeff;
                }
                return true;
            }
        }

        if (is_linear_arg) {
            Scalar h_coeff = arg_coeff;
            result.assign(degree + 1, Scalar(0));

            if (func_name == "exp") {
                Scalar exp_center = mymath::precise128::exp(Scalar(center));
                result[0] = exp_center;
                for (int n = 1; n <= degree; ++n) {
                    result[n] = result[n - 1] * Scalar(h_coeff) / Scalar(n);
                }
                return true;
            }

            if (func_name == "sin") {
                Scalar sin_a = mymath::precise128::sin(Scalar(center));
                Scalar cos_a = mymath::precise128::cos(Scalar(center));
                for (int n = 0; n <= degree; ++n) {
                    Scalar h_pow = mymath::precise128::pow(Scalar(h_coeff), Scalar(n));
                    if (n % 2 == 0) {
                        int sign = (n % 4 == 0) ? 1 : -1;
                        result[n] = sin_a * Scalar(sign) * h_pow / Scalar(prob::factorial(n));
                    } else {
                        int sign = (n % 4 == 1) ? 1 : -1;
                        result[n] = cos_a * Scalar(sign) * h_pow / Scalar(prob::factorial(n));
                    }
                }
                return true;
            }

            if (func_name == "cos") {
                Scalar sin_a = mymath::precise128::sin(Scalar(center));
                Scalar cos_a = mymath::precise128::cos(Scalar(center));
                for (int n = 0; n <= degree; ++n) {
                    Scalar h_pow = mymath::precise128::pow(Scalar(h_coeff), Scalar(n));
                    if (n % 2 == 0) {
                        int sign = (n % 4 == 0) ? 1 : -1;
                        result[n] = cos_a * Scalar(sign) * h_pow / Scalar(prob::factorial(n));
                    } else {
                        int sign = (n % 4 == 1) ? -1 : 1;
                        result[n] = sin_a * Scalar(sign) * h_pow / Scalar(prob::factorial(n));
                    }
                }
                return true;
            }
        }
    }

    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exponent(expr.node_->right);

        Scalar exp_val = 0.0L;
        if (base.node_->type == NodeType::kVariable && base.node_->text == var_name &&
            exponent.is_number(&exp_val) && mymath::is_near_zero(center, 1e-10)) {
            if (exp_val > 0 && mymath::is_integer(exp_val, 1e-10)) {
                int p = static_cast<int>(exp_val + 0.5);
                result.assign(degree + 1, Scalar(0));
                if (p <= degree) {
                    result[p] = Scalar(1);
                }
                return true;
            }
        }
    }

    return false;
}

}  // namespace predefined_series

bool compute_taylor_coefficients_ad(const SeriesContext& ctx,
                                    const SymbolicExpression& expression,
                                    const std::string& variable_name,
                                    Scalar center,
                                    int degree,
                                    std::vector<Scalar>& result) {
    result.clear();
    result.reserve(degree + 1);

    auto optimal_step = [](int order) -> Scalar {
        const Scalar eps = Scalar(1e-34L);
        Scalar exponent = Scalar(1) / Scalar(order + 2);
        Scalar ln_eps = Scalar(-78.8L);
        Scalar ln_h = exponent * ln_eps;
        if (ln_h > Scalar(-40)) {
            return mymath::precise128::exp(ln_h);
        } else {
            return mymath::precise128::pow(eps, exponent);
        }
    };

    Scalar f0 = Scalar(0);
    try {
        f0 = Scalar(ctx.evaluate_at(expression, variable_name, center));
    } catch (...) {
        return false;
    }
    if (!mymath::precise128::isfinite(f0)) return false;
    result.push_back(f0);

    auto compute_derivative = [&](int order, Scalar h) -> Scalar {
        if (order == 0) return f0;

        int n_points = 2 * order + 1;
        std::vector<Scalar> points;
        points.reserve(n_points);

        for (int k = -order; k <= order; ++k) {
            Scalar x = Scalar(center) + Scalar(k) * h;
            Scalar val = Scalar(0);
            try {
                val = Scalar(ctx.evaluate_at(expression, variable_name, x.to_long_double()));
            } catch (...) {
                return Scalar(mymath::quiet_nan());
            }
            if (!mymath::precise128::isfinite(val)) {
                return Scalar(mymath::quiet_nan());
            }
            points.push_back(val);
        }

        std::vector<std::vector<Scalar>> diff_table;
        diff_table.push_back(points);

        for (int d = 1; d <= order; ++d) {
            std::vector<Scalar> next_diff;
            next_diff.reserve(diff_table[d - 1].size() - 1);
            for (std::size_t i = 0; i + 1 < diff_table[d - 1].size(); ++i) {
                next_diff.push_back(diff_table[d - 1][i + 1] - diff_table[d - 1][i]);
            }
            diff_table.push_back(next_diff);
        }

        const auto& nth_diff = diff_table[order];
        if (nth_diff.empty()) {
            return Scalar(mymath::quiet_nan());
        }

        Scalar coeff = nth_diff[nth_diff.size() / 2];
        return coeff / mymath::precise128::pow(h, Scalar(order));
    };

    auto richardson_derivative = [&](int order) -> Scalar {
        Scalar base_h = optimal_step(order);
        Scalar h1 = base_h;
        Scalar h2 = base_h / Scalar(2);

        Scalar d1 = compute_derivative(order, h1);
        Scalar d2 = compute_derivative(order, h2);

        if (!mymath::precise128::isfinite(d1) || !mymath::precise128::isfinite(d2)) {
            return compute_derivative(order, base_h);
        }

        Scalar factor = mymath::precise128::pow(Scalar(2), Scalar(order)) - Scalar(1);
        if (mymath::precise128::abs(factor) < Scalar(1e-10L)) {
            return d2;
        }
        return d2 + (d2 - d1) / factor;
    };

    for (int order = 1; order <= degree; ++order) {
        Scalar deriv = richardson_derivative(order);

        if (!mymath::precise128::isfinite(deriv)) {
            Scalar fallback_h = optimal_step(order) / Scalar(10);
            deriv = compute_derivative(order, fallback_h);
        }

        if (!mymath::precise128::isfinite(deriv)) {
            return result.size() > 1;
        }

        result.push_back(deriv / Scalar(prob::factorial(order)));
    }

    return true;
}

std::vector<Scalar> build_taylor_coefficients(
    const SeriesContext& ctx,
    const SymbolicExpression& expression,
    const std::string& variable_name,
    Scalar center,
    int degree) {
    struct TaylorDerivativeCacheEntry {
        SymbolicExpression derivative;
        Scalar value = Scalar(0);
        bool has_value = false;
    };
    static thread_local std::map<std::string, TaylorDerivativeCacheEntry> derivative_cache;
    static constexpr std::size_t kMaxTaylorDerivativeCacheSize = 256;

    std::vector<Scalar> predefined_result;
    if (predefined_series::try_predefined_taylor(expression, variable_name, center, degree, predefined_result)) {
        return predefined_result;
    }

    std::vector<Scalar> psa_result;
    if (internal::evaluate_psa(expression.simplify(), variable_name, center, degree, psa_result, ctx)) {
        return psa_result;
    }

    std::vector<Scalar> ad_result;
    if (degree >= 5 && compute_taylor_coefficients_ad(ctx, expression, variable_name, center, degree, ad_result)) {
        return ad_result;
    }

    const std::string base_key =
        variable_name + "|" + format_symbolic_scalar(center) + "|" +
        expression.simplify().to_string();
    std::vector<Scalar> coefficients;
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
                Scalar(ctx.evaluate_at(found->second.derivative, variable_name, center));
            found->second.has_value = true;
        }
        const Scalar derivative_value = found->second.value;
        if (!mymath::precise128::isfinite(derivative_value)) {
            std::vector<Scalar> ad_fallback;
            if (compute_taylor_coefficients_ad(ctx, expression, variable_name, center, degree, ad_fallback)) {
                return ad_fallback;
            }
            throw std::runtime_error(
                "taylor expansion produced a non-finite coefficient");
        }
        coefficients.push_back(derivative_value / Scalar(prob::factorial(order)));
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

std::string taylor(const SeriesContext& ctx,
                   const std::string& expr,
                   Scalar center,
                   int degree) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

    const std::vector<Scalar> coefficients =
        build_taylor_coefficients(ctx, expression, variable_name, center, degree);
    std::vector<Scalar> ld_coefficients;
    for (const auto& c : coefficients) {
        ld_coefficients.push_back(c.to_long_double());
    }
    return taylor_series_to_string(ld_coefficients, variable_name, center);
}

}  // namespace taylor
}  // namespace series_ops

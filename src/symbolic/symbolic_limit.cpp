// ============================================================================
// 符号极限模块
// ============================================================================
//
// 实现符号极限计算，支持：
// 1. 直接代入
// 2. 不定型检测与处理
// 3. L'Hôpital 法则
// 4. 已知极限模式
// 5. 无穷远点极限
//
// ============================================================================

#include "symbolic/symbolic_limit.h"
#include "symbolic/symbolic_expression_internal.h"
#include "math/mymath.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace symbolic_limit {

namespace {

using namespace symbolic_expression_internal;

// 检查表达式是否为数值常量
bool is_numeric_constant(const SymbolicExpression& expr) {
    return expr.node_->type == NodeType::kNumber ||
           expr.node_->type == NodeType::kPi ||
           expr.node_->type == NodeType::kE;
}

// 获取表达式的数值（如果是常量）
std::optional<double> get_numeric_value(const SymbolicExpression& expr) {
    if (expr.node_->type == NodeType::kNumber) {
        return expr.node_->number_value;
    }
    if (expr.node_->type == NodeType::kPi) {
        return mymath::kPi;
    }
    if (expr.node_->type == NodeType::kE) {
        return mymath::kE;
    }
    return std::nullopt;
}

// 判断表达式在变量趋于某值时是否趋于 0
bool tends_to_zero(const SymbolicExpression& expr, const std::string& var, const BoundArgument& point) {
    // 简单情况：表达式不包含变量
    if (expr.to_string().find(var) == std::string::npos) {
        auto val = get_numeric_value(expr);
        return val.has_value() && mymath::is_near_zero(*val, 1e-12);
    }
    // 变量本身趋于某值
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        if (point.is_finite()) {
            return mymath::is_near_zero(point.value, 1e-12);
        }
        return false;
    }
    return false;
}

// 判断表达式在变量趋于某值时是否趋于无穷
bool tends_to_infinity(const SymbolicExpression& expr, const std::string& var, const BoundArgument& point) {
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        return point.is_infinite();
    }
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);
        // 如果分母趋于 0 且分子不趋于 0
        if (tends_to_zero(den, var, point) && !tends_to_zero(num, var, point)) {
            return true;
        }
    }
    return false;
}

// 提取分式的分子和分母
bool extract_numerator_denominator(const SymbolicExpression& expr,
                                   SymbolicExpression* numerator,
                                   SymbolicExpression* denominator) {
    if (expr.node_->type == NodeType::kDivide) {
        *numerator = SymbolicExpression(expr.node_->left);
        *denominator = SymbolicExpression(expr.node_->right);
        return true;
    }
    *numerator = expr;
    *denominator = SymbolicExpression::number(1.0);
    return false;
}

// 计算表达式在某点的值（用于直接代入）
std::optional<double> evaluate_at_point(const SymbolicExpression& expr,
                                        const std::string& var,
                                        double point) {
    SymbolicExpression substituted = expr.substitute(var, SymbolicExpression::number(point));
    substituted = substituted.simplify();
    double val = 0.0;
    if (substituted.is_number(&val)) {
        return val;
    }
    return std::nullopt;
}

}  // namespace

// ============================================================================
// SymbolicLimitEngine 实现
// ============================================================================

LimitResult SymbolicLimitEngine::compute_limit(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point,
    int direction) {

    // 策略 1: 直接代入（仅对有限点）
    if (point.is_finite()) {
        auto direct_val = try_direct_substitution(expr, var, point.value);
        if (direct_val.has_value() && mymath::isfinite(*direct_val)) {
            return LimitResult::elementary(
                SymbolicExpression::number(*direct_val),
                "direct_substitution");
        }
    }

    // 策略 2: 无穷远点极限
    if (point.is_infinite()) {
        return limit_at_infinity(expr, var, direction);
    }

    // 策略 3: 检测不定型
    IndeterminateForm form = detect_indeterminate_form(expr, var, point);

    // 策略 4: 尝试已知模式
    LimitResult pattern_result;
    if (try_known_pattern(expr, var, point, direction, &pattern_result)) {
        return pattern_result;
    }

    // 策略 5: 根据不定型处理
    switch (form) {
        case IndeterminateForm::kZeroOverZero:
        case IndeterminateForm::kInfOverInf: {
            // 应用 L'Hôpital 法则
            SymbolicExpression num, den;
            extract_numerator_denominator(expr, &num, &den);
            LimitResult lhopital_result;
            if (apply_lhopital(num, den, var, point, direction, &lhopital_result)) {
                return lhopital_result;
            }
            break;
        }
        case IndeterminateForm::kZeroTimesInf: {
            // 转换为 0/0 型: a * b → a / (1/b)
            if (expr.node_->type == NodeType::kMultiply) {
                SymbolicExpression left(expr.node_->left);
                SymbolicExpression right(expr.node_->right);
                // 尝试转换为分式
                SymbolicExpression converted = make_divide(left, make_power(right, SymbolicExpression::number(-1.0)));
                return compute_limit(converted.simplify(), var, point, direction);
            }
            break;
        }
        case IndeterminateForm::kInfMinusInf: {
            // 尝试合并为单个分式
            if (expr.node_->type == NodeType::kSubtract) {
                SymbolicExpression left(expr.node_->left);
                SymbolicExpression right(expr.node_->right);
                // 尝试提取公分母
                SymbolicExpression combined = (left - right).simplify();
                return compute_limit(combined, var, point, direction);
            }
            break;
        }
        case IndeterminateForm::kZeroToZero: {
            // x^0 → 1 (当 x → 0 时，需要考虑左右极限)
            if (expr.node_->type == NodeType::kPower) {
                SymbolicExpression base(expr.node_->left);
                SymbolicExpression exp(expr.node_->right);
                // 如果底数趋于 0，指数趋于 0
                // 使用 exp(exp * ln(base)) 变换
                SymbolicExpression transformed = make_function("exp",
                    make_multiply(exp, make_function("ln", base)));
                return compute_limit(transformed.simplify(), var, point, direction);
            }
            break;
        }
        case IndeterminateForm::kOneToInf: {
            // 1^inf 型: 使用 exp(inf * ln(1)) = exp(inf * 0)
            if (expr.node_->type == NodeType::kPower) {
                SymbolicExpression base(expr.node_->left);
                SymbolicExpression exp(expr.node_->right);
                LimitResult result;
                if (handle_one_to_infinity(base, exp, var, point, &result)) {
                    return result;
                }
            }
            break;
        }
        default:
            break;
    }

    // 策略 6: 泰勒展开（对于有限点）
    if (point.is_finite()) {
        // 尝试使用级数展开
        SymbolicExpression series = expr;  // 简化版本，实际应调用泰勒展开
        // 这里可以集成现有的 PSA 功能
    }

    return LimitResult::unknown();
}

IndeterminateForm SymbolicLimitEngine::detect_indeterminate_form(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point) {

    // 检查除法形式
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        bool num_zero = is_zero_at_point(num, var, point);
        bool num_inf = is_infinite_at_point(num, var, point);
        bool den_zero = is_zero_at_point(den, var, point);
        bool den_inf = is_infinite_at_point(den, var, point);

        if (num_zero && den_zero) return IndeterminateForm::kZeroOverZero;
        if (num_inf && den_inf) return IndeterminateForm::kInfOverInf;
    }

    // 检查乘法形式
    if (expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        bool left_zero = is_zero_at_point(left, var, point);
        bool right_inf = is_infinite_at_point(right, var, point);
        bool right_zero = is_zero_at_point(right, var, point);
        bool left_inf = is_infinite_at_point(left, var, point);

        if ((left_zero && right_inf) || (left_inf && right_zero)) {
            return IndeterminateForm::kZeroTimesInf;
        }
    }

    // 检查减法形式
    if (expr.node_->type == NodeType::kSubtract) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        if (is_infinite_at_point(left, var, point) &&
            is_infinite_at_point(right, var, point)) {
            return IndeterminateForm::kInfMinusInf;
        }
    }

    // 检查幂形式
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        SymbolicExpression exp(expr.node_->right);

        bool base_zero = is_zero_at_point(base, var, point);
        bool base_one = false;  // 需要检查是否趋于 1
        bool exp_zero = is_zero_at_point(exp, var, point);
        bool exp_inf = is_infinite_at_point(exp, var, point);

        // 检查 base 是否趋于 1
        if (point.is_finite()) {
            auto base_val = evaluate_at_point(base, var, point.value);
            if (base_val.has_value() && mymath::is_near_zero(*base_val - 1.0, 1e-10)) {
                base_one = true;
            }
        }

        if (base_zero && exp_zero) return IndeterminateForm::kZeroToZero;
        if (is_infinite_at_point(base, var, point) && exp_zero) return IndeterminateForm::kInfToZero;
        if (base_one && exp_inf) return IndeterminateForm::kOneToInf;
    }

    return IndeterminateForm::kNone;
}

bool SymbolicLimitEngine::apply_lhopital(
    const SymbolicExpression& numerator,
    const SymbolicExpression& denominator,
    const std::string& var,
    const BoundArgument& point,
    int direction,
    LimitResult* result) {

    // 计算导数
    SymbolicExpression num_deriv = numerator.derivative(var).simplify();
    SymbolicExpression den_deriv = denominator.derivative(var).simplify();

    // 检查导数是否有效
    if (expr_is_zero(den_deriv)) {
        return false;
    }

    // 递归计算极限
    SymbolicExpression ratio = make_divide(num_deriv, den_deriv).simplify();
    LimitResult new_result = compute_limit(ratio, var, point, direction);

    if (new_result.is_definite) {
        new_result.method_used = "lhopital";
        *result = new_result;
        return true;
    }

    return false;
}

bool SymbolicLimitEngine::try_known_pattern(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point,
    int direction,
    LimitResult* result) {

    // 获取表达式的字符串形式用于模式匹配
    std::string expr_str = expr.to_string();

    // 模式 1: sin(x)/x → 1 (x → 0)
    if (point.is_finite() && mymath::is_near_zero(point.value, 1e-12)) {
        // 检查 sin(var)/var 形式
        if (expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expr.node_->left);
            SymbolicExpression den(expr.node_->right);

            if (num.node_->type == NodeType::kFunction && num.node_->text == "sin") {
                SymbolicExpression arg(num.node_->left);
                // 检查参数是否为 var 或 k*var
                if (arg.node_->type == NodeType::kVariable && arg.node_->text == var) {
                    if (den.node_->type == NodeType::kVariable && den.node_->text == var) {
                        *result = LimitResult::elementary(SymbolicExpression::number(1.0), "known_pattern: sin(x)/x");
                        return true;
                    }
                }
                // sin(kx)/(kx) → 1
                if (arg.node_->type == NodeType::kMultiply) {
                    double k = 0.0;
                    SymbolicExpression other;
                    if (SymbolicExpression(arg.node_->left).is_number(&k) &&
                        SymbolicExpression(arg.node_->right).node_->type == NodeType::kVariable &&
                        SymbolicExpression(arg.node_->right).node_->text == var) {
                        SymbolicExpression expected_den = make_multiply(SymbolicExpression::number(k), SymbolicExpression::variable(var));
                        if (expressions_match(den, expected_den)) {
                            *result = LimitResult::elementary(SymbolicExpression::number(1.0), "known_pattern: sin(kx)/(kx)");
                            return true;
                        }
                    }
                }
            }
        }

        // 模式 2: tan(x)/x → 1 (x → 0)
        if (expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expr.node_->left);
            SymbolicExpression den(expr.node_->right);

            if (num.node_->type == NodeType::kFunction && num.node_->text == "tan") {
                SymbolicExpression arg(num.node_->left);
                if (arg.node_->type == NodeType::kVariable && arg.node_->text == var) {
                    if (den.node_->type == NodeType::kVariable && den.node_->text == var) {
                        *result = LimitResult::elementary(SymbolicExpression::number(1.0), "known_pattern: tan(x)/x");
                        return true;
                    }
                }
            }
        }

        // 模式 3: (exp(x) - 1)/x → 1 (x → 0)
        if (expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expr.node_->left);
            SymbolicExpression den(expr.node_->right);

            if (num.node_->type == NodeType::kSubtract) {
                SymbolicExpression left(num.node_->left);
                SymbolicExpression right(num.node_->right);

                if (left.node_->type == NodeType::kFunction && left.node_->text == "exp") {
                    SymbolicExpression arg(left.node_->left);
                    if (arg.node_->type == NodeType::kVariable && arg.node_->text == var &&
                        right.is_number(nullptr) && mymath::is_near_zero(right.node_->number_value - 1.0, 1e-12)) {
                        if (den.node_->type == NodeType::kVariable && den.node_->text == var) {
                            *result = LimitResult::elementary(SymbolicExpression::number(1.0), "known_pattern: (exp(x)-1)/x");
                            return true;
                        }
                    }
                }
            }
        }

        // 模式 4: ln(1+x)/x → 1 (x → 0)
        if (expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expr.node_->left);
            SymbolicExpression den(expr.node_->right);

            if (num.node_->type == NodeType::kFunction && num.node_->text == "ln") {
                SymbolicExpression arg(num.node_->left);
                if (arg.node_->type == NodeType::kAdd) {
                    SymbolicExpression left(arg.node_->left);
                    SymbolicExpression right(arg.node_->right);
                    double one = 0.0;
                    if (left.is_number(&one) && mymath::is_near_zero(one - 1.0, 1e-12) &&
                        right.node_->type == NodeType::kVariable && right.node_->text == var) {
                        if (den.node_->type == NodeType::kVariable && den.node_->text == var) {
                            *result = LimitResult::elementary(SymbolicExpression::number(1.0), "known_pattern: ln(1+x)/x");
                            return true;
                        }
                    }
                }
            }
        }
    }

    // 模式 5: (1 + 1/n)^n → e (n → inf)
    if (point.is_infinite()) {
        if (expr.node_->type == NodeType::kPower) {
            SymbolicExpression base(expr.node_->left);
            SymbolicExpression exp(expr.node_->right);

            if (base.node_->type == NodeType::kAdd && exp.node_->type == NodeType::kVariable && exp.node_->text == var) {
                SymbolicExpression left(base.node_->left);
                SymbolicExpression right(base.node_->right);

                double one = 0.0;
                if (left.is_number(&one) && mymath::is_near_zero(one - 1.0, 1e-12)) {
                    // 检查 right 是否为 1/var
                    if (right.node_->type == NodeType::kDivide) {
                        SymbolicExpression rnum(right.node_->left);
                        SymbolicExpression rden(right.node_->right);
                        double one_check = 0.0;
                        if (rnum.is_number(&one_check) && mymath::is_near_zero(one_check - 1.0, 1e-12) &&
                            rden.node_->type == NodeType::kVariable && rden.node_->text == var) {
                            *result = LimitResult::elementary(
                                SymbolicExpression::parse("e"),
                                "known_pattern: (1+1/n)^n");
                            return true;
                        }
                    }
                }
            }
        }

        // 模式 6: (1 + a/n)^n → e^a (n → inf)
        if (expr.node_->type == NodeType::kPower) {
            SymbolicExpression base(expr.node_->left);
            SymbolicExpression exp(expr.node_->right);

            if (base.node_->type == NodeType::kAdd && exp.node_->type == NodeType::kVariable && exp.node_->text == var) {
                SymbolicExpression left(base.node_->left);
                SymbolicExpression right(base.node_->right);

                double one = 0.0;
                if (left.is_number(&one) && mymath::is_near_zero(one - 1.0, 1e-12)) {
                    // 检查 right 是否为 a/var
                    if (right.node_->type == NodeType::kDivide) {
                        SymbolicExpression rnum(right.node_->left);
                        SymbolicExpression rden(right.node_->right);
                        if (rden.node_->type == NodeType::kVariable && rden.node_->text == var) {
                            double a = 0.0;
                            if (rnum.is_number(&a)) {
                                *result = LimitResult::elementary(
                                    make_function("exp", SymbolicExpression::number(a)),
                                    "known_pattern: (1+a/n)^n");
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }

    // 模式 7: x^x → 1 (x → 0+)
    if (point.is_finite() && mymath::is_near_zero(point.value, 1e-12) && direction >= 0) {
        if (expr.node_->type == NodeType::kPower) {
            SymbolicExpression base(expr.node_->left);
            SymbolicExpression exp(expr.node_->right);

            if (base.node_->type == NodeType::kVariable && base.node_->text == var &&
                exp.node_->type == NodeType::kVariable && exp.node_->text == var) {
                *result = LimitResult::elementary(SymbolicExpression::number(1.0), "known_pattern: x^x (right limit)");
                return true;
            }
        }
    }

    return false;
}

std::optional<double> SymbolicLimitEngine::try_direct_substitution(
    const SymbolicExpression& expr,
    const std::string& var,
    double point) {

    return evaluate_at_point(expr, var, point);
}

LimitResult SymbolicLimitEngine::limit_at_infinity(
    const SymbolicExpression& expr,
    const std::string& var,
    int direction) {

    // 对于无穷远点，分析主导项
    // 简化版本：检查多项式主导项

    // 尝试提取多项式系数
    std::vector<double> coeffs;
    if (expr.polynomial_coefficients(var, &coeffs)) {
        if (coeffs.empty()) {
            return LimitResult::elementary(SymbolicExpression::number(0.0), "polynomial_at_infinity");
        }

        // 找到最高次项
        int degree = static_cast<int>(coeffs.size()) - 1;
        while (degree >= 0 && mymath::is_near_zero(coeffs[degree], 1e-15)) {
            degree--;
        }

        if (degree < 0) {
            return LimitResult::elementary(SymbolicExpression::number(0.0), "zero_polynomial");
        }

        // 多项式在无穷远的行为由最高次项决定
        if (degree == 0) {
            return LimitResult::elementary(SymbolicExpression::number(coeffs[0]), "constant_polynomial");
        }

        // 高次多项式趋于无穷
        double leading_coeff = coeffs[degree];
        bool positive_inf = (leading_coeff > 0) == (degree % 2 == 0 || direction >= 0);
        return LimitResult::infinite(positive_inf, "polynomial_dominant_term");
    }

    // 检查 1/x^n 形式
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        double num_val = 0.0;
        if (num.is_number(&num_val) && !mymath::is_near_zero(num_val, 1e-15)) {
            // 检查分母是否为 x 的幂
            if (den.node_->type == NodeType::kVariable && den.node_->text == var) {
                return LimitResult::elementary(SymbolicExpression::number(0.0), "1/x_at_infinity");
            }
            if (den.node_->type == NodeType::kPower) {
                SymbolicExpression base(den.node_->left);
                double exp = 0.0;
                if (base.node_->type == NodeType::kVariable && base.node_->text == var &&
                    SymbolicExpression(den.node_->right).is_number(&exp) && exp > 0) {
                    return LimitResult::elementary(SymbolicExpression::number(0.0), "1/x^n_at_infinity");
                }
            }
        }
    }

    // 检查指数函数
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
        SymbolicExpression arg(expr.node_->left);
        // exp(-x) → 0 as x → +inf
        if (arg.node_->type == NodeType::kNegate) {
            SymbolicExpression inner(arg.node_->left);
            if (inner.node_->type == NodeType::kVariable && inner.node_->text == var) {
                return LimitResult::elementary(SymbolicExpression::number(0.0), "exp(-x)_at_infinity");
            }
        }
        // exp(x) → inf as x → +inf
        if (arg.node_->type == NodeType::kVariable && arg.node_->text == var) {
            return LimitResult::infinite(true, "exp(x)_at_infinity");
        }
    }

    return LimitResult::unknown();
}

bool SymbolicLimitEngine::handle_one_to_infinity(
    const SymbolicExpression& base,
    const SymbolicExpression& exponent,
    const std::string& var,
    const BoundArgument& point,
    LimitResult* result) {

    // 1^inf 型极限: 使用 exp(inf * ln(1)) = exp(inf * 0)
    // 需要计算 (base - 1) * exponent 的极限

    // 计算 base - 1
    SymbolicExpression base_minus_one = make_subtract(base, SymbolicExpression::number(1.0));

    // 计算 (base - 1) * exponent
    SymbolicExpression product = make_multiply(base_minus_one, exponent);

    // 计算 product 的极限
    LimitResult product_limit = compute_limit(product.simplify(), var, point, 0);

    if (product_limit.is_definite && product_limit.is_elementary) {
        // 结果是 exp(product_limit)
        SymbolicExpression final_result = make_function("exp", product_limit.value).simplify();
        *result = LimitResult::elementary(final_result, "one_to_infinity_transform");
        return true;
    }

    return false;
}

bool SymbolicLimitEngine::is_infinite_at_point(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point) {

    if (expr.node_->type == NodeType::kInfinity) {
        return true;
    }

    // 检查 1/0 形式
    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        if (is_zero_at_point(den, var, point) && !is_zero_at_point(num, var, point)) {
            return true;
        }
    }

    return false;
}

bool SymbolicLimitEngine::is_zero_at_point(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point) {

    // 常量零
    double val = 0.0;
    if (expr.is_number(&val) && mymath::is_near_zero(val, 1e-15)) {
        return true;
    }

    // 变量在趋于零点时
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        if (point.is_finite() && mymath::is_near_zero(point.value, 1e-15)) {
            return true;
        }
    }

    return false;
}

// ============================================================================
// 辅助函数实现
// ============================================================================

std::optional<SymbolicExpression> SymbolicLimitEngine::limit(
    const std::string& expr_str,
    const std::string& var,
    const BoundArgument& point,
    int direction) {

    try {
        SymbolicExpression expr = SymbolicExpression::parse(expr_str);
        SymbolicLimitEngine engine;
        LimitResult result = engine.compute_limit(expr, var, point, direction);

        if (result.is_definite) {
            return result.value;
        }
        return std::nullopt;
    } catch (...) {
        return std::nullopt;
    }
}

bool parse_limit_arguments(
    const std::vector<std::string>& args,
    SymbolicExpression* expr,
    std::string* var,
    BoundArgument* point,
    int* direction) {

    if (args.size() < 3) {
        return false;
    }

    try {
        *expr = SymbolicExpression::parse(args[0]);
        *var = args[1];

        // 解析极限点
        std::string point_str = args[2];
        if (point_str == "inf" || point_str == "infinity" || point_str == "oo" || point_str == "+inf") {
            *point = BoundArgument::pos_inf();
        } else if (point_str == "-inf" || point_str == "-infinity" || point_str == "-oo") {
            *point = BoundArgument::neg_inf();
        } else {
            double p = std::stod(point_str);
            *point = BoundArgument::finite(p);
        }

        // 解析方向
        *direction = 0;  // 默认双侧
        if (args.size() > 3) {
            if (args[3] == "left" || args[3] == "-") {
                *direction = -1;
            } else if (args[3] == "right" || args[3] == "+") {
                *direction = 1;
            }
        }

        return true;
    } catch (...) {
        return false;
    }
}

}  // namespace symbolic_limit

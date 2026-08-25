// ============================================================================
// 符号极限模块
// ============================================================================
//
// 实现符号极限计算，支持：
// 1. 直接代入（连续点）
// 2. 幂级数展开（Taylor / PSA 分析）
// 3. 不定型检测与自适应转换（0/0, inf/inf, 0*inf, inf-inf, 0^0, inf^0, 1^inf）
// 4. 有理分式渐近阶分析
// 5. 共轭根式有理化
// 6. 夹逼准则与有界量乘无穷小
// 7. L'Hôpital 法则
// 8. 扩展经典极限模式与等价无穷小
// 9. 无穷远点极限
//
// ============================================================================

#include "symbolic/calculus/limit/symbolic_limit.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "analysis/series/psa_engine.h"
#include "analysis/modules/series_module.h"
#include "core/execution_context.h"

#include "types/scalar_type.h"
#include "math/mymath.h"
#include "math/numeric/precision/tolerances.h"

#include <algorithm>
#include <functional>
#include <sstream>

namespace symbolic_limit {

namespace {

using namespace symbolic_expression_internal;
using Scalar = mymath::Scalar;

// 提取分式的分子和分母
bool extract_numerator_denominator(const SymbolicExpression& expr,
                                   SymbolicExpression* numerator,
                                   SymbolicExpression* denominator) {
    if (expr.node_ == nullptr) {
        *numerator = expr;
        *denominator = SymbolicExpression::number(Scalar(1.0L));
        return false;
    }
    if (expr.node_->type == NodeType::kDivide) {
        *numerator = SymbolicExpression(expr.node_->left);
        *denominator = SymbolicExpression(expr.node_->right);
        return true;
    }
    *numerator = expr;
    *denominator = SymbolicExpression::number(Scalar(1.0L));
    return false;
}

// 计算表达式在某点的值（用于直接代入）
std::optional<Scalar> evaluate_at_point(const SymbolicExpression& expr,
                                        const std::string& var,
                                        Scalar point) {
    try {
        SymbolicExpression substituted = expr.substitute(var, SymbolicExpression::number(point));
        Scalar val = Scalar(0.0L);
        if (substituted.is_number(&val) && mymath::isfinite(val)) {
            return val;
        }
        core::ExecutionContext ctx;
        StoredValue sv = substituted.evalf(ctx);
        if (sv.is_scalar() && mymath::isfinite(sv.get_decimal())) {
            return sv.get_decimal();
        }
    } catch (...) {
    }
    return std::nullopt;
}

bool extract_linear_term(const SymbolicExpression& expr, const std::string& var, Scalar* k) {
    if (expr.node_ == nullptr) return false;
    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        if (k) *k = Scalar(1.0L);
        return true;
    }
    if (expr.node_->type == NodeType::kMultiply) {
        Scalar num = Scalar(0.0L);
        if (SymbolicExpression(expr.node_->left).is_number(&num) &&
            SymbolicExpression(expr.node_->right).node_->type == NodeType::kVariable &&
            SymbolicExpression(expr.node_->right).node_->text == var) {
            if (k) *k = num;
            return true;
        }
        if (SymbolicExpression(expr.node_->right).is_number(&num) &&
            SymbolicExpression(expr.node_->left).node_->type == NodeType::kVariable &&
            SymbolicExpression(expr.node_->left).node_->text == var) {
            if (k) *k = num;
            return true;
        }
    }
    if (expr.node_->type == NodeType::kNegate) {
        Scalar inner_k = Scalar(0.0L);
        if (extract_linear_term(SymbolicExpression(expr.node_->left), var, &inner_k)) {
            if (k) *k = -inner_k;
            return true;
        }
    }
    return false;
}

bool extract_quadratic_term(const SymbolicExpression& expr, const std::string& var, Scalar* m) {
    if (expr.node_ == nullptr) return false;
    if (expr.node_->type == NodeType::kPower) {
        SymbolicExpression base(expr.node_->left);
        Scalar exp = Scalar(0);
        if (base.is_variable_named(var) && SymbolicExpression(expr.node_->right).is_number(&exp) &&
            mymath::is_near_zero(exp - Scalar(2.0L), Scalar(1e-9L))) {
            if (m) *m = Scalar(1.0L);
            return true;
        }
    }
    if (expr.node_->type == NodeType::kMultiply) {
        Scalar num = Scalar(0.0L);
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);
        if (left.is_number(&num)) {
            Scalar inner_m = Scalar(0);
            if (extract_quadratic_term(right, var, &inner_m)) {
                if (m) *m = num * inner_m;
                return true;
            }
        }
        if (right.is_number(&num)) {
            Scalar inner_m = Scalar(0);
            if (extract_quadratic_term(left, var, &inner_m)) {
                if (m) *m = num * inner_m;
                return true;
            }
        }
    }
    return false;
}

bool extract_one_plus_linear(const SymbolicExpression& expr, const std::string& var, Scalar* k) {
    if (expr.node_ == nullptr || expr.node_->type != NodeType::kAdd) return false;
    SymbolicExpression left(expr.node_->left);
    SymbolicExpression right(expr.node_->right);
    Scalar one = Scalar(0.0L);
    if (left.is_number(&one) && mymath::is_near_zero(one - Scalar(1.0L), Scalar(1e-9L))) {
        return extract_linear_term(right, var, k);
    }
    if (right.is_number(&one) && mymath::is_near_zero(one - Scalar(1.0L), Scalar(1e-9L))) {
        return extract_linear_term(left, var, k);
    }
    return false;
}

bool extract_exp_minus_one(const SymbolicExpression& expr, const std::string& var, Scalar* k) {
    if (expr.node_ == nullptr) return false;
    if (expr.node_->type == NodeType::kSubtract) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);
        Scalar one = Scalar(0.0L);
        if (right.is_number(&one) && mymath::is_near_zero(one - Scalar(1.0L), Scalar(1e-9L)) &&
            left.node_->type == NodeType::kFunction && left.node_->text == "exp") {
            return extract_linear_term(SymbolicExpression(left.node_->left), var, k);
        }
    }
    if (expr.node_->type == NodeType::kAdd) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);
        Scalar neg_one = Scalar(0.0L);
        if (left.is_number(&neg_one) && mymath::is_near_zero(neg_one + Scalar(1.0L), Scalar(1e-9L)) &&
            right.node_->type == NodeType::kFunction && right.node_->text == "exp") {
            return extract_linear_term(SymbolicExpression(right.node_->left), var, k);
        }
        if (right.is_number(&neg_one) && mymath::is_near_zero(neg_one + Scalar(1.0L), Scalar(1e-9L)) &&
            left.node_->type == NodeType::kFunction && left.node_->text == "exp") {
            return extract_linear_term(SymbolicExpression(left.node_->left), var, k);
        }
    }
    return false;
}

bool is_sqrt_node(const SymbolicExpression& expr, SymbolicExpression* inner) {
    if (expr.node_ == nullptr) return false;
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "sqrt") {
        if (inner) *inner = SymbolicExpression(expr.node_->left);
        return true;
    }
    if (expr.node_->type == NodeType::kPower) {
        Scalar exp = Scalar(0);
        if (SymbolicExpression(expr.node_->right).is_number(&exp) &&
            mymath::is_near_zero(exp - Scalar(0.5L), Scalar(1e-9L))) {
            if (inner) *inner = SymbolicExpression(expr.node_->left);
            return true;
        }
    }
    return false;
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

    static thread_local int limit_depth = 0;
    if (limit_depth >= 6) {
        return LimitResult::unknown();
    }
    ++limit_depth;
    struct LimitDepthGuard {
        ~LimitDepthGuard() { --limit_depth; }
    } limit_guard;

    // 0. 特殊情况：如果表达式不含目标变量，极限即为常数表达式自身
    auto all_vars = expr.identifier_variables();
    if (std::find(all_vars.begin(), all_vars.end(), var) == all_vars.end()) {
        return LimitResult::elementary(expr, "constant_expression");
    }

    // 1. 双侧有限点极限：当且仅当左右极限均存在且相等
    if (point.is_finite() && direction == 0) {
        // 先尝试通过解析级数直接判断（若在中心点解析连续）
        LimitResult series_res;
        if (apply_series_expansion(expr, var, point, 0, &series_res)) {
            return series_res;
        }

        const LimitResult left = compute_limit(expr, var, point, -1);
        const LimitResult right = compute_limit(expr, var, point, 1);
        if (!left.is_definite || !right.is_definite) {
            if (left.is_definite && !right.is_definite) return left;
            if (!left.is_definite && right.is_definite) return right;
            return LimitResult::unknown();
        }

        Scalar left_value = Scalar(0), right_value = Scalar(0);
        if (left.value.is_number(&left_value) && right.value.is_number(&right_value)) {
            if (!mymath::is_near_zero(left_value - right_value, Scalar(1e-9L))) {
                return LimitResult::does_not_exist("left and right limits differ");
            }
        } else if (left.value.simplify().to_string() != right.value.simplify().to_string()) {
            return LimitResult::does_not_exist("left and right limits differ");
        }
        return left;
    }

    // 2. 策略 1: 直接代入（仅对有限点且可能连续的函数）
    if (point.is_finite()) {
        bool is_discontinuous = false;
        std::function<void(const SymbolicExpression&)> check_disc = [&](const SymbolicExpression& e) {
            if (!e.node_ || is_discontinuous) return;
            if (e.node_->type == NodeType::kFunction) {
                if (e.node_->text == "floor" || e.node_->text == "ceil" || e.node_->text == "sign" ||
                    e.node_->text == "abs") {
                    is_discontinuous = true;
                    return;
                }
            }
            if (e.node_->left) check_disc(SymbolicExpression(e.node_->left));
            if (e.node_->right) check_disc(SymbolicExpression(e.node_->right));
            for (const auto& child : e.node_->children) {
                check_disc(SymbolicExpression(child));
            }
        };
        check_disc(expr);

        IndeterminateForm form_check = detect_indeterminate_form(expr, var, point);
        if (!is_discontinuous && form_check == IndeterminateForm::kNone) {
            auto direct_val = try_direct_substitution(expr, var, point.value);
            if (direct_val.has_value() && mymath::isfinite(*direct_val)) {
                return LimitResult::elementary(
                    SymbolicExpression::number(*direct_val),
                    "direct_substitution");
            }
        }
    }

    // 3. 策略 2: 无穷远点极限
    if (point.is_infinite()) {
        return limit_at_infinity(expr, var, point, direction);
    }

    // 4. 策略 3: 夹逼准则与有界量乘无穷小
    LimitResult squeeze_result;
    if (apply_squeeze_or_bounded(expr, var, point, direction, &squeeze_result)) {
        return squeeze_result;
    }

    // 5. 策略 4: 尝试已知模式与等价无穷小
    LimitResult pattern_result;
    if (try_known_pattern(expr, var, point, direction, &pattern_result)) {
        return pattern_result;
    }

    // 6. 策略 5: 有理分式多项式极限
    SymbolicExpression num, den;
    if (extract_numerator_denominator(expr, &num, &den)) {
        LimitResult rational_res;
        if (try_rational_limit(num, den, var, point, direction, &rational_res)) {
            return rational_res;
        }
    }

    // 7. 策略 6: 共轭根式有理化
    LimitResult conjugate_res;
    if (try_conjugate_rationalization(expr, var, point, direction, &conjugate_res)) {
        return conjugate_res;
    }

    // 8. 策略 7: 泰勒/幂级数 (PSA) 展开（严格符号级数）
    LimitResult series_result;
    if (apply_series_expansion(expr, var, point, direction, &series_result)) {
        return series_result;
    }

    // 9. 策略 8: 根据不定型分类转化处理
    IndeterminateForm form = detect_indeterminate_form(expr, var, point);
    switch (form) {
        case IndeterminateForm::kZeroOverZero:
        case IndeterminateForm::kInfOverInf: {
            SymbolicExpression n, d;
            extract_numerator_denominator(expr, &n, &d);
            LimitResult lhopital_result;
            if (apply_lhopital(n, d, var, point, direction, &lhopital_result)) {
                return lhopital_result;
            }
            break;
        }
        case IndeterminateForm::kZeroTimesInf: {
            if (expr.node_ && expr.node_->type == NodeType::kMultiply) {
                SymbolicExpression left(expr.node_->left);
                SymbolicExpression right(expr.node_->right);

                // 启发式：如果其中一项包含对数，保留对数在分子，将另一项下放到分母
                bool right_is_log = (right.node_ && right.node_->type == NodeType::kFunction &&
                                     (right.node_->text == "ln" || right.node_->text == "log"));
                bool left_is_log = (left.node_ && left.node_->type == NodeType::kFunction &&
                                    (left.node_->text == "ln" || left.node_->text == "log"));

                SymbolicExpression num_part, den_part;
                if (right_is_log && !left_is_log) {
                    num_part = right;
                    den_part = make_divide(SymbolicExpression::number(1.0L), left);
                } else if (left_is_log && !right_is_log) {
                    num_part = left;
                    den_part = make_divide(SymbolicExpression::number(1.0L), right);
                } else {
                    num_part = left;
                    den_part = make_divide(SymbolicExpression::number(1.0L), right);
                }
                LimitResult lhop_res;
                if (apply_lhopital(num_part, den_part, var, point, direction, &lhop_res)) {
                    return lhop_res;
                }
            }
            break;
        }
        case IndeterminateForm::kInfMinusInf: {
            if (expr.node_ && expr.node_->type == NodeType::kSubtract) {
                SymbolicExpression left(expr.node_->left);
                SymbolicExpression right(expr.node_->right);

                if (left.node_ && left.node_->type == NodeType::kDivide &&
                    right.node_ && right.node_->type == NodeType::kDivide) {
                    SymbolicExpression A(left.node_->left), B(left.node_->right);
                    SymbolicExpression C(right.node_->left), D(right.node_->right);
                    SymbolicExpression combined = make_divide(
                        make_subtract(make_multiply(A, D), make_multiply(B, C)),
                        make_multiply(B, D)
                    ).simplify();
                    return compute_limit(combined, var, point, direction);
                }

                SymbolicExpression combined = (left - right).simplify();
                if (combined.to_string() != expr.to_string()) {
                    return compute_limit(combined, var, point, direction);
                }
            }
            break;
        }
        case IndeterminateForm::kZeroToZero:
        case IndeterminateForm::kInfToZero: {
            if (expr.node_ && expr.node_->type == NodeType::kPower) {
                SymbolicExpression base(expr.node_->left);
                SymbolicExpression exp(expr.node_->right);
                LimitResult result;
                if (handle_power_indeterminate(base, exp, var, point, direction, &result)) {
                    return result;
                }
            }
            break;
        }
        case IndeterminateForm::kOneToInf: {
            if (expr.node_ && expr.node_->type == NodeType::kPower) {
                SymbolicExpression base(expr.node_->left);
                SymbolicExpression exp(expr.node_->right);
                LimitResult result;
                if (handle_one_to_infinity(base, exp, var, point, direction, &result)) {
                    return result;
                }
            }
            break;
        }
        default:
            break;
    }

    // 10. 策略 9: Gruntz 采样探针后备
    LimitResult gruntz_result;
    if (apply_gruntz(expr, var, point, direction, &gruntz_result)) {
        return gruntz_result;
    }

    return LimitResult::unknown();
}

IndeterminateForm SymbolicLimitEngine::detect_indeterminate_form(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point) {

    if (expr.node_ == nullptr) return IndeterminateForm::kNone;

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
        bool base_inf = is_infinite_at_point(base, var, point);
        bool base_one = false;
        bool exp_zero = is_zero_at_point(exp, var, point);
        bool exp_inf = is_infinite_at_point(exp, var, point);

        if (point.is_finite()) {
            auto base_val = evaluate_at_point(base, var, point.value);
            if (base_val.has_value() && mymath::is_near_zero(Scalar(*base_val - 1.0L), Scalar(1e-10))) {
                base_one = true;
            }
        } else {
            LimitResult base_lim = limit_at_infinity(base, var, point, 0);
            Scalar bv = Scalar(0);
            if (base_lim.is_definite && base_lim.value.is_number(&bv) &&
                mymath::is_near_zero(bv - Scalar(1.0L), Scalar(1e-10))) {
                base_one = true;
            }
        }

        if (base_zero && exp_zero) return IndeterminateForm::kZeroToZero;
        if (base_inf && exp_zero) return IndeterminateForm::kInfToZero;
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

    static thread_local int lhopital_depth = 0;
    constexpr int kMaxLhopitalDepth = 8;

    if (lhopital_depth >= kMaxLhopitalDepth) {
        return false;
    }

    ++lhopital_depth;
    struct DepthGuard {
        ~DepthGuard() { --lhopital_depth; }
    } guard;

    SymbolicExpression num_deriv = numerator.derivative(var).simplify();
    SymbolicExpression den_deriv = denominator.derivative(var).simplify();

    if (expr_is_zero(den_deriv)) {
        return false;
    }

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

    // 模式 1：x → 0 经典极限定理
    if (point.is_finite() && mymath::is_near_zero(Scalar(point.value), Scalar(1e-12))) {
        if (expr.node_ && expr.node_->type == NodeType::kDivide) {
            SymbolicExpression num(expr.node_->left);
            SymbolicExpression den(expr.node_->right);

            Scalar k = Scalar(0.0L), m = Scalar(0.0L);

            // sin(kx) / (mx) -> k/m
            if (num.node_ && num.node_->type == NodeType::kFunction && num.node_->text == "sin" &&
                extract_linear_term(SymbolicExpression(num.node_->left), var, &k) &&
                extract_linear_term(den, var, &m) && !mymath::is_near_zero(m, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(k / m), "known_pattern: sin(kx)/(mx)");
                return true;
            }
            if (den.node_ && den.node_->type == NodeType::kFunction && den.node_->text == "sin" &&
                extract_linear_term(SymbolicExpression(den.node_->left), var, &k) &&
                extract_linear_term(num, var, &m) && !mymath::is_near_zero(k, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(m / k), "known_pattern: mx/sin(kx)");
                return true;
            }

            // tan(kx) / (mx) -> k/m
            if (num.node_ && num.node_->type == NodeType::kFunction && num.node_->text == "tan" &&
                extract_linear_term(SymbolicExpression(num.node_->left), var, &k) &&
                extract_linear_term(den, var, &m) && !mymath::is_near_zero(m, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(k / m), "known_pattern: tan(kx)/(mx)");
                return true;
            }
            if (den.node_ && den.node_->type == NodeType::kFunction && den.node_->text == "tan" &&
                extract_linear_term(SymbolicExpression(den.node_->left), var, &k) &&
                extract_linear_term(num, var, &m) && !mymath::is_near_zero(k, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(m / k), "known_pattern: mx/tan(kx)");
                return true;
            }

            // arcsin/asin(kx) / (mx) -> k/m
            if (num.node_ && num.node_->type == NodeType::kFunction &&
                (num.node_->text == "asin" || num.node_->text == "arcsin") &&
                extract_linear_term(SymbolicExpression(num.node_->left), var, &k) &&
                extract_linear_term(den, var, &m) && !mymath::is_near_zero(m, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(k / m), "known_pattern: asin(kx)/(mx)");
                return true;
            }

            // arctan/atan(kx) / (mx) -> k/m
            if (num.node_ && num.node_->type == NodeType::kFunction &&
                (num.node_->text == "atan" || num.node_->text == "arctan") &&
                extract_linear_term(SymbolicExpression(num.node_->left), var, &k) &&
                extract_linear_term(den, var, &m) && !mymath::is_near_zero(m, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(k / m), "known_pattern: atan(kx)/(mx)");
                return true;
            }

            // sinh(kx) / (mx) -> k/m
            if (num.node_ && num.node_->type == NodeType::kFunction && num.node_->text == "sinh" &&
                extract_linear_term(SymbolicExpression(num.node_->left), var, &k) &&
                extract_linear_term(den, var, &m) && !mymath::is_near_zero(m, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(k / m), "known_pattern: sinh(kx)/(mx)");
                return true;
            }

            // tanh(kx) / (mx) -> k/m
            if (num.node_ && num.node_->type == NodeType::kFunction && num.node_->text == "tanh" &&
                extract_linear_term(SymbolicExpression(num.node_->left), var, &k) &&
                extract_linear_term(den, var, &m) && !mymath::is_near_zero(m, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(k / m), "known_pattern: tanh(kx)/(mx)");
                return true;
            }

            // (exp(kx) - 1) / (mx) -> k/m
            if (extract_exp_minus_one(num, var, &k) &&
                extract_linear_term(den, var, &m) && !mymath::is_near_zero(m, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(k / m), "known_pattern: (exp(kx)-1)/(mx)");
                return true;
            }
            if (extract_exp_minus_one(den, var, &k) &&
                extract_linear_term(num, var, &m) && !mymath::is_near_zero(k, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(m / k), "known_pattern: mx/(exp(kx)-1)");
                return true;
            }

            // ln(1 + kx) / (mx) -> k/m
            if (num.node_ && num.node_->type == NodeType::kFunction && num.node_->text == "ln" &&
                extract_one_plus_linear(SymbolicExpression(num.node_->left), var, &k) &&
                extract_linear_term(den, var, &m) && !mymath::is_near_zero(m, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(k / m), "known_pattern: ln(1+kx)/(mx)");
                return true;
            }
            if (den.node_ && den.node_->type == NodeType::kFunction && den.node_->text == "ln" &&
                extract_one_plus_linear(SymbolicExpression(den.node_->left), var, &k) &&
                extract_linear_term(num, var, &m) && !mymath::is_near_zero(k, Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::number(m / k), "known_pattern: mx/ln(1+kx)");
                return true;
            }

            // (1 - cos(kx)) / (m x^2) -> k^2 / (2m)
            if (num.node_ && num.node_->type == NodeType::kSubtract &&
                SymbolicExpression(num.node_->left).is_number(&k) && mymath::is_near_zero(k - Scalar(1.0L), Scalar(1e-9L)) &&
                SymbolicExpression(num.node_->right).node_ &&
                SymbolicExpression(num.node_->right).node_->type == NodeType::kFunction &&
                SymbolicExpression(num.node_->right).node_->text == "cos" &&
                extract_linear_term(SymbolicExpression(SymbolicExpression(num.node_->right).node_->left), var, &k)) {
                Scalar quad_coeff = Scalar(0);
                if (extract_quadratic_term(den, var, &quad_coeff) && !mymath::is_near_zero(quad_coeff, Scalar(1e-9L))) {
                    *result = LimitResult::elementary(
                        SymbolicExpression::number(k * k / (Scalar(2.0L) * quad_coeff)),
                        "known_pattern: (1-cos(kx))/(m*x^2)");
                    return true;
                }
            }
        }

        // (1 + kx)^(m/x) -> exp(k*m) (x -> 0)
        if (expr.node_ && expr.node_->type == NodeType::kPower) {
            SymbolicExpression base(expr.node_->left);
            SymbolicExpression exp(expr.node_->right);

            Scalar k = Scalar(0.0L);
            if (extract_one_plus_linear(base, var, &k)) {
                if (exp.node_ && exp.node_->type == NodeType::kDivide) {
                    SymbolicExpression exp_num(exp.node_->left);
                    SymbolicExpression exp_den(exp.node_->right);
                    Scalar m_val = Scalar(0.0L);
                    if (exp_den.is_variable_named(var) && exp_num.is_number(&m_val)) {
                        Scalar exponent_val = k * m_val;
                        if (mymath::is_near_zero(exponent_val - Scalar(1.0L), Scalar(1e-9L))) {
                            *result = LimitResult::elementary(SymbolicExpression::parse("e"), "known_pattern: (1+kx)^(1/x)");
                        } else {
                            *result = LimitResult::elementary(
                                make_function("exp", SymbolicExpression::number(exponent_val)),
                                "known_pattern: (1+kx)^(m/x)");
                        }
                        return true;
                    }
                }
            }
        }
    }

    // 模式 2：n → inf 经典自然底数极限 (1 + a/n)^(b*n) → e^(a*b)
    if (point.is_infinite()) {
        // x * sin(1/x) -> 1, x * tan(1/x) -> 1
        if (expr.node_ && expr.node_->type == NodeType::kMultiply) {
            SymbolicExpression left(expr.node_->left);
            SymbolicExpression right(expr.node_->right);
            auto check_x_trig_inv = [&](const SymbolicExpression& var_term,
                                        const SymbolicExpression& fn_term) -> bool {
                if (var_term.is_variable_named(var) && fn_term.node_ && fn_term.node_->type == NodeType::kFunction) {
                    if (fn_term.node_->text == "sin" || fn_term.node_->text == "tan" ||
                        fn_term.node_->text == "atan" || fn_term.node_->text == "asin") {
                        SymbolicExpression arg(fn_term.node_->left);
                        if (arg.node_ && arg.node_->type == NodeType::kDivide) {
                            SymbolicExpression n(arg.node_->left);
                            SymbolicExpression d(arg.node_->right);
                            Scalar n_val = Scalar(0);
                            if (n.is_number(&n_val) && d.is_variable_named(var)) {
                                *result = LimitResult::elementary(SymbolicExpression::number(n_val), "known_pattern: x*trig(a/x)");
                                return true;
                            }
                        }
                    }
                }
                return false;
            };
            if (check_x_trig_inv(left, right) || check_x_trig_inv(right, left)) {
                return true;
            }
        }

        if (expr.node_ && expr.node_->type == NodeType::kPower) {
            SymbolicExpression base(expr.node_->left);
            SymbolicExpression exp(expr.node_->right);

            Scalar b_coeff = Scalar(0.0L);
            if (base.node_ && base.node_->type == NodeType::kAdd && extract_linear_term(exp, var, &b_coeff)) {
                SymbolicExpression left(base.node_->left);
                SymbolicExpression right(base.node_->right);

                auto check_one_plus_a_over_n = [&](const SymbolicExpression& const_part,
                                                   const SymbolicExpression& frac_part,
                                                   Scalar* a_val) -> bool {
                    Scalar one = Scalar(0.0L);
                    if (!const_part.is_number(&one) || !mymath::is_near_zero(one - Scalar(1.0L), Scalar(1e-9L))) {
                        return false;
                    }
                    if (frac_part.node_ && frac_part.node_->type == NodeType::kDivide) {
                        SymbolicExpression rnum(frac_part.node_->left);
                        SymbolicExpression rden(frac_part.node_->right);
                        if (rden.node_->type == NodeType::kVariable && rden.node_->text == var) {
                            return rnum.is_number(a_val);
                        }
                    }
                    return false;
                };

                Scalar a_val = Scalar(0.0L);
                if (check_one_plus_a_over_n(left, right, &a_val) || check_one_plus_a_over_n(right, left, &a_val)) {
                    Scalar exponent_val = a_val * b_coeff;
                    if (mymath::is_near_zero(exponent_val - Scalar(1.0L), Scalar(1e-9L))) {
                        *result = LimitResult::elementary(SymbolicExpression::parse("e"), "known_pattern: (1+1/n)^n");
                    } else {
                        *result = LimitResult::elementary(
                            make_function("exp", SymbolicExpression::number(exponent_val)),
                            "known_pattern: (1+a/n)^(b*n)");
                    }
                    return true;
                }
            }
        }
    }

    // 模式 3：x^x → 1 (x → 0+)
    if (point.is_finite() && mymath::is_near_zero(Scalar(point.value), Scalar(1e-12)) && direction >= 0) {
        if (expr.node_ && expr.node_->type == NodeType::kPower) {
            SymbolicExpression base(expr.node_->left);
            SymbolicExpression exp(expr.node_->right);

            if (base.is_variable_named(var) && exp.is_variable_named(var)) {
                *result = LimitResult::elementary(SymbolicExpression::number(Scalar(1.0L)), "known_pattern: x^x");
                return true;
            }
        }
    }

    return false;
}

std::optional<Scalar> SymbolicLimitEngine::try_direct_substitution(
    const SymbolicExpression& expr,
    const std::string& var,
    Scalar point) {

    return evaluate_at_point(expr, var, point);
}

bool SymbolicLimitEngine::try_rational_limit(
    const SymbolicExpression& num,
    const SymbolicExpression& den,
    const std::string& var,
    const BoundArgument& point,
    int direction,
    LimitResult* result) {

    static thread_local int rational_depth = 0;
    if (rational_depth >= 6) return false;
    ++rational_depth;
    struct RationalDepthGuard { ~RationalDepthGuard() { --rational_depth; } } rat_guard;

    std::vector<Scalar> num_coeffs, den_coeffs;
    if (!num.polynomial_coefficients(var, &num_coeffs) || !den.polynomial_coefficients(var, &den_coeffs)) {
        return false;
    }

    int deg_num = static_cast<int>(num_coeffs.size()) - 1;
    while (deg_num >= 0 && mymath::is_near_zero(num_coeffs[deg_num], Scalar(1e-15L))) deg_num--;

    int deg_den = static_cast<int>(den_coeffs.size()) - 1;
    while (deg_den >= 0 && mymath::is_near_zero(den_coeffs[deg_den], Scalar(1e-15L))) deg_den--;

    if (deg_den < 0) {
        return false;
    }
    if (deg_num < 0) {
        *result = LimitResult::elementary(SymbolicExpression::number(Scalar(0.0L)), "rational_zero_numerator");
        return true;
    }

    if (point.is_infinite()) {
        Scalar c_num = num_coeffs[deg_num];
        Scalar c_den = den_coeffs[deg_den];
        Scalar ratio = c_num / c_den;

        if (deg_num < deg_den) {
            *result = LimitResult::elementary(SymbolicExpression::number(Scalar(0.0L)), "rational_infinity_deg_less");
            return true;
        } else if (deg_num == deg_den) {
            *result = LimitResult::elementary(SymbolicExpression::number(ratio), "rational_infinity_deg_equal");
            return true;
        } else {
            int deg_diff = deg_num - deg_den;
            bool is_pos = (ratio > Scalar(0.0L));
            if (point.is_neg_inf() && deg_diff % 2 != 0) {
                is_pos = !is_pos;
            }
            *result = LimitResult::infinite(is_pos, "rational_infinity_deg_greater");
            return true;
        }
    } else {
        Scalar x0 = point.value;
        auto eval_poly = [](const std::vector<Scalar>& c, int deg, Scalar x) -> Scalar {
            Scalar val = Scalar(0);
            Scalar power = Scalar(1);
            for (int i = 0; i <= deg; ++i) {
                val += c[i] * power;
                power *= x;
            }
            return val;
        };

        Scalar num_val = eval_poly(num_coeffs, deg_num, x0);
        Scalar den_val = eval_poly(den_coeffs, deg_den, x0);

        if (!mymath::is_near_zero(den_val, Scalar(1e-12L))) {
            *result = LimitResult::elementary(SymbolicExpression::number(num_val / den_val), "rational_finite_regular");
            return true;
        }

        if (mymath::is_near_zero(num_val, Scalar(1e-12L))) {
            SymbolicExpression num_d = num.derivative(var).simplify();
            SymbolicExpression den_d = den.derivative(var).simplify();
            return try_rational_limit(num_d, den_d, var, point, direction, result);
        } else {
            int order = 0;
            SymbolicExpression cur_den = den;
            while (order < 10) {
                order++;
                cur_den = cur_den.derivative(var).simplify();
                std::vector<Scalar> d_coeffs;
                if (cur_den.polynomial_coefficients(var, &d_coeffs)) {
                    int d_deg = static_cast<int>(d_coeffs.size()) - 1;
                    while (d_deg >= 0 && mymath::is_near_zero(d_coeffs[d_deg], Scalar(1e-15L))) d_deg--;
                    Scalar d_val = eval_poly(d_coeffs, d_deg, x0);
                    if (!mymath::is_near_zero(d_val, Scalar(1e-12L))) {
                        bool sign_pos = (num_val / d_val) > Scalar(0);
                        if (direction == 0) {
                            if (order % 2 == 0) {
                                *result = LimitResult::infinite(sign_pos, "rational_even_pole");
                                return true;
                            } else {
                                *result = LimitResult::does_not_exist("left and right limits differ (odd pole)");
                                return true;
                            }
                        } else if (direction == 1) {
                            *result = LimitResult::infinite(sign_pos, "rational_right_pole");
                            return true;
                        } else {
                            if (order % 2 != 0) sign_pos = !sign_pos;
                            *result = LimitResult::infinite(sign_pos, "rational_left_pole");
                            return true;
                        }
                    }
                } else {
                    break;
                }
            }
        }
    }

    return false;
}

bool SymbolicLimitEngine::try_conjugate_rationalization(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point,
    int direction,
    LimitResult* result) {

    static thread_local int conjugate_depth = 0;
    if (conjugate_depth >= 3) return false;
    ++conjugate_depth;
    struct DepthGuard { ~DepthGuard() { --conjugate_depth; } } guard;

    if (expr.node_ && expr.node_->type == NodeType::kSubtract) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        SymbolicExpression u, v;
        bool left_is_sqrt = is_sqrt_node(left, &u);
        bool right_is_sqrt = is_sqrt_node(right, &v);

        if (left_is_sqrt || right_is_sqrt) {
            SymbolicExpression new_num;
            SymbolicExpression new_den;

            if (left_is_sqrt && right_is_sqrt) {
                new_num = make_subtract(u, v);
                new_den = make_add(left, right);
            } else if (left_is_sqrt && !right_is_sqrt) {
                new_num = make_subtract(u, make_power(right, SymbolicExpression::number(Scalar(2.0L))));
                new_den = make_add(left, right);
            } else {
                new_num = make_subtract(make_power(left, SymbolicExpression::number(Scalar(2.0L))), v);
                new_den = make_add(left, right);
            }

            SymbolicExpression transformed = make_divide(new_num.simplify(), new_den.simplify()).simplify();
            LimitResult res = compute_limit(transformed, var, point, direction);
            if (res.is_definite) {
                res.method_used = "conjugate_rationalization";
                *result = res;
                return true;
            }
        }
    }

    if (expr.node_ && expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        if (num.node_ && num.node_->type == NodeType::kSubtract) {
            SymbolicExpression left(num.node_->left);
            SymbolicExpression right(num.node_->right);

            SymbolicExpression u, v;
            bool left_is_sqrt = is_sqrt_node(left, &u);
            bool right_is_sqrt = is_sqrt_node(right, &v);

            if (left_is_sqrt || right_is_sqrt) {
                SymbolicExpression new_num;
                SymbolicExpression conj_part;

                if (left_is_sqrt && right_is_sqrt) {
                    new_num = make_subtract(u, v);
                    conj_part = make_add(left, right);
                } else if (left_is_sqrt && !right_is_sqrt) {
                    new_num = make_subtract(u, make_power(right, SymbolicExpression::number(Scalar(2.0L))));
                    conj_part = make_add(left, right);
                } else {
                    new_num = make_subtract(make_power(left, SymbolicExpression::number(Scalar(2.0L))), v);
                    conj_part = make_add(left, right);
                }

                SymbolicExpression new_den = make_multiply(den, conj_part);
                SymbolicExpression transformed = make_divide(new_num.simplify(), new_den.simplify()).simplify();
                LimitResult res = compute_limit(transformed, var, point, direction);
                if (res.is_definite) {
                    res.method_used = "fraction_conjugate_rationalization";
                    *result = res;
                    return true;
                }
            }
        }
    }

    return false;
}

bool SymbolicLimitEngine::apply_series_expansion(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point,
    int direction,
    LimitResult* result) {

    series_ops::SeriesContext ctx;
    ctx.evaluate_at = [](const SymbolicExpression& e, const std::string&, Scalar) {
        Scalar val = 0.0L;
        if (e.is_number(&val)) return val;
        if (e.node_type() == NodeType::kPi) return mymath::pi();
        if (e.node_type() == NodeType::kE) return mymath::e();
        return Scalar(0.0L);
    };

    if (point.is_finite()) {
        std::vector<Scalar> coeffs;
        try {
            if (series_ops::internal::evaluate_psa(expr, var, point.value, 4, coeffs, ctx)) {
                if (!coeffs.empty() && mymath::isfinite(coeffs[0])) {
                    *result = LimitResult::elementary(SymbolicExpression::number(coeffs[0]), "psa_taylor_series");
                    return true;
                }
            }
        } catch (...) {
            return false;
        }
    } else {
        SymbolicExpression t_var = SymbolicExpression::variable("__t_psa_inf");
        SymbolicExpression inv_t = point.is_pos_inf()
            ? SymbolicExpression::number(Scalar(1.0L)) / t_var
            : SymbolicExpression::number(Scalar(-1.0L)) / t_var;
        SymbolicExpression substituted = expr.substitute(var, inv_t).simplify();

        std::vector<Scalar> coeffs;
        try {
            if (series_ops::internal::evaluate_psa(substituted, "__t_psa_inf", Scalar(0), 4, coeffs, ctx)) {
                if (!coeffs.empty() && mymath::isfinite(coeffs[0])) {
                    *result = LimitResult::elementary(SymbolicExpression::number(coeffs[0]), "psa_at_infinity");
                    return true;
                }
            }
        } catch (...) {
            return false;
        }
    }

    return false;
}

bool SymbolicLimitEngine::apply_squeeze_or_bounded(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point,
    int /*direction*/,
    LimitResult* result) {

    if (point.is_infinite()) {
        if (expr.node_ && expr.node_->type == NodeType::kFunction) {
            if (expr.node_->text == "sin" || expr.node_->text == "cos") {
                SymbolicExpression arg(expr.node_->left);
                Scalar k = Scalar(0);
                if (extract_linear_term(arg, var, &k) && !mymath::is_near_zero(k, Scalar(1e-9L))) {
                    *result = LimitResult::does_not_exist("oscillates at infinity");
                    return true;
                }
            }
        }
    }

    if (expr.node_ && expr.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(expr.node_->left);
        SymbolicExpression right(expr.node_->right);

        auto is_bounded_oscillating = [&](const SymbolicExpression& e) -> bool {
            if (e.node_ && e.node_->type == NodeType::kFunction) {
                if (e.node_->text == "sin" || e.node_->text == "cos") {
                    return true;
                }
            }
            return false;
        };

        if (is_bounded_oscillating(right) && is_zero_at_point(left, var, point)) {
            *result = LimitResult::elementary(SymbolicExpression::number(Scalar(0.0L)), "squeeze_bounded_times_infinitesimal");
            return true;
        }
        if (is_bounded_oscillating(left) && is_zero_at_point(right, var, point)) {
            *result = LimitResult::elementary(SymbolicExpression::number(Scalar(0.0L)), "squeeze_bounded_times_infinitesimal");
            return true;
        }
    }

    if (point.is_infinite() && expr.node_ && expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);
        if (num.node_ && num.node_->type == NodeType::kFunction &&
            (num.node_->text == "sin" || num.node_->text == "cos") &&
            is_infinite_at_point(den, var, point)) {
            *result = LimitResult::elementary(SymbolicExpression::number(Scalar(0.0L)), "squeeze_sin_over_infinity");
            return true;
        }
    }

    return false;
}

LimitResult SymbolicLimitEngine::limit_at_infinity(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point,
    int direction) {

    LimitResult squeeze_res;
    if (apply_squeeze_or_bounded(expr, var, point, direction, &squeeze_res)) {
        return squeeze_res;
    }

    LimitResult pattern_res;
    if (try_known_pattern(expr, var, point, direction, &pattern_res)) {
        return pattern_res;
    }

    SymbolicExpression num, den;
    extract_numerator_denominator(expr, &num, &den);
    LimitResult rational_res;
    if (try_rational_limit(num, den, var, point, direction, &rational_res)) {
        return rational_res;
    }

    LimitResult conjugate_res;
    if (try_conjugate_rationalization(expr, var, point, direction, &conjugate_res)) {
        return conjugate_res;
    }

    // 5. 对数分式: ln(x) / x^p -> 0 (x -> +inf, p > 0)
    if (expr.node_ && expr.node_->type == NodeType::kDivide) {
        SymbolicExpression n(expr.node_->left);
        SymbolicExpression d(expr.node_->right);
        if (n.node_ && n.node_->type == NodeType::kFunction && n.node_->text == "ln" && point.is_pos_inf()) {
            if (d.is_variable_named(var)) {
                return LimitResult::elementary(SymbolicExpression::number(Scalar(0.0L)), "ln_over_x_at_infinity");
            }
            if (d.node_ && d.node_->type == NodeType::kPower &&
                SymbolicExpression(d.node_->left).is_variable_named(var)) {
                Scalar p = Scalar(0);
                if (SymbolicExpression(d.node_->right).is_number(&p) && p > Scalar(0)) {
                    return LimitResult::elementary(SymbolicExpression::number(Scalar(0.0L)), "ln_over_x_p_at_infinity");
                }
            }
        }
    }

    // 6. 指数函数特征分析
    if (expr.node_ && expr.node_->type == NodeType::kFunction && expr.node_->text == "exp") {
        SymbolicExpression arg(expr.node_->left);
        const bool at_negative_infinity = point.is_neg_inf();
        if (arg.node_ && arg.node_->type == NodeType::kNegate) {
            SymbolicExpression inner(arg.node_->left);
            if (inner.node_ && inner.node_->type == NodeType::kVariable && inner.node_->text == var) {
                return at_negative_infinity
                    ? LimitResult::infinite(true, "exp(-x)_at_negative_infinity")
                    : LimitResult::elementary(SymbolicExpression::number(0.0L), "exp(-x)_at_infinity");
            }
        }
        if (arg.node_ && arg.node_->type == NodeType::kVariable && arg.node_->text == var) {
            return at_negative_infinity
                ? LimitResult::elementary(SymbolicExpression::number(0.0L), "exp(x)_at_negative_infinity")
                : LimitResult::infinite(true, "exp(x)_at_infinity");
        }
    }

    // 7. PSA 展开 (t = 1/x -> 0+)
    LimitResult series_res;
    if (apply_series_expansion(expr, var, point, direction, &series_res)) {
        return series_res;
    }

    // 8. Gruntz 后备探针
    LimitResult gruntz_res;
    if (apply_gruntz(expr, var, point, direction, &gruntz_res)) {
        return gruntz_res;
    }

    return LimitResult::unknown();
}

bool SymbolicLimitEngine::handle_one_to_infinity(
    const SymbolicExpression& base,
    const SymbolicExpression& exponent,
    const std::string& var,
    const BoundArgument& point,
    int direction,
    LimitResult* result) {

    static thread_local int one_to_inf_depth = 0;
    if (one_to_inf_depth >= 3) return false;
    ++one_to_inf_depth;
    struct DepthGuard { ~DepthGuard() { --one_to_inf_depth; } } guard;

    SymbolicExpression base_minus_one = make_subtract(base, SymbolicExpression::number(Scalar(1.0L)));
    SymbolicExpression product = make_multiply(base_minus_one, exponent);

    LimitResult product_limit = compute_limit(product.simplify(), var, point, direction);

    if (product_limit.is_definite && product_limit.is_elementary) {
        Scalar p_val = Scalar(0);
        if (product_limit.value.is_number(&p_val)) {
            if (mymath::is_near_zero(p_val - Scalar(1.0L), Scalar(1e-9L))) {
                *result = LimitResult::elementary(SymbolicExpression::parse("e"), "one_to_infinity_e");
                return true;
            }
        }
        SymbolicExpression final_result = make_function("exp", product_limit.value).simplify();
        *result = LimitResult::elementary(final_result, "one_to_infinity_transform");
        return true;
    }

    return false;
}

bool SymbolicLimitEngine::handle_power_indeterminate(
    const SymbolicExpression& base,
    const SymbolicExpression& exponent,
    const std::string& var,
    const BoundArgument& point,
    int direction,
    LimitResult* result) {

    static thread_local int power_depth = 0;
    if (power_depth >= 3) return false;
    ++power_depth;
    struct DepthGuard { ~DepthGuard() { --power_depth; } } guard;

    SymbolicExpression transformed_exp = make_multiply(exponent, make_function("ln", base));
    LimitResult exp_limit = compute_limit(transformed_exp.simplify(), var, point, direction);

    if (exp_limit.is_definite && exp_limit.is_elementary) {
        Scalar val = Scalar(0);
        if (exp_limit.value.is_number(&val)) {
            *result = LimitResult::elementary(SymbolicExpression::number(mymath::exp(val)), "power_indeterminate_transform");
            return true;
        }
        SymbolicExpression final_result = make_function("exp", exp_limit.value).simplify();
        *result = LimitResult::elementary(final_result, "power_indeterminate_transform");
        return true;
    }

    return false;
}

bool SymbolicLimitEngine::is_infinite_at_point(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point) {

    if (expr.node_ == nullptr) return false;

    if (expr.node_->type == NodeType::kInfinity) {
        return true;
    }

    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        return point.is_infinite();
    }

    if (expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);

        if (is_zero_at_point(den, var, point) && !is_zero_at_point(num, var, point)) {
            return true;
        }
    }

    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "ln") {
        SymbolicExpression arg(expr.node_->left);
        if (is_zero_at_point(arg, var, point)) {
            return true;
        }
    }

    if (point.is_infinite()) {
        std::vector<Scalar> coeffs;
        if (expr.polynomial_coefficients(var, &coeffs)) {
            int deg = static_cast<int>(coeffs.size()) - 1;
            while (deg >= 0 && mymath::is_near_zero(coeffs[deg], Scalar(1e-15L))) deg--;
            if (deg > 0) return true;
        }
    }

    return false;
}

bool SymbolicLimitEngine::is_zero_at_point(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point) {

    if (expr.node_ == nullptr) return false;

    if (point.is_finite()) {
        const auto value = evaluate_at_point(expr, var, point.value);
        if (value.has_value() && mymath::is_near_zero(*value, Scalar(1e-10L))) {
            return true;
        }
    }

    Scalar val = Scalar(0);
    if (expr.is_number(&val) && mymath::is_near_zero(val, Scalar(1e-15L))) {
        return true;
    }

    if (expr.node_->type == NodeType::kVariable && expr.node_->text == var) {
        if (point.is_finite() && mymath::is_near_zero(Scalar(point.value), Scalar(1e-15L))) {
            return true;
        }
    }

    if (point.is_infinite() && expr.node_->type == NodeType::kDivide) {
        SymbolicExpression num(expr.node_->left);
        SymbolicExpression den(expr.node_->right);
        if (!is_infinite_at_point(num, var, point) && is_infinite_at_point(den, var, point)) {
            return true;
        }
    }

    return false;
}

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

bool SymbolicLimitEngine::apply_gruntz(
    const SymbolicExpression& expr,
    const std::string& var,
    const BoundArgument& point,
    int direction,
    LimitResult* result) {

    std::optional<Scalar> val1, val2, val3;
    if (point.is_infinite()) {
        const Scalar sign = (direction < 0 || point.is_neg_inf()) ? Scalar(-1.0L) : Scalar(1.0L);
        val1 = evaluate_at_point(expr, var, sign * Scalar(1e4L));
        val2 = evaluate_at_point(expr, var, sign * Scalar(1e6L));
        val3 = evaluate_at_point(expr, var, sign * Scalar(1e8L));
    } else {
        const Scalar sign = (direction < 0) ? Scalar(-1.0L) : Scalar(1.0L);
        val1 = evaluate_at_point(expr, var, point.value + sign * Scalar(1e-4L));
        val2 = evaluate_at_point(expr, var, point.value + sign * Scalar(1e-6L));
        val3 = evaluate_at_point(expr, var, point.value + sign * Scalar(1e-8L));
    }

    if (val1.has_value() && val2.has_value() && val3.has_value()) {
        Scalar v1 = *val1, v2 = *val2, v3 = *val3;
        if (mymath::isfinite(v1) && mymath::isfinite(v2) && mymath::isfinite(v3)) {
            if (mymath::abs(v3 - v2) < Scalar(1e-4L) * (Scalar(1) + mymath::abs(v3))) {
                // Snap to simple rationals if very close
                Scalar snapped = v3;
                if (mymath::abs(v3) < Scalar(1e-5L)) {
                    snapped = Scalar(0.0L);
                } else {
                    for (int den = 1; den <= 12; ++den) {
                        Scalar num = mymath::round(v3 * Scalar(den));
                        if (mymath::abs(v3 - num / Scalar(den)) < Scalar(1e-4L)) {
                            snapped = num / Scalar(den);
                            break;
                        }
                    }
                }
                *result = LimitResult::elementary(SymbolicExpression::number(snapped), "gruntz_dominant_scale");
                return true;
            }
            if (mymath::abs(v3) < Scalar(1e-8L) && mymath::abs(v3) <= mymath::abs(v2)) {
                *result = LimitResult::elementary(SymbolicExpression::number(0.0L), "gruntz_infinitesimal");
                return true;
            }
            if (mymath::abs(v3) > Scalar(1e6L) && mymath::abs(v3) > mymath::abs(v2)) {
                bool is_pos = (v3 > Scalar(0));
                *result = LimitResult::infinite(is_pos, "gruntz_divergent");
                return true;
            }
        }
    }

    return false;
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

        std::string point_str = args[2];
        if (point_str == "inf" || point_str == "infinity" || point_str == "oo" || point_str == "+inf") {
            *point = BoundArgument::pos_inf();
        } else if (point_str == "-inf" || point_str == "-infinity" || point_str == "-oo") {
            *point = BoundArgument::neg_inf();
        } else {
            Scalar p = Scalar(point_str);
            *point = BoundArgument::finite(p);
        }

        *direction = 0;
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


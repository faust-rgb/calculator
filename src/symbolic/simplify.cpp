// ============================================================================
// 符号表达式简化模块
// ============================================================================
//
// 本文件实现符号表达式的代数简化规则。简化是符号计算的核心，应用于：
// - 解析后的表达式规范化
// - 微积分运算后的结果整理
// - 用户交互时的输出优化
//
// 简化策略采用启发式规则集，按表达式类型分派：
// 1. 常数折叠：数值运算的直接求值
// 2. 恒等式消去：x+0→x, x*1→x, x^0→1 等
// 3. 三角恒等式：sin²(x)+cos²(x)→1, sin(0)→0 等
// 4. 幂运算简化：x^a * x^b → x^(a+b), (x^a)^b → x^(a*b)
// 5. 有理式约分：多项式 GCD 约分、整除化简
// 6. 公因子提取：2x+2y → 2(x+y)
// 7. 多项式规范化：单变量多项式按降幂排列
//
// 简化通过多轮迭代进行，直到表达式结构不再变化（最多 16 轮）。
// 结果通过 LRU 缓存记忆，避免重复计算。
// ============================================================================

#include "symbolic_expression_internal.h"

#include "mymath.h"

#include <string>
#include <vector>

namespace symbolic_expression_internal {

// ============================================================================
// 单轮简化主函数
// ============================================================================

/**
 * @brief 对表达式应用一轮简化规则
 *
 * 这是简化的核心调度函数，根据节点类型分派到相应的规则：
 * - kNumber/kVariable: 原样返回
 * - kFunction: 简化参数，应用函数特定规则
 * - kNegate: 简化操作数，处理双重取负
 * - 二元运算: 简化两操作数，应用相应代数规则
 *
 * 注意：此函数只应用规则一次，不保证达到最优形式。
 * 完整简化通过 simplify_impl() 多轮调用此函数实现。
 */
SymbolicExpression simplify_once(const SymbolicExpression& expression) {
    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
        case NodeType::kVariable:
            return expression;
        // ========================================================================
    // 函数节点简化
    // ========================================================================
    // 策略：先简化参数，然后应用函数特定的规则：
    // 1. 数值参数：直接调用数值函数求值
    // 2. 特殊恒等式：ln(e)→1, exp(ln(x))→x（需正性检查）, sqrt(x²)→abs(x)
    // 3. 奇偶性处理：sin(-x)→-sin(x), cos(-x)→cos(x)
    // 4. π 的特殊角度：sin(π/6)→1/2, cos(π/3)→1/2 等
    case NodeType::kFunction: {
            const SymbolicExpression argument = SymbolicExpression(node->left).simplify();
            double numeric = 0.0;
            if (argument.is_number(&numeric)) {
                if (node->text == "asin") {
                    return SymbolicExpression::number(mymath::asin(numeric));
                }
                if (node->text == "acos") {
                    return SymbolicExpression::number(mymath::acos(numeric));
                }
                if (node->text == "atan") {
                    return SymbolicExpression::number(mymath::atan(numeric));
                }
                if (node->text == "sin") {
                    return SymbolicExpression::number(mymath::sin(numeric));
                }
                if (node->text == "cos") {
                    return SymbolicExpression::number(mymath::cos(numeric));
                }
                if (node->text == "tan") {
                    return SymbolicExpression::number(mymath::tan(numeric));
                }
                if (node->text == "exp") {
                    return SymbolicExpression::number(mymath::exp(numeric));
                }
                if (node->text == "sinh") {
                    return SymbolicExpression::number(mymath::sinh(numeric));
                }
                if (node->text == "cosh") {
                    return SymbolicExpression::number(mymath::cosh(numeric));
                }
                if (node->text == "tanh") {
                    return SymbolicExpression::number(mymath::tanh(numeric));
                }
                if (node->text == "ln") {
                    return SymbolicExpression::number(mymath::ln(numeric));
                }
                if (node->text == "sqrt") {
                    return SymbolicExpression::number(mymath::sqrt(numeric));
                }
                if (node->text == "abs") {
                    return SymbolicExpression::number(mymath::abs(numeric));
                }
                if (node->text == "floor") {
                    return SymbolicExpression::number(static_cast<double>(
                        static_cast<long long>(numeric < 0.0 && static_cast<double>(static_cast<long long>(numeric)) != numeric
                                                   ? numeric - 1.0
                                                   : numeric)));
                }
                if (node->text == "ceil") {
                    long long truncated = static_cast<long long>(numeric);
                    if (numeric > 0.0 && static_cast<double>(truncated) != numeric) {
                        ++truncated;
                    }
                    return SymbolicExpression::number(static_cast<double>(truncated));
                }
                if (node->text == "cbrt") {
                    return SymbolicExpression::number(mymath::cbrt(numeric));
                }
                if (node->text == "sign") {
                    if (mymath::is_near_zero(numeric, kFormatEps)) {
                        return SymbolicExpression::number(0.0);
                    }
                    return SymbolicExpression::number(numeric > 0.0 ? 1.0 : -1.0);
                }
            }

            if (node->text == "ln" && expr_is_variable(argument, "e")) {
                return SymbolicExpression::number(1.0);
            }
            if (node->text == "exp" &&
                argument.node_->type == NodeType::kFunction &&
                argument.node_->text == "ln") {
                const SymbolicExpression inner(argument.node_->left);
                if (is_known_positive_expression(inner.simplify())) {
                    return inner.simplify();
                }
                return make_function(node->text, argument);
            }
            if (node->text == "ln" &&
                argument.node_->type == NodeType::kFunction &&
                argument.node_->text == "exp") {
                return SymbolicExpression(argument.node_->left).simplify();
            }
            if (node->text == "sqrt" &&
                argument.node_->type == NodeType::kPower) {
                double exponent = 0.0;
                if (SymbolicExpression(argument.node_->right).is_number(&exponent) &&
                    mymath::is_near_zero(exponent - 2.0, kFormatEps)) {
                    return make_function("abs",
                                         SymbolicExpression(argument.node_->left))
                        .simplify();
                }
            }
            if (node->text == "abs" &&
                argument.node_->type == NodeType::kFunction &&
                (argument.node_->text == "abs" || argument.node_->text == "sqrt")) {
                return argument;
            }
            if ((node->text == "sin" || node->text == "tan" ||
                 node->text == "sinh" || node->text == "tanh") &&
                argument.node_->type == NodeType::kNegate) {
                return make_negate(
                           make_function(node->text,
                                         SymbolicExpression(argument.node_->left)))
                    .simplify();
            }
            if ((node->text == "cos" || node->text == "cosh" ||
                 node->text == "abs") &&
                argument.node_->type == NodeType::kNegate) {
                return make_function(node->text,
                                     SymbolicExpression(argument.node_->left))
                    .simplify();
            }

            double pi_multiple = 0.0;
            if (decompose_numeric_multiple_of_symbol(argument, "pi", &pi_multiple)) {
                if (node->text == "sin") {
                    if (numeric_matches_any(pi_multiple, {0.0, 1.0, -1.0})) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (numeric_matches_any(pi_multiple, {0.5})) {
                        return SymbolicExpression::number(1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {-0.5})) {
                        return SymbolicExpression::number(-1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, 5.0 / 6.0})) {
                        return half_symbol();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 6.0, -5.0 / 6.0})) {
                        return make_negate(half_symbol()).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, 2.0 / 3.0})) {
                        return make_divide(sqrt3_symbol(),
                                           SymbolicExpression::number(2.0)).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 3.0, -2.0 / 3.0})) {
                        return make_negate(
                            make_divide(sqrt3_symbol(),
                                        SymbolicExpression::number(2.0))).simplify();
                    }
                }
                if (node->text == "cos") {
                    if (numeric_matches_any(pi_multiple, {0.0})) {
                        return SymbolicExpression::number(1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0, -1.0})) {
                        return SymbolicExpression::number(-1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {0.5, -0.5})) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, -1.0 / 3.0})) {
                        return half_symbol();
                    }
                    if (numeric_matches_any(pi_multiple, {2.0 / 3.0, -2.0 / 3.0})) {
                        return make_negate(half_symbol()).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, -1.0 / 6.0})) {
                        return make_divide(sqrt3_symbol(),
                                           SymbolicExpression::number(2.0)).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {5.0 / 6.0, -5.0 / 6.0})) {
                        return make_negate(
                            make_divide(sqrt3_symbol(),
                                        SymbolicExpression::number(2.0))).simplify();
                    }
                }
                if (node->text == "tan") {
                    if (numeric_matches_any(pi_multiple, {0.0, 1.0, -1.0})) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (numeric_matches_any(pi_multiple, {0.25})) {
                        return SymbolicExpression::number(1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {-0.25})) {
                        return SymbolicExpression::number(-1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, -5.0 / 6.0})) {
                        return make_divide(SymbolicExpression::number(1.0),
                                           sqrt3_symbol()).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 6.0, 5.0 / 6.0})) {
                        return make_negate(
                            make_divide(SymbolicExpression::number(1.0),
                                        sqrt3_symbol())).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, -2.0 / 3.0})) {
                        return sqrt3_symbol();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 3.0, 2.0 / 3.0})) {
                        return make_negate(sqrt3_symbol()).simplify();
                    }
                }
            }

            return make_function(node->text, argument);
        }
        // ========================================================================
    // 取负节点简化
    // ========================================================================
    // 规则：
    // 1. 数值取负：直接计算
    // 2. 双重取负：-(-x) → x
    case NodeType::kNegate: {
            const SymbolicExpression operand = SymbolicExpression(node->left).simplify();
            double value = 0.0;
            if (operand.is_number(&value)) {
                return SymbolicExpression::number(-value);
            }
            if (operand.node_->type == NodeType::kNegate) {
                return SymbolicExpression(operand.node_->left).simplify();
            }
            return make_negate(operand);
        }
        // ========================================================================
    // 二元运算节点简化
    // ========================================================================
    // 首先简化左右操作数，然后根据运算类型应用规则
    case NodeType::kAdd:
    case NodeType::kSubtract:
    case NodeType::kMultiply:
    case NodeType::kDivide:
    case NodeType::kPower:
        break;
    }

    // 简化后的左右操作数
    const SymbolicExpression left = SymbolicExpression(node->left).simplify();
    const SymbolicExpression right = SymbolicExpression(node->right).simplify();
    double left_value = 0.0;
    double right_value = 0.0;

    switch (node->type) {
        // ====================================================================
        // 加法简化
        // ====================================================================
        // 规则优先级：
        // 1. 常数折叠：数值直接相加
        // 2. 零元消去：x+0 → x
        // 3. 三角恒等式：sin²(x)+cos²(x) → 1
        // 4. 相似项合并：2x+3x → 5x
        // 5. 多项式规范化
        // 6. 公因子提取：ax+ay → a(x+y)
        // 7. 项排序：按结构键字典序排列
        case NodeType::kAdd:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value + right_value);
            }
            if (expr_is_zero(left)) {
                return right;
            }
            if (expr_is_zero(right)) {
                return left;
            }
            {
                std::string left_argument;
                std::string right_argument;
                if ((is_squared_function(left, "sin", &left_argument) &&
                     is_squared_function(right, "cos", &right_argument) &&
                     left_argument == right_argument) ||
                    (is_squared_function(left, "cos", &left_argument) &&
                     is_squared_function(right, "sin", &right_argument) &&
                     left_argument == right_argument)) {
                    return SymbolicExpression::number(1.0);
                }
            }
            {
                SymbolicExpression combined;
                if (try_combine_like_terms(left, right, 1.0, &combined)) {
                    return combined;
                }
            }
            {
                const SymbolicExpression sum = make_add(left, right);
                if (is_single_variable_polynomial(sum)) {
                    return maybe_canonicalize_polynomial(sum);
                }
            }
            {
                SymbolicExpression factored;
                if (try_factor_common_terms(left, right, 1.0, &factored)) {
                    return factored;
                }
            }
            {
                std::vector<SymbolicExpression> terms;
                collect_additive_expressions(make_add(left, right), &terms);
                const SymbolicExpression sorted = make_sorted_sum(terms);
                if (node_structural_key(sorted.node_) !=
                    node_structural_key(make_add(left, right).node_)) {
                    return sorted;
                }
            }
            return make_add(left, right);
        // ====================================================================
        // 减法简化
        // ====================================================================
        // 规则：
        // 1. 常数折叠
        // 2. 零元消去：x-0 → x, 0-x → -x
        // 3. 负负得正：x-(-y) → x+y
        // 4. 相似项合并
        // 5. 多项式规范化
        // 6. 公因子提取
        case NodeType::kSubtract:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value - right_value);
            }
            if (expr_is_zero(right)) {
                return left;
            }
            if (expr_is_zero(left)) {
                return make_negate(right).simplify();
            }
            if (right.node_->type == NodeType::kNegate) {
                return make_add(left, SymbolicExpression(right.node_->left)).simplify();
            }
            if (right.is_number(&right_value) && right_value < 0.0) {
                return make_add(left, SymbolicExpression::number(-right_value)).simplify();
            }
            {
                SymbolicExpression combined;
                if (try_combine_like_terms(left, right, -1.0, &combined)) {
                    return combined;
                }
            }
            {
                const SymbolicExpression difference = make_subtract(left, right);
                if (is_single_variable_polynomial(difference)) {
                    return maybe_canonicalize_polynomial(difference);
                }
            }
            {
                SymbolicExpression factored;
                if (try_factor_common_terms(left, right, -1.0, &factored)) {
                    return factored;
                }
            }
            return make_subtract(left, right);
        // ====================================================================
        // 乘法简化
        // ====================================================================
        // 规则：
        // 1. 常数折叠
        // 2. 零元消去：x*0 → 0
        // 3. 单位元消去：x*1 → x
        // 4. 负号处理：x*(-1) → -x
        // 5. 幂次合并：x^a * x^b → x^(a+b)
        // 6. 乘积规范化：收集因子、排序、合并幂次
        case NodeType::kMultiply:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value * right_value);
            }
            if (expr_is_zero(left) || expr_is_zero(right)) {
                return SymbolicExpression::number(0.0);
            }
            if (expr_is_one(left)) {
                return right;
            }
            if (expr_is_one(right)) {
                return left;
            }
            if (expr_is_minus_one(left)) {
                return make_negate(right).simplify();
            }
            if (expr_is_minus_one(right)) {
                return make_negate(left).simplify();
            }
            {
                SymbolicExpression left_base;
                SymbolicExpression right_base;
                double left_exponent = 0.0;
                double right_exponent = 0.0;
                decompose_power_factor(left, &left_base, &left_exponent);
                decompose_power_factor(right, &right_base, &right_exponent);
                if (expressions_match(left_base, right_base)) {
                    return rebuild_power_difference(left_base,
                                                    left_exponent + right_exponent);
                }
            }
            {
                double numeric_factor = 1.0;
                std::vector<SymbolicExpression> symbolic_factors;
                collect_multiplicative_terms(left, &numeric_factor, &symbolic_factors);
                collect_multiplicative_terms(right, &numeric_factor, &symbolic_factors);
                return make_sorted_product(numeric_factor, symbolic_factors);
            }
            // ====================================================================
            // 除法简化
            // ====================================================================
            // 规则：
            // 1. 常数折叠
            // 2. 零分子：0/x → 0
            // 3. 单位分母：x/1 → x
            // 4. 负号提取：(-a)/b → -(a/b)
            // 5. 多项式整除约分
            // 6. 多项式 GCD 约分
            // 7. 幂次约分：x^a / x^b → x^(a-b)
            // 8. 因子约分：展开幂次后对消相同因子
            case NodeType::kDivide:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value / right_value);
            }
            if (expr_is_zero(left)) {
                return SymbolicExpression::number(0.0);
            }
            if (expr_is_one(right)) {
                return left;
            }
            if (left.is_number(&left_value) && left_value < 0.0) {
                return make_negate(make_divide(SymbolicExpression::number(-left_value),
                                               right));
            }
            {
                SymbolicExpression reduced;
                if (try_reduce_polynomial_quotient(left, right, &reduced)) {
                    return reduced;
                }
            }
            {
                SymbolicExpression reduced;
                if (try_reduce_polynomial_gcd_quotient(left, right, &reduced)) {
                    return reduced;
                }
            }
            {
                SymbolicExpression left_base;
                SymbolicExpression right_base;
                double left_exponent = 0.0;
                double right_exponent = 0.0;
                decompose_power_factor(left, &left_base, &left_exponent);
                decompose_power_factor(right, &right_base, &right_exponent);
                if (expressions_match(left_base, right_base)) {
                    return rebuild_power_difference(left_base,
                                                    left_exponent - right_exponent);
                }
            }
            {
                SymbolicExpression quotient;
                if (try_canonical_factor_quotient(left, right, &quotient)) {
                    return quotient;
                }
            }
            {
                double numerator_coefficient = 1.0;
                double denominator_coefficient = 1.0;
                std::vector<SymbolicExpression> numerator_factors;
                std::vector<SymbolicExpression> denominator_factors;
                collect_division_factors(left, &numerator_coefficient, &numerator_factors);
                collect_division_factors(right, &denominator_coefficient, &denominator_factors);
                if (!mymath::is_near_zero(denominator_coefficient, kFormatEps)) {
                    std::vector<bool> denominator_used(denominator_factors.size(), false);
                    std::vector<SymbolicExpression> reduced_numerator_factors;
                    for (const SymbolicExpression& numerator_factor : numerator_factors) {
                        const std::string numerator_key =
                            node_structural_key(numerator_factor.node_);
                        bool canceled = false;
                        for (std::size_t i = 0; i < denominator_factors.size(); ++i) {
                            if (denominator_used[i]) {
                                continue;
                            }
                            if (numerator_key ==
                                node_structural_key(denominator_factors[i].node_)) {
                                denominator_used[i] = true;
                                canceled = true;
                                break;
                            }
                        }
                        if (!canceled) {
                            reduced_numerator_factors.push_back(numerator_factor);
                        }
                    }

                    std::vector<SymbolicExpression> reduced_denominator_factors;
                    for (std::size_t i = 0; i < denominator_factors.size(); ++i) {
                        if (!denominator_used[i]) {
                            reduced_denominator_factors.push_back(denominator_factors[i]);
                        }
                    }

                    const bool symbolic_cancellation_happened =
                        reduced_numerator_factors.size() != numerator_factors.size() ||
                        reduced_denominator_factors.size() != denominator_factors.size();

                    const double reduced_coefficient =
                        numerator_coefficient / denominator_coefficient;
                    if (!symbolic_cancellation_happened &&
                        reduced_denominator_factors.empty() &&
                        mymath::is_near_zero(numerator_coefficient - 1.0, kFormatEps) &&
                        !mymath::is_near_zero(denominator_coefficient - 1.0, kFormatEps)) {
                        return make_divide(left,
                                           SymbolicExpression::number(denominator_coefficient));
                    }
                    SymbolicExpression numerator_expression =
                        rebuild_product_expression(reduced_coefficient,
                                                   reduced_numerator_factors);
                    SymbolicExpression denominator_expression =
                        rebuild_product_expression(1.0,
                                                   reduced_denominator_factors);

                    if (expr_is_zero(numerator_expression)) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (expr_is_one(denominator_expression)) {
                        return numerator_expression;
                    }
                    if (numerator_expression.is_number(&left_value) &&
                        denominator_expression.is_number(&right_value)) {
                        return SymbolicExpression::number(left_value / right_value);
                    }
                    return make_divide(numerator_expression, denominator_expression);
                }
            }
            return make_divide(left, right);
        // ====================================================================
        // 幂运算简化
        // ====================================================================
        // 规则：
        // 1. 底数为零或一：0^x → 0, 1^x → 1
        // 2. 指数为零或一：x^0 → 1, x^1 → x
        // 3. 幂的幂：(x^a)^b → x^(a*b)
        // 4. 常数幂：数值直接计算
        case NodeType::kPower:
            if (left.is_number(&left_value)) {
                if (mymath::is_near_zero(left_value, kFormatEps)) {
                    return SymbolicExpression::number(0.0);
                }
                if (mymath::is_near_zero(left_value - 1.0, kFormatEps)) {
                    return SymbolicExpression::number(1.0);
                }
            }
            if (right.is_number(&right_value)) {
                if (mymath::is_near_zero(right_value, kFormatEps)) {
                    return SymbolicExpression::number(1.0);
                }
                if (mymath::is_near_zero(right_value - 1.0, kFormatEps)) {
                    return left;
                }
            }
            if (left.node_->type == NodeType::kPower) {
                double inner_exponent = 0.0;
                if (SymbolicExpression(left.node_->right).is_number(&inner_exponent) &&
                    right.is_number(&right_value)) {
                    return make_power(SymbolicExpression(left.node_->left).simplify(),
                                      SymbolicExpression::number(inner_exponent * right_value))
                        .simplify();
                }
            }
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(mymath::pow(left_value, right_value));
            }
            return make_power(left, right);
        // 已在上面处理的节点类型，此处仅为消除编译器警告
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kNegate:
        case NodeType::kFunction:
            break;
    }
    return expression;
}

// ============================================================================
// 多轮简化实现
// ============================================================================

/**
 * @brief 简化的完整实现（多轮迭代）
 *
 * 简化通过迭代进行，直到表达式结构不再变化。
 * 最多进行 kMaxSimplifyPasses (16) 轮迭代。
 *
 * 使用结构键检测变化，避免不必要的字符串比较。
 *
 * 注意：递归调用时使用 simplify_once 而非 simplify_impl，
 * 以避免无限循环和重复的缓存查询。
 */
SymbolicExpression simplify_impl(const SymbolicExpression& expression) {
    // 简化深度计数器，用于检测递归调用
    static thread_local int simplify_depth = 0;

    // RAII 守卫，确保深度计数器正确递减
    struct SimplifyDepthGuard {
        int* depth;
        explicit SimplifyDepthGuard(int* value) : depth(value) {
            ++(*depth);
        }
        ~SimplifyDepthGuard() {
            --(*depth);
        }
    };

    // 如果已经在递归调用中，直接执行单轮简化
    // 这避免了嵌套的多轮迭代
    if (simplify_depth > 0) {
        SimplifyDepthGuard guard(&simplify_depth);
        return simplify_once(expression);
    }

    SimplifyDepthGuard guard(&simplify_depth);
    SymbolicExpression current = expression;

    // 最多 16 轮迭代，通常 2-4 轮即可收敛
    constexpr int kMaxSimplifyPasses = 16;
    for (int pass = 0; pass < kMaxSimplifyPasses; ++pass) {
        // 使用结构键检测是否还有变化
        const std::string current_key = node_structural_key(current.node_);
        SymbolicExpression next = simplify_once(current);
        const std::string next_key = node_structural_key(next.node_);

        // 结构键不变表示已收敛
        if (next_key == current_key) {
            return next;
        }
        current = next;
    }

    // 达到最大迭代次数，返回当前结果
    return current;
}


}  // namespace symbolic_expression_internal

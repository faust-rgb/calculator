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

#include "symbolic/symbolic_expression_internal.h"

#include "math/mymath.h"

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>

namespace symbolic_expression_internal {

// ============================================================================
// Expression size monitoring
// ============================================================================

std::size_t count_nodes(const std::shared_ptr<SymbolicExpression::Node>& node) {
    if (!node) {
        return 0;
    }

    // Use structural key for memoization
    static thread_local std::unordered_map<std::string, std::size_t> cache;
    const std::string key = node_structural_key(node);
    auto it = cache.find(key);
    if (it != cache.end()) {
        return it->second;
    }

    std::size_t count = 1;  // This node
    count += count_nodes(node->left);
    count += count_nodes(node->right);

    cache[key] = count;
    return count;
}

bool combine_all_like_additive_terms(const SymbolicExpression& expression,
                                     SymbolicExpression* combined_expression) {
    std::vector<SymbolicExpression> terms;
    collect_additive_expressions(expression, &terms);
    if (terms.size() < 3) {
        return false;
    }

    std::vector<SymbolicExpression> combined_terms;
    bool changed = false;
    for (const SymbolicExpression& term : terms) {
        bool merged = false;
        for (SymbolicExpression& existing : combined_terms) {
            SymbolicExpression combined;
            if (try_combine_like_terms(existing, term, 1.0, &combined)) {
                existing = combined;
                changed = true;
                merged = true;
                break;
            }
        }
        if (!merged) {
            combined_terms.push_back(term);
        }
    }

    if (!changed) {
        return false;
    }

    std::vector<SymbolicExpression> nonzero_terms;
    for (const SymbolicExpression& term : combined_terms) {
        if (!expr_is_zero(term)) {
            nonzero_terms.push_back(term);
        }
    }
    *combined_expression = make_sorted_sum(nonzero_terms);
    return true;
}

bool try_extract_common_symbolic_factor(const SymbolicExpression& expression,
                                        SymbolicExpression* common_factor,
                                        SymbolicExpression* remaining) {
    if (expression.node_->type != NodeType::kAdd && expression.node_->type != NodeType::kSubtract) {
        return false;
    }
    
    std::vector<SymbolicExpression> terms;
    collect_additive_expressions(expression, &terms);
    if (terms.size() < 2) return false;
    
    double dummy_num = 1.0;
    std::vector<SymbolicExpression> common_factors;
    collect_multiplicative_terms(terms[0], &dummy_num, &common_factors);
    
    if (common_factors.empty()) return false;
    
    for (size_t i = 1; i < terms.size(); ++i) {
        double t_num = 1.0;
        std::vector<SymbolicExpression> t_factors;
        collect_multiplicative_terms(terms[i], &t_num, &t_factors);
        
        std::vector<SymbolicExpression> new_common;
        std::vector<bool> used(t_factors.size(), false);
        for (const auto& cf : common_factors) {
            const std::string c_key = node_structural_key(cf.node_);
            for (size_t j = 0; j < t_factors.size(); ++j) {
                if (!used[j] && c_key == node_structural_key(t_factors[j].node_)) {
                    new_common.push_back(cf);
                    used[j] = true;
                    break;
                }
            }
        }
        common_factors = new_common;
        if (common_factors.empty()) break;
    }
    
    if (common_factors.empty()) return false;
    
    *common_factor = make_sorted_product(1.0, common_factors).simplify();
    
    std::vector<SymbolicExpression> new_terms;
    for (const auto& term : terms) {
        new_terms.push_back(make_divide(term, *common_factor).simplify());
    }
    *remaining = make_sorted_sum(new_terms).simplify();
    return true;
}

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
    static thread_local int depth = 0;
    static constexpr int kMaxDepth = 512;
    if (++depth > kMaxDepth) {
        --depth;
        return expression; // 达到最大深度，不再继续简化
    }

    struct DepthGuard {
        int* d;
        ~DepthGuard() { if (d) (*d)--; }
    } guard{&depth};

    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
            return expression;
        
        case NodeType::kFunction: {
            const SymbolicExpression argument = SymbolicExpression(node->left).simplify();
            double numeric = 0.0;
            if (argument.is_number(&numeric)) {
                if (node->text == "asin") return SymbolicExpression::number(mymath::asin(numeric));
                if (node->text == "acos") return SymbolicExpression::number(mymath::acos(numeric));
                if (node->text == "atan") return SymbolicExpression::number(mymath::atan(numeric));
                if (node->text == "sin") return SymbolicExpression::number(mymath::sin(numeric));
                if (node->text == "cos") return SymbolicExpression::number(mymath::cos(numeric));
                if (node->text == "tan") return SymbolicExpression::number(mymath::tan(numeric));
                if (node->text == "exp") return SymbolicExpression::number(mymath::exp(numeric));
                if (node->text == "sinh") return SymbolicExpression::number(mymath::sinh(numeric));
                if (node->text == "cosh") return SymbolicExpression::number(mymath::cosh(numeric));
                if (node->text == "tanh") return SymbolicExpression::number(mymath::tanh(numeric));
                if (node->text == "ln") {
                    if (numeric <= 0.0) return make_function(node->text, argument);
                    return SymbolicExpression::number(mymath::ln(numeric));
                }
                if (node->text == "sqrt") {
                    if (numeric < 0.0) return make_function(node->text, argument);
                    double root = mymath::sqrt(numeric);
                    if (mymath::is_near_zero(root * root - numeric, 1e-12) && mymath::is_integer(root, 1e-10)) {
                        return SymbolicExpression::number(root);
                    }
                    if (mymath::is_integer(numeric, 1e-10) && numeric > 0.0) {
                        long long val = static_cast<long long>(numeric + 0.5);
                        long long extracted = 1;
                        for (long long i = 2; i * i <= val; ++i) {
                            while (val % (i * i) == 0) {
                                extracted *= i;
                                val /= (i * i);
                            }
                        }
                        if (extracted > 1) {
                            return make_multiply(SymbolicExpression::number(static_cast<double>(extracted)),
                                                 make_function("sqrt", SymbolicExpression::number(static_cast<double>(val)))).simplify();
                        }
                    }
                    return make_function("sqrt", SymbolicExpression::number(numeric));
                }
                if (node->text == "erf") return SymbolicExpression::number(mymath::erf(numeric));
                if (node->text == "erfc") return SymbolicExpression::number(mymath::erfc(numeric));
                if (node->text == "gamma") {
                    if (mymath::is_integer(numeric, 1e-10) && numeric <= 0) return make_function(node->text, argument);
                    return SymbolicExpression::number(mymath::gamma(numeric));
                }
                if (node->text == "abs") return SymbolicExpression::number(mymath::abs(numeric));
                if (node->text == "floor") return SymbolicExpression::number(mymath::floor(numeric));
                if (node->text == "ceil") return SymbolicExpression::number(mymath::ceil(numeric));
                if (node->text == "cbrt") return SymbolicExpression::number(mymath::cbrt(numeric));
                if (node->text == "sign") {
                    if (mymath::is_near_zero(numeric, kFormatEps)) return SymbolicExpression::number(0.0);
                    return SymbolicExpression::number(numeric > 0.0 ? 1.0 : -1.0);
                }
            }

            if (node->text == "ln" && expr_is_variable(argument, "e")) return SymbolicExpression::number(1.0);
            if (node->text == "exp" && argument.node_->type == NodeType::kFunction && argument.node_->text == "ln") {
                const SymbolicExpression inner(argument.node_->left);
                if (is_known_positive_expression(inner.simplify())) return inner.simplify();
            }
            if (node->text == "exp" && argument.node_->type == NodeType::kAdd) {
                SymbolicExpression left(argument.node_->left);
                SymbolicExpression right(argument.node_->right);
                if (left.node_->type == NodeType::kFunction && left.node_->text == "ln") {
                    return (SymbolicExpression(left.node_->left) * make_function("exp", right)).simplify();
                }
                if (right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
                    return (SymbolicExpression(right.node_->left) * make_function("exp", left)).simplify();
                }
            }
            if (node->text == "exp" && argument.node_->type == NodeType::kSubtract) {
                SymbolicExpression left(argument.node_->left);
                SymbolicExpression right(argument.node_->right);
                if (left.node_->type == NodeType::kFunction && left.node_->text == "ln") {
                    return (SymbolicExpression(left.node_->left) / make_function("exp", right)).simplify();
                }
                if (right.node_->type == NodeType::kFunction && right.node_->text == "ln") {
                    return (make_function("exp", left) / SymbolicExpression(right.node_->left)).simplify();
                }
            }
            if (node->text == "ln" && argument.node_->type == NodeType::kFunction && argument.node_->text == "exp") {
                return SymbolicExpression(argument.node_->left).simplify();
            }
            if (node->text == "ln" && argument.node_->type == NodeType::kPower) {
                double exponent = 0.0;
                if (SymbolicExpression(argument.node_->right).is_number(&exponent) &&
                    is_known_positive_expression(SymbolicExpression(argument.node_->left))) {
                    return make_multiply(SymbolicExpression::number(exponent),
                                         make_function("ln", SymbolicExpression(argument.node_->left))).simplify();
                }
            }
            if (node->text == "ln" && argument.node_->type == NodeType::kMultiply) {
                SymbolicExpression mult_left(argument.node_->left);
                SymbolicExpression mult_right(argument.node_->right);
                if (is_known_positive_expression(mult_left) && is_known_positive_expression(mult_right)) {
                    return (make_function("ln", mult_left) + make_function("ln", mult_right)).simplify();
                }
            }
            if (node->text == "sqrt" && argument.node_->type == NodeType::kPower) {
                double exponent = 0.0;
                if (SymbolicExpression(argument.node_->right).is_number(&exponent) && mymath::is_near_zero(exponent - 2.0, kFormatEps)) {
                    return make_function("abs", SymbolicExpression(argument.node_->left)).simplify();
                }
            }
            if (node->text == "abs" && argument.node_->type == NodeType::kFunction && (argument.node_->text == "abs" || argument.node_->text == "sqrt")) {
                return argument;
            }
            if ((node->text == "sin" || node->text == "tan" || node->text == "sinh" || node->text == "tanh") && argument.node_->type == NodeType::kNegate) {
                return make_negate(make_function(node->text, SymbolicExpression(argument.node_->left))).simplify();
            }
            if ((node->text == "cos" || node->text == "cosh" || node->text == "abs") && argument.node_->type == NodeType::kNegate) {
                return make_function(node->text, SymbolicExpression(argument.node_->left)).simplify();
            }

            double pi_multiple = 0.0;
            if (decompose_numeric_multiple_of_symbol(argument, "pi", &pi_multiple)) {
                if (node->text == "sin") {
                    if (numeric_matches_any(pi_multiple, {0.0, 1.0, -1.0})) return SymbolicExpression::number(0.0);
                    if (numeric_matches_any(pi_multiple, {0.5})) return SymbolicExpression::number(1.0);
                    if (numeric_matches_any(pi_multiple, {-0.5})) return SymbolicExpression::number(-1.0);
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, 5.0 / 6.0})) return half_symbol();
                    if (numeric_matches_any(pi_multiple, {-1.0 / 6.0, -5.0 / 6.0})) return make_negate(half_symbol()).simplify();
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, 2.0 / 3.0})) return make_divide(sqrt3_symbol(), SymbolicExpression::number(2.0)).simplify();
                    if (numeric_matches_any(pi_multiple, {-1.0 / 3.0, -2.0 / 3.0})) return make_negate(make_divide(sqrt3_symbol(), SymbolicExpression::number(2.0))).simplify();
                }
                if (node->text == "cos") {
                    if (numeric_matches_any(pi_multiple, {0.0})) return SymbolicExpression::number(1.0);
                    if (numeric_matches_any(pi_multiple, {1.0, -1.0})) return SymbolicExpression::number(-1.0);
                    if (numeric_matches_any(pi_multiple, {0.5, -0.5})) return SymbolicExpression::number(0.0);
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, -1.0 / 3.0})) return half_symbol();
                    if (numeric_matches_any(pi_multiple, {2.0 / 3.0, -2.0 / 3.0})) return make_negate(half_symbol()).simplify();
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, -1.0 / 6.0})) return make_divide(sqrt3_symbol(), SymbolicExpression::number(2.0)).simplify();
                    if (numeric_matches_any(pi_multiple, {5.0 / 6.0, -5.0 / 6.0})) return make_negate(make_divide(sqrt3_symbol(), SymbolicExpression::number(2.0))).simplify();
                }
                if (node->text == "tan") {
                    if (numeric_matches_any(pi_multiple, {0.0, 1.0, -1.0})) return SymbolicExpression::number(0.0);
                    if (numeric_matches_any(pi_multiple, {0.25})) return SymbolicExpression::number(1.0);
                    if (numeric_matches_any(pi_multiple, {-0.25})) return SymbolicExpression::number(-1.0);
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, -5.0 / 6.0})) return make_divide(SymbolicExpression::number(1.0), sqrt3_symbol()).simplify();
                    if (numeric_matches_any(pi_multiple, {-1.0 / 6.0, 5.0 / 6.0})) return make_negate(make_divide(SymbolicExpression::number(1.0), sqrt3_symbol())).simplify();
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, -2.0 / 3.0})) return sqrt3_symbol();
                    if (numeric_matches_any(pi_multiple, {-1.0 / 3.0, 2.0 / 3.0})) return make_negate(sqrt3_symbol()).simplify();
                }
            }
            return make_function(node->text, argument);
        }

        case NodeType::kNegate: {
            const SymbolicExpression operand = SymbolicExpression(node->left).simplify();
            double value = 0.0;
            if (operand.is_number(&value)) return SymbolicExpression::number(-value);
            if (operand.node_->type == NodeType::kNegate) return SymbolicExpression(operand.node_->left).simplify();
            return make_negate(operand);
        }

        default: break;
    }

    const SymbolicExpression left = SymbolicExpression(node->left).simplify();
    const SymbolicExpression right = SymbolicExpression(node->right).simplify();
    double left_value = 0.0;
    double right_value = 0.0;

    switch (node->type) {
        case NodeType::kAdd:
            if (left.is_number(&left_value) && right.is_number(&right_value)) return SymbolicExpression::number(left_value + right_value);
            if (expr_is_zero(left)) return right;
            if (expr_is_zero(right)) return left;
            {
                std::string left_arg, right_arg;
                if ((is_squared_function(left, "sin", &left_arg) && is_squared_function(right, "cos", &right_arg) && left_arg == right_arg) ||
                    (is_squared_function(left, "cos", &left_arg) && is_squared_function(right, "sin", &right_arg) && left_arg == right_arg)) {
                    return SymbolicExpression::number(1.0);
                }
            }
            {
                SymbolicExpression combined;
                if (try_combine_like_terms(left, right, 1.0, &combined)) return combined;
            }
            {
                const SymbolicExpression sum = make_add(left, right);
                if (is_single_variable_polynomial(sum)) return maybe_canonicalize_polynomial(sum);
            }
            {
                SymbolicExpression factored;
                if (try_factor_common_terms(left, right, 1.0, &factored)) return factored;
            }
            {
                std::vector<SymbolicExpression> terms;
                collect_additive_expressions(make_add(left, right), &terms);
                SymbolicExpression combined_terms;
                if (combine_all_like_additive_terms(make_add(left, right), &combined_terms)) return combined_terms;
                return make_sorted_sum(terms);
            }

        case NodeType::kSubtract:
            if (left.is_number(&left_value) && right.is_number(&right_value)) return SymbolicExpression::number(left_value - right_value);
            if (expr_is_zero(right)) return left;
            if (expr_is_zero(left)) return make_negate(right).simplify();
            if (right.node_->type == NodeType::kNegate) return make_add(left, SymbolicExpression(right.node_->left)).simplify();
            {
                std::string left_arg, right_arg;
                if (is_squared_function(left, "sec", &left_arg) && is_squared_function(right, "tan", &right_arg) && left_arg == right_arg) return SymbolicExpression::number(1.0);
                if (is_squared_function(left, "csc", &left_arg) && is_squared_function(right, "cot", &right_arg) && left_arg == right_arg) return SymbolicExpression::number(1.0);
            }
            {
                SymbolicExpression combined;
                if (try_combine_like_terms(left, right, -1.0, &combined)) return combined;
            }
            {
                SymbolicExpression factored;
                if (try_factor_common_terms(left, right, -1.0, &factored)) return factored;
            }
            return make_subtract(left, right);

        case NodeType::kMultiply:
            if (left.is_number(&left_value) && right.is_number(&right_value)) return SymbolicExpression::number(left_value * right_value);
            if (expr_is_zero(left) || expr_is_zero(right)) return SymbolicExpression::number(0.0);
            if (expr_is_one(left)) return right;
            if (expr_is_one(right)) return left;
            if (left.is_variable_named("i") && right.is_variable_named("i")) return SymbolicExpression::number(-1.0);
            if (expr_is_minus_one(left)) return make_negate(right).simplify();
            if (expr_is_minus_one(right)) return make_negate(left).simplify();
            if (left.node_->type == NodeType::kFunction && right.node_->type == NodeType::kFunction && left.node_->text == "exp" && right.node_->text == "exp") {
                return make_function("exp", (SymbolicExpression(left.node_->left) + SymbolicExpression(right.node_->left))).simplify();
            }
            {
                SymbolicExpression left_base, right_base;
                double left_exp, right_exp;
                decompose_power_factor(left, &left_base, &left_exp);
                decompose_power_factor(right, &right_base, &right_exp);
                if (expressions_match(left_base, right_base)) return rebuild_power_difference(left_base, left_exp + right_exp);
            }
            {
                double numeric_factor = 1.0;
                std::vector<SymbolicExpression> symbolic_factors;
                collect_multiplicative_terms(left, &numeric_factor, &symbolic_factors);
                collect_multiplicative_terms(right, &numeric_factor, &symbolic_factors);
                return make_sorted_product(numeric_factor, symbolic_factors);
            }

        case NodeType::kDivide:
            if (left.is_number(&left_value) && right.is_number(&right_value)) return SymbolicExpression::number(left_value / right_value);
            if (expr_is_zero(left)) return SymbolicExpression::number(0.0);
            if (expr_is_one(right)) return left;
            {
                SymbolicExpression reduced;
                if (try_reduce_polynomial_quotient(left, right, &reduced)) return reduced;
                if (try_reduce_polynomial_gcd_quotient(left, right, &reduced)) return reduced;
            }
            {
                SymbolicExpression left_base, right_base;
                double left_exp, right_exp;
                decompose_power_factor(left, &left_base, &left_exp);
                decompose_power_factor(right, &right_base, &right_exp);
                if (expressions_match(left_base, right_base)) return rebuild_power_difference(left_base, left_exp - right_exp);
            }
            {
                SymbolicExpression common, rest;
                if (try_extract_common_symbolic_factor(left, &common, &rest)) return make_multiply(common, make_divide(rest, right)).simplify();
                if (try_extract_common_symbolic_factor(right, &common, &rest)) return make_divide(left, make_multiply(common, rest)).simplify();
            }
            {
                double num_c = 1.0, den_c = 1.0;
                std::vector<SymbolicExpression> num_f, den_f;
                collect_division_factors(left, &num_c, &num_f);
                collect_division_factors(right, &den_c, &den_f);
                if (!mymath::is_near_zero(den_c, kFormatEps)) {
                    std::vector<bool> den_used(den_f.size(), false);
                    std::vector<SymbolicExpression> red_num;
                    for (const auto& nf : num_f) {
                        const std::string nk = node_structural_key(nf.node_);
                        bool canceled = false;
                        for (size_t i = 0; i < den_f.size(); ++i) {
                            if (!den_used[i] && nk == node_structural_key(den_f[i].node_)) {
                                den_used[i] = true;
                                canceled = true;
                                break;
                            }
                        }
                        if (!canceled) red_num.push_back(nf);
                    }
                    std::vector<SymbolicExpression> red_den;
                    for (size_t i = 0; i < den_f.size(); ++i) if (!den_used[i]) red_den.push_back(den_f[i]);
                    SymbolicExpression num_expr = rebuild_product_expression(num_c / den_c, red_num);
                    SymbolicExpression den_expr = rebuild_product_expression(1.0, red_den);
                    if (expr_is_one(den_expr)) return num_expr;
                    return make_divide(num_expr, den_expr);
                }
            }
            return make_divide(left, right);

        case NodeType::kPower:
            if (left.node_->type == NodeType::kE) return make_function("exp", right).simplify();
            if (left.is_variable_named("i")) {
                double exp_v;
                if (right.is_number(&exp_v)) {
                    int exp_i = static_cast<int>(mymath::round(exp_v));
                    if (mymath::abs(exp_v - exp_i) < 1e-10) {
                        int m4 = ((exp_i % 4) + 4) % 4;
                        if (m4 == 0) return SymbolicExpression::number(1.0);
                        if (m4 == 1) return left;
                        if (m4 == 2) return SymbolicExpression::number(-1.0);
                        if (m4 == 3) return make_negate(left).simplify();
                    }
                }
            }
            if (left.is_number(&left_value)) {
                if (mymath::is_near_zero(left_value, kFormatEps)) {
                    if (right.is_number(&right_value)) return mymath::is_near_zero(right_value, kFormatEps) ? SymbolicExpression::number(1.0) : SymbolicExpression::number(0.0);
                }
                if (mymath::is_near_zero(left_value - 1.0, kFormatEps)) return SymbolicExpression::number(1.0);
            }
            if (right.is_number(&right_value)) {
                if (mymath::is_near_zero(right_value, kFormatEps)) return SymbolicExpression::number(1.0);
                if (mymath::is_near_zero(right_value - 1.0, kFormatEps)) return left;
            }
            if (left.node_->type == NodeType::kPower) {
                double inner_exp;
                if (SymbolicExpression(left.node_->right).is_number(&inner_exp) && right.is_number(&right_value)) {
                    return make_power(SymbolicExpression(left.node_->left).simplify(), SymbolicExpression::number(inner_exp * right_value)).simplify();
                }
            }
            if (left.is_number(&left_value) && right.is_number(&right_value)) return SymbolicExpression::number(mymath::pow(left_value, right_value));
            return make_power(left, right);

        case NodeType::kVector: {
            std::vector<SymbolicExpression> components;
            for (const auto& child : node->children) components.push_back(simplify_once(SymbolicExpression(child)).simplify());
            return SymbolicExpression::vector(components);
        }

        case NodeType::kTensor: {
            std::vector<std::vector<SymbolicExpression>> rows;
            for (const auto& row_node : node->children) {
                std::vector<SymbolicExpression> row;
                if (row_node->type == NodeType::kVector) {
                    for (const auto& comp : row_node->children) row.push_back(simplify_once(SymbolicExpression(comp)).simplify());
                }
                rows.push_back(row);
            }
            return SymbolicExpression::tensor(rows);
        }

        case NodeType::kDifferentialOp: {
            SymbolicExpression op = simplify_once(SymbolicExpression(node->left));
            if (node->text == "div" && op.node_->type == NodeType::kDifferentialOp && op.node_->text == "grad") return make_function("laplacian", simplify_once(SymbolicExpression(op.node_->left)));
            if (node->text == "div" && op.node_->type == NodeType::kDifferentialOp && op.node_->text == "curl") return SymbolicExpression::number(0.0);
            if (node->text == "curl" && op.node_->type == NodeType::kDifferentialOp && op.node_->text == "grad") return SymbolicExpression::number(0.0);
            return SymbolicExpression(make_function(node->text, op).node_);
        }

        default: break;
    }
    return expression;
}

SymbolicExpression simplify_impl(const SymbolicExpression& expression) {
    static thread_local int simplify_depth = 0;
    struct Guard { int* d; Guard(int* v) : d(v) { (*d)++; } ~Guard() { (*d)--; } };
    if (simplify_depth > 0) return simplify_once(expression);
    Guard guard(&simplify_depth);
    SymbolicExpression current = expression;
    constexpr int kMaxPasses = 24;
    for (int p = 0; p < kMaxPasses; ++p) {
        const std::string ck = node_structural_key(current.node_);
        SymbolicExpression next = simplify_once(current);
        if (node_structural_key(next.node_) == ck) return next;
        if (count_nodes(next.node_) > 10000) return current;
        current = next;
    }
    return current;
}

SymbolicExpression simplify_with_budget_impl(const SymbolicExpression& expression, std::size_t max_nodes) {
    static thread_local int simplify_depth = 0;
    struct Guard { int* d; Guard(int* v) : d(v) { (*d)++; } ~Guard() { (*d)--; } };
    if (simplify_depth > 0) return simplify_once(expression);
    Guard guard(&simplify_depth);
    SymbolicExpression current = expression;
    constexpr int kMaxPasses = 24;
    for (int p = 0; p < kMaxPasses; ++p) {
        const std::string ck = node_structural_key(current.node_);
        SymbolicExpression next = simplify_once(current);
        if (node_structural_key(next.node_) == ck) return next;
        if (count_nodes(next.node_) > max_nodes) return current;
        current = next;
    }
    return current;
}

SymbolicExpression expand_impl(const SymbolicExpression& expression) {
    const std::shared_ptr<SymbolicExpression::Node>& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
            return expression;
        case NodeType::kFunction:
            return make_function(node->text, expand_impl(SymbolicExpression(node->left))).simplify();
        case NodeType::kNegate:
            return make_negate(expand_impl(SymbolicExpression(node->left))).simplify();
        case NodeType::kAdd:
            return make_add(expand_impl(SymbolicExpression(node->left)),
                            expand_impl(SymbolicExpression(node->right))).simplify();
        case NodeType::kSubtract:
            return make_subtract(expand_impl(SymbolicExpression(node->left)),
                                 expand_impl(SymbolicExpression(node->right))).simplify();
        case NodeType::kMultiply: {
            SymbolicExpression left = expand_impl(SymbolicExpression(node->left));
            SymbolicExpression right = expand_impl(SymbolicExpression(node->right));
            if (left.node_->type == NodeType::kAdd) return make_add(expand_impl(make_multiply(SymbolicExpression(left.node_->left), right)), expand_impl(make_multiply(SymbolicExpression(left.node_->right), right))).simplify();
            if (left.node_->type == NodeType::kSubtract) return make_subtract(expand_impl(make_multiply(SymbolicExpression(left.node_->left), right)), expand_impl(make_multiply(SymbolicExpression(left.node_->right), right))).simplify();
            if (right.node_->type == NodeType::kAdd) return make_add(expand_impl(make_multiply(left, SymbolicExpression(right.node_->left))), expand_impl(make_multiply(left, SymbolicExpression(right.node_->right)))).simplify();
            if (right.node_->type == NodeType::kSubtract) return make_subtract(expand_impl(make_multiply(left, SymbolicExpression(right.node_->left))), expand_impl(make_multiply(left, SymbolicExpression(right.node_->right)))).simplify();
            return make_multiply(left, right).simplify();
        }
        case NodeType::kDivide: {
            SymbolicExpression left = expand_impl(SymbolicExpression(node->left));
            SymbolicExpression right = expand_impl(SymbolicExpression(node->right));
            if (left.node_->type == NodeType::kAdd) return make_add(expand_impl(make_divide(SymbolicExpression(left.node_->left), right)), expand_impl(make_divide(SymbolicExpression(left.node_->right), right))).simplify();
            if (left.node_->type == NodeType::kSubtract) return make_subtract(expand_impl(make_divide(SymbolicExpression(left.node_->left), right)), expand_impl(make_divide(SymbolicExpression(left.node_->right), right))).simplify();
            return make_divide(left, right).simplify();
        }
        case NodeType::kPower: {
            SymbolicExpression left = expand_impl(SymbolicExpression(node->left));
            SymbolicExpression right = expand_impl(SymbolicExpression(node->right));
            double n = 0.0;
            if (right.is_number(&n) && n > 1.0 && mymath::is_integer(n, 1e-10)) {
                if (left.node_->type == NodeType::kAdd || left.node_->type == NodeType::kSubtract) {
                    SymbolicExpression rest = make_power(left, SymbolicExpression::number(n - 1.0));
                    return expand_impl(make_multiply(left, rest));
                }
            }
            return make_power(left, right).simplify();
        }
        case NodeType::kVector: {
            std::vector<SymbolicExpression> components;
            for (const auto& child : node->children) components.push_back(expand_impl(SymbolicExpression(child)));
            return SymbolicExpression::vector(components).simplify();
        }
        case NodeType::kTensor: {
            std::vector<std::vector<SymbolicExpression>> rows;
            for (const auto& row_node : node->children) {
                std::vector<SymbolicExpression> row;
                for (const auto& child : row_node->children) row.push_back(expand_impl(SymbolicExpression(child)));
                rows.push_back(std::move(row));
            }
            return SymbolicExpression::tensor(rows).simplify();
        }
        case NodeType::kDifferentialOp:
            return expression;
        case NodeType::kRootOf:
            // RootOf 节点不需要展开，保持原样
            return expression;
    }
    return expression;
}


}  // namespace symbolic_expression_internal

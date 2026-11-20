// ============================================================================
// 符号微分实现
// ============================================================================

#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/algebra/polynomial/symbolic_polynomial.h"

#include "math/mymath.h"
#include "polynomial/polynomial.h"

#include <algorithm>
#include <vector>

namespace symbolic_expression_internal {

// ============================================================================
// ============================================================================

/**
 * @brief 符号微分的内部实现（无缓存）
 *
 * 递归应用微分规则：
 * - 常数：导数为 0
 * - 变量：导数为 1（对目标变量）或 0
 * - 加减法：线性性
 * - 乘法：乘积法则
 * - 除法：商法则
 * - 幂函数：幂法则 + 对数微分
 * - 三角/反三角函数：标准公式
 * - 指数/对数函数：标准公式
 * - 特殊函数（erf, sign, step, delta）
 */
SymbolicExpression derivative_uncached(const SymbolicExpression& expression,
                                       const std::string& variable_name) {
    const std::shared_ptr<SymbolicExpression::Node>& node_ = expression.node_;
    switch (node_->type) {
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
            return SymbolicExpression::number(0.0L);
        case NodeType::kVariable:
            return SymbolicExpression::number(node_->text == variable_name ? 1.0L : 0.0L);
        case NodeType::kNegate:
            return make_negate(SymbolicExpression(node_->left).derivative(variable_name)).simplify();
        case NodeType::kAdd:
            return make_add(SymbolicExpression(node_->left).derivative(variable_name),
                            SymbolicExpression(node_->right).derivative(variable_name)).simplify();
        case NodeType::kSubtract:
            return make_subtract(SymbolicExpression(node_->left).derivative(variable_name),
                                 SymbolicExpression(node_->right).derivative(variable_name)).simplify();
        case NodeType::kMultiply: {
            const SymbolicExpression left(node_->left);
            const SymbolicExpression right(node_->right);
            return make_add(make_multiply(left.derivative(variable_name), right),
                            make_multiply(left, right.derivative(variable_name)))
                .simplify();
        }
        case NodeType::kDivide: {
            const SymbolicExpression left(node_->left);
            const SymbolicExpression right(node_->right);
            return make_divide(
                       make_subtract(make_multiply(left.derivative(variable_name), right),
                                     make_multiply(left, right.derivative(variable_name))),
                       make_power(right, SymbolicExpression::number(2.0)))
                .simplify();
        }
        case NodeType::kPower: {
            const SymbolicExpression base(node_->left);
            const SymbolicExpression exponent(node_->right);
            if (exponent.is_constant(variable_name)) {
                return make_multiply(
                           make_multiply(exponent,
                                         make_power(base, make_subtract(exponent, SymbolicExpression::number(1.0L)))),
                           base.derivative(variable_name))
                    .simplify();
            }
            if (base.is_constant(variable_name)) {
                return make_multiply(
                           make_multiply(expression, make_function("ln", base)),
                           exponent.derivative(variable_name))
                    .simplify();
            }
            return make_multiply(
                       expression,
                       make_add(
                           make_multiply(exponent.derivative(variable_name), make_function("ln", base)),
                           make_multiply(exponent,
                                         make_divide(base.derivative(variable_name), base))))
                .simplify();
        }
        case NodeType::kFunction: {
            if ((node_->text == "besselj" || node_->text == "bessely" ||
                 node_->text == "besseli" || node_->text == "besselk") &&
                node_->children.size() >= 2) {
                const SymbolicExpression order(node_->children[0]);
                const SymbolicExpression argument(node_->children[1]);
                if (!order.is_constant(variable_name)) {
                    throw std::runtime_error(
                        "symbolic derivative does not support variable-order " + node_->text);
                }
                const SymbolicExpression inner = argument.derivative(variable_name);
                Scalar order_val = Scalar(0.0L);
                if (order.is_number(&order_val) && mymath::is_near_zero(order_val)) {
                    const SymbolicExpression first = SymbolicExpression::function(
                        node_->text, std::vector<SymbolicExpression>{SymbolicExpression::number(1.0L), argument});
                    if (node_->text == "besseli") {
                        return make_multiply(first, inner).simplify();
                    } else {
                        return make_multiply(make_negate(first), inner).simplify();
                    }
                }
                const SymbolicExpression lower = SymbolicExpression::function(
                    node_->text, std::vector<SymbolicExpression>{
                        (order - SymbolicExpression::number(1.0L)).simplify(), argument});
                const SymbolicExpression upper = SymbolicExpression::function(
                    node_->text, std::vector<SymbolicExpression>{
                        (order + SymbolicExpression::number(1.0L)).simplify(), argument});
                if (node_->text == "besselj" || node_->text == "bessely") {
                    return make_multiply(
                        make_multiply(SymbolicExpression::number(0.5L),
                                      make_subtract(lower, upper)),
                        inner).simplify();
                } else if (node_->text == "besseli") {
                    return make_multiply(
                        make_multiply(SymbolicExpression::number(0.5L),
                                      make_add(lower, upper)),
                        inner).simplify();
                } else { // besselk
                    return make_multiply(
                        make_multiply(SymbolicExpression::number(-0.5L),
                                      make_add(lower, upper)),
                        inner).simplify();
                }
            }
            if (node_->text == "polygamma" && node_->children.size() >= 2) {
                const SymbolicExpression order(node_->children[0]);
                const SymbolicExpression argument(node_->children[1]);
                if (!order.is_constant(variable_name)) {
                    throw std::runtime_error("symbolic derivative does not support variable-order polygamma");
                }
                const SymbolicExpression inner = argument.derivative(variable_name);
                const SymbolicExpression next_order = (order + SymbolicExpression::number(1.0L)).simplify();
                const SymbolicExpression next_poly = SymbolicExpression::function(
                    "polygamma", std::vector<SymbolicExpression>{next_order, argument});
                return make_multiply(next_poly, inner).simplify();
            }
            if (node_->text == "atan2" && node_->children.size() >= 2) {
                const SymbolicExpression y(node_->children[0]);
                const SymbolicExpression x(node_->children[1]);
                const SymbolicExpression dy = y.derivative(variable_name);
                const SymbolicExpression dx = x.derivative(variable_name);
                const SymbolicExpression denom = make_add(make_power(x, SymbolicExpression::number(2.0L)),
                                                          make_power(y, SymbolicExpression::number(2.0L)));
                const SymbolicExpression numer = make_subtract(make_multiply(x, dy), make_multiply(y, dx));
                return make_divide(numer, denom).simplify();
            }
            if ((node_->text == "min" || node_->text == "max") && node_->children.size() >= 2) {
                const SymbolicExpression u(node_->children[0]);
                const SymbolicExpression v(node_->children[1]);
                const SymbolicExpression du = u.derivative(variable_name);
                const SymbolicExpression dv = v.derivative(variable_name);
                if (node_->text == "min") {
                    return make_add(make_multiply(du, make_function("step", make_subtract(v, u))),
                                    make_multiply(dv, make_function("step", make_subtract(u, v)))).simplify();
                } else {
                    return make_add(make_multiply(du, make_function("step", make_subtract(u, v))),
                                    make_multiply(dv, make_function("step", make_subtract(v, u)))).simplify();
                }
            }
            if (node_->text == "pow" && node_->children.size() >= 2) {
                const SymbolicExpression base(node_->children[0]);
                const SymbolicExpression exponent(node_->children[1]);
                return make_power(base, exponent).derivative(variable_name).simplify();
            }
            if (node_->text == "log" && node_->children.size() >= 2) {
                // log(b, x) = ln(x) / ln(b)
                const SymbolicExpression base(node_->children[0]);
                const SymbolicExpression arg(node_->children[1]);
                return make_divide(make_function("ln", arg), make_function("ln", base)).derivative(variable_name).simplify();
            }
            if (node_->children.size() > 1) {
                if (node_->text == "Integral" && node_->children.size() >= 2) {
                    SymbolicExpression integrand(node_->children[0]);
                    SymbolicExpression int_var(node_->children[1]);
                    if (int_var.is_variable_named(variable_name)) {
                        return integrand.simplify();
                    }
                    throw std::runtime_error("symbolic derivative: differentiation variable does not match Integral variable");
                }
                throw std::runtime_error(
                    "symbolic derivative does not support multi-argument function: " + node_->text);
            }
            if (node_->text == "Integral") {
                throw std::runtime_error("symbolic derivative does not support unparameterized Integral placeholder");
            }
            const SymbolicExpression argument(node_->left ? node_->left : node_->children[0]);
            const SymbolicExpression inner = argument.derivative(variable_name);
            if (node_->text == "asin") {
                return make_divide(
                           inner,
                           make_function("sqrt",
                                         make_subtract(SymbolicExpression::number(Scalar(1.0L)),
                                                       make_power(argument, SymbolicExpression::number(2.0)))))
                    .simplify();
            }
            if (node_->text == "acos") {
                return make_negate(
                           make_divide(inner,
                                       make_function("sqrt",
                                                     make_subtract(SymbolicExpression::number(Scalar(1.0L)),
                                                                   make_power(argument, SymbolicExpression::number(2.0))))))
                    .simplify();
            }
            if (node_->text == "atan") {
                return make_divide(inner,
                                   make_add(SymbolicExpression::number(Scalar(1.0L)), make_power(argument, SymbolicExpression::number(2.0))))
                    .simplify();
            }
            if (node_->text == "sin") {
                return make_multiply(make_function("cos", argument), inner).simplify();
            }
            if (node_->text == "cos") {
                return make_multiply(make_negate(make_function("sin", argument)), inner).simplify();
            }
            if (node_->text == "tan") {
                return make_multiply(make_divide(SymbolicExpression::number(Scalar(1.0L)),
                                                 make_power(make_function("cos", argument),
                                                            SymbolicExpression::number(2.0))),
                                     inner)
                    .simplify();
            }
            if (node_->text == "sec") {
                return make_multiply(
                           make_multiply(make_function("sec", argument),
                                         make_function("tan", argument)),
                           inner)
                    .simplify();
            }
            if (node_->text == "csc") {
                return make_multiply(
                           make_negate(make_multiply(make_function("csc", argument),
                                                     make_function("cot", argument))),
                           inner)
                    .simplify();
            }
            if (node_->text == "cot") {
                return make_multiply(
                           make_negate(make_power(make_function("csc", argument),
                                                  SymbolicExpression::number(2.0))),
                           inner)
                    .simplify();
            }
            if (node_->text == "exp") {
                return make_multiply(make_function("exp", argument), inner).simplify();
            }
            if (node_->text == "sinh") {
                return make_multiply(make_function("cosh", argument), inner).simplify();
            }
            if (node_->text == "cosh") {
                return make_multiply(make_function("sinh", argument), inner).simplify();
            }
            if (node_->text == "tanh") {
                return make_divide(inner,
                                   make_power(make_function("cosh", argument),
                                              SymbolicExpression::number(2.0)))
                    .simplify();
            }
            if (node_->text == "asinh") {
                return make_divide(
                           inner,
                           make_function("sqrt",
                                         make_add(make_power(argument, SymbolicExpression::number(2.0)),
                                                  SymbolicExpression::number(Scalar(1.0L)))))
                    .simplify();
            }
            if (node_->text == "acosh") {
                return make_divide(
                           inner,
                           make_function("sqrt",
                                         make_subtract(make_power(argument, SymbolicExpression::number(2.0)),
                                                       SymbolicExpression::number(Scalar(1.0L)))))
                    .simplify();
            }
            if (node_->text == "atanh") {
                return make_divide(
                           inner,
                           make_subtract(SymbolicExpression::number(Scalar(1.0L)),
                                         make_power(argument, SymbolicExpression::number(2.0))))
                    .simplify();
            }
            if (node_->text == "coth") {
                return make_divide(
                           make_negate(inner),
                           make_power(make_function("sinh", argument), SymbolicExpression::number(2.0)))
                    .simplify();
            }
            if (node_->text == "sech") {
                return make_multiply(
                           make_negate(make_multiply(make_function("sech", argument),
                                                     make_function("tanh", argument))),
                           inner)
                    .simplify();
            }
            if (node_->text == "csch") {
                return make_multiply(
                           make_negate(make_multiply(make_function("csch", argument),
                                                     make_function("coth", argument))),
                           inner)
                    .simplify();
            }
            if (node_->text == "ln" || node_->text == "log") {
                if (argument.node_->type == NodeType::kFunction && argument.node_->text == "abs") {
                    auto inner_node = argument.node_->left ? argument.node_->left : (!argument.node_->children.empty() ? argument.node_->children[0] : nullptr);
                    if (inner_node) {
                        SymbolicExpression inner_abs(inner_node);
                        return make_divide(inner_abs.derivative(variable_name), inner_abs).simplify();
                    }
                }
                return make_divide(inner, argument).simplify();
            }
            if (node_->text == "sqrt") {
                return make_divide(inner,
                                   make_multiply(SymbolicExpression::number(2.0), make_function("sqrt", argument)))
                    .simplify();
            }
            if (node_->text == "cbrt") {
                return make_divide(inner,
                                   make_multiply(SymbolicExpression::number(3.0),
                                                 make_power(make_function("cbrt", argument),
                                                            SymbolicExpression::number(2.0))))
                    .simplify();
            }
            if (node_->text == "erf") {
                const SymbolicExpression pi_val(make_unary(NodeType::kPi, nullptr));
                const SymbolicExpression exp_part = make_function("exp", make_negate(make_power(argument, SymbolicExpression::number(2.0))));
                const SymbolicExpression factor = make_divide(SymbolicExpression::number(2.0), make_function("sqrt", pi_val));
                return make_multiply(make_multiply(factor, exp_part), inner).simplify();
            }
            if (node_->text == "erfc") {
                const SymbolicExpression pi_val(make_unary(NodeType::kPi, nullptr));
                const SymbolicExpression exp_part = make_function("exp", make_negate(make_power(argument, SymbolicExpression::number(2.0))));
                const SymbolicExpression factor = make_divide(SymbolicExpression::number(-2.0), make_function("sqrt", pi_val));
                return make_multiply(make_multiply(factor, exp_part), inner).simplify();
            }
            if (node_->text == "erfi") {
                const SymbolicExpression pi_val(make_unary(NodeType::kPi, nullptr));
                const SymbolicExpression exp_part = make_function("exp", make_power(argument, SymbolicExpression::number(2.0)));
                const SymbolicExpression factor = make_divide(SymbolicExpression::number(2.0), make_function("sqrt", pi_val));
                return make_multiply(make_multiply(factor, exp_part), inner).simplify();
            }
            if (node_->text == "Ei") {
                return make_multiply(make_divide(make_function("exp", argument), argument), inner).simplify();
            }
            if (node_->text == "Si") {
                return make_multiply(make_divide(make_function("sin", argument), argument), inner).simplify();
            }
            if (node_->text == "Ci") {
                return make_multiply(make_divide(make_function("cos", argument), argument), inner).simplify();
            }
            if (node_->text == "abs") {
                return make_multiply(make_function("sign", argument), inner).simplify();
            }
            if (node_->text == "sign") {
                return make_multiply(
                           make_multiply(SymbolicExpression::number(2.0),
                                         make_function("delta", argument)),
                           inner)
                    .simplify();
            }
            if (node_->text == "step" || node_->text == "heaviside" || node_->text == "Heaviside") {
                return make_multiply(make_function("delta", argument), inner).simplify();
            }
            if (node_->text == "gamma") {
                // d/dx gamma(x) = gamma(x) * digamma(x)
                return make_multiply(make_multiply(expression, make_function("digamma", argument)), inner).simplify();
            }
            if (node_->text == "digamma") {
                // d/dx digamma(x) = polygamma(1, x)
                const SymbolicExpression poly = SymbolicExpression::function(
                    "polygamma", std::vector<SymbolicExpression>{SymbolicExpression::number(1.0L), argument});
                return make_multiply(poly, inner).simplify();
            }
            if (node_->text == "lambertw" || node_->text == "LambertW") {
                // d/dx W(x) = W(x) / (x * (1 + W(x)))
                const SymbolicExpression denom = make_multiply(argument, make_add(SymbolicExpression::number(1.0L), expression));
                return make_multiply(make_divide(expression, denom), inner).simplify();
            }
            if (node_->text == "sinc") {
                // d/dx sinc(x) = (cos(x) - sinc(x)) / x
                const SymbolicExpression numer = make_subtract(make_function("cos", argument), expression);
                return make_multiply(make_divide(numer, argument), inner).simplify();
            }
            if (node_->text == "besselj") {
                // d/dx J_0(x) = -J_1(x)
                const SymbolicExpression j1 = SymbolicExpression::function(
                    "besselj", std::vector<SymbolicExpression>{SymbolicExpression::number(1.0L), argument});
                return make_multiply(make_negate(j1), inner).simplify();
            }
            if (node_->text == "bessely") {
                const SymbolicExpression y1 = SymbolicExpression::function(
                    "bessely", std::vector<SymbolicExpression>{SymbolicExpression::number(1.0L), argument});
                return make_multiply(make_negate(y1), inner).simplify();
            }
            if (node_->text == "besseli") {
                const SymbolicExpression i1 = SymbolicExpression::function(
                    "besseli", std::vector<SymbolicExpression>{SymbolicExpression::number(1.0L), argument});
                return make_multiply(i1, inner).simplify();
            }
            if (node_->text == "besselk") {
                const SymbolicExpression k1 = SymbolicExpression::function(
                    "besselk", std::vector<SymbolicExpression>{SymbolicExpression::number(1.0L), argument});
                return make_multiply(make_negate(k1), inner).simplify();
            }
            throw std::runtime_error("symbolic derivative does not support function: " + node_->text);
        }
        case NodeType::kVector: {
            std::vector<SymbolicExpression> components;
            for (const auto& child : node_->children) {
                components.push_back(SymbolicExpression(child).derivative(variable_name));
            }
            return SymbolicExpression::vector(components).simplify();
        }
        case NodeType::kTensor: {
            std::vector<std::vector<SymbolicExpression>> rows;
            for (const auto& row_node : node_->children) {
                std::vector<SymbolicExpression> row;
                for (const auto& child : row_node->children) {
                    row.push_back(SymbolicExpression(child).derivative(variable_name));
                }
                rows.push_back(std::move(row));
            }
            return SymbolicExpression::tensor(rows).simplify();
        }
        case NodeType::kDifferentialOp:
            throw std::runtime_error("symbolic derivative does not support differential operator: " + node_->text);
        case NodeType::kRootOf:
            // RootOf 表示代数数，是常数，导数为 0
            return SymbolicExpression::number(0.0L);
    }
    throw std::runtime_error("unsupported symbolic derivative");
}

}  // namespace symbolic_expression_internal

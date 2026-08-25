/**
 * @file symbolic_commands_misc.cpp
 * @brief 其他符号命令实现
 *
 * 本文件实现了其他符号计算命令：
 * - solve: 方程求解（代数方程、微分方程）
 * - sum: 符号求和
 * - ODE 求解器：一阶线性、二阶常系数等
 *
 * 这些命令扩展了符号计算的功能范围。
 */

#include "symbolic/modules/commands/symbolic_commands_internal.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/solver/symbolic_solver.h"
#include "symbolic/calculus/sum/symbolic_sum.h"
#include "core/services/string_utils.h"
#include <vector>
#include <regex>

namespace symbolic_commands {

namespace {
using namespace symbolic_expression_internal;

// Detect ODE type from expression
enum class ODEType {
    UNKNOWN,
    FIRST_ORDER_SEPARABLE,    // y' = f(x)g(y) -> separable
    FIRST_ORDER_LINEAR,       // y' + p(x)y = q(x) -> integrating factor
    SECOND_ORDER_CONSTANT_COEFF, // ay'' + by' + cy = f(x)
    HOMOGENEOUS               // y' = f(y/x)
};

static SymbolicExpression replace_diff_nodes_for_ode(const SymbolicExpression& expr, const std::string& y_var, const std::string& /*x_var*/) {
    if (!expr.node_) return expr;
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "diff") {
        bool is_order_2 = false;
        if (expr.node_->children.size() >= 3) {
            Scalar ord;
            if (SymbolicExpression(expr.node_->children[2]).is_number(&ord)) {
                if (mymath::abs(ord - Scalar(2.0L)) < Scalar(1e-5L)) is_order_2 = true;
            } else if (expr.node_->children[2]->text == "2") {
                is_order_2 = true;
            }
        }
        if (is_order_2) {
            return SymbolicExpression::variable("_y2");
        } else {
            return SymbolicExpression::variable("_y1");
        }
    }
    if (expr.node_->type == NodeType::kVariable) {
        if (expr.node_->text == y_var + "''") return SymbolicExpression::variable("_y2");
        if (expr.node_->text == y_var + "'") return SymbolicExpression::variable("_y1");
        if (expr.node_->text == y_var) return SymbolicExpression::variable("_y0");
    }
    if (!expr.node_->children.empty()) {
        std::vector<SymbolicExpression> new_children;
        for (const auto& ch : expr.node_->children) {
            new_children.push_back(replace_diff_nodes_for_ode(SymbolicExpression(ch), y_var, ""));
        }
        return SymbolicExpression::function(expr.node_->text, new_children);
    }
    if (expr.node_->left && expr.node_->right) {
        auto l = replace_diff_nodes_for_ode(SymbolicExpression(expr.node_->left), y_var, "");
        auto r = replace_diff_nodes_for_ode(SymbolicExpression(expr.node_->right), y_var, "");
        return SymbolicExpression(symbolic_expression_internal::make_binary(expr.node_->type, l.node_, r.node_));
    }
    if (expr.node_->left) {
        auto l = replace_diff_nodes_for_ode(SymbolicExpression(expr.node_->left), y_var, "");
        return SymbolicExpression(symbolic_expression_internal::make_unary(expr.node_->type, l.node_, expr.node_->text));
    }
    return expr;
}

bool solve_first_order_linear(const SymbolicExpression& eq, const std::string& x_var,
                              const std::string& y_var, std::string* output) {
    try {
        SymbolicExpression x = SymbolicExpression::variable(x_var);
        SymbolicExpression y = SymbolicExpression::variable(y_var);
        SymbolicExpression replaced = replace_diff_nodes_for_ode(eq, y_var, x_var).simplify();

        SymbolicExpression coeff_y1 = replaced.derivative("_y1").simplify();
        SymbolicExpression coeff_y0 = replaced.derivative("_y0").simplify();

        SymbolicExpression p, q;

        if (!expr_is_zero(coeff_y1)) {
            // Form: coeff_y1 * y' + coeff_y0 * y + const_term = 0
            // y' + (coeff_y0 / coeff_y1) * y = -const_term / coeff_y1
            SymbolicExpression rem = (replaced - (coeff_y1 * SymbolicExpression::variable("_y1")) - (coeff_y0 * SymbolicExpression::variable("_y0"))).simplify();
            if (!rem.is_constant(y_var) && !rem.is_constant("_y0") && !rem.is_constant("_y1")) {
                return false;
            }
            p = (coeff_y0 / coeff_y1).simplify();
            q = (make_negate(rem) / coeff_y1).simplify();
        } else {
            // Form: y' = rhs(x, y) -> y' - p(x)*y = q(x)
            SymbolicExpression rhs = eq.simplify();
            p = rhs.derivative(y_var).simplify();
            SymbolicExpression p_dy = p.derivative(y_var).simplify();
            if (!expr_is_zero(p_dy)) {
                return false;
            }
            q = (rhs - (p * y)).simplify();
            SymbolicExpression q_dy = q.derivative(y_var).simplify();
            if (!expr_is_zero(q_dy)) {
                return false;
            }
            p = make_negate(p).simplify();
        }

        if (!p.is_constant(y_var) || !q.is_constant(y_var)) {
            return false;
        }

        // Integrating factor: μ(x) = e^(∫p dx)
        SymbolicExpression integral_p = p.integral(x_var).simplify();
        SymbolicExpression mu = make_function("exp", integral_p).simplify();
        SymbolicExpression inv_mu = make_function("exp", make_negate(integral_p).simplify()).simplify();

        // ∫q*μ dx
        SymbolicExpression q_mu = (q * mu).simplify();
        SymbolicExpression integral_q_mu = q_mu.integral(x_var).simplify();

        SymbolicExpression c_const = SymbolicExpression::variable("C");
        SymbolicExpression y_solution = (inv_mu * (integral_q_mu + c_const)).simplify();

        *output = y_var + "(" + x_var + ") = " + y_solution.to_string();
        return true;
    } catch (...) {
        return false;
    }
}

bool solve_second_order_constant_coeff(
    const SymbolicExpression& eq,
    const std::string& x_var,
    const std::string& y_var,
    std::string* output) {

    try {
        SymbolicExpression x = SymbolicExpression::variable(x_var);
        SymbolicExpression replaced = replace_diff_nodes_for_ode(eq, y_var, x_var).simplify();

        // Extract a, b, c
        SymbolicExpression a = replaced.derivative("_y2").simplify();
        SymbolicExpression b = replaced.derivative("_y1").simplify();
        SymbolicExpression c = replaced.derivative("_y0").simplify();

        if (expr_is_zero(a)) {
            return false;
        }

        if (!a.is_constant(x_var) || !b.is_constant(x_var) || !c.is_constant(x_var)) {
            return false; // Not constant coefficients
        }

        // a*r^2 + b*r + c = 0
        SymbolicExpression two_a = (SymbolicExpression::number(2.0L) * a).simplify();
        SymbolicExpression four_ac = (SymbolicExpression::number(4.0L) * a * c).simplify();
        SymbolicExpression discriminant = (b * b - four_ac).simplify();

        SymbolicExpression c1 = SymbolicExpression::variable("C1");
        SymbolicExpression c2 = SymbolicExpression::variable("C2");

        Scalar disc_val = Scalar(0.0L);
        if (discriminant.is_number(&disc_val)) {
            if (mymath::is_near_zero(disc_val, Scalar(1e-12L))) {
                // Repeated root: r = -b / (2a)
                SymbolicExpression r = (make_negate(b) / two_a).simplify();
                SymbolicExpression exp_rx = make_function("exp", (r * x).simplify()).simplify();
                SymbolicExpression y_sol = ((c1 + c2 * x) * exp_rx).simplify();
                *output = y_var + "(" + x_var + ") = " + y_sol.to_string();
                return true;
            } else if (disc_val > Scalar(0.0L)) {
                // Two distinct real roots
                SymbolicExpression sqrt_disc = make_function("sqrt", discriminant).simplify();
                SymbolicExpression r1 = ((make_negate(b) + sqrt_disc) / two_a).simplify();
                SymbolicExpression r2 = ((make_negate(b) - sqrt_disc) / two_a).simplify();
                SymbolicExpression exp1 = make_function("exp", (r1 * x).simplify()).simplify();
                SymbolicExpression exp2 = make_function("exp", (r2 * x).simplify()).simplify();
                SymbolicExpression y_sol = (c1 * exp1 + c2 * exp2).simplify();
                *output = y_var + "(" + x_var + ") = " + y_sol.to_string();
                return true;
            } else {
                // Complex conjugate roots: alpha +- beta * i
                SymbolicExpression abs_disc = make_function("sqrt", make_negate(discriminant).simplify()).simplify();
                SymbolicExpression alpha = (make_negate(b) / two_a).simplify();
                SymbolicExpression beta = (abs_disc / two_a).simplify();

                SymbolicExpression exp_ax = make_function("exp", (alpha * x).simplify()).simplify();
                SymbolicExpression cos_bx = make_function("cos", (beta * x).simplify()).simplify();
                SymbolicExpression sin_bx = make_function("sin", (beta * x).simplify()).simplify();

                SymbolicExpression y_sol = (exp_ax * (c1 * cos_bx + c2 * sin_bx)).simplify();
                *output = y_var + "(" + x_var + ") = " + y_sol.to_string();
                return true;
            }
        } else {
            // General symbolic discriminant
            SymbolicExpression sqrt_disc = make_function("sqrt", discriminant).simplify();
            SymbolicExpression r1 = ((make_negate(b) + sqrt_disc) / two_a).simplify();
            SymbolicExpression r2 = ((make_negate(b) - sqrt_disc) / two_a).simplify();
            SymbolicExpression exp1 = make_function("exp", (r1 * x).simplify()).simplify();
            SymbolicExpression exp2 = make_function("exp", (r2 * x).simplify()).simplify();
            SymbolicExpression y_sol = (c1 * exp1 + c2 * exp2).simplify();
            *output = y_var + "(" + x_var + ") = " + y_sol.to_string();
            return true;
        }
    } catch (...) {
        return false;
    }
}

} // namespace

bool handle_misc_commands(const SymbolicCommandContext& ctx,
                         const std::string& command,
                         const std::string& /*inside*/,
                         const std::vector<std::string>& arguments,
                         std::string* output) {
    if (command == "dsolve") {
        if (arguments.size() < 1) throw std::runtime_error("dsolve expects equation");

        std::string var; SymbolicExpression eq;
        ctx.resolve_symbolic(arguments[0], false, &var, &eq);

        std::string y_var = "y";
        std::string x_var = "x";
        if (arguments.size() >= 3) {
            y_var = trim_copy(arguments[1]);
            x_var = trim_copy(arguments[2]);
        } else if (arguments.size() == 2) {
            std::string arg1 = trim_copy(arguments[1]);
            if (arg1 == "y") {
                y_var = "y";
                x_var = "x";
            } else if (arg1 == "x" || arg1 == "t") {
                y_var = "y";
                x_var = arg1;
            } else {
                y_var = arg1;
                x_var = "x";
            }
        }

        try {
            SymbolicExpression rhs = eq.simplify();
            SymbolicExpression y = SymbolicExpression::variable(y_var);
            SymbolicExpression x = SymbolicExpression::variable(x_var);

            // Try second-order ODE first
            if (solve_second_order_constant_coeff(rhs, x_var, y_var, output)) {
                return true;
            }

            // Try first-order linear ODE
            if (solve_first_order_linear(rhs, x_var, y_var, output)) {
                return true;
            }

            // Fallback: simple separable case y' = f(x)
            // Check if rhs is independent of y
            SymbolicExpression rhs_dy = rhs.derivative(y_var).simplify();
            if (expr_is_zero(rhs_dy)) {
                // y' = f(x) -> y = ∫f(x)dx + C
                SymbolicExpression sol = rhs.integral(x_var).simplify();
                *output = y_var + "(" + x_var + ") = " + sol.to_string() + " + C";
                return true;
            }

            // Check for separable form: y' = f(x)*g(y)
            // Try to separate: dy/g(y) = f(x)dx
            SymbolicExpression rhs_dx = rhs.derivative(x_var).simplify();
            SymbolicExpression rhs_dx_dy = rhs_dx.derivative(y_var).simplify();

            // If ∂²(rhs)/(∂x∂y) = 0, then rhs = f(x) + g(y), not separable product
            // For separable: rhs = f(x)*g(y), check if rhs/(∂rhs/∂x) is independent of x
            if (!expr_is_zero(rhs_dx)) {
                SymbolicExpression ratio = (rhs / rhs_dx).simplify();
                SymbolicExpression ratio_dx = ratio.derivative(x_var).simplify();
                if (expr_is_zero(ratio_dx)) {
                    // Separable! g(y) = rhs/f(x), so ∫dy/g(y) = ∫f(x)dx
                    // ratio = g(y), so 1/g(y) dy = f(x)/g(y) dx = rhs/g(y)² dx
                    // Actually: rhs = f(x)*g(y), rhs_dx = f'(x)*g(y)
                    // ratio = f(x)/f'(x) which should be independent of x for separable

                    // For now, try direct integration approach
                    SymbolicExpression one_over_g = (SymbolicExpression::number(1.0L) / ratio).simplify();
                    SymbolicExpression integral_dy = one_over_g.integral(y_var).simplify();

                    // f(x) = rhs / g(y) = rhs * one_over_g
                    SymbolicExpression f_x = (rhs * one_over_g).simplify();
                    SymbolicExpression integral_dx = f_x.integral(x_var).simplify();

                    *output = "Separable ODE: ∫dy/g(y) = ∫f(x)dx\n" +
                              std::string("∫") + one_over_g.to_string() + " dy = " +
                              integral_dy.to_string() + "\n" +
                              std::string("∫") + f_x.to_string() + " dx = " +
                              integral_dx.to_string() + "\n" +
                              "Solution: " + integral_dy.to_string() + " = " + integral_dx.to_string() + " + C";
                    return true;
                }
            }

            // Unable to solve symbolically
            *output = "Unable to solve ODE symbolically. The equation may require numerical methods or is not in a supported form.\n"
                      "Supported forms:\n"
                      "  - y' = f(x) (direct integration)\n"
                      "  - y' + p(x)y = q(x) (first-order linear, integrating factor)\n"
                      "  - y' = f(x)g(y) (separable)\n"
                      "  - ay'' + by' + cy = 0 (second-order constant coefficients)";
            return true;

        } catch (const std::exception& e) {
            *output = std::string("Error solving ODE: ") + e.what();
            return true;
        } catch (...) {
            return false;
        }
    }

    if (command == "sum") {
        if (arguments.size() < 4) throw std::runtime_error("sum expects expr, var, lower, upper");
        
        SymbolicExpression term;
        std::string var;
        BoundArgument lower, upper;
        
        if (symbolic_sum::parse_sum_arguments(arguments, &term, &var, &lower, &upper)) {
            symbolic_sum::SymbolicSumEngine engine;
            symbolic_sum::SumResult res = engine.compute_sum(term, var, lower, upper);
            *output = symbolic_sum::format_sum_result(res);
            return true;
        }
        return false;
    }

    return false;
}

} // namespace symbolic_commands

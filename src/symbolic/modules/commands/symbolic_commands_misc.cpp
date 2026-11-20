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
// Detect ODE type from expression
enum class ODEType {
    UNKNOWN,
    FIRST_ORDER_SEPARABLE,
    FIRST_ORDER_LINEAR,
    FIRST_ORDER_BERNOULLI,
    SECOND_ORDER_CONSTANT_COEFF,
    HOMOGENEOUS
};

static SymbolicExpression replace_diff_nodes_for_ode(const SymbolicExpression& expr,
                                                     const std::string& y_var,
                                                     const std::string& x_var) {
    if (!expr.node_) return expr;
    if (expr.node_->type == NodeType::kFunction && expr.node_->text == "diff") {
        bool is_y_diff = false;
        if (!expr.node_->children.empty()) {
            SymbolicExpression target(expr.node_->children[0]);
            if (target.is_variable_named(y_var) ||
                (target.node_ && target.node_->type == NodeType::kFunction && target.node_->text == y_var)) {
                is_y_diff = true;
            }
        } else if (expr.node_->left) {
            SymbolicExpression target(expr.node_->left);
            if (target.is_variable_named(y_var) ||
                (target.node_ && target.node_->type == NodeType::kFunction && target.node_->text == y_var)) {
                is_y_diff = true;
            }
        }
        if (is_y_diff) {
            bool is_order_2 = false;
            if (expr.node_->children.size() >= 3) {
                Scalar ord;
                if (SymbolicExpression(expr.node_->children[2]).is_number(&ord)) {
                    if (mymath::abs(ord - Scalar(2.0L)) < Scalar(1e-5L)) is_order_2 = true;
                } else if (expr.node_->children[2]->text == "2") {
                    is_order_2 = true;
                }
            }
            return is_order_2 ? SymbolicExpression::variable("_y2") : SymbolicExpression::variable("_y1");
        }
    }
    if (expr.node_->type == NodeType::kFunction) {
        if (expr.node_->text == y_var + "''") return SymbolicExpression::variable("_y2");
        if (expr.node_->text == y_var + "'") return SymbolicExpression::variable("_y1");
        if (expr.node_->text == y_var) return SymbolicExpression::variable("_y0");
    }
    if (expr.node_->type == NodeType::kVariable) {
        if (expr.node_->text == y_var + "''") return SymbolicExpression::variable("_y2");
        if (expr.node_->text == y_var + "'") return SymbolicExpression::variable("_y1");
        if (expr.node_->text == y_var) return SymbolicExpression::variable("_y0");
    }
    if (!expr.node_->children.empty()) {
        std::vector<SymbolicExpression> new_children;
        for (const auto& ch : expr.node_->children) {
            new_children.push_back(replace_diff_nodes_for_ode(SymbolicExpression(ch), y_var, x_var));
        }
        return SymbolicExpression::function(expr.node_->text, new_children);
    }
    if (expr.node_->left && expr.node_->right) {
        auto l = replace_diff_nodes_for_ode(SymbolicExpression(expr.node_->left), y_var, x_var);
        auto r = replace_diff_nodes_for_ode(SymbolicExpression(expr.node_->right), y_var, x_var);
        return SymbolicExpression(symbolic_expression_internal::make_binary(expr.node_->type, l.node_, r.node_));
    }
    if (expr.node_->left) {
        auto l = replace_diff_nodes_for_ode(SymbolicExpression(expr.node_->left), y_var, x_var);
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
    } catch (const std::exception&) {
        return false;
    }
}

static inline bool is_same_expr(const SymbolicExpression& a, const SymbolicExpression& b) {
    if (!a.node_ && !b.node_) return true;
    if (!a.node_ || !b.node_) return false;
    return expr_is_zero((a - b).simplify());
}

bool solve_bernoulli(const SymbolicExpression& rhs,
                     const std::string& x_var,
                     const std::string& y_var,
                     std::string* output) {
    // Form: y' = A(x)*y + B(x)*y^n (n != 0, 1)
    // => y' - A(x)*y = B(x)*y^n
    try {
        SymbolicExpression y = SymbolicExpression::variable(y_var);
        SymbolicExpression x = SymbolicExpression::variable(x_var);

        if (rhs.node_->type != NodeType::kAdd && rhs.node_->type != NodeType::kSubtract) {
            return false;
        }

        SymbolicExpression term1(rhs.node_->left);
        SymbolicExpression term2(rhs.node_->right);
        if (rhs.node_->type == NodeType::kSubtract) {
            term2 = make_negate(term2).simplify();
        }

        auto extract_power_of_y = [&](const SymbolicExpression& term, Scalar* exp_val, SymbolicExpression* coeff) -> bool {
            if (is_same_expr(term, y)) {
                *exp_val = 1.0L;
                *coeff = SymbolicExpression::number(1.0L);
                return true;
            }
            if (term.is_constant(y_var)) {
                *exp_val = 0.0L;
                *coeff = term;
                return true;
            }
            if (term.node_->type == NodeType::kPower &&
                is_same_expr(SymbolicExpression(term.node_->left), y)) {
                if (SymbolicExpression(term.node_->right).is_number(exp_val)) {
                    *coeff = SymbolicExpression::number(1.0L);
                    return true;
                }
            }
            if (term.node_->type == NodeType::kMultiply) {
                SymbolicExpression l(term.node_->left);
                SymbolicExpression r(term.node_->right);
                if (l.is_constant(y_var)) {
                    if (is_same_expr(r, y)) {
                        *exp_val = 1.0L; *coeff = l; return true;
                    }
                    if (r.node_->type == NodeType::kPower &&
                        is_same_expr(SymbolicExpression(r.node_->left), y) &&
                        SymbolicExpression(r.node_->right).is_number(exp_val)) {
                        *coeff = l; return true;
                    }
                } else if (r.is_constant(y_var)) {
                    if (is_same_expr(l, y)) {
                        *exp_val = 1.0L; *coeff = r; return true;
                    }
                    if (l.node_->type == NodeType::kPower &&
                        is_same_expr(SymbolicExpression(l.node_->left), y) &&
                        SymbolicExpression(l.node_->right).is_number(exp_val)) {
                        *coeff = r; return true;
                    }
                }
            }
            return false;
        };

        Scalar exp1 = 0, exp2 = 0;
        SymbolicExpression coeff1, coeff2;
        if (!extract_power_of_y(term1, &exp1, &coeff1) || !extract_power_of_y(term2, &exp2, &coeff2)) {
            return false;
        }

        SymbolicExpression P_x, Q_x;
        Scalar n = 0;
        if (mymath::is_near_zero(exp1 - 1.0L) && !mymath::is_near_zero(exp2 - 1.0L) && !mymath::is_near_zero(exp2)) {
            P_x = make_negate(coeff1).simplify();
            Q_x = coeff2;
            n = exp2;
        } else if (mymath::is_near_zero(exp2 - 1.0L) && !mymath::is_near_zero(exp1 - 1.0L) && !mymath::is_near_zero(exp1)) {
            P_x = make_negate(coeff2).simplify();
            Q_x = coeff1;
            n = exp1;
        } else {
            return false;
        }

        // v = y^(1-n) => v' + (1-n)*P(x)*v = (1-n)*Q(x)
        Scalar one_minus_n = 1.0L - n;
        SymbolicExpression p_v = (SymbolicExpression::number(one_minus_n) * P_x).simplify();
        SymbolicExpression q_v = (SymbolicExpression::number(one_minus_n) * Q_x).simplify();

        SymbolicExpression int_p_v = p_v.integral(x_var).simplify();
        SymbolicExpression mu_v = make_function("exp", int_p_v).simplify();
        SymbolicExpression inv_mu_v = make_function("exp", make_negate(int_p_v).simplify()).simplify();
        SymbolicExpression int_q_mu = (q_v * mu_v).simplify().integral(x_var).simplify();

        SymbolicExpression c_const = SymbolicExpression::variable("C");
        SymbolicExpression v_sol = (inv_mu_v * (int_q_mu + c_const)).simplify();

        // y = v^(1/(1-n))
        SymbolicExpression y_sol = make_power(v_sol, SymbolicExpression::number(1.0L / one_minus_n)).simplify();
        *output = y_var + "(" + x_var + ") = " + y_sol.to_string();
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool solve_separable_ode(const SymbolicExpression& rhs,
                         const std::string& x_var,
                         const std::string& y_var,
                         std::string* output) {
    try {
        SymbolicExpression x = SymbolicExpression::variable(x_var);
        SymbolicExpression y = SymbolicExpression::variable(y_var);

        // Case 1: rhs independent of y: y' = f(x)
        SymbolicExpression rhs_dy = rhs.derivative(y_var).simplify();
        if (expr_is_zero(rhs_dy)) {
            SymbolicExpression sol = rhs.integral(x_var).simplify();
            *output = y_var + "(" + x_var + ") = " + sol.to_string() + " + C";
            return true;
        }

        // Case 2: rhs independent of x: y' = g(y)
        SymbolicExpression rhs_dx = rhs.derivative(x_var).simplify();
        if (expr_is_zero(rhs_dx)) {
            SymbolicExpression one_over_g = (SymbolicExpression::number(1.0L) / rhs).simplify();
            SymbolicExpression int_y = one_over_g.integral(y_var).simplify();
            if (is_same_expr(rhs, y)) {
                SymbolicExpression sol = make_multiply(SymbolicExpression::variable("C"),
                                                       make_function("exp", x)).simplify();
                *output = y_var + "(" + x_var + ") = " + sol.to_string();
                return true;
            }
            *output = int_y.to_string() + " = " + x_var + " + C";
            return true;
        }

        // Case 3: Check product/quotient separation rhs = f(x) * g(y)
        SymbolicExpression f_part, g_part;
        if (rhs.node_->type == NodeType::kMultiply) {
            SymbolicExpression l(rhs.node_->left);
            SymbolicExpression r(rhs.node_->right);
            if (l.is_constant(y_var) && r.is_constant(x_var)) {
                f_part = l; g_part = r;
            } else if (l.is_constant(x_var) && r.is_constant(y_var)) {
                f_part = r; g_part = l;
            }
        } else if (rhs.node_->type == NodeType::kDivide) {
            SymbolicExpression num(rhs.node_->left);
            SymbolicExpression den(rhs.node_->right);
            if (num.is_constant(y_var) && den.is_constant(x_var)) {
                f_part = num; g_part = (SymbolicExpression::number(1.0L) / den).simplify();
            } else if (num.is_constant(x_var) && den.is_constant(y_var)) {
                f_part = (SymbolicExpression::number(1.0L) / den).simplify(); g_part = num;
            }
        }

        if (!f_part.node_ || !g_part.node_) {
            SymbolicExpression log_deriv_x = (rhs_dx / rhs).simplify();
            SymbolicExpression log_deriv_x_dy = log_deriv_x.derivative(y_var).simplify();
            if (expr_is_zero(log_deriv_x_dy)) {
                SymbolicExpression int_log = log_deriv_x.integral(x_var).simplify();
                f_part = make_function("exp", int_log).simplify();
                g_part = (rhs / f_part).simplify();
            }
        }

        if (f_part.node_ && g_part.node_ && f_part.is_constant(y_var) && g_part.is_constant(x_var)) {
            SymbolicExpression one_over_g = (SymbolicExpression::number(1.0L) / g_part).simplify();
            SymbolicExpression int_y = one_over_g.integral(y_var).simplify();
            SymbolicExpression int_x = f_part.integral(x_var).simplify();

            if (is_same_expr(g_part, y)) {
                SymbolicExpression sol = make_multiply(SymbolicExpression::variable("C"),
                                                       make_function("exp", int_x)).simplify();
                *output = y_var + "(" + x_var + ") = " + sol.to_string();
                return true;
            }

            *output = int_y.to_string() + " = " + int_x.to_string() + " + C";
            return true;
        }

        return false;
    } catch (const std::exception&) {
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

        SymbolicExpression rem = (replaced - (a * SymbolicExpression::variable("_y2"))
                                           - (b * SymbolicExpression::variable("_y1"))
                                           - (c * SymbolicExpression::variable("_y0"))).simplify();
        if (!rem.is_constant("_y0") || !rem.is_constant("_y1") || !rem.is_constant("_y2")) {
            return false;
        }

        SymbolicExpression f_x = (make_negate(rem) / a).simplify();

        // a*r^2 + b*r + c = 0
        SymbolicExpression two_a = (SymbolicExpression::number(2.0L) * a).simplify();
        SymbolicExpression four_ac = (SymbolicExpression::number(4.0L) * a * c).simplify();
        SymbolicExpression discriminant = (b * b - four_ac).simplify();

        SymbolicExpression c1 = SymbolicExpression::variable("C1");
        SymbolicExpression c2 = SymbolicExpression::variable("C2");

        SymbolicExpression y1, y2, wronskian, y_h;

        Scalar disc_val = Scalar(0.0L);
        if (discriminant.is_number(&disc_val)) {
            if (mymath::is_near_zero(disc_val, Scalar(1e-12L))) {
                // Repeated root: r = -b / (2a)
                SymbolicExpression r = (make_negate(b) / two_a).simplify();
                SymbolicExpression exp_rx = make_function("exp", (r * x).simplify()).simplify();
                y1 = exp_rx;
                y2 = (x * exp_rx).simplify();
                wronskian = make_function("exp", ((SymbolicExpression::number(2.0L) * r).simplify() * x).simplify()).simplify();
                y_h = ((c1 + c2 * x) * exp_rx).simplify();
            } else if (disc_val > Scalar(0.0L)) {
                // Two distinct real roots
                SymbolicExpression sqrt_disc = make_function("sqrt", discriminant).simplify();
                SymbolicExpression r1 = ((make_negate(b) + sqrt_disc) / two_a).simplify();
                SymbolicExpression r2 = ((make_negate(b) - sqrt_disc) / two_a).simplify();
                y1 = make_function("exp", (r1 * x).simplify()).simplify();
                y2 = make_function("exp", (r2 * x).simplify()).simplify();
                wronskian = ((r2 - r1).simplify() * make_function("exp", ((r1 + r2).simplify() * x).simplify())).simplify();
                y_h = (c1 * y1 + c2 * y2).simplify();
            } else {
                // Complex conjugate roots: alpha +- beta * i
                SymbolicExpression abs_disc = make_function("sqrt", make_negate(discriminant).simplify()).simplify();
                SymbolicExpression alpha = (make_negate(b) / two_a).simplify();
                SymbolicExpression beta = (abs_disc / two_a).simplify();

                SymbolicExpression exp_ax = make_function("exp", (alpha * x).simplify()).simplify();
                SymbolicExpression cos_bx = make_function("cos", (beta * x).simplify()).simplify();
                SymbolicExpression sin_bx = make_function("sin", (beta * x).simplify()).simplify();
                y1 = (exp_ax * cos_bx).simplify();
                y2 = (exp_ax * sin_bx).simplify();
                wronskian = (beta * make_function("exp", ((SymbolicExpression::number(2.0L) * alpha).simplify() * x).simplify())).simplify();

                y_h = (exp_ax * (c1 * cos_bx + c2 * sin_bx)).simplify();
            }
        } else {
            SymbolicExpression sqrt_disc = make_function("sqrt", discriminant).simplify();
            SymbolicExpression r1 = ((make_negate(b) + sqrt_disc) / two_a).simplify();
            SymbolicExpression r2 = ((make_negate(b) - sqrt_disc) / two_a).simplify();
            y1 = make_function("exp", (r1 * x).simplify()).simplify();
            y2 = make_function("exp", (r2 * x).simplify()).simplify();
            wronskian = ((r2 - r1).simplify() * make_function("exp", ((r1 + r2).simplify() * x).simplify())).simplify();
            y_h = (c1 * y1 + c2 * y2).simplify();
        }

        if (expr_is_zero(f_x)) {
            *output = y_var + "(" + x_var + ") = " + y_h.to_string();
            return true;
        }

        // Non-homogeneous: Variation of parameters
        try {
            SymbolicExpression integrand1 = (make_negate(y2 * f_x) / wronskian).simplify();
            SymbolicExpression integrand2 = ((y1 * f_x) / wronskian).simplify();
            SymbolicExpression u1 = integrand1.integral(x_var).simplify();
            SymbolicExpression u2 = integrand2.integral(x_var).simplify();
            SymbolicExpression y_p = (u1 * y1 + u2 * y2).simplify();
            SymbolicExpression y_total = (y_h + y_p).simplify();
            *output = y_var + "(" + x_var + ") = " + y_total.to_string();
            return true;
        } catch (const std::exception&) {
            *output = y_var + "(" + x_var + ") = " + y_h.to_string() + " + y_p(" + x_var + ")";
            return true;
        }
    } catch (const std::exception&) {
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

            // 1. Try second-order ODE
            if (solve_second_order_constant_coeff(rhs, x_var, y_var, output)) {
                return true;
            }

            // 2. Try first-order linear ODE
            if (solve_first_order_linear(rhs, x_var, y_var, output)) {
                return true;
            }

            // 3. Try Bernoulli ODE
            if (solve_bernoulli(rhs, x_var, y_var, output)) {
                return true;
            }

            // 4. Try separable ODE
            if (solve_separable_ode(rhs, x_var, y_var, output)) {
                return true;
            }

            // Unable to solve symbolically
            *output = "Unable to solve ODE symbolically. The equation may require numerical methods or is not in a supported form.\n"
                      "Supported forms:\n"
                      "  - y' = f(x) (direct integration)\n"
                      "  - y' + p(x)y = q(x) (first-order linear, integrating factor)\n"
                      "  - y' + p(x)y = q(x)y^n (Bernoulli ODE)\n"
                      "  - y' = f(x)g(y) (separable)\n"
                      "  - ay'' + by' + cy = f(x) (second-order constant coefficients, homogeneous/non-homogeneous)";
            return true;

        } catch (const std::exception& e) {
            *output = std::string("Error solving ODE: ") + e.what();
            return true;
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

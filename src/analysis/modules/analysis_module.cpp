// ============================================================================
// 函数分析模块命令处理
// ============================================================================
//
// 本文件是函数分析模块的入口，负责：
// - 解析命令参数
// - 调用 calculus/ 子目录中的核心算法
// - 格式化输出结果
//
// 核心计算已拆分到：
// - function_analysis.cpp: 函数分析器
// - limit_solver.cpp: 极限计算
// - critical_point_solver.cpp: 临界点求解
// - analysis_command_helpers.cpp: 辅助函数

#include "analysis/modules/analysis_module.h"
#include "execution/engine/script_context.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "analysis/base/precision_constants.h"
#include "analysis/calculus/analysis_command_helpers.h"
#include "analysis/calculus/critical_point_solver.h"
#include "symbolic/modules/symbolic_module.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/algebra/groebner/groebner_basis.h"
#include "symbolic/solver/symbolic_solver.h"
#include "analysis/calculus/function_analysis.h"
#include "parser/grammars/unified_expression_parser.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "math/mymath.h"
#include "math/helpers/integer_helpers.h"
#include <algorithm>
#include <sstream>
#include <iterator>

namespace analysis_cmds {

using namespace utils;
using Scalar = mymath::Scalar;

bool is_analysis_command(const std::string& command) {
    return command == "limit" || command == "critical" || command == "extrema" || command == "lagrange";
}

bool handle_analysis_command(const AnalysisContext& ctx,
                             const std::string& command,
                             const std::vector<std::string>& arguments,
                             std::string* output) {
    if (command == "limit") {
        if (arguments.size() < 2) throw std::runtime_error("limit expects expr, point");
        std::string point_arg; size_t dir_idx; FunctionAnalysis analysis;
        bool explicit_var = arguments.size() >= 3 && is_identifier_text(utils::trim_copy(arguments[1])) && !is_infinity_literal(arguments[1]);
        if (explicit_var) {
            std::string inf; SymbolicExpression expr; ctx.resolve_symbolic(arguments[0], false, &inf, &expr);
            analysis = FunctionAnalysis(utils::trim_copy(arguments[1])); analysis.define(expr.to_string());
            point_arg = arguments[2]; dir_idx = 3;
        } else {
            analysis = ctx.build_analysis(arguments[0]); point_arg = arguments[1]; dir_idx = 2;
        }
        const std::string trimmed_point = utils::trim_copy(point_arg);
        const bool infinite_point = is_infinity_literal(trimmed_point);
        const bool negative_infinity =
            infinite_point && !trimmed_point.empty() && trimmed_point.front() == '-';
        int dir = 0;
        if (arguments.size() > dir_idx) {
            dir = static_cast<int>(round_to_long_long(ctx.parse_decimal(arguments[dir_idx])));
        } else if (infinite_point) {
            dir = negative_infinity ? -1 : 1;
        }
        const Scalar point_value = infinite_point
            ? (negative_infinity ? -Scalar::infinity() : Scalar::infinity())
            : ctx.parse_decimal(point_arg);
        Scalar limit_value = 0.0L;
        try {
            limit_value = analysis.limit(point_value, dir);
        } catch (const std::runtime_error& ex) {
            const std::string message = ex.what();
            if (message.find("limit does not exist") != std::string::npos) {
                throw std::runtime_error("limit did not converge");
            }
            throw;
        }
        
        if (mymath::isnan(limit_value)) {
            throw std::runtime_error("limit did not converge");
        }
        
        if (!mymath::isfinite(limit_value)) {
            *output = (limit_value > 0.0L) ? "inf" : "-inf";
        } else {
            *output = format_decimal(ctx.normalize_result(limit_value));
        }
        return true;
    }

    if (command == "extrema") {
        if (arguments.size() != 3) throw std::runtime_error("extrema expects expression, a, b");
        FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
        const std::vector<ExtremumPoint> points = analysis.solve_extrema(ctx.parse_decimal(arguments[1]), ctx.parse_decimal(arguments[2]));
        if (points.empty()) { *output = "No extrema found."; return true; }
        std::ostringstream out;
        for (size_t i = 0; i < points.size(); ++i) {
            if (i != 0) out << '\n';
            out << (points[i].is_maximum ? "max" : "min") << ": x = " << format_decimal(points[i].x) << ", f(x) = " << format_decimal(points[i].value);
        }
        *output = out.str(); return true;
    }

    if (command == "critical") {
        if (arguments.empty()) throw std::runtime_error("critical expects expression and optional variables");
        std::string var;
        SymbolicExpression expr;
        ctx.resolve_symbolic(arguments[0], false, &var, &expr);
        std::vector<std::string> variables = ctx.parse_symbolic_variable_arguments(arguments, 1, expr.identifier_variables());
        if (variables.empty()) variables = {var};

        std::vector<SymbolicExpression> gradient;
        for (const auto& v : variables) {
            gradient.push_back(expr.derivative(v).simplify());
        }

        std::vector<std::map<std::string, Scalar>> critical_points;

        auto eval_ast = [](auto& self, const std::shared_ptr<SymbolicExpression::Node>& n, const std::map<std::string, Scalar>& pt) -> Scalar {
            if (!n) return 0.0L;
            switch (n->type) {
                case NodeType::kNumber: return Scalar(n->number_value);
                case NodeType::kVariable: {
                    auto it = pt.find(n->text);
                    if (it != pt.end()) return it->second;
                    throw std::runtime_error("unknown variable: " + n->text);
                }
                case NodeType::kPi: return Scalar(mymath::pi());
                case NodeType::kE: return Scalar(mymath::e());
                case NodeType::kAdd: return self(self, n->left, pt) + self(self, n->right, pt);
                case NodeType::kSubtract: return self(self, n->left, pt) - self(self, n->right, pt);
                case NodeType::kMultiply: return self(self, n->left, pt) * self(self, n->right, pt);
                case NodeType::kDivide: return self(self, n->left, pt) / self(self, n->right, pt);
                case NodeType::kPower: return mymath::pow(self(self, n->left, pt), self(self, n->right, pt));
                case NodeType::kNegate: return -self(self, n->left, pt);
                case NodeType::kFunction: {
                    Scalar arg = self(self, n->left, pt);
                    if (n->text == "sin") return mymath::sin(arg);
                    if (n->text == "cos") return mymath::cos(arg);
                    if (n->text == "tan") return mymath::tan(arg);
                    if (n->text == "exp") return mymath::exp(arg);
                    if (n->text == "ln" || n->text == "log") return mymath::log(arg);
                    if (n->text == "sqrt") return mymath::sqrt(arg);
                    throw std::runtime_error("unsupported function in numeric eval: " + n->text);
                }
                default: throw std::runtime_error("unsupported node type in numeric eval");
            }
        };

        if (variables.size() == 1) {
            const std::string& variable = variables[0];
            const SymbolicExpression derivative = gradient[0];
            auto eval_derivative = [&](Scalar x) {
                try {
                    return eval_ast(eval_ast, derivative.node_, {{variable, x}});
                } catch (...) {
                    throw std::runtime_error("derivative is not numeric at this point");
                }
            };
            auto add_point = [&](Scalar x) {
                for (const auto& existing : critical_points) {
                    const auto it = existing.find(variable);
                    if (it != existing.end() && mymath::abs(it->second - x) < Scalar(1e-5L)) {
                        return;
                    }
                }
                critical_points.push_back({{variable, Scalar(ctx.normalize_result(static_cast<double>(x)))}});
            };

            const Scalar scan_min = Scalar(-100);
            const Scalar scan_max = Scalar(100);
            const int coarse_segments = 512;

            SymbolicExpression second_deriv = derivative.derivative(variable).simplify();
            auto eval_second = [&](Scalar x) {
                try {
                    return eval_ast(eval_ast, second_deriv.node_, {{variable, x}});
                } catch (...) {
                    return Scalar(0.0L);
                }
                return Scalar(0);
            };

            Scalar previous_x = scan_min;
            Scalar previous_value = eval_derivative(previous_x);
            for (int i = 1; i <= coarse_segments; ++i) {
                const Scalar current_x = scan_min + (scan_max - scan_min) * Scalar(i) / Scalar(coarse_segments);
                const Scalar current_value = eval_derivative(current_x);
                if ((previous_value < Scalar(0) && current_value > Scalar(0)) ||
                    (previous_value > Scalar(0) && current_value < Scalar(0))) {
                    Scalar left = previous_x;
                    Scalar right = current_x;
                    Scalar left_value = previous_value;
                    for (int iter = 0; iter < 80; ++iter) {
                        const Scalar mid = (left + right) * Scalar(0.5L);
                        const Scalar mid_value = eval_derivative(mid);
                        if (mymath::isfinite(mid_value) && mymath::abs(mid_value) < precision::newton_tolerance<Scalar>()) {
                            left = right = mid;
                            break;
                        }
                        if ((left_value < Scalar(0) && mid_value > Scalar(0)) ||
                            (left_value > Scalar(0) && mid_value < Scalar(0))) {
                            right = mid;
                        } else {
                            left = mid;
                            left_value = mid_value;
                        }
                    }
                    add_point((left + right) * Scalar(0.5L));
                }
                previous_x = current_x;
                previous_value = current_value;
            }

            for (int i = 0; i <= coarse_segments; ++i) {
                const Scalar x = scan_min + (scan_max - scan_min) * Scalar(i) / Scalar(coarse_segments);
                const Scalar deriv_val = eval_derivative(x);
                const Scalar second_val = eval_second(x);

                if (mymath::abs(deriv_val) < precision::gradient_convergence_threshold<Scalar>() && mymath::abs(second_val) > precision::sqrt_epsilon<Scalar>()) {
                    Scalar refined_x = x;
                    for (int iter = 0; iter < 20; ++iter) {
                        const Scalar f_prime = eval_derivative(refined_x);
                        const Scalar f_double_prime = eval_second(refined_x);
                        if (mymath::isfinite(f_prime) && mymath::abs(f_prime) < precision::newton_tolerance<Scalar>()) break;
                        if (mymath::abs(f_double_prime) < precision::gradient_convergence_threshold<Scalar>()) break;
                        refined_x = refined_x - f_prime / f_double_prime;
                    }
                    if (mymath::isfinite(eval_derivative(refined_x)) && mymath::abs(eval_derivative(refined_x)) < precision::newton_tolerance<Scalar>() &&
                        mymath::abs(eval_second(refined_x)) > precision::sqrt_epsilon<Scalar>()) {
                        add_point(refined_x);
                    }
                }
            }
        } else {
            auto eval_gradient_at = [&](const SymbolicExpression& g,
                                        const std::map<std::string, Scalar>& point) {
                return eval_ast(eval_ast, g.node_, point);
            };

            auto gradient_norm_at = [&](const std::map<std::string, Scalar>& point) {
                Scalar norm = Scalar(0);
                for (const auto& g : gradient) {
                    const Scalar val = eval_gradient_at(g, point);
                    norm = norm + val * val;
                }
                return mymath::sqrt(norm);
            };

            auto add_critical_point = [&](const std::map<std::string, Scalar>& point) {
                bool duplicate = false;
                for (const auto& existing : critical_points) {
                    bool same = true;
                    for (const auto& v : variables) {
                        const auto it_existing = existing.find(v);
                        const auto it_current = point.find(v);
                        if (it_existing == existing.end() || it_current == point.end() ||
                            mymath::abs(it_existing->second - it_current->second) > Scalar(1e-2L)) {
                            same = false;
                            break;
                        }
                    }
                    if (same) { duplicate = true; break; }
                }
                if (!duplicate) {
                    std::map<std::string, Scalar> normalized;
                    for (const auto& [k, v] : point) {
                        normalized[k] = Scalar(ctx.normalize_result(static_cast<double>(v)));
                    }
                    critical_points.push_back(normalized);
                }
            };

            std::vector<std::map<std::string, Scalar>> starting_points;

            auto is_polynomial_node = [&](const std::shared_ptr<SymbolicExpression::Node>& node) {
                std::function<bool(const std::shared_ptr<SymbolicExpression::Node>&)> check = [&](const std::shared_ptr<SymbolicExpression::Node>& n) {
                    if (!n) return true;
                    switch (n->type) {
                        case NodeType::kNumber:
                        case NodeType::kPi:
                        case NodeType::kE:
                        case NodeType::kVariable:
                            return true;
                        case NodeType::kAdd:
                        case NodeType::kSubtract:
                        case NodeType::kMultiply:
                            return check(n->left) && check(n->right);
                        case NodeType::kNegate:
                            return check(n->left);
                        case NodeType::kPower: {
                            Scalar val;
                            if (SymbolicExpression(n->right).is_number(&val)) {
                                if (val >= Scalar(0) && mymath::abs(val - mymath::floor(val)) < precision::epsilon<Scalar>() * Scalar(100)) {
                                    return check(n->left);
                                }
                            }
                            return false;
                        }
                        default: return false;
                    }
                };
                return check(node);
            };

            bool is_poly_system = true;
            for (const auto& g : gradient) {
                if (!is_polynomial_node(g.node_)) {
                    is_poly_system = false; break;
                }
            }

            if (is_poly_system && variables.size() > 1) {
                try {
                    auto basis = symbolic_groebner::compute_groebner_basis(gradient, variables);
                    if (!basis.empty()) {
                        const SymbolicExpression& last_poly = basis.back();
                        symbolic_solver::SymbolicSolver solver;
                        auto sol = solver.solve(last_poly, variables.back());
                        for (const auto& val_expr : sol.values) {
                            Scalar val;
                            if (val_expr.is_number(&val)) {
                                std::map<std::string, Scalar> pt;
                                pt[variables.back()] = val;
                                for (size_t i = 0; i < variables.size() - 1; ++i) pt[variables[i]] = Scalar(0);
                                starting_points.push_back(pt);
                                std::map<std::string, Scalar> pt1;
                                pt1[variables.back()] = val;
                                for (size_t i = 0; i < variables.size() - 1; ++i) pt1[variables[i]] = Scalar(1);
                                starting_points.push_back(pt1);
                            }
                        }
                    }
                } catch (...) {
                }
            }

            std::map<std::string, Scalar> origin;
            for (const auto& v : variables) origin[v] = Scalar(0);
            starting_points.push_back(origin);

            const std::vector<Scalar> grid_values = {Scalar(-10), Scalar(-1), Scalar(1), Scalar(10)};
            if (variables.size() <= 3) {
                for (Scalar v0 : grid_values) {
                    for (Scalar v1 : grid_values) {
                        std::map<std::string, Scalar> pt;
                        pt[variables[0]] = v0;
                        pt[variables[1]] = v1;
                        if (variables.size() == 3) {
                            for (Scalar v2 : grid_values) {
                                pt[variables[2]] = v2;
                                starting_points.push_back(pt);
                            }
                        } else {
                            starting_points.push_back(pt);
                        }
                    }
                }
            }

            for (int i = 0; i < 5; ++i) {
                std::map<std::string, Scalar> pt;
                for (const auto& v : variables) {
                    pt[v] = Scalar(static_cast<long double>(rand()) / RAND_MAX * 20.0L - 10.0L);
                }
                starting_points.push_back(pt);
            }

            std::vector<std::vector<SymbolicExpression>> symbolic_hessian(variables.size(),
                std::vector<SymbolicExpression>(variables.size()));
            for (std::size_t i = 0; i < variables.size(); ++i) {
                for (std::size_t j = 0; j < variables.size(); ++j) {
                    symbolic_hessian[i][j] = gradient[i].derivative(variables[j]).simplify();
                }
            }

            for (const auto& start : starting_points) {
                try {
                    std::map<std::string, Scalar> current = start;
                    if (gradient_norm_at(current) < precision::gradient_convergence_threshold<Scalar>()) {
                        add_critical_point(current);
                        continue;
                    }

                    for (int iter = 0; iter < 100; ++iter) {
                        std::vector<Scalar> rhs(variables.size(), Scalar(0));
                        std::vector<std::vector<Scalar>> jac(variables.size(),
                            std::vector<Scalar>(variables.size(), Scalar(0)));

                        Scalar current_norm = gradient_norm_at(current);
                        if (current_norm < precision::newton_tolerance<Scalar>()) break;

                        bool eval_ok = true;
                        for (std::size_t row = 0; row < variables.size(); ++row) {
                            rhs[row] = Scalar(0) - eval_gradient_at(gradient[row], current);
                            for (std::size_t col = 0; col < variables.size(); ++col) {
                                try {
                                    jac[row][col] = eval_ast(eval_ast, symbolic_hessian[row][col].node_, current);
                                } catch (...) {
                                    eval_ok = false;
                                    break;
                                }
                            }
                            if (!eval_ok) break;
                        }

                        if (!eval_ok) break;

                        const std::vector<Scalar> delta = solve_linear_system(jac, rhs);

                        Scalar alpha = 1.0L;
                        bool line_search_success = false;
                        for (int ls_iter = 0; ls_iter < 10; ++ls_iter) {
                            std::map<std::string, Scalar> next = current;
                            for (std::size_t i = 0; i < variables.size(); ++i) {
                                next[variables[i]] = next[variables[i]] + alpha * delta[i];
                            }

                            if (gradient_norm_at(next) < current_norm) {
                                current = next;
                                line_search_success = true;
                                break;
                            }
                            alpha *= 0.5L;
                        }

                        if (!line_search_success) {
                            for (std::size_t i = 0; i < variables.size(); ++i) {
                                current[variables[i]] = current[variables[i]] + alpha * delta[i];
                            }
                            if (alpha < precision::line_search_min_step<Scalar>()) break;
                        }

                        Scalar max_delta = 0.0L;
                        for (auto d : delta) max_delta = mymath::fmax(max_delta, mymath::abs(d * alpha));
                        if (max_delta < precision::line_search_min_step<Scalar>()) break;
                    }

                    const Scalar grad_norm = gradient_norm_at(current);
                    if (grad_norm < precision::gradient_convergence_threshold<Scalar>()) {
                        add_critical_point(current);
                    }
                } catch (const std::exception&) {
                }
            }
        }

        if (critical_points.empty()) {
            *output = "No critical points found.";
            return true;
        }

        std::ostringstream out;
        for (size_t i = 0; i < critical_points.size(); ++i) {
            const auto& pt = critical_points[i];
            if (i > 0) out << "\n";

            out << "[";
            bool first = true;
            for (const auto& v : variables) {
                if (!first) out << ", ";
                first = false;
                auto it = pt.find(v);
                if (it != pt.end()) {
                    out << v << " = " << format_decimal((it->second));
                }
            }
            out << "]";

            if (variables.size() == 1) {
                SymbolicExpression second_deriv = expr.derivative(variables[0]).derivative(variables[0]).simplify();
                Scalar second_val = 0.0L;
                try {
                    second_val = eval_ast(eval_ast, second_deriv.node_, pt);
                } catch (...) {}
                const Scalar hessian_threshold = precision::positive_definite_threshold<Scalar>();
                if (second_val > hessian_threshold) out << " (local min)";
                else if (second_val < -hessian_threshold) out << " (local max)";
                else out << " (inflection)";
            } else {
                auto hessian = expr.hessian(variables);

                std::vector<std::vector<Scalar>> hessian_values(hessian.size(), std::vector<Scalar>(hessian.size(), Scalar(0)));
                bool hessian_evaluable = true;
                for (size_t r = 0; r < hessian.size(); ++r) {
                    for (size_t c = 0; c < hessian[r].size(); ++c) {
                        try {
                            hessian_values[r][c] = eval_ast(eval_ast, hessian[r][c].node_, pt);
                        } catch (...) {
                            hessian_evaluable = false;
                            break;
                        }
                    }
                }

                if (hessian_evaluable && !hessian_values.empty()) {
                    bool positive_definite = true;
                    bool negative_definite = true;

                    const Scalar hessian_threshold = precision::positive_definite_threshold<Scalar>();
                    for (size_t i = 0; i < hessian_values.size(); ++i) {
                        if (hessian_values[i][i] <= hessian_threshold) positive_definite = false;
                        if (hessian_values[i][i] >= -hessian_threshold) negative_definite = false;
                    }

                    if (hessian_values.size() == 2) {
                        Scalar det = hessian_values[0][0] * hessian_values[1][1] - hessian_values[0][1] * hessian_values[1][0];
                        if (det <= hessian_threshold) positive_definite = false;
                        if (det <= hessian_threshold) negative_definite = false;
                    }

                    if (positive_definite) {
                        out << " (local min)";
                    } else if (negative_definite) {
                        out << " (local max)";
                    } else {
                        bool all_zero = true;
                        const Scalar zero_threshold = precision::positive_definite_threshold<Scalar>();
                        for (const auto& row : hessian_values) {
                            for (const Scalar& val : row) {
                                if (mymath::abs(val) > zero_threshold) {
                                    all_zero = false;
                                    break;
                                }
                            }
                        }
                        if (all_zero) {
                            out << " (degenerate)";
                        } else {
                            out << " (saddle)";
                        }
                    }
                } else {
                    out << " (degenerate)";
                }
            }
        }
        *output = out.str();
        return true;
    }

    if (command == "lagrange") {
        if (arguments.size() < 2) throw std::runtime_error("lagrange expects f, [g1, g2, ...], [vars]");
        SymbolicExpression f = SymbolicExpression::parse(arguments[0]);
        std::vector<SymbolicExpression> constraints = symbolic_commands::parse_symbolic_expression_list(arguments[1], nullptr);
        std::vector<std::string> variables = ctx.parse_symbolic_variable_arguments(arguments, 2, f.identifier_variables());
        SymbolicExpression lagrangian = f;
        std::vector<std::string> all_vars = variables;
        for (std::size_t i = 0; i < constraints.size(); ++i) {
            std::string lambda_var = "L" + std::to_string(i + 1);
            all_vars.push_back(lambda_var);
            lagrangian = (lagrangian - SymbolicExpression::variable(lambda_var) * constraints[i]).simplify();
        }
        std::string all_vars_str;
        for (std::size_t i = 0; i < all_vars.size(); ++i) { if (i > 0) all_vars_str += ", "; all_vars_str += all_vars[i]; }
        return handle_analysis_command(ctx, "critical", {lagrangian.to_string(), all_vars_str}, output);
    }
    
    return false;
}


std::string AnalysisModule::execute_args(const std::string& command,
                                        const std::vector<std::string>& args,
                                        ::ServiceLocator& locator) {
    auto engine = locator.resolve<IEvaluationEngine>();
    AnalysisContext ctx;
    ctx.resolve_symbolic = [&locator](const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) {
        locator.resolve<IExecutionContext>()->services().symbolic.resolve_symbolic(arg, req, var, expr);
    };
    ctx.parse_symbolic_variable_arguments = [engine](const std::vector<std::string>& arguments, std::size_t start_index, const std::vector<std::string>& defaults) {
        return engine->parse_symbolic_vars(arguments, start_index, defaults);
    };
    ctx.parse_decimal = [engine](const std::string& expr) {
        return engine->parse_decimal(expr);
    };
    ctx.normalize_result = [engine](Scalar value) {
        return engine->normalize_result(value);
    };
    ctx.build_analysis = [&locator, engine](const std::string& expression) {
        SymbolicExpression expr; std::string var;
        locator.resolve<IExecutionContext>()->services().symbolic.resolve_symbolic(expression, false, &var, &expr);
        FunctionAnalysis analysis(var);
        analysis.define(expr.to_string());
        analysis.set_evaluator(engine->build_scoped_evaluator(expr.to_string()));
        return analysis;
    };
    std::string out;
    if (handle_analysis_command(ctx, command, args, &out)) return out;
    throw std::runtime_error("Unknown analysis command: " + command);
}

std::vector<std::string> AnalysisModule::get_commands() const {
    return {"limit", "extrema", "critical", "lagrange"};
}

std::string AnalysisModule::get_help_snippet(const std::string& topic) const {
    if (topic == "analysis") {
        return "Function Analysis:\n"
               "  limit(f, x0, [dir])    Limit of f at x=x0\n"
               "  critical(f, [vars])    Find and classify critical points\n"
               "  extrema(f, a, b)       Global and local extrema in [a, b]\n"
               "  lagrange(f, g, [vars]) Constrained optimization";
    }
    return "";
}

}  // namespace analysis_cmds

#include "module/module_registration.h"
REGISTER_CALCULATOR_MODULE(analysis_cmds::AnalysisModule)

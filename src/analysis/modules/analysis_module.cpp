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
#include "math/numeric/precision/tolerances.h"
#include "analysis/calculus/analysis_command_helpers.h"
#include "analysis/calculus/critical_point_solver.h"
#include "symbolic/modules/symbolic_module.h"
#include "symbolic/core/symbolic_expression.h"
#include "symbolic/public/symbolic_node_types.h"
#include "symbolic/algebra/groebner/groebner_basis.h"
#include "symbolic/solver/symbolic_solver.h"
#include "symbolic/calculus/limit/symbolic_limit.h"
#include "analysis/calculus/function_analysis.h"
#include "parser/grammars/unified_expression_parser.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "math/mymath.h"
#include "math/functions/integer/integer_helpers.h"
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
            std::string dir_str = utils::trim_copy(arguments[dir_idx]);
            if (dir_str == "+" || dir_str == "right" || dir_str == ">" || dir_str == "1") {
                dir = 1;
            } else if (dir_str == "-" || dir_str == "left" || dir_str == "<" || dir_str == "-1") {
                dir = -1;
            } else if (dir_str == "0" || dir_str == "both") {
                dir = 0;
            } else {
                try {
                    dir = static_cast<int>(round_to_long_long(ctx.parse_decimal(dir_str)));
                } catch (...) {
                    dir = 0;
                }
            }
        } else if (infinite_point) {
            dir = negative_infinity ? -1 : 1;
        }
        const Scalar point_value = infinite_point
            ? (negative_infinity ? -Scalar::infinity() : Scalar::infinity())
            : ctx.parse_decimal(point_arg);

        // Prefer the symbolic engine when it can prove a numeric result or
        // prove that the two-sided limit does not exist.  Keep the numerical
        // analyzer as a fallback for symbolic results such as `e` and for
        // expressions outside the symbolic engine's supported rules.
        bool symbolic_nonexistent = false;
        try {
            std::string symbolic_default_var;
            SymbolicExpression symbolic_expr;
            ctx.resolve_symbolic(arguments[0], false, &symbolic_default_var, &symbolic_expr);
            const std::string symbolic_var = explicit_var
                ? utils::trim_copy(arguments[1])
                : analysis.variable_name();
            const BoundArgument bound = infinite_point
                ? (negative_infinity ? BoundArgument::neg_inf() : BoundArgument::pos_inf())
                : BoundArgument::finite(point_value);
            symbolic_limit::SymbolicLimitEngine symbolic_engine;
            const auto symbolic_result = symbolic_engine.compute_limit(
                symbolic_expr, symbolic_var, bound, dir);
            if (!symbolic_result.is_definite) {
                if (symbolic_result.is_oscillating) {
                    symbolic_nonexistent = true;
                }
            } else if (symbolic_result.is_infinite) {
                *output = symbolic_result.value.to_string();
                return true;
            } else {
                Scalar sym_num = Scalar(0);
                if (symbolic_result.value.is_number(&sym_num)) {
                    *output = format_decimal(ctx.normalize_result(sym_num));
                    return true;
                } else if (symbolic_result.value.node_type() == NodeType::kE) {
                    *output = format_decimal(ctx.normalize_result(mymath::e()));
                    return true;
                } else if (symbolic_result.value.node_type() == NodeType::kPi) {
                    *output = format_decimal(ctx.normalize_result(mymath::pi()));
                    return true;
                } else if (symbolic_result.value.has_node()) {
                    *output = symbolic_result.value.to_string();
                    return true;
                }
            }
        } catch (...) {
            // Fall through to the numerical analyzer.
        }
        if (symbolic_nonexistent) {
            throw std::runtime_error("limit did not converge");
        }

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

        if (variables.size() == 1) {
            const std::string& variable = variables[0];
            const SymbolicExpression derivative = gradient[0];
            std::vector<Scalar> pts = analysis::find_univariate_critical_points(derivative, variable, ctx.normalize_result);
            for (Scalar p : pts) {
                critical_points.push_back({{variable, p}});
            }
        } else {
            critical_points = analysis::find_multivariate_critical_points(gradient, variables, ctx.normalize_result);
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
                    second_val = analysis::evaluate_ast_node(second_deriv.node_, pt);
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
                            hessian_values[r][c] = analysis::evaluate_ast_node(hessian[r][c].node_, pt);
                        } catch (...) {
                            hessian_evaluable = false;
                            break;
                        }
                    }
                }

                if (hessian_evaluable && !hessian_values.empty()) {
                    bool positive_definite = true;
                    bool negative_definite = true;
                    const size_t n_dim = hessian_values.size();

                    auto compute_minor = [](const std::vector<std::vector<Scalar>>& mat, size_t k) -> Scalar {
                        std::vector<std::vector<Scalar>> sub(k, std::vector<Scalar>(k));
                        for (size_t r = 0; r < k; ++r) {
                            for (size_t c = 0; c < k; ++c) {
                                sub[r][c] = mat[r][c];
                            }
                        }
                        Scalar det = Scalar(1);
                        for (size_t col = 0; col < k; ++col) {
                            size_t pivot = col;
                            for (size_t r = col + 1; r < k; ++r) {
                                if (mymath::abs(sub[r][col]) > mymath::abs(sub[pivot][col])) {
                                    pivot = r;
                                }
                            }
                            if (mymath::abs(sub[pivot][col]) < Scalar(1e-15L)) return Scalar(0);
                            if (pivot != col) {
                                std::swap(sub[pivot], sub[col]);
                                det = -det;
                            }
                            det *= sub[col][col];
                            for (size_t r = col + 1; r < k; ++r) {
                                Scalar factor = sub[r][col] / sub[col][col];
                                for (size_t c = col; c < k; ++c) {
                                    sub[r][c] -= factor * sub[col][c];
                                }
                            }
                        }
                        return det;
                    };

                    const Scalar hessian_threshold = precision::positive_definite_threshold<Scalar>();
                    for (size_t k = 1; k <= n_dim; ++k) {
                        Scalar minor_k = compute_minor(hessian_values, k);
                        if (minor_k <= hessian_threshold) positive_definite = false;
                        Scalar sign_k = (k % 2 == 1) ? Scalar(-1.0L) : Scalar(1.0L);
                        if (sign_k * minor_k <= hessian_threshold) negative_definite = false;
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
        std::vector<SymbolicExpression> constraints = symbolic_commands::parse_symbolic_expression_list(arguments[1], [](const std::string& s) { return s; });
        std::vector<std::string> variables = ctx.parse_symbolic_variable_arguments(arguments, 2, f.identifier_variables());
        SymbolicExpression lagrangian = f;
        std::vector<std::string> all_vars = variables;
        for (std::size_t i = 0; i < constraints.size(); ++i) {
            std::string lambda_var = "L" + std::to_string(i + 1);
            all_vars.push_back(lambda_var);
            lagrangian = (lagrangian - SymbolicExpression::variable(lambda_var) * constraints[i]).simplify();
        }
        std::vector<std::string> critical_args;
        critical_args.push_back(lagrangian.to_string());
        for (const auto& v : all_vars) critical_args.push_back(v);
        return handle_analysis_command(ctx, "critical", critical_args, output);
    }
    
    return false;
}


std::string AnalysisModule::execute_args_view(std::string_view command,
                                             const std::vector<std::string_view>& args,
                                             ::ServiceLocator& locator) {
    using namespace module_helpers;
    auto services = locator.resolve<CoreServices>();
    AnalysisContext ctx;
    ctx.resolve_symbolic = [&locator](const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) {
        locator.resolve<CoreServices>()->symbolic.resolve_symbolic(arg, req, var, expr);
    };
    ctx.parse_symbolic_variable_arguments = [services](const std::vector<std::string>& arguments, std::size_t start_index, const std::vector<std::string>& defaults) {
        return services->parse_symbolic_vars(arguments, start_index, defaults);
    };
    ctx.parse_decimal = [services](const std::string& expr) {
        return services->evaluation.parse_decimal(expr);
    };
    ctx.normalize_result = [services](Scalar value) {
        return services->evaluation.normalize_result(value);
    };
    ctx.build_analysis = [&locator, services](const std::string& expression) {
        SymbolicExpression expr; std::string var;
        locator.resolve<CoreServices>()->symbolic.resolve_symbolic(expression, false, &var, &expr);
        FunctionAnalysis analysis(var);
        analysis.define(expr.to_string());
        analysis.set_evaluator(services->evaluation.build_decimal_evaluator(expr.to_string()));
        return analysis;
    };

    // 转换参数为 string 向量用于 handle_analysis_command
    std::vector<std::string> string_args;
    for (const auto& arg : args) string_args.emplace_back(arg);

    std::string out;
    if (handle_analysis_command(ctx, std::string(command), string_args, &out)) return out;
    throw std::runtime_error("Unknown analysis command: " + std::string(command));
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

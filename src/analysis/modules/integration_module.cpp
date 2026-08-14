// ============================================================================
// 积分模块命令处理
// ============================================================================
//
// 本文件是积分模块的入口，负责：
// - 解析命令参数
// - 调用 integration/ 子目录中的核心算法
// - 格式化输出结果
//
// 核心计算已拆分到：
// - integration_commands.cpp: 线积分、曲面积分、重积分
// - vector_field_theorems.cpp: 向量场定理
// - multivariable_integrator.cpp: 多重积分器
// - multidim_integration.cpp: 多维积分

#include "analysis/modules/integration_module.h"
#include "execution/engine/script_context.h"
#include "symbolic/core/symbolic_expression.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "math/base/precision_constants.h"
#include "analysis/calculus/function_analysis.h"
#include "analysis/integration/integration_engine.h"
#include "analysis/integration/integration_commands.h"
#include "analysis/integration/multidim_integration.h"
#include "analysis/integration/multivariable_integrator.h"
#include "analysis/integration/vector_field_theorems.h"
#include "core/types/module_types.h"
#include "parser/grammars/unified_expression_parser.h"
#include "math/helpers/integer_helpers.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "types/scalar_type.h"

#include <stdexcept>
#include <vector>
#include <sstream>

namespace integration_ops {

using namespace utils;
using Scalar = mymath::Scalar;

namespace {

std::vector<int> parse_subdivisions(const IntegrationContext& ctx,
                                    const std::vector<std::string>& arguments,
                                    std::size_t offset,
                                    const std::vector<int>& defaults) {
    std::vector<int> subdivisions = defaults;
    if (arguments.size() <= offset) return subdivisions;
    std::size_t count = std::min(arguments.size() - offset, defaults.size());
    for (std::size_t i = 0; i < count; ++i) {
        const Scalar value = ctx.parse_decimal(arguments[offset + i]);
        if (!is_integer_double(value) || value <= 0.0L) throw std::runtime_error("integration subdivision counts must be positive integers");
        subdivisions[i] = static_cast<int>(round_to_long_long(value));
    }
    return subdivisions;
}

}  // namespace

bool is_integration_command(const std::string& command) {
    return command == "double_integral" || command == "double_integral_cyl" ||
           command == "triple_integral" || command == "triple_integral_sph" ||
           command == "line_integral" || command == "surface_integral" ||
           command == "integrate_region" || command == "greens_theorem" ||
           command == "divergence_theorem" || command == "stokes_theorem";
}

bool handle_integration_command(const IntegrationContext& ctx,
                                const std::string& command,
                                const std::vector<std::string>& arguments_orig,
                                std::string* output) {
    std::vector<std::string> arguments = arguments_orig;
    if (arguments.empty()) return false;
    if (command == "integral") {
        if (arguments.size() < 3) throw std::runtime_error("integral expects expr, x0, x1");
        Scalar x0, x1;
        if (arguments.size() >= 4 && is_identifier_text(utils::trim_copy(arguments[1]))) {
            // integral(f, x, a, b) - 显式指定变量
            std::string x_var = utils::trim_copy(arguments[1]);
            x0 = ctx.parse_decimal(arguments[2]);
            x1 = ctx.parse_decimal(arguments[3]);
            // 使用 FunctionAnalysis 进行自适应 Gauss-Kronrod 积分
            if (ctx.build_analysis) {
                FunctionAnalysis analysis(x_var);
                analysis.set_evaluator(ctx.build_scoped_evaluator(arguments[0]));
                *output = format_decimal(ctx.normalize_result(analysis.definite_integral(x0, x1)));
            } else {
                // 回退到 Simpson 方法
                auto f = ctx.build_scoped_evaluator(arguments[0]);
                MultivariableIntegrator integrator([&f, x_var](const std::vector<Scalar>& pt) { return f({{ x_var, pt[0] }}); });
                std::vector<MultivariableIntegrator::BoundFunc> bounds = { [x0, x1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {x0, x1}; } };
                *output = format_decimal(ctx.normalize_result(integrator.integrate(bounds, {1024})));
            }
            return true;
        } else {
            // integral(f, a, b) - 默认变量 x
            x0 = ctx.parse_decimal(arguments[1]);
            x1 = ctx.parse_decimal(arguments[2]);
            // 使用 FunctionAnalysis 进行自适应 Gauss-Kronrod 积分
            if (ctx.build_analysis) {
                FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
                *output = format_decimal(ctx.normalize_result(analysis.definite_integral(x0, x1)));
            } else {
                // 回退到 Simpson 方法
                auto f = ctx.build_scoped_evaluator(arguments[0]);
                MultivariableIntegrator integrator([&f](const std::vector<Scalar>& pt) { return f({{ "x", pt[0] }}); });
                std::vector<MultivariableIntegrator::BoundFunc> bounds = { [x0, x1](const std::vector<Scalar>&) -> std::pair<Scalar, Scalar> { return {x0, x1}; } };
                *output = format_decimal(ctx.normalize_result(integrator.integrate(bounds, {1024})));
            }
            return true;
        }
    }
    std::string coord_system = "rect";
    std::string method = "simpson";
    Scalar tol = 1e-6;
    bool has_options = true;
    while (arguments.size() >= 2 && has_options) {
        has_options = false;
        std::string last = utils::trim_copy(arguments.back());
        if (last.front() == '"' && last.back() == '"') last = last.substr(1, last.size() - 2);

        if (last == "polar" || last == "cyl" || last == "sph" || last == "rect") {
            coord_system = last;
            arguments.pop_back();
            has_options = true;
        } else if (last.find("method=") == 0) {
            method = last.substr(7);
            if (method.front() == '"' && method.back() == '"') method = method.substr(1, method.size() - 2);
            arguments.pop_back();
            has_options = true;
        } else if (last.find("tol=") == 0) {
            try { tol = ctx.parse_decimal(last.substr(4)); } catch (...) {}
            arguments.pop_back();
            has_options = true;
        } else if (last == "adaptive" || last == "monte_carlo" || last == "sparse_grid" || last == "simpson" || last == "tensor_product") {
            method = last;
            arguments.pop_back();
            has_options = true;
        }
    }

    if (command == "line_integral") {
        if (arguments.size() < 6) throw std::runtime_error("line_integral expects expr, t_var, t0, t1, x_expr, y_expr, [z_expr], [subdivides]");
        std::string t_var = utils::trim_copy(arguments[1]);
        Scalar t0 = ctx.parse_decimal(arguments[2]);
        Scalar t1 = ctx.parse_decimal(arguments[3]);
        std::string x_expr = arguments[4];
        std::string y_expr = arguments[5];
        std::string z_expr = "";
        int subdivides = 1024;
        if (arguments.size() >= 7) {
            try {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[6]));
            } catch (...) {
                z_expr = arguments[6];
                if (arguments.size() >= 8) {
                    subdivides = static_cast<int>(ctx.parse_decimal(arguments[7]));
                }
            }
        }
        *output = format_decimal(ctx.normalize_result(line_integral(ctx, arguments[0], t_var, t0, t1, x_expr, y_expr, z_expr, subdivides)));
        return true;
    }

    if (command == "surface_integral") {
        if (arguments.size() < 10) throw std::runtime_error("surface_integral expects expr, u_var, u0, u1, v_var, v0, v1, x_expr, y_expr, z_expr, [nu], [nv]");
        std::string u_var = utils::trim_copy(arguments[1]);
        Scalar u0 = ctx.parse_decimal(arguments[2]);
        Scalar u1 = ctx.parse_decimal(arguments[3]);
        std::string v_var = utils::trim_copy(arguments[4]);
        Scalar v0 = ctx.parse_decimal(arguments[5]);
        Scalar v1 = ctx.parse_decimal(arguments[6]);
        std::string x_expr = arguments[7];
        std::string y_expr = arguments[8];
        std::string z_expr = arguments[9];
        auto s = parse_subdivisions(ctx, arguments, 10, {256, 256});
        *output = format_decimal(ctx.normalize_result(surface_integral(ctx, arguments[0], u_var, u0, u1, v_var, v0, v1, x_expr, y_expr, z_expr, s[0], s[1])));
        return true;
    }
    if (command == "double_integral") {
        if (coord_system == "polar" || coord_system == "cyl") {
            if (arguments.size() >= 7 && is_identifier_text(utils::trim_copy(arguments[1]))) {
                auto subs = parse_subdivisions(ctx, arguments, 7, {64, 64});
                *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], subs[0], subs[1], method, tol)));
            } else {
                auto subs = parse_subdivisions(ctx, arguments, 5, {64, 64});
                *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], subs[1], subs[0], method, tol)));
            }
        } else {
            std::string xv = "x", yv = "y"; Scalar x0, x1; std::string y0e, y1e; size_t next;
            if (is_identifier_text(utils::trim_copy(arguments[1])) && arguments.size() >= 7) { xv = utils::trim_copy(arguments[1]); x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]); yv = utils::trim_copy(arguments[4]); y0e = arguments[5]; y1e = arguments[6]; next = 7; }
            else { x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]); y0e = arguments[3]; y1e = arguments[4]; next = 5; }
            auto subs = parse_subdivisions(ctx, arguments, next, {32, 32});
            *output = format_decimal(ctx.normalize_result(double_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, subs[0], subs[1], method, tol)));
        }
        return true;
    }
    if (command == "double_integral_cyl") {
        if (arguments.size() >= 7 && is_identifier_text(utils::trim_copy(arguments[1]))) {
            auto subs = parse_subdivisions(ctx, arguments, 7, {64, 64});
            *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], subs[0], subs[1], method, tol)));
        } else {
            auto subs = parse_subdivisions(ctx, arguments, 5, {64, 64});
            *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], subs[1], subs[0], method, tol)));
        }
        return true;
    }
    if (command == "triple_integral") {
        if (coord_system == "cyl") {
            if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[1]))) {
                auto s = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_cyl(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], arguments[7], arguments[8], arguments[9], s[0], s[1], s[2], method, tol)));
            } else {
                auto s = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_cyl(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], "z", arguments[5], arguments[6], s[1], s[0], s[2], method, tol)));
            }
        } else if (coord_system == "sph") {
            if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[1]))) {
                auto s = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), arguments[7], arguments[8], arguments[9], s[0], s[1], s[2], method, tol)));
            } else {
                auto s = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "phi", ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), "rho", arguments[1], arguments[2], s[1], s[2], s[0], method, tol)));
            }
        } else {
            std::string xv = "x", yv = "y", zv = "z"; Scalar x0, x1; std::string y0e, y1e, z0e, z1e; size_t next;
            if (is_identifier_text(utils::trim_copy(arguments[1])) && arguments.size() >= 10) { xv = utils::trim_copy(arguments[1]); x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]); yv = utils::trim_copy(arguments[4]); y0e = arguments[5]; y1e = arguments[6]; zv = utils::trim_copy(arguments[7]); z0e = arguments[8]; z1e = arguments[9]; next = 10; }
            else { x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]); y0e = arguments[3]; y1e = arguments[4]; z0e = arguments[5]; z1e = arguments[6]; next = 7; }
            auto s = parse_subdivisions(ctx, arguments, next, {16, 16, 16});
            *output = format_decimal(ctx.normalize_result(triple_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, zv, z0e, z1e, s[0], s[1], s[2], method, tol)));
        }
        return true;
    }
    if (command == "triple_integral_sph") {
        if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[1]))) {
            auto s = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
            *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), arguments[7], arguments[8], arguments[9], s[0], s[1], s[2], method, tol)));
        } else {
            auto s = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
            *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "phi", ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), "rho", arguments[1], arguments[2], s[1], s[2], s[0], method, tol)));
        }
        return true;
    }

    if (command == "integrate_region") {
        if (arguments.size() < 6) {
            throw std::runtime_error("integrate_region expects f, constraint, x0, x1, y0, y1, [z0, z1], [method], [samples]");
        }

        std::string f_expr = arguments[0];
        std::string constraint_expr = arguments[1];

        Scalar x0 = ctx.parse_decimal(arguments[2]);
        Scalar x1 = ctx.parse_decimal(arguments[3]);
        Scalar y0 = ctx.parse_decimal(arguments[4]);
        Scalar y1 = ctx.parse_decimal(arguments[5]);

        bool is_3d = arguments.size() >= 8 && is_identifier_text(utils::trim_copy(arguments[6]));
        Scalar z0 = 0.0L, z1 = 0.0L;
        std::string z_var = "z";

        if (is_3d) {
            z_var = utils::trim_copy(arguments[6]);
            z0 = ctx.parse_decimal(arguments[7]);
            z1 = ctx.parse_decimal(arguments[8]);
        }

        std::string region_method = "quasi_monte_carlo";
        int num_samples = 10000;

        size_t opt_offset = is_3d ? 9 : 6;
        if (arguments.size() > opt_offset) {
            std::string opt = utils::trim_copy(arguments[opt_offset]);
            if (opt == "monte_carlo" || opt == "quasi_monte_carlo" || opt == "adaptive") {
                region_method = opt;
                opt_offset++;
            }
        }
        if (arguments.size() > opt_offset) {
            num_samples = static_cast<int>(ctx.parse_decimal(arguments[opt_offset]));
        }

        auto f_eval = ctx.build_scoped_evaluator(f_expr);
        auto constraint_eval = ctx.build_scoped_evaluator(constraint_expr);

        if (is_3d) {
            multidim::MultidimFunction integrand = [f_eval, z_var](const std::vector<Scalar>& pt) {
                return f_eval({{"x", pt[0]}, {"y", pt[1]}, {z_var, pt[2]}});
            };
            multidim::RegionConstraint constraint = [constraint_eval, z_var](const std::vector<Scalar>& pt) {
                return constraint_eval({{"x", pt[0]}, {"y", pt[1]}, {z_var, pt[2]}});
            };

            std::vector<multidim::IntegrationBounds> bounds = {
                {x0, x1}, {y0, y1}, {z0, z1}
            };

            multidim::IntegrationOptions opts;
            opts.relative_tolerance = tol;
            opts.max_evaluations = num_samples;

            if (region_method == "monte_carlo") {
                auto result = multidim::integrate_monte_carlo(integrand, bounds, constraint, num_samples);
                *output = format_decimal(ctx.normalize_result(result.value)) + " (error: " +
                          format_decimal(result.error_estimate) + ")";
            } else {
                auto result = multidim::integrate_quasi_monte_carlo(integrand, bounds, constraint, num_samples);
                *output = format_decimal(ctx.normalize_result(result.value)) + " (error: " +
                          format_decimal(result.error_estimate) + ")";
            }
        } else {
            multidim::MultidimFunction integrand = [f_eval](const std::vector<Scalar>& pt) {
                return f_eval({{"x", pt[0]}, {"y", pt[1]}});
            };
            multidim::RegionConstraint constraint = [constraint_eval](const std::vector<Scalar>& pt) {
                return constraint_eval({{"x", pt[0]}, {"y", pt[1]}});
            };

            std::vector<multidim::IntegrationBounds> bounds = {{x0, x1}, {y0, y1}};

            if (region_method == "monte_carlo") {
                auto result = multidim::integrate_monte_carlo(integrand, bounds, constraint, num_samples);
                *output = format_decimal(ctx.normalize_result(result.value)) + " (error: " +
                          format_decimal(result.error_estimate) + ")";
            } else {
                auto result = multidim::integrate_quasi_monte_carlo(integrand, bounds, constraint, num_samples);
                *output = format_decimal(ctx.normalize_result(result.value)) + " (error: " +
                          format_decimal(result.error_estimate) + ")";
            }
        }
        return true;
    }

    if (command == "greens_theorem") {
        if (arguments.size() < 7) {
            throw std::runtime_error("greens_theorem expects P, Q, curve_x, curve_y, t, t0, t1, [subdivides]");
        }

        std::string P = arguments[0];
        std::string Q = arguments[1];

        if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[5]))) {
            std::string x_var = utils::trim_copy(arguments[2]);
            Scalar x0 = ctx.parse_decimal(arguments[3]);
            Scalar x1 = ctx.parse_decimal(arguments[4]);
            std::string y_var = utils::trim_copy(arguments[5]);
            std::string y0_expr = arguments[6];
            std::string y1_expr = arguments[7];

            int subdivides = 32;
            if (arguments.size() > 8) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[8]));
            }

            auto result = greens_theorem_area(ctx, P, Q, x_var, x0, x1, y_var, y0_expr, y1_expr, subdivides);
            *output = format_decimal(ctx.normalize_result(result.value)) + " (method: " + result.method_used + ")";
        } else {
            std::string curve_x = arguments[2];
            std::string curve_y = arguments[3];
            std::string t_var = utils::trim_copy(arguments[4]);
            Scalar t0 = ctx.parse_decimal(arguments[5]);
            Scalar t1 = ctx.parse_decimal(arguments[6]);

            int subdivides = 64;
            if (arguments.size() > 7) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[7]));
            }

            auto result = greens_theorem(ctx, P, Q, curve_x, curve_y, t_var, t0, t1, subdivides);
            *output = format_decimal(ctx.normalize_result(result.value)) + " (method: " + result.method_used + ")";
        }
        return true;
    }

    if (command == "divergence_theorem") {
        if (arguments.size() < 9) {
            throw std::runtime_error("divergence_theorem expects Fx, Fy, Fz, surface/region params");
        }

        std::string Fx = arguments[0];
        std::string Fy = arguments[1];
        std::string Fz = arguments[2];

        bool is_surface = arguments.size() >= 12 && is_identifier_text(utils::trim_copy(arguments[6]));

        if (is_surface) {
            std::string surface_x = arguments[3];
            std::string surface_y = arguments[4];
            std::string surface_z = arguments[5];
            std::string u_var = utils::trim_copy(arguments[6]);
            Scalar u0 = ctx.parse_decimal(arguments[7]);
            Scalar u1 = ctx.parse_decimal(arguments[8]);
            std::string v_var = utils::trim_copy(arguments[9]);
            Scalar v0 = ctx.parse_decimal(arguments[10]);
            Scalar v1 = ctx.parse_decimal(arguments[11]);

            std::string orientation = "outward";
            int subdivides = 32;
            size_t opt_offset = 12;

            if (arguments.size() > opt_offset) {
                std::string opt = utils::trim_copy(arguments[opt_offset]);
                if (opt == "outward" || opt == "inward") {
                    orientation = opt;
                    opt_offset++;
                }
            }
            if (arguments.size() > opt_offset) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[opt_offset]));
            }

            auto result = divergence_theorem(ctx, Fx, Fy, Fz, surface_x, surface_y, surface_z,
                                             u_var, u0, u1, v_var, v0, v1, orientation, subdivides);
            *output = format_decimal(ctx.normalize_result(result.value)) + " (method: " + result.method_used + ")";
        } else {
            std::string x_var = utils::trim_copy(arguments[3]);
            Scalar x0 = ctx.parse_decimal(arguments[4]);
            Scalar x1 = ctx.parse_decimal(arguments[5]);
            std::string y_var = utils::trim_copy(arguments[6]);
            std::string y0_expr = arguments[7];
            std::string y1_expr = arguments[8];
            std::string z_var = utils::trim_copy(arguments[9]);
            std::string z0_expr = arguments[10];
            std::string z1_expr = arguments[11];

            int subdivides = 16;
            if (arguments.size() > 12) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[12]));
            }

            auto result = divergence_theorem_volume(ctx, Fx, Fy, Fz,
                                                    x_var, x0, x1, y_var, y0_expr, y1_expr,
                                                    z_var, z0_expr, z1_expr, subdivides);
            *output = format_decimal(ctx.normalize_result(result.value)) + " (method: " + result.method_used + ")";
        }
        return true;
    }

    if (command == "stokes_theorem") {
        if (arguments.size() < 9) {
            throw std::runtime_error("stokes_theorem expects Fx, Fy, Fz, curve params");
        }

        std::string Fx = arguments[0];
        std::string Fy = arguments[1];
        std::string Fz = arguments[2];
        std::string curve_x = arguments[3];
        std::string curve_y = arguments[4];
        std::string curve_z = arguments[5];
        std::string t_var = utils::trim_copy(arguments[6]);
        Scalar t0 = ctx.parse_decimal(arguments[7]);
        Scalar t1 = ctx.parse_decimal(arguments[8]);

        int subdivides = 64;
        if (arguments.size() > 9 && !is_identifier_text(utils::trim_copy(arguments[9]))) {
            subdivides = static_cast<int>(ctx.parse_decimal(arguments[9]));
            *output = format_decimal(ctx.normalize_result(
                stokes_theorem_line(ctx, Fx, Fy, Fz, curve_x, curve_y, curve_z, t_var, t0, t1, subdivides).value));
            return true;
        }

        if (arguments.size() >= 18) {
            std::string surface_x = arguments[9];
            std::string surface_y = arguments[10];
            std::string surface_z = arguments[11];
            std::string u_var = utils::trim_copy(arguments[12]);
            Scalar u0 = ctx.parse_decimal(arguments[13]);
            Scalar u1 = ctx.parse_decimal(arguments[14]);
            std::string v_var = utils::trim_copy(arguments[15]);
            Scalar v0 = ctx.parse_decimal(arguments[16]);
            Scalar v1 = ctx.parse_decimal(arguments[17]);

            std::string orientation = "outward";
            subdivides = 32;
            size_t opt_offset = 18;

            if (arguments.size() > opt_offset) {
                std::string opt = utils::trim_copy(arguments[opt_offset]);
                if (opt == "outward" || opt == "inward") {
                    orientation = opt;
                    opt_offset++;
                }
            }
            if (arguments.size() > opt_offset) {
                subdivides = static_cast<int>(ctx.parse_decimal(arguments[opt_offset]));
            }

            auto result = stokes_theorem(ctx, Fx, Fy, Fz, curve_x, curve_y, curve_z, t_var, t0, t1,
                                         surface_x, surface_y, surface_z, u_var, u0, u1, v_var, v0, v1,
                                         orientation, subdivides);
            std::ostringstream oss;
            oss << format_decimal(ctx.normalize_result(result.value)) << " (method: " << result.method_used;
            if (result.verified) {
                oss << ", verified, diff: " << format_decimal(result.verification_diff);
            }
            oss << ")";
            *output = oss.str();
            return true;
        }

        *output = format_decimal(ctx.normalize_result(
            stokes_theorem_line(ctx, Fx, Fy, Fz, curve_x, curve_y, curve_z, t_var, t0, t1, subdivides).value));
        return true;
    }

    return false;
}

bool handle_integration_command(const IntegrationContext& ctx, const std::string& command, const std::string& inside, std::string* output) {
    return handle_integration_command(ctx, command, split_top_level_arguments(inside), output);
}

std::string IntegrationModule::execute_args(const std::string& command,
                                            const std::vector<std::string>& args,
                                            ServiceLocator& locator) {
    auto services = locator.resolve<CoreServices>();
    IntegrationContext ctx;
    ctx.parse_decimal = [services](const std::string& expr) { return services->evaluation.parse_decimal(expr); };
    ctx.build_scoped_evaluator = [services](const std::string& expression) { return services->evaluation.build_decimal_evaluator(expression); };
    ctx.normalize_result = [services](Scalar value) { return services->evaluation.normalize_result(value); };
    ctx.build_analysis = [&locator, services](const std::string& expression) {
        SymbolicExpression expr; std::string var;
        locator.resolve<CoreServices>()->symbolic.resolve_symbolic(expression, false, &var, &expr);
        FunctionAnalysis analysis(var);
        analysis.define(expr.to_string());
        analysis.set_evaluator(services->evaluation.build_decimal_evaluator(expr.to_string()));
        return analysis;
    };

    std::string output;
    if (handle_integration_command(ctx, command, args, &output)) {
        return output;
    }
    throw std::runtime_error("Unknown integration command: " + command);
}

std::string IntegrationModule::get_help_snippet(const std::string& topic) const {
    if (topic == "analysis") {
        return "Numerical Integration:\n"
               "  double_integral(f, x0, x1, y0, y1, [\"polar\"]) 2D numerical integration\n"
               "  triple_integral(f, x0, x1, y0, y1, z0, z1, [\"cyl\"|\"sph\"]) 3D numerical integration\n"
               "  line_integral(f, t, t0, t1, x(t), y(t), [z(t)]) Line integral\n"
               "  surface_integral(f, u, u0, u1, v, v0, v1, x(u,v), y(u,v), z(u,v)) Surface integral\n"
               "  integrate_region(f, constraint, x0, x1, y0, y1, [z0, z1]) Implicit region integral\n"
               "Vector Field Theorems:\n"
               "  greens_theorem(P, Q, curve_x, curve_y, t, t0, t1) Green's theorem\n"
               "  divergence_theorem(Fx, Fy, Fz, surface...) Divergence theorem\n"
               "  stokes_theorem(Fx, Fy, Fz, curve..., surface...) Stokes' theorem";
    }
    return "";
}

std::vector<std::string> IntegrationModule::get_commands() const {
    return {"double_integral", "double_integral_cyl", "triple_integral", "triple_integral_sph",
            "integrate_region", "greens_theorem", "divergence_theorem", "stokes_theorem"};
}

}  // namespace integration_ops

REGISTER_CALCULATOR_MODULE(integration_ops::IntegrationModule)

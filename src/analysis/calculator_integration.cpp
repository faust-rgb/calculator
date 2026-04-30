// ============================================================================
// 多重积分命令实现
// ============================================================================
//
// 本文件实现了多重积分命令的数值计算，包括：
// - double_integral: 二重积分（直角坐标）
// - double_integral_polar / double_integral_cyl: 极坐标二重积分
// - triple_integral: 三重积分（直角坐标）
// - triple_integral_cyl: 柱坐标三重积分
// - triple_integral_sph: 球坐标三重积分

#include "calculator_integration.h"
#include "calculator_internal_types.h"

#include "mymath.h"
#include "multidim_integration.h"

#include <stdexcept>
#include <vector>

namespace integration_ops {

namespace {

/**
 * @brief 解析积分分割数参数
 *
 * @param ctx 积分上下文
 * @param arguments 参数列表
 * @param offset 起始索引
 * @param defaults 默认分割数
 * @return 实际分割数
 */
std::vector<int> parse_subdivisions(const IntegrationContext& ctx,
                                    const std::vector<std::string>& arguments,
                                    std::size_t offset,
                                    const std::vector<int>& defaults) {
    std::vector<int> subdivisions = defaults;
    if (arguments.size() <= offset) {
        return subdivisions;
    }

    std::size_t count = arguments.size() - offset;
    if (count > defaults.size()) {
        count = defaults.size();
    }

    for (std::size_t i = 0; i < count; ++i) {
        const double value = ctx.parse_decimal(arguments[offset + i]);
        if (!is_integer_double(value) || value <= 0.0) {
            throw std::runtime_error(
                "integration subdivision counts must be positive integers");
        }
        subdivisions[i] = static_cast<int>(round_to_long_long(value));
    }
    return subdivisions;
}

/**
 * @brief 创建积分边界函数
 *
 * 边界可以是常数或表达式。如果是表达式，则返回一个依赖于外层变量的函数。
 */
std::function<double(const std::vector<double>&)> make_scalar_bound_func(
    const IntegrationContext& ctx,
    const std::string& bound_expr,
    const std::vector<std::string>& var_names,
    std::size_t current_dim) {

    // 首先尝试解析为常数
    try {
        double val = ctx.parse_decimal(bound_expr);
        return [val](const std::vector<double>&) { return val; };
    } catch (...) {
        // 不是常数，回退到表达式求值
    }

    // 创建表达式求值器
    auto evaluator = ctx.build_scoped_evaluator(bound_expr);
    return [evaluator, var_names, current_dim](const std::vector<double>& point) {
        std::vector<std::pair<std::string, double>> scope;
        scope.reserve(current_dim);
        for (std::size_t i = 0; i < current_dim; ++i) {
            scope.push_back({var_names[i], point[i]});
        }
        return evaluator(scope);
    };
}

}  // namespace

/**
 * @brief 计算二重积分（直角坐标）
 *
 * 计算 ∫∫ f(x,y) dx dy，支持变边界。
 */
double double_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& x_var, double x0, double x1,
    const std::string& y_var, const std::string& y0_expr, const std::string& y1_expr,
    int nx, int ny) {

    // 优化点：如果边界是常数，优先使用自适应 Gauss-Kronrod 积分以提高精度
    bool y_bounds_are_constant = false;
    double y0_val = 0, y1_val = 0;
    try {
        y0_val = ctx.parse_decimal(y0_expr);
        y1_val = ctx.parse_decimal(y1_expr);
        y_bounds_are_constant = true;
    } catch (...) {
        y_bounds_are_constant = false;
    }

    if (y_bounds_are_constant) {
        const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
        auto result = multidim::integrate_2d_adaptive(
            [evaluate_expression, x_var, y_var](double x, double y) {
                return evaluate_expression({{x_var, x}, {y_var, y}});
            },
            multidim::IntegrationBounds(x0, x1),
            multidim::IntegrationBounds(y0_val, y1_val));
        if (result.converged) {
            return result.value;
        }
    }

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression, x_var, y_var](const std::vector<double>& point) {
            return evaluate_expression({{x_var, point[0]}, {y_var, point[1]}});
        });

    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([x0, x1](const std::vector<double>&) -> std::pair<double, double> { 
        return std::make_pair(x0, x1); 
    });
    
    auto y0_f = make_scalar_bound_func(ctx, y0_expr, {x_var}, 1);
    auto y1_f = make_scalar_bound_func(ctx, y1_expr, {x_var}, 1);
    bounds.push_back([y0_f, y1_f](const std::vector<double>& pt) -> std::pair<double, double> {
        return std::make_pair(y0_f(pt), y1_f(pt));
    });

    return integrator.integrate(bounds, {nx, ny});
}

/**
 * @brief 计算极坐标二重积分
 *
 * 计算 ∫∫ f(r,theta) * r dr dtheta。
 * 支持变边界：r = r(theta)。
 */
double double_integral_polar(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& theta_var, double theta0, double theta1,
    const std::string& r_var, const std::string& r0_expr, const std::string& r1_expr,
    int ntheta, int nr) {

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression, theta_var, r_var](const std::vector<double>& point) {
            const double theta = point[0];
            const double r = point[1];
            const double x = r * mymath::cos(theta);
            const double y = r * mymath::sin(theta);
            return evaluate_expression(
                       {{r_var, r}, {theta_var, theta}, {"x", x}, {"y", y}}) * r;
        });

    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    // theta 为外层变量
    bounds.push_back([theta0, theta1](const std::vector<double>&) -> std::pair<double, double> { 
        return std::make_pair(theta0, theta1); 
    });
    
    // r 为内层变量，可依赖于 theta
    auto r0_f = make_scalar_bound_func(ctx, r0_expr, {theta_var}, 1);
    auto r1_f = make_scalar_bound_func(ctx, r1_expr, {theta_var}, 1);
    bounds.push_back([r0_f, r1_f](const std::vector<double>& pt) -> std::pair<double, double> {
        return std::make_pair(r0_f(pt), r1_f(pt));
    });

    return integrator.integrate(bounds, {ntheta, nr});
}

/**
 * @brief 计算三重积分（直角坐标）
 *
 * 计算 ∫∫∫ f(x,y,z) dx dy dz，支持变边界。
 */
double triple_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& x_var, double x0, double x1,
    const std::string& y_var, const std::string& y0_expr, const std::string& y1_expr,
    const std::string& z_var, const std::string& z0_expr, const std::string& z1_expr,
    int nx, int ny, int nz) {

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression, x_var, y_var, z_var](const std::vector<double>& point) {
            return evaluate_expression({{x_var, point[0]}, {y_var, point[1]}, {z_var, point[2]}});
        });

    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([x0, x1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(x0, x1); });
    
    auto y0_f = make_scalar_bound_func(ctx, y0_expr, {x_var}, 1);
    auto y1_f = make_scalar_bound_func(ctx, y1_expr, {x_var}, 1);
    bounds.push_back([y0_f, y1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(y0_f(pt), y1_f(pt)); });

    auto z0_f = make_scalar_bound_func(ctx, z0_expr, {x_var, y_var}, 2);
    auto z1_f = make_scalar_bound_func(ctx, z1_expr, {x_var, y_var}, 2);
    bounds.push_back([z0_f, z1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(z0_f(pt), z1_f(pt)); });

    return integrator.integrate(bounds, {nx, ny, nz});
}

/**
 * @brief 计算柱坐标三重积分
 *
 * 计算 ∫∫∫ f(r,theta,z) * r dr dtheta dz。
 * 支持变边界：r(theta), z(r, theta)。
 */
double triple_integral_cyl(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& theta_var, double theta0, double theta1,
    const std::string& r_var, const std::string& r0_expr, const std::string& r1_expr,
    const std::string& z_var, const std::string& z0_expr, const std::string& z1_expr,
    int ntheta, int nr, int nz) {

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression, theta_var, r_var, z_var](const std::vector<double>& point) {
            const double theta = point[0];
            const double r = point[1];
            const double z = point[2];
            const double x = r * mymath::cos(theta);
            const double y = r * mymath::sin(theta);
            return evaluate_expression(
                       {{r_var, r}, {theta_var, theta}, {z_var, z}, {"x", x}, {"y", y}}) * r;
        });

    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    // theta 为最外层
    bounds.push_back([theta0, theta1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(theta0, theta1); });
    
    // r(theta)
    auto r0_f = make_scalar_bound_func(ctx, r0_expr, {theta_var}, 1);
    auto r1_f = make_scalar_bound_func(ctx, r1_expr, {theta_var}, 1);
    bounds.push_back([r0_f, r1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(r0_f(pt), r1_f(pt)); });

    // z(theta, r)
    auto z0_f = make_scalar_bound_func(ctx, z0_expr, {theta_var, r_var}, 2);
    auto z1_f = make_scalar_bound_func(ctx, z1_expr, {theta_var, r_var}, 2);
    bounds.push_back([z0_f, z1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(z0_f(pt), z1_f(pt)); });

    return integrator.integrate(bounds, {ntheta, nr, nz});
}

/**
 * @brief 计算球坐标三重积分
 *
 * 计算 ∫∫∫ f(rho,theta,phi) * rho²sin(phi) drho dtheta dphi。
 * 支持变边界。
 */
double triple_integral_sph(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& theta_var, double theta0, double theta1,
    const std::string& phi_var, double phi0, double phi1,
    const std::string& rho_var, const std::string& rho0_expr, const std::string& rho1_expr,
    int ntheta, int nphi, int nrho) {

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression, theta_var, phi_var, rho_var](const std::vector<double>& point) {
            const double theta = point[0];
            const double phi = point[1];
            const double rho = point[2];
            const double sin_phi = mymath::sin(phi);
            const double cos_phi = mymath::cos(phi);
            const double cos_theta = mymath::cos(theta);
            const double sin_theta = mymath::sin(theta);
            const double x = rho * sin_phi * cos_theta;
            const double y = rho * sin_phi * sin_theta;
            const double z = rho * cos_phi;
            return evaluate_expression(
                       {{rho_var, rho}, {theta_var, theta}, {phi_var, phi},
                        {"x", x}, {"y", y}, {"z", z}}) *
                   rho * rho * sin_phi;
        });

    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    // theta, phi 为外层
    bounds.push_back([theta0, theta1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(theta0, theta1); });
    bounds.push_back([phi0, phi1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(phi0, phi1); });
    
    // rho(theta, phi)
    auto rho0_f = make_scalar_bound_func(ctx, rho0_expr, {theta_var, phi_var}, 2);
    auto rho1_f = make_scalar_bound_func(ctx, rho1_expr, {theta_var, phi_var}, 2);
    bounds.push_back([rho0_f, rho1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(rho0_f(pt), rho1_f(pt)); });

    return integrator.integrate(bounds, {ntheta, nphi, nrho});
}

/**
 * @brief 检查是否为多重积分命令
 */
bool is_integration_command(const std::string& command) {
    return command == "double_integral" ||
           command == "triple_integral";
}

/**
 * @brief 处理多重积分命令
 *
 * 解析参数并调用相应的积分函数。
 */
bool handle_integration_command(const IntegrationContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output) {
    std::vector<std::string> arguments = split_top_level_arguments(inside);
    if (arguments.empty()) return false;

    // 检查最后一个参数是否为坐标系指定
    std::string coord_system = "rect";
    bool has_coord_arg = false;
    if (arguments.size() >= 2) {
        std::string last_arg = trim_copy(arguments.back());
        if (last_arg == "\"polar\"" || last_arg == "polar" ||
            last_arg == "\"cyl\"" || last_arg == "cyl" ||
            last_arg == "\"sph\"" || last_arg == "sph" ||
            last_arg == "\"rect\"" || last_arg == "rect") {
            coord_system = last_arg;
            if (coord_system.front() == '"') {
                coord_system = coord_system.substr(1, coord_system.size() - 2);
            }
            arguments.pop_back();
            has_coord_arg = true;
        }
    }

    if (command == "double_integral") {
        if (coord_system == "polar" || coord_system == "cyl") {
            // 极坐标处理逻辑
            bool is_new_style = false;
            if (arguments.size() >= 7 && is_identifier_text(arguments[1]) && is_identifier_text(arguments[4])) {
                is_new_style = true;
            }

            if (is_new_style) {
                const std::string theta_var = arguments[1];
                const double theta0 = ctx.parse_decimal(arguments[2]);
                const double theta1 = ctx.parse_decimal(arguments[3]);
                const std::string r_var = arguments[4];
                const std::string r0_e = arguments[5];
                const std::string r1_e = arguments[6];
                const std::vector<int> subs = parse_subdivisions(ctx, arguments, 7, {64, 64});
                
                double result = double_integral_polar(ctx, arguments[0], theta_var, theta0, theta1, r_var, r0_e, r1_e, subs[0], subs[1]);
                *output = format_decimal(ctx.normalize_result(result));
            } else {
                if (arguments.size() != 5 && arguments.size() != 7) {
                    throw std::runtime_error("double_integral (polar) expects expr, r0, r1, theta0, theta1, and optional nr, ntheta");
                }
                const double r0 = ctx.parse_decimal(arguments[1]);
                const double r1 = ctx.parse_decimal(arguments[2]);
                const double t0 = ctx.parse_decimal(arguments[3]);
                const double t1 = ctx.parse_decimal(arguments[4]);
                const std::vector<int> subs = parse_subdivisions(ctx, arguments, 5, {64, 64});
                double result = double_integral_polar(ctx, arguments[0], "theta", t0, t1, "r", std::to_string(r0), std::to_string(r1), subs[1], subs[0]);
                *output = format_decimal(ctx.normalize_result(result));
            }
            return true;
        } else {
            // 直角坐标处理逻辑
            if (arguments.size() < 5) {
                throw std::runtime_error("double_integral expects expr, [x,] x0, x1, [y,] y0, y1, [nx, ny]");
            }
            
            std::string x_var = "x", y_var = "y";
            double x0, x1;
            std::string y0_expr, y1_expr;
            std::size_t next_idx = 0;
            
            bool custom_vars = false;
            try {
                (void)ctx.parse_decimal(arguments[1]);
            } catch (...) {
                if (is_identifier_text(arguments[1])) custom_vars = true;
            }

            if (custom_vars) {
                if (arguments.size() != 7 && arguments.size() != 9) {
                    throw std::runtime_error("double_integral with variable names expects expr, x, x0, x1, y, y0, y1, and optional nx, ny");
                }
                x_var = arguments[1]; x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]);
                y_var = arguments[4]; y0_expr = arguments[5]; y1_expr = arguments[6];
                next_idx = 7;
            } else {
                if (arguments.size() != 5 && arguments.size() != 7) {
                    throw std::runtime_error("double_integral expects expr, x0, x1, y0, y1, and optional nx, ny");
                }
                x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]);
                y0_expr = arguments[3]; y1_expr = arguments[4];
                next_idx = 5;
            }

            const std::vector<int> subdivisions = parse_subdivisions(ctx, arguments, next_idx, {32, 32});
            double result = double_integral(ctx, arguments[0], x_var, x0, x1, y_var, y0_expr, y1_expr, subdivisions[0], subdivisions[1]);
            *output = format_decimal(ctx.normalize_result(result));
            return true;
        }
    }

    if (command == "triple_integral") {
        if (coord_system == "cyl") {
            // 柱坐标处理逻辑
            bool is_new = (arguments.size() >= 10 && is_identifier_text(arguments[1]));
            if (is_new) {
                const std::vector<int> subs = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
                double result = triple_integral_cyl(ctx, arguments[0], 
                    arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]),
                    arguments[4], arguments[5], arguments[6],
                    arguments[7], arguments[8], arguments[9],
                    subs[0], subs[1], subs[2]);
                *output = format_decimal(ctx.normalize_result(result));
            } else {
                if (arguments.size() != 7 && arguments.size() != 10) {
                    throw std::runtime_error("triple_integral (cyl) expects expr, r0, r1, theta0, theta1, z0, z1, and optional nr, ntheta, nz");
                }
                const std::vector<int> subs = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
                double result = triple_integral_cyl(ctx, arguments[0],
                    "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]),
                    "r", arguments[1], arguments[2], "z", arguments[5], arguments[6],
                    subs[1], subs[0], subs[2]);
                *output = format_decimal(ctx.normalize_result(result));
            }
            return true;
        } else if (coord_system == "sph") {
            // 球坐标处理逻辑
            bool is_new = (arguments.size() >= 10 && is_identifier_text(arguments[1]));
            if (is_new) {
                const std::vector<int> subs = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
                double result = triple_integral_sph(ctx, arguments[0],
                    arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]),
                    arguments[4], ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]),
                    arguments[7], arguments[8], arguments[9],
                    subs[0], subs[1], subs[2]);
                *output = format_decimal(ctx.normalize_result(result));
            } else {
                if (arguments.size() != 7 && arguments.size() != 10) {
                    throw std::runtime_error("triple_integral (sph) expects expr, rho0, rho1, theta0, theta1, phi0, phi1, and optional nrho, ntheta, nphi");
                }
                const std::vector<int> subs = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
                double result = triple_integral_sph(ctx, arguments[0],
                    "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]),
                    "phi", ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]),
                    "rho", arguments[1], arguments[2],
                    subs[1], subs[2], subs[0]);
                *output = format_decimal(ctx.normalize_result(result));
            }
            return true;
        } else {
            // 直角坐标处理逻辑
            if (arguments.size() < 7) {
                throw std::runtime_error("triple_integral expects expr, [x, x0, x1, y, y0, y1, z, z0, z1] or [x0, x1, y0, y1, z0, z1]");
            }

            std::string x_var = "x", y_var = "y", z_var = "z";
            double x0, x1;
            std::string y0_e, y1_e, z0_e, z1_e;
            std::size_t next_idx = 0;

            bool custom_vars = false;
            try {
                (void)ctx.parse_decimal(arguments[1]);
            } catch (...) {
                if (is_identifier_text(arguments[1])) custom_vars = true;
            }

            if (custom_vars) {
                if (arguments.size() != 10 && arguments.size() != 13) {
                    throw std::runtime_error("triple_integral with variable names expects expr, x, x0, x1, y, y0, y1, z, z0, z1, and optional nx, ny, nz");
                }
                x_var = arguments[1]; x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]);
                y_var = arguments[4]; y0_e = arguments[5]; y1_e = arguments[6];
                z_var = arguments[7]; z0_e = arguments[8]; z1_e = arguments[9];
                next_idx = 10;
            } else {
                if (arguments.size() != 7 && arguments.size() != 10) {
                    throw std::runtime_error("triple_integral expects expr, x0, x1, y0, y1, z0, z1, and optional nx, ny, nz");
                }
                x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]);
                y0_e = arguments[3]; y1_e = arguments[4]; z0_e = arguments[5]; z1_e = arguments[6];
                next_idx = 7;
            }

            const std::vector<int> subdivisions = parse_subdivisions(ctx, arguments, next_idx, {16, 16, 16});
            double result = triple_integral(ctx, arguments[0], x_var, x0, x1, y_var, y0_e, y1_e, z_var, z0_e, z1_e, subdivisions[0], subdivisions[1], subdivisions[2]);
            *output = format_decimal(ctx.normalize_result(result));
            return true;
        }
    }

    return false;
}

}  // namespace integration_ops

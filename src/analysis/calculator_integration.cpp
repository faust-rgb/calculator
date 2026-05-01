// ============================================================================
// 多重积分命令实现
// ============================================================================

#include "calculator_integration.h"
#include "calculator_internal_types.h"
#include "mymath.h"
#include "multidim_integration.h"
#include "multivariable_integrator.h"
#include "utils.h"

#include <stdexcept>
#include <vector>
#include <sstream>

namespace integration_ops {

using namespace utils;

namespace {

std::vector<int> parse_subdivisions(const IntegrationContext& ctx,
                                    const std::vector<std::string>& arguments,
                                    std::size_t offset,
                                    const std::vector<int>& defaults) {
    std::vector<int> subdivisions = defaults;
    if (arguments.size() <= offset) return subdivisions;
    std::size_t count = std::min(arguments.size() - offset, defaults.size());
    for (std::size_t i = 0; i < count; ++i) {
        const double value = ctx.parse_decimal(arguments[offset + i]);
        if (!is_integer_double(value) || value <= 0.0) throw std::runtime_error("integration subdivision counts must be positive integers");
        subdivisions[i] = static_cast<int>(round_to_long_long(value));
    }
    return subdivisions;
}

std::function<double(const std::vector<double>&)> make_scalar_bound_func(
    const IntegrationContext& ctx,
    const std::string& bound_expr,
    const std::vector<std::string>& var_names,
    std::size_t current_dim) {
    try {
        double val = ctx.parse_decimal(bound_expr);
        return [val](const std::vector<double>&) { return val; };
    } catch (...) {}
    auto evaluator = ctx.build_scoped_evaluator(bound_expr);
    return [evaluator, var_names, current_dim](const std::vector<double>& point) {
        std::vector<std::pair<std::string, double>> scope;
        for (std::size_t i = 0; i < current_dim; ++i) scope.push_back({var_names[i], point[i]});
        return evaluator(scope);
    };
}

}  // namespace

double double_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& x_var, double x0, double x1,
    const std::string& y_var, const std::string& y0_expr, const std::string& y1_expr,
    int nx, int ny) {
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression, x_var, y_var](const std::vector<double>& point) {
            return evaluate_expression({{x_var, point[0]}, {y_var, point[1]}});
        });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([x0, x1](const std::vector<double>&) -> std::pair<double, double> { return {x0, x1}; });
    auto y0_f = make_scalar_bound_func(ctx, y0_expr, {x_var}, 1);
    auto y1_f = make_scalar_bound_func(ctx, y1_expr, {x_var}, 1);
    bounds.push_back([y0_f, y1_f](const std::vector<double>& pt) -> std::pair<double, double> { return {y0_f(pt), y1_f(pt)}; });
    return integrator.integrate(bounds, {nx, ny});
}

double double_integral_polar(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& theta_var, double theta0, double theta1,
    const std::string& r_var, const std::string& r0_expr, const std::string& r1_expr,
    int ntheta, int nr) {
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression, theta_var, r_var](const std::vector<double>& point) {
            const double theta = point[0], r = point[1];
            return evaluate_expression({{r_var, r}, {theta_var, theta}, {"x", r * mymath::cos(theta)}, {"y", r * mymath::sin(theta)}}) * r;
        });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([theta0, theta1](const std::vector<double>&) -> std::pair<double, double> { return {theta0, theta1}; });
    auto r0_f = make_scalar_bound_func(ctx, r0_expr, {theta_var}, 1);
    auto r1_f = make_scalar_bound_func(ctx, r1_expr, {theta_var}, 1);
    bounds.push_back([r0_f, r1_f](const std::vector<double>& pt) -> std::pair<double, double> { return {r0_f(pt), r1_f(pt)}; });
    return integrator.integrate(bounds, {ntheta, nr});
}

double triple_integral(const IntegrationContext& ctx, const std::string& expr, const std::string& x_v, double x0, double x1, const std::string& y_v, const std::string& y0_e, const std::string& y1_e, const std::string& z_v, const std::string& z0_e, const std::string& z1_e, int nx, int ny, int nz) {
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator([evaluate_expression, x_v, y_v, z_v](const std::vector<double>& pt) { return evaluate_expression({{x_v, pt[0]}, {y_v, pt[1]}, {z_v, pt[2]}}); });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([x0, x1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(x0, x1); });
    auto y0_f = make_scalar_bound_func(ctx, y0_e, {x_v}, 1);
    auto y1_f = make_scalar_bound_func(ctx, y1_e, {x_v}, 1);
    bounds.push_back([y0_f, y1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(y0_f(pt), y1_f(pt)); });
    auto z0_f = make_scalar_bound_func(ctx, z0_e, {x_v, y_v}, 2);
    auto z1_f = make_scalar_bound_func(ctx, z1_e, {x_v, y_v}, 2);
    bounds.push_back([z0_f, z1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(z0_f(pt), z1_f(pt)); });
    return integrator.integrate(bounds, {nx, ny, nz});
}

double triple_integral_cyl(const IntegrationContext& ctx, const std::string& expr, const std::string& t_v, double t0, double t1, const std::string& r_v, const std::string& r0_e, const std::string& r1_e, const std::string& z_v, const std::string& z0_e, const std::string& z1_e, int nt, int nr, int nz) {
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator([evaluate_expression, t_v, r_v, z_v](const std::vector<double>& pt) { 
        double t = pt[0], r = pt[1], z = pt[2];
        return evaluate_expression({{r_v, r}, {t_v, t}, {z_v, z}, {"x", r * mymath::cos(t)}, {"y", r * mymath::sin(t)}}) * r;
    });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([t0, t1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(t0, t1); });
    auto r0_f = make_scalar_bound_func(ctx, r0_e, {t_v}, 1);
    auto r1_f = make_scalar_bound_func(ctx, r1_e, {t_v}, 1);
    bounds.push_back([r0_f, r1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(r0_f(pt), r1_f(pt)); });
    auto z0_f = make_scalar_bound_func(ctx, z0_e, {t_v, r_v}, 2);
    auto z1_f = make_scalar_bound_func(ctx, z1_e, {t_v, r_v}, 2);
    bounds.push_back([z0_f, z1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(z0_f(pt), z1_f(pt)); });
    return integrator.integrate(bounds, {nt, nr, nz});
}

double triple_integral_sph(const IntegrationContext& ctx, const std::string& expr, const std::string& t_v, double t0, double t1, const std::string& p_v, double p0, double p1, const std::string& r_v, const std::string& r0_e, const std::string& r1_e, int nt, int np, int nr) {
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator([evaluate_expression, t_v, p_v, r_v](const std::vector<double>& pt) {
        double t = pt[0], p = pt[1], r = pt[2], sp = mymath::sin(p);
        return evaluate_expression({{r_v, r}, {t_v, t}, {p_v, p}, {"x", r * sp * mymath::cos(t)}, {"y", r * sp * mymath::sin(t)}, {"z", r * mymath::cos(p)}}) * r * r * sp;
    });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([t0, t1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(t0, t1); });
    bounds.push_back([p0, p1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(p0, p1); });
    auto r0_f = make_scalar_bound_func(ctx, r0_e, {t_v, p_v}, 2);
    auto r1_f = make_scalar_bound_func(ctx, r1_e, {t_v, p_v}, 2);
    bounds.push_back([r0_f, r1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(r0_f(pt), r1_f(pt)); });
    return integrator.integrate(bounds, {nt, np, nr});
}

bool is_integration_command(const std::string& command) {
    return command == "double_integral" || command == "triple_integral";
}

bool handle_integration_command(const IntegrationContext& ctx,
                                const std::string& command,
                                const std::vector<std::string>& arguments_orig,
                                std::string* output) {
    std::vector<std::string> arguments = arguments_orig;
    if (arguments.empty()) return false;
    if (command == "integral") {
        if (arguments.size() < 3) throw std::runtime_error("integral expects expr, x0, x1");
        std::string x_var = "x"; double x0, x1;
        if (arguments.size() >= 4 && is_identifier_text(utils::trim_copy(arguments[1]))) {
            x_var = utils::trim_copy(arguments[1]); x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]);
        } else {
            x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]);
        }
        auto f = ctx.build_scoped_evaluator(arguments[0]);
        MultivariableIntegrator integrator([&f, x_var](const std::vector<double>& pt) { return f({{ x_var, pt[0] }}); });
        std::vector<MultivariableIntegrator::BoundFunc> bounds = { [x0, x1](const std::vector<double>&) -> std::pair<double, double> { return {x0, x1}; } };
        *output = format_decimal(ctx.normalize_result(integrator.integrate(bounds, {1024})));
        return true;
    }
    std::string coord_system = "rect";
    if (arguments.size() >= 2) {
        std::string last = utils::trim_copy(arguments.back());
        if (last == "\"polar\"" || last == "polar" || last == "\"cyl\"" || last == "cyl" || last == "\"sph\"" || last == "sph" || last == "\"rect\"" || last == "rect") {
            coord_system = last; if (coord_system.front() == '"') coord_system = coord_system.substr(1, coord_system.size() - 2);
            arguments.pop_back();
        }
    }
    if (command == "double_integral") {
        if (coord_system == "polar" || coord_system == "cyl") {
            if (arguments.size() >= 7 && is_identifier_text(utils::trim_copy(arguments[1]))) {
                auto subs = parse_subdivisions(ctx, arguments, 7, {64, 64});
                *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], subs[0], subs[1])));
            } else {
                auto subs = parse_subdivisions(ctx, arguments, 5, {64, 64});
                *output = format_decimal(ctx.normalize_result(double_integral_polar(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], subs[1], subs[0])));
            }
        } else {
            std::string xv = "x", yv = "y"; double x0, x1; std::string y0e, y1e; size_t next;
            if (is_identifier_text(utils::trim_copy(arguments[1])) && arguments.size() >= 7) { xv = utils::trim_copy(arguments[1]); x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]); yv = utils::trim_copy(arguments[4]); y0e = arguments[5]; y1e = arguments[6]; next = 7; }
            else { x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]); y0e = arguments[3]; y1e = arguments[4]; next = 5; }
            auto subs = parse_subdivisions(ctx, arguments, next, {32, 32});
            *output = format_decimal(ctx.normalize_result(double_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, subs[0], subs[1])));
        }
        return true;
    }
    if (command == "triple_integral") {
        if (coord_system == "cyl") {
            if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[1]))) {
                auto s = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_cyl(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], arguments[7], arguments[8], arguments[9], s[0], s[1], s[2])));
            } else {
                auto s = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_cyl(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], "z", arguments[5], arguments[6], s[1], s[0], s[2])));
            }
        } else if (coord_system == "sph") {
            if (arguments.size() >= 10 && is_identifier_text(utils::trim_copy(arguments[1]))) {
                auto s = parse_subdivisions(ctx, arguments, 10, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), arguments[7], arguments[8], arguments[9], s[0], s[1], s[2])));
            } else {
                auto s = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
                *output = format_decimal(ctx.normalize_result(triple_integral_sph(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "phi", ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), "rho", arguments[1], arguments[2], s[1], s[2], s[0])));
            }
        } else {
            std::string xv = "x", yv = "y", zv = "z"; double x0, x1; std::string y0e, y1e, z0e, z1e; size_t next;
            if (is_identifier_text(utils::trim_copy(arguments[1])) && arguments.size() >= 10) { xv = utils::trim_copy(arguments[1]); x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]); yv = utils::trim_copy(arguments[4]); y0e = arguments[5]; y1e = arguments[6]; zv = utils::trim_copy(arguments[7]); z0e = arguments[8]; z1e = arguments[9]; next = 10; }
            else { x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]); y0e = arguments[3]; y1e = arguments[4]; z0e = arguments[5]; z1e = arguments[6]; next = 7; }
            auto s = parse_subdivisions(ctx, arguments, next, {16, 16, 16});
            *output = format_decimal(ctx.normalize_result(triple_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, zv, z0e, z1e, s[0], s[1], s[2])));
        }
        return true;
    }
    return false;
}

bool handle_integration_command(const IntegrationContext& ctx, const std::string& command, const std::string& inside, std::string* output) {
    return handle_integration_command(ctx, command, split_top_level_arguments(inside), output);
}

bool IntegrationModule::can_handle(const std::string& command) const { return is_integration_command(command); }

std::string IntegrationModule::execute_args(const std::string& command, const std::vector<std::string>& args, const CoreServices& services) {
    IntegrationContext ctx;
    ctx.parse_decimal = services.evaluation.parse_decimal;
    ctx.build_scoped_evaluator = services.evaluation.build_decimal_evaluator;
    ctx.normalize_result = services.evaluation.normalize_result;
    std::string out;
    if (handle_integration_command(ctx, command, args, &out)) return out;
    throw std::runtime_error("Unknown integration command: " + command);
}

std::string IntegrationModule::get_help_snippet(const std::string& topic) const {
    if (topic == "analysis") {
        return "Numerical Integration:\n"
               "  double_integral(f, x0, x1, y0, y1, [\"polar\"]) 2D numerical integration\n"
               "  triple_integral(f, x0, x1, y0, y1, z0, z1, [\"cyl\"|\"sph\"]) 3D numerical integration";
    }
    return "";
}

std::vector<std::string> IntegrationModule::get_commands() const {
    return {"double_integral", "double_integral_cyl", "triple_integral", "triple_integral_sph"};
}

}  // namespace integration_ops

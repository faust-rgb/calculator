// ============================================================================
// 多重积分命令实现
// ============================================================================

#include "calculator_integration.h"

#include "mymath.h"

#include <stdexcept>
#include <vector>

namespace integration_ops {

namespace {

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

// Fixed: helper function to always return a double instead of potentially ambiguous lambda
std::function<double(const std::vector<double>&)> make_scalar_bound_func(
    const IntegrationContext& ctx,
    const std::string& bound_expr,
    const std::vector<std::string>& var_names,
    std::size_t current_dim) {
    
    try {
        double val = ctx.parse_decimal(bound_expr);
        return [val](const std::vector<double>&) { return val; };
    } catch (...) {
        // Fallback to evaluator
    }

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

double double_integral_polar(
    const IntegrationContext& ctx,
    const std::string& expr,
    double r0, double r1,
    double theta0, double theta1,
    int nr, int ntheta) {

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression](const std::vector<double>& point) {
            const double r = point[0];
            const double theta = point[1];
            const double x = r * mymath::cos(theta);
            const double y = r * mymath::sin(theta);
            return evaluate_expression(
                       {{"r", r}, {"theta", theta}, {"x", x}, {"y", y}}) *
                   r;
        });
    return integrator.integrate({
        [r0, r1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(r0, r1); },
        [theta0, theta1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(theta0, theta1); }
    }, {nr, ntheta});
}

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

double triple_integral_cyl(
    const IntegrationContext& ctx,
    const std::string& expr,
    double r0, double r1,
    double theta0, double theta1,
    double z0, double z1,
    int nr, int ntheta, int nz) {

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression](const std::vector<double>& point) {
            const double r = point[0];
            const double theta = point[1];
            const double z = point[2];
            const double x = r * mymath::cos(theta);
            const double y = r * mymath::sin(theta);
            return evaluate_expression(
                       {{"r", r}, {"theta", theta}, {"z", z}, {"x", x}, {"y", y}}) *
                   r;
        });
    return integrator.integrate({
        [r0, r1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(r0, r1); },
        [theta0, theta1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(theta0, theta1); },
        [z0, z1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(z0, z1); }
    }, {nr, ntheta, nz});
}

double triple_integral_sph(
    const IntegrationContext& ctx,
    const std::string& expr,
    double rho0, double rho1,
    double theta0, double theta1,
    double phi0, double phi1,
    int nrho, int ntheta, int nphi) {

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression](const std::vector<double>& point) {
            const double rho = point[0];
            const double theta = point[1];
            const double phi = point[2];
            const double sin_phi = mymath::sin(phi);
            const double cos_phi = mymath::cos(phi);
            const double cos_theta = mymath::cos(theta);
            const double sin_theta = mymath::sin(theta);
            const double x = rho * sin_phi * cos_theta;
            const double y = rho * sin_phi * sin_theta;
            const double z = rho * cos_phi;
            return evaluate_expression(
                       {{"rho", rho}, {"theta", theta}, {"phi", phi},
                        {"x", x}, {"y", y}, {"z", z}}) *
                   rho * rho * sin_phi;
        });
    return integrator.integrate({
        [rho0, rho1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(rho0, rho1); },
        [theta0, theta1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(theta0, theta1); },
        [phi0, phi1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(phi0, phi1); }
    }, {nrho, ntheta, nphi});
}

bool is_integration_command(const std::string& command) {
    return command == "double_integral" ||
           command == "double_integral_cyl" ||
           command == "double_integral_polar" ||
           command == "triple_integral" ||
           command == "triple_integral_cyl" ||
           command == "triple_integral_sph";
}

bool handle_integration_command(const IntegrationContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "double_integral") {
        if (arguments.size() < 5) {
            throw std::runtime_error(
                "double_integral expects expr, [x,] x0, x1, [y,] y0, y1, [nx, ny]");
        }
        
        std::string x_var = "x", y_var = "y";
        double x0, x1;
        std::string y0_expr, y1_expr;
        std::size_t next_idx = 0;
        
        bool custom_vars = false;
        try {
            (void)ctx.parse_decimal(arguments[1]);
        } catch (...) {
            if (is_identifier_text(arguments[1])) {
                custom_vars = true;
            }
        }

        if (custom_vars) {
            x_var = arguments[1];
            x0 = ctx.parse_decimal(arguments[2]);
            x1 = ctx.parse_decimal(arguments[3]);
            y_var = arguments[4];
            y0_expr = arguments[5];
            y1_expr = arguments[6];
            next_idx = 7;
        } else {
            x0 = ctx.parse_decimal(arguments[1]);
            x1 = ctx.parse_decimal(arguments[2]);
            y0_expr = arguments[3];
            y1_expr = arguments[4];
            next_idx = 5;
        }

        const std::vector<int> subdivisions = parse_subdivisions(ctx, arguments, next_idx, {32, 32});
        double result = double_integral(
            ctx, arguments[0],
            x_var, x0, x1,
            y_var, y0_expr, y1_expr,
            subdivisions[0], subdivisions[1]);
        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    if (command == "double_integral_cyl" || command == "double_integral_polar") {
        if (arguments.size() != 5 && arguments.size() != 7) {
            throw std::runtime_error(
                command + " expects expr, r0, r1, theta0, theta1, and optional nr, ntheta");
        }
        const std::vector<int> subdivisions = parse_subdivisions(ctx, arguments, 5, {32, 32});
        double result = double_integral_polar(
            ctx, arguments[0],
            ctx.parse_decimal(arguments[1]), ctx.parse_decimal(arguments[2]),
            ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]),
            subdivisions[0], subdivisions[1]);
        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    if (command == "triple_integral") {
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
            if (is_identifier_text(arguments[1])) {
                custom_vars = true;
            }
        }

        if (custom_vars) {
            x_var = arguments[1];
            x0 = ctx.parse_decimal(arguments[2]); x1 = ctx.parse_decimal(arguments[3]);
            y_var = arguments[4];
            y0_e = arguments[5]; y1_e = arguments[6];
            z_var = arguments[7];
            z0_e = arguments[8]; z1_e = arguments[9];
            next_idx = 10;
        } else {
            x0 = ctx.parse_decimal(arguments[1]); x1 = ctx.parse_decimal(arguments[2]);
            y0_e = arguments[3]; y1_e = arguments[4];
            z0_e = arguments[5]; z1_e = arguments[6];
            next_idx = 7;
        }

        const std::vector<int> subdivisions = parse_subdivisions(ctx, arguments, next_idx, {16, 16, 16});
        double result = triple_integral(
            ctx, arguments[0],
            x_var, x0, x1,
            y_var, y0_e, y1_e,
            z_var, z0_e, z1_e,
            subdivisions[0], subdivisions[1], subdivisions[2]);
        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    if (command == "triple_integral_cyl") {
        if (arguments.size() != 7 && arguments.size() != 10) {
            throw std::runtime_error(
                "triple_integral_cyl expects expr, r0, r1, theta0, theta1, z0, z1, and optional nr, ntheta, nz");
        }
        const std::vector<int> subdivisions = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
        double result = triple_integral_cyl(
            ctx, arguments[0],
            ctx.parse_decimal(arguments[1]), ctx.parse_decimal(arguments[2]),
            ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]),
            ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]),
            subdivisions[0], subdivisions[1], subdivisions[2]);
        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    if (command == "triple_integral_sph") {
        if (arguments.size() != 7 && arguments.size() != 10) {
            throw std::runtime_error(
                "triple_integral_sph expects expr, rho0, rho1, theta0, theta1, phi0, phi1, and optional nrho, ntheta, nphi");
        }
        const std::vector<int> subdivisions = parse_subdivisions(ctx, arguments, 7, {16, 16, 16});
        double result = triple_integral_sph(
            ctx, arguments[0],
            ctx.parse_decimal(arguments[1]), ctx.parse_decimal(arguments[2]),
            ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]),
            ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]),
            subdivisions[0], subdivisions[1], subdivisions[2]);
        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    return false;
}

}  // namespace integration_ops

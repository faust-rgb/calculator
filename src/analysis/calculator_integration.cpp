// ============================================================================
// 多重积分命令实现
// ============================================================================

#include "calculator_integration.h"

#include "mymath.h"

#include <vector>

namespace integration_ops {

double double_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    double x0, double x1,
    double y0, double y1,
    int nx, int ny) {

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression](const std::vector<double>& point) {
            return evaluate_expression({{"x", point[0]}, {"y", point[1]}});
        });
    return integrator.integrate({{x0, x1}, {y0, y1}}, {nx, ny});
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
    return integrator.integrate({{r0, r1}, {theta0, theta1}}, {nr, ntheta});
}

double triple_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    double x0, double x1,
    double y0, double y1,
    double z0, double z1,
    int nx, int ny, int nz) {

    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    const MultivariableIntegrator integrator(
        [evaluate_expression](const std::vector<double>& point) {
            return evaluate_expression(
                {{"x", point[0]}, {"y", point[1]}, {"z", point[2]}});
        });
    return integrator.integrate({{x0, x1}, {y0, y1}, {z0, z1}}, {nx, ny, nz});
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
    return integrator.integrate({{r0, r1}, {theta0, theta1}, {z0, z1}}, {nr, ntheta, nz});
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
    return integrator.integrate({{rho0, rho1}, {theta0, theta1}, {phi0, phi1}}, {nrho, ntheta, nphi});
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
        if (arguments.size() != 5 && arguments.size() != 7) {
            throw std::runtime_error(
                "double_integral expects expr, x0, x1, y0, y1, and optional nx, ny");
        }
        const std::vector<int> subdivisions = ctx.parse_subdivisions(arguments, 5, {32, 32});
        double result = double_integral(
            ctx, arguments[0],
            ctx.parse_decimal(arguments[1]), ctx.parse_decimal(arguments[2]),
            ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]),
            subdivisions[0], subdivisions[1]);
        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    if (command == "double_integral_cyl" || command == "double_integral_polar") {
        if (arguments.size() != 5 && arguments.size() != 7) {
            throw std::runtime_error(
                command + " expects expr, r0, r1, theta0, theta1, and optional nr, ntheta");
        }
        const std::vector<int> subdivisions = ctx.parse_subdivisions(arguments, 5, {32, 32});
        double result = double_integral_polar(
            ctx, arguments[0],
            ctx.parse_decimal(arguments[1]), ctx.parse_decimal(arguments[2]),
            ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]),
            subdivisions[0], subdivisions[1]);
        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    if (command == "triple_integral") {
        if (arguments.size() != 7 && arguments.size() != 10) {
            throw std::runtime_error(
                "triple_integral expects expr, x0, x1, y0, y1, z0, z1, and optional nx, ny, nz");
        }
        const std::vector<int> subdivisions = ctx.parse_subdivisions(arguments, 7, {16, 16, 16});
        double result = triple_integral(
            ctx, arguments[0],
            ctx.parse_decimal(arguments[1]), ctx.parse_decimal(arguments[2]),
            ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]),
            ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]),
            subdivisions[0], subdivisions[1], subdivisions[2]);
        *output = format_decimal(ctx.normalize_result(result));
        return true;
    }

    if (command == "triple_integral_cyl") {
        if (arguments.size() != 7 && arguments.size() != 10) {
            throw std::runtime_error(
                "triple_integral_cyl expects expr, r0, r1, theta0, theta1, z0, z1, and optional nr, ntheta, nz");
        }
        const std::vector<int> subdivisions = ctx.parse_subdivisions(arguments, 7, {16, 16, 16});
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
        const std::vector<int> subdivisions = ctx.parse_subdivisions(arguments, 7, {16, 16, 16});
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

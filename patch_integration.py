import re

header_file = "src/analysis/calculator_integration.h"
with open(header_file, "r") as f:
    header_content = f.read()

# Add line_integral and surface_integral
header_content = header_content.replace(
    "// Double integral functions",
    "// Line and Surface integral functions\n"
    "double line_integral(const IntegrationContext& ctx, const std::string& expr,\n"
    "                     const std::string& t_var, double t0, double t1,\n"
    "                     const std::string& x_expr, const std::string& y_expr, const std::string& z_expr, int subdivides);\n\n"
    "double surface_integral(const IntegrationContext& ctx, const std::string& expr,\n"
    "                        const std::string& u_var, double u0, double u1,\n"
    "                        const std::string& v_var, double v0, double v1,\n"
    "                        const std::string& x_expr, const std::string& y_expr, const std::string& z_expr, int nu, int nv);\n\n"
    "// Double integral functions"
)

with open(header_file, "w") as f:
    f.write(header_content)

cpp_file = "src/analysis/calculator_integration.cpp"
with open(cpp_file, "r") as f:
    cpp_content = f.read()

# Insert the line_integral and surface_integral implementation right after make_scalar_bound_func
line_surf_impl = """
double line_integral(const IntegrationContext& ctx, const std::string& expr,
                     const std::string& t_var, double t0, double t1,
                     const std::string& x_expr, const std::string& y_expr, const std::string& z_expr,
                     int subdivides) {
    auto f_eval = ctx.build_scoped_evaluator(expr);
    auto x_eval = ctx.build_scoped_evaluator(x_expr);
    auto y_eval = ctx.build_scoped_evaluator(y_expr);
    bool has_z = !z_expr.empty() && z_expr != "0" && z_expr != "\"\"";
    std::function<double(const std::vector<std::pair<std::string, double>>&)> z_eval;
    if (has_z) z_eval = ctx.build_scoped_evaluator(z_expr);

    auto integrand = [&](const std::vector<double>& pt) {
        double t = pt[0];
        double h = 1e-6;
        double x_val = x_eval({{t_var, t}});
        double y_val = y_eval({{t_var, t}});
        double z_val = has_z ? z_eval({{t_var, t}}) : 0.0;
        
        double dx = (x_eval({{t_var, t + h}}) - x_eval({{t_var, t - h}})) / (2 * h);
        double dy = (y_eval({{t_var, t + h}}) - y_eval({{t_var, t - h}})) / (2 * h);
        double dz = has_z ? (z_eval({{t_var, t + h}}) - z_eval({{t_var, t - h}})) / (2 * h) : 0.0;
        
        double ds = mymath::sqrt(dx*dx + dy*dy + dz*dz);
        
        std::vector<std::pair<std::string, double>> scope = {
            {t_var, t}, {"x", x_val}, {"y", y_val}
        };
        if (has_z) scope.push_back({"z", z_val});
        
        return f_eval(scope) * ds;
    };

    MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds = {
        [t0, t1](const std::vector<double>&) -> std::pair<double, double> { return {t0, t1}; }
    };
    return integrator.integrate(bounds, {subdivides});
}

double surface_integral(const IntegrationContext& ctx, const std::string& expr,
                        const std::string& u_var, double u0, double u1,
                        const std::string& v_var, double v0, double v1,
                        const std::string& x_expr, const std::string& y_expr, const std::string& z_expr,
                        int nu, int nv) {
    auto f_eval = ctx.build_scoped_evaluator(expr);
    auto x_eval = ctx.build_scoped_evaluator(x_expr);
    auto y_eval = ctx.build_scoped_evaluator(y_expr);
    auto z_eval = ctx.build_scoped_evaluator(z_expr);

    auto integrand = [&](const std::vector<double>& pt) {
        double u = pt[0];
        double v = pt[1];
        double h = 1e-5;
        
        double x_val = x_eval({{u_var, u}, {v_var, v}});
        double y_val = y_eval({{u_var, u}, {v_var, v}});
        double z_val = z_eval({{u_var, u}, {v_var, v}});
        
        double xu = (x_eval({{u_var, u + h}, {v_var, v}}) - x_eval({{u_var, u - h}, {v_var, v}})) / (2 * h);
        double yu = (y_eval({{u_var, u + h}, {v_var, v}}) - y_eval({{u_var, u - h}, {v_var, v}})) / (2 * h);
        double zu = (z_eval({{u_var, u + h}, {v_var, v}}) - z_eval({{u_var, u - h}, {v_var, v}})) / (2 * h);
        
        double xv = (x_eval({{u_var, u}, {v_var, v + h}}) - x_eval({{u_var, u}, {v_var, v - h}})) / (2 * h);
        double yv = (y_eval({{u_var, u}, {v_var, v + h}}) - y_eval({{u_var, u}, {v_var, v - h}})) / (2 * h);
        double zv = (z_eval({{u_var, u}, {v_var, v + h}}) - z_eval({{u_var, u}, {v_var, v - h}})) / (2 * h);
        
        double cx = yu * zv - zu * yv;
        double cy = zu * xv - xu * zv;
        double cz = xu * yv - yu * xv;
        
        double dS = mymath::sqrt(cx*cx + cy*cy + cz*cz);
        
        return f_eval({{u_var, u}, {v_var, v}, {"x", x_val}, {"y", y_val}, {"z", z_val}}) * dS;
    };

    MultivariableIntegrator integrator(integrand);
    std::vector<MultivariableIntegrator::BoundFunc> bounds = {
        [u0, u1](const std::vector<double>&) -> std::pair<double, double> { return {u0, u1}; },
        [v0, v1](const std::vector<double>&) -> std::pair<double, double> { return {v0, v1}; }
    };
    return integrator.integrate(bounds, {nu, nv});
}

"""

cpp_content = cpp_content.replace(
    "}  // namespace\n\ndouble double_integral(",
    "}  // namespace\n\n" + line_surf_impl + "double double_integral("
)

# Update is_integration_command
cpp_content = cpp_content.replace(
    'return command == "double_integral" || command == "triple_integral";',
    'return command == "double_integral" || command == "triple_integral" || command == "line_integral" || command == "surface_integral";'
)

# Insert handlers for line_integral and surface_integral
handlers = """
    if (command == "line_integral") {
        if (arguments.size() < 6) throw std::runtime_error("line_integral expects expr, t_var, t0, t1, x_expr, y_expr, [z_expr], [subdivides]");
        std::string t_var = utils::trim_copy(arguments[1]);
        double t0 = ctx.parse_decimal(arguments[2]);
        double t1 = ctx.parse_decimal(arguments[3]);
        std::string x_expr = arguments[4];
        std::string y_expr = arguments[5];
        std::string z_expr = "";
        int subdivides = 1024;
        if (arguments.size() >= 7) {
            // Check if arguments[6] is a number (subdivides) or z_expr
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
        double u0 = ctx.parse_decimal(arguments[2]);
        double u1 = ctx.parse_decimal(arguments[3]);
        std::string v_var = utils::trim_copy(arguments[4]);
        double v0 = ctx.parse_decimal(arguments[5]);
        double v1 = ctx.parse_decimal(arguments[6]);
        std::string x_expr = arguments[7];
        std::string y_expr = arguments[8];
        std::string z_expr = arguments[9];
        auto s = parse_subdivisions(ctx, arguments, 10, {32, 32});
        *output = format_decimal(ctx.normalize_result(surface_integral(ctx, arguments[0], u_var, u0, u1, v_var, v0, v1, x_expr, y_expr, z_expr, s[0], s[1])));
        return true;
    }
"""

cpp_content = cpp_content.replace(
    'if (command == "double_integral") {',
    handlers + '    if (command == "double_integral") {'
)

# Update help snippet
cpp_content = cpp_content.replace(
    '               "  triple_integral(f, x0, x1, y0, y1, z0, z1, [\\"cyl\\"|\\"sph\\"]) 3D numerical integration";',
    '               "  triple_integral(f, x0, x1, y0, y1, z0, z1, [\\"cyl\\"|\\"sph\\"]) 3D numerical integration\\n"'
    '               "  line_integral(f, t, t0, t1, x(t), y(t), [z(t)]) Line integral\\n"'
    '               "  surface_integral(f, u, u0, u1, v, v0, v1, x(u,v), y(u,v), z(u,v)) Surface integral";'
)

# Update get_commands
cpp_content = cpp_content.replace(
    'return {"double_integral", "double_integral_cyl", "triple_integral", "triple_integral_sph"};',
    'return {"double_integral", "double_integral_cyl", "triple_integral", "triple_integral_sph", "line_integral", "surface_integral"};'
)

with open(cpp_file, "w") as f:
    f.write(cpp_content)


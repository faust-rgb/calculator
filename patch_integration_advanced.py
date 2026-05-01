import re

header_file = "src/analysis/calculator_integration.h"
with open(header_file, "r") as f:
    header_content = f.read()

# Add method and tol to double_integral
header_content = header_content.replace(
    "                       int nx, int ny);",
    "                       int nx, int ny, const std::string& method = \"simpson\", double tol = 1e-6);"
)

# Add method and tol to double_integral_polar
header_content = header_content.replace(
    "                              int ntheta, int nr);",
    "                              int ntheta, int nr, const std::string& method = \"simpson\", double tol = 1e-6);"
)

# Add method and tol to triple_integral
header_content = header_content.replace(
    "                       int nx, int ny, int nz);",
    "                       int nx, int ny, int nz, const std::string& method = \"simpson\", double tol = 1e-6);"
)

# Add method and tol to triple_integral_cyl
header_content = header_content.replace(
    "                           int ntheta, int nr, int nz);",
    "                           int ntheta, int nr, int nz, const std::string& method = \"simpson\", double tol = 1e-6);"
)

# Add method and tol to triple_integral_sph
header_content = header_content.replace(
    "                           int ntheta, int nphi, int nr);",
    "                           int ntheta, int nphi, int nr, const std::string& method = \"simpson\", double tol = 1e-6);"
)

with open(header_file, "w") as f:
    f.write(header_content)


cpp_file = "src/analysis/calculator_integration.cpp"
with open(cpp_file, "r") as f:
    cpp_content = f.read()

old_parse = """    std::string coord_system = "rect";
    if (arguments.size() >= 2) {
        std::string last = utils::trim_copy(arguments.back());
        if (last == "\\"polar\\"" || last == "polar" || last == "\\"cyl\\"" || last == "cyl" || last == "\\"sph\\"" || last == "sph" || last == "\\"rect\\"" || last == "rect") {
            coord_system = last; if (coord_system.front() == '"') coord_system = coord_system.substr(1, coord_system.size() - 2);
            arguments.pop_back();
        }
    }"""

new_parse = """    std::string coord_system = "rect";
    std::string method = "simpson";
    double tol = 1e-6;
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
    }"""

cpp_content = cpp_content.replace(old_parse, new_parse)

# Update handler calls
cpp_content = cpp_content.replace(
    'double_integral_polar(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], subs[0], subs[1])',
    'double_integral_polar(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], subs[0], subs[1], method, tol)'
)
cpp_content = cpp_content.replace(
    'double_integral_polar(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], subs[1], subs[0])',
    'double_integral_polar(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], subs[1], subs[0], method, tol)'
)
cpp_content = cpp_content.replace(
    'double_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, subs[0], subs[1])',
    'double_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, subs[0], subs[1], method, tol)'
)

cpp_content = cpp_content.replace(
    'triple_integral_cyl(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], arguments[7], arguments[8], arguments[9], s[0], s[1], s[2])',
    'triple_integral_cyl(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], arguments[5], arguments[6], arguments[7], arguments[8], arguments[9], s[0], s[1], s[2], method, tol)'
)
cpp_content = cpp_content.replace(
    'triple_integral_cyl(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], "z", arguments[5], arguments[6], s[1], s[0], s[2])',
    'triple_integral_cyl(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "r", arguments[1], arguments[2], "z", arguments[5], arguments[6], s[1], s[0], s[2], method, tol)'
)
cpp_content = cpp_content.replace(
    'triple_integral_sph(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), arguments[7], arguments[8], arguments[9], s[0], s[1], s[2])',
    'triple_integral_sph(ctx, arguments[0], arguments[1], ctx.parse_decimal(arguments[2]), ctx.parse_decimal(arguments[3]), arguments[4], ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), arguments[7], arguments[8], arguments[9], s[0], s[1], s[2], method, tol)'
)
cpp_content = cpp_content.replace(
    'triple_integral_sph(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "phi", ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), "rho", arguments[1], arguments[2], s[1], s[2], s[0])',
    'triple_integral_sph(ctx, arguments[0], "theta", ctx.parse_decimal(arguments[3]), ctx.parse_decimal(arguments[4]), "phi", ctx.parse_decimal(arguments[5]), ctx.parse_decimal(arguments[6]), "rho", arguments[1], arguments[2], s[1], s[2], s[0], method, tol)'
)
cpp_content = cpp_content.replace(
    'triple_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, zv, z0e, z1e, s[0], s[1], s[2])',
    'triple_integral(ctx, arguments[0], xv, x0, x1, yv, y0e, y1e, zv, z0e, z1e, s[0], s[1], s[2], method, tol)'
)


# Now update the implementations to use method and tol.
# Let's replace the double_integral body.

def replace_impl(func_name, old_args, inner_block):
    global cpp_content
    pattern = r"double " + func_name + r"\([^)]+\)\s*\{[^{]*\{[^{]*\}[^}]*\}[^{]*\{[^{]*\}[^{]*\}[^}]*\}"
    # Use a generic regex to replace the function body
    pass

# We will just use string replacement on the exact function definition
old_double_integral = """double double_integral(
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
}"""

new_double_integral = """double double_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    const std::string& x_var, double x0, double x1,
    const std::string& y_var, const std::string& y0_expr, const std::string& y1_expr,
    int nx, int ny, const std::string& method, double tol) {
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    auto y0_f = make_scalar_bound_func(ctx, y0_expr, {x_var}, 1);
    auto y1_f = make_scalar_bound_func(ctx, y1_expr, {x_var}, 1);
    
    bool is_constant_bounds = false;
    double y0_c = 0.0, y1_c = 0.0;
    try {
        y0_c = ctx.parse_decimal(y0_expr);
        y1_c = ctx.parse_decimal(y1_expr);
        is_constant_bounds = true;
    } catch (...) {}
    
    if (is_constant_bounds && method != "simpson") {
        auto integrand = [evaluate_expression, x_var, y_var](const std::vector<double>& pt) {
            return evaluate_expression({{x_var, pt[0]}, {y_var, pt[1]}});
        };
        std::vector<multidim::IntegrationBounds> rect_bounds = {{x0, x1}, {y0_c, y1_c}};
        multidim::IntegrationOptions opts;
        opts.relative_tolerance = tol;
        opts.absolute_tolerance = 1e-12;
        opts.max_evaluations = 500000;
        if (method == "adaptive") opts.method = multidim::IntegrationMethod::Adaptive;
        else if (method == "monte_carlo") opts.method = multidim::IntegrationMethod::MonteCarlo;
        else if (method == "sparse_grid") opts.method = multidim::IntegrationMethod::SparseGrid;
        else if (method == "tensor_product") opts.method = multidim::IntegrationMethod::TensorProduct;
        return multidim::integrate_rectangular(integrand, rect_bounds, opts).value;
    }

    const MultivariableIntegrator integrator(
        [evaluate_expression, x_var, y_var](const std::vector<double>& point) {
            return evaluate_expression({{x_var, point[0]}, {y_var, point[1]}});
        });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([x0, x1](const std::vector<double>&) -> std::pair<double, double> { return {x0, x1}; });
    bounds.push_back([y0_f, y1_f](const std::vector<double>& pt) -> std::pair<double, double> { return {y0_f(pt), y1_f(pt)}; });
    return integrator.integrate(bounds, {nx, ny});
}"""
cpp_content = cpp_content.replace(old_double_integral, new_double_integral)

old_triple_integral = """double triple_integral(const IntegrationContext& ctx, const std::string& expr, const std::string& x_v, double x0, double x1, const std::string& y_v, const std::string& y0_e, const std::string& y1_e, const std::string& z_v, const std::string& z0_e, const std::string& z1_e, int nx, int ny, int nz) {
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
}"""

new_triple_integral = """double triple_integral(const IntegrationContext& ctx, const std::string& expr, const std::string& x_v, double x0, double x1, const std::string& y_v, const std::string& y0_e, const std::string& y1_e, const std::string& z_v, const std::string& z0_e, const std::string& z1_e, int nx, int ny, int nz, const std::string& method, double tol) {
    const auto evaluate_expression = ctx.build_scoped_evaluator(expr);
    auto y0_f = make_scalar_bound_func(ctx, y0_e, {x_v}, 1);
    auto y1_f = make_scalar_bound_func(ctx, y1_e, {x_v}, 1);
    auto z0_f = make_scalar_bound_func(ctx, z0_e, {x_v, y_v}, 2);
    auto z1_f = make_scalar_bound_func(ctx, z1_e, {x_v, y_v}, 2);
    
    bool is_constant_bounds = false;
    double y0_c = 0.0, y1_c = 0.0, z0_c = 0.0, z1_c = 0.0;
    try {
        y0_c = ctx.parse_decimal(y0_e); y1_c = ctx.parse_decimal(y1_e);
        z0_c = ctx.parse_decimal(z0_e); z1_c = ctx.parse_decimal(z1_e);
        is_constant_bounds = true;
    } catch (...) {}
    
    if (is_constant_bounds && method != "simpson") {
        auto integrand = [evaluate_expression, x_v, y_v, z_v](const std::vector<double>& pt) {
            return evaluate_expression({{x_v, pt[0]}, {y_v, pt[1]}, {z_v, pt[2]}});
        };
        std::vector<multidim::IntegrationBounds> rect_bounds = {{x0, x1}, {y0_c, y1_c}, {z0_c, z1_c}};
        multidim::IntegrationOptions opts;
        opts.relative_tolerance = tol;
        opts.absolute_tolerance = 1e-12;
        opts.max_evaluations = 1000000;
        if (method == "adaptive") opts.method = multidim::IntegrationMethod::Adaptive;
        else if (method == "monte_carlo") opts.method = multidim::IntegrationMethod::MonteCarlo;
        else if (method == "sparse_grid") opts.method = multidim::IntegrationMethod::SparseGrid;
        else if (method == "tensor_product") opts.method = multidim::IntegrationMethod::TensorProduct;
        return multidim::integrate_rectangular(integrand, rect_bounds, opts).value;
    }

    const MultivariableIntegrator integrator([evaluate_expression, x_v, y_v, z_v](const std::vector<double>& pt) { return evaluate_expression({{x_v, pt[0]}, {y_v, pt[1]}, {z_v, pt[2]}}); });
    std::vector<MultivariableIntegrator::BoundFunc> bounds;
    bounds.push_back([x0, x1](const std::vector<double>&) -> std::pair<double, double> { return std::make_pair(x0, x1); });
    bounds.push_back([y0_f, y1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(y0_f(pt), y1_f(pt)); });
    bounds.push_back([z0_f, z1_f](const std::vector<double>& pt) -> std::pair<double, double> { return std::make_pair(z0_f(pt), z1_f(pt)); });
    return integrator.integrate(bounds, {nx, ny, nz});
}"""
cpp_content = cpp_content.replace(old_triple_integral, new_triple_integral)

# Also update polar, cyl and sph since we changed the signatures
cpp_content = cpp_content.replace(
    'int ntheta, int nr) {',
    'int ntheta, int nr, const std::string& method, double tol) {'
)
cpp_content = cpp_content.replace(
    'int nt, int nr, int nz) {',
    'int nt, int nr, int nz, const std::string& method, double tol) {'
)
cpp_content = cpp_content.replace(
    'int nt, int np, int nr) {',
    'int nt, int np, int nr, const std::string& method, double tol) {'
)

with open(cpp_file, "w") as f:
    f.write(cpp_content)


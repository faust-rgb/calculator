#ifndef CALCULATOR_INTEGRATION_H
#define CALCULATOR_INTEGRATION_H

#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"

namespace integration_ops {

/**
 * @class IntegrationModule
 * @brief 提供数值积分功能（一重、二重、三重积分）的模块
 */
class IntegrationModule : public CalculatorModule {
public:
    std::string name() const override { return "Integration"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct IntegrationContext {
    std::function<long double(const std::string&)> parse_decimal;
    std::function<std::function<long double(const std::vector<std::pair<std::string, long double>>&)>(const std::string&)> build_scoped_evaluator;
    std::function<long double(double)> normalize_result;
};

// Line and Surface integral functions
long double line_integral(const IntegrationContext& ctx, const std::string& expr,
                     const std::string& t_var, long double t0, long double t1,
                     const std::string& x_expr, const std::string& y_expr, const std::string& z_expr, int subdivides);

long double surface_integral(const IntegrationContext& ctx, const std::string& expr,
                        const std::string& u_var, long double u0, long double u1,
                        const std::string& v_var, long double v0, long double v1,
                        const std::string& x_expr, const std::string& y_expr, const std::string& z_expr, int nu, int nv);

// long double integral functions
long double double_integral(const IntegrationContext& ctx, const std::string& expr,
                       const std::string& x_v, long double x0, long double x1,
                       const std::string& y_v, const std::string& y0_e, const std::string& y1_e,
                       int nx, int ny, const std::string& method = "simpson", long double tol = 1e-6);

long double double_integral_polar(const IntegrationContext& ctx, const std::string& expr,
                              const std::string& theta_v, long double theta0, long double theta1,
                              const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                              int ntheta, int nr, const std::string& method = "simpson", long double tol = 1e-6);

// Triple integral functions
long double triple_integral(const IntegrationContext& ctx, const std::string& expr,
                       const std::string& x_v, long double x0, long double x1,
                       const std::string& y_v, const std::string& y0_e, const std::string& y1_e,
                       const std::string& z_v, const std::string& z0_e, const std::string& z1_e,
                       int nx, int ny, int nz, const std::string& method = "simpson", long double tol = 1e-6);

long double triple_integral_cyl(const IntegrationContext& ctx, const std::string& expr,
                           const std::string& theta_v, long double theta0, long double theta1,
                           const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                           const std::string& z_v, const std::string& z0_e, const std::string& z1_e,
                           int ntheta, int nr, int nz, const std::string& method = "simpson", long double tol = 1e-6);

long double triple_integral_sph(const IntegrationContext& ctx, const std::string& expr,
                           const std::string& theta_v, long double theta0, long double theta1,
                           const std::string& phi_v, long double phi0, long double phi1,
                           const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                           int ntheta, int nphi, int nr, const std::string& method = "simpson", long double tol = 1e-6);

bool is_integration_command(const std::string& command);

bool handle_integration_command(const IntegrationContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output);

}  // namespace integration_ops

#endif  // CALCULATOR_INTEGRATION_H

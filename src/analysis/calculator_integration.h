#ifndef CALCULATOR_INTEGRATION_H
#define CALCULATOR_INTEGRATION_H

#include <string>
#include <functional>
#include <vector>
#include "../core/calculator_module.h"

namespace integration_ops {

/**
 * @class IntegrationModule
 * @brief 提供数值积分功能（一重、二重、三重积分）的模块
 */
class IntegrationModule : public CalculatorModule {
public:
    std::string name() const override { return "Integration"; }
    std::vector<std::string> get_commands() const override;
    bool can_handle(const std::string& command) const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct IntegrationContext {
    std::function<double(const std::string&)> parse_decimal;
    std::function<std::function<double(const std::vector<std::pair<std::string, double>>&)>(const std::string&)> build_scoped_evaluator;
    std::function<double(double)> normalize_result;
};

// Double integral functions
double double_integral(const IntegrationContext& ctx, const std::string& expr,
                       const std::string& x_v, double x0, double x1,
                       const std::string& y_v, const std::string& y0_e, const std::string& y1_e,
                       int nx, int ny);

double double_integral_polar(const IntegrationContext& ctx, const std::string& expr,
                              const std::string& theta_v, double theta0, double theta1,
                              const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                              int ntheta, int nr);

// Triple integral functions
double triple_integral(const IntegrationContext& ctx, const std::string& expr,
                       const std::string& x_v, double x0, double x1,
                       const std::string& y_v, const std::string& y0_e, const std::string& y1_e,
                       const std::string& z_v, const std::string& z0_e, const std::string& z1_e,
                       int nx, int ny, int nz);

double triple_integral_cyl(const IntegrationContext& ctx, const std::string& expr,
                           const std::string& theta_v, double theta0, double theta1,
                           const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                           const std::string& z_v, const std::string& z0_e, const std::string& z1_e,
                           int ntheta, int nr, int nz);

double triple_integral_sph(const IntegrationContext& ctx, const std::string& expr,
                           const std::string& theta_v, double theta0, double theta1,
                           const std::string& phi_v, double phi0, double phi1,
                           const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                           int ntheta, int nphi, int nr);

bool is_integration_command(const std::string& command);

bool handle_integration_command(const IntegrationContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output);

}  // namespace integration_ops

#endif  // CALCULATOR_INTEGRATION_H

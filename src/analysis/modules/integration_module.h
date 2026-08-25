/**
 * @file integration_module.h
 * @brief 数值积分模块
 *
 * 本文件定义了数值积分模块：
 * - 一重积分：Gauss-Kronrod 自适应积分
 * - 二重积分：嵌套的一重积分
 * - 三重积分：嵌套的二重积分
 * - 多重积分：蒙特卡洛方法
 */

#ifndef INTEGRATION_MODULE_H
#define INTEGRATION_MODULE_H

#include <string>
#include <functional>
#include <vector>
#include "module/calculator_module.h"

// 前向声明
class FunctionAnalysis;
class ServiceLocator;
class SymbolicExpression;

namespace integration_ops {

/**
 * @class IntegrationModule
 * @brief 提供数值积分功能（一重、二重、三重积分）的模块
 */
class IntegrationModule : public CommandModuleBase {
public:
    std::string name() const override { return "Integration"; }
    std::vector<std::string> get_commands() const override;
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ServiceLocator& locator) override;
    std::string get_help_snippet(const std::string& topic) const override;
};

struct IntegrationContext {
    std::function<Scalar(const std::string&)> parse_decimal;
    std::function<std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)>(const std::string&)> build_scoped_evaluator;
    std::function<Scalar(Scalar)> normalize_result;
    std::function<FunctionAnalysis(const std::string&)> build_analysis;
    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)> resolve_symbolic;
};

// Line and Surface integral functions
Scalar line_integral(const IntegrationContext& ctx, const std::string& expr,
                     const std::string& t_var, Scalar t0, Scalar t1,
                     const std::string& x_expr, const std::string& y_expr, const std::string& z_expr, int subdivides);

Scalar surface_integral(const IntegrationContext& ctx, const std::string& expr,
                        const std::string& u_var, Scalar u0, Scalar u1,
                        const std::string& v_var, Scalar v0, Scalar v1,
                        const std::string& x_expr, const std::string& y_expr, const std::string& z_expr, int nu, int nv);

// Scalar integral functions
Scalar double_integral(const IntegrationContext& ctx, const std::string& expr,
                       const std::string& x_v, Scalar x0, Scalar x1,
                       const std::string& y_v, const std::string& y0_e, const std::string& y1_e,
                       int nx, int ny, const std::string& method = "simpson", Scalar tol = 1e-6);

Scalar double_integral_polar(const IntegrationContext& ctx, const std::string& expr,
                              const std::string& theta_v, Scalar theta0, Scalar theta1,
                              const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                              int ntheta, int nr, const std::string& method = "simpson", Scalar tol = 1e-6);

// Triple integral functions
Scalar triple_integral(const IntegrationContext& ctx, const std::string& expr,
                       const std::string& x_v, Scalar x0, Scalar x1,
                       const std::string& y_v, const std::string& y0_e, const std::string& y1_e,
                       const std::string& z_v, const std::string& z0_e, const std::string& z1_e,
                       int nx, int ny, int nz, const std::string& method = "simpson", Scalar tol = 1e-6);

Scalar triple_integral_cyl(const IntegrationContext& ctx, const std::string& expr,
                           const std::string& theta_v, Scalar theta0, Scalar theta1,
                           const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                           const std::string& z_v, const std::string& z0_e, const std::string& z1_e,
                           int ntheta, int nr, int nz, const std::string& method = "simpson", Scalar tol = 1e-6);

Scalar triple_integral_sph(const IntegrationContext& ctx, const std::string& expr,
                           const std::string& theta_v, Scalar theta0, Scalar theta1,
                           const std::string& phi_v, Scalar phi0, Scalar phi1,
                           const std::string& r_v, const std::string& r0_e, const std::string& r1_e,
                           int ntheta, int nphi, int nr, const std::string& method = "simpson", Scalar tol = 1e-6);

bool is_integration_command(const std::string& command);

bool handle_integration_command(const IntegrationContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output);

}  // namespace integration_ops

#endif  // INTEGRATION_MODULE_H

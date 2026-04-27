// ============================================================================
// 多重积分命令
// ============================================================================
//
// 提供多重积分命令的计算逻辑，包括：
// - 二重积分 (double_integral)
// - 极坐标二重积分 (double_integral_polar, double_integral_cyl)
// - 三重积分 (triple_integral)
// - 柱坐标三重积分 (triple_integral_cyl)
// - 球坐标三重积分 (triple_integral_sph)

#ifndef CALCULATOR_INTEGRATION_H
#define CALCULATOR_INTEGRATION_H

#include "calculator_internal_types.h"

#include "multivariable_integrator.h"

#include <string>
#include <functional>
#include <vector>

namespace integration_ops {

// ============================================================================
// 积分上下文
// ============================================================================

/**
 * @brief 多重积分计算上下文
 */
struct IntegrationContext {
    // 数值求值
    std::function<double(const std::string&)> parse_decimal;

    // 带作用域的求值器
    std::function<std::function<double(const std::vector<std::pair<std::string, double>>&)>(const std::string&)>
        build_scoped_evaluator;

    // 结果归一化
    std::function<double(double)> normalize_result;

    // 解析分割数
    std::function<std::vector<int>(const std::vector<std::string>&, std::size_t, const std::vector<int>&)>
        parse_subdivisions;
};

// ============================================================================
// 积分函数
// ============================================================================

/**
 * @brief 二重积分
 */
double double_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    double x0, double x1,
    double y0, double y1,
    int nx, int ny);

/**
 * @brief 极坐标二重积分
 */
double double_integral_polar(
    const IntegrationContext& ctx,
    const std::string& expr,
    double r0, double r1,
    double theta0, double theta1,
    int nr, int ntheta);

/**
 * @brief 三重积分
 */
double triple_integral(
    const IntegrationContext& ctx,
    const std::string& expr,
    double x0, double x1,
    double y0, double y1,
    double z0, double z1,
    int nx, int ny, int nz);

/**
 * @brief 柱坐标三重积分
 */
double triple_integral_cyl(
    const IntegrationContext& ctx,
    const std::string& expr,
    double r0, double r1,
    double theta0, double theta1,
    double z0, double z1,
    int nr, int ntheta, int nz);

/**
 * @brief 球坐标三重积分
 */
double triple_integral_sph(
    const IntegrationContext& ctx,
    const std::string& expr,
    double rho0, double rho1,
    double theta0, double theta1,
    double phi0, double phi1,
    int nrho, int ntheta, int nphi);

// ============================================================================
// 命令处理
// ============================================================================

/**
 * @brief 检查是否为多重积分命令
 */
bool is_integration_command(const std::string& command);

/**
 * @brief 处理多重积分命令
 */
bool handle_integration_command(const IntegrationContext& ctx,
                                const std::string& command,
                                const std::string& inside,
                                std::string* output);

}  // namespace integration_ops

#endif  // CALCULATOR_INTEGRATION_H

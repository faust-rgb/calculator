// ============================================================================
// 多维数值积分
// ============================================================================
//
// 本文件实现了多维数值积分功能，支持：
// - 矩形区域上的张量积积分
// - 非矩形区域上的积分（通过变换或蒙特卡洛方法）
// - 自适应积分策略
// - 高维问题的稀疏网格方法 (Smolyak)
// - 隐式区域积分（通过区域约束）
// - 误差估计和自动精度控制
// - 准蒙特卡洛方法 (Sobol 序列)
//
// 算法选择策略：
// - 低维 (d ≤ 2): 自适应 Gauss-Kronrod
// - 中维 (d ≤ 4): 张量积 Gauss 积分
// - 高维 (d > 4): 蒙特卡洛或稀疏网格

#ifndef MULTIDIM_INTEGRATION_H
#define MULTIDIM_INTEGRATION_H

#include <functional>
#include <vector>

namespace multidim {

// ============================================================================
// 数据结构
// ============================================================================

/**
 * @struct IntegrationBounds
 * @brief 积分区域边界
 */
struct IntegrationBounds {
    double lower = 0.0;  ///< 下界
    double upper = 0.0;  ///< 上界

    IntegrationBounds(double l, double u) : lower(l), upper(u) {}
};

/**
 * @struct IntegrationResult
 * @brief 积分结果
 */
struct IntegrationResult {
    double value = 0.0;         ///< 积分值
    double error_estimate = 0.0; ///< 误差估计
    int function_evaluations = 0; ///< 函数求值次数
    bool converged = false;     ///< 是否收敛
    int subdivisions_used = 0;  ///< 使用的细分数
};

/**
 * @enum IntegrationMethod
 * @brief 积分方法
 */
enum class IntegrationMethod {
    Adaptive,      ///< 自适应 Gauss-Kronrod (低维)
    TensorProduct, ///< 张量积 Gauss 积分 (中维)
    MonteCarlo,    ///< 蒙特卡洛积分 (高维)
    SparseGrid,    ///< 稀疏网格积分 (高维)
    QuasiMonteCarlo ///< 准蒙特卡洛 (Sobol 序列)
};

/**
 * @struct IntegrationOptions
 * @brief 积分选项
 */
struct IntegrationOptions {
    double relative_tolerance = 1e-8;  ///< 相对容差
    double absolute_tolerance = 1e-10; ///< 绝对容差
    int max_evaluations = 100000;      ///< 最大函数求值次数
    int max_depth = 20;                ///< 自适应细分最大深度
    int initial_subdivisions = 16;     ///< 初始分割数
    IntegrationMethod method = IntegrationMethod::Adaptive; ///< 积分方法
    bool return_error_estimate = true; ///< 是否返回误差估计
};

// ============================================================================
// 多维函数类型
// ============================================================================

/// 多维函数类型: f(x1, x2, ..., xd)
using MultidimFunction = std::function<double(const std::vector<double>&)>;

/// 区域约束函数类型: g(x1, x2, ..., xd) <= 0 表示在区域内
using RegionConstraint = std::function<double(const std::vector<double>&)>;

// ============================================================================
// 矩形区域积分
// ============================================================================

/**
 * @brief 矩形区域上的多维积分
 *
 * 使用自适应 Gauss-Kronrod 方法（低维）或张量积 Gauss 方法（中维）。
 *
 * @param function 多维函数
 * @param bounds 各维度的边界列表
 * @param relative_tolerance 相对容差
 * @param absolute_tolerance 绝对容差
 * @param max_evaluations 最大函数求值次数
 * @return 积分结果
 */
IntegrationResult integrate_rectangular(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    double relative_tolerance = 1e-8,
    double absolute_tolerance = 1e-10,
    int max_evaluations = 100000);

/**
 * @brief 矩形区域上的多维积分（带选项）
 */
IntegrationResult integrate_rectangular(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const IntegrationOptions& options);

/**
 * @brief 矩形区域上的张量积 Gauss 积分
 *
 * 使用 Gauss-Legendre 积分点的张量积。
 * 适用于中等维度 (d ≤ 4)。
 *
 * @param function 多维函数
 * @param bounds 各维度的边界
 * @param points_per_dimension 每维度的积分点数
 * @return 积分值
 */
double integrate_tensor_product(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    int points_per_dimension = 15);

// ============================================================================
// 非矩形区域积分
// ============================================================================

/**
 * @brief 非矩形区域上的积分（变换法）
 *
 * 通过变量变换将非矩形区域映射到矩形区域。
 * 例如：极坐标变换、球坐标变换等。
 *
 * @param function 多维函数
 * @param transform 变换函数: (u1,...,ud) -> (x1,...,xd)
 * @param jacobian Jacobian 行列式函数
 * @param bounds 变换后变量的边界（矩形）
 * @param points_per_dimension 积分点数
 * @return 积分值
 */
double integrate_with_transform(
    const MultidimFunction& function,
    const std::function<std::vector<double>(const std::vector<double>&)>& transform,
    const std::function<double(const std::vector<double>&)>& jacobian,
    const std::vector<IntegrationBounds>& bounds,
    int points_per_dimension = 15);

/**
 * @brief 非矩形区域上的积分（蒙特卡洛方法）
 *
 * 使用蒙特卡洛方法估计积分值。
 * 适用于高维问题或复杂边界。
 *
 * @param function 多维函数
 * @param bounds 各维度的边界（定义采样区域）
 * @param constraint 区域约束函数（可选）
 * @param num_samples 采样点数
 * @return 积分结果（包含误差估计）
 */
IntegrationResult integrate_monte_carlo(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const RegionConstraint& constraint = RegionConstraint(),
    int num_samples = 10000);

/**
 * @brief 非矩形区域上的积分（重要性采样）
 *
 * 使用重要性采样改进蒙特卡洛效率。
 *
 * @param function 多维函数
 * @param bounds 边界
 * @param constraint 区域约束
 * @param importance_density 重要性密度函数
 * @param num_samples 采样点数
 * @return 积分结果
 */
IntegrationResult integrate_importance_sampling(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const RegionConstraint& constraint,
    const MultidimFunction& importance_density,
    int num_samples = 10000);

/**
 * @brief 准蒙特卡洛积分（使用 Sobol 序列）
 *
 * 使用低差异序列提高收敛速度。
 * 收敛速度 O(1/N) 优于普通蒙特卡洛的 O(1/sqrt(N))。
 *
 * @param function 多维函数
 * @param bounds 边界
 * @param constraint 区域约束（可选）
 * @param num_samples 采样点数
 * @return 积分结果
 */
IntegrationResult integrate_quasi_monte_carlo(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const RegionConstraint& constraint = RegionConstraint(),
    int num_samples = 10000);

// ============================================================================
// 隐式区域积分
// ============================================================================

/**
 * @brief 隐式区域上的积分
 *
 * 计算由不等式 g(x) <= 0 定义的隐式区域上的积分。
 * 使用自适应拒绝采样或水平集方法。
 *
 * @param function 多维函数
 * @param bounds 采样边界框
 * @param constraint 区域约束函数（g(x) <= 0 表示在区域内）
 * @param options 积分选项
 * @return 积分结果
 */
IntegrationResult integrate_implicit_region(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    const RegionConstraint& constraint,
    const IntegrationOptions& options = IntegrationOptions());

// ============================================================================
// 特殊区域积分
// ============================================================================

/**
 * @brief 圆形区域上的积分
 *
 * 使用极坐标变换计算圆形区域上的积分。
 *
 * @param function 二维函数 f(x, y)
 * @param center_x 圆心 x 坐标
 * @param center_y 圆心 y 坐标
 * @param radius 半径
 * @param points_per_dimension 积分点数
 * @return 积分值
 */
double integrate_over_circle(
    const std::function<double(double, double)>& function,
    double center_x,
    double center_y,
    double radius,
    int points_per_dimension = 15);

/**
 * @brief 球形区域上的积分
 *
 * 使用球坐标变换计算球形区域上的积分。
 *
 * @param function 三维函数 f(x, y, z)
 * @param center 球心坐标
 * @param radius 半径
 * @param points_per_dimension 积分点数
 * @return 积分值
 */
double integrate_over_sphere(
    const std::function<double(double, double, double)>& function,
    const std::vector<double>& center,
    double radius,
    int points_per_dimension = 15);

/**
 * @brief 三角形区域上的积分
 *
 * 使用重心坐标变换计算三角形区域上的积分。
 *
 * @param function 二维函数 f(x, y)
 * @param vertices 三角形顶点坐标 {{x1,y1}, {x2,y2}, {x3,y3}}
 * @param points_per_dimension 积分点数
 * @return 积分值
 */
double integrate_over_triangle(
    const std::function<double(double, double)>& function,
    const std::vector<std::vector<double>>& vertices,
    int points_per_dimension = 15);

/**
 * @brief 多边形区域上的积分
 *
 * 将多边形分解为三角形后计算积分。
 *
 * @param function 二维函数 f(x, y)
 * @param vertices 多边形顶点坐标（按顺序）
 * @param points_per_dimension 积分点数
 * @return 积分值
 */
double integrate_over_polygon(
    const std::function<double(double, double)>& function,
    const std::vector<std::vector<double>>& vertices,
    int points_per_dimension = 15);

// ============================================================================
// 稀疏网格积分 (Smolyak)
// ============================================================================

/**
 * @brief 稀疏网格积分
 *
 * 使用 Smolyak 算法构造稀疏网格。
 * 相比张量积网格，显著减少高维积分的求值次数。
 *
 * @param function 多维函数
 * @param bounds 边界
 * @param level 网格层级（控制精度）
 * @return 积分值
 */
double integrate_sparse_grid(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    int level = 5);

/**
 * @brief 稀疏网格积分（带误差估计）
 */
IntegrationResult integrate_sparse_grid_with_error(
    const MultidimFunction& function,
    const std::vector<IntegrationBounds>& bounds,
    int level = 5);

// ============================================================================
// 自适应积分
// ============================================================================

/**
 * @brief 二维自适应积分
 *
 * 使用自适应细分策略计算二维积分。
 *
 * @param function 二维函数
 * @param x_bounds x 方向边界
 * @param y_bounds y 方向边界
 * @param relative_tolerance 相对容差
 * @param max_depth 最大细分深度
 * @return 积分结果
 */
IntegrationResult integrate_2d_adaptive(
    const std::function<double(double, double)>& function,
    const IntegrationBounds& x_bounds,
    const IntegrationBounds& y_bounds,
    double relative_tolerance = 1e-8,
    int max_depth = 15);

/**
 * @brief 三维自适应积分
 */
IntegrationResult integrate_3d_adaptive(
    const std::function<double(double, double, double)>& function,
    const IntegrationBounds& x_bounds,
    const IntegrationBounds& y_bounds,
    const IntegrationBounds& z_bounds,
    double relative_tolerance = 1e-6,
    int max_depth = 10);

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 估计积分误差
 *
 * 使用 Richardson 外推法估计积分误差。
 *
 * @param fine 高精度积分值
 * @param coarse 低精度积分值
 * @param order 方法阶数
 * @return 误差估计
 */
double estimate_error(double fine, double coarse, int order = 4);

/**
 * @brief 自动选择积分方法
 *
 * 根据维度和边界特征自动选择最优积分方法。
 *
 * @param dimension 积分维度
 * @param has_variable_bounds 是否有变边界
 * @param has_constraint 是否有区域约束
 * @return 推荐的积分方法
 */
IntegrationMethod select_method(int dimension, bool has_variable_bounds, bool has_constraint);

}  // namespace multidim

#endif  // MULTIDIM_INTEGRATION_H

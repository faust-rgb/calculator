// ============================================================================
// 常微分方程求解器
// ============================================================================
//
// 本文件实现了常微分方程（ODE）的数值求解器。
// 使用 RKF45（Runge-Kutta-Fehlberg）方法，具有自适应步长控制。
//
// 支持功能：
// - 单方程求解
// - 方程组求解
// - 事件检测（Event Detection）

#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <functional>
#include <utility>
#include <vector>

// ============================================================================
// 数据结构
// ============================================================================

/**
 * @struct ODEPoint
 * @brief 单方程 ODE 解的点
 */
struct ODEPoint {
    double x = 0.0;  ///< 自变量值
    double y = 0.0;  ///< 因变量值
};

/**
 * @struct ODESystemPoint
 * @brief 方程组 ODE 解的点
 */
struct ODESystemPoint {
    double x = 0.0;           ///< 自变量值
    std::vector<double> y;    ///< 因变量向量
};

// ============================================================================
// 单方程求解器
// ============================================================================

/**
 * @class ODESolver
 * @brief 单方程常微分方程求解器
 *
 * 求解形如 dy/dx = f(x, y) 的一阶常微分方程。
 * 使用 RKF45 方法，具有自适应步长控制和事件检测功能。
 */
class ODESolver {
public:
    /// 右端函数类型：f(x, y)
    using RHSFunction = std::function<double(double, double)>;

    /// 事件函数类型：g(x, y)，当 g 变号时触发事件
    using EventFunction = std::function<double(double, double)>;

    /**
     * @brief 构造函数
     * @param rhs 右端函数 f(x, y)
     * @param event 事件函数（可选）
     * @param relative_tolerance 相对容差
     * @param absolute_tolerance 绝对容差
     */
    explicit ODESolver(RHSFunction rhs,
                       EventFunction event = EventFunction(),
                       double relative_tolerance = 1e-9,
                       double absolute_tolerance = 1e-12);

    /**
     * @brief 求解 ODE，返回终点值
     * @param x0 初始自变量
     * @param y0 初始因变量
     * @param x1 终点自变量
     * @param steps 步数（用于确定初始步长）
     * @return 终点处的 y 值
     */
    double solve(double x0, double y0, double x1, int steps = 100) const;

    /**
     * @brief 求解 ODE，返回整个轨迹
     * @param x0 初始自变量
     * @param y0 初始因变量
     * @param x1 终点自变量
     * @param steps 步数
     * @return 解轨迹点列表
     */
    std::vector<ODEPoint> solve_trajectory(double x0,
                                           double y0,
                                           double x1,
                                           int steps = 100) const;

private:
    /// 使用自适应步长积分一个区间
    double integrate_segment(double x0, double y0, double x1) const;

    /// 带事件检测的区间积分
    ODEPoint integrate_segment_with_event(double x0, double y0, double x1, bool* stopped) const;

    /// RKF45 单步
    std::pair<double, double> rkf45_step(double x, double y, double h) const;

    RHSFunction rhs_;                  ///< 右端函数
    EventFunction event_;              ///< 事件函数
    double relative_tolerance_ = 1e-9; ///< 相对容差
    double absolute_tolerance_ = 1e-12;///< 绝对容差
};

// ============================================================================
// 方程组求解器
// ============================================================================

/**
 * @class ODESystemSolver
 * @brief 方程组常微分方程求解器
 *
 * 求解形如 dy/dx = f(x, y) 的一阶常微分方程组。
 * 其中 y 是向量，f 是向量值函数。
 * 使用 RKF45 方法，具有自适应步长控制和事件检测功能。
 */
class ODESystemSolver {
public:
    /// 右端函数类型：f(x, y) 返回向量
    using RHSFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

    /// 事件函数类型：g(x, y) 返回标量
    using EventFunction = std::function<double(double, const std::vector<double>&)>;

    /**
     * @brief 构造函数
     * @param rhs 右端函数
     * @param event 事件函数（可选）
     * @param relative_tolerance 相对容差
     * @param absolute_tolerance 绝对容差
     */
    explicit ODESystemSolver(RHSFunction rhs,
                             EventFunction event = EventFunction(),
                             double relative_tolerance = 1e-9,
                             double absolute_tolerance = 1e-12);

    /**
     * @brief 求解 ODE 方程组，返回终点值
     * @param x0 初始自变量
     * @param y0 初始状态向量
     * @param x1 终点自变量
     * @param steps 步数
     * @return 终点处的状态向量
     */
    std::vector<double> solve(double x0,
                              const std::vector<double>& y0,
                              double x1,
                              int steps = 100) const;

    /**
     * @brief 求解 ODE 方程组，返回整个轨迹
     * @param x0 初始自变量
     * @param y0 初始状态向量
     * @param x1 终点自变量
     * @param steps 步数
     * @return 解轨迹点列表
     */
    std::vector<ODESystemPoint> solve_trajectory(double x0,
                                                 const std::vector<double>& y0,
                                                 double x1,
                                                 int steps = 100) const;

private:
    /// 使用自适应步长积分一个区间
    std::vector<double> integrate_segment(double x0,
                                          const std::vector<double>& y0,
                                          double x1) const;

    /// 带事件检测的区间积分
    ODESystemPoint integrate_segment_with_event(double x0,
                                                const std::vector<double>& y0,
                                                double x1,
                                                bool* stopped) const;

    /// 经典 4 阶 Runge-Kutta 单步
    std::vector<double> rk4_step(double x,
                                 const std::vector<double>& y,
                                 double h) const;

    /// RKF45 单步
    std::pair<std::vector<double>, double> rkf45_step(double x,
                                                       const std::vector<double>& y,
                                                       double h) const;

    RHSFunction rhs_;                  ///< 右端函数
    EventFunction event_;              ///< 事件函数
    double relative_tolerance_ = 1e-9; ///< 相对容差
    double absolute_tolerance_ = 1e-12;///< 绝对容差
};

// ============================================================================
// 刚性 ODE 求解器 (BDF 方法)
// ============================================================================

/**
 * @class StiffODESolver
 * @brief 刚性常微分方程求解器
 *
 * 使用 BDF (Backward Differentiation Formula) 方法求解刚性 ODE。
 * BDF 方法是隐式多步法，对刚性问题具有优良的稳定性。
 *
 * 支持阶数 1-5 的 BDF 方法，自动选择最优阶数。
 * 使用 Newton 迭代求解隐式方程。
 */
class StiffODESolver {
public:
    /// 右端函数类型：f(x, y)
    using RHSFunction = std::function<double(double, double)>;

    /// Jacobian 函数类型：∂f/∂y (可选，若不提供则使用数值近似)
    using JacobianFunction = std::function<double(double, double)>;

    /**
     * @brief 构造函数
     * @param rhs 右端函数 f(x, y)
     * @param jacobian Jacobian 函数 ∂f/∂y (可选)
     * @param relative_tolerance 相对容差
     * @param absolute_tolerance 绝对容差
     */
    explicit StiffODESolver(RHSFunction rhs,
                            JacobianFunction jacobian = JacobianFunction(),
                            double relative_tolerance = 1e-6,
                            double absolute_tolerance = 1e-9);

    /**
     * @brief 求解 ODE，返回终点值
     * @param x0 初始自变量
     * @param y0 初始因变量
     * @param x1 终点自变量
     * @param steps 步数（用于确定初始步长）
     * @return 终点处的 y 值
     */
    double solve(double x0, double y0, double x1, int steps = 100) const;

    /**
     * @brief 求解 ODE，返回整个轨迹
     */
    std::vector<ODEPoint> solve_trajectory(double x0, double y0, double x1, int steps = 100) const;

private:
    /// BDF 单步 (隐式求解)
    double bdf_step(double x, double y, double h, int order,
                    const std::vector<double>& prev_y,
                    const std::vector<double>& prev_h) const;

    /// Newton 迭代求解隐式方程
    double newton_implicit(double x, double y_pred, double h,
                           double gamma, double rhs_val) const;

    /// 数值 Jacobian 近似
    double numerical_jacobian(double x, double y) const;

    RHSFunction rhs_;
    JacobianFunction jacobian_;
    double relative_tolerance_;
    double absolute_tolerance_;
};

/**
 * @class StiffODESystemSolver
 * @brief 刚性常微分方程组求解器
 *
 * 使用 BDF 方法求解刚性 ODE 方程组。
 * 使用 Newton-Krylov 方法求解隐式方程组。
 */
class StiffODESystemSolver {
public:
    using RHSFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;
    using JacobianFunction = std::function<std::vector<std::vector<double>>(double, const std::vector<double>&)>;

    explicit StiffODESystemSolver(RHSFunction rhs,
                                  JacobianFunction jacobian = JacobianFunction(),
                                  double relative_tolerance = 1e-6,
                                  double absolute_tolerance = 1e-9);

    std::vector<double> solve(double x0,
                              const std::vector<double>& y0,
                              double x1,
                              int steps = 100) const;

    std::vector<ODESystemPoint> solve_trajectory(double x0,
                                                  const std::vector<double>& y0,
                                                  double x1,
                                                  int steps = 100) const;

private:
    std::vector<double> bdf_step(double x, const std::vector<double>& y,
                                  double h, int order,
                                  const std::vector<std::vector<double>>& prev_y,
                                  const std::vector<double>& prev_h) const;

    std::vector<double> newton_implicit_system(
        double x, const std::vector<double>& y_pred, double h,
        double gamma, const std::vector<double>& rhs_val) const;

    std::vector<std::vector<double>> numerical_jacobian_matrix(
        double x, const std::vector<double>& y) const;

    RHSFunction rhs_;
    JacobianFunction jacobian_;
    double relative_tolerance_;
    double absolute_tolerance_;
};

#endif

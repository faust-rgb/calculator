/**
 * @file pde_solver.h
 * @brief 偏微分方程 (PDE) 数值求解器头文件
 *
 * 实现了三类经典二阶 PDE 数值求解算法：
 * 1. 抛物型 (Parabolic) - 热传导/扩散方程: 1D/2D Crank-Nicolson 隐式与 ADI 求解器
 * 2. 双曲型 (Hyperbolic) - 波动方程: 1D/2D 时域有限差分法 (FDTD) 显式求解器
 * 3. 椭圆型 (Elliptic) - 泊松/拉普拉斯方程: 五点差分格式 + 逐次超松弛迭代法 (SOR)
 */

#ifndef ANALYSIS_PDE_SOLVER_H
#define ANALYSIS_PDE_SOLVER_H

#include "types/scalar_type.h"
#include "matrix/matrix.h"
#include <functional>
#include <vector>

namespace analysis::pde {

using Scalar = mymath::Scalar;
using Matrix = matrix::Matrix;

// 边界条件类型
enum class BoundaryType {
    kDirichlet,
    kNeumann,
    kPeriodic
};

struct BoundaryCondition1D {
    BoundaryType type_left = BoundaryType::kDirichlet;
    BoundaryType type_right = BoundaryType::kDirichlet;
    std::function<Scalar(Scalar)> val_left = [](Scalar) { return Scalar(0); };
    std::function<Scalar(Scalar)> val_right = [](Scalar) { return Scalar(0); };
};

struct BoundaryCondition2D {
    std::function<Scalar(Scalar, Scalar)> boundary_fn; // 边界上的真值
};

// ============================================================================
// 1. 抛物型方程：热传导 / 扩散方程求解器 (Heat Equation)
// ============================================================================

/**
 * @brief 一维热传导方程 u_t = alpha * u_xx + f(x, t)
 * 采用 Crank-Nicolson 无条件稳定隐式差分法求解
 */
class HeatSolver1D {
public:
    HeatSolver1D(Scalar alpha, Scalar x_min, Scalar x_max, std::size_t nx);

    /**
     * @brief 求解 1D 热传导方程
     * @param u0 初始条件 u(x, 0)
     * @param t_end 演化截止时间
     * @param nt 时间步数
     * @param bc 边界条件
     * @param source 源项 f(x, t)，默认为 0
     * @return 返回大小为 nt x nx 的矩阵，每行对应一个时间步的空间分布
     */
    Matrix solve(std::function<Scalar(Scalar)> u0,
                 Scalar t_end,
                 std::size_t nt,
                 BoundaryCondition1D bc = {},
                 std::function<Scalar(Scalar, Scalar)> source = nullptr) const;

private:
    Scalar alpha_;
    Scalar x_min_, x_max_;
    std::size_t nx_;
};

/**
 * @brief 二维热传导方程 u_t = alpha * (u_xx + u_yy)
 * 采用交替方向隐式法 (ADI) 求解，无条件稳定且时间复杂度 O(Nx * Ny)
 */
class HeatSolver2D {
public:
    HeatSolver2D(Scalar alpha, Scalar x_min, Scalar x_max, std::size_t nx,
                 Scalar y_min, Scalar y_max, std::size_t ny);

    /**
     * @brief 求解 2D 热传导方程到指定时间
     * @param u0 初始条件 u(x, y, 0)
     * @param t_end 演化截止时间
     * @param nt 时间步数
     * @param bc 边界条件
     * @return 返回最终时刻 nx x ny 的空间温度矩阵
     */
    Matrix solve(std::function<Scalar(Scalar, Scalar)> u0,
                 Scalar t_end,
                 std::size_t nt,
                 BoundaryCondition2D bc) const;

private:
    Scalar alpha_;
    Scalar x_min_, x_max_;
    Scalar y_min_, y_max_;
    std::size_t nx_, ny_;
};

// ============================================================================
// 2. 双曲型方程：波动方程求解器 (Wave Equation)
// ============================================================================

/**
 * @brief 一维波动方程 u_tt = c^2 * u_xx
 * 采用时域有限差分法 (FDTD Leapfrog)
 */
class WaveSolver1D {
public:
    WaveSolver1D(Scalar c, Scalar x_min, Scalar x_max, std::size_t nx);

    /**
     * @brief 求解 1D 波动方程
     * @param u0 初始位移 u(x, 0)
     * @param v0 初始速度 u_t(x, 0)
     * @param t_end 演化截止时间
     * @param nt 时间步数
     * @param bc 边界条件
     * @return 返回大小为 nt x nx 的波形演化矩阵
     */
    Matrix solve(std::function<Scalar(Scalar)> u0,
                 std::function<Scalar(Scalar)> v0,
                 Scalar t_end,
                 std::size_t nt,
                 BoundaryCondition1D bc = {}) const;

private:
    Scalar c_;
    Scalar x_min_, x_max_;
    std::size_t nx_;
};

/**
 * @brief 二维波动方程 u_tt = c^2 * (u_xx + u_yy)
 */
class WaveSolver2D {
public:
    WaveSolver2D(Scalar c, Scalar x_min, Scalar x_max, std::size_t nx,
                 Scalar y_min, Scalar y_max, std::size_t ny);

    Matrix solve(std::function<Scalar(Scalar, Scalar)> u0,
                 std::function<Scalar(Scalar, Scalar)> v0,
                 Scalar t_end,
                 std::size_t nt,
                 BoundaryCondition2D bc) const;

private:
    Scalar c_;
    Scalar x_min_, x_max_, y_min_, y_max_;
    std::size_t nx_, ny_;
};

// ============================================================================
// 3. 椭圆型方程：泊松 / 拉普拉斯方程求解器 (Poisson / Laplace Equation)
// ============================================================================

/**
 * @brief 二维泊松方程 u_xx + u_yy = f(x, y) (当 f=0 时为拉普拉斯方程)
 * 采用五点差分格式 + 逐次超松弛迭代法 (SOR)
 */
class PoissonSolver2D {
public:
    PoissonSolver2D(Scalar x_min, Scalar x_max, std::size_t nx,
                    Scalar y_min, Scalar y_max, std::size_t ny);

    /**
     * @brief 求解 2D 泊松/拉普拉斯方程
     * @param source 源项函数 f(x, y)
     * @param bc 边界条件
     * @param tol 相对收敛容差，默认 1e-8
     * @param max_iter 最大迭代步数
     * @return 返回 nx x ny 的电位/稳态温度场解矩阵
     */
    Matrix solve(std::function<Scalar(Scalar, Scalar)> source,
                 BoundaryCondition2D bc,
                 Scalar tol = Scalar(1e-8L),
                 std::size_t max_iter = 5000) const;

private:
    Scalar x_min_, x_max_;
    Scalar y_min_, y_max_;
    std::size_t nx_, ny_;
};

} // namespace analysis::pde

#endif // ANALYSIS_PDE_SOLVER_H

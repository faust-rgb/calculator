/**
 * @file pde_solver.cpp
 * @brief 偏微分方程数值求解器实现
 */

#include "pde_solver.h"
#include "math/mymath.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace analysis::pde {

namespace {

// 三对角线性方程组追赶法求解器 (Thomas Algorithm) - O(N) 复杂度
std::vector<double> solve_tridiagonal_d(const std::vector<double>& a,
                                        const std::vector<double>& b,
                                        const std::vector<double>& c,
                                        std::vector<double> d) {
    const std::size_t n = d.size();
    if (n == 0) return {};
    if (n == 1) return {d[0] / b[0]};

    std::vector<double> c_prime(n, 0.0);
    std::vector<double> d_prime(n, 0.0);

    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (std::size_t i = 1; i < n; ++i) {
        double denom = b[i] - a[i] * c_prime[i - 1];
        if (std::abs(denom) < 1e-15) {
            denom = 1e-15;
        }
        if (i < n - 1) {
            c_prime[i] = c[i] / denom;
        }
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom;
    }

    std::vector<double> x(n, 0.0);
    x[n - 1] = d_prime[n - 1];
    for (std::size_t i = n - 1; i > 0; --i) {
        x[i - 1] = d_prime[i - 1] - c_prime[i - 1] * x[i];
    }
    return x;
}

} // namespace

// ============================================================================
// 1. 抛物型方程：热传导 / 扩散方程求解器实现
// ============================================================================

HeatSolver1D::HeatSolver1D(Scalar alpha, Scalar x_min, Scalar x_max, std::size_t nx)
    : alpha_(alpha), x_min_(x_min), x_max_(x_max), nx_(nx) {
    if (nx < 3) throw std::runtime_error("HeatSolver1D requires at least 3 spatial points");
    if (alpha <= Scalar(0)) throw std::runtime_error("HeatSolver1D: alpha must be positive");
    if (x_max <= x_min) throw std::runtime_error("HeatSolver1D: x_max must be greater than x_min");
}

Matrix HeatSolver1D::solve(std::function<Scalar(Scalar)> u0,
                           Scalar t_end,
                           std::size_t nt,
                           BoundaryCondition1D bc,
                           std::function<Scalar(Scalar, Scalar)> source) const {
    if (nt < 2) throw std::runtime_error("HeatSolver1D requires at least 2 time steps");
    if (t_end <= Scalar(0)) throw std::runtime_error("HeatSolver1D: t_end must be positive");

    const double x_min_d = static_cast<double>(x_min_);
    const double x_max_d = static_cast<double>(x_max_);
    const double t_end_d = static_cast<double>(t_end);
    const double alpha_d = static_cast<double>(alpha_);

    const Scalar dx_s = (x_max_ - x_min_) / Scalar(static_cast<long long>(nx_ - 1));
    const Scalar dt_s = t_end / Scalar(static_cast<long long>(nt - 1));

    const double dx = (x_max_d - x_min_d) / static_cast<double>(nx_ - 1);
    const double dt = t_end_d / static_cast<double>(nt - 1);
    const double r = (alpha_d * dt) / (dx * dx);
    const double r_half = r / 2.0;

    Matrix result(nt, nx_, Scalar(0));

    // 初始化 t = 0 的解
    std::vector<double> u_curr(nx_, 0.0);
    for (std::size_t i = 0; i < nx_; ++i) {
        Scalar x_s = x_min_ + Scalar(static_cast<long long>(i)) * dx_s;
        Scalar u_val = u0(x_s);
        u_curr[i] = static_cast<double>(u_val);
        result.at(0, i) = Scalar(u_curr[i]);
    }
    double left_0 = static_cast<double>(bc.val_left(Scalar(0)));
    double right_0 = static_cast<double>(bc.val_right(Scalar(0)));
    u_curr[0] = left_0;
    u_curr[nx_ - 1] = right_0;
    result.at(0, 0) = Scalar(left_0);
    result.at(0, nx_ - 1) = Scalar(right_0);

    const std::size_t inner_n = nx_ - 2;
    std::vector<double> a(inner_n, -r_half);
    std::vector<double> b(inner_n, 1.0 + r);
    std::vector<double> c(inner_n, -r_half);

    for (std::size_t t_step = 1; t_step < nt; ++t_step) {
        Scalar t_curr_s = Scalar(static_cast<long long>(t_step)) * dt_s;
        Scalar t_prev_s = t_curr_s - dt_s;

        double left_curr = static_cast<double>(bc.val_left(t_curr_s));
        double right_curr = static_cast<double>(bc.val_right(t_curr_s));

        std::vector<double> d(inner_n, 0.0);
        for (std::size_t i = 0; i < inner_n; ++i) {
            std::size_t orig_i = i + 1;
            double u_prev = u_curr[orig_i];
            double u_prev_l = u_curr[orig_i - 1];
            double u_prev_r = u_curr[orig_i + 1];

            d[i] = r_half * u_prev_l + (1.0 - r) * u_prev + r_half * u_prev_r;

            if (source) {
                Scalar x_s = x_min_ + Scalar(static_cast<long long>(orig_i)) * dx_s;
                double s_prev = static_cast<double>(source(x_s, t_prev_s));
                double s_curr = static_cast<double>(source(x_s, t_curr_s));
                d[i] += (dt / 2.0) * (s_prev + s_curr);
            }
        }

        d[0] += r_half * left_curr;
        d[inner_n - 1] += r_half * right_curr;

        std::vector<double> inner_sol = solve_tridiagonal_d(a, b, c, d);

        u_curr[0] = left_curr;
        result.at(t_step, 0) = Scalar(left_curr);
        for (std::size_t i = 0; i < inner_n; ++i) {
            u_curr[i + 1] = inner_sol[i];
            result.at(t_step, i + 1) = Scalar(inner_sol[i]);
        }
        u_curr[nx_ - 1] = right_curr;
        result.at(t_step, nx_ - 1) = Scalar(right_curr);
    }

    return result;
}

HeatSolver2D::HeatSolver2D(Scalar alpha, Scalar x_min, Scalar x_max, std::size_t nx,
                           Scalar y_min, Scalar y_max, std::size_t ny)
    : alpha_(alpha), x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max),
      nx_(nx), ny_(ny) {
    if (nx < 3 || ny < 3) throw std::runtime_error("HeatSolver2D requires at least 3x3 grid points");
    if (alpha <= Scalar(0)) throw std::runtime_error("HeatSolver2D: alpha must be positive");
    if (x_max <= x_min || y_max <= y_min) throw std::runtime_error("HeatSolver2D: spatial bounds invalid (max must be > min)");
}

Matrix HeatSolver2D::solve(std::function<Scalar(Scalar, Scalar)> u0,
                           Scalar t_end,
                           std::size_t nt,
                           BoundaryCondition2D bc) const {
    if (nt < 2) throw std::runtime_error("HeatSolver2D requires at least 2 time steps");
    if (t_end <= Scalar(0)) throw std::runtime_error("HeatSolver2D: t_end must be positive");
    const double x_min_d = static_cast<double>(x_min_);
    const double x_max_d = static_cast<double>(x_max_);
    const double y_min_d = static_cast<double>(y_min_);
    const double y_max_d = static_cast<double>(y_max_);
    const double alpha_d = static_cast<double>(alpha_);
    const double t_end_d = static_cast<double>(t_end);

    const Scalar dx_s = (x_max_ - x_min_) / Scalar(static_cast<long long>(nx_ - 1));
    const Scalar dy_s = (y_max_ - y_min_) / Scalar(static_cast<long long>(ny_ - 1));

    const double dx = (x_max_d - x_min_d) / static_cast<double>(nx_ - 1);
    const double dy = (y_max_d - y_min_d) / static_cast<double>(ny_ - 1);
    const double dt = t_end_d / static_cast<double>(nt - 1);

    const double rx = (alpha_d * (dt / 2.0)) / (dx * dx);
    const double ry = (alpha_d * (dt / 2.0)) / (dy * dy);

    std::vector<std::vector<double>> u(nx_, std::vector<double>(ny_, 0.0));

    // 初始化网格
    for (std::size_t i = 0; i < nx_; ++i) {
        Scalar x_s = x_min_ + Scalar(static_cast<long long>(i)) * dx_s;
        for (std::size_t j = 0; j < ny_; ++j) {
            Scalar y_s = y_min_ + Scalar(static_cast<long long>(j)) * dy_s;
            if (i == 0 || i == nx_ - 1 || j == 0 || j == ny_ - 1) {
                u[i][j] = bc.boundary_fn ? static_cast<double>(bc.boundary_fn(x_s, y_s)) : 0.0;
            } else {
                u[i][j] = static_cast<double>(u0(x_s, y_s));
            }
        }
    }

    // ADI 迭代步
    const std::size_t in_x = nx_ - 2;
    const std::size_t in_y = ny_ - 2;

    std::vector<double> ax(in_x, -rx);
    std::vector<double> bx(in_x, 1.0 + 2.0 * rx);
    std::vector<double> cx(in_x, -rx);

    std::vector<double> ay(in_y, -ry);
    std::vector<double> by(in_y, 1.0 + 2.0 * ry);
    std::vector<double> cy(in_y, -ry);

    std::vector<std::vector<double>> u_half = u;

    for (std::size_t t = 0; t < nt; ++t) {
        // 第一半步：x 方向隐式求解
        for (std::size_t j = 1; j <= in_y; ++j) {
            std::vector<double> d(in_x, 0.0);
            for (std::size_t i = 1; i <= in_x; ++i) {
                d[i - 1] = (1.0 - 2.0 * ry) * u[i][j] + ry * (u[i][j - 1] + u[i][j + 1]);
            }
            d[0] += rx * u[0][j];
            d[in_x - 1] += rx * u[nx_ - 1][j];

            std::vector<double> sol = solve_tridiagonal_d(ax, bx, cx, d);
            for (std::size_t i = 1; i <= in_x; ++i) {
                u_half[i][j] = sol[i - 1];
            }
        }

        // 第二半步：y 方向隐式求解
        for (std::size_t i = 1; i <= in_x; ++i) {
            std::vector<double> d(in_y, 0.0);
            for (std::size_t j = 1; j <= in_y; ++j) {
                d[j - 1] = (1.0 - 2.0 * rx) * u_half[i][j] + rx * (u_half[i - 1][j] + u_half[i + 1][j]);
            }
            d[0] += ry * u_half[i][0];
            d[in_y - 1] += ry * u_half[i][ny_ - 1];

            std::vector<double> sol = solve_tridiagonal_d(ay, by, cy, d);
            for (std::size_t j = 1; j <= in_y; ++j) {
                u[i][j] = sol[j - 1];
            }
        }
    }

    Matrix res(nx_, ny_, Scalar(0));
    for (std::size_t i = 0; i < nx_; ++i) {
        for (std::size_t j = 0; j < ny_; ++j) {
            res.at(i, j) = Scalar(u[i][j]);
        }
    }
    return res;
}

// ============================================================================
// 2. 双曲型方程：波动方程求解器实现
// ============================================================================

WaveSolver1D::WaveSolver1D(Scalar c, Scalar x_min, Scalar x_max, std::size_t nx)
    : c_(c), x_min_(x_min), x_max_(x_max), nx_(nx) {
    if (nx < 3) throw std::runtime_error("WaveSolver1D requires at least 3 spatial points");
    if (c <= Scalar(0)) throw std::runtime_error("WaveSolver1D: wave speed c must be positive");
    if (x_max <= x_min) throw std::runtime_error("WaveSolver1D: x_max must be greater than x_min");
}

Matrix WaveSolver1D::solve(std::function<Scalar(Scalar)> u0,
                           std::function<Scalar(Scalar)> v0,
                           Scalar t_end,
                           std::size_t nt,
                           BoundaryCondition1D bc) const {
    if (nt < 2) throw std::runtime_error("WaveSolver1D requires at least 2 time steps");
    if (t_end <= Scalar(0)) throw std::runtime_error("WaveSolver1D: t_end must be positive");
    const double x_min_d = static_cast<double>(x_min_);
    const double x_max_d = static_cast<double>(x_max_);
    const double c_d = static_cast<double>(c_);
    const double t_end_d = static_cast<double>(t_end);

    const Scalar dx_s = (x_max_ - x_min_) / Scalar(static_cast<long long>(nx_ - 1));
    const Scalar dt_s = t_end / Scalar(static_cast<long long>(nt - 1));

    const double dx = (x_max_d - x_min_d) / static_cast<double>(nx_ - 1);
    const double dt = t_end_d / static_cast<double>(nt - 1);
    const double r = (c_d * dt) / dx;
    const double r2 = r * r;

    Matrix result(nt, nx_, Scalar(0));
    std::vector<std::vector<double>> u(nt, std::vector<double>(nx_, 0.0));

    // t = 0
    for (std::size_t i = 0; i < nx_; ++i) {
        Scalar x_s = x_min_ + Scalar(static_cast<long long>(i)) * dx_s;
        u[0][i] = static_cast<double>(u0(x_s));
        result.at(0, i) = Scalar(u[0][i]);
    }
    double l0 = static_cast<double>(bc.val_left(Scalar(0)));
    double r0 = static_cast<double>(bc.val_right(Scalar(0)));
    u[0][0] = l0; u[0][nx_ - 1] = r0;
    result.at(0, 0) = Scalar(l0); result.at(0, nx_ - 1) = Scalar(r0);

    // t = 1
    double t1_d = dt;
    double l1 = static_cast<double>(bc.val_left(dt_s));
    double r1 = static_cast<double>(bc.val_right(dt_s));
    u[1][0] = l1; u[1][nx_ - 1] = r1;
    result.at(1, 0) = Scalar(l1); result.at(1, nx_ - 1) = Scalar(r1);

    for (std::size_t i = 1; i < nx_ - 1; ++i) {
        Scalar x_s = x_min_ + Scalar(static_cast<long long>(i)) * dx_s;
        double v = v0 ? static_cast<double>(v0(x_s)) : 0.0;
        u[1][i] = u[0][i] + t1_d * v + (r2 / 2.0) * (u[0][i + 1] - 2.0 * u[0][i] + u[0][i - 1]);
        result.at(1, i) = Scalar(u[1][i]);
    }

    // 显式 Leapfrog 迭代
    for (std::size_t t = 2; t < nt; ++t) {
        Scalar t_curr_s = Scalar(static_cast<long long>(t)) * dt_s;
        double lt = static_cast<double>(bc.val_left(t_curr_s));
        double rt = static_cast<double>(bc.val_right(t_curr_s));
        u[t][0] = lt; u[t][nx_ - 1] = rt;
        result.at(t, 0) = Scalar(lt); result.at(t, nx_ - 1) = Scalar(rt);

        for (std::size_t i = 1; i < nx_ - 1; ++i) {
            u[t][i] = 2.0 * u[t - 1][i] - u[t - 2][i] +
                      r2 * (u[t - 1][i + 1] - 2.0 * u[t - 1][i] + u[t - 1][i - 1]);
            result.at(t, i) = Scalar(u[t][i]);
        }
    }

    return result;
}

WaveSolver2D::WaveSolver2D(Scalar c, Scalar x_min, Scalar x_max, std::size_t nx,
                           Scalar y_min, Scalar y_max, std::size_t ny)
    : c_(c), x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max), nx_(nx), ny_(ny) {
    if (nx < 3 || ny < 3) throw std::runtime_error("WaveSolver2D requires at least 3x3 grid points");
    if (c <= Scalar(0)) throw std::runtime_error("WaveSolver2D: wave speed c must be positive");
    if (x_max <= x_min || y_max <= y_min) throw std::runtime_error("WaveSolver2D: spatial bounds invalid (max must be > min)");
}

Matrix WaveSolver2D::solve(std::function<Scalar(Scalar, Scalar)> u0,
                           std::function<Scalar(Scalar, Scalar)> v0,
                           Scalar t_end,
                           std::size_t nt,
                           BoundaryCondition2D bc) const {
    if (nt < 2) throw std::runtime_error("WaveSolver2D requires at least 2 time steps");
    if (t_end <= Scalar(0)) throw std::runtime_error("WaveSolver2D: t_end must be positive");
    const double x_min_d = static_cast<double>(x_min_);
    const double x_max_d = static_cast<double>(x_max_);
    const double y_min_d = static_cast<double>(y_min_);
    const double y_max_d = static_cast<double>(y_max_);
    const double c_d = static_cast<double>(c_);
    const double t_end_d = static_cast<double>(t_end);

    const Scalar dx_s = (x_max_ - x_min_) / Scalar(static_cast<long long>(nx_ - 1));
    const Scalar dy_s = (y_max_ - y_min_) / Scalar(static_cast<long long>(ny_ - 1));

    const double dx = (x_max_d - x_min_d) / static_cast<double>(nx_ - 1);
    const double dy = (y_max_d - y_min_d) / static_cast<double>(ny_ - 1);
    const double dt = t_end_d / static_cast<double>(nt - 1);

    const double rx2 = (c_d * dt / dx) * (c_d * dt / dx);
    const double ry2 = (c_d * dt / dy) * (c_d * dt / dy);

    std::vector<std::vector<double>> u_prev(nx_, std::vector<double>(ny_, 0.0));
    std::vector<std::vector<double>> u_curr(nx_, std::vector<double>(ny_, 0.0));
    std::vector<std::vector<double>> u_next(nx_, std::vector<double>(ny_, 0.0));

    for (std::size_t i = 0; i < nx_; ++i) {
        Scalar x_s = x_min_ + Scalar(static_cast<long long>(i)) * dx_s;
        for (std::size_t j = 0; j < ny_; ++j) {
            Scalar y_s = y_min_ + Scalar(static_cast<long long>(j)) * dy_s;
            if (i == 0 || i == nx_ - 1 || j == 0 || j == ny_ - 1) {
                u_prev[i][j] = bc.boundary_fn ? static_cast<double>(bc.boundary_fn(x_s, y_s)) : 0.0;
                u_curr[i][j] = u_prev[i][j];
            } else {
                u_prev[i][j] = static_cast<double>(u0(x_s, y_s));
            }
        }
    }

    for (std::size_t i = 1; i < nx_ - 1; ++i) {
        Scalar x_s = x_min_ + Scalar(static_cast<long long>(i)) * dx_s;
        for (std::size_t j = 1; j < ny_ - 1; ++j) {
            Scalar y_s = y_min_ + Scalar(static_cast<long long>(j)) * dy_s;
            double v = v0 ? static_cast<double>(v0(x_s, y_s)) : 0.0;
            double lap = rx2 * (u_prev[i + 1][j] - 2.0 * u_prev[i][j] + u_prev[i - 1][j]) +
                         ry2 * (u_prev[i][j + 1] - 2.0 * u_prev[i][j] + u_prev[i][j - 1]);
            u_curr[i][j] = u_prev[i][j] + dt * v + 0.5 * lap;
        }
    }

    for (std::size_t t = 1; t < nt; ++t) {
        for (std::size_t i = 1; i < nx_ - 1; ++i) {
            for (std::size_t j = 1; j < ny_ - 1; ++j) {
                double lap = rx2 * (u_curr[i + 1][j] - 2.0 * u_curr[i][j] + u_curr[i - 1][j]) +
                             ry2 * (u_curr[i][j + 1] - 2.0 * u_curr[i][j] + u_curr[i][j - 1]);
                u_next[i][j] = 2.0 * u_curr[i][j] - u_prev[i][j] + lap;
            }
        }
        u_prev = u_curr;
        u_curr = u_next;
    }

    Matrix res(nx_, ny_, Scalar(0));
    for (std::size_t i = 0; i < nx_; ++i) {
        for (std::size_t j = 0; j < ny_; ++j) {
            res.at(i, j) = Scalar(u_curr[i][j]);
        }
    }
    return res;
}

// ============================================================================
// 3. 椭圆型方程：泊松 / 拉普拉斯方程求解器实现 (SOR 迭代)
// ============================================================================

PoissonSolver2D::PoissonSolver2D(Scalar x_min, Scalar x_max, std::size_t nx,
                                 Scalar y_min, Scalar y_max, std::size_t ny)
    : x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max), nx_(nx), ny_(ny) {
    if (nx < 3 || ny < 3) throw std::runtime_error("PoissonSolver2D requires at least 3x3 grid");
    if (x_max <= x_min || y_max <= y_min) throw std::runtime_error("PoissonSolver2D: spatial bounds invalid (max must be > min)");
}

Matrix PoissonSolver2D::solve(std::function<Scalar(Scalar, Scalar)> source,
                             BoundaryCondition2D bc,
                             Scalar tol,
                             std::size_t max_iter) const {
    if (tol <= Scalar(0)) throw std::runtime_error("PoissonSolver2D: tol must be positive");
    if (max_iter == 0) throw std::runtime_error("PoissonSolver2D: max_iter must be greater than 0");
    const double x_min_d = static_cast<double>(x_min_);
    const double x_max_d = static_cast<double>(x_max_);
    const double y_min_d = static_cast<double>(y_min_);
    const double y_max_d = static_cast<double>(y_max_);
    const double tol_d = static_cast<double>(tol);

    const Scalar dx_s = (x_max_ - x_min_) / Scalar(static_cast<long long>(nx_ - 1));
    const Scalar dy_s = (y_max_ - y_min_) / Scalar(static_cast<long long>(ny_ - 1));

    const double dx = (x_max_d - x_min_d) / static_cast<double>(nx_ - 1);
    const double dy = (y_max_d - y_min_d) / static_cast<double>(ny_ - 1);
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const double factor = 2.0 * (1.0 / dx2 + 1.0 / dy2);

    const double pi = 3.14159265358979323846;
    const double rho = (std::cos(pi / static_cast<double>(nx_)) +
                        std::cos(pi / static_cast<double>(ny_))) / 2.0;
    const double rad = std::max(0.0, 1.0 - rho * rho);
    const double omega = std::clamp(2.0 / (1.0 + std::sqrt(rad)), 1.0, 1.999);

    std::vector<std::vector<double>> u(nx_, std::vector<double>(ny_, 0.0));

    // 初始化边界
    for (std::size_t i = 0; i < nx_; ++i) {
        Scalar x_s = x_min_ + Scalar(static_cast<long long>(i)) * dx_s;
        for (std::size_t j = 0; j < ny_; ++j) {
            Scalar y_s = y_min_ + Scalar(static_cast<long long>(j)) * dy_s;
            if (i == 0 || i == nx_ - 1 || j == 0 || j == ny_ - 1) {
                u[i][j] = bc.boundary_fn ? static_cast<double>(bc.boundary_fn(x_s, y_s)) : 0.0;
            }
        }
    }

    // SOR 迭代
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        double max_diff = 0.0;

        for (std::size_t i = 1; i < nx_ - 1; ++i) {
            Scalar x_s = x_min_ + Scalar(static_cast<long long>(i)) * dx_s;
            double dx_sq_inv = 1.0 / dx2;
            double dy_sq_inv = 1.0 / dy2;

            for (std::size_t j = 1; j < ny_ - 1; ++j) {
                Scalar y_s = y_min_ + Scalar(static_cast<long long>(j)) * dy_s;
                double f_val = source ? static_cast<double>(source(x_s, y_s)) : 0.0;

                double u_star = ((u[i + 1][j] + u[i - 1][j]) * dx_sq_inv +
                                 (u[i][j + 1] + u[i][j - 1]) * dy_sq_inv - f_val) / factor;

                double new_val = (1.0 - omega) * u[i][j] + omega * u_star;
                double diff = std::abs(new_val - u[i][j]);
                if (diff > max_diff) max_diff = diff;

                u[i][j] = new_val;
            }
        }

        if (max_diff < tol_d) {
            break;
        }
    }

    Matrix res(nx_, ny_, Scalar(0));
    for (std::size_t i = 0; i < nx_; ++i) {
        for (std::size_t j = 0; j < ny_; ++j) {
            res.at(i, j) = Scalar(u[i][j]);
        }
    }
    return res;
}

} // namespace analysis::pde

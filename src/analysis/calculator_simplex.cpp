// ============================================================================
// 单纯形法求解线性规划问题
// ============================================================================

#include "calculator_simplex.h"

#include "mymath.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace simplex {

namespace {

// 使用 Sherman-Morrison 公式更新基矩阵之逆 B_inv
// B_new = B + (a_entering - B_leaving) * e_leaving'
// B_inv_new = B_inv - (B_inv * (a_entering - B_leaving) * e_leaving' * B_inv) / (1 + e_leaving' * B_inv * (a_entering - B_leaving))
// 简化版：令 u = a_entering, v = e_leaving, 则 B_inv_new = B_inv - (B_inv * u - e_leaving) * (row_leaving_of_B_inv) / (e_leaving' * B_inv * u)
bool update_basis_inverse(
    matrix::Matrix& B_inv,
    const std::vector<double>& a_enter,
    std::size_t leaving_idx,
    double eps) {

    const std::size_t m = B_inv.rows;
    
    // 计算 w = B_inv * a_enter
    std::vector<double> w(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            w[i] += B_inv.at(i, j) * a_enter[j];
        }
    }

    const double pivot = w[leaving_idx];
    if (std::abs(pivot) < eps * 1e-3) {
        return false; // 数值不稳定
    }

    // 更新 B_inv
    for (std::size_t i = 0; i < m; ++i) {
        if (i == leaving_idx) continue;
        const double factor = w[i] / pivot;
        for (std::size_t j = 0; j < m; ++j) {
            B_inv.at(i, j) -= factor * B_inv.at(leaving_idx, j);
        }
    }
    
    // 处理出基行
    for (std::size_t j = 0; j < m; ++j) {
        B_inv.at(leaving_idx, j) /= pivot;
    }

    return true;
}

// 重新从头计算基矩阵之逆
bool reinvert_basis(
    matrix::Matrix& B_inv,
    const std::vector<std::size_t>& basis_curr,
    const std::vector<std::vector<double>>& A_full,
    std::size_t m) {

    matrix::Matrix B(m, m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t k = 0; k < m; ++k) {
            B.at(i, k) = A_full[i][basis_curr[k]];
        }
    }
    try {
        B_inv = matrix::inverse(B);
        return true;
    } catch (...) {
        return false;
    }
}

// 单纯形迭代核心算法 (优化版: O(m^2) 迭代)
bool simplex_iterate(
    std::vector<double>& x_curr,
    std::vector<std::size_t>& basis_curr,
    const std::vector<double>& c_obj,
    const std::vector<std::vector<double>>& A_full,
    const std::vector<double>& lb_full,
    const std::vector<double>& ub_full,
    std::size_t m_total,
    std::size_t n_full,
    bool minimize,
    std::size_t max_iters,
    double eps) {

    matrix::Matrix B_inv(m_total, m_total, 0.0);
    if (!reinvert_basis(B_inv, basis_curr, A_full, m_total)) {
        return false;
    }

    for (std::size_t iter = 0; iter < max_iters; ++iter) {
        // 每 50 次迭代重新求逆以保持数值稳定性
        if (iter > 0 && iter % 50 == 0) {
            if (!reinvert_basis(B_inv, basis_curr, A_full, m_total)) return false;
        }

        // 计算对偶变量 y: y' = c_B' * B_inv
        std::vector<double> y(m_total, 0.0);
        for (std::size_t j = 0; j < m_total; ++j) {
            const double cB_val = c_obj[basis_curr[j]];
            for (std::size_t i = 0; i < m_total; ++i) {
                y[i] += cB_val * B_inv.at(j, i);
            }
        }

        // 找入基变量（定价）
        std::size_t entering = n_full;
        double best_rc = 0.0;
        for (std::size_t j = 0; j < n_full; ++j) {
            bool is_basic = false;
            for (std::size_t k = 0; k < m_total; ++k) {
                if (basis_curr[k] == j) { is_basic = true; break; }
            }
            if (is_basic) continue;

            double rc = c_obj[j];
            for (std::size_t i = 0; i < m_total; ++i) {
                rc -= y[i] * A_full[i][j];
            }

            double xj = x_curr[j];
            bool at_lower = std::abs(xj - lb_full[j]) <= eps;
            bool at_upper = (ub_full[j] < std::numeric_limits<double>::infinity()) &&
                            std::abs(xj - ub_full[j]) <= eps;

            if (minimize) {
                if (at_lower && rc < -eps) {
                    if (rc < best_rc) { best_rc = rc; entering = j; }
                } else if (at_upper && rc > eps) {
                    if (-rc < best_rc) { best_rc = -rc; entering = j; }
                }
            } else {
                if (at_lower && rc > eps) {
                    if (-rc < best_rc) { best_rc = -rc; entering = j; }
                } else if (at_upper && rc < -eps) {
                    if (rc < best_rc) { best_rc = rc; entering = j; }
                }
            }
        }

        if (entering >= n_full || std::abs(best_rc) <= eps) {
            return true; // 找到最优解
        }

        // 计算搜索方向 d = B_inv * A_entering
        std::vector<double> a_enter(m_total);
        for (std::size_t i = 0; i < m_total; ++i) a_enter[i] = A_full[i][entering];

        std::vector<double> d(m_total, 0.0);
        for (std::size_t i = 0; i < m_total; ++i) {
            for (std::size_t j = 0; j < m_total; ++j) {
                d[i] += B_inv.at(i, j) * a_enter[j];
            }
        }

        // 确定步长 (Ratio Test)
        double theta = std::numeric_limits<double>::infinity();
        std::size_t leaving = m_total;

        // 如果变量从上界减少
        bool decreasing_entering = false;
        if (ub_full[entering] < std::numeric_limits<double>::infinity() &&
            x_curr[entering] >= ub_full[entering] - eps) {
            decreasing_entering = true;
        }

        double max_theta = (ub_full[entering] < std::numeric_limits<double>::infinity())
                            ? (ub_full[entering] - lb_full[entering])
                            : std::numeric_limits<double>::infinity();

        theta = max_theta;

        for (std::size_t i = 0; i < m_total; ++i) {
            std::size_t j = basis_curr[i];
            double di = decreasing_entering ? -d[i] : d[i];
            if (std::abs(di) <= eps) continue;

            if (di > 0) {
                double ratio = (x_curr[j] - lb_full[j]) / di;
                if (ratio < theta) { theta = ratio; leaving = i; }
            } else {
                if (ub_full[j] < std::numeric_limits<double>::infinity()) {
                    double ratio = (ub_full[j] - x_curr[j]) / (-di);
                    if (ratio < theta) { theta = ratio; leaving = i; }
                }
            }
        }

        if (theta >= std::numeric_limits<double>::infinity() - 1e9) {
            return true; // 无界或已达到极值
        }

        // 更新当前解
        double shift = decreasing_entering ? -theta : theta;
        x_curr[entering] += shift;
        for (std::size_t i = 0; i < m_total; ++i) {
            x_curr[basis_curr[i]] -= shift * d[i];
        }

        // 换基与 B_inv 更新
        if (leaving < m_total) {
            if (!update_basis_inverse(B_inv, a_enter, leaving, eps)) {
                // 如果秩 1 更新失败，尝试重新从头计算
                if (!reinvert_basis(B_inv, basis_curr, A_full, m_total)) return false;
            }
            basis_curr[leaving] = entering;
        }
    }

    return true;
}

}  // namespace

bool solve_linear_box_problem(
    const std::vector<double>& objective,
    const matrix::Matrix& inequality_matrix,
    const std::vector<double>& inequality_rhs,
    const matrix::Matrix& equality_matrix,
    const std::vector<double>& equality_rhs,
    const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds,
    double tolerance,
    std::vector<double>* solution,
    double* objective_value,
    std::string* diagnostic) {

    const std::size_t n = objective.size();
    if (inequality_matrix.cols != n ||
        inequality_rhs.size() != inequality_matrix.rows ||
        equality_matrix.cols != n ||
        equality_rhs.size() != equality_matrix.rows ||
        lower_bounds.size() != n ||
        upper_bounds.size() != n) {
        throw std::runtime_error("planning dimension mismatch");
    }

    if (n == 0) {
        *solution = {};
        *objective_value = 0.0;
        return true;
    }

    const std::size_t m_ineq = inequality_matrix.rows;
    const std::size_t m_eq = equality_matrix.rows;
    const std::size_t m_total = m_ineq + m_eq;

    if (m_total == 0) {
        std::vector<double> best(n, 0.0);
        double best_val = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (objective[j] >= 0) {
                best[j] = lower_bounds[j];
            } else {
                best[j] = upper_bounds[j];
            }
            best_val += objective[j] * best[j];
        }
        *solution = best;
        *objective_value = best_val;
        return true;
    }

    const std::size_t n_slack = m_ineq;
    const std::size_t n_total = n + n_slack;

    std::vector<double> c(n_total, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
        c[j] = objective[j];
    }

    std::vector<double> lb(n_total, 0.0);
    std::vector<double> ub(n_total, std::numeric_limits<double>::infinity());
    for (std::size_t j = 0; j < n; ++j) {
        lb[j] = lower_bounds[j];
        ub[j] = upper_bounds[j];
    }

    std::vector<std::vector<double>> A(m_total, std::vector<double>(n_total, 0.0));
    std::vector<double> b(m_total, 0.0);

    for (std::size_t i = 0; i < m_ineq; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            A[i][j] = inequality_matrix.at(i, j);
        }
        A[i][n + i] = 1.0;
        b[i] = inequality_rhs[i];
    }

    for (std::size_t i = 0; i < m_eq; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            A[m_ineq + i][j] = equality_matrix.at(i, j);
        }
        b[m_ineq + i] = equality_rhs[i];
    }

    for (std::size_t i = 0; i < m_total; ++i) {
        if (b[i] < 0) {
            b[i] = -b[i];
            for (std::size_t j = 0; j < n_total; ++j) {
                A[i][j] = -A[i][j];
            }
        }
    }

    const std::size_t n_art = m_total;
    const std::size_t n_full = n_total + n_art;

    std::vector<std::vector<double>> A_full(m_total, std::vector<double>(n_full, 0.0));
    for (std::size_t i = 0; i < m_total; ++i) {
        for (std::size_t j = 0; j < n_total; ++j) {
            A_full[i][j] = A[i][j];
        }
        A_full[i][n_total + i] = 1.0;
    }

    std::vector<double> lb_full(n_full, 0.0);
    std::vector<double> ub_full(n_full, std::numeric_limits<double>::infinity());
    for (std::size_t j = 0; j < n_total; ++j) {
        lb_full[j] = lb[j];
        ub_full[j] = ub[j];
    }

    std::vector<double> c_phase1(n_full, 0.0);
    for (std::size_t j = n_total; j < n_full; ++j) {
        c_phase1[j] = 1.0;
    }

    std::vector<std::size_t> basis(m_total);
    for (std::size_t i = 0; i < m_total; ++i) {
        basis[i] = n_total + i;
    }

    std::vector<double> x(n_full, 0.0);
    for (std::size_t j = 0; j < n_total; ++j) {
        x[j] = lb_full[j];
    }
    for (std::size_t i = 0; i < m_total; ++i) {
        x[n_total + i] = b[i];
    }

    constexpr std::size_t kMaxSimplexIters = 10000;
    bool phase1_ok = simplex_iterate(x, basis, c_phase1, A_full, lb_full, ub_full,
                                      m_total, n_full, true, kMaxSimplexIters, tolerance);

    if (!phase1_ok) {
        if (diagnostic) *diagnostic = "phase 1 simplex failed";
        return false;
    }

    double art_sum = 0.0;
    for (std::size_t j = n_total; j < n_full; ++j) {
        art_sum += std::abs(x[j]);
    }
    if (art_sum > tolerance * static_cast<double>(m_total)) {
        if (diagnostic) *diagnostic = "no feasible solution found";
        return false;
    }

    for (std::size_t i = 0; i < m_total; ++i) {
        if (basis[i] >= n_total) {
            for (std::size_t j = 0; j < n_total; ++j) {
                bool in_basis = false;
                for (std::size_t k = 0; k < m_total; ++k) {
                    if (basis[k] == j) { in_basis = true; break; }
                }
                if (!in_basis) {
                    bool can_enter = false;
                    for (std::size_t row = 0; row < m_total; ++row) {
                        if (std::abs(A_full[row][j]) > tolerance) {
                            can_enter = true;
                            break;
                        }
                    }
                    if (can_enter) {
                        basis[i] = j;
                        break;
                    }
                }
            }
        }
    }

    std::vector<double> c_phase2(n_full, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
        c_phase2[j] = -objective[j];
    }

    bool phase2_ok = simplex_iterate(x, basis, c_phase2, A_full, lb_full, ub_full,
                                      m_total, n_full, true, kMaxSimplexIters, tolerance);

    if (!phase2_ok) {
        if (diagnostic) *diagnostic = "phase 2 simplex failed";
        return false;
    }

    std::vector<double> result(n, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
        result[j] = x[j];
        if (result[j] < lower_bounds[j] - tolerance ||
            result[j] > upper_bounds[j] + tolerance) {
            result[j] = std::max(lower_bounds[j], std::min(upper_bounds[j], result[j]));
        }
    }

    for (std::size_t i = 0; i < m_ineq; ++i) {
        double lhs = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            lhs += inequality_matrix.at(i, j) * result[j];
        }
        if (lhs > inequality_rhs[i] + tolerance) {
            if (diagnostic) *diagnostic = "solution violates inequality constraint";
            return false;
        }
    }
    for (std::size_t i = 0; i < m_eq; ++i) {
        double lhs = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            lhs += equality_matrix.at(i, j) * result[j];
        }
        if (std::abs(lhs - equality_rhs[i]) > tolerance) {
            if (diagnostic) *diagnostic = "solution violates equality constraint";
            return false;
        }
    }

    *solution = result;
    *objective_value = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        *objective_value += objective[j] * result[j];
    }
    return true;
}

}  // namespace simplex

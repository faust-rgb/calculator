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

// 单纯形迭代核心算法
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

    for (std::size_t iter = 0; iter < max_iters; ++iter) {
        // 计算对偶变量 y: B' * y = c_B
        std::vector<double> y(m_total, 0.0);
        {
            matrix::Matrix B(m_total, m_total, 0.0);
            matrix::Matrix cB(m_total, 1, 0.0);
            for (std::size_t i = 0; i < m_total; ++i) {
                cB.at(i, 0) = c_obj[basis_curr[i]];
                for (std::size_t k = 0; k < m_total; ++k) {
                    B.at(i, k) = A_full[i][basis_curr[k]];
                }
            }
            matrix::Matrix BT(m_total, m_total, 0.0);
            for (std::size_t i = 0; i < m_total; ++i) {
                for (std::size_t k = 0; k < m_total; ++k) {
                    BT.at(i, k) = B.at(k, i);
                }
            }
            try {
                matrix::Matrix y_mat = matrix::solve(BT, cB);
                for (std::size_t i = 0; i < m_total; ++i) {
                    y[i] = y_mat.at(i, 0);
                }
            } catch (...) {
                return false;
            }
        }

        // 找入基变量（最负检验数）
        std::size_t entering = n_full;
        double best_rc = 0.0;
        for (std::size_t j = 0; j < n_full; ++j) {
            bool is_basic = false;
            for (std::size_t k = 0; k < m_total; ++k) {
                if (basis_curr[k] == j) {
                    is_basic = true;
                    break;
                }
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
                if (at_lower && rc < best_rc) {
                    best_rc = rc;
                    entering = j;
                } else if (at_upper && rc > -best_rc) {
                    best_rc = -rc;
                    entering = j;
                } else if (!at_lower && !at_upper && rc < best_rc) {
                    best_rc = rc;
                    entering = j;
                }
            } else {
                if (at_lower && rc > -best_rc) {
                    best_rc = -rc;
                    entering = j;
                } else if (at_upper && rc < best_rc) {
                    best_rc = rc;
                    entering = j;
                } else if (!at_lower && !at_upper && rc > -best_rc) {
                    best_rc = -rc;
                    entering = j;
                }
            }
        }

        if (entering >= n_full || std::abs(best_rc) <= eps) {
            return true;
        }

        // 计算方向: B * d = A_entering
        std::vector<double> a_enter(m_total, 0.0);
        for (std::size_t i = 0; i < m_total; ++i) {
            a_enter[i] = A_full[i][entering];
        }

        std::vector<double> d(m_total, 0.0);
        {
            matrix::Matrix B(m_total, m_total, 0.0);
            for (std::size_t i = 0; i < m_total; ++i) {
                for (std::size_t k = 0; k < m_total; ++k) {
                    B.at(i, k) = A_full[i][basis_curr[k]];
                }
            }
            matrix::Matrix a_mat(m_total, 1, 0.0);
            for (std::size_t i = 0; i < m_total; ++i) {
                a_mat.at(i, 0) = a_enter[i];
            }
            try {
                matrix::Matrix d_mat = matrix::solve(B, a_mat);
                for (std::size_t i = 0; i < m_total; ++i) {
                    d[i] = d_mat.at(i, 0);
                }
            } catch (...) {
                return false;
            }
        }

        // 比值检验
        double theta = std::numeric_limits<double>::infinity();
        std::size_t leaving = n_full;
        bool from_upper = false;

        if (x_curr[entering] <= lb_full[entering] + eps) {
            theta = ub_full[entering] - lb_full[entering];
        } else if (ub_full[entering] < std::numeric_limits<double>::infinity() &&
                   x_curr[entering] >= ub_full[entering] - eps) {
            theta = x_curr[entering] - lb_full[entering];
            from_upper = true;
        }

        for (std::size_t i = 0; i < m_total; ++i) {
            std::size_t j = basis_curr[i];
            if (std::abs(d[i]) <= eps) continue;

            double ratio;
            if (d[i] > 0) {
                ratio = (x_curr[j] - lb_full[j]) / d[i];
                if (ratio < theta) {
                    theta = ratio;
                    leaving = i;
                    from_upper = false;
                }
            } else {
                if (ub_full[j] < std::numeric_limits<double>::infinity()) {
                    ratio = (ub_full[j] - x_curr[j]) / (-d[i]);
                    if (ratio < theta) {
                        theta = ratio;
                        leaving = i;
                        from_upper = true;
                    }
                }
            }
        }

        if (theta >= std::numeric_limits<double>::infinity() - 1) {
            return true;
        }

        x_curr[entering] += theta;
        for (std::size_t i = 0; i < m_total; ++i) {
            x_curr[basis_curr[i]] -= theta * d[i];
        }

        if (leaving < m_total) {
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

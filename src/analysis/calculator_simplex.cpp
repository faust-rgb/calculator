// ============================================================================
// 单纯形法求解线性规划问题
// ============================================================================
//
// 本文件实现了两阶段修正单纯形法，用于求解线性规划问题。
// 支持有界变量和等式/不等式约束。
//
// 算法特点：
// - 两阶段法：第一阶段找可行解，第二阶段优化
// - 修正单纯形：使用矩阵求逆更新，O(m²) 每次迭代
// - 有界变量处理：支持上下界约束

#include "calculator_simplex.h"

#include "mymath.h"

#include <algorithm>
#include <stdexcept>

namespace simplex {

namespace {

/**
 * @brief 使用 Sherman-Morrison 公式更新基矩阵之逆
 *
 * 当换基时，避免重新计算逆矩阵，而是使用秩 1 更新。
 * B_new = B + (a_entering - B_leaving) * e_leaving'
 * B_inv_new = B_inv - (B_inv * (a_entering - B_leaving) * e_leaving' * B_inv) / (1 + e_leaving' * B_inv * (a_entering - B_leaving))
 *
 * @param B_inv 当前基矩阵之逆（会被更新）
 * @param a_enter 入基变量的列向量
 * @param leaving_idx 出基变量在基中的索引
 * @param eps 数值容差
 * @return 是否成功更新
 */
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
    if (mymath::abs(pivot) < eps * 1e-3) {
        return false; // 数值不稳定
    }

    // 更新 B_inv 的非出基行
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

/**
 * @brief 重新从头计算基矩阵之逆
 *
 * 当秩 1 更新数值不稳定时，重新计算逆矩阵。
 */
bool reinvert_basis(
    matrix::Matrix& B_inv,
    const std::vector<std::size_t>& basis_curr,
    const std::vector<std::vector<double>>& A_full,
    std::size_t m) {

    // 构建基矩阵
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

/**
 * @brief 单纯形迭代核心算法
 *
 * 执行修正单纯形法的迭代过程：
 * 1. 计算对偶变量 y = c_B' * B_inv
 * 2. 选择入基变量（定价）
 * 3. 计算搜索方向 d = B_inv * a_enter
 * 4. 确定步长（比值检验）
 * 5. 更新解和基
 *
 * @param x_curr 当前解（会被更新）
 * @param basis_curr 当前基（会被更新）
 * @param c_obj 目标函数系数
 * @param A_full 约束矩阵（含松弛变量）
 * @param lb_full 变量下界
 * @param ub_full 变量上界
 * @param m_total 约束数量
 * @param n_full 变量数量（含松弛和人工变量）
 * @param minimize 是否为最小化
 * @param max_iters 最大迭代次数
 * @param eps 数值容差
 * @return 是否找到最优解
 */
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

    // 初始化基矩阵之逆
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

        // 寻找入基变量（定价阶段）
        std::size_t entering = n_full;
        double best_rc = 0.0;
        for (std::size_t j = 0; j < n_full; ++j) {
            // 检查是否已在基中
            bool is_basic = false;
            for (std::size_t k = 0; k < m_total; ++k) {
                if (basis_curr[k] == j) { is_basic = true; break; }
            }
            if (is_basic) continue;

            // 计算检验数（reduced cost）
            double rc = c_obj[j];
            for (std::size_t i = 0; i < m_total; ++i) {
                rc -= y[i] * A_full[i][j];
            }

            // 根据变量当前状态选择入基
            double xj = x_curr[j];
            bool at_lower = mymath::abs(xj - lb_full[j]) <= eps;
            bool at_upper = (ub_full[j] < mymath::infinity()) &&
                            mymath::abs(xj - ub_full[j]) <= eps;

            if (minimize) {
                // 最小化：检验数为负时可入基
                if (at_lower && rc < -eps) {
                    if (rc < best_rc) { best_rc = rc; entering = j; }
                } else if (at_upper && rc > eps) {
                    if (-rc < best_rc) { best_rc = -rc; entering = j; }
                }
            } else {
                // 最大化：检验数为正时可入基
                if (at_lower && rc > eps) {
                    if (-rc < best_rc) { best_rc = -rc; entering = j; }
                } else if (at_upper && rc < -eps) {
                    if (rc < best_rc) { best_rc = rc; entering = j; }
                }
            }
        }

        // 检查最优性
        if (entering >= n_full || mymath::abs(best_rc) <= eps) {
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

        // 比值检验确定步长
        double theta = mymath::infinity();
        std::size_t leaving = m_total;

        // 检查入基变量是否从上界减少
        bool decreasing_entering = false;
        if (ub_full[entering] < mymath::infinity() &&
            x_curr[entering] >= ub_full[entering] - eps) {
            decreasing_entering = true;
        }

        // 计算最大步长
        double max_theta = (ub_full[entering] < mymath::infinity())
                            ? (ub_full[entering] - lb_full[entering])
                            : mymath::infinity();

        theta = max_theta;

        // 对每个基变量进行比值检验
        for (std::size_t i = 0; i < m_total; ++i) {
            std::size_t j = basis_curr[i];
            double di = decreasing_entering ? -d[i] : d[i];
            if (mymath::abs(di) <= eps) continue;

            if (di > 0) {
                // 基变量可能到达下界
                double ratio = (x_curr[j] - lb_full[j]) / di;
                if (ratio < theta) { theta = ratio; leaving = i; }
            } else {
                // 基变量可能到达上界
                if (ub_full[j] < mymath::infinity()) {
                    double ratio = (ub_full[j] - x_curr[j]) / (-di);
                    if (ratio < theta) { theta = ratio; leaving = i; }
                }
            }
        }

        if (theta >= mymath::infinity() - 1e9) {
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
                // 秩 1 更新失败，尝试重新计算
                if (!reinvert_basis(B_inv, basis_curr, A_full, m_total)) return false;
            }
            basis_curr[leaving] = entering;
        }
    }

    return true;
}

}  // namespace

/**
 * @brief 使用两阶段单纯形法求解线性规划问题
 *
 * 求解问题:
 *   maximize c'x
 *   subject to: A_ineq * x <= b_ineq
 *               A_eq * x = b_eq
 *               lower <= x <= upper
 *
 * 第一阶段：添加人工变量，最小化人工变量之和，找到基本可行解
 * 第二阶段：删除人工变量，优化原目标函数
 */
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
    std::vector<double> ub(n_total, mymath::infinity());
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
    std::vector<double> ub_full(n_full, mymath::infinity());
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
        art_sum += mymath::abs(x[j]);
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
                        if (mymath::abs(A_full[row][j]) > tolerance) {
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
        if (mymath::abs(lhs - equality_rhs[i]) > tolerance) {
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

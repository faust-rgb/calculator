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

#include "analysis/optimization/simplex_engine.h"

#include "math/base/precision_constants.h"
#include "app/scalar_type.h"
#include "math/mymath.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace simplex {

namespace {

using Scalar = mymath::Scalar;

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
    const std::vector<Scalar>& a_enter,
    std::size_t leaving_idx,
    Scalar eps) {

    const std::size_t m = B_inv.rows;

    // 计算 w = B_inv * a_enter
    std::vector<Scalar> w(m, Scalar(0));
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            w[i] += Scalar(B_inv.at(i, j)) * a_enter[j];
        }
    }

    const Scalar pivot = w[leaving_idx];
    // 使用精度感知的阈值检测数值不稳定
    const Scalar pivot_threshold = eps * precision::sqrt_epsilon<Scalar>();
    if (mymath::abs(pivot) < pivot_threshold) {
        return false; // 数值不稳定
    }

    // 更新 B_inv 的非出基行
    for (std::size_t i = 0; i < m; ++i) {
        if (i == leaving_idx) continue;
        const Scalar factor = w[i] / pivot;
        for (std::size_t j = 0; j < m; ++j) {
            B_inv.at(i, j) = (Scalar(B_inv.at(i, j)) - factor * Scalar(B_inv.at(leaving_idx, j)));
        }
    }

    // 处理出基行
    for (std::size_t j = 0; j < m; ++j) {
        B_inv.at(leaving_idx, j) = (Scalar(B_inv.at(leaving_idx, j)) / pivot);
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
    const std::vector<std::vector<Scalar>>& A_full,
    std::size_t m) {

    // 构建基矩阵
    matrix::Matrix B(m, m, 0.0L);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t k = 0; k < m; ++k) {
            B.at(i, k) = (A_full[i][basis_curr[k]]);
        }
    }
    try {
        B_inv = matrix::inverse(B);
        return true;
    } catch (...) {
        return false;
    }
}

bool solve_square_system(std::vector<std::vector<Scalar>> matrix,
                         std::vector<Scalar> rhs,
                         Scalar eps,
                         std::vector<Scalar>* solution) {
    const std::size_t n = rhs.size();
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        Scalar best_abs = mymath::abs(matrix[col][col]);
        for (std::size_t row = col + 1; row < n; ++row) {
            const Scalar candidate_abs = mymath::abs(matrix[row][col]);
            if (candidate_abs > best_abs) {
                best_abs = candidate_abs;
                pivot = row;
            }
        }

        if (best_abs <= eps) return false;

        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(rhs[pivot], rhs[col]);
        }

        const Scalar pivot_value = matrix[col][col];
        for (std::size_t j = col; j < n; ++j) matrix[col][j] /= pivot_value;
        rhs[col] /= pivot_value;

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;
            const Scalar factor = matrix[row][col];
            if (mymath::abs(factor) <= eps) continue;
            for (std::size_t j = col; j < n; ++j) {
                matrix[row][j] -= factor * matrix[col][j];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    *solution = rhs;
    return true;
}

struct LinearConstraint {
    std::vector<Scalar> coeffs;
    Scalar rhs;
};

bool satisfies_equalities(const std::vector<Scalar>& point,
                          const std::vector<LinearConstraint>& equalities,
                          Scalar eps) {
    for (const LinearConstraint& constraint : equalities) {
        Scalar lhs = Scalar(0);
        for (std::size_t j = 0; j < point.size(); ++j) {
            lhs += constraint.coeffs[j] * point[j];
        }
        if (mymath::abs(lhs - constraint.rhs) > eps) return false;
    }
    return true;
}

bool satisfies_inequalities(const std::vector<Scalar>& point,
                            const std::vector<LinearConstraint>& inequalities,
                            Scalar eps) {
    for (const LinearConstraint& constraint : inequalities) {
        Scalar lhs = Scalar(0);
        for (std::size_t j = 0; j < point.size(); ++j) {
            lhs += constraint.coeffs[j] * point[j];
        }
        if (lhs > constraint.rhs + eps) return false;
    }
    return true;
}

bool solve_by_vertex_enumeration(
    const std::vector<Scalar>& objective,
    const matrix::Matrix& inequality_matrix,
    const std::vector<Scalar>& inequality_rhs,
    const matrix::Matrix& equality_matrix,
    const std::vector<Scalar>& equality_rhs,
    const std::vector<Scalar>& lower_bounds,
    const std::vector<Scalar>& upper_bounds,
    Scalar tolerance,
    std::vector<Scalar>* solution,
    Scalar* objective_value,
    std::string* diagnostic) {

    const std::size_t n = objective.size();
    const Scalar inf_val = Scalar(mymath::infinity());
    // 使用精度感知的内部容差
    const Scalar eps = std::max(tolerance * Scalar(100), precision::default_absolute_tolerance<Scalar>());

    std::vector<LinearConstraint> equalities;
    equalities.reserve(equality_matrix.rows);
    for (std::size_t i = 0; i < equality_matrix.rows; ++i) {
        LinearConstraint constraint{{}, Scalar(equality_rhs[i])};
        constraint.coeffs.resize(n, Scalar(0));
        for (std::size_t j = 0; j < n; ++j) constraint.coeffs[j] = Scalar(equality_matrix.at(i, j));
        equalities.push_back(constraint);
    }

    std::vector<LinearConstraint> inequalities;
    inequalities.reserve(inequality_matrix.rows + n * 2);
    for (std::size_t i = 0; i < inequality_matrix.rows; ++i) {
        LinearConstraint constraint{{}, Scalar(inequality_rhs[i])};
        constraint.coeffs.resize(n, Scalar(0));
        for (std::size_t j = 0; j < n; ++j) constraint.coeffs[j] = Scalar(inequality_matrix.at(i, j));
        inequalities.push_back(constraint);
    }

    for (std::size_t j = 0; j < n; ++j) {
        if (tolerance < 0) return false;
        if (lower_bounds[j] > upper_bounds[j] + eps) {
            if (diagnostic) *diagnostic = "inconsistent variable bounds";
            return false;
        }

        LinearConstraint lower{{}, -Scalar(lower_bounds[j])};
        lower.coeffs.resize(n, Scalar(0));
        lower.coeffs[j] = Scalar(-1);
        inequalities.push_back(lower);

        if (upper_bounds[j] < inf_val) {
            LinearConstraint upper{{}, Scalar(upper_bounds[j])};
            upper.coeffs.resize(n, Scalar(0));
            upper.coeffs[j] = Scalar(1);
            inequalities.push_back(upper);
        } else if (objective[j] > eps) {
            bool bounded_by_constraint = false;
            for (const LinearConstraint& constraint : inequalities) {
                if (constraint.coeffs[j] > eps) {
                    bounded_by_constraint = true;
                    break;
                }
            }
            if (!bounded_by_constraint && equality_matrix.rows == 0) {
                if (diagnostic) *diagnostic = "linear objective is unbounded";
                return false;
            }
        }
    }

    bool have_best = false;
    std::vector<Scalar> best_solution(n, Scalar(0));
    Scalar best_value = -inf_val;

    auto consider = [&](const std::vector<Scalar>& candidate) {
        if (!satisfies_equalities(candidate, equalities, eps) ||
            !satisfies_inequalities(candidate, inequalities, eps)) {
            return;
        }
        Scalar value = Scalar(0);
        for (std::size_t j = 0; j < n; ++j) value += objective[j] * candidate[j];
        if (!have_best || value > best_value + eps) {
            have_best = true;
            best_value = value;
            best_solution = candidate;
        }
    };

    auto solve_active_set = [&](const std::vector<std::size_t>& active) {
        std::vector<std::vector<Scalar>> system;
        std::vector<Scalar> rhs;
        system.reserve(equalities.size() + active.size());
        rhs.reserve(equalities.size() + active.size());

        for (const LinearConstraint& equality : equalities) {
            system.push_back(equality.coeffs);
            rhs.push_back(equality.rhs);
        }
        for (std::size_t index : active) {
            system.push_back(inequalities[index].coeffs);
            rhs.push_back(inequalities[index].rhs);
        }
        if (system.size() != n) return;

        std::vector<Scalar> candidate;
        if (solve_square_system(system, rhs, eps, &candidate)) {
            consider(candidate);
        }
    };

    if (n == 0) {
        *solution = {};
        *objective_value = Scalar(0);
        return true;
    }

    if (equalities.size() > n) {
        if (diagnostic) *diagnostic = "too many equality constraints for vertex enumeration";
        return false;
    }

    const std::size_t active_needed = n - equalities.size();
    std::vector<std::size_t> active;
    auto choose_active = [&](auto&& self, std::size_t start) -> void {
        if (active.size() == active_needed) {
            solve_active_set(active);
            return;
        }
        const std::size_t remaining = active_needed - active.size();
        for (std::size_t i = start; i + remaining <= inequalities.size(); ++i) {
            active.push_back(i);
            self(self, i + 1);
            active.pop_back();
        }
    };

    choose_active(choose_active, 0);

    if (!have_best) {
        if (diagnostic) *diagnostic = "no feasible solution found";
        return false;
    }

    for (std::size_t j = 0; j < n; ++j) {
        if (mymath::abs(best_solution[j]) <= eps) best_solution[j] = Scalar(0);
        if (mymath::abs(best_solution[j] - lower_bounds[j]) <= eps) {
            best_solution[j] = lower_bounds[j];
        }
        if (upper_bounds[j] < inf_val &&
            mymath::abs(best_solution[j] - upper_bounds[j]) <= eps) {
            best_solution[j] = upper_bounds[j];
        }
    }

    *solution = best_solution;
    *objective_value = best_value;
    return true;
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
    std::vector<Scalar>& x_curr,
    std::vector<std::size_t>& basis_curr,
    const std::vector<Scalar>& c_obj,
    const std::vector<std::vector<Scalar>>& A_full,
    const std::vector<Scalar>& lb_full,
    const std::vector<Scalar>& ub_full,
    std::size_t m_total,
    std::size_t n_full,
    bool minimize,
    std::size_t max_iters,
    Scalar eps) {

    const Scalar inf_val = Scalar(mymath::infinity());

    // 初始化基矩阵之逆
    matrix::Matrix B_inv(m_total, m_total, 0.0L);
    if (!reinvert_basis(B_inv, basis_curr, A_full, m_total)) {
        return false;
    }

    for (std::size_t iter = 0; iter < max_iters; ++iter) {
        // 每 50 次迭代重新求逆以保持数值稳定性
        if (iter > 0 && iter % 50 == 0) {
            if (!reinvert_basis(B_inv, basis_curr, A_full, m_total)) return false;
        }

        // 快速检查变量是否在基中，将定价复杂度降为 O(n)
        std::vector<bool> is_basic_mask(n_full, false);
        for (std::size_t k = 0; k < m_total; ++k) is_basic_mask[basis_curr[k]] = true;

        // 计算对偶变量 y: y' = c_B' * B_inv
        std::vector<Scalar> y(m_total, Scalar(0));
        const Scalar cB_threshold = eps * precision::epsilon<Scalar>();
        for (std::size_t j = 0; j < m_total; ++j) {
            const Scalar cB_val = c_obj[basis_curr[j]];
            if (mymath::abs(cB_val) < cB_threshold) continue;
            for (std::size_t i = 0; i < m_total; ++i) {
                y[i] += cB_val * Scalar(B_inv.at(j, i));
            }
        }

        // 寻找入基变量（定价阶段）
        // 使用 Bland 规则：选择满足条件的最小索引，防止循环
        std::size_t entering = n_full;
        Scalar best_rc = Scalar(0);
        
        for (std::size_t j = 0; j < n_full; ++j) {
            if (is_basic_mask[j]) continue;

            // 计算检验数（reduced cost）
            Scalar rc = c_obj[j];
            for (std::size_t i = 0; i < m_total; ++i) {
                rc -= y[i] * A_full[i][j];
            }

            // 根据变量当前状态选择入基
            Scalar xj = x_curr[j];
            bool at_lower = mymath::abs(xj - lb_full[j]) <= eps;
            bool at_upper = (ub_full[j] < inf_val) &&
                            mymath::abs(xj - ub_full[j]) <= eps;

            if (minimize) {
                // 最小化：检验数为负时可入基
                if ((at_lower && rc < -eps) || (at_upper && rc > eps)) {
                    if (entering == n_full || j < entering) { 
                        best_rc = (at_lower ? rc : -rc);
                        entering = j; 
                    }
                }
            } else {
                // 最大化：检验数为正时可入基
                if ((at_lower && rc > eps) || (at_upper && rc < -eps)) {
                    if (entering == n_full || j < entering) { 
                        best_rc = (at_lower ? -rc : rc);
                        entering = j; 
                    }
                }
            }
        }

        // 检查最优性
        if (entering >= n_full || mymath::abs(best_rc) <= eps) {
            return true; // 找到最优解
        }

        // 计算搜索方向 d = B_inv * A_entering
        std::vector<Scalar> a_enter(m_total);
        for (std::size_t i = 0; i < m_total; ++i) a_enter[i] = A_full[i][entering];

        std::vector<Scalar> d(m_total, Scalar(0));
        for (std::size_t i = 0; i < m_total; ++i) {
            for (std::size_t j = 0; j < m_total; ++j) {
                d[i] += Scalar(B_inv.at(i, j)) * a_enter[j];
            }
        }

        // 比值检验确定步长
        Scalar theta = inf_val;
        std::size_t leaving = m_total;

        // 检查入基变量是否从上界减少
        bool decreasing_entering = false;
        if (ub_full[entering] < inf_val &&
            x_curr[entering] >= ub_full[entering] - eps) {
            decreasing_entering = true;
        }

        // 计算最大允许步长（入基变量本身到达另一边界）
        theta = (ub_full[entering] < inf_val)
                ? (ub_full[entering] - lb_full[entering])
                : inf_val;

        // 对每个基变量进行比值检验
        const Scalar di_threshold = eps * precision::sqrt_epsilon<Scalar>();
        for (std::size_t i = 0; i < m_total; ++i) {
            std::size_t j = basis_curr[i];
            Scalar di = decreasing_entering ? -d[i] : d[i];

            if (di > di_threshold) {
                // 基变量可能到达下界
                Scalar ratio = (x_curr[j] - lb_full[j]) / di;
                if (ratio < 0) ratio = 0;
                if (ratio < theta - eps) {
                    theta = ratio; leaving = i;
                } else if (mymath::abs(ratio - theta) <= eps && leaving < m_total) {
                    // Bland 规则：比值相同时选择索引最小的出基变量
                    if (basis_curr[i] < basis_curr[leaving]) leaving = i;
                }
            } else if (di < -di_threshold) {
                // 基变量可能到达上界
                if (ub_full[j] < inf_val) {
                    Scalar ratio = (ub_full[j] - x_curr[j]) / (-di);
                    if (ratio < 0) ratio = 0;
                    if (ratio < theta - eps) {
                        theta = ratio; leaving = i;
                    } else if (mymath::abs(ratio - theta) <= eps && leaving < m_total) {
                        if (basis_curr[i] < basis_curr[leaving]) leaving = i;
                    }
                }
            }
        }

        if (theta >= inf_val - Scalar(1e9L)) {
            // 如果目标函数系数在该方向上有利且无边界限制，则问题无界
            return false;
        }

        // 更新当前解
        Scalar shift = decreasing_entering ? -theta : theta;
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
    const std::vector<Scalar>& objective,
    const matrix::Matrix& inequality_matrix,
    const std::vector<Scalar>& inequality_rhs,
    const matrix::Matrix& equality_matrix,
    const std::vector<Scalar>& equality_rhs,
    const std::vector<Scalar>& lower_bounds,
    const std::vector<Scalar>& upper_bounds,
    Scalar tolerance,
    std::vector<Scalar>* solution,
    Scalar* objective_value,
    std::string* diagnostic) {

    const Scalar inf_val = Scalar(mymath::infinity());
    const Scalar eps = Scalar(tolerance);
    // 使用精度感知的内部容差
    const Scalar internal_eps = precision::default_absolute_tolerance<Scalar>();

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
        *objective_value = 0.0L;
        return true;
    }

    const std::size_t m_ineq = inequality_matrix.rows;
    const std::size_t m_eq = equality_matrix.rows;
    const std::size_t m_total = m_ineq + m_eq;

    if (m_total == 0) {
        std::vector<Scalar> best(n, 0.0L);
        Scalar best_val = 0.0L;
        for (std::size_t j = 0; j < n; ++j) {
            if (objective[j] >= 0) {
                best[j] = upper_bounds[j];
            } else {
                best[j] = lower_bounds[j];
            }
            best_val += objective[j] * best[j];
        }
        *solution = best;
        *objective_value = best_val;
        return true;
    }

    const std::size_t n_slack = m_ineq;
    const std::size_t n_total = n + n_slack;

    std::vector<Scalar> c(n_total, Scalar(0));
    for (std::size_t j = 0; j < n; ++j) {
        c[j] = Scalar(objective[j]);
    }

    std::vector<Scalar> lb(n_total, Scalar(0));
    std::vector<Scalar> ub(n_total, inf_val);
    for (std::size_t j = 0; j < n; ++j) {
        lb[j] = Scalar(lower_bounds[j]);
        ub[j] = Scalar(upper_bounds[j]);
    }

    std::vector<std::vector<Scalar>> A(m_total, std::vector<Scalar>(n_total, Scalar(0)));
    std::vector<Scalar> b(m_total, Scalar(0));

    for (std::size_t i = 0; i < m_ineq; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            A[i][j] = Scalar(inequality_matrix.at(i, j));
        }
        A[i][n + i] = Scalar(1);
        b[i] = Scalar(inequality_rhs[i]);
    }

    for (std::size_t i = 0; i < m_eq; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            A[m_ineq + i][j] = Scalar(equality_matrix.at(i, j));
        }
        b[m_ineq + i] = Scalar(equality_rhs[i]);
    }

    // 调整 b 以反映 x = lb 时的剩余量，并确保 b >= 0
    for (std::size_t i = 0; i < m_total; ++i) {
        Scalar current_lhs = 0;
        for (std::size_t j = 0; j < n_total; ++j) {
            current_lhs += A[i][j] * lb[j];
        }
        b[i] -= current_lhs;
        
        if (b[i] < -eps) {
            b[i] = -b[i];
            for (std::size_t j = 0; j < n_total; ++j) {
                A[i][j] = -A[i][j];
            }
        } else if (b[i] < 0) {
            b[i] = 0;
        }
    }

    const std::size_t n_art = m_total;
    const std::size_t n_full = n_total + n_art;

    std::vector<std::vector<Scalar>> A_full(m_total, std::vector<Scalar>(n_full, Scalar(0)));
    for (std::size_t i = 0; i < m_total; ++i) {
        for (std::size_t j = 0; j < n_total; ++j) {
            A_full[i][j] = A[i][j];
        }
        A_full[i][n_total + i] = Scalar(1);
    }

    std::vector<Scalar> lb_full(n_full, Scalar(0));
    std::vector<Scalar> ub_full(n_full, inf_val);
    for (std::size_t j = 0; j < n_total; ++j) {
        lb_full[j] = lb[j];
        ub_full[j] = ub[j];
    }

    std::vector<Scalar> c_phase1(n_full, Scalar(0));
    for (std::size_t j = n_total; j < n_full; ++j) {
        c_phase1[j] = Scalar(1);
    }

    std::vector<std::size_t> basis(m_total);
    for (std::size_t i = 0; i < m_total; ++i) {
        basis[i] = n_total + i;
    }

    std::vector<Scalar> x(n_full, Scalar(0));
    for (std::size_t j = 0; j < n_total; ++j) {
        x[j] = lb_full[j];
    }
    for (std::size_t i = 0; i < m_total; ++i) {
        x[n_total + i] = b[i];
    }

    constexpr std::size_t kMaxSimplexIters = 10000;
    bool phase1_ok = simplex_iterate(x, basis, c_phase1, A_full, lb_full, ub_full,
                                      m_total, n_full, true, kMaxSimplexIters, internal_eps);

    auto try_fallback = [&](std::string* diag) {
        if (n + m_total > 15) {
            if (diag) *diag = "problem infeasible or numerical failure, and too large for vertex fallback";
            return false;
        }
        return solve_by_vertex_enumeration(objective,
                                           inequality_matrix,
                                           inequality_rhs,
                                           equality_matrix,
                                           equality_rhs,
                                           lower_bounds,
                                           upper_bounds,
                                           tolerance,
                                           solution,
                                           objective_value,
                                           diag);
    };

    if (!phase1_ok) {
        return try_fallback(diagnostic);
    }

    Scalar art_sum = Scalar(0);
    for (std::size_t j = n_total; j < n_full; ++j) {
        art_sum += mymath::abs(x[j]);
    }
    if (art_sum > eps * Scalar(static_cast<long long>(m_total))) {
        return try_fallback(diagnostic);
    }

    for (std::size_t i = 0; i < m_total; ++i) {
        if (basis[i] >= n_total) {
            // 需要将人工变量从基中移出
            matrix::Matrix B_inv_temp(m_total, m_total, 0.0L);
            if (!reinvert_basis(B_inv_temp, basis, A_full, m_total)) break;

            for (std::size_t j = 0; j < n_total; ++j) {
                // 检查变量 j 是否可以入基替代人工变量 i
                bool in_basis = false;
                for (std::size_t k = 0; k < m_total; ++k) {
                    if (basis[k] == j) { in_basis = true; break; }
                }
                if (in_basis) continue;

                // 计算该列在单纯形表中的第 i 行的值: (B_inv * A_j)[i]
                Scalar tableau_val = 0;
                for (std::size_t k = 0; k < m_total; ++k) {
                    tableau_val += Scalar(B_inv_temp.at(i, k)) * A_full[k][j];
                }

                if (mymath::abs(tableau_val) > eps * precision::sqrt_epsilon<Scalar>()) {
                    basis[i] = j;
                    break;
                }
            }
        }
    }

    std::vector<Scalar> c_phase2(n_full, Scalar(0));
    for (std::size_t j = 0; j < n; ++j) {
        c_phase2[j] = Scalar(-objective[j]);
    }

    bool phase2_ok = simplex_iterate(x, basis, c_phase2, A_full, lb_full, ub_full,
                                      m_total, n_full, true, kMaxSimplexIters, internal_eps);

    if (!phase2_ok) {
        if (diagnostic) *diagnostic = "problem is unbounded or numerical failure in phase 2";
        return false;
    }

    std::vector<Scalar> result(n, 0.0L);
    for (std::size_t j = 0; j < n; ++j) {
        result[j] = (x[j]);
        if (result[j] < lower_bounds[j] - tolerance ||
            result[j] > upper_bounds[j] + tolerance) {
            result[j] = std::max(lower_bounds[j], std::min(upper_bounds[j], result[j]));
        }
    }

    for (std::size_t i = 0; i < m_ineq; ++i) {
        Scalar lhs = Scalar(0);
        for (std::size_t j = 0; j < n; ++j) {
            lhs += Scalar(inequality_matrix.at(i, j)) * Scalar(result[j]);
        }
        if (lhs > Scalar(inequality_rhs[i]) + eps) {
            if (diagnostic) *diagnostic = "solution violates inequality constraint";
            return false;
        }
    }
    for (std::size_t i = 0; i < m_eq; ++i) {
        Scalar lhs = Scalar(0);
        for (std::size_t j = 0; j < n; ++j) {
            lhs += Scalar(equality_matrix.at(i, j)) * Scalar(result[j]);
        }
        if (mymath::abs(lhs - Scalar(equality_rhs[i])) > eps) {
            if (diagnostic) *diagnostic = "solution violates equality constraint";
            return false;
        }
    }

    *solution = result;
    Scalar obj_val = Scalar(0);
    for (std::size_t j = 0; j < n; ++j) {
        obj_val += Scalar(objective[j]) * Scalar(result[j]);
    }
    *objective_value = (obj_val);
    return true;
}

}  // namespace simplex

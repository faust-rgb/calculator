// ============================================================================
// 优化辅助函数实现
// ============================================================================

#include "analysis/optimization/optimization_helpers.h"
#include "analysis/optimization/simplex_engine.h"

#include "math/mymath.h"

#include <algorithm>
#include <queue>
#include <sstream>
#include <stdexcept>

namespace optimization_helpers {

// ============================================================================
// 向量运算
// ============================================================================

Scalar dot_product(const std::vector<Scalar>& lhs, const std::vector<Scalar>& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("vector dimension mismatch");
    }
    Scalar total = Scalar(0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        total += Scalar(lhs[i]) * Scalar(rhs[i]);
    }
    return (total);
}

std::string format_planning_result(const std::vector<Scalar>& solution, Scalar objective) {
    std::ostringstream out;
    out << "x = " << matrix::Matrix::vector(solution).to_string()
        << "\nobjective = " << objective;
    return out.str();
}

// ============================================================================
// 整数规划分支定界
// ============================================================================

namespace {

struct Node {
    std::vector<Scalar> lower;
    std::vector<Scalar> upper;
    Scalar estimated_value;

    // 用于优先队列 (Best-Bound Search)
    // 假设是最大化问题：估计值越大越优先
    bool operator<(const Node& other) const {
        // Tie-breaker: 优先选择深度较大的节点，更早找到可行解
        if (mymath::abs(estimated_value - other.estimated_value) < 1e-6L) {
            return depth < other.depth;
        }
        return estimated_value < other.estimated_value;
    }
    int depth = 0;
};

bool is_integer_val(Scalar val, Scalar eps) {
    return mymath::abs(val - mymath::round(val)) <= eps;
}

Scalar get_infinity() {
    return Scalar(mymath::infinity());
}

}  // namespace

void search_integer_branch_and_bound(IntegerSearchContext& ctx,
                                      const std::vector<Scalar>& initial_lower,
                                      const std::vector<Scalar>& initial_upper) {

    // 预求解阶段：收紧整数变量的边界 (Presolve Bound Tightening)
    std::vector<Scalar> lower = initial_lower;
    std::vector<Scalar> upper = initial_upper;
    Scalar tolerance_scalar = Scalar(ctx.tolerance);

    for (std::size_t idx : *ctx.integer_indices) {
        if (mymath::isfinite(lower[idx])) {
            lower[idx] = mymath::ceil(lower[idx] - tolerance_scalar);
        }
        if (mymath::isfinite(upper[idx])) {
            upper[idx] = mymath::floor(upper[idx] + tolerance_scalar);
        }
        if (lower[idx] > upper[idx]) {
            // 整数变量下界大于上界，直接判定不可行
            return;
        }
    }

    std::priority_queue<Node> nodes;

    // 根节点：初始 LP 松弛
    // 初始估计值设为无穷大，确保根节点首先被探索
    nodes.push({lower, upper, get_infinity(), 0});

    std::vector<Scalar> pseudo_costs_up(ctx.variable_count, 1.0L);
    std::vector<Scalar> pseudo_costs_down(ctx.variable_count, 1.0L);

    while (!nodes.empty()) {
        Node current = nodes.top();
        nodes.pop();

        ++(*ctx.visited_nodes);
        if (*ctx.visited_nodes > ctx.max_nodes) {
            throw std::runtime_error(
                *ctx.command_name + " integer search node limit exceeded after " +
                std::to_string(ctx.max_nodes) + " nodes");
        }

        std::vector<Scalar> sol;
        Scalar obj_val = 0.0L;
        std::string diag;

        // 求解当前节点的 LP 松弛子问题
        bool feasible = simplex::solve_linear_box_problem(
            *ctx.objective, *ctx.inequality_matrix, *ctx.inequality_rhs,
            *ctx.equality_matrix, *ctx.equality_rhs,
            current.lower, current.upper, ctx.tolerance,
            &sol, &obj_val, &diag);

        if (!feasible) {
            continue; // 剪枝：LP 不可行
        }

        Scalar obj_val_scalar = Scalar(obj_val);

        // 边界剪枝：当 LP 上界严格小于已知最优解时安全剪枝
        if (*ctx.found && obj_val_scalar < Scalar(*ctx.best_value) - tolerance_scalar) {
            continue;
        }

        // 检查所有应为整数的变量
        std::size_t branch_var = ctx.variable_count;
        Scalar best_score = Scalar(-1);

        for (std::size_t idx : *ctx.integer_indices) {
            Scalar val = Scalar(sol[idx]);
            if (!is_integer_val(val, tolerance_scalar)) {
                // 分支策略：伪成本启发式 (Pseudo-cost heuristic)
                Scalar frac = val - mymath::floor(val);
                Scalar score = std::min(pseudo_costs_down[idx] * frac, pseudo_costs_up[idx] * (1.0L - frac));
                if (score > best_score) {
                    best_score = score;
                    branch_var = idx;
                }
            }
        }

        if (branch_var == ctx.variable_count) {
            // 所有整数变量确实都取了整数值，这是一个改进的整数可行解
            *ctx.found = true;
            *ctx.best_value = obj_val;
            *ctx.best_solution = sol;
        } else {
            // 需要分支
            Scalar val = Scalar(sol[branch_var]);
            Scalar floor_val = mymath::floor(val + tolerance_scalar);
            Scalar ceil_val = mymath::ceil(val - tolerance_scalar);

            // 更新伪成本估计 (简单启发式，通常使用强分支初始化)
            if (current.estimated_value < get_infinity()) {
                Scalar obj_diff = current.estimated_value - obj_val_scalar;
                Scalar frac_down = val - floor_val;
                Scalar frac_up = ceil_val - val;
                if (frac_down > 1e-4L) pseudo_costs_down[branch_var] = 0.8L * pseudo_costs_down[branch_var] + 0.2L * (obj_diff / frac_down);
                if (frac_up > 1e-4L) pseudo_costs_up[branch_var] = 0.8L * pseudo_costs_up[branch_var] + 0.2L * (obj_diff / frac_up);
            }

            // 下分支节点: x_i <= floor(v)
            Node left = current;
            left.upper[branch_var] = floor_val;
            left.estimated_value = obj_val_scalar;
            left.depth = current.depth + 1;
            if (Scalar(left.upper[branch_var]) >= Scalar(left.lower[branch_var])) {
                nodes.push(left);
            }

            // 上分支节点: x_i >= ceil(v)
            Node right = current;
            right.lower[branch_var] = ceil_val;
            right.estimated_value = obj_val_scalar;
            right.depth = current.depth + 1;
            if (Scalar(right.lower[branch_var]) <= Scalar(right.upper[branch_var])) {
                nodes.push(right);
            }
        }
    }
}

}  // namespace optimization_helpers

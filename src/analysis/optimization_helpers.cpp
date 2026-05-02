// ============================================================================
// 优化辅助函数实现
// ============================================================================

#include "analysis/optimization_helpers.h"
#include "analysis/calculator_simplex.h"

#include "math/mymath.h"

#include <algorithm>
#include <queue>
#include <sstream>
#include <stdexcept>

namespace optimization_helpers {

// ============================================================================
// 向量运算
// ============================================================================

double dot_product(const std::vector<double>& lhs, const std::vector<double>& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("vector dimension mismatch");
    }
    long double total = 0.0L;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        total += static_cast<long double>(lhs[i]) * static_cast<long double>(rhs[i]);
    }
    return static_cast<double>(total);
}

std::string format_planning_result(const std::vector<double>& solution, double objective) {
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
    std::vector<double> lower;
    std::vector<double> upper;
    double estimated_value;

    // 用于优先队列 (Best-First Search)
    // 假设是最大化问题：估计值越大越优先
    bool operator<(const Node& other) const {
        return estimated_value < other.estimated_value;
    }
};

bool is_integer_val(double val, double eps) {
    return mymath::abs(val - mymath::round(val)) <= eps;
}

}  // namespace

void search_integer_branch_and_bound(IntegerSearchContext& ctx,
                                      const std::vector<double>& initial_lower,
                                      const std::vector<double>& initial_upper) {

    std::priority_queue<Node> nodes;
    
    // 根节点：初始 LP 松弛
    // 初始估计值设为无穷大，确保根节点首先被探索
    nodes.push({initial_lower, initial_upper, mymath::infinity()});

    while (!nodes.empty()) {
        Node current = nodes.top();
        nodes.pop();

        ++(*ctx.visited_nodes);
        if (*ctx.visited_nodes > ctx.max_nodes) {
            throw std::runtime_error(
                *ctx.command_name + " integer search node limit exceeded after " +
                std::to_string(ctx.max_nodes) + " nodes");
        }

        std::vector<double> sol;
        double obj_val = 0.0;
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

        // 边界剪枝
        if (*ctx.found && obj_val <= *ctx.best_value + ctx.tolerance) {
            continue; 
        }

        // 检查所有应为整数的变量
        std::size_t branch_var = ctx.variable_count;
        double max_fractionality = -1.0;

        for (std::size_t idx : *ctx.integer_indices) {
            double val = sol[idx];
            if (!is_integer_val(val, ctx.tolerance)) {
                // 分支策略：选取最接近 0.5 的变量（Most fractional）
                double fractionality = mymath::abs(val - mymath::round(val));
                if (fractionality > max_fractionality) {
                    max_fractionality = fractionality;
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
            double val = sol[branch_var];
            
            // 下分支节点: x_i <= floor(v)
            Node left = current;
            left.upper[branch_var] = mymath::floor(val + ctx.tolerance);
            left.estimated_value = obj_val;
            if (left.upper[branch_var] >= left.lower[branch_var]) {
                nodes.push(left);
            }

            // 上分支节点: x_i >= ceil(v)
            Node right = current;
            right.lower[branch_var] = mymath::ceil(val - ctx.tolerance);
            right.estimated_value = obj_val;
            if (right.lower[branch_var] <= right.upper[branch_var]) {
                nodes.push(right);
            }
        }
    }
}

}  // namespace optimization_helpers

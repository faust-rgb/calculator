// ============================================================================
// 优化辅助函数实现
// ============================================================================

#include "optimization_helpers.h"

#include "mymath.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

// 检查约束可行性剪枝
bool check_constraint_feasibility(
    const IntegerSearchContext& ctx,
    std::size_t depth,
    const std::vector<long long>& current_values) {

    const double eps = ctx.tolerance;

    // 检查不等式约束
    for (std::size_t row = 0; row < ctx.inequality_matrix->rows; ++row) {
        long double assigned_total = 0.0L;
        for (std::size_t d = 0; d < depth; ++d) {
            const std::size_t col = (*ctx.integer_indices)[d];
            assigned_total += static_cast<long double>(ctx.inequality_matrix->at(row, col)) *
                              static_cast<long double>(current_values[col]);
        }

        long double minimum_possible = assigned_total;
        for (std::size_t remaining = depth; remaining < ctx.integer_indices->size(); ++remaining) {
            const std::size_t col = (*ctx.integer_indices)[remaining];
            const long double coeff = static_cast<long double>(ctx.inequality_matrix->at(row, col));
            minimum_possible += coeff >= 0.0L
                ? coeff * static_cast<long double>((*ctx.integer_lower)[col])
                : coeff * static_cast<long double>((*ctx.integer_upper)[col]);
        }
        for (std::size_t col : *ctx.continuous_indices) {
            const long double coeff = static_cast<long double>(ctx.inequality_matrix->at(row, col));
            minimum_possible += coeff >= 0.0L
                ? coeff * static_cast<long double>((*ctx.lower_bounds)[col])
                : coeff * static_cast<long double>((*ctx.upper_bounds)[col]);
        }

        if (minimum_possible > static_cast<long double>((*ctx.inequality_rhs)[row]) + eps) {
            return false;
        }
    }

    // 检查等式约束
    for (std::size_t row = 0; row < ctx.equality_matrix->rows; ++row) {
        long double assigned_total = 0.0L;
        for (std::size_t d = 0; d < depth; ++d) {
            const std::size_t col = (*ctx.integer_indices)[d];
            assigned_total += static_cast<long double>(ctx.equality_matrix->at(row, col)) *
                              static_cast<long double>(current_values[col]);
        }

        long double minimum_possible = assigned_total;
        long double maximum_possible = assigned_total;
        for (std::size_t remaining = depth; remaining < ctx.integer_indices->size(); ++remaining) {
            const std::size_t col = (*ctx.integer_indices)[remaining];
            const long double coeff = static_cast<long double>(ctx.equality_matrix->at(row, col));
            minimum_possible += coeff >= 0.0L
                ? coeff * static_cast<long double>((*ctx.integer_lower)[col])
                : coeff * static_cast<long double>((*ctx.integer_upper)[col]);
            maximum_possible += coeff >= 0.0L
                ? coeff * static_cast<long double>((*ctx.integer_upper)[col])
                : coeff * static_cast<long double>((*ctx.integer_lower)[col]);
        }
        for (std::size_t col : *ctx.continuous_indices) {
            const long double coeff = static_cast<long double>(ctx.equality_matrix->at(row, col));
            minimum_possible += coeff >= 0.0L
                ? coeff * static_cast<long double>((*ctx.lower_bounds)[col])
                : coeff * static_cast<long double>((*ctx.upper_bounds)[col]);
            maximum_possible += coeff >= 0.0L
                ? coeff * static_cast<long double>((*ctx.upper_bounds)[col])
                : coeff * static_cast<long double>((*ctx.lower_bounds)[col]);
        }

        const long double target = static_cast<long double>((*ctx.equality_rhs)[row]);
        if (target < minimum_possible - eps || target > maximum_possible + eps) {
            return false;
        }
    }

    return true;
}

// 计算乐观目标值
long double compute_optimistic_objective(
    const IntegerSearchContext& ctx,
    std::size_t depth,
    long double current_objective) {

    long double optimistic = current_objective;
    for (std::size_t remaining = depth; remaining < ctx.integer_indices->size(); ++remaining) {
        const std::size_t col = (*ctx.integer_indices)[remaining];
        const long double coeff = static_cast<long double>((*ctx.objective)[col]);
        optimistic += coeff >= 0.0L
            ? coeff * static_cast<long double>((*ctx.integer_upper)[col])
            : coeff * static_cast<long double>((*ctx.integer_lower)[col]);
    }
    for (std::size_t col : *ctx.continuous_indices) {
        const long double coeff = static_cast<long double>((*ctx.objective)[col]);
        optimistic += coeff >= 0.0L
            ? coeff * static_cast<long double>((*ctx.upper_bounds)[col])
            : coeff * static_cast<long double>((*ctx.lower_bounds)[col]);
    }
    return optimistic;
}

}  // namespace

void search_integer_branch_and_bound(IntegerSearchContext& ctx,
                                      std::size_t depth,
                                      long double current_objective) {
    ++(*ctx.visited_nodes);
    if (*ctx.visited_nodes > ctx.max_nodes) {
        throw std::runtime_error(
            *ctx.command_name + " integer search node limit exceeded after " +
            std::to_string(ctx.max_nodes) + " nodes");
    }

    // 约束可行性剪枝
    if (!check_constraint_feasibility(ctx, depth, *ctx.current_integer_values)) {
        return;
    }

    // 目标值剪枝
    const long double optimistic = compute_optimistic_objective(ctx, depth, current_objective);
    if (*ctx.found && optimistic <= static_cast<long double>(*ctx.best_value) + ctx.tolerance) {
        return;
    }

    // 到达叶子节点
    if (depth == ctx.integer_indices->size()) {
        // 构建候选解
        std::vector<double> candidate(ctx.variable_count, 0.0);
        for (std::size_t col = 0; col < ctx.variable_count; ++col) {
            candidate[col] = (*ctx.lower_bounds)[col];
        }
        for (std::size_t col : *ctx.integer_indices) {
            candidate[col] = static_cast<double>((*ctx.current_integer_values)[col]);
        }

        // 检查可行性
        const double eps = ctx.tolerance;
        for (std::size_t row = 0; row < ctx.inequality_matrix->rows; ++row) {
            long double total = 0.0L;
            for (std::size_t col = 0; col < ctx.variable_count; ++col) {
                total += static_cast<long double>(ctx.inequality_matrix->at(row, col)) *
                         static_cast<long double>(candidate[col]);
            }
            if (total > static_cast<long double>((*ctx.inequality_rhs)[row]) + eps) {
                return;
            }
        }
        for (std::size_t row = 0; row < ctx.equality_matrix->rows; ++row) {
            long double total = 0.0L;
            for (std::size_t col = 0; col < ctx.variable_count; ++col) {
                total += static_cast<long double>(ctx.equality_matrix->at(row, col)) *
                         static_cast<long double>(candidate[col]);
            }
            if (mymath::abs(static_cast<double>(total - (*ctx.equality_rhs)[row])) > eps) {
                return;
            }
        }

        // 更新最优解
        const double obj_value = dot_product(*ctx.objective, candidate);
        if (!*ctx.found || obj_value > *ctx.best_value + eps) {
            *ctx.found = true;
            *ctx.best_value = obj_value;
            *ctx.best_solution = candidate;
        }
        return;
    }

    // 分支搜索
    const std::size_t current_col = (*ctx.integer_indices)[depth];
    const bool descending = (*ctx.objective)[current_col] >= 0.0;

    if (descending) {
        for (long long value = (*ctx.integer_upper)[current_col];
             value >= (*ctx.integer_lower)[current_col]; --value) {
            (*ctx.current_integer_values)[current_col] = value;
            search_integer_branch_and_bound(
                ctx, depth + 1,
                current_objective + static_cast<long double>((*ctx.objective)[current_col]) *
                                    static_cast<long double>(value));
        }
    } else {
        for (long long value = (*ctx.integer_lower)[current_col];
             value <= (*ctx.integer_upper)[current_col]; ++value) {
            (*ctx.current_integer_values)[current_col] = value;
            search_integer_branch_and_bound(
                ctx, depth + 1,
                current_objective + static_cast<long double>((*ctx.objective)[current_col]) *
                                    static_cast<long double>(value));
        }
    }
}

}  // namespace optimization_helpers

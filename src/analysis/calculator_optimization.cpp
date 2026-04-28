// ============================================================================
// 线性规划优化命令实现
// ============================================================================

#include "calculator_optimization.h"

#include "mymath.h"
#include "optimization_helpers.h"
#include "calculator_simplex.h"

#include <algorithm>
#include <vector>

namespace optimization {

namespace {

std::vector<double> matrix_to_vector_values(const matrix::Matrix& value,
                                            const std::string& context) {
    if (!value.is_vector()) {
        throw std::runtime_error(context + " expects vector arguments");
    }
    const std::size_t size = value.rows == 1 ? value.cols : value.rows;
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < size; ++i) {
        result[i] = value.rows == 1 ? value.at(0, i) : value.at(i, 0);
    }
    return result;
}

std::string format_planning_result(const OptimizationContext& ctx,
                                   const std::vector<double>& solution,
                                   double objective) {
    return "x = " + matrix::Matrix::vector(solution).to_string() +
           "\nobjective = " + format_decimal(ctx.normalize_result(objective));
}

}  // namespace

ProblemType get_problem_type(const std::string& command) {
    if (command == "lp_max" || command == "lp_min") {
        return ProblemType::LP;
    }
    if (command == "ilp_max" || command == "ilp_min") {
        return ProblemType::ILP;
    }
    if (command == "milp_max" || command == "milp_min") {
        return ProblemType::MILP;
    }
    if (command == "bip_max" || command == "bip_min" ||
        command == "binary_max" || command == "binary_min") {
        return ProblemType::BIP;
    }
    return ProblemType::LP;  // 默认
}

bool is_maximize(const std::string& command) {
    return command == "lp_max" || command == "ilp_max" ||
           command == "milp_max" || command == "bip_max" ||
           command == "binary_max";
}

bool is_optimization_command(const std::string& command) {
    return command == "lp_max" || command == "lp_min" ||
           command == "ilp_max" || command == "ilp_min" ||
           command == "milp_max" || command == "milp_min" ||
           command == "bip_max" || command == "bip_min" ||
           command == "binary_max" || command == "binary_min";
}

std::string normalize_optimization_command(const std::string& command) {
    if (command == "binary_max") return "bip_max";
    if (command == "binary_min") return "bip_min";
    return command;
}

bool solve_linear_box_problem(const std::vector<double>& objective,
                              const matrix::Matrix& inequality_matrix,
                              const std::vector<double>& inequality_rhs,
                              const matrix::Matrix& equality_matrix,
                              const std::vector<double>& equality_rhs,
                              const std::vector<double>& lower_bounds,
                              const std::vector<double>& upper_bounds,
                              double planning_tolerance,
                              std::vector<double>* solution,
                              double* objective_value,
                              std::string* diagnostic) {
    // 委托给专业的单纯形法实现，告别基枚举组合爆炸
    return simplex::solve_linear_box_problem(
        objective, inequality_matrix, inequality_rhs,
        equality_matrix, equality_rhs,
        lower_bounds, upper_bounds, planning_tolerance,
        solution, objective_value, diagnostic);
}

bool handle_optimization_command(const OptimizationContext& ctx,
                                 const std::string& command,
                                 const std::string& inside,
                                 std::string* output) {
    const std::string normalized = normalize_optimization_command(command);
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    const ProblemType problem_type = get_problem_type(normalized);
    const bool maximize = is_maximize(normalized);
    const double planning_tolerance = 1e-8;

    // 解析目标函数
    const std::vector<double> objective = matrix_to_vector_values(
        ctx.parse_matrix_argument(arguments[0], normalized), normalized);
    const std::size_t variable_count = objective.size();

    // 解析约束
    std::size_t argument_index = 1;
    const matrix::Matrix inequality_matrix =
        ctx.parse_matrix_argument(arguments[argument_index++], normalized);
    const std::vector<double> inequality_rhs = matrix_to_vector_values(
        ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);

    matrix::Matrix equality_matrix(0, variable_count, 0.0);
    std::vector<double> equality_rhs;
    std::vector<double> lower_bounds(variable_count, 0.0);
    std::vector<double> upper_bounds(variable_count, 1e20); // 默认无上限
    std::vector<double> integrality(variable_count, 0.0);

    if (problem_type == ProblemType::BIP) {
        if (arguments.size() != 3 && arguments.size() != 5) {
            throw std::runtime_error(normalized + " expects objective, A, b [, Aeq, beq]");
        }
        if (arguments.size() == 5) {
            equality_matrix = ctx.parse_matrix_argument(arguments[argument_index++], normalized);
            equality_rhs = matrix_to_vector_values(ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        }
        std::fill(lower_bounds.begin(), lower_bounds.end(), 0.0);
        std::fill(upper_bounds.begin(), upper_bounds.end(), 1.0);
        std::fill(integrality.begin(), integrality.end(), 1.0);
    } else if (problem_type == ProblemType::MILP) {
        if (arguments.size() != 6 && arguments.size() != 8) {
            throw std::runtime_error(normalized + " expects objective, A, b, lb, ub, integrality [, Aeq, beq]");
        }
        if (arguments.size() == 8) {
            equality_matrix = ctx.parse_matrix_argument(arguments[argument_index++], normalized);
            equality_rhs = matrix_to_vector_values(ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        }
        lower_bounds = matrix_to_vector_values(ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        upper_bounds = matrix_to_vector_values(ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        integrality = matrix_to_vector_values(ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
    } else {
        if (arguments.size() != 5 && arguments.size() != 7) {
            throw std::runtime_error(normalized + " expects objective, A, b, lb, ub [, Aeq, beq]");
        }
        if (arguments.size() == 7) {
            equality_matrix = ctx.parse_matrix_argument(arguments[argument_index++], normalized);
            equality_rhs = matrix_to_vector_values(ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        }
        lower_bounds = matrix_to_vector_values(ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        upper_bounds = matrix_to_vector_values(ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        if (problem_type == ProblemType::ILP) {
            std::fill(integrality.begin(), integrality.end(), 1.0);
        }
    }

    // 转换最小化为最大化
    std::vector<double> simplex_objective = objective;
    if (!maximize) {
        for (double& val : simplex_objective) val = -val;
    }

    // 识别整数索引
    std::vector<std::size_t> integer_indices;
    for (std::size_t i = 0; i < variable_count; ++i) {
        if (std::abs(integrality[i]) > planning_tolerance) integer_indices.push_back(i);
    }

    std::vector<double> final_solution;
    double final_obj = 0.0;

    if (integer_indices.empty()) {
        // 纯 LP
        std::string diag;
        if (!simplex::solve_linear_box_problem(simplex_objective, inequality_matrix, inequality_rhs,
                                               equality_matrix, equality_rhs, lower_bounds, upper_bounds,
                                               planning_tolerance, &final_solution, &final_obj, &diag)) {
            throw std::runtime_error(normalized + " failed: " + diag);
        }
    } else {
        // 真正的大型分支定界
        bool found = false;
        double best_val = -std::numeric_limits<double>::infinity();
        std::vector<double> best_sol;
        std::size_t visited = 0;

        optimization_helpers::IntegerSearchContext bb_ctx;
        bb_ctx.variable_count = variable_count;
        bb_ctx.objective = &simplex_objective;
        bb_ctx.inequality_matrix = &inequality_matrix;
        bb_ctx.inequality_rhs = &inequality_rhs;
        bb_ctx.equality_matrix = &equality_matrix;
        bb_ctx.equality_rhs = &equality_rhs;
        bb_ctx.integer_indices = &integer_indices;
        bb_ctx.tolerance = planning_tolerance;
        bb_ctx.found = &found;
        bb_ctx.best_value = &best_val;
        bb_ctx.best_solution = &best_sol;
        bb_ctx.visited_nodes = &visited;
        bb_ctx.max_nodes = 50000; // 调高节点上限
        bb_ctx.command_name = &normalized;

        optimization_helpers::search_integer_branch_and_bound(bb_ctx, lower_bounds, upper_bounds);

        if (!found) {
            throw std::runtime_error(normalized + " found no feasible integer solution");
        }
        final_solution = best_sol;
        final_obj = best_val;
    }

    *output = format_planning_result(ctx, final_solution, 
                                     optimization_helpers::dot_product(objective, final_solution));
    return true;
}

}  // namespace optimization

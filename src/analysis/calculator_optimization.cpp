// ============================================================================
// 线性规划优化命令实现
// ============================================================================
//
// 本文件实现了线性规划优化命令的数值计算，包括：
// - lp_max / lp_min: 线性规划
// - ilp_max / ilp_min: 整数规划
// - milp_max / milp_min: 混合整数规划
// - bip_max / bip_min: 二进制规划

#include "calculator_optimization.h"

#include "mymath.h"
#include "optimization_helpers.h"
#include "calculator_simplex.h"

#include <algorithm>
#include <vector>

namespace optimization {

enum class ProblemType { LP, ILP, MILP, BIP };

namespace {

/**
 * @brief 将矩阵转换为向量
 */
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

/**
 * @brief 格式化规划问题结果
 */
std::string format_planning_result(const OptimizationContext& ctx,
                                   const std::vector<double>& solution,
                                   double objective) {
    return "x = " + matrix::Matrix::vector(solution).to_string() +
           "\nobjective = " + format_decimal(ctx.normalize_result(objective));
}

}  // namespace

/**
 * @brief 解析问题类型
 */
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

/**
 * @brief 检查是否为最大化问题
 */
bool is_maximize(const std::string& command) {
    return command == "lp_max" || command == "ilp_max" ||
           command == "milp_max" || command == "bip_max" ||
           command == "binary_max";
}

/**
 * @brief 检查是否为优化命令
 */
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

/**
 * @brief 求解带上下界和线性约束的连续线性规划
 *
 * 委托给单纯形法实现。
 */
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
    // 委托给专业的单纯形法实现
    return simplex::solve_linear_box_problem(
        objective, inequality_matrix, inequality_rhs,
        equality_matrix, equality_rhs,
        lower_bounds, upper_bounds, planning_tolerance,
        solution, objective_value, diagnostic);
}

/**
 * @brief 处理优化命令
 *
 * 解析参数并调用相应的求解器：
 * - 纯线性规划：直接使用单纯形法
 * - 整数/混合整数规划：使用分支定界法
 */
bool handle_optimization_command(const OptimizationContext& ctx,
                                 const std::string& command,
                                 const std::string& inside,
                                 std::string* output) {
    const std::string normalized = normalize_optimization_command(command);
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    const ProblemType problem_type = get_problem_type(normalized);
    const bool maximize = is_maximize(normalized);
    const double planning_tolerance = 1e-8;

    // 解析目标函数系数
    const std::vector<double> objective = matrix_to_vector_values(
        ctx.parse_matrix_argument(arguments[0], normalized), normalized);
    const std::size_t variable_count = objective.size();

    // 解析约束条件
    std::size_t argument_index = 1;
    const matrix::Matrix inequality_matrix =
        ctx.parse_matrix_argument(arguments[argument_index++], normalized);
    const std::vector<double> inequality_rhs = matrix_to_vector_values(
        ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);

    // 初始化可选约束
    matrix::Matrix equality_matrix(0, variable_count, 0.0);
    std::vector<double> equality_rhs;
    std::vector<double> lower_bounds(variable_count, 0.0);
    std::vector<double> upper_bounds(variable_count, 1e20); // 默认无上限
    std::vector<double> integrality(variable_count, 0.0);

    // 根据问题类型解析不同参数
    if (problem_type == ProblemType::BIP) {
        // 二进制规划：变量必须为 0 或 1
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
        // 混合整数规划：部分变量必须为整数
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
        // 纯线性规划或整数规划
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
            // 整数规划：所有变量必须为整数
            std::fill(integrality.begin(), integrality.end(), 1.0);
        }
    }

    // 转换最小化为最大化（改变目标函数符号）
    std::vector<double> simplex_objective = objective;
    if (!maximize) {
        for (double& val : simplex_objective) val = -val;
    }

    // 识别整数变量索引
    std::vector<std::size_t> integer_indices;
    for (std::size_t i = 0; i < variable_count; ++i) {
        if (mymath::abs(integrality[i]) > planning_tolerance) integer_indices.push_back(i);
    }

    std::vector<double> final_solution;
    double final_obj = 0.0;

    if (integer_indices.empty()) {
        // 纯线性规划：直接使用单纯形法
        std::string diag;
        if (!simplex::solve_linear_box_problem(simplex_objective, inequality_matrix, inequality_rhs,
                                               equality_matrix, equality_rhs, lower_bounds, upper_bounds,
                                               planning_tolerance, &final_solution, &final_obj, &diag)) {
            throw std::runtime_error(normalized + " failed: " + diag);
        }
    } else {
        // 整数/混合整数规划：使用分支定界法
        bool found = false;
        double best_val = -mymath::infinity();
        std::vector<double> best_sol;
        std::size_t visited = 0;

        // 设置分支定界上下文
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
        bb_ctx.max_nodes = 50000; // 节点上限
        bb_ctx.command_name = &normalized;

        // 执行分支定界搜索
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


std::string OptimizationModule::execute_args(const std::string& command,
                                            const std::vector<std::string>& args,
                                            const CoreServices& services) {
    OptimizationContext ctx;
    ctx.parse_matrix_argument = services.parse_matrix_argument;
    ctx.normalize_result = services.evaluation.normalize_result;
    ctx.is_integer_double = services.is_integer_double;
    ctx.round_to_long_long = services.round_to_long_long;

    std::string inside;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += args[i];
    }

    std::string output;
    if (handle_optimization_command(ctx, command, inside, &output)) {
        return output;
    }
    throw std::runtime_error("Optimization command failed: " + command);
}

std::vector<std::string> OptimizationModule::get_commands() const {
    return {"lp_max", "lp_min", "ilp_max", "ilp_min", "milp_max", "milp_min", "bip_max", "bip_min", "binary_max", "binary_min"};
}

std::string OptimizationModule::get_help_snippet(const std::string& topic) const {
    if (topic == "planning") {
        return "Optimization:\n"
               "  lp_max(c, A, b, lb, ub)       Linear programming (max)\n"
               "  lp_min(c, A, b, lb, ub)       Linear programming (min)\n"
               "  ilp_max(c, A, b, lb, ub)      Integer linear programming\n"
               "  milp_max(c, A, b, lb, ub, int) Mixed-integer programming\n"
               "  bip_max(c, A, b)              Binary integer programming";
    }
    return "";
}

}  // namespace optimization

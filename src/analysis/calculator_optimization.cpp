// ============================================================================
// 线性规划优化命令实现
// ============================================================================

#include "calculator_optimization.h"

#include "mymath.h"
#include "optimization_helpers.h"

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
    const std::size_t variable_count = objective.size();
    if (inequality_matrix.cols != variable_count ||
        inequality_rhs.size() != inequality_matrix.rows ||
        equality_matrix.cols != variable_count ||
        equality_rhs.size() != equality_matrix.rows ||
        lower_bounds.size() != variable_count ||
        upper_bounds.size() != variable_count) {
        throw std::runtime_error("planning dimension mismatch");
    }

    auto feasible_solution = [&](const std::vector<double>& x) {
        if (x.size() != variable_count) {
            return false;
        }
        for (std::size_t col = 0; col < variable_count; ++col) {
            if (x[col] < lower_bounds[col] - planning_tolerance ||
                x[col] > upper_bounds[col] + planning_tolerance) {
                return false;
            }
        }
        for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
            long double total = 0.0L;
            for (std::size_t col = 0; col < variable_count; ++col) {
                total += static_cast<long double>(inequality_matrix.at(row, col)) *
                         static_cast<long double>(x[col]);
            }
            if (total >
                static_cast<long double>(inequality_rhs[row]) + planning_tolerance) {
                return false;
            }
        }
        for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
            long double total = 0.0L;
            for (std::size_t col = 0; col < variable_count; ++col) {
                total += static_cast<long double>(equality_matrix.at(row, col)) *
                         static_cast<long double>(x[col]);
            }
            if (mymath::abs(static_cast<double>(total - equality_rhs[row])) >
                planning_tolerance) {
                return false;
            }
        }
        return true;
    };

    if (variable_count == 0) {
        const std::vector<double> empty_solution;
        if (!feasible_solution(empty_solution)) {
            return false;
        }
        *solution = empty_solution;
        *objective_value = 0.0;
        return true;
    }

    static constexpr std::size_t kMaxPlanningBasisCount = 250000;
    std::vector<std::vector<double>> all_rows;
    std::vector<double> all_rhs;
    all_rows.reserve(inequality_matrix.rows + 2 * equality_matrix.rows +
                     2 * variable_count);
    all_rhs.reserve(inequality_matrix.rows + 2 * equality_matrix.rows +
                    2 * variable_count);

    for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
        std::vector<double> coefficients(variable_count, 0.0);
        for (std::size_t col = 0; col < variable_count; ++col) {
            coefficients[col] = inequality_matrix.at(row, col);
        }
        all_rows.push_back(coefficients);
        all_rhs.push_back(inequality_rhs[row]);
    }
    for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
        std::vector<double> positive_row(variable_count, 0.0);
        std::vector<double> negative_row(variable_count, 0.0);
        for (std::size_t col = 0; col < variable_count; ++col) {
            positive_row[col] = equality_matrix.at(row, col);
            negative_row[col] = -equality_matrix.at(row, col);
        }
        all_rows.push_back(positive_row);
        all_rhs.push_back(equality_rhs[row]);
        all_rows.push_back(negative_row);
        all_rhs.push_back(-equality_rhs[row]);
    }
    for (std::size_t col = 0; col < variable_count; ++col) {
        std::vector<double> upper_row(variable_count, 0.0);
        upper_row[col] = 1.0;
        all_rows.push_back(upper_row);
        all_rhs.push_back(upper_bounds[col]);

        std::vector<double> lower_row(variable_count, 0.0);
        lower_row[col] = -1.0;
        all_rows.push_back(lower_row);
        all_rhs.push_back(-lower_bounds[col]);
    }

    auto capped_choose = [](std::size_t n, std::size_t k, std::size_t cap) {
        if (k > n) {
            return std::size_t{0};
        }
        if (k > n - k) {
            k = n - k;
        }

        long double value = 1.0L;
        for (std::size_t i = 1; i <= k; ++i) {
            value *= static_cast<long double>(n - k + i);
            value /= static_cast<long double>(i);
            if (value > static_cast<long double>(cap)) {
                return cap + 1;
            }
        }
        return static_cast<std::size_t>(value + 0.5L);
    };

    const std::size_t basis_count_estimate =
        capped_choose(all_rows.size(), variable_count, kMaxPlanningBasisCount);
    if (basis_count_estimate > kMaxPlanningBasisCount) {
        throw std::runtime_error(
            "planning basis enumeration limit exceeded: " +
            std::to_string(basis_count_estimate) +
            "+ candidate bases for " + std::to_string(variable_count) +
            " variables and " + std::to_string(all_rows.size()) +
            " active constraints");
    }

    bool found = false;
    std::vector<double> best_solution(variable_count, 0.0);
    double best_value = 0.0;
    std::vector<std::size_t> selection(variable_count, 0);
    std::size_t checked_bases = 0;
    std::size_t singular_bases = 0;
    std::size_t infeasible_bases = 0;

    std::function<void(std::size_t, std::size_t)> enumerate_bases =
        [&](std::size_t start, std::size_t depth) {
            if (depth == variable_count) {
                ++checked_bases;
                matrix::Matrix basis(variable_count, variable_count, 0.0);
                matrix::Matrix basis_rhs(variable_count, 1, 0.0);
                for (std::size_t row = 0; row < variable_count; ++row) {
                    basis_rhs.at(row, 0) = all_rhs[selection[row]];
                    for (std::size_t col = 0; col < variable_count; ++col) {
                        basis.at(row, col) = all_rows[selection[row]][col];
                    }
                }

                try {
                    const matrix::Matrix solved = matrix::solve(basis, basis_rhs);
                    std::vector<double> candidate(variable_count, 0.0);
                    for (std::size_t i = 0; i < variable_count; ++i) {
                        const double value = solved.at(i, 0);
                        candidate[i] =
                            mymath::abs(value) <= planning_tolerance ? 0.0 : value;
                    }
                    if (!feasible_solution(candidate)) {
                        ++infeasible_bases;
                        return;
                    }

                    const double candidate_objective =
                        optimization_helpers::dot_product(objective, candidate);
                    if (!found ||
                        candidate_objective > best_value + planning_tolerance) {
                        found = true;
                        best_value = candidate_objective;
                        best_solution = candidate;
                    }
                } catch (const std::exception&) {
                    ++singular_bases;
                }
                return;
            }

            const std::size_t remaining = variable_count - depth;
            for (std::size_t index = start;
                 index + remaining <= all_rows.size();
                 ++index) {
                selection[depth] = index;
                enumerate_bases(index + 1, depth + 1);
            }
        };

    enumerate_bases(0, 0);
    if (!found) {
        if (diagnostic != nullptr) {
            *diagnostic = "checked " + std::to_string(checked_bases) +
                          " bases, skipped " + std::to_string(singular_bases) +
                          " singular bases and " +
                          std::to_string(infeasible_bases) +
                          " infeasible candidates";
        }
        return false;
    }

    *solution = best_solution;
    *objective_value = best_value;
    return true;
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

    // 解析不等式约束
    std::size_t argument_index = 1;
    const matrix::Matrix inequality_matrix =
        ctx.parse_matrix_argument(arguments[argument_index++], normalized);
    const std::vector<double> inequality_rhs = matrix_to_vector_values(
        ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);

    // 初始化默认值
    matrix::Matrix equality_matrix(0, variable_count, 0.0);
    std::vector<double> equality_rhs;
    std::vector<double> lower_bounds(variable_count, 0.0);
    std::vector<double> upper_bounds(variable_count, 0.0);
    std::vector<double> integrality(variable_count, 0.0);

    // 根据问题类型解析参数
    if (problem_type == ProblemType::BIP) {
        if (arguments.size() != 3 && arguments.size() != 5) {
            throw std::runtime_error(
                normalized +
                " expects objective_vector, A, b, and optional Aeq, beq");
        }
        if (arguments.size() == 5) {
            equality_matrix =
                ctx.parse_matrix_argument(arguments[argument_index++], normalized);
            equality_rhs = matrix_to_vector_values(
                ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        }
        std::fill(lower_bounds.begin(), lower_bounds.end(), 0.0);
        std::fill(upper_bounds.begin(), upper_bounds.end(), 1.0);
        std::fill(integrality.begin(), integrality.end(), 1.0);
    } else if (problem_type == ProblemType::MILP) {
        if (arguments.size() != 6 && arguments.size() != 8) {
            throw std::runtime_error(
                normalized +
                " expects objective_vector, A, b, lower_bounds, upper_bounds, integrality, and optional Aeq, beq");
        }
        if (arguments.size() == 8) {
            equality_matrix =
                ctx.parse_matrix_argument(arguments[argument_index++], normalized);
            equality_rhs = matrix_to_vector_values(
                ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        }
        lower_bounds = matrix_to_vector_values(
            ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        upper_bounds = matrix_to_vector_values(
            ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        integrality = matrix_to_vector_values(
            ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
    } else {
        // LP or ILP
        if (arguments.size() != 5 && arguments.size() != 7) {
            throw std::runtime_error(
                normalized +
                " expects objective_vector, A, b, lower_bounds, upper_bounds, and optional Aeq, beq");
        }
        if (arguments.size() == 7) {
            equality_matrix =
                ctx.parse_matrix_argument(arguments[argument_index++], normalized);
            equality_rhs = matrix_to_vector_values(
                ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        }
        lower_bounds = matrix_to_vector_values(
            ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        upper_bounds = matrix_to_vector_values(
            ctx.parse_matrix_argument(arguments[argument_index++], normalized), normalized);
        if (problem_type == ProblemType::ILP) {
            std::fill(integrality.begin(), integrality.end(), 1.0);
        }
    }

    if (argument_index != arguments.size()) {
        throw std::runtime_error(normalized + " received an invalid argument count");
    }

    // 维度检查
    if (inequality_matrix.cols != variable_count ||
        inequality_rhs.size() != inequality_matrix.rows ||
        equality_matrix.cols != variable_count ||
        equality_rhs.size() != equality_matrix.rows ||
        lower_bounds.size() != variable_count ||
        upper_bounds.size() != variable_count ||
        integrality.size() != variable_count) {
        throw std::runtime_error(normalized + " dimension mismatch");
    }

    // 边界检查
    for (std::size_t i = 0; i < variable_count; ++i) {
        if (lower_bounds[i] > upper_bounds[i]) {
            throw std::runtime_error(normalized + " requires lower_bounds <= upper_bounds");
        }
        if (integrality[i] != 0.0 && !is_integer_double(lower_bounds[i])) {
            throw std::runtime_error(normalized + " requires integer lower bounds for integer variables");
        }
        if (integrality[i] != 0.0 && !is_integer_double(upper_bounds[i])) {
            throw std::runtime_error(normalized + " requires integer upper bounds for integer variables");
        }
    }

    // 转换目标函数（最小化转为最大化）
    std::vector<double> transformed_objective = objective;
    if (!maximize) {
        for (double& value : transformed_objective) {
            value = -value;
        }
    }

    // 检查是否有整数变量
    std::vector<std::size_t> integer_indices;
    for (std::size_t i = 0; i < variable_count; ++i) {
        if (!mymath::is_near_zero(integrality[i], planning_tolerance)) {
            integer_indices.push_back(i);
        }
    }

    // 纯线性规划（无整数变量）
    if (integer_indices.empty()) {
        std::vector<double> best_solution;
        double best_value = 0.0;
        std::string planning_diagnostic;
        if (!solve_linear_box_problem(transformed_objective,
                                      inequality_matrix,
                                      inequality_rhs,
                                      equality_matrix,
                                      equality_rhs,
                                      lower_bounds,
                                      upper_bounds,
                                      planning_tolerance,
                                      &best_solution,
                                      &best_value,
                                      &planning_diagnostic)) {
            std::string message =
                normalized + " found no feasible bounded solution";
            if (!planning_diagnostic.empty()) {
                message += " (" + planning_diagnostic + ")";
            }
            throw std::runtime_error(message);
        }
        *output = format_planning_result(
            ctx,
            best_solution,
            optimization_helpers::dot_product(objective, best_solution));
        return true;
    }

    // 整数规划 - 分支定界算法
    static constexpr std::size_t kMaxIntegerAssignments = 1000000;
    static constexpr std::size_t kMaxIntegerSearchNodes = 2000000;

    std::vector<std::size_t> continuous_indices;
    for (std::size_t i = 0; i < variable_count; ++i) {
        if (mymath::is_near_zero(integrality[i], planning_tolerance)) {
            continuous_indices.push_back(i);
        }
    }

    std::size_t estimated_integer_assignments = 1;
    std::vector<long long> integer_lower(variable_count, 0);
    std::vector<long long> integer_upper(variable_count, 0);
    for (std::size_t index : integer_indices) {
        integer_lower[index] = ctx.round_to_long_long(lower_bounds[index]);
        integer_upper[index] = ctx.round_to_long_long(upper_bounds[index]);
        const long double width =
            static_cast<long double>(integer_upper[index]) -
            static_cast<long double>(integer_lower[index]) + 1.0L;
        if (width <= 0.0L ||
            width > static_cast<long double>(kMaxIntegerAssignments) ||
            estimated_integer_assignments >
                kMaxIntegerAssignments / static_cast<std::size_t>(width)) {
            estimated_integer_assignments = kMaxIntegerAssignments + 1;
        } else {
            estimated_integer_assignments *= static_cast<std::size_t>(width);
        }
    }
    if (estimated_integer_assignments > kMaxIntegerAssignments) {
        throw std::runtime_error(
            normalized + " integer search limit exceeded: more than " +
            std::to_string(kMaxIntegerAssignments) +
            " possible integer assignments");
    }

    bool found = false;
    double best_value = 0.0;
    std::vector<double> best_solution(variable_count, 0.0);
    std::vector<long long> current_integer_values(variable_count, 0);
    std::size_t visited_integer_nodes = 0;

    std::function<void(std::size_t, long double)> search_integer =
        [&](std::size_t depth, long double current_objective) {
            ++visited_integer_nodes;
            if (visited_integer_nodes > kMaxIntegerSearchNodes) {
                throw std::runtime_error(
                    normalized + " integer search node limit exceeded after " +
                    std::to_string(kMaxIntegerSearchNodes) + " nodes");
            }

            // 检查不等式约束的可行性剪枝
            for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                long double assigned_total = 0.0L;
                for (std::size_t assigned_depth = 0; assigned_depth < depth; ++assigned_depth) {
                    const std::size_t col = integer_indices[assigned_depth];
                    assigned_total += static_cast<long double>(inequality_matrix.at(row, col)) *
                                      static_cast<long double>(current_integer_values[col]);
                }

                long double minimum_possible = assigned_total;
                for (std::size_t remaining_depth = depth;
                     remaining_depth < integer_indices.size();
                     ++remaining_depth) {
                    const std::size_t col = integer_indices[remaining_depth];
                    const long double coefficient =
                        static_cast<long double>(inequality_matrix.at(row, col));
                    minimum_possible +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(integer_lower[col])
                            : coefficient * static_cast<long double>(integer_upper[col]);
                }
                for (std::size_t col : continuous_indices) {
                    const long double coefficient =
                        static_cast<long double>(inequality_matrix.at(row, col));
                    minimum_possible +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(lower_bounds[col])
                            : coefficient * static_cast<long double>(upper_bounds[col]);
                }
                if (minimum_possible >
                    static_cast<long double>(inequality_rhs[row]) + planning_tolerance) {
                    return;
                }
            }

            // 检查等式约束的可行性剪枝
            for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                long double assigned_total = 0.0L;
                for (std::size_t assigned_depth = 0; assigned_depth < depth; ++assigned_depth) {
                    const std::size_t col = integer_indices[assigned_depth];
                    assigned_total += static_cast<long double>(equality_matrix.at(row, col)) *
                                      static_cast<long double>(current_integer_values[col]);
                }

                long double minimum_possible = assigned_total;
                long double maximum_possible = assigned_total;
                for (std::size_t remaining_depth = depth;
                     remaining_depth < integer_indices.size();
                     ++remaining_depth) {
                    const std::size_t col = integer_indices[remaining_depth];
                    const long double coefficient =
                        static_cast<long double>(equality_matrix.at(row, col));
                    minimum_possible +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(integer_lower[col])
                            : coefficient * static_cast<long double>(integer_upper[col]);
                    maximum_possible +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(integer_upper[col])
                            : coefficient * static_cast<long double>(integer_lower[col]);
                }
                for (std::size_t col : continuous_indices) {
                    const long double coefficient =
                        static_cast<long double>(equality_matrix.at(row, col));
                    minimum_possible +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(lower_bounds[col])
                            : coefficient * static_cast<long double>(upper_bounds[col]);
                    maximum_possible +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(upper_bounds[col])
                            : coefficient * static_cast<long double>(lower_bounds[col]);
                }

                const long double target = static_cast<long double>(equality_rhs[row]);
                if (target < minimum_possible - planning_tolerance ||
                    target > maximum_possible + planning_tolerance) {
                    return;
                }
            }

            // 目标函数上界剪枝
            long double optimistic_objective = current_objective;
            for (std::size_t remaining_depth = depth;
                 remaining_depth < integer_indices.size();
                 ++remaining_depth) {
                const std::size_t col = integer_indices[remaining_depth];
                const long double coefficient =
                    static_cast<long double>(transformed_objective[col]);
                optimistic_objective +=
                    coefficient >= 0.0L
                        ? coefficient * static_cast<long double>(integer_upper[col])
                        : coefficient * static_cast<long double>(integer_lower[col]);
            }
            for (std::size_t col : continuous_indices) {
                const long double coefficient =
                    static_cast<long double>(transformed_objective[col]);
                optimistic_objective +=
                    coefficient >= 0.0L
                        ? coefficient * static_cast<long double>(upper_bounds[col])
                        : coefficient * static_cast<long double>(lower_bounds[col]);
            }
            if (found &&
                optimistic_objective <=
                    static_cast<long double>(best_value) + planning_tolerance) {
                return;
            }

            // 到达叶子节点
            if (depth == integer_indices.size()) {
                std::vector<double> candidate(variable_count, 0.0);
                for (std::size_t col = 0; col < variable_count; ++col) {
                    candidate[col] = lower_bounds[col];
                }
                for (std::size_t col : integer_indices) {
                    candidate[col] = static_cast<double>(current_integer_values[col]);
                }

                if (continuous_indices.empty()) {
                    // 纯整数规划 - 直接检查可行性
                    bool feasible = true;
                    for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                        long double total = 0.0L;
                        for (std::size_t col = 0; col < variable_count; ++col) {
                            total += static_cast<long double>(inequality_matrix.at(row, col)) *
                                     static_cast<long double>(candidate[col]);
                        }
                        if (total >
                            static_cast<long double>(inequality_rhs[row]) + planning_tolerance) {
                            feasible = false;
                            break;
                        }
                    }
                    if (!feasible) {
                        return;
                    }
                    for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                        long double total = 0.0L;
                        for (std::size_t col = 0; col < variable_count; ++col) {
                            total += static_cast<long double>(equality_matrix.at(row, col)) *
                                     static_cast<long double>(candidate[col]);
                        }
                        if (mymath::abs(static_cast<double>(total - equality_rhs[row])) >
                            planning_tolerance) {
                            feasible = false;
                            break;
                        }
                    }
                    if (!feasible) {
                        return;
                    }

                    const double objective_value =
                        optimization_helpers::dot_product(transformed_objective, candidate);
                    if (!found || objective_value > best_value + planning_tolerance) {
                        found = true;
                        best_value = objective_value;
                        best_solution = candidate;
                    }
                    return;
                }

                // 混合整数规划 - 求解连续变量子问题
                matrix::Matrix reduced_inequality(inequality_matrix.rows,
                                                  continuous_indices.size(),
                                                  0.0);
                std::vector<double> reduced_inequality_rhs(inequality_rhs.size(), 0.0);
                for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                    long double rhs_adjustment = static_cast<long double>(inequality_rhs[row]);
                    for (std::size_t col : integer_indices) {
                        rhs_adjustment -=
                            static_cast<long double>(inequality_matrix.at(row, col)) *
                            static_cast<long double>(candidate[col]);
                    }
                    reduced_inequality_rhs[row] = static_cast<double>(rhs_adjustment);
                    for (std::size_t reduced_col = 0;
                         reduced_col < continuous_indices.size();
                         ++reduced_col) {
                        reduced_inequality.at(row, reduced_col) =
                            inequality_matrix.at(row, continuous_indices[reduced_col]);
                    }
                }

                matrix::Matrix reduced_equality(equality_matrix.rows,
                                                continuous_indices.size(),
                                                0.0);
                std::vector<double> reduced_equality_rhs(equality_rhs.size(), 0.0);
                for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                    long double rhs_adjustment = static_cast<long double>(equality_rhs[row]);
                    for (std::size_t col : integer_indices) {
                        rhs_adjustment -=
                            static_cast<long double>(equality_matrix.at(row, col)) *
                            static_cast<long double>(candidate[col]);
                    }
                    reduced_equality_rhs[row] = static_cast<double>(rhs_adjustment);
                    for (std::size_t reduced_col = 0;
                         reduced_col < continuous_indices.size();
                         ++reduced_col) {
                        reduced_equality.at(row, reduced_col) =
                            equality_matrix.at(row, continuous_indices[reduced_col]);
                    }
                }

                std::vector<double> reduced_objective(continuous_indices.size(), 0.0);
                std::vector<double> reduced_lower(continuous_indices.size(), 0.0);
                std::vector<double> reduced_upper(continuous_indices.size(), 0.0);
                for (std::size_t reduced_col = 0;
                     reduced_col < continuous_indices.size();
                     ++reduced_col) {
                    const std::size_t original_col = continuous_indices[reduced_col];
                    reduced_objective[reduced_col] = transformed_objective[original_col];
                    reduced_lower[reduced_col] = lower_bounds[original_col];
                    reduced_upper[reduced_col] = upper_bounds[original_col];
                }

                std::vector<double> reduced_solution;
                double reduced_objective_value = 0.0;
                std::string reduced_diagnostic;
                if (!solve_linear_box_problem(reduced_objective,
                                              reduced_inequality,
                                              reduced_inequality_rhs,
                                              reduced_equality,
                                              reduced_equality_rhs,
                                              reduced_lower,
                                              reduced_upper,
                                              planning_tolerance,
                                              &reduced_solution,
                                              &reduced_objective_value,
                                              &reduced_diagnostic)) {
                    return;
                }

                for (std::size_t reduced_col = 0;
                     reduced_col < continuous_indices.size();
                     ++reduced_col) {
                    candidate[continuous_indices[reduced_col]] = reduced_solution[reduced_col];
                }

                const double objective_value =
                    optimization_helpers::dot_product(transformed_objective, candidate);
                if (!found || objective_value > best_value + planning_tolerance) {
                    found = true;
                    best_value = objective_value;
                    best_solution = candidate;
                }
                return;
            }

            // 递归搜索
            const std::size_t current_col = integer_indices[depth];
            const bool descending = transformed_objective[current_col] >= 0.0;
            if (descending) {
                for (long long value = integer_upper[current_col];
                     value >= integer_lower[current_col];
                     --value) {
                    current_integer_values[current_col] = value;
                    search_integer(
                        depth + 1,
                        current_objective +
                            static_cast<long double>(transformed_objective[current_col]) *
                                static_cast<long double>(value));
                }
            } else {
                for (long long value = integer_lower[current_col];
                     value <= integer_upper[current_col];
                     ++value) {
                    current_integer_values[current_col] = value;
                    search_integer(
                        depth + 1,
                        current_objective +
                            static_cast<long double>(transformed_objective[current_col]) *
                                static_cast<long double>(value));
                }
            }
        };

    search_integer(0, 0.0L);
    if (!found) {
        throw std::runtime_error(normalized + " found no feasible mixed-integer solution");
    }

        *output = format_planning_result(
            ctx,
            best_solution,
            optimization_helpers::dot_product(objective, best_solution));
    return true;
}

}  // namespace optimization

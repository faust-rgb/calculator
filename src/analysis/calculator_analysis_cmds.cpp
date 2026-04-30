// ============================================================================
// 函数分析命令实现
// ============================================================================
//
// 本文件实现了函数分析相关的命令处理逻辑，包括：
// - limit: 极限计算（支持双侧极限和单侧极限）
// - critical: 临界点分析（求解梯度为零的点）
// - extrema: 极值点查找（在给定区间内寻找极值点）
// - lagrange: 拉格朗日乘数法求解约束优化问题

#include "calculator_analysis_cmds.h"

#include "mymath.h"

#include <algorithm>
#include <sstream>
#include <vector>

namespace analysis_cmds {

namespace {

/**
 * @brief 在指定点求值符号表达式
 *
 * @param expression 符号表达式
 * @param variables 变量名列表
 * @param values 变量值列表
 * @param result 输出求值结果
 * @return 是否成功求值为数值
 */
bool evaluate_symbolic_at_point(SymbolicExpression expression,
                                const std::vector<std::string>& variables,
                                const std::vector<double>& values,
                                double* result) {
    // 将所有变量替换为对应的数值
    for (std::size_t i = 0; i < variables.size(); ++i) {
        expression = expression.substitute(
            variables[i], SymbolicExpression::number(values[i]));
    }
    // 简化并检查是否为数值
    return expression.simplify().is_number(result);
}

bool is_infinity_literal(const std::string& text) {
    std::string value = trim_copy(text);
    if (!value.empty() && value.front() == '+') {
        value = trim_copy(value.substr(1));
    } else if (!value.empty() && value.front() == '-') {
        value = trim_copy(value.substr(1));
    }
    return value == "inf" || value == "infinity" || value == "oo";
}

}  // namespace

/**
 * @brief 分类临界点类型
 *
 * 通过计算 Hessian 矩阵的特征值来判断临界点的类型：
 * - 所有特征值 > 0：局部极小值
 * - 所有特征值 < 0：局部极大值
 * - 特征值有正有负：鞍点
 * - 存在零特征值：退化点
 *
 * @param hessian Hessian 矩阵的符号表达式
 * @param variables 变量名列表
 * @param values 临界点坐标
 * @return 分类结果字符串
 */
std::string classify_critical_point(
    const std::vector<std::vector<SymbolicExpression>>& hessian,
    const std::vector<std::string>& variables,
    const std::vector<double>& values) {

    const std::size_t n = variables.size();

    // 构建 Hessian 矩阵的数值矩阵
    matrix::Matrix h_mat(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            double val = 0.0;
            if (!evaluate_symbolic_at_point(
                    hessian[row][col], variables, values, &val)) {
                return "unclassified";
            }
            h_mat.at(row, col) = val;
        }
    }

    // 计算特征值
    const matrix::Matrix ev_result = matrix::eigenvalues(h_mat);
    std::vector<double> evals;
    evals.reserve(n);

    // 从特征值结果中提取特征值（处理不同的返回格式）
    if (ev_result.cols == 1 || ev_result.rows == 1) {
        // 向量结果
        std::size_t len = ev_result.rows * ev_result.cols;
        for (std::size_t i = 0; i < len; ++i) {
            evals.push_back(ev_result.data[i]);
        }
    } else if (ev_result.rows == 2 && ev_result.cols == 2) {
        // 2x2 复数结果格式: [real1, imag1; real2, imag2]
        evals.push_back(ev_result.at(0, 0));
        evals.push_back(ev_result.at(1, 0));
    } else {
        // 一般对角矩阵结果
        for (std::size_t i = 0; i < std::min(ev_result.rows, ev_result.cols); ++i) {
            evals.push_back(ev_result.at(i, i));
        }
    }

    // 分析特征值符号
    bool positive = false;
    bool negative = false;
    bool zero = false;
    constexpr double eps = 1e-5;

    for (double ev : evals) {
        if (ev > eps) {
            positive = true;
        } else if (ev < -eps) {
            negative = true;
        } else {
            zero = true;
        }
    }

    // 根据特征值符号判断临界点类型
    if (zero) {
        return "degenerate";  // 退化点：存在零特征值
    }
    if (positive && !negative) {
        return "local min";   // 局部极小值：所有特征值 > 0
    }
    if (negative && !positive) {
        return "local max";   // 局部极大值：所有特征值 < 0
    }
    if (positive && negative) {
        return "saddle";      // 鞍点：特征值有正有负
    }

    return "unclassified";
}

/**
 * @brief 检查是否为分析命令
 */
bool is_analysis_command(const std::string& command) {
    return command == "limit" || command == "critical" || command == "extrema" || command == "lagrange";
}

/**
 * @brief 处理分析命令的统一入口
 *
 * 根据命令类型分发到相应的处理函数：
 * - lagrange: 使用拉格朗日乘数法处理约束优化
 * - limit: 计算函数极限
 * - critical: 查找临界点
 * - extrema: 查找极值点
 *
 * @param ctx 分析上下文
 * @param command 命令名
 * @param inside 括号内的参数字符串
 * @param output 输出结果
 * @return 是否成功处理
 */
bool handle_analysis_command(const AnalysisContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    // ==================== 拉格朗日乘数法 ====================
    if (command == "lagrange") {
        if (arguments.size() < 2) {
            throw std::runtime_error("lagrange expects (expr, constraints, [vars...])");
        }

        // 解析目标函数 f
        std::string variable_name;
        SymbolicExpression f;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &f);

        // 解析约束条件 [g1, g2, ...]
        std::vector<SymbolicExpression> constraints;
        if (arguments[1].front() == '[' && arguments[1].back() == ']') {
             // 多个约束：解析列表形式
             std::vector<std::string> const_strs = split_top_level_arguments(arguments[1].substr(1, arguments[1].size() - 2));
             for (const auto& s : const_strs) constraints.push_back(SymbolicExpression::parse(s));
        } else {
             // 单个约束
             constraints.push_back(SymbolicExpression::parse(arguments[1]));
        }

        // 解析变量列表
        std::vector<std::string> variables = ctx.parse_symbolic_variable_arguments(arguments, 2, f.identifier_variables());

        // 构造拉格朗日函数 L = f - sum(lambda_i * gi)
        SymbolicExpression lagrangian = f;
        std::vector<std::string> all_vars = variables;
        for (std::size_t i = 0; i < constraints.size(); ++i) {
            std::string lambda_var = "L" + std::to_string(i + 1); // 使用 L1, L2, ... 作为拉格朗日乘子（避免 l 和 1 混淆）
            all_vars.push_back(lambda_var);
            lagrangian = (lagrangian - SymbolicExpression::variable(lambda_var) * constraints[i]).simplify();
        }

        // 转换为 critical 命令求解
        std::string lagrangian_str = lagrangian.to_string();
        std::string all_vars_str;
        for (std::size_t i = 0; i < all_vars.size(); ++i) {
            if (i > 0) all_vars_str += ", ";
            all_vars_str += all_vars[i];
        }

        return handle_analysis_command(ctx, "critical", lagrangian_str + ", " + all_vars_str, output);
    }

    // ==================== 极限计算 ====================
    if (command == "limit") {
        if (arguments.size() != 2 && arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "limit expects (expr, point[, direction]) or (expr, variable, point[, direction])");
        }

        std::string point_argument;
        std::size_t direction_index = 0;
        FunctionAnalysis analysis;

        const bool explicit_variable =
            arguments.size() >= 3 &&
            is_identifier_text(trim_copy(arguments[1])) &&
            !is_infinity_literal(arguments[1]);
        if (explicit_variable) {
            const std::string variable_name = trim_copy(arguments[1]);
            std::string inferred_variable;
            SymbolicExpression expression;
            ctx.resolve_symbolic(arguments[0], false, &inferred_variable, &expression);
            analysis = FunctionAnalysis(variable_name);
            analysis.define(expression.to_string());
            point_argument = arguments[2];
            direction_index = 3;
        } else {
            analysis = ctx.build_analysis(arguments[0]);
            point_argument = arguments[1];
            direction_index = 2;
        }

        int direction = 0;  // 默认双侧极限
        if (arguments.size() > direction_index) {
            const double direction_value = ctx.parse_decimal(arguments[direction_index]);
            if (!is_integer_double(direction_value)) {
                throw std::runtime_error("limit direction must be -1, 0, or 1");
            }
            direction = static_cast<int>(round_to_long_long(direction_value));
        }

        *output = format_decimal(ctx.normalize_result(
            analysis.limit(ctx.parse_decimal(point_argument), direction)));
        return true;
    }

    // ==================== 临界点分析 ====================
    if (command == "critical") {
        if (arguments.empty()) {
            throw std::runtime_error(
                "critical expects a symbolic expression and optional variable names");
        }

        // 解析符号表达式
        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            ctx.parse_symbolic_variable_arguments(arguments,
                                                  1,
                                                  expression.identifier_variables());

        // 计算梯度和 Hessian 矩阵
        const std::vector<SymbolicExpression> gradient =
            expression.gradient(variables);
        const std::vector<std::vector<SymbolicExpression>> hessian =
            expression.hessian(variables);

        // 格式化临界点解的 lambda 函数
        auto format_critical_solution = [&](const std::vector<double>& values) {
            std::ostringstream out;
            out << "[";
            for (std::size_t i = 0; i < variables.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << variables[i] << " = "
                    << format_decimal(ctx.normalize_result(values[i]));
            }
            out << "]";
            out << " (" << classify_critical_point(
                hessian, variables, values) << ")";
            return out.str();
        };

        // 求解线性方程组的 lambda 函数（用于仿射梯度）
        auto solve_linear_system =
            [](std::vector<std::vector<double>> matrix,
               std::vector<double> rhs,
               std::vector<double>* solution) {
                const std::size_t size = rhs.size();
                // 高斯消元法求解
                for (std::size_t col = 0; col < size; ++col) {
                    // 选主元
                    std::size_t pivot = col;
                    for (std::size_t row = col + 1; row < size; ++row) {
                        if (mymath::abs(matrix[row][col]) > mymath::abs(matrix[pivot][col])) {
                            pivot = row;
                        }
                    }
                    if (mymath::abs(matrix[pivot][col]) < 1e-12) {
                        return false;  // 矩阵奇异
                    }
                    if (pivot != col) {
                        std::swap(matrix[pivot], matrix[col]);
                        std::swap(rhs[pivot], rhs[col]);
                    }

                    // 消元
                    const double pivot_value = matrix[col][col];
                    for (std::size_t j = col; j < size; ++j) {
                        matrix[col][j] /= pivot_value;
                    }
                    rhs[col] /= pivot_value;

                    for (std::size_t row = 0; row < size; ++row) {
                        if (row == col) {
                            continue;
                        }
                        const double factor = matrix[row][col];
                        if (mymath::abs(factor) < 1e-12) {
                            continue;
                        }
                        for (std::size_t j = col; j < size; ++j) {
                            matrix[row][j] -= factor * matrix[col][j];
                        }
                        rhs[row] -= factor * rhs[col];
                    }
                }
                *solution = rhs;
                return true;
            };

        // 检查梯度是否为仿射函数（线性 + 常数）
        // 如果是，可以直接用线性方程组求解
        std::vector<std::vector<double>> coefficients(
            variables.size(), std::vector<double>(variables.size(), 0.0));
        std::vector<double> rhs(variables.size(), 0.0);
        const std::vector<double> zeros(variables.size(), 0.0);
        bool affine_gradient = true;

        // 提取梯度的线性系数
        for (std::size_t row = 0; row < gradient.size(); ++row) {
            double constant = 0.0;
            if (!evaluate_symbolic_at_point(
                    gradient[row], variables, zeros, &constant)) {
                affine_gradient = false;
                break;
            }
            rhs[row] = -constant;

            // 计算每个变量的系数
            for (std::size_t col = 0; col < variables.size(); ++col) {
                std::vector<double> sample = zeros;
                sample[col] = 1.0;
                double value = 0.0;
                if (!evaluate_symbolic_at_point(
                        gradient[row], variables, sample, &value)) {
                    affine_gradient = false;
                    break;
                }
                coefficients[row][col] = value - constant;
            }
            if (!affine_gradient) {
                break;
            }

            // 验证梯度确实是仿射的（通过多点采样）
            std::vector<std::vector<double>> validation_samples;
            validation_samples.push_back(std::vector<double>(variables.size(), 1.0));
            for (std::size_t col = 0; col < variables.size(); ++col) {
                std::vector<double> sample = zeros;
                sample[col] = 2.0;
                validation_samples.push_back(sample);
            }
            for (const std::vector<double>& sample : validation_samples) {
                double actual = 0.0;
                if (!evaluate_symbolic_at_point(
                        gradient[row], variables, sample, &actual)) {
                    affine_gradient = false;
                    break;
                }
                double predicted = constant;
                for (std::size_t col = 0; col < variables.size(); ++col) {
                    predicted += coefficients[row][col] * sample[col];
                }
                if (!mymath::is_near_zero(actual - predicted, 1e-8)) {
                    affine_gradient = false;
                    break;
                }
            }
            if (!affine_gradient) {
                break;
            }
        }

        // 如果梯度是仿射的，直接求解线性方程组
        if (affine_gradient) {
            std::vector<double> solution;
            if (!solve_linear_system(coefficients, rhs, &solution)) {
                *output = "No isolated critical point.";
                return true;
            }
            *output = format_critical_solution(solution);
            return true;
        }

        // 非仿射梯度：使用牛顿迭代法求解非线性方程组
        if (variables.size() > 3) {
            throw std::runtime_error(
                "critical nonlinear search supports up to 3 variables");
        }

        // 生成多个起始点以避免局部收敛
        std::vector<std::vector<double>> starts = {
            std::vector<double>(variables.size(), 0.0)};
        const std::vector<double> seeds = {-2.0, -1.0, 1.0, 2.0};
        for (std::size_t dimension = 0; dimension < variables.size(); ++dimension) {
            std::vector<std::vector<double>> next = starts;
            for (const std::vector<double>& start : starts) {
                for (double seed : seeds) {
                    std::vector<double> candidate = start;
                    candidate[dimension] = seed;
                    next.push_back(candidate);
                }
            }
            starts.swap(next);
        }

        std::vector<std::vector<double>> solutions;
        // 从每个起始点开始牛顿迭代
        for (std::vector<double> current : starts) {
            bool converged = false;
            for (int iteration = 0; iteration < 40; ++iteration) {
                // 计算当前点的梯度值
                std::vector<double> gradient_values(variables.size(), 0.0);
                double gradient_norm = 0.0;
                bool numeric_ok = true;
                for (std::size_t row = 0; row < variables.size(); ++row) {
                    if (!evaluate_symbolic_at_point(
                            gradient[row], variables, current, &gradient_values[row])) {
                        numeric_ok = false;
                        break;
                    }
                    gradient_norm += gradient_values[row] * gradient_values[row];
                }
                if (!numeric_ok) {
                    break;
                }
                if (gradient_norm < 1e-24) {
                    converged = true;  // 梯度足够小，认为收敛
                    break;
                }

                // 计算 Jacobian 矩阵（即 Hessian 矩阵）
                std::vector<std::vector<double>> jacobian(
                    variables.size(), std::vector<double>(variables.size(), 0.0));
                for (std::size_t row = 0; row < variables.size() && numeric_ok; ++row) {
                    for (std::size_t col = 0; col < variables.size(); ++col) {
                        if (!evaluate_symbolic_at_point(
                                hessian[row][col], variables, current, &jacobian[row][col])) {
                            numeric_ok = false;
                            break;
                        }
                    }
                }
                if (!numeric_ok) {
                    break;
                }

                // 牛顿法：求解 J * delta = -gradient
                for (double& value : gradient_values) {
                    value = -value;
                }
                std::vector<double> step;
                if (!solve_linear_system(jacobian, gradient_values, &step)) {
                    break;  // Jacobian 奇异，无法继续
                }
                double step_norm = 0.0;
                for (std::size_t i = 0; i < current.size(); ++i) {
                    current[i] += step[i];
                    step_norm += step[i] * step[i];
                }
                if (step_norm < 1e-24) {
                    converged = true;  // 步长足够小，认为收敛
                    break;
                }
            }

            if (!converged) {
                continue;
            }

            // 处理接近原点的退化情况
            double current_norm = 0.0;
            for (double value : current) {
                current_norm += value * value;
            }
            if (current_norm < 1e-2 &&
                (classify_critical_point(hessian, variables, current) == "degenerate" ||
                 classify_critical_point(hessian, variables, current) == "local min")) {
                for (double& value : current) {
                    value = 0.0;
                }
            }

            // 检查是否与已找到的解重复
            bool duplicate = false;
            for (const std::vector<double>& existing : solutions) {
                double distance = 0.0;
                for (std::size_t i = 0; i < existing.size(); ++i) {
                    const double diff = existing[i] - current[i];
                    distance += diff * diff;
                }
                if (distance < 1e-10) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                solutions.push_back(current);
            }
        }

        // 输出结果
        if (solutions.empty()) {
            *output = "No isolated critical point.";
            return true;
        }
        std::sort(solutions.begin(), solutions.end());
        std::ostringstream out;
        if (solutions.size() == 1) {
            *output = format_critical_solution(solutions.front());
            return true;
        }
        // 多个解：输出列表
        out << "[";
        for (std::size_t i = 0; i < solutions.size(); ++i) {
            if (i != 0) {
                out << ", ";
            }
            out << format_critical_solution(solutions[i]);
        }
        out << "]";
        *output = out.str();
        return true;
    }

    // ==================== 极值点查找 ====================
    if (command == "extrema") {
        if (arguments.size() != 3) {
            throw std::runtime_error("extrema expects expression, a, b");
        }

        FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
        double left = ctx.parse_decimal(arguments[1]);
        double right = ctx.parse_decimal(arguments[2]);

        // 在给定区间内查找极值点
        const std::vector<ExtremumPoint> points = analysis.solve_extrema(left, right);

        if (points.empty()) {
            *output = "No extrema found in the given interval.";
            return true;
        }

        // 格式化输出
        std::ostringstream out;
        auto display_value = [](double value) {
            if (is_integer_double(value, 1e-6)) {
                return format_decimal(static_cast<double>(round_to_long_long(value)));
            }
            return format_decimal(value);
        };
        for (std::size_t i = 0; i < points.size(); ++i) {
            if (i != 0) {
                out << '\n';
            }
            out << (points[i].is_maximum ? "max" : "min")
                << ": x = " << display_value(points[i].x)
                << ", f(x) = " << display_value(points[i].value);
        }
        *output = out.str();
        return true;
    }

    return false;
}

}  // namespace analysis_cmds

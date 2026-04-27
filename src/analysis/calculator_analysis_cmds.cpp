// ============================================================================
// 函数分析命令实现
// ============================================================================

#include "calculator_analysis_cmds.h"

#include "mymath.h"

#include <sstream>
#include <vector>

namespace analysis_cmds {

std::string classify_critical_point(
    const std::vector<std::vector<SymbolicExpression>>& hessian,
    const std::vector<std::string>& variables,
    const std::vector<double>& values,
    const std::function<bool(SymbolicExpression, const std::vector<std::string>&, const std::vector<double>&, double*)>& evaluate_at_point) {

    const std::size_t n = variables.size();
    std::vector<std::vector<double>> evaluated(n, std::vector<double>(n, 0.0));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            if (!evaluate_at_point(hessian[row][col], variables, values, &evaluated[row][col])) {
                return "unclassified";
            }
        }
    }

    if (n == 1) {
        if (evaluated[0][0] > 1e-8) {
            return "local min";
        }
        if (evaluated[0][0] < -1e-8) {
            return "local max";
        }
        return "degenerate";
    }

    if (n == 2) {
        const double det =
            evaluated[0][0] * evaluated[1][1] -
            evaluated[0][1] * evaluated[1][0];
        if (det > 1e-8 && evaluated[0][0] > 1e-8) {
            return "local min";
        }
        if (det > 1e-8 && evaluated[0][0] < -1e-8) {
            return "local max";
        }
        if (det < -1e-8) {
            return "saddle";
        }
        return "degenerate";
    }

    if (n == 3) {
        const double d1 = evaluated[0][0];
        const double d2 =
            evaluated[0][0] * evaluated[1][1] -
            evaluated[0][1] * evaluated[1][0];
        const double d3 =
            evaluated[0][0] *
                (evaluated[1][1] * evaluated[2][2] -
                 evaluated[1][2] * evaluated[2][1]) -
            evaluated[0][1] *
                (evaluated[1][0] * evaluated[2][2] -
                 evaluated[1][2] * evaluated[2][0]) +
            evaluated[0][2] *
                (evaluated[1][0] * evaluated[2][1] -
                 evaluated[1][1] * evaluated[2][0]);
        if (d1 > 1e-8 && d2 > 1e-8 && d3 > 1e-8) {
            return "local min";
        }
        if (d1 < -1e-8 && d2 > 1e-8 && d3 < -1e-8) {
            return "local max";
        }
        if (mymath::abs(d1) <= 1e-8 ||
            mymath::abs(d2) <= 1e-8 ||
            mymath::abs(d3) <= 1e-8) {
            return "degenerate";
        }
        return "saddle";
    }

    return "unclassified";
}

bool is_analysis_command(const std::string& command) {
    return command == "critical" || command == "extrema";
}

bool handle_analysis_command(const AnalysisContext& ctx,
                             const std::string& command,
                             const std::string& inside,
                             std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "critical") {
        if (arguments.empty()) {
            throw std::runtime_error(
                "critical expects a symbolic expression and optional variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        ctx.resolve_symbolic(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            ctx.parse_symbolic_variable_arguments(arguments,
                                                  1,
                                                  expression.identifier_variables());
        const std::vector<SymbolicExpression> gradient =
            expression.gradient(variables);
        const std::vector<std::vector<SymbolicExpression>> hessian =
            expression.hessian(variables);

        // 格式化临界点解
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
            return out.str();
        };

        // 求解临界点（简化版本，实际实现更复杂）
        // 这里返回 false 让 calculator_commands.cpp 处理完整逻辑
        return false;
    }

    if (command == "extrema") {
        if (arguments.size() != 3) {
            throw std::runtime_error("extrema expects expression, a, b");
        }

        FunctionAnalysis analysis = ctx.build_analysis(arguments[0]);
        double left = ctx.parse_decimal(arguments[1]);
        double right = ctx.parse_decimal(arguments[2]);

        const std::vector<ExtremumPoint> points = analysis.solve_extrema(left, right);

        if (points.empty()) {
            *output = "No extrema found in the given interval.";
            return true;
        }

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

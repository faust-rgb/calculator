// ============================================================================
// 函数分析命令辅助函数实现
// ============================================================================

#include "analysis/calculus/analysis_command_helpers.h"
#include "math/numeric/precision/tolerances.h"
#include "math/mymath.h"
#include "symbolic/core/symbolic_expression.h"

#include <algorithm>
#include <stdexcept>

namespace analysis_cmds {


// ============================================================================
// 临界点分类
// ============================================================================

std::string classify_critical_point(
    const std::vector<std::vector<SymbolicExpression>>& hessian,
    const std::vector<std::string>& variables,
    const std::vector<Scalar>& values) {
    std::vector<std::vector<Scalar>> numeric_hessian(variables.size(), std::vector<Scalar>(variables.size(), Scalar(0)));
    for (std::size_t i = 0; i < variables.size(); ++i) {
        for (std::size_t j = 0; j < variables.size(); ++j) {
            SymbolicExpression current = hessian[i][j];
            for (std::size_t k = 0; k < variables.size(); ++k) {
                current = current.substitute(variables[k], SymbolicExpression::number((values[k]))).simplify();
            }
            if (current.node_type() == NodeType::kNumber) numeric_hessian[i][j] = Scalar(current.node_numeric_value());
            else return "unknown";
        }
    }

    if (variables.size() == 1) {
        Scalar d2f = numeric_hessian[0][0];
        if (mymath::isfinite(d2f) && mymath::abs(d2f) < precision::sqrt_epsilon<Scalar>() * Scalar(10)) return "degenerate";
        return d2f > Scalar(0) ? "local min" : "local max";
    }

    if (variables.size() == 2) {
        Scalar A = numeric_hessian[0][0], B = numeric_hessian[0][1], C = numeric_hessian[1][1];
        Scalar D = A * C - B * B;
        if (mymath::isfinite(D) && mymath::abs(D) < precision::sqrt_epsilon<Scalar>() * Scalar(10)) return "degenerate";
        if (D < Scalar(0)) return "saddle point";
        return A > Scalar(0) ? "local min" : "local max";
    }

    return "higher order";
}

// ============================================================================
// 线性系统求解
// ============================================================================

std::vector<Scalar> solve_linear_system(
    std::vector<std::vector<Scalar>> matrix,
    std::vector<Scalar> rhs) {
    const std::size_t n = rhs.size();
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < n; ++row) {
            if (mymath::abs(matrix[row][col]) > mymath::abs(matrix[pivot][col])) {
                pivot = row;
            }
        }
        if (mymath::isfinite(matrix[pivot][col]) && mymath::abs(matrix[pivot][col]) < precision::singular_value_threshold<Scalar>()) {
            throw std::runtime_error("singular critical point system");
        }
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(rhs[pivot], rhs[col]);
        }
        const Scalar divisor = matrix[col][col];
        for (std::size_t c = col; c < n; ++c) matrix[col][c] = matrix[col][c] / divisor;
        rhs[col] = rhs[col] / divisor;
        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;
            const Scalar factor = matrix[row][col];
            if (mymath::isfinite(factor) && mymath::abs(factor) < precision::epsilon<Scalar>() * Scalar(100)) continue;
            for (std::size_t c = col; c < n; ++c) matrix[row][c] = matrix[row][c] - factor * matrix[col][c];
            rhs[row] = rhs[row] - factor * rhs[col];
        }
    }
    return rhs;
}

// ============================================================================
// 无穷大检测
// ============================================================================

bool is_infinity_literal(const std::string& text) {
    std::string value = text;
    // 去除前后空格
    size_t start = value.find_first_not_of(" \t");
    size_t end = value.find_last_not_of(" \t");
    if (start == std::string::npos) return false;
    value = value.substr(start, end - start + 1);

    // 去除符号
    if (!value.empty() && value.front() == '+') value = value.substr(1);
    else if (!value.empty() && value.front() == '-') value = value.substr(1);

    return value == "inf" || value == "infinity" || value == "oo";
}

}  // namespace analysis_cmds

// ============================================================================
// symbolic_system_solver.cpp - 符号方程组解析求解器实现
// ============================================================================

#include "symbolic/solver/symbolic_system_solver.h"
#include "symbolic/algebra/polynomial/symbolic_polynomial.h"
#include "core/services/string_utils.h"
#include <stdexcept>
#include <iostream>

namespace symbolic_solvers {

namespace {

SymbolicExpression parse_single_equation(const std::string& eq_str) {
    std::string trimmed = trim_copy(eq_str);
    std::size_t eq_pos = trimmed.find("==");
    if (eq_pos != std::string::npos) {
        SymbolicExpression lhs = SymbolicExpression::parse(trimmed.substr(0, eq_pos));
        SymbolicExpression rhs = SymbolicExpression::parse(trimmed.substr(eq_pos + 2));
        return (lhs - rhs).simplify();
    }
    eq_pos = trimmed.find('=');
    if (eq_pos != std::string::npos) {
        SymbolicExpression lhs = SymbolicExpression::parse(trimmed.substr(0, eq_pos));
        SymbolicExpression rhs = SymbolicExpression::parse(trimmed.substr(eq_pos + 1));
        return (lhs - rhs).simplify();
    }
    return SymbolicExpression::parse(trimmed).simplify();
}

} // namespace

std::map<std::string, SymbolicExpression> solve_linear_system_symbolic(
    const std::vector<SymbolicExpression>& equations,
    const std::vector<std::string>& vars) {

    const std::size_t n = vars.size();
    if (equations.size() < n) {
        throw std::runtime_error("underdetermined system: fewer equations than variables");
    }

    // 构建系数矩阵 A[n][n] 和常数向量 B[n]
    std::vector<std::vector<SymbolicExpression>> A(n, std::vector<SymbolicExpression>(n, SymbolicExpression::number(0)));
    std::vector<SymbolicExpression> B(n, SymbolicExpression::number(0));

    for (std::size_t i = 0; i < n; ++i) {
        SymbolicExpression eq = equations[i].simplify();
        SymbolicExpression linear_terms_sum = SymbolicExpression::number(0);

        for (std::size_t j = 0; j < n; ++j) {
            // A[i][j] = d(eq)/d(vars[j])
            SymbolicExpression coeff = eq.derivative(vars[j]).simplify();
            A[i][j] = coeff;
            linear_terms_sum = (linear_terms_sum + coeff * SymbolicExpression::variable(vars[j])).simplify();
        }

        // B[i] = linear_terms_sum - eq
        B[i] = (linear_terms_sum - eq).simplify();
    }

    // 高斯-若尔当消元法 (Gauss-Jordan Elimination)
    for (std::size_t k = 0; k < n; ++k) {
        // 寻找主元
        std::size_t pivot = k;
        bool found = false;
        for (std::size_t r = k; r < n; ++r) {
            SymbolicExpression pivot_val = A[r][k].simplify();
            Scalar num;
            if (pivot_val.is_number(&num) && num == Scalar(0)) {
                continue;
            }
            pivot = r;
            found = true;
            break;
        }

        if (!found) {
            throw std::runtime_error("singular matrix or degenerate system in solve_linear_system_symbolic");
        }

        if (pivot != k) {
            std::swap(A[k], A[pivot]);
            std::swap(B[k], B[pivot]);
        }

        // 归一化第 k 行
        SymbolicExpression diag = A[k][k];
        for (std::size_t j = k; j < n; ++j) {
            A[k][j] = (A[k][j] / diag).simplify();
        }
        B[k] = (B[k] / diag).simplify();

        // 消去其他行的第 k 列
        for (std::size_t r = 0; r < n; ++r) {
            if (r != k) {
                SymbolicExpression factor = A[r][k];
                if (factor.to_string() != "0") {
                    for (std::size_t j = k; j < n; ++j) {
                        A[r][j] = (A[r][j] - factor * A[k][j]).simplify();
                    }
                    B[r] = (B[r] - factor * B[k]).simplify();
                }
            }
        }
    }

    std::map<std::string, SymbolicExpression> result;
    for (std::size_t i = 0; i < n; ++i) {
        result[vars[i]] = B[i].simplify();
    }
    return result;
}

std::map<std::string, SymbolicExpression> solve_linear_system_symbolic(
    const std::vector<std::string>& equations,
    const std::vector<std::string>& vars) {

    std::vector<SymbolicExpression> parsed_eqs;
    parsed_eqs.reserve(equations.size());
    for (const auto& s : equations) {
        parsed_eqs.push_back(parse_single_equation(s));
    }
    return solve_linear_system_symbolic(parsed_eqs, vars);
}

} // namespace symbolic_solvers

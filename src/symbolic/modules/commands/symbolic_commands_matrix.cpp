/**
 * @file symbolic_commands_matrix.cpp
 * @brief 符号矩阵命令实现
 *
 * 本文件实现了符号矩阵相关的命令处理：
 * - det: 符号行列式计算（使用 Bareiss 算法）
 * - inverse: 符号矩阵求逆
 * - eigenvalues: 符号特征值计算
 *
 * Bareiss 算法是无分数的高斯消元，适合符号计算。
 */

#include "symbolic/modules/commands/symbolic_commands_internal.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "symbolic/solver/symbolic_solver.h"
#include "symbolic/matrix/symbolic_matrix.h"
#include "core/services/string_utils.h"
#include <vector>
#include <sstream>

namespace symbolic_commands {

namespace {
using namespace symbolic_expression_internal;

// Bareiss algorithm for symbolic determinant - O(n^3) complexity
// Based on fraction-free Gaussian elimination
SymbolicExpression compute_det_bareiss(std::vector<std::vector<SymbolicExpression>> m) {
    std::size_t n = m.size();
    if (n == 0) return SymbolicExpression::number(1.0L);
    if (n == 1) return m[0][0];
    if (n == 2) return (m[0][0] * m[1][1] - m[0][1] * m[1][0]).simplify();

    SymbolicExpression det = SymbolicExpression::number(1.0L);
    int sign = 1;

    for (std::size_t k = 0; k < n; ++k) {
        // Find pivot
        std::size_t pivot_row = k;
        while (pivot_row < n && expr_is_zero(m[pivot_row][k])) {
            ++pivot_row;
        }

        if (pivot_row == n) return SymbolicExpression::number(0.0L);

        // Swap rows if needed
        if (pivot_row != k) {
            std::swap(m[k], m[pivot_row]);
            sign = -sign;
        }

        // Bareiss elimination
        for (std::size_t i = k + 1; i < n; ++i) {
            for (std::size_t j = k + 1; j < n; ++j) {
                // m[i][j] = (m[k][k] * m[i][j] - m[i][k] * m[k][j]) / det
                SymbolicExpression numerator = (m[k][k] * m[i][j] - m[i][k] * m[k][j]).simplify();
                if (k > 0) {
                    m[i][j] = (numerator / det).simplify();
                } else {
                    m[i][j] = numerator;
                }
            }
            m[i][k] = SymbolicExpression::number(0.0L);
        }
        det = m[k][k];
    }

    SymbolicExpression result = m[n-1][n-1];
    if (sign < 0) result = make_negate(result).simplify();
    return result;
}

// Compute nullspace from RREF matrix
std::vector<std::vector<SymbolicExpression>> compute_nullspace(std::vector<std::vector<SymbolicExpression>> rref_mat,
                                                                std::size_t pivot_count) {
    std::size_t m_rows = rref_mat.size();
    if (m_rows == 0) return {};
    std::size_t n_cols = rref_mat[0].size();

    // Identify pivot columns
    std::vector<bool> is_pivot(n_cols, false);
    std::size_t row = 0;
    for (std::size_t col = 0; col < n_cols && row < pivot_count; ++col) {
        if (!expr_is_zero(rref_mat[row][col])) {
            is_pivot[col] = true;
            ++row;
        }
    }

    // Free columns form nullspace basis
    std::vector<std::vector<SymbolicExpression>> nullspace_basis;
    for (std::size_t col = 0; col < n_cols; ++col) {
        if (!is_pivot[col]) {
            // Create basis vector for this free variable
            std::vector<SymbolicExpression> basis_vec(n_cols, SymbolicExpression::number(0.0L));
            basis_vec[col] = SymbolicExpression::number(1.0L);

            // Fill in pivot variable values
            std::size_t pivot_row = 0;
            for (std::size_t j = 0; j < col && pivot_row < pivot_count; ++j) {
                if (is_pivot[j]) {
                    basis_vec[j] = make_negate(rref_mat[pivot_row][col]).simplify();
                    ++pivot_row;
                }
            }
            nullspace_basis.push_back(basis_vec);
        }
    }

    return nullspace_basis;
}

// Compute characteristic polynomial coefficients for eigenvalue computation
// Returns coefficients of det(A - λI) as polynomial in λ
std::vector<SymbolicExpression> characteristic_poly_coeffs(const std::vector<std::vector<SymbolicExpression>>& mat) {
    std::size_t n = mat.size();
    if (n == 0) return {};
    if (n == 1) return {make_negate(mat[0][0]).simplify(), SymbolicExpression::number(1.0L)};

    // Build A - λI matrix
    SymbolicExpression lambda = SymbolicExpression::variable("λ");
    std::vector<std::vector<SymbolicExpression>> char_mat(n, std::vector<SymbolicExpression>(n));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) {
                char_mat[i][j] = (mat[i][j] - lambda).simplify();
            } else {
                char_mat[i][j] = mat[i][j];
            }
        }
    }

    // Compute determinant symbolically - this gives characteristic polynomial
    SymbolicExpression det = compute_det_bareiss(char_mat);

    // Return the determinant expression - caller will handle it
    return {det};
}

}

bool handle_matrix_commands(const SymbolicCommandContext& ctx,
                           const std::string& command,
                           const std::string& /*inside*/,
                           const std::vector<std::string>& arguments,
                           std::string* output) {
    if (command == "det" || command == "inverse" || command == "transpose" || command == "rref" || command == "rank" || command == "nullspace" || command == "eigenvalues") {
        if (arguments.size() != 1) throw std::runtime_error(command + " expects exactly one argument");
        std::string var; SymbolicExpression expr;
        ctx.resolve_symbolic(arguments[0], false, &var, &expr);

        if (!expr.is_tensor() && !expr.is_vector()) throw std::runtime_error(command + " expects a matrix or vector");

        if (command == "transpose") {
            if (expr.is_vector()) {
                auto components = expr.vector_components();
                std::vector<std::vector<SymbolicExpression>> rows;
                for (const auto& comp : components) rows.push_back({comp});
                *output = SymbolicExpression::tensor(rows).simplify().to_string();
            } else {
                auto rows = expr.tensor_rows();
                if (rows.empty()) { *output = expr.to_string(); return true; }
                std::size_t old_rows = rows.size(), old_cols = rows[0].size();
                std::vector<std::vector<SymbolicExpression>> new_rows(old_cols, std::vector<SymbolicExpression>(old_rows));
                for (size_t i = 0; i < old_rows; ++i)
                    for (size_t j = 0; j < old_cols; ++j) new_rows[j][i] = rows[i][j];
                *output = SymbolicExpression::tensor(new_rows).simplify().to_string();
            }
            return true;
        }

        if (!expr.is_tensor()) throw std::runtime_error(command + " expects a matrix");
        auto rows = expr.tensor_rows();
        std::size_t m_rows = rows.size();
        if (m_rows == 0) throw std::runtime_error("empty matrix");
        std::size_t n_cols = rows[0].size();

        if (command == "det") {
            if (m_rows != n_cols) throw std::runtime_error("det expects a square matrix");
            *output = compute_det_bareiss(rows).to_string();
            return true;
        }

        if (command == "inverse") {
            if (m_rows != n_cols) throw std::runtime_error("inverse expects a square matrix");
            SymbolicExpression det = compute_det_bareiss(rows);
            if (expr_is_zero(det)) throw std::runtime_error("matrix is singular");

            // Use augmented matrix [A|I] and compute RREF
            std::vector<std::vector<SymbolicExpression>> aug(m_rows, std::vector<SymbolicExpression>(2 * m_rows));
            for (std::size_t i = 0; i < m_rows; ++i) {
                for (std::size_t j = 0; j < m_rows; ++j) {
                    aug[i][j] = rows[i][j];
                    aug[i][j + m_rows] = (i == j) ? SymbolicExpression::number(1.0L) : SymbolicExpression::number(0.0L);
                }
            }

            // Gaussian elimination on augmented matrix
            std::size_t pivot_row = 0;
            for (std::size_t col = 0; col < m_rows && pivot_row < m_rows; ++col) {
                // Find pivot
                std::size_t best_row = pivot_row;
                for (std::size_t i = pivot_row; i < m_rows; ++i) {
                    if (!expr_is_zero(aug[i][col])) {
                        best_row = i;
                        break;
                    }
                }
                if (expr_is_zero(aug[best_row][col])) continue;

                std::swap(aug[pivot_row], aug[best_row]);
                SymbolicExpression pivot = aug[pivot_row][col];

                // Normalize pivot row
                for (std::size_t j = 0; j < 2 * m_rows; ++j) {
                    aug[pivot_row][j] = (aug[pivot_row][j] / pivot).simplify();
                }

                // Eliminate column
                for (std::size_t i = 0; i < m_rows; ++i) {
                    if (i != pivot_row && !expr_is_zero(aug[i][col])) {
                        SymbolicExpression factor = aug[i][col];
                        for (std::size_t j = 0; j < 2 * m_rows; ++j) {
                            aug[i][j] = (aug[i][j] - factor * aug[pivot_row][j]).simplify();
                        }
                    }
                }
                ++pivot_row;
            }

            // Extract inverse from right half
            std::vector<std::vector<SymbolicExpression>> inv(m_rows, std::vector<SymbolicExpression>(m_rows));
            for (std::size_t i = 0; i < m_rows; ++i) {
                for (std::size_t j = 0; j < m_rows; ++j) {
                    inv[i][j] = aug[i][j + m_rows];
                }
            }
            *output = SymbolicExpression::tensor(inv).to_string();
            return true;
        }

        if (command == "rref" || command == "rank") {
            // Symbolic Gaussian Elimination for RREF
            std::vector<std::vector<SymbolicExpression>> mat = rows;
            size_t pivot_row = 0;
            size_t pivot_col = 0;
            while (pivot_row < m_rows && pivot_col < n_cols) {
                size_t best_row = pivot_row;
                bool found = false;
                for (size_t i = pivot_row; i < m_rows; ++i) {
                    if (!expr_is_zero(mat[i][pivot_col])) {
                        best_row = i;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    pivot_col++;
                    continue;
                }
                std::swap(mat[pivot_row], mat[best_row]);
                SymbolicExpression pivot_val = mat[pivot_row][pivot_col];
                for (size_t j = pivot_col; j < n_cols; ++j) {
                    mat[pivot_row][j] = (mat[pivot_row][j] / pivot_val).simplify();
                }
                for (size_t i = 0; i < m_rows; ++i) {
                    if (i != pivot_row && !expr_is_zero(mat[i][pivot_col])) {
                        SymbolicExpression factor = mat[i][pivot_col];
                        for (size_t j = pivot_col; j < n_cols; ++j) {
                            mat[i][j] = (mat[i][j] - factor * mat[pivot_row][j]).simplify();
                        }
                    }
                }
                pivot_row++;
                pivot_col++;
            }
            if (command == "rank") {
                *output = std::to_string(pivot_row);
            } else {
                *output = SymbolicExpression::tensor(mat).to_string();
            }
            return true;
        }

        if (command == "nullspace") {
            // Compute RREF first
            std::vector<std::vector<SymbolicExpression>> mat = rows;
            size_t pivot_row = 0;
            size_t pivot_col = 0;
            std::vector<std::size_t> pivot_cols;
            while (pivot_row < m_rows && pivot_col < n_cols) {
                size_t best_row = pivot_row;
                bool found = false;
                for (size_t i = pivot_row; i < m_rows; ++i) {
                    if (!expr_is_zero(mat[i][pivot_col])) {
                        best_row = i;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    pivot_col++;
                    continue;
                }
                std::swap(mat[pivot_row], mat[best_row]);
                pivot_cols.push_back(pivot_col);
                SymbolicExpression pivot_val = mat[pivot_row][pivot_col];
                for (size_t j = pivot_col; j < n_cols; ++j) {
                    mat[pivot_row][j] = (mat[pivot_row][j] / pivot_val).simplify();
                }
                for (size_t i = 0; i < m_rows; ++i) {
                    if (i != pivot_row && !expr_is_zero(mat[i][pivot_col])) {
                        SymbolicExpression factor = mat[i][pivot_col];
                        for (size_t j = pivot_col; j < n_cols; ++j) {
                            mat[i][j] = (mat[i][j] - factor * mat[pivot_row][j]).simplify();
                        }
                    }
                }
                pivot_row++;
                pivot_col++;
            }

            // Compute nullspace basis
            auto nullspace_basis = compute_nullspace(mat, pivot_row);

            if (nullspace_basis.empty()) {
                *output = "{}";
            } else {
                std::ostringstream oss;
                oss << "{";
                for (size_t i = 0; i < nullspace_basis.size(); ++i) {
                    if (i > 0) oss << ", ";
                    oss << SymbolicExpression::tensor({nullspace_basis[i]}).to_string();
                }
                oss << "}";
                *output = oss.str();
            }
            return true;
        }

        if (command == "eigenvalues") {
            if (m_rows != n_cols) throw std::runtime_error("eigenvalues expects a square matrix");

            // Compute characteristic polynomial
            auto coeffs = characteristic_poly_coeffs(rows);

            if (coeffs.size() == 1) {
                // Full symbolic characteristic polynomial
                *output = "char_poly(λ) = " + coeffs[0].to_string();
            } else {
                // Multiple coefficients - format as polynomial
                std::ostringstream oss;
                oss << "characteristic polynomial: ";
                for (size_t i = coeffs.size(); i-- > 0; ) {
                    if (i < coeffs.size() - 1) oss << " + ";
                    if (i == 0) oss << coeffs[i].to_string();
                    else if (i == 1) oss << coeffs[i].to_string() << "*λ";
                    else oss << coeffs[i].to_string() << "*λ^" << i;
                }

                // Try to find eigenvalues symbolically
                SymbolicExpression poly = coeffs[0];
                SymbolicExpression lambda = SymbolicExpression::variable("λ");
                for (size_t i = 1; i < coeffs.size(); ++i) {
                    poly = (poly + coeffs[i] * make_power(lambda, SymbolicExpression::number(static_cast<long long>(i)))).simplify();
                }

                symbolic_solver::SymbolicSolver solver;
                symbolic_solver::Solution sol = solver.solve(poly, "λ");
                if (sol.is_complete || !sol.values.empty() || sol.uses_root_of) {
                    oss << "\neigenvalues: " << symbolic_solver::format_solution(sol);
                }

                *output = oss.str();
            }
            return true;
        }
    }

    if (command == "smith" || command == "snf" ||
        command == "hermite" || command == "hnf" ||
        command == "symbolic_lu" || command == "sym_lu" ||
        command == "symbolic_qr" || command == "sym_qr" ||
        command == "symbolic_jordan" || command == "sym_jordan" || command == "jordan" ||
        command == "symbolic_expm" || command == "sym_expm") {
        if (arguments.empty()) throw std::runtime_error(command + " expects matrix argument");
        std::string var; SymbolicExpression expr;
        ctx.resolve_symbolic(arguments[0], false, &var, &expr);
        if (!expr.is_tensor()) throw std::runtime_error(command + " expects a matrix");
        auto rows = expr.tensor_rows();

        if (command == "smith" || command == "snf") {
            std::string pvar = (arguments.size() > 1) ? arguments[1] : "x";
            auto res = symbolic_matrix::smith_normal_form(rows, pvar);
            *output = "S = " + symbolic_matrix::matrix_to_string(res.S) +
                      ", P = " + symbolic_matrix::matrix_to_string(res.P) +
                      ", Q = " + symbolic_matrix::matrix_to_string(res.Q);
            return true;
        }

        if (command == "hermite" || command == "hnf") {
            std::string pvar = (arguments.size() > 1) ? arguments[1] : "x";
            auto res = symbolic_matrix::hermite_normal_form(rows, pvar);
            *output = "H = " + symbolic_matrix::matrix_to_string(res.H) +
                      ", P = " + symbolic_matrix::matrix_to_string(res.P);
            return true;
        }

        if (command == "symbolic_lu" || command == "sym_lu") {
            auto res = symbolic_matrix::symbolic_lu(rows);
            *output = "L = " + symbolic_matrix::matrix_to_string(res.L) +
                      ", U = " + symbolic_matrix::matrix_to_string(res.U) +
                      ", P = " + symbolic_matrix::matrix_to_string(res.P);
            return true;
        }

        if (command == "symbolic_qr" || command == "sym_qr") {
            auto res = symbolic_matrix::symbolic_qr(rows);
            *output = "Q = " + symbolic_matrix::matrix_to_string(res.Q) +
                      ", R = " + symbolic_matrix::matrix_to_string(res.R);
            return true;
        }

        if (command == "symbolic_jordan" || command == "sym_jordan" || command == "jordan") {
            auto res = symbolic_matrix::symbolic_jordan(rows);
            *output = "J = " + symbolic_matrix::matrix_to_string(res.J) +
                      ", P = " + symbolic_matrix::matrix_to_string(res.P);
            return true;
        }

        if (command == "symbolic_expm" || command == "sym_expm") {
            std::string tvar = (arguments.size() > 1) ? arguments[1] : "t";
            auto res = symbolic_matrix::symbolic_matrix_exponential(rows, tvar);
            *output = symbolic_matrix::matrix_to_string(res);
            return true;
        }
    }

    if (command == "jacobian" || command == "hessian") {
        if (command == "jacobian") {
            if (arguments.size() < 1) throw std::runtime_error("jacobian expects expression list");
            auto exprs = ctx.parse_symbolic_expression_list(arguments[0]);
            auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {});
            auto jac = SymbolicExpression::jacobian(exprs, vars);
            *output = symbolic_matrix_to_string(jac);
        } else {
            if (arguments.size() < 1) throw std::runtime_error("hessian expects expression");
            std::string v; SymbolicExpression e; ctx.resolve_symbolic(arguments[0], false, &v, &e);
            auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {v});
            auto hess = e.hessian(vars);
            *output = symbolic_matrix_to_string(hess);
        }
        return true;
    }

    return false;
}

} // namespace symbolic_commands

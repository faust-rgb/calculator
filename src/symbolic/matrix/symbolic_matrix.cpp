/**
 * @file symbolic_matrix.cpp
 * @brief 纯符号线性代数算法实现
 */

#include "symbolic_matrix.h"
#include "symbolic/solver/symbolic_solver.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "math/mymath.h"
#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace symbolic_matrix {

using namespace symbolic_expression_internal;

namespace {
inline SymbolicExpression make_sqrt(const SymbolicExpression& expr) {
    return SymbolicExpression::function("sqrt", expr);
}
inline SymbolicExpression make_exp(const SymbolicExpression& expr) {
    return SymbolicExpression::function("exp", expr);
}
} // namespace

// ============================================================================
// 辅助矩阵工具
// ============================================================================

SymbolicMat identity_matrix(std::size_t n) {
    SymbolicMat I(n, std::vector<SymbolicExpression>(n, SymbolicExpression::number(0.0L)));
    for (std::size_t i = 0; i < n; ++i) {
        I[i][i] = SymbolicExpression::number(1.0L);
    }
    return I;
}

SymbolicMat zero_matrix(std::size_t rows, std::size_t cols) {
    return SymbolicMat(rows, std::vector<SymbolicExpression>(cols, SymbolicExpression::number(0.0L)));
}

SymbolicMat matrix_multiply(const SymbolicMat& A, const SymbolicMat& B) {
    if (A.empty() || B.empty()) return {};
    std::size_t m = A.size();
    std::size_t k = A[0].size();
    std::size_t n = B[0].size();
    if (B.size() != k) {
        throw std::runtime_error("matrix_multiply: incompatible dimensions");
    }

    SymbolicMat C = zero_matrix(m, n);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            SymbolicExpression sum = SymbolicExpression::number(0.0L);
            for (std::size_t p = 0; p < k; ++p) {
                sum = (sum + A[i][p] * B[p][j]).simplify();
            }
            C[i][j] = sum;
        }
    }
    return C;
}

SymbolicMat matrix_add(const SymbolicMat& A, const SymbolicMat& B) {
    std::size_t m = A.size();
    std::size_t n = A[0].size();
    SymbolicMat C(m, std::vector<SymbolicExpression>(n));
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            C[i][j] = (A[i][j] + B[i][j]).simplify();
        }
    }
    return C;
}

SymbolicMat matrix_scale(const SymbolicMat& A, const SymbolicExpression& factor) {
    std::size_t m = A.size();
    std::size_t n = A[0].size();
    SymbolicMat C(m, std::vector<SymbolicExpression>(n));
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            C[i][j] = (A[i][j] * factor).simplify();
        }
    }
    return C;
}

std::string matrix_to_string(const SymbolicMat& A) {
    if (A.empty()) return "[]";
    std::ostringstream oss;
    oss << "[";
    for (std::size_t i = 0; i < A.size(); ++i) {
        if (i > 0) oss << "; ";
        for (std::size_t j = 0; j < A[i].size(); ++j) {
            if (j > 0) oss << ", ";
            oss << A[i][j].to_string();
        }
    }
    oss << "]";
    return oss.str();
}

std::string SymbolicMatrixResult::to_string() const {
    return matrix_to_string(matrix);
}

// 符号行列式 (Bareiss 算法)
static SymbolicExpression compute_bareiss_det(SymbolicMat m) {
    std::size_t n = m.size();
    if (n == 0) return SymbolicExpression::number(1.0L);
    if (n == 1) return m[0][0];
    if (n == 2) return (m[0][0] * m[1][1] - m[0][1] * m[1][0]).simplify();

    SymbolicExpression det = SymbolicExpression::number(1.0L);
    int sign = 1;

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot_row = k;
        while (pivot_row < n && expr_is_zero(m[pivot_row][k])) {
            ++pivot_row;
        }
        if (pivot_row == n) return SymbolicExpression::number(0.0L);

        if (pivot_row != k) {
            std::swap(m[k], m[pivot_row]);
            sign = -sign;
        }

        for (std::size_t i = k + 1; i < n; ++i) {
            for (std::size_t j = k + 1; j < n; ++j) {
                SymbolicExpression num = (m[k][k] * m[i][j] - m[i][k] * m[k][j]).simplify();
                if (k > 0) m[i][j] = (num / det).simplify();
                else m[i][j] = num;
            }
            m[i][k] = SymbolicExpression::number(0.0L);
        }
        det = m[k][k];
    }

    SymbolicExpression result = m[n - 1][n - 1];
    if (sign < 0) result = make_negate(result).simplify();
    return result;
}

// ============================================================================
// 1. Smith 标准型 (Smith Normal Form, SNF)
// ============================================================================

SmithNormalFormResult smith_normal_form(const SymbolicMat& A, const std::string& /*var*/) {
    if (A.empty()) return {};
    std::size_t m = A.size();
    std::size_t n = A[0].size();

    SymbolicMat S = A;
    SymbolicMat P = identity_matrix(m);
    SymbolicMat Q = identity_matrix(n);

    auto row_swap = [&](std::size_t r1, std::size_t r2) {
        if (r1 == r2) return;
        std::swap(S[r1], S[r2]);
        std::swap(P[r1], P[r2]);
    };

    auto col_swap = [&](std::size_t c1, std::size_t c2) {
        if (c1 == c2) return;
        for (std::size_t i = 0; i < m; ++i) std::swap(S[i][c1], S[i][c2]);
        for (std::size_t i = 0; i < n; ++i) std::swap(Q[i][c1], Q[i][c2]);
    };

    auto row_add = [&](std::size_t target, std::size_t src, const SymbolicExpression& factor) {
        for (std::size_t j = 0; j < n; ++j) {
            S[target][j] = (S[target][j] + factor * S[src][j]).simplify();
        }
        for (std::size_t j = 0; j < m; ++j) {
            P[target][j] = (P[target][j] + factor * P[src][j]).simplify();
        }
    };

    auto col_add = [&](std::size_t target, std::size_t src, const SymbolicExpression& factor) {
        for (std::size_t i = 0; i < m; ++i) {
            S[i][target] = (S[i][target] + factor * S[i][src]).simplify();
        }
        for (std::size_t i = 0; i < n; ++i) {
            Q[i][target] = (Q[i][target] + factor * Q[i][src]).simplify();
        }
    };

    std::size_t min_dim = std::min(m, n);

    for (std::size_t k = 0; k < min_dim; ++k) {
        bool progress = true;
        while (progress) {
            progress = false;

            // 寻找主元
            std::size_t p_row = k, p_col = k;
            bool found = false;
            for (std::size_t i = k; i < m && !found; ++i) {
                for (std::size_t j = k; j < n && !found; ++j) {
                    if (!expr_is_zero(S[i][j])) {
                        p_row = i;
                        p_col = j;
                        found = true;
                    }
                }
            }
            if (!found) break;

            row_swap(k, p_row);
            col_swap(k, p_col);

            // 消去第 k 列
            for (std::size_t i = k + 1; i < m; ++i) {
                if (!expr_is_zero(S[i][k])) {
                    SymbolicExpression factor = make_negate(S[i][k] / S[k][k]).simplify();
                    row_add(i, k, factor);
                }
            }

            // 消去第 k 行
            for (std::size_t j = k + 1; j < n; ++j) {
                if (!expr_is_zero(S[k][j])) {
                    SymbolicExpression factor = make_negate(S[k][j] / S[k][k]).simplify();
                    col_add(j, k, factor);
                    progress = true; // 列操作可能破坏之前的零，需复查
                }
            }
        }

        // 归一化主元 (变为首一多项式/常数)
        if (!expr_is_zero(S[k][k])) {
            SymbolicExpression leading = S[k][k];
            for (std::size_t j = 0; j < n; ++j) {
                S[k][j] = (S[k][j] / leading).simplify();
            }
            for (std::size_t j = 0; j < m; ++j) {
                P[k][j] = (P[k][j] / leading).simplify();
            }
        }
    }

    return {S, P, Q};
}

// ============================================================================
// 2. Hermite 标准型 (Hermite Normal Form, HNF)
// ============================================================================

HermiteNormalFormResult hermite_normal_form(const SymbolicMat& A, const std::string& /*var*/) {
    if (A.empty()) return {};
    std::size_t m = A.size();
    std::size_t n = A[0].size();

    SymbolicMat H = A;
    SymbolicMat P = identity_matrix(m);

    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < n && pivot_row < m; ++col) {
        // 找到该列非零元
        std::size_t best_row = pivot_row;
        bool found = false;
        for (std::size_t i = pivot_row; i < m; ++i) {
            if (!expr_is_zero(H[i][col])) {
                best_row = i;
                found = true;
                break;
            }
        }
        if (!found) continue;

        // 行交换
        if (best_row != pivot_row) {
            std::swap(H[pivot_row], H[best_row]);
            std::swap(P[pivot_row], P[best_row]);
        }

        // 归一化主元行
        SymbolicExpression pivot_val = H[pivot_row][col];
        for (std::size_t j = 0; j < n; ++j) {
            H[pivot_row][j] = (H[pivot_row][j] / pivot_val).simplify();
        }
        for (std::size_t j = 0; j < m; ++j) {
            P[pivot_row][j] = (P[pivot_row][j] / pivot_val).simplify();
        }

        // 消去其它所有行的该列元素（包括上方和下方）
        for (std::size_t i = 0; i < m; ++i) {
            if (i != pivot_row && !expr_is_zero(H[i][col])) {
                SymbolicExpression factor = H[i][col];
                for (std::size_t j = 0; j < n; ++j) {
                    H[i][j] = (H[i][j] - factor * H[pivot_row][j]).simplify();
                }
                for (std::size_t j = 0; j < m; ++j) {
                    P[i][j] = (P[i][j] - factor * P[pivot_row][j]).simplify();
                }
            }
        }

        ++pivot_row;
    }

    return {H, P};
}

// ============================================================================
// 3. 符号特征多项式 (Symbolic Charpoly)
// ============================================================================

SymbolicExpression symbolic_charpoly(const SymbolicMat& A, const std::string& lambda_var) {
    std::size_t n = A.size();
    if (n == 0) return SymbolicExpression::number(1.0L);
    if (n != A[0].size()) {
        throw std::runtime_error("symbolic_charpoly requires a square matrix");
    }

    SymbolicExpression lambda = SymbolicExpression::variable(lambda_var);
    SymbolicMat char_mat = A;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) {
                char_mat[i][j] = (lambda - A[i][j]).simplify();
            } else {
                char_mat[i][j] = make_negate(A[i][j]).simplify();
            }
        }
    }

    return compute_bareiss_det(char_mat).simplify();
}

// ============================================================================
// 4. 符号若尔当标准型 (Jordan Canonical Form)
// ============================================================================

SymbolicJordanResult symbolic_jordan(const SymbolicMat& A) {
    std::size_t n = A.size();
    if (n == 0 || n != A[0].size()) {
        throw std::runtime_error("symbolic_jordan requires a square matrix");
    }

    if (n == 1) {
        return {A, identity_matrix(1), identity_matrix(1)};
    }

    if (n == 2) {
        // 2x2 特征多项式求根
        SymbolicExpression tr = (A[0][0] + A[1][1]).simplify();
        SymbolicExpression det = (A[0][0] * A[1][1] - A[0][1] * A[1][0]).simplify();
        SymbolicExpression disc = (tr * tr - SymbolicExpression::number(4.0L) * det).simplify();

        SymbolicExpression sqrt_disc = make_sqrt(disc).simplify();
        SymbolicExpression lam1 = ((tr + sqrt_disc) / SymbolicExpression::number(2.0L)).simplify();
        SymbolicExpression lam2 = ((tr - sqrt_disc) / SymbolicExpression::number(2.0L)).simplify();

        if (expr_is_zero(disc)) {
            // 重特征值 lambda
            SymbolicExpression lam = (tr / SymbolicExpression::number(2.0L)).simplify();
            SymbolicMat A_minus_lam = A;
            A_minus_lam[0][0] = (A[0][0] - lam).simplify();
            A_minus_lam[1][1] = (A[1][1] - lam).simplify();

            if (expr_is_zero(A_minus_lam[0][0]) && expr_is_zero(A_minus_lam[0][1]) &&
                expr_is_zero(A_minus_lam[1][0]) && expr_is_zero(A_minus_lam[1][1])) {
                // 对角阵
                SymbolicMat J = zero_matrix(2, 2);
                J[0][0] = lam; J[1][1] = lam;
                return {J, identity_matrix(2), identity_matrix(2)};
            }

            // 单个 2x2 若尔当块 [lam, 1; 0, lam]
            SymbolicMat J = zero_matrix(2, 2);
            J[0][0] = lam; J[0][1] = SymbolicExpression::number(1.0L);
            J[1][1] = lam;

            // 构造广义特征向量
            SymbolicMat P = identity_matrix(2);
            if (!expr_is_zero(A[0][1])) {
                P[0][0] = A[0][1];
                P[1][0] = (lam - A[0][0]).simplify();
                P[0][1] = SymbolicExpression::number(0.0L);
                P[1][1] = SymbolicExpression::number(1.0L);
            } else if (!expr_is_zero(A[1][0])) {
                P[0][0] = (lam - A[1][1]).simplify();
                P[1][0] = A[1][0];
                P[0][1] = SymbolicExpression::number(1.0L);
                P[1][1] = SymbolicExpression::number(0.0L);
            }
            return {J, P, identity_matrix(2)};
        }

        // 两个互异特征值，对角形
        SymbolicMat J = zero_matrix(2, 2);
        J[0][0] = lam1;
        J[1][1] = lam2;

        SymbolicMat P = zero_matrix(2, 2);
        if (!expr_is_zero(A[0][1])) {
            P[0][0] = A[0][1]; P[1][0] = (lam1 - A[0][0]).simplify();
            P[0][1] = A[0][1]; P[1][1] = (lam2 - A[0][0]).simplify();
        } else if (!expr_is_zero(A[1][0])) {
            P[0][0] = (lam1 - A[1][1]).simplify(); P[1][0] = A[1][0];
            P[0][1] = (lam2 - A[1][1]).simplify(); P[1][1] = A[1][0];
        } else {
            P = identity_matrix(2);
        }

        return {J, P, identity_matrix(2)};
    }

    // 默认回退特征系统
    SymbolicExpression poly = symbolic_charpoly(A, "lambda");
    symbolic_solver::SymbolicSolver solver;
    auto sol = solver.solve(poly, "lambda");

    SymbolicMat J = zero_matrix(n, n);
    for (std::size_t i = 0; i < std::min(n, sol.values.size()); ++i) {
        J[i][i] = sol.values[i];
    }
    return {J, identity_matrix(n), identity_matrix(n)};
}

// ============================================================================
// 5. 符号矩阵指数 exp(A * t) 解析解 (Putzer 算法)
// ============================================================================

SymbolicMat symbolic_matrix_exponential(const SymbolicMat& A, const std::string& t_var) {
    std::size_t n = A.size();
    if (n == 0 || n != A[0].size()) {
        throw std::runtime_error("symbolic_matrix_exponential requires a square matrix");
    }

    SymbolicExpression t = SymbolicExpression::variable(t_var);

    if (n == 1) {
        SymbolicMat res = zero_matrix(1, 1);
        res[0][0] = make_exp(A[0][0] * t).simplify();
        return res;
    }

    // 针对 2x2 矩阵的精确解析式
    if (n == 2) {
        SymbolicExpression tr = (A[0][0] + A[1][1]).simplify();
        SymbolicExpression det = (A[0][0] * A[1][1] - A[0][1] * A[1][0]).simplify();
        SymbolicExpression disc = (tr * tr - SymbolicExpression::number(4.0L) * det).simplify();

        SymbolicExpression half_tr = (tr / SymbolicExpression::number(2.0L)).simplify();
        SymbolicExpression exp_half_tr = make_exp(half_tr * t).simplify();

        if (expr_is_zero(disc)) {
            // 重根 lambda = tr/2
            // exp(At) = e^(lambda*t) * [I + (A - lambda*I)*t]
            SymbolicMat I = identity_matrix(2);
            SymbolicMat A_minus_lam = A;
            A_minus_lam[0][0] = (A[0][0] - half_tr).simplify();
            A_minus_lam[1][1] = (A[1][1] - half_tr).simplify();

            SymbolicMat term = matrix_add(I, matrix_scale(A_minus_lam, t));
            return matrix_scale(term, exp_half_tr);
        }

        // 互异特征值或复共轭特征值
        SymbolicExpression delta = (make_sqrt(disc) / SymbolicExpression::number(2.0L)).simplify();

        // r1(t) = e^(lam1 * t), r2(t) = (e^(lam1 * t) - e^(lam2 * t)) / (lam1 - lam2)
        // Putzer: exp(At) = e^(lam1*t)*I + (e^(lam2*t) - e^(lam1*t))/(lam2 - lam1) * (A - lam1*I)
        SymbolicExpression lam1 = (half_tr + delta).simplify();
        SymbolicExpression lam2 = (half_tr - delta).simplify();

        SymbolicExpression exp1 = make_exp(lam1 * t).simplify();
        SymbolicExpression exp2 = make_exp(lam2 * t).simplify();

        SymbolicExpression coeff2 = ((exp2 - exp1) / (lam2 - lam1)).simplify();

        SymbolicMat I = identity_matrix(2);
        SymbolicMat A_minus_lam1 = A;
        A_minus_lam1[0][0] = (A[0][0] - lam1).simplify();
        A_minus_lam1[1][1] = (A[1][1] - lam1).simplify();

        SymbolicMat res = matrix_add(matrix_scale(I, exp1), matrix_scale(A_minus_lam1, coeff2));
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                res[i][j] = res[i][j].simplify();
            }
        }
        return res;
    }

    // 针对任意维度的 Putzer's Algorithm
    // 特征值列表
    SymbolicExpression poly = symbolic_charpoly(A, "lambda");
    symbolic_solver::SymbolicSolver solver;
    auto sol = solver.solve(poly, "lambda");

    std::vector<SymbolicExpression> lams = sol.values;
    while (lams.size() < n) lams.push_back(SymbolicExpression::number(0.0L));

    std::vector<SymbolicMat> P(n + 1);
    P[0] = identity_matrix(n);

    for (std::size_t k = 1; k < n; ++k) {
        SymbolicMat A_minus_lam = A;
        for (std::size_t i = 0; i < n; ++i) {
            A_minus_lam[i][i] = (A[i][i] - lams[k - 1]).simplify();
        }
        P[k] = matrix_multiply(A_minus_lam, P[k - 1]);
    }

    // r_k(t) via Putzer's formula
    std::vector<SymbolicExpression> r(n);
    for (std::size_t k = 0; k < n; ++k) {
        // Check if all eigenvalues up to k are identical
        bool all_same = true;
        for (std::size_t j = 1; j <= k; ++j) {
            if (!expressions_match(lams[j], lams[0])) {
                all_same = false;
                break;
            }
        }

        if (all_same) {
            // r_{k+1}(t) = t^k / k! * e^(lambda*t)
            SymbolicExpression exp_lam = make_exp(lams[0] * t).simplify();
            if (k == 0) {
                r[0] = exp_lam;
            } else {
                SymbolicExpression t_pow = make_power(t, SymbolicExpression::number(Scalar(static_cast<long long>(k))));
                SymbolicExpression fact = SymbolicExpression::number(Scalar(static_cast<long long>(factorial_double(static_cast<int>(k)))));
                r[k] = ((t_pow / fact) * exp_lam).simplify();
            }
        } else {
            // General distinct eigenvalues formula:
            // r_{k+1}(t) = sum_{j=0}^k [ e^(lam_j * t) / prod_{m != j, m=0}^k (lam_j - lam_m) ]
            SymbolicExpression sum_r = SymbolicExpression::number(Scalar(0.0L));
            for (std::size_t j = 0; j <= k; ++j) {
                SymbolicExpression denom = SymbolicExpression::number(Scalar(1.0L));
                for (std::size_t m = 0; m <= k; ++m) {
                    if (m != j) {
                        SymbolicExpression diff = (lams[j] - lams[m]).simplify();
                        if (!expr_is_zero(diff)) {
                            denom = (denom * diff).simplify();
                        }
                    }
                }
                SymbolicExpression term = (make_exp(lams[j] * t) / denom).simplify();
                sum_r = (sum_r + term).simplify();
            }
            r[k] = sum_r;
        }
    }

    SymbolicMat exp_At = zero_matrix(n, n);
    for (std::size_t k = 0; k < n; ++k) {
        exp_At = matrix_add(exp_At, matrix_scale(P[k], r[k]));
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            exp_At[i][j] = exp_At[i][j].simplify();
        }
    }
    return exp_At;
}

// ============================================================================
// 6. 符号矩阵分解 (Symbolic LU & QR)
// ============================================================================

SymbolicLUResult symbolic_lu(const SymbolicMat& A) {
    std::size_t n = A.size();
    if (n == 0 || n != A[0].size()) {
        throw std::runtime_error("symbolic_lu requires a square matrix");
    }

    SymbolicMat L = identity_matrix(n);
    SymbolicMat U = zero_matrix(n, n);
    SymbolicMat P = identity_matrix(n);
    SymbolicMat A_work = A;

    for (std::size_t k = 0; k < n; ++k) {
        // 主元选取
        std::size_t pivot = k;
        while (pivot < n && expr_is_zero(A_work[pivot][k])) {
            ++pivot;
        }
        if (pivot < n && pivot != k) {
            std::swap(A_work[k], A_work[pivot]);
            std::swap(P[k], P[pivot]);
            for (std::size_t j = 0; j < k; ++j) {
                std::swap(L[k][j], L[pivot][j]);
            }
        }

        // 计算 U 的第 k 行
        for (std::size_t j = k; j < n; ++j) {
            SymbolicExpression sum = SymbolicExpression::number(0.0L);
            for (std::size_t m = 0; m < k; ++m) {
                sum = (sum + L[k][m] * U[m][j]).simplify();
            }
            U[k][j] = (A_work[k][j] - sum).simplify();
        }

        // 计算 L 的第 k 列
        if (!expr_is_zero(U[k][k])) {
            for (std::size_t i = k + 1; i < n; ++i) {
                SymbolicExpression sum = SymbolicExpression::number(0.0L);
                for (std::size_t m = 0; m < k; ++m) {
                    sum = (sum + L[i][m] * U[m][k]).simplify();
                }
                L[i][k] = ((A_work[i][k] - sum) / U[k][k]).simplify();
            }
        }
    }

    return {L, U, P};
}

SymbolicQRResult symbolic_qr(const SymbolicMat& A) {
    std::size_t m = A.size();
    if (m == 0) return {};
    std::size_t n = A[0].size();

    SymbolicMat Q = zero_matrix(m, n);
    SymbolicMat R = zero_matrix(n, n);

    // 提取列向量
    std::vector<std::vector<SymbolicExpression>> a_cols(n, std::vector<SymbolicExpression>(m));
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < m; ++i) {
            a_cols[j][i] = A[i][j];
        }
    }

    std::vector<std::vector<SymbolicExpression>> q_cols(n, std::vector<SymbolicExpression>(m));

    auto dot_prod = [&](const std::vector<SymbolicExpression>& u, const std::vector<SymbolicExpression>& v) {
        SymbolicExpression sum = SymbolicExpression::number(0.0L);
        for (std::size_t i = 0; i < m; ++i) {
            sum = (sum + u[i] * v[i]).simplify();
        }
        return sum;
    };

    // Modified Gram-Schmidt
    for (std::size_t j = 0; j < n; ++j) {
        std::vector<SymbolicExpression> v = a_cols[j];
        for (std::size_t i = 0; i < j; ++i) {
            R[i][j] = dot_prod(q_cols[i], v);
            for (std::size_t p = 0; p < m; ++p) {
                v[p] = (v[p] - R[i][j] * q_cols[i][p]).simplify();
            }
        }

        SymbolicExpression norm_sq = dot_prod(v, v);
        R[j][j] = make_sqrt(norm_sq).simplify();

        if (!expr_is_zero(R[j][j])) {
            for (std::size_t p = 0; p < m; ++p) {
                q_cols[j][p] = (v[p] / R[j][j]).simplify();
            }
        } else {
            for (std::size_t p = 0; p < m; ++p) {
                q_cols[j][p] = SymbolicExpression::number(0.0L);
            }
        }
    }

    // 组装 Q
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < m; ++i) {
            Q[i][j] = q_cols[j][i];
        }
    }

    return {Q, R};
}

} // namespace symbolic_matrix

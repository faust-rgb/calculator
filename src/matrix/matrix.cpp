/**
 * @file matrix.cpp
 * @brief 矩阵核心类实现
 *
 * 本文件包含 Matrix 类的核心方法：构造函数、索引访问、大小调整等。
 * 其他矩阵运算、线性代数、统计等函数已移至独立的模块文件中。
 */

#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include <stdexcept>
#include <utility>
#include <vector>

namespace matrix {

using namespace internal;

std::pair<Matrix, Matrix> qr_decompose(const Matrix& matrix);

Matrix::Matrix(std::size_t row_count, std::size_t col_count, double fill_value)
    : rows(row_count),
      cols(col_count),
      data(row_count * col_count, fill_value) {}

Matrix Matrix::vector(const std::vector<double>& values) {
    Matrix matrix(1, values.size(), 0.0);
    matrix.data = values;
    return matrix;
}

Matrix Matrix::zero(std::size_t row_count, std::size_t col_count) {
    return Matrix(row_count, col_count, 0.0);
}

Matrix Matrix::identity(std::size_t size) {
    Matrix matrix(size, size, 0.0);
    for (std::size_t i = 0; i < size; ++i) {
        matrix.at(i, i) = 1.0;
    }
    return matrix;
}

bool Matrix::is_vector() const {
    return rows == 1 || cols == 1;
}

bool Matrix::is_square() const {
    return rows == cols;
}

double& Matrix::at(std::size_t row, std::size_t col) {
    if (row >= rows || col >= cols) {
        throw std::out_of_range("matrix index out of range");
    }
    return data[row * cols + col];
}

double Matrix::at(std::size_t row, std::size_t col) const {
    if (row >= rows || col >= cols) {
        throw std::out_of_range("matrix index out of range");
    }
    return data[row * cols + col];
}

Matrix& Matrix::operator+=(const Matrix& rhs) {
    if (rows != rhs.rows || cols != rhs.cols) {
        throw std::runtime_error("matrix addition requires same shape");
    }
    for (std::size_t i = 0; i < data.size(); ++i) {
        data[i] += rhs.data[i];
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs) {
    if (rows != rhs.rows || cols != rhs.cols) {
        throw std::runtime_error("matrix subtraction requires same shape");
    }
    for (std::size_t i = 0; i < data.size(); ++i) {
        data[i] -= rhs.data[i];
    }
    return *this;
}

Matrix& Matrix::operator*=(double scalar) {
    for (double& val : data) {
        val *= scalar;
    }
    return *this;
}

Matrix& Matrix::operator/=(double scalar) {
    if (mymath::is_near_zero(scalar)) {
        throw std::runtime_error("division by zero");
    }
    const double inv_scalar = 1.0 / scalar;
    for (double& val : data) {
        val *= inv_scalar;
    }
    return *this;
}

// Value 辅助构造函数
Value Value::from_scalar(double scalar_value) {
    Value value;
    value.is_matrix = false;
    value.scalar = scalar_value;
    return value;
}

Value Value::from_complex(double real, double imag) {
    return from_complex(ComplexNumber{real, imag});
}

Value Value::from_complex(ComplexNumber complex_value) {
    Value value;
    value.is_complex = true;
    value.complex = internal::normalize_complex(complex_value);
    return value;
}

Value Value::from_matrix(const Matrix& matrix_value) {
    Value value;
    value.is_matrix = true;
    value.matrix = matrix_value;
    return value;
}

Value Value::from_matrix(Matrix&& matrix_value) {
    Value value;
    value.is_matrix = true;
    value.matrix = std::move(matrix_value);
    return value;
}

namespace internal {

std::vector<double> matrix_column(const std::vector<std::vector<double>>& matrix, std::size_t col) {
    if (matrix.empty()) return {};
    std::vector<double> result(matrix.size());
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        if (col < matrix[i].size()) {
            result[i] = matrix[i][col];
        }
    }
    return result;
}

// ==================== 特征值计算（简化实现） ====================

EigenResult eigenvalues_with_vectors(const Matrix& matrix) {
    if (!matrix.is_square()) {
        throw std::runtime_error("eigenvalues requires a square matrix");
    }

    EigenResult result;
    const std::size_t n = matrix.rows;

    if (n == 0) return result;

    // 检测矩阵是否对称
    bool is_sym = true;
    for (std::size_t i = 0; i < n && is_sym; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (mymath::abs(matrix.at(i, j) - matrix.at(j, i)) > 1e-10) {
                is_sym = false;
                break;
            }
        }
    }

    Matrix A = matrix;
    Matrix V = Matrix::identity(n);

    if (is_sym && n > 1) {
        // 对称矩阵使用 Jacobi 旋转法，非常稳定且能得到精确的正交特征向量
        const double tolerance = 1e-12;
        const int max_sweeps = 100;

        for (int sweep = 0; sweep < max_sweeps; ++sweep) {
            double max_off_diag = 0.0;
            for (std::size_t p = 0; p < n - 1; ++p) {
                for (std::size_t q = p + 1; q < n; ++q) {
                    double off_diag = mymath::abs(A.at(p, q));
                    if (off_diag > max_off_diag) max_off_diag = off_diag;
                }
            }

            if (max_off_diag < tolerance) break;

            for (std::size_t p = 0; p < n - 1; ++p) {
                for (std::size_t q = p + 1; q < n; ++q) {
                    double app = A.at(p, p);
                    double aqq = A.at(q, q);
                    double apq = A.at(p, q);

                    if (mymath::abs(apq) > tolerance * 1e-3) {
                        double theta = (aqq - app) / (2.0 * apq);
                        double t;
                        if (mymath::abs(theta) > 1e15) {
                            t = 1.0 / (2.0 * theta);
                        } else {
                            double sign_theta = theta >= 0.0 ? 1.0 : -1.0;
                            t = sign_theta / (mymath::abs(theta) + mymath::sqrt(1.0 + theta * theta));
                        }
                        double c = 1.0 / mymath::sqrt(1.0 + t * t);
                        double s = c * t;

                        // 更新 A 矩阵（只更新受旋转影响的行和列）
                        A.at(p, p) = app - t * apq;
                        A.at(q, q) = aqq + t * apq;
                        A.at(p, q) = 0.0;
                        A.at(q, p) = 0.0;

                        for (std::size_t j = 0; j < n; ++j) {
                            if (j != p && j != q) {
                                double apj = A.at(p, j);
                                double aqj = A.at(q, j);
                                A.at(p, j) = A.at(j, p) = c * apj - s * aqj;
                                A.at(q, j) = A.at(j, q) = s * apj + c * aqj;
                            }
                        }

                        // 更新特征向量矩阵 V
                        for (std::size_t i = 0; i < n; ++i) {
                            double vip = V.at(i, p);
                            double viq = V.at(i, q);
                            V.at(i, p) = c * vip - s * viq;
                            V.at(i, q) = s * vip + c * viq;
                        }
                    }
                }
            }
        }
    } else if (n > 1) {
        // 非对称矩阵或小矩阵使用通用 QR 迭代
        const int max_iter = n <= 3 ? 200 : 500;
        for (int iter = 0; iter < max_iter; ++iter) {
            auto [Q, R] = matrix::qr_decompose(A);
            A = multiply(R, Q);
            V = multiply(V, Q);

            // 检查收敛：下三角部分是否接近于零（Schur 分解）
            double max_lower = 0.0;
            for (std::size_t i = 1; i < n; ++i) {
                for (std::size_t j = 0; j < i; ++j) {
                    double val = mymath::abs(A.at(i, j));
                    if (val > max_lower) max_lower = val;
                }
            }
            if (max_lower < 1e-10) break;
        }
    }

    // 提取结果并按特征值降序排序
    std::vector<std::pair<double, std::vector<double>>> sorted_results(n);
    for (std::size_t i = 0; i < n; ++i) {
        sorted_results[i] = {A.at(i, i), matrix_column(V, i)};
    }
    std::sort(sorted_results.begin(), sorted_results.end(), [](const auto& a, const auto& b) {
        return a.first > b.first;
    });

    result.values.reserve(n);
    result.vectors.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        result.values.push_back(sorted_results[i].first);
        result.vectors.push_back(sorted_results[i].second);
    }

    return result;
}

}  // namespace internal

}  // namespace matrix

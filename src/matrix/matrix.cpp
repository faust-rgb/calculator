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

    // 对于小矩阵，使用幂迭代法
    if (n <= 3) {
        // 简化实现：使用 QR 算法的简化版本
        Matrix A = matrix;
        Matrix Q_total = Matrix::identity(n);

        for (int iter = 0; iter < 100; ++iter) {
            auto [Q, R] = matrix::qr_decompose(A);
            A = multiply(R, Q);
            Q_total = multiply(Q_total, Q);
        }

        // 提取特征值（对角元素）
        result.values.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            result.values.push_back(A.at(i, i));
        }

        // 特征向量是 Q_total 的列
        result.vectors.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            result.vectors.push_back(matrix_column(Q_total, i));
        }
    } else {
        // 对于大矩阵，返回占位符
        result.values.assign(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            result.vectors.push_back(standard_basis_vector(n, i));
        }
    }

    return result;
}

}  // namespace internal

}  // namespace matrix

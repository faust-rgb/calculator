#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include <stdexcept>

namespace matrix {

using namespace internal;

Matrix add(const Matrix& lhs, const Matrix& rhs) {
    Matrix result = lhs;
    result += rhs;
    return result;
}

Matrix add(Matrix&& lhs, const Matrix& rhs) {
    lhs += rhs;
    return std::move(lhs);
}

Matrix add(const Matrix& lhs, Matrix&& rhs) {
    rhs += lhs;
    return std::move(rhs);
}

Matrix add(const Matrix& lhs, double scalar) {
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] + scalar;
    }
    return result;
}

Matrix add(Matrix&& lhs, double scalar) {
    for (double& val : lhs.data) {
        val += scalar;
    }
    return std::move(lhs);
}

Matrix subtract(const Matrix& lhs, const Matrix& rhs) {
    Matrix result = lhs;
    result -= rhs;
    return result;
}

Matrix subtract(Matrix&& lhs, const Matrix& rhs) {
    lhs -= rhs;
    return std::move(lhs);
}

Matrix subtract(const Matrix& lhs, Matrix&& rhs) {
    for (std::size_t i = 0; i < rhs.data.size(); ++i) {
        rhs.data[i] = lhs.data[i] - rhs.data[i];
    }
    return std::move(rhs);
}

Matrix subtract(const Matrix& lhs, double scalar) {
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] - scalar;
    }
    return result;
}

Matrix subtract(Matrix&& lhs, double scalar) {
    for (double& val : lhs.data) {
        val -= scalar;
    }
    return std::move(lhs);
}

Matrix multiply(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.cols != rhs.rows) {
        throw std::runtime_error("matrix multiplication requires lhs.cols == rhs.rows");
    }

    // Optimization: check if one is identity
    if (lhs.is_square() && lhs.rows == rhs.rows) {
        bool is_lhs_identity = true;
        for (std::size_t i = 0; i < lhs.rows; ++i) {
            for (std::size_t j = 0; j < lhs.cols; ++j) {
                const double val = lhs.at(i, j);
                if (i == j) {
                    if (!mymath::is_near_zero(val - 1.0)) { is_lhs_identity = false; break; }
                } else {
                    if (!mymath::is_near_zero(val)) { is_lhs_identity = false; break; }
                }
            }
            if (!is_lhs_identity) break;
        }
        if (is_lhs_identity) return rhs;
    }
    if (rhs.is_square() && rhs.rows == lhs.cols) {
        bool is_rhs_identity = true;
        for (std::size_t i = 0; i < rhs.rows; ++i) {
            for (std::size_t j = 0; j < rhs.cols; ++j) {
                const double val = rhs.at(i, j);
                if (i == j) {
                    if (!mymath::is_near_zero(val - 1.0)) { is_rhs_identity = false; break; }
                } else {
                    if (!mymath::is_near_zero(val)) { is_rhs_identity = false; break; }
                }
            }
            if (!is_rhs_identity) break;
        }
        if (is_rhs_identity) return lhs;
    }

    Matrix result(lhs.rows, rhs.cols, 0.0);
    std::vector<long double> sums(lhs.rows * rhs.cols, 0.0L);
    
    for (std::size_t i = 0; i < lhs.rows; ++i) {
        for (std::size_t k = 0; k < lhs.cols; ++k) {
            const long double lhs_val = static_cast<long double>(lhs.at(i, k));
            if (mymath::is_near_zero(static_cast<double>(lhs_val))) continue;
            
            const std::size_t row_offset = i * rhs.cols;
            const std::size_t rhs_offset = k * rhs.cols;
            for (std::size_t j = 0; j < rhs.cols; ++j) {
                sums[row_offset + j] += lhs_val * static_cast<long double>(rhs.data[rhs_offset + j]);
            }
        }
    }

    for (std::size_t i = 0; i < sums.size(); ++i) {
        result.data[i] = static_cast<double>(sums[i]);
    }
    return result;
}

Matrix multiply(const Matrix& lhs, double scalar) {
    Matrix result = lhs;
    result *= scalar;
    return result;
}

Matrix multiply(Matrix&& lhs, double scalar) {
    lhs *= scalar;
    return std::move(lhs);
}

Matrix divide(const Matrix& lhs, double scalar) {
    Matrix result = lhs;
    result /= scalar;
    return result;
}

Matrix divide(Matrix&& lhs, double scalar) {
    lhs /= scalar;
    return std::move(lhs);
}

double dot(const Matrix& lhs, const Matrix& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "dot");
    const std::size_t rhs_size = vector_length(rhs, "dot");
    if (lhs_size != rhs_size) {
        throw std::runtime_error("dot requires vectors of the same length");
    }

    long double sum = 0.0L;
    long double compensation = 0.0L;
    for (std::size_t i = 0; i < lhs_size; ++i) {
        const long double term =
            static_cast<long double>(vector_entry(lhs, i)) *
            static_cast<long double>(vector_entry(rhs, i));
        const long double adjusted = term - compensation;
        const long double next = sum + adjusted;
        compensation = (next - sum) - adjusted;
        sum = next;
    }
    return static_cast<double>(sum);
}

Matrix outer(const Matrix& lhs, const Matrix& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "outer");
    const std::size_t rhs_size = vector_length(rhs, "outer");

    Matrix result(lhs_size, rhs_size, 0.0);
    for (std::size_t row = 0; row < lhs_size; ++row) {
        for (std::size_t col = 0; col < rhs_size; ++col) {
            result.at(row, col) = vector_entry(lhs, row) * vector_entry(rhs, col);
        }
    }
    return result;
}

Matrix cross(const Matrix& lhs, const Matrix& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "cross");
    const std::size_t rhs_size = vector_length(rhs, "cross");

    if (lhs_size != 3 || rhs_size != 3) {
        throw std::runtime_error("cross product is only defined for 3D vectors");
    }

    const double x1 = vector_entry(lhs, 0);
    const double y1 = vector_entry(lhs, 1);
    const double z1 = vector_entry(lhs, 2);
    const double x2 = vector_entry(rhs, 0);
    const double y2 = vector_entry(rhs, 1);
    const double z2 = vector_entry(rhs, 2);

    std::vector<double> result_data = {
        y1 * z2 - z1 * y2,
        z1 * x2 - x1 * z2,
        x1 * y2 - y1 * x2
    };

    if (lhs.rows == 1) {
        Matrix res(1, 3);
        res.data = result_data;
        return res;
    } else {
        return Matrix::vector(result_data);
    }
}

Matrix project(const Matrix& lhs, const Matrix& rhs) {
    const double rhs_norm_sq = dot(rhs, rhs);
    if (mymath::is_near_zero(rhs_norm_sq)) {
        throw std::runtime_error("cannot project onto a zero vector");
    }
    return multiply(rhs, dot(lhs, rhs) / rhs_norm_sq);
}

Matrix normalize(const Matrix& matrix) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("normalize only works on vectors");
    }
    const double v_norm = norm(matrix);
    if (mymath::is_near_zero(v_norm)) {
        throw std::runtime_error("cannot normalize a zero vector");
    }
    return divide(matrix, v_norm);
}

Matrix kronecker(const Matrix& lhs, const Matrix& rhs) {
    Matrix result(lhs.rows * rhs.rows, lhs.cols * rhs.cols, 0.0);
    for (std::size_t row = 0; row < lhs.rows; ++row) {
        for (std::size_t col = 0; col < lhs.cols; ++col) {
            const double scale = lhs.at(row, col);
            for (std::size_t rhs_row = 0; rhs_row < rhs.rows; ++rhs_row) {
                for (std::size_t rhs_col = 0; rhs_col < rhs.cols; ++rhs_col) {
                    result.at(row * rhs.rows + rhs_row,
                              col * rhs.cols + rhs_col) =
                        scale * rhs.at(rhs_row, rhs_col);
                }
            }
        }
    }
    return result;
}

Matrix hadamard(const Matrix& lhs, const Matrix& rhs) {
    require_same_shape(lhs, rhs, "hadamard");
    Matrix result(lhs.rows, lhs.cols, 0.0);
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] * rhs.data[i];
    }
    return result;
}

} // namespace matrix

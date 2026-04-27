#include "number_matrix.h"
#include "matrix.h"
#include "conversion.h"
#include "functions.h"

#include <sstream>

namespace matrix {

NumberMatrix::NumberMatrix(std::size_t row_count, std::size_t col_count,
                           const numeric::Number& fill_value,
                           const numeric::PrecisionContext& ctx)
    : rows(row_count), cols(col_count), context(ctx) {
    data.resize(rows * cols, fill_value);
}

NumberMatrix NumberMatrix::from_double_matrix(const Matrix& m,
                                              const numeric::PrecisionContext& ctx) {
    NumberMatrix result(m.rows, m.cols, numeric::Number(numeric::BigInt(0)), ctx);
    for (std::size_t i = 0; i < m.rows; ++i) {
        for (std::size_t j = 0; j < m.cols; ++j) {
            result.at(i, j) = m.at(i, j).with_context(ctx);
        }
    }
    return result;
}

Matrix NumberMatrix::to_double_matrix() const {
    Matrix result(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            result.at(i, j) = numeric::to_double(at(i, j));
        }
    }
    return result;
}

NumberMatrix NumberMatrix::zero(std::size_t row_count, std::size_t col_count,
                                const numeric::PrecisionContext& ctx) {
    return NumberMatrix(row_count, col_count, numeric::Number(numeric::BigInt(0)), ctx);
}

NumberMatrix NumberMatrix::identity(std::size_t size,
                                    const numeric::PrecisionContext& ctx) {
    NumberMatrix result(size, size, numeric::Number(numeric::BigInt(0)), ctx);
    for (std::size_t i = 0; i < size; ++i) {
        result.at(i, i) = numeric::Number(numeric::BigInt(1)).with_context(ctx);
    }
    return result;
}

bool NumberMatrix::is_vector() const {
    return rows == 1 || cols == 1;
}

bool NumberMatrix::is_square() const {
    return rows == cols;
}

numeric::Number& NumberMatrix::at(std::size_t row, std::size_t col) {
    return data[row * cols + col];
}

const numeric::Number& NumberMatrix::at(std::size_t row, std::size_t col) const {
    return data[row * cols + col];
}

void NumberMatrix::resize(std::size_t new_rows, std::size_t new_cols,
                          const numeric::Number& fill_value) {
    std::vector<numeric::Number> new_data(new_rows * new_cols, fill_value);
    for (std::size_t i = 0; i < rows && i < new_rows; ++i) {
        for (std::size_t j = 0; j < cols && j < new_cols; ++j) {
            new_data[i * new_cols + j] = at(i, j);
        }
    }
    data = std::move(new_data);
    rows = new_rows;
    cols = new_cols;
}

std::string NumberMatrix::to_string() const {
    std::ostringstream out;
    out << "[";
    for (std::size_t i = 0; i < rows; ++i) {
        if (i > 0) out << "; ";
        for (std::size_t j = 0; j < cols; ++j) {
            if (j > 0) out << ", ";
            out << at(i, j).to_string();
        }
    }
    out << "]";
    return out.str();
}

NumberMatrix add(const NumberMatrix& lhs, const NumberMatrix& rhs) {
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) {
        return NumberMatrix();
    }
    NumberMatrix result(lhs.rows, lhs.cols, numeric::Number(numeric::BigInt(0)), lhs.context);
    for (std::size_t i = 0; i < lhs.rows; ++i) {
        for (std::size_t j = 0; j < lhs.cols; ++j) {
            result.at(i, j) = numeric::add(lhs.at(i, j), rhs.at(i, j));
        }
    }
    return result;
}

NumberMatrix subtract(const NumberMatrix& lhs, const NumberMatrix& rhs) {
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) {
        return NumberMatrix();
    }
    NumberMatrix result(lhs.rows, lhs.cols, numeric::Number(numeric::BigInt(0)), lhs.context);
    for (std::size_t i = 0; i < lhs.rows; ++i) {
        for (std::size_t j = 0; j < lhs.cols; ++j) {
            result.at(i, j) = numeric::subtract(lhs.at(i, j), rhs.at(i, j));
        }
    }
    return result;
}

NumberMatrix multiply(const NumberMatrix& lhs, const NumberMatrix& rhs) {
    if (lhs.cols != rhs.rows) {
        return NumberMatrix();
    }
    NumberMatrix result(lhs.rows, rhs.cols, numeric::Number(numeric::BigInt(0)), lhs.context);
    for (std::size_t i = 0; i < lhs.rows; ++i) {
        for (std::size_t j = 0; j < rhs.cols; ++j) {
            numeric::Number sum(numeric::BigInt(0));
            for (std::size_t k = 0; k < lhs.cols; ++k) {
                sum = numeric::add(sum, numeric::multiply(lhs.at(i, k), rhs.at(k, j)));
            }
            result.at(i, j) = sum;
        }
    }
    return result;
}

NumberMatrix multiply(const NumberMatrix& matrix, const numeric::Number& scalar) {
    NumberMatrix result(matrix.rows, matrix.cols, numeric::Number(numeric::BigInt(0)), matrix.context);
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        for (std::size_t j = 0; j < matrix.cols; ++j) {
            result.at(i, j) = numeric::multiply(matrix.at(i, j), scalar);
        }
    }
    return result;
}

NumberMatrix divide(const NumberMatrix& matrix, const numeric::Number& scalar) {
    NumberMatrix result(matrix.rows, matrix.cols, numeric::Number(numeric::BigInt(0)), matrix.context);
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        for (std::size_t j = 0; j < matrix.cols; ++j) {
            result.at(i, j) = numeric::divide(matrix.at(i, j), scalar);
        }
    }
    return result;
}

NumberMatrix transpose(const NumberMatrix& matrix) {
    NumberMatrix result(matrix.cols, matrix.rows, numeric::Number(numeric::BigInt(0)), matrix.context);
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        for (std::size_t j = 0; j < matrix.cols; ++j) {
            result.at(j, i) = matrix.at(i, j);
        }
    }
    return result;
}

numeric::Number determinant(const NumberMatrix& matrix) {
    if (!matrix.is_square() || matrix.rows == 0) {
        return numeric::Number(numeric::BigInt(0));
    }

    std::size_t n = matrix.rows;
    if (n == 1) {
        return matrix.at(0, 0);
    }
    if (n == 2) {
        return numeric::subtract(
            numeric::multiply(matrix.at(0, 0), matrix.at(1, 1)),
            numeric::multiply(matrix.at(0, 1), matrix.at(1, 0)));
    }

    // LU decomposition for determinant
    NumberMatrix m = matrix;
    numeric::Number det(numeric::BigInt(1));

    for (std::size_t col = 0; col < n; ++col) {
        // Find pivot
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < n; ++row) {
            numeric::Number pivot_val = m.at(pivot, col);
            numeric::Number row_val = m.at(row, col);
            if (numeric::abs(row_val).compare(numeric::abs(pivot_val)) > 0) {
                pivot = row;
            }
        }

        if (pivot != col) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(m.at(col, j), m.at(pivot, j));
            }
            det = numeric::negate(det);
        }

        const numeric::Number& pivot_val = m.at(col, col);
        if (numeric::is_zero(pivot_val)) {
            return numeric::Number(numeric::BigInt(0));
        }

        det = numeric::multiply(det, pivot_val);

        for (std::size_t row = col + 1; row < n; ++row) {
            numeric::Number factor = numeric::divide(m.at(row, col), pivot_val);
            for (std::size_t j = col + 1; j < n; ++j) {
                m.at(row, j) = numeric::subtract(
                    m.at(row, j),
                    numeric::multiply(factor, m.at(col, j)));
            }
        }
    }

    return det;
}

NumberMatrix inverse(const NumberMatrix& matrix) {
    if (!matrix.is_square() || matrix.rows == 0) {
        return NumberMatrix();
    }

    std::size_t n = matrix.rows;
    NumberMatrix m = matrix;
    NumberMatrix result = NumberMatrix::identity(n, matrix.context);

    // Forward elimination
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < n; ++row) {
            if (numeric::abs(m.at(row, col)).compare(numeric::abs(m.at(pivot, col))) > 0) {
                pivot = row;
            }
        }

        if (pivot != col) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(m.at(col, j), m.at(pivot, j));
                std::swap(result.at(col, j), result.at(pivot, j));
            }
        }

        numeric::Number pivot_val = m.at(col, col);
        if (numeric::is_zero(pivot_val)) {
            return NumberMatrix();
        }

        for (std::size_t j = 0; j < n; ++j) {
            m.at(col, j) = numeric::divide(m.at(col, j), pivot_val);
            result.at(col, j) = numeric::divide(result.at(col, j), pivot_val);
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;
            numeric::Number factor = m.at(row, col);
            for (std::size_t j = 0; j < n; ++j) {
                m.at(row, j) = numeric::subtract(
                    m.at(row, j),
                    numeric::multiply(factor, m.at(col, j)));
                result.at(row, j) = numeric::subtract(
                    result.at(row, j),
                    numeric::multiply(factor, result.at(col, j)));
            }
        }
    }

    return result;
}

NumberMatrix solve(const NumberMatrix& coefficients, const NumberMatrix& rhs) {
    if (coefficients.rows != coefficients.cols || coefficients.rows != rhs.rows) {
        return NumberMatrix();
    }

    NumberMatrix inv = inverse(coefficients);
    if (inv.rows == 0) {
        return NumberMatrix();
    }

    return multiply(inv, rhs);
}

numeric::Number norm(const NumberMatrix& matrix) {
    numeric::Number sum(numeric::BigInt(0));
    for (const auto& elem : matrix.data) {
        sum = numeric::add(sum, numeric::multiply(elem, elem));
    }
    return numeric::sqrt(sum);
}

numeric::Number trace(const NumberMatrix& matrix) {
    if (!matrix.is_square()) {
        return numeric::Number(numeric::BigInt(0));
    }
    numeric::Number sum(numeric::BigInt(0));
    for (std::size_t i = 0; i < matrix.rows; ++i) {
        sum = numeric::add(sum, matrix.at(i, i));
    }
    return sum;
}

}  // namespace matrix

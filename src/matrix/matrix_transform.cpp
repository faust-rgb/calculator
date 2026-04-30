#include "matrix.h"
#include "matrix_internal.h"
#include <sstream>
#include <stdexcept>

namespace matrix {

using namespace internal;

Matrix reshape(const Matrix& matrix, std::size_t rows, std::size_t cols) {
    if (rows * cols != matrix.rows * matrix.cols) {
        throw std::runtime_error("reshape requires the element count to stay unchanged");
    }
    Matrix result(rows, cols, 0.0);
    result.data = matrix.data;
    return result;
}

Matrix vectorize(const Matrix& matrix) {
    if (matrix.is_vector()) {
        return matrix.rows == 1 ? transpose(matrix) : matrix;
    }
    Matrix result(matrix.rows * matrix.cols, 1, 0.0);
    std::size_t index = 0;
    for (std::size_t col = 0; col < matrix.cols; ++col) {
        for (std::size_t row = 0; row < matrix.rows; ++row) {
            result.at(index, 0) = matrix.at(row, col);
            ++index;
        }
    }
    return result;
}

Matrix diag(const Matrix& matrix) {
    if (matrix.is_vector()) {
        const std::vector<double> values = as_vector_values(matrix, "diag");
        Matrix result(values.size(), values.size(), 0.0);
        for (std::size_t i = 0; i < values.size(); ++i) {
            result.at(i, i) = values[i];
        }
        return result;
    }

    const std::size_t diagonal = matrix.rows < matrix.cols ? matrix.rows : matrix.cols;
    Matrix result(diagonal, 1, 0.0);
    for (std::size_t i = 0; i < diagonal; ++i) {
        result.at(i, 0) = matrix.at(i, i);
    }
    return result;
}

Matrix transpose(const Matrix& matrix) {
    Matrix result(matrix.cols, matrix.rows, 0.0);
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            result.at(col, row) = matrix.at(row, col);
        }
    }
    return result;
}

void Matrix::resize(std::size_t new_rows,
                    std::size_t new_cols,
                    double fill_value) {
    std::vector<double> resized(new_rows * new_cols, fill_value);
    const std::size_t shared_rows = rows < new_rows ? rows : new_rows;
    const std::size_t shared_cols = cols < new_cols ? cols : new_cols;

    for (std::size_t row = 0; row < shared_rows; ++row) {
        for (std::size_t col = 0; col < shared_cols; ++col) {
            resized[row * new_cols + col] = at(row, col);
        }
    }

    rows = new_rows;
    cols = new_cols;
    data.swap(resized);
}

void Matrix::append_row(const std::vector<double>& values) {
    if (cols == 0) {
        rows = 1;
        cols = values.size();
        data = values;
        return;
    }

    if (values.size() > cols) {
        resize(rows, values.size(), 0.0);
    }

    data.reserve(data.size() + cols);
    for (std::size_t col = 0; col < cols; ++col) {
        data.push_back(col < values.size() ? values[col] : 0.0);
    }
    ++rows;
}

void Matrix::append_col(const std::vector<double>& values) {
    if (rows == 0) {
        rows = values.size();
        cols = 1;
        data = values;
        return;
    }

    if (values.size() > rows) {
        resize(values.size(), cols, 0.0);
    }

    std::vector<double> resized;
    resized.reserve(rows * (cols + 1));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
            resized.push_back(at(row, col));
        }
        resized.push_back(row < values.size() ? values[row] : 0.0);
    }

    ++cols;
    data.swap(resized);
}

std::string Matrix::to_string() const {
    if (rows == 1) {
        std::ostringstream out;
        out << "[";
        for (std::size_t col = 0; col < cols; ++col) {
            if (col != 0) {
                out << ", ";
            }
            out << format_number(at(0, col));
        }
        out << "]";
        return out.str();
    }

    std::ostringstream out;
    out << "[";
    for (std::size_t row = 0; row < rows; ++row) {
        if (row != 0) {
            out << ", ";
        }
        out << "[";
        for (std::size_t col = 0; col < cols; ++col) {
            if (col != 0) {
                out << ", ";
            }
            out << format_number(at(row, col));
        }
        out << "]";
    }
    out << "]";
    return out.str();
}

} // namespace matrix

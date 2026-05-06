#include "matrix.h"
#include "matrix_internal.h"
#include "precise/precise_decimal.h"
#include <sstream>
#include <stdexcept>

namespace matrix {

using namespace internal;

template <typename T>
TMatrix<T> reshape(const TMatrix<T>& matrix, std::size_t rows, std::size_t cols) {
    if (rows * cols != matrix.rows * matrix.cols) {
        throw std::runtime_error("reshape requires the element count to stay unchanged");
    }
    TMatrix<T> result(rows, cols, T(static_cast<long long>(0)));
    result.data = matrix.data;
    return result;
}

template <typename T>
TMatrix<T> vectorize(const TMatrix<T>& matrix) {
    if (matrix.is_vector()) {
        return matrix.rows == 1 ? matrix::transpose<T>(matrix) : matrix;
    }
    TMatrix<T> result(matrix.rows * matrix.cols, 1, T(static_cast<long long>(0)));
    std::size_t index = 0;
    for (std::size_t col = 0; col < matrix.cols; ++col) {
        for (std::size_t row = 0; row < matrix.rows; ++row) {
            result.at(index, 0) = matrix.at(row, col);
            ++index;
        }
    }
    return result;
}

template <typename T>
TMatrix<T> diag(const TMatrix<T>& matrix) {
    if (matrix.is_vector()) {
        const std::vector<T> values = internal::as_vector_values<T>(matrix, "diag");
        TMatrix<T> result(values.size(), values.size(), T(static_cast<long long>(0)));
        for (std::size_t i = 0; i < values.size(); ++i) {
            result.at(i, i) = values[i];
        }
        return result;
    }

    const std::size_t diagonal = matrix.rows < matrix.cols ? matrix.rows : matrix.cols;
    TMatrix<T> result(diagonal, 1, T(static_cast<long long>(0)));
    for (std::size_t i = 0; i < diagonal; ++i) {
        result.at(i, 0) = matrix.at(i, i);
    }
    return result;
}

template <typename T>
TMatrix<T> transpose(const TMatrix<T>& matrix) {
    TMatrix<T> result(matrix.cols, matrix.rows, T(static_cast<long long>(0)));
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            result.at(col, row) = matrix.at(row, col);
        }
    }
    return result;
}

template <typename T>
void TMatrix<T>::resize(std::size_t new_rows,
                        std::size_t new_cols,
                        T fill_value) {
    std::vector<T> resized(new_rows * new_cols, fill_value);
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

template <typename T>
void TMatrix<T>::append_row(const std::vector<T>& values) {
    if (cols == 0) {
        rows = 1;
        cols = values.size();
        data = values;
        return;
    }

    if (values.size() > cols) {
        resize(rows, values.size(), T(static_cast<long long>(0)));
    }

    data.reserve(data.size() + cols);
    for (std::size_t col = 0; col < cols; ++col) {
        data.push_back(col < values.size() ? values[col] : T(static_cast<long long>(0)));
    }
    ++rows;
}

template <typename T>
void TMatrix<T>::append_col(const std::vector<T>& values) {
    if (rows == 0) {
        rows = values.size();
        cols = 1;
        data = values;
        return;
    }

    if (values.size() > rows) {
        resize(values.size(), cols, T(static_cast<long long>(0)));
    }

    std::vector<T> resized;
    resized.reserve(rows * (cols + 1));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
            resized.push_back(at(row, col));
        }
        resized.push_back(row < values.size() ? values[row] : T(static_cast<long long>(0)));
    }

    ++cols;
    data.swap(resized);
}

template <typename T>
std::string TMatrix<T>::to_string() const {
    if (rows == 1) {
        std::ostringstream out;
        out << "[";
        for (std::size_t col = 0; col < cols; ++col) {
            if (col != 0) {
                out << ", ";
            }
            out << internal::format_number<T>(at(0, col));
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
            out << internal::format_number<T>(at(row, col));
        }
        out << "]";
    }
    out << "]";
    return out.str();
}

// Explicit template instantiations
template TMatrix<double> reshape<double>(const TMatrix<double>&, std::size_t, std::size_t);
template TMatrix<PreciseDecimal> reshape<PreciseDecimal>(const TMatrix<PreciseDecimal>&, std::size_t, std::size_t);

template TMatrix<double> vectorize<double>(const TMatrix<double>&);
template TMatrix<PreciseDecimal> vectorize<PreciseDecimal>(const TMatrix<PreciseDecimal>&);

template TMatrix<double> diag<double>(const TMatrix<double>&);
template TMatrix<PreciseDecimal> diag<PreciseDecimal>(const TMatrix<PreciseDecimal>&);

template TMatrix<double> transpose<double>(const TMatrix<double>&);
template TMatrix<PreciseDecimal> transpose<PreciseDecimal>(const TMatrix<PreciseDecimal>&);

template void TMatrix<double>::resize(std::size_t, std::size_t, double);
template void TMatrix<PreciseDecimal>::resize(std::size_t, std::size_t, PreciseDecimal);

template void TMatrix<double>::append_row(const std::vector<double>&);
template void TMatrix<PreciseDecimal>::append_row(const std::vector<PreciseDecimal>&);

template void TMatrix<double>::append_col(const std::vector<double>&);
template void TMatrix<PreciseDecimal>::append_col(const std::vector<PreciseDecimal>&);

template std::string TMatrix<double>::to_string() const;
template std::string TMatrix<PreciseDecimal>::to_string() const;

} // namespace matrix

#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include "precise/precise_decimal.h"
#include <stdexcept>
#include <utility>
#include <vector>
#include <algorithm>

namespace matrix {

using namespace internal;

// ============================================================================
// TMatrix 核心方法实现 (基础构造与访问)
// ============================================================================

template <typename T>
TMatrix<T>::TMatrix(std::size_t row_count, std::size_t col_count, T fill_value)
    : rows(row_count),
      cols(col_count),
      data(row_count * col_count, fill_value) {}

template <typename T>
TMatrix<T> TMatrix<T>::vector(const std::vector<T>& values) {
    TMatrix<T> m(1, values.size(), T(0));
    m.data = values;
    return m;
}

template <typename T>
TMatrix<T> TMatrix<T>::zero(std::size_t row_count, std::size_t col_count) {
    return TMatrix<T>(row_count, col_count, T(0));
}

template <typename T>
TMatrix<T> TMatrix<T>::identity(std::size_t size) {
    TMatrix<T> m(size, size, T(0));
    for (std::size_t i = 0; i < size; ++i) m.at(i, i) = T(1);
    return m;
}

template <typename T> bool TMatrix<T>::is_vector() const { return rows == 1 || cols == 1; }
template <typename T> bool TMatrix<T>::is_square() const { return rows == cols; }

template <typename T>
T& TMatrix<T>::at(std::size_t row, std::size_t col) {
    if (row >= rows || col >= cols) throw std::out_of_range("matrix index out of range");
    return data[row * cols + col];
}

template <typename T>
T TMatrix<T>::at(std::size_t row, std::size_t col) const {
    if (row >= rows || col >= cols) throw std::out_of_range("matrix index out of range");
    return data[row * cols + col];
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator+=(const TMatrix<T>& rhs) {
    if (rows != rhs.rows || cols != rhs.cols) throw std::runtime_error("shape mismatch");
    for (std::size_t i = 0; i < data.size(); ++i) data[i] += rhs.data[i];
    return *this;
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator-=(const TMatrix<T>& rhs) {
    if (rows != rhs.rows || cols != rhs.cols) throw std::runtime_error("shape mismatch");
    for (std::size_t i = 0; i < data.size(); ++i) data[i] -= rhs.data[i];
    return *this;
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator*=(T scalar) {
    for (T& val : data) val *= scalar;
    return *this;
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator/=(T scalar) {
    if (scalar == T(0)) throw std::runtime_error("division by zero");
    for (T& val : data) val /= scalar;
    return *this;
}

// ============================================================================
// TValue 方法实现
// ============================================================================

template <typename T> TValue<T> TValue<T>::from_scalar(T v) { TValue<T> r; r.scalar = v; return r; }
template <typename T> TValue<T> TValue<T>::from_complex(T re, T im) { TValue<T> r; r.is_complex = true; r.complex.real = re; r.complex.imag = im; return r; }
template <typename T> TValue<T> TValue<T>::from_complex(TComplex<T> c) { TValue<T> r; r.is_complex = true; r.complex = c; return r; }
template <typename T> TValue<T> TValue<T>::from_matrix(const TMatrix<T>& m) { TValue<T> r; r.is_matrix = true; r.matrix = m; return r; }
template <typename T> TValue<T> TValue<T>::from_matrix(TMatrix<T>&& m) { TValue<T> r; r.is_matrix = true; r.matrix = std::move(m); return r; }

} // namespace matrix

// 显式实例化 (Core)
template class matrix::TMatrix<double>;
template class matrix::TMatrix<PreciseDecimal>;
template class matrix::TValue<double>;
template class matrix::TValue<PreciseDecimal>;

#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include "precise/precise_decimal.h"
#include <stdexcept>

namespace matrix {

using namespace internal;

template <typename T>
TMatrix<T> add(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    TMatrix<T> result = lhs;
    result += rhs;
    return result;
}

template <typename T>
TMatrix<T> add(TMatrix<T>&& lhs, const TMatrix<T>& rhs) {
    lhs += rhs;
    return std::move(lhs);
}

template <typename T>
TMatrix<T> add(const TMatrix<T>& lhs, TMatrix<T>&& rhs) {
    rhs += lhs;
    return std::move(rhs);
}

template <typename T>
TMatrix<T> add(const TMatrix<T>& lhs, T scalar) {
    TMatrix<T> result(lhs.rows, lhs.cols, T(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] + scalar;
    }
    return result;
}

template <typename T>
TMatrix<T> add(TMatrix<T>&& lhs, T scalar) {
    for (T& val : lhs.data) {
        val += scalar;
    }
    return std::move(lhs);
}

template <typename T>
TMatrix<T> subtract(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    TMatrix<T> result = lhs;
    result -= rhs;
    return result;
}

template <typename T>
TMatrix<T> subtract(TMatrix<T>&& lhs, const TMatrix<T>& rhs) {
    lhs -= rhs;
    return std::move(lhs);
}

template <typename T>
TMatrix<T> subtract(const TMatrix<T>& lhs, TMatrix<T>&& rhs) {
    require_same_shape(lhs, rhs, "subtract");
    for (std::size_t i = 0; i < rhs.data.size(); ++i) {
        rhs.data[i] = lhs.data[i] - rhs.data[i];
    }
    return std::move(rhs);
}

template <typename T>
TMatrix<T> subtract(const TMatrix<T>& lhs, T scalar) {
    TMatrix<T> result(lhs.rows, lhs.cols, T(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] - scalar;
    }
    return result;
}

template <typename T>
TMatrix<T> subtract(TMatrix<T>&& lhs, T scalar) {
    for (T& val : lhs.data) {
        val -= scalar;
    }
    return std::move(lhs);
}

template <typename T>
TMatrix<T> multiply(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    if (lhs.cols != rhs.rows) {
        throw std::runtime_error("matrix multiplication requires lhs.cols == rhs.rows");
    }

    TMatrix<T> result(lhs.rows, rhs.cols, T(0));
    for (std::size_t i = 0; i < lhs.rows; ++i) {
        for (std::size_t k = 0; k < lhs.cols; ++k) {
            const T& lhs_val = lhs.at(i, k);
            if (lhs_val == T(0)) continue;
            
            const std::size_t row_offset = i * rhs.cols;
            const std::size_t rhs_offset = k * rhs.cols;
            for (std::size_t j = 0; j < rhs.cols; ++j) {
                result.data[row_offset + j] += lhs_val * rhs.data[rhs_offset + j];
            }
        }
    }
    return result;
}

// 针对 double 的优化版本 (使用 long double 中间累加)
template <>
TMatrix<double> multiply<double>(const TMatrix<double>& lhs, const TMatrix<double>& rhs) {
    if (lhs.cols != rhs.rows) {
        throw std::runtime_error("matrix multiplication requires lhs.cols == rhs.rows");
    }
    TMatrix<double> result(lhs.rows, rhs.cols, 0.0);
    std::vector<long double> sums(lhs.rows * rhs.cols, 0.0L);
    for (std::size_t i = 0; i < lhs.rows; ++i) {
        for (std::size_t k = 0; k < lhs.cols; ++k) {
            const long double lhs_val = static_cast<long double>(lhs.at(i, k));
            if (lhs_val == 0.0L) continue;
            const std::size_t row_offset = i * rhs.cols;
            const std::size_t rhs_offset = k * rhs.cols;
            for (std::size_t j = 0; j < rhs.cols; ++j) {
                sums[row_offset + j] += lhs_val * static_cast<long double>(rhs.data[rhs_offset + j]);
            }
        }
    }
    for (std::size_t i = 0; i < sums.size(); ++i) result.data[i] = static_cast<double>(sums[i]);
    return result;
}

template <typename T>
TMatrix<T> multiply(const TMatrix<T>& lhs, T scalar) {
    TMatrix<T> result = lhs;
    result *= scalar;
    return result;
}

template <typename T>
TMatrix<T> multiply(TMatrix<T>&& lhs, T scalar) {
    lhs *= scalar;
    return std::move(lhs);
}

template <typename T>
TMatrix<T> divide(const TMatrix<T>& lhs, T scalar) {
    TMatrix<T> result = lhs;
    result /= scalar;
    return result;
}

template <typename T>
TMatrix<T> divide(TMatrix<T>&& lhs, T scalar) {
    lhs /= scalar;
    return std::move(lhs);
}

template <typename T>
T dot(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "dot");
    const std::size_t rhs_size = vector_length(rhs, "dot");
    if (lhs_size != rhs_size) {
        throw std::runtime_error("dot requires vectors of the same length");
    }
    T sum = T(0);
    for (std::size_t i = 0; i < lhs_size; ++i) {
        sum += vector_entry(lhs, i) * vector_entry(rhs, i);
    }
    return sum;
}

// 针对 double 的优化版 dot (Kahan 求和)
template <>
double dot<double>(const TMatrix<double>& lhs, const TMatrix<double>& rhs) {
    const std::size_t n = vector_length(lhs, "dot");
    const std::size_t m = vector_length(rhs, "dot");
    if (n != m) throw std::runtime_error("dot size mismatch");
    long double sum = 0.0L;
    long double c = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
        long double y = static_cast<long double>(vector_entry(lhs, i)) * static_cast<long double>(vector_entry(rhs, i)) - c;
        long double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return static_cast<double>(sum);
}

template <typename T>
TMatrix<T> outer(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "outer");
    const std::size_t rhs_size = vector_length(rhs, "outer");
    TMatrix<T> result(lhs_size, rhs_size, T(0));
    for (std::size_t row = 0; row < lhs_size; ++row) {
        for (std::size_t col = 0; col < rhs_size; ++col) {
            result.at(row, col) = vector_entry(lhs, row) * vector_entry(rhs, col);
        }
    }
    return result;
}

template <typename T>
TMatrix<T> cross(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    const std::size_t lhs_size = vector_length(lhs, "cross");
    const std::size_t rhs_size = vector_length(rhs, "cross");
    if (lhs_size != 3 || rhs_size != 3) {
        throw std::runtime_error("cross product is only defined for 3D vectors");
    }
    const T x1 = vector_entry(lhs, 0);
    const T y1 = vector_entry(lhs, 1);
    const T z1 = vector_entry(lhs, 2);
    const T x2 = vector_entry(rhs, 0);
    const T y2 = vector_entry(rhs, 1);
    const T z2 = vector_entry(rhs, 2);
    std::vector<T> res_data = {
        y1 * z2 - z1 * y2,
        z1 * x2 - x1 * z2,
        x1 * y2 - y1 * x2
    };
    if (lhs.rows == 1) {
        TMatrix<T> res(1, 3);
        res.data = res_data;
        return res;
    } else {
        return TMatrix<T>::vector(res_data);
    }
}

template <typename T>
TMatrix<T> project(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    const T rhs_norm_sq = dot(rhs, rhs);
    if (rhs_norm_sq == T(0)) {
        throw std::runtime_error("cannot project onto a zero vector");
    }
    return multiply(rhs, dot(lhs, rhs) / rhs_norm_sq);
}

template <typename T>
TMatrix<T> normalize(const TMatrix<T>& matrix) {
    if (!matrix.is_vector()) {
        throw std::runtime_error("normalize only works on vectors");
    }
    const T v_norm = norm(matrix);
    if (v_norm == T(0)) {
        throw std::runtime_error("cannot normalize a zero vector");
    }
    return divide(matrix, v_norm);
}

template <typename T>
TMatrix<T> kronecker(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    TMatrix<T> result(lhs.rows * rhs.rows, lhs.cols * rhs.cols, T(0));
    for (std::size_t row = 0; row < lhs.rows; ++row) {
        for (std::size_t col = 0; col < lhs.cols; ++col) {
            const T& scale = lhs.at(row, col);
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

template <typename T>
TMatrix<T> hadamard(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    require_same_shape(lhs, rhs, "hadamard");
    TMatrix<T> result(lhs.rows, lhs.cols, T(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] * rhs.data[i];
    }
    return result;
}

// 显式实例化
#define INSTANTIATE_OPS(TYPE) \
    template TMatrix<TYPE> add(const TMatrix<TYPE>&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> add(TMatrix<TYPE>&&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> add(const TMatrix<TYPE>&, TMatrix<TYPE>&&); \
    template TMatrix<TYPE> add(const TMatrix<TYPE>&, TYPE); \
    template TMatrix<TYPE> add(TMatrix<TYPE>&&, TYPE); \
    template TMatrix<TYPE> subtract(const TMatrix<TYPE>&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> subtract(TMatrix<TYPE>&&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> subtract(const TMatrix<TYPE>&, TMatrix<TYPE>&&); \
    template TMatrix<TYPE> subtract(const TMatrix<TYPE>&, TYPE); \
    template TMatrix<TYPE> subtract(TMatrix<TYPE>&&, TYPE); \
    template TMatrix<TYPE> multiply(const TMatrix<TYPE>&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> multiply(const TMatrix<TYPE>&, TYPE); \
    template TMatrix<TYPE> multiply(TMatrix<TYPE>&&, TYPE); \
    template TMatrix<TYPE> divide(const TMatrix<TYPE>&, TYPE); \
    template TMatrix<TYPE> divide(TMatrix<TYPE>&&, TYPE); \
    template TYPE dot(const TMatrix<TYPE>&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> outer(const TMatrix<TYPE>&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> cross(const TMatrix<TYPE>&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> project(const TMatrix<TYPE>&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> normalize(const TMatrix<TYPE>&); \
    template TMatrix<TYPE> kronecker(const TMatrix<TYPE>&, const TMatrix<TYPE>&); \
    template TMatrix<TYPE> hadamard(const TMatrix<TYPE>&, const TMatrix<TYPE>&);

INSTANTIATE_OPS(double)
INSTANTIATE_OPS(PreciseDecimal)

} // namespace matrix

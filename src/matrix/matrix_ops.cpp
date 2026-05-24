/**
 * @file matrix_ops.cpp
 * @brief 矩阵基本运算实现
 *
 * 本文件实现了矩阵的基本算术运算，包括：
 * - 加法运算：矩阵+矩阵、矩阵+标量（支持移动语义优化）
 * - 减法运算：矩阵-矩阵、矩阵-标量
 * - 乘法运算：矩阵乘法、标量乘法（包含 long double 和 PreciseDecimal 的优化版本）
 * - 除法运算：矩阵/标量
 * - 向量运算：点积（使用 Kahan 求和优化）、外积、叉积、投影、归一化
 * - 特殊积：Kronecker 积、Hadamard 积
 *
 * 所有函数都支持泛型模板，可以处理 long double 和 PreciseDecimal 类型。
 *
 * @author Calculator Team
 * @date 2024
 */

#include "matrix.h"
#include "matrix_internal.h"
#include "math/mymath.h"
#include "math/precise/precise_decimal.h"
#include "math/types/float128.h"
#include <stdexcept>
#include <type_traits>

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

// ============================================================================
// 矩阵乘法实现
// ============================================================================

/**
 * @brief 矩阵乘法（延迟规范化版本）
 *
 * 使用延迟规范化优化：在累加过程中不规范化，
 * 最后统一规范化，减少规范化开销。
 * 对于浮点类型，规范化为空操作。
 *
 * @param lhs 左操作数矩阵 (m×k)
 * @param rhs 右操作数矩阵 (k×n)
 * @return 乘积矩阵 (m×n)
 * @throws 如果 lhs.cols != rhs.rows 则抛出异常
 */
template <typename T>
TMatrix<T> multiply(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    if (lhs.cols != rhs.rows) {
        throw std::runtime_error("matrix multiplication requires lhs.cols == rhs.rows");
    }

    TMatrix<T> result(lhs.rows, rhs.cols, T(0));

    if constexpr (std::is_same_v<T, PreciseDecimal>) {
        for (std::size_t i = 0; i < lhs.rows; ++i) {
            for (std::size_t j = 0; j < rhs.cols; ++j) {
                T sum = T(0);
                for (std::size_t k = 0; k < lhs.cols; ++k) {
                    const T& lhs_val = lhs.at(i, k);
                    const T& rhs_val = rhs.at(k, j);
                    if (lhs_val == T(0) || rhs_val == T(0)) continue;
                    sum += lhs_val * rhs_val;
                }
                result.at(i, j) = sum;
            }
        }
        return result;
    }

    for (std::size_t i = 0; i < lhs.rows; ++i) {
        for (std::size_t k = 0; k < lhs.cols; ++k) {
            const T& lhs_val = lhs.at(i, k);
            if (lhs_val == T(0)) continue;

            const std::size_t row_offset = i * rhs.cols;
            const std::size_t rhs_offset = k * rhs.cols;
            for (std::size_t j = 0; j < rhs.cols; ++j) {
                multiply_add_without_normalize(result.data[row_offset + j], lhs_val, rhs.data[rhs_offset + j]);
            }
        }
    }

    // 最终对结果矩阵全量规范化。矩阵乘法的 PreciseDecimal 快速累加会
    // 延迟规范化中间值，收尾时需要确保规范化没有被外层抑制状态跳过。
    precise::ScopedNormalizationEnable enable_normalization;
    for (auto& val : result.data) normalize(val);

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

// ============================================================================
// 向量运算实现
// ============================================================================

/**
 * @brief 向量点积（泛型版本）
 *
 * 计算两个向量的点积（内积）。支持行向量和列向量。
 *
 * @param lhs 第一个向量
 * @param rhs 第二个向量
 * @return 点积值
 * @throws 如果向量长度不匹配则抛出异常
 */
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
#if 0
/**
 * @brief 向量点积 - long double 特化版本
 *
 * 使用 Kahan 补偿求和算法，减少浮点累加误差。
 * 对于大量元素的点积运算，精度损失显著降低。
 */
// 针对 long double 的优化版 dot (Kahan 求和)
template <>
long double dot<long double>(const TMatrix<long double>& lhs, const TMatrix<long double>& rhs) {
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
    return static_cast<long double>(sum);
}
#endif
/**
 * @brief 向量外积
 *
 * 计算两个向量的外积（张量积）。结果是一个矩阵。
 * 例如，[a, b] 与 [c, d] 的外积为 [[ac, ad], [bc, bd]]。
 *
 * @param lhs 第一个向量（结果矩阵的行数由其长度决定）
 * @param rhs 第二个向量（结果矩阵的列数由其长度决定）
 * @return 外积矩阵
 */
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

/**
 * @brief 向量叉积
 *
 * 计算两个三维向量的叉积。仅支持三维向量。
 *
 * @param lhs 第一个三维向量
 * @param rhs 第二个三维向量
 * @return 叉积向量（方向遵循右手定则）
 * @throws 如果输入不是三维向量则抛出异常
 */
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

/**
 * @brief 向量投影
 *
 * 计算向量 lhs 在向量 rhs 方向上的投影向量。
 * 投影向量为 (lhs·rhs / rhs·rhs) * rhs。
 *
 * @param lhs 被投影的向量
 * @param rhs 投影方向向量
 * @return 投影向量
 * @throws 如果 rhs 是零向量则抛出异常
 */
template <typename T>
TMatrix<T> project(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    const T rhs_norm_sq = dot(rhs, rhs);
    if (rhs_norm_sq == T(0)) {
        throw std::runtime_error("cannot project onto a zero vector");
    }
    return multiply(rhs, dot(lhs, rhs) / rhs_norm_sq);
}

/**
 * @brief 向量归一化
 *
 * 将向量归一化为单位向量（长度为1）。
 *
 * @param matrix 输入向量
 * @return 单位向量
 * @throws 如果输入不是向量或为零向量则抛出异常
 */
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

// ============================================================================
// 特殊矩阵积实现
// ============================================================================

/**
 * @brief Kronecker 积
 *
 * 计算两个矩阵的 Kronecker 积（张量积）。
 * 若 A 为 m×n，B 为 p×q，则结果为 mp×nq 矩阵。
 * 结果矩阵的每个元素 a[i][j] 被替换为 a[i][j] * B。
 *
 * @param lhs 左操作数矩阵
 * @param rhs 右操作数矩阵
 * @return Kronecker 积矩阵
 */
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

/**
 * @brief Hadamard 积（逐元素乘积）
 *
 * 计算两个相同形状矩阵的逐元素乘积。
 * 要求两个矩阵的行列数完全相同。
 *
 * @param lhs 左操作数矩阵
 * @param rhs 右操作数矩阵
 * @return 逐元素乘积矩阵
 * @throws 如果矩阵形状不匹配则抛出异常
 */
template <typename T>
TMatrix<T> hadamard(const TMatrix<T>& lhs, const TMatrix<T>& rhs) {
    require_same_shape(lhs, rhs, "hadamard");
    TMatrix<T> result(lhs.rows, lhs.cols, T(0));
    for (std::size_t i = 0; i < lhs.data.size(); ++i) {
        result.data[i] = lhs.data[i] * rhs.data[i];
    }
    return result;
}

// ============================================================================
// 显式模板实例化
// ============================================================================

/**
 * @brief 模板实例化宏
 *
 * 为指定类型实例化所有矩阵运算函数模板。
 * 支持 long double 和 PreciseDecimal 两种类型。
 */
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

//INSTANTIATE_OPS(long double)
//INSTANTIATE_OPS(PreciseDecimal)
INSTANTIATE_OPS(mymath::Scalar)

} // namespace matrix

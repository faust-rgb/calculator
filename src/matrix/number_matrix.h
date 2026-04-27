#ifndef MATRIX_NUMBER_MATRIX_H
#define MATRIX_NUMBER_MATRIX_H

#include "number.h"
#include "precision_context.h"

#include <cstddef>
#include <string>
#include <vector>

namespace matrix {

/**
 * @struct NumberMatrix
 * @brief 使用任意精度 Number 的矩阵数据结构
 *
 * 与 Matrix 类似，但使用 numeric::Number 存储元素，
 * 支持任意精度运算。
 */
struct NumberMatrix {
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<numeric::Number> data;
    numeric::PrecisionContext context;

    NumberMatrix() = default;
    NumberMatrix(std::size_t row_count, std::size_t col_count,
                 const numeric::Number& fill_value = numeric::Number(numeric::BigInt(0)),
                 const numeric::PrecisionContext& ctx = numeric::PrecisionContext{});

    static NumberMatrix from_double_matrix(const struct Matrix& m,
                                           const numeric::PrecisionContext& ctx = numeric::PrecisionContext{});
    struct Matrix to_double_matrix() const;

    static NumberMatrix zero(std::size_t row_count, std::size_t col_count,
                             const numeric::PrecisionContext& ctx = numeric::PrecisionContext{});
    static NumberMatrix identity(std::size_t size,
                                 const numeric::PrecisionContext& ctx = numeric::PrecisionContext{});

    bool is_vector() const;
    bool is_square() const;

    numeric::Number& at(std::size_t row, std::size_t col);
    const numeric::Number& at(std::size_t row, std::size_t col) const;

    void resize(std::size_t new_rows, std::size_t new_cols,
                const numeric::Number& fill_value = numeric::Number(numeric::BigInt(0)));

    std::string to_string() const;
};

// 基本运算
NumberMatrix add(const NumberMatrix& lhs, const NumberMatrix& rhs);
NumberMatrix subtract(const NumberMatrix& lhs, const NumberMatrix& rhs);
NumberMatrix multiply(const NumberMatrix& lhs, const NumberMatrix& rhs);
NumberMatrix multiply(const NumberMatrix& matrix, const numeric::Number& scalar);
NumberMatrix divide(const NumberMatrix& matrix, const numeric::Number& scalar);
NumberMatrix transpose(const NumberMatrix& matrix);

// 线性代数
NumberMatrix inverse(const NumberMatrix& matrix);
numeric::Number determinant(const NumberMatrix& matrix);
NumberMatrix solve(const NumberMatrix& coefficients, const NumberMatrix& rhs);

// 范数和迹
numeric::Number norm(const NumberMatrix& matrix);
numeric::Number trace(const NumberMatrix& matrix);

}  // namespace matrix

#endif

#ifndef RUNTIME_MATRIX_VALUE_H
#define RUNTIME_MATRIX_VALUE_H

#include "expr.h"

#include <cstddef>
#include <string>
#include <vector>

namespace runtime {

class MatrixValue {
public:
    MatrixValue();
    MatrixValue(std::size_t rows, std::size_t cols, const expression::Expr& fill);
    explicit MatrixValue(const std::vector<std::vector<expression::Expr>>& entries);

    std::size_t rows() const;
    std::size_t cols() const;
    bool empty() const;

    const expression::Expr& at(std::size_t row, std::size_t col) const;
    expression::Expr& at(std::size_t row, std::size_t col);
    const std::vector<std::vector<expression::Expr>>& entries() const;

    std::string to_string() const;

private:
    std::vector<std::vector<expression::Expr>> entries_;
};

}  // namespace runtime

#endif

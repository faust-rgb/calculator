#include "matrix_value.h"

#include "printer.h"

#include <sstream>
#include <stdexcept>

namespace runtime {

MatrixValue::MatrixValue() = default;

MatrixValue::MatrixValue(std::size_t rows,
                         std::size_t cols,
                         const expression::Expr& fill)
    : entries_(rows, std::vector<expression::Expr>(cols, fill)) {}

MatrixValue::MatrixValue(const std::vector<std::vector<expression::Expr>>& entries)
    : entries_(entries) {
    if (entries_.empty()) {
        return;
    }
    const std::size_t width = entries_[0].size();
    if (width == 0) {
        throw std::runtime_error("matrix rows cannot be empty");
    }
    for (const auto& row : entries_) {
        if (row.size() != width) {
            throw std::runtime_error("matrix rows must have the same length");
        }
    }
}

std::size_t MatrixValue::rows() const {
    return entries_.size();
}

std::size_t MatrixValue::cols() const {
    return entries_.empty() ? 0 : entries_[0].size();
}

bool MatrixValue::empty() const {
    return rows() == 0 || cols() == 0;
}

const expression::Expr& MatrixValue::at(std::size_t row, std::size_t col) const {
    return entries_.at(row).at(col);
}

expression::Expr& MatrixValue::at(std::size_t row, std::size_t col) {
    return entries_.at(row).at(col);
}

const std::vector<std::vector<expression::Expr>>& MatrixValue::entries() const {
    return entries_;
}

std::string MatrixValue::to_string() const {
    std::ostringstream out;
    out << "[";
    for (std::size_t row = 0; row < rows(); ++row) {
        if (row != 0) {
            out << ", ";
        }
        out << "[";
        for (std::size_t col = 0; col < cols(); ++col) {
            if (col != 0) {
                out << ", ";
            }
            out << expression::print(at(row, col));
        }
        out << "]";
    }
    out << "]";
    return out.str();
}

}  // namespace runtime

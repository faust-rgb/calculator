// ============================================================================
// 矩阵命令实现
// ============================================================================

#include "calculator_matrix_commands.h"

#include "calculator_internal_types.h"
#include "parser/unified_expression_parser.h"
#include "matrix_internal.h"
#include "mymath.h"

#include <sstream>
#include <stdexcept>

namespace matrix_commands {
namespace {

std::string format_eigenvalue_matrix(const matrix::Matrix& values) {
    if (values.rows <= 1 || values.cols != 2) {
        return values.to_string();
    }
    std::ostringstream out;
    out << "[";
    for (std::size_t row = 0; row < values.rows; ++row) {
        if (row != 0) {
            out << ", ";
        }
        out << matrix::internal::format_complex<double>({values.at(row, 0),
                                                           values.at(row, 1)});    }
    out << "]";
    return out.str();
}

}  // namespace

bool is_matrix_command(const std::string& command) {
    return command == "eig" || command == "svd" || command == "lu_p";
}

bool handle_matrix_command(const MatrixCommandContext& ctx,
                           const std::string& command,
                           const std::string& inside,
                           std::string* output) {
    if (!is_matrix_command(command)) {
        return false;
    }

    const std::vector<std::string> arguments = split_top_level_arguments(inside);
    if (arguments.size() != 1 || !ctx.is_matrix_argument(arguments[0])) {
        throw std::runtime_error(command + " expects exactly one matrix argument");
    }

    const matrix::Matrix matrix_value =
        ctx.parse_matrix_argument(arguments[0], command);
    if (command == "svd") {
        *output = "U: " + matrix::svd_u(matrix_value).to_string() +
                  "\nS: " + matrix::svd_s(matrix_value).to_string() +
                  "\nVt: " + matrix::svd_vt(matrix_value).to_string();
        return true;
    }

    if (command == "lu_p") {
        *output = matrix::lu_p(matrix_value).to_string();
        return true;
    }

    try {
        const matrix::Matrix values = matrix::eigenvalues(matrix_value);
        if (values.rows > 1 && values.cols == 2) {
            *output = "values: " + format_eigenvalue_matrix(values) +
                      "\nvectors: unavailable for complex eigenvalues";
            return true;
        }
        *output = "values: " + values.to_string() +
                  "\nvectors: " + matrix::eigenvectors(matrix_value).to_string();
        return true;
    } catch (const std::exception&) {
        if (matrix_value.rows == 2 && matrix_value.cols == 2) {
            const double trace = matrix_value.at(0, 0) + matrix_value.at(1, 1);
            const double det = matrix::determinant(matrix_value);
            const double discriminant = trace * trace - 4.0 * det;
            if (discriminant < 0.0) {
                const double real = trace * 0.5;
                const double imag = mymath::sqrt(-discriminant) * 0.5;
                std::ostringstream out;
                out << "values: [complex(" << format_decimal(real) << ", "
                    << format_decimal(imag) << "), complex("
                    << format_decimal(real) << ", " << format_decimal(-imag)
                    << ")]\nvectors: unavailable for complex eigenvalues";
                *output = out.str();
                return true;
            }
        }
        throw;
    }
}

}  // namespace matrix_commands

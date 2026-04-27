// ============================================================================
// 矩阵命令
// ============================================================================

#ifndef CALCULATOR_MATRIX_COMMANDS_H
#define CALCULATOR_MATRIX_COMMANDS_H

#include "matrix.h"

#include <functional>
#include <string>

namespace matrix_commands {

struct MatrixCommandContext {
    std::function<bool(const std::string&)> is_matrix_argument;
    std::function<matrix::Matrix(const std::string&, const std::string&)>
        parse_matrix_argument;
};

bool is_matrix_command(const std::string& command);

bool handle_matrix_command(const MatrixCommandContext& ctx,
                           const std::string& command,
                           const std::string& inside,
                           std::string* output);

}  // namespace matrix_commands

#endif  // CALCULATOR_MATRIX_COMMANDS_H

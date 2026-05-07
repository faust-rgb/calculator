#include "symbolic/commands/symbolic_commands_internal.h"
#include "symbolic/symbolic_expression_internal.h"
#include "core/string_utils.h"
#include <vector>

namespace symbolic_commands {

namespace {
using namespace symbolic_expression_internal;

SymbolicExpression compute_det_recursive(const std::vector<std::vector<SymbolicExpression>>& m) {
    std::size_t size = m.size();
    if (size == 1) return m[0][0];
    if (size == 2) return (m[0][0] * m[1][1] - m[0][1] * m[1][0]).simplify();
    
    SymbolicExpression d = SymbolicExpression::number(0.0L);
    for (std::size_t j = 0; j < size; ++j) {
        std::vector<std::vector<SymbolicExpression>> minor;
        for (std::size_t i = 1; i < size; ++i) {
            std::vector<SymbolicExpression> row;
            for (std::size_t k = 0; k < size; ++k) {
                if (k != j) row.push_back(m[i][k]);
            }
            minor.push_back(row);
        }
        SymbolicExpression term = (m[0][j] * compute_det_recursive(minor)).simplify();
        if (j % 2 == 1) d = (d - term).simplify();
        else d = (d + term).simplify();
    }
    return d;
}

std::string symbolic_matrix_to_semicolon_string(
    const std::vector<std::vector<SymbolicExpression>>& values) {
    std::string out = "[";
    for (std::size_t row = 0; row < values.size(); ++row) {
        if (row != 0) out += "; ";
        out += "[";
        for (std::size_t col = 0; col < values[row].size(); ++col) {
            if (col != 0) out += ", ";
            out += values[row][col].simplify().to_string();
        }
        out += "]";
    }
    out += "]";
    return out;
}
}

bool handle_matrix_commands(const SymbolicCommandContext& ctx,
                           const std::string& command,
                           const std::string& /*inside*/,
                           const std::vector<std::string>& arguments,
                           std::string* output) {
    if (command == "det" || command == "inverse" || command == "transpose") {
        if (arguments.size() != 1) throw std::runtime_error(command + " expects exactly one argument");
        std::string var; SymbolicExpression expr;
        ctx.resolve_symbolic(arguments[0], false, &var, &expr);

        if (!expr.is_tensor() && !expr.is_vector()) throw std::runtime_error(command + " expects a matrix or vector");

        if (command == "transpose") {
            if (expr.is_vector()) {
                auto components = expr.vector_components();
                std::vector<std::vector<SymbolicExpression>> rows;
                for (const auto& comp : components) rows.push_back({comp});
                *output = SymbolicExpression::tensor(rows).simplify().to_string();
            } else {
                auto rows = expr.tensor_rows();
                if (rows.empty()) { *output = expr.to_string(); return true; }
                std::size_t old_rows = rows.size(), old_cols = rows[0].size();
                std::vector<std::vector<SymbolicExpression>> new_rows(old_cols, std::vector<SymbolicExpression>(old_rows));
                for (size_t i = 0; i < old_rows; ++i)
                    for (size_t j = 0; j < old_cols; ++j) new_rows[j][i] = rows[i][j];
                *output = SymbolicExpression::tensor(new_rows).simplify().to_string();
            }
            return true;
        }

        if (!expr.is_tensor()) throw std::runtime_error(command + " expects a matrix");
        auto rows = expr.tensor_rows();
        std::size_t n = rows.size();
        if (n == 0 || rows[0].size() != n) throw std::runtime_error(command + " expects a square matrix");

        if (command == "det") {
            *output = compute_det_recursive(rows).to_string();
            return true;
        }

        if (command == "inverse") {
            SymbolicExpression det = compute_det_recursive(rows);
            if (expr_is_zero(det)) throw std::runtime_error("matrix is singular");
            std::vector<std::vector<SymbolicExpression>> adj(n, std::vector<SymbolicExpression>(n));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    std::vector<std::vector<SymbolicExpression>> minor;
                    for (size_t row = 0; row < n; ++row) {
                        if (row == i) continue;
                        std::vector<SymbolicExpression> minor_row;
                        for (size_t col = 0; col < n; ++col) if (col != j) minor_row.push_back(rows[row][col]);
                        minor.push_back(minor_row);
                    }
                    SymbolicExpression cofactor = compute_det_recursive(minor);
                    if ((i + j) % 2 == 1) cofactor = make_negate(cofactor).simplify();
                    adj[j][i] = cofactor;
                }
            }
            std::vector<std::vector<SymbolicExpression>> inv(n, std::vector<SymbolicExpression>(n));
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j) inv[i][j] = (adj[i][j] / det).simplify();
            *output = SymbolicExpression::tensor(inv).to_string();
            return true;
        }
    }

    if (command == "jacobian" || command == "hessian") {
        if (command == "jacobian") {
            if (arguments.size() < 1) throw std::runtime_error("jacobian expects expression list");
            auto exprs = ctx.parse_symbolic_expression_list(arguments[0]);
            auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {});
            auto jac = SymbolicExpression::jacobian(exprs, vars);
            *output = symbolic_matrix_to_string(jac);
        } else {
            if (arguments.size() < 1) throw std::runtime_error("hessian expects expression");
            std::string v; SymbolicExpression e; ctx.resolve_symbolic(arguments[0], false, &v, &e);
            auto vars = ctx.parse_symbolic_variable_arguments(arguments, 1, {v});
            auto hess = e.hessian(vars);
            *output = symbolic_matrix_to_semicolon_string(hess);
        }
        return true;
    }

    return false;
}

} // namespace symbolic_commands

#include "evaluator.h"

#include "differentiate.h"
#include "function_registry.h"
#include "integrate.h"
#include "simplify.h"

#include <cstddef>
#include <stdexcept>

namespace expression {
namespace {

bool integer_exponent(const numeric::Number& number, int* exponent) {
    const std::string text = number.to_string();
    if (text.empty() || text.find('/') != std::string::npos ||
        text.find('.') != std::string::npos || text.find('i') != std::string::npos) {
        return false;
    }
    long long value = 0;
    std::size_t pos = 0;
    bool negative = false;
    if (text[pos] == '-') {
        negative = true;
        ++pos;
    }
    for (; pos < text.size(); ++pos) {
        if (text[pos] < '0' || text[pos] > '9') {
            return false;
        }
        value = value * 10 + (text[pos] - '0');
        if (value > 1000000) {
            return false;
        }
    }
    *exponent = static_cast<int>(negative ? -value : value);
    return true;
}

runtime::Value make_symbolic(const Expr& expr) {
    return runtime::Value(expr);
}

Expr value_to_expr(const runtime::Value& value) {
    if (value.is_expr()) {
        return value.as_expr();
    }
    if (value.is_number()) {
        return Expr::number(value.as_number());
    }
    if (value.is_matrix()) {
        throw std::runtime_error("matrix value cannot be used as a scalar expression");
    }
    throw std::runtime_error("string value cannot be used as an expression");
}

bool expr_is_number(const Expr& expr) {
    return expr.kind() == ExprKind::Number;
}

bool is_zero_number(const numeric::Number& number) {
    return number.to_string() == "0";
}

std::size_t size_from_number(const numeric::Number& number, const std::string& context) {
    if (!number.is_integer()) {
        throw std::runtime_error(context + " requires integer dimensions");
    }
    const std::string text = number.as_integer().to_string();
    if (text.empty() || text[0] == '-') {
        throw std::runtime_error(context + " requires non-negative dimensions");
    }
    std::size_t value = 0;
    for (char ch : text) {
        value = value * 10 + static_cast<std::size_t>(ch - '0');
        if (value > 10000) {
            throw std::runtime_error(context + " dimension is too large");
        }
    }
    return value;
}

Expr add_expr(const Expr& lhs, const Expr& rhs) {
    if (expr_is_number(lhs) && expr_is_number(rhs)) {
        return Expr::number(numeric::add(lhs.number_value(), rhs.number_value()));
    }
    return cas::simplify(Expr::add(lhs, rhs));
}

Expr negate_expr(const Expr& expr) {
    return add_expr(Expr::number(numeric::Number(0)),
                    Expr::mul(Expr::number(numeric::Number(-1)), expr));
}

Expr subtract_expr(const Expr& lhs, const Expr& rhs) {
    if (expr_is_number(lhs) && expr_is_number(rhs)) {
        return Expr::number(numeric::subtract(lhs.number_value(), rhs.number_value()));
    }
    return add_expr(lhs, Expr::mul(Expr::number(numeric::Number(-1)), rhs));
}

Expr multiply_expr(const Expr& lhs, const Expr& rhs) {
    if (expr_is_number(lhs) && expr_is_number(rhs)) {
        return Expr::number(numeric::multiply(lhs.number_value(), rhs.number_value()));
    }
    return cas::simplify(Expr::mul(lhs, rhs));
}

runtime::MatrixValue matrix_from_values(const std::vector<runtime::Value>& args,
                                        std::size_t start,
                                        std::size_t rows,
                                        std::size_t cols) {
    if (args.size() - start != rows * cols) {
        throw std::runtime_error("mat entry count does not match dimensions");
    }
    std::vector<std::vector<Expr>> entries(rows, std::vector<Expr>(cols));
    std::size_t index = start;
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
            entries[row][col] = value_to_expr(args[index++]);
        }
    }
    return runtime::MatrixValue(entries);
}

runtime::MatrixValue add_matrix(const runtime::MatrixValue& lhs,
                                const runtime::MatrixValue& rhs) {
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
        throw std::runtime_error("matrix dimensions must agree");
    }
    runtime::MatrixValue result(lhs.rows(), lhs.cols(), Expr::number(numeric::Number(0)));
    for (std::size_t row = 0; row < lhs.rows(); ++row) {
        for (std::size_t col = 0; col < lhs.cols(); ++col) {
            result.at(row, col) = add_expr(lhs.at(row, col), rhs.at(row, col));
        }
    }
    return result;
}

runtime::MatrixValue multiply_matrix(const runtime::MatrixValue& lhs,
                                     const runtime::MatrixValue& rhs) {
    if (lhs.cols() != rhs.rows()) {
        throw std::runtime_error("matrix inner dimensions must agree");
    }
    runtime::MatrixValue result(lhs.rows(), rhs.cols(), Expr::number(numeric::Number(0)));
    for (std::size_t row = 0; row < lhs.rows(); ++row) {
        for (std::size_t col = 0; col < rhs.cols(); ++col) {
            Expr total = Expr::number(numeric::Number(0));
            for (std::size_t k = 0; k < lhs.cols(); ++k) {
                total = add_expr(total,
                                 multiply_expr(lhs.at(row, k), rhs.at(k, col)));
            }
            result.at(row, col) = total;
        }
    }
    return result;
}

runtime::MatrixValue scale_matrix(const Expr& scalar,
                                  const runtime::MatrixValue& matrix) {
    runtime::MatrixValue result(matrix.rows(), matrix.cols(), Expr::number(numeric::Number(0)));
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            result.at(row, col) = multiply_expr(scalar, matrix.at(row, col));
        }
    }
    return result;
}

Expr determinant(const runtime::MatrixValue& matrix) {
    if (matrix.rows() != matrix.cols()) {
        throw std::runtime_error("det expects a square matrix");
    }
    if (matrix.rows() == 0) {
        return Expr::number(numeric::Number(1));
    }
    if (matrix.rows() == 1) {
        return matrix.at(0, 0);
    }
    if (matrix.rows() == 2) {
        return subtract_expr(multiply_expr(matrix.at(0, 0), matrix.at(1, 1)),
                             multiply_expr(matrix.at(0, 1), matrix.at(1, 0)));
    }
    Expr total = Expr::number(numeric::Number(0));
    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        std::vector<std::vector<Expr>> minor_entries;
        minor_entries.reserve(matrix.rows() - 1);
        for (std::size_t row = 1; row < matrix.rows(); ++row) {
            std::vector<Expr> minor_row;
            minor_row.reserve(matrix.cols() - 1);
            for (std::size_t minor_col = 0; minor_col < matrix.cols(); ++minor_col) {
                if (minor_col != col) {
                    minor_row.push_back(matrix.at(row, minor_col));
                }
            }
            minor_entries.push_back(minor_row);
        }
        Expr term = multiply_expr(matrix.at(0, col),
                                  determinant(runtime::MatrixValue(minor_entries)));
        if (col % 2 == 1) {
            term = negate_expr(term);
        }
        total = add_expr(total, term);
    }
    return total;
}

runtime::MatrixValue inverse_matrix(const runtime::MatrixValue& matrix) {
    if (matrix.rows() != matrix.cols()) {
        throw std::runtime_error("inverse expects a square matrix");
    }
    const std::size_t n = matrix.rows();
    std::vector<std::vector<numeric::Number>> augmented(
        n, std::vector<numeric::Number>(2 * n, numeric::Number(0)));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            if (!expr_is_number(matrix.at(row, col))) {
                throw std::runtime_error("inverse currently requires numeric matrix entries");
            }
            augmented[row][col] = matrix.at(row, col).number_value();
        }
        augmented[row][n + row] = numeric::Number(1);
    }
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        while (pivot < n && is_zero_number(augmented[pivot][col])) {
            ++pivot;
        }
        if (pivot == n) {
            throw std::runtime_error("matrix is singular");
        }
        if (pivot != col) {
            std::swap(augmented[pivot], augmented[col]);
        }
        const numeric::Number pivot_value = augmented[col][col];
        for (std::size_t j = 0; j < 2 * n; ++j) {
            augmented[col][j] = numeric::divide(augmented[col][j], pivot_value);
        }
        for (std::size_t row = 0; row < n; ++row) {
            if (row == col || is_zero_number(augmented[row][col])) {
                continue;
            }
            const numeric::Number factor = augmented[row][col];
            for (std::size_t j = 0; j < 2 * n; ++j) {
                augmented[row][j] = numeric::subtract(
                    augmented[row][j],
                    numeric::multiply(factor, augmented[col][j]));
            }
        }
    }
    std::vector<std::vector<Expr>> entries(n, std::vector<Expr>(n));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            entries[row][col] = Expr::number(augmented[row][n + col]);
        }
    }
    return runtime::MatrixValue(entries);
}

runtime::MatrixValue rref_matrix(const runtime::MatrixValue& matrix) {
    std::vector<std::vector<Expr>> entries = matrix.entries();
    std::size_t lead = 0;
    for (std::size_t row = 0; row < matrix.rows() && lead < matrix.cols(); ++row) {
        std::size_t pivot = row;
        while (pivot < matrix.rows()) {
            if (!expr_is_number(entries[pivot][lead])) {
                throw std::runtime_error("rref currently requires numeric matrix entries");
            }
            if (!is_zero_number(entries[pivot][lead].number_value())) {
                break;
            }
            ++pivot;
        }
        if (pivot == matrix.rows()) {
            ++lead;
            --row;
            continue;
        }
        if (pivot != row) {
            std::swap(entries[pivot], entries[row]);
        }
        const numeric::Number pivot_value = entries[row][lead].number_value();
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            entries[row][col] = Expr::number(numeric::divide(entries[row][col].number_value(),
                                                             pivot_value));
        }
        for (std::size_t other = 0; other < matrix.rows(); ++other) {
            if (other == row) {
                continue;
            }
            const numeric::Number factor = entries[other][lead].number_value();
            if (is_zero_number(factor)) {
                continue;
            }
            for (std::size_t col = 0; col < matrix.cols(); ++col) {
                entries[other][col] = Expr::number(numeric::subtract(
                    entries[other][col].number_value(),
                    numeric::multiply(factor, entries[row][col].number_value())));
            }
        }
        ++lead;
    }
    return runtime::MatrixValue(entries);
}

std::vector<std::string> variables_from_list_expr(const Expr& expr) {
    if (expr.kind() != ExprKind::List) {
        throw std::runtime_error("expected a variable list");
    }
    std::vector<std::string> variables;
    variables.reserve(expr.children().size());
    for (const Expr& child : expr.children()) {
        if (child.kind() != ExprKind::Symbol) {
            throw std::runtime_error("variable list entries must be symbols");
        }
        variables.push_back(child.symbol_name());
    }
    return variables;
}

}  // namespace

runtime::Value evaluate(const Expr& expr, runtime::Environment& env) {
    switch (expr.kind()) {
        case ExprKind::Number:
            return runtime::Value(expr.number_value());
        case ExprKind::Symbol: {
            runtime::Value value;
            if (env.lookup(expr.symbol_name(), &value)) {
                return value;
            }
            if (expr.symbol_name() == "i") {
                return runtime::Value(numeric::Number(numeric::Complex(
                    numeric::BigDecimal(), numeric::BigDecimal(numeric::BigInt(1)))));
            }
            if (expr.symbol_name() == "pi" || expr.symbol_name() == "e") {
                return make_symbolic(expr);
            }
            return make_symbolic(expr);
        }
        case ExprKind::Add: {
            const runtime::Value lhs = evaluate(expr.children()[0], env);
            const runtime::Value rhs = evaluate(expr.children()[1], env);
            if (lhs.is_number() && rhs.is_number()) {
                return runtime::Value(numeric::add(lhs.as_number(), rhs.as_number()));
            }
            if (lhs.is_matrix() && rhs.is_matrix()) {
                return runtime::Value(add_matrix(lhs.as_matrix(), rhs.as_matrix()));
            }
            if (lhs.is_matrix() || rhs.is_matrix()) {
                throw std::runtime_error("matrix addition requires two matrices");
            }
            return make_symbolic(Expr::add(lhs.is_expr() ? lhs.as_expr()
                                                         : Expr::number(lhs.as_number()),
                                          rhs.is_expr() ? rhs.as_expr()
                                                        : Expr::number(rhs.as_number())));
        }
        case ExprKind::Mul: {
            const runtime::Value lhs = evaluate(expr.children()[0], env);
            const runtime::Value rhs = evaluate(expr.children()[1], env);
            if (lhs.is_number() && rhs.is_number()) {
                return runtime::Value(numeric::multiply(lhs.as_number(), rhs.as_number()));
            }
            if (lhs.is_matrix() && rhs.is_matrix()) {
                return runtime::Value(multiply_matrix(lhs.as_matrix(), rhs.as_matrix()));
            }
            if (lhs.is_number() && rhs.is_matrix()) {
                return runtime::Value(scale_matrix(Expr::number(lhs.as_number()),
                                                   rhs.as_matrix()));
            }
            if (lhs.is_expr() && rhs.is_matrix()) {
                return runtime::Value(scale_matrix(lhs.as_expr(), rhs.as_matrix()));
            }
            if (lhs.is_matrix() && rhs.is_number()) {
                return runtime::Value(scale_matrix(Expr::number(rhs.as_number()),
                                                   lhs.as_matrix()));
            }
            if (lhs.is_matrix() && rhs.is_expr()) {
                return runtime::Value(scale_matrix(rhs.as_expr(), lhs.as_matrix()));
            }
            if (lhs.is_matrix() || rhs.is_matrix()) {
                throw std::runtime_error("unsupported matrix multiplication operands");
            }
            return make_symbolic(Expr::mul(lhs.is_expr() ? lhs.as_expr()
                                                         : Expr::number(lhs.as_number()),
                                          rhs.is_expr() ? rhs.as_expr()
                                                        : Expr::number(rhs.as_number())));
        }
        case ExprKind::Pow: {
            const runtime::Value base = evaluate(expr.children()[0], env);
            const runtime::Value exponent = evaluate(expr.children()[1], env);
            if (base.is_number() && exponent.is_number()) {
                int exp = 0;
                if (integer_exponent(exponent.as_number(), &exp)) {
                    numeric::Number result(1);
                    const int count = exp < 0 ? -exp : exp;
                    for (int i = 0; i < count; ++i) {
                        result = numeric::multiply(result, base.as_number());
                    }
                    if (exp < 0) {
                        result = numeric::divide(numeric::Number(1), result);
                    }
                    return runtime::Value(result);
                }
            }
            return make_symbolic(Expr::pow(base.is_expr() ? base.as_expr()
                                                          : Expr::number(base.as_number()),
                                          exponent.is_expr() ? exponent.as_expr()
                                                             : Expr::number(exponent.as_number())));
        }
        case ExprKind::Function:
        {
            std::vector<runtime::Value> evaluated_args;
            std::vector<numeric::Number> numeric_args;
            bool all_numeric = true;
            for (const Expr& child : expr.children()) {
                runtime::Value value = evaluate(child, env);
                if (!value.is_number()) {
                    all_numeric = false;
                } else {
                    numeric_args.push_back(value.as_number());
                }
                evaluated_args.push_back(value);
            }
            if (expr.function_name() == "mat") {
                if (evaluated_args.size() < 2 ||
                    !evaluated_args[0].is_number() ||
                    !evaluated_args[1].is_number()) {
                    throw std::runtime_error("mat expects row count, column count, and entries");
                }
                const std::size_t rows =
                    size_from_number(evaluated_args[0].as_number(), "mat");
                const std::size_t cols =
                    size_from_number(evaluated_args[1].as_number(), "mat");
                return runtime::Value(matrix_from_values(evaluated_args, 2, rows, cols));
            }
            if (expr.function_name() == "det") {
                if (evaluated_args.size() != 1 || !evaluated_args[0].is_matrix()) {
                    throw std::runtime_error("det expects exactly one matrix argument");
                }
                const Expr result = determinant(evaluated_args[0].as_matrix());
                if (result.kind() == ExprKind::Number) {
                    return runtime::Value(result.number_value());
                }
                return make_symbolic(result);
            }
            if (expr.function_name() == "inverse") {
                if (evaluated_args.size() != 1 || !evaluated_args[0].is_matrix()) {
                    throw std::runtime_error("inverse expects exactly one matrix argument");
                }
                return runtime::Value(inverse_matrix(evaluated_args[0].as_matrix()));
            }
            if (expr.function_name() == "rref") {
                if (evaluated_args.size() != 1 || !evaluated_args[0].is_matrix()) {
                    throw std::runtime_error("rref expects exactly one matrix argument");
                }
                return runtime::Value(rref_matrix(evaluated_args[0].as_matrix()));
            }
            if (expr.function_name() == "simplify") {
                if (evaluated_args.size() != 1) {
                    throw std::runtime_error("simplify expects exactly one argument");
                }
                return make_symbolic(cas::simplify(value_to_expr(evaluated_args[0])));
            }
            if (expr.function_name() == "diff") {
                if (evaluated_args.size() != 2) {
                    throw std::runtime_error("diff expects exactly two arguments");
                }
                const Expr variable = value_to_expr(evaluated_args[1]);
                if (variable.kind() != ExprKind::Symbol) {
                    throw std::runtime_error("diff variable must be a symbol");
                }
                return make_symbolic(cas::differentiate(value_to_expr(evaluated_args[0]),
                                                        variable.symbol_name()));
            }
            if (expr.function_name() == "integrate") {
                if (evaluated_args.size() != 2) {
                    throw std::runtime_error("integrate expects exactly two arguments");
                }
                const Expr variable = value_to_expr(evaluated_args[1]);
                if (variable.kind() != ExprKind::Symbol) {
                    throw std::runtime_error("integrate variable must be a symbol");
                }
                return make_symbolic(cas::integrate(value_to_expr(evaluated_args[0]),
                                                    variable.symbol_name()));
            }
            if (expr.function_name() == "gradient") {
                if (evaluated_args.size() != 2) {
                    throw std::runtime_error("gradient expects expression and variable list");
                }
                const Expr body = value_to_expr(evaluated_args[0]);
                const std::vector<std::string> variables =
                    variables_from_list_expr(value_to_expr(evaluated_args[1]));
                std::vector<Expr> entries;
                entries.reserve(variables.size());
                for (const std::string& variable : variables) {
                    entries.push_back(cas::differentiate(body, variable));
                }
                return make_symbolic(Expr::list(entries));
            }
            if (expr.function_name() == "jacobian") {
                if (evaluated_args.size() != 2) {
                    throw std::runtime_error("jacobian expects expression list and variable list");
                }
                const Expr functions = value_to_expr(evaluated_args[0]);
                if (functions.kind() != ExprKind::List) {
                    throw std::runtime_error("jacobian first argument must be a list");
                }
                const std::vector<std::string> variables =
                    variables_from_list_expr(value_to_expr(evaluated_args[1]));
                std::vector<Expr> rows;
                rows.reserve(functions.children().size());
                for (const Expr& function : functions.children()) {
                    std::vector<Expr> row;
                    row.reserve(variables.size());
                    for (const std::string& variable : variables) {
                        row.push_back(cas::differentiate(function, variable));
                    }
                    rows.push_back(Expr::list(row));
                }
                return make_symbolic(Expr::list(rows));
            }
            if (expr.function_name() == "hessian") {
                if (evaluated_args.size() != 2) {
                    throw std::runtime_error("hessian expects expression and variable list");
                }
                const Expr body = value_to_expr(evaluated_args[0]);
                const std::vector<std::string> variables =
                    variables_from_list_expr(value_to_expr(evaluated_args[1]));
                std::vector<Expr> rows;
                rows.reserve(variables.size());
                for (const std::string& row_variable : variables) {
                    std::vector<Expr> row;
                    row.reserve(variables.size());
                    const Expr first = cas::differentiate(body, row_variable);
                    for (const std::string& col_variable : variables) {
                        row.push_back(cas::differentiate(first, col_variable));
                    }
                    rows.push_back(Expr::list(row));
                }
                return make_symbolic(Expr::list(rows));
            }
            if (all_numeric && runtime::FunctionRegistry::builtins().find(expr.function_name())) {
                try {
                    return runtime::Value(runtime::FunctionRegistry::builtins().evaluate_numeric(
                        expr.function_name(), numeric_args, env.precision()));
                } catch (const std::runtime_error&) {
                    // Unsupported exact cases intentionally remain symbolic while
                    // the 2.0 function surface is being migrated.
                }
            }
            std::vector<Expr> symbolic_args;
            symbolic_args.reserve(evaluated_args.size());
            for (const runtime::Value& value : evaluated_args) {
                if (value.is_matrix()) {
                    throw std::runtime_error("function does not accept matrix arguments yet");
                }
                symbolic_args.push_back(value.is_expr() ? value.as_expr()
                                                        : Expr::number(value.as_number()));
            }
            return make_symbolic(Expr::function(expr.function_name(), symbolic_args));
        }
        case ExprKind::List:
        {
            std::vector<runtime::Value> values;
            values.reserve(expr.children().size());
            for (const Expr& child : expr.children()) {
                values.push_back(evaluate(child, env));
            }
            bool matrix_literal = !values.empty();
            std::size_t matrix_cols = 0;
            for (const runtime::Value& value : values) {
                if (!value.is_expr() || value.as_expr().kind() != ExprKind::List) {
                    matrix_literal = false;
                    break;
                }
                if (matrix_cols == 0) {
                    matrix_cols = value.as_expr().children().size();
                } else if (value.as_expr().children().size() != matrix_cols) {
                    matrix_literal = false;
                    break;
                }
            }
            if (matrix_literal && matrix_cols != 0) {
                std::vector<std::vector<Expr>> entries;
                entries.reserve(values.size());
                for (const runtime::Value& value : values) {
                    entries.push_back(value.as_expr().children());
                }
                return runtime::Value(runtime::MatrixValue(entries));
            }
            std::vector<Expr> items;
            items.reserve(values.size());
            for (const runtime::Value& value : values) {
                items.push_back(value_to_expr(value));
            }
            return make_symbolic(Expr::list(items));
        }
        case ExprKind::Integral:
            return make_symbolic(expr);
    }
    throw std::runtime_error("unknown expression kind");
}

}  // namespace expression

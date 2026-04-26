#include "evaluator.h"

#include "differentiate.h"
#include "function_registry.h"
#include "integrate.h"
#include "simplify.h"

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
    throw std::runtime_error("string value cannot be used as an expression");
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
                return runtime::Value(numeric::add(lhs.as_number(), rhs.as_number(),
                                                   env.precision()));
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
                return runtime::Value(numeric::multiply(lhs.as_number(), rhs.as_number(),
                                                        env.precision()));
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
                        result = numeric::multiply(result, base.as_number(), env.precision());
                    }
                    if (exp < 0) {
                        result = numeric::divide(numeric::Number(1), result, env.precision());
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
                symbolic_args.push_back(value.is_expr() ? value.as_expr()
                                                        : Expr::number(value.as_number()));
            }
            return make_symbolic(Expr::function(expr.function_name(), symbolic_args));
        }
        case ExprKind::List:
        {
            std::vector<Expr> items;
            items.reserve(expr.children().size());
            for (const Expr& child : expr.children()) {
                items.push_back(value_to_expr(evaluate(child, env)));
            }
            return make_symbolic(Expr::list(items));
        }
        case ExprKind::Integral:
            return make_symbolic(expr);
    }
    throw std::runtime_error("unknown expression kind");
}

}  // namespace expression

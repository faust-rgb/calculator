#include "differentiate.h"

#include "printer.h"
#include "simplify.h"

namespace cas {
namespace {

expression::Expr number(long long value) {
    return expression::Expr::number(numeric::Number(value));
}

expression::Expr add(const expression::Expr& lhs, const expression::Expr& rhs) {
    return expression::Expr::add(lhs, rhs);
}

expression::Expr sub(const expression::Expr& lhs, const expression::Expr& rhs) {
    return add(lhs, expression::Expr::mul(number(-1), rhs));
}

expression::Expr mul(const expression::Expr& lhs, const expression::Expr& rhs) {
    return expression::Expr::mul(lhs, rhs);
}

expression::Expr div(const expression::Expr& lhs, const expression::Expr& rhs) {
    return mul(lhs, expression::Expr::pow(rhs, number(-1)));
}

expression::Expr fn(const std::string& name, const expression::Expr& arg) {
    return expression::Expr::function(name, {arg});
}

bool integer_number(const expression::Expr& expr, int* value) {
    if (expr.kind() != expression::ExprKind::Number) {
        return false;
    }
    const std::string text = expr.number_value().to_string();
    if (text.empty() || text.find('/') != std::string::npos ||
        text.find('.') != std::string::npos || text.find('i') != std::string::npos) {
        return false;
    }
    int parsed = 0;
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
        parsed = parsed * 10 + (text[pos] - '0');
    }
    *value = negative ? -parsed : parsed;
    return true;
}

}  // namespace

expression::Expr differentiate(const expression::Expr& expr, const std::string& variable_name) {
    switch (expr.kind()) {
        case expression::ExprKind::Number:
            return number(0);
        case expression::ExprKind::Symbol:
            return number(expr.symbol_name() == variable_name ? 1 : 0);
        case expression::ExprKind::Add:
            return simplify(add(differentiate(expr.children()[0], variable_name),
                                differentiate(expr.children()[1], variable_name)));
        case expression::ExprKind::Mul: {
            const expression::Expr lhs = expr.children()[0];
            const expression::Expr rhs = expr.children()[1];
            return simplify(add(mul(differentiate(lhs, variable_name), rhs),
                                mul(lhs, differentiate(rhs, variable_name))));
        }
        case expression::ExprKind::Pow: {
            const expression::Expr base = expr.children()[0];
            const expression::Expr exponent = expr.children()[1];
            int n = 0;
            if (integer_number(exponent, &n)) {
                return simplify(mul(mul(number(n),
                                        expression::Expr::pow(base, number(n - 1))),
                                    differentiate(base, variable_name)));
            }
            // General power rule: (a^b)' = a^b * (b' ln(a) + b a'/a).
            return simplify(mul(expression::Expr::pow(base, exponent),
                                add(mul(differentiate(exponent, variable_name), fn("ln", base)),
                                    mul(exponent,
                                        div(differentiate(base, variable_name), base)))));
        }
        case expression::ExprKind::Function: {
            if (expr.children().empty()) {
                return expression::Expr::function("Derivative",
                                                 {expr, expression::Expr::symbol(variable_name)});
            }
            const expression::Expr arg = expr.children()[0];
            const expression::Expr d_arg = differentiate(arg, variable_name);
            const std::string& name = expr.function_name();
            if (name == "sin") {
                return simplify(mul(fn("cos", arg), d_arg));
            }
            if (name == "cos") {
                return simplify(mul(mul(number(-1), fn("sin", arg)), d_arg));
            }
            if (name == "tan") {
                return simplify(mul(expression::Expr::pow(fn("sec", arg), number(2)), d_arg));
            }
            if (name == "cot") {
                return simplify(mul(mul(number(-1),
                                        expression::Expr::pow(fn("csc", arg), number(2))),
                                    d_arg));
            }
            if (name == "sec") {
                return simplify(mul(mul(fn("sec", arg), fn("tan", arg)), d_arg));
            }
            if (name == "csc") {
                return simplify(mul(mul(number(-1),
                                        mul(fn("csc", arg), fn("cot", arg))),
                                    d_arg));
            }
            if (name == "exp") {
                return simplify(mul(fn("exp", arg), d_arg));
            }
            if (name == "ln" || name == "log") {
                if (arg.kind() == expression::ExprKind::Function &&
                    arg.function_name() == "sin" &&
                    arg.children().size() == 1) {
                    return simplify(mul(fn("cot", arg.children()[0]),
                                        differentiate(arg.children()[0], variable_name)));
                }
                return simplify(div(d_arg, arg));
            }
            if (name == "sqrt") {
                return simplify(div(d_arg, mul(number(2), fn("sqrt", arg))));
            }
            if (name == "sinh") {
                return simplify(mul(fn("cosh", arg), d_arg));
            }
            if (name == "cosh") {
                return simplify(mul(fn("sinh", arg), d_arg));
            }
            if (name == "tanh") {
                return simplify(mul(expression::Expr::pow(fn("sech", arg), number(2)), d_arg));
            }
            if (name == "asinh") {
                return simplify(div(d_arg, fn("sqrt", add(expression::Expr::pow(arg, number(2)),
                                                          number(1)))));
            }
            if (name == "acosh") {
                return simplify(div(d_arg,
                                    mul(fn("sqrt", sub(arg, number(1))),
                                        fn("sqrt", add(arg, number(1))))));
            }
            if (name == "atanh") {
                return simplify(div(d_arg, sub(number(1), expression::Expr::pow(arg, number(2)))));
            }
            if (name == "asin") {
                return simplify(div(d_arg, fn("sqrt", sub(number(1), expression::Expr::pow(arg, number(2))))));
            }
            if (name == "acos") {
                return simplify(mul(number(-1),
                                    div(d_arg,
                                        fn("sqrt", sub(number(1),
                                                       expression::Expr::pow(arg, number(2)))))));
            }
            if (name == "atan") {
                return simplify(div(d_arg, add(number(1), expression::Expr::pow(arg, number(2)))));
            }
            if (name == "abs") {
                return simplify(mul(fn("sign", arg), d_arg));
            }
            if (name == "step" || name == "heaviside") {
                return simplify(mul(fn("delta", arg), d_arg));
            }
            if (name == "sign") {
                return expression::Expr::function("Derivative",
                                                 {expr, expression::Expr::symbol(variable_name)});
            }
            return expression::Expr::function("Derivative",
                                             {expr, expression::Expr::symbol(variable_name)});
        }
        case expression::ExprKind::List: {
            std::vector<expression::Expr> items;
            items.reserve(expr.children().size());
            for (const expression::Expr& child : expr.children()) {
                items.push_back(differentiate(child, variable_name));
            }
            return expression::Expr::list(items);
        }
        case expression::ExprKind::Integral:
            return expression::Expr::function("Derivative",
                                             {expr, expression::Expr::symbol(variable_name)});
    }
    return number(0);
}

}  // namespace cas

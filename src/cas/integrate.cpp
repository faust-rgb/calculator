#include "integrate.h"

#include "differentiate.h"
#include "printer.h"
#include "simplify.h"

#include <optional>
#include <string>
#include <vector>

namespace cas {
namespace {

using expression::Expr;
using expression::ExprKind;

Expr number(long long value) {
    return Expr::number(numeric::Number(value));
}

Expr rational(long long numerator, long long denominator) {
    return Expr::number(numeric::Number(numeric::Rational(numerator, denominator)));
}

Expr symbol(const std::string& name) {
    return Expr::symbol(name);
}

Expr add(const Expr& lhs, const Expr& rhs) {
    return Expr::add(lhs, rhs);
}

Expr sub(const Expr& lhs, const Expr& rhs) {
    return add(lhs, Expr::mul(number(-1), rhs));
}

Expr mul(const Expr& lhs, const Expr& rhs) {
    return Expr::mul(lhs, rhs);
}

Expr div(const Expr& lhs, const Expr& rhs) {
    return mul(lhs, Expr::pow(rhs, number(-1)));
}

Expr pow(const Expr& base, const Expr& exponent) {
    return Expr::pow(base, exponent);
}

Expr fn(const std::string& name, const Expr& arg) {
    return Expr::function(name, {arg});
}

bool structurally_equal(const Expr& lhs, const Expr& rhs) {
    return expression::print(simplify(lhs)) == expression::print(simplify(rhs));
}

bool is_number_text(const Expr& expr, const std::string& text) {
    return expr.kind() == ExprKind::Number && expr.number_value().to_string() == text;
}

bool is_zero(const Expr& expr) {
    return is_number_text(expr, "0");
}

bool is_one(const Expr& expr) {
    return is_number_text(expr, "1");
}

bool is_minus_one(const Expr& expr) {
    return is_number_text(expr, "-1");
}

Expr expr_from_number(const numeric::Number& value) {
    return Expr::number(value);
}

bool contains_variable(const Expr& expr, const std::string& variable_name) {
    switch (expr.kind()) {
        case ExprKind::Number:
            return false;
        case ExprKind::Symbol:
            return expr.symbol_name() == variable_name;
        case ExprKind::Add:
        case ExprKind::Mul:
        case ExprKind::Pow:
        case ExprKind::Function:
        case ExprKind::Integral:
        case ExprKind::List:
            for (const Expr& child : expr.children()) {
                if (contains_variable(child, variable_name)) {
                    return true;
                }
            }
            return false;
    }
    return false;
}

bool integer_number(const Expr& expr, int* value) {
    if (expr.kind() != ExprKind::Number) {
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

void append_mul_factors(const Expr& expr, std::vector<Expr>* factors) {
    if (expr.kind() == ExprKind::Mul && expr.children().size() == 2) {
        append_mul_factors(expr.children()[0], factors);
        append_mul_factors(expr.children()[1], factors);
        return;
    }
    factors->push_back(expr);
}

std::optional<Expr> integrate_internal(const Expr& expr, const std::string& variable_name);
bool is_variable(const Expr& expr, const std::string& variable_name);

Expr monomial(const std::string& variable_name, int power_value) {
    if (power_value == 0) {
        return number(1);
    }
    if (power_value == 1) {
        return symbol(variable_name);
    }
    return pow(symbol(variable_name), number(power_value));
}

std::optional<int> monomial_power(const Expr& expr, const std::string& variable_name) {
    if (is_variable(expr, variable_name)) {
        return 1;
    }
    if (expr.kind() == ExprKind::Pow &&
        expr.children().size() == 2 &&
        is_variable(expr.children()[0], variable_name)) {
        int power_value = 0;
        if (integer_number(expr.children()[1], &power_value) && power_value >= 0) {
            return power_value;
        }
    }
    return std::nullopt;
}

Expr integrate_monomial_times_function(const std::string& function_name,
                                       int power_value,
                                       const std::string& variable_name);

Expr integrate_monomial_times_exp(int power_value, const std::string& variable_name) {
    const Expr x = symbol(variable_name);
    if (power_value == 0) {
        return fn("exp", x);
    }
    return sub(mul(monomial(variable_name, power_value), fn("exp", x)),
               mul(number(power_value),
                   integrate_monomial_times_exp(power_value - 1, variable_name)));
}

Expr integrate_monomial_times_sin(int power_value, const std::string& variable_name) {
    const Expr x = symbol(variable_name);
    if (power_value == 0) {
        return mul(number(-1), fn("cos", x));
    }
    return add(mul(number(-1), mul(monomial(variable_name, power_value), fn("cos", x))),
               mul(number(power_value),
                   integrate_monomial_times_function("cos", power_value - 1, variable_name)));
}

Expr integrate_monomial_times_cos(int power_value, const std::string& variable_name) {
    const Expr x = symbol(variable_name);
    if (power_value == 0) {
        return fn("sin", x);
    }
    return sub(mul(monomial(variable_name, power_value), fn("sin", x)),
               mul(number(power_value),
                   integrate_monomial_times_function("sin", power_value - 1, variable_name)));
}

Expr integrate_monomial_times_function(const std::string& function_name,
                                       int power_value,
                                       const std::string& variable_name) {
    if (function_name == "exp") {
        return integrate_monomial_times_exp(power_value, variable_name);
    }
    if (function_name == "sin") {
        return integrate_monomial_times_sin(power_value, variable_name);
    }
    return integrate_monomial_times_cos(power_value, variable_name);
}

std::optional<Expr> integrate_polynomial_times_known_function(
    const Expr& expr,
    const std::string& variable_name) {
    std::vector<Expr> factors;
    append_mul_factors(expr, &factors);
    numeric::Number coefficient(1);
    std::optional<int> power_value;
    std::optional<std::string> function_name;

    for (const Expr& factor : factors) {
        if (factor.kind() == ExprKind::Number) {
            coefficient = numeric::multiply(coefficient, factor.number_value());
            continue;
        }
        if (const std::optional<int> n = monomial_power(factor, variable_name)) {
            if (power_value) {
                return std::nullopt;
            }
            power_value = *n;
            continue;
        }
        if (factor.kind() == ExprKind::Function &&
            factor.children().size() == 1 &&
            is_variable(factor.children()[0], variable_name) &&
            (factor.function_name() == "exp" ||
             factor.function_name() == "sin" ||
             factor.function_name() == "cos")) {
            if (function_name) {
                return std::nullopt;
            }
            function_name = factor.function_name();
            continue;
        }
        return std::nullopt;
    }

    if (!power_value || !function_name || *power_value == 0) {
        return std::nullopt;
    }
    Expr result = integrate_monomial_times_function(*function_name,
                                                    *power_value,
                                                    variable_name);
    if (coefficient.to_string() != "1") {
        result = mul(expr_from_number(coefficient), result);
    }
    return result;
}

std::optional<Expr> integrate_chain_rule(const Expr& expr,
                                         const std::string& variable_name) {
    std::vector<Expr> factors;
    append_mul_factors(expr, &factors);
    if (factors.size() < 2) {
        return std::nullopt;
    }

    for (std::size_t i = 0; i < factors.size(); ++i) {
        const Expr candidate = factors[i];
        if (candidate.kind() != ExprKind::Function || candidate.children().size() != 1) {
            continue;
        }
        const Expr inner = candidate.children()[0];
        const Expr derivative = simplify(differentiate(inner, variable_name));
        std::vector<Expr> rest;
        for (std::size_t j = 0; j < factors.size(); ++j) {
            if (j != i) {
                rest.push_back(factors[j]);
            }
        }
        Expr rest_expr = rest.empty() ? number(1) : rest[0];
        for (std::size_t j = 1; j < rest.size(); ++j) {
            rest_expr = mul(rest_expr, rest[j]);
        }
        if (!structurally_equal(rest_expr, derivative)) {
            continue;
        }
        if (candidate.function_name() == "cos") {
            return fn("sin", inner);
        }
        if (candidate.function_name() == "sin") {
            return mul(number(-1), fn("cos", inner));
        }
        if (candidate.function_name() == "exp") {
            return fn("exp", inner);
        }
    }
    return std::nullopt;
}

bool is_variable(const Expr& expr, const std::string& variable_name) {
    return expr.kind() == ExprKind::Symbol && expr.symbol_name() == variable_name;
}

bool is_x_squared(const Expr& expr, const std::string& variable_name) {
    return expr.kind() == ExprKind::Pow &&
           expr.children().size() == 2 &&
           is_variable(expr.children()[0], variable_name) &&
           is_number_text(expr.children()[1], "2");
}

bool is_negative_x_squared(const Expr& expr, const std::string& variable_name) {
    if (expr.kind() != ExprKind::Mul || expr.children().size() != 2) {
        return false;
    }
    return (is_minus_one(expr.children()[0]) && is_x_squared(expr.children()[1], variable_name)) ||
           (is_minus_one(expr.children()[1]) && is_x_squared(expr.children()[0], variable_name));
}

std::optional<Expr> constant_derivative(const Expr& arg,
                                        const std::string& variable_name) {
    const Expr derivative = simplify(differentiate(arg, variable_name));
    if (!contains_variable(derivative, variable_name) && !is_zero(derivative)) {
        return derivative;
    }
    return std::nullopt;
}

std::optional<Expr> integrate_function(const Expr& expr, const std::string& variable_name) {
    if (expr.children().size() != 1) {
        return std::nullopt;
    }
    const std::string& name = expr.function_name();
    const Expr arg = expr.children()[0];
    if (is_variable(arg, variable_name)) {
        if (name == "sin") {
            return mul(number(-1), fn("cos", arg));
        }
        if (name == "cos") {
            return fn("sin", arg);
        }
        if (name == "sec") {
            return fn("ln", add(fn("sec", arg), fn("tan", arg)));
        }
        if (name == "csc") {
            return mul(number(-1), fn("ln", add(fn("csc", arg), fn("cot", arg))));
        }
        if (name == "exp") {
            return fn("exp", arg);
        }
        if (name == "tan") {
            return mul(number(-1), fn("ln", fn("cos", arg)));
        }
        if (name == "cot") {
            return fn("ln", fn("sin", arg));
        }
        if (name == "ln" || name == "log") {
            return sub(mul(arg, fn("ln", arg)), arg);
        }
    }
    if (const std::optional<Expr> scale = constant_derivative(arg, variable_name)) {
        if (name == "sin") {
            return div(mul(number(-1), fn("cos", arg)), *scale);
        }
        if (name == "cos") {
            return div(fn("sin", arg), *scale);
        }
        if (name == "exp") {
            return div(fn("exp", arg), *scale);
        }
    }
    if (name == "exp" && is_negative_x_squared(arg, variable_name)) {
        return mul(div(fn("sqrt", symbol("pi")), number(2)), fn("erf", symbol(variable_name)));
    }
    if (name == "exp" && is_x_squared(arg, variable_name)) {
        return mul(div(fn("sqrt", symbol("pi")), number(2)), fn("erfi", symbol(variable_name)));
    }
    return std::nullopt;
}

std::optional<Expr> integrate_power(const Expr& base,
                                    const Expr& exponent,
                                    const std::string& variable_name) {
    if (!is_variable(base, variable_name)) {
        return std::nullopt;
    }
    int n = 0;
    if (!integer_number(exponent, &n)) {
        return std::nullopt;
    }
    if (n == -1) {
        return fn("ln", symbol(variable_name));
    }
    return mul(rational(1, n + 1), pow(symbol(variable_name), number(n + 1)));
}

std::optional<Expr> integrate_rational_pattern(const Expr& expr,
                                               const std::string& variable_name) {
    if (expr.kind() != ExprKind::Pow || expr.children().size() != 2 ||
        !is_minus_one(expr.children()[1])) {
        return std::nullopt;
    }
    const Expr denominator = simplify(expr.children()[0]);
    if (const std::optional<Expr> derivative = constant_derivative(denominator, variable_name)) {
        return div(fn("ln", denominator), *derivative);
    }
    if (denominator.kind() == ExprKind::Add && denominator.children().size() == 2) {
        const Expr lhs = denominator.children()[0];
        const Expr rhs = denominator.children()[1];
        if ((is_x_squared(lhs, variable_name) && is_one(rhs)) ||
            (is_x_squared(rhs, variable_name) && is_one(lhs))) {
            return fn("atan", symbol(variable_name));
        }
    }
    return std::nullopt;
}

std::optional<Expr> integrate_log_derivative_pattern(const Expr& expr,
                                                     const std::string& variable_name) {
    std::vector<Expr> factors;
    append_mul_factors(expr, &factors);
    if (factors.size() < 2) {
        return std::nullopt;
    }

    for (std::size_t i = 0; i < factors.size(); ++i) {
        const Expr factor = factors[i];
        if (factor.kind() != ExprKind::Pow ||
            factor.children().size() != 2 ||
            !is_minus_one(factor.children()[1])) {
            continue;
        }

        const Expr denominator = simplify(factor.children()[0]);
        std::vector<Expr> rest;
        for (std::size_t j = 0; j < factors.size(); ++j) {
            if (j != i) {
                rest.push_back(factors[j]);
            }
        }
        Expr numerator = rest.empty() ? number(1) : rest[0];
        for (std::size_t j = 1; j < rest.size(); ++j) {
            numerator = mul(numerator, rest[j]);
        }
        numerator = simplify(numerator);

        const Expr derivative = simplify(differentiate(denominator, variable_name));
        if (structurally_equal(numerator, derivative)) {
            return fn("ln", denominator);
        }
        if (!contains_variable(numerator, variable_name) &&
            !contains_variable(derivative, variable_name) &&
            derivative.kind() == ExprKind::Number) {
            return mul(div(numerator, derivative), fn("ln", denominator));
        }
    }
    return std::nullopt;
}

std::optional<Expr> integrate_by_parts_pattern(const Expr& expr,
                                               const std::string& variable_name) {
    std::vector<Expr> factors;
    append_mul_factors(expr, &factors);
    if (factors.size() != 2) {
        return std::nullopt;
    }
    const Expr x = symbol(variable_name);
    for (int order = 0; order < 2; ++order) {
        const Expr lhs = factors[order == 0 ? 0 : 1];
        const Expr rhs = factors[order == 0 ? 1 : 0];
        if (!is_variable(lhs, variable_name) ||
            rhs.kind() != ExprKind::Function ||
            rhs.children().size() != 1 ||
            !is_variable(rhs.children()[0], variable_name)) {
            continue;
        }
        if (rhs.function_name() == "exp") {
            return mul(fn("exp", x), sub(x, number(1)));
        }
        if (rhs.function_name() == "sin") {
            return add(mul(number(-1), mul(x, fn("cos", x))), fn("sin", x));
        }
        if (rhs.function_name() == "cos") {
            return add(mul(x, fn("sin", x)), fn("cos", x));
        }
        if (rhs.function_name() == "ln") {
            return sub(mul(div(pow(x, number(2)), number(2)), fn("ln", x)),
                       div(pow(x, number(2)), number(4)));
        }
    }
    return std::nullopt;
}

std::optional<Expr> integrate_internal(const Expr& raw_expr, const std::string& variable_name) {
    const Expr expr = simplify(raw_expr);

    if (!contains_variable(expr, variable_name)) {
        return mul(expr, symbol(variable_name));
    }
    if (is_variable(expr, variable_name)) {
        return div(pow(symbol(variable_name), number(2)), number(2));
    }

    if (expr.kind() == ExprKind::Add && expr.children().size() == 2) {
        const std::optional<Expr> lhs = integrate_internal(expr.children()[0], variable_name);
        const std::optional<Expr> rhs = integrate_internal(expr.children()[1], variable_name);
        const Expr lhs_result = lhs ? *lhs : Expr::integral(expr.children()[0], symbol(variable_name));
        const Expr rhs_result = rhs ? *rhs : Expr::integral(expr.children()[1], symbol(variable_name));
        return add(lhs_result, rhs_result);
    }

    if (expr.kind() == ExprKind::Mul && expr.children().size() == 2) {
        if (const std::optional<Expr> log_derivative =
                integrate_log_derivative_pattern(expr, variable_name)) {
            return log_derivative;
        }
        if (const std::optional<Expr> polynomial_product =
                integrate_polynomial_times_known_function(expr, variable_name)) {
            return polynomial_product;
        }
        if (!contains_variable(expr.children()[0], variable_name)) {
            const std::optional<Expr> rhs = integrate_internal(expr.children()[1], variable_name);
            const Expr rhs_result = rhs ? *rhs : Expr::integral(expr.children()[1], symbol(variable_name));
            return mul(expr.children()[0], rhs_result);
        }
        if (!contains_variable(expr.children()[1], variable_name)) {
            const std::optional<Expr> lhs = integrate_internal(expr.children()[0], variable_name);
            const Expr lhs_result = lhs ? *lhs : Expr::integral(expr.children()[0], symbol(variable_name));
            return mul(expr.children()[1], lhs_result);
        }
        if (const std::optional<Expr> chain = integrate_chain_rule(expr, variable_name)) {
            return chain;
        }
        if (const std::optional<Expr> parts = integrate_by_parts_pattern(expr, variable_name)) {
            return parts;
        }
    }

    if (expr.kind() == ExprKind::Pow && expr.children().size() == 2) {
        if (const std::optional<Expr> rational_pattern =
                integrate_rational_pattern(expr, variable_name)) {
            return rational_pattern;
        }
        if (const std::optional<Expr> power_result =
                integrate_power(expr.children()[0], expr.children()[1], variable_name)) {
            return power_result;
        }
    }

    if (expr.kind() == ExprKind::Function) {
        if (const std::optional<Expr> function_result =
                integrate_function(expr, variable_name)) {
            return function_result;
        }
    }

    return std::nullopt;
}

}  // namespace

Expr integrate(const Expr& expr, const std::string& variable_name) {
    if (const std::optional<Expr> result = integrate_internal(expr, variable_name)) {
        return simplify(*result);
    }
    return Expr::integral(simplify(expr), symbol(variable_name));
}

}  // namespace cas

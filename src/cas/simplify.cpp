#include "simplify.h"

#include "printer.h"

#include <algorithm>
#include <map>

namespace cas {
namespace {

expression::Expr number(long long value) {
    return expression::Expr::number(numeric::Number(value));
}

bool is_number_text(const expression::Expr& expr, const std::string& text) {
    return expr.kind() == expression::ExprKind::Number &&
           expr.number_value().to_string() == text;
}

bool is_zero(const expression::Expr& expr) {
    return is_number_text(expr, "0");
}

bool is_one(const expression::Expr& expr) {
    return is_number_text(expr, "1");
}

bool is_two(const expression::Expr& expr) {
    return is_number_text(expr, "2");
}

bool structurally_equal(const expression::Expr& lhs, const expression::Expr& rhs) {
    return expression::print(lhs) == expression::print(rhs);
}

bool power_parts(const expression::Expr& expr,
                 expression::Expr* base,
                 expression::Expr* exponent) {
    if (expr.kind() == expression::ExprKind::Pow &&
        expr.children().size() == 2) {
        *base = expr.children()[0];
        *exponent = expr.children()[1];
        return true;
    }
    *base = expr;
    *exponent = number(1);
    return true;
}

bool is_square_of_function(const expression::Expr& expr,
                           const std::string& function_name,
                           expression::Expr* argument) {
    if (expr.kind() != expression::ExprKind::Pow ||
        expr.children().size() != 2 ||
        !is_two(expr.children()[1])) {
        return false;
    }
    const expression::Expr& base = expr.children()[0];
    if (base.kind() != expression::ExprKind::Function ||
        base.function_name() != function_name ||
        base.children().size() != 1) {
        return false;
    }
    *argument = base.children()[0];
    return true;
}

bool trig_pythagorean_identity(const expression::Expr& lhs,
                               const expression::Expr& rhs) {
    expression::Expr sin_arg;
    expression::Expr cos_arg;
    if (is_square_of_function(lhs, "sin", &sin_arg) &&
        is_square_of_function(rhs, "cos", &cos_arg)) {
        return structurally_equal(sin_arg, cos_arg);
    }
    if (is_square_of_function(lhs, "cos", &cos_arg) &&
        is_square_of_function(rhs, "sin", &sin_arg)) {
        return structurally_equal(sin_arg, cos_arg);
    }
    return false;
}

bool integer_exponent(const numeric::Number& value, int* exponent) {
    const std::string text = value.to_string();
    if (text.empty() || text.find('/') != std::string::npos ||
        text.find('.') != std::string::npos || text.find('i') != std::string::npos) {
        return false;
    }
    long long parsed = 0;
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
        if (parsed > 100000) {
            return false;
        }
    }
    *exponent = static_cast<int>(negative ? -parsed : parsed);
    return true;
}

expression::Expr fold_power(const expression::Expr& base,
                            const expression::Expr& exponent) {
    if (base.kind() != expression::ExprKind::Number ||
        exponent.kind() != expression::ExprKind::Number) {
        return expression::Expr::pow(base, exponent);
    }
    int exp = 0;
    if (!integer_exponent(exponent.number_value(), &exp)) {
        return expression::Expr::pow(base, exponent);
    }
    numeric::PrecisionContext context;
    numeric::Number result(1);
    const int count = exp < 0 ? -exp : exp;
    for (int i = 0; i < count; ++i) {
        result = numeric::multiply(result, base.number_value(), context);
    }
    if (exp < 0) {
        result = numeric::divide(numeric::Number(1), result, context);
    }
    return expression::Expr::number(result);
}

bool compare_expr(const expression::Expr& lhs, const expression::Expr& rhs) {
    return expression::print(lhs) < expression::print(rhs);
}

void append_add_terms(const expression::Expr& expr, std::vector<expression::Expr>* terms) {
    if (expr.kind() == expression::ExprKind::Add && expr.children().size() == 2) {
        append_add_terms(expr.children()[0], terms);
        append_add_terms(expr.children()[1], terms);
        return;
    }
    terms->push_back(expr);
}

void append_mul_factors(const expression::Expr& expr, std::vector<expression::Expr>* factors) {
    if (expr.kind() == expression::ExprKind::Mul && expr.children().size() == 2) {
        append_mul_factors(expr.children()[0], factors);
        append_mul_factors(expr.children()[1], factors);
        return;
    }
    factors->push_back(expr);
}

expression::Expr build_add(std::vector<expression::Expr> terms) {
    if (terms.empty()) {
        return number(0);
    }
    std::sort(terms.begin(), terms.end(), compare_expr);
    expression::Expr result = terms[0];
    for (std::size_t i = 1; i < terms.size(); ++i) {
        result = expression::Expr::add(result, terms[i]);
    }
    return result;
}

expression::Expr build_mul(std::vector<expression::Expr> factors) {
    if (factors.empty()) {
        return number(1);
    }
    std::sort(factors.begin(), factors.end(), compare_expr);
    expression::Expr result = factors[0];
    for (std::size_t i = 1; i < factors.size(); ++i) {
        result = expression::Expr::mul(result, factors[i]);
    }
    return result;
}

expression::Expr canonicalize_mul(const std::vector<expression::Expr>& raw_factors) {
    numeric::PrecisionContext context;
    numeric::Number coefficient(1);
    std::map<std::string, std::pair<expression::Expr, numeric::Number>> powers;

    for (const expression::Expr& factor : raw_factors) {
        std::vector<expression::Expr> nested;
        append_mul_factors(factor, &nested);
        for (const expression::Expr& nested_factor : nested) {
            if (nested_factor.kind() == expression::ExprKind::Number) {
                coefficient = numeric::multiply(coefficient, nested_factor.number_value(), context);
                continue;
            }
            expression::Expr base;
            expression::Expr exponent;
            power_parts(nested_factor, &base, &exponent);
            if (exponent.kind() != expression::ExprKind::Number) {
                const std::string key = expression::print(nested_factor);
                powers.emplace(key, std::make_pair(nested_factor, numeric::Number(1)));
                continue;
            }
            const std::string key = expression::print(base);
            auto it = powers.find(key);
            if (it == powers.end()) {
                powers.emplace(key, std::make_pair(base, exponent.number_value()));
            } else {
                it->second.second = numeric::add(it->second.second,
                                                 exponent.number_value(),
                                                 context);
            }
        }
    }

    if (coefficient.to_string() == "0") {
        return number(0);
    }

    std::vector<expression::Expr> factors;
    if (coefficient.to_string() != "1") {
        factors.push_back(expression::Expr::number(coefficient));
    }
    for (const auto& [_, entry] : powers) {
        const expression::Expr& base = entry.first;
        const numeric::Number& exponent = entry.second;
        if (exponent.to_string() == "0") {
            continue;
        }
        if (exponent.to_string() == "1") {
            factors.push_back(base);
        } else {
            factors.push_back(expression::Expr::pow(base, expression::Expr::number(exponent)));
        }
    }
    return build_mul(factors);
}

struct AddTerm {
    numeric::Number coefficient;
    expression::Expr base;
};

AddTerm split_add_term(const expression::Expr& term) {
    numeric::PrecisionContext context;
    std::vector<expression::Expr> factors;
    append_mul_factors(term, &factors);
    numeric::Number coefficient(1);
    std::vector<expression::Expr> symbolic_factors;
    for (const expression::Expr& factor : factors) {
        if (factor.kind() == expression::ExprKind::Number) {
            coefficient = numeric::multiply(coefficient, factor.number_value(), context);
        } else {
            symbolic_factors.push_back(factor);
        }
    }
    return {coefficient, build_mul(symbolic_factors)};
}

expression::Expr canonicalize_add(const std::vector<expression::Expr>& raw_terms) {
    numeric::PrecisionContext context;
    numeric::Number constant(0);
    std::map<std::string, AddTerm> terms;

    for (const expression::Expr& term : raw_terms) {
        std::vector<expression::Expr> nested;
        append_add_terms(term, &nested);
        for (const expression::Expr& nested_term : nested) {
            if (nested_term.kind() == expression::ExprKind::Number) {
                constant = numeric::add(constant, nested_term.number_value(), context);
                continue;
            }
            const AddTerm split = split_add_term(nested_term);
            if (split.base.kind() == expression::ExprKind::Number && is_one(split.base)) {
                constant = numeric::add(constant, split.coefficient, context);
                continue;
            }
            const std::string key = expression::print(split.base);
            auto it = terms.find(key);
            if (it == terms.end()) {
                terms.emplace(key, split);
            } else {
                it->second.coefficient = numeric::add(it->second.coefficient,
                                                      split.coefficient,
                                                      context);
            }
        }
    }

    std::vector<expression::Expr> rebuilt;
    for (const auto& [_, term] : terms) {
        if (term.coefficient.to_string() == "0") {
            continue;
        }
        if (term.coefficient.to_string() == "1") {
            rebuilt.push_back(term.base);
        } else {
            rebuilt.push_back(canonicalize_mul({expression::Expr::number(term.coefficient),
                                                term.base}));
        }
    }
    if (constant.to_string() != "0") {
        rebuilt.push_back(expression::Expr::number(constant));
    }

    if (rebuilt.size() == 2 && trig_pythagorean_identity(rebuilt[0], rebuilt[1])) {
        return number(1);
    }
    return build_add(rebuilt);
}

}  // namespace

expression::Expr simplify(const expression::Expr& expr) {
    switch (expr.kind()) {
        case expression::ExprKind::Number:
        case expression::ExprKind::Symbol:
            return expr;
        case expression::ExprKind::Add: {
            const expression::Expr lhs = simplify(expr.children()[0]);
            const expression::Expr rhs = simplify(expr.children()[1]);
            return canonicalize_add({lhs, rhs});
        }
        case expression::ExprKind::Mul: {
            const expression::Expr lhs = simplify(expr.children()[0]);
            const expression::Expr rhs = simplify(expr.children()[1]);
            return canonicalize_mul({lhs, rhs});
        }
        case expression::ExprKind::Pow: {
            const expression::Expr base = simplify(expr.children()[0]);
            const expression::Expr exponent = simplify(expr.children()[1]);
            if (is_zero(exponent)) {
                return number(1);
            }
            if (is_one(exponent)) {
                return base;
            }
            if (is_one(base)) {
                return number(1);
            }
            if (is_zero(base)) {
                return number(0);
            }
            if (base.kind() == expression::ExprKind::Pow &&
                base.children().size() == 2 &&
                base.children()[1].kind() == expression::ExprKind::Number &&
                exponent.kind() == expression::ExprKind::Number) {
                numeric::PrecisionContext context;
                return simplify(expression::Expr::pow(base.children()[0],
                    expression::Expr::number(numeric::multiply(base.children()[1].number_value(),
                                                               exponent.number_value(),
                                                               context))));
            }
            return fold_power(base, exponent);
        }
        case expression::ExprKind::Function: {
            std::vector<expression::Expr> args;
            args.reserve(expr.children().size());
            for (const expression::Expr& child : expr.children()) {
                args.push_back(simplify(child));
            }
            if (expr.function_name() == "exp" && args.size() == 1 &&
                args[0].kind() == expression::ExprKind::Function &&
                args[0].function_name() == "ln" &&
                args[0].children().size() == 1) {
                return args[0].children()[0];
            }
            if (expr.function_name() == "ln" && args.size() == 1 &&
                args[0].kind() == expression::ExprKind::Function &&
                args[0].function_name() == "exp" &&
                args[0].children().size() == 1) {
                return args[0].children()[0];
            }
            return expression::Expr::function(expr.function_name(), args);
        }
        case expression::ExprKind::Integral:
            return expression::Expr::integral(simplify(expr.children()[0]), expr.children()[1]);
        case expression::ExprKind::List: {
            std::vector<expression::Expr> items;
            items.reserve(expr.children().size());
            for (const expression::Expr& child : expr.children()) {
                items.push_back(simplify(child));
            }
            return expression::Expr::list(items);
        }
    }
    return expr;
}

}  // namespace cas

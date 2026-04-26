#include "value.h"

#include "printer.h"

#include <stdexcept>

namespace runtime {

Value::Value() : value_(numeric::Number()) {}
Value::Value(const numeric::Number& number) : value_(number) {}
Value::Value(const expression::Expr& expr) : value_(expr) {}
Value::Value(const std::string& text) : value_(text) {}

bool Value::is_number() const {
    return std::holds_alternative<numeric::Number>(value_);
}

bool Value::is_expr() const {
    return std::holds_alternative<expression::Expr>(value_);
}

bool Value::is_string() const {
    return std::holds_alternative<std::string>(value_);
}

const numeric::Number& Value::as_number() const {
    if (!is_number()) {
        throw std::runtime_error("value is not a number");
    }
    return std::get<numeric::Number>(value_);
}

const expression::Expr& Value::as_expr() const {
    if (!is_expr()) {
        throw std::runtime_error("value is not an expression");
    }
    return std::get<expression::Expr>(value_);
}

const std::string& Value::as_string() const {
    if (!is_string()) {
        throw std::runtime_error("value is not a string");
    }
    return std::get<std::string>(value_);
}

std::string Value::to_string() const {
    if (is_number()) {
        return as_number().to_string();
    }
    if (is_expr()) {
        return expression::print(as_expr());
    }
    return as_string();
}

}  // namespace runtime

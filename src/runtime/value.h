#ifndef RUNTIME_VALUE_H
#define RUNTIME_VALUE_H

#include "expr.h"
#include "number.h"

#include <string>
#include <variant>

namespace runtime {

class Value {
public:
    Value();
    Value(const numeric::Number& number);
    Value(const expression::Expr& expr);
    Value(const std::string& text);

    bool is_number() const;
    bool is_expr() const;
    bool is_string() const;

    const numeric::Number& as_number() const;
    const expression::Expr& as_expr() const;
    const std::string& as_string() const;
    std::string to_string() const;

private:
    std::variant<numeric::Number, expression::Expr, std::string> value_;
};

}  // namespace runtime

#endif

#ifndef EXPRESSION_EXPR_H
#define EXPRESSION_EXPR_H

#include "number.h"

#include <memory>
#include <string>
#include <vector>

namespace expression {

enum class ExprKind {
    Number,
    Symbol,
    Add,
    Mul,
    Pow,
    Function,
    Integral,
    List,
};

class Expr {
public:
    struct Node;

    Expr();

    static Expr number(const numeric::Number& value);
    static Expr symbol(const std::string& name);
    static Expr add(const Expr& lhs, const Expr& rhs);
    static Expr mul(const Expr& lhs, const Expr& rhs);
    static Expr pow(const Expr& base, const Expr& exponent);
    static Expr function(const std::string& name, const std::vector<Expr>& args);
    static Expr integral(const Expr& integrand, const Expr& variable);
    static Expr list(const std::vector<Expr>& items);

    ExprKind kind() const;
    const numeric::Number& number_value() const;
    const std::string& symbol_name() const;
    const std::string& function_name() const;
    const std::vector<Expr>& children() const;

private:
    explicit Expr(std::shared_ptr<Node> node);

    std::shared_ptr<Node> node_;
};

}  // namespace expression

#endif

#include "expr.h"

#include <stdexcept>
#include <utility>

namespace expression {

struct Expr::Node {
    ExprKind kind = ExprKind::Number;
    numeric::Number number;
    std::string text;
    std::vector<Expr> children;
};

Expr::Expr() : node_(std::make_shared<Node>()) {}

Expr::Expr(std::shared_ptr<Node> node) : node_(std::move(node)) {}

Expr Expr::number(const numeric::Number& value) {
    auto node = std::make_shared<Node>();
    node->kind = ExprKind::Number;
    node->number = value;
    return Expr(node);
}

Expr Expr::symbol(const std::string& name) {
    auto node = std::make_shared<Node>();
    node->kind = ExprKind::Symbol;
    node->text = name;
    return Expr(node);
}

Expr Expr::add(const Expr& lhs, const Expr& rhs) {
    auto node = std::make_shared<Node>();
    node->kind = ExprKind::Add;
    node->children = {lhs, rhs};
    return Expr(node);
}

Expr Expr::mul(const Expr& lhs, const Expr& rhs) {
    auto node = std::make_shared<Node>();
    node->kind = ExprKind::Mul;
    node->children = {lhs, rhs};
    return Expr(node);
}

Expr Expr::pow(const Expr& base, const Expr& exponent) {
    auto node = std::make_shared<Node>();
    node->kind = ExprKind::Pow;
    node->children = {base, exponent};
    return Expr(node);
}

Expr Expr::function(const std::string& name, const std::vector<Expr>& args) {
    auto node = std::make_shared<Node>();
    node->kind = ExprKind::Function;
    node->text = name;
    node->children = args;
    return Expr(node);
}

Expr Expr::integral(const Expr& integrand, const Expr& variable) {
    auto node = std::make_shared<Node>();
    node->kind = ExprKind::Integral;
    node->children = {integrand, variable};
    return Expr(node);
}

Expr Expr::list(const std::vector<Expr>& items) {
    auto node = std::make_shared<Node>();
    node->kind = ExprKind::List;
    node->children = items;
    return Expr(node);
}

ExprKind Expr::kind() const {
    return node_->kind;
}

const numeric::Number& Expr::number_value() const {
    if (kind() != ExprKind::Number) {
        throw std::runtime_error("expression is not a number");
    }
    return node_->number;
}

const std::string& Expr::symbol_name() const {
    if (kind() != ExprKind::Symbol) {
        throw std::runtime_error("expression is not a symbol");
    }
    return node_->text;
}

const std::string& Expr::function_name() const {
    if (kind() != ExprKind::Function) {
        throw std::runtime_error("expression is not a function");
    }
    return node_->text;
}

const std::vector<Expr>& Expr::children() const {
    return node_->children;
}

}  // namespace expression

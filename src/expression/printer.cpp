#include "printer.h"

#include <sstream>

namespace expression {
namespace {

int precedence(ExprKind kind) {
    switch (kind) {
        case ExprKind::Add:
            return 1;
        case ExprKind::Mul:
            return 2;
        case ExprKind::Pow:
            return 3;
        default:
            return 4;
    }
}

std::string print_with_parent(const Expr& expr, int parent_precedence) {
    std::ostringstream out;
    switch (expr.kind()) {
        case ExprKind::Number:
            out << expr.number_value().to_string();
            break;
        case ExprKind::Symbol:
            out << expr.symbol_name();
            break;
        case ExprKind::Add:
            out << print_with_parent(expr.children()[0], precedence(expr.kind()))
                << " + "
                << print_with_parent(expr.children()[1], precedence(expr.kind()));
            break;
        case ExprKind::Mul:
            out << print_with_parent(expr.children()[0], precedence(expr.kind()))
                << " * "
                << print_with_parent(expr.children()[1], precedence(expr.kind()));
            break;
        case ExprKind::Pow:
            out << print_with_parent(expr.children()[0], precedence(expr.kind()))
                << " ^ "
                << print_with_parent(expr.children()[1], precedence(expr.kind()));
            break;
        case ExprKind::Function:
            out << expr.function_name() << "(";
            for (std::size_t i = 0; i < expr.children().size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << print(expr.children()[i]);
            }
            out << ")";
            break;
        case ExprKind::Integral:
            out << "integral(" << print(expr.children()[0]) << ", "
                << print(expr.children()[1]) << ")";
            break;
        case ExprKind::List:
            out << "[";
            for (std::size_t i = 0; i < expr.children().size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << print(expr.children()[i]);
            }
            out << "]";
            break;
    }
    const std::string text = out.str();
    if (precedence(expr.kind()) < parent_precedence) {
        return "(" + text + ")";
    }
    return text;
}

}  // namespace

std::string print(const Expr& expr) {
    return print_with_parent(expr, 0);
}

}  // namespace expression

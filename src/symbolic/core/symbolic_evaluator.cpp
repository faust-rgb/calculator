/**
 * @file symbolic_evaluator.cpp
 * @brief 显式符号数值转换器与代数解耦实现
 */

#include "symbolic/core/symbolic_evaluator.h"
#include "math/mymath.h"
#include "math/helpers/integer_helpers.h"
#include "matrix.h"
#include <cmath>
#include <stdexcept>

namespace symbolic {

StoredValue SymbolicEvaluator::evalf(const SymbolicExpression& expr, core::ExecutionContext& ctx) {
    if (!expr.has_node()) {
        return StoredValue(Scalar(0.0L));
    }

    NodeType type = expr.node_type();

    switch (type) {
    case NodeType::kNumber: {
        return StoredValue(expr.node_numeric_value());
    }
    case NodeType::kPi: {
        return StoredValue(mymath::pi());
    }
    case NodeType::kE: {
        return StoredValue(mymath::e());
    }
    case NodeType::kVariable: {
        const std::string& name = expr.node_text();
        const auto* slot = ctx.scope().find(name);
        if (slot) {
            return slot->value;
        }
        if (name == "pi" || name == "PI") return StoredValue(mymath::pi());
        if (name == "e" || name == "E") return StoredValue(mymath::e());
        throw std::runtime_error("Undefined symbolic variable in evalf: " + name);
    }
    case NodeType::kAdd: {
        StoredValue l = evalf(expr.left_child(), ctx);
        StoredValue r = evalf(expr.right_child(), ctx);
        if (l.is_scalar() && r.is_scalar()) {
            return StoredValue(l.get_decimal() + r.get_decimal());
        }
        if (l.is_matrix() && r.is_matrix()) {
            return StoredValue(matrix::add(l.as_matrix(), r.as_matrix()));
        }
        throw std::runtime_error("Type error in symbolic add evalf");
    }
    case NodeType::kSubtract: {
        StoredValue l = evalf(expr.left_child(), ctx);
        StoredValue r = evalf(expr.right_child(), ctx);
        if (l.is_scalar() && r.is_scalar()) {
            return StoredValue(l.get_decimal() - r.get_decimal());
        }
        if (l.is_matrix() && r.is_matrix()) {
            return StoredValue(matrix::subtract(l.as_matrix(), r.as_matrix()));
        }
        throw std::runtime_error("Type error in symbolic sub evalf");
    }
    case NodeType::kMultiply: {
        StoredValue l = evalf(expr.left_child(), ctx);
        StoredValue r = evalf(expr.right_child(), ctx);
        if (l.is_scalar() && r.is_scalar()) {
            return StoredValue(l.get_decimal() * r.get_decimal());
        }
        if (l.is_scalar() && r.is_matrix()) {
            return StoredValue(matrix::multiply(r.as_matrix(), l.get_decimal()));
        }
        if (l.is_matrix() && r.is_scalar()) {
            return StoredValue(matrix::multiply(l.as_matrix(), r.get_decimal()));
        }
        if (l.is_matrix() && r.is_matrix()) {
            return StoredValue(matrix::multiply(l.as_matrix(), r.as_matrix()));
        }
        throw std::runtime_error("Type error in symbolic mul evalf");
    }
    case NodeType::kDivide: {
        StoredValue l = evalf(expr.left_child(), ctx);
        StoredValue r = evalf(expr.right_child(), ctx);
        if (l.is_scalar() && r.is_scalar()) {
            Scalar denom = r.get_decimal();
            if (denom == Scalar(0.0L)) throw std::runtime_error("Division by zero in evalf");
            return StoredValue(l.get_decimal() / denom);
        }
        throw std::runtime_error("Type error in symbolic div evalf");
    }
    case NodeType::kPower: {
        StoredValue l = evalf(expr.left_child(), ctx);
        StoredValue r = evalf(expr.right_child(), ctx);
        if (l.is_scalar() && r.is_scalar()) {
            return StoredValue(mymath::pow(l.get_decimal(), r.get_decimal()));
        }
        throw std::runtime_error("Type error in symbolic power evalf");
    }
    case NodeType::kNegate: {
        StoredValue child = evalf(expr.left_child(), ctx);
        if (child.is_scalar()) {
            return StoredValue(-child.get_decimal());
        }
        if (child.is_matrix()) {
            return StoredValue(matrix::multiply(child.as_matrix(), Scalar(-1.0L)));
        }
        throw std::runtime_error("Type error in symbolic negate evalf");
    }
    case NodeType::kFunction: {
        const std::string& func_name = expr.node_text();
        std::vector<StoredValue> args;
        std::size_t count = expr.child_count();
        for (std::size_t i = 0; i < count; ++i) {
            args.push_back(evalf(expr.child_at(i), ctx));
        }

        const auto* native_func = ctx.functions().find_function(func_name);
        if (native_func) {
            return (*native_func)(args, ctx);
        }

        // 内置标量函数 fallback
        if (args.size() == 1 && args[0].is_scalar()) {
            Scalar val = args[0].get_decimal();
            if (func_name == "sin") return StoredValue(mymath::sin(val));
            if (func_name == "cos") return StoredValue(mymath::cos(val));
            if (func_name == "tan") return StoredValue(mymath::tan(val));
            if (func_name == "exp") return StoredValue(mymath::exp(val));
            if (func_name == "ln" || func_name == "log") return StoredValue(mymath::ln(val));
            if (func_name == "sqrt") return StoredValue(mymath::sqrt(val));
        }
        throw std::runtime_error("Unknown function in symbolic evalf: " + func_name);
    }
    default:
        throw std::runtime_error("Unsupported node type in symbolic evalf");
    }
}

bool SymbolicEvaluator::is_exact_algebraic(const SymbolicExpression& expr) {
    if (!expr.has_node()) return true;
    NodeType type = expr.node_type();
    if (type == NodeType::kNumber) {
        Scalar val = expr.node_numeric_value();
        return is_integer_double(val);
    }
    if (type == NodeType::kRootOf) return true;
    if (type == NodeType::kPi || type == NodeType::kE ||
        type == NodeType::kInfinity || type == NodeType::kVariable ||
        type == NodeType::kFunction || type == NodeType::kDifferentialOp) {
        return false;
    }
    for (std::size_t i = 0; i < expr.child_count(); ++i) {
        if (!is_exact_algebraic(expr.child_at(i))) return false;
    }
    return true;
}

} // namespace symbolic

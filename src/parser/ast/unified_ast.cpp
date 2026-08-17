/**
 * @file unified_ast.cpp
 * @brief 统一 AST 节点与求值实现
 */

#include "parser/ast/unified_ast.h"
#include "parser/ast/expression_ast.h"
#include "matrix.h"
#include "math/mymath.h"
#include <cmath>

namespace core {

StoredValue VariableNode::evaluate(ExecutionContext& ctx) const {
    const auto* slot = ctx.scope().find(name_);
    if (slot) {
        return slot->value;
    }
    if (ctx.external_variable_lookup()) {
        StoredValue value;
        if (ctx.external_variable_lookup()(name_, &value)) return value;
    }
    // 特殊常量预设
    if (name_ == "pi" || name_ == "PI") {
        return StoredValue(Scalar(mymath::kPi));
    }
    if (name_ == "e" || name_ == "E") {
        return StoredValue(Scalar(mymath::kE));
    }
    if (name_ == "i") {
        return StoredValue(mymath::complex<Scalar>(0.0L, 1.0L));
    }
    throw std::runtime_error("Undefined variable: " + name_);
}

StoredValue UnaryOpNode::evaluate(ExecutionContext& ctx) const {
    StoredValue child_val = child_->evaluate(ctx);
    if (op_ == '-') {
        if (child_val.is_scalar()) {
            return StoredValue(-child_val.get_decimal());
        } else if (child_val.is_complex()) {
            auto c = child_val.get_complex();
            return StoredValue(mymath::complex<Scalar>(-c.real(), -c.imag()));
        } else if (child_val.is_matrix()) {
            auto m = child_val.as_matrix();
            return StoredValue(matrix::multiply(m, Scalar(-1.0L)));
        }
    } else if (op_ == '+') {
        return child_val;
    } else if (op_ == '!') {
        if (child_val.is_scalar()) {
            bool b = static_cast<long double>(child_val.get_decimal()) != 0.0L;
            return StoredValue(Scalar(!b ? 1.0L : 0.0L));
        }
    }
    return child_val;
}

StoredValue BinaryOpNode::evaluate(ExecutionContext& ctx) const {
    StoredValue l = left_->evaluate(ctx);
    StoredValue r = right_->evaluate(ctx);

    if (op_ == "&&" || op_ == "||") {
        const bool lv = l.get_decimal() != Scalar(0);
        if (op_ == "&&" && !lv) return StoredValue(Scalar(0));
        if (op_ == "||" && lv) return StoredValue(Scalar(1));
        return StoredValue(Scalar(r.get_decimal() != Scalar(0) ? 1 : 0));
    }
    if (op_ == "==" || op_ == "!=" || op_ == "<" || op_ == ">" || op_ == "<=" || op_ == ">=") {
        const Scalar a = l.get_decimal(), b = r.get_decimal();
        bool result = false;
        if (op_ == "==") result = a == b;
        else if (op_ == "!=") result = a != b;
        else if (op_ == "<") result = a < b;
        else if (op_ == ">") result = a > b;
        else if (op_ == "<=") result = a <= b;
        else result = a >= b;
        return StoredValue(Scalar(result ? 1 : 0));
    }

    if (op_ == "+") {
        if (l.is_scalar() && r.is_scalar()) {
            return StoredValue(l.get_decimal() + r.get_decimal());
        }
        if (l.is_matrix() && r.is_matrix()) {
            return StoredValue(matrix::add(l.as_matrix(), r.as_matrix()));
        }
        if (l.is_complex() || r.is_complex()) {
            return StoredValue(l.get_complex() + r.get_complex());
        }
        if (l.is_string() || r.is_string()) {
            return StoredValue(l.as_string() + r.as_string());
        }
    } else if (op_ == "-") {
        if (l.is_scalar() && r.is_scalar()) {
            return StoredValue(l.get_decimal() - r.get_decimal());
        }
        if (l.is_matrix() && r.is_matrix()) {
            return StoredValue(matrix::subtract(l.as_matrix(), r.as_matrix()));
        }
        if (l.is_complex() || r.is_complex()) {
            return StoredValue(l.get_complex() - r.get_complex());
        }
    } else if (op_ == "*") {
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
        if (l.is_complex() || r.is_complex()) {
            return StoredValue(l.get_complex() * r.get_complex());
        }
    } else if (op_ == "/") {
        if (l.is_scalar() && r.is_scalar()) {
            Scalar denom = r.get_decimal();
            if (denom == Scalar(0.0L)) {
                throw std::runtime_error("Division by zero");
            }
            return StoredValue(l.get_decimal() / denom);
        }
        if (l.is_complex() || r.is_complex()) {
            return StoredValue(l.get_complex() / r.get_complex());
        }
    } else if (op_ == "^") {
        if (l.is_scalar() && r.is_scalar()) {
            return StoredValue(mymath::pow(l.get_decimal(), r.get_decimal()));
        }
    }
    throw std::runtime_error("Unsupported binary operator '" + op_ + "' for given types");
}

StoredValue ConditionalNode::evaluate(ExecutionContext& ctx) const {
    return condition_->evaluate(ctx).get_decimal() != Scalar(0)
        ? then_->evaluate(ctx) : else_->evaluate(ctx);
}

StoredValue FunctionCallNode::evaluate(ExecutionContext& ctx) const {
    std::vector<StoredValue> arg_values;
    arg_values.reserve(args_.size());
    for (const auto& arg : args_) {
        arg_values.push_back(arg->evaluate(ctx));
    }

    // 优先查 FunctionRegistry
    const auto* func = ctx.functions().find_function(name_);
    if (func) {
        return (*func)(arg_values, ctx);
    }

    throw std::runtime_error("Unknown function: " + name_);
}

StoredValue MatrixLiteralNode::evaluate(ExecutionContext& ctx) const {
    if (rows_.empty()) {
        return StoredValue(matrix::Matrix(0, 0));
    }
    std::size_t num_rows = rows_.size();
    std::size_t num_cols = rows_[0].size();
    matrix::Matrix mat(num_rows, num_cols);

    for (std::size_t r = 0; r < num_rows; ++r) {
        if (rows_[r].size() != num_cols) {
            throw std::runtime_error("Inconsistent matrix row lengths");
        }
        for (std::size_t c = 0; c < num_cols; ++c) {
            StoredValue elem = rows_[r][c]->evaluate(ctx);
            mat.at(r, c) = elem.get_decimal();
        }
    }
    return StoredValue(mat);
}

StoredValue IndexAccessNode::evaluate(ExecutionContext& ctx) const {
    if (parts_.empty()) throw std::runtime_error("empty index expression");
    StoredValue base = parts_[0]->evaluate(ctx);
    if (!base.is_matrix()) throw std::runtime_error("index access requires a matrix");
    if (parts_.size() != 2 && parts_.size() != 3)
        throw std::runtime_error("matrix index expects row and column");
    auto index = [&](std::size_t n) -> std::size_t {
        StoredValue v = parts_[n]->evaluate(ctx);
        if (!v.is_scalar()) throw std::runtime_error("matrix index must be scalar");
        auto i = static_cast<long long>(v.get_decimal());
        if (i < 0) throw std::runtime_error("matrix index must be non-negative");
        return static_cast<std::size_t>(i);
    };
    const auto& m = base.as_matrix();
    std::size_t row = index(1);
    if (parts_.size() == 2) {
        if (row >= m.rows || m.cols != 1) throw std::runtime_error("matrix index out of range");
        return StoredValue(m.at(row, 0));
    }
    std::size_t col = index(2);
    if (row >= m.rows || col >= m.cols) throw std::runtime_error("matrix index out of range");
    return StoredValue(m.at(row, col));
}

StoredValue ListLiteralNode::evaluate(ExecutionContext& ctx) const {
    std::vector<StoredValue> values;
    values.reserve(elements_.size());
    for (const auto& elem : elements_) {
        values.push_back(elem->evaluate(ctx));
    }
    StoredValue res;
    res.data = std::make_shared<std::vector<StoredValue>>(std::move(values));
    return res;
}

StoredValue DictLiteralNode::evaluate(ExecutionContext& ctx) const {
    std::map<std::string, StoredValue> values;
    for (const auto& [k_node, v_node] : pairs_) {
        StoredValue k_val = k_node->evaluate(ctx);
        StoredValue v_val = v_node->evaluate(ctx);
        std::string key_str = k_val.is_string() ? k_val.as_string() : (k_val.exact ? k_val.rational.to_string() : std::to_string(static_cast<double>(k_val.as_scalar())));
        values[key_str] = std::move(v_val);
    }
    StoredValue res;
    res.data = std::make_shared<std::map<std::string, StoredValue>>(std::move(values));
    return res;
}

StoredValue SliceNode::evaluate(ExecutionContext& ctx) const {
    StoredValue start_v = start_ ? start_->evaluate(ctx) : StoredValue(Scalar(0.0L));
    StoredValue stop_v = stop_ ? stop_->evaluate(ctx) : StoredValue(Scalar(0.0L));
    StoredValue step_v = step_ ? step_->evaluate(ctx) : StoredValue(Scalar(1.0L));

    std::vector<StoredValue> slice_tuple;
    slice_tuple.push_back(start_v);
    slice_tuple.push_back(stop_v);
    slice_tuple.push_back(step_v);

    StoredValue res;
    res.is_list = true;
    res.list_value = std::make_shared<std::vector<StoredValue>>(std::move(slice_tuple));
    return res;
}

std::unique_ptr<ASTNode> build_unified_ast(const ExpressionAST* ast) {
    if (!ast) return nullptr;
    switch (ast->kind) {
    case ExprKind::kNumber:
        return std::make_unique<ScalarNode>(StoredValue(ast->number_value));
    case ExprKind::kString:
        return std::make_unique<StringNode>(ast->string_value);
    case ExprKind::kVariable:
        return std::make_unique<VariableNode>(ast->identifier);
    case ExprKind::kImaginary:
        return std::make_unique<ScalarNode>(StoredValue(mymath::complex<Scalar>(0.0L, 1.0L)));
    case ExprKind::kMatrixLiteral: {
        std::vector<std::vector<std::unique_ptr<ASTNode>>> rows;
        for (const auto& row : ast->matrix_rows) {
            std::vector<std::unique_ptr<ASTNode>> out;
            for (const auto& cell : row) out.push_back(build_unified_ast(cell.get()));
            rows.push_back(std::move(out));
        }
        return std::make_unique<MatrixLiteralNode>(std::move(rows));
    }
    case ExprKind::kIndexAccess: {
        std::vector<std::unique_ptr<ASTNode>> parts;
        for (const auto& child : ast->children) parts.push_back(build_unified_ast(child.get()));
        return std::make_unique<IndexAccessNode>(std::move(parts));
    }
    case ExprKind::kPostfixOp:
        if (!ast->children.empty()) return std::make_unique<UnaryOpNode>(ast->op_char, build_unified_ast(ast->children[0].get()));
        break;
    case ExprKind::kUnaryOp:
        if (!ast->children.empty()) {
            return std::make_unique<UnaryOpNode>(ast->op_char, build_unified_ast(ast->children[0].get()));
        }
        break;
    case ExprKind::kBinaryOp:
        if (ast->children.size() >= 2) {
            return std::make_unique<BinaryOpNode>(
                std::string(1, ast->op_char),
                build_unified_ast(ast->children[0].get()),
                build_unified_ast(ast->children[1].get())
            );
        }
        break;
    case ExprKind::kComparison:
        if (ast->children.size() >= 2) {
            return std::make_unique<BinaryOpNode>(
                ast->comparison_op,
                build_unified_ast(ast->children[0].get()),
                build_unified_ast(ast->children[1].get())
            );
        }
        break;
    case ExprKind::kLogicalOp:
        if (ast->children.size() >= 2) {
            return std::make_unique<BinaryOpNode>(
                ast->comparison_op,
                build_unified_ast(ast->children[0].get()),
                build_unified_ast(ast->children[1].get()));
        }
        break;
    case ExprKind::kConditional:
        if (ast->children.size() == 3) {
            return std::make_unique<ConditionalNode>(
                build_unified_ast(ast->children[0].get()),
                build_unified_ast(ast->children[1].get()),
                build_unified_ast(ast->children[2].get()));
        }
        break;
        break;
    case ExprKind::kFunctionCall: {
        std::vector<std::unique_ptr<ASTNode>> args;
        for (const auto& child : ast->children) {
            args.push_back(build_unified_ast(child.get()));
        }
        return std::make_unique<FunctionCallNode>(ast->identifier, std::move(args));
    }
    default:
        break;
    }
    return std::make_unique<ScalarNode>(StoredValue(Scalar(0.0L)));
}

StoredValue evaluate_unified_ast(const ExpressionAST* ast, ExecutionContext& ctx) {
    auto node = build_unified_ast(ast);
    if (!node) {
        return StoredValue(Scalar(0.0L));
    }
    return ASTEvaluator::evaluate(*node, ctx);
}

} // namespace core

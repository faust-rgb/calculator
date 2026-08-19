/**
 * @file unified_ast.cpp
 * @brief 统一 AST 节点与求值实现
 */

#include "parser/ast/unified_ast.h"
#include "parser/ast/expression_ast.h"
#include "matrix.h"
#include "math/mymath.h"
#include <algorithm>
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
        if (!l.is_scalar() || !r.is_scalar()) {
            if (op_ == "==" || op_ == "!=") {
                if (l.is_complex() && r.is_complex()) {
                    const bool equal = l.as_complex() == r.as_complex();
                    return StoredValue(Scalar((op_ == "==") == equal ? 1 : 0));
                }
            }
            throw std::runtime_error("comparisons require scalar values");
        }
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
        if (l.is_matrix() && r.is_scalar()) {
            return StoredValue(matrix::add(l.as_matrix(), r.get_decimal()));
        }
        if (l.is_scalar() && r.is_matrix()) {
            return StoredValue(matrix::add(r.as_matrix(), l.get_decimal()));
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
        if (l.is_matrix() && r.is_scalar()) {
            return StoredValue(matrix::subtract(l.as_matrix(), r.get_decimal()));
        }
        if (l.is_scalar() && r.is_matrix()) {
            return StoredValue(matrix::add(matrix::multiply(r.as_matrix(), Scalar(-1.0L)), l.get_decimal()));
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
        if (l.is_matrix() && r.is_scalar()) {
            if (r.get_decimal() == Scalar(0.0L)) throw std::runtime_error("Division by zero");
            return StoredValue(matrix::divide(l.as_matrix(), r.get_decimal()));
        }
        if (l.is_complex() || r.is_complex()) {
            return StoredValue(l.get_complex() / r.get_complex());
        }
    } else if (op_ == "^") {
        if (l.is_scalar() && r.is_scalar()) {
            return StoredValue(mymath::pow(l.get_decimal(), r.get_decimal()));
        }
        if (l.is_matrix() && r.is_scalar()) {
            const Scalar exponent = r.get_decimal();
            if (!mymath::isfinite(exponent) || mymath::floor(exponent) != exponent) {
                throw std::runtime_error("matrix exponent must be an integer");
            }
            return StoredValue(matrix::power(l.as_matrix(), static_cast<long long>(exponent)));
        }
    }
    throw std::runtime_error("Unsupported binary operator '" + op_ + "' for given types");
}

StoredValue ConditionalNode::evaluate(ExecutionContext& ctx) const {
    return condition_->evaluate(ctx).get_decimal() != Scalar(0)
        ? then_->evaluate(ctx) : else_->evaluate(ctx);
}

StoredValue PostfixOpNode::evaluate(ExecutionContext& ctx) const {
    StoredValue value = child_->evaluate(ctx);
    if (op_ == "'") {
        if (!value.is_matrix()) throw std::runtime_error("transpose requires a matrix");
        return StoredValue(matrix::transpose(value.as_matrix()));
    }
    if (!value.is_scalar()) throw std::runtime_error("postfix operator requires a scalar");
    const Scalar n = value.get_decimal();
    if (op_ == "%") return StoredValue(n / Scalar(100));
    if (op_ != "!" && op_ != "!!") throw std::runtime_error("unknown postfix operator: " + op_);
    if (n < 0 || mymath::floor(n) != n) throw std::runtime_error("factorial requires a non-negative integer");
    const long long step = op_ == "!!" ? 2 : 1;
    const long long start = op_ == "!!" && static_cast<long long>(n) % 2 == 0 ? 2 : 1;
    Scalar result = 1;
    for (long long i = static_cast<long long>(n); i >= start; i -= step) result *= Scalar(i);
    return StoredValue(result);
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
    std::size_t num_cols = 0;
    for (const auto& row : rows_) num_cols = std::max(num_cols, row.size());
    matrix::Matrix mat(num_rows, num_cols);

    for (std::size_t r = 0; r < num_rows; ++r) {
        for (std::size_t c = 0; c < num_cols; ++c) {
            if (c >= rows_[r].size()) continue;
            StoredValue elem = rows_[r][c]->evaluate(ctx);
            if (!elem.is_scalar()) {
                throw std::runtime_error("matrix literal entries must be scalar values");
            }
            mat.at(r, c) = elem.get_decimal();
        }
    }
    return StoredValue(mat);
}

StoredValue IndexAccessNode::evaluate(ExecutionContext& ctx) const {
    if (parts_.empty()) throw std::runtime_error("empty index expression");
    StoredValue base = parts_[0]->evaluate(ctx);
    if (parts_.size() == 2 && parts_[1]->type() == ASTNodeType::kSlice) {
        const StoredValue bounds = parts_[1]->evaluate(ctx);
        if (!base.is_list() || !bounds.is_list() || !bounds.list_value || bounds.list_value->size() != 3) {
            throw std::runtime_error("slice access requires a list");
        }
        const auto& list = *base.list_value;
        const auto& b = *bounds.list_value;
        auto bound = [&](const StoredValue& value, long long fallback) {
            if (value.is_nil()) return fallback;
            if (!value.is_scalar() || mymath::floor(value.get_decimal()) != value.get_decimal())
                throw std::runtime_error("slice bound must be an integer");
            return static_cast<long long>(value.get_decimal());
        };
        long long start = bound(b[0], 0);
        long long stop = bound(b[1], static_cast<long long>(list.size()));
        long long step = bound(b[2], 1);
        if (step == 0) throw std::runtime_error("slice step cannot be zero");
        auto normalize = [&](long long index) {
            if (index < 0) index += static_cast<long long>(list.size());
            return std::max(0LL, std::min(index, static_cast<long long>(list.size())));
        };
        start = normalize(start); stop = normalize(stop);
        std::vector<StoredValue> result;
        if (step > 0) for (long long i = start; i < stop; i += step) result.push_back(list[static_cast<std::size_t>(i)]);
        else for (long long i = start; i > stop; i += step) result.push_back(list[static_cast<std::size_t>(i)]);
        StoredValue output; output.data = std::make_shared<StoredValue::ListType>(std::move(result)); return output;
    }
    if (parts_.size() == 2 && base.is_list()) {
        const StoredValue index_value = parts_[1]->evaluate(ctx);
        if (!index_value.is_scalar() || mymath::floor(index_value.get_decimal()) != index_value.get_decimal())
            throw std::runtime_error("list index must be an integer");
        long long index = static_cast<long long>(index_value.get_decimal());
        if (index < 0) index += static_cast<long long>(base.list_value->size());
        if (index < 0 || index >= static_cast<long long>(base.list_value->size())) throw std::runtime_error("list index out of range");
        return (*base.list_value)[static_cast<std::size_t>(index)];
    }
    if (parts_.size() == 2 && base.is_dict()) {
        const StoredValue key = parts_[1]->evaluate(ctx);
        const std::string key_text = key.is_string() ? key.as_string() :
            std::to_string(static_cast<double>(key.get_decimal()));
        const auto it = base.dict_value->find(key_text);
        if (it == base.dict_value->end()) throw std::runtime_error("dictionary key not found: " + key_text);
        return it->second;
    }
    if (!base.is_matrix()) throw std::runtime_error("index access requires a matrix");
    if (parts_.size() != 2 && parts_.size() != 3)
        throw std::runtime_error("matrix index expects row and column");
    auto index = [&](std::size_t n) -> long long {
        StoredValue v = parts_[n]->evaluate(ctx);
        if (!v.is_scalar()) throw std::runtime_error("matrix index must be scalar");
        const Scalar raw = v.get_decimal();
        if (!mymath::isfinite(raw) || mymath::floor(raw) != raw) {
            throw std::runtime_error("matrix index must be an integer");
        }
        return static_cast<long long>(raw);
    };
    const auto& m = base.as_matrix();
    long long row = index(1);
    if (parts_.size() == 2) {
        if (row < 0) row += static_cast<long long>(m.rows);
        if (row >= static_cast<long long>(m.rows) || m.cols != 1) throw std::runtime_error("matrix index out of range");
        if (row < 0) throw std::runtime_error("matrix index out of range");
        return StoredValue(m.at(static_cast<std::size_t>(row), 0));
    }
    long long col = index(2);
    if (row < 0) row += static_cast<long long>(m.rows);
    if (col < 0) col += static_cast<long long>(m.cols);
    if (row >= static_cast<long long>(m.rows) || col >= static_cast<long long>(m.cols)) throw std::runtime_error("matrix index out of range");
    if (row < 0 || col < 0) throw std::runtime_error("matrix index out of range");
    return StoredValue(m.at(static_cast<std::size_t>(row), static_cast<std::size_t>(col)));
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
    StoredValue start_v = start_ ? start_->evaluate(ctx) : StoredValue();
    StoredValue stop_v = stop_ ? stop_->evaluate(ctx) : StoredValue();
    StoredValue step_v = step_ ? step_->evaluate(ctx) : StoredValue();

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
    case ExprKind::kListLiteral: {
        std::vector<std::unique_ptr<ASTNode>> elements;
        for (const auto& child : ast->children) elements.push_back(build_unified_ast(child.get()));
        return std::make_unique<ListLiteralNode>(std::move(elements));
    }
    case ExprKind::kIndexAccess: {
        std::vector<std::unique_ptr<ASTNode>> parts;
        for (const auto& child : ast->children) parts.push_back(build_unified_ast(child.get()));
        return std::make_unique<IndexAccessNode>(std::move(parts));
    }
    case ExprKind::kPostfixOp:
        if (!ast->children.empty()) return std::make_unique<PostfixOpNode>(
            ast->postfix_op, build_unified_ast(ast->children[0].get()));
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
    case ExprKind::kDictLiteral: {
        std::vector<std::pair<std::unique_ptr<ASTNode>, std::unique_ptr<ASTNode>>> entries;
        for (const auto& [key, value] : ast->dict_entries) {
            entries.emplace_back(build_unified_ast(key.get()), build_unified_ast(value.get()));
        }
        return std::make_unique<DictLiteralNode>(std::move(entries));
    }
    case ExprKind::kSlice:
        if (ast->children.size() == 3) return std::make_unique<SliceNode>(
            build_unified_ast(ast->children[0].get()), build_unified_ast(ast->children[1].get()),
            build_unified_ast(ast->children[2].get()));
        break;
    default:
        throw std::runtime_error("unsupported expression AST node at position " + std::to_string(ast->position));
    }
    throw std::runtime_error("invalid expression AST node at position " + std::to_string(ast->position));
}

StoredValue evaluate_unified_ast(const ExpressionAST* ast, ExecutionContext& ctx) {
    auto node = build_unified_ast(ast);
    if (!node) throw std::runtime_error("cannot build unified expression AST");
    return ASTEvaluator::evaluate(*node, ctx);
}

} // namespace core

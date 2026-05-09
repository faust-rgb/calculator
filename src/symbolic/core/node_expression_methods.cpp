// ============================================================================
// 向量/张量表达式实现
// ============================================================================

#include "symbolic/core/symbolic_expression_internal.h"

#include <algorithm>
#include <cctype>
#include <memory>
#include <vector>

using namespace symbolic_expression_internal;

// ============================================================================
// 向量/张量表达式实现
// ============================================================================

SymbolicExpression SymbolicExpression::vector(
    const std::vector<SymbolicExpression>& components) {
    auto node = std::make_shared<Node>();
    node->type = NodeType::kVector;
    node->shape = {components.size()};
    for (const SymbolicExpression& comp : components) {
        node->children.push_back(comp.node_);
    }
    return SymbolicExpression(node);
}

SymbolicExpression SymbolicExpression::tensor(
    const std::vector<std::vector<SymbolicExpression>>& rows) {
    auto node = std::make_shared<Node>();
    node->type = NodeType::kTensor;
    node->shape = {rows.size(), rows.empty() ? 0 : rows[0].size()};
    for (const auto& row : rows) {
        auto row_node = std::make_shared<Node>();
        row_node->type = NodeType::kVector;
        row_node->shape = {row.size()};
        for (const SymbolicExpression& comp : row) {
            row_node->children.push_back(comp.node_);
        }
        node->children.push_back(row_node);
    }
    return SymbolicExpression(node);
}

bool SymbolicExpression::is_vector() const {
    return node_->type == NodeType::kVector;
}

bool SymbolicExpression::is_tensor() const {
    return node_->type == NodeType::kTensor;
}

std::vector<SymbolicExpression> SymbolicExpression::vector_components() const {
    if (!is_vector()) {
        return {};
    }
    std::vector<SymbolicExpression> components;
    for (const auto& child : node_->children) {
        components.push_back(SymbolicExpression(child));
    }
    return components;
}

std::vector<std::vector<SymbolicExpression>> SymbolicExpression::tensor_rows() const {
    if (!is_tensor()) {
        return {};
    }
    std::vector<std::vector<SymbolicExpression>> rows;
    for (const auto& row_node : node_->children) {
        std::vector<SymbolicExpression> row;
        if (row_node->type == NodeType::kVector) {
            for (const auto& comp : row_node->children) {
                row.push_back(SymbolicExpression(comp));
            }
        }
        rows.push_back(row);
    }
    return rows;
}

std::vector<std::size_t> SymbolicExpression::get_shape() const {
    return node_->shape;
}

SymbolicExpression SymbolicExpression::substitute(
    const std::string& variable_name,
    const SymbolicExpression& replacement) const {
    if (!is_identifier_variable_name(variable_name) ||
        variable_name == "pi" || variable_name == "e" || variable_name == "i") {
        throw std::runtime_error(
            "symbolic substitution variable must be a non-reserved identifier");
    }
    return substitute_impl(*this, variable_name, replacement).simplify();
}

SymbolicExpression SymbolicExpression::substitute_expression(
    const SymbolicExpression& target,
    const SymbolicExpression& replacement) const {
    return substitute_expression_impl(*this, target, replacement).simplify();
}
// ============================================================================
// 边界参数解析实现
// ============================================================================

Scalar BoundArgument::to_scalar() const {
    switch (kind) {
        case BoundKind::kFinite:
            return value;
        case BoundKind::kPosInf:
            return mymath::precise128::infinity();
        case BoundKind::kNegInf:
            return -mymath::precise128::infinity();
    }
    return value;
}

BoundArgument BoundArgument::finite(Scalar v) {
    BoundArgument result;
    result.kind = BoundKind::kFinite;
    result.value = v;
    return result;
}

BoundArgument BoundArgument::pos_inf() {
    BoundArgument result;
    result.kind = BoundKind::kPosInf;
    return result;
}

BoundArgument BoundArgument::neg_inf() {
    BoundArgument result;
    result.kind = BoundKind::kNegInf;
    return result;
}

bool is_infinity_literal(const std::string& text) {
    std::string value = text;
    // 去除空白
    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.front()))) {
        value.erase(0, 1);
    }
    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.back()))) {
        value.pop_back();
    }
    // 去除正号
    if (!value.empty() && value.front() == '+') {
        value.erase(0, 1);
    }
    return value == "inf" || value == "infinity" || value == "oo";
}

BoundArgument parse_bound_argument(const std::string& text) {
    std::string value = text;
    // 去除空白
    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.front()))) {
        value.erase(0, 1);
    }
    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.back()))) {
        value.pop_back();
    }

    if (value.empty()) {
        throw std::runtime_error("empty bound argument");
    }

    // 检查符号
    bool negative = false;
    if (value.front() == '-') {
        negative = true;
        value.erase(0, 1);
        // 再次去除空白
        while (!value.empty() && std::isspace(static_cast<unsigned char>(value.front()))) {
            value.erase(0, 1);
        }
    } else if (value.front() == '+') {
        value.erase(0, 1);
        while (!value.empty() && std::isspace(static_cast<unsigned char>(value.front()))) {
            value.erase(0, 1);
        }
    }

    // 检查无穷大字面量
    if (value == "inf" || value == "infinity" || value == "oo") {
        return negative ? BoundArgument::neg_inf() : BoundArgument::pos_inf();
    }

    // 尝试解析为数值
    try {
        // 使用简单的数值解析
        Scalar num = Scalar(0.0L);
        bool has_digit = false;
        std::size_t i = 0;

        // 解析整数部分
        while (i < value.size() && std::isdigit(static_cast<unsigned char>(value[i]))) {
            num = num * 10.0L + (value[i] - '0');
            has_digit = true;
            ++i;
        }

        // 解析小数部分
        if (i < value.size() && value[i] == '.') {
            ++i;
            Scalar place = 0.1;
            while (i < value.size() && std::isdigit(static_cast<unsigned char>(value[i]))) {
                num += (value[i] - '0') * place;
                place *= 0.1;
                has_digit = true;
                ++i;
            }
        }

        // 解析科学计数法
        if (i < value.size() && (value[i] == 'e' || value[i] == 'E')) {
            ++i;
            bool exp_neg = false;
            if (i < value.size() && value[i] == '-') {
                exp_neg = true;
                ++i;
            } else if (i < value.size() && value[i] == '+') {
                ++i;
            }
            int exponent = 0;
            while (i < value.size() && std::isdigit(static_cast<unsigned char>(value[i]))) {
                exponent = exponent * 10 + (value[i] - '0');
                ++i;
            }
            if (exp_neg) {
                num /= mymath::pow(10.0L, exponent);
            } else {
                num *= mymath::pow(10.0L, exponent);
            }
        }

        if (!has_digit || i != value.size()) {
            throw std::runtime_error("invalid bound argument: " + text);
        }

        return BoundArgument::finite(negative ? -num : num);
    } catch (...) {
        throw std::runtime_error("invalid bound argument: " + text);
    }
}

bool expr_is_infinity(const SymbolicExpression& expression, bool* positive) {
    const auto& node = expression.node_;
    if (node->type == NodeType::kInfinity) {
        if (positive) {
            *positive = node->number_value > 0;
        }
        return true;
    }
    // 也检查数值无穷大
    if (node->type == NodeType::kNumber) {
        if (mymath::isinf(node->number_value)) {
            if (positive) {
                *positive = node->number_value > 0;
            }
            return true;
        }
    }
    return false;
}
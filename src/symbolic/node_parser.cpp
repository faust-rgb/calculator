#include "symbolic_expression_internal.h"

#include "mymath.h"

#include <algorithm>
#include <cctype>
#include <list>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace symbolic_expression_internal {

std::shared_ptr<SymbolicExpression::Node> intern_node(
    std::shared_ptr<SymbolicExpression::Node> node) {
    static thread_local std::unordered_map<std::string,
                                           std::weak_ptr<SymbolicExpression::Node>>
        interned_nodes;
    static constexpr std::size_t kMaxInternedNodes = 8192;

    const std::string key = node_structural_key(node);
    const auto found = interned_nodes.find(key);
    if (found != interned_nodes.end()) {
        if (std::shared_ptr<SymbolicExpression::Node> existing = found->second.lock()) {
            return existing;
        }
    }

    if (interned_nodes.size() >= kMaxInternedNodes) {
        for (auto it = interned_nodes.begin(); it != interned_nodes.end();) {
            if (it->second.expired()) {
                it = interned_nodes.erase(it);
            } else {
                ++it;
            }
        }
        if (interned_nodes.size() >= kMaxInternedNodes) {
            interned_nodes.clear();
        }
    }

    interned_nodes[key] = node;
    return node;
}

class SymbolicExpressionLruCache {
public:
    explicit SymbolicExpressionLruCache(std::size_t capacity)
        : capacity_(capacity) {}

    bool get(const std::string& key, SymbolicExpression* value) {
        const auto found = index_.find(key);
        if (found == index_.end()) {
            return false;
        }
        entries_.splice(entries_.begin(), entries_, found->second);
        *value = found->second->second;
        return true;
    }

    void put(const std::string& key, const SymbolicExpression& value) {
        const auto found = index_.find(key);
        if (found != index_.end()) {
            found->second->second = value;
            entries_.splice(entries_.begin(), entries_, found->second);
            return;
        }

        entries_.push_front({key, value});
        index_[key] = entries_.begin();
        while (entries_.size() > capacity_) {
            index_.erase(entries_.back().first);
            entries_.pop_back();
        }
    }

private:
    std::size_t capacity_ = 0;
    std::list<std::pair<std::string, SymbolicExpression>> entries_;
    std::unordered_map<std::string,
                       std::list<std::pair<std::string, SymbolicExpression>>::iterator>
        index_;
};

std::string format_number(double value) {
    // 符号输出会频繁生成中间常数。
    // 这里统一做“接近 0 归零”和“接近整数按整数打印”，
    // 避免输出里出现 -0、2.00000000001 这类噪声。
    if (mymath::is_near_zero(value, kFormatEps)) {
        return "0";
    }
    if (mymath::is_integer(value, 1e-10)) {
        long long rounded = static_cast<long long>(value >= 0.0 ? value + 0.5 : value - 0.5);
        return std::to_string(rounded);
    }

    long long numerator = 0;
    long long denominator = 1;
    if (mymath::approximate_fraction(value,
                                     &numerator,
                                     &denominator,
                                     999,
                                     1e-10)) {
        if (value < 0.0) {
            numerator = -numerator;
        }
        if (denominator == 1) {
            return std::to_string(numerator);
        }
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }

    std::ostringstream out;
    out.precision(12);
    out << value;
    return out.str();
}

std::shared_ptr<SymbolicExpression::Node> make_number(double value) {
    return intern_node(std::make_shared<SymbolicExpression::Node>(value));
}

std::shared_ptr<SymbolicExpression::Node> make_variable(const std::string& name) {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    node->type = NodeType::kVariable;
    node->text = name;
    return intern_node(node);
}

std::shared_ptr<SymbolicExpression::Node> make_unary(NodeType type,
                                                     std::shared_ptr<SymbolicExpression::Node> operand,
                                                     const std::string& text) {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    node->type = type;
    node->left = std::move(operand);
    node->text = text;
    return intern_node(node);
}

std::shared_ptr<SymbolicExpression::Node> make_binary(NodeType type,
                                                      std::shared_ptr<SymbolicExpression::Node> left,
                                                      std::shared_ptr<SymbolicExpression::Node> right) {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    node->type = type;
    node->left = std::move(left);
    node->right = std::move(right);
    return intern_node(node);
}

int precedence(const std::shared_ptr<SymbolicExpression::Node>& node) {
    // 这个优先级表直接决定 to_string() 时哪些子表达式要补括号。
    // 目标不是做“最少括号竞赛”，而是保证输出稳定、易读、可再次解析。
    switch (node->type) {
        case NodeType::kAdd:
        case NodeType::kSubtract:
            return 1;
        case NodeType::kMultiply:
        case NodeType::kDivide:
            return 2;
        case NodeType::kPower:
            return 3;
        case NodeType::kNegate:
            return 4;
        case NodeType::kFunction:
        case NodeType::kNumber:
        case NodeType::kVariable:
            return 5;
    }
    return 5;
}

std::string to_string_impl(const std::shared_ptr<SymbolicExpression::Node>& node, int parent_precedence) {
    std::string text;
    switch (node->type) {
        case NodeType::kNumber:
            text = format_number(node->number_value);
            break;
        case NodeType::kVariable:
            text = node->text;
            break;
        case NodeType::kNegate:
            text = "-" + to_string_impl(node->left, precedence(node));
            break;
        case NodeType::kFunction:
            text = node->text + "(" + to_string_impl(node->left, 0) + ")";
            break;
        case NodeType::kAdd:
            if (node->right->type == NodeType::kNegate) {
                text = to_string_impl(node->left, precedence(node)) + " - " +
                       to_string_impl(node->right->left, precedence(node) + 1);
            } else {
                text = to_string_impl(node->left, precedence(node)) + " + " +
                       to_string_impl(node->right, precedence(node));
            }
            break;
        case NodeType::kSubtract:
            text = to_string_impl(node->left, precedence(node)) + " - " +
                   to_string_impl(node->right, precedence(node) + 1);
            break;
        case NodeType::kMultiply:
            text = to_string_impl(node->left, precedence(node)) + " * " +
                   to_string_impl(node->right, precedence(node));
            break;
        case NodeType::kDivide:
            text = to_string_impl(node->left, precedence(node)) + " / " +
                   to_string_impl(node->right, precedence(node) + 1);
            break;
        case NodeType::kPower:
            text = to_string_impl(node->left, precedence(node)) + " ^ " +
                   to_string_impl(node->right, precedence(node));
            break;
    }

    if (node->type == NodeType::kNumber &&
        text.find('/') != std::string::npos &&
        parent_precedence >= 3) {
        return "(" + text + ")";
    }

    if (precedence(node) < parent_precedence) {
        return "(" + text + ")";
    }
    return text;
}

std::string node_structural_key(const std::shared_ptr<SymbolicExpression::Node>& node) {
    if (!node->structural_key_cache.empty()) {
        return node->structural_key_cache;
    }

    std::string key;
    switch (node->type) {
        case NodeType::kNumber:
            key = "N(" + format_number(node->number_value) + ")";
            break;
        case NodeType::kVariable:
            key = "V(" + node->text + ")";
            break;
        case NodeType::kNegate:
            key = "NEG(" + node_structural_key(node->left) + ")";
            break;
        case NodeType::kFunction:
            key = "F(" + node->text + ":" + node_structural_key(node->left) + ")";
            break;
        case NodeType::kAdd:
            key = "ADD(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
        case NodeType::kSubtract:
            key = "SUB(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
        case NodeType::kMultiply:
            key = "MUL(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
        case NodeType::kDivide:
            key = "DIV(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
        case NodeType::kPower:
            key = "POW(" + node_structural_key(node->left) + "," +
                  node_structural_key(node->right) + ")";
            break;
    }
    node->structural_key_cache = key;
    return node->structural_key_cache;
}

class Parser {
public:
    explicit Parser(std::string source) : source_(std::move(source)) {}

    SymbolicExpression parse() {
        SymbolicExpression expression = parse_expression();
        skip_spaces();
        if (pos_ != source_.size()) {
            throw std::runtime_error("unexpected token near: " + source_.substr(pos_, 1));
        }
        return expression;
    }

private:
    SymbolicExpression parse_expression() {
        // 解析器采用经典的递归下降结构：
        // expression -> term ((+|-) term)*
        // term       -> unary ((*|/) unary)*
        // unary      -> (+|-) unary | power
        // power      -> primary (^ unary)?
        //
        // 它和 src/core/decimal_parser.cpp 里的数值解析器保持同样的运算优先级，
        // 这样符号路径和数值路径看到的是一致的表达式结构。
        SymbolicExpression value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = SymbolicExpression(make_binary(NodeType::kAdd, value.simplify().node_, parse_term().simplify().node_));
            } else if (match('-')) {
                value = SymbolicExpression(make_binary(NodeType::kSubtract, value.simplify().node_, parse_term().simplify().node_));
            } else {
                return value;
            }
        }
    }

    SymbolicExpression parse_term() {
        SymbolicExpression value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = SymbolicExpression(make_binary(NodeType::kMultiply, value.node_, parse_unary().node_));
            } else if (match('/')) {
                value = SymbolicExpression(make_binary(NodeType::kDivide, value.node_, parse_unary().node_));
            } else {
                return value;
            }
        }
    }

    SymbolicExpression parse_power() {
        SymbolicExpression value = parse_primary();
        skip_spaces();
        if (match('^')) {
            return SymbolicExpression(make_binary(NodeType::kPower, value.node_, parse_unary().node_));
        }
        return value;
    }

    SymbolicExpression parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            return SymbolicExpression(make_unary(NodeType::kNegate, parse_unary().node_));
        }
        return parse_power();
    }

    SymbolicExpression parse_primary() {
        skip_spaces();
        if (match('(')) {
            SymbolicExpression value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        if (peek_is_alpha()) {
            std::string identifier = parse_identifier();
            if (identifier == "pi") {
                return SymbolicExpression(make_variable("pi"));
            }
            if (identifier == "e") {
                return SymbolicExpression(make_variable("e"));
            }

            skip_spaces();
            if (match('(')) {
                if (identifier == "u" || identifier == "heaviside") {
                    identifier = "step";
                } else if (identifier == "impulse") {
                    identifier = "delta";
                }
                SymbolicExpression argument = parse_expression();
                skip_spaces();
                expect(')');
                return SymbolicExpression(make_unary(NodeType::kFunction, argument.node_, identifier));
            }
            return SymbolicExpression(make_variable(identifier));
        }

        return parse_number();
    }

    SymbolicExpression parse_number() {
        skip_spaces();
        const std::size_t start = pos_;
        bool has_digit = false;
        bool seen_dot = false;

        while (pos_ < source_.size()) {
            const char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                has_digit = true;
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else {
                break;
            }
        }

        if (!has_digit) {
            throw std::runtime_error("expected number");
        }

        double value = 0.0;
        std::size_t index = start;
        while (index < pos_ && source_[index] != '.') {
            value = value * 10.0 + static_cast<double>(source_[index] - '0');
            ++index;
        }
        if (index < pos_ && source_[index] == '.') {
            ++index;
            double place = 0.1;
            while (index < pos_) {
                value += static_cast<double>(source_[index] - '0') * place;
                place *= 0.1;
                ++index;
            }
        }

        return SymbolicExpression(make_number(value));
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (pos_ < source_.size()) {
            const char ch = source_[pos_];
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
                ++pos_;
            } else {
                break;
            }
        }
        return source_.substr(start, pos_ - start);
    }

    bool peek_is_alpha() const {
        return pos_ < source_.size() &&
               std::isalpha(static_cast<unsigned char>(source_[pos_]));
    }

    bool match(char ch) {
        if (pos_ >= source_.size() || source_[pos_] != ch) {
            return false;
        }
        ++pos_;
        return true;
    }

    void expect(char ch) {
        if (!match(ch)) {
            throw std::runtime_error(std::string("expected '") + ch + "'");
        }
    }

    void skip_spaces() {
        while (pos_ < source_.size() &&
               std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }

    std::string source_;
    std::size_t pos_ = 0;
};

bool expr_is_number(const SymbolicExpression& expression, double* value);

SymbolicExpression simplify_impl(const SymbolicExpression& expression);
SymbolicExpression simplify_once(const SymbolicExpression& expression);

SymbolicExpression substitute_impl(const SymbolicExpression& expression,
                                  const std::string& variable_name,
                                  const SymbolicExpression& replacement) {
    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
            return expression;
        case NodeType::kVariable:
            if (node->text == variable_name) {
                return replacement;
            }
            return expression;
        case NodeType::kNegate:
            return SymbolicExpression(
                       make_unary(NodeType::kNegate,
                                  substitute_impl(SymbolicExpression(node->left),
                                                  variable_name,
                                                  replacement)
                                      .node_))
                .simplify();
        case NodeType::kFunction:
            return SymbolicExpression(
                       make_unary(NodeType::kFunction,
                                  substitute_impl(SymbolicExpression(node->left),
                                                  variable_name,
                                                  replacement)
                                      .node_,
                                  node->text))
                .simplify();
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            return SymbolicExpression(
                       make_binary(node->type,
                                   substitute_impl(SymbolicExpression(node->left),
                                                   variable_name,
                                                   replacement)
                                       .node_,
                                   substitute_impl(SymbolicExpression(node->right),
                                                   variable_name,
                                                   replacement)
                                       .node_))
                .simplify();
    }
    throw std::runtime_error("unsupported symbolic substitution");
}

bool try_evaluate_numeric_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                               double* value) {
    switch (node->type) {
        case NodeType::kNumber:
            *value = node->number_value;
            return true;
        case NodeType::kVariable:
            if (node->text == "1 / 2") {
                *value = 0.5;
                return true;
            }
            return false;
        case NodeType::kNegate: {
            double operand = 0.0;
            if (!try_evaluate_numeric_node(node->left, &operand)) {
                return false;
            }
            *value = -operand;
            return true;
        }
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower: {
            double left = 0.0;
            double right = 0.0;
            if (!try_evaluate_numeric_node(node->left, &left) ||
                !try_evaluate_numeric_node(node->right, &right)) {
                return false;
            }

            switch (node->type) {
                case NodeType::kAdd:
                    *value = left + right;
                    return true;
                case NodeType::kSubtract:
                    *value = left - right;
                    return true;
                case NodeType::kMultiply:
                    *value = left * right;
                    return true;
                case NodeType::kDivide:
                    *value = left / right;
                    return true;
                case NodeType::kPower:
                    *value = mymath::pow(left, right);
                    return true;
                case NodeType::kNumber:
                case NodeType::kVariable:
                case NodeType::kNegate:
                case NodeType::kFunction:
                    break;
            }
            return false;
        }
        case NodeType::kFunction: {
            double argument = 0.0;
            if (!try_evaluate_numeric_node(node->left, &argument)) {
                return false;
            }
            if (node->text == "asin") {
                *value = mymath::asin(argument);
                return true;
            }
            if (node->text == "acos") {
                *value = mymath::acos(argument);
                return true;
            }
            if (node->text == "atan") {
                *value = mymath::atan(argument);
                return true;
            }
            if (node->text == "sin") {
                *value = mymath::sin(argument);
                return true;
            }
            if (node->text == "cos") {
                *value = mymath::cos(argument);
                return true;
            }
            if (node->text == "tan") {
                *value = mymath::tan(argument);
                return true;
            }
            if (node->text == "exp") {
                *value = mymath::exp(argument);
                return true;
            }
            if (node->text == "sinh") {
                *value = mymath::sinh(argument);
                return true;
            }
            if (node->text == "cosh") {
                *value = mymath::cosh(argument);
                return true;
            }
            if (node->text == "tanh") {
                *value = mymath::tanh(argument);
                return true;
            }
            if (node->text == "ln") {
                *value = mymath::ln(argument);
                return true;
            }
            if (node->text == "sqrt") {
                *value = mymath::sqrt(argument);
                return true;
            }
            if (node->text == "abs") {
                *value = mymath::abs(argument);
                return true;
            }
            if (node->text == "floor") {
                *value = static_cast<double>(
                    static_cast<long long>(argument < 0.0 &&
                                                   static_cast<double>(static_cast<long long>(argument)) != argument
                                               ? argument - 1.0
                                               : argument));
                return true;
            }
            if (node->text == "ceil") {
                long long truncated = static_cast<long long>(argument);
                if (argument > 0.0 && static_cast<double>(truncated) != argument) {
                    ++truncated;
                }
                *value = static_cast<double>(truncated);
                return true;
            }
            if (node->text == "cbrt") {
                *value = mymath::cbrt(argument);
                return true;
            }
            if (node->text == "sign") {
                if (mymath::is_near_zero(argument, kFormatEps)) {
                    *value = 0.0;
                } else {
                    *value = argument > 0.0 ? 1.0 : -1.0;
                }
                return true;
            }
            if (node->text == "step") {
                *value = argument >= 0.0 ? 1.0 : 0.0;
                return true;
            }
            if (node->text == "delta") {
                *value = mymath::is_near_zero(argument, kFormatEps) ? 1.0 : 0.0;
                return true;
            }
            return false;
        }
    }
    return false;
}

bool expr_is_variable(const SymbolicExpression& expression, const std::string& name) {
    return expression.node_->type == NodeType::kVariable && expression.node_->text == name;
}

}  // namespace symbolic_expression_internal

using namespace symbolic_expression_internal;

SymbolicExpression::SymbolicExpression() : node_(make_number(0.0)) {}

SymbolicExpression::SymbolicExpression(std::shared_ptr<Node> node)
    : node_(std::move(node)) {}

SymbolicExpression SymbolicExpression::parse(const std::string& text) {
    Parser parser(text);
    return parser.parse().simplify();
}

SymbolicExpression SymbolicExpression::number(double value) {
    return SymbolicExpression(make_number(value));
}

SymbolicExpression SymbolicExpression::variable(const std::string& name) {
    return SymbolicExpression(make_variable(name));
}

std::string SymbolicExpression::to_string() const {
    return to_string_impl(simplify().node_, 0);
}

bool SymbolicExpression::is_constant(const std::string& variable_name) const {
    switch (node_->type) {
        case NodeType::kNumber:
            return true;
        case NodeType::kVariable:
            return node_->text != variable_name;
        case NodeType::kNegate:
        case NodeType::kFunction:
            return SymbolicExpression(node_->left).is_constant(variable_name);
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            return SymbolicExpression(node_->left).is_constant(variable_name) &&
                   SymbolicExpression(node_->right).is_constant(variable_name);
    }
    return false;
}

bool SymbolicExpression::is_number(double* value) const {
    double numeric = 0.0;
    if (!try_evaluate_numeric_node(node_, &numeric)) {
        return false;
    }
    if (value != nullptr) {
        *value = numeric;
    }
    return true;
}

bool SymbolicExpression::is_variable_named(const std::string& variable_name) const {
    return node_->type == NodeType::kVariable && node_->text == variable_name;
}

bool SymbolicExpression::polynomial_coefficients(
    const std::string& variable_name,
    std::vector<double>* coefficients) const {
    const SymbolicExpression simplified = simplify();
    return polynomial_coefficients_from_simplified(simplified, variable_name, coefficients);
}

std::vector<std::string> SymbolicExpression::identifier_variables() const {
    std::vector<std::string> names;
    collect_identifier_variables(simplify(), &names);
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
}

SymbolicExpression SymbolicExpression::simplify() const {
    static constexpr std::size_t kMaxSimplifyCacheSize = 4096;
    static thread_local SymbolicExpressionLruCache cache(kMaxSimplifyCacheSize);

    const std::string key = node_structural_key(node_);
    SymbolicExpression cached;
    if (cache.get(key, &cached)) {
        return cached;
    }

    SymbolicExpression simplified = simplify_impl(*this);
    cache.put(key, simplified);
    return simplified;
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

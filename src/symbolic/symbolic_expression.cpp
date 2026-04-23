/**
 * @file symbolic_expression.cpp
 * @brief 符号表达式实现
 *
 * 实现符号表达式的解析、操作和简化。
 * 支持符号微分、积分、代数简化和多项式分析。
 */

#include "symbolic_expression.h"

#include "mymath.h"
#include "polynomial.h"

#include <cctype>
#include <memory>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/** @brief 表达式节点类型枚举 */
enum class NodeType {
    kNumber,     ///< 数字常量
    kVariable,   ///< 变量
    kAdd,        ///< 加法
    kSubtract,   ///< 减法
    kMultiply,   ///< 乘法
    kDivide,     ///< 除法
    kPower,      ///< 幂运算
    kNegate,     ///< 取负
    kFunction,   ///< 函数调用
};

struct SymbolicExpression::Node {
    // 这是一个非常轻量的表达式树节点：
    // - number_value 用于常数
    // - text 用于变量名或函数名
    // - left/right 用于一元或二元表达式的子树
    //
    // 这里不单独拆很多不同节点类型，而是统一由 NodeType 驱动，
    // 这样实现求导、积分、化简时可以直接按 kind 分派。
    NodeType type = NodeType::kNumber;
    double number_value = 0.0;
    std::string text;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;

    Node() = default;
    explicit Node(double value) : type(NodeType::kNumber), number_value(value) {}
};

namespace {

constexpr double kFormatEps = 1e-12;

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
    return std::make_shared<SymbolicExpression::Node>(value);
}

std::shared_ptr<SymbolicExpression::Node> make_variable(const std::string& name) {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    node->type = NodeType::kVariable;
    node->text = name;
    return node;
}

std::shared_ptr<SymbolicExpression::Node> make_unary(NodeType type,
                                                     std::shared_ptr<SymbolicExpression::Node> operand,
                                                     const std::string& text = "") {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    node->type = type;
    node->left = std::move(operand);
    node->text = text;
    return node;
}

std::shared_ptr<SymbolicExpression::Node> make_binary(NodeType type,
                                                      std::shared_ptr<SymbolicExpression::Node> left,
                                                      std::shared_ptr<SymbolicExpression::Node> right) {
    std::shared_ptr<SymbolicExpression::Node> node =
        std::make_shared<SymbolicExpression::Node>();
    node->type = type;
    node->left = std::move(left);
    node->right = std::move(right);
    return node;
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
            text = to_string_impl(node->left, precedence(node)) + " + " +
                   to_string_impl(node->right, precedence(node));
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
    switch (node->type) {
        case NodeType::kNumber:
            return "N(" + format_number(node->number_value) + ")";
        case NodeType::kVariable:
            return "V(" + node->text + ")";
        case NodeType::kNegate:
            return "NEG(" + node_structural_key(node->left) + ")";
        case NodeType::kFunction:
            return "F(" + node->text + ":" + node_structural_key(node->left) + ")";
        case NodeType::kAdd:
            return "ADD(" + node_structural_key(node->left) + "," +
                   node_structural_key(node->right) + ")";
        case NodeType::kSubtract:
            return "SUB(" + node_structural_key(node->left) + "," +
                   node_structural_key(node->right) + ")";
        case NodeType::kMultiply:
            return "MUL(" + node_structural_key(node->left) + "," +
                   node_structural_key(node->right) + ")";
        case NodeType::kDivide:
            return "DIV(" + node_structural_key(node->left) + "," +
                   node_structural_key(node->right) + ")";
        case NodeType::kPower:
            return "POW(" + node_structural_key(node->left) + "," +
                   node_structural_key(node->right) + ")";
    }
    return "";
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
        // 它和 calculator.cpp 里的数值解析器保持同样的运算优先级，
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

bool expr_is_number(const SymbolicExpression& expression, double* value = nullptr);

SymbolicExpression simplify_impl(const SymbolicExpression& expression);

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

bool expr_is_zero(const SymbolicExpression& expression) {
    double value = 0.0;
    return expr_is_number(expression, &value) && mymath::is_near_zero(value, kFormatEps);
}

bool expr_is_one(const SymbolicExpression& expression) {
    double value = 0.0;
    return expr_is_number(expression, &value) && mymath::is_near_zero(value - 1.0, kFormatEps);
}

bool expr_is_minus_one(const SymbolicExpression& expression) {
    double value = 0.0;
    return expr_is_number(expression, &value) && mymath::is_near_zero(value + 1.0, kFormatEps);
}

bool expr_is_number(const SymbolicExpression& expression, double* value) {
    return expression.is_number(value);
}

bool decompose_numeric_multiple_of_symbol(const SymbolicExpression& expression,
                                          const std::string& symbol_name,
                                          double* coefficient) {
    const SymbolicExpression simplified = expression.simplify();
    if (expr_is_variable(simplified, symbol_name)) {
        *coefficient = 1.0;
        return true;
    }

    const auto& node = simplified.node_;
    if (node->type == NodeType::kNegate) {
        double nested = 0.0;
        if (decompose_numeric_multiple_of_symbol(SymbolicExpression(node->left),
                                                 symbol_name,
                                                 &nested)) {
            *coefficient = -nested;
            return true;
        }
        return false;
    }

    if (node->type == NodeType::kMultiply) {
        double numeric = 0.0;
        if (SymbolicExpression(node->left).is_number(&numeric) &&
            decompose_numeric_multiple_of_symbol(SymbolicExpression(node->right),
                                                 symbol_name,
                                                 coefficient)) {
            *coefficient *= numeric;
            return true;
        }
        if (SymbolicExpression(node->right).is_number(&numeric) &&
            decompose_numeric_multiple_of_symbol(SymbolicExpression(node->left),
                                                 symbol_name,
                                                 coefficient)) {
            *coefficient *= numeric;
            return true;
        }
        return false;
    }

    if (node->type == NodeType::kDivide) {
        double divisor = 0.0;
        if (!SymbolicExpression(node->right).is_number(&divisor) ||
            mymath::is_near_zero(divisor, kFormatEps)) {
            return false;
        }
        if (decompose_numeric_multiple_of_symbol(SymbolicExpression(node->left),
                                                 symbol_name,
                                                 coefficient)) {
            *coefficient /= divisor;
            return true;
        }
        return false;
    }

    return false;
}

bool numeric_matches_any(double value, const std::initializer_list<double>& candidates) {
    for (double candidate : candidates) {
        if (mymath::is_near_zero(value - candidate, kFormatEps)) {
            return true;
        }
    }
    return false;
}

SymbolicExpression sqrt3_symbol() {
    return SymbolicExpression::variable("sqrt(3)");
}

SymbolicExpression half_symbol() {
    return SymbolicExpression::variable("1 / 2");
}

bool decompose_constant_times_expression(const SymbolicExpression& expression,
                                         const std::string& variable_name,
                                         double* constant,
                                         SymbolicExpression* rest) {
    double numeric = 0.0;
    if (expression.is_constant(variable_name) && expression.is_number(&numeric)) {
        *constant = numeric;
        *rest = SymbolicExpression::number(1.0);
        return true;
    }

    const auto& node = expression.simplify().node_;
    if (node->type != NodeType::kMultiply) {
        return false;
    }

    SymbolicExpression left(node->left);
    SymbolicExpression right(node->right);
    if (left.is_constant(variable_name) && left.is_number(constant)) {
        *rest = right;
        return true;
    }
    if (right.is_constant(variable_name) && right.is_number(constant)) {
        *rest = left;
        return true;
    }
    return false;
}

void collect_multiplicative_terms(const SymbolicExpression& expression,
                                  double* numeric_factor,
                                  std::vector<SymbolicExpression>* symbolic_factors) {
    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        *numeric_factor *= numeric;
        return;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kNegate) {
        *numeric_factor *= -1.0;
        collect_multiplicative_terms(SymbolicExpression(node->left).simplify(),
                                     numeric_factor,
                                     symbolic_factors);
        return;
    }

    if (node->type == NodeType::kMultiply) {
        collect_multiplicative_terms(SymbolicExpression(node->left).simplify(),
                                     numeric_factor,
                                     symbolic_factors);
        collect_multiplicative_terms(SymbolicExpression(node->right).simplify(),
                                     numeric_factor,
                                     symbolic_factors);
        return;
    }

    symbolic_factors->push_back(expression);
}

void collect_division_factors(const SymbolicExpression& expression,
                              double* numeric_factor,
                              std::vector<SymbolicExpression>* symbolic_factors) {
    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        *numeric_factor *= numeric;
        return;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kNegate) {
        *numeric_factor *= -1.0;
        collect_division_factors(SymbolicExpression(node->left).simplify(),
                                 numeric_factor,
                                 symbolic_factors);
        return;
    }
    if (node->type == NodeType::kMultiply) {
        collect_division_factors(SymbolicExpression(node->left).simplify(),
                                 numeric_factor,
                                 symbolic_factors);
        collect_division_factors(SymbolicExpression(node->right).simplify(),
                                 numeric_factor,
                                 symbolic_factors);
        return;
    }
    if (node->type == NodeType::kPower) {
        double exponent = 0.0;
        if (SymbolicExpression(node->right).is_number(&exponent) &&
            mymath::is_integer(exponent, 1e-10) &&
            exponent > 0.0) {
            const int count = static_cast<int>(exponent + 0.5);
            const SymbolicExpression base = SymbolicExpression(node->left).simplify();
            for (int i = 0; i < count; ++i) {
                symbolic_factors->push_back(base);
            }
            return;
        }
    }

    symbolic_factors->push_back(expression);
}

SymbolicExpression rebuild_product_expression(double numeric_factor,
                                              const std::vector<SymbolicExpression>& factors) {
    if (mymath::is_near_zero(numeric_factor, kFormatEps)) {
        return SymbolicExpression::number(0.0);
    }

    SymbolicExpression combined;
    bool has_combined = false;
    for (const SymbolicExpression& factor : factors) {
        if (!has_combined) {
            combined = factor;
            has_combined = true;
        } else {
            combined = SymbolicExpression(
                make_binary(NodeType::kMultiply, combined.node_, factor.node_));
        }
    }

    if (!has_combined) {
        return SymbolicExpression::number(numeric_factor);
    }
    if (mymath::is_near_zero(numeric_factor - 1.0, kFormatEps)) {
        return combined;
    }
    if (mymath::is_near_zero(numeric_factor + 1.0, kFormatEps)) {
        return SymbolicExpression(make_unary(NodeType::kNegate, combined.node_)).simplify();
    }
    return SymbolicExpression(
               make_binary(NodeType::kMultiply,
                           SymbolicExpression::number(numeric_factor).node_,
                           combined.node_))
        .simplify();
}

bool decompose_numeric_factor(const SymbolicExpression& expression,
                              double* coefficient,
                              SymbolicExpression* rest) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric_factor = 1.0;
    std::vector<SymbolicExpression> symbolic_factors;
    collect_multiplicative_terms(simplified, &numeric_factor, &symbolic_factors);

    *coefficient = numeric_factor;
    if (symbolic_factors.empty()) {
        *rest = SymbolicExpression::number(1.0);
        return true;
    }

    SymbolicExpression combined = symbolic_factors.front();
    for (std::size_t i = 1; i < symbolic_factors.size(); ++i) {
        combined = SymbolicExpression(
            make_binary(NodeType::kMultiply, combined.node_, symbolic_factors[i].node_));
    }
    *rest = combined.simplify();
    return true;
}

std::string canonical_multiplicative_key(const SymbolicExpression& expression) {
    double ignored = 1.0;
    const SymbolicExpression simplified = expression.simplify();
    std::vector<SymbolicExpression> symbolic_factors;
    collect_multiplicative_terms(simplified, &ignored, &symbolic_factors);
    std::vector<std::string> parts;
    parts.reserve(symbolic_factors.size());
    for (const SymbolicExpression& factor : symbolic_factors) {
        parts.push_back(node_structural_key(factor.node_));
    }
    std::sort(parts.begin(), parts.end());

    std::ostringstream out;
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i != 0) {
            out << " * ";
        }
        out << parts[i];
    }
    return out.str();
}

void collect_additive_terms(const SymbolicExpression& expression,
                            std::vector<std::string>* parts) {
    const auto& node = expression.node_;
    if (node->type == NodeType::kAdd) {
        collect_additive_terms(SymbolicExpression(node->left).simplify(), parts);
        collect_additive_terms(SymbolicExpression(node->right).simplify(), parts);
        return;
    }

    parts->push_back(node_structural_key(node));
}

std::string canonical_expression_key(const SymbolicExpression& expression) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type == NodeType::kAdd) {
        std::vector<std::string> parts;
        collect_additive_terms(simplified, &parts);
        std::sort(parts.begin(), parts.end());
        std::ostringstream out;
        for (std::size_t i = 0; i < parts.size(); ++i) {
            if (i != 0) {
                out << " + ";
            }
            out << parts[i];
        }
        return out.str();
    }
    return canonical_multiplicative_key(simplified);
}

bool try_combine_like_terms(const SymbolicExpression& left,
                            const SymbolicExpression& right,
                            double right_sign,
                            SymbolicExpression* combined) {
    double left_coefficient = 0.0;
    double right_coefficient = 0.0;
    SymbolicExpression left_rest;
    SymbolicExpression right_rest;
    if (!decompose_numeric_factor(left, &left_coefficient, &left_rest) ||
        !decompose_numeric_factor(right, &right_coefficient, &right_rest)) {
        return false;
    }

    if (canonical_expression_key(left_rest) !=
        canonical_expression_key(right_rest)) {
        return false;
    }

    const double result_coefficient =
        left_coefficient + right_sign * right_coefficient;
    if (mymath::is_near_zero(result_coefficient, kFormatEps)) {
        *combined = SymbolicExpression::number(0.0);
        return true;
    }
    if (left_rest.is_number()) {
        *combined = SymbolicExpression::number(result_coefficient);
        return true;
    }
    if (mymath::is_near_zero(result_coefficient - 1.0, kFormatEps)) {
        *combined = left_rest;
        return true;
    }
    if (mymath::is_near_zero(result_coefficient + 1.0, kFormatEps)) {
        *combined = SymbolicExpression(
                        make_unary(NodeType::kNegate, left_rest.node_))
                        .simplify();
        return true;
    }
    *combined = SymbolicExpression(
                    make_binary(NodeType::kMultiply,
                                SymbolicExpression::number(result_coefficient).node_,
                                left_rest.node_))
                    .simplify();
    return true;
}

bool decompose_linear(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      double* coefficient,
                      double* intercept) {
    // 把表达式尝试识别成 a*x + b 的形式。
    // 这是符号积分里处理 sin(ax+b)、cos(ax+b)、exp(ax+b)、1/(ax+b)
    // 这些常见模式的基础工具函数。
    const SymbolicExpression simplified = expression.simplify();
    double number = 0.0;
    if (simplified.is_variable_named(variable_name)) {
        *coefficient = 1.0;
        *intercept = 0.0;
        return true;
    }
    if (simplified.is_number(&number)) {
        *coefficient = 0.0;
        *intercept = number;
        return true;
    }

    const auto& node = simplified.node_;
    if (node->type == NodeType::kNegate) {
        if (decompose_linear(SymbolicExpression(node->left), variable_name, coefficient, intercept)) {
            *coefficient = -*coefficient;
            *intercept = -*intercept;
            return true;
        }
    }
    if (node->type == NodeType::kAdd || node->type == NodeType::kSubtract) {
        double left_a = 0.0;
        double left_b = 0.0;
        double right_a = 0.0;
        double right_b = 0.0;
        if (!decompose_linear(SymbolicExpression(node->left), variable_name, &left_a, &left_b) ||
            !decompose_linear(SymbolicExpression(node->right), variable_name, &right_a, &right_b)) {
            return false;
        }
        *coefficient = left_a + (node->type == NodeType::kAdd ? right_a : -right_a);
        *intercept = left_b + (node->type == NodeType::kAdd ? right_b : -right_b);
        return true;
    }
    if (node->type == NodeType::kMultiply) {
        double factor = 0.0;
        SymbolicExpression rest;
        if (decompose_constant_times_expression(simplified, variable_name, &factor, &rest) &&
            rest.is_variable_named(variable_name)) {
            *coefficient = factor;
            *intercept = 0.0;
            return true;
        }
    }

    return false;
}

void trim_polynomial_coefficients(std::vector<double>* coefficients) {
    while (coefficients->size() > 1 &&
           mymath::is_near_zero(coefficients->back(), kFormatEps)) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(0.0);
    }
}

std::vector<double> polynomial_add_impl(const std::vector<double>& lhs,
                                        const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] += lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] += rhs[i];
    }
    trim_polynomial_coefficients(&result);
    return result;
}

std::vector<double> polynomial_subtract_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] += lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] -= rhs[i];
    }
    trim_polynomial_coefficients(&result);
    return result;
}

std::vector<double> polynomial_multiply_impl(const std::vector<double>& lhs,
                                             const std::vector<double>& rhs) {
    std::vector<double> result(lhs.size() + rhs.size() - 1, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result[i + j] += lhs[i] * rhs[j];
        }
    }
    trim_polynomial_coefficients(&result);
    return result;
}

SymbolicExpression make_add(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kAdd, left.node_, right.node_));
}

SymbolicExpression make_subtract(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kSubtract, left.node_, right.node_));
}

SymbolicExpression make_multiply(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kMultiply, left.node_, right.node_));
}

SymbolicExpression make_divide(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kDivide, left.node_, right.node_));
}

SymbolicExpression make_power(SymbolicExpression left, SymbolicExpression right) {
    return SymbolicExpression(make_binary(NodeType::kPower, left.node_, right.node_));
}

SymbolicExpression make_negate(SymbolicExpression operand) {
    return SymbolicExpression(make_unary(NodeType::kNegate, operand.node_));
}

SymbolicExpression make_function(const std::string& name, SymbolicExpression argument) {
    return SymbolicExpression(make_unary(NodeType::kFunction, argument.node_, name));
}

SymbolicExpression build_polynomial_expression_from_coefficients(
    const std::vector<double>& coefficients,
    const std::string& variable_name) {
    std::vector<double> normalized = coefficients;
    trim_polynomial_coefficients(&normalized);

    SymbolicExpression result = SymbolicExpression::number(0.0);
    bool has_term = false;
    for (std::size_t index = normalized.size(); index > 0; --index) {
        const std::size_t degree = index - 1;
        const double coefficient = normalized[degree];
        if (mymath::is_near_zero(coefficient, kFormatEps)) {
            continue;
        }

        const bool negative = coefficient < 0.0;
        const double magnitude = negative ? -coefficient : coefficient;

        SymbolicExpression term;
        if (degree == 0) {
            term = SymbolicExpression::number(magnitude);
        } else {
            term = degree == 1
                       ? SymbolicExpression::variable(variable_name)
                       : make_power(SymbolicExpression::variable(variable_name),
                                    SymbolicExpression::number(
                                        static_cast<double>(degree)));
            if (!mymath::is_near_zero(magnitude - 1.0, kFormatEps)) {
                term = make_multiply(SymbolicExpression::number(magnitude), term);
            }
        }

        if (!has_term) {
            result = negative ? make_negate(term) : term;
            has_term = true;
            continue;
        }

        result = negative ? make_subtract(result, term) : make_add(result, term);
    }

    return has_term ? result : SymbolicExpression::number(0.0);
}

bool expressions_match(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return node_structural_key(lhs.node_) == node_structural_key(rhs.node_);
}

bool is_known_positive_expression(const SymbolicExpression& expression) {
    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        return numeric > 0.0;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kVariable) {
        return node->text == "e" || node->text == "pi";
    }
    if (node->type == NodeType::kFunction) {
        return node->text == "exp" || node->text == "sqrt" || node->text == "abs";
    }
    return false;
}

bool polynomial_expression(const SymbolicExpression& expression,
                           const std::string& variable_name,
                           SymbolicExpression* polynomial) {
    std::vector<double> coefficients;
    if (!expression.polynomial_coefficients(variable_name, &coefficients)) {
        return false;
    }
    *polynomial = build_polynomial_expression_from_coefficients(coefficients,
                                                                variable_name);
    return true;
}

bool is_linear_function_argument(const SymbolicExpression& argument,
                                 const std::string& variable_name,
                                 double* a) {
    double b = 0.0;
    return decompose_linear(argument, variable_name, a, &b) &&
           !mymath::is_near_zero(*a, kFormatEps);
}

bool integrate_polynomial_times_function(const SymbolicExpression& polynomial,
                                         const std::string& function_name,
                                         const SymbolicExpression& argument,
                                         const std::string& variable_name,
                                         SymbolicExpression* integrated) {
    double a = 0.0;
    if (!is_linear_function_argument(argument, variable_name, &a)) {
        return false;
    }

    double constant = 0.0;
    if (polynomial.is_number(&constant)) {
        if (function_name == "exp") {
            *integrated = make_multiply(
                              SymbolicExpression::number(constant),
                              make_divide(make_function("exp", argument),
                                          SymbolicExpression::number(a)))
                              .simplify();
            return true;
        }
        if (function_name == "sin") {
            *integrated = make_multiply(
                              SymbolicExpression::number(constant),
                              make_divide(make_negate(make_function("cos", argument)),
                                          SymbolicExpression::number(a)))
                              .simplify();
            return true;
        }
        if (function_name == "cos") {
            *integrated = make_multiply(
                              SymbolicExpression::number(constant),
                              make_divide(make_function("sin", argument),
                                          SymbolicExpression::number(a)))
                              .simplify();
            return true;
        }
        return false;
    }

    const SymbolicExpression derivative = polynomial.derivative(variable_name).simplify();
    SymbolicExpression recursive;
    if (!integrate_polynomial_times_function(derivative,
                                             function_name,
                                             argument,
                                             variable_name,
                                             &recursive)) {
        return false;
    }

    if (function_name == "exp") {
        *integrated = make_subtract(
                          make_divide(make_multiply(polynomial,
                                                    make_function("exp", argument)),
                                      SymbolicExpression::number(a)),
                          make_divide(recursive, SymbolicExpression::number(a)))
                          .simplify();
        return true;
    }
    if (function_name == "sin") {
        *integrated = make_add(
                          make_divide(make_negate(make_multiply(polynomial,
                                                                make_function("cos", argument))),
                                      SymbolicExpression::number(a)),
                          make_divide(recursive, SymbolicExpression::number(a)))
                          .simplify();
        return true;
    }
    if (function_name == "cos") {
        *integrated = make_subtract(
                          make_divide(make_multiply(polynomial,
                                                    make_function("sin", argument)),
                                      SymbolicExpression::number(a)),
                          make_divide(recursive, SymbolicExpression::number(a)))
                          .simplify();
        return true;
    }
    return false;
}

bool decompose_power_factor(const SymbolicExpression& expression,
                            SymbolicExpression* base,
                            double* exponent) {
    if (expression.node_->type == NodeType::kPower &&
        SymbolicExpression(expression.node_->right).is_number(exponent)) {
        *base = SymbolicExpression(expression.node_->left).simplify();
        return true;
    }

    *base = expression;
    *exponent = 1.0;
    return true;
}

SymbolicExpression rebuild_power_difference(const SymbolicExpression& base, double exponent) {
    if (mymath::is_near_zero(exponent, kFormatEps)) {
        return SymbolicExpression::number(1.0);
    }
    if (mymath::is_near_zero(exponent - 1.0, kFormatEps)) {
        return base;
    }
    if (mymath::is_near_zero(exponent + 1.0, kFormatEps)) {
        return make_divide(SymbolicExpression::number(1.0), base).simplify();
    }
    if (exponent < 0.0) {
        return make_divide(SymbolicExpression::number(1.0),
                           make_power(base, SymbolicExpression::number(-exponent)))
            .simplify();
    }
    return make_power(base, SymbolicExpression::number(exponent)).simplify();
}

double common_numeric_factor(double lhs, double rhs) {
    const double lhs_abs = mymath::abs(lhs);
    const double rhs_abs = mymath::abs(rhs);
    if (mymath::is_near_zero(lhs_abs, kFormatEps) ||
        mymath::is_near_zero(rhs_abs, kFormatEps)) {
        return 0.0;
    }
    if (mymath::is_integer(lhs_abs, 1e-10) && mymath::is_integer(rhs_abs, 1e-10)) {
        long long a = static_cast<long long>(lhs_abs + 0.5);
        long long b = static_cast<long long>(rhs_abs + 0.5);
        while (b != 0) {
            const long long next = a % b;
            a = b;
            b = next;
        }
        return static_cast<double>(a);
    }
    if (mymath::is_near_zero(lhs_abs - rhs_abs, kFormatEps)) {
        return lhs_abs;
    }
    return 1.0;
}

bool try_factor_common_terms(const SymbolicExpression& left,
                             const SymbolicExpression& right,
                             double right_sign,
                             SymbolicExpression* combined) {
    double left_coefficient = 1.0;
    double right_coefficient = 1.0;
    std::vector<SymbolicExpression> left_factors;
    std::vector<SymbolicExpression> right_factors;
    collect_multiplicative_terms(left, &left_coefficient, &left_factors);
    collect_multiplicative_terms(right, &right_coefficient, &right_factors);

    std::vector<bool> right_used(right_factors.size(), false);
    std::vector<SymbolicExpression> common_factors;
    std::vector<SymbolicExpression> left_remaining;
    for (const SymbolicExpression& left_factor : left_factors) {
        bool matched = false;
        for (std::size_t i = 0; i < right_factors.size(); ++i) {
            if (right_used[i] || !expressions_match(left_factor, right_factors[i])) {
                continue;
            }
            right_used[i] = true;
            common_factors.push_back(left_factor);
            matched = true;
            break;
        }
        if (!matched) {
            left_remaining.push_back(left_factor);
        }
    }

    std::vector<SymbolicExpression> right_remaining;
    for (std::size_t i = 0; i < right_factors.size(); ++i) {
        if (!right_used[i]) {
            right_remaining.push_back(right_factors[i]);
        }
    }

    const double numeric_factor = common_numeric_factor(left_coefficient, right_coefficient);
    if (common_factors.empty() &&
        mymath::is_near_zero(numeric_factor - 1.0, kFormatEps)) {
        return false;
    }

    if (mymath::is_near_zero(numeric_factor, kFormatEps)) {
        return false;
    }

    const SymbolicExpression outer =
        rebuild_product_expression(numeric_factor, common_factors);
    const SymbolicExpression left_inner =
        rebuild_product_expression(left_coefficient / numeric_factor, left_remaining);
    const SymbolicExpression right_inner =
        rebuild_product_expression(right_coefficient / numeric_factor, right_remaining);
    const SymbolicExpression inner =
        right_sign > 0.0 ? make_add(left_inner, right_inner).simplify()
                         : make_subtract(left_inner, right_inner).simplify();

    if (expr_is_one(outer)) {
        *combined = inner;
        return true;
    }

    *combined = make_multiply(outer, inner).simplify();
    return true;
}

bool is_squared_function(const SymbolicExpression& expression,
                         const std::string& function_name,
                         std::string* argument_key) {
    if (expression.node_->type != NodeType::kPower) {
        return false;
    }

    double exponent = 0.0;
    if (!SymbolicExpression(expression.node_->right).is_number(&exponent) ||
        !mymath::is_near_zero(exponent - 2.0, kFormatEps)) {
        return false;
    }

    const SymbolicExpression base(expression.node_->left);
    if (base.node_->type != NodeType::kFunction || base.node_->text != function_name) {
        return false;
    }

    *argument_key = node_structural_key(base.node_->left);
    return true;
}

bool is_identifier_variable_name(const std::string& name) {
    if (name.empty() ||
        !std::isalpha(static_cast<unsigned char>(name.front()))) {
        return false;
    }

    for (char ch : name) {
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') {
            return false;
        }
    }
    return true;
}

void collect_identifier_variables(const SymbolicExpression& expression,
                                  std::vector<std::string>* names) {
    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kVariable:
            if (node->text != "pi" && node->text != "e" && node->text != "i" &&
                is_identifier_variable_name(node->text)) {
                names->push_back(node->text);
            }
            return;
        case NodeType::kNumber:
            return;
        case NodeType::kNegate:
        case NodeType::kFunction:
            collect_identifier_variables(SymbolicExpression(node->left), names);
            return;
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            collect_identifier_variables(SymbolicExpression(node->left), names);
            collect_identifier_variables(SymbolicExpression(node->right), names);
            return;
    }
}

std::string unique_identifier_variable(const SymbolicExpression& expression) {
    std::vector<std::string> names;
    collect_identifier_variables(expression, &names);
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names.size() == 1 ? names.front() : "";
}

bool polynomial_coefficients_from_simplified(const SymbolicExpression& expression,
                                             const std::string& variable_name,
                                             std::vector<double>* coefficients) {
    double numeric = 0.0;
    if (expression.is_number(&numeric)) {
        *coefficients = {numeric};
        return true;
    }
    if (expression.is_variable_named(variable_name)) {
        *coefficients = {0.0, 1.0};
        return true;
    }
    if (expression.is_constant(variable_name)) {
        if (expression.is_number(&numeric)) {
            *coefficients = {numeric};
            return true;
        }
        return false;
    }

    const auto& node = expression.node_;
    if (node->type == NodeType::kNegate) {
        if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                     variable_name,
                                                     coefficients)) {
            return false;
        }
        for (double& value : *coefficients) {
            value = -value;
        }
        trim_polynomial_coefficients(coefficients);
        return true;
    }

    std::vector<double> left;
    std::vector<double> right;
    switch (node->type) {
        case NodeType::kAdd:
            return polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                           variable_name,
                                                           &left) &&
                   polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                           variable_name,
                                                           &right) &&
                   ((*coefficients = polynomial_add_impl(left, right)), true);
        case NodeType::kSubtract:
            return polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                           variable_name,
                                                           &left) &&
                   polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                           variable_name,
                                                           &right) &&
                   ((*coefficients = polynomial_subtract_impl(left, right)), true);
        case NodeType::kMultiply:
            return polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                           variable_name,
                                                           &left) &&
                   polynomial_coefficients_from_simplified(SymbolicExpression(node->right),
                                                           variable_name,
                                                           &right) &&
                   ((*coefficients = polynomial_multiply_impl(left, right)), true);
        case NodeType::kDivide: {
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left)) {
                return false;
            }
            double divisor = 0.0;
            if (!SymbolicExpression(node->right).is_number(&divisor) ||
                mymath::is_near_zero(divisor, kFormatEps)) {
                return false;
            }
            for (double& value : left) {
                value /= divisor;
            }
            trim_polynomial_coefficients(&left);
            *coefficients = left;
            return true;
        }
        case NodeType::kPower: {
            double exponent = 0.0;
            if (!SymbolicExpression(node->right).is_number(&exponent) ||
                !mymath::is_integer(exponent, 1e-10) || exponent < 0.0) {
                return false;
            }
            if (!polynomial_coefficients_from_simplified(SymbolicExpression(node->left),
                                                         variable_name,
                                                         &left)) {
                return false;
            }
            std::vector<double> result = {1.0};
            for (int i = 0; i < static_cast<int>(exponent + 0.5); ++i) {
                result = polynomial_multiply_impl(result, left);
            }
            *coefficients = result;
            return true;
        }
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kFunction:
        case NodeType::kNegate:
            return false;
    }
    return false;
}

bool polynomial_is_zero_remainder(const std::vector<double>& coefficients) {
    for (double coefficient : coefficients) {
        if (!mymath::is_near_zero(coefficient, kFormatEps)) {
            return false;
        }
    }
    return true;
}

bool try_reduce_polynomial_quotient(const SymbolicExpression& left,
                                    const SymbolicExpression& right,
                                    SymbolicExpression* reduced) {
    const std::string variable_name = unique_identifier_variable(make_add(left, right));
    if (variable_name.empty()) {
        return false;
    }

    std::vector<double> numerator;
    std::vector<double> denominator;
    if (!polynomial_coefficients_from_simplified(left, variable_name, &numerator) ||
        !polynomial_coefficients_from_simplified(right, variable_name, &denominator)) {
        return false;
    }
    trim_polynomial_coefficients(&denominator);
    if (denominator.size() <= 1) {
        return false;
    }

    const PolynomialDivisionResult division = polynomial_divide(numerator, denominator);
    if (!polynomial_is_zero_remainder(division.remainder)) {
        return false;
    }

    *reduced = build_polynomial_expression_from_coefficients(division.quotient,
                                                             variable_name);
    return true;
}

bool is_single_variable_polynomial(const SymbolicExpression& expression) {
    const std::string variable_name = unique_identifier_variable(expression);
    if (variable_name.empty()) {
        return false;
    }

    std::vector<double> coefficients;
    return polynomial_coefficients_from_simplified(expression, variable_name, &coefficients);
}

SymbolicExpression maybe_canonicalize_polynomial(const SymbolicExpression& expression) {
    const std::string variable_name = unique_identifier_variable(expression);
    if (variable_name.empty()) {
        return expression;
    }

    std::vector<double> coefficients;
    if (!polynomial_coefficients_from_simplified(expression, variable_name, &coefficients)) {
        return expression;
    }

    const std::string canonical = polynomial_to_string(coefficients, variable_name);
    if (canonical == to_string_impl(expression.node_, 0)) {
        return expression;
    }
    return build_polynomial_expression_from_coefficients(coefficients, variable_name);
}

SymbolicExpression simplify_impl(const SymbolicExpression& expression) {
    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
        case NodeType::kVariable:
            return expression;
        case NodeType::kFunction: {
            const SymbolicExpression argument = SymbolicExpression(node->left).simplify();
            double numeric = 0.0;
            if (argument.is_number(&numeric)) {
                if (node->text == "asin") {
                    return SymbolicExpression::number(mymath::asin(numeric));
                }
                if (node->text == "acos") {
                    return SymbolicExpression::number(mymath::acos(numeric));
                }
                if (node->text == "atan") {
                    return SymbolicExpression::number(mymath::atan(numeric));
                }
                if (node->text == "sin") {
                    return SymbolicExpression::number(mymath::sin(numeric));
                }
                if (node->text == "cos") {
                    return SymbolicExpression::number(mymath::cos(numeric));
                }
                if (node->text == "tan") {
                    return SymbolicExpression::number(mymath::tan(numeric));
                }
                if (node->text == "exp") {
                    return SymbolicExpression::number(mymath::exp(numeric));
                }
                if (node->text == "sinh") {
                    return SymbolicExpression::number(mymath::sinh(numeric));
                }
                if (node->text == "cosh") {
                    return SymbolicExpression::number(mymath::cosh(numeric));
                }
                if (node->text == "tanh") {
                    return SymbolicExpression::number(mymath::tanh(numeric));
                }
                if (node->text == "ln") {
                    return SymbolicExpression::number(mymath::ln(numeric));
                }
                if (node->text == "sqrt") {
                    return SymbolicExpression::number(mymath::sqrt(numeric));
                }
                if (node->text == "abs") {
                    return SymbolicExpression::number(mymath::abs(numeric));
                }
                if (node->text == "floor") {
                    return SymbolicExpression::number(static_cast<double>(
                        static_cast<long long>(numeric < 0.0 && static_cast<double>(static_cast<long long>(numeric)) != numeric
                                                   ? numeric - 1.0
                                                   : numeric)));
                }
                if (node->text == "ceil") {
                    long long truncated = static_cast<long long>(numeric);
                    if (numeric > 0.0 && static_cast<double>(truncated) != numeric) {
                        ++truncated;
                    }
                    return SymbolicExpression::number(static_cast<double>(truncated));
                }
                if (node->text == "cbrt") {
                    return SymbolicExpression::number(mymath::cbrt(numeric));
                }
                if (node->text == "sign") {
                    if (mymath::is_near_zero(numeric, kFormatEps)) {
                        return SymbolicExpression::number(0.0);
                    }
                    return SymbolicExpression::number(numeric > 0.0 ? 1.0 : -1.0);
                }
            }

            if (node->text == "ln" && expr_is_variable(argument, "e")) {
                return SymbolicExpression::number(1.0);
            }
            if (node->text == "exp" &&
                argument.node_->type == NodeType::kFunction &&
                argument.node_->text == "ln" &&
                is_known_positive_expression(
                    SymbolicExpression(argument.node_->left).simplify())) {
                return SymbolicExpression(argument.node_->left).simplify();
            }
            if (node->text == "ln" &&
                argument.node_->type == NodeType::kFunction &&
                argument.node_->text == "exp") {
                return SymbolicExpression(argument.node_->left).simplify();
            }

            double pi_multiple = 0.0;
            if (decompose_numeric_multiple_of_symbol(argument, "pi", &pi_multiple)) {
                if (node->text == "sin") {
                    if (numeric_matches_any(pi_multiple, {0.0, 1.0, -1.0})) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (numeric_matches_any(pi_multiple, {0.5})) {
                        return SymbolicExpression::number(1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {-0.5})) {
                        return SymbolicExpression::number(-1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, 5.0 / 6.0})) {
                        return half_symbol();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 6.0, -5.0 / 6.0})) {
                        return make_negate(half_symbol()).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, 2.0 / 3.0})) {
                        return make_divide(sqrt3_symbol(),
                                           SymbolicExpression::number(2.0)).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 3.0, -2.0 / 3.0})) {
                        return make_negate(
                            make_divide(sqrt3_symbol(),
                                        SymbolicExpression::number(2.0))).simplify();
                    }
                }
                if (node->text == "cos") {
                    if (numeric_matches_any(pi_multiple, {0.0})) {
                        return SymbolicExpression::number(1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0, -1.0})) {
                        return SymbolicExpression::number(-1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {0.5, -0.5})) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, -1.0 / 3.0})) {
                        return half_symbol();
                    }
                    if (numeric_matches_any(pi_multiple, {2.0 / 3.0, -2.0 / 3.0})) {
                        return make_negate(half_symbol()).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, -1.0 / 6.0})) {
                        return make_divide(sqrt3_symbol(),
                                           SymbolicExpression::number(2.0)).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {5.0 / 6.0, -5.0 / 6.0})) {
                        return make_negate(
                            make_divide(sqrt3_symbol(),
                                        SymbolicExpression::number(2.0))).simplify();
                    }
                }
                if (node->text == "tan") {
                    if (numeric_matches_any(pi_multiple, {0.0, 1.0, -1.0})) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (numeric_matches_any(pi_multiple, {0.25})) {
                        return SymbolicExpression::number(1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {-0.25})) {
                        return SymbolicExpression::number(-1.0);
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 6.0, -5.0 / 6.0})) {
                        return make_divide(SymbolicExpression::number(1.0),
                                           sqrt3_symbol()).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 6.0, 5.0 / 6.0})) {
                        return make_negate(
                            make_divide(SymbolicExpression::number(1.0),
                                        sqrt3_symbol())).simplify();
                    }
                    if (numeric_matches_any(pi_multiple, {1.0 / 3.0, -2.0 / 3.0})) {
                        return sqrt3_symbol();
                    }
                    if (numeric_matches_any(pi_multiple, {-1.0 / 3.0, 2.0 / 3.0})) {
                        return make_negate(sqrt3_symbol()).simplify();
                    }
                }
            }

            return make_function(node->text, argument);
        }
        case NodeType::kNegate: {
            const SymbolicExpression operand = SymbolicExpression(node->left).simplify();
            double value = 0.0;
            if (operand.is_number(&value)) {
                return SymbolicExpression::number(-value);
            }
            if (operand.node_->type == NodeType::kNegate) {
                return SymbolicExpression(operand.node_->left).simplify();
            }
            return make_negate(operand);
        }
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower:
            break;
    }

    const SymbolicExpression left = SymbolicExpression(node->left).simplify();
    const SymbolicExpression right = SymbolicExpression(node->right).simplify();
    double left_value = 0.0;
    double right_value = 0.0;

    switch (node->type) {
        case NodeType::kAdd:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value + right_value);
            }
            if (expr_is_zero(left)) {
                return right;
            }
            if (expr_is_zero(right)) {
                return left;
            }
            {
                std::string left_argument;
                std::string right_argument;
                if ((is_squared_function(left, "sin", &left_argument) &&
                     is_squared_function(right, "cos", &right_argument) &&
                     left_argument == right_argument) ||
                    (is_squared_function(left, "cos", &left_argument) &&
                     is_squared_function(right, "sin", &right_argument) &&
                     left_argument == right_argument)) {
                    return SymbolicExpression::number(1.0);
                }
            }
            {
                SymbolicExpression combined;
                if (try_combine_like_terms(left, right, 1.0, &combined)) {
                    return combined;
                }
            }
            {
                const SymbolicExpression sum = make_add(left, right);
                if (is_single_variable_polynomial(sum)) {
                    return maybe_canonicalize_polynomial(sum);
                }
            }
            {
                SymbolicExpression factored;
                if (try_factor_common_terms(left, right, 1.0, &factored)) {
                    return factored;
                }
            }
            return make_add(left, right);
        case NodeType::kSubtract:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value - right_value);
            }
            if (expr_is_zero(right)) {
                return left;
            }
            if (expr_is_zero(left)) {
                return make_negate(right).simplify();
            }
            {
                SymbolicExpression combined;
                if (try_combine_like_terms(left, right, -1.0, &combined)) {
                    return combined;
                }
            }
            {
                const SymbolicExpression difference = make_subtract(left, right);
                if (is_single_variable_polynomial(difference)) {
                    return maybe_canonicalize_polynomial(difference);
                }
            }
            {
                SymbolicExpression factored;
                if (try_factor_common_terms(left, right, -1.0, &factored)) {
                    return factored;
                }
            }
            return make_subtract(left, right);
        case NodeType::kMultiply:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value * right_value);
            }
            if (expr_is_zero(left) || expr_is_zero(right)) {
                return SymbolicExpression::number(0.0);
            }
            if (expr_is_one(left)) {
                return right;
            }
            if (expr_is_one(right)) {
                return left;
            }
            if (expr_is_minus_one(left)) {
                return make_negate(right).simplify();
            }
            if (expr_is_minus_one(right)) {
                return make_negate(left).simplify();
            }
            {
                SymbolicExpression left_base;
                SymbolicExpression right_base;
                double left_exponent = 0.0;
                double right_exponent = 0.0;
                decompose_power_factor(left, &left_base, &left_exponent);
                decompose_power_factor(right, &right_base, &right_exponent);
                if (expressions_match(left_base, right_base)) {
                    return rebuild_power_difference(left_base,
                                                    left_exponent + right_exponent);
                }
            }
            {
                double numeric_factor = 1.0;
                std::vector<SymbolicExpression> symbolic_factors;
                collect_multiplicative_terms(left, &numeric_factor, &symbolic_factors);
                collect_multiplicative_terms(right, &numeric_factor, &symbolic_factors);

                if (mymath::is_near_zero(numeric_factor, kFormatEps)) {
                    return SymbolicExpression::number(0.0);
                }

                SymbolicExpression combined;
                bool has_combined = false;
                for (const SymbolicExpression& factor : symbolic_factors) {
                    if (!has_combined) {
                        combined = factor;
                        has_combined = true;
                    } else {
                        combined = make_multiply(combined, factor);
                    }
                }

                if (!has_combined) {
                    return SymbolicExpression::number(numeric_factor);
                }
                if (mymath::is_near_zero(numeric_factor - 1.0, kFormatEps)) {
                    return combined;
                }
                if (mymath::is_near_zero(numeric_factor + 1.0, kFormatEps)) {
                    return make_negate(combined).simplify();
                }
                return make_multiply(SymbolicExpression::number(numeric_factor), combined);
            }
        case NodeType::kDivide:
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(left_value / right_value);
            }
            if (expr_is_zero(left)) {
                return SymbolicExpression::number(0.0);
            }
            if (expr_is_one(right)) {
                return left;
            }
            {
                SymbolicExpression reduced;
                if (try_reduce_polynomial_quotient(left, right, &reduced)) {
                    return reduced;
                }
            }
            {
                SymbolicExpression left_base;
                SymbolicExpression right_base;
                double left_exponent = 0.0;
                double right_exponent = 0.0;
                decompose_power_factor(left, &left_base, &left_exponent);
                decompose_power_factor(right, &right_base, &right_exponent);
                if (expressions_match(left_base, right_base)) {
                    return rebuild_power_difference(left_base,
                                                    left_exponent - right_exponent);
                }
            }
            {
                double numerator_coefficient = 1.0;
                double denominator_coefficient = 1.0;
                std::vector<SymbolicExpression> numerator_factors;
                std::vector<SymbolicExpression> denominator_factors;
                collect_division_factors(left, &numerator_coefficient, &numerator_factors);
                collect_division_factors(right, &denominator_coefficient, &denominator_factors);
                if (!mymath::is_near_zero(denominator_coefficient, kFormatEps)) {
                    std::vector<bool> denominator_used(denominator_factors.size(), false);
                    std::vector<SymbolicExpression> reduced_numerator_factors;
                    for (const SymbolicExpression& numerator_factor : numerator_factors) {
                        const std::string numerator_key =
                            node_structural_key(numerator_factor.node_);
                        bool canceled = false;
                        for (std::size_t i = 0; i < denominator_factors.size(); ++i) {
                            if (denominator_used[i]) {
                                continue;
                            }
                            if (numerator_key ==
                                node_structural_key(denominator_factors[i].node_)) {
                                denominator_used[i] = true;
                                canceled = true;
                                break;
                            }
                        }
                        if (!canceled) {
                            reduced_numerator_factors.push_back(numerator_factor);
                        }
                    }

                    std::vector<SymbolicExpression> reduced_denominator_factors;
                    for (std::size_t i = 0; i < denominator_factors.size(); ++i) {
                        if (!denominator_used[i]) {
                            reduced_denominator_factors.push_back(denominator_factors[i]);
                        }
                    }

                    const bool symbolic_cancellation_happened =
                        reduced_numerator_factors.size() != numerator_factors.size() ||
                        reduced_denominator_factors.size() != denominator_factors.size();

                    const double reduced_coefficient =
                        numerator_coefficient / denominator_coefficient;
                    if (!symbolic_cancellation_happened &&
                        reduced_denominator_factors.empty() &&
                        mymath::is_near_zero(numerator_coefficient - 1.0, kFormatEps) &&
                        !mymath::is_near_zero(denominator_coefficient - 1.0, kFormatEps)) {
                        return make_divide(left,
                                           SymbolicExpression::number(denominator_coefficient));
                    }
                    SymbolicExpression numerator_expression =
                        rebuild_product_expression(reduced_coefficient,
                                                   reduced_numerator_factors);
                    SymbolicExpression denominator_expression =
                        rebuild_product_expression(1.0,
                                                   reduced_denominator_factors);

                    if (expr_is_zero(numerator_expression)) {
                        return SymbolicExpression::number(0.0);
                    }
                    if (expr_is_one(denominator_expression)) {
                        return numerator_expression;
                    }
                    if (numerator_expression.is_number(&left_value) &&
                        denominator_expression.is_number(&right_value)) {
                        return SymbolicExpression::number(left_value / right_value);
                    }
                    return make_divide(numerator_expression, denominator_expression);
                }
            }
            return make_divide(left, right);
        case NodeType::kPower:
            if (right.is_number(&right_value)) {
                if (mymath::is_near_zero(right_value, kFormatEps)) {
                    return SymbolicExpression::number(1.0);
                }
                if (mymath::is_near_zero(right_value - 1.0, kFormatEps)) {
                    return left;
                }
            }
            if (left.node_->type == NodeType::kPower) {
                double inner_exponent = 0.0;
                if (SymbolicExpression(left.node_->right).is_number(&inner_exponent) &&
                    right.is_number(&right_value)) {
                    return make_power(SymbolicExpression(left.node_->left).simplify(),
                                      SymbolicExpression::number(inner_exponent * right_value))
                        .simplify();
                }
            }
            if (left.is_number(&left_value) && right.is_number(&right_value)) {
                return SymbolicExpression::number(mymath::pow(left_value, right_value));
            }
            return make_power(left, right);
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kNegate:
        case NodeType::kFunction:
            break;
    }
    return expression;
}

double factorial_double(int exponent) {
    double value = 1.0;
    for (int i = 2; i <= exponent; ++i) {
        value *= static_cast<double>(i);
    }
    return value;
}

bool is_function_named(const SymbolicExpression& expression,
                       const std::string& function_name,
                       SymbolicExpression* argument = nullptr) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kFunction ||
        simplified.node_->text != function_name) {
        return false;
    }
    if (argument != nullptr) {
        *argument = SymbolicExpression(simplified.node_->left).simplify();
    }
    return true;
}

bool is_i_variable(const SymbolicExpression& expression) {
    return expr_is_variable(expression.simplify(), "i");
}

bool decompose_i_times_variable(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                double* coefficient) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric_factor = 1.0;
    std::vector<SymbolicExpression> symbolic_factors;
    collect_multiplicative_terms(simplified, &numeric_factor, &symbolic_factors);

    bool has_i = false;
    bool has_variable = false;
    for (const SymbolicExpression& factor : symbolic_factors) {
        if (is_i_variable(factor)) {
            if (has_i) {
                return false;
            }
            has_i = true;
            continue;
        }
        if (expr_is_variable(factor.simplify(), variable_name)) {
            if (has_variable) {
                return false;
            }
            has_variable = true;
            continue;
        }
        return false;
    }

    if (!has_i || !has_variable) {
        return false;
    }

    *coefficient = numeric_factor;
    return true;
}

bool decompose_unit_slope_shift(const SymbolicExpression& expression,
                                const std::string& variable_name,
                                double* shift) {
    double coefficient = 0.0;
    double intercept = 0.0;
    if (!decompose_linear(expression.simplify(),
                          variable_name,
                          &coefficient,
                          &intercept) ||
        !mymath::is_near_zero(coefficient - 1.0, kFormatEps)) {
        return false;
    }
    *shift = -intercept;
    return true;
}

bool match_step_shift(const SymbolicExpression& expression,
                      const std::string& variable_name,
                      double* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "step", &argument) &&
           decompose_unit_slope_shift(argument, variable_name, shift);
}

bool match_delta_shift(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       double* shift) {
    SymbolicExpression argument;
    return is_function_named(expression, "delta", &argument) &&
           decompose_unit_slope_shift(argument, variable_name, shift);
}

bool match_exponential_linear(const SymbolicExpression& expression,
                              const std::string& variable_name,
                              double* coefficient,
                              double* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "exp", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_sine_linear(const SymbolicExpression& expression,
                       const std::string& variable_name,
                       double* coefficient,
                       double* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "sin", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_cosine_linear(const SymbolicExpression& expression,
                         const std::string& variable_name,
                         double* coefficient,
                         double* intercept) {
    SymbolicExpression argument;
    return is_function_named(expression, "cos", &argument) &&
           decompose_linear(argument, variable_name, coefficient, intercept);
}

bool match_constant_power_sequence(const SymbolicExpression& expression,
                                   const std::string& index_variable,
                                   double* base) {
    const SymbolicExpression simplified = expression.simplify();
    if (simplified.node_->type != NodeType::kPower ||
        !SymbolicExpression(simplified.node_->left).is_number(base) ||
        !SymbolicExpression(simplified.node_->right).is_variable_named(index_variable)) {
        return false;
    }
    return true;
}

bool match_non_negative_integer_shift(const SymbolicExpression& expression,
                                      const std::string& variable_name,
                                      int* shift) {
    double raw_shift = 0.0;
    if (!decompose_unit_slope_shift(expression, variable_name, &raw_shift) ||
        !mymath::is_integer(raw_shift, 1e-10) ||
        raw_shift < 0.0) {
        return false;
    }
    *shift = static_cast<int>(raw_shift + 0.5);
    return true;
}

SymbolicExpression make_shifted_variable_expression(const std::string& variable_name,
                                                    double shift) {
    if (mymath::is_near_zero(shift, kFormatEps)) {
        return SymbolicExpression::variable(variable_name);
    }
    if (shift > 0.0) {
        return make_subtract(SymbolicExpression::variable(variable_name),
                             SymbolicExpression::number(shift)).simplify();
    }
    return make_add(SymbolicExpression::variable(variable_name),
                    SymbolicExpression::number(-shift)).simplify();
}

SymbolicExpression make_step_expression(const std::string& variable_name, double shift) {
    return make_function("step",
                         make_shifted_variable_expression(variable_name, shift)).simplify();
}

SymbolicExpression make_delta_expression(const std::string& variable_name, double shift) {
    return make_function("delta",
                         make_shifted_variable_expression(variable_name, shift)).simplify();
}

SymbolicExpression make_complex_phase(double coefficient,
                                      const std::string& variable_name) {
    SymbolicExpression phase =
        make_multiply(SymbolicExpression::variable("i"),
                      SymbolicExpression::variable(variable_name));
    if (!mymath::is_near_zero(coefficient - 1.0, kFormatEps)) {
        phase = make_multiply(SymbolicExpression::number(coefficient), phase);
    }
    return make_function("exp", phase).simplify();
}

SymbolicExpression make_z_shift_term(const std::string& z_variable, int shift) {
    if (shift == 0) {
        return SymbolicExpression::number(1.0);
    }
    return make_power(SymbolicExpression::variable(z_variable),
                      SymbolicExpression::number(-static_cast<double>(shift)))
        .simplify();
}

bool match_i_frequency_minus_constant(const SymbolicExpression& expression,
                                      const std::string& frequency_variable,
                                      double* constant) {
    const SymbolicExpression simplified = expression.simplify();
    const auto& node = simplified.node_;
    if (node->type == NodeType::kSubtract) {
        double coefficient = 0.0;
        if (decompose_i_times_variable(SymbolicExpression(node->left),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - 1.0, kFormatEps) &&
            SymbolicExpression(node->right).is_number(constant)) {
            return true;
        }
    }
    if (node->type == NodeType::kAdd) {
        double coefficient = 0.0;
        double numeric = 0.0;
        if (decompose_i_times_variable(SymbolicExpression(node->left),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - 1.0, kFormatEps) &&
            SymbolicExpression(node->right).is_number(&numeric)) {
            *constant = -numeric;
            return true;
        }
        if (decompose_i_times_variable(SymbolicExpression(node->right),
                                       frequency_variable,
                                       &coefficient) &&
            mymath::is_near_zero(coefficient - 1.0, kFormatEps) &&
            SymbolicExpression(node->left).is_number(&numeric)) {
            *constant = -numeric;
            return true;
        }
    }
    return false;
}

SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable);
SymbolicExpression inverse_laplace_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& transform_variable,
                                                  const std::string& time_variable);
SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable);
SymbolicExpression inverse_fourier_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& frequency_variable,
                                                  const std::string& time_variable);
SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable);
SymbolicExpression inverse_z_transform_impl(const SymbolicExpression& expression,
                                            const std::string& transform_variable,
                                            const std::string& index_variable);

SymbolicExpression laplace_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& transform_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_divide(SymbolicExpression::number(numeric),
                           SymbolicExpression::variable(transform_variable))
            .simplify();
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              transform_variable),
                       laplace_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              transform_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              transform_variable),
                       laplace_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              transform_variable))
                .simplify();
        case NodeType::kMultiply: {
            double constant = 0.0;
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   time_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     laplace_transform_impl(rest,
                                                            time_variable,
                                                            transform_variable))
                    .simplify();
            }

            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);
            double shift = 0.0;
            if (match_step_shift(left, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
                return laplace_transform_impl(right, time_variable, transform_variable);
            }
            if (match_step_shift(right, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
                return laplace_transform_impl(left, time_variable, transform_variable);
            }
            break;
        }
        case NodeType::kPower: {
            if (simplified.is_variable_named(time_variable)) {
                return make_divide(SymbolicExpression::number(1.0),
                                   make_power(SymbolicExpression::variable(transform_variable),
                                              SymbolicExpression::number(2.0)))
                    .simplify();
            }

            const SymbolicExpression base(simplified.node_->left);
            double exponent = 0.0;
            if (base.is_variable_named(time_variable) &&
                SymbolicExpression(simplified.node_->right).is_number(&exponent) &&
                mymath::is_integer(exponent, 1e-10) &&
                exponent >= 0.0) {
                const int order = static_cast<int>(exponent + 0.5);
                return make_divide(
                           SymbolicExpression::number(factorial_double(order)),
                           make_power(SymbolicExpression::variable(transform_variable),
                                      SymbolicExpression::number(
                                          static_cast<double>(order + 1))))
                    .simplify();
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kVariable:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                          time_variable,
                                          transform_variable))
            .simplify();
    }

    if (simplified.is_variable_named(time_variable)) {
        return make_divide(SymbolicExpression::number(1.0),
                           make_power(SymbolicExpression::variable(transform_variable),
                                      SymbolicExpression::number(2.0)))
            .simplify();
    }

    double shift = 0.0;
    if (match_step_shift(simplified, time_variable, &shift)) {
        if (shift < -kFormatEps) {
            throw std::runtime_error("laplace(step(t + a)) is not supported for negative shifts");
        }
        SymbolicExpression result = make_divide(SymbolicExpression::number(1.0),
                                                SymbolicExpression::variable(transform_variable));
        if (!mymath::is_near_zero(shift, kFormatEps)) {
            result = make_multiply(
                         make_function("exp",
                                       make_negate(
                                           make_multiply(SymbolicExpression::number(shift),
                                                         SymbolicExpression::variable(
                                                             transform_variable)))),
                         result)
                         .simplify();
        }
        return result.simplify();
    }

    if (match_delta_shift(simplified, time_variable, &shift)) {
        if (shift < -kFormatEps) {
            throw std::runtime_error("laplace(delta(t + a)) is not supported for negative shifts");
        }
        if (mymath::is_near_zero(shift, kFormatEps)) {
            return SymbolicExpression::number(1.0);
        }
        return make_function(
                   "exp",
                   make_negate(make_multiply(SymbolicExpression::number(shift),
                                             SymbolicExpression::variable(
                                                 transform_variable))))
            .simplify();
    }

    double linear_coefficient = 0.0;
    double linear_intercept = 0.0;
    if (match_exponential_linear(simplified,
                                 time_variable,
                                 &linear_coefficient,
                                 &linear_intercept) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps)) {
        SymbolicExpression denominator =
            make_subtract(SymbolicExpression::variable(transform_variable),
                          SymbolicExpression::number(linear_coefficient))
                .simplify();
        SymbolicExpression result =
            make_divide(SymbolicExpression::number(1.0), denominator).simplify();
        if (!mymath::is_near_zero(linear_intercept, kFormatEps)) {
            result = make_multiply(SymbolicExpression::number(mymath::exp(linear_intercept)),
                                   result)
                         .simplify();
        }
        return result;
    }

    if (match_sine_linear(simplified,
                          time_variable,
                          &linear_coefficient,
                          &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps)) {
        return make_divide(
                   SymbolicExpression::number(linear_coefficient),
                   make_add(make_power(SymbolicExpression::variable(transform_variable),
                                       SymbolicExpression::number(2.0)),
                            SymbolicExpression::number(
                                linear_coefficient * linear_coefficient)))
            .simplify();
    }

    if (match_cosine_linear(simplified,
                            time_variable,
                            &linear_coefficient,
                            &linear_intercept) &&
        mymath::is_near_zero(linear_intercept, kFormatEps) &&
        !mymath::is_near_zero(linear_coefficient, kFormatEps)) {
        return make_divide(
                   SymbolicExpression::variable(transform_variable),
                   make_add(make_power(SymbolicExpression::variable(transform_variable),
                                       SymbolicExpression::number(2.0)),
                            SymbolicExpression::number(
                                linear_coefficient * linear_coefficient)))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic Laplace transform");
}

SymbolicExpression inverse_laplace_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& transform_variable,
                                                  const std::string& time_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, 0.0))
            .simplify();
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           transform_variable,
                           time_variable),
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           transform_variable,
                           time_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           transform_variable,
                           time_variable),
                       inverse_laplace_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           transform_variable,
                           time_variable))
                .simplify();
        case NodeType::kMultiply: {
            double constant = 0.0;
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   transform_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     inverse_laplace_transform_impl(rest,
                                                                    transform_variable,
                                                                    time_variable))
                    .simplify();
            }

            double shift = 0.0;
            double a = 0.0;
            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);
            if (match_exponential_linear(left, transform_variable, &a, &shift) &&
                mymath::is_near_zero(a + shift, kFormatEps) &&
                a < 0.0) {
                const double delay = -a;
                if (right.simplify().is_variable_named(transform_variable)) {
                    return make_step_expression(time_variable, delay);
                }
                const SymbolicExpression shifted =
                    inverse_laplace_transform_impl(right, transform_variable, time_variable);
                if (right.simplify().node_->type == NodeType::kDivide &&
                    SymbolicExpression(right.simplify().node_->left).is_number(&numeric) &&
                    mymath::is_near_zero(numeric - 1.0, kFormatEps)) {
                    double pole = 0.0;
                    if (SymbolicExpression(right.simplify().node_->right)
                            .simplify()
                            .node_->type == NodeType::kSubtract &&
                        SymbolicExpression(
                            SymbolicExpression(right.simplify().node_->right)
                                .simplify()
                                .node_->left)
                            .is_variable_named(transform_variable) &&
                        SymbolicExpression(
                            SymbolicExpression(right.simplify().node_->right)
                                .simplify()
                                .node_->right)
                            .is_number(&pole)) {
                        return make_multiply(
                                   make_function(
                                       "exp",
                                       make_multiply(SymbolicExpression::number(pole),
                                                     make_shifted_variable_expression(
                                                         time_variable,
                                                         delay))),
                                   make_step_expression(time_variable, delay))
                            .simplify();
                    }
                }
                return shifted;
            }
            break;
        }
        case NodeType::kDivide:
        case NodeType::kFunction:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   inverse_laplace_transform_impl(SymbolicExpression(simplified.node_->left),
                                                  transform_variable,
                                                  time_variable))
            .simplify();
    }

    if (simplified.is_variable_named(transform_variable)) {
        return make_delta_expression(time_variable, 0.0);
    }

    if (simplified.node_->type == NodeType::kDivide) {
        const SymbolicExpression numerator(simplified.node_->left);
        const SymbolicExpression denominator = SymbolicExpression(simplified.node_->right).simplify();

        if (numerator.is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0, kFormatEps)) {
            if (denominator.is_variable_named(transform_variable)) {
                return make_step_expression(time_variable, 0.0);
            }

            if (denominator.node_->type == NodeType::kPower &&
                SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable)) {
                double exponent = 0.0;
                if (SymbolicExpression(denominator.node_->right).is_number(&exponent) &&
                    mymath::is_integer(exponent, 1e-10) &&
                    exponent >= 1.0) {
                    const int order = static_cast<int>(exponent + 0.5) - 1;
                    SymbolicExpression result;
                    if (order == 0) {
                        result = SymbolicExpression::number(1.0);
                    } else {
                        result = make_divide(
                                     make_power(SymbolicExpression::variable(time_variable),
                                                SymbolicExpression::number(
                                                    static_cast<double>(order))),
                                     SymbolicExpression::number(
                                         factorial_double(order)))
                                     .simplify();
                    }
                    return make_multiply(result,
                                         make_step_expression(time_variable, 0.0))
                        .simplify();
                }
            }

            if (denominator.node_->type == NodeType::kSubtract &&
                SymbolicExpression(denominator.node_->left)
                    .is_variable_named(transform_variable) &&
                SymbolicExpression(denominator.node_->right).is_number(&numeric)) {
                return make_multiply(
                           make_function(
                               "exp",
                               make_multiply(SymbolicExpression::number(numeric),
                                             SymbolicExpression::variable(time_variable))),
                           make_step_expression(time_variable, 0.0))
                    .simplify();
            }
        }

        double sine_frequency = 0.0;
        if (denominator.node_->type == NodeType::kAdd &&
            SymbolicExpression(denominator.node_->left).node_->type == NodeType::kPower &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->left)
                .is_variable_named(transform_variable) &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->right)
                .is_number(&numeric) &&
            mymath::is_near_zero(numeric - 2.0, kFormatEps) &&
            SymbolicExpression(denominator.node_->right).is_number(&sine_frequency) &&
            sine_frequency > 0.0) {
            const double frequency = mymath::sqrt(sine_frequency);
            if (numerator.is_variable_named(transform_variable)) {
                return make_multiply(make_function("cos",
                                                   make_multiply(
                                                       SymbolicExpression::number(frequency),
                                                       SymbolicExpression::variable(
                                                           time_variable))),
                                     make_step_expression(time_variable, 0.0))
                    .simplify();
            }
            if (numerator.is_number(&numeric) &&
                mymath::is_near_zero(numeric - frequency, 1e-10)) {
                return make_multiply(make_function("sin",
                                                   make_multiply(
                                                       SymbolicExpression::number(frequency),
                                                       SymbolicExpression::variable(
                                                           time_variable))),
                                     make_step_expression(time_variable, 0.0))
                    .simplify();
            }
        }
    }

    double shift = 0.0;
    if (match_exponential_linear(simplified,
                                 transform_variable,
                                 &numeric,
                                 &shift) &&
        mymath::is_near_zero(shift, kFormatEps) &&
        numeric < 0.0) {
        return make_delta_expression(time_variable, -numeric);
    }

    throw std::runtime_error("unsupported symbolic inverse Laplace transform");
}

SymbolicExpression fourier_transform_impl(const SymbolicExpression& expression,
                                          const std::string& time_variable,
                                          const std::string& frequency_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_multiply(
                   SymbolicExpression::number(2.0 * mymath::kPi * numeric),
                   make_delta_expression(frequency_variable, 0.0))
            .simplify();
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              frequency_variable),
                       fourier_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              frequency_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                              time_variable,
                                              frequency_variable),
                       fourier_transform_impl(SymbolicExpression(simplified.node_->right),
                                              time_variable,
                                              frequency_variable))
                .simplify();
        case NodeType::kMultiply: {
            double constant = 0.0;
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   time_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     fourier_transform_impl(rest,
                                                            time_variable,
                                                            frequency_variable))
                    .simplify();
            }

            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);
            double shift = 0.0;
            if (match_step_shift(left, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
                double exponent = 0.0;
                double intercept = 0.0;
                if (match_exponential_linear(right,
                                             time_variable,
                                             &exponent,
                                             &intercept) &&
                    !mymath::is_near_zero(exponent, kFormatEps)) {
                    return make_divide(
                               SymbolicExpression::number(mymath::exp(intercept)),
                               make_subtract(
                                   make_multiply(SymbolicExpression::variable("i"),
                                                 SymbolicExpression::variable(
                                                     frequency_variable)),
                                   SymbolicExpression::number(exponent)))
                        .simplify();
                }
            }
            if (match_step_shift(right, time_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
                double exponent = 0.0;
                double intercept = 0.0;
                if (match_exponential_linear(left,
                                             time_variable,
                                             &exponent,
                                             &intercept) &&
                    !mymath::is_near_zero(exponent, kFormatEps)) {
                    return make_divide(
                               SymbolicExpression::number(mymath::exp(intercept)),
                               make_subtract(
                                   make_multiply(SymbolicExpression::variable("i"),
                                                 SymbolicExpression::variable(
                                                     frequency_variable)),
                                   SymbolicExpression::number(exponent)))
                        .simplify();
                }
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                          time_variable,
                                          frequency_variable))
            .simplify();
    }

    double shift = 0.0;
    if (match_delta_shift(simplified, time_variable, &shift)) {
        if (mymath::is_near_zero(shift, kFormatEps)) {
            return SymbolicExpression::number(1.0);
        }
        return make_complex_phase(-shift, frequency_variable);
    }

    if (match_step_shift(simplified, time_variable, &shift)) {
        if (!mymath::is_near_zero(shift, kFormatEps)) {
            return make_multiply(make_complex_phase(-shift, frequency_variable),
                                 fourier_transform_impl(make_step_expression(time_variable, 0.0),
                                                        time_variable,
                                                        frequency_variable))
                .simplify();
        }
        return make_add(
                   make_multiply(SymbolicExpression::number(mymath::kPi),
                                 make_delta_expression(frequency_variable, 0.0)),
                   make_divide(
                       SymbolicExpression::number(1.0),
                       make_multiply(SymbolicExpression::variable("i"),
                                     SymbolicExpression::variable(frequency_variable))))
            .simplify();
    }

    double coefficient = 0.0;
    double intercept = 0.0;
    if (match_cosine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps) &&
        !mymath::is_near_zero(coefficient, kFormatEps)) {
        return make_multiply(
                   SymbolicExpression::number(mymath::kPi),
                   make_add(make_delta_expression(frequency_variable, coefficient),
                            make_delta_expression(frequency_variable, -coefficient)))
            .simplify();
    }

    if (match_sine_linear(simplified, time_variable, &coefficient, &intercept) &&
        mymath::is_near_zero(intercept, kFormatEps) &&
        !mymath::is_near_zero(coefficient, kFormatEps)) {
        return make_divide(
                   make_multiply(
                       SymbolicExpression::number(mymath::kPi),
                       make_subtract(make_delta_expression(frequency_variable, coefficient),
                                     make_delta_expression(frequency_variable, -coefficient))),
                   SymbolicExpression::variable("i"))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic Fourier transform");
}

SymbolicExpression inverse_fourier_transform_impl(const SymbolicExpression& expression,
                                                  const std::string& frequency_variable,
                                                  const std::string& time_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(time_variable, 0.0))
            .simplify();
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           frequency_variable,
                           time_variable),
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           frequency_variable,
                           time_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->left),
                           frequency_variable,
                           time_variable),
                       inverse_fourier_transform_impl(
                           SymbolicExpression(simplified.node_->right),
                           frequency_variable,
                           time_variable))
                .simplify();
        case NodeType::kMultiply: {
            double constant = 0.0;
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   frequency_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     inverse_fourier_transform_impl(rest,
                                                                    frequency_variable,
                                                                    time_variable))
                    .simplify();
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   inverse_fourier_transform_impl(SymbolicExpression(simplified.node_->left),
                                                  frequency_variable,
                                                  time_variable))
            .simplify();
    }

    double shift = 0.0;
    if (match_delta_shift(simplified, frequency_variable, &shift)) {
        return make_multiply(
                   SymbolicExpression::number(1.0 / (2.0 * mymath::kPi)),
                   make_complex_phase(shift, time_variable))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kFunction &&
        simplified.node_->text == "exp") {
        double phase = 0.0;
        if (decompose_i_times_variable(SymbolicExpression(simplified.node_->left),
                                       frequency_variable,
                                       &phase)) {
            return make_delta_expression(time_variable, -phase);
        }
    }

    double constant = 0.0;
    if (simplified.node_->type == NodeType::kDivide &&
        SymbolicExpression(simplified.node_->left).is_number(&constant) &&
        match_i_frequency_minus_constant(SymbolicExpression(simplified.node_->right),
                                         frequency_variable,
                                         &shift)) {
        return make_multiply(
                   make_multiply(SymbolicExpression::number(constant),
                                 make_function(
                                     "exp",
                                     make_multiply(SymbolicExpression::number(shift),
                                                   SymbolicExpression::variable(
                                                       time_variable)))),
                   make_step_expression(time_variable, 0.0))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic inverse Fourier transform");
}

SymbolicExpression z_transform_impl(const SymbolicExpression& expression,
                                    const std::string& index_variable,
                                    const std::string& transform_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_divide(
                   make_multiply(SymbolicExpression::number(numeric),
                                 SymbolicExpression::variable(transform_variable)),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 SymbolicExpression::number(1.0)))
            .simplify();
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       z_transform_impl(SymbolicExpression(simplified.node_->left),
                                        index_variable,
                                        transform_variable),
                       z_transform_impl(SymbolicExpression(simplified.node_->right),
                                        index_variable,
                                        transform_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       z_transform_impl(SymbolicExpression(simplified.node_->left),
                                        index_variable,
                                        transform_variable),
                       z_transform_impl(SymbolicExpression(simplified.node_->right),
                                        index_variable,
                                        transform_variable))
                .simplify();
        case NodeType::kMultiply: {
            double constant = 0.0;
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   index_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     z_transform_impl(rest,
                                                      index_variable,
                                                      transform_variable))
                    .simplify();
            }

            const SymbolicExpression left(simplified.node_->left);
            const SymbolicExpression right(simplified.node_->right);
            double shift = 0.0;
            if (match_step_shift(left, index_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
                return z_transform_impl(right, index_variable, transform_variable);
            }
            if (match_step_shift(right, index_variable, &shift) &&
                mymath::is_near_zero(shift, kFormatEps)) {
                return z_transform_impl(left, index_variable, transform_variable);
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   z_transform_impl(SymbolicExpression(simplified.node_->left),
                                    index_variable,
                                    transform_variable))
            .simplify();
    }

    int shift = 0;
    SymbolicExpression argument;
    if (is_function_named(simplified, "delta", &argument) &&
        match_non_negative_integer_shift(argument, index_variable, &shift)) {
        return make_z_shift_term(transform_variable, shift);
    }

    if (is_function_named(simplified, "step", &argument) &&
        match_non_negative_integer_shift(argument, index_variable, &shift)) {
        SymbolicExpression numerator = SymbolicExpression::variable(transform_variable);
        if (shift != 0) {
            numerator = make_multiply(
                            make_power(SymbolicExpression::variable(transform_variable),
                                       SymbolicExpression::number(
                                           static_cast<double>(1 - shift))),
                            SymbolicExpression::number(1.0))
                            .simplify();
        }
        return make_divide(
                   shift == 0
                       ? numerator
                       : make_power(SymbolicExpression::variable(transform_variable),
                                    SymbolicExpression::number(
                                        static_cast<double>(1 - shift))),
                   make_subtract(SymbolicExpression::variable(transform_variable),
                                 SymbolicExpression::number(1.0)))
            .simplify();
    }

    if (simplified.is_variable_named(index_variable)) {
        return make_divide(SymbolicExpression::variable(transform_variable),
                           make_power(
                               make_subtract(SymbolicExpression::variable(transform_variable),
                                             SymbolicExpression::number(1.0)),
                               SymbolicExpression::number(2.0)))
            .simplify();
    }

    double base = 0.0;
    if (match_constant_power_sequence(simplified, index_variable, &base)) {
        return make_divide(SymbolicExpression::variable(transform_variable),
                           make_subtract(SymbolicExpression::variable(transform_variable),
                                         SymbolicExpression::number(base)))
            .simplify();
    }

    throw std::runtime_error("unsupported symbolic z transform");
}

SymbolicExpression inverse_z_transform_impl(const SymbolicExpression& expression,
                                            const std::string& transform_variable,
                                            const std::string& index_variable) {
    const SymbolicExpression simplified = expression.simplify();
    double numeric = 0.0;
    if (simplified.is_number(&numeric)) {
        return make_multiply(SymbolicExpression::number(numeric),
                             make_delta_expression(index_variable, 0.0))
            .simplify();
    }

    switch (simplified.node_->type) {
        case NodeType::kAdd:
            return make_add(
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->left),
                                                transform_variable,
                                                index_variable),
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->right),
                                                transform_variable,
                                                index_variable))
                .simplify();
        case NodeType::kSubtract:
            return make_subtract(
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->left),
                                                transform_variable,
                                                index_variable),
                       inverse_z_transform_impl(SymbolicExpression(simplified.node_->right),
                                                transform_variable,
                                                index_variable))
                .simplify();
        case NodeType::kMultiply: {
            double constant = 0.0;
            SymbolicExpression rest;
            if (decompose_constant_times_expression(simplified,
                                                   transform_variable,
                                                   &constant,
                                                   &rest)) {
                return make_multiply(SymbolicExpression::number(constant),
                                     inverse_z_transform_impl(rest,
                                                              transform_variable,
                                                              index_variable))
                    .simplify();
            }
            break;
        }
        case NodeType::kFunction:
        case NodeType::kDivide:
        case NodeType::kNegate:
        case NodeType::kNumber:
        case NodeType::kVariable:
        case NodeType::kPower:
            break;
    }

    if (simplified.node_->type == NodeType::kNegate) {
        return make_negate(
                   inverse_z_transform_impl(SymbolicExpression(simplified.node_->left),
                                            transform_variable,
                                            index_variable))
            .simplify();
    }

    if (simplified.node_->type == NodeType::kPower &&
        SymbolicExpression(simplified.node_->left).is_variable_named(transform_variable) &&
        SymbolicExpression(simplified.node_->right).is_number(&numeric) &&
        mymath::is_integer(numeric, 1e-10) &&
        numeric <= 0.0) {
        return make_delta_expression(index_variable, static_cast<int>(-numeric + 0.5));
    }

    if (simplified.node_->type == NodeType::kDivide) {
        const SymbolicExpression numerator = SymbolicExpression(simplified.node_->left).simplify();
        const SymbolicExpression denominator = SymbolicExpression(simplified.node_->right).simplify();

        if (numerator.is_variable_named(transform_variable) &&
            denominator.node_->type == NodeType::kSubtract &&
            SymbolicExpression(denominator.node_->left).is_variable_named(transform_variable) &&
            SymbolicExpression(denominator.node_->right).is_number(&numeric)) {
            return make_multiply(
                       make_power(SymbolicExpression::number(numeric),
                                  SymbolicExpression::variable(index_variable)),
                       make_step_expression(index_variable, 0.0))
                .simplify();
        }

        if (numerator.is_variable_named(transform_variable) &&
            denominator.node_->type == NodeType::kPower &&
            SymbolicExpression(denominator.node_->left).node_->type == NodeType::kSubtract &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->left)
                .is_variable_named(transform_variable) &&
            SymbolicExpression(
                SymbolicExpression(denominator.node_->left).node_->right)
                .is_number(&numeric) &&
            mymath::is_near_zero(numeric - 1.0, kFormatEps) &&
            SymbolicExpression(denominator.node_->right).is_number(&numeric) &&
            mymath::is_near_zero(numeric - 2.0, kFormatEps)) {
            return make_multiply(SymbolicExpression::variable(index_variable),
                                 make_step_expression(index_variable, 0.0))
                .simplify();
        }
    }

    throw std::runtime_error("unsupported symbolic inverse z transform");
}

}  // namespace

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
    return simplify_impl(*this);
}

SymbolicExpression SymbolicExpression::derivative(const std::string& variable_name) const {
    // 这里实现的是规则式符号求导。
    // 每种节点都直接对应一条微分规则：
    // - 和差法则
    // - 乘法法则
    // - 商法则
    // - 幂函数求导
    // - 链式法则（通过 inner derivative 接进去）
    switch (node_->type) {
        case NodeType::kNumber:
            return number(0.0);
        case NodeType::kVariable:
            return number(node_->text == variable_name ? 1.0 : 0.0);
        case NodeType::kNegate:
            return make_negate(SymbolicExpression(node_->left).derivative(variable_name)).simplify();
        case NodeType::kAdd:
            return make_add(SymbolicExpression(node_->left).derivative(variable_name),
                            SymbolicExpression(node_->right).derivative(variable_name)).simplify();
        case NodeType::kSubtract:
            return make_subtract(SymbolicExpression(node_->left).derivative(variable_name),
                                 SymbolicExpression(node_->right).derivative(variable_name)).simplify();
        case NodeType::kMultiply: {
            const SymbolicExpression left(node_->left);
            const SymbolicExpression right(node_->right);
            return make_add(
                       make_multiply(left.derivative(variable_name), right),
                       make_multiply(left, right.derivative(variable_name)))
                .simplify();
        }
        case NodeType::kDivide: {
            const SymbolicExpression left(node_->left);
            const SymbolicExpression right(node_->right);
            return make_divide(
                       make_subtract(make_multiply(left.derivative(variable_name), right),
                                     make_multiply(left, right.derivative(variable_name))),
                       make_power(right, number(2.0)))
                .simplify();
        }
        case NodeType::kPower: {
            const SymbolicExpression base(node_->left);
            const SymbolicExpression exponent(node_->right);
            double exponent_value = 0.0;
            if (exponent.is_number(&exponent_value)) {
                return make_multiply(
                           make_multiply(number(exponent_value),
                                         make_power(base, number(exponent_value - 1.0))),
                           base.derivative(variable_name))
                    .simplify();
            }
            if (base.is_constant(variable_name)) {
                return make_multiply(
                           make_multiply(*this, make_function("ln", base)),
                           exponent.derivative(variable_name))
                    .simplify();
            }
            return make_multiply(
                       *this,
                       make_add(
                           make_multiply(exponent.derivative(variable_name), make_function("ln", base)),
                           make_multiply(exponent,
                                         make_divide(base.derivative(variable_name), base))))
                .simplify();
        }
        case NodeType::kFunction: {
            const SymbolicExpression argument(node_->left);
            const SymbolicExpression inner = argument.derivative(variable_name);
            if (node_->text == "asin") {
                return make_divide(
                           inner,
                           make_function("sqrt",
                                         make_subtract(number(1.0),
                                                       make_power(argument, number(2.0)))))
                    .simplify();
            }
            if (node_->text == "acos") {
                return make_negate(
                           make_divide(
                               inner,
                               make_function("sqrt",
                                             make_subtract(number(1.0),
                                                           make_power(argument, number(2.0))))))
                    .simplify();
            }
            if (node_->text == "atan") {
                return make_divide(
                           inner,
                           make_add(number(1.0), make_power(argument, number(2.0))))
                    .simplify();
            }
            if (node_->text == "sin") {
                return make_multiply(make_function("cos", argument), inner).simplify();
            }
            if (node_->text == "cos") {
                return make_multiply(make_negate(make_function("sin", argument)), inner).simplify();
            }
            if (node_->text == "tan") {
                return make_multiply(make_divide(number(1.0),
                                                 make_power(make_function("cos", argument),
                                                            number(2.0))),
                                     inner)
                    .simplify();
            }
            if (node_->text == "exp") {
                return make_multiply(make_function("exp", argument), inner).simplify();
            }
            if (node_->text == "sinh") {
                return make_multiply(make_function("cosh", argument), inner).simplify();
            }
            if (node_->text == "cosh") {
                return make_multiply(make_function("sinh", argument), inner).simplify();
            }
            if (node_->text == "tanh") {
                return make_divide(inner,
                                   make_power(make_function("cosh", argument),
                                              number(2.0)))
                    .simplify();
            }
            if (node_->text == "ln") {
                return make_divide(inner, argument).simplify();
            }
            if (node_->text == "sqrt") {
                return make_divide(inner,
                                   make_multiply(number(2.0), make_function("sqrt", argument)))
                    .simplify();
            }
            if (node_->text == "cbrt") {
                return make_divide(inner,
                                   make_multiply(number(3.0),
                                                 make_power(make_function("cbrt", argument),
                                                            number(2.0))))
                    .simplify();
            }
            if (node_->text == "abs") {
                return make_multiply(make_function("sign", argument), inner).simplify();
            }
            if (node_->text == "step") {
                return make_multiply(make_function("delta", argument), inner).simplify();
            }
            throw std::runtime_error("symbolic derivative does not support function: " + node_->text);
        }
    }
    throw std::runtime_error("unsupported symbolic derivative");
}

SymbolicExpression SymbolicExpression::integral(const std::string& variable_name) const {
    // 这里的符号积分是“有限规则覆盖”，不是完整 CAS。
    // 当前主要支持：
    // - 常数、多项式
    // - 常数倍
    // - 和差
    // - sin(ax+b), cos(ax+b), exp(ax+b), 1/(ax+b)
    //
    // 不支持的结构会明确报错，而不是静默返回错误结果。
    double numeric_value = 0.0;
    if (is_constant(variable_name) && is_number(&numeric_value)) {
        return make_multiply(number(numeric_value), variable(variable_name)).simplify();
    }

    switch (node_->type) {
        case NodeType::kNumber:
            return make_multiply(number(node_->number_value), variable(variable_name)).simplify();
        case NodeType::kVariable:
            if (node_->text == variable_name) {
                return make_divide(make_power(variable(variable_name), number(2.0)),
                                   number(2.0))
                    .simplify();
            }
            return make_multiply(variable(node_->text), variable(variable_name)).simplify();
        case NodeType::kNegate:
            return make_negate(SymbolicExpression(node_->left).integral(variable_name)).simplify();
        case NodeType::kAdd:
            return make_add(SymbolicExpression(node_->left).integral(variable_name),
                            SymbolicExpression(node_->right).integral(variable_name)).simplify();
        case NodeType::kSubtract:
            return make_subtract(SymbolicExpression(node_->left).integral(variable_name),
                                 SymbolicExpression(node_->right).integral(variable_name)).simplify();
        case NodeType::kMultiply: {
            double constant = 0.0;
            SymbolicExpression rest;
            if (decompose_constant_times_expression(*this, variable_name, &constant, &rest)) {
                return make_multiply(number(constant), rest.integral(variable_name)).simplify();
            }
            const SymbolicExpression left(node_->left);
            const SymbolicExpression right(node_->right);
            SymbolicExpression polynomial;
            SymbolicExpression integrated;
            if (polynomial_expression(left, variable_name, &polynomial) &&
                right.node_->type == NodeType::kFunction &&
                integrate_polynomial_times_function(polynomial,
                                                    right.node_->text,
                                                    SymbolicExpression(right.node_->left),
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            if (polynomial_expression(right, variable_name, &polynomial) &&
                left.node_->type == NodeType::kFunction &&
                integrate_polynomial_times_function(polynomial,
                                                    left.node_->text,
                                                    SymbolicExpression(left.node_->left),
                                                    variable_name,
                                                    &integrated)) {
                return integrated.simplify();
            }
            throw std::runtime_error("symbolic integral does not support this product");
        }
        case NodeType::kPower:
        case NodeType::kFunction:
        case NodeType::kDivide:
            break;
    }

    if (node_->type == NodeType::kDivide) {
        const SymbolicExpression left(node_->left);
        const SymbolicExpression right(node_->right);
        if (left.is_number(&numeric_value) && mymath::is_near_zero(numeric_value - 1.0, kFormatEps)) {
            double a = 0.0;
            double b = 0.0;
            if (decompose_linear(right, variable_name, &a, &b) &&
                !mymath::is_near_zero(a, kFormatEps)) {
                return make_divide(make_function("ln", make_function("abs", right)),
                                   number(a))
                    .simplify();
            }
        }
        throw std::runtime_error("symbolic integral does not support this quotient");
    }

    if (node_->type == NodeType::kPower) {
        const SymbolicExpression base(node_->left);
        const SymbolicExpression exponent(node_->right);
        double exponent_value = 0.0;
        double a = 0.0;
        double b = 0.0;
        if (exponent.is_number(&exponent_value) &&
            decompose_linear(base, variable_name, &a, &b) &&
            !mymath::is_near_zero(a, kFormatEps)) {
            if (mymath::is_near_zero(exponent_value + 1.0, kFormatEps)) {
                return make_divide(make_function("ln", make_function("abs", base)),
                                   number(a))
                    .simplify();
            }
            return make_divide(make_power(base, number(exponent_value + 1.0)),
                               number(a * (exponent_value + 1.0)))
                .simplify();
        }
        throw std::runtime_error("symbolic integral only supports powers of the integration variable");
    }

    if (node_->type == NodeType::kFunction) {
        const SymbolicExpression argument(node_->left);
        double a = 0.0;
        double b = 0.0;
        const bool linear = decompose_linear(argument, variable_name, &a, &b) &&
                            !mymath::is_near_zero(a, kFormatEps);
        if (node_->text == "sin" && linear) {
            return make_divide(make_negate(make_function("cos", argument)),
                               number(a))
                .simplify();
        }
        if (node_->text == "cos" && linear) {
            return make_divide(make_function("sin", argument), number(a)).simplify();
        }
        if (node_->text == "exp" && linear) {
            return make_divide(make_function("exp", argument), number(a)).simplify();
        }
        if (node_->text == "sqrt" && linear) {
            return make_divide(make_multiply(number(2.0),
                                             make_power(make_function("sqrt", argument),
                                                        number(3.0))),
                               number(3.0 * a))
                .simplify();
        }
        if (node_->text == "cbrt" && linear) {
            return make_divide(make_multiply(number(3.0),
                                             make_power(make_function("cbrt", argument),
                                                        number(4.0))),
                               number(4.0 * a))
                .simplify();
        }
        if (node_->text == "tan" && linear) {
            return make_divide(make_negate(make_function("ln",
                                                         make_function("abs",
                                                                       make_function("cos",
                                                                                     argument)))),
                               number(a))
                .simplify();
        }
        if (node_->text == "delta" &&
            argument.is_variable_named(variable_name)) {
            return make_step_expression(variable_name, 0.0);
        }
        throw std::runtime_error("symbolic integral does not support function: " + node_->text);
    }

    throw std::runtime_error("unsupported symbolic integral");
}

SymbolicExpression SymbolicExpression::fourier_transform(
    const std::string& time_variable,
    const std::string& frequency_variable) const {
    return fourier_transform_impl(*this, time_variable, frequency_variable).simplify();
}

SymbolicExpression SymbolicExpression::inverse_fourier_transform(
    const std::string& frequency_variable,
    const std::string& time_variable) const {
    return inverse_fourier_transform_impl(*this, frequency_variable, time_variable).simplify();
}

SymbolicExpression SymbolicExpression::laplace_transform(
    const std::string& time_variable,
    const std::string& transform_variable) const {
    return laplace_transform_impl(*this, time_variable, transform_variable).simplify();
}

SymbolicExpression SymbolicExpression::inverse_laplace_transform(
    const std::string& transform_variable,
    const std::string& time_variable) const {
    return inverse_laplace_transform_impl(*this, transform_variable, time_variable).simplify();
}

SymbolicExpression SymbolicExpression::z_transform(
    const std::string& index_variable,
    const std::string& transform_variable) const {
    return z_transform_impl(*this, index_variable, transform_variable).simplify();
}

SymbolicExpression SymbolicExpression::inverse_z_transform(
    const std::string& transform_variable,
    const std::string& index_variable) const {
    return inverse_z_transform_impl(*this, transform_variable, index_variable).simplify();
}

SymbolicExpression SymbolicExpression::substitute(
    const std::string& variable_name,
    const SymbolicExpression& replacement) const {
    return substitute_impl(*this, variable_name, replacement).simplify();
}

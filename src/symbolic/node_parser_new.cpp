// ============================================================================
// 表达式解析器
// ============================================================================

#include "symbolic/symbolic_expression_internal.h"

#include "core/format_utils.h"
#include "math/mymath.h"
#include "math/mymath_dual.h"

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

#include <mutex>

#include "parser/base_parser.h"

namespace symbolic_expression_internal {

// ============================================================================
// 表达式解析器
// ============================================================================

/**
 * @class Parser
 * @brief 递归下降表达式解析器
 *
 * 解析数学表达式字符串，构建表达式树。
 *
 * 文法规则：
 * ```
 * expression -> term (('+' | '-') term)*
 * term       -> unary (('*' | '/') unary)*
 * unary      -> ('+' | '-') unary | power
 * power      -> primary ('^' unary)?
 * primary    -> '(' expression ')' | identifier '(' expression ')' | identifier | number
 * ```
 *
 * 运算优先级（从低到高）：
 * 1. 加法、减法
 * 2. 乘法、除法
 * 3. 幂运算（右结合）
 * 4. 一元正负
 * 5. 函数调用、原子表达式
 */

class Parser : public BaseParser {
public:
    explicit Parser(std::string_view source) : BaseParser(source) {}

    /**
     * @brief 解析表达式
     * @return 解析后的符号表达式
     * @throw std::runtime_error 当语法错误时
     */
    SymbolicExpression parse() {
        SymbolicExpression expression = parse_expression();
        skip_spaces();
        if (pos_ != source_.size()) {
            throw SyntaxError("unexpected token near: " + std::string(source_.substr(pos_, 1)));
        }
        return expression;
    }

private:
    int depth_ = 0;
    static constexpr int kMaxDepth = 256;

    struct DepthGuard {
        int& depth_;
        explicit DepthGuard(int& depth) : depth_(depth) {
            if (++depth_ > kMaxDepth) {
                throw std::runtime_error("expression too complex: maximum parsing depth exceeded");
            }
        }
        ~DepthGuard() { --depth_; }
    };

    // ========================================================================
    // 解析规则实现
    // ========================================================================

    /**
     * @brief 解析加减法表达式
     *
     * expression -> term (('+' | '-') term)*
     */
    SymbolicExpression parse_expression() {
        DepthGuard guard(depth_);
        SymbolicExpression value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                // 加法在解析时简化，减少后续工作量
                value = SymbolicExpression(make_binary(NodeType::kAdd, value.simplify().node_, parse_term().simplify().node_));
            } else if (match('-')) {
                value = SymbolicExpression(make_binary(NodeType::kSubtract, value.simplify().node_, parse_term().simplify().node_));
            } else {
                return value;
            }
        }
    }

    /**
     * @brief 解析乘除法表达式
     *
     * term -> unary (('*' | '/') unary)*
     */
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

    /**
     * @brief 解析幂运算表达式
     *
     * power -> primary ('^' unary)?
     * 注意：幂运算是右结合的
     */
    SymbolicExpression parse_power() {
        SymbolicExpression value = parse_primary();
        skip_spaces();
        if (match('^')) {
            // 幂运算右结合：x^y^z = x^(y^z)
            return SymbolicExpression(make_binary(NodeType::kPower, value.node_, parse_unary().node_));
        }
        return value;
    }

    /**
     * @brief 解析一元表达式
     *
     * unary -> ('+' | '-') unary | power
     */
    SymbolicExpression parse_unary() {
        skip_spaces();
        if (match('+')) {
            // 正号可忽略
            return parse_unary();
        }
        if (match('-')) {
            // 负号创建取负节点
            return SymbolicExpression(make_unary(NodeType::kNegate, parse_unary().node_));
        }
        return parse_power();
    }

    /**
     * @brief 解析基本表达式
     *
     * primary -> '(' expression ')' | identifier '(' expression ')' | identifier | number
     */
    SymbolicExpression parse_primary() {
        skip_spaces();

        // 括号表达式
        if (match('(')) {
            SymbolicExpression value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        // 标识符（变量或函数）
        if (peek_is_alpha()) {
            std::string identifier(parse_identifier());

            // 特殊常量
            if (identifier == "pi") {
                return SymbolicExpression(make_unary(NodeType::kPi, nullptr));
            }
            if (identifier == "e") {
                return SymbolicExpression(make_unary(NodeType::kE, nullptr));
            }
            if (identifier == "inf" || identifier == "infinity" || identifier == "oo") {
                return SymbolicExpression(make_infinity(true));
            }
            if (identifier == "nan") {
                return SymbolicExpression::number(0.0L / 0.0L);
            }

            skip_spaces();

            // 函数调用
            if (match('(')) {
                // 别名规范化
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

            // 变量
            return SymbolicExpression(make_variable(identifier));
        }

        // 数值
        return parse_number();
    }

    /**
     * @brief 解析数值
     *
     * 支持整数和小数，不支持科学计数法。
     */
    SymbolicExpression parse_number() {
        skip_spaces();
        const std::size_t start = pos_;
        bool has_digit = false;
        bool seen_dot = false;

        // 收集数字字符
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
            throw SyntaxError("expected number");
        }

        // 手动解析数值（避免 strtod 的依赖和区域设置问题）
        long double value = 0.0L;
        std::size_t index = start;
        while (index < pos_ && source_[index] != '.') {
            value = value * 10.0L + static_cast<long double>(source_[index] - '0');
            ++index;
        }
        if (index < pos_ && source_[index] == '.') {
            ++index;
            long double place = 0.1;
            while (index < pos_) {
                value += static_cast<long double>(source_[index] - '0') * place;
                place *= 0.1;
                ++index;
            }
        }

        // 解析科学计数法
        if (pos_ < source_.size() && (source_[pos_] == 'e' || source_[pos_] == 'E')) {
            std::size_t exp_pos = pos_ + 1;
            bool exp_negative = false;
            if (exp_pos < source_.size() && source_[exp_pos] == '+') {
                ++exp_pos;
            } else if (exp_pos < source_.size() && source_[exp_pos] == '-') {
                exp_negative = true;
                ++exp_pos;
            }
            
            if (exp_pos < source_.size() && std::isdigit(static_cast<unsigned char>(source_[exp_pos]))) {
                long double exponent = 0.0L;
                while (exp_pos < source_.size() && std::isdigit(static_cast<unsigned char>(source_[exp_pos]))) {
                    exponent = exponent * 10.0L + static_cast<long double>(source_[exp_pos] - '0');
                    ++exp_pos;
                }
                pos_ = exp_pos;
                if (exp_negative) {
                    value /= mymath::pow(10.0L, exponent);
                } else {
                    value *= mymath::pow(10.0L, exponent);
                }
            }
        }

        return SymbolicExpression(make_number(value));
    }
};

bool expr_is_number(const SymbolicExpression& expression, long double* value);

SymbolicExpression simplify_impl(const SymbolicExpression& expression);
SymbolicExpression simplify_once(const SymbolicExpression& expression);

SymbolicExpression substitute_impl(const SymbolicExpression& expression,
                                  const std::string& variable_name,
                                  const SymbolicExpression& replacement) {
    static thread_local int depth = 0;
    static constexpr int kMaxDepth = 512;
    if (++depth > kMaxDepth) {
        --depth;
        throw std::runtime_error("symbolic substitution too complex: maximum depth exceeded");
    }

    struct DepthGuard {
        int* d;
        ~DepthGuard() { if (d) (*d)--; }
    } guard{&depth};

    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
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
        case NodeType::kVector: {
            std::vector<SymbolicExpression> components;
            for (const auto& child : node->children) {
                components.push_back(substitute_impl(SymbolicExpression(child),
                                                     variable_name,
                                                     replacement));
            }
            return SymbolicExpression::vector(components).simplify();
        }
        case NodeType::kTensor: {
            std::vector<std::vector<SymbolicExpression>> rows;
            for (const auto& row_node : node->children) {
                std::vector<SymbolicExpression> row;
                for (const auto& child : row_node->children) {
                    row.push_back(substitute_impl(SymbolicExpression(child),
                                                  variable_name,
                                                  replacement));
                }
                rows.push_back(std::move(row));
            }
            return SymbolicExpression::tensor(rows).simplify();
        }
        case NodeType::kDifferentialOp:
            return SymbolicExpression(
                       make_unary(NodeType::kDifferentialOp,
                                  substitute_impl(SymbolicExpression(node->left),
                                                  variable_name,
                                                  replacement)
                                      .node_,
                                  node->text))
                .simplify();
        case NodeType::kRootOf:
            return expression;
    }
    throw std::runtime_error("unsupported symbolic substitution");
}

SymbolicExpression substitute_expression_impl(const SymbolicExpression& expression,
                                              const SymbolicExpression& target,
                                              const SymbolicExpression& replacement) {
    if (expressions_match(expression, target)) {
        return replacement;
    }
    const auto& node = expression.node_;
    switch (node->type) {
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
        case NodeType::kVariable:
            return expression;
        case NodeType::kNegate:
            return SymbolicExpression(
                       make_unary(NodeType::kNegate,
                                  substitute_expression_impl(SymbolicExpression(node->left),
                                                             target,
                                                             replacement)
                                      .node_))
                .simplify();
        case NodeType::kFunction:
            return SymbolicExpression(
                       make_unary(NodeType::kFunction,
                                  substitute_expression_impl(SymbolicExpression(node->left),
                                                             target,
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
                                   substitute_expression_impl(SymbolicExpression(node->left),
                                                              target,
                                                              replacement)
                                       .node_,
                                   substitute_expression_impl(SymbolicExpression(node->right),
                                                              target,
                                                              replacement)
                                       .node_))
                .simplify();
        case NodeType::kVector: {
            std::vector<SymbolicExpression> components;
            for (const auto& child : node->children) {
                components.push_back(substitute_expression_impl(SymbolicExpression(child),
                                                                target,
                                                                replacement));
            }
            return SymbolicExpression::vector(components).simplify();
        }
        case NodeType::kTensor: {
            std::vector<std::vector<SymbolicExpression>> rows;
            for (const auto& row_node : node->children) {
                std::vector<SymbolicExpression> row;
                for (const auto& child : row_node->children) {
                    row.push_back(substitute_expression_impl(SymbolicExpression(child),
                                                             target,
                                                             replacement));
                }
                rows.push_back(std::move(row));
            }
            return SymbolicExpression::tensor(rows).simplify();
        }
        case NodeType::kDifferentialOp:
            return SymbolicExpression(
                       make_unary(NodeType::kDifferentialOp,
                                  substitute_expression_impl(SymbolicExpression(node->left),
                                                             target,
                                                             replacement)
                                      .node_,
                                  node->text))
                .simplify();
        case NodeType::kRootOf:
            return expression;
    }
    throw std::runtime_error("unsupported symbolic substitution");
}

bool try_evaluate_numeric_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                               long double* value) {
    switch (node->type) {
        case NodeType::kNumber:
            *value = node->number_value;
            return true;
        case NodeType::kPi:
        case NodeType::kE:
            return false;
        case NodeType::kInfinity:
            *value = (node->number_value > 0) ?
                mymath::infinity() :
                -mymath::infinity();
            return true;
        case NodeType::kVariable:
            if (node->text == "1 / 2") {
                *value = 0.5;
                return true;
            }
            return false;
        case NodeType::kNegate: {
            long double operand = 0.0L;
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
            long double left = 0.0L;
            long double right = 0.0L;
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
                case NodeType::kPi:
                case NodeType::kE:
                case NodeType::kInfinity:
                case NodeType::kVariable:
                case NodeType::kNegate:
                case NodeType::kFunction:
                case NodeType::kVector:
                case NodeType::kTensor:
                case NodeType::kDifferentialOp:
                case NodeType::kRootOf:
                    break;
            }
            return false;
        }
        case NodeType::kFunction: {
            long double argument = 0.0L;
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
                if (argument <= 0.0L) {
                    return false;
                }
                *value = mymath::ln(argument);
                return true;
            }
            if (node->text == "sqrt") {
                if (argument < 0.0L) {
                    return false;
                }
                long double root = mymath::sqrt(argument);
                // Only evaluate to number if it's a perfect square
                if (mymath::is_near_zero(root * root - argument, 1e-12) && mymath::is_integer(root, 1e-10)) {
                    *value = root;
                    return true;
                }
                return false;
            }
            if (node->text == "erf") {
                *value = mymath::erf(argument);
                return true;
            }
            if (node->text == "erfc") {
                *value = mymath::erfc(argument);
                return true;
            }
            if (node->text == "gamma") {
                *value = mymath::gamma(argument);
                return true;
            }
            if (node->text == "abs") {
                *value = mymath::abs(argument);
                return true;
            }
            if (node->text == "floor") {
                *value = static_cast<long double>(
                    static_cast<long long>(argument < 0.0L &&
                                                   static_cast<long double>(static_cast<long long>(argument)) != argument
                                               ? argument - 1.0L
                                               : argument));
                return true;
            }
            if (node->text == "ceil") {
                long long truncated = static_cast<long long>(argument);
                if (argument > 0.0L && static_cast<long double>(truncated) != argument) {
                    ++truncated;
                }
                *value = static_cast<long double>(truncated);
                return true;
            }
            if (node->text == "cbrt") {
                *value = mymath::cbrt(argument);
                return true;
            }
            if (node->text == "sign") {
                if (mymath::is_near_zero(argument, kFormatEps)) {
                    *value = 0.0L;
                } else {
                    *value = argument > 0.0L ? 1.0L : -1.0L;
                }
                return true;
            }
            if (node->text == "step") {
                *value = argument >= 0.0L ? 1.0L : 0.0L;
                return true;
            }
            if (node->text == "delta") {
                *value = mymath::is_near_zero(argument, kFormatEps) ? 1.0L : 0.0L;
                return true;
            }
            return false;
        }
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
            return false;
    }
    return false;
}

/**
 * @brief Evaluate expression node with dual numbers for automatic differentiation.
 *
 * For the differentiation variable, returns dual(point_value, 1.0L).
 * For other variables, returns dual(their_value, 0.0L).
 * All operations propagate derivatives using dual arithmetic.
 */
bool try_evaluate_dual_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                            const DualEvaluationContext& ctx,
                            mymath::dual<long double>* result) {
    switch (node->type) {
        case NodeType::kNumber:
            *result = mymath::dual<long double>(node->number_value, 0.0L);
            return true;
        case NodeType::kPi:
            *result = mymath::dual<long double>(mymath::kPi, 0.0L);
            return true;
        case NodeType::kE:
            *result = mymath::dual<long double>(mymath::kE, 0.0L);
            return true;
        case NodeType::kInfinity:
            *result = mymath::dual<long double>(
                (node->number_value > 0) ? mymath::infinity() : -mymath::infinity(), 0.0L);
            return true;
        case NodeType::kVariable: {
            if (node->text == ctx.differentiation_variable) {
                *result = mymath::dual<long double>(ctx.point_value, 1.0L);
                return true;
            }
            auto it = ctx.other_values.find(node->text);
            if (it != ctx.other_values.end()) {
                *result = mymath::dual<long double>(it->second, 0.0L);
                return true;
            }
            if (node->text == "1 / 2") {
                *result = mymath::dual<long double>(0.5, 0.0L);
                return true;
            }
            return false;
        }
        case NodeType::kNegate: {
            mymath::dual<long double> operand;
            if (!try_evaluate_dual_node(node->left, ctx, &operand)) {
                return false;
            }
            *result = -operand;
            return true;
        }
        case NodeType::kAdd:
        case NodeType::kSubtract:
        case NodeType::kMultiply:
        case NodeType::kDivide:
        case NodeType::kPower: {
            mymath::dual<long double> left, right;
            if (!try_evaluate_dual_node(node->left, ctx, &left) ||
                !try_evaluate_dual_node(node->right, ctx, &right)) {
                return false;
            }

            switch (node->type) {
                case NodeType::kAdd:
                    *result = left + right;
                    return true;
                case NodeType::kSubtract:
                    *result = left - right;
                    return true;
                case NodeType::kMultiply:
                    *result = left * right;
                    return true;
                case NodeType::kDivide:
                    *result = left / right;
                    return true;
                case NodeType::kPower:
                    *result = mymath::pow(left, right);
                    return true;
                case NodeType::kNumber:
                case NodeType::kPi:
                case NodeType::kE:
                case NodeType::kInfinity:
                case NodeType::kVariable:
                case NodeType::kNegate:
                case NodeType::kFunction:
                case NodeType::kVector:
                case NodeType::kTensor:
                case NodeType::kDifferentialOp:
                case NodeType::kRootOf:
                    break;
            }
            return false;
        }
        case NodeType::kFunction: {
            mymath::dual<long double> argument;
            if (!try_evaluate_dual_node(node->left, ctx, &argument)) {
                return false;
            }

            if (node->text == "sin") {
                *result = mymath::sin(argument);
                return true;
            }
            if (node->text == "cos") {
                *result = mymath::cos(argument);
                return true;
            }
            if (node->text == "tan") {
                *result = mymath::tan(argument);
                return true;
            }
            if (node->text == "sec") {
                *result = mymath::sec(argument);
                return true;
            }
            if (node->text == "csc") {
                *result = mymath::csc(argument);
                return true;
            }
            if (node->text == "cot") {
                *result = mymath::cot(argument);
                return true;
            }
            if (node->text == "asin") {
                *result = mymath::asin(argument);
                return true;
            }
            if (node->text == "acos") {
                *result = mymath::acos(argument);
                return true;
            }
            if (node->text == "atan") {
                *result = mymath::atan(argument);
                return true;
            }
            if (node->text == "sinh") {
                *result = mymath::sinh(argument);
                return true;
            }
            if (node->text == "cosh") {
                *result = mymath::cosh(argument);
                return true;
            }
            if (node->text == "tanh") {
                *result = mymath::tanh(argument);
                return true;
            }
            if (node->text == "asinh") {
                *result = mymath::asinh(argument);
                return true;
            }
            if (node->text == "acosh") {
                *result = mymath::acosh(argument);
                return true;
            }
            if (node->text == "atanh") {
                *result = mymath::atanh(argument);
                return true;
            }
            if (node->text == "exp") {
                *result = mymath::exp(argument);
                return true;
            }
            if (node->text == "ln" || node->text == "log") {
                *result = mymath::ln(argument);
                return true;
            }
            if (node->text == "log10") {
                *result = mymath::log10(argument);
                return true;
            }
            if (node->text == "sqrt") {
                *result = mymath::sqrt(argument);
                return true;
            }
            if (node->text == "cbrt") {
                *result = mymath::cbrt(argument);
                return true;
            }
            if (node->text == "abs") {
                *result = mymath::abs(argument);
                return true;
            }
            if (node->text == "erf") {
                *result = mymath::erf(argument);
                return true;
            }
            if (node->text == "erfc") {
                *result = mymath::erfc(argument);
                return true;
            }
            if (node->text == "sign") {
                *result = mymath::sign(argument);
                return true;
            }
            if (node->text == "floor") {
                *result = mymath::floor(argument);
                return true;
            }
            if (node->text == "ceil") {
                *result = mymath::ceil(argument);
                return true;
            }
            if (node->text == "round") {
                *result = mymath::round(argument);
                return true;
            }
            if (node->text == "trunc") {
                *result = mymath::trunc(argument);
                return true;
            }
            return false;
        }
        case NodeType::kVector:
        case NodeType::kTensor:
        case NodeType::kDifferentialOp:
        case NodeType::kRootOf:
            return false;
    }
    return false;
}

bool expr_is_variable(const SymbolicExpression& expression, const std::string& name) {
    if (name == "pi") return expression.node_->type == NodeType::kPi;
    if (name == "e") return expression.node_->type == NodeType::kE;
    return expression.node_->type == NodeType::kVariable && expression.node_->text == name;
}

}  // namespace symbolic_expression_internal

using namespace symbolic_expression_internal;

SymbolicExpression::SymbolicExpression() : node_(make_number(0.0L)) {}

SymbolicExpression::SymbolicExpression(std::shared_ptr<Node> node)
    : node_(std::move(node)) {}

SymbolicExpression SymbolicExpression::parse(const std::string& text) {
    Parser parser(text);
    return parser.parse().simplify();
}

SymbolicExpression SymbolicExpression::number(long double value) {
    return SymbolicExpression(make_number(value));
}

void SymbolicExpression::set_display_precision(int precision) {
    mutable_display_precision() = clamp_display_precision(precision);
}

SymbolicExpression operator+(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_add(lhs, rhs);
}

SymbolicExpression operator-(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_subtract(lhs, rhs);
}

SymbolicExpression operator*(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_multiply(lhs, rhs);
}

SymbolicExpression operator/(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_divide(lhs, rhs);
}

SymbolicExpression operator^(const SymbolicExpression& lhs, const SymbolicExpression& rhs) {
    return make_power(lhs, rhs);
}

SymbolicExpression operator-(const SymbolicExpression& expr) {
    return make_negate(expr);
}

SymbolicExpression SymbolicExpression::variable(const std::string& name) {
    return SymbolicExpression(make_variable(name));
}

std::string to_latex_impl(const std::shared_ptr<SymbolicExpression::Node>& node) {
    if (!node) return "";
    switch (node->type) {
        case NodeType::kNumber: {
            long double val = 0.0L;
            try_evaluate_numeric_node(node, &val);
            return format_decimal(val);
        }
        case NodeType::kPi: return "\\pi";
        case NodeType::kE: return "e";
        case NodeType::kInfinity: return "\\infty";
        case NodeType::kVariable: return node->text;
        case NodeType::kNegate: return "-" + to_latex_impl(node->left);
        case NodeType::kAdd: return to_latex_impl(node->left) + " + " + to_latex_impl(node->right);
        case NodeType::kSubtract: return to_latex_impl(node->left) + " - " + to_latex_impl(node->right);
        case NodeType::kMultiply: {
            std::string left_str = to_latex_impl(node->left);
            std::string right_str = to_latex_impl(node->right);
            // 简单逻辑：如果右侧是变量，可以省略乘号
            if (node->right->type == NodeType::kVariable) return left_str + right_str;
            return left_str + " \\cdot " + right_str;
        }
        case NodeType::kDivide: return "\\frac{" + to_latex_impl(node->left) + "}{" + to_latex_impl(node->right) + "}";
        case NodeType::kPower: {
            // e^x -> e^x
            if (node->left->type == NodeType::kE) return "e^{" + to_latex_impl(node->right) + "}";
            return "{" + to_latex_impl(node->left) + "}^{" + to_latex_impl(node->right) + "}";
        }
        case NodeType::kVector: {
            std::string res = "\\begin{pmatrix} ";
            for (size_t i = 0; i < node->children.size(); ++i) {
                if (i > 0) res += " \\\\ ";
                res += to_latex_impl(node->children[i]);
            }
            res += " \\end{pmatrix}";
            return res;
        }
        case NodeType::kTensor: {
            std::string res = "\\begin{pmatrix} ";
            for (size_t i = 0; i < node->children.size(); ++i) {
                if (i > 0) res += " \\\\ ";
                // 每行是一个 kVector
                if (node->children[i]->type == NodeType::kVector) {
                    for (size_t j = 0; j < node->children[i]->children.size(); ++j) {
                        if (j > 0) res += " & ";
                        res += to_latex_impl(node->children[i]->children[j]);
                    }
                }
            }
            res += " \\end{pmatrix}";
            return res;
        }
        case NodeType::kFunction: {
            if (node->text == "sqrt") return "\\sqrt{" + to_latex_impl(node->left) + "}";
            if (node->text == "sin") return "\\sin(" + to_latex_impl(node->left) + ")";
            if (node->text == "cos") return "\\cos(" + to_latex_impl(node->left) + ")";
            return "\\" + node->text + "(" + to_latex_impl(node->left) + ")";
        }
        default: return "???";
    }
}

std::string SymbolicExpression::to_latex() const {
    return to_latex_impl(simplify().node_);
}

std::string SymbolicExpression::to_string() const {
    return to_string_impl(simplify().node_, 0);
}

bool SymbolicExpression::is_constant(const std::string& variable_name) const {
    switch (node_->type) {
        case NodeType::kNumber:
        case NodeType::kPi:
        case NodeType::kE:
        case NodeType::kInfinity:
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
        case NodeType::kVector:
        case NodeType::kTensor:
            for (const auto& child : node_->children) {
                if (!SymbolicExpression(child).is_constant(variable_name)) {
                    return false;
                }
            }
            return true;
        case NodeType::kDifferentialOp:
            return SymbolicExpression(node_->left).is_constant(variable_name);
        case NodeType::kRootOf:
            return true;
    }
    return false;
}

bool SymbolicExpression::is_number(long double* value) const {
    long double numeric = 0.0L;
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
    std::vector<long double>* coefficients) const {
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

namespace {
void collect_subexpression_counts(const std::shared_ptr<SymbolicExpression::Node>& node,
                                  std::unordered_map<std::string, std::pair<std::shared_ptr<SymbolicExpression::Node>, int>>& counts) {
    if (!node) return;

    // 忽略叶子节点（数字、变量、常数、无穷大），因为它们作为 CSE 提取没有意义
    if (node->type != NodeType::kNumber && node->type != NodeType::kVariable &&
        node->type != NodeType::kPi && node->type != NodeType::kE &&
        node->type != NodeType::kInfinity) {
        const std::string key = node_structural_key(node);
        auto& entry = counts[key];
        if (entry.second == 0) {
            entry.first = node;
        }
        entry.second++;
    }

    collect_subexpression_counts(node->left, counts);
    collect_subexpression_counts(node->right, counts);
}
} // namespace

std::vector<std::pair<SymbolicExpression, int>> SymbolicExpression::common_subexpressions() const {
    std::unordered_map<std::string, std::pair<std::shared_ptr<Node>, int>> counts;
    collect_subexpression_counts(node_, counts);

    std::vector<std::pair<SymbolicExpression, int>> result;
    for (auto& [key, entry] : counts) {
        if (entry.second > 1) {
            result.push_back({SymbolicExpression(entry.first), entry.second});
        }
    }

    // 按出现次数降序排列，次数相同时按长度（通过字符串表示的长度估计）降序
    std::sort(result.begin(), result.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first.to_string().size() > b.first.to_string().size();
    });

    return result;
}

SymbolicExpression SymbolicExpression::cse_rewrite(std::vector<std::pair<std::string, SymbolicExpression>>& bindings, const std::string& prefix) const {
    auto subs = common_subexpressions();
    if (subs.empty()) {
        return *this;
    }
    
    // 按节点数（树大小）降序排序，优先替换大表达式
    // 使得提取出的临时变量更少但覆盖更大
    std::sort(subs.begin(), subs.end(), [](const auto& a, const auto& b) {
        std::size_t a_size = a.first.node_count();
        std::size_t b_size = b.first.node_count();
        if (a_size != b_size) return a_size > b_size;
        return a.second > b.second;
    });

    SymbolicExpression current = *this;
    int t_index = 0;
    for (const auto& sub : subs) {
        std::string t_name = prefix + std::to_string(t_index++);
        bindings.push_back({t_name, sub.first});
        current = current.substitute_expression(sub.first, SymbolicExpression::variable(t_name)).simplify();
        
        // 更新已经提取的表达式（它们可能包含刚刚提取的小表达式）
        for (std::size_t i = 0; i < bindings.size() - 1; ++i) {
            bindings[i].second = bindings[i].second.substitute_expression(sub.first, SymbolicExpression::variable(t_name)).simplify();
        }
    }
    return current;
}

SymbolicExpression SymbolicExpression::power(const SymbolicExpression& exponent) const {
    return make_power(*this, exponent);
}

SymbolicExpression SymbolicExpression::simplify() const {
    static constexpr std::size_t kMaxSimplifyCacheSize = 4096;
    static thread_local SymbolicExpressionLruCache cache(kMaxSimplifyCacheSize);

    const std::string key = node_structural_key(node_);
    SymbolicExpression cached;
    if (cache.get(key, &cached)) {
        return cached;
    }

    // 1. 检查表达式膨胀
    std::size_t count = node_count();
    if (count > 500) {
        // 自动触发 CSE 以控制膨胀
        std::vector<std::pair<std::string, SymbolicExpression>> bindings;
        SymbolicExpression cse_expr = cse_rewrite(bindings, "_auto_t");
        if (bindings.size() > 5) {
             SymbolicExpression simplified = simplify_with_budget(200); 
             cache.put(key, simplified);
             return simplified;
        }
    }

    SymbolicExpression simplified = simplify_impl(*this);
    cache.put(key, simplified);
    return simplified;
}

SymbolicExpression SymbolicExpression::simplify_with_budget(std::size_t max_nodes) const {
    static constexpr std::size_t kMaxSimplifyCacheSize = 4096;
    static thread_local SymbolicExpressionLruCache cache(kMaxSimplifyCacheSize);

    const std::string key = node_structural_key(node_) + "_budget_" + std::to_string(max_nodes);
    SymbolicExpression cached;
    if (cache.get(key, &cached)) {
        return cached;
    }

    SymbolicExpression simplified = simplify_with_budget_impl(*this, max_nodes);
    cache.put(key, simplified);
    return simplified;
}

std::size_t SymbolicExpression::node_count() const {
    return count_nodes(node_);
}

SymbolicExpression SymbolicExpression::expand() const {
    return expand_impl(*this);
}


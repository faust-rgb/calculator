#include "symbolic_expression_internal.h"

#include "mymath.h"

#include <stdexcept>

using namespace symbolic_expression_internal;

namespace {

bool is_one_plus_variable_squared(const SymbolicExpression& expression,
                                  const std::string& variable_name) {
    std::vector<double> coefficients;
    return expression.polynomial_coefficients(variable_name, &coefficients) &&
           coefficients.size() == 3 &&
           mymath::is_near_zero(coefficients[0] - 1.0, kFormatEps) &&
           mymath::is_near_zero(coefficients[1], kFormatEps) &&
           mymath::is_near_zero(coefficients[2] - 1.0, kFormatEps);
}

bool is_one_minus_variable_squared(const SymbolicExpression& expression,
                                   const std::string& variable_name) {
    std::vector<double> coefficients;
    return expression.polynomial_coefficients(variable_name, &coefficients) &&
           coefficients.size() == 3 &&
           mymath::is_near_zero(coefficients[0] - 1.0, kFormatEps) &&
           mymath::is_near_zero(coefficients[1], kFormatEps) &&
           mymath::is_near_zero(coefficients[2] + 1.0, kFormatEps);
}

bool is_sqrt_one_minus_variable_squared(const SymbolicExpression& expression,
                                        const std::string& variable_name) {
    if (expression.node_->type != NodeType::kFunction ||
        expression.node_->text != "sqrt") {
        return false;
    }
    return is_one_minus_variable_squared(SymbolicExpression(expression.node_->left),
                                         variable_name);
}

}  // namespace

SymbolicExpression SymbolicExpression::derivative(const std::string& variable_name) const {
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
            return make_add(make_multiply(left.derivative(variable_name), right),
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
                           make_divide(inner,
                                       make_function("sqrt",
                                                     make_subtract(number(1.0),
                                                                   make_power(argument, number(2.0))))))
                    .simplify();
            }
            if (node_->text == "atan") {
                return make_divide(inner,
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

std::vector<SymbolicExpression> SymbolicExpression::gradient(
    const std::vector<std::string>& variable_names) const {
    std::vector<SymbolicExpression> result;
    result.reserve(variable_names.size());
    for (const std::string& variable_name : variable_names) {
        result.push_back(derivative(variable_name).simplify());
    }
    return result;
}

std::vector<std::vector<SymbolicExpression>> SymbolicExpression::hessian(
    const std::vector<std::string>& variable_names) const {
    std::vector<std::vector<SymbolicExpression>> result;
    result.reserve(variable_names.size());
    for (const std::string& outer_variable : variable_names) {
        std::vector<SymbolicExpression> row;
        row.reserve(variable_names.size());
        const SymbolicExpression outer_derivative =
            derivative(outer_variable).simplify();
        for (const std::string& inner_variable : variable_names) {
            row.push_back(outer_derivative.derivative(inner_variable).simplify());
        }
        result.push_back(row);
    }
    return result;
}

std::vector<std::vector<SymbolicExpression>> SymbolicExpression::jacobian(
    const std::vector<SymbolicExpression>& expressions,
    const std::vector<std::string>& variable_names) {
    std::vector<std::vector<SymbolicExpression>> result;
    result.reserve(expressions.size());
    for (const SymbolicExpression& expression : expressions) {
        result.push_back(expression.gradient(variable_names));
    }
    return result;
}

SymbolicExpression SymbolicExpression::integral(const std::string& variable_name) const {
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
            if (is_one_plus_variable_squared(right, variable_name)) {
                return make_function("atan", variable(variable_name)).simplify();
            }
            if (is_sqrt_one_minus_variable_squared(right, variable_name)) {
                return make_function("asin", variable(variable_name)).simplify();
            }
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
        if (node_->text == "sqrt" &&
            is_one_minus_variable_squared(argument, variable_name)) {
            const SymbolicExpression x = variable(variable_name);
            return make_divide(
                       make_add(make_multiply(x, make_function("sqrt", argument)),
                                make_function("asin", x)),
                       number(2.0))
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

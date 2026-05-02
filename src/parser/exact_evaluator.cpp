// ============================================================================
// 精确有理数 AST 求值器
// ============================================================================
//
// 提供表达式 AST 的精确有理数求值功能。
// ============================================================================

#include "parser/unified_expression_parser.h"
#include "command/expression_ast.h"
#include "calculator_exceptions.h"
#include "variable_resolver.h"
#include "types/function.h"
#include "precise/rational.h"
#include "math/helpers/base_conversions.h"
#include "math/helpers/bitwise_helpers.h"
#include "math/helpers/combinatorics.h"
#include "math/helpers/integer_helpers.h"
#include "core/calculator_internal_types.h"
#include "mymath.h"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace {

template <typename ExceptionType = std::runtime_error>
[[noreturn]] void throw_ast_error(const std::string& message, std::size_t pos) {
    std::ostringstream oss;
    oss << message << " at position " << pos;
    throw ExceptionType(oss.str());
}

Rational parse_rational_literal(const std::string& token) {
    std::string significand = token;
    long long exponent_adjust = 0;
    const std::size_t exponent_pos = token.find_first_of("eE");
    if (exponent_pos != std::string::npos) {
        significand = token.substr(0, exponent_pos);
        exponent_adjust = std::stoll(token.substr(exponent_pos + 1));
    }

    long long numerator = 0;
    long long denominator = 1;
    std::size_t idx = 0;

    // Handle optional sign
    int sign = 1;
    if (idx < significand.size() && significand[idx] == '-') {
        sign = -1;
        ++idx;
    } else if (idx < significand.size() && significand[idx] == '+') {
        ++idx;
    }

    while (idx < significand.size() && significand[idx] != '.') {
        numerator =
            numerator * 10 + static_cast<long long>(significand[idx] - '0');
        ++idx;
    }

    if (idx < significand.size() && significand[idx] == '.') {
        ++idx;
        while (idx < significand.size()) {
            numerator =
                numerator * 10 + static_cast<long long>(significand[idx] - '0');
            denominator *= 10;
            ++idx;
        }
    }

    numerator *= sign;

    while (exponent_adjust > 0) {
        numerator *= 10;
        --exponent_adjust;
    }
    while (exponent_adjust < 0) {
        denominator *= 10;
        ++exponent_adjust;
    }

    return Rational(numerator, denominator);
}

Rational lookup_variable_exact(const std::string& name, const VariableResolver& variables, std::size_t pos) {
    const StoredValue* found = variables.lookup(name);
    if (!found) {
        double constant_value = 0.0;
        if (lookup_builtin_constant(name, &constant_value)) {
            throw_ast_error<ExactModeUnsupported>("built-in constants are not rational", pos);
        }
        throw_ast_error<std::runtime_error>("unknown variable: " + name, pos);
    }
    if (found->is_matrix || found->is_complex) {
        throw_ast_error<ExactModeUnsupported>("matrix or complex variable " + name + " cannot be used exactly", pos);
    }
    if (found->is_string) {
        throw_ast_error<ExactModeUnsupported>("string variable " + name + " cannot be used exactly", pos);
    }
    if (!found->exact) {
        throw_ast_error<ExactModeUnsupported>("variable " + name + " is only stored approximately", pos);
    }
    return found->rational;
}

Rational apply_function_exact(const std::string& name,
                        const std::vector<Rational>& arguments,
                        const std::map<std::string, CustomFunction>* functions,
                        HasScriptFunctionCallback has_script_function,
                        std::size_t pos) {
    if (functions && functions->find(name) != functions->end()) {
        throw_ast_error<ExactModeUnsupported>("custom function " + name +
                                   " is not supported exactly", pos);
    }
    if (has_script_function && has_script_function(name)) {
        throw_ast_error<ExactModeUnsupported>("script function " + name +
                                   " is not supported exactly", pos);
    }
    if (name == "pow") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("pow expects exactly two arguments", pos);
        }
        if (!arguments[1].is_integer()) {
            throw_ast_error<ExactModeUnsupported>("exact rational mode does not support non-integer exponents", pos);
        }
        return pow_rational(arguments[0], arguments[1].numerator);
    }
    if (name == "abs") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("abs expects exactly one argument", pos);
        }
        return abs_rational(arguments[0]);
    }
    if (name == "step" || name == "u" || name == "heaviside") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("step expects exactly one argument", pos);
        }
        return Rational(arguments[0].numerator >= 0 ? 1 : 0, 1);
    }
    if (name == "delta" || name == "impulse") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("delta expects exactly one argument", pos);
        }
        return Rational(arguments[0].numerator == 0 ? 1 : 0, 1);
    }
    if (name == "not") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("not expects exactly one argument", pos);
        }
        if (!arguments[0].is_integer()) {
            throw_ast_error<std::runtime_error>("not only accepts integers", pos);
        }
        return Rational(~arguments[0].numerator, 1);
    }
    if (name == "sign") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("sign expects exactly one argument", pos);
        }
        if (arguments[0].numerator == 0) {
            return Rational(0, 1);
        }
        return Rational(arguments[0].numerator > 0 ? 1 : -1, 1);
    }
    if (name == "gcd") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("gcd expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("gcd only accepts integers", pos);
        }
        return Rational(gcd_ll(arguments[0].numerator, arguments[1].numerator), 1);
    }
    if (name == "lcm") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("lcm expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("lcm only accepts integers", pos);
        }
        return Rational(lcm_ll(arguments[0].numerator, arguments[1].numerator), 1);
    }
    if (name == "mod") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("mod expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("mod only accepts integers", pos);
        }
        if (arguments[1].numerator == 0) {
            throw_ast_error<std::runtime_error>("mod divisor cannot be zero", pos);
        }
        return Rational(arguments[0].numerator % arguments[1].numerator, 1);
    }
    if (name == "rol") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("rol expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("rol only accepts integers", pos);
        }
        const unsigned count = normalize_rotation_count(arguments[1].numerator);
        return Rational(
            from_unsigned_bits(rotate_left_bits(
                to_unsigned_bits(arguments[0].numerator), count)),
            1);
    }
    if (name == "ror") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("ror expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("ror only accepts integers", pos);
        }
        const unsigned count = normalize_rotation_count(arguments[1].numerator);
        return Rational(
            from_unsigned_bits(rotate_right_bits(
                to_unsigned_bits(arguments[0].numerator), count)),
            1);
    }
    if (name == "floor") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("floor expects exactly one argument", pos);
        }
        return Rational(floor_to_long_long(rational_to_double(arguments[0])), 1);
    }
    if (name == "ceil") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("ceil expects exactly one argument", pos);
        }
        return Rational(ceil_to_long_long(rational_to_double(arguments[0])), 1);
    }
    if (name == "round") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("round expects exactly one argument", pos);
        }
        return Rational(round_to_long_long(rational_to_double(arguments[0])), 1);
    }
    if (name == "trunc") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("trunc expects exactly one argument", pos);
        }
        return Rational(trunc_to_long_long(rational_to_double(arguments[0])), 1);
    }
    if (name == "min") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("min expects exactly two arguments", pos);
        }
        return rational_to_double(arguments[0]) < rational_to_double(arguments[1])
                   ? arguments[0]
                   : arguments[1];
    }
    if (name == "max") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("max expects exactly two arguments", pos);
        }
        return rational_to_double(arguments[0]) > rational_to_double(arguments[1])
                   ? arguments[0]
                   : arguments[1];
    }
    if (name == "clamp") {
        if (arguments.size() != 3) {
            throw_ast_error<std::runtime_error>("clamp expects exactly three arguments", pos);
        }
        Rational lower = arguments[1];
        Rational upper = arguments[2];
        if (rational_to_double(lower) > rational_to_double(upper)) {
            std::swap(lower, upper);
        }
        if (rational_to_double(arguments[0]) < rational_to_double(lower)) {
            return lower;
        }
        if (rational_to_double(arguments[0]) > rational_to_double(upper)) {
            return upper;
        }
        return arguments[0];
    }
    if (name == "sum") {
        if (arguments.empty()) {
            throw_ast_error<std::runtime_error>("sum expects at least one argument", pos);
        }
        Rational total(0, 1);
        for (const Rational& value : arguments) {
            total = total + value;
        }
        return total;
    }
    if (name == "avg" || name == "mean") {
        if (arguments.empty()) {
            throw_ast_error<std::runtime_error>("mean expects at least one argument", pos);
        }
        Rational total(0, 1);
        for (const Rational& value : arguments) {
            total = total + value;
        }
        return total / Rational(static_cast<long long>(arguments.size()), 1);
    }
    if (name == "median") {
        if (arguments.empty()) {
            throw_ast_error<std::runtime_error>("median expects at least one argument", pos);
        }
        std::vector<Rational> sorted = arguments;
        std::sort(sorted.begin(), sorted.end(),
                  [](const Rational& lhs, const Rational& rhs) {
                      return rational_to_double(lhs) < rational_to_double(rhs);
                  });
        const std::size_t middle = sorted.size() / 2;
        if (sorted.size() % 2 == 1) {
            return sorted[middle];
        }
        return (sorted[middle - 1] + sorted[middle]) / Rational(2, 1);
    }
    if (name == "factorial") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("factorial expects exactly one argument", pos);
        }
        if (!arguments[0].is_integer()) {
            throw_ast_error<std::runtime_error>("factorial only accepts integers", pos);
        }
        return factorial_rational(arguments[0].numerator);
    }
    if (name == "nCr" || name == "binom") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("nCr expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("nCr only accepts integers", pos);
        }
        return combination_rational(arguments[0].numerator, arguments[1].numerator);
    }
    if (name == "nPr") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("nPr expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("nPr only accepts integers", pos);
        }
        return permutation_rational(arguments[0].numerator, arguments[1].numerator);
    }
    if (name == "popcount") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("popcount expects exactly one argument", pos);
        }
        if (!arguments[0].is_integer()) {
            throw_ast_error<std::runtime_error>("popcount only accepts integers", pos);
        }
        return Rational(popcount_bits(to_unsigned_bits(arguments[0].numerator)), 1);
    }
    if (name == "bitlen") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("bitlen expects exactly one argument", pos);
        }
        if (!arguments[0].is_integer()) {
            throw_ast_error<std::runtime_error>("bitlen only accepts integers", pos);
        }
        return Rational(bit_length_bits(to_unsigned_bits(arguments[0].numerator)), 1);
    }
    if (name == "ctz") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("ctz expects exactly one argument", pos);
        }
        if (!arguments[0].is_integer()) {
            throw_ast_error<std::runtime_error>("ctz only accepts integers", pos);
        }
        return Rational(trailing_zero_count_bits(to_unsigned_bits(arguments[0].numerator)), 1);
    }
    if (name == "clz") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("clz expects exactly one argument", pos);
        }
        if (!arguments[0].is_integer()) {
            throw_ast_error<std::runtime_error>("clz only accepts integers", pos);
        }
        return Rational(leading_zero_count_bits(to_unsigned_bits(arguments[0].numerator)), 1);
    }
    if (name == "parity") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("parity expects exactly one argument", pos);
        }
        if (!arguments[0].is_integer()) {
            throw_ast_error<std::runtime_error>("parity only accepts integers", pos);
        }
        return Rational(parity_bits(to_unsigned_bits(arguments[0].numerator)), 1);
    }
    if (name == "reverse_bits") {
        if (arguments.size() != 1) {
            throw_ast_error<std::runtime_error>("reverse_bits expects exactly one argument", pos);
        }
        if (!arguments[0].is_integer()) {
            throw_ast_error<std::runtime_error>("reverse_bits only accepts integers", pos);
        }
        return Rational(
            from_unsigned_bits(reverse_bits(to_unsigned_bits(arguments[0].numerator))),
            1);
    }
    if (name == "and") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("and expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("and only accepts integers", pos);
        }
        return Rational(arguments[0].numerator & arguments[1].numerator, 1);
    }
    if (name == "or") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("or expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("or only accepts integers", pos);
        }
        return Rational(arguments[0].numerator | arguments[1].numerator, 1);
    }
    if (name == "xor") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("xor expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("xor only accepts integers", pos);
        }
        return Rational(arguments[0].numerator ^ arguments[1].numerator, 1);
    }
    if (name == "shl") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("shl expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("shl only accepts integers", pos);
        }
        if (arguments[1].numerator < 0) {
            throw_ast_error<std::runtime_error>("shift count cannot be negative", pos);
        }
        return Rational(arguments[0].numerator << arguments[1].numerator, 1);
    }
    if (name == "shr") {
        if (arguments.size() != 2) {
            throw_ast_error<std::runtime_error>("shr expects exactly two arguments", pos);
        }
        if (!arguments[0].is_integer() || !arguments[1].is_integer()) {
            throw_ast_error<std::runtime_error>("shr only accepts integers", pos);
        }
        if (arguments[1].numerator < 0) {
            throw_ast_error<std::runtime_error>("shift count cannot be negative", pos);
        }
        return Rational(arguments[0].numerator >> arguments[1].numerator, 1);
    }

    throw_ast_error<ExactModeUnsupported>("function " + name + " is not supported exactly", pos);
}

} // namespace

// ============================================================================
// 公开的精确 AST 求值函数
// ============================================================================

Rational evaluate_ast_exact(const ExpressionAST* ast,
                            const VariableResolver& variables,
                            const std::map<std::string, CustomFunction>* functions,
                            HasScriptFunctionCallback has_script_function) {
    if (!ast) {
        throw MathError("null AST node");
    }

    switch (ast->kind) {
        case ExprKind::kNumber: {
            if (ast->string_value.empty()) {
                throw_ast_error<std::runtime_error>("empty number literal", ast->position);
            }
            if (ast->string_value.size() > 2 && ast->string_value[0] == '0' &&
                std::isalpha(ast->string_value[1])) {
                // handle 0x, 0b, 0o
                return Rational(parse_prefixed_integer_token(ast->string_value), 1);
            }
            return parse_rational_literal(ast->string_value);
        }

        case ExprKind::kVariable:
            return lookup_variable_exact(ast->identifier, variables, ast->position);

        case ExprKind::kBinaryOp: {
            if (ast->children.size() != 2) {
                throw_ast_error<std::runtime_error>("invalid binary operation", ast->position);
            }
            Rational left = evaluate_ast_exact(ast->children[0].get(), variables, functions, has_script_function);
            Rational right = evaluate_ast_exact(ast->children[1].get(), variables, functions, has_script_function);

            switch (ast->op_char) {
                case '+': return left + right;
                case '-': return left - right;
                case '*': return left * right;
                case '/':
                    if (right.numerator == 0) throw_ast_error<std::runtime_error>("division by zero", ast->position);
                    return left / right;
                case '%':
                    if (right.numerator == 0) throw_ast_error<std::runtime_error>("modulo by zero", ast->position);
                    if (!left.is_integer() || !right.is_integer()) {
                        throw_ast_error<ExactModeUnsupported>("exact rational mode does not support non-integer modulo", ast->position);
                    }
                    return Rational(left.numerator % right.numerator, 1);
                case '^':
                    if (!right.is_integer()) {
                        throw_ast_error<ExactModeUnsupported>("exact rational mode does not support non-integer exponents", ast->position);
                    }
                    return pow_rational(left, right.numerator);
                default:
                    throw_ast_error<std::runtime_error>("unknown operator", ast->position);
            }
        }

        case ExprKind::kUnaryOp: {
            if (ast->children.size() != 1) {
                throw_ast_error<std::runtime_error>("invalid unary operation", ast->position);
            }
            Rational operand = evaluate_ast_exact(ast->children[0].get(), variables, functions, has_script_function);
            switch (ast->op_char) {
                case '-': return Rational(-operand.numerator, operand.denominator);
                case '+': return operand;
                default:
                    throw_ast_error<std::runtime_error>("unknown unary operator", ast->position);
            }
        }

        case ExprKind::kComparison: {
            if (ast->children.size() != 2) {
                throw_ast_error<std::runtime_error>("invalid comparison", ast->position);
            }
            Rational left = evaluate_ast_exact(ast->children[0].get(), variables, functions, has_script_function);
            Rational right = evaluate_ast_exact(ast->children[1].get(), variables, functions, has_script_function);

            if (ast->comparison_op == "==") return Rational(left.numerator * right.denominator == right.numerator * left.denominator ? 1 : 0, 1);
            if (ast->comparison_op == "!=") return Rational(left.numerator * right.denominator != right.numerator * left.denominator ? 1 : 0, 1);
            if (ast->comparison_op == "<") return Rational(rational_to_double(left) < rational_to_double(right) ? 1 : 0, 1);
            if (ast->comparison_op == ">") return Rational(rational_to_double(left) > rational_to_double(right) ? 1 : 0, 1);
            if (ast->comparison_op == "<=") return Rational(rational_to_double(left) <= rational_to_double(right) ? 1 : 0, 1);
            if (ast->comparison_op == ">=") return Rational(rational_to_double(left) >= rational_to_double(right) ? 1 : 0, 1);

            throw_ast_error<std::runtime_error>("unknown comparison operator", ast->position);
        }

        case ExprKind::kFunctionCall: {
            std::vector<Rational> args;
            args.reserve(ast->children.size());
            for (const auto& child : ast->children) {
                args.push_back(evaluate_ast_exact(child.get(), variables, functions, has_script_function));
            }
            return apply_function_exact(ast->identifier, args, functions, has_script_function, ast->position);
        }

        case ExprKind::kConditional: {
            if (ast->children.size() != 3) {
                throw_ast_error<std::runtime_error>("invalid conditional", ast->position);
            }
            Rational cond = evaluate_ast_exact(ast->children[0].get(), variables, functions, has_script_function);
            if (cond.numerator != 0) {
                return evaluate_ast_exact(ast->children[1].get(), variables, functions, has_script_function);
            }
            return evaluate_ast_exact(ast->children[2].get(), variables, functions, has_script_function);
        }

        default:
            throw_ast_error<std::runtime_error>("unknown AST node kind", ast->position);
    }
}

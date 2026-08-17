/**
 * @file calculator_polynomial.cpp
 * @brief 多项式操作命令实现
 *
 * 本文件实现了计算器中的多项式操作命令，作为多项式运算库与计算器用户界面之间的桥梁。
 * 主要功能包括：
 * - 多项式四则运算（加、减、乘、除）
 * - 多项式求根（实根和复根）
 * - 嵌套多项式表达式的解析和计算
 *
 * 支持的输入形式：
 * 1. 直接的符号表达式（如 "x^2 + 2*x + 1"）
 * 2. 自定义函数名（如 "f"，其中 f = x^2 - 1）
 * 3. 嵌套的多项式操作（如 "poly_add(poly_mul(p, q), r)"）
 */

#include "execution/engine/script_context.h"
#include "calculator_polynomial.h"

#include "types/scalar_type.h"
#include "parser/grammars/unified_expression_parser.h"
#include "matrix/public/matrix_format.h"
#include "polynomial.h"
#include "mymath.h"
#include "math/functions/integer/integer_helpers.h"
#include "parser/grammars/command_parser.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"

#include <sstream>

namespace polynomial_ops {

namespace {

using Scalar = mymath::Scalar;

// ============================================================================
// 内部辅助函数
// ============================================================================

bool try_parse_vector_coefficients(const std::string& input, std::vector<Scalar>* coeffs) {
    std::string s = trim_copy(input);
    if (s.empty()) return false;
    if (s.front() == '[' && s.back() == ']') {
        s = s.substr(1, s.size() - 2);
    } else if (s.size() >= 5 && s.substr(0, 4) == "vec(" && s.back() == ')') {
        s = s.substr(4, s.size() - 5);
    } else {
        return false;
    }
    for (char& c : s) {
        if (c == ',' || c == ';') c = ' ';
    }
    std::istringstream iss(s);
    std::string token;
    std::vector<Scalar> parsed;
    while (iss >> token) {
        SymbolicExpression expr = SymbolicExpression::parse(token);
        Scalar v = 0;
        if (expr.is_number(&v)) {
            parsed.push_back(v);
        } else {
            return false;
        }
    }
    if (!parsed.empty()) {
        *coeffs = parsed;
        return true;
    }
    return false;
}

/**
 * @brief 递归构建多项式系数
 *
 * 从参数字符串解析多项式，支持嵌套的多项式操作。
 *
 * @param ctx 多项式构建上下文，提供符号解析等功能
 * @param argument 参数字符串（表达式、函数名或嵌套调用）
 * @param variable_name 输出：多项式变量名
 * @param coefficients 输出：多项式系数向量（低次到高次）
 *
 * @throw std::runtime_error 当参数无法解析为多项式时抛出
 */
void build_polynomial_recursive(
    const PolynomialContext& ctx,
    const std::string& argument,
    std::string* variable_name,
    std::vector<Scalar>* coefficients) {

    const std::string trimmed_argument = trim_copy(argument);

    // 检查是否为向量字面量，如 [1, 2, 1] 或 vec(1, 2, 1)
    if (try_parse_vector_coefficients(trimmed_argument, coefficients)) {
        if (variable_name->empty()) {
            *variable_name = "x";
        }
        return;
    }

    CommandASTNode ast = parse_command(trimmed_argument);

    // 检查是否为嵌套的多项式操作
    if (ast.kind == CommandKind::kFunctionCall) {
        const auto* call = ast.as_function_call();
        if (call->name == "poly_add" || call->name == "poly_sub" ||
            call->name == "poly_mul" || call->name == "poly_div" ||
            call->name == "poly_deriv" || call->name == "poly_integ" ||
            call->name == "poly_compose" || call->name == "poly_gcd") {

            if (call->name == "poly_deriv" || call->name == "poly_integ") {
                if (call->arguments.empty() || call->arguments.size() > 2) {
                    throw std::runtime_error(std::string(call->name) + " expects 1 or 2 arguments");
                }
                std::string var;
                std::vector<Scalar> inner_coeffs;
                build_polynomial_recursive(ctx, std::string(call->arguments[0].text),
                                           &var, &inner_coeffs);
                *variable_name = var;
                if (call->name == "poly_deriv") {
                    *coefficients = polynomial_derivative(inner_coeffs);
                } else {
                    *coefficients = polynomial_integral(inner_coeffs);
                }
                return;
            }

            if (call->arguments.size() != 2) {
                throw std::runtime_error(
                    "polynomial binary operations expect exactly two arguments");
            }

            std::string lhs_variable;
            std::string rhs_variable;
            std::vector<Scalar> lhs_coefficients;
            std::vector<Scalar> rhs_coefficients;
            build_polynomial_recursive(ctx, std::string(call->arguments[0].text),
                                       &lhs_variable, &lhs_coefficients);
            build_polynomial_recursive(ctx, std::string(call->arguments[1].text),
                                       &rhs_variable, &rhs_coefficients);

            *variable_name = lhs_variable.empty() ? rhs_variable : lhs_variable;
            if (call->name == "poly_add") {
                *coefficients = polynomial_add(lhs_coefficients, rhs_coefficients);
                return;
            }
            if (call->name == "poly_sub") {
                *coefficients = polynomial_subtract(lhs_coefficients, rhs_coefficients);
                return;
            }
            if (call->name == "poly_compose") {
                *coefficients = polynomial_compose(lhs_coefficients, rhs_coefficients);
                return;
            }
            if (call->name == "poly_gcd") {
                *coefficients = polynomial_gcd(lhs_coefficients, rhs_coefficients);
                return;
            }
            if (call->name == "poly_div") {
                const PolynomialDivisionResult division =
                    polynomial_divide(lhs_coefficients, rhs_coefficients);
                bool zero_remainder = true;
                for (const Scalar& coefficient : division.remainder) {
                    if (mymath::abs(coefficient) >= Scalar(1e-10L)) {
                        zero_remainder = false;
                        break;
                    }
                }
                if (!zero_remainder) {
                    throw std::runtime_error(
                        "nested poly_div requires zero remainder");
                }
                *coefficients = division.quotient;
                return;
            }

            *coefficients = polynomial_multiply(lhs_coefficients, rhs_coefficients);
            return;
        }
    }

    // 从符号表达式构建多项式
    SymbolicExpression expression =
        ctx.resolve_symbolic(trimmed_argument, variable_name);
    if (!expression.polynomial_coefficients(*variable_name, coefficients)) {
        throw std::runtime_error("custom function " + trimmed_argument +
                                  " is not a polynomial");
    }
}

}  // namespace

// ============================================================================
// 公共接口实现
// ============================================================================

/**
 * @brief 从参数构建多项式数据
 *
 * 这是构建多项式的主入口函数，内部调用递归解析函数。
 *
 * @param ctx 多项式构建上下文
 * @param argument 参数字符串
 * @return 多项式数据结构，包含变量名和系数向量
 */
PolynomialData build_polynomial(const PolynomialContext& ctx,
                                const std::string& argument) {
    PolynomialData result;
    build_polynomial_recursive(ctx, argument,
                               &result.variable_name, &result.coefficients);
    return result;
}

// ============================================================================
// 多项式运算函数实现
// ============================================================================

/**
 * @brief 多项式加法运算
 *
 * 将两个多项式相加，返回格式化的字符串结果。
 *
 * @param lhs 左操作数（加数）
 * @param rhs 右操作数（加数）
 * @return 和的多项式字符串表示
 */
std::string poly_add(const PolynomialData& lhs, const PolynomialData& rhs) {
    return polynomial_to_string(
        polynomial_add(lhs.coefficients, rhs.coefficients),
        lhs.variable_name.empty() ? rhs.variable_name : lhs.variable_name);
}

/**
 * @brief 多项式减法运算
 *
 * 从第一个多项式中减去第二个多项式，返回格式化的字符串结果。
 *
 * @param lhs 左操作数（被减数）
 * @param rhs 右操作数（减数）
 * @return 差的多项式字符串表示
 */
std::string poly_sub(const PolynomialData& lhs, const PolynomialData& rhs) {
    return polynomial_to_string(
        polynomial_subtract(lhs.coefficients, rhs.coefficients),
        lhs.variable_name.empty() ? rhs.variable_name : lhs.variable_name);
}

/**
 * @brief 多项式乘法运算
 *
 * 将两个多项式相乘，返回格式化的字符串结果。
 *
 * @param lhs 左操作数（乘数）
 * @param rhs 右操作数（乘数）
 * @return 积的多项式字符串表示
 */
std::string poly_mul(const PolynomialData& lhs, const PolynomialData& rhs) {
    return polynomial_to_string(
        polynomial_multiply(lhs.coefficients, rhs.coefficients),
        lhs.variable_name.empty() ? rhs.variable_name : lhs.variable_name);
}

/**
 * @brief 多项式除法运算
 *
 * 使用多项式长除法算法，返回商和余数。
 *
 * @param lhs 被除数多项式
 * @param rhs 除数多项式
 * @return 格式化的商和余数字符串，格式为 "quotient: <商>, remainder: <余数>"
 */
std::string poly_div(const PolynomialData& lhs, const PolynomialData& rhs) {
    const PolynomialDivisionResult division =
        polynomial_divide(lhs.coefficients, rhs.coefficients);
    const std::string var = lhs.variable_name.empty() ? rhs.variable_name : lhs.variable_name;
    return "quotient: " +
           polynomial_to_string(division.quotient, var) +
           ", remainder: " +
           polynomial_to_string(division.remainder, var);
}

/**
 * @brief 计算多项式的所有根
 *
 * 计算多项式的全部实根和复根，返回格式化的字符串。
 *
 * @param poly 输入多项式
 * @return 格式化的根字符串，多个根用逗号分隔
 */
std::string roots(const PolynomialData& poly) {
    const std::vector<mymath::complex<Scalar>> roots =
        polynomial_complex_roots(poly.coefficients);
    if (roots.empty()) {
        return "No roots.";
    }

    std::ostringstream out;
    bool wrote_root = false;
    mymath::complex<Scalar> previous_root(0.0L, 0.0L);
    for (std::size_t i = 0; i < roots.size(); ++i) {
        const Scalar real_diff = roots[i].real() - previous_root.real();
        const Scalar imag_diff = roots[i].imag() - previous_root.imag();
        if (wrote_root &&
            mymath::abs(real_diff) <= Scalar(1e-7L) &&
            mymath::abs(imag_diff) <= Scalar(1e-7L)) {
            continue;
        }
        if (wrote_root) {
            out << ", ";
        }
        Scalar real = roots[i].real();
        Scalar imag = roots[i].imag();
        if (is_integer_double(real, 1e-6)) {
            real = static_cast<Scalar>(round_to_long_long(real));
        }
        if (is_integer_double(imag, 1e-6)) {
            imag = static_cast<Scalar>(round_to_long_long(imag));
        }
        if (is_complex_effectively_real(std::complex<long double>(real.to_long_double(), imag.to_long_double()))) {
            out << format_symbolic_scalar(real);
        } else {
            out << matrix_utils::format_complex<Scalar>({real, imag});
        }
        previous_root = {real, imag};
        wrote_root = true;
    }
    return out.str();
}

// ============================================================================
// 命令处理函数实现
// ============================================================================

/**
 * @brief 检查命令是否为多项式命令
 *
 * 判断给定命令字符串是否属于多项式操作命令。
 *
 * @param command 命令字符串
 * @return 如果是多项式命令返回 true，否则返回 false
 */
bool is_polynomial_command(const std::string& command) {
    return command == "poly_add" ||
           command == "poly_sub" ||
           command == "poly_mul" ||
           command == "poly_div" ||
           command == "roots" ||
           command == "poly_eval" ||
           command == "poly_deriv" ||
           command == "poly_integ" ||
           command == "poly_compose" ||
           command == "poly_gcd";
}

/**
 * @brief 处理多项式命令
 *
 * 解析命令参数并执行相应的多项式操作。支持嵌套的多项式表达式。
 *
 * @param ctx 多项式构建上下文，提供符号解析功能
 * @param command 命令名称
 * @param arguments 参数字符串列表
 * @param output 输出结果字符串指针
 * @return 是否成功处理命令
 */
bool handle_polynomial_command(const PolynomialContext& ctx,
                               const std::string& command,
                               const std::vector<std::string>& arguments,
                               std::string* output) {
    if (command == "roots") {
        if (arguments.size() != 1) {
            throw std::runtime_error("roots expects exactly one argument");
        }
        PolynomialData poly = build_polynomial(ctx, arguments[0]);
        *output = roots(poly);
        return true;
    }

    if (command == "poly_eval") {
        if (arguments.size() != 2) {
            throw std::runtime_error("poly_eval expects polynomial expression and x value");
        }
        PolynomialData poly = build_polynomial(ctx, arguments[0]);
        Scalar x = 0;
        SymbolicExpression x_expr = SymbolicExpression::parse(arguments[1]);
        if (!x_expr.is_number(&x)) {
            std::string temp_var;
            SymbolicExpression resolved = ctx.resolve_symbolic(arguments[1], &temp_var);
            if (!resolved.is_number(&x)) {
                throw std::runtime_error("poly_eval expects a numeric evaluation point");
            }
        }
        Scalar val = polynomial_evaluate(poly.coefficients, x);
        *output = format_decimal(val);
        return true;
    }

    if (command == "poly_deriv") {
        if (arguments.empty() || arguments.size() > 2) {
            throw std::runtime_error("poly_deriv expects polynomial and optional variable");
        }
        PolynomialData poly = build_polynomial(ctx, arguments[0]);
        *output = polynomial_to_string(polynomial_derivative(poly.coefficients),
                                       poly.variable_name.empty() ? "x" : poly.variable_name);
        return true;
    }

    if (command == "poly_integ") {
        if (arguments.empty() || arguments.size() > 2) {
            throw std::runtime_error("poly_integ expects polynomial and optional variable");
        }
        PolynomialData poly = build_polynomial(ctx, arguments[0]);
        *output = polynomial_to_string(polynomial_integral(poly.coefficients),
                                       poly.variable_name.empty() ? "x" : poly.variable_name);
        return true;
    }

    if (command == "poly_compose") {
        if (arguments.size() != 2) {
            throw std::runtime_error("poly_compose expects exactly two arguments");
        }
        PolynomialData lhs = build_polynomial(ctx, arguments[0]);
        PolynomialData rhs = build_polynomial(ctx, arguments[1]);
        *output = polynomial_to_string(polynomial_compose(lhs.coefficients, rhs.coefficients),
                                       lhs.variable_name.empty() ? rhs.variable_name : lhs.variable_name);
        return true;
    }

    if (command == "poly_gcd") {
        if (arguments.size() != 2) {
            throw std::runtime_error("poly_gcd expects exactly two arguments");
        }
        PolynomialData lhs = build_polynomial(ctx, arguments[0]);
        PolynomialData rhs = build_polynomial(ctx, arguments[1]);
        *output = polynomial_to_string(polynomial_gcd(lhs.coefficients, rhs.coefficients),
                                       lhs.variable_name.empty() ? rhs.variable_name : lhs.variable_name);
        return true;
    }

    if (command == "poly_fit" || command == "polynomial_fit") {
        if (arguments.size() != 3) {
            throw std::runtime_error("poly_fit expects x_samples, y_samples, and degree");
        }
        std::vector<Scalar> x_samples;
        std::vector<Scalar> y_samples;
        if (!try_parse_vector_coefficients(arguments[0], &x_samples)) {
            throw std::runtime_error("poly_fit expects x_samples vector");
        }
        if (!try_parse_vector_coefficients(arguments[1], &y_samples)) {
            throw std::runtime_error("poly_fit expects y_samples vector");
        }
        int degree = std::stoi(trim_copy(arguments[2]));
        std::vector<Scalar> fit_coeffs = polynomial_fit(x_samples, y_samples, degree);
        *output = polynomial_to_string(fit_coeffs, "x");
        return true;
    }

    // poly_add, poly_sub, poly_mul, poly_div
    if (arguments.size() != 2) {
        throw std::runtime_error("polynomial operations expect exactly two arguments");
    }

    PolynomialData lhs = build_polynomial(ctx, arguments[0]);
    PolynomialData rhs = build_polynomial(ctx, arguments[1]);

    if (command == "poly_add") {
        *output = poly_add(lhs, rhs);
    } else if (command == "poly_sub") {
        *output = poly_sub(lhs, rhs);
    } else if (command == "poly_mul") {
        *output = poly_mul(lhs, rhs);
    } else if (command == "poly_div") {
        *output = poly_div(lhs, rhs);
    } else {
        return false;
    }
    return true;
}

// ============================================================================
// PolynomialModule 类实现
// ============================================================================

/**
 * @brief 执行多项式命令
 *
 * 实现 CalculatorModule 接口，将命令分发给相应的处理函数。
 *
 * @param command 命令名称
 * @param args 命令参数列表
 * @param services 核心服务接口，提供符号解析等功能
 * @return 命令执行结果字符串
 * @throw std::runtime_error 当命令未知或执行失败时抛出
 */
std::string PolynomialModule::execute_args(const std::string& command,
                                          const std::vector<std::string>& args,
                                          ServiceLocator& locator) {
    auto services = locator.resolve<CoreServices>();
    PolynomialContext ctx;
    ctx.functions = nullptr;
    ctx.resolve_symbolic = [&](const std::string& name, std::string* var) {
        SymbolicExpression expr;
        locator.resolve<CoreServices>()->symbolic.resolve_symbolic(name, false, var, &expr);
        return expr;
    };

    std::string output;
    if (handle_polynomial_command(ctx, command, args, &output)) {
        return output;
    }
    throw std::runtime_error("Unknown polynomial command: " + command);
}

/**
 * @brief 获取帮助信息片段
 *
 * 返回多项式命令的帮助文本，用于计算器的帮助系统。
 *
 * @param topic 帮助主题（如 "functions"）
 * @return 对应主题的帮助文本，如果主题不匹配则返回空字符串
 */
std::string PolynomialModule::get_help_snippet(const std::string& topic) const {
    if (topic == "functions") {
        return "Polynomials:\n"
               "  poly_add(p, q), poly_sub(p, q), poly_mul(p, q), poly_div(p, q)\n"
               "  poly_deriv(p), poly_integ(p), poly_compose(p, q), poly_gcd(p, q)\n"
               "  poly_eval(p, x), roots(p)   Real and complex roots";
    }
    return "";
}

}  // namespace polynomial_ops

REGISTER_CALCULATOR_MODULE(polynomial_ops::PolynomialModule)

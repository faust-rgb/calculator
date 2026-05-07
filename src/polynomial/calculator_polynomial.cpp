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

#include "calculator_polynomial.h"

#include "parser/unified_expression_parser.h"
#include "matrix_internal.h"
#include "polynomial.h"
#include "mymath.h"
#include "math/helpers/integer_helpers.h"
#include "parser/command_parser.h"
#include "core/string_utils.h"

#include <sstream>

namespace polynomial_ops {

namespace {

// ============================================================================
// 内部辅助函数
// ============================================================================

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
 *
 * 支持的输入形式：
 * 1. 直接的符号表达式（如 "x^2 + 2*x + 1"）
 * 2. 嵌套的多项式操作：
 *    - poly_add(p, q): 多项式加法
 *    - poly_sub(p, q): 多项式减法
 *    - poly_mul(p, q): 多项式乘法
 *    - poly_div(p, q): 多项式除法（要求余数为零）
 *
 * 注意：嵌套的 poly_div 要求余数必须为零，否则抛出异常。
 */
void build_polynomial_recursive(
    const PolynomialContext& ctx,
    const std::string& argument,
    std::string* variable_name,
    std::vector<long double>* coefficients) {

    const std::string trimmed_argument = trim_copy(argument);
    CommandASTNode ast = parse_command(trimmed_argument);

    // 检查是否为嵌套的多项式操作
    if (ast.kind == CommandKind::kFunctionCall) {
        const auto* call = ast.as_function_call();
        if (call->name == "poly_add" || call->name == "poly_sub" ||
            call->name == "poly_mul" || call->name == "poly_div") {

            if (call->arguments.size() != 2) {
                throw std::runtime_error(
                    "polynomial operations expect exactly two arguments");
            }

            std::string lhs_variable;
            std::string rhs_variable;
            std::vector<long double> lhs_coefficients;
            std::vector<long double> rhs_coefficients;
            build_polynomial_recursive(ctx, std::string(call->arguments[0].text),
                                       &lhs_variable, &lhs_coefficients);
            build_polynomial_recursive(ctx, std::string(call->arguments[1].text),
                                       &rhs_variable, &rhs_coefficients);

            *variable_name = lhs_variable;
            if (call->name == "poly_add") {
                *coefficients = polynomial_add(lhs_coefficients, rhs_coefficients);
                return;
            }
            if (call->name == "poly_sub") {
                *coefficients = polynomial_subtract(lhs_coefficients, rhs_coefficients);
                return;
            }
            if (call->name == "poly_div") {
                const PolynomialDivisionResult division =
                    polynomial_divide(lhs_coefficients, rhs_coefficients);
                bool zero_remainder = true;
                for (long double coefficient : division.remainder) {
                    if (!mymath::is_near_zero(coefficient, 1e-10)) {
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
        lhs.variable_name);
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
        lhs.variable_name);
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
        lhs.variable_name);
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
    return "quotient: " +
           polynomial_to_string(division.quotient, lhs.variable_name) +
           ", remainder: " +
           polynomial_to_string(division.remainder, lhs.variable_name);
}

/**
 * @brief 计算多项式的所有根
 *
 * 计算多项式的全部实根和复根，返回格式化的字符串。
 *
 * @param poly 输入多项式
 * @return 格式化的根字符串，多个根用逗号分隔
 *
 * 算法：
 * 1. 调用 polynomial_complex_roots 获取所有复根
 * 2. 过滤重复的根（在容差范围内）
 * 3. 对接近整数的实部和虚部进行取整
 * 4. 实根直接输出数值，复根以 a + bi 形式输出
 */
std::string roots(const PolynomialData& poly) {
    const std::vector<mymath::complex<long double>> roots =
        polynomial_complex_roots(poly.coefficients);
    if (roots.empty()) {
        return "No roots.";
    }

    std::ostringstream out;
    bool wrote_root = false;
    mymath::complex<long double> previous_root(0.0L, 0.0L);
    for (std::size_t i = 0; i < roots.size(); ++i) {
        if (wrote_root &&
            mymath::abs(roots[i].real() - previous_root.real()) <= 1e-7 &&
            mymath::abs(roots[i].imag() - previous_root.imag()) <= 1e-7) {
            continue;
        }
        if (wrote_root) {
            out << ", ";
        }
        long double real = roots[i].real();
        long double imag = roots[i].imag();
        if (is_integer_double(real, 1e-6)) {
            real = static_cast<long double>(round_to_long_long(real));
        }
        if (is_integer_double(imag, 1e-6)) {
            imag = static_cast<long double>(round_to_long_long(imag));
        }
        if (mymath::is_near_zero(imag, 1e-8)) {
            out << format_symbolic_scalar(real);
        } else {
            out << matrix::internal::format_complex<long double>({real, imag});
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
 *
 * 支持的命令：poly_add, poly_sub, poly_mul, poly_div, roots
 */
bool is_polynomial_command(const std::string& command) {
    return command == "poly_add" ||
           command == "poly_sub" ||
           command == "poly_mul" ||
           command == "poly_div" ||
           command == "roots";
}

/**
 * @brief 处理多项式命令
 *
 * 解析命令参数并执行相应的多项式操作。支持嵌套的多项式表达式。
 *
 * @param ctx 多项式构建上下文，提供符号解析功能
 * @param command 命令名称（poly_add, poly_sub, poly_mul, poly_div, roots）
 * @param inside 括号内的参数字符串
 * @param output 输出结果字符串指针
 * @return 是否成功处理命令
 *
 * @throw std::runtime_error 当参数数量或格式错误时抛出
 *
 * 支持的命令：
 * - roots(p): 计算多项式 p 的所有根（实根和复根）
 * - poly_add(p, q): 多项式加法 p + q
 * - poly_sub(p, q): 多项式减法 p - q
 * - poly_mul(p, q): 多项式乘法 p * q
 * - poly_div(p, q): 多项式除法 p / q，返回商和余数
 *
 * @code
 * // 示例用法
 * PolynomialContext ctx = {...};
 * std::string result;
 * handle_polynomial_command(ctx, "roots", "x^2 - 1", &result);
 * // result = "-1, 1"
 * handle_polynomial_command(ctx, "poly_add", "x + 1, x - 1", &result);
 * // result = "2 * x"
 * @endcode
 */
bool handle_polynomial_command(const PolynomialContext& ctx,
                               const std::string& command,
                               const std::string& inside,
                               std::string* output) {
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (command == "roots") {
        if (arguments.size() != 1) {
            throw std::runtime_error("roots expects exactly one argument");
        }
        PolynomialData poly = build_polynomial(ctx, arguments[0]);
        *output = roots(poly);
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
                                          const CoreServices& services) {
    PolynomialContext ctx;
    ctx.functions = nullptr;
    ctx.resolve_symbolic = [&](const std::string& name, std::string* var) {
        SymbolicExpression expr;
        services.symbolic.resolve_symbolic(name, false, var, &expr);
        return expr;
    };

    std::string inside;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i != 0) inside += ", ";
        inside += args[i];
    }

    std::string output;
    if (handle_polynomial_command(ctx, command, inside, &output)) {
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
               "  roots(p)            Real and complex roots";
    }
    return "";
}

}  // namespace polynomial_ops

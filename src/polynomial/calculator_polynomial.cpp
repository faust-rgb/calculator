// ============================================================================
// 多项式操作命令实现
// ============================================================================
//
// 本文件实现了计算器中的多项式操作命令，包括：
// - 多项式四则运算（加、减、乘、除）
// - 多项式求根（实根和复根）
// - 嵌套多项式表达式的解析和计算

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
// 多项式构建辅助函数
// ============================================================================

/**
 * @brief 递归构建多项式系数
 * @param ctx 多项式构建上下文
 * @param argument 参数字符串（函数名或嵌套调用）
 * @param variable_name 输出：多项式变量名
 * @param coefficients 输出：多项式系数（低次到高次）
 *
 * 支持的输入形式：
 * 1. 直接的符号表达式（如 "x^2 + 2*x + 1"）
 * 2. 嵌套的多项式操作（如 "poly_add(p, q)"）
 *
 * 对于除法，要求余数为零。
 */
void build_polynomial_recursive(
    const PolynomialContext& ctx,
    const std::string& argument,
    std::string* variable_name,
    std::vector<double>* coefficients) {

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
            std::vector<double> lhs_coefficients;
            std::vector<double> rhs_coefficients;
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
                for (double coefficient : division.remainder) {
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
 * @param ctx 多项式构建上下文
 * @param argument 参数字符串
 * @return 多项式数据（包含变量名和系数）
 */
PolynomialData build_polynomial(const PolynomialContext& ctx,
                                const std::string& argument) {
    PolynomialData result;
    build_polynomial_recursive(ctx, argument,
                               &result.variable_name, &result.coefficients);
    return result;
}

// ============================================================================
// 多项式运算函数
// ============================================================================

/**
 * @brief 多项式加法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 和的多项式字符串
 */
std::string poly_add(const PolynomialData& lhs, const PolynomialData& rhs) {
    return polynomial_to_string(
        polynomial_add(lhs.coefficients, rhs.coefficients),
        lhs.variable_name);
}

/**
 * @brief 多项式减法
 * @param lhs 左操作数（被减数）
 * @param rhs 右操作数（减数）
 * @return 差的多项式字符串
 */
std::string poly_sub(const PolynomialData& lhs, const PolynomialData& rhs) {
    return polynomial_to_string(
        polynomial_subtract(lhs.coefficients, rhs.coefficients),
        lhs.variable_name);
}

/**
 * @brief 多项式乘法
 * @param lhs 左操作数
 * @param rhs 右操作数
 * @return 积的多项式字符串
 */
std::string poly_mul(const PolynomialData& lhs, const PolynomialData& rhs) {
    return polynomial_to_string(
        polynomial_multiply(lhs.coefficients, rhs.coefficients),
        lhs.variable_name);
}

/**
 * @brief 多项式除法
 * @param lhs 被除数
 * @param rhs 除数
 * @return 包含商和余数的结果字符串
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
 * @param poly 多项式数据
 * @return 格式化的根字符串（实根或复根）
 *
 * 显示所有不同的根，复根以 a + bi 形式输出。
 * 对接近整数的实部和虚部进行取整处理。
 */
std::string roots(const PolynomialData& poly) {
    const std::vector<mymath::complex<double>> roots =
        polynomial_complex_roots(poly.coefficients);
    if (roots.empty()) {
        return "No roots.";
    }

    std::ostringstream out;
    bool wrote_root = false;
    mymath::complex<double> previous_root(0.0, 0.0);
    for (std::size_t i = 0; i < roots.size(); ++i) {
        if (wrote_root &&
            mymath::abs(roots[i].real() - previous_root.real()) <= 1e-7 &&
            mymath::abs(roots[i].imag() - previous_root.imag()) <= 1e-7) {
            continue;
        }
        if (wrote_root) {
            out << ", ";
        }
        double real = roots[i].real();
        double imag = roots[i].imag();
        if (is_integer_double(real, 1e-6)) {
            real = static_cast<double>(round_to_long_long(real));
        }
        if (is_integer_double(imag, 1e-6)) {
            imag = static_cast<double>(round_to_long_long(imag));
        }
        if (mymath::is_near_zero(imag, 1e-8)) {
            out << format_symbolic_scalar(real);
        } else {
            out << matrix::internal::format_complex<double>({real, imag});
        }
        previous_root = {real, imag};
        wrote_root = true;
    }
    return out.str();
}

// ============================================================================
// 命令处理函数
// ============================================================================

/**
 * @brief 检查是否为多项式命令
 * @param command 命令字符串
 * @return 如果是多項式命令返回 true
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
 * @param ctx 多项式构建上下文
 * @param command 命令名称
 * @param inside 括号内的参数字符串
 * @param output 输出结果字符串
 * @return 是否成功处理
 *
 * 支持的命令：
 * - poly_add(p, q): 多项式加法
 * - poly_sub(p, q): 多项式减法
 * - poly_mul(p, q): 多项式乘法
 * - poly_div(p, q): 多项式除法
 * - roots(p): 多项式求根
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
// PolynomialModule 实现
// ============================================================================


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

std::string PolynomialModule::get_help_snippet(const std::string& topic) const {
    if (topic == "functions") {
        return "Polynomials:\n"
               "  poly_add(p, q), poly_sub(p, q), poly_mul(p, q), poly_div(p, q)\n"
               "  roots(p)            Real and complex roots";
    }
    return "";
}

}  // namespace polynomial_ops

/**
 * @file calculator_polynomial.h
 * @brief 多项式操作命令模块
 *
 * 本文件定义了计算器中多项式运算的命令接口，作为多项式运算库 (polynomial.h)
 * 与计算器核心之间的桥梁。提供以下功能：
 * - 多项式四则运算 (poly_add, poly_sub, poly_mul, poly_div)
 * - 多项式求根 (roots)
 * - 多项式求值 (poly_eval)
 * - 多项式微分 (poly_deriv)
 * - 多项式积分 (poly_integ)
 * - 多项式拟合 (poly_fit)
 * - 多项式复合 (poly_compose)
 * - 多项式最大公因式 (poly_gcd)
 *
 * 支持嵌套的多项式表达式解析，如 poly_add(poly_mul(p, q), r)。
 */

#ifndef CALCULATOR_POLYNOMIAL_H
#define CALCULATOR_POLYNOMIAL_H

#include "core/types/module_types.h"

#include "symbolic_expression.h"

#include <string>
#include <vector>
#include <map>
#include <functional>

#include "module/calculator_module.h"

/**
 * @namespace polynomial_ops
 * @brief 多项式操作命令命名空间
 *
 * 封装所有多项式相关的命令处理逻辑，与计算器核心解耦。
 */
class ServiceLocator;

namespace polynomial_ops {

/**
 * @class PolynomialModule
 * @brief 提供多项式操作的计算器模块
 *
 * 该模块继承自 CalculatorModule，注册并处理各种多项式命令：
 * - poly_add: 多项式加法
 * - poly_sub: 多项式减法
 * - poly_mul: 多项式乘法
 * - poly_div: 多项式除法（返回商和余数）
 * - roots: 多项式求根（实根和复根）
 * - poly_eval: 多项式求值
 * - poly_deriv: 多项式求导
 * - poly_integ: 多项式积分
 * - poly_fit: 多项式拟合
 * - poly_compose: 多项式复合
 * - poly_gcd: 多项式最大公因式
 */
class PolynomialModule : public CalculatorModule {
public:
    /**
     * @brief 获取模块名称
     * @return 模块名称字符串 "Polynomial"
     */
    std::string name() const override { return "Polynomial"; }

    /**
     * @brief 获取模块提供的所有命令
     * @return 命令名称列表
     */
    std::vector<std::string> get_commands() const override {
        return {"poly_add", "poly_sub", "poly_mul", "poly_div", "roots",
                "poly_eval", "poly_deriv", "poly_integ", "poly_fit", "poly_compose", "poly_gcd"};
    }

    /**
     * @brief 执行多项式命令
     * @param command 命令名称
     * @param args 命令参数列表
     * @param locator 服务定位器
     * @return 命令执行结果字符串
     */
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ServiceLocator& locator) override;


    /**
     * @brief 获取帮助信息片段
     * @param topic 帮助主题
     * @return 对应主题的帮助文本
     */
    std::string get_help_snippet(const std::string& topic) const override;
};

// ============================================================================
// 多项式数据结构
// ============================================================================

/**
 * @struct PolynomialData
 * @brief 多项式数据结构
 *
 * 存储多项式的完整信息，包括变量名和系数向量。
 * 系数按低次到高次的顺序存储，例如 x^2 + 2x + 1 表示为 [1, 2, 1]。
 */
struct PolynomialData {
    std::string variable_name;           ///< 多项式变量名（如 "x"）
    std::vector<Scalar> coefficients; ///< 多项式系数（低次到高次）
};

/**
 * @struct PolynomialContext
 * @brief 多项式构建上下文
 *
 * 封装构建多项式所需的外部依赖，避免直接依赖 Calculator 类。
 * 用于解析符号表达式和访问自定义函数。
 */
struct PolynomialContext {
    const std::map<std::string, CustomFunction>* functions; ///< 自定义函数映射表
    /**
     * @brief 符号表达式解析函数
     * @param name 表达式字符串
     * @param error_output 错误信息输出
     * @return 解析后的符号表达式
     */
    std::function<SymbolicExpression(const std::string&, std::string*)> resolve_symbolic;
};

/**
 * @brief 从参数构建多项式系数
 *
 * 支持嵌套的多项式操作（如 poly_add(poly_mul(p, q), r)）。
 * 是多项式命令的核心解析函数。
 *
 * @param ctx 多项式构建上下文，提供符号解析等功能
 * @param argument 参数字符串（函数名、表达式或嵌套调用）
 * @return 多项式数据（包含变量名和系数向量）
 * @throw std::runtime_error 当参数无法解析为多项式时抛出
 */
PolynomialData build_polynomial(const PolynomialContext& ctx, const std::string& argument);

// ============================================================================
// 多项式运算函数
// ============================================================================

/**
 * @brief 多项式加法运算
 * @param lhs 左操作数（加数）
 * @param rhs 右操作数（加数）
 * @return 和的多项式字符串表示
 */
std::string poly_add(const PolynomialData& lhs, const PolynomialData& rhs);

/**
 * @brief 多项式减法运算
 * @param lhs 左操作数（被减数）
 * @param rhs 右操作数（减数）
 * @return 差的多项式字符串表示
 */
std::string poly_sub(const PolynomialData& lhs, const PolynomialData& rhs);

/**
 * @brief 多项式乘法运算
 * @param lhs 左操作数（乘数）
 * @param rhs 右操作数（乘数）
 * @return 积的多项式字符串表示
 */
std::string poly_mul(const PolynomialData& lhs, const PolynomialData& rhs);

/**
 * @brief 多项式除法运算
 * @param lhs 被除数
 * @param rhs 除数
 * @return 格式化的商和余数字符串，格式为 "quotient: ..., remainder: ..."
 */
std::string poly_div(const PolynomialData& lhs, const PolynomialData& rhs);

/**
 * @brief 多项式求根
 * @param poly 输入多项式
 * @return 格式化的根字符串，包含所有实根和复根
 *
 * 复根以 a + bi 形式输出，接近整数的实部和虚部会被取整。
 */
std::string roots(const PolynomialData& poly);

// ============================================================================
// 命令处理接口
// ============================================================================

/**
 * @brief 检查命令是否为多项式命令
 * @param command 命令字符串
 * @return 如果是多项式命令（poly_add, poly_sub, poly_mul, poly_div, roots）返回 true
 */
bool is_polynomial_command(const std::string& command);

/**
 * @brief 处理多项式命令
 *
 * 解析并执行多项式命令，支持嵌套的多项式表达式。
 *
 * @param ctx 多项式构建上下文
 * @param command 命令名称（poly_add, poly_sub, poly_mul, poly_div, roots）
 * @param inside 括号内的参数字符串
 * @param output 输出结果字符串指针
 * @return 是否成功处理命令
 * @throw std::runtime_error 当参数数量或格式错误时抛出
 *
 * @code
 * // 示例用法
 * PolynomialContext ctx = {...};
 * std::string result;
 * handle_polynomial_command(ctx, "roots", "x^2 - 1", &result);
 * // result = "-1, 1"
 * @endcode
 */
bool handle_polynomial_command(const PolynomialContext& ctx,
                               const std::string& command,
                               const std::vector<std::string>& arguments,
                               std::string* output);

}  // namespace polynomial_ops

#endif  // CALCULATOR_POLYNOMIAL_H

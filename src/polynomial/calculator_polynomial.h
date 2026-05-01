// ============================================================================
// 多项式操作命令
// ============================================================================
//
// 提供多项式运算命令的计算逻辑，包括：
// - 多项式四则运算 (poly_add, poly_sub, poly_mul, poly_div)
// - 多项式求根 (roots)

#ifndef CALCULATOR_POLYNOMIAL_H
#define CALCULATOR_POLYNOMIAL_H

#include "calculator_internal_types.h"

#include "symbolic_expression.h"

#include <string>
#include <vector>
#include <map>
#include <functional>

#include "../core/calculator_module.h"

namespace polynomial_ops {

/**
 * @class PolynomialModule
 * @brief 提供多项式运算命令的模块
 */
class PolynomialModule : public CalculatorModule {
public:
    std::string name() const override { return "Polynomial"; }
    
    std::vector<std::string> get_commands() const override {
        return {"poly_add", "poly_sub", "poly_mul", "poly_div", "roots", 
                "poly_eval", "poly_deriv", "poly_integ", "poly_fit", "poly_compose", "poly_gcd"};
    }


    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;

    std::string get_help_snippet(const std::string& topic) const override;
};

// ============================================================================
// 多项式构建
// ============================================================================

/**
 * @brief 多项式构建结果
 */
struct PolynomialData {
    std::string variable_name;
    std::vector<double> coefficients;
};

/**
 * @brief 多项式构建上下文
 *
 * 封装构建多项式所需的依赖，避免直接依赖 Calculator 类。
 */
struct PolynomialContext {
    const std::map<std::string, CustomFunction>* functions;
    std::function<SymbolicExpression(const std::string&, std::string*)> resolve_symbolic;
};

/**
 * @brief 从参数构建多项式系数
 *
 * 支持嵌套的多项式操作（如 poly_add(poly_mul(p, q), r)）。
 *
 * @param ctx 多项式构建上下文
 * @param argument 参数字符串（函数名或嵌套调用）
 * @return 多项式数据（变量名和系数）
 */
PolynomialData build_polynomial(const PolynomialContext& ctx, const std::string& argument);

// ============================================================================
// 多项式运算
// ============================================================================

/**
 * @brief 多项式加法
 */
std::string poly_add(const PolynomialData& lhs, const PolynomialData& rhs);

/**
 * @brief 多项式减法
 */
std::string poly_sub(const PolynomialData& lhs, const PolynomialData& rhs);

/**
 * @brief 多项式乘法
 */
std::string poly_mul(const PolynomialData& lhs, const PolynomialData& rhs);

/**
 * @brief 多项式除法
 * @return 格式化的商和余数字符串
 */
std::string poly_div(const PolynomialData& lhs, const PolynomialData& rhs);

/**
 * @brief 多项式求根
 * @return 格式化的实根字符串
 */
std::string roots(const PolynomialData& poly);

// ============================================================================
// 命令处理
// ============================================================================

/**
 * @brief 检查是否为多项式命令
 */
bool is_polynomial_command(const std::string& command);

/**
 * @brief 处理多项式命令
 *
 * @param ctx 多项式构建上下文
 * @param command 命令名（poly_add, poly_sub, poly_mul, poly_div, roots）
 * @param inside 括号内的参数字符串
 * @param output 输出字符串
 * @return 是否成功处理
 */
bool handle_polynomial_command(const PolynomialContext& ctx,
                               const std::string& command,
                               const std::string& inside,
                               std::string* output);

}  // namespace polynomial_ops

#endif  // CALCULATOR_POLYNOMIAL_H

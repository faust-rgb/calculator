/**
 * @file symbolic_expression_integral_helpers.h
 * @brief 符号积分辅助函数
 *
 * 本文件定义了符号积分的辅助函数：
 * - 换元积分：识别并应用换元法
 * - 三角积分：处理三角函数乘积
 * - 特殊形式：处理 sec^n, csc^n 等形式
 */

#ifndef SYMBOLIC_EXPRESSION_INTEGRAL_HELPERS_H
#define SYMBOLIC_EXPRESSION_INTEGRAL_HELPERS_H

#include "symbolic/core/symbolic_expression_internal.h"
#include <string>

namespace symbolic_expression_internal {

bool try_integrate_substitution_product(const SymbolicExpression& derivative_factor,
                                        const SymbolicExpression& function_factor,
                                        const std::string& variable_name,
                                        SymbolicExpression* result);

bool try_integrate_trig_product_identity(const SymbolicExpression& left,
                                         const SymbolicExpression& right,
                                         const std::string& variable_name,
                                         SymbolicExpression* result);

bool try_integrate_sec_csc_power_product(const SymbolicExpression& left,
                                         const SymbolicExpression& right,
                                         const std::string& variable_name,
                                         SymbolicExpression* result);

bool try_integrate_by_parts(const SymbolicExpression& left,
                            const SymbolicExpression& right,
                            const std::string& variable_name,
                            SymbolicExpression* result);

bool is_pure_quadratic(const SymbolicExpression& expr,
                       const std::string& var,
                       SymbolicExpression* c_term,
                       SymbolicExpression* x2_coeff);

bool same_simplified_expression(const SymbolicExpression& a, const SymbolicExpression& b);

bool try_integrate_symbolic_two_linear_factors(const SymbolicExpression& numerator,
                                               const SymbolicExpression& denominator,
                                               const std::string& var,
                                               SymbolicExpression* result);

bool try_integrate_symbolic_repeated_linear_and_linear(const SymbolicExpression& numerator,
                                                       const SymbolicExpression& denominator,
                                                       const std::string& var,
                                                       SymbolicExpression* result);

bool try_integrate_symbolic_two_linear_times_pure_quadratic(const SymbolicExpression& numerator,
                                                            const SymbolicExpression& denominator,
                                                            const std::string& var,
                                                            SymbolicExpression* result);

bool try_integrate_symbolic_linear_times_pure_quadratic(const SymbolicExpression& numerator,
                                                        const SymbolicExpression& denominator,
                                                        const std::string& var,
                                                        SymbolicExpression* result);

bool try_integrate_symbolic_repeated_linear_times_pure_quadratic(const SymbolicExpression& numerator,
                                                                 const SymbolicExpression& denominator,
                                                                 const std::string& var,
                                                                 SymbolicExpression* result);

bool try_integrate_symbolic_repeated_pure_quadratic(const SymbolicExpression& numerator,
                                                    const SymbolicExpression& denominator,
                                                    const std::string& var,
                                                    SymbolicExpression* result);

bool try_integrate_symbolic_two_pure_quadratics(const SymbolicExpression& numerator,
                                                const SymbolicExpression& denominator,
                                                const std::string& var,
                                                SymbolicExpression* result);

bool try_integrate_polynomial_quotient(const SymbolicExpression& numerator,
                                       const SymbolicExpression& denominator,
                                       const std::string& var,
                                       SymbolicExpression* result);

bool try_integrate_weierstrass_substitution(const SymbolicExpression& expr,
                                            const std::string& var,
                                            SymbolicExpression* result);

bool try_integrate_trig_power_identity(const SymbolicExpression& expr,
                                       Scalar exponent_value,
                                       const std::string& var,
                                       SymbolicExpression* result);

bool primitive_for_outer_function(const std::string& function_name,
                                  const SymbolicExpression& argument,
                                  SymbolicExpression* primitive);

bool is_general_quadratic(const SymbolicExpression& expr,
                          const std::string& var,
                          SymbolicExpression* a,
                          SymbolicExpression* b,
                          SymbolicExpression* c);

} // namespace symbolic_expression_internal

#endif // SYMBOLIC_EXPRESSION_INTEGRAL_HELPERS_H

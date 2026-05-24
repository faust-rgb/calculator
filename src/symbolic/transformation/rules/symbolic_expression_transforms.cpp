// ============================================================================
// 积分变换实现模块
// ============================================================================
//
// 本文件为 SymbolicExpression 类提供积分变换的公共接口：
//
// 1. Fourier 变换
//    - fourier_transform: 正向变换
//    - inverse_fourier_transform: 逆向变换
//
// 2. Laplace 变换
//    - laplace_transform: 正向变换
//    - inverse_laplace_transform: 逆向变换
//
// 3. Z 变换
//    - z_transform: 正向变换
//    - inverse_z_transform: 逆向变换
//
// 具体实现在 transforms.cpp 中，本文件仅提供包装接口。
// ============================================================================

#include "symbolic_expression_internal.h"

using namespace symbolic_expression_internal;

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

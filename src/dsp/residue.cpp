/**
 * @file residue.cpp
 * @brief 留数计算实现
 *
 * 本文件实现了有理函数留数的计算功能：
 * - 解析有理函数表达式
 * - 提取分子分母多项式系数
 * - 计算指定极点处的留数
 *
 * 留数计算在复分析和信号处理中有重要应用，
 * 如部分分式分解和拉普拉斯逆变换。
 */

#include "residue.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "parser/grammars/unified_expression_parser.h"
#include "core/services/string_utils.h"
#include "core/types/module_types.h"
#include "math/mymath.h"
#include "math/types/complex.h"
#include "symbolic_expression.h"
#include "symbolic/public/symbolic_node_types.h"
#include "matrix.h"
#include "polynomial.h"
#include "calculator_exceptions.h"
#include "math/base/precision_constants.h"
#include "module/calculator_module.h"

namespace dsp_ops {

std::string handle_residue_command(const std::string& command,
                                   const std::vector<std::string>& arguments,
                                   ServiceLocator& locator) {
    (void)command;
    auto services = locator.resolve<CoreServices>();
    if (arguments.size() != 3) {
        throw DimensionError("residue(expression, variable, point) expects 3 arguments");
    }

    const std::string variable_name = trim_copy(arguments[1]);
    if (!is_identifier_text(variable_name)) {
        throw SyntaxError("residue variable must be an identifier");
    }

    const SymbolicExpression expression =
        SymbolicExpression::parse(
            trim_copy(services->symbolic.expand_inline(arguments[0])))
            .simplify();
    SymbolicExpression numerator = expression;
    SymbolicExpression denominator = SymbolicExpression::number(1.0L);
    if (expression.node_type() == NodeType::kDivide) {
        numerator = expression.left_child().simplify();
        denominator = expression.right_child().simplify();
    }

    std::vector<Scalar> numerator_coefficients;
    std::vector<Scalar> denominator_coefficients;
    if (!numerator.polynomial_coefficients(variable_name,
                                           &numerator_coefficients) ||
        !denominator.polynomial_coefficients(variable_name,
                                             &denominator_coefficients)) {
        throw MathError("residue currently supports rational polynomial expressions");
    }

    StoredValue point_value = services->evaluation.evaluate_value(arguments[2], false);

    mymath::complex<Scalar> point(point_value.get_decimal(), 0.0L);
    if (point_value.is_matrix && point_value.matrix_ptr) {
        const matrix::Matrix& point_matrix = *point_value.matrix_ptr;
        if (!point_matrix.is_vector() ||
            point_matrix.rows * point_matrix.cols != 2) {
            throw DimensionError("residue point must be scalar or complex(real, imag)");
        }
        const Scalar real = point_matrix.rows == 1 ? point_matrix.at(0, 0)
                                                   : point_matrix.at(0, 0);
        const Scalar imag = point_matrix.rows == 1 ? point_matrix.at(0, 1)
                                                   : point_matrix.at(1, 0);
        point = mymath::complex<Scalar>(real, imag);
    } else if (point_value.is_complex) {
        point = point_value.complex;
    }

    auto evaluate_polynomial_complex =
        [](const std::vector<Scalar>& coefficients,
           mymath::complex<Scalar> value) {
            mymath::complex<Scalar> result(0.0L, 0.0L);
            for (std::size_t i = coefficients.size(); i > 0; --i) {
                result = result * value + coefficients[i - 1];
            }
            return result;
        };

    // 计算分母在 point 处的阶数 m (泰勒级数首个非零项)
    std::vector<std::vector<Scalar>> d_derivs;
    d_derivs.push_back(denominator_coefficients);
    while (d_derivs.back().size() > 1) {
        d_derivs.push_back(polynomial_derivative(d_derivs.back()));
    }

    std::vector<mymath::complex<Scalar>> q_coeffs; // Q(w + point) 的展开系数
    Scalar fact = Scalar(1.0L);
    for (std::size_t k = 0; k < d_derivs.size(); ++k) {
        if (k > 0) fact *= Scalar(static_cast<long long>(k));
        mymath::complex<Scalar> val = evaluate_polynomial_complex(d_derivs[k], point) / fact;
        q_coeffs.push_back(val);
    }

    std::size_t m = 0;
    while (m < q_coeffs.size() && mymath::abs(q_coeffs[m]) <= Scalar(1e-9L)) {
        ++m;
    }

    // 如果分母在 point 处不为 0 (m == 0)，则不是极点，留数为 0
    if (m == 0) {
        return matrix::Matrix::vector({0.0L, 0.0L}).to_string();
    }

    // 计算分子在 point 处的泰勒展开系数 P(w + point)
    std::vector<std::vector<Scalar>> n_derivs;
    n_derivs.push_back(numerator_coefficients);
    while (n_derivs.back().size() > 1) {
        n_derivs.push_back(polynomial_derivative(n_derivs.back()));
    }

    std::vector<mymath::complex<Scalar>> p_coeffs;
    fact = Scalar(1.0L);
    for (std::size_t k = 0; k < m; ++k) {
        if (k > 0) fact *= Scalar(static_cast<long long>(k));
        mymath::complex<Scalar> val(0.0L, 0.0L);
        if (k < n_derivs.size()) {
            val = evaluate_polynomial_complex(n_derivs[k], point) / fact;
        }
        p_coeffs.push_back(val);
    }

    // 幂级数长除法：P(w) / \tilde{Q}(w)，求到 w^{m-1} 项的系数
    // \tilde{Q}(w) = q_coeffs[m] + q_coeffs[m+1]*w + ...
    mymath::complex<Scalar> q0 = q_coeffs[m];
    if (mymath::abs(q0) <= Scalar(1e-18L)) {
        throw MathError("Failed to evaluate denominator derivative at pole");
    }

    std::vector<mymath::complex<Scalar>> c_coeffs(m, mymath::complex<Scalar>(0.0L, 0.0L));
    for (std::size_t k = 0; k < m; ++k) {
        mymath::complex<Scalar> sum = p_coeffs[k];
        for (std::size_t j = 1; j <= k; ++j) {
            if (m + j < q_coeffs.size()) {
                sum = sum - q_coeffs[m + j] * c_coeffs[k - j];
            }
        }
        c_coeffs[k] = sum / q0;
    }

    const mymath::complex<Scalar> residue = c_coeffs[m - 1];

    // 规范化结果
    auto normalize = [](Scalar x) -> Scalar {
        if (!mymath::isfinite(x)) return x;
        if (mymath::is_near_zero(x, precision::epsilon<Scalar>() * Scalar(100))) return Scalar(0);
        return x;
    };

    return matrix::Matrix::vector(
                  {normalize(residue.real()),
                   normalize(residue.imag())})
                  .to_string();
}

} // namespace dsp_ops

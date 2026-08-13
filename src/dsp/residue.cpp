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

    const std::vector<Scalar> denominator_derivative =
        polynomial_derivative(denominator_coefficients);
    const mymath::complex<Scalar> denominator_value =
        evaluate_polynomial_complex(denominator_coefficients, point);
    if (mymath::abs(denominator_value) > 1e-8) {
        return matrix::Matrix::vector({0.0L, 0.0L}).to_string();
    }
    const mymath::complex<Scalar> denominator_prime =
        evaluate_polynomial_complex(denominator_derivative, point);
    if (mymath::abs(denominator_prime) <= 1e-10) {
        throw MathError("residue currently supports only simple poles");
    }

    const mymath::complex<Scalar> residue =
        evaluate_polynomial_complex(numerator_coefficients, point) /
        denominator_prime;

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

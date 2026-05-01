#include "residue.h"
#include "calculator_internal_types.h"
#include "../math/mymath.h"
#include "mymath_complex.h"
#include "symbolic_expression.h"
#include "symbolic_expression_internal.h"
#include "matrix.h"
#include "polynomial.h"
#include "calculator_exceptions.h"
#include "../core/calculator_module.h"

namespace dsp_ops {

std::string handle_residue_command(const std::string& command,
                                   const std::string& inside,
                                   const CoreServices& svc) {
    std::vector<std::string> arguments = split_top_level_arguments(inside);
    if (arguments.size() != 3) {
        throw DimensionError("residue(expression, variable, point) expects 3 arguments");
    }

    const std::string variable_name = trim_copy(arguments[1]);
    if (!is_identifier_text(variable_name)) {
        throw SyntaxError("residue variable must be an identifier");
    }

    const SymbolicExpression expression =
        SymbolicExpression::parse(
            trim_copy(svc.symbolic.expand_inline(arguments[0])))
            .simplify();
    SymbolicExpression numerator = expression;
    SymbolicExpression denominator = SymbolicExpression::number(1.0);
    if (expression.node_->type == NodeType::kDivide) {
        numerator = SymbolicExpression(expression.node_->left).simplify();
        denominator = SymbolicExpression(expression.node_->right).simplify();
    }

    std::vector<double> numerator_coefficients;
    std::vector<double> denominator_coefficients;
    if (!numerator.polynomial_coefficients(variable_name,
                                           &numerator_coefficients) ||
        !denominator.polynomial_coefficients(variable_name,
                                             &denominator_coefficients)) {
        throw MathError("residue currently supports rational polynomial expressions");
    }

    StoredValue point_value = svc.evaluation.evaluate_value(arguments[2], false);

    mymath::complex<double> point(point_value.exact
                                   ? rational_to_double(point_value.rational)
                                   : point_value.decimal,
                               0.0);
    if (point_value.is_matrix) {
        const matrix::Matrix& point_matrix = point_value.matrix;
        if (!point_matrix.is_vector() ||
            point_matrix.rows * point_matrix.cols != 2) {
            throw DimensionError("residue point must be scalar or complex(real, imag)");
        }
        const double real = point_matrix.rows == 1 ? point_matrix.at(0, 0)
                                                   : point_matrix.at(0, 0);
        const double imag = point_matrix.rows == 1 ? point_matrix.at(0, 1)
                                                   : point_matrix.at(1, 0);
        point = {real, imag};
    } else if (point_value.is_complex) {
        point = {point_value.complex.real, point_value.complex.imag};
    }

    auto evaluate_polynomial_complex =
        [](const std::vector<double>& coefficients,
           mymath::complex<double> value) {
            mymath::complex<double> result(0.0, 0.0);
            for (std::size_t i = coefficients.size(); i > 0; --i) {
                result = result * value + coefficients[i - 1];
            }
            return result;
        };

    const std::vector<double> denominator_derivative =
        polynomial_derivative(denominator_coefficients);
    const mymath::complex<double> denominator_value =
        evaluate_polynomial_complex(denominator_coefficients, point);
    if (mymath::abs(denominator_value) > 1e-8) {
        return matrix::Matrix::vector({0.0, 0.0}).to_string();
    }
    const mymath::complex<double> denominator_prime =
        evaluate_polynomial_complex(denominator_derivative, point);
    if (mymath::abs(denominator_prime) <= 1e-10) {
        throw MathError("residue currently supports only simple poles");
    }

    const mymath::complex<double> residue =
        evaluate_polynomial_complex(numerator_coefficients, point) /
        denominator_prime;

    // 规范化结果
    auto normalize = [](double x) -> double {
        if (!mymath::isfinite(x)) return x;
        if (mymath::is_near_zero(x, 1e-10)) return 0.0;
        return x;
    };

    return matrix::Matrix::vector(
                  {normalize(residue.real()),
                   normalize(residue.imag())})
                  .to_string();
}

} // namespace dsp_ops

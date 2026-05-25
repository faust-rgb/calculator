/**
 * @file pade_approximation.cpp
 * @brief Padé 逼近实现
 *
 * 本文件实现了 Padé 逼近算法：
 * - 有理逼近：用有理函数逼近幂级数
 * - 分子分母计算：求解 Padé 系数的线性方程组
 * - 收敛加速：Padé 逼近常比 Taylor 级数收敛更快
 *
 * Padé 逼近在物理和工程中有广泛应用。
 */

#include "analysis/series/pade_approximation.h"
#include "analysis/series/taylor_series.h"
#include "analysis/modules/series_module.h"
#include "app/default_precision.h"
#include "app/scalar_type.h"
#include "symbolic/core/symbolic_expression.h"
#include "math/mymath.h"
#include "expression_utils.h"
#include "string_utils.h"
#include "polynomial/polynomial.h"
#include <sstream>
#include <algorithm>

namespace series_ops {
namespace pade {

using internal::PoleException;

std::string format_pade_result(const std::vector<Scalar>& numerator,
                                const std::vector<Scalar>& denominator) {
    const std::string base = "x";

    Scalar scale = denominator[0];
    if (mymath::abs(scale) < Scalar(app::algebraic_tolerance())) {
        scale = Scalar(1);
    }

    std::vector<Scalar> num_normalized(numerator.size());
    std::vector<Scalar> den_normalized(denominator.size());

    for (std::size_t i = 0; i < numerator.size(); ++i) {
        num_normalized[i] = (numerator[i] / scale).to_long_double();
    }
    for (std::size_t i = 0; i < denominator.size(); ++i) {
        den_normalized[i] = (denominator[i] / scale).to_long_double();
    }

    std::string num_text = polynomial_to_string(num_normalized, base);
    std::string den_text = polynomial_to_string(den_normalized, base);

    if (den_text == "1" || den_text == "1.0L") {
        return SymbolicExpression::parse(num_text).simplify().to_string();
    }

    return SymbolicExpression::parse("(" + num_text + ") / (" + den_text + ")").simplify().to_string();
}

std::string format_simple_pade(const std::vector<Scalar>& numerator,
                                const std::vector<Scalar>& denominator) {
    return format_pade_result(numerator, denominator);
}

bool solve_tohplitz_stable(std::function<Scalar(int)> c, int n, std::vector<Scalar>& q) {
    if (n == 0) return true;

    const Scalar singular_threshold = Scalar(app::summation_tolerance());

    std::vector<Scalar> f(n + 1, Scalar(0));
    std::vector<Scalar> b(n + 1, Scalar(0));
    std::vector<Scalar> q_new(n + 1, Scalar(0));
    std::vector<Scalar> q_128(n + 1, Scalar(0));

    f[0] = Scalar(1);
    b[0] = Scalar(1);
    q_128[0] = Scalar(1);

    Scalar ef = c(1);

    if (mymath::abs(ef) < singular_threshold) {
        ef = singular_threshold;
    }

    for (int k = 0; k < n; ++k) {
        Scalar denom = ef;
        if (mymath::abs(denom) < singular_threshold) {
            denom = singular_threshold * (denom < Scalar(0) ? Scalar(-1) : Scalar(1));
        }
        Scalar kappa = Scalar(0);

        Scalar sum_f = Scalar(0), sum_b = Scalar(0);
        for (int i = 0; i <= k; ++i) {
            sum_f = sum_f + f[i] * c(k + 2 - i);
            sum_b = sum_b + b[i] * c(i + 1);
        }

        if (mymath::abs(sum_b) > singular_threshold) {
            kappa = -sum_f / sum_b;
        }

        if (mymath::abs(kappa) > Scalar(1)) {
            kappa = kappa / mymath::abs(kappa) * Scalar(0.99L);
        }

        std::vector<Scalar> f_new = f;
        for (int i = 0; i <= k; ++i) {
            f_new[i + 1] = f_new[i + 1] + kappa * b[k - i];
        }
        for (int i = 0; i <= k; ++i) {
            b[i + 1] = b[i] + kappa * f[k - i];
        }
        f = f_new;

        ef = ef * (Scalar(1) - kappa * kappa);
        if (mymath::abs(ef) < singular_threshold) {
            ef = singular_threshold;
        }

        Scalar q_sum = Scalar(0);
        for (int i = 0; i <= k; ++i) {
            q_sum = q_sum + q_128[i] * c(k + 1 - i);
        }
        Scalar delta = -q_sum / ef;

        for (int i = 0; i <= k + 1; ++i) {
            q_new[i] = q_128[i] + delta * f[i];
        }
        q_128 = q_new;
    }

    q.assign(n + 1, Scalar(0));
    for (int i = 0; i <= n; ++i) {
        q[i] = q_128[i];
    }

    return true;
}

bool solve_pade_denominator(std::function<Scalar(int)> c,
                             int numerator_degree,
                             int denominator_degree,
                             std::vector<Scalar>& q) {
    if (denominator_degree == 0) return true;

    const Scalar singular_threshold = Scalar(1e-30L);

    std::vector<std::vector<Scalar>> matrix(
        static_cast<std::size_t>(denominator_degree),
        std::vector<Scalar>(static_cast<std::size_t>(denominator_degree), Scalar(0)));
    std::vector<Scalar> rhs(static_cast<std::size_t>(denominator_degree), Scalar(0));

    for (int row = 0; row < denominator_degree; ++row) {
        const int k = numerator_degree + 1 + row;
        for (int col = 0; col < denominator_degree; ++col) {
            matrix[row][col] = c(k - (col + 1));
        }
        rhs[row] = -c(k);
    }

    std::vector<Scalar> diag_scale(denominator_degree, Scalar(1));
    for (int i = 0; i < denominator_degree; ++i) {
        Scalar max_row = Scalar(0);
        for (int j = 0; j < denominator_degree; ++j) {
            max_row = mymath::fmax(max_row, mymath::abs(matrix[i][j]));
        }
        if (max_row > Scalar(0)) {
            diag_scale[i] = Scalar(1) / max_row;
            for (int j = 0; j < denominator_degree; ++j) {
                matrix[i][j] = matrix[i][j] * diag_scale[i];
            }
            rhs[i] = rhs[i] * diag_scale[i];
        }
    }

    for (int col = 0; col < denominator_degree; ++col) {
        int pivot = col;
        Scalar pivot_abs = mymath::abs(matrix[col][col]);
        for (int row = col + 1; row < denominator_degree; ++row) {
            const Scalar candidate = mymath::abs(matrix[row][col]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot = row;
            }
        }

        if (pivot_abs < singular_threshold) {
            return solve_tohplitz_stable(c, denominator_degree, q);
        }

        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
            std::swap(rhs[pivot], rhs[col]);
            std::swap(diag_scale[pivot], diag_scale[col]);
        }

        const Scalar divisor = matrix[col][col];
        for (int c_col = col; c_col < denominator_degree; ++c_col) {
            matrix[col][c_col] = matrix[col][c_col] / divisor;
        }
        rhs[col] = rhs[col] / divisor;

        for (int row = 0; row < denominator_degree; ++row) {
            if (row == col) continue;
            const Scalar factor = matrix[row][col];
            if (mymath::abs(factor) < singular_threshold) continue;
            for (int c_col = col; c_col < denominator_degree; ++c_col) {
                matrix[row][c_col] = matrix[row][c_col] - factor * matrix[col][c_col];
            }
            rhs[row] = rhs[row] - factor * rhs[col];
        }
    }

    q.assign(static_cast<std::size_t>(denominator_degree + 1), Scalar(0));
    q[0] = Scalar(1);
    for (int i = 0; i < denominator_degree; ++i) {
        q[i + 1] = rhs[i];
    }
    return true;
}

std::string pade_from_coeffs(const std::vector<Scalar>& coefficients,
                              int numerator_degree,
                              int denominator_degree) {
    //const int total_degree = numerator_degree + denominator_degree;

    auto c = [&](int k) -> Scalar {
        if (k < 0 || k >= static_cast<int>(coefficients.size())) return Scalar(0);
        return coefficients[static_cast<std::size_t>(k)];
    };

    if (denominator_degree == 0) {
        std::vector<Scalar> num(numerator_degree + 1, 0.0L);
        for (int i = 0; i <= numerator_degree && i < static_cast<int>(coefficients.size()); ++i) {
            num[i] = coefficients[i].to_long_double();
        }
        return polynomial_to_string(num, "x");
    }

    if (numerator_degree == 0) {
        std::vector<Scalar> q_coeffs(denominator_degree + 1, Scalar(0));
        q_coeffs[0] = Scalar(1);

        std::vector<std::vector<Scalar>> H(
            static_cast<std::size_t>(denominator_degree),
            std::vector<Scalar>(static_cast<std::size_t>(denominator_degree + 1), Scalar(0)));

        for (int i = 0; i < denominator_degree; ++i) {
            for (int j = 0; j < denominator_degree; ++j) {
                H[i][j] = c(i - j + 1);
            }
            H[i][denominator_degree] = -c(i + 1);
        }

        std::vector<Scalar> q(denominator_degree + 1, Scalar(0));
        q[0] = Scalar(1);
        if (!solve_tohplitz_stable(c, denominator_degree, q)) {
            return format_simple_pade({c(0)}, std::vector<Scalar>(denominator_degree + 1, Scalar(1)));
        }

        Scalar p0 = c(0);
        for (int j = 1; j <= denominator_degree; ++j) {
            p0 = p0 + q[j] * c(-j);
        }

        std::vector<Scalar> num_vec = {p0};
        std::vector<Scalar> den_vec(denominator_degree + 1);
        for (int i = 0; i <= denominator_degree; ++i) {
            den_vec[i] = q[i];
        }
        return format_pade_result(num_vec, den_vec);
    }

    std::vector<Scalar> p_coeffs(numerator_degree + 1, Scalar(0));
    std::vector<Scalar> q_coeffs(denominator_degree + 1, Scalar(0));
    q_coeffs[0] = Scalar(1);

    if (denominator_degree > 0) {
        if (!solve_pade_denominator(c, numerator_degree, denominator_degree, q_coeffs)) {
            throw std::runtime_error("pade denominator system is singular");
        }
    }

    for (int i = 0; i <= numerator_degree; ++i) {
        Scalar sum = Scalar(0);
        for (int j = 0; j <= denominator_degree && j <= i; ++j) {
            sum = sum + q_coeffs[j] * c(i - j);
        }
        p_coeffs[i] = sum;
    }

    return format_pade_result(p_coeffs, q_coeffs);
}

std::string pade(const SeriesContext& ctx,
                 const std::string& expr,
                 Scalar center,
                 int numerator_degree,
                 int denominator_degree) {
    if (numerator_degree == 0 && denominator_degree == 0) {
        throw std::runtime_error("pade requires at least one non-zero degree");
    }

    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, true, &variable_name, &expression);

    const int total_degree = numerator_degree + denominator_degree;
    const std::vector<Scalar> coefficients = taylor::build_taylor_coefficients(
        ctx, expression, variable_name, center, total_degree);

    if (coefficients.empty() || mymath::abs(coefficients[0]) < Scalar(1e-15L)) {
        int first_nonzero = 0;
        for (int i = 0; i < static_cast<int>(coefficients.size()); ++i) {
            if (mymath::abs(coefficients[i]) >= Scalar(1e-15L)) {
                first_nonzero = i;
                break;
            }
        }
        if (first_nonzero > 0) {
            std::vector<Scalar> shifted_coeffs(coefficients.begin() + first_nonzero, coefficients.end());
            std::string inner_result = pade_from_coeffs(shifted_coeffs, numerator_degree, denominator_degree);
            const std::string base = shifted_series_base(variable_name, center);
            if (first_nonzero == 1) {
                return SymbolicExpression::parse(base + " * (" + inner_result + ")").simplify().to_string();
            } else {
                return SymbolicExpression::parse(base + "^" + std::to_string(first_nonzero) + " * (" + inner_result + ")").simplify().to_string();
            }
        }
    }

    return pade_from_coeffs(coefficients, numerator_degree, denominator_degree);
}

}  // namespace pade
}  // namespace series_ops

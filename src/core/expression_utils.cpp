// ============================================================================
// 表达式工具函数实现
// ============================================================================
//
// 提供级数格式化和数值容差计算功能。
// ============================================================================

#include "expression_utils.h"
#include "format_utils.h"
#include "math/mymath.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

// ============================================================================
// 数值容差
// ============================================================================

double root_position_tolerance(double value) {
    return 1e-10 * std::max(1.0, mymath::abs(value));
}

double root_function_tolerance(double value) {
    return 1e-10 * std::max(1.0, mymath::abs(value));
}

double root_derivative_step(double value) {
    return 1e-6 * std::max(1.0, mymath::abs(value));
}

// ============================================================================
// 级数格式化
// ============================================================================

std::string shifted_series_base(const std::string& variable_name, double center) {
    if (mymath::is_near_zero(center, 1e-12)) {
        return variable_name;
    }
    return "(" + variable_name + signed_center_text(center) + ")";
}

std::string generalized_series_to_string(const std::vector<double>& coefficients,
                                         const std::string& variable_name,
                                         double center,
                                         int denominator) {
    if (denominator <= 0) {
        throw std::runtime_error("series denominator must be positive");
    }

    const std::string base = shifted_series_base(variable_name, center);
    std::vector<std::string> terms;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const double coefficient = coefficients[i];
        if (mymath::is_near_zero(coefficient, 1e-12)) {
            continue;
        }
        const std::string factor =
            power_term(base, static_cast<int>(i), denominator);
        terms.push_back(format_term(coefficient, factor));
    }

    if (terms.empty()) {
        return "0";
    }

    std::ostringstream out;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        if (i == 0) {
            out << terms[i];
        } else if (!terms[i].empty() && terms[i][0] == '-') {
            out << " - " << terms[i].substr(1);
        } else {
            out << " + " << terms[i];
        }
    }
    std::string result = out.str();
    if (result.rfind("1 * ", 0) == 0) {
        result.erase(0, 4);
    }
    return result;
}

std::string taylor_series_to_string(const std::vector<double>& coefficients,
                                    const std::string& variable_name,
                                    double center) {
    return generalized_series_to_string(coefficients, variable_name, center, 1);
}

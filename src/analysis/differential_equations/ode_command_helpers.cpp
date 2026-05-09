// ============================================================================
// ODE 命令辅助函数实现
// ============================================================================

#include "analysis/differential_equations/ode_command_helpers.h"
#include "analysis/modules/ode_module.h"
#include "symbolic/core/symbolic_expression_internal.h"
#include "core/common/scalar_type.h"

#include <algorithm>
#include <sstream>

namespace ode_ops {

using Scalar = mymath::Scalar;

// ============================================================================
// 导数阶数解析
// ============================================================================

int get_derivative_order(const std::string& var) {
    if (var.empty() || var[0] != 'y') return 0;
    if (var == "y") return 0;
    int order = 0;
    for (std::size_t i = 1; i < var.size(); ++i) {
        if (var[i] == '\'') {
            order++;
        } else {
            return 0;
        }
    }
    return order;
}

std::string order_to_var(int order) {
    std::string s = "y";
    for (int i = 0; i < order; ++i) s += "'";
    return s;
}

// ============================================================================
// ODE 分析
// ============================================================================

ODEInfo analyze_ode_expression(const std::string& expr_str) {
    using namespace symbolic_expression_internal;

    SymbolicExpression expr = SymbolicExpression::parse(expr_str);
    std::vector<std::string> vars = expr.identifier_variables();
    int max_order = 0;
    for (const std::string& v : vars) {
        max_order = std::max(max_order, get_derivative_order(v));
    }

    ODEInfo info;
    if (max_order <= 1 && expr_str.find("y'") == std::string::npos) {
        info.is_high_order = false;
        info.order = 1;
        info.rhs = expr;
        return info;
    }

    info.is_high_order = true;
    info.order = std::max(1, max_order);
    const std::string highest_var = order_to_var(info.order);

    SymbolicExpression coeff_A = expr.derivative(highest_var).simplify();
    if (coeff_A.is_constant(highest_var) && !expr_is_zero(coeff_A)) {
        SymbolicExpression term_B = expr.substitute(highest_var, SymbolicExpression::number(0.0L)).simplify();
        info.rhs = ((-term_B) / coeff_A).simplify();
    } else {
        throw std::runtime_error("Could not solve for highest derivative " + highest_var + ". The equation must be linear in the highest derivative.");
    }

    return info;
}

// ============================================================================
// 高阶 ODE 转换
// ============================================================================

std::string replace_derivative_notation(
    const std::string& expr_str,
    int max_order) {
    std::string result = expr_str;

    auto replace_identifier = [](std::string& s, const std::string& id, const std::string& replacement) {
        std::string res;
        std::string current_id;
        auto flush = [&]() {
            if (current_id == id) res += replacement;
            else res += current_id;
            current_id.clear();
        };
        for (char c : s) {
            if (std::isalnum(c) || c == '_' || c == '\'') {
                current_id += c;
            } else {
                flush();
                res += c;
            }
        }
        flush();
        s = res;
    };

    for (int i = max_order - 1; i >= 0; --i) {
        std::string from = "y";
        for (int j = 0; j < i; ++j) from += "'";
        std::string to = "y" + std::to_string(i + 1);
        replace_identifier(result, from, to);
    }

    return result;
}

std::string convert_high_order_to_system(
    const ODEInfo& info,
    std::vector<Scalar>& initial_state) {
    std::vector<std::string> system_exprs;
    for (int i = 1; i < info.order; ++i) {
        system_exprs.push_back("y" + std::to_string(i + 1));
    }

    std::string rhs_str = info.rhs.to_string();
    rhs_str = replace_derivative_notation(rhs_str, info.order);
    system_exprs.push_back(rhs_str);

    // Ensure initial state matches order
    if (initial_state.size() < static_cast<std::size_t>(info.order)) {
        initial_state.resize(info.order, Scalar(0));
    } else if (initial_state.size() > static_cast<std::size_t>(info.order)) {
        initial_state.resize(info.order);
    }

    std::string system_arg = "[";
    for (std::size_t i = 0; i < system_exprs.size(); ++i) {
        if (i > 0) system_arg += "; ";
        system_arg += system_exprs[i];
    }
    system_arg += "]";

    return system_arg;
}

// ============================================================================
// 矩阵转换辅助
// ============================================================================

std::vector<Scalar> matrix_to_vector_values(
    const matrix::Matrix& value,
    const std::string& context) {
    if (!value.is_vector()) {
        throw std::runtime_error(context + " expects vector arguments");
    }
    const std::size_t size = value.rows == 1 ? value.cols : value.rows;
    std::vector<Scalar> result(size, Scalar(0));
    for (std::size_t i = 0; i < size; ++i) {
        result[i] = Scalar(value.rows == 1 ? value.at(0, i) : value.at(i, 0));
    }
    return result;
}

matrix::Matrix vector_to_column_matrix(
    const ODEContext& ctx,
    const std::vector<Scalar>& values) {
    matrix::Matrix result(values.size(), 1, 0.0L);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result.at(i, 0) = ctx.normalize_result((values[i]));
    }
    return result;
}

// ============================================================================
// 存储值辅助
// ============================================================================

StoredValue make_scalar_stored(const ODEContext& ctx, Scalar value) {
    StoredValue stored;
    stored.decimal = ctx.normalize_result((value));
    return stored;
}

void append_parameter_assignments(
    const ODEContext& ctx,
    const StoredValue& parameter_value,
    std::vector<std::pair<std::string, StoredValue>>* assignments) {
    assignments->push_back({"p", parameter_value});
    if (!parameter_value.is_matrix || !parameter_value.matrix.is_vector()) {
        return;
    }

    const std::size_t size =
        parameter_value.matrix.rows == 1
            ? parameter_value.matrix.cols
            : parameter_value.matrix.rows;
    for (std::size_t i = 0; i < size; ++i) {
        const Scalar component_value =
            parameter_value.matrix.rows == 1
                ? Scalar(parameter_value.matrix.at(0, i))
                : Scalar(parameter_value.matrix.at(i, 0));
        assignments->push_back({"p" + std::to_string(i + 1),
                                make_scalar_stored(ctx, component_value)});
    }
}

}  // namespace ode_ops

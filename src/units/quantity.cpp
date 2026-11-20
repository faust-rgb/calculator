// ============================================================================
// quantity.cpp - 物理量实现
// ============================================================================

#include "units/quantity.h"
#include "core/services/format_utils.h"
#include "core/services/string_utils.h"
#include <sstream>
#include <stdexcept>
#include <cctype>

namespace units {

Quantity::Quantity(Scalar val, const std::string& unit_expr)
    : value(val), unit_symbol(trim_copy(unit_expr)) {
    if (!unit_symbol.empty()) {
        if (!UnitRegistry::instance().parse_unit_expression(unit_symbol, &dimension, &scale_to_si, &offset_to_si)) {
            throw std::runtime_error("unknown or invalid unit: " + unit_symbol);
        }
    } else {
        dimension = Dimension::Dimensionless();
        scale_to_si = 1.0L;
        offset_to_si = 0.0L;
    }
}

Quantity::Quantity(Scalar val, Dimension dim, Scalar scale, const std::string& symbol, Scalar offset)
    : value(val), unit_symbol(symbol), dimension(dim), scale_to_si(scale), offset_to_si(offset) {}

Quantity Quantity::convert_to(const std::string& target_unit) const {
    std::string trimmed = trim_copy(target_unit);
    Dimension t_dim;
    Scalar t_scale = 1.0L;
    Scalar t_offset = 0.0L;

    if (!UnitRegistry::instance().parse_unit_expression(trimmed, &t_dim, &t_scale, &t_offset)) {
        throw std::runtime_error("unknown target unit: " + target_unit);
    }

    if (dimension != t_dim) {
        throw std::runtime_error("dimension mismatch in unit conversion: " +
                                 dimension.to_string() + " cannot convert to " + t_dim.to_string());
    }

    Scalar si = to_si();
    Scalar new_val = (si - t_offset) / t_scale;
    return Quantity(new_val, t_dim, t_scale, trimmed, t_offset);
}

Quantity Quantity::operator+(const Quantity& other) const {
    if (dimension != other.dimension) {
        throw std::runtime_error("cannot add quantities with different dimensions: " +
                                 dimension.to_string() + " vs " + other.dimension.to_string());
    }
    Quantity other_conv = other.convert_to(unit_symbol);
    return Quantity(value + other_conv.value, dimension, scale_to_si, unit_symbol, offset_to_si);
}

Quantity Quantity::operator-(const Quantity& other) const {
    if (dimension != other.dimension) {
        throw std::runtime_error("cannot subtract quantities with different dimensions: " +
                                 dimension.to_string() + " vs " + other.dimension.to_string());
    }
    Quantity other_conv = other.convert_to(unit_symbol);
    return Quantity(value - other_conv.value, dimension, scale_to_si, unit_symbol, offset_to_si);
}

Quantity Quantity::operator*(const Quantity& other) const {
    Dimension new_dim = dimension * other.dimension;
    Scalar new_si = to_si() * other.to_si();

    if (new_dim.is_dimensionless()) {
        return Quantity(new_si, new_dim, 1.0L, "");
    }

    std::string simplified_unit = UnitRegistry::instance().simplify_dimension(new_dim);
    Dimension s_dim;
    Scalar s_scale = 1.0L;
    if (UnitRegistry::instance().parse_unit_expression(simplified_unit, &s_dim, &s_scale) && s_dim == new_dim) {
        return Quantity(new_si / s_scale, new_dim, s_scale, simplified_unit);
    }

    return Quantity(new_si, new_dim, 1.0L, new_dim.to_string());
}

Quantity Quantity::operator/(const Quantity& other) const {
    if (other.to_si() == Scalar(0)) {
        throw std::runtime_error("division by zero in quantity division");
    }

    Dimension new_dim = dimension / other.dimension;
    Scalar new_si = to_si() / other.to_si();

    if (new_dim.is_dimensionless()) {
        return Quantity(new_si, new_dim, 1.0L, "");
    }

    std::string simplified_unit = UnitRegistry::instance().simplify_dimension(new_dim);
    Dimension s_dim;
    Scalar s_scale = 1.0L;
    if (UnitRegistry::instance().parse_unit_expression(simplified_unit, &s_dim, &s_scale) && s_dim == new_dim) {
        return Quantity(new_si / s_scale, new_dim, s_scale, simplified_unit);
    }

    return Quantity(new_si, new_dim, 1.0L, new_dim.to_string());
}

Quantity Quantity::operator*(Scalar scalar) const {
    return Quantity(value * scalar, dimension, scale_to_si, unit_symbol, offset_to_si * scalar);
}

Quantity Quantity::operator/(Scalar scalar) const {
    if (scalar == Scalar(0)) {
        throw std::runtime_error("division by zero in quantity division");
    }
    return Quantity(value / scalar, dimension, scale_to_si, unit_symbol, offset_to_si / scalar);
}

Quantity Quantity::pow(int exponent) const {
    if (exponent == 0) {
        return Quantity(1.0L, Dimension::Dimensionless(), 1.0L, "");
    }

    Dimension new_dim = dimension.pow(exponent);
    Scalar new_val = value;
    Scalar v = 1.0L;
    if (exponent > 0) {
        for (int i = 0; i < exponent; ++i) v = v * new_val;
    } else {
        for (int i = 0; i < -exponent; ++i) v = v / new_val;
    }

    std::string new_sym = unit_symbol;
    if (!new_sym.empty() && exponent != 1) {
        new_sym = new_sym + "^" + std::to_string(exponent);
    }

    return Quantity(v, new_dim, 1.0L, new_sym);
}

std::string Quantity::to_string() const {
    std::string val_str = format_decimal(value);
    if (unit_symbol.empty()) {
        return val_str;
    }
    return val_str + " " + unit_symbol;
}

std::string Quantity::to_latex() const {
    std::string val_str = format_decimal(value);
    if (unit_symbol.empty()) {
        return val_str;
    }
    return val_str + " \\,\\mathrm{" + unit_symbol + "}";
}

bool Quantity::try_parse(const std::string& text, Quantity* out) {
    std::string trimmed = trim_copy(text);
    if (trimmed.empty()) return false;

    // 分离开头的数值部分和后续的单位部分
    std::size_t i = 0;
    if (trimmed[i] == '+' || trimmed[i] == '-') ++i;

    bool has_digits = false;
    bool has_dot = false;
    while (i < trimmed.size() && (std::isdigit(static_cast<unsigned char>(trimmed[i])) || trimmed[i] == '.')) {
        if (trimmed[i] == '.') {
            if (has_dot) break;
            has_dot = true;
        } else {
            has_digits = true;
        }
        ++i;
    }

    // 支持科学计数法 e/E
    if (has_digits && i < trimmed.size() && (trimmed[i] == 'e' || trimmed[i] == 'E')) {
        std::size_t exp_pos = i;
        ++i;
        if (i < trimmed.size() && (trimmed[i] == '+' || trimmed[i] == '-')) ++i;
        bool exp_digits = false;
        while (i < trimmed.size() && std::isdigit(static_cast<unsigned char>(trimmed[i]))) {
            exp_digits = true;
            ++i;
        }
        if (!exp_digits) {
            i = exp_pos; // 回退
        }
    }

    if (!has_digits) return false;

    std::string num_str = trimmed.substr(0, i);
    std::string unit_part = trim_copy(trimmed.substr(i));

    // 如果包含方括号包装，如 "[m/s]"，去除括号
    if (unit_part.size() >= 2 && unit_part.front() == '[' && unit_part.back() == ']') {
        unit_part = unit_part.substr(1, unit_part.size() - 2);
        unit_part = trim_copy(unit_part);
    }

    try {
        Scalar val(num_str);
        if (out) {
            *out = Quantity(val, unit_part);
        }
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace units

// ============================================================================
// quantity.h - 带单位物理量类
// ============================================================================

#ifndef UNITS_QUANTITY_H
#define UNITS_QUANTITY_H

#include "units/dimension.h"
#include "units/unit_registry.h"
#include "types/scalar_type.h"
#include <string>

namespace units {

class Quantity {
public:
    Scalar value{0.0L};
    std::string unit_symbol;
    Dimension dimension;
    Scalar scale_to_si{1.0L};
    Scalar offset_to_si{0.0L};

    Quantity() = default;
    Quantity(Scalar val, const std::string& unit_expr);
    Quantity(Scalar val, Dimension dim, Scalar scale, const std::string& symbol, Scalar offset = 0.0L);

    // 转换为 SI 基础数值
    Scalar to_si() const {
        return value * scale_to_si + offset_to_si;
    }

    // 单位转换
    Quantity convert_to(const std::string& target_unit) const;

    // 物理量四则运算
    Quantity operator+(const Quantity& other) const;
    Quantity operator-(const Quantity& other) const;
    Quantity operator*(const Quantity& other) const;
    Quantity operator/(const Quantity& other) const;
    Quantity operator*(Scalar scalar) const;
    Quantity operator/(Scalar scalar) const;
    Quantity pow(int exponent) const;

    // 格式化输出
    std::string to_string() const;
    std::string to_latex() const;

    // 解析形如 "50 km/h" 或 "9.8 m/s^2" 或 "25 degC" 的字符串
    static bool try_parse(const std::string& text, Quantity* out);
};

} // namespace units

#endif // UNITS_QUANTITY_H

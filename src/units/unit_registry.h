// ============================================================================
// unit_registry.h - 单位注册表与单位表达式解析器
// ============================================================================

#ifndef UNITS_UNIT_REGISTRY_H
#define UNITS_UNIT_REGISTRY_H

#include "units/dimension.h"
#include "types/scalar_type.h"
#include <string>
#include <vector>
#include <map>
#include <optional>

namespace units {

struct UnitInfo {
    std::string name;
    std::string symbol;
    std::string category;
    Dimension dimension;
    Scalar scale_to_si{1.0L};
    Scalar offset_to_si{0.0L};
    std::string description;
};

class UnitRegistry {
public:
    UnitRegistry();

    static UnitRegistry& instance();

    // 查找单一已知单位（如 "km", "m/s", "lbf", "degC" 等）
    const UnitInfo* find_unit(const std::string& symbol) const;

    // 解析复合单位表达式，如 "kg*m/s^2", "kW*h", "J/(mol*K)"
    bool parse_unit_expression(const std::string& expr,
                               Dimension* out_dimension,
                               Scalar* out_scale_to_si,
                               Scalar* out_offset_to_si = nullptr) const;

    // 根据量纲推荐简化的具名单位（如 [M L T^-2] -> "N", [M L^2 T^-2] -> "J"）
    std::string simplify_dimension(const Dimension& dim) const;

    // 获取所有注册单位
    const std::vector<UnitInfo>& get_all_units() const { return units_list_; }

    // 注册新单位
    void register_unit(UnitInfo unit);

private:
    void initialize_standard_units();

    std::vector<UnitInfo> units_list_;
    std::map<std::string, std::size_t> symbol_map_;
    std::map<std::string, std::string> dimension_canonical_unit_;
};

} // namespace units

#endif // UNITS_UNIT_REGISTRY_H

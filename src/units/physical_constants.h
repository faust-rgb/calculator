// ============================================================================
// physical_constants.h - 物理常数库
// ============================================================================

#ifndef UNITS_PHYSICAL_CONSTANTS_H
#define UNITS_PHYSICAL_CONSTANTS_H

#include "units/quantity.h"
#include <string>
#include <vector>

namespace units {

struct PhysicalConstant {
    std::string name;
    std::string symbol;
    Quantity quantity;
    std::string description;
};

class PhysicalConstants {
public:
    static const std::vector<PhysicalConstant>& get_all();
    static const PhysicalConstant* find(const std::string& name_or_symbol);
};

} // namespace units

#endif // UNITS_PHYSICAL_CONSTANTS_H

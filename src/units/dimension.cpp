// ============================================================================
// dimension.cpp - 物理量纲实现
// ============================================================================

#include "units/dimension.h"
#include <sstream>

namespace units {

bool Dimension::is_dimensionless() const {
    for (int8_t exp : exponents) {
        if (exp != 0) return false;
    }
    return true;
}

Dimension Dimension::operator*(const Dimension& other) const {
    Dimension result;
    for (std::size_t i = 0; i < kNumDimensions; ++i) {
        result.exponents[i] = static_cast<int8_t>(exponents[i] + other.exponents[i]);
    }
    return result;
}

Dimension Dimension::operator/(const Dimension& other) const {
    Dimension result;
    for (std::size_t i = 0; i < kNumDimensions; ++i) {
        result.exponents[i] = static_cast<int8_t>(exponents[i] - other.exponents[i]);
    }
    return result;
}

Dimension Dimension::pow(int power) const {
    Dimension result;
    for (std::size_t i = 0; i < kNumDimensions; ++i) {
        result.exponents[i] = static_cast<int8_t>(exponents[i] * power);
    }
    return result;
}

bool Dimension::operator==(const Dimension& other) const {
    return exponents == other.exponents;
}

bool Dimension::operator!=(const Dimension& other) const {
    return exponents != other.exponents;
}

std::string Dimension::to_string() const {
    if (is_dimensionless()) {
        return "[1]";
    }

    static const char* kSymbols[kNumDimensions] = {
        "L", "M", "T", "I", "Theta", "N", "J", "rad"
    };

    std::ostringstream oss;
    oss << "[";
    bool first = true;
    for (std::size_t i = 0; i < kNumDimensions; ++i) {
        if (exponents[i] != 0) {
            if (!first) oss << " ";
            first = false;
            oss << kSymbols[i];
            if (exponents[i] != 1) {
                oss << "^" << static_cast<int>(exponents[i]);
            }
        }
    }
    oss << "]";
    return oss.str();
}

} // namespace units

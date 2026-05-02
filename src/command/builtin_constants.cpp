// ============================================================================
// 内置常量实现
// ============================================================================
//
// 提供内置数学常量的查找功能。
// 从 core/utils.cpp 移出，职责更清晰。
// ============================================================================

#include "builtin_constants.h"
#include "math/mymath.h"

bool lookup_builtin_constant(const std::string& name, double* value) {
    if (name == "pi") {
        *value = mymath::kPi;
        return true;
    }
    if (name == "e") {
        *value = mymath::kE;
        return true;
    }
    if (name == "c") {
        *value = mymath::kSpeedOfLight;
        return true;
    }
    if (name == "G") {
        *value = mymath::kGravitationalConstant;
        return true;
    }
    if (name == "h") {
        *value = mymath::kPlanckConstant;
        return true;
    }
    if (name == "k") {
        *value = mymath::kBoltzmannConstant;
        return true;
    }
    if (name == "NA") {
        *value = mymath::kAvogadroNumber;
        return true;
    }
    if (name == "inf" || name == "infinity" || name == "oo") {
        *value = mymath::infinity();
        return true;
    }
    return false;
}

bool is_builtin_constant_name(const std::string& name) {
    return name == "pi" || name == "e" || name == "c" || name == "G" ||
           name == "h" || name == "k" || name == "NA" ||
           name == "inf" || name == "infinity" || name == "oo";
}

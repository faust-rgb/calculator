// ============================================================================
// 内置常量实现文件
// ============================================================================
//
// 本文件提供内置数学常量的查找功能实现。
// 支持的常量包括：
//   - 数学常量：pi, e
//   - 物理常量：c (光速), G (引力常数), h (普朗克常数), k (玻尔兹曼常数), NA (阿伏伽德罗常数)
//   - 特殊值：inf, infinity, oo (无穷大)
//
// 从 core/utils.cpp 移出，职责更清晰。
// ============================================================================

#include "builtin_constants.h"
#include "math/mymath.h"

// ============================================================================
// 常量查找函数实现
// ============================================================================

/**
 * @brief 查找内置常量的值
 * @param name 常量名称（如 "pi", "e", "c" 等）
 * @param value 输出参数，用于存储找到的常量值
 * @return 如果找到常量返回 true，否则返回 false
 */
bool lookup_builtin_constant(const std::string& name, Scalar* value) {
    // 数学常量 pi
    if (name == "pi") {
        *value = mymath::pi();
        return true;
    }
    // 数学常量 e (自然对数的底)
    if (name == "e") {
        *value = mymath::e();
        return true;
    }
    // 物理常量：光速
    if (name == "c") {
        *value = mymath::kSpeedOfLight;
        return true;
    }
    // 物理常量：万有引力常数
    if (name == "G") {
        *value = mymath::kGravitationalConstant;
        return true;
    }
    // 物理常量：普朗克常数
    if (name == "h") {
        *value = mymath::kPlanckConstant;
        return true;
    }
    // 物理常量：玻尔兹曼常数
    if (name == "k") {
        *value = mymath::kBoltzmannConstant;
        return true;
    }
    // 物理常量：阿伏伽德罗常数
    if (name == "NA") {
        *value = mymath::kAvogadroNumber;
        return true;
    }
    // 无穷大（支持多种表示方式）
    if (name == "inf" || name == "infinity" || name == "oo") {
        *value = mymath::infinity();
        return true;
    }
    return false;
}

/**
 * @brief 检查给定名称是否为内置常量
 * @param name 要检查的名称
 * @return 如果是内置常量名称返回 true，否则返回 false
 */
bool is_builtin_constant_name(const std::string& name) {
    return name == "pi" || name == "e" || name == "c" || name == "G" ||
           name == "h" || name == "k" || name == "NA" ||
           name == "inf" || name == "infinity" || name == "oo";
}

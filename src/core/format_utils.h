// ============================================================================
// 格式化工具函数
// ============================================================================
//
// 提供数值和存储值的格式化函数。
// 从 utils.h 拆分以便更好地组织代码。
// ============================================================================

#ifndef CORE_FORMAT_UTILS_H
#define CORE_FORMAT_UTILS_H

#include <string>

struct StoredValue;
struct Rational;

// ============================================================================
// 有理数近似
// ============================================================================

/**
 * @brief 尝试将 double 转换为简单的 Rational
 */
bool try_make_simple_rational(double value, int max_denominator, Rational* rational);


// ============================================================================
// 显示精度
// ============================================================================

/**
 * @brief 获取/设置显示精度
 */
int process_display_precision();
void set_process_display_precision(int precision);

/**
 * @brief 规范化显示的小数值（处理接近零和接近整数的情况）
 */
double normalize_display_decimal(double value);

// ============================================================================
// 数值格式化
// ============================================================================

/**
 * @brief 格式化小数值
 */
std::string format_decimal(double value);
std::string format_decimal(double value, int precision);

/**
 * @brief 格式化符号数值（如 pi, e）
 */
std::string format_symbolic_number(double value);
std::string format_symbolic_scalar(double value);

/**
 * @brief 格式化幂次项
 */
std::string power_term(const std::string& base, int numerator, int denominator);

/**
 * @brief 格式化带符号的中心文本
 */
std::string signed_center_text(double center);

/**
 * @brief 格式化项
 */
std::string format_term(double coefficient, const std::string& factor);

// ============================================================================
// 存储值格式化
// ============================================================================

/**
 * @brief 格式化 StoredValue 用于显示
 */
std::string format_stored_value(const StoredValue& value, bool symbolic_constants_mode);

/**
 * @brief 格式化 StoredValue 用于 print 命令
 */
std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode);

#endif // CORE_FORMAT_UTILS_H

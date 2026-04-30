// ============================================================================
// 存储值类型
// ============================================================================
//
// 计算器中存储的值，支持多种类型：
// - 标量（double 或有理数）
// - 矩阵
// - 字符串
// - 符号表达式文本
// - 精确小数文本

#ifndef TYPES_STORED_VALUE_H
#define TYPES_STORED_VALUE_H

#include "rational.h"
#include "matrix.h"

#include <string>

/**
 * @struct StoredValue
 * @brief 计算器中存储的值
 *
 * 通过标志位区分值类型，实现统一的存储接口。
 */
struct StoredValue {
    bool is_matrix = false;              ///< 是否为矩阵
    bool is_complex = false;             ///< 是否为复数标量
    bool is_string = false;              ///< 是否为字符串
    bool has_symbolic_text = false;      ///< 是否有符号表达式文本
    bool has_precise_decimal_text = false; ///< 是否有精确小数文本
    bool exact = false;                  ///< 是否在精确模式下创建

    Rational rational;                   ///< 有理数值
    double decimal = 0.0;                ///< 浮点数值
    matrix::ComplexNumber complex;       ///< 复数值
    std::string string_value;            ///< 字符串值
    std::string symbolic_text;           ///< 符号表达式文本
    std::string precise_decimal_text;    ///< 精确小数文本
    matrix::Matrix matrix;               ///< 矩阵值
};

/**
 * @brief 获取存储值的精确小数文本表示
 * @param value 存储值
 * @return 精确小数文本
 */
std::string stored_value_precise_decimal_text(const StoredValue& value);

#endif // TYPES_STORED_VALUE_H
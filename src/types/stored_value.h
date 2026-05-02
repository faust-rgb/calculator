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
 * 支持延迟符号计算：符号文本仅在需要时才计算。
 */
struct StoredValue {
    bool is_matrix = false;              ///< 是否为矩阵
    bool is_complex = false;             ///< 是否为复数标量
    bool is_string = false;              ///< 是否为字符串
    mutable bool has_symbolic_text = false;      ///< 是否有符号表达式文本（mutable 支持延迟计算）
    bool has_precise_decimal_text = false; ///< 是否有精确小数文本
    bool exact = false;                  ///< 是否在精确模式下创建

    // 延迟符号计算支持
    mutable bool symbolic_computed = false;  ///< 符号文本是否已计算
    mutable std::string source_expression;   ///< 源表达式（用于延迟符号计算）

    Rational rational;                   ///< 有理数值
    double decimal = 0.0;                ///< 浮点数值
    matrix::ComplexNumber complex;       ///< 复数值
    std::string string_value;            ///< 字符串值
    mutable std::string symbolic_text;   ///< 符号表达式文本（mutable 支持延迟计算）
    std::string precise_decimal_text;    ///< 精确小数文本
    matrix::Matrix matrix;               ///< 矩阵值

    /**
     * @brief 获取符号文本（延迟计算）
     * @param need_symbolic 是否需要符号结果
     * @return 符号文本，如果不需要或无法计算则返回空
     */
    const std::string& get_symbolic_text(bool need_symbolic) const {
        if (!need_symbolic || is_matrix || is_complex || is_string) {
            static const std::string empty;
            return empty;
        }
        if (!symbolic_computed && !source_expression.empty()) {
            // 标记为已计算（实际计算在 format 时进行）
            symbolic_computed = true;
            if (symbolic_text.empty()) {
                symbolic_text = source_expression;
            }
            has_symbolic_text = true;
        }
        return symbolic_text;
    }

    /**
     * @brief 设置源表达式（用于延迟符号计算）
     * @param expr 源表达式
     */
    void set_source_expression(const std::string& expr) {
        source_expression = expr;
        symbolic_computed = false;
    }
};

/**
 * @brief 获取存储值的精确小数文本表示
 * @param value 存储值
 * @return 精确小数文本
 */
std::string stored_value_precise_decimal_text(const StoredValue& value);

#endif // TYPES_STORED_VALUE_H
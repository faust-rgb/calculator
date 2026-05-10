// ============================================================================
// 精确小数解析器接口
// ============================================================================

#ifndef PRECISE_PRECISE_PARSER_H
#define PRECISE_PRECISE_PARSER_H

#include "precise_decimal.h"
#include "types/stored_value.h"

#include <map>
#include <string>

/**
 * @brief 解析精确小数表达式
 * @param expression 表达式字符串
 * @param variables 变量表指针
 * @return 解析结果
 */
PreciseDecimal parse_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables);

#endif // PRECISE_PRECISE_PARSER_H

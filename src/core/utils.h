#ifndef CALCULATOR_UTILS_H
#define CALCULATOR_UTILS_H

#include <string>

namespace utils {

/**
 * @brief 去除字符串首尾空白
 * @param text 输入字符串
 * @return 处理后的字符串
 */
std::string trim_copy(const std::string& text);

/**
 * @brief 检查字符串是否为合法的标识符
 * @param name 字符串
 * @return 是否合法
 */
bool is_valid_identifier(const std::string& name);

} // namespace utils

#endif // CALCULATOR_UTILS_H

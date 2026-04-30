#ifndef CALCULATOR_RESIDUE_H
#define CALCULATOR_RESIDUE_H

#include <string>
#include <vector>
#include "calculator.h"
#include "calculator_internal_types.h"

namespace dsp_ops {

/**
 * @brief 处理 residue 命令，计算有理分式在特定点的留数
 * @param calculator 计算器实例
 * @param impl 计算器实现指针
 * @param arguments 命令参数列表
 * @return 留数的字符串表示（复数向量）
 */
std::string handle_residue_command(Calculator* calculator,
                                   Calculator::Impl* impl,
                                   const std::vector<std::string>& arguments);

} // namespace dsp_ops

#endif

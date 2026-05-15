#ifndef CALCULATOR_RESIDUE_H
#define CALCULATOR_RESIDUE_H

#include <string>
#include <vector>
#include "module/calculator_module.h"

class ServiceLocator;

namespace dsp_ops {

/**
 * @brief 处理 residue 命令，计算有理分式在特定点的留数
 * @param command 命令名称
 * @param arguments 已解析的参数列表
 * @param locator 服务定位器
 * @return 留数的字符串表示（复数向量）
 */
std::string handle_residue_command(const std::string& command,
                                   const std::vector<std::string>& arguments,
                                   ServiceLocator& locator);

} // namespace dsp_ops

#endif

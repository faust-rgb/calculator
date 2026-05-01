#ifndef CALCULATOR_RESIDUE_H
#define CALCULATOR_RESIDUE_H

#include <string>
#include <vector>
#include "../core/calculator_module.h"

namespace dsp_ops {

/**
 * @brief 处理 residue 命令，计算有理分式在特定点的留数
 * @param command 命令名称
 * @param inside 括号内的参数字符串
 * @param svc 核心服务接口
 * @return 留数的字符串表示（复数向量）
 */
std::string handle_residue_command(const std::string& command,
                                   const std::string& inside,
                                   const CoreServices& svc);

} // namespace dsp_ops

#endif

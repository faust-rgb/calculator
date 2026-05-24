/**
 * @file io_module.h
 * @brief 输入输出模块头文件
 *
 * 本文件定义了 IoModule 类，提供文件读写功能的模块接口。
 * 支持文件的打开、关闭、读取、写入等基本操作，以及 CSV 和 JSON 格式的数据读写。
 *
 * @author Calculator Team
 * @date 2024
 */

#ifndef IO_MODULE_H
#define IO_MODULE_H

#include "module/calculator_module.h"
#include "types/stored_value.h"

#include <fstream>
#include <map>
#include <memory>

class ServiceLocator;

/**
 * @class IoModule
 * @brief 提供文件读写功能的模块
 *
 * 该模块继承自 CalculatorModule，为计算器提供完整的文件输入输出能力。
 * 所有打开的文件通过文件描述符（fd）进行管理，支持多种文件打开模式。
 *
 * 支持的函数：
 * - open(path, mode) - 打开文件，返回文件描述符
 * - close(fd) - 关闭文件
 * - read(fd) - 读取整个文件内容
 * - write(fd, text) - 写入文本到文件
 * - read_lines(fd) - 读取所有行到列表
 * - readline(fd) - 读取单行
 * - seek(fd, pos) - 定位文件指针
 * - tell(fd) - 获取当前位置
 * - exists(path) - 检查文件是否存在
 * - delete(path) - 删除文件
 * - read_csv(path) - 从 CSV 文件读取矩阵
 * - write_csv(path, matrix) - 将矩阵写入 CSV 文件
 * - read_json(path) - 从 JSON 文件读取数据
 * - write_json(path, data) - 将数据写入 JSON 文件
 */
class IoModule : public CalculatorModule {
public:
    /**
     * @brief 获取模块名称
     * @return 返回模块名称 "IO"
     */
    std::string name() const override { return "IO"; }

    /**
     * @brief 获取模块支持的所有命令列表
     * @return 返回命令名称字符串向量
     */
    std::vector<std::string> get_commands() const override {
        return {"open", "close", "read", "write", "read_lines", "readline",
                "seek", "tell", "exists", "delete", "read_csv", "write_csv",
                "read_json", "write_json"};
    }

    /**
     * @brief 获取模块提供的原生函数映射
     * @return 返回函数名到函数实现的映射表
     */
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_native_functions() const override;

    /**
     * @brief 执行命令行参数形式的命令
     * @param command 命令名称
     * @param args 参数列表
     * @param locator 服务定位器
     * @return 返回命令执行结果的字符串表示
     */
    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             ServiceLocator& locator) override;

    /**
     * @brief 获取指定主题的帮助信息
     * @param topic 帮助主题名称
     * @return 返回该主题的帮助文本
     */
    std::string get_help_snippet(const std::string& topic) const override;

private:
    /// 文件描述符到文件流的映射表，存储所有打开的文件
    mutable std::map<int, std::shared_ptr<std::fstream>> files_;
    /// 下一个可用的文件描述符，从1开始递增
    mutable int next_fd_ = 1;
};

#endif // IO_MODULE_H

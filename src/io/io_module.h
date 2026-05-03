#ifndef IO_MODULE_H
#define IO_MODULE_H

#include "module/calculator_module.h"
#include "types/stored_value.h"

#include <fstream>
#include <map>
#include <memory>

/**
 * @class IoModule
 * @brief 提供文件读写功能的模块
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
    std::string name() const override { return "IO"; }

    std::vector<std::string> get_commands() const override {
        return {"open", "close", "read", "write", "read_lines", "readline",
                "seek", "tell", "exists", "delete", "read_csv", "write_csv",
                "read_json", "write_json"};
    }

    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_native_functions() const override;

    std::string execute_args(const std::string& command,
                             const std::vector<std::string>& args,
                             const CoreServices& services) override;

    std::string get_help_snippet(const std::string& topic) const override;

private:
    mutable std::map<int, std::shared_ptr<std::fstream>> files_;
    mutable int next_fd_ = 1;
};

#endif // IO_MODULE_H

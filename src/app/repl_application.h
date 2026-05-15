// ============================================================================
// repl_application.h - 交互式命令行应用程序
// ============================================================================

#ifndef APP_REPL_APPLICATION_H
#define APP_REPL_APPLICATION_H

#include "core/api/calculator.h"
#include <string>
#include <vector>

/**
 * @class ReplApplication
 * @brief 负责处理交互式命令行界面的逻辑
 */
class ReplApplication {
public:
    explicit ReplApplication(Calculator& calculator);
    
    /**
     * @brief 运行 REPL 循环
     * @param plain_mode 是否为简洁模式（禁用提示符和彩色输出）
     * @return 退出码
     */
    int run(bool plain_mode);

private:
    std::string read_line(const std::string& prompt, bool plain_mode);
    std::string execute_line(const std::string& line, bool* should_exit);
    
    Calculator& calculator_;
    std::vector<std::string> history_;
    bool exact_mode_ = false;
};

#endif // APP_REPL_APPLICATION_H

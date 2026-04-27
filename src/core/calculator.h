#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <memory>
#include <string>
#include <vector>

/**
 * @class Calculator
 * @brief 高级科学计算器的核心类，提供表达式求值、变量管理、脚本执行等功能
 *
 * 采用 Pimpl (Pointer to Implementation) 设计模式隐藏内部实现细节，
 * 减少编译依赖并允许在不修改头文件的情况下更改内部数据结构。
 *
 * 主要功能：
 * - 数学表达式求值（支持标准数学函数、变量、常量）
 * - 变量存储与管理（支持标量和矩阵）
 * - 脚本执行（支持自定义函数、控制流）
 * - 进制转换（二进制、八进制、十进制、十六进制）
 * - 整数因式分解
 * - 状态保存与加载
 */
class Calculator {
public:
    /** @brief 前向声明的实现类，用于 Pimpl 模式 */
    struct Impl;

    /** @brief 默认构造函数，初始化计算器状态 */
    Calculator();

    /** @brief 析构函数，释放实现类资源 */
    ~Calculator();

    /**
     * @brief 计算数学表达式的数值结果
     * @param expression 数学表达式字符串，如 "2 + 3 * sin(pi/4)"
     * @return 表达式的数值结果
     * @throw std::runtime_error 表达式语法错误或计算错误时抛出
     */
    double evaluate(const std::string& expression);

    /**
     * @brief 计算表达式的原始数值结果，不做显示层归整
     * @param expression 数学表达式字符串
     * @return 原始浮点结果
     */
    double evaluate_raw(const std::string& expression);

    /**
     * @brief 计算表达式并以字符串形式返回，支持精确模式
     * @param expression 数学表达式字符串
     * @param exact_mode 是否使用精确模式（尝试保留分数形式）
     * @return 格式化后的结果字符串
     */
    std::string evaluate_for_display(const std::string& expression, bool exact_mode);

    /**
     * @brief 处理单行输入（命令或表达式）
     * @param expression 输入字符串
     * @param exact_mode 是否使用精确模式
     * @return 处理结果字符串
     *
     * 支持命令：help, vars, clear, save, load, exit 等
     */
    std::string process_line(const std::string& expression, bool exact_mode);

    /**
     * @brief 执行脚本代码
     * @param source 脚本源代码
     * @param exact_mode 是否使用精确模式
     * @return 脚本执行结果或最后表达式的值
     */
    std::string execute_script(const std::string& source, bool exact_mode);

    /**
     * @brief 获取帮助文本
     * @return 包含所有可用命令和函数的说明文本
     */
    std::string help_text() const;

    /**
     * @brief 获取特定主题的帮助
     * @param topic 帮助主题，如 "functions", "variables", "matrices"
     * @return 该主题的详细说明
     */
    std::string help_topic(const std::string& topic) const;

    /**
     * @brief 列出所有已定义的变量
     * @return 变量列表字符串
     */
    std::string list_variables() const;

    /**
     * @brief 清除指定变量
     * @param name 变量名
     * @return 操作结果消息
     */
    std::string clear_variable(const std::string& name);

    /**
     * @brief 清除所有变量
     * @return 操作结果消息
     */
    std::string clear_all_variables();

    /**
     * @brief 对整数表达式进行因式分解
     * @param expression 整数表达式，如 "120"
     * @return 因式分解结果，如 "2^3 * 3 * 5"
     */
    std::string factor_expression(const std::string& expression) const;

    /**
     * @brief 执行进制转换
     * @param expression 转换表达式，如 "0b1010 to hex"
     * @return 转换结果
     */
    std::string base_conversion_expression(const std::string& expression) const;

    /**
     * @brief 尝试处理函数分析相关命令
     * @param expression 输入表达式
     * @param output 输出结果存储位置
     * @return 是否成功识别并处理为函数命令
     *
     * 支持的命令：limit, derivative, integral, extrema 等
     */
    bool try_process_function_command(const std::string& expression, std::string* output);

    /**
     * @brief 保存计算器状态到文件
     * @param path 文件路径
     * @return 操作结果消息
     */
    std::string save_state(const std::string& path) const;

    /**
     * @brief 从文件加载计算器状态
     * @param path 文件路径
     * @return 操作结果消息
     */
    std::string load_state(const std::string& path);

    /**
     * @brief 设置十六进制输出是否带 0x/0X 前缀
     * @param enabled 是否启用前缀
     * @return 操作结果消息
     */
    std::string set_hex_prefix_mode(bool enabled);

    /**
     * @brief 查询十六进制输出前缀模式
     * @return true 如果启用了前缀
     */
    bool hex_prefix_mode() const;

    /**
     * @brief 设置十六进制字母是否使用大写
     * @param enabled true 表示大写，false 表示小写
     * @return 操作结果消息
     */
    std::string set_hex_uppercase_mode(bool enabled);

    /**
     * @brief 查询十六进制大小写模式
     * @return true 如果使用大写字母
     */
    bool hex_uppercase_mode() const;

    /**
     * @brief 设置符号常量模式
     * @param enabled 是否启用符号常量（如 pi, e 以符号形式保留）
     * @return 操作结果消息
     */
    std::string set_symbolic_constants_mode(bool enabled);

    /**
     * @brief 检查符号常量模式是否启用
     * @return true 如果符号常量模式已启用
     */
    bool symbolic_constants_mode() const;

    /**
     * @brief 设置十进制输出显示有效位数
     * @param precision 有效位数，范围 1..17
     * @return 操作结果消息
     */
    std::string set_display_precision(int precision);

    /**
     * @brief 查询当前十进制输出显示有效位数
     * @return 当前有效位数
     */
    int display_precision() const;

    /**
     * @brief 列出当前会话可补全的变量名
     * @return 变量名列表
     */
    std::vector<std::string> variable_names() const;

    /**
     * @brief 列出当前会话可补全的自定义函数名
     * @return 函数名列表
     */
    std::vector<std::string> custom_function_names() const;

private:
    /** @brief Pimpl 模式的实现指针 */
    std::unique_ptr<Impl> impl_;

    /**
     * @brief 规范化计算结果（处理 -0 等边界情况）
     * @param value 原始计算结果
     * @return 规范化后的结果
     */
    static double normalize_result(double value);
};

#endif

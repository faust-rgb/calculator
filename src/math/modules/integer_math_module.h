/**
 * @file integer_math_module.h
 * @brief 整数数学模块头文件
 *
 * 定义 IntegerMathModule 类，提供整数数学、数论和进制转换功能。
 */

#ifndef INTEGER_MATH_MODULE_H
#define INTEGER_MATH_MODULE_H

#include "module/calculator_module.h"

class ServiceLocator;

/**
 * @class IntegerMathModule
 * @brief 提供整数数学、数论和进制转换功能的模块
 *
 * 该模块实现了以下功能：
 * - 因式分解（factor）
 * - 进制转换（bin、oct、hex、base）
 * - 位运算（and、or、xor、not、shl、shr、rol、ror）
 * - 位统计（popcount、bitlen、ctz、clz、parity、reverse_bits）
 * - 数论函数（gcd、lcm、mod、factorial、nCr、nPr、fib）
 * - 素数操作（is_prime、next_prime、prev_prime）
 * - 数论函数（euler_phi、mobius、prime_pi、egcd）
 */
class IntegerMathModule : public CalculatorModule {
public:
    /// 获取模块名称
    std::string name() const override { return "IntegerMath"; }

    /// 获取模块提供的命令列表
    std::vector<std::string> get_commands() const override {
        return {"factor", "bin", "oct", "hex", "base"};
    }

    /**
     * @brief 执行命令式操作
     * @param command 命令名称
     * @param args 参数列表
     * @param locator 服务定位器
     * @return 执行结果的字符串表示
     */
    std::string execute_args_view(std::string_view command,
                                  const std::vector<std::string_view>& args,
                                  ServiceLocator& locator) override;

    /** @brief 返回支持的帮助主题 */
    std::vector<std::string> get_help_topics() const override;

    /**
     * @brief 获取模块提供的标量函数映射
     * @return 函数名称到函数实现的映射
     */
    std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>> get_scalar_functions() const override;

    /**
     * @brief 获取模块提供的所有函数名称列表
     * @return 函数名称列表
     */
    std::vector<std::string> get_functions() const override;

    /**
     * @brief 获取帮助信息片段
     * @param topic 主题名称
     * @return 对应主题的帮助文本
     */
    std::string get_help_snippet(const std::string& topic) const override;
};

#endif

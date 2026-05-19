/**
 * @file statistics_module.h
 * @brief 统计模块头文件
 *
 * 本文件定义了 StatisticsModule 类，作为计算器统计功能的模块化接口。
 * 该模块提供了统计摘要、描述性统计量计算等功能，可被计算器主程序动态加载。
 */

#ifndef STATISTICS_MODULE_H
#define STATISTICS_MODULE_H

#include "module/calculator_module.h"
#include "parser/grammars/unified_expression_parser.h"
#include "app/scalar_type.h"
#include "statistics/calculator_statistics.h"
#include "statistics/statistics.h"
#include "statistics/probability.h"
#include "math/mymath.h"
#include <numeric>
#include <sstream>
#include <span>

class ServiceLocator;

/**
 * @brief 统计模块类
 *
 * 继承自 CalculatorModule 基类，提供统计相关的命令和函数。
 * 支持的命令包括：
 * - describe / stat_summary: 生成完整的统计摘要
 *
 * 支持的原生函数包括：
 * - mean / avg: 平均值
 * - median: 中位数
 * - std: 标准差
 * - var: 方差
 * - t_test2: 双样本 T 检验
 * - chi2_test: 卡方检验
 */
class StatisticsModule : public CalculatorModule {
public:
    using Scalar = mymath::Scalar;
    /**
     * @brief 获取模块元数据
     * @return 模块元数据
     */
    ModuleMetadata get_metadata() const override {
        return ModuleMetadata(
            "Statistics",
            "1.0.0",
            "Statistical analysis module",
            "Calculator Team",
            {}  // 无依赖
        );
    }

    /**
     * @brief 执行统计命令
     *
     * 处理统计相关的命令请求，如统计摘要生成等。
     */
    std::string execute_command(const CommandASTNode& node, ServiceLocator& locator) override;

    /**
     * @brief 获取模块提供的原生函数映射
     */
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> get_native_functions() const override;

    /**
     * @brief 获取帮助信息片段
     *
     * 根据主题返回相关的帮助文本。
     *
     * @param topic 帮助主题
     * @return 帮助文本
     */
    std::string get_help_snippet(const std::string& topic) const override {
        if (topic == "analysis") {
            return "Statistics:\n"
                   "  describe(data)      Full statistical summary\n"
                   "  stat_summary(data)  Full statistical summary";
        }
        return "";
    }

};

#endif
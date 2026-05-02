// ============================================================================
// 高精度计算模块
// ============================================================================

#ifndef PRECISE_PRECISE_MODULE_H
#define PRECISE_PRECISE_MODULE_H

#include "module/calculator_module.h"
#include "types/stored_value.h"

#include <map>
#include <string>

/**
 * @class PreciseModule
 * @brief 高精度计算模块
 *
 * 提供精确小数计算功能，避免浮点误差。
 * 当表达式满足特定条件时自动尝试使用精确模式计算。
 */
class PreciseModule : public CalculatorModule {
public:
    std::string name() const override { return "PreciseDecimal"; }

    bool wants_implicit_evaluation() const override { return true; }

    std::string get_implicit_trigger_chars() const override;

    /**
     * @brief 尝试使用精确模式隐式求值
     *
     * @param expression 表达式
     * @param out 输出值
     * @param variables 变量表
     * @return true 如果成功使用精确模式计算
     */
    bool try_evaluate_implicit(const std::string& expression,
                              StoredValue* out,
                              const std::map<std::string, StoredValue>& variables) const override;

private:
    /**
     * @brief 判断是否应该尝试精确小数计算
     */
    bool should_try_precise_decimal_expression(
        const std::string& expression,
        const std::map<std::string, StoredValue>& variables) const;
};

#endif // PRECISE_PRECISE_MODULE_H

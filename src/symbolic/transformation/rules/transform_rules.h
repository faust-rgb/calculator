// ============================================================================
// 变换规则框架
// ============================================================================
//
// 本文件实现可扩展的变换规则框架，支持：
// - Laplace 变换
// - Fourier 变换
// - Z 变换
// - 用户自定义变换规则
//
// 规则通过规则表注册，支持优先级排序和模式匹配。

#ifndef TRANSFORM_RULES_H
#define TRANSFORM_RULES_H

#include "symbolic/core/symbolic_expression.h"

#include <string>
#include <vector>
#include <map>
#include <functional>
#include <memory>

namespace transform_rules {

/**
 * @struct TransformRule
 * @brief 变换规则结构
 *
 * 表示一条变换规则，包含模式匹配器和变换函数。
 */
struct TransformRule {
    std::string name;                    ///< 规则名称
    std::string transform_type;          ///< 变换类型: "laplace", "fourier", "z"

    /**
     * @brief 模式匹配函数
     * @param expr 输入表达式
     * @param input_var 输入变量 (如 t)
     * @return 是否匹配此规则
     */
    std::function<bool(const SymbolicExpression& expr,
                       const std::string& input_var)> matcher;

    /**
     * @brief 变换函数
     * @param expr 输入表达式
     * @param input_var 输入变量
     * @param output_var 输出变量 (如 s)
     * @return 变换结果
     */
    std::function<SymbolicExpression(const SymbolicExpression& expr,
                                     const std::string& input_var,
                                     const std::string& output_var)> transformer;

    int priority = 0;                    ///< 优先级 (高优先级先尝试)
    std::string description;             ///< 规则描述
    std::string formula;                 ///< 数学公式 (文档用)
};

/**
 * @class TransformRuleRegistry
 * @brief 变换规则注册表
 *
 * 单例类，管理所有变换规则。支持：
 * - 注册内置规则
 * - 添加用户自定义规则
 * - 按优先级查询规则
 */
class TransformRuleRegistry {
public:
    static TransformRuleRegistry& instance();

    /**
     * @brief 注册变换规则
     * @param rule 规则对象
     */
    void register_rule(const TransformRule& rule);

    /**
     * @brief 获取指定类型的所有规则 (按优先级排序)
     * @param transform_type 变换类型
     * @return 规则列表
     */
    std::vector<TransformRule> get_rules(const std::string& transform_type) const;

    /**
     * @brief 添加用户自定义规则
     * @param rule 规则对象
     */
    void add_user_rule(const TransformRule& rule);

    /**
     * @brief 清除所有用户规则
     */
    void clear_user_rules();

    /**
     * @brief 初始化内置规则
     */
    void initialize_builtin_rules();

private:
    TransformRuleRegistry();
    ~TransformRuleRegistry() = default;

    std::map<std::string, std::vector<TransformRule>> rules_;
    std::vector<TransformRule> user_rules_;
    bool initialized_ = false;
};

// ============================================================================
// 内置规则辅助函数
// ============================================================================

/**
 * @brief 检查表达式是否为常数
 */
bool is_constant_expression(const SymbolicExpression& expr, const std::string& var);

/**
 * @brief 检查表达式是否为变量本身
 */
bool is_variable_expression(const SymbolicExpression& expr, const std::string& var);

/**
 * @brief 检查表达式是否为变量的幂次
 */
bool is_power_of_variable(const SymbolicExpression& expr, const std::string& var, int* power);

/**
 * @brief 检查表达式是否为指数函数 exp(a*x + b)
 */
bool is_exponential_form(const SymbolicExpression& expr, const std::string& var,
                         SymbolicExpression* coefficient, SymbolicExpression* constant);

/**
 * @brief 检查表达式是否为三角函数 sin(a*x + b) 或 cos(a*x + b)
 */
bool is_trig_form(const SymbolicExpression& expr, const std::string& func_name,
                  const std::string& var, SymbolicExpression* coefficient,
                  SymbolicExpression* constant);

} // namespace transform_rules

#endif // TRANSFORM_RULES_H

// ============================================================================
// service_access.h - 模块服务访问助手
// ============================================================================
//
// 提供统一的服务访问模式，简化模块对核心服务的使用。
// 包括：
// - ServiceAccess 类：封装常用服务访问
// - 快捷方法：evaluate, parse_decimal 等
//
// 使用示例：
//   ServiceAccess access(locator);
//   Scalar value = access.evaluate_decimal("1+2");
//   Matrix m = access.parse_matrix("[1,2;3,4]");
// ============================================================================

#ifndef MODULE_SERVICE_ACCESS_H
#define MODULE_SERVICE_ACCESS_H

#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "core/services/service_interfaces.h"
#include "core/types/module_types.h"
#include "app/scalar_type.h"

#include <string>
#include <functional>

namespace module {

/**
 * @class ServiceAccess
 * @brief 模块服务访问助手，提供统一的服务访问模式
 *
 * 封装 ServiceLocator 的常用解析操作，简化模块代码。
 */
class ServiceAccess {
public:
    /**
     * @brief 构造服务访问助手
     * @param locator 服务定位器
     */
    explicit ServiceAccess(ServiceLocator& locator) : locator_(locator) {}

    // ========================================================================
    // 主要服务访问
    // ========================================================================

    /**
     * @brief 获取求值引擎
     * @return 求值引擎引用
     */
    IEvaluationEngine& engine() {
        return *locator_.resolve<IEvaluationEngine>();
    }

    /**
     * @brief 获取执行上下文
     * @return 执行上下文引用
     */
    IExecutionContext& context() {
        return *locator_.resolve<IExecutionContext>();
    }

    /**
     * @brief 获取核心服务接口
     * @return 核心服务引用
     */
    const CoreServices& services() {
        return context().services();
    }

    // ========================================================================
    // 变量和函数管理
    // ========================================================================

    /**
     * @brief 获取变量管理器
     * @return 变量管理器引用
     */
    IVariableManager& variables() {
        return *locator_.resolve<IVariableManager>();
    }

    /**
     * @brief 获取函数管理器
     * @return 函数管理器引用
     */
    IFunctionManager& functions() {
        return *locator_.resolve<IFunctionManager>();
    }

    /**
     * @brief 获取配置管理器
     * @return 配置管理器引用
     */
    IConfigManager& config() {
        return *locator_.resolve<IConfigManager>();
    }

    // ========================================================================
    // 便捷求值方法
    // ========================================================================

    /**
     * @brief 解析小数表达式
     * @param expr 表达式字符串
     * @return 解析后的数值
     */
    Scalar parse_decimal(const std::string& expr) {
        return engine().parse_decimal(expr);
    }

    /**
     * @brief 求值表达式
     * @param expr 表达式字符串
     * @return 求值结果
     */
    Scalar evaluate(const std::string& expr) {
        return engine().evaluate(expr);
    }

    /**
     * @brief 求值表达式为存储值
     * @param expr 表达式字符串
     * @param exact_mode 是否精确模式
     * @return 存储值结果
     */
    StoredValue evaluate_value(const std::string& expr, bool exact_mode = false) {
        return engine().evaluate_expression_value(expr, exact_mode);
    }

    /**
     * @brief 检查参数是否为矩阵
     * @param arg 参数字符串
     * @return true 如果是矩阵格式
     */
    bool is_matrix_argument(const std::string& arg) {
        return engine().is_matrix_argument(arg);
    }

    /**
     * @brief 解析矩阵参数
     * @param arg 参数字符串
     * @param context 上下文名称（用于错误信息）
     * @return 解析后的矩阵
     */
    matrix::Matrix parse_matrix(const std::string& arg, const std::string& context = "") {
        return engine().parse_matrix_argument(arg, context);
    }

    /**
     * @brief 规范化结果值
     * @param value 输入值
     * @return 规范化后的值
     */
    Scalar normalize_result(Scalar value) {
        return engine().normalize_result(value);
    }

private:
    ServiceLocator& locator_;
};

} // namespace module

#endif // MODULE_SERVICE_ACCESS_H
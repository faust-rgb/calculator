// ============================================================================
// Calculator 实现类
// ============================================================================
//
// Calculator 类的内部实现，使用 Pimpl 模式隐藏细节。
// 存储计算器的所有状态：变量、函数、模块、显示选项等。
//
// 注意：此文件仅供 Calculator 内部使用，不应被外部模块直接引用。
// ============================================================================

#ifndef CORE_CALCULATOR_IMPL_H
#define CORE_CALCULATOR_IMPL_H

#include "core/api/calculator.h"
#include "core/environment/scope.h"
#include "core/execution_context.h"
#include "types/scalar_type.h"
#include "math/runtime/precision/default_precision.h"
#include "execution/registry/command_registry.h"
#include "execution/functions/user_function.h"

#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// 前向声明
class CalculatorModule;
struct CoreServices;

// ============================================================================
// 显示精度常量（使用 app 命名空间中的统一定义）
// ============================================================================

namespace core {

using Scalar = mymath::Scalar;

/** @brief 判断数值是否为零的显示阈值 */
inline Scalar kDisplayZeroEpsScalar() { return Scalar(app::display_zero_threshold()); }

/** @brief 判断数值是否为整数的显示阈值 */
inline Scalar kDisplayIntegerEpsScalar() { return Scalar(app::display_integer_threshold()); }

} // namespace core

/** @brief 判断数值是否为零的显示阈值 */
inline Scalar kDisplayZeroEps() { return Scalar(app::display_zero_threshold()); }

/** @brief 判断数值是否为整数的显示阈值 */
inline Scalar kDisplayIntegerEps() { return Scalar(app::display_integer_threshold()); }

// 使用 app 命名空间中的统一定义
using app::kDefaultDisplayPrecision;
using app::kMinDisplayPrecision;
using app::kMaxDisplayPrecision;

// ============================================================================
// Calculator 实现类
// ============================================================================

#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "execution/engine/script_context.h"

// ============================================================================
// Calculator 实现类
// ============================================================================

/**
 * @struct Calculator::Impl
 * @brief Calculator 的内部实现
 *
 * 存储计算器的所有状态：
 * - 变量表（全局和局部作用域）
 * - 函数表（简单函数和脚本函数）
 * - 模块注册
 * - 显示选项
 */
struct Calculator::Impl : public IExecutionContext {
    Calculator* parent = nullptr;
    ServiceLocator locator;
    core::ExecutionContext execution_ctx;

    // 快捷访问服务接口（延迟初始化）
    std::shared_ptr<IVariableManager> variables_ptr;
    std::shared_ptr<IFunctionManager> functions_ptr;
    std::shared_ptr<IConfigManager> config_ptr;
    std::shared_ptr<ICommandRegistry> commands_ptr;
    std::shared_ptr<IModuleManager> modules;
    std::shared_ptr<IStatePersistence> persistence;

    // IExecutionContext implementation
    core::ExecutionContext& core_context() override { return execution_ctx; }
    const core::ExecutionContext& core_context() const override { return execution_ctx; }

    IVariableManager& variables() override { return *variables_ptr; }
    const IVariableManager& variables() const override { return *variables_ptr; }
    
    IFunctionManager& functions() override { return *functions_ptr; }
    const IFunctionManager& functions() const override { return *functions_ptr; }
    
    IConfigManager& config() override { return *config_ptr; }
    const IConfigManager& config() const override { return *config_ptr; }
    
    ICommandRegistry& commands() override { return *commands_ptr; }
    const ICommandRegistry& commands() const override { return *commands_ptr; }
    
    const CoreServices& services() const override;
    StoredValue evaluate(const std::string& expression, bool exact_mode) override;
    bool try_evaluate_implicit(const std::string& expression, 
                               StoredValue* value, 
                               const std::map<std::string, StoredValue>& vars) override;
    std::string expand_inline(const std::string& expression) override;
    std::string execute_script_file(const std::string& path, 
                                   bool exact_mode, 
                                   bool create_scope) override;
    bool try_process_function_command(const std::string& expression, 
                                     std::string* output, 
                                     bool exact_mode = false) override;

    std::vector<std::string> module_commands;
    std::vector<std::string> module_functions;

    std::map<std::string, std::vector<std::shared_ptr<CalculatorModule>>> help_topic_to_modules;
    std::vector<std::shared_ptr<CalculatorModule>> implicit_evaluation_modules;

    // 优化的隐式求值：触发字符到模块的映射
    std::array<std::vector<std::shared_ptr<CalculatorModule>>, 256> trigger_char_to_modules;

    // 核心服务缓存
    std::unique_ptr<CoreServices> core_services;

    // 脚本执行状态
    ScriptContext script_context;
};

// ============================================================================
// 辅助函数声明
// ============================================================================

// 显示精度
void set_process_display_precision(int precision);
int process_display_precision();

// 字符串处理
bool is_valid_variable_name(std::string_view name);
bool is_identifier_text(std::string_view text);
bool is_string_literal(std::string_view text);
std::string parse_string_literal_value(std::string_view text);
bool is_reserved_user_function_name(IExecutionContext* ctx, std::string_view name);

// 状态持久化
std::string encode_state_field(const std::string& text);
std::string decode_state_field(const std::string& text);

// 值格式化
void apply_calculator_display_precision(const Calculator::Impl* impl);

// 线性方程组
std::vector<Scalar> solve_dense_linear_system(
    std::vector<std::vector<Scalar>> matrix,
    std::vector<Scalar> rhs,
    const std::string& context);

// 模块注册
void register_standard_modules(Calculator* calculator);

#endif // CORE_CALCULATOR_IMPL_H

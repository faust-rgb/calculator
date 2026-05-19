#include "execution/engine/inline_expander.h"
// ============================================================================
// Calculator 核心实现
// ============================================================================
//
// 本文件合并了以下原本分离的实现：
// - calculator_lifecycle.cpp: 生命周期、设置、模块注册
// - calculator_commands.cpp: 命令处理与分发
// - calculator_basic_commands.cpp: 基础求值与显示
// - calculator_state_persistence.cpp: 状态保存与加载
// ============================================================================

#include "core/api/calculator_internal_types.h"
#include "app/scalar_type.h"
#include "parser/grammars/unified_expression_parser.h"
#include "parser/grammars/command_parser.h"
#include "analysis/calculus/function_analysis.h"
#include "matrix/matrix.h"
#include "math/mymath.h"
#include "symbolic/core/symbolic_expression.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "parser/grammars/command_parser.h"
#include "execution/engine/script_runtime.h"
#include "parser/grammars/script_parser.h"
#include "module/module_registration.h"
#include "module/calculator_module.h"
#include "math/helpers/integer_helpers.h"
#include "plot/calculator_plot.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "core/services/core_managers.h"
#include "execution/registry/command_registry.h"
#include "core/services/evaluation_engine_impl.h"
#include "core/services/state_persistence.h"

// ============================================================================
// 辅助函数
// ============================================================================

namespace {

// 检查表达式是否包含触发字符
bool has_trigger_char(std::string_view expression, const std::array<bool, 256>& table) {
    for (char c : expression) {
        if (table[static_cast<unsigned char>(c)]) {
            return true;
        }
    }
    return false;
}

} // namespace

void apply_calculator_display_precision(const Calculator::Impl* impl) {
    const int precision = (impl == nullptr || !impl->config_ptr) ? kDefaultDisplayPrecision : impl->config_ptr->get_display_precision();
    set_process_display_precision(precision);
    matrix::set_display_precision(precision);
    SymbolicExpression::set_display_precision(precision);
}

void broadcast_settings(Calculator* calculator, Calculator::Impl* impl) {
    (void)calculator;
    if (!impl->config_ptr || !impl->modules) return;

    CalculatorSettings settings;
    settings.display_precision = impl->config_ptr->get_display_precision();
    settings.exact_mode = impl->config_ptr->is_exact_mode();
    settings.symbolic_constants_mode = impl->config_ptr->is_symbolic_constants_mode();
    settings.hex_prefix_mode = impl->config_ptr->is_hex_prefix_mode();
    settings.hex_uppercase_mode = impl->config_ptr->is_hex_uppercase_mode();

    for (auto& module : impl->modules->get_all_modules()) {
        module->on_settings_changed(settings);
    }
}

// ============================================================================
// 生命周期与模块注册
// ============================================================================

Calculator::Calculator() : impl_(new Impl()) {
    impl_->parent = this;
    // 1. 注册核心管理服务到 Locator
    impl_->variables_ptr = std::make_shared<VariableManager>();
    impl_->functions_ptr = std::make_shared<FunctionManager>();
    impl_->config_ptr = std::make_shared<ConfigManager>();
    impl_->commands_ptr = std::make_shared<CommandRegistry>();
    impl_->modules = std::make_shared<ModuleManager>(impl_->locator);

    // 创建状态持久化服务
    auto persistence = std::make_shared<StatePersistenceService>(impl_->variables_ptr, impl_->functions_ptr);
    persistence->set_script_executor([this](const std::string& source) {
        execute_script(source, false);
    });
    impl_->persistence = persistence;

    impl_->locator.register_service<IVariableManager>(impl_->variables_ptr);
    impl_->locator.register_service<IFunctionManager>(impl_->functions_ptr);
    impl_->locator.register_service<IConfigManager>(impl_->config_ptr);
    impl_->locator.register_service<ICommandRegistry>(impl_->commands_ptr);
    impl_->locator.register_service<IModuleManager>(impl_->modules);
    impl_->locator.register_service<IStatePersistence>(impl_->persistence);

    // 注册执行上下文
    auto ctx_ptr = std::shared_ptr<IExecutionContext>(impl_.get(), [](IExecutionContext*){});
    impl_->locator.register_service<IExecutionContext>(ctx_ptr);


    // 2. 注册逻辑引擎服务
    auto engine = std::make_shared<EvaluationEngineImpl>(this, impl_.get());
    impl_->locator.register_service<IEvaluationEngine>(engine);

    apply_calculator_display_precision(impl_.get());

    // 3. 注册所有标准模块（按依赖拓扑顺序）
    register_standard_modules(this);
    broadcast_settings(this, impl_.get());
}

Calculator::~Calculator() = default;

void Calculator::register_module(std::shared_ptr<CalculatorModule> module) {
    if (!module) return;

    // 隐式求值路由优化
    if (module->wants_implicit_evaluation()) {
        impl_->implicit_evaluation_modules.push_back(module);
        const auto* trigger_table = module->get_cached_trigger_table();
        if (trigger_table) {
            for (int c = 0; c < 256; ++c) {
                if ((*trigger_table)[c]) {
                    impl_->trigger_char_to_modules[c].push_back(module);
                }
            }
        }
    }

    // 委托给 ModuleManager 处理核心注册逻辑
    impl_->modules->register_module(module);

    // 记录命令和函数名用于补全
    auto specs = module->get_command_specs();
    for (const auto& spec : specs) {
        impl_->module_commands.push_back(command_key_display(spec.key));
    }
    auto funcs = module->get_functions();
    impl_->module_functions.insert(impl_->module_functions.end(), funcs.begin(), funcs.end());
}

bool is_reserved_user_function_name(IExecutionContext* ctx, std::string_view name) {
    if (utils::is_reserved_function_name(name)) {
        return true;
    }
    if (ctx == nullptr) {
        return false;
    }
    const std::string name_text(name);
    // 检查是否是已注册的命令（使用 CommandRegistry）
    if (ctx->commands().has_command(name_text)) {
        return true;
    }
    return ctx->functions().get_scalar_functions()->count(name_text) > 0 ||
           ctx->functions().get_matrix_functions()->count(name_text) > 0 ||
           ctx->functions().get_value_functions()->count(name_text) > 0;
}


// ============================================================================
// 变量与设置
// ============================================================================

std::string Calculator::clear_variable(const std::string& name) {
    if (!impl_->variables_ptr->has(name)) throw std::runtime_error("unknown variable: " + name);
    impl_->variables_ptr->remove(name);
    return "Cleared variable: " + name;
}

std::string Calculator::clear_all_variables() {
    impl_->variables_ptr->clear_all();
    return "Cleared all variables.";
}

std::string Calculator::set_hex_prefix_mode(bool enabled) {
    impl_->config_ptr->set_hex_prefix_mode(enabled);
    broadcast_settings(this, impl_.get());
    return std::string("Hex prefix mode: ") + (enabled ? "ON" : "OFF");
}

bool Calculator::hex_prefix_mode() const { return impl_->config_ptr->is_hex_prefix_mode(); }

std::string Calculator::set_hex_uppercase_mode(bool enabled) {
    impl_->config_ptr->set_hex_uppercase_mode(enabled);
    broadcast_settings(this, impl_.get());
    return std::string("Hex letter case: ") + (enabled ? "UPPER" : "LOWER");
}

bool Calculator::hex_uppercase_mode() const { return impl_->config_ptr->is_hex_uppercase_mode(); }

std::string Calculator::set_symbolic_constants_mode(bool enabled) {
    impl_->config_ptr->set_symbolic_constants_mode(enabled);
    broadcast_settings(this, impl_.get());
    return std::string("Symbolic constants mode: ") + (enabled ? "ON" : "OFF");
}

bool Calculator::symbolic_constants_mode() const { return impl_->config_ptr->is_symbolic_constants_mode(); }

std::string Calculator::set_display_precision(int precision) {
    if (precision < kMinDisplayPrecision || precision > kMaxDisplayPrecision) {
        throw std::runtime_error("display precision must be in range " +
                                 std::to_string(kMinDisplayPrecision) + ".." +
                                 std::to_string(kMaxDisplayPrecision));
    }
    impl_->config_ptr->set_display_precision(precision);
    apply_calculator_display_precision(impl_.get());
    broadcast_settings(this, impl_.get());
    return "Display precision: " + std::to_string(precision);
}

int Calculator::display_precision() const { return impl_->config_ptr->get_display_precision(); }
std::vector<std::string> Calculator::module_command_names() const { return impl_->module_commands; }
std::vector<std::string> Calculator::module_function_names() const { return impl_->module_functions; }
std::vector<std::string> Calculator::variable_names() const {
    return impl_->variables_ptr->get_all_names();
}
std::vector<std::string> Calculator::custom_function_names() const {
    return impl_->functions_ptr->get_custom_names();
}

Scalar Calculator::normalize_result(Scalar value) {
    using Scalar = mymath::Scalar;
    if (!mymath::isfinite(value)) return value;
    Scalar v(value);
    if constexpr (mymath::is_scalar_float128) {
        if (mymath::abs(v) < Scalar(1e-70L)) return Scalar(0.0L);
    } else if constexpr (!mymath::is_scalar_precise_decimal) {
        if (mymath::abs(v) < kDisplayZeroEps()) return Scalar(0.0);
    }
    // Only round to integer if within long long range
    const Scalar kMaxLL(static_cast<long double>(std::numeric_limits<long long>::max()));
    const Scalar kMinLL(static_cast<long double>(std::numeric_limits<long long>::min()));
    if (mymath::abs(v) > kDisplayIntegerEps() &&
        is_integer_double(value, 1e-9) &&
        v < kMaxLL && v > kMinLL) {
        return (round_to_long_long(value));
    }
    return value;
}

// ============================================================================
// 命令处理与分发
// ============================================================================

bool Calculator::try_process_function_command(const std::string& expression,
                                              std::string* output, bool exact_mode) const {
    auto is_command = [this](std::string_view name) {
        return impl_->commands_ptr->has_command(std::string(name)) ||
               impl_->commands_ptr->has_command(":" + std::string(name));
    };

    CommandASTNode root = parse_command(expression, is_command);
    if (root.kind == CommandKind::kEmpty) return false;

    *output = execute_command_ast(impl_.get(), root, exact_mode);
    return true;
}

bool Calculator::try_evaluate_implicit(std::string_view expression,
                                       StoredValue* output,
                                       const std::map<std::string, StoredValue>& vars) const {
    if (expression.empty()) return false;

    // 优化：使用触发字符到模块的映射，直接找到相关模块
    // 首先收集表达式中出现的所有触发字符对应的模块
    std::set<std::shared_ptr<CalculatorModule>> candidate_modules;

    for (char c : expression) {
        const auto& modules_for_char = impl_->trigger_char_to_modules[static_cast<unsigned char>(c)];
        for (const auto& mod : modules_for_char) {
            candidate_modules.insert(mod);
        }
    }

    // 如果没有找到任何候选模块，回退到遍历所有隐式求值模块
    if (candidate_modules.empty() && !impl_->implicit_evaluation_modules.empty()) {
        for (const auto& module : impl_->implicit_evaluation_modules) {
            if (module->try_evaluate_implicit(std::string(expression), output, vars)) {
                return true;
            }
        }
        return false;
    }

    // 尝试候选模块
    for (const auto& module : candidate_modules) {
        if (module->try_evaluate_implicit(std::string(expression), output, vars)) {
            return true;
        }
    }
    return false;
}

// ============================================================================
// 基础求值与显示
// ============================================================================

Scalar Calculator::evaluate(const std::string& expression) {
    return normalize_result(evaluate_raw(expression));
}

Scalar Calculator::evaluate_raw(const std::string& expression) {
    const StoredValue value = evaluate_expression_value(impl_.get(), expression, false);
    if (value.is_matrix || value.is_complex) {
        throw std::runtime_error("matrix or complex expression cannot be used as a scalar");
    }
    return value.decimal;
}

std::string Calculator::evaluate_for_display(const std::string& expression, bool exact_mode) {
    apply_calculator_display_precision(impl_.get());

    std::string converted;
    if (try_base_conversion_expression(expression,
                                       impl_->variables_ptr->create_resolver(),
                                       impl_->functions_ptr->get_custom_functions_map(),
                                       {impl_->config_ptr->is_hex_prefix_mode(), impl_->config_ptr->is_hex_uppercase_mode()},
                                       &converted)) {
        return converted;
    }

    const StoredValue value = evaluate_expression_value(impl_.get(), expression, exact_mode);
    return format_stored_value(value, impl_->config_ptr->is_symbolic_constants_mode());
}

std::string Calculator::process_line(const std::string& expression, bool exact_mode) {
    std::string output;
    if (try_process_function_command(expression, &output, exact_mode)) {
        return output;
    }
    return evaluate_for_display(expression, exact_mode);
}

namespace {

std::string execute_script_source(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const std::string& source,
                                  bool exact_mode,
                                  bool suppress_implicit_output) {
    auto is_command = [impl](std::string_view name) {
        return impl->commands_ptr->has_command(std::string(name)) ||
               impl->commands_ptr->has_command(":" + std::string(name));
    };
    script::Program program = script::parse_program(source, is_command);
    std::string accumulated_output;
    std::string last_statement_output;
    for (const auto& statement : program.statements) {
        std::string current_output;
        const ScriptSignal signal =
            execute_script_statement(impl, *statement, exact_mode, &current_output, false);

        // Accumulate output if:
        // 1. It's a print call (recognized by the script engine)
        // 2. It's a "command" (like diff, integral) that produced output
        // 3. It's the very last statement in the program
        bool should_accumulate = false;
        if (statement->kind == script::Statement::Kind::kSimple) {
            const auto& simple = static_cast<const script::SimpleStatement&>(*statement);
            const std::string trimmed = trim_copy(simple.text);
            // Explicit print call
            if (trimmed.compare(0, 6, "print(") == 0 || trimmed == "print") {
                should_accumulate = true;
            }
            // Commands that are NOT assignments but produced output should usually be shown in scripts
            else if (trimmed.find('=') == std::string::npos && !current_output.empty()) {
                should_accumulate = true;
            }
        } else if (statement->kind == script::Statement::Kind::kImport && !current_output.empty()) {
            should_accumulate = true;
        }

        bool is_last = (statement == program.statements.back());
        if (is_last && !suppress_implicit_output) {
            should_accumulate = true;
        }

        if (should_accumulate && !current_output.empty()) {
            if (!accumulated_output.empty()) {
                accumulated_output += "\n";
            }
            accumulated_output += current_output;
        }
        last_statement_output = current_output;

        if (signal.kind == ScriptSignal::Kind::kReturn) {
            std::string return_val = signal.has_value ? format_stored_value(signal.value, impl->config_ptr->is_symbolic_constants_mode())
                                                      : (last_statement_output.empty() ? "OK" : last_statement_output);
            if (accumulated_output.empty()) {
                return return_val;
            } else {
                return accumulated_output + "\n" + return_val;
            }
        }
        if (signal.kind == ScriptSignal::Kind::kBreak ||
            signal.kind == ScriptSignal::Kind::kContinue) {
            throw std::runtime_error("break/continue can only be used inside loops");
        }
    }
    if (!accumulated_output.empty()) return accumulated_output;
    if (suppress_implicit_output) return "";
    return last_statement_output.empty() ? "OK" : last_statement_output;
}

std::filesystem::path resolve_script_file_path(const Calculator::Impl* impl, const std::string& path_text) {
    std::filesystem::path path(path_text);
    if (path.is_relative()) {
        const std::filesystem::path base =
            impl->script_context.file_stack.empty()
                ? std::filesystem::current_path()
                : impl->script_context.file_stack.back().parent_path();
        path = base / path;
    }
    return std::filesystem::weakly_canonical(path);
}

std::string read_script_file(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("unable to open script file: " + path.string());
    }
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

} // namespace

std::string Calculator::execute_script(const std::string& source, bool exact_mode) {
    apply_calculator_display_precision(impl_.get());
    return execute_script_source(this, impl_.get(), source, exact_mode, false);
}

std::string Calculator::execute_script_file(const std::string& path,
                                            bool exact_mode,
                                            bool suppress_implicit_output) {
    apply_calculator_display_precision(impl_.get());

    const std::filesystem::path resolved = resolve_script_file_path(impl_.get(), path);
    if (impl_->script_context.importing_files.find(resolved) != impl_->script_context.importing_files.end()) {
        throw std::runtime_error("circular script import detected: " + resolved.string());
    }

    const std::string source = read_script_file(resolved);
    impl_->script_context.importing_files.insert(resolved);
    impl_->script_context.file_stack.push_back(resolved);
    try {
        const std::string output =
            execute_script_source(this, impl_.get(), source, exact_mode, suppress_implicit_output);
        impl_->script_context.file_stack.pop_back();
        impl_->script_context.importing_files.erase(resolved);
        return output;
    } catch (...) {
        impl_->script_context.file_stack.pop_back();
        impl_->script_context.importing_files.erase(resolved);
        throw;
    }
}

std::string Calculator::list_variables() const {
    apply_calculator_display_precision(impl_.get());

    auto all_vars = impl_->variables_ptr->get_all_globals();
    if (all_vars.empty()) {
        return "No variables defined.";
    }

    std::ostringstream out;
    bool first = true;
    for (const auto& [name, value] : all_vars) {
        // std::map 保证变量按名字稳定排序，便于人读和测试断言。
        if (!first) {
            out << '\n';
        }
        first = false;
        out << name << " = " << format_stored_value(value, impl_->config_ptr->is_symbolic_constants_mode());
    }
    return out.str();
}

std::string Calculator::factor_expression(const std::string& expression) const {
    std::string output;
    if (try_process_function_command(expression, &output, false)) {
        return output;
    }
    throw std::runtime_error("factor command failed or not available");
}

std::string Calculator::plot_expression(const std::string& expression) const {
    std::string output;
    if (try_process_function_command(expression, &output, false)) {
        return output;
    }
    throw std::runtime_error("plot command failed or not available");
}

std::string Calculator::export_variable(const std::string& line) const {
    std::string output;
    if (try_process_function_command(line, &output, false)) {
        return output;
    }
    throw std::runtime_error("export command failed or not available");
}

std::string Calculator::base_conversion_expression(const std::string& expression) const {
    std::string output;
    if (try_process_function_command(expression, &output, false)) {
        return output;
    }

    // 回退到旧实现（处理非函数形式的进制转换）
    std::string converted;
    if (try_base_conversion_expression(expression,
                                        impl_->variables_ptr->create_resolver(),
                                        impl_->functions_ptr->get_custom_functions_map(),
                                        {impl_->config_ptr->is_hex_prefix_mode(), impl_->config_ptr->is_hex_uppercase_mode()},
                                        &converted)) {
        return converted;
    }
    throw std::runtime_error("base conversion failed or not available");
}

// ============================================================================
// 状态持久化
// ============================================================================

std::string Calculator::save_state(const std::string& path) const {
    return impl_->persistence->save_state(path);
}

std::string Calculator::load_state(const std::string& path) {
    return impl_->persistence->load_state(path);
}

// ============================================================================
// IExecutionContext 实现
// ============================================================================

StoredValue Calculator::Impl::evaluate(const std::string& expression, bool exact_mode) {
    auto engine = locator.resolve<IEvaluationEngine>();
    return engine->evaluate_expression_value(expression, exact_mode);
}

bool Calculator::Impl::try_evaluate_implicit(const std::string& expression,
                                            StoredValue* value,
                                            const std::map<std::string, StoredValue>& vars) {
    return parent->try_evaluate_implicit(expression, value, vars);
}

std::string Calculator::Impl::expand_inline(const std::string& expression) {
    return expand_inline_function_commands(this, expression);
}

std::string Calculator::Impl::execute_script_file(const std::string& path,
                                                bool exact_mode,
                                                bool create_scope) {
    return parent->execute_script_file(path, exact_mode, create_scope);
}

bool Calculator::Impl::try_process_function_command(const std::string& expression,
                                                 std::string* output,
                                                 bool exact_mode) const {
    return parent->try_process_function_command(expression, output, exact_mode);
}

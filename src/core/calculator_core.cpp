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

#include "core/calculator_internal_types.h"
#include "parser/unified_expression_parser.h"
#include "parser/command_parser.h"
#include "analysis/function_analysis.h"
#include "matrix/matrix.h"
#include "math/mymath.h"
#include "symbolic/symbolic_expression.h"
#include "core/string_utils.h"
#include "core/format_utils.h"
#include "parser/command_parser.h"
#include "core/calculator_service_factory.h"
#include "script/script_runtime.h"
#include "script/script_parser.h"
#include "module/module_registration.h"
#include "math/helpers/integer_helpers.h"
#include "plot/calculator_plot.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

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
    const int precision = impl == nullptr ? kDefaultDisplayPrecision : impl->display_precision;
    set_process_display_precision(precision);
    matrix::set_display_precision(precision);
    SymbolicExpression::set_display_precision(precision);
}

void broadcast_settings(Calculator* calculator, Calculator::Impl* impl) {
    (void)calculator;
    CalculatorSettings settings;
    settings.display_precision = impl->display_precision;
    settings.exact_mode = false;
    settings.symbolic_constants_mode = impl->symbolic_constants_mode;
    settings.hex_prefix_mode = impl->hex_prefix_mode;
    settings.hex_uppercase_mode = impl->hex_uppercase_mode;

    for (auto& module : impl->registered_modules) {
        module->on_settings_changed(settings);
    }
}

// ============================================================================
// 生命周期与模块注册
// ============================================================================

Calculator::Calculator() : impl_(new Impl()) {
    apply_calculator_display_precision(impl_.get());
    
    // 初始化核心服务缓存
    impl_->core_services = std::make_unique<CoreServices>(core::build_core_services(this, impl_.get()));
    
    register_standard_modules(this);
    broadcast_settings(this, impl_.get());
}

const CoreServices& Calculator::get_core_services() const {
    return *impl_->core_services;
}

Calculator::~Calculator() = default;

void Calculator::register_module(std::shared_ptr<CalculatorModule> module) {
    if (!module) return;

    // 合并函数映射
    auto new_scalar_funcs = module->get_scalar_functions();
    for (auto& [name, func] : new_scalar_funcs) impl_->scalar_functions[name] = std::move(func);
    auto new_matrix_funcs = module->get_matrix_functions();
    for (auto& [name, func] : new_matrix_funcs) impl_->matrix_functions[name] = std::move(func);
    auto new_value_funcs = module->get_value_functions();
    for (auto& [name, func] : new_value_funcs) impl_->value_functions[name] = std::move(func);
    auto new_native_funcs = module->get_native_functions();
    for (auto& [name, func] : new_native_funcs) {
        impl_->native_functions[name] = std::move(func);
        impl_->module_functions.push_back(name);
        impl_->help_topic_to_modules[name].push_back(module);
    }

    // 隐式求值路由优化
    if (module->wants_implicit_evaluation()) {
        impl_->implicit_evaluation_modules.push_back(module);
    }

    // 收集元数据并注册命令到 CommandRegistry
    auto specs = module->get_command_specs();
    for (const auto& spec : specs) {
        std::string cmd_name = command_key_display(spec.key);

        // 检查重复注册
        if (impl_->command_registry.has_command(cmd_name)) {
            throw std::runtime_error("duplicate command registration: " + cmd_name);
        }

        impl_->module_commands.push_back(cmd_name);

        // 保存模块引用用于命令分发
        impl_->command_to_module[spec.key] = {module, spec.dispatch_name};

        // 注册到 CommandRegistry
        // 使用 weak_ptr 避免循环引用
        std::weak_ptr<CalculatorModule> weak_module = module;
        std::string dispatch_name = spec.dispatch_name;
        CommandKey key = spec.key;

        impl_->command_registry.register_command(
            cmd_name,
            [weak_module, dispatch_name, key, this](
                const std::string& /*input*/,
                const std::vector<std::string_view>& args,
                std::string* output,
                bool exact_mode,
                const CoreServices& services) -> bool {
                auto mod = weak_module.lock();
                if (!mod) return false;

                // 使用模块的 execute_args_view 方法处理命令
                *output = mod->execute_args_view(dispatch_name, args, services);
                return true;
            },
            module->get_help_snippet("commands")
        );
    }
    auto funcs = module->get_functions();
    impl_->module_functions.insert(impl_->module_functions.end(), funcs.begin(), funcs.end());

    // 建立帮助索引
    static const std::vector<std::string> topics = {
        "commands", "functions", "matrix", "symbolic", "analysis", "planning",
        "examples", "exact", "variables", "persistence", "programmer"
    };
    for (const auto& topic : topics) {
        if (!module->get_help_snippet(topic).empty()) {
            impl_->help_topic_to_modules[topic].push_back(module);
        }
    }

    impl_->registered_modules.push_back(module);
}

bool is_reserved_user_function_name(const Calculator::Impl* impl, std::string_view name) {
    if (utils::is_reserved_function_name(name)) {
        return true;
    }
    if (impl == nullptr) {
        return false;
    }
    const std::string name_text(name);
    if (impl->command_to_module.find(call_command_key(name)) != impl->command_to_module.end()) {
        return true;
    }
    return impl->scalar_functions.find(name_text) != impl->scalar_functions.end() ||
           impl->matrix_functions.find(name_text) != impl->matrix_functions.end() ||
           impl->value_functions.find(name_text) != impl->value_functions.end();
}

// ============================================================================
// 变量与设置
// ============================================================================

std::string Calculator::clear_variable(const std::string& name) {
    const auto it = impl_->variables.find(name);
    if (it == impl_->variables.end()) throw std::runtime_error("unknown variable: " + name);
    impl_->variables.erase(it);
    return "Cleared variable: " + name;
}

std::string Calculator::clear_all_variables() {
    impl_->variables.clear();
    return "Cleared all variables.";
}

std::string Calculator::set_hex_prefix_mode(bool enabled) {
    impl_->hex_prefix_mode = enabled;
    broadcast_settings(this, impl_.get());
    return std::string("Hex prefix mode: ") + (enabled ? "ON" : "OFF");
}

bool Calculator::hex_prefix_mode() const { return impl_->hex_prefix_mode; }

std::string Calculator::set_hex_uppercase_mode(bool enabled) {
    impl_->hex_uppercase_mode = enabled;
    broadcast_settings(this, impl_.get());
    return std::string("Hex letter case: ") + (enabled ? "UPPER" : "LOWER");
}

bool Calculator::hex_uppercase_mode() const { return impl_->hex_uppercase_mode; }

std::string Calculator::set_symbolic_constants_mode(bool enabled) {
    impl_->symbolic_constants_mode = enabled;
    broadcast_settings(this, impl_.get());
    return std::string("Symbolic constants mode: ") + (enabled ? "ON" : "OFF");
}

bool Calculator::symbolic_constants_mode() const { return impl_->symbolic_constants_mode; }

std::string Calculator::set_display_precision(int precision) {
    if (precision < kMinDisplayPrecision || precision > kMaxDisplayPrecision) {
        throw std::runtime_error("display precision must be in range 1..17");
    }
    impl_->display_precision = precision;
    apply_calculator_display_precision(impl_.get());
    broadcast_settings(this, impl_.get());
    return "Display precision: " + std::to_string(precision);
}

int Calculator::display_precision() const { return impl_->display_precision; }
std::vector<std::string> Calculator::module_command_names() const { return impl_->module_commands; }
std::vector<std::string> Calculator::module_function_names() const { return impl_->module_functions; }
std::vector<std::string> Calculator::variable_names() const {
    std::vector<std::string> names;
    for (const auto& [name, _] : impl_->variables) names.push_back(name);
    return names;
}
std::vector<std::string> Calculator::custom_function_names() const {
    std::vector<std::string> names;
    for (const auto& [name, _] : impl_->functions) names.push_back(name);
    return names;
}

double Calculator::normalize_result(double value) {
    if (!mymath::isfinite(value)) return value;
    if (mymath::abs(value) < kDisplayZeroEps) return 0.0;
    if (mymath::abs(value) > kDisplayIntegerEps && is_integer_double(value, kDisplayIntegerEps)) {
        return static_cast<double>(round_to_long_long(value));
    }
    return value;
}

// ============================================================================
// 命令处理与分发
// ============================================================================

bool Calculator::try_process_function_command(const std::string& expression,
                                              std::string* output, bool exact_mode) {
    auto is_command = [this](std::string_view name) {
        return impl_->command_registry.has_command(std::string(name)) ||
               impl_->command_registry.has_command(":" + std::string(name));
    };

    CommandASTNode root = parse_command(expression, is_command);
    if (root.kind == CommandKind::kEmpty) return false;

    *output = execute_command_ast(this, impl_.get(), root, exact_mode);
    return true;
}

bool Calculator::try_evaluate_implicit(std::string_view expression,
                                       StoredValue* output,
                                       const std::map<std::string, StoredValue>& vars) const {
    if (expression.empty()) return false;

    for (const auto& module : impl_->implicit_evaluation_modules) {
        // 使用模块缓存的触发字符表，避免重复构建
        const auto* trigger_table = module->get_cached_trigger_table();
        if (trigger_table && !has_trigger_char(expression, *trigger_table)) {
            continue;
        }

        if (module->try_evaluate_implicit(std::string(expression), output, vars)) {
            return true;
        }
    }
    return false;
}

// ============================================================================
// 基础求值与显示
// ============================================================================

double Calculator::evaluate(const std::string& expression) {
    return normalize_result(evaluate_raw(expression));
}

double Calculator::evaluate_raw(const std::string& expression) {
    const StoredValue value = evaluate_expression_value(this, impl_.get(), expression, false);
    if (value.is_matrix || value.is_complex) {
        throw std::runtime_error("matrix or complex expression cannot be used as a scalar");
    }
    return value.decimal;
}

std::string Calculator::evaluate_for_display(const std::string& expression, bool exact_mode) {
    apply_calculator_display_precision(impl_.get());

    std::string converted;
    if (try_base_conversion_expression(expression,
                                       VariableResolver(&impl_->variables, nullptr),
                                       &impl_->functions,
                                       {impl_->hex_prefix_mode, impl_->hex_uppercase_mode},
                                       &converted)) {
        return converted;
    }

    const StoredValue value = evaluate_expression_value(this, impl_.get(), expression, exact_mode);
    return format_stored_value(value, impl_->symbolic_constants_mode);
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
        return impl->command_registry.has_command(std::string(name)) ||
               impl->command_registry.has_command(":" + std::string(name));
    };
    script::Program program = script::parse_program(source, is_command);
    std::string accumulated_output;
    std::string last_statement_output;
    for (const auto& statement : program.statements) {
        std::string current_output;
        const ScriptSignal signal =
            execute_script_statement(calculator, impl, *statement, exact_mode, &current_output, false);

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
            std::string return_val = signal.has_value ? format_stored_value(signal.value, impl->symbolic_constants_mode)
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
            impl->script_file_stack.empty()
                ? std::filesystem::current_path()
                : impl->script_file_stack.back().parent_path();
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
    if (impl_->importing_script_files.find(resolved) != impl_->importing_script_files.end()) {
        throw std::runtime_error("circular script import detected: " + resolved.string());
    }

    const std::string source = read_script_file(resolved);
    impl_->importing_script_files.insert(resolved);
    impl_->script_file_stack.push_back(resolved);
    try {
        const std::string output =
            execute_script_source(this, impl_.get(), source, exact_mode, suppress_implicit_output);
        impl_->script_file_stack.pop_back();
        impl_->importing_script_files.erase(resolved);
        return output;
    } catch (...) {
        impl_->script_file_stack.pop_back();
        impl_->importing_script_files.erase(resolved);
        throw;
    }
}

std::string Calculator::list_variables() const {
    apply_calculator_display_precision(impl_.get());

    if (impl_->variables.empty()) {
        return "No variables defined.";
    }

    std::ostringstream out;
    bool first = true;
    for (const auto& [name, value] : impl_->variables) {
        // std::map 保证变量按名字稳定排序，便于人读和测试断言。
        if (!first) {
            out << '\n';
        }
        first = false;
        out << name << " = " << format_stored_value(value, impl_->symbolic_constants_mode);
    }
    return out.str();
}

std::string Calculator::factor_expression(const std::string& expression) const {
    auto is_command = [this](std::string_view name) {
        return impl_->command_registry.has_command(std::string(name));
    };
    CommandASTNode ast = parse_command(expression, is_command);
    const auto* call = ast.as_function_call();
    if (!call || call->name != "factor" || call->arguments.size() != 1) {
        throw std::runtime_error("expected factor(expression)");
    }

    // 先允许 inside 是一个普通表达式或变量，再检查最终值是否为整数。
    const double value = parse_decimal_expression(std::string(call->arguments[0].text), VariableResolver(&impl_->variables, nullptr), &impl_->functions, &impl_->scalar_functions);
    if (!is_integer_double(value)) {
        throw std::runtime_error("factor only accepts integers");
    }

    return factor_integer(round_to_long_long(value));
}

std::string Calculator::plot_expression(const std::string& expression) const {
    CommandASTNode ast = parse_command(expression);
    const auto* call = ast.as_function_call();
    if (!call || call->name != "plot") {
        throw std::runtime_error("expected plot(...)");
    }

    std::vector<std::string> arguments;
    for (const auto& arg : call->arguments) {
        arguments.emplace_back(arg.text);
    }

    plot::PlotContext ctx;
    ctx.variables = visible_variables(impl_.get());
    ctx.functions = &impl_->functions;
    ctx.scalar_functions = &impl_->scalar_functions;
    ctx.has_script_function = [this](const std::string& name) {
        return has_visible_script_function(impl_.get(), name);
    };
    ctx.invoke_script_function = [this](const std::string& name, const std::vector<double>& args) {
        return invoke_script_function_decimal(const_cast<Calculator*>(this), impl_.get(), name, args);
    };

    return plot::handle_plot_command(ctx, arguments);
}

std::string Calculator::export_variable(const std::string& line) const {
    plot::PlotContext ctx;
    ctx.variables = visible_variables(impl_.get());
    ctx.functions = &impl_->functions;
    ctx.scalar_functions = &impl_->scalar_functions;
    ctx.has_script_function = [this](const std::string& name) {
        return has_visible_script_function(impl_.get(), name);
    };
    ctx.invoke_script_function = [this](const std::string& name, const std::vector<double>& args) {
        return invoke_script_function_decimal(const_cast<Calculator*>(this), impl_.get(), name, args);
    };
    return plot::handle_export_command(ctx, line);
}

std::string Calculator::base_conversion_expression(const std::string& expression) const {
    auto is_command = [this](std::string_view name) {
        return impl_->command_registry.has_command(std::string(name));
    };
    CommandASTNode ast = parse_command(expression, is_command);
    const auto* call = ast.as_function_call();
    if (!call) {
        throw std::runtime_error("expected bin(...), oct(...), hex(...), or base(value, base)");
    }

    std::string converted;
    if (!try_base_conversion_expression(expression,
                                        VariableResolver(&impl_->variables, nullptr),
                                        &impl_->functions,
                                        {impl_->hex_prefix_mode, impl_->hex_uppercase_mode},
                                        &converted)) {
        throw std::runtime_error("expected bin(...), oct(...), hex(...), or base(value, base)");
    }
    return converted;
}

// ============================================================================
// 状态持久化
// ============================================================================

std::string Calculator::save_state(const std::string& path) const {
    const std::filesystem::path target_path(path);
    const std::filesystem::path temp_path =
        target_path.parent_path() /
        (target_path.filename().string() + ".tmp-save");
    std::ofstream out(temp_path);
    if (!out) {
        throw std::runtime_error("unable to open file for writing: " + path);
    }

    out << "STATE_V5\n";

    for (const auto& [name, value] : impl_->variables) {
        if (value.is_matrix) {
            out << "VAR\t" << encode_state_field(name)
                << "\tMATRIX\t" << value.matrix.rows
                << '\t' << value.matrix.cols;
            for (double element : value.matrix.data) {
                out << '\t' << std::setprecision(17) << element;
            }
            out << '\n';
        } else if (value.is_complex) {
            out << "VAR\t" << encode_state_field(name)
                << "\tCOMPLEX\t" << std::setprecision(17) << value.complex.real
                << '\t' << std::setprecision(17) << value.complex.imag << '\n';
        } else if (value.is_string) {
            out << "VAR\t" << encode_state_field(name)
                << "\tSTRING\t" << encode_state_field(value.string_value) << '\n';
        } else if (value.exact) {
            out << "VAR\t" << encode_state_field(name)
                << "\tEXACT\t" << value.rational.numerator
                << '\t' << value.rational.denominator
                << '\t' << std::setprecision(17) << value.decimal << '\n';
        } else {
            out << "VAR\t" << encode_state_field(name)
                << "\tDECIMAL\t" << std::setprecision(17) << value.decimal << '\n';
        }
        if (value.has_precise_decimal_text) {
            out << "PRECISE\t" << encode_state_field(name)
                << '\t' << encode_state_field(value.precise_decimal_text) << '\n';
        }
        if (value.has_symbolic_text) {
            out << "SYMBOLIC\t" << encode_state_field(name)
                << '\t' << encode_state_field(value.symbolic_text) << '\n';
        }
    }

    for (const auto& [name, function] : impl_->functions) {
        std::string params_str;
        for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
            if (i != 0) params_str += ", ";
            params_str += function.parameter_names[i];
        }
        out << "EXPRFUNC\t"
            << encode_state_field(name + "(" + params_str + ") = " +
                                  function.expression)
            << '\n';
    }

    for (const auto& [name, function] : impl_->script_functions) {
        std::ostringstream source;
        source << "fn " << name << "(";
        for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
            if (i != 0) {
                source << ", ";
            }
            source << function.parameter_names[i];
        }
        source << ") " << render_script_block(*function.body, 0);
        out << "SCRIPT\t" << encode_state_field(source.str()) << '\n';
    }

    out.close();
    if (!out) {
        std::error_code remove_error;
        std::filesystem::remove(temp_path, remove_error);
        throw std::runtime_error("unable to finish writing state file: " + path);
    }

    std::error_code rename_error;
    std::filesystem::rename(temp_path, target_path, rename_error);
    if (rename_error) {
        std::error_code remove_existing_error;
        std::filesystem::remove(target_path, remove_existing_error);
        rename_error.clear();
        std::filesystem::rename(temp_path, target_path, rename_error);
    }
    if (rename_error) {
        std::error_code remove_error;
        std::filesystem::remove(temp_path, remove_error);
        throw std::runtime_error("unable to replace state file: " + path);
    }

    return "Saved variables to: " + path;
}

std::string Calculator::load_state(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("unable to open file for reading: " + path);
    }

    std::map<std::string, StoredValue> loaded;
    std::map<std::string, CustomFunction> loaded_functions;
    std::map<std::string, ScriptFunction> loaded_script_functions;
    std::string line;
    int state_version = 1;

    auto split_tab_fields = [](const std::string& row_text) {
        std::vector<std::string> parts;
        std::size_t start = 0;
        for (std::size_t i = 0; i <= row_text.size(); ++i) {
            if (i == row_text.size() || row_text[i] == '\t') {
                parts.push_back(row_text.substr(start, i - start));
                start = i + 1;
            }
        }
        return parts;
    };

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line == "STATE_V2") {
            state_version = 2;
            continue;
        }
        if (line == "STATE_V3") {
            state_version = 3;
            continue;
        }
        if (line == "STATE_V4") {
            state_version = 4;
            continue;
        }
        if (line == "STATE_V5") {
            state_version = 5;
            continue;
        }

        const std::vector<std::string> parts = split_tab_fields(line);
        if (state_version >= 2) {
            if (parts.empty()) {
                continue;
            }
            if (parts[0] == "VAR") {
                if (parts.size() < 4) {
                    throw std::runtime_error("invalid save file format");
                }

                StoredValue value;
                const std::string name = decode_state_field(parts[1]);
                if (parts[2] == "STRING") {
                    if (parts.size() != 4) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.is_string = true;
                    value.string_value = decode_state_field(parts[3]);
                } else if (parts[2] == "MATRIX" && state_version >= 4) {
                    if (parts.size() < 5) {
                        throw std::runtime_error("invalid save file format");
                    }
                    const std::size_t rows =
                        static_cast<std::size_t>(std::stoull(parts[3]));
                    const std::size_t cols =
                        static_cast<std::size_t>(std::stoull(parts[4]));
                    if (parts.size() != rows * cols + 5) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.is_matrix = true;
                    value.matrix = matrix::Matrix(rows, cols, 0.0);
                    for (std::size_t i = 0; i < value.matrix.data.size(); ++i) {
                        value.matrix.data[i] = std::stod(parts[i + 5]);
                    }
                } else if (parts[2] == "COMPLEX" && state_version >= 5) {
                    if (parts.size() != 5) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.is_complex = true;
                    value.complex.real = std::stod(parts[3]);
                    value.complex.imag = std::stod(parts[4]);
                } else if (parts[2] == "EXACT") {
                    if (parts.size() != 6) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.exact = true;
                    value.rational = Rational(std::stoll(parts[3]), std::stoll(parts[4]));
                    value.decimal = std::stod(parts[5]);
                } else if (parts[2] == "DECIMAL") {
                    if ((state_version == 2 && parts.size() != 4 && parts.size() != 5) ||
                        (state_version >= 3 && parts.size() != 4)) {
                        throw std::runtime_error("invalid save file format");
                    }
                    value.decimal = std::stod(parts[3]);
                    if (state_version == 2 && parts.size() == 5) {
                        value.has_precise_decimal_text = true;
                        value.precise_decimal_text = decode_state_field(parts[4]);
                    }
                } else {
                    throw std::runtime_error("invalid save file format");
                }
                loaded[name] = value;
                continue;
            }

            if (state_version >= 3 && parts[0] == "PRECISE") {
                if (parts.size() != 3) {
                    throw std::runtime_error("invalid save file format");
                }
                const std::string name = decode_state_field(parts[1]);
                auto it = loaded.find(name);
                if (it == loaded.end() || it->second.is_matrix || it->second.is_complex ||
                    it->second.is_string) {
                    throw std::runtime_error("invalid save file format");
                }
                it->second.has_precise_decimal_text = true;
                it->second.precise_decimal_text = decode_state_field(parts[2]);
                continue;
            }

            if (state_version >= 3 && parts[0] == "SYMBOLIC") {
                if (parts.size() != 3) {
                    throw std::runtime_error("invalid save file format");
                }
                const std::string name = decode_state_field(parts[1]);
                auto it = loaded.find(name);
                if (it == loaded.end() || it->second.is_matrix || it->second.is_complex ||
                    it->second.is_string) {
                    throw std::runtime_error("invalid save file format");
                }
                it->second.has_symbolic_text = true;
                it->second.symbolic_text = decode_state_field(parts[2]);
                continue;
            }

            if (parts[0] == "EXPRFUNC") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                const std::string definition = decode_state_field(parts[1]);
                // 使用 CommandParser 解析函数定义
                CommandASTNode ast = parse_command(definition);
                if (ast.kind == CommandKind::kFunctionDefinition) {
                    const FunctionDefinitionInfo* def = ast.as_function_definition();
                    if (def) {
                        std::vector<std::string> params;
                        for (auto p : def->parameters) {
                            params.emplace_back(p);
                        }
                        loaded_functions[std::string(def->name)] = {params, std::string(def->body.text)};
                    }
                } else {
                    throw std::runtime_error("invalid save file format");
                }
                continue;
            }

            if (parts[0] == "SCRIPT") {
                if (parts.size() != 2) {
                    throw std::runtime_error("invalid save file format");
                }
                Calculator temp;
                temp.execute_script(decode_state_field(parts[1]), false);
                for (const auto& [name, function] : temp.impl_->script_functions) {
                    loaded_script_functions[name] = function;
                }
                continue;
            }

            throw std::runtime_error("invalid save file format");
        }

        if (parts.size() != 5) {
            throw std::runtime_error("invalid save file format");
        }

        StoredValue value;
        value.exact = std::stoi(parts[1]) != 0;
        value.rational = Rational(std::stoll(parts[2]), std::stoll(parts[3]));
        value.decimal = std::stod(parts[4]);
        loaded[parts[0]] = value;
    }

    impl_->variables = loaded;
    if (state_version >= 2) {
        impl_->functions = loaded_functions;
        impl_->script_functions = loaded_script_functions;
    }
    return "Loaded variables from: " + path;
}

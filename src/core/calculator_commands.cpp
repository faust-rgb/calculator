#include "calculator_internal_types.h"
#include "calculator_simplex.h"
#include "calculator_polynomial.h"
#include "calculator_series.h"
#include "calculator_transforms.h"
#include "calculator_rootfinding.h"
#include "calculator_integration.h"
#include "calculator_ode.h"
#include "calculator_optimization.h"
#include "calculator_analysis_cmds.h"
#include "calculator_signal_commands.h"

#include "function_analysis.h"
#include "calculator_matrix_commands.h"
#include "calculator_symbolic_commands.h"
#include "plot/calculator_plot.h"
#include "matrix.h"
#include "multivariable_integrator.h"
#include "mymath.h"
#include "ode_solver.h"
#include "optimization_helpers.h"
#include "polynomial.h"
#include "symbolic_expression.h"
#include "symbolic_expression_internal.h"
#include "system_module.h"
#include "utils.h"

#include "dsp/dsp_module.h"
#include <algorithm>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "../precise/precise_module.h"
#include "../statistics/statistics_module.h"
#include "../math/standard_math_module.h"
#include "../math/integer_math_module.h"
#include "../matrix/matrix_module.h"
#include "../analysis/calculator_series.h"
#include "../analysis/calculator_integration.h"
#include "../analysis/calculator_rootfinding.h"
#include "../analysis/calculator_optimization.h"
#include "../analysis/calculator_analysis_cmds.h"
#include "../analysis/calculator_ode.h"
#include "../symbolic/calculator_symbolic_commands.h"
#include "../symbolic/calculator_transforms.h"
#include "calculator_module.h"

namespace {

} // namespace

void register_standard_modules(Calculator* calculator) {
    // 注册标准数学函数模块
    calculator->register_module(std::make_shared<StandardMathModule>());

    // 注册矩阵函数模块（包含 eig, svd, lu_p 命令）
    calculator->register_module(std::make_shared<MatrixModule>());

    // 注册 DSP 模块（包含 residue 命令）
    calculator->register_module(std::make_shared<DspModule>());

    // 注册系统管理模块
    // 注册核心功能模块
    calculator->register_module(std::make_shared<SystemModule>());
    calculator->register_module(std::make_shared<polynomial_ops::PolynomialModule>());
    calculator->register_module(std::make_shared<series_ops::SeriesModule>());
    calculator->register_module(std::make_shared<transforms::TransformModule>());
    calculator->register_module(std::make_shared<integration_ops::IntegrationModule>());
    calculator->register_module(std::make_shared<rootfinding::RootfindingModule>());
    calculator->register_module(std::make_shared<optimization::OptimizationModule>());
    calculator->register_module(std::make_shared<symbolic_commands::SymbolicModule>());
    calculator->register_module(std::make_shared<analysis_cmds::AnalysisModule>());
    calculator->register_module(std::make_shared<ode_ops::ODEModule>());

    // 注册其他辅助模块
    calculator->register_module(std::make_shared<StatisticsModule>());
    calculator->register_module(std::make_shared<IntegerMathModule>());
    calculator->register_module(std::make_shared<PreciseModule>());
}
bool Calculator::try_process_function_command(const std::string& expression,
                                              std::string* output, bool exact_mode) {
    // 1. 处理自定义函数定义 f(x) = ...
    std::string function_name;
    std::string parameter_name;
    std::string body;
    if (split_function_definition(expression, &function_name, &parameter_name, &body)) {
        if (is_reserved_function_name(function_name)) {
            throw std::runtime_error("function name is reserved: " + function_name);
        }
        impl_->functions[function_name] = {parameter_name, body};
        *output = function_name + "(" + parameter_name + ") = " + body;
        return true;
    }

    const std::string trimmed = trim_copy(expression);
    if (trimmed.empty()) return false;

    // 3. 尝试解析为命令名和括号内容
    std::string command_name;
    std::string inside;
    const std::size_t open_paren = trimmed.find('(');
    bool is_command_style = (open_paren != std::string::npos && trimmed.back() == ')');

    // 检查元命令
    bool is_meta_command = (trimmed.front() == ':');

    if (is_command_style || is_meta_command) {
        if (is_command_style) {
            command_name = trim_copy(trimmed.substr(0, open_paren));
            inside = trimmed.substr(open_paren + 1, trimmed.size() - open_paren - 2);
        } else {
            const std::size_t space = trimmed.find(' ');
            if (space != std::string::npos) {
                command_name = trimmed.substr(0, space);
                inside = trimmed.substr(space + 1);
            } else {
                command_name = trimmed;
                inside = "";
            }
        }

        // 统一参数拆分
        std::vector<std::string> arguments = split_top_level_arguments(inside);

        // 指令别名解析
        if (command_name == ":plot") {
            command_name = "plot";
            arguments.insert(arguments.begin(), "__gnuplot__");
        }

        // 5. 构造 CoreServices (分项构造)
        auto services_factory = [&]() {
            CoreServices s;
            
            // Evaluation Service
            s.evaluation.parse_decimal = [this](const std::string& arg) {
                DecimalParser parser(arg, VariableResolver(&impl_->variables, nullptr), &impl_->functions, &impl_->scalar_functions);
                return parser.parse();
            };
            s.evaluation.evaluate_value = [this](const std::string& arg, bool exact) {
                return evaluate_expression_value(this, impl_.get(), arg, exact);
            };
            s.evaluation.normalize_result = [](double v) { return Calculator::normalize_result(v); };
            s.evaluation.build_decimal_evaluator = [this, variables = VariableResolver::make_owned(visible_variables(impl_.get()))](const std::string& arg) {
                const std::string scoped_expression = trim_copy(expand_inline_function_commands(this, arg));
                return [this, scoped_expression, variables](const std::vector<std::pair<std::string, double>>& assignments) {
                    std::map<std::string, StoredValue> override_vars;
                    for (const auto& [name, value] : assignments) {
                        StoredValue stored;
                        stored.decimal = normalize_display_decimal(value);
                        stored.exact = false;
                        override_vars[name] = stored;
                    }
                    const HasScriptFunctionCallback has_script_function = [this](const std::string& name) { return has_visible_script_function(impl_.get(), name); };
                    const InvokeScriptFunctionDecimalCallback invoke_script_function = [this](const std::string& name, const std::vector<double>& args) { return invoke_script_function_decimal(this, impl_.get(), name, args); };
                    VariableResolver chained_resolver(nullptr, nullptr, &override_vars, &variables);
                    DecimalParser parser(scoped_expression, chained_resolver, &impl_->functions, &impl_->scalar_functions, has_script_function, invoke_script_function);
                    return Calculator::normalize_result(parser.parse());
                };
            };
            s.evaluation.build_scalar_evaluator = [this](const std::string& arg) {
                const std::string scoped_expression = trim_copy(expand_inline_function_commands(this, arg));
                return [this, scoped_expression](const std::vector<std::pair<std::string, StoredValue>>& assignments) {
                    std::map<std::string, StoredValue> frame;
                    for (const auto& [name, value] : assignments) frame[name] = value;
                    impl_->local_scopes.push_back(frame);
                    try {
                        const StoredValue value = evaluate_expression_value(this, impl_.get(), scoped_expression, false);
                        impl_->local_scopes.pop_back();
                        if (value.is_matrix || value.is_complex || value.is_string) throw std::runtime_error("expected a scalar-valued expression");
                        return Calculator::normalize_result(value.exact ? rational_to_double(value.rational) : value.decimal);
                    } catch (...) { impl_->local_scopes.pop_back(); throw; }
                };
            };
            s.evaluation.build_matrix_evaluator = [this](const std::string& arg) {
                const std::string scoped_expression = trim_copy(expand_inline_function_commands(this, arg));
                return [this, scoped_expression](const std::vector<std::pair<std::string, StoredValue>>& assignments) {
                    std::map<std::string, StoredValue> scoped_variables = visible_variables(impl_.get()).snapshot();
                    for (const auto& [name, value] : assignments) scoped_variables[name] = value;
                    const HasScriptFunctionCallback has_script_function = [this](const std::string& name) { return has_visible_script_function(impl_.get(), name); };
                    const InvokeScriptFunctionDecimalCallback invoke_script_function = [this](const std::string& name, const std::vector<double>& args) { return invoke_script_function_decimal(this, impl_.get(), name, args); };
                    matrix::Value val;
                    if (!try_evaluate_matrix_expression(scoped_expression, VariableResolver(&scoped_variables, nullptr), &impl_->functions, &impl_->scalar_functions, &impl_->matrix_functions, &impl_->value_functions, has_script_function, invoke_script_function, &val) || !val.is_matrix)
                        throw std::runtime_error("expected a matrix-valued expression");
                    return val.matrix;
                };
            };

            // Symbolic Service
            s.symbolic.resolve_symbolic = [this](const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) {
                symbolic_commands::SymbolicResolverContext symbolic_resolver_ctx;
                symbolic_resolver_ctx.resolve_custom_function = [this](const std::string& name, std::string* v) {
                    const auto it = impl_->functions.find(name);
                    if (it == impl_->functions.end()) throw std::runtime_error("unknown custom function: " + name);
                    *v = it->second.parameter_name;
                    return SymbolicExpression::parse(it->second.expression);
                };
                symbolic_resolver_ctx.has_custom_function = [this](const std::string& name) {
                    return impl_->functions.find(name) != impl_->functions.end();
                };
                symbolic_resolver_ctx.expand_inline = [this](const std::string& a) {
                    return expand_inline_function_commands(this, a);
                };
                symbolic_commands::resolve_symbolic_expression(symbolic_resolver_ctx, arg, req, var, expr);
            };
            s.symbolic.expand_inline = [this](const std::string& arg) { return expand_inline_function_commands(this, arg); };
            s.symbolic.simplify_symbolic = [](const std::string& text) { return SymbolicExpression::parse(text).simplify().to_string(); };
            s.symbolic.evaluate_symbolic_at = [this](const SymbolicExpression& expr, const std::string& var, double p) {
                const auto existing = impl_->variables.find(var);
                const bool had_existing = existing != impl_->variables.end();
                StoredValue backup;
                if (had_existing) backup = existing->second;
                StoredValue temporary;
                temporary.decimal = p;
                temporary.exact = false;
                impl_->variables[var] = temporary;
                auto cleanup = [&]() {
                    if (had_existing) impl_->variables[var] = backup;
                    else impl_->variables.erase(var);
                };
                try {
                    const double value = this->evaluate(expr.to_string());
                    cleanup();
                    if (!mymath::isfinite(value)) throw std::runtime_error("Non-finite value");
                    return value;
                } catch (...) {
                    cleanup();
                    try {
                        FunctionAnalysis analysis(var);
                        analysis.define(expr.to_string());
                        return analysis.limit(p, 1);
                    } catch (...) { return mymath::quiet_nan(); }
                }
            };
            s.symbolic.parse_symbolic_expr_list = [this](const std::string& arg) { return symbolic_commands::parse_symbolic_expression_list(arg, [this](const std::string& a) { return expand_inline_function_commands(this, a); }); };

            // Environment Service
            s.env.has_variable = [this](const std::string& name) { return impl_->variables.find(name) != impl_->variables.end(); };
            s.env.has_function = [this](const std::string& name) { return impl_->functions.find(name) != impl_->functions.end(); };
            s.env.list_variables = [this]() { return this->list_variables(); };
            s.env.list_functions = [this]() {
                if (impl_->functions.empty() && impl_->script_functions.empty()) return std::string("No custom functions defined.");
                std::ostringstream out;
                bool first = true;
                for (const auto& [name, function] : impl_->functions) {
                    if (!first) out << '\n';
                    first = false;
                    out << name << "(" << function.parameter_name << ") = " << function.expression;
                }
                for (const auto& [name, function] : impl_->script_functions) {
                    if (!first) out << '\n';
                    first = false;
                    out << name << "(";
                    for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
                        if (i != 0) out << ", ";
                        out << function.parameter_names[i];
                    }
                    out << ") = { ... }";
                }
                return out.str();
            };
            s.env.clear_variable = [this](const std::string& name) { return this->clear_variable(name); };
            s.env.clear_function = [this](const std::string& name) {
                const auto simple_it = impl_->functions.find(name);
                if (simple_it != impl_->functions.end()) {
                    impl_->functions.erase(simple_it);
                    return std::string("Cleared custom function: ") + name;
                }
                const auto script_it = impl_->script_functions.find(name);
                if (script_it != impl_->script_functions.end()) {
                    impl_->script_functions.erase(script_it);
                    return std::string("Cleared custom function: ") + name;
                }
                throw std::runtime_error("unknown custom function: " + name);
            };
            s.env.clear_all_variables = [this]() { return this->clear_all_variables(); };
            s.env.save_state = [this](const std::string& p) { return this->save_state(p); };
            s.env.load_state = [this](const std::string& p) { return this->load_state(p); };
            s.env.export_variable = [this](const std::string& p) { return this->export_variable(p); };
            s.env.execute_script = [this](const std::string& c, bool e) { return this->execute_script(c, e); };
            s.env.clear_all_functions = [this]() {
                impl_->functions.clear();
                impl_->script_functions.clear();
                return std::string("Cleared all custom functions.");
            };
            s.env.set_exact_mode = [this](bool m) { return std::string("Exact mode: ") + (m ? "ON" : "OFF"); };
            s.env.set_symbolic_mode = [this](bool m) { return this->set_symbolic_constants_mode(m); };
            s.env.set_precision = [this](int p) { return this->set_display_precision(p); };
            s.env.set_hex_prefix = [this](bool m) { return this->set_hex_prefix_mode(m); };
            s.env.set_hex_uppercase = [this](bool m) { return this->set_hex_uppercase_mode(m); };

            // Global Helpers
            s.parse_symbolic_vars = [](const std::vector<std::string>& args, std::size_t start, const std::vector<std::string>& fallback) { return symbolic_commands::parse_symbolic_variable_arguments(args, start, fallback); };
            s.is_matrix_argument = [this](const std::string& arg) {
                const VariableResolver visible = visible_variables(impl_.get());
                const HasScriptFunctionCallback has_script_function = [this](const std::string& name) { return has_visible_script_function(impl_.get(), name); };
                const InvokeScriptFunctionDecimalCallback invoke_script_function = [this](const std::string& name, const std::vector<double>& args) { return invoke_script_function_decimal(this, impl_.get(), name, args); };
                matrix::Value value;
                return try_evaluate_matrix_expression(trim_copy(arg), visible, &impl_->functions, &impl_->scalar_functions, &impl_->matrix_functions, &impl_->value_functions, has_script_function, invoke_script_function, &value) && value.is_matrix;
            };
            s.parse_matrix_argument = [this](const std::string& arg, const std::string& c) {
                const StoredValue value = evaluate_expression_value(this, impl_.get(), arg, false);
                if (!value.is_matrix) throw std::runtime_error(c + " expects a matrix or vector argument");
                return value.matrix;
            };
            s.is_integer_double = [](double x, double eps) { return is_integer_double(x, eps); };
            s.round_to_long_long = [](double x) { return round_to_long_long(x); };

            return s;
        };

        CoreServices svc = services_factory();

        // 5.1 设置需要引用 svc 自身的成员
        svc.symbolic.build_analysis = [this, &svc](const std::string& argument) {
            std::string variable_name;
            SymbolicExpression expression;
            symbolic_commands::SymbolicResolverContext symbolic_resolver_ctx;
            symbolic_resolver_ctx.resolve_custom_function = [this](const std::string& name, std::string* v) {
                const auto it = impl_->functions.find(name);
                if (it == impl_->functions.end()) throw std::runtime_error("unknown custom function: " + name);
                *v = it->second.parameter_name;
                return SymbolicExpression::parse(it->second.expression);
            };
            symbolic_resolver_ctx.has_custom_function = [this](const std::string& name) { return impl_->functions.find(name) != impl_->functions.end(); };
            symbolic_resolver_ctx.expand_inline = [this](const std::string& a) { return expand_inline_function_commands(this, a); };
            symbolic_commands::resolve_symbolic_expression(symbolic_resolver_ctx, argument, true, &variable_name, &expression);
            
            FunctionAnalysis analysis(variable_name);
            analysis.define(expression.to_string());
            analysis.set_evaluator(svc.evaluation.build_decimal_evaluator(expression.to_string()));
            return analysis;
        };

        // 5.2 确保模块已初始化 (简单起见每次调用都注入最新 svc)
        for (auto& module : impl_->registered_modules) {
            module->initialize(svc);
        }

        // 6. 优先通过命令映射快速分发
        const auto cmd_it = impl_->command_to_module.find(command_name);
        if (cmd_it != impl_->command_to_module.end()) {
            *output = cmd_it->second->execute_args(command_name, arguments, svc);
            return true;
        }

        // 7. 轮询注册模块 (作为回退)
        for (const auto& module : impl_->registered_modules) {
            if (module->can_handle(command_name)) {
                *output = module->execute_args(command_name, arguments, svc);
                return true;
            }
        }
    }

    // 7. 回退到标准表达式求值
    *output = this->process_line(trimmed, exact_mode);
    return true;
}

bool Calculator::try_evaluate_implicit(const std::string& expression,
                                       StoredValue* output,
                                       const std::map<std::string, StoredValue>& vars) const {
    if (expression.empty()) return false;
    
    for (const auto& module : impl_->implicit_evaluation_modules) {
        const std::string triggers = module->get_implicit_trigger_chars();
        if (!triggers.empty()) {
            bool matched = false;
            for (char c : expression) {
                if (triggers.find(c) != std::string::npos) {
                    matched = true;
                    break;
                }
            }
            if (!matched) continue;
        }

        if (module->try_evaluate_implicit(expression, output, vars)) {
            return true;
        }
    }
    return false;
}

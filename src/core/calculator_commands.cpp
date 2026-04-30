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

#include "dsp/residue.h"
#include <algorithm>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "../precise/precise_module.h"
#include "../statistics/statistics_module.h"
#include "calculator_module.h"
#include <set>

namespace {

/**
 * @class BridgeModule
 * @brief 桥接模块，用于将现有的模块处理逻辑包装成 CalculatorModule 接口
 */
class BridgeModule : public CalculatorModule {
public:
    using Matcher = std::function<bool(const std::string&)>;
    using Handler = std::function<bool(const std::string&, const std::string&, std::string*, const CoreServices&)>;

    BridgeModule(std::string name, Matcher matcher, Handler handler)
        : name_(std::move(name)), matcher_(std::move(matcher)), handler_(std::move(handler)) {}

    std::string name() const override { return name_; }
    bool can_handle(const std::string& command) const override { return matcher_(command); }
    std::string execute(const std::string& command, const std::string& inside, const CoreServices& services) override {
        std::string output;
        if (handler_(command, inside, &output, services)) {
            return output;
        }
        throw std::runtime_error("Command failed in module: " + name_);
    }

private:
    std::string name_;
    Matcher matcher_;
    Handler handler_;
};

} // namespace

void register_standard_modules(Calculator* calculator, Calculator::Impl* impl) {
    // 注册各个领域的 BridgeModule
    calculator->register_module(std::make_shared<BridgeModule>("Polynomial", 
        polynomial_ops::is_polynomial_command,
        [impl](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices&) {
            polynomial_ops::PolynomialContext ctx;
            ctx.functions = &impl->functions;
            ctx.resolve_symbolic = [&](const std::string& name, std::string* var) {
                const auto it = impl->functions.find(name);
                if (it == impl->functions.end()) throw std::runtime_error("unknown custom function: " + name);
                *var = it->second.parameter_name;
                return SymbolicExpression::parse(it->second.expression);
            };
            return polynomial_ops::handle_polynomial_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("Series",
        series_ops::is_series_command,
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            series_ops::SeriesContext ctx;
            ctx.resolve_symbolic = svc.resolve_symbolic;
            ctx.parse_decimal = svc.parse_decimal;
            ctx.evaluate_at = svc.evaluate_symbolic_at;
            ctx.simplify_symbolic = svc.simplify_symbolic;
            ctx.expand_inline = svc.expand_inline;
            return series_ops::handle_series_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("Transforms",
        transforms::is_transform_command,
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            transforms::TransformContext ctx;
            ctx.resolve_symbolic = svc.resolve_symbolic;
            return transforms::handle_transform_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("Integration",
        integration_ops::is_integration_command,
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            integration_ops::IntegrationContext ctx;
            ctx.parse_decimal = svc.parse_decimal;
            ctx.build_scoped_evaluator = svc.build_decimal_evaluator;
            ctx.normalize_result = svc.normalize_result;
            return integration_ops::handle_integration_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("Rootfinding",
        rootfinding::is_rootfinding_command,
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            rootfinding::RootfindingContext ctx;
            ctx.parse_decimal = svc.parse_decimal;
            ctx.build_scoped_evaluator = svc.build_decimal_evaluator;
            ctx.get_derivative_expression = [&](const std::string& expr_str, const std::string& var_name) {
                try {
                    std::string var;
                    SymbolicExpression expr;
                    svc.resolve_symbolic(expr_str, false, &var, &expr);
                    if (expr.node_) return expr.derivative(var_name).simplify().to_string();
                } catch (...) {}
                return std::string();
            };
            ctx.is_matrix_argument = svc.is_matrix_argument;
            ctx.parse_matrix_argument = svc.parse_matrix_argument;
            ctx.normalize_result = svc.normalize_result;
            return rootfinding::handle_rootfinding_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("Optimization",
        optimization::is_optimization_command,
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            optimization::OptimizationContext ctx;
            ctx.parse_matrix_argument = svc.parse_matrix_argument;
            ctx.normalize_result = svc.normalize_result;
            ctx.is_integer_double = svc.is_integer_double;
            ctx.round_to_long_long = svc.round_to_long_long;
            return optimization::handle_optimization_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("Symbolic",
        symbolic_commands::is_symbolic_command,
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            symbolic_commands::SymbolicCommandContext ctx;
            ctx.resolve_symbolic = svc.resolve_symbolic;
            ctx.parse_symbolic_variable_arguments = svc.parse_symbolic_vars;
            ctx.parse_symbolic_expression_list = svc.parse_symbolic_expr_list;
            ctx.build_analysis = svc.build_analysis;
            ctx.build_scoped_evaluator = svc.build_decimal_evaluator;
            ctx.parse_decimal = svc.parse_decimal;
            ctx.normalize_result = svc.normalize_result;
            return symbolic_commands::handle_symbolic_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("Analysis",
        analysis_cmds::is_analysis_command,
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            analysis_cmds::AnalysisContext ctx;
            ctx.resolve_symbolic = svc.resolve_symbolic;
            ctx.parse_symbolic_variable_arguments = svc.parse_symbolic_vars;
            ctx.parse_decimal = svc.parse_decimal;
            ctx.normalize_result = svc.normalize_result;
            ctx.build_analysis = svc.build_analysis;
            return analysis_cmds::handle_analysis_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("ODE",
        ode_ops::is_ode_command,
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            ode_ops::ODEContext ctx;
            ctx.parse_decimal = svc.parse_decimal;
            ctx.build_scoped_scalar_evaluator = svc.build_scalar_evaluator;
            ctx.build_scoped_matrix_evaluator = svc.build_matrix_evaluator;
            ctx.is_matrix_argument = svc.is_matrix_argument;
            ctx.parse_matrix_argument = svc.parse_matrix_argument;
            ctx.evaluate_expression_value = svc.evaluate_value;
            ctx.normalize_result = svc.normalize_result;
            return ode_ops::handle_ode_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("Matrix",
        matrix_commands::is_matrix_command,
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            matrix_commands::MatrixCommandContext ctx;
            ctx.is_matrix_argument = svc.is_matrix_argument;
            ctx.parse_matrix_argument = svc.parse_matrix_argument;
            return matrix_commands::handle_matrix_command(ctx, cmd, inside, out);
        }
    ));

    calculator->register_module(std::make_shared<BridgeModule>("Signal",
        signal_cmds::is_signal_command,
        [impl](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            signal_cmds::SignalContext ctx;
            ctx.functions = &impl->functions;
            ctx.resolve_scalar = [&](const std::string& arg, std::string* err) {
                try { return svc.parse_decimal(arg); }
                catch (...) { *err = "Failed to parse scalar: " + arg; return 0.0; }
            };
            ctx.resolve_signal = [&](const std::string& arg, std::string* err) {
                signal_cmds::SignalData data;
                std::string trimmed = trim_copy(arg);
                if (trimmed.front() == '[') {
                    std::string inner = trimmed.substr(1, trimmed.size() - 2);
                    std::stringstream ss(inner);
                    std::string token;
                    while (std::getline(ss, token, ',')) {
                        token = trim_copy(token);
                        if (!token.empty()) {
                            try { data.samples.push_back(std::stod(token)); }
                            catch (...) { *err = "Failed to parse vector element: " + token; }
                        }
                    }
                }
                return data;
            };
            return signal_cmds::handle_signal_command(ctx, cmd, inside, out);
        }
    ));

    // 注册绘图模块
    calculator->register_module(std::make_shared<BridgeModule>("Plot",
        [](const std::string& cmd) {
            return cmd == "plot" || cmd == "imshow" || cmd == "bar" || cmd == "hist";
        },
        [calculator, impl](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices&) {
            std::vector<std::string> arguments = split_top_level_arguments(inside);
            plot::PlotContext plot_ctx;
            plot_ctx.variables = visible_variables(impl);
            plot_ctx.functions = &impl->functions;
            plot_ctx.has_script_function = [impl](const std::string& name) {
                return has_visible_script_function(impl, name);
            };
            plot_ctx.invoke_script_function = [calculator, impl](const std::string& name, const std::vector<double>& args) {
                return invoke_script_function_decimal(calculator, impl, name, args);
            };

            if (cmd == "plot") {
                *out = plot::handle_plot_command(plot_ctx, arguments);
            } else if (cmd == "imshow") {
                *out = plot::handle_imshow_command(plot_ctx, arguments);
            } else if (cmd == "bar") {
                *out = plot::handle_bar_command(plot_ctx, arguments);
            } else if (cmd == "hist") {
                *out = plot::handle_hist_command(plot_ctx, arguments);
            } else {
                return false;
            }
            return true;
        }
    ));

    // 注册全新的独立统计模块
    calculator->register_module(std::make_shared<StatisticsModule>());
    // 注册基础数学模块
    calculator->register_module(std::make_shared<BridgeModule>("Basic",
        [](const std::string& cmd) {
            return cmd == "factor";
        },
        [](const std::string& cmd, const std::string& inside, std::string* out, const CoreServices& svc) {
            if (cmd == "factor") {
                double val = svc.parse_decimal(inside);
                if (!mymath::is_near_zero(val - std::round(val), 1e-10)) {
                    throw std::runtime_error("factor only accepts integers");
                }
                *out = factor_integer(static_cast<long long>(std::round(val)));
                return true;
            }
            return false;
        }
    ));

    calculator->register_module(std::make_shared<PreciseModule>());
}
bool Calculator::try_process_function_command(const std::string& expression,
                                              std::string* output) {
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

    // 2. 处理元命令
    if (trimmed == ":funcs") {
        if (impl_->functions.empty() && impl_->script_functions.empty()) {
            *output = "No custom functions defined.";
            return true;
        }
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
        *output = out.str();
        return true;
    }

    if (trimmed == ":clearfuncs") {
        impl_->functions.clear();
        impl_->script_functions.clear();
        *output = "Cleared all custom functions.";
        return true;
    }

    if (trimmed.rfind(":clearfunc ", 0) == 0) {
        const std::string name = trim_copy(trimmed.substr(11));
        if (impl_->functions.erase(name) > 0 || impl_->script_functions.erase(name) > 0) {
            *output = "Cleared custom function: " + name;
            return true;
        }
        throw std::runtime_error("unknown custom function: " + name);
    }

    // 3. 解析命令名和括号内容
    std::string command_name;
    std::string inside;
    const std::size_t open_paren = trimmed.find('(');
    if (open_paren == std::string::npos || trimmed.back() != ')') return false;
    
    command_name = trim_copy(trimmed.substr(0, open_paren));
    inside = trimmed.substr(open_paren + 1, trimmed.size() - open_paren - 2);

    if (command_name == "residue") {
        *output = dsp_ops::handle_residue_command(this, impl_.get(), split_top_level_arguments(inside));
        return true;
    }

    // 5. 构造 CoreServices (按需构造以提高效率，或在此处构造一次)
    auto services_factory = [&]() {
        CoreServices s;
        s.parse_decimal = [this](const std::string& arg) {
            DecimalParser parser(arg, VariableResolver(&impl_->variables, nullptr), &impl_->functions);
            return parser.parse();
        };
        s.evaluate_value = [this](const std::string& arg, bool exact) {
            return evaluate_expression_value(this, impl_.get(), arg, exact);
        };
        s.normalize_result = [](double v) { return Calculator::normalize_result(v); };
        s.resolve_symbolic = [this](const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) {
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
        s.expand_inline = [this](const std::string& arg) { return expand_inline_function_commands(this, arg); };
        s.simplify_symbolic = [](const std::string& text) { return SymbolicExpression::parse(text).simplify().to_string(); };
        s.evaluate_symbolic_at = [this](const SymbolicExpression& expr, const std::string& var, double p) {
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
        s.build_decimal_evaluator = [this, variables = VariableResolver::make_owned(visible_variables(impl_.get()))](const std::string& arg) {
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
                
                // Use the captured 'variables' as parent
                VariableResolver chained_resolver(nullptr, nullptr, &override_vars, &variables);
                
                DecimalParser parser(scoped_expression, chained_resolver, &impl_->functions, has_script_function, invoke_script_function);
                return Calculator::normalize_result(parser.parse());
            };
        };
        s.build_scalar_evaluator = [this](const std::string& arg) {
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
        s.build_matrix_evaluator = [this](const std::string& arg) {
            const std::string scoped_expression = trim_copy(expand_inline_function_commands(this, arg));
            return [this, scoped_expression](const std::vector<std::pair<std::string, StoredValue>>& assignments) {
                std::map<std::string, StoredValue> scoped_variables = visible_variables(impl_.get()).snapshot();
                for (const auto& [name, value] : assignments) scoped_variables[name] = value;
                const HasScriptFunctionCallback has_script_function = [this](const std::string& name) { return has_visible_script_function(impl_.get(), name); };
                const InvokeScriptFunctionDecimalCallback invoke_script_function = [this](const std::string& name, const std::vector<double>& args) { return invoke_script_function_decimal(this, impl_.get(), name, args); };
                matrix::Value val;
                if (!try_evaluate_matrix_expression(scoped_expression, VariableResolver(&scoped_variables, nullptr), &impl_->functions, has_script_function, invoke_script_function, &val) || !val.is_matrix)
                    throw std::runtime_error("expected a matrix-valued expression");
                return val.matrix;
            };
        };
        s.parse_symbolic_vars = [](const std::vector<std::string>& args, std::size_t start, const std::vector<std::string>& fallback) { return symbolic_commands::parse_symbolic_variable_arguments(args, start, fallback); };
        s.parse_symbolic_expr_list = [this](const std::string& arg) { return symbolic_commands::parse_symbolic_expression_list(arg, [this](const std::string& a) { return expand_inline_function_commands(this, a); }); };
        s.is_matrix_argument = [this](const std::string& arg) {
            const VariableResolver visible = visible_variables(impl_.get());
            const HasScriptFunctionCallback has_script_function = [this](const std::string& name) { return has_visible_script_function(impl_.get(), name); };
            const InvokeScriptFunctionDecimalCallback invoke_script_function = [this](const std::string& name, const std::vector<double>& args) { return invoke_script_function_decimal(this, impl_.get(), name, args); };
            matrix::Value value;
            return try_evaluate_matrix_expression(trim_copy(arg), visible, &impl_->functions, has_script_function, invoke_script_function, &value) && value.is_matrix;
        };
        s.parse_matrix_argument = [this](const std::string& arg, const std::string& c) {
            const StoredValue value = evaluate_expression_value(this, impl_.get(), arg, false);
            if (!value.is_matrix) throw std::runtime_error(c + " expects a matrix or vector argument");
            return value.matrix;
        };
        s.is_integer_double = [](double x, double eps) { return is_integer_double(x, eps); };
        s.round_to_long_long = [](double x) { return round_to_long_long(x); };
        s.has_variable = [this](const std::string& name) { return impl_->variables.find(name) != impl_->variables.end(); };
        s.has_function = [this](const std::string& name) { return impl_->functions.find(name) != impl_->functions.end(); };
        return s;
    };

    CoreServices svc = services_factory();

    // 5.1 设置需要引用 svc 自身的成员
    svc.build_analysis = [this, &svc](const std::string& argument) {
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
        analysis.set_evaluator(svc.build_decimal_evaluator(expression.to_string()));
        return analysis;
    };

    // 6. 轮询注册模块
    for (const auto& module : impl_->registered_modules) {
        if (module->can_handle(command_name)) {
            *output = module->execute(command_name, inside, svc);
            return true;
        }
    }

    return false;
}

// ============================================================================
// 核心服务工厂实现
// ============================================================================
//
// 构建计算器核心提供的所有服务接口，包括：
// - 求值服务（表达式求值、作用域求值器构建）
// - 符号服务（符号表达式解析、简化、求值）
// - 环境服务（变量/函数管理、状态持久化）
// ============================================================================

#include "core/calculator_service_factory.h"
#include "core/calculator_internal_types.h"
#include "math/helpers/integer_helpers.h"
#include "math/helpers/combinatorics.h"
#include "math/helpers/bitwise_helpers.h"
#include "math/helpers/unit_conversions.h"
#include "math/helpers/base_conversions.h"
#include "command/variable_resolver.h"
#include "parser/unified_expression_parser.h"
#include "script/script_runtime.h"
#include "symbolic/calculator_symbolic_commands.h"
#include "symbolic/symbolic_expression.h"
#include "analysis/function_analysis.h"
#include "plot/calculator_plot.h"
#include "core/string_utils.h"
#include "core/format_utils.h"
#include "statistics/statistics.h"
#include "statistics/probability.h"
#include <sstream>

namespace core {

/// 构建核心服务对象
CoreServices build_core_services(Calculator* calculator, Calculator::Impl* impl) {
    CoreServices s;

    // Evaluation Service
    s.evaluation.parse_decimal = [calculator, impl](const std::string& arg) {
        return parse_decimal_expression(arg, VariableResolver(&impl->variables, nullptr), &impl->functions, &impl->scalar_functions);
    };
    s.evaluation.evaluate_value = [calculator, impl](const std::string& arg, bool exact) {
        return evaluate_expression_value(calculator, impl, arg, exact);
    };
    s.evaluation.normalize_result = [](double v) { return Calculator::normalize_result(v); };
    
    s.evaluation.build_decimal_evaluator = [calculator, impl, variables = VariableResolver::make_owned(visible_variables(impl))](const std::string& arg) {
        const std::string scoped_expression = trim_copy(expand_inline_function_commands(calculator, arg));
        return [calculator, impl, scoped_expression, variables](const std::vector<std::pair<std::string, double>>& assignments) {
            std::map<std::string, StoredValue> override_vars;
            for (const auto& [name, value] : assignments) {
                StoredValue stored;
                stored.decimal = normalize_display_decimal(value);
                stored.exact = false;
                override_vars[name] = stored;
            }
            const HasScriptFunctionCallback has_script_function = [calculator, impl](const std::string& name) { return has_visible_script_function(impl, name); };
            const InvokeScriptFunctionDecimalCallback invoke_script_function = [calculator, impl](const std::string& name, const std::vector<double>& args) { return invoke_script_function_decimal(calculator, impl, name, args); };
            VariableResolver chained_resolver(nullptr, nullptr, &override_vars, &variables);
            return parse_decimal_expression(scoped_expression, chained_resolver, &impl->functions, &impl->scalar_functions, has_script_function, invoke_script_function);
        };
    };
    
    s.evaluation.build_scalar_evaluator = [calculator, impl](const std::string& arg) {
        const std::string scoped_expression = trim_copy(expand_inline_function_commands(calculator, arg));
        return [calculator, impl, scoped_expression](const std::vector<std::pair<std::string, StoredValue>>& assignments) {
            impl->flat_scopes.push_scope();
            for (const auto& [name, value] : assignments) {
                impl->flat_scopes.set(name, value);
            }
            try {
                const StoredValue value = evaluate_expression_value(calculator, impl, scoped_expression, false);
                impl->flat_scopes.pop_scope();
                if (value.is_matrix || value.is_complex || value.is_string) throw std::runtime_error("expected a scalar-valued expression");
                return Calculator::normalize_result(value.exact ? rational_to_double(value.rational) : value.decimal);
            } catch (...) {
                impl->flat_scopes.pop_scope();
                throw;
            }
        };
    };
    
    s.evaluation.build_matrix_evaluator = [calculator, impl](const std::string& arg) {
        const std::string scoped_expression = trim_copy(expand_inline_function_commands(calculator, arg));
        return [calculator, impl, scoped_expression](const std::vector<std::pair<std::string, StoredValue>>& assignments) {
            std::map<std::string, StoredValue> scoped_variables = visible_variables(impl).snapshot();
            for (const auto& [name, value] : assignments) scoped_variables[name] = value;
            const HasScriptFunctionCallback has_script_function = [calculator, impl](const std::string& name) { return has_visible_script_function(impl, name); };
            const InvokeScriptFunctionDecimalCallback invoke_script_function = [calculator, impl](const std::string& name, const std::vector<double>& args) { return invoke_script_function_decimal(calculator, impl, name, args); };
            matrix::Value val;
            if (!try_evaluate_matrix_expression(scoped_expression, VariableResolver(&scoped_variables, nullptr), &impl->functions, &impl->scalar_functions, &impl->matrix_functions, &impl->value_functions, has_script_function, invoke_script_function, &val) || !val.is_matrix)
                throw std::runtime_error("expected a matrix-valued expression");
            return val.matrix;
        };
    };

    // Symbolic Service
    s.symbolic.resolve_symbolic = [calculator, impl](const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) {
        symbolic_commands::SymbolicResolverContext symbolic_resolver_ctx;
        symbolic_resolver_ctx.resolve_custom_function = [calculator, impl](const std::string& name, std::string* v) {
            const auto it = impl->functions.find(name);
            if (it == impl->functions.end()) throw std::runtime_error("unknown custom function: " + name);
            
            // 为兼容符号解析接口，合并参数名（注：符号引擎目前主要支持单变量解析）
            std::string params;
            for (std::size_t i = 0; i < it->second.parameter_names.size(); ++i) {
                params += it->second.parameter_names[i];
                if (i + 1 < it->second.parameter_names.size()) params += ",";
            }
            *v = params;
            return SymbolicExpression::parse(it->second.expression);
        };
        symbolic_resolver_ctx.has_custom_function = [impl](const std::string& name) {
            return impl->functions.find(name) != impl->functions.end();
        };
        symbolic_resolver_ctx.expand_inline = [calculator](const std::string& a) {
            return expand_inline_function_commands(calculator, a);
        };
        symbolic_commands::resolve_symbolic_expression(symbolic_resolver_ctx, arg, req, var, expr);
    };
    s.symbolic.expand_inline = [calculator](const std::string& arg) { return expand_inline_function_commands(calculator, arg); };
    s.symbolic.simplify_symbolic = [](const std::string& text) { return SymbolicExpression::parse(text).simplify().to_string(); };
    s.symbolic.evaluate_symbolic_at = [calculator, impl](const SymbolicExpression& expr, const std::string& var, double p) {
        const auto existing = impl->variables.find(var);
        const bool had_existing = existing != impl->variables.end();
        StoredValue backup;
        if (had_existing) backup = existing->second;
        StoredValue temporary;
        temporary.decimal = p;
        temporary.exact = false;
        impl->variables[var] = temporary;
        auto cleanup = [&]() {
            if (had_existing) impl->variables[var] = backup;
            else impl->variables.erase(var);
        };
        try {
            const double value = calculator->evaluate(expr.to_string());
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
    s.symbolic.parse_symbolic_expr_list = [calculator](const std::string& arg) { return symbolic_commands::parse_symbolic_expression_list(arg, [calculator](const std::string& a) { return expand_inline_function_commands(calculator, a); }); };
    
    // Symbolic build_analysis is a special case that needs svc reference, 
    // it will be set by the caller if needed or we can wrap it here.
    s.symbolic.build_analysis = [calculator, impl, s_capture = s](const std::string& argument) mutable {
        const std::string trimmed_argument = trim_copy(argument);
        const auto direct_function = impl->functions.find(trimmed_argument);
        if (direct_function != impl->functions.end() &&
            direct_function->second.parameter_names.size() == 1) {
            const std::string variable_name = direct_function->second.parameter_names.front();
            const std::string expression = direct_function->second.expression;
            FunctionAnalysis analysis(variable_name);
            analysis.define(expression);
            analysis.set_evaluator(s_capture.evaluation.build_decimal_evaluator(expression));
            return analysis;
        }

        std::string variable_name;
        SymbolicExpression expression;
        symbolic_commands::SymbolicResolverContext symbolic_resolver_ctx;
        symbolic_resolver_ctx.resolve_custom_function = [calculator, impl](const std::string& name, std::string* v) {
            const auto it = impl->functions.find(name);
            if (it == impl->functions.end()) throw std::runtime_error("unknown custom function: " + name);
            
            // 为兼容符号解析接口，合并参数名（注：符号引擎目前主要支持单变量解析）
            std::string params;
            for (std::size_t i = 0; i < it->second.parameter_names.size(); ++i) {
                params += it->second.parameter_names[i];
                if (i + 1 < it->second.parameter_names.size()) params += ",";
            }
            *v = params;
            return SymbolicExpression::parse(it->second.expression);
        };
        symbolic_resolver_ctx.has_custom_function = [impl](const std::string& name) { return impl->functions.find(name) != impl->functions.end(); };
        symbolic_resolver_ctx.expand_inline = [calculator](const std::string& a) { return expand_inline_function_commands(calculator, a); };
        symbolic_commands::resolve_symbolic_expression(symbolic_resolver_ctx, argument, true, &variable_name, &expression);

        FunctionAnalysis analysis(variable_name);
        analysis.define(expression.to_string());
        analysis.set_evaluator(s_capture.evaluation.build_decimal_evaluator(expression.to_string()));
        return analysis;
    };

    // Environment Service
    s.env.has_variable = [impl](const std::string& name) { return impl->variables.find(name) != impl->variables.end(); };
    s.env.has_function = [impl](const std::string& name) { return impl->functions.find(name) != impl->functions.end(); };
    s.env.list_variables = [calculator]() { return calculator->list_variables(); };
    s.env.list_functions = [impl]() {
        if (impl->functions.empty() && impl->script_functions.empty()) return std::string("No custom functions defined.");
        std::ostringstream out;
        bool first = true;
        for (const auto& [name, function] : impl->functions) {
            if (!first) out << '\n';
            first = false;
            out << name << "(";
            for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
                if (i != 0) out << ", ";
                out << function.parameter_names[i];
            }
            out << ") = " << function.expression;
        }
        for (const auto& [name, function] : impl->script_functions) {
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
    s.env.clear_variable = [calculator](const std::string& name) { return calculator->clear_variable(name); };
    s.env.clear_function = [impl](const std::string& name) {
        const auto simple_it = impl->functions.find(name);
        if (simple_it != impl->functions.end()) {
            impl->functions.erase(simple_it);
            return std::string("Cleared custom function: ") + name;
        }
        const auto script_it = impl->script_functions.find(name);
        if (script_it != impl->script_functions.end()) {
            impl->script_functions.erase(script_it);
            return std::string("Cleared custom function: ") + name;
        }
        throw std::runtime_error("unknown custom function: " + name);
    };
    s.env.clear_all_variables = [calculator]() { return calculator->clear_all_variables(); };
    s.env.save_state = [calculator](const std::string& p) { return calculator->save_state(p); };
    s.env.load_state = [calculator](const std::string& p) { return calculator->load_state(p); };
    s.env.export_variable = [calculator](const std::string& p) { return calculator->export_variable(p); };
    s.env.execute_script = [calculator](const std::string& c, bool e) { return calculator->execute_script(c, e); };
    s.env.execute_script_file = [calculator](const std::string& p, bool e) { return calculator->execute_script_file(p, e); };
    s.env.clear_all_functions = [impl]() {
        impl->functions.clear();
        impl->script_functions.clear();
        return std::string("Cleared all custom functions.");
    };
    s.env.set_exact_mode = [](bool m) { return std::string("Exact mode: ") + (m ? "ON" : "OFF"); };
    s.env.set_symbolic_mode = [calculator](bool m) { return calculator->set_symbolic_constants_mode(m); };
    s.env.set_precision = [calculator](int p) { return calculator->set_display_precision(p); };
    s.env.set_hex_prefix = [calculator](bool m) { return calculator->set_hex_prefix_mode(m); };
    s.env.set_hex_uppercase = [calculator](bool m) { return calculator->set_hex_uppercase_mode(m); };

    // Global Helpers
    s.parse_symbolic_vars = [](const std::vector<std::string>& args, std::size_t start, const std::vector<std::string>& fallback) { return symbolic_commands::parse_symbolic_variable_arguments(args, start, fallback); };
    s.is_matrix_argument = [calculator, impl](const std::string& arg) {
        const VariableResolver visible = visible_variables(impl);
        const HasScriptFunctionCallback has_script_function = [calculator, impl](const std::string& name) { return has_visible_script_function(impl, name); };
        const InvokeScriptFunctionDecimalCallback invoke_script_function = [calculator, impl](const std::string& name, const std::vector<double>& args) { return invoke_script_function_decimal(calculator, impl, name, args); };
        matrix::Value value;
        return try_evaluate_matrix_expression(trim_copy(arg), visible, &impl->functions, &impl->scalar_functions, &impl->matrix_functions, &impl->value_functions, has_script_function, invoke_script_function, &value) && value.is_matrix;
    };
    s.parse_matrix_argument = [calculator, impl](const std::string& arg, const std::string& c) {
        const StoredValue value = evaluate_expression_value(calculator, impl, arg, false);
        if (!value.is_matrix) throw std::runtime_error(c + " expects a matrix or vector argument");
        return value.matrix;
    };
    s.render_plot = [calculator, impl](const std::vector<std::string>& args, bool gnuplot) {
        plot::PlotContext ctx;
        ctx.variables = visible_variables(impl);
        ctx.functions = &impl->functions;
        ctx.scalar_functions = &impl->scalar_functions;
        ctx.has_script_function = [impl](const std::string& name) {
            return has_visible_script_function(impl, name);
        };
        ctx.invoke_script_function = [calculator, impl](const std::string& name, const std::vector<double>& call_args) {
            return invoke_script_function_decimal(calculator, impl, name, call_args);
        };
        return gnuplot ? plot::handle_gnuplot_command(ctx, args)
                       : plot::handle_plot_command(ctx, args);
    };
    s.is_integer_double = [](double x, double eps) { return is_integer_double(x, eps); };
    s.round_to_long_long = [](double x) { return round_to_long_long(x); };

    return s;
}

} // namespace core

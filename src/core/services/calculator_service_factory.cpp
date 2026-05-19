// ============================================================================
// 核心服务工厂实现
// ============================================================================
//
// 构建计算器核心提供的所有服务接口，包括：
// - 求值服务（表达式求值、作用域求值器构建）
// - 符号服务（符号表达式解析、简化、求值）
// - 环境服务（变量/函数管理、状态持久化）
// ============================================================================

#include "core/services/calculator_service_factory.h"
#include "core/services/core_managers.h"
#include "core/services/service_registry.h"
#include "core/api/calculator_internal_types.h"
#include "math/helpers/integer_helpers.h"
#include "math/helpers/combinatorics.h"
#include "math/helpers/bitwise_helpers.h"
#include "math/helpers/unit_conversions.h"
#include "math/helpers/base_conversions.h"
#include "execution/resolver/variable_resolver.h"
#include "parser/grammars/unified_expression_parser.h"
#include "execution/engine/script_runtime.h"
#include "symbolic/modules/symbolic_module.h"
#include "symbolic/core/symbolic_expression.h"
#include "analysis/calculus/function_analysis.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "statistics/statistics.h"
#include "statistics/probability.h"
#include "execution/engine/inline_expander.h"
#include "math/mymath.h"
#include <sstream>

namespace core {

/// 构建核心服务对象
CoreServices build_core_services(Calculator* calculator, Calculator::Impl* impl) {
    CoreServices s;


    // Evaluation Service
    s.evaluation.parse_decimal = [impl](const std::string& arg) {
        return parse_decimal_expression(arg, impl->variables_ptr->create_resolver(),
            impl->functions_ptr->get_custom_functions_map(),
            impl->functions_ptr->get_scalar_functions());
    };

    s.evaluation.evaluate_value = [impl](const std::string& arg, bool exact) {
        return evaluate_expression_value(impl, arg, exact);
    };

    s.evaluation.normalize_result = [](Scalar v) { return Calculator::normalize_result(v); };

    s.evaluation.build_decimal_evaluator = [impl](const std::string& arg) {
        const std::string scoped_expression = trim_copy(impl->expand_inline(arg));
        auto base_resolver = impl->variables_ptr->create_resolver();
        return [impl, scoped_expression, base_resolver](const std::vector<std::pair<std::string, Scalar>>& assignments) {
            std::map<std::string, StoredValue> override_vars;
            for (const auto& [name, value] : assignments) {
                StoredValue stored;
                stored.decimal = value;
                stored.exact = false;
                override_vars[name] = stored;
            }
            const HasScriptFunctionCallback has_script_function = [impl](const std::string& name) {
                return impl->functions_ptr->get_script(name) != nullptr;
            };
            const InvokeScriptFunctionDecimalCallback invoke_script_function = [impl](const std::string& name, const std::vector<Scalar>& args) {
                return invoke_script_function_decimal(impl, name, args);
            };
            VariableResolver chained_resolver(nullptr, nullptr, &override_vars, &base_resolver);
            return parse_decimal_expression(scoped_expression, chained_resolver,
                impl->functions_ptr->get_custom_functions_map(),
                impl->functions_ptr->get_scalar_functions(),
                has_script_function, invoke_script_function);
        };
    };

    s.evaluation.build_scalar_evaluator = [impl](const std::string& arg) {
        const std::string scoped_expression = trim_copy(impl->expand_inline(arg));
        return [impl, scoped_expression](const std::vector<std::pair<std::string, StoredValue>>& assignments) {
            impl->variables_ptr->push_scope();
            for (const auto& [name, value] : assignments) {
                impl->variables_ptr->set_local(name, value);
            }
            try {
                const StoredValue value = evaluate_expression_value(impl, scoped_expression, false);
                impl->variables_ptr->pop_scope();
                if (value.is_matrix || value.is_complex || value.is_string)
                    throw std::runtime_error("expected a scalar-valued expression");
                return Calculator::normalize_result(value.exact ? rational_to_double(value.rational) : value.decimal);
            } catch (...) {
                impl->variables_ptr->pop_scope();
                throw;
            }
        };
    };

    s.evaluation.build_matrix_evaluator = [impl](const std::string& arg) {
        const std::string scoped_expression = trim_copy(impl->expand_inline(arg));
        return [impl, scoped_expression](const std::vector<std::pair<std::string, StoredValue>>& assignments) {
            std::map<std::string, StoredValue> scoped_variables = impl->variables_ptr->get_all_globals();
            for (const auto& [name, value] : assignments) {
                scoped_variables[name] = value;
            }
            const HasScriptFunctionCallback has_script_function = [impl](const std::string& name) {
                return impl->functions_ptr->get_script(name) != nullptr;
            };
            const InvokeScriptFunctionDecimalCallback invoke_script_function = [impl](const std::string& name, const std::vector<Scalar>& args) {
                return invoke_script_function_decimal(impl, name, args);
            };
            matrix::Value val;
            if (!try_evaluate_matrix_expression(scoped_expression,
                    VariableResolver(&scoped_variables, nullptr),
                    impl->functions_ptr->get_custom_functions_map(),
                    impl->functions_ptr->get_scalar_functions(),
                    impl->functions_ptr->get_matrix_functions(),
                    impl->functions_ptr->get_value_functions(),
                    has_script_function, invoke_script_function, &val) || !val.is_matrix) {
                throw std::runtime_error("expected a matrix-valued expression");
            }
            return val.matrix;
        };
    };

    // Environment Service
    s.env.has_variable = [impl](const std::string& name) { return impl->variables_ptr->has(name); };
    s.env.has_function = [impl](const std::string& name) { return impl->functions_ptr->has_function(name); };
    s.env.list_variables = [impl]() { return impl->parent->list_variables(); };
    s.env.list_functions = [impl]() {
        const auto custom_names = impl->functions_ptr->get_custom_names();
        const auto script_names = impl->functions_ptr->get_script_names();
        if (custom_names.empty() && script_names.empty()) return std::string("No custom functions defined.");
        std::ostringstream out;
        bool first = true;
        for (const auto& name : custom_names) {
            if (!first) out << '\n';
            first = false;
            const CustomFunction* func = impl->functions_ptr->get_custom(name);
            out << name << "(";
            for (std::size_t i = 0; i < func->parameter_names.size(); ++i) {
                if (i != 0) out << ", ";
                out << func->parameter_names[i];
            }
            out << ") = " << func->expression;
        }
        for (const auto& name : script_names) {
            if (!first) out << '\n';
            first = false;
            const ScriptFunction* func = impl->functions_ptr->get_script(name);
            out << name << "(";
            for (std::size_t i = 0; i < func->parameter_names.size(); ++i) {
                if (i != 0) out << ", ";
                out << func->parameter_names[i];
            }
            out << ") = { ... }";
        }
        return out.str();
    };
    s.env.clear_variable = [impl](const std::string& name) { return impl->parent->clear_variable(name); };
    s.env.clear_function = [impl](const std::string& name) {
        if (impl->functions_ptr->get_custom(name)) {
            impl->functions_ptr->remove_function(name);
            return std::string("Cleared custom function: ") + name;
        }
        if (impl->functions_ptr->get_script(name)) {
            impl->functions_ptr->remove_function(name);
            return std::string("Cleared custom function: ") + name;
        }
        throw std::runtime_error("unknown custom function: " + name);
    };
    s.env.clear_all_variables = [impl]() { return impl->parent->clear_all_variables(); };
    s.env.save_state = [impl](const std::string& p) { return impl->parent->save_state(p); };
    s.env.load_state = [impl](const std::string& p) { return impl->parent->load_state(p); };
    s.env.export_variable = [impl](const std::string& p) { return impl->parent->export_variable(p); };
    s.env.execute_script = [impl](const std::string& c, bool e) { return impl->parent->execute_script(c, e); };
    s.env.execute_script_file = [impl](const std::string& p, bool e) { return impl->parent->execute_script_file(p, e); };
    s.env.clear_all_functions = [impl]() {
        impl->functions_ptr->clear_all();
        return std::string("Cleared all custom functions.");
    };
    s.env.set_exact_mode = [](bool m) { return std::string("Exact mode: ") + (m ? "ON" : "OFF"); };
    s.env.set_symbolic_mode = [impl](bool m) { return impl->parent->set_symbolic_constants_mode(m); };
    s.env.set_precision = [impl](int p) { return impl->parent->set_display_precision(p); };
    s.env.set_hex_prefix = [impl](bool m) { return impl->parent->set_hex_prefix_mode(m); };
    s.env.set_hex_uppercase = [impl](bool m) { return impl->parent->set_hex_uppercase_mode(m); };

    // Help Services
    s.env.help_text = [impl]() { return impl->parent->help_text(); };
    s.env.help_topic = [impl](const std::string& t) { return impl->parent->help_topic(t); };

    // Global Helpers
    s.parse_symbolic_vars = [](const std::vector<std::string>& args, std::size_t start, const std::vector<std::string>& fallback) {
        return symbolic_commands::parse_symbolic_variable_arguments(args, start, fallback);
    };

    s.is_matrix_argument = [impl](const std::string& arg) {
        const VariableResolver visible = impl->variables_ptr->create_resolver();
        const HasScriptFunctionCallback has_script_function = [impl](const std::string& name) {
            return impl->functions_ptr->get_script(name) != nullptr;
        };
        const InvokeScriptFunctionDecimalCallback invoke_script_function = [impl](const std::string& name, const std::vector<Scalar>& args) {
            return invoke_script_function_decimal(impl, name, args);
        };
        matrix::Value value;
        return try_evaluate_matrix_expression(trim_copy(arg), visible,
            impl->functions_ptr->get_custom_functions_map(),
            impl->functions_ptr->get_scalar_functions(),
            impl->functions_ptr->get_matrix_functions(),
            impl->functions_ptr->get_value_functions(),
            has_script_function, invoke_script_function, &value) && value.is_matrix;
    };

    s.parse_matrix_argument = [impl](const std::string& arg, const std::string& c) {
        const StoredValue value = evaluate_expression_value(impl, arg, false);
        if (!value.is_matrix) throw std::runtime_error(c + " expects a matrix or vector argument");
        return value.matrix;
    };

    s.is_integer_double = [](Scalar x, Scalar eps) { return is_integer_double(x, eps); };
    s.round_to_long_long = [](Scalar x) { return round_to_long_long(x); };

    // 动态加载模块服务
    global_service_registry().build_all(s, calculator, impl);

    return s;
}

} // namespace core

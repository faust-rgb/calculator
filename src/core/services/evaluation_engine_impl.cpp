// ============================================================================
// evaluation_engine_impl.cpp - IEvaluationEngine 的具体实现
// ============================================================================

#include "core/services/evaluation_engine_impl.h"
#include "core/api/calculator_impl.h"
#include "parser/grammars/unified_expression_parser.h"
#include "matrix/matrix_internal.h"
#include "symbolic/base/assumptions.h"
#include "execution/resolver/variable_resolver.h"
#include "core/services/string_utils.h"
#include "execution/engine/script_runtime.h"
#include "plot/calculator_plot.h"

bool EvaluationEngineImpl::is_matrix_argument(const std::string& arg) {
    matrix::Value value;
    return try_evaluate_matrix_expression(arg,
        impl_->variables_ptr->create_resolver(),
        impl_->functions_ptr->get_custom_functions_map(),
        impl_->functions_ptr->get_scalar_functions(),
        impl_->functions_ptr->get_matrix_functions(),
        impl_->functions_ptr->get_value_functions(),
        {}, {}, &value) && value.is_matrix;
}

matrix::Matrix EvaluationEngineImpl::parse_matrix_argument(const std::string& arg, const std::string& command) {
    const StoredValue value = evaluate_expression_value(arg, false);
    if (!value.is_matrix) {
        throw std::runtime_error(command + " expects a matrix or vector argument");
    }
    return value.matrix;
}

void EvaluationEngineImpl::resolve_symbolic(const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) {
    symbolic_commands::SymbolicResolverContext ctx;
    ctx.resolve_custom_function = [this](const std::string& name, std::string* v) {
        const CustomFunction* func = impl_->functions_ptr->get_custom(name);
        if (!func) throw std::runtime_error("unknown custom function: " + name);

        std::string params;
        for (std::size_t i = 0; i < func->parameter_names.size(); ++i) {
            params += func->parameter_names[i];
            if (i + 1 < func->parameter_names.size()) params += ",";
        }
        *v = params;
        return SymbolicExpression::parse(func->expression);
    };
    ctx.has_custom_function = [this](const std::string& name) {
        return impl_->functions_ptr->get_custom(name) != nullptr;
    };
    ctx.expand_inline = [this](const std::string& a) {
        return expand_inline_function_commands(impl_, a);
    };
    symbolic_commands::resolve_symbolic_expression(ctx, arg, req, var, expr);
}

std::string EvaluationEngineImpl::render_plot(const std::vector<std::string>& args, bool gnuplot) {
    plot::PlotContext ctx;
    ctx.variables = impl_->variables_ptr->create_resolver();
    ctx.functions = impl_->functions_ptr->get_custom_functions_map();
    ctx.scalar_functions = impl_->functions_ptr->get_scalar_functions();
    ctx.has_script_function = [this](const std::string& name) {
        return has_visible_script_function(impl_, name);
    };
    ctx.invoke_script_function = [this](const std::string& name, const std::vector<Scalar>& arguments) {
        return invoke_script_function_decimal(impl_, name, arguments);
    };

    if (gnuplot) {
        return plot::handle_gnuplot_command(ctx, args);
    }
    return plot::handle_plot_command(ctx, args);
}

Scalar EvaluationEngineImpl::parse_decimal(const std::string& expr) {
    return parse_decimal_expression(expr,
        impl_->variables_ptr->create_resolver(),
        impl_->functions_ptr->get_custom_functions_map(),
        impl_->functions_ptr->get_scalar_functions());
}

std::function<Scalar(const std::vector<std::pair<std::string, Scalar>>&)> EvaluationEngineImpl::build_scoped_evaluator(const std::string& expression) {
    const std::string scoped_expression = trim_copy(expand_inline_function_commands(impl_, expression));
    auto base_resolver = impl_->variables_ptr->create_resolver();
    return [this, scoped_expression, base_resolver](const std::vector<std::pair<std::string, Scalar>>& assignments) {
        std::map<std::string, StoredValue> override_vars;
        for (const auto& [name, value] : assignments) {
            StoredValue stored;
            stored.decimal = value;
            stored.exact = false;
            override_vars[name] = stored;
        }
        VariableResolver chained_resolver(nullptr, nullptr, &override_vars, &base_resolver);
        return parse_decimal_expression(scoped_expression, chained_resolver,
            impl_->functions_ptr->get_custom_functions_map(),
            impl_->functions_ptr->get_scalar_functions());
    };
}

std::function<Scalar(const std::vector<std::pair<std::string, StoredValue>>&)> EvaluationEngineImpl::build_scoped_scalar_evaluator(const std::string& expression) {
    const std::string scoped_expression = trim_copy(expand_inline_function_commands(impl_, expression));
    return [this, scoped_expression](const std::vector<std::pair<std::string, StoredValue>>& assignments) {
        impl_->variables_ptr->push_scope();
        for (const auto& [name, value] : assignments) {
            impl_->variables_ptr->set_local(name, value);
        }
        try {
            const StoredValue value = evaluate_expression_value(scoped_expression, false);
            impl_->variables_ptr->pop_scope();
            if (value.is_matrix || value.is_complex || value.is_string)
                throw std::runtime_error("expected a scalar-valued expression");
            return Calculator::normalize_result(value.exact ? rational_to_double(value.rational) : value.decimal);
        } catch (...) {
            impl_->variables_ptr->pop_scope();
            throw;
        }
    };
}

std::function<matrix::Matrix(const std::vector<std::pair<std::string, StoredValue>>&)> EvaluationEngineImpl::build_scoped_matrix_evaluator(const std::string& expression) {
    const std::string scoped_expression = trim_copy(expand_inline_function_commands(impl_, expression));
    return [this, scoped_expression](const std::vector<std::pair<std::string, StoredValue>>& assignments) {
        std::map<std::string, StoredValue> scoped_variables = impl_->variables_ptr->get_all_globals();
        for (const auto& [name, value] : assignments) {
            scoped_variables[name] = value;
        }
        matrix::Value val;
        if (!try_evaluate_matrix_expression(scoped_expression,
                VariableResolver(&scoped_variables, nullptr),
                impl_->functions_ptr->get_custom_functions_map(),
                impl_->functions_ptr->get_scalar_functions(),
                impl_->functions_ptr->get_matrix_functions(),
                impl_->functions_ptr->get_value_functions(),
                {}, {}, &val) || !val.is_matrix) {
            throw std::runtime_error("expected a matrix-valued expression");
        }
        return val.matrix;
    };
}

StoredValue EvaluationEngineImpl::evaluate_expression_value(const std::string& arg, bool exact) {
    // 直接调用求值函数，避免通过 impl_->evaluate 造成循环调用
    auto ctx = impl_->locator.resolve<IExecutionContext>();
    if (!ctx) throw std::runtime_error("Invalid execution context");
    return ::evaluate_expression_value(ctx.get(), arg, exact);
}

FunctionAnalysis EvaluationEngineImpl::build_analysis(const std::string& expression) {
    const std::string trimmed_argument = trim_copy(expression);
    const CustomFunction* direct_function = impl_->functions_ptr->get_custom(trimmed_argument);
    if (direct_function && direct_function->parameter_names.size() == 1) {
        const std::string variable_name = direct_function->parameter_names.front();
        const std::string expr = direct_function->expression;
        FunctionAnalysis analysis(variable_name);
        analysis.define(expr);
        analysis.set_evaluator(build_scoped_evaluator(expr));
        return analysis;
    }

    std::string variable_name;
    SymbolicExpression expr;
    resolve_symbolic(trimmed_argument, true, &variable_name, &expr);

    FunctionAnalysis analysis(variable_name);
    analysis.define(expr.to_string());
    analysis.set_evaluator(build_scoped_evaluator(expr.to_string()));
    return analysis;
}

std::vector<std::string> EvaluationEngineImpl::parse_symbolic_vars(const std::vector<std::string>& arguments, std::size_t start_index, const std::vector<std::string>& defaults) {
    return symbolic_commands::parse_symbolic_variable_arguments(arguments, start_index, defaults);
}

std::string EvaluationEngineImpl::expand_inline(const std::string& expression) {
    return expand_inline_function_commands(impl_, expression);
}

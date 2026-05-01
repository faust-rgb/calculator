#include "calculator_internal_types.h"
#include "plot/calculator_plot.h"

#include "script_parser.h"

#include <sstream>
#include <stdexcept>

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

    // 显示型功能优先于普通数值/分数显示，例如 hex(255) 应直接得到 "FF"。
    const VariableResolver variables = visible_variables(impl_.get());
    std::string converted;
    if (try_base_conversion_expression(expression,
                                       variables,
                                       &impl_->functions,
                                       {impl_->hex_prefix_mode, impl_->hex_uppercase_mode},
                                       &converted)) {
        return converted;
    }

    if (impl_->symbolic_constants_mode) {
        std::string symbolic_output;
        if (try_symbolic_constant_expression(expression,
                                             variables,
                                             &impl_->functions,
                                             &symbolic_output)) {
            return symbolic_output;
        }
    }

    return format_stored_value(
        evaluate_expression_value(this, impl_.get(), expression, exact_mode),
        impl_->symbolic_constants_mode);
}

std::string Calculator::process_line(const std::string& expression, bool exact_mode) {
    apply_calculator_display_precision(impl_.get());

    std::string_view lhs;
    std::string_view rhs;
    if (!split_assignment(expression, &lhs, &rhs)) {
        return evaluate_for_display(expression, exact_mode);
    }

    if (!is_valid_variable_name(lhs)) {
        throw std::runtime_error("invalid variable name: " + std::string(lhs));
    }
    if (rhs.empty()) {
        throw std::runtime_error("assignment requires a value");
    }

    const StoredValue stored = evaluate_expression_value(this, impl_.get(), std::string(rhs), exact_mode);
    assign_visible_variable(impl_.get(), std::string(lhs), stored);
    return std::string(lhs) + " = " + format_stored_value(stored, impl_->symbolic_constants_mode);
}

std::string Calculator::execute_script(const std::string& source, bool exact_mode) {
    apply_calculator_display_precision(impl_.get());

    script::Program program = script::parse_program(source);
    std::string accumulated_output;
    std::string last_statement_output;
    for (const auto& statement : program.statements) {
        std::string current_output;
        const ScriptSignal signal =
            execute_script_statement(this, impl_.get(), *statement, exact_mode, &current_output, false);
        
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
        }
        
        bool is_last = (statement == program.statements.back());
        if (is_last) {
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
            std::string return_val = signal.has_value ? format_stored_value(signal.value, impl_->symbolic_constants_mode)
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
    return accumulated_output.empty() ? (last_statement_output.empty() ? "OK" : last_statement_output) 
                                      : accumulated_output;
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
    std::string_view inside;
    if (!split_named_call(expression, "factor", &inside)) {
        throw std::runtime_error("expected factor(expression)");
    }

    // 先允许 inside 是一个普通表达式或变量，再检查最终值是否为整数。
    DecimalParser parser(inside, VariableResolver(&impl_->variables, nullptr), &impl_->functions);
    const double value = normalize_result(parser.parse());
    if (!is_integer_double(value)) {
        throw std::runtime_error("factor only accepts integers");
    }

    return factor_integer(round_to_long_long(value));
}

std::string Calculator::plot_expression(const std::string& expression) const {
    std::vector<std::string_view> argument_views;
    if (!split_named_call_with_arguments(expression, "plot", &argument_views)) {
        throw std::runtime_error("expected plot(...)");
    }

    std::vector<std::string> arguments;
    for (auto arg : argument_views) {
        arguments.emplace_back(arg);
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

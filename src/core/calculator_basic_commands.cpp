#include "calculator_internal_types.h"

#include "script_parser.h"

#include <sstream>
#include <stdexcept>

double Calculator::evaluate(const std::string& expression) {
    return normalize_result(evaluate_raw(expression));
}

double Calculator::evaluate_raw(const std::string& expression) {
    const StoredValue value = evaluate_expression_value(this, impl_.get(), expression, false);
    if (value.is_matrix) {
        throw std::runtime_error("matrix expression cannot be used as a scalar");
    }
    return value.decimal;
}

std::string Calculator::evaluate_for_display(const std::string& expression, bool exact_mode) {
    // 显示型功能优先于普通数值/分数显示，例如 hex(255) 应直接得到 "FF"。
    const std::map<std::string, StoredValue> variables = visible_variables(impl_.get());
    std::string converted;
    if (try_base_conversion_expression(expression,
                                       &variables,
                                       &impl_->functions,
                                       {impl_->hex_prefix_mode, impl_->hex_uppercase_mode},
                                       &converted)) {
        return converted;
    }

    if (impl_->symbolic_constants_mode) {
        std::string symbolic_output;
        if (try_symbolic_constant_expression(expression,
                                             &variables,
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
    std::string lhs;
    std::string rhs;
    if (!split_assignment(expression, &lhs, &rhs)) {
        return evaluate_for_display(expression, exact_mode);
    }

    if (!is_valid_variable_name(lhs)) {
        throw std::runtime_error("invalid variable name: " + lhs);
    }
    if (rhs.empty()) {
        throw std::runtime_error("assignment requires a value");
    }

    const StoredValue stored = evaluate_expression_value(this, impl_.get(), rhs, exact_mode);
    assign_visible_variable(impl_.get(), lhs, stored);
    return lhs + " = " + format_stored_value(stored, impl_->symbolic_constants_mode);
}

std::string Calculator::execute_script(const std::string& source, bool exact_mode) {
    script::Program program = script::parse_program(source);
    std::string last_output;
    for (const auto& statement : program.statements) {
        const ScriptSignal signal =
            execute_script_statement(this, impl_.get(), *statement, exact_mode, &last_output, false);
        if (signal.kind == ScriptSignal::Kind::kReturn) {
            return signal.has_value ? format_stored_value(signal.value, impl_->symbolic_constants_mode)
                                    : (last_output.empty() ? "OK" : last_output);
        }
        if (signal.kind == ScriptSignal::Kind::kBreak ||
            signal.kind == ScriptSignal::Kind::kContinue) {
            throw std::runtime_error("break/continue can only be used inside loops");
        }
    }
    return last_output.empty() ? "OK" : last_output;
}

std::string Calculator::list_variables() const {
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
    std::string inside;
    if (!split_named_call(expression, "factor", &inside)) {
        throw std::runtime_error("expected factor(expression)");
    }

    // 先允许 inside 是一个普通表达式或变量，再检查最终值是否为整数。
    DecimalParser parser(inside, &impl_->variables, &impl_->functions);
    const double value = normalize_result(parser.parse());
    if (!is_integer_double(value)) {
        throw std::runtime_error("factor only accepts integers");
    }

    return factor_integer(round_to_long_long(value));
}

std::string Calculator::base_conversion_expression(const std::string& expression) const {
    std::string converted;
    if (!try_base_conversion_expression(expression,
                                        &impl_->variables,
                                        &impl_->functions,
                                        {impl_->hex_prefix_mode, impl_->hex_uppercase_mode},
                                        &converted)) {
        throw std::runtime_error("expected bin(...), oct(...), hex(...), or base(value, base)");
    }
    return converted;
}

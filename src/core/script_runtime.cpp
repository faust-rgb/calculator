#include "calculator_internal_types.h"

#include "mymath.h"
#include "script_parser.h"

#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

std::map<std::string, StoredValue> visible_variables(const Calculator::Impl* impl) {
    std::map<std::string, StoredValue> merged = impl->variables;
    for (const auto& scope : impl->local_scopes) {
        for (const auto& [name, value] : scope) {
            merged[name] = value;
        }
    }
    return merged;
}

bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name) {
    return impl->script_functions.find(name) != impl->script_functions.end();
}

void assign_visible_variable(Calculator::Impl* impl,
                             const std::string& name,
                             const StoredValue& value) {
    for (auto it = impl->local_scopes.rbegin(); it != impl->local_scopes.rend(); ++it) {
        const auto existing = it->find(name);
        if (existing != it->end()) {
            (*it)[name] = value;
            return;
        }
    }

    const auto global_existing = impl->variables.find(name);
    if (global_existing != impl->variables.end()) {
        impl->variables[name] = value;
        return;
    }

    if (!impl->local_scopes.empty()) {
        impl->local_scopes.back()[name] = value;
        return;
    }

    impl->variables[name] = value;
}

ScriptSignal ScriptSignal::make_return(const StoredValue& return_value) {
    ScriptSignal signal;
    signal.kind = Kind::kReturn;
    signal.has_value = true;
    signal.value = return_value;
    return signal;
}

ScriptSignal ScriptSignal::make_break() {
    ScriptSignal signal;
    signal.kind = Kind::kBreak;
    return signal;
}

ScriptSignal ScriptSignal::make_continue() {
    ScriptSignal signal;
    signal.kind = Kind::kContinue;
    return signal;
}

StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode);
ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope);
ScriptSignal execute_script_block(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope);

script::StatementPtr clone_statement(const script::Statement& statement);
std::string render_script_statement(const script::Statement& statement, int indent);
std::string render_script_block(const script::BlockStatement& block, int indent);

std::string indent_text(int indent) {
    return std::string(static_cast<std::size_t>(indent) * 2, ' ');
}

std::unique_ptr<script::BlockStatement> clone_block_statement(const script::BlockStatement& block) {
    auto clone = std::make_unique<script::BlockStatement>();
    for (const auto& statement : block.statements) {
        clone->statements.push_back(clone_statement(*statement));
    }
    return clone;
}

std::string render_script_block(const script::BlockStatement& block, int indent) {
    std::ostringstream out;
    out << "{\n";
    for (const auto& statement : block.statements) {
        out << render_script_statement(*statement, indent + 1);
    }
    out << indent_text(indent) << "}";
    return out.str();
}

std::string render_script_statement(const script::Statement& statement, int indent) {
    const std::string prefix = indent_text(indent);
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return prefix + render_script_block(
                                static_cast<const script::BlockStatement&>(statement),
                                indent) + "\n";
        case script::Statement::Kind::kSimple:
            return prefix + static_cast<const script::SimpleStatement&>(statement).text + ";\n";
        case script::Statement::Kind::kIf: {
            const auto& source = static_cast<const script::IfStatement&>(statement);
            std::string rendered =
                prefix + "if (" + source.condition + ") " +
                render_script_statement(*source.then_branch, indent).substr(prefix.size());
            if (source.else_branch) {
                rendered.pop_back();
                rendered += prefix + "else " +
                            render_script_statement(*source.else_branch, indent).substr(prefix.size());
            }
            return rendered;
        }
        case script::Statement::Kind::kWhile: {
            const auto& source = static_cast<const script::WhileStatement&>(statement);
            return prefix + "while (" + source.condition + ") " +
                   render_script_statement(*source.body, indent).substr(prefix.size());
        }
        case script::Statement::Kind::kFor: {
            const auto& source = static_cast<const script::ForStatement&>(statement);
            return prefix + "for (" + source.initializer + "; " + source.condition + "; " +
                   source.step + ") " +
                   render_script_statement(*source.body, indent).substr(prefix.size());
        }
        case script::Statement::Kind::kFunction: {
            const auto& source = static_cast<const script::FunctionStatement&>(statement);
            std::ostringstream out;
            out << prefix << "fn " << source.name << "(";
            for (std::size_t i = 0; i < source.parameters.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << source.parameters[i];
            }
            out << ") " << render_script_block(*source.body, indent) << '\n';
            return out.str();
        }
        case script::Statement::Kind::kReturn: {
            const auto& source = static_cast<const script::ReturnStatement&>(statement);
            return prefix + "return" +
                   (source.has_expression ? " " + source.expression : "") +
                   ";\n";
        }
        case script::Statement::Kind::kBreak:
            return prefix + "break;\n";
        case script::Statement::Kind::kContinue:
            return prefix + "continue;\n";
    }

    throw std::runtime_error("unknown script statement kind");
}

script::StatementPtr clone_statement(const script::Statement& statement) {
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return clone_block_statement(static_cast<const script::BlockStatement&>(statement));
        case script::Statement::Kind::kSimple: {
            auto clone = std::make_unique<script::SimpleStatement>();
            clone->text = static_cast<const script::SimpleStatement&>(statement).text;
            return clone;
        }
        case script::Statement::Kind::kIf: {
            const auto& source = static_cast<const script::IfStatement&>(statement);
            auto clone = std::make_unique<script::IfStatement>();
            clone->condition = source.condition;
            clone->then_branch = clone_statement(*source.then_branch);
            if (source.else_branch) {
                clone->else_branch = clone_statement(*source.else_branch);
            }
            return clone;
        }
        case script::Statement::Kind::kWhile: {
            const auto& source = static_cast<const script::WhileStatement&>(statement);
            auto clone = std::make_unique<script::WhileStatement>();
            clone->condition = source.condition;
            clone->body = clone_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kFor: {
            const auto& source = static_cast<const script::ForStatement&>(statement);
            auto clone = std::make_unique<script::ForStatement>();
            clone->initializer = source.initializer;
            clone->condition = source.condition;
            clone->step = source.step;
            clone->body = clone_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kFunction: {
            const auto& source = static_cast<const script::FunctionStatement&>(statement);
            auto clone = std::make_unique<script::FunctionStatement>();
            clone->name = source.name;
            clone->parameters = source.parameters;
            clone->body = clone_block_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kReturn: {
            const auto& source = static_cast<const script::ReturnStatement&>(statement);
            auto clone = std::make_unique<script::ReturnStatement>();
            clone->has_expression = source.has_expression;
            clone->expression = source.expression;
            return clone;
        }
        case script::Statement::Kind::kBreak:
            return std::make_unique<script::BreakStatement>();
        case script::Statement::Kind::kContinue:
            return std::make_unique<script::ContinueStatement>();
    }

    throw std::runtime_error("unknown script statement kind");
}

bool truthy_value(const StoredValue& value) {
    if (value.is_matrix) {
        throw std::runtime_error("matrix values cannot be used as script conditions");
    }
    return !mymath::is_near_zero(value.exact
                                     ? rational_to_double(value.rational)
                                     : value.decimal,
                                 1e-10);
}

double invoke_script_function_decimal(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& name,
                                      const std::vector<double>& arguments) {
    auto it = impl->script_functions.find(name);
    if (it == impl->script_functions.end()) {
        throw std::runtime_error("unknown function: " + name);
    }

    const ScriptFunction& function = it->second;
    if (arguments.size() != function.parameter_names.size()) {
        throw std::runtime_error("script function " + name + " expects " +
                                 std::to_string(function.parameter_names.size()) +
                                 " arguments");
    }

    std::map<std::string, StoredValue> frame;
    for (std::size_t i = 0; i < arguments.size(); ++i) {
        StoredValue parameter_value;
        parameter_value.decimal = arguments[i];
        parameter_value.exact = false;
        frame[function.parameter_names[i]] = parameter_value;
    }

    impl->local_scopes.push_back(frame);
    std::string ignored_output;
    try {
        const ScriptSignal signal =
            execute_script_block(calculator, impl, *function.body, false, &ignored_output, false);
        impl->local_scopes.pop_back();

        if (signal.kind != ScriptSignal::Kind::kReturn || !signal.has_value) {
            throw std::runtime_error("script function " + name + " must return a value");
        }
        if (signal.value.is_matrix) {
            throw std::runtime_error("script function " + name +
                                     " cannot be used as a scalar expression");
        }
        if (signal.value.is_string) {
            throw std::runtime_error("script function " + name +
                                     " cannot be used as a numeric expression");
        }
        return signal.value.exact
                   ? rational_to_double(signal.value.rational)
                   : signal.value.decimal;
    } catch (...) {
        impl->local_scopes.pop_back();
        throw;
    }
}

StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode) {
    const std::string expanded_expression =
        expand_inline_function_commands(calculator, expression);
    const std::string trimmed = trim_copy(expanded_expression);
    const std::map<std::string, StoredValue> variables = visible_variables(impl);
    if (is_string_literal(trimmed)) {
        StoredValue stored;
        stored.is_string = true;
        stored.string_value = parse_string_literal_value(trimmed);
        return stored;
    }
    if (is_identifier_text(trimmed)) {
        const auto it = variables.find(trimmed);
        if (it != variables.end() && it->second.is_string) {
            return it->second;
        }
    }

    const HasScriptFunctionCallback has_script_function =
        [impl](const std::string& name) {
            return has_visible_script_function(impl, name);
        };
    const InvokeScriptFunctionDecimalCallback invoke_script_function =
        [calculator, impl](const std::string& name, const std::vector<double>& arguments) {
            return invoke_script_function_decimal(calculator, impl, name, arguments);
        };

    StoredValue stored;
    std::vector<std::string> rational_arguments;
    if (split_named_call_with_arguments(trimmed, "rat", &rational_arguments)) {
        if (rational_arguments.size() != 1 && rational_arguments.size() != 2) {
            throw std::runtime_error(
                "rat expects one argument or expression plus max_denominator");
        }

        const StoredValue value =
            evaluate_expression_value(calculator, impl, rational_arguments[0], false);
        if (value.is_matrix) {
            throw std::runtime_error("rat cannot approximate a matrix value");
        }
        if (value.is_string) {
            throw std::runtime_error("rat cannot approximate a string value");
        }

        long long max_denominator = 999;
        if (rational_arguments.size() == 2) {
            const StoredValue max_denominator_value =
                evaluate_expression_value(calculator, impl, rational_arguments[1], false);
            if (max_denominator_value.is_matrix || max_denominator_value.is_string) {
                throw std::runtime_error("rat max_denominator must be a positive integer");
            }
            const double scalar =
                max_denominator_value.exact
                    ? rational_to_double(max_denominator_value.rational)
                    : max_denominator_value.decimal;
            if (!is_integer_double(scalar) || scalar <= 0.0) {
                throw std::runtime_error("rat max_denominator must be a positive integer");
            }
            max_denominator = round_to_long_long(scalar);
        }

        if (value.exact && value.rational.denominator <= max_denominator) {
            return value;
        }

        const double decimal_value = value.exact
                                         ? rational_to_double(value.rational)
                                         : value.decimal;
        long long numerator = 0;
        long long denominator = 1;
        if (!mymath::best_rational_approximation(decimal_value,
                                                 &numerator,
                                                 &denominator,
                                                 max_denominator)) {
            throw std::runtime_error("rat could not compute a rational approximation");
        }

        stored.exact = true;
        stored.rational = Rational(numerator, denominator);
        stored.decimal = rational_to_double(stored.rational);
        return stored;
    }

    if (exact_mode) {
        try {
            stored.rational = parse_exact_expression(trimmed,
                                                     &variables,
                                                     &impl->functions,
                                                     has_script_function);
            stored.exact = true;
            stored.decimal = rational_to_double(stored.rational);
            return stored;
        } catch (const ExactModeUnsupported&) {
        }
    }

    matrix::Value matrix_value;
    if (try_evaluate_matrix_expression(trimmed,
                                       &variables,
                                       &impl->functions,
                                       has_script_function,
                                       invoke_script_function,
                                       &matrix_value)) {
        if (matrix_value.is_matrix) {
            stored.is_matrix = true;
            stored.matrix = matrix_value.matrix;
        } else {
            stored.decimal = matrix_value.scalar;
            stored.exact = false;
        }
        return stored;
    }

    if (!exact_mode) {
        try {
            const PreciseDecimal precise_value =
                parse_precise_decimal_expression(trimmed, &variables);
            stored.decimal = precise_value.to_double();
            stored.exact = false;
            stored.has_precise_decimal_text = true;
            stored.precise_decimal_text = precise_value.to_string();
            return stored;
        } catch (const PreciseDecimalUnsupported&) {
        }
    }

    DecimalParser parser(trimmed,
                         &variables,
                         &impl->functions,
                         has_script_function,
                         invoke_script_function);
    const double parsed_value = parser.parse();
    stored.decimal = parsed_value;
    stored.exact = false;
    if (impl->symbolic_constants_mode) {
        std::string symbolic_output;
        if (try_symbolic_constant_expression(trimmed,
                                             &variables,
                                             &impl->functions,
                                             &symbolic_output)) {
            stored.has_symbolic_text = true;
            stored.symbolic_text = symbolic_output;
        }
    }
    return stored;
}

std::string execute_simple_script_line(Calculator* calculator,
                                       Calculator::Impl* impl,
                                       const std::string& text,
                                       bool exact_mode) {
    std::vector<std::string> print_arguments;
    if (split_named_call_with_arguments(text, "print", &print_arguments)) {
        if (print_arguments.empty()) {
            throw std::runtime_error("print expects at least one argument");
        }
        std::ostringstream out;
        for (std::size_t i = 0; i < print_arguments.size(); ++i) {
            if (i != 0) {
                out << ' ';
            }
            out << format_print_value(
                evaluate_expression_value(calculator, impl, print_arguments[i], exact_mode),
                impl->symbolic_constants_mode);
        }
        return out.str();
    }

    std::string function_output;
    if (calculator->try_process_function_command(text, &function_output)) {
        return function_output;
    }

    std::string lhs;
    std::string rhs;
    if (split_assignment(text, &lhs, &rhs)) {
        if (!is_valid_variable_name(lhs)) {
            throw std::runtime_error("invalid variable name: " + lhs);
        }
        if (rhs.empty()) {
            throw std::runtime_error("assignment requires a value");
        }
        const StoredValue stored = evaluate_expression_value(calculator, impl, rhs, exact_mode);
        assign_visible_variable(impl, lhs, stored);
        return lhs + " = " + format_stored_value(stored, impl->symbolic_constants_mode);
    }

    return format_stored_value(evaluate_expression_value(calculator, impl, text, exact_mode),
                               impl->symbolic_constants_mode);
}

ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope) {
    switch (statement.kind) {
        case script::Statement::Kind::kBlock:
            return execute_script_block(calculator,
                                        impl,
                                        static_cast<const script::BlockStatement&>(statement),
                                        exact_mode,
                                        last_output,
                                        create_scope);
        case script::Statement::Kind::kSimple: {
            const auto& simple = static_cast<const script::SimpleStatement&>(statement);
            *last_output = execute_simple_script_line(calculator, impl, simple.text, exact_mode);
            return {};
        }
        case script::Statement::Kind::kIf: {
            const auto& if_statement = static_cast<const script::IfStatement&>(statement);
            if (truthy_value(evaluate_expression_value(calculator, impl, if_statement.condition, false))) {
                return execute_script_statement(calculator,
                                                impl,
                                                *if_statement.then_branch,
                                                exact_mode,
                                                last_output,
                                                true);
            }
            if (if_statement.else_branch) {
                return execute_script_statement(calculator,
                                                impl,
                                                *if_statement.else_branch,
                                                exact_mode,
                                                last_output,
                                                true);
            }
            return {};
        }
        case script::Statement::Kind::kWhile: {
            const auto& while_statement = static_cast<const script::WhileStatement&>(statement);
            while (truthy_value(evaluate_expression_value(calculator,
                                                          impl,
                                                          while_statement.condition,
                                                          false))) {
                const ScriptSignal signal =
                    execute_script_statement(calculator,
                                             impl,
                                             *while_statement.body,
                                             exact_mode,
                                             last_output,
                                             true);
                if (signal.kind == ScriptSignal::Kind::kReturn) {
                    return signal;
                }
                if (signal.kind == ScriptSignal::Kind::kBreak) {
                    break;
                }
                if (signal.kind == ScriptSignal::Kind::kContinue) {
                    continue;
                }
            }
            return {};
        }
        case script::Statement::Kind::kFor: {
            const auto& for_statement = static_cast<const script::ForStatement&>(statement);
            impl->local_scopes.push_back({});
            try {
                if (!for_statement.initializer.empty()) {
                    (void)execute_simple_script_line(calculator,
                                                     impl,
                                                     for_statement.initializer,
                                                     exact_mode);
                }
                while (for_statement.condition.empty() ||
                       truthy_value(evaluate_expression_value(calculator,
                                                              impl,
                                                              for_statement.condition,
                                                              false))) {
                    const ScriptSignal signal =
                        execute_script_statement(calculator,
                                                 impl,
                                                 *for_statement.body,
                                                 exact_mode,
                                                 last_output,
                                                 true);
                    if (signal.kind == ScriptSignal::Kind::kReturn) {
                        impl->local_scopes.pop_back();
                        return signal;
                    }
                    if (signal.kind == ScriptSignal::Kind::kBreak) {
                        break;
                    }
                    if (!for_statement.step.empty()) {
                        (void)execute_simple_script_line(calculator,
                                                         impl,
                                                         for_statement.step,
                                                         exact_mode);
                    }
                    if (signal.kind == ScriptSignal::Kind::kContinue) {
                        continue;
                    }
                }
                impl->local_scopes.pop_back();
                return {};
            } catch (...) {
                impl->local_scopes.pop_back();
                throw;
            }
        }
        case script::Statement::Kind::kFunction: {
            const auto& function_statement =
                static_cast<const script::FunctionStatement&>(statement);
            if (!is_valid_variable_name(function_statement.name)) {
                throw std::runtime_error("invalid function name: " + function_statement.name);
            }
            if (is_reserved_function_name(function_statement.name)) {
                throw std::runtime_error("function name is reserved: " +
                                         function_statement.name);
            }
            for (const std::string& parameter_name : function_statement.parameters) {
                if (!is_valid_variable_name(parameter_name)) {
                    throw std::runtime_error("invalid parameter name: " + parameter_name);
                }
            }
            ScriptFunction function;
            function.parameter_names = function_statement.parameters;
            function.body = std::shared_ptr<const script::BlockStatement>(
                clone_block_statement(*function_statement.body).release());
            impl->script_functions[function_statement.name] = function;
            *last_output = function_statement.name + "(...)";
            return {};
        }
        case script::Statement::Kind::kReturn: {
            const auto& return_statement =
                static_cast<const script::ReturnStatement&>(statement);
            if (!return_statement.has_expression) {
                ScriptSignal signal;
                signal.kind = ScriptSignal::Kind::kReturn;
                return signal;
            }
            return ScriptSignal::make_return(
                evaluate_expression_value(calculator, impl, return_statement.expression, exact_mode));
        }
        case script::Statement::Kind::kBreak:
            return ScriptSignal::make_break();
        case script::Statement::Kind::kContinue:
            return ScriptSignal::make_continue();
    }

    throw std::runtime_error("unknown script statement kind");
}

ScriptSignal execute_script_block(Calculator* calculator,
                                  Calculator::Impl* impl,
                                  const script::BlockStatement& block,
                                  bool exact_mode,
                                  std::string* last_output,
                                  bool create_scope) {
    if (create_scope) {
        impl->local_scopes.push_back({});
    }

    try {
        for (const auto& statement : block.statements) {
            const ScriptSignal signal =
                execute_script_statement(calculator,
                                         impl,
                                         *statement,
                                         exact_mode,
                                         last_output,
                                         true);
            if (signal.kind != ScriptSignal::Kind::kNone) {
                if (create_scope) {
                    impl->local_scopes.pop_back();
                }
                return signal;
            }
        }
        if (create_scope) {
            impl->local_scopes.pop_back();
        }
        return {};
    } catch (...) {
        if (create_scope) {
            impl->local_scopes.pop_back();
        }
        throw;
    }
}

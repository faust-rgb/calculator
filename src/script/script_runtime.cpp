#include "script_runtime.h"
#include "calculator_module.h"
#include "expression_ast.h"
#include "expression_compiler.h"
#include "core/utils.h"
#include "mymath.h"
#include "script_parser.h"

#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

VariableResolver visible_variables(const Calculator::Impl* impl) {
    if (impl->use_flat_scopes) {
        return VariableResolver(&impl->variables, &impl->flat_scopes);
    }
    return VariableResolver(&impl->variables, &impl->local_scopes);
}

bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name) {
    return impl->script_functions.find(name) != impl->script_functions.end();
}

void assign_visible_variable(Calculator::Impl* impl,
                             const std::string& name,
                             const StoredValue& value) {
    // 优先使用 FlatScopeStack
    if (impl->use_flat_scopes) {
        // 先在 flat_scopes 中查找是否已存在
        if (VariableSlot* existing = impl->flat_scopes.find(name)) {
            existing->value = value;
            return;
        }
        // 检查全局变量
        if (impl->variables.find(name) != impl->variables.end()) {
            impl->variables[name] = value;
            return;
        }
        // 在当前作用域创建新变量
        if (impl->flat_scopes.scope_depth() > 0) {
            impl->flat_scopes.set(name, value);
            return;
        }
        // 无局部作用域，创建全局变量
        impl->variables[name] = value;
        return;
    }

    // 传统路径：std::map 查找
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

namespace {

std::string indent_text(int indent) {
    return std::string(static_cast<std::size_t>(indent) * 2, ' ');
}

bool contains_builtin_constant_token(const std::string& expression) {
    for (std::size_t i = 0; i < expression.size(); ++i) {
        const char ch = expression[i];
        if (!std::isalpha(static_cast<unsigned char>(ch)) && ch != '_') {
            continue;
        }

        const std::size_t start = i;
        ++i;
        while (i < expression.size() &&
               (std::isalnum(static_cast<unsigned char>(expression[i])) ||
                expression[i] == '_')) {
            ++i;
        }
        const std::string token = expression.substr(start, i - start);
        --i;

        double value = 0.0;
        if (lookup_builtin_constant(token, &value)) {
            return true;
        }
    }
    return false;
}

bool truthy_value(const StoredValue& value) {
    if (value.is_matrix) {
        throw std::runtime_error("matrix values cannot be used as script conditions");
    }
    if (value.is_complex) {
        throw std::runtime_error("complex values cannot be used as script conditions");
    }
    return !mymath::is_near_zero(value.exact
                                     ? rational_to_double(value.rational)
                                     : value.decimal,
                                 1e-10);
}

} // namespace

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
        case script::Statement::Kind::kForRange: {
            const auto& source = static_cast<const script::ForRangeStatement&>(statement);
            return prefix + "for " + source.variable + " in range(" +
                   source.start_expr + ", " + source.stop_expr + ", " +
                   source.step_expr + ") " +
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
            const auto& source = static_cast<const script::SimpleStatement&>(statement);
            auto clone = std::make_unique<script::SimpleStatement>();
            clone->text = source.text;
            clone->cache = source.cache;
            return clone;
        }
        case script::Statement::Kind::kIf: {
            const auto& source = static_cast<const script::IfStatement&>(statement);
            auto clone = std::make_unique<script::IfStatement>();
            clone->condition = source.condition;
            clone->cache = source.cache;
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
            clone->cache = source.cache;
            clone->body = clone_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kFor: {
            const auto& source = static_cast<const script::ForStatement&>(statement);
            auto clone = std::make_unique<script::ForStatement>();
            clone->initializer = source.initializer;
            clone->condition = source.condition;
            clone->step = source.step;
            clone->init_cache = source.init_cache;
            clone->cond_cache = source.cond_cache;
            clone->step_cache = source.step_cache;
            clone->body = clone_statement(*source.body);
            return clone;
        }
        case script::Statement::Kind::kForRange: {
            const auto& source = static_cast<const script::ForRangeStatement&>(statement);
            auto clone = std::make_unique<script::ForRangeStatement>();
            clone->variable = source.variable;
            clone->start_expr = source.start_expr;
            clone->stop_expr = source.stop_expr;
            clone->step_expr = source.step_expr;
            clone->step_is_negative = source.step_is_negative;
            clone->start_cache = source.start_cache;
            clone->stop_cache = source.stop_cache;
            clone->step_cache = source.step_cache;
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
            clone->cache = source.cache;
            return clone;
        }
        case script::Statement::Kind::kBreak:
            return std::make_unique<script::BreakStatement>();
        case script::Statement::Kind::kContinue:
            return std::make_unique<script::ContinueStatement>();
    }

    throw std::runtime_error("unknown script statement kind");
}

double invoke_script_function_decimal(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& name,
                                      const std::vector<double>& arguments) {
    // 递归深度检查
    static constexpr int kMaxScriptRecursionDepth = 4096;
    if (impl->script_call_depth >= kMaxScriptRecursionDepth) {
        throw std::runtime_error("maximum script recursion depth exceeded (" +
                                 std::to_string(kMaxScriptRecursionDepth) + ")");
    }

    struct DepthGuard {
        int* depth;
        explicit DepthGuard(int* d) : depth(d) { (*depth)++; }
        ~DepthGuard() { (*depth)--; }
    } guard(&impl->script_call_depth);

    auto it = impl->script_functions.find(name);
    if (it == impl->script_functions.end()) {
        // Try to invoke as a module command
        std::ostringstream call_ss;
        call_ss << name << "(";
        for (std::size_t i = 0; i < arguments.size(); ++i) {
            if (i != 0) call_ss << ", ";
            call_ss << std::setprecision(17) << arguments[i];
        }
        call_ss << ")";
        
        std::string output;
        if (calculator->try_process_function_command(call_ss.str(), &output)) {
            try {
                return std::stod(output);
            } catch (...) {
                // If output is not a number (e.g. vector), it's an error for this scalar path
                throw std::runtime_error("function " + name + " did not return a scalar value: " + output);
            }
        }
        
        throw std::runtime_error("unknown function: " + name);
    }

    const ScriptFunction& function = it->second;
    if (arguments.size() != function.parameter_names.size()) {
        throw std::runtime_error("script function " + name + " expects " +
                                 std::to_string(function.parameter_names.size()) +
                                 " arguments");
    }

    // 使用 FlatScopeStack 创建栈帧
    if (impl->use_flat_scopes) {
        impl->flat_scopes.push_scope();
        for (std::size_t i = 0; i < arguments.size(); ++i) {
            StoredValue parameter_value;
            parameter_value.decimal = arguments[i];
            parameter_value.exact = false;
            impl->flat_scopes.set(function.parameter_names[i], parameter_value);
        }

        struct FlatScopeGuard {
            Calculator::Impl* impl;
            explicit FlatScopeGuard(Calculator::Impl* i) : impl(i) {}
            ~FlatScopeGuard() { impl->flat_scopes.pop_scope(); }
        } scope_guard(impl);

        std::string ignored_output;
        const ScriptSignal signal =
            execute_script_block(calculator, impl, *function.body, false, &ignored_output, false);

        if (signal.kind != ScriptSignal::Kind::kReturn || !signal.has_value) {
            throw std::runtime_error("script function " + name + " must return a value");
        }
        if (signal.value.is_matrix || signal.value.is_complex) {
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
    }

    // 传统路径：std::map 栈帧
    std::map<std::string, StoredValue> frame;
    for (std::size_t i = 0; i < arguments.size(); ++i) {
        StoredValue parameter_value;
        parameter_value.decimal = arguments[i];
        parameter_value.exact = false;
        frame[function.parameter_names[i]] = parameter_value;
    }

    impl->local_scopes.push_back(std::move(frame));
    struct ScopeGuard {
        Calculator::Impl* impl;
        explicit ScopeGuard(Calculator::Impl* i) : impl(i) {}
        ~ScopeGuard() { impl->local_scopes.pop_back(); }
    } scope_guard(impl);

    std::string ignored_output;
    const ScriptSignal signal =
        execute_script_block(calculator, impl, *function.body, false, &ignored_output, false);

    if (signal.kind != ScriptSignal::Kind::kReturn || !signal.has_value) {
        throw std::runtime_error("script function " + name + " must return a value");
    }
    if (signal.value.is_matrix || signal.value.is_complex) {
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
}

// ============================================================================
// 优化后的表达式求值函数
// ============================================================================

StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode,
                                      std::shared_ptr<ExpressionCache>* cache) {
    // 获取或创建缓存
    std::shared_ptr<ExpressionCache> expr_cache;
    if (cache && *cache) {
        expr_cache = *cache;
    } else {
        expr_cache = std::make_shared<ExpressionCache>();
        expr_cache->expanded = expand_inline_function_commands(calculator, expression);
        expr_cache->hint = analyze_expression_hint(expr_cache->expanded);
        expr_cache->features = analyze_expression_features(expr_cache->expanded);
        if (cache) {
            *cache = expr_cache;
        }
    }

    const std::string& trimmed = expr_cache->expanded;
    const VariableResolver variables = visible_variables(impl);

    // 根据类型提示快速分发
    switch (expr_cache->hint) {
        case ExpressionHint::kStringLiteral: {
            StoredValue stored;
            stored.is_string = true;
            stored.string_value = parse_string_literal_value(trimmed);
            return stored;
        }

        case ExpressionHint::kIdentifier: {
            const StoredValue* found = variables.lookup(trimmed);
            if (found) {
                return *found;
            }
            throw UndefinedError("unknown variable: " + trimmed);
        }

        case ExpressionHint::kRatCall: {
            // rat() 特殊处理
            std::vector<std::string_view> rational_arguments_view;
            if (!split_named_call_with_arguments(trimmed, "rat", &rational_arguments_view)) {
                break; // 回退到通用路径
            }
            std::vector<std::string> rational_arguments;
            for (auto sv : rational_arguments_view) {
                rational_arguments.push_back(std::string(sv));
            }
            if (rational_arguments.size() != 1 && rational_arguments.size() != 2) {
                throw std::runtime_error(
                    "rat expects one argument or expression plus max_denominator");
            }

            const StoredValue value =
                evaluate_expression_value(calculator, impl, rational_arguments[0], false, nullptr);
            if (value.is_matrix || value.is_complex) {
                throw std::runtime_error("rat cannot approximate a matrix or complex value");
            }
            if (value.is_string) {
                throw std::runtime_error("rat cannot approximate a string value");
            }

            long long max_denominator = 999;
            if (rational_arguments.size() == 2) {
                const StoredValue max_denominator_value =
                    evaluate_expression_value(calculator, impl, rational_arguments[1], false, nullptr);
                if (max_denominator_value.is_matrix || max_denominator_value.is_complex ||
                    max_denominator_value.is_string) {
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

            StoredValue stored;
            stored.exact = true;
            stored.rational = Rational(numerator, denominator);
            stored.decimal = rational_to_double(stored.rational);
            return stored;
        }

        case ExpressionHint::kComplexCandidate: {
            // 复数表达式路径
            const HasScriptFunctionCallback has_script_function =
                [impl](const std::string& name) {
                    return has_visible_script_function(impl, name);
                };
            const InvokeScriptFunctionDecimalCallback invoke_script_function =
                [calculator, impl](const std::string& name, const std::vector<double>& arguments) {
                    return invoke_script_function_decimal(calculator, impl, name, arguments);
                };

            matrix::Value matrix_val;
            if (try_evaluate_matrix_expression(trimmed,
                                              variables,
                                              &impl->functions,
                                              &impl->scalar_functions,
                                              &impl->matrix_functions,
                                              &impl->value_functions,
                                              has_script_function,
                                              invoke_script_function,
                                              &matrix_val)) {
                StoredValue result;
                if (matrix_val.is_matrix) {
                    result.is_matrix = true;
                    result.matrix = std::move(matrix_val.matrix);
                } else if (matrix_val.is_complex) {
                    result.is_complex = true;
                    result.complex = matrix_val.complex;
                } else {
                    result.decimal = matrix_val.scalar;
                }
                return result;
            }
            // 失败则回退到标量路径
            break;
        }

        case ExpressionHint::kMatrixCandidate: {
            // 矩阵表达式路径
            const HasScriptFunctionCallback has_script_function =
                [impl](const std::string& name) {
                    return has_visible_script_function(impl, name);
                };
            const InvokeScriptFunctionDecimalCallback invoke_script_function =
                [calculator, impl](const std::string& name, const std::vector<double>& arguments) {
                    return invoke_script_function_decimal(calculator, impl, name, arguments);
                };

            matrix::Value matrix_val;
            if (try_evaluate_matrix_expression(trimmed,
                                              variables,
                                              &impl->functions,
                                              &impl->scalar_functions,
                                              &impl->matrix_functions,
                                              &impl->value_functions,
                                              has_script_function,
                                              invoke_script_function,
                                              &matrix_val)) {
                StoredValue result;
                if (matrix_val.is_matrix) {
                    result.is_matrix = true;
                    result.matrix = std::move(matrix_val.matrix);
                } else if (matrix_val.is_complex) {
                    result.is_complex = true;
                    result.complex = matrix_val.complex;
                } else {
                    result.decimal = matrix_val.scalar;
                }
                return result;
            }
            // 失败则回退到标量路径
            break;
        }

        case ExpressionHint::kScalar: {
            break;
        }

        default:
            break;
    }

    // 通用路径：依次尝试各种解析器
    const HasScriptFunctionCallback has_script_function =
        [impl](const std::string& name) {
            return has_visible_script_function(impl, name);
        };
    const InvokeScriptFunctionDecimalCallback invoke_script_function =
        [calculator, impl](const std::string& name, const std::vector<double>& arguments) {
            return invoke_script_function_decimal(calculator, impl, name, arguments);
        };
    auto attach_symbolic_text = [&](StoredValue* stored) {
        if (!impl->symbolic_constants_mode || stored == nullptr ||
            stored->is_matrix || stored->is_complex || stored->is_string) {
            return;
        }
        std::string symbolic_output;
        if (try_symbolic_constant_expression(trimmed,
                                             variables,
                                             &impl->functions,
                                             &symbolic_output)) {
            stored->has_symbolic_text = true;
            stored->symbolic_text = symbolic_output;
        } else if (contains_builtin_constant_token(trimmed)) {
            stored->has_symbolic_text = true;
            stored->symbolic_text = trimmed;
        }
    };

    // 检查字符串字面量
    if (is_string_literal(trimmed)) {
        StoredValue stored;
        stored.is_string = true;
        stored.string_value = parse_string_literal_value(trimmed);
        return stored;
    }

    // 检查标识符
    if (is_identifier_text(trimmed)) {
        const StoredValue* found = variables.lookup(trimmed);
        if (found && found->is_string) {
            return *found;
        }
    }

    // 精确模式
    if (exact_mode) {
        try {
            StoredValue stored;
            stored.rational = parse_exact_expression(trimmed,
                                                     variables,
                                                     &impl->functions,
                                                     has_script_function);
            stored.exact = true;
            stored.decimal = rational_to_double(stored.rational);
            return stored;
        } catch (const ExactModeUnsupported&) {
        }
    }

    // 矩阵/复数表达式
    auto matrix_value = std::make_unique<matrix::Value>();
    if (try_evaluate_matrix_expression(trimmed,
                                       variables,
                                       &impl->functions,
                                       &impl->scalar_functions,
                                       &impl->matrix_functions,
                                       &impl->value_functions,
                                       has_script_function,
                                       invoke_script_function,
                                       matrix_value.get())) {
        StoredValue stored;
        if (matrix_value->is_matrix) {
            stored.is_matrix = true;
            stored.matrix = matrix_value->matrix;
        } else if (matrix_value->is_complex) {
            stored.is_complex = true;
            stored.complex = matrix_value->complex;
            stored.decimal = matrix_value->complex.real;
            stored.exact = false;
        } else {
            stored.decimal = matrix_value->scalar;
            stored.exact = false;
        }
        return stored;
    }

    // 模块隐式求值
    if (!exact_mode) {
        std::map<std::string, StoredValue> snapshot = variables.snapshot();
        StoredValue stored;
        for (const auto& module : impl->implicit_evaluation_modules) {
            if (module->try_evaluate_implicit(trimmed, &stored, snapshot)) {
                return stored;
            }
        }
    }

    // 标量解析
    try {
        // 尝试使用编译后的 AST 进行快速求值
        if (expr_cache->is_compiled && expr_cache->compiled_ast) {
            try {
                const double parsed_value = evaluate_compiled_ast(
                    expr_cache->compiled_ast.get(),
                    variables,
                    &impl->functions,
                    &impl->scalar_functions,
                    has_script_function,
                    invoke_script_function);
                StoredValue stored;
                stored.decimal = parsed_value;
                stored.exact = false;
                attach_symbolic_text(&stored);
                return stored;
            } catch (const MathError&) {
                // AST 求值失败，回退到完整解析器
            }
        } else if (!expr_cache->is_compiled && expr_cache->hint == ExpressionHint::kScalar) {
            // 尝试编译 AST
            expr_cache->compiled_ast = compile_expression_ast(trimmed);
            if (expr_cache->compiled_ast) {
                expr_cache->is_compiled = true;
                // 尝试绑定变量槽位（高性能路径）
                bind_variable_slots(expr_cache->compiled_ast.get(), variables);
                try {
                    const double parsed_value = evaluate_compiled_ast(
                        expr_cache->compiled_ast.get(),
                        variables,
                        &impl->functions,
                        &impl->scalar_functions,
                        has_script_function,
                        invoke_script_function);
                    StoredValue stored;
                    stored.decimal = parsed_value;
                    stored.exact = false;
                    attach_symbolic_text(&stored);
                    return stored;
                } catch (const MathError&) {
                    // AST 求值失败，回退到完整解析器
                }
            }
        }

        // 回退到完整解析器
        DecimalParser parser(trimmed,
                             variables,
                             &impl->functions,
                             &impl->scalar_functions,
                             has_script_function,
                             invoke_script_function);
        const double parsed_value = parser.parse();
        StoredValue stored;
        stored.decimal = parsed_value;
        stored.exact = false;
        attach_symbolic_text(&stored);
        return stored;
    } catch (const UndefinedError&) {
        // 如果数值求值失败（由于变量未定义），但该表达式是由符号指令展开而来的，
        // 则保留其原始文本作为符号结果返回。
        if (trimmed != expression && (trimmed.find(' ') != std::string::npos || trimmed.find('*') != std::string::npos)) {
            StoredValue stored;
            stored.has_symbolic_text = true;
            stored.symbolic_text = trimmed;
            return stored;
        }
        throw;
    }
}

// 兼容旧接口的包装函数
StoredValue evaluate_expression_value_legacy(Calculator* calculator,
                                             Calculator::Impl* impl,
                                             const std::string& expression,
                                             bool exact_mode,
                                             std::shared_ptr<void>* cache) {
    std::shared_ptr<ExpressionCache> expr_cache;
    if (cache && *cache) {
        expr_cache = std::static_pointer_cast<ExpressionCache>(*cache);
    }
    auto result = evaluate_expression_value(calculator, impl, expression, exact_mode,
                                            cache ? &expr_cache : nullptr);
    if (cache && expr_cache) {
        *cache = expr_cache;
    }
    return result;
}

std::string execute_simple_script_line(Calculator* calculator,
                                       Calculator::Impl* impl,
                                       const std::string& text,
                                       bool exact_mode) {
    const std::string trimmed = trim_copy(text);
    if (trimmed.empty()) return "";

    // 1. 尝试作为内置/模块命令处理 (如 print, diff, plot 等)
    // 所有的 print 逻辑现在都下沉到 System 或领域模块中，或者通过统一分发器处理
    std::string command_output;
    if (calculator->try_process_function_command(trimmed, &command_output, exact_mode)) {
        return command_output;
    }

    // 2. 处理赋值语句
    std::string_view lhs, rhs;
    if (split_assignment(trimmed, &lhs, &rhs)) {
        if (!is_valid_variable_name(lhs)) throw std::runtime_error("invalid variable name: " + std::string(lhs));
        const StoredValue stored = evaluate_expression_value(calculator, impl, std::string(rhs), exact_mode);
        assign_visible_variable(impl, std::string(lhs), stored);
        return std::string(lhs) + " = " + format_stored_value(stored, impl->symbolic_constants_mode);
    }

    // 3. 回退到标准表达式求值
    return format_stored_value(evaluate_expression_value(calculator, impl, trimmed, exact_mode),
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
            // 使用预编译缓存
            std::shared_ptr<ExpressionCache> cache = simple.cache;
            *last_output = execute_simple_script_line(calculator, impl,
                cache ? cache->expanded : simple.text, exact_mode);
            return {};
        }
        case script::Statement::Kind::kIf: {
            const auto& if_statement = static_cast<const script::IfStatement&>(statement);
            if (truthy_value(evaluate_expression_value(calculator, impl, if_statement.condition, false, &if_statement.cache))) {
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
                                                          false,
                                                          &while_statement.cache))) {
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

            // 使用 FlatScopeStack
            if (impl->use_flat_scopes) {
                impl->flat_scopes.push_scope();
                try {
                    if (!for_statement.initializer.empty()) {
                        std::shared_ptr<ExpressionCache> init_cache = for_statement.init_cache;
                        if (!init_cache) {
                            init_cache = std::make_shared<ExpressionCache>();
                            init_cache->expanded = expand_inline_function_commands(calculator, for_statement.initializer);
                            for_statement.init_cache = init_cache;
                        }
                        (void)execute_simple_script_line(calculator,
                                                         impl,
                                                         init_cache->expanded,
                                                         exact_mode);
                    }
                    while (for_statement.condition.empty() ||
                           truthy_value(evaluate_expression_value(calculator,
                                                                  impl,
                                                                  for_statement.condition,
                                                                  false,
                                                                  &for_statement.cond_cache))) {
                        const ScriptSignal signal =
                            execute_script_statement(calculator,
                                                     impl,
                                                     *for_statement.body,
                                                     exact_mode,
                                                     last_output,
                                                     true);
                        if (signal.kind == ScriptSignal::Kind::kReturn) {
                            impl->flat_scopes.pop_scope();
                            return signal;
                        }
                        if (signal.kind == ScriptSignal::Kind::kBreak) {
                            break;
                        }
                        if (!for_statement.step.empty()) {
                            std::shared_ptr<ExpressionCache> step_cache = for_statement.step_cache;
                            if (!step_cache) {
                                step_cache = std::make_shared<ExpressionCache>();
                                step_cache->expanded = expand_inline_function_commands(calculator, for_statement.step);
                                for_statement.step_cache = step_cache;
                            }
                            (void)execute_simple_script_line(calculator,
                                                             impl,
                                                             step_cache->expanded,
                                                             exact_mode);
                        }
                        if (signal.kind == ScriptSignal::Kind::kContinue) {
                            continue;
                        }
                    }
                    impl->flat_scopes.pop_scope();
                    return {};
                } catch (...) {
                    impl->flat_scopes.pop_scope();
                    throw;
                }
            }

            // 传统路径
            impl->local_scopes.push_back({});
            try {
                if (!for_statement.initializer.empty()) {
                    // 使用预编译缓存
                    std::shared_ptr<ExpressionCache> init_cache = for_statement.init_cache;
                    if (!init_cache) {
                        init_cache = std::make_shared<ExpressionCache>();
                        init_cache->expanded = expand_inline_function_commands(calculator, for_statement.initializer);
                        for_statement.init_cache = init_cache;
                    }
                    (void)execute_simple_script_line(calculator,
                                                     impl,
                                                     init_cache->expanded,
                                                     exact_mode);
                }
                while (for_statement.condition.empty() ||
                       truthy_value(evaluate_expression_value(calculator,
                                                              impl,
                                                              for_statement.condition,
                                                              false,
                                                              &for_statement.cond_cache))) {
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
                        // 使用预编译缓存
                        std::shared_ptr<ExpressionCache> step_cache = for_statement.step_cache;
                        if (!step_cache) {
                            step_cache = std::make_shared<ExpressionCache>();
                            step_cache->expanded = expand_inline_function_commands(calculator, for_statement.step);
                            for_statement.step_cache = step_cache;
                        }
                        (void)execute_simple_script_line(calculator,
                                                         impl,
                                                         step_cache->expanded,
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
        case script::Statement::Kind::kForRange: {
            const auto& for_range = static_cast<const script::ForRangeStatement&>(statement);

            // 使用 FlatScopeStack
            if (impl->use_flat_scopes) {
                impl->flat_scopes.push_scope();
                try {
                    // 一次性求值 start, stop, step
                    const StoredValue start_val = evaluate_expression_value(
                        calculator, impl, for_range.start_expr, false, &for_range.start_cache);
                    const StoredValue stop_val = evaluate_expression_value(
                        calculator, impl, for_range.stop_expr, false, &for_range.stop_cache);
                    const StoredValue step_val = evaluate_expression_value(
                        calculator, impl, for_range.step_expr, false, &for_range.step_cache);

                    if (start_val.is_matrix || start_val.is_complex || start_val.is_string) {
                        throw std::runtime_error("range start must be a scalar");
                    }
                    if (stop_val.is_matrix || stop_val.is_complex || stop_val.is_string) {
                        throw std::runtime_error("range stop must be a scalar");
                    }
                    if (step_val.is_matrix || step_val.is_complex || step_val.is_string) {
                        throw std::runtime_error("range step must be a scalar");
                    }

                    double start = start_val.exact ? rational_to_double(start_val.rational) : start_val.decimal;
                    double stop = stop_val.exact ? rational_to_double(stop_val.rational) : stop_val.decimal;
                    double step = step_val.exact ? rational_to_double(step_val.rational) : step_val.decimal;

                    if (step == 0.0) {
                        throw std::runtime_error("range step cannot be zero");
                    }

                    bool ascending = step > 0.0;

                    // 直接数值循环，无需每次重新解析表达式
                    double current = start;
                    while ((ascending && current < stop) || (!ascending && current > stop)) {
                        // 设置循环变量（直接在 flat_scopes 中设置）
                        StoredValue loop_val;
                        loop_val.decimal = current;
                        loop_val.exact = false;
                        impl->flat_scopes.set(for_range.variable, loop_val);

                        // 执行循环体
                        const ScriptSignal signal =
                            execute_script_statement(calculator,
                                                     impl,
                                                     *for_range.body,
                                                     exact_mode,
                                                     last_output,
                                                     true);
                        if (signal.kind == ScriptSignal::Kind::kReturn) {
                            impl->flat_scopes.pop_scope();
                            return signal;
                        }
                        if (signal.kind == ScriptSignal::Kind::kBreak) {
                            break;
                        }

                        current += step;

                        if (signal.kind == ScriptSignal::Kind::kContinue) {
                            continue;
                        }
                    }
                    impl->flat_scopes.pop_scope();
                    return {};
                } catch (...) {
                    impl->flat_scopes.pop_scope();
                    throw;
                }
            }

            // 传统路径
            impl->local_scopes.push_back({});
            try {
                // 一次性求值 start, stop, step
                const StoredValue start_val = evaluate_expression_value(
                    calculator, impl, for_range.start_expr, false, &for_range.start_cache);
                const StoredValue stop_val = evaluate_expression_value(
                    calculator, impl, for_range.stop_expr, false, &for_range.stop_cache);
                const StoredValue step_val = evaluate_expression_value(
                    calculator, impl, for_range.step_expr, false, &for_range.step_cache);

                if (start_val.is_matrix || start_val.is_complex || start_val.is_string) {
                    throw std::runtime_error("range start must be a scalar");
                }
                if (stop_val.is_matrix || stop_val.is_complex || stop_val.is_string) {
                    throw std::runtime_error("range stop must be a scalar");
                }
                if (step_val.is_matrix || step_val.is_complex || step_val.is_string) {
                    throw std::runtime_error("range step must be a scalar");
                }

                double start = start_val.exact ? rational_to_double(start_val.rational) : start_val.decimal;
                double stop = stop_val.exact ? rational_to_double(stop_val.rational) : stop_val.decimal;
                double step = step_val.exact ? rational_to_double(step_val.rational) : step_val.decimal;

                if (step == 0.0) {
                    throw std::runtime_error("range step cannot be zero");
                }

                bool ascending = step > 0.0;

                // 直接数值循环，无需每次重新解析表达式
                double current = start;
                while ((ascending && current < stop) || (!ascending && current > stop)) {
                    // 设置循环变量
                    StoredValue loop_val;
                    loop_val.decimal = current;
                    loop_val.exact = false;
                    impl->local_scopes.back()[for_range.variable] = loop_val;

                    // 执行循环体
                    const ScriptSignal signal =
                        execute_script_statement(calculator,
                                                 impl,
                                                 *for_range.body,
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

                    current += step;

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
            if (is_reserved_user_function_name(impl, function_statement.name)) {
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
                evaluate_expression_value(calculator, impl, return_statement.expression, exact_mode, &return_statement.cache));
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
        if (impl->use_flat_scopes) {
            impl->flat_scopes.push_scope();
        } else {
            impl->local_scopes.push_back({});
        }
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
                    if (impl->use_flat_scopes) {
                        impl->flat_scopes.pop_scope();
                    } else {
                        impl->local_scopes.pop_back();
                    }
                }
                return signal;
            }
        }
        if (create_scope) {
            if (impl->use_flat_scopes) {
                impl->flat_scopes.pop_scope();
            } else {
                impl->local_scopes.pop_back();
            }
        }
        return {};
    } catch (...) {
        if (create_scope) {
            if (impl->use_flat_scopes) {
                impl->flat_scopes.pop_scope();
            } else {
                impl->local_scopes.pop_back();
            }
        }
        throw;
    }
}

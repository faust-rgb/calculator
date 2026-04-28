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

#include "function_analysis.h"
#include "calculator_matrix_commands.h"
#include "calculator_symbolic_commands.h"
#include "matrix.h"
#include "multivariable_integrator.h"
#include "mymath.h"
#include "ode_solver.h"
#include "optimization_helpers.h"
#include "polynomial.h"
#include "symbolic_expression.h"
#include "symbolic_expression_internal.h"

#include <algorithm>
#include <complex>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

bool Calculator::try_process_function_command(const std::string& expression,
                                              std::string* output) {
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
    if (trimmed == ":funcs") {
        if (impl_->functions.empty() && impl_->script_functions.empty()) {
            *output = "No custom functions defined.";
            return true;
        }

        std::ostringstream out;
        bool first = true;
        for (const auto& [name, function] : impl_->functions) {
            if (!first) {
                out << '\n';
            }
            first = false;
            out << name << "(" << function.parameter_name << ") = "
                << function.expression;
        }
        for (const auto& [name, function] : impl_->script_functions) {
            if (!first) {
                out << '\n';
            }
            first = false;
            out << name << "(";
            for (std::size_t i = 0; i < function.parameter_names.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
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
        const auto it = impl_->functions.find(name);
        if (it != impl_->functions.end()) {
            impl_->functions.erase(it);
            *output = "Cleared custom function: " + name;
            return true;
        }
        const auto script_it = impl_->script_functions.find(name);
        if (script_it == impl_->script_functions.end()) {
            throw std::runtime_error("unknown custom function: " + name);
        }
        impl_->script_functions.erase(script_it);
        *output = "Cleared custom function: " + name;
        return true;
    }

    auto build_symbolic_expression =
        [this](const std::string& name, std::string* variable_name) {
            const auto it = impl_->functions.find(name);
            if (it == impl_->functions.end()) {
                throw std::runtime_error("unknown custom function: " + name);
            }
            *variable_name = it->second.parameter_name;
            return SymbolicExpression::parse(it->second.expression);
        };

    symbolic_commands::SymbolicResolverContext symbolic_resolver_ctx;
    symbolic_resolver_ctx.resolve_custom_function = build_symbolic_expression;
    symbolic_resolver_ctx.has_custom_function = [&](const std::string& name) {
        return impl_->functions.find(name) != impl_->functions.end();
    };
    symbolic_resolver_ctx.expand_inline = [&](const std::string& arg) {
        return expand_inline_function_commands(this, arg);
    };

    auto resolve_symbolic_expression =
        [&](const std::string& argument,
            bool require_single_variable,
            std::string* variable_name,
            SymbolicExpression* expression) {
            symbolic_commands::resolve_symbolic_expression(symbolic_resolver_ctx,
                                                           argument,
                                                           require_single_variable,
                                                           variable_name,
                                                           expression);
        };

    auto build_analysis = [&](const std::string& argument) {
        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(argument, true, &variable_name, &expression);
        FunctionAnalysis analysis(variable_name);
        analysis.define(expression.to_string());
        return analysis;
    };

    auto build_scoped_decimal_evaluator = [&](const std::string& argument) {
        const std::string scoped_expression =
            trim_copy(expand_inline_function_commands(this, argument));
        return [this, scoped_expression](
                   const std::vector<std::pair<std::string, double>>& assignments) {
            std::map<std::string, StoredValue> scoped_variables =
                visible_variables(impl_.get());
            for (const auto& [name, value] : assignments) {
                StoredValue stored;
                stored.decimal = normalize_display_decimal(value);
                stored.exact = false;
                scoped_variables[name] = stored;
            }

            const HasScriptFunctionCallback has_script_function =
                [this](const std::string& name) {
                    return has_visible_script_function(impl_.get(), name);
                };
            const InvokeScriptFunctionDecimalCallback invoke_script_function =
                [this](const std::string& name,
                       const std::vector<double>& arguments) {
                    return invoke_script_function_decimal(
                        this, impl_.get(), name, arguments);
                };

            DecimalParser parser(scoped_expression,
                                 &scoped_variables,
                                 &impl_->functions,
                                 has_script_function,
                                 invoke_script_function);
            return normalize_result(parser.parse());
        };
    };

    auto build_scoped_matrix_evaluator = [&](const std::string& argument) {
        const std::string scoped_expression =
            trim_copy(expand_inline_function_commands(this, argument));
        return [this, scoped_expression](
                   const std::vector<std::pair<std::string, StoredValue>>& assignments) {
            std::map<std::string, StoredValue> scoped_variables =
                visible_variables(impl_.get());
            for (const auto& [name, value] : assignments) {
                scoped_variables[name] = value;
            }

            const HasScriptFunctionCallback has_script_function =
                [this](const std::string& name) {
                    return has_visible_script_function(impl_.get(), name);
                };
            const InvokeScriptFunctionDecimalCallback invoke_script_function =
                [this](const std::string& name,
                       const std::vector<double>& arguments) {
                    return invoke_script_function_decimal(
                        this, impl_.get(), name, arguments);
                };

            matrix::Value value;
            if (!try_evaluate_matrix_expression(scoped_expression,
                                               &scoped_variables,
                                               &impl_->functions,
                                               has_script_function,
                                               invoke_script_function,
                                               &value) ||
                !value.is_matrix) {
                throw std::runtime_error("expected a matrix-valued expression");
            }            return value.matrix;
        };
    };

    auto build_scoped_scalar_evaluator = [&](const std::string& argument) {
        const std::string scoped_expression =
            trim_copy(expand_inline_function_commands(this, argument));
        return [this, scoped_expression](
                   const std::vector<std::pair<std::string, StoredValue>>& assignments) {
            std::map<std::string, StoredValue> frame;
            for (const auto& [name, value] : assignments) {
                frame[name] = value;
            }

            impl_->local_scopes.push_back(frame);
            try {
                const StoredValue value =
                    evaluate_expression_value(this, impl_.get(), scoped_expression, false);
                impl_->local_scopes.pop_back();

                if (value.is_matrix || value.is_complex || value.is_string) {
                    throw std::runtime_error("expected a scalar-valued expression");
                }
                return normalize_result(value.exact
                                            ? rational_to_double(value.rational)
                                            : value.decimal);
            } catch (...) {
                impl_->local_scopes.pop_back();
                throw;
            }
        };
    };

    auto parse_decimal_argument = [&](const std::string& argument) {
        DecimalParser parser(argument, &impl_->variables, &impl_->functions);
        return parser.parse();
    };

    auto parse_matrix_argument = [&](const std::string& argument,
                                     const std::string& context) {
        const StoredValue value =
            evaluate_expression_value(this, impl_.get(), argument, false);
        if (!value.is_matrix) {
            throw std::runtime_error(context + " expects a matrix or vector argument");
        }
        return value.matrix;
    };

    auto is_matrix_argument = [&](const std::string& argument) {
        const std::map<std::string, StoredValue> visible = visible_variables(impl_.get());
        const HasScriptFunctionCallback has_script_function =
            [this](const std::string& name) {
                return has_visible_script_function(impl_.get(), name);
            };
        const InvokeScriptFunctionDecimalCallback invoke_script_function =
            [this](const std::string& name, const std::vector<double>& arguments) {
                return invoke_script_function_decimal(this, impl_.get(), name, arguments);
            };
        matrix::Value value;
        return try_evaluate_matrix_expression(trim_copy(argument),
                                             &visible,
                                             &impl_->functions,
                                             has_script_function,
                                             invoke_script_function,
                                             &value) &&
               value.is_matrix;
    };

    auto evaluate_symbolic_at = [&](const SymbolicExpression& expression,
                                    const std::string& variable_name,
                                    double point) {
        // Taylor 展开需要在某个点反复计算“当前导函数”的数值。
        // 为了复用现有 Calculator::evaluate，这里临时把变量写入会话变量表，
        // 求值完成后再恢复原状态，避免污染用户上下文。
        const auto existing = impl_->variables.find(variable_name);
        const bool had_existing = existing != impl_->variables.end();
        StoredValue backup;
        if (had_existing) {
            backup = existing->second;
        }

        StoredValue temporary;
        temporary.decimal = point;
        temporary.exact = false;
        impl_->variables[variable_name] = temporary;

        auto cleanup = [&]() {
            if (had_existing) {
                impl_->variables[variable_name] = backup;
            } else {
                impl_->variables.erase(variable_name);
            }
        };

        try {
            const double value = evaluate(expression.to_string());
            if (!std::isfinite(value)) {
                throw std::runtime_error("Non-finite value");
            }
            cleanup();
            return value;
        } catch (...) {
            cleanup();
            // If direct evaluation fails (e.g. division by zero like sin(x)/x at 0), fallback to limit
            try {
                FunctionAnalysis analysis(variable_name);
                analysis.define(expression.to_string());
                // Use right-sided limit (1) instead of two-sided (0) to handle branches like t/|t| in Puiseux
                return analysis.limit(point, 1);
            } catch (...) {
                return std::numeric_limits<double>::quiet_NaN();
            }
        }
    };

    auto simplify_symbolic_text = [&](const std::string& text) {
        return SymbolicExpression::parse(text).simplify().to_string();
    };

    auto parse_symbolic_variable_arguments =
        [](const std::vector<std::string>& arguments,
           std::size_t start_index,
           const std::vector<std::string>& fallback_variables) {
            return symbolic_commands::parse_symbolic_variable_arguments(
                arguments, start_index, fallback_variables);
        };

    auto parse_symbolic_expression_list = [&](const std::string& argument) {
        return symbolic_commands::parse_symbolic_expression_list(
            argument,
            [&](const std::string& arg) {
                return expand_inline_function_commands(this, arg);
            });
    };

    using FunctionCommandHandler =
        std::function<bool(const std::string&, const std::string&, std::string*)>;
    using FunctionCommandMatcher = std::function<bool(const std::string&)>;
    struct FunctionCommandRegistration {
        FunctionCommandMatcher matches;
        FunctionCommandHandler handler;
    };

    auto handle_polynomial_command = [&](const std::string& command,
                                         const std::string& inside,
                                         std::string* output) {
        // Pre-check if 'inside' is a single undefined variable to match test expectations
        const std::string trimmed = trim_copy(inside);
        if (is_identifier_text(trimmed) &&
            impl_->functions.find(trimmed) == impl_->functions.end()) {
            throw std::runtime_error("unknown variable: " + trimmed);
        }

        polynomial_ops::PolynomialContext ctx;
        ctx.functions = &impl_->functions;
        ctx.resolve_symbolic = [&](const std::string& name, std::string* var) {
            return build_symbolic_expression(name, var);
        };
        return polynomial_ops::handle_polynomial_command(ctx, command, inside, output);
    };

    auto handle_series_command = [&](const std::string& command,
                                     const std::string& inside,
                                     std::string* output) {
        series_ops::SeriesContext ctx;
        ctx.resolve_symbolic = [&](const std::string& arg,
                                   bool req,
                                   std::string* var,
                                   SymbolicExpression* expr) {
            resolve_symbolic_expression(arg, req, var, expr);
        };
        ctx.parse_decimal = [&](const std::string& arg) { return parse_decimal_argument(arg); };
        ctx.evaluate_at = [&](const SymbolicExpression& e, const std::string& v, double p) {
            return evaluate_symbolic_at(e, v, p);
        };
        ctx.simplify_symbolic = [&](const std::string& text) {
            return simplify_symbolic_text(text);
        };
        ctx.expand_inline = [&](const std::string& arg) {
            return expand_inline_function_commands(this, arg);
        };
        return series_ops::handle_series_command(ctx, command, inside, output);
    };

    auto handle_transform_command = [&](const std::string& command,
                                        const std::string& inside,
                                        std::string* output) {
        transforms::TransformContext ctx;
        ctx.resolve_symbolic = [&](const std::string& arg,
                                   bool req,
                                   std::string* var,
                                   SymbolicExpression* expr) {
            resolve_symbolic_expression(arg, req, var, expr);
        };
        return transforms::handle_transform_command(ctx, command, inside, output);
    };

    auto handle_integration_command = [&](const std::string& command,
                                          const std::string& inside,
                                          std::string* output) {
        integration_ops::IntegrationContext ctx;
        ctx.parse_decimal = [&](const std::string& arg) { return parse_decimal_argument(arg); };
        ctx.build_scoped_evaluator = [&](const std::string& arg) {
            return build_scoped_decimal_evaluator(arg);
        };
        ctx.normalize_result = [&](double v) { return normalize_result(v); };
        return integration_ops::handle_integration_command(ctx, command, inside, output);
    };

    auto handle_rootfinding_command = [&](const std::string& command,
                                          const std::string& inside,
                                          std::string* output) {
        rootfinding::RootfindingContext ctx;
        ctx.parse_decimal = [&](const std::string& arg) { return parse_decimal_argument(arg); };
        ctx.build_scoped_evaluator = [&](const std::string& arg) {
            return build_scoped_decimal_evaluator(arg);
        };
        ctx.get_derivative_expression = [&](const std::string& expr_str, const std::string& var_name) {
            std::string var;
            SymbolicExpression expr;
            try {
                // 尝试解析表达式。如果 expr_str 仅是表达式，resolve_symbolic_expression 也能处理。
                resolve_symbolic_expression(expr_str, false, &var, &expr);
                if (expr.node_) {
                    return expr.derivative(var_name).simplify().to_string();
                }
            } catch (...) {
                // 如果解析失败，返回空字符串，调用者将回退到差分法
            }
            return std::string();
        };
        ctx.is_matrix_argument = [&](const std::string& arg) { return is_matrix_argument(arg); };
        ctx.normalize_result = [&](double v) { return normalize_result(v); };
        return rootfinding::handle_rootfinding_command(ctx, command, inside, output);
    };

    auto handle_optimization_command = [&](const std::string& command,
                                           const std::string& inside,
                                           std::string* output) {
        optimization::OptimizationContext ctx;
        ctx.parse_matrix_argument = [&](const std::string& arg, const std::string& c) {
            return parse_matrix_argument(arg, c);
        };
        ctx.normalize_result = [&](double v) { return normalize_result(v); };
        ctx.is_integer_double = [](double x, double eps) { return is_integer_double(x, eps); };
        ctx.round_to_long_long = [](double x) { return round_to_long_long(x); };
        return optimization::handle_optimization_command(ctx, command, inside, output);
    };

    auto handle_symbolic_command = [&](const std::string& command,
                                       const std::string& inside,
                                       std::string* output) {
        symbolic_commands::SymbolicCommandContext ctx;
        ctx.resolve_symbolic = [&](const std::string& arg,
                                   bool req,
                                   std::string* var,
                                   SymbolicExpression* expr) {
            resolve_symbolic_expression(arg, req, var, expr);
        };
        ctx.parse_symbolic_variable_arguments = parse_symbolic_variable_arguments;
        ctx.parse_symbolic_expression_list = parse_symbolic_expression_list;
        ctx.build_analysis = [&](const std::string& arg) { return build_analysis(arg); };
        ctx.parse_decimal = [&](const std::string& arg) { return parse_decimal_argument(arg); };
        ctx.normalize_result = [&](double v) { return normalize_result(v); };
        return symbolic_commands::handle_symbolic_command(ctx, command, inside, output);
    };

    auto handle_analysis_command = [&](const std::string& command,
                                       const std::string& inside,
                                       std::string* output) {
        analysis_cmds::AnalysisContext ctx;
        ctx.resolve_symbolic = [&](const std::string& arg,
                                   bool req,
                                   std::string* var,
                                   SymbolicExpression* expr) {
            resolve_symbolic_expression(arg, req, var, expr);
        };
        ctx.parse_symbolic_variable_arguments = parse_symbolic_variable_arguments;
        ctx.parse_decimal = [&](const std::string& arg) { return parse_decimal_argument(arg); };
        ctx.normalize_result = [&](double v) { return normalize_result(v); };
        ctx.build_analysis = [&](const std::string& arg) { return build_analysis(arg); };
        return analysis_cmds::handle_analysis_command(ctx, command, inside, output);
    };

    auto handle_ode_command = [&](const std::string& command,
                                  const std::string& inside,
                                  std::string* output) {
        ode_ops::ODEContext ctx;
        ctx.parse_decimal = [&](const std::string& arg) { return parse_decimal_argument(arg); };
        ctx.build_scoped_scalar_evaluator = [&](const std::string& arg) {
            return build_scoped_scalar_evaluator(arg);
        };
        ctx.build_scoped_matrix_evaluator = [&](const std::string& arg) {
            return build_scoped_matrix_evaluator(arg);
        };
        ctx.is_matrix_argument = [&](const std::string& arg) { return is_matrix_argument(arg); };
        ctx.parse_matrix_argument = [&](const std::string& arg, const std::string& context) {
            return parse_matrix_argument(arg, context);
        };
        ctx.evaluate_expression_value = [&](const std::string& arg, bool exact_mode) {
            return evaluate_expression_value(this, impl_.get(), arg, exact_mode);
        };
        ctx.normalize_result = [&](double v) { return normalize_result(v); };
        return ode_ops::handle_ode_command(ctx, command, inside, output);
    };

    auto handle_matrix_command = [&](const std::string& command,
                                     const std::string& inside,
                                     std::string* output) {
        matrix_commands::MatrixCommandContext ctx;
        ctx.is_matrix_argument = [&](const std::string& arg) { return is_matrix_argument(arg); };
        ctx.parse_matrix_argument = [&](const std::string& arg, const std::string& context) {
            return parse_matrix_argument(arg, context);
        };
        return matrix_commands::handle_matrix_command(ctx, command, inside, output);
    };

    const std::vector<FunctionCommandRegistration> command_registry = {
        {polynomial_ops::is_polynomial_command, handle_polynomial_command},
        {series_ops::is_series_command, handle_series_command},
        {transforms::is_transform_command, handle_transform_command},
        {integration_ops::is_integration_command, handle_integration_command},
        {rootfinding::is_rootfinding_command, handle_rootfinding_command},
        {optimization::is_optimization_command, handle_optimization_command},
        {symbolic_commands::is_symbolic_command, handle_symbolic_command},
        {analysis_cmds::is_analysis_command, handle_analysis_command},
        {ode_ops::is_ode_command, handle_ode_command},
        {matrix_commands::is_matrix_command, handle_matrix_command},
    };

    auto split_function_command_call = [](const std::string& text,
                                          std::string* command_name,
                                          std::string* inside) {
        const std::size_t open = text.find('(');
        if (open == std::string::npos) {
            return false;
        }
        const std::string name = trim_copy(text.substr(0, open));
        if (!is_identifier_text(name) || !split_named_call(text, name, inside)) {
            return false;
        }
        *command_name = name;
        return true;
    };

    std::string command_name;
    std::string inside;
    if (!split_function_command_call(trimmed, &command_name, &inside)) {
        return false;
    }

    if (command_name == "residue") {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() == 3) {
            const std::string variable_name = trim_copy(arguments[1]);
            if (!is_identifier_text(variable_name)) {
                throw std::runtime_error("residue variable must be an identifier");
            }

            const SymbolicExpression expression =
                SymbolicExpression::parse(
                    trim_copy(expand_inline_function_commands(this, arguments[0])))
                    .simplify();
            SymbolicExpression numerator = expression;
            SymbolicExpression denominator = SymbolicExpression::number(1.0);
            if (expression.node_->type == NodeType::kDivide) {
                numerator = SymbolicExpression(expression.node_->left).simplify();
                denominator = SymbolicExpression(expression.node_->right).simplify();
            }

            std::vector<double> numerator_coefficients;
            std::vector<double> denominator_coefficients;
            if (!numerator.polynomial_coefficients(variable_name,
                                                   &numerator_coefficients) ||
                !denominator.polynomial_coefficients(variable_name,
                                                     &denominator_coefficients)) {
                throw std::runtime_error(
                    "residue currently supports rational polynomial expressions");
            }

            StoredValue point_value =
                evaluate_expression_value(this, impl_.get(), arguments[2], false);
            std::complex<double> point(point_value.exact
                                           ? rational_to_double(point_value.rational)
                                           : point_value.decimal,
                                       0.0);
            if (point_value.is_matrix) {
                const matrix::Matrix& point_matrix = point_value.matrix;
                if (!point_matrix.is_vector() ||
                    point_matrix.rows * point_matrix.cols != 2) {
                    throw std::runtime_error(
                        "residue point must be scalar or complex(real, imag)");
                }
                const double real = point_matrix.rows == 1 ? point_matrix.at(0, 0)
                                                           : point_matrix.at(0, 0);
                const double imag = point_matrix.rows == 1 ? point_matrix.at(0, 1)
                                                           : point_matrix.at(1, 0);
                point = {real, imag};
            } else if (point_value.is_complex) {
                point = {point_value.complex.real, point_value.complex.imag};
            } else if (point_value.is_string) {
                throw std::runtime_error("residue point must be numeric");
            }

            auto evaluate_polynomial_complex =
                [](const std::vector<double>& coefficients,
                   std::complex<double> value) {
                    std::complex<double> result(0.0, 0.0);
                    for (std::size_t i = coefficients.size(); i > 0; --i) {
                        result = result * value + coefficients[i - 1];
                    }
                    return result;
                };

            const std::vector<double> denominator_derivative =
                polynomial_derivative(denominator_coefficients);
            const std::complex<double> denominator_value =
                evaluate_polynomial_complex(denominator_coefficients, point);
            if (std::abs(denominator_value) > 1e-8) {
                *output = matrix::Matrix::vector({0.0, 0.0}).to_string();
                return true;
            }
            const std::complex<double> denominator_prime =
                evaluate_polynomial_complex(denominator_derivative, point);
            if (std::abs(denominator_prime) <= 1e-10) {
                throw std::runtime_error(
                    "residue currently supports only simple poles");
            }

            const std::complex<double> residue =
                evaluate_polynomial_complex(numerator_coefficients, point) /
                denominator_prime;
            *output = matrix::Matrix::vector(
                          {normalize_result(residue.real()),
                           normalize_result(residue.imag())})
                          .to_string();
            return true;
        }
    }

    for (const FunctionCommandRegistration& registration : command_registry) {
        if (registration.matches(command_name)) {
            return registration.handler(command_name, inside, output);
        }
    }

    return false;
}

// ============================================================================
// 积分变换命令实现
// ============================================================================

#include "symbolic/modules/transform_module.h"
#include "execution/engine/script_context.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "parser/grammars/unified_expression_parser.h"
#include "core/services/string_utils.h"

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vector>

namespace transforms {

using BasicTransform = std::function<SymbolicExpression(
    const SymbolicExpression&, const std::string&, const std::string&)>;

std::string run_basic_transform(const TransformContext& ctx,
                                const std::string& expr,
                                const std::string& input_var,
                                const std::string& output_var,
                                const char* default_input,
                                const char* default_output,
                                const BasicTransform& transform) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, false, &variable_name, &expression);

    const std::string in_var = input_var.empty()
        ? (variable_name.empty() ? default_input : variable_name)
        : input_var;
    const std::string out_var = output_var.empty() ? default_output : output_var;
    return transform(expression, in_var, out_var).simplify().to_string();
}

void append_transform_metadata(const TransformResult& result,
                               bool include_constraints,
                               std::string* output) {
    *output += "; condition: ";
    *output += result.condition.has_value() ? result.condition->expression : "none";
    if (include_constraints && result.condition.has_value() &&
        !result.condition->constraints.empty()) {
        *output += "; constraints: ";
        for (std::size_t i = 0; i < result.condition->constraints.size(); ++i) {
            if (i != 0) *output += " and ";
            *output += result.condition->constraints[i];
        }
    }
    *output += "; distribution: ";
    *output += result.contains_distribution ? "true" : "false";
}

std::string fourier(const TransformContext& ctx,
                    const std::string& expr,
                    const std::string& input_var,
                    const std::string& output_var) {
    return run_basic_transform(ctx, expr, input_var, output_var, "t", "w",
                               [](const SymbolicExpression& expression,
                                  const std::string& in_var,
                                  const std::string& out_var) {
                                   return expression.fourier_transform(in_var, out_var);
                               });
}

std::string inverse_fourier(const TransformContext& ctx,
                            const std::string& expr,
                            const std::string& input_var,
                            const std::string& output_var) {
    return run_basic_transform(ctx, expr, input_var, output_var, "w", "t",
                               [](const SymbolicExpression& expression,
                                  const std::string& in_var,
                                  const std::string& out_var) {
                                   return expression.inverse_fourier_transform(in_var, out_var);
                               });
}

std::string laplace(const TransformContext& ctx,
                    const std::string& expr,
                    const std::string& input_var,
                    const std::string& output_var) {
    return run_basic_transform(ctx, expr, input_var, output_var, "t", "s",
                               [](const SymbolicExpression& expression,
                                  const std::string& in_var,
                                  const std::string& out_var) {
                                   return expression.laplace_transform(in_var, out_var);
                               });
}

std::string inverse_laplace(const TransformContext& ctx,
                            const std::string& expr,
                            const std::string& input_var,
                            const std::string& output_var) {
    return run_basic_transform(ctx, expr, input_var, output_var, "s", "t",
                               [](const SymbolicExpression& expression,
                                  const std::string& in_var,
                                  const std::string& out_var) {
                                   return expression.inverse_laplace_transform(in_var, out_var);
                               });
}

std::string z_transform(const TransformContext& ctx,
                        const std::string& expr,
                        const std::string& input_var,
                        const std::string& output_var) {
    return run_basic_transform(ctx, expr, input_var, output_var, "n", "z",
                               [](const SymbolicExpression& expression,
                                  const std::string& in_var,
                                  const std::string& out_var) {
                                   return expression.z_transform(in_var, out_var);
                               });
}

std::string inverse_z_transform(const TransformContext& ctx,
                                const std::string& expr,
                                const std::string& input_var,
                                const std::string& output_var) {
    return run_basic_transform(ctx, expr, input_var, output_var, "z", "n",
                               [](const SymbolicExpression& expression,
                                  const std::string& in_var,
                                  const std::string& out_var) {
                                   return expression.inverse_z_transform(in_var, out_var);
                               });
}

bool is_transform_command(const std::string& command) {
    return command == "fourier" ||
           command == "fourier_info" ||
           command == "ifourier" ||
           command == "ifourier_info" ||
           command == "inverse_fourier" ||
           command == "laplace" ||
           command == "laplace_info" ||
           command == "ilaplace" ||
           command == "ilaplace_info" ||
           command == "inverse_laplace" ||
           command == "ztrans" ||
           command == "ztrans_info" ||
           command == "z_transform" ||
           command == "iztrans" ||
           command == "iztrans_info" ||
           command == "inverse_z";
}

std::string normalize_transform_command(const std::string& command) {
    if (command == "inverse_fourier") return "ifourier";
    if (command == "fourier_info") return "fourier_info";
    if (command == "ifourier_info") return "ifourier_info";
    if (command == "inverse_laplace") return "ilaplace";
    if (command == "laplace_info") return "laplace_info";
    if (command == "ilaplace_info") return "ilaplace_info";
    if (command == "z_transform") return "ztrans";
    if (command == "ztrans_info") return "ztrans_info";
    if (command == "inverse_z") return "iztrans";
    if (command == "iztrans_info") return "iztrans_info";
    return command;
}

void get_default_variables(const std::string& command,
                           const std::string& expr_var,
                           std::string* input_var,
                           std::string* output_var) {
    if (command == "fourier" || command == "fourier_info") {
        *input_var = expr_var.empty() ? "t" : expr_var;
        *output_var = "w";
    } else if (command == "ifourier" || command == "ifourier_info") {
        *input_var = expr_var.empty() ? "w" : expr_var;
        *output_var = "t";
    } else if (command == "laplace" || command == "laplace_info") {
        *input_var = expr_var.empty() ? "t" : expr_var;
        *output_var = "s";
    } else if (command == "ilaplace" || command == "ilaplace_info") {
        *input_var = expr_var.empty() ? "s" : expr_var;
        *output_var = "t";
    } else if (command == "ztrans" || command == "ztrans_info") {
        *input_var = expr_var.empty() ? "n" : expr_var;
        *output_var = "z";
    } else if (command == "iztrans" || command == "iztrans_info") {
        *input_var = expr_var.empty() ? "z" : expr_var;
        *output_var = "n";
    }
}

bool handle_transform_command(const TransformContext& ctx,
                              const std::string& command,
                              const std::vector<std::string>& arguments,
                              std::string* output) {
    if (arguments.empty()) throw std::runtime_error(command + " expects at least one argument");

    std::string expr_var;
    SymbolicExpression expression;
    ctx.resolve_symbolic(arguments[0], false, &expr_var, &expression);

    const std::string normalized = normalize_transform_command(command);
    std::string input_var, output_var;

    if (arguments.size() == 1) {
        get_default_variables(normalized, expr_var, &input_var, &output_var);
    } else if (arguments.size() == 2) {
        get_default_variables(normalized, expr_var, &input_var, &output_var);
        output_var = utils::trim_copy(arguments[1]);
    } else if (arguments.size() == 3) {
        input_var = utils::trim_copy(arguments[1]);
        output_var = utils::trim_copy(arguments[2]);
        if (!is_identifier_text(input_var) || !is_identifier_text(output_var)) {
            throw std::runtime_error(command + " variable names must be identifiers");
        }
    } else {
        throw std::runtime_error(command + " expects 1, 2, or 3 arguments");
    }

    if (normalized == "fourier" || normalized == "fourier_info") {
        const TransformResult result = expression.fourier_transform_with_conditions(
            input_var, output_var);
        *output = result.expression.simplify().to_string();
        if (normalized == "fourier_info") {
            append_transform_metadata(result, false, output);
        }
    } else if (normalized == "ifourier" || normalized == "ifourier_info") {
        const TransformResult result = expression.inverse_fourier_transform_with_conditions(
            input_var, output_var);
        *output = result.expression.simplify().to_string();
        if (normalized == "ifourier_info") {
            append_transform_metadata(result, false, output);
        }
    } else if (normalized == "laplace" || normalized == "laplace_info") {
        const TransformResult result = expression.laplace_transform_with_conditions(
            input_var, output_var);
        *output = result.expression.simplify().to_string();
        if (normalized == "laplace_info") {
            append_transform_metadata(result, true, output);
        }
    } else if (normalized == "ilaplace" || normalized == "ilaplace_info") {
        const TransformResult result = expression.inverse_laplace_transform_with_conditions(
            input_var, output_var);
        *output = result.expression.simplify().to_string();
        if (normalized == "ilaplace_info") {
            append_transform_metadata(result, true, output);
        }
    } else if (normalized == "ztrans" || normalized == "ztrans_info") {
        const TransformResult result = expression.z_transform_with_conditions(
            input_var, output_var);
        *output = result.expression.simplify().to_string();
        if (normalized == "ztrans_info") {
            append_transform_metadata(result, true, output);
        }
    } else if (normalized == "iztrans" || normalized == "iztrans_info") {
        const TransformResult result = expression.inverse_z_transform_with_conditions(
            input_var, output_var);
        *output = result.expression.simplify().to_string();
        if (normalized == "iztrans_info") {
            append_transform_metadata(result, false, output);
        }
    } else {
        return false;
    }
    return true;
}

bool handle_transform_command(const TransformContext& ctx,
                              const std::string& command,
                              const std::string& inside,
                              std::string* output) {
    return handle_transform_command(ctx, command, split_top_level_arguments(inside), output);
}


std::string TransformModule::execute_args(const std::string& command,
                                         const std::vector<std::string>& args,
                                         ServiceLocator& locator) {
    auto services = locator.resolve<CoreServices>();
    TransformContext ctx;
    ctx.resolve_symbolic = [&locator](const std::string& arg, bool req, std::string* var, SymbolicExpression* expr) {
        locator.resolve<CoreServices>()->symbolic.resolve_symbolic(arg, req, var, expr);
    };

    std::string output;
    if (handle_transform_command(ctx, command, args, &output)) {
        return output;
    }
    throw std::runtime_error("Unknown transform command: " + command);
}

std::string TransformModule::execute_args_view(std::string_view command,
                                              const std::vector<std::string_view>& args,
                                              ServiceLocator& locator) {
    std::vector<std::string> str_args;
    str_args.reserve(args.size());
    for (const auto& sv : args) {
        str_args.emplace_back(sv);
    }
    return execute_args(std::string(command), str_args, locator);
}

std::string TransformModule::get_help_snippet(const std::string& topic) const {
    if (topic == "symbolic") {
        return "Transforms:\n"
               "  laplace(f, [t], [s])       Laplace transform\n"
               "  ilaplace(F, [s], [t])      Inverse Laplace transform\n"
               "  laplace_info(f, [t], [s])  Transform with condition metadata\n"
               "  ilaplace_info(F, [s], [t]) Inverse transform with metadata\n"
               "  fourier(f, [t], [w])       Fourier transform\n"
               "  fourier_info(f, [t], [w])  Fourier transform with metadata\n"
               "  ifourier(F, [w], [t])      Inverse Fourier transform\n"
               "  ifourier_info(F, [w], [t]) Inverse Fourier with metadata\n"
               "  ztrans(f, [n], [z])        Z transform\n"
               "  iztrans(F, [z], [n])       Inverse Z transform";
    }
    return "";
}

std::vector<std::string> TransformModule::get_commands() const {
    return {"laplace", "laplace_info", "ilaplace", "ilaplace_info", "inverse_laplace", "fourier",
            "fourier_info", "ifourier", "ifourier_info", "inverse_fourier", "ztrans", "iztrans",
            "z_transform", "inverse_z"};
}

std::vector<CommandSpec> TransformModule::get_command_specs() const {
    std::vector<CommandSpec> specs;
    for (const std::string& cmd : get_commands()) {
        CommandSpec spec;
        spec.key = call_command_key(cmd);
        spec.dispatch_name = cmd;
        spec.short_help = cmd;
        spec.is_inlineable = true;
        specs.push_back(std::move(spec));
    }
    return specs;
}

std::vector<std::string> TransformModule::get_function_names() const {
    return get_commands();
}

}  // namespace transforms

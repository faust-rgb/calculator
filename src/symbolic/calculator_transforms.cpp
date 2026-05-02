// ============================================================================
// 积分变换命令实现
// ============================================================================

#include "symbolic/calculator_transforms.h"
#include "core/string_utils.h"

#include <vector>
#include <algorithm>

namespace transforms {

std::string fourier(const TransformContext& ctx,
                    const std::string& expr,
                    const std::string& input_var,
                    const std::string& output_var) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, false, &variable_name, &expression);

    std::string in_var = input_var.empty() ? (variable_name.empty() ? "t" : variable_name) : input_var;
    std::string out_var = output_var.empty() ? "w" : output_var;

    return expression.fourier_transform(in_var, out_var)
               .simplify()
               .to_string();
}

std::string inverse_fourier(const TransformContext& ctx,
                            const std::string& expr,
                            const std::string& input_var,
                            const std::string& output_var) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, false, &variable_name, &expression);

    std::string in_var = input_var.empty() ? (variable_name.empty() ? "w" : variable_name) : input_var;
    std::string out_var = output_var.empty() ? "t" : output_var;

    return expression.inverse_fourier_transform(in_var, out_var)
               .simplify()
               .to_string();
}

std::string laplace(const TransformContext& ctx,
                    const std::string& expr,
                    const std::string& input_var,
                    const std::string& output_var) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, false, &variable_name, &expression);

    std::string in_var = input_var.empty() ? (variable_name.empty() ? "t" : variable_name) : input_var;
    std::string out_var = output_var.empty() ? "s" : output_var;

    return expression.laplace_transform(in_var, out_var)
               .simplify()
               .to_string();
}

std::string inverse_laplace(const TransformContext& ctx,
                            const std::string& expr,
                            const std::string& input_var,
                            const std::string& output_var) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, false, &variable_name, &expression);

    std::string in_var = input_var.empty() ? (variable_name.empty() ? "s" : variable_name) : input_var;
    std::string out_var = output_var.empty() ? "t" : output_var;

    return expression.inverse_laplace_transform(in_var, out_var)
               .simplify()
               .to_string();
}

std::string z_transform(const TransformContext& ctx,
                        const std::string& expr,
                        const std::string& input_var,
                        const std::string& output_var) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, false, &variable_name, &expression);

    std::string in_var = input_var.empty() ? (variable_name.empty() ? "n" : variable_name) : input_var;
    std::string out_var = output_var.empty() ? "z" : output_var;

    return expression.z_transform(in_var, out_var)
               .simplify()
               .to_string();
}

std::string inverse_z_transform(const TransformContext& ctx,
                                const std::string& expr,
                                const std::string& input_var,
                                const std::string& output_var) {
    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(expr, false, &variable_name, &expression);

    std::string in_var = input_var.empty() ? (variable_name.empty() ? "z" : variable_name) : input_var;
    std::string out_var = output_var.empty() ? "n" : output_var;

    return expression.inverse_z_transform(in_var, out_var)
               .simplify()
               .to_string();
}

bool is_transform_command(const std::string& command) {
    return command == "fourier" ||
           command == "ifourier" ||
           command == "inverse_fourier" ||
           command == "laplace" ||
           command == "ilaplace" ||
           command == "inverse_laplace" ||
           command == "ztrans" ||
           command == "z_transform" ||
           command == "iztrans" ||
           command == "inverse_z";
}

std::string normalize_transform_command(const std::string& command) {
    if (command == "inverse_fourier") return "ifourier";
    if (command == "inverse_laplace") return "ilaplace";
    if (command == "z_transform") return "ztrans";
    if (command == "inverse_z") return "iztrans";
    return command;
}

void get_default_variables(const std::string& command,
                           const std::string& expr_var,
                           std::string* input_var,
                           std::string* output_var) {
    if (command == "fourier") {
        *input_var = expr_var.empty() ? "t" : expr_var;
        *output_var = "w";
    } else if (command == "ifourier") {
        *input_var = expr_var.empty() ? "w" : expr_var;
        *output_var = "t";
    } else if (command == "laplace") {
        *input_var = expr_var.empty() ? "t" : expr_var;
        *output_var = "s";
    } else if (command == "ilaplace") {
        *input_var = expr_var.empty() ? "s" : expr_var;
        *output_var = "t";
    } else if (command == "ztrans") {
        *input_var = expr_var.empty() ? "n" : expr_var;
        *output_var = "z";
    } else if (command == "iztrans") {
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

    if (normalized == "fourier") {
        *output = expression.fourier_transform(input_var, output_var).simplify().to_string();
    } else if (normalized == "ifourier") {
        *output = expression.inverse_fourier_transform(input_var, output_var).simplify().to_string();
    } else if (normalized == "laplace") {
        *output = expression.laplace_transform(input_var, output_var).simplify().to_string();
    } else if (normalized == "ilaplace") {
        *output = expression.inverse_laplace_transform(input_var, output_var).simplify().to_string();
    } else if (normalized == "ztrans") {
        *output = expression.z_transform(input_var, output_var).simplify().to_string();
    } else if (normalized == "iztrans") {
        *output = expression.inverse_z_transform(input_var, output_var).simplify().to_string();
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
                                         const CoreServices& services) {
    TransformContext ctx;
    ctx.resolve_symbolic = services.symbolic.resolve_symbolic;

    std::string output;
    if (handle_transform_command(ctx, command, args, &output)) {
        return output;
    }
    throw std::runtime_error("Unknown transform command: " + command);
}

std::string TransformModule::get_help_snippet(const std::string& topic) const {
    if (topic == "symbolic") {
        return "Transforms:\n"
               "  laplace(f, [s], [t])       Laplace transform\n"
               "  ilaplace(F, [t], [s])      Inverse Laplace transform\n"
               "  fourier(f, [w], [t])       Fourier transform\n"
               "  ifourier(F, [t], [w])      Inverse Fourier transform\n"
               "  ztrans(f, [z], [n])        Z transform\n"
               "  iztrans(F, [n], [z])       Inverse Z transform";
    }
    return "";
}

std::vector<std::string> TransformModule::get_commands() const {
    return {"laplace", "ilaplace", "inverse_laplace", "fourier", "ifourier",
            "inverse_fourier", "ztrans", "iztrans", "z_transform", "inverse_z"};
}

}  // namespace transforms

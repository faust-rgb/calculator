// ============================================================================
// 积分变换命令实现
// ============================================================================
//
// 本文件实现各种积分变换的符号计算：
//
// 1. Fourier 变换族
//    - fourier: 正向 Fourier 变换
//    - inverse_fourier: 逆向 Fourier 变换
//
// 2. Laplace 变换族
//    - laplace: 正向 Laplace 变换
//    - inverse_laplace: 逆向 Laplace 变换
//
// 3. Z 变换族
//    - z_transform: 正向 Z 变换
//    - inverse_z_transform: 逆向 Z 变换
//
// 各变换函数支持指定输入/输出变量名，
// 默认变量根据变换类型自动选择（如 Fourier: t→w, Laplace: t→s, Z: n→z）。
// ============================================================================

#include "calculator_transforms.h"

#include <vector>

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
                              const std::string& inside,
                              std::string* output) {
    const std::string normalized = normalize_transform_command(command);
    const std::vector<std::string> arguments = split_top_level_arguments(inside);

    if (arguments.size() != 1 && arguments.size() != 3) {
        throw std::runtime_error(
            command +
            " expects either one symbolic expression or expression plus input/output variable names");
    }

    std::string variable_name;
    SymbolicExpression expression;
    ctx.resolve_symbolic(arguments[0], false, &variable_name, &expression);

    std::string input_var;
    std::string output_var;
    if (arguments.size() == 3) {
        input_var = trim_copy(arguments[1]);
        output_var = trim_copy(arguments[2]);
        if (!is_identifier_text(input_var) || !is_identifier_text(output_var)) {
            throw std::runtime_error(command + " variable names must be identifiers");
        }
    } else {
        get_default_variables(normalized, variable_name, &input_var, &output_var);
    }

    if (normalized == "fourier") {
        *output = expression.fourier_transform(input_var, output_var)
                      .simplify()
                      .to_string();
    } else if (normalized == "ifourier") {
        *output = expression.inverse_fourier_transform(input_var, output_var)
                      .simplify()
                      .to_string();
    } else if (normalized == "laplace") {
        *output = expression.laplace_transform(input_var, output_var)
                      .simplify()
                      .to_string();
    } else if (normalized == "ilaplace") {
        *output = expression.inverse_laplace_transform(input_var, output_var)
                      .simplify()
                      .to_string();
    } else if (normalized == "ztrans") {
        *output = expression.z_transform(input_var, output_var)
                      .simplify()
                      .to_string();
    } else if (normalized == "iztrans") {
        *output = expression.inverse_z_transform(input_var, output_var)
                      .simplify()
                      .to_string();
    } else {
        return false;
    }
    return true;
}

}  // namespace transforms

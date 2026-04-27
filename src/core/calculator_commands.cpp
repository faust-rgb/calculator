#include "calculator_internal_types.h"
#include "calculator_command_helpers.h"
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
#include "matrix.h"
#include "multivariable_integrator.h"
#include "mymath.h"
#include "ode_solver.h"
#include "optimization_helpers.h"
#include "polynomial.h"
#include "symbolic_expression.h"

#include <algorithm>
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

    std::function<void(const std::string&, std::string*, std::vector<double>*)>
        build_polynomial = [&](const std::string& argument,
                               std::string* variable_name,
                               std::vector<double>* coefficients) {
            const std::string trimmed_argument = trim_copy(argument);
            std::string nested_inside;
            if (split_named_call(trimmed_argument, "poly_add", &nested_inside) ||
                split_named_call(trimmed_argument, "poly_sub", &nested_inside) ||
                split_named_call(trimmed_argument, "poly_mul", &nested_inside) ||
                split_named_call(trimmed_argument, "poly_div", &nested_inside)) {
                const std::vector<std::string> nested_arguments =
                    split_top_level_arguments(nested_inside);
                if (nested_arguments.size() != 2) {
                    throw std::runtime_error(
                        "polynomial operations expect exactly two arguments");
                }

                std::string lhs_variable;
                std::string rhs_variable;
                std::vector<double> lhs_coefficients;
                std::vector<double> rhs_coefficients;
                build_polynomial(nested_arguments[0], &lhs_variable, &lhs_coefficients);
                build_polynomial(nested_arguments[1], &rhs_variable, &rhs_coefficients);

                *variable_name = lhs_variable;
                if (trimmed_argument.rfind("poly_add", 0) == 0) {
                    *coefficients = polynomial_add(lhs_coefficients, rhs_coefficients);
                    return;
                }
                if (trimmed_argument.rfind("poly_sub", 0) == 0) {
                    *coefficients = polynomial_subtract(lhs_coefficients, rhs_coefficients);
                    return;
                }
                if (trimmed_argument.rfind("poly_div", 0) == 0) {
                    const PolynomialDivisionResult division =
                        polynomial_divide(lhs_coefficients, rhs_coefficients);
                    bool zero_remainder = true;
                    for (double coefficient : division.remainder) {
                        if (!mymath::is_near_zero(coefficient, 1e-10)) {
                            zero_remainder = false;
                            break;
                        }
                    }
                    if (!zero_remainder) {
                        throw std::runtime_error(
                            "nested poly_div requires zero remainder");
                    }
                    *coefficients = division.quotient;
                    return;
                }

                *coefficients = polynomial_multiply(lhs_coefficients, rhs_coefficients);
                return;
            }

            // 命令层不直接处理原始字符串，而是统一走：
            // 自定义函数名 -> 符号表达式 -> 多项式系数
            //
            // 这样 poly_add/poly_mul/roots 等功能只依赖一个稳定的中间表示，
            // 不需要各自重复做“它是不是多项式”的判断。
            SymbolicExpression expression =
                build_symbolic_expression(trimmed_argument, variable_name);
            if (!expression.polynomial_coefficients(*variable_name, coefficients)) {
                throw std::runtime_error("custom function " + trimmed_argument +
                                         " is not a polynomial");
            }
        };

    std::function<void(const std::string&, bool, std::string*, SymbolicExpression*)>
        resolve_symbolic_expression =
            [&](const std::string& argument,
                bool require_single_variable,
                std::string* variable_name,
                SymbolicExpression* expression) {
                const std::string trimmed_argument = trim_copy(argument);
                std::string nested_inside;
                if (split_named_call(trimmed_argument, "diff", &nested_inside)) {
                    const std::vector<std::string> nested_arguments =
                        split_top_level_arguments(nested_inside);
                    if (nested_arguments.empty()) {
                        throw std::runtime_error(
                            "nested symbolic diff expects at least one argument");
                    }
                    SymbolicExpression nested_expression;
                    resolve_symbolic_expression(nested_arguments[0],
                                                nested_arguments.size() == 1,
                                                variable_name,
                                                &nested_expression);
                    if (nested_arguments.size() == 1) {
                        *expression =
                            nested_expression.derivative(*variable_name).simplify();
                    } else {
                        SymbolicExpression differentiated = nested_expression;
                        for (std::size_t i = 1; i < nested_arguments.size(); ++i) {
                            const std::string derivative_variable =
                                trim_copy(nested_arguments[i]);
                            if (!is_identifier_text(derivative_variable)) {
                                throw std::runtime_error(
                                    "nested symbolic diff variable arguments must be identifiers");
                            }
                            differentiated =
                                differentiated.derivative(derivative_variable).simplify();
                        }
                        *variable_name = trim_copy(nested_arguments[1]);
                        *expression = differentiated;
                    }
                    return;
                }
                if (split_named_call(trimmed_argument, "integral", &nested_inside)) {
                    const std::vector<std::string> nested_arguments =
                        split_top_level_arguments(nested_inside);
                    if (nested_arguments.size() != 1 && nested_arguments.size() != 2) {
                        throw std::runtime_error(
                            "nested symbolic integral expects expression and optional variable");
                    }
                    SymbolicExpression nested_expression;
                    resolve_symbolic_expression(nested_arguments[0],
                                                nested_arguments.size() == 1,
                                                variable_name,
                                                &nested_expression);
                    if (nested_arguments.size() == 2) {
                        const std::string integral_variable =
                            trim_copy(nested_arguments[1]);
                        if (!is_identifier_text(integral_variable)) {
                            throw std::runtime_error(
                                "nested symbolic integral variable must be an identifier");
                        }
                        *variable_name = integral_variable;
                    }
                    *expression =
                        nested_expression.integral(*variable_name).simplify();
                    return;
                }
                if (split_named_call(trimmed_argument, "poly_add", &nested_inside) ||
                    split_named_call(trimmed_argument, "poly_sub", &nested_inside) ||
                    split_named_call(trimmed_argument, "poly_mul", &nested_inside) ||
                    split_named_call(trimmed_argument, "poly_div", &nested_inside)) {
                    std::vector<double> coefficients;
                    build_polynomial(trimmed_argument, variable_name, &coefficients);
                    *expression = SymbolicExpression::parse(
                        polynomial_to_string(coefficients, *variable_name));
                    return;
                }
                const auto function_it = impl_->functions.find(trimmed_argument);
                if (function_it != impl_->functions.end()) {
                    *expression = build_symbolic_expression(trimmed_argument, variable_name);
                    return;
                }
                *expression = SymbolicExpression::parse(
                    expand_inline_function_commands(this, trimmed_argument));
                const std::vector<std::string> identifiers =
                    expression->identifier_variables();
                if (identifiers.size() == 1) {
                    *variable_name = identifiers[0];
                } else if (identifiers.empty()) {
                    *variable_name = "x";
                } else if (require_single_variable) {
                    throw std::runtime_error(
                        "symbolic expressions with multiple variables must use a custom function");
                } else {
                    *variable_name = identifiers.front();
                }
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
            }
            return value.matrix;
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

                if (value.is_matrix || value.is_string) {
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

    auto matrix_to_vector_values = [&](const matrix::Matrix& value,
                                       const std::string& context) {
        if (!value.is_vector()) {
            throw std::runtime_error(context + " expects vector arguments");
        }
        const std::size_t size = value.rows == 1 ? value.cols : value.rows;
        std::vector<double> result(size, 0.0);
        for (std::size_t i = 0; i < size; ++i) {
            result[i] = value.rows == 1 ? value.at(0, i) : value.at(i, 0);
        }
        return result;
    };

    auto vector_to_column_matrix = [&](const std::vector<double>& values) {
        matrix::Matrix result(values.size(), 1, 0.0);
        for (std::size_t i = 0; i < values.size(); ++i) {
            result.at(i, 0) = normalize_result(values[i]);
        }
        return result;
    };

    auto make_scalar_stored = [&](double value) {
        StoredValue stored;
        stored.decimal = normalize_result(value);
        return stored;
    };

    auto append_parameter_assignments =
        [&](const StoredValue& parameter_value,
            std::vector<std::pair<std::string, StoredValue>>* assignments) {
            assignments->push_back({"p", parameter_value});
            if (!parameter_value.is_matrix || !parameter_value.matrix.is_vector()) {
                return;
            }

            const std::size_t size =
                parameter_value.matrix.rows == 1
                    ? parameter_value.matrix.cols
                    : parameter_value.matrix.rows;
            for (std::size_t i = 0; i < size; ++i) {
                const double component_value =
                    parameter_value.matrix.rows == 1
                        ? parameter_value.matrix.at(0, i)
                        : parameter_value.matrix.at(i, 0);
                assignments->push_back({"p" + std::to_string(i + 1),
                                        make_scalar_stored(component_value)});
            }
        };

    auto try_parse_positive_step_argument = [&](const std::string& argument,
                                                int* steps) {
        try {
            const double value = parse_decimal_argument(argument);
            if (!is_integer_double(value) || value <= 0.0) {
                return false;
            }
            *steps = static_cast<int>(round_to_long_long(value));
            return true;
        } catch (const std::exception&) {
            return false;
        }
    };

    auto dot_product = [](const std::vector<double>& lhs,
                          const std::vector<double>& rhs) {
        return optimization_helpers::dot_product(lhs, rhs);
    };

    auto format_planning_result = [&](const std::vector<double>& solution,
                                      double objective) {
        return "x = " + matrix::Matrix::vector(solution).to_string() +
               "\nobjective = " + format_decimal(normalize_result(objective));
    };

    auto solve_linear_box_problem =
        [&](const std::vector<double>& objective,
            const matrix::Matrix& inequality_matrix,
            const std::vector<double>& inequality_rhs,
            const matrix::Matrix& equality_matrix,
            const std::vector<double>& equality_rhs,
            const std::vector<double>& lower_bounds,
            const std::vector<double>& upper_bounds,
            double planning_tolerance,
            std::vector<double>* solution,
            double* objective_value,
            std::string* diagnostic) {
            const std::size_t variable_count = objective.size();
            if (inequality_matrix.cols != variable_count ||
                inequality_rhs.size() != inequality_matrix.rows ||
                equality_matrix.cols != variable_count ||
                equality_rhs.size() != equality_matrix.rows ||
                lower_bounds.size() != variable_count ||
                upper_bounds.size() != variable_count) {
                throw std::runtime_error("planning dimension mismatch");
            }

            auto feasible_solution = [&](const std::vector<double>& x) {
                if (x.size() != variable_count) {
                    return false;
                }
                for (std::size_t col = 0; col < variable_count; ++col) {
                    if (x[col] < lower_bounds[col] - planning_tolerance ||
                        x[col] > upper_bounds[col] + planning_tolerance) {
                        return false;
                    }
                }
                for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                    long double total = 0.0L;
                    for (std::size_t col = 0; col < variable_count; ++col) {
                        total += static_cast<long double>(inequality_matrix.at(row, col)) *
                                 static_cast<long double>(x[col]);
                    }
                    if (total >
                        static_cast<long double>(inequality_rhs[row]) + planning_tolerance) {
                        return false;
                    }
                }
                for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                    long double total = 0.0L;
                    for (std::size_t col = 0; col < variable_count; ++col) {
                        total += static_cast<long double>(equality_matrix.at(row, col)) *
                                 static_cast<long double>(x[col]);
                    }
                    if (mymath::abs(static_cast<double>(total - equality_rhs[row])) >
                        planning_tolerance) {
                        return false;
                    }
                }
                return true;
            };

            if (variable_count == 0) {
                const std::vector<double> empty_solution;
                if (!feasible_solution(empty_solution)) {
                    return false;
                }
                *solution = empty_solution;
                *objective_value = 0.0;
                return true;
            }

            static constexpr std::size_t kMaxPlanningBasisCount = 250000;
            std::vector<std::vector<double>> all_rows;
            std::vector<double> all_rhs;
            all_rows.reserve(inequality_matrix.rows + 2 * equality_matrix.rows +
                             2 * variable_count);
            all_rhs.reserve(inequality_matrix.rows + 2 * equality_matrix.rows +
                            2 * variable_count);

            for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                std::vector<double> coefficients(variable_count, 0.0);
                for (std::size_t col = 0; col < variable_count; ++col) {
                    coefficients[col] = inequality_matrix.at(row, col);
                }
                all_rows.push_back(coefficients);
                all_rhs.push_back(inequality_rhs[row]);
            }
            for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                std::vector<double> positive_row(variable_count, 0.0);
                std::vector<double> negative_row(variable_count, 0.0);
                for (std::size_t col = 0; col < variable_count; ++col) {
                    positive_row[col] = equality_matrix.at(row, col);
                    negative_row[col] = -equality_matrix.at(row, col);
                }
                all_rows.push_back(positive_row);
                all_rhs.push_back(equality_rhs[row]);
                all_rows.push_back(negative_row);
                all_rhs.push_back(-equality_rhs[row]);
            }
            for (std::size_t col = 0; col < variable_count; ++col) {
                std::vector<double> upper_row(variable_count, 0.0);
                upper_row[col] = 1.0;
                all_rows.push_back(upper_row);
                all_rhs.push_back(upper_bounds[col]);

                std::vector<double> lower_row(variable_count, 0.0);
                lower_row[col] = -1.0;
                all_rows.push_back(lower_row);
                all_rhs.push_back(-lower_bounds[col]);
            }

            auto capped_choose = [](std::size_t n, std::size_t k, std::size_t cap) {
                if (k > n) {
                    return std::size_t{0};
                }
                if (k > n - k) {
                    k = n - k;
                }

                long double value = 1.0L;
                for (std::size_t i = 1; i <= k; ++i) {
                    value *= static_cast<long double>(n - k + i);
                    value /= static_cast<long double>(i);
                    if (value > static_cast<long double>(cap)) {
                        return cap + 1;
                    }
                }
                return static_cast<std::size_t>(value + 0.5L);
            };

            const std::size_t basis_count_estimate =
                capped_choose(all_rows.size(), variable_count, kMaxPlanningBasisCount);
            if (basis_count_estimate > kMaxPlanningBasisCount) {
                throw std::runtime_error(
                    "planning basis enumeration limit exceeded: " +
                    std::to_string(basis_count_estimate) +
                    "+ candidate bases for " + std::to_string(variable_count) +
                    " variables and " + std::to_string(all_rows.size()) +
                    " active constraints");
            }

            bool found = false;
            std::vector<double> best_solution(variable_count, 0.0);
            double best_value = 0.0;
            std::vector<std::size_t> selection(variable_count, 0);
            std::size_t checked_bases = 0;
            std::size_t singular_bases = 0;
            std::size_t infeasible_bases = 0;

            std::function<void(std::size_t, std::size_t)> enumerate_bases =
                [&](std::size_t start, std::size_t depth) {
                    if (depth == variable_count) {
                        ++checked_bases;
                        matrix::Matrix basis(variable_count, variable_count, 0.0);
                        matrix::Matrix basis_rhs(variable_count, 1, 0.0);
                        for (std::size_t row = 0; row < variable_count; ++row) {
                            basis_rhs.at(row, 0) = all_rhs[selection[row]];
                            for (std::size_t col = 0; col < variable_count; ++col) {
                                basis.at(row, col) = all_rows[selection[row]][col];
                            }
                        }

                        try {
                            const matrix::Matrix solved = matrix::solve(basis, basis_rhs);
                            std::vector<double> candidate(variable_count, 0.0);
                            for (std::size_t i = 0; i < variable_count; ++i) {
                                const double value = solved.at(i, 0);
                                candidate[i] =
                                    mymath::abs(value) <= planning_tolerance ? 0.0 : value;
                            }
                            if (!feasible_solution(candidate)) {
                                ++infeasible_bases;
                                return;
                            }

                            const double candidate_objective =
                                dot_product(objective, candidate);
                            if (!found ||
                                candidate_objective > best_value + planning_tolerance) {
                                found = true;
                                best_value = candidate_objective;
                                best_solution = candidate;
                            }
                        } catch (const std::exception&) {
                            ++singular_bases;
                        }
                        return;
                    }

                    const std::size_t remaining = variable_count - depth;
                    for (std::size_t index = start;
                         index + remaining <= all_rows.size();
                         ++index) {
                        selection[depth] = index;
                        enumerate_bases(index + 1, depth + 1);
                    }
                };

            enumerate_bases(0, 0);
            if (!found) {
                if (diagnostic != nullptr) {
                    *diagnostic = "checked " + std::to_string(checked_bases) +
                                  " bases, skipped " + std::to_string(singular_bases) +
                                  " singular bases and " +
                                  std::to_string(infeasible_bases) +
                                  " infeasible candidates";
                }
                return false;
            }

            *solution = best_solution;
            *objective_value = best_value;
            return true;
        };

    auto parse_subdivisions = [&](const std::vector<std::string>& arguments,
                                  std::size_t offset,
                                  const std::vector<int>& defaults) {
        std::vector<int> subdivisions = defaults;
        if (arguments.size() == offset) {
            return subdivisions;
        }
        if (arguments.size() != offset + defaults.size()) {
            throw std::runtime_error("unexpected subdivision argument count");
        }

        for (std::size_t i = 0; i < defaults.size(); ++i) {
            const double value = parse_decimal_argument(arguments[offset + i]);
            if (!is_integer_double(value) || value <= 0.0) {
                throw std::runtime_error(
                    "integration subdivision counts must be positive integers");
            }
            subdivisions[i] = static_cast<int>(round_to_long_long(value));
        }
        return subdivisions;
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

        try {
            const double value = evaluate(expression.to_string());
            if (had_existing) {
                impl_->variables[variable_name] = backup;
            } else {
                impl_->variables.erase(variable_name);
            }
            return value;
        } catch (...) {
            if (had_existing) {
                impl_->variables[variable_name] = backup;
            } else {
                impl_->variables.erase(variable_name);
            }
            throw;
        }
    };

    auto build_taylor_coefficients = [&](const SymbolicExpression& expression,
                                         const std::string& variable_name,
                                         double center,
                                         int degree) {
        struct TaylorDerivativeCacheEntry {
            SymbolicExpression derivative;
            double value = 0.0;
            bool has_value = false;
        };
        static thread_local std::map<std::string, TaylorDerivativeCacheEntry> derivative_cache;
        static constexpr std::size_t kMaxTaylorDerivativeCacheSize = 256;

        const std::string base_key =
            variable_name + "|" + format_symbolic_scalar(center) + "|" +
            expression.simplify().to_string();
        std::vector<double> coefficients;
        coefficients.reserve(static_cast<std::size_t>(degree + 1));
        SymbolicExpression current = expression;
        for (int order = 0; order <= degree; ++order) {
            const std::string order_key = base_key + "|" + std::to_string(order);
            auto found = derivative_cache.find(order_key);
            if (found == derivative_cache.end()) {
                if (derivative_cache.size() >= kMaxTaylorDerivativeCacheSize) {
                    derivative_cache.clear();
                }
                TaylorDerivativeCacheEntry entry;
                entry.derivative = current.simplify();
                found = derivative_cache.emplace(order_key, entry).first;
            } else {
                current = found->second.derivative;
            }

            if (!found->second.has_value) {
                found->second.value =
                    evaluate_symbolic_at(found->second.derivative,
                                         variable_name,
                                         center);
                found->second.has_value = true;
            }
            const double derivative_value = found->second.value;
            coefficients.push_back(derivative_value / factorial_value(order));
            if (order != degree) {
                const std::string next_key = base_key + "|" + std::to_string(order + 1);
                auto next_found = derivative_cache.find(next_key);
                if (next_found != derivative_cache.end()) {
                    current = next_found->second.derivative;
                } else {
                    current = found->second.derivative.derivative(variable_name).simplify();
                }
            }
        }
        return coefficients;
    };

    auto simplify_symbolic_text = [&](const std::string& text) {
        return SymbolicExpression::parse(text).simplify().to_string();
    };

    auto symbolic_vector_to_string =
        [](const std::vector<SymbolicExpression>& values) {
            std::ostringstream out;
            out << "[";
            for (std::size_t i = 0; i < values.size(); ++i) {
                if (i != 0) {
                    out << ", ";
                }
                out << values[i].simplify().to_string();
            }
            out << "]";
            return out.str();
        };

    auto symbolic_matrix_to_string =
        [](const std::vector<std::vector<SymbolicExpression>>& values) {
            std::ostringstream out;
            out << "[";
            for (std::size_t row = 0; row < values.size(); ++row) {
                if (row != 0) {
                    out << ", ";
                }
                out << "[";
                for (std::size_t col = 0; col < values[row].size(); ++col) {
                    if (col != 0) {
                        out << ", ";
                    }
                    out << values[row][col].simplify().to_string();
                }
                out << "]";
            }
            out << "]";
            return out.str();
        };

    auto evaluate_symbolic_at_point =
        [](SymbolicExpression expression,
           const std::vector<std::string>& variables,
           const std::vector<double>& values,
           double* result) {
            for (std::size_t i = 0; i < variables.size(); ++i) {
                expression = expression.substitute(
                    variables[i], SymbolicExpression::number(values[i]));
            }
            return expression.simplify().is_number(result);
        };

    auto solve_linear_system =
        [](std::vector<std::vector<double>> matrix,
           std::vector<double> rhs,
           std::vector<double>* solution) {
            const std::size_t size = rhs.size();
            for (std::size_t col = 0; col < size; ++col) {
                std::size_t pivot = col;
                for (std::size_t row = col + 1; row < size; ++row) {
                    if (std::abs(matrix[row][col]) > std::abs(matrix[pivot][col])) {
                        pivot = row;
                    }
                }
                if (std::abs(matrix[pivot][col]) < 1e-12) {
                    return false;
                }
                if (pivot != col) {
                    std::swap(matrix[pivot], matrix[col]);
                    std::swap(rhs[pivot], rhs[col]);
                }

                const double pivot_value = matrix[col][col];
                for (std::size_t j = col; j < size; ++j) {
                    matrix[col][j] /= pivot_value;
                }
                rhs[col] /= pivot_value;

                for (std::size_t row = 0; row < size; ++row) {
                    if (row == col) {
                        continue;
                    }
                    const double factor = matrix[row][col];
                    if (std::abs(factor) < 1e-12) {
                        continue;
                    }
                    for (std::size_t j = col; j < size; ++j) {
                        matrix[row][j] -= factor * matrix[col][j];
                    }
                    rhs[row] -= factor * rhs[col];
                }
            }
            *solution = rhs;
            return true;
        };

    auto parse_symbolic_variable_arguments =
        [](const std::vector<std::string>& arguments,
           std::size_t start_index,
           const std::vector<std::string>& fallback_variables) {
            std::vector<std::string> variables;
            for (std::size_t i = start_index; i < arguments.size(); ++i) {
                const std::string variable = trim_copy(arguments[i]);
                if (!is_identifier_text(variable)) {
                    throw std::runtime_error(
                        "symbolic variable arguments must be identifiers");
                }
                variables.push_back(variable);
            }
            if (variables.empty()) {
                variables = fallback_variables;
            }
            if (variables.empty()) {
                variables.push_back("x");
            }
            return variables;
        };

    auto parse_symbolic_expression_list = [&](const std::string& argument) {
        std::string text = trim_copy(argument);
        if (text.size() < 2 || text.front() != '[' || text.back() != ']') {
            throw std::runtime_error(
                "jacobian expects its first argument to be a bracketed expression list");
        }
        text = trim_copy(text.substr(1, text.size() - 2));

        std::vector<std::string> expression_texts;
        int paren_depth = 0;
        int bracket_depth = 0;
        std::size_t start = 0;
        for (std::size_t i = 0; i < text.size(); ++i) {
            const char ch = text[i];
            if (ch == '(') {
                ++paren_depth;
            } else if (ch == ')') {
                --paren_depth;
            } else if (ch == '[') {
                ++bracket_depth;
            } else if (ch == ']') {
                --bracket_depth;
            } else if ((ch == ';' || ch == ',') &&
                       paren_depth == 0 &&
                       bracket_depth == 0) {
                expression_texts.push_back(trim_copy(text.substr(start, i - start)));
                start = i + 1;
            }
        }
        if (!text.empty()) {
            expression_texts.push_back(trim_copy(text.substr(start)));
        }
        if (expression_texts.empty()) {
            throw std::runtime_error("jacobian expression list cannot be empty");
        }

        std::vector<SymbolicExpression> expressions;
        expressions.reserve(expression_texts.size());
        for (const std::string& expression_text : expression_texts) {
            if (expression_text.empty()) {
                throw std::runtime_error("jacobian expression list cannot contain empty items");
            }
            expressions.push_back(SymbolicExpression::parse(
                expand_inline_function_commands(this, expression_text)));
        }
        return expressions;
    };

    using FunctionCommandHandler =
        std::function<bool(const std::string&, const std::string&, std::string*)>;
    struct FunctionCommandRegistration {
        std::vector<std::string> names;
        FunctionCommandHandler handler;
    };

    auto handle_polynomial_command = [&](const std::string& command,
                                         const std::string& inside,
                                         std::string* output) {
        polynomial_ops::PolynomialContext ctx;
        ctx.functions = &impl_->functions;
        ctx.resolve_symbolic = [&](const std::string& name, std::string* var) {
            return build_symbolic_expression(name, var);
        };
        return polynomial_ops::handle_polynomial_command(ctx, command, inside, output);
    };

    auto handle_series_module_command = [&](const std::string& command,
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
        ctx.build_taylor_coefficients = [&](const SymbolicExpression& e,
                                            const std::string& v,
                                            double c,
                                            int d) {
            return build_taylor_coefficients(e, v, c, d);
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
        ctx.parse_subdivisions = [&](const std::vector<std::string>& args,
                                     std::size_t off,
                                     const std::vector<int>& def) {
            return parse_subdivisions(args, off, def);
        };
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
        ctx.matrix_to_vector_values = [&](const matrix::Matrix& m, const std::string& c) {
            return matrix_to_vector_values(m, c);
        };
        ctx.dot_product = [&](const std::vector<double>& lhs, const std::vector<double>& rhs) {
            return dot_product(lhs, rhs);
        };
        ctx.format_planning_result = [&](const std::vector<double>& sol, double obj) {
            return format_planning_result(sol, obj);
        };
        ctx.solve_linear_box_problem = [&](const std::vector<double>& obj,
                                           const matrix::Matrix& im,
                                           const std::vector<double>& ir,
                                           const matrix::Matrix& em,
                                           const std::vector<double>& er,
                                           const std::vector<double>& lb,
                                           const std::vector<double>& ub,
                                           double tol,
                                           std::vector<double>* sol,
                                           double* val,
                                           std::string* diag) {
            return solve_linear_box_problem(obj, im, ir, em, er, lb, ub, tol, sol, val, diag);
        };
        ctx.is_integer_double = [](double x, double eps) { return is_integer_double(x, eps); };
        ctx.round_to_long_long = [](double x) { return round_to_long_long(x); };
        return optimization::handle_optimization_command(ctx, command, inside, output);
    };

    auto handle_series_sum_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 4) {
            throw std::runtime_error(
                "series_sum expects expr, index, lower, upper");
        }

        const std::string index_name = trim_copy(arguments[1]);
        if (!is_identifier_text(index_name)) {
            throw std::runtime_error("series_sum index must be an identifier");
        }

        SymbolicExpression summand =
            SymbolicExpression::parse(expand_inline_function_commands(this, arguments[0]));
        SymbolicExpression upper_expression;
        SymbolicExpression lower_expression;
        const std::string upper_text = trim_copy(arguments[3]);
        const bool upper_is_infinite =
            upper_text == "inf" || upper_text == "oo" ||
            upper_text == "infinity";
        if (!upper_is_infinite) {
            upper_expression = SymbolicExpression::parse(
                expand_inline_function_commands(this, upper_text));
        }
        lower_expression = SymbolicExpression::parse(
            expand_inline_function_commands(this, arguments[2]));

        auto make_polynomial_sum_primitive =
            [&](const std::vector<double>& coefficients) {
                if (coefficients.size() > 4) {
                    throw std::runtime_error(
                        "series_sum polynomial summands are currently supported up to degree 3");
                }

                std::vector<std::string> pieces;
                if (coefficients.size() >= 1 &&
                    !mymath::is_near_zero(coefficients[0], 1e-10)) {
                    pieces.push_back("(" + format_symbolic_scalar(coefficients[0]) +
                                     ") * (" + index_name + " + 1)");
                }
                if (coefficients.size() >= 2 &&
                    !mymath::is_near_zero(coefficients[1], 1e-10)) {
                    pieces.push_back("(" + format_symbolic_scalar(coefficients[1]) +
                                     ") * (" + index_name + " * (" + index_name +
                                     " + 1) / 2)");
                }
                if (coefficients.size() >= 3 &&
                    !mymath::is_near_zero(coefficients[2], 1e-10)) {
                    pieces.push_back("(" + format_symbolic_scalar(coefficients[2]) +
                                     ") * (" + index_name + " * (" + index_name +
                                     " + 1) * (2 * " + index_name + " + 1) / 6)");
                }
                if (coefficients.size() >= 4 &&
                    !mymath::is_near_zero(coefficients[3], 1e-10)) {
                    pieces.push_back("(" + format_symbolic_scalar(coefficients[3]) +
                                     ") * ((" + index_name + " * (" + index_name +
                                     " + 1) / 2) ^ 2)");
                }

                if (pieces.empty()) {
                    return SymbolicExpression::number(0.0);
                }
                std::ostringstream out;
                for (std::size_t i = 0; i < pieces.size(); ++i) {
                    if (i != 0) {
                        out << " + ";
                    }
                    out << pieces[i];
                }
                return SymbolicExpression::parse(out.str()).simplify();
            };

        auto finite_sum_from_primitive =
            [&](const SymbolicExpression& primitive) {
                const SymbolicExpression lower_minus_one =
                    SymbolicExpression::parse("(" + arguments[2] + ") - 1").simplify();
                return SymbolicExpression::parse(
                           "(" +
                           primitive.substitute(index_name, upper_expression)
                               .to_string() +
                           ") - (" +
                           primitive.substitute(index_name, lower_minus_one)
                               .to_string() +
                           ")")
                    .simplify()
                    .to_string();
            };

        std::vector<double> polynomial_coefficients;
        if (summand.polynomial_coefficients(index_name, &polynomial_coefficients)) {
            if (upper_is_infinite) {
                bool all_zero = true;
                for (double coefficient : polynomial_coefficients) {
                    if (!mymath::is_near_zero(coefficient, 1e-10)) {
                        all_zero = false;
                        break;
                    }
                }
                if (!all_zero) {
                    throw std::runtime_error(
                        "series_sum does not support infinite polynomial sums");
                }
                *output = "0";
                return true;
            }

            const SymbolicExpression primitive =
                make_polynomial_sum_primitive(polynomial_coefficients);
            *output = finite_sum_from_primitive(primitive);
            return true;
        }

        auto geometric_ratio = [&](double* coefficient, double* ratio) {
            const double s0 = evaluate_symbolic_at(summand, index_name, 0.0);
            const double s1 = evaluate_symbolic_at(summand, index_name, 1.0);
            const double s2 = evaluate_symbolic_at(summand, index_name, 2.0);
            const double s3 = evaluate_symbolic_at(summand, index_name, 3.0);
            if (mymath::is_near_zero(s0, 1e-10)) {
                return false;
            }
            const double candidate = s1 / s0;
            if (!mymath::is_near_zero(s2 - s1 * candidate, 1e-8) ||
                !mymath::is_near_zero(s3 - s2 * candidate, 1e-8)) {
                return false;
            }
            *coefficient = s0;
            *ratio = candidate;
            return true;
        };

        double geometric_coefficient = 0.0;
        double geometric_ratio_value = 0.0;
        if (!geometric_ratio(&geometric_coefficient, &geometric_ratio_value)) {
            throw std::runtime_error(
                "series_sum currently supports polynomial summands up to degree 3 and common geometric series");
        }

        const std::string coefficient_text =
            format_symbolic_scalar(geometric_coefficient);
        const std::string ratio_text =
            format_symbolic_scalar(geometric_ratio_value);

        if (upper_is_infinite) {
            if (mymath::abs(geometric_ratio_value) >= 1.0 - 1e-10) {
                throw std::runtime_error(
                    "series_sum infinite geometric series requires |r| < 1");
            }
            if (mymath::is_near_zero(geometric_ratio_value - 1.0, 1e-10)) {
                throw std::runtime_error(
                    "series_sum infinite geometric series diverges for r = 1");
            }
            *output = simplify_symbolic_text(
                "(" + coefficient_text + ") * (" + ratio_text + ") ^ (" +
                arguments[2] + ") / (1 - (" + ratio_text + "))");
            return true;
        }

        const std::string geometric_primitive_text =
            mymath::is_near_zero(geometric_ratio_value - 1.0, 1e-10)
                ? "(" + coefficient_text + ") * (" + index_name + " + 1)"
                : "(" + coefficient_text + ") * (1 - (" + ratio_text +
                      ") ^ (" + index_name + " + 1)) / (1 - (" +
                      ratio_text + "))";
        const SymbolicExpression primitive =
            SymbolicExpression::parse(geometric_primitive_text).simplify();
        *output = finite_sum_from_primitive(primitive);
        return true;
    };

    auto handle_simplify_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::string argument = trim_copy(inside);
        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(argument, false, &variable_name, &expression);
        *output = expression.simplify().to_string();
        return true;
    };

    auto handle_gradient_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.empty()) {
            throw std::runtime_error(
                "gradient expects a symbolic expression and optional variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            parse_symbolic_variable_arguments(arguments,
                                              1,
                                              expression.identifier_variables());
        *output = symbolic_vector_to_string(expression.gradient(variables));
        return true;
    };

    auto handle_hessian_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.empty()) {
            throw std::runtime_error(
                "hessian expects a symbolic expression and optional variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            parse_symbolic_variable_arguments(arguments,
                                              1,
                                              expression.identifier_variables());
        *output = symbolic_matrix_to_string(expression.hessian(variables));
        return true;
    };

    auto handle_critical_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.empty()) {
            throw std::runtime_error(
                "critical expects a symbolic expression and optional variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], false, &variable_name, &expression);
        const std::vector<std::string> variables =
            parse_symbolic_variable_arguments(arguments,
                                              1,
                                              expression.identifier_variables());
        const std::vector<SymbolicExpression> gradient =
            expression.gradient(variables);
        const std::vector<std::vector<SymbolicExpression>> hessian =
            expression.hessian(variables);
        auto classify_critical_point =
            [&](const std::vector<double>& values) {
                const std::size_t n = variables.size();
                std::vector<std::vector<double>> evaluated(
                    n, std::vector<double>(n, 0.0));
                for (std::size_t row = 0; row < n; ++row) {
                    for (std::size_t col = 0; col < n; ++col) {
                        if (!evaluate_symbolic_at_point(hessian[row][col],
                                                        variables,
                                                        values,
                                                        &evaluated[row][col])) {
                            return std::string("unclassified");
                        }
                    }
                }

                if (n == 1) {
                    if (evaluated[0][0] > 1e-8) {
                        return std::string("local min");
                    }
                    if (evaluated[0][0] < -1e-8) {
                        return std::string("local max");
                    }
                    return std::string("degenerate");
                }
                if (n == 2) {
                    const double det =
                        evaluated[0][0] * evaluated[1][1] -
                        evaluated[0][1] * evaluated[1][0];
                    if (det > 1e-8 && evaluated[0][0] > 1e-8) {
                        return std::string("local min");
                    }
                    if (det > 1e-8 && evaluated[0][0] < -1e-8) {
                        return std::string("local max");
                    }
                    if (det < -1e-8) {
                        return std::string("saddle");
                    }
                    return std::string("degenerate");
                }
                if (n == 3) {
                    const double d1 = evaluated[0][0];
                    const double d2 =
                        evaluated[0][0] * evaluated[1][1] -
                        evaluated[0][1] * evaluated[1][0];
                    const double d3 =
                        evaluated[0][0] *
                            (evaluated[1][1] * evaluated[2][2] -
                             evaluated[1][2] * evaluated[2][1]) -
                        evaluated[0][1] *
                            (evaluated[1][0] * evaluated[2][2] -
                             evaluated[1][2] * evaluated[2][0]) +
                        evaluated[0][2] *
                            (evaluated[1][0] * evaluated[2][1] -
                             evaluated[1][1] * evaluated[2][0]);
                    if (d1 > 1e-8 && d2 > 1e-8 && d3 > 1e-8) {
                        return std::string("local min");
                    }
                    if (d1 < -1e-8 && d2 > 1e-8 && d3 < -1e-8) {
                        return std::string("local max");
                    }
                    if (mymath::abs(d1) <= 1e-8 ||
                        mymath::abs(d2) <= 1e-8 ||
                        mymath::abs(d3) <= 1e-8) {
                        return std::string("degenerate");
                    }
                    return std::string("saddle");
                }
                return std::string("unclassified");
            };
        auto format_critical_solution =
            [&](const std::vector<double>& values) {
                std::ostringstream out;
                out << "[";
                for (std::size_t i = 0; i < variables.size(); ++i) {
                    if (i != 0) {
                        out << ", ";
                    }
                    out << variables[i] << " = "
                        << format_decimal(normalize_result(values[i]));
                }
                out << "]";
                out << " (" << classify_critical_point(values) << ")";
                return out.str();
            };
        std::vector<std::vector<double>> coefficients(
            variables.size(), std::vector<double>(variables.size(), 0.0));
        std::vector<double> rhs(variables.size(), 0.0);
        const std::vector<double> zeros(variables.size(), 0.0);
        bool affine_gradient = true;

        for (std::size_t row = 0; row < gradient.size(); ++row) {
            double constant = 0.0;
            if (!evaluate_symbolic_at_point(gradient[row], variables, zeros, &constant)) {
                affine_gradient = false;
                break;
            }
            rhs[row] = -constant;

            for (std::size_t col = 0; col < variables.size(); ++col) {
                std::vector<double> sample = zeros;
                sample[col] = 1.0;
                double value = 0.0;
                if (!evaluate_symbolic_at_point(gradient[row],
                                                variables,
                                                sample,
                                                &value)) {
                    affine_gradient = false;
                    break;
                }
                coefficients[row][col] = value - constant;
            }
            if (!affine_gradient) {
                break;
            }

            std::vector<std::vector<double>> validation_samples;
            validation_samples.push_back(std::vector<double>(variables.size(), 1.0));
            for (std::size_t col = 0; col < variables.size(); ++col) {
                std::vector<double> sample = zeros;
                sample[col] = 2.0;
                validation_samples.push_back(sample);
            }
            for (const std::vector<double>& sample : validation_samples) {
                double actual = 0.0;
                if (!evaluate_symbolic_at_point(gradient[row],
                                                variables,
                                                sample,
                                                &actual)) {
                    affine_gradient = false;
                    break;
                }
                double predicted = constant;
                for (std::size_t col = 0; col < variables.size(); ++col) {
                    predicted += coefficients[row][col] * sample[col];
                }
                if (!mymath::is_near_zero(actual - predicted, 1e-8)) {
                    affine_gradient = false;
                    break;
                }
            }
            if (!affine_gradient) {
                break;
            }
        }

        if (affine_gradient) {
            std::vector<double> solution;
            if (!solve_linear_system(coefficients, rhs, &solution)) {
                *output = "No isolated critical point.";
                return true;
            }
            *output = format_critical_solution(solution);
            return true;
        }

        if (variables.size() > 3) {
            throw std::runtime_error(
                "critical nonlinear search supports up to 3 variables");
        }

        std::vector<std::vector<double>> starts = {std::vector<double>(variables.size(), 0.0)};
        const std::vector<double> seeds = {-2.0, -1.0, 1.0, 2.0};
        for (std::size_t dimension = 0; dimension < variables.size(); ++dimension) {
            std::vector<std::vector<double>> next = starts;
            for (const std::vector<double>& start : starts) {
                for (double seed : seeds) {
                    std::vector<double> candidate = start;
                    candidate[dimension] = seed;
                    next.push_back(candidate);
                }
            }
            starts.swap(next);
        }

        std::vector<std::vector<double>> solutions;
        for (std::vector<double> current : starts) {
            bool converged = false;
            for (int iteration = 0; iteration < 40; ++iteration) {
                std::vector<double> gradient_values(variables.size(), 0.0);
                double gradient_norm = 0.0;
                bool numeric_ok = true;
                for (std::size_t row = 0; row < variables.size(); ++row) {
                    if (!evaluate_symbolic_at_point(gradient[row],
                                                    variables,
                                                    current,
                                                    &gradient_values[row])) {
                        numeric_ok = false;
                        break;
                    }
                    gradient_norm += gradient_values[row] * gradient_values[row];
                }
                if (!numeric_ok) {
                    break;
                }
                if (gradient_norm < 1e-16) {
                    converged = true;
                    break;
                }

                std::vector<std::vector<double>> jacobian(
                    variables.size(), std::vector<double>(variables.size(), 0.0));
                for (std::size_t row = 0; row < variables.size() && numeric_ok; ++row) {
                    for (std::size_t col = 0; col < variables.size(); ++col) {
                        if (!evaluate_symbolic_at_point(hessian[row][col],
                                                        variables,
                                                        current,
                                                        &jacobian[row][col])) {
                            numeric_ok = false;
                            break;
                        }
                    }
                }
                if (!numeric_ok) {
                    break;
                }

                for (double& value : gradient_values) {
                    value = -value;
                }
                std::vector<double> step;
                if (!solve_linear_system(jacobian, gradient_values, &step)) {
                    break;
                }
                double step_norm = 0.0;
                for (std::size_t i = 0; i < current.size(); ++i) {
                    current[i] += step[i];
                    step_norm += step[i] * step[i];
                }
                if (step_norm < 1e-18) {
                    converged = true;
                    break;
                }
            }

            if (!converged) {
                continue;
            }

            bool duplicate = false;
            for (const std::vector<double>& existing : solutions) {
                double distance = 0.0;
                for (std::size_t i = 0; i < existing.size(); ++i) {
                    const double diff = existing[i] - current[i];
                    distance += diff * diff;
                }
                if (distance < 1e-10) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                solutions.push_back(current);
            }
        }

        if (solutions.empty()) {
            *output = "No isolated critical point.";
            return true;
        }
        std::sort(solutions.begin(), solutions.end());
        std::ostringstream out;
        if (solutions.size() == 1) {
            *output = format_critical_solution(solutions.front());
            return true;
        }
        out << "[";
        for (std::size_t i = 0; i < solutions.size(); ++i) {
            if (i != 0) {
                out << ", ";
            }
            out << format_critical_solution(solutions[i]);
        }
        out << "]";
        *output = out.str();
        return true;
    };

    auto handle_jacobian_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() < 2) {
            throw std::runtime_error(
                "jacobian expects a bracketed expression list and variable names");
        }

        const std::vector<SymbolicExpression> expressions =
            parse_symbolic_expression_list(arguments[0]);
        std::vector<std::string> fallback_variables;
        for (const SymbolicExpression& expression_item : expressions) {
            const std::vector<std::string> identifiers =
                expression_item.identifier_variables();
            fallback_variables.insert(fallback_variables.end(),
                                      identifiers.begin(),
                                      identifiers.end());
        }
        std::sort(fallback_variables.begin(), fallback_variables.end());
        fallback_variables.erase(std::unique(fallback_variables.begin(),
                                             fallback_variables.end()),
                                 fallback_variables.end());
        const std::vector<std::string> variables =
            parse_symbolic_variable_arguments(arguments, 1, fallback_variables);
        *output = symbolic_matrix_to_string(
            SymbolicExpression::jacobian(expressions, variables));
        return true;
    };

    auto handle_diff_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.empty()) {
            throw std::runtime_error(
                "diff expects a symbolic expression and optional variable names");
        }
        bool symbolic_derivative = arguments.size() == 1;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            symbolic_derivative =
                symbolic_derivative || is_identifier_text(trim_copy(arguments[i]));
        }
        if (symbolic_derivative) {
            std::string variable_name;
            SymbolicExpression expression;
            resolve_symbolic_expression(arguments[0],
                                        arguments.size() == 1,
                                        &variable_name,
                                        &expression);
            SymbolicExpression differentiated = expression;
            if (arguments.size() == 1) {
                differentiated = differentiated.derivative(variable_name).simplify();
            } else {
                for (std::size_t i = 1; i < arguments.size(); ++i) {
                    const std::string derivative_variable = trim_copy(arguments[i]);
                    if (!is_identifier_text(derivative_variable)) {
                        throw std::runtime_error(
                            "diff variable arguments must be identifiers");
                    }
                    differentiated =
                        differentiated.derivative(derivative_variable).simplify();
                }
            }
            *output = differentiated.simplify().to_string();
            return true;
        }
        if (arguments.size() != 2) {
            throw std::runtime_error(
                "diff numeric evaluation expects exactly two arguments");
        }
        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        DecimalParser parser(arguments[1], &impl_->variables, &impl_->functions);
        *output = format_decimal(normalize_result(analysis.derivative(parser.parse())));
        return true;
    };

    auto handle_limit_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 2 && arguments.size() != 3) {
            throw std::runtime_error(
                "limit expects 2 arguments for a two-sided limit or 3 with direction");
        }

        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        DecimalParser point_parser(arguments[1], &impl_->variables, &impl_->functions);
        int direction = 0;
        if (arguments.size() == 3) {
            DecimalParser direction_parser(arguments[2], &impl_->variables, &impl_->functions);
            const double direction_value = direction_parser.parse();
            if (!is_integer_double(direction_value)) {
                throw std::runtime_error("limit direction must be -1, 0, or 1");
            }
            direction = static_cast<int>(round_to_long_long(direction_value));
        }

        *output = format_decimal(normalize_result(
            analysis.limit(point_parser.parse(), direction)));
        return true;
    };

    auto handle_integral_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 1 && arguments.size() != 2 &&
            arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "integral expects 1 argument for symbolic indefinite integral, "
                "identifier arguments for chained symbolic integrals, "
                "2 arguments for indefinite value, "
                "3 for definite integral, or 4 for anchor and constant");
        }

        bool symbolic_integral = true;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            symbolic_integral =
                symbolic_integral && is_identifier_text(trim_copy(arguments[i]));
        }
        if (symbolic_integral) {
            std::string variable_name;
            SymbolicExpression expression;
            resolve_symbolic_expression(arguments[0],
                                        arguments.size() == 1,
                                        &variable_name,
                                        &expression);
            if (arguments.size() == 2) {
                variable_name = trim_copy(arguments[1]);
            }
            SymbolicExpression integrated = expression;
            if (arguments.size() == 1) {
                integrated = integrated.integral(variable_name).simplify();
            } else {
                for (std::size_t i = 1; i < arguments.size(); ++i) {
                    integrated =
                        integrated.integral(trim_copy(arguments[i])).simplify();
                }
            }
            *output = integrated.simplify().to_string() + " + C";
            return true;
        }
        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        if (arguments.size() == 2) {
            DecimalParser x_parser(arguments[1], &impl_->variables, &impl_->functions);
            *output = format_decimal(
                normalize_result(analysis.indefinite_integral_at(x_parser.parse())));
            return true;
        }
        if (arguments.size() == 3) {
            DecimalParser left_parser(arguments[1], &impl_->variables, &impl_->functions);
            DecimalParser right_parser(arguments[2], &impl_->variables, &impl_->functions);
            *output = format_decimal(normalize_result(
                analysis.definite_integral(left_parser.parse(), right_parser.parse())));
            return true;
        }

        DecimalParser x_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser anchor_parser(arguments[2], &impl_->variables, &impl_->functions);
        DecimalParser constant_parser(arguments[3], &impl_->variables, &impl_->functions);
        *output = format_decimal(normalize_result(
            analysis.indefinite_integral_at(x_parser.parse(),
                                            anchor_parser.parse(),
                                            constant_parser.parse())));
        return true;
    };

    auto handle_ode_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                command +
                " expects rhs, x0, y0, x1, optional steps, optional event, and optional params");
        }

        DecimalParser x0_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser y0_parser(arguments[2], &impl_->variables, &impl_->functions);
        DecimalParser x1_parser(arguments[3], &impl_->variables, &impl_->functions);
        int steps = command == "ode" ? 100 : 10;

        std::size_t optional_index = 4;
        int parsed_steps = steps;
        if (optional_index < arguments.size() &&
            try_parse_positive_step_argument(arguments[optional_index], &parsed_steps)) {
            steps = parsed_steps;
            ++optional_index;
        }

        std::string event_expression;
        bool has_event = false;
        StoredValue parameter_value;
        bool has_parameter = false;
        if (optional_index < arguments.size()) {
            if (optional_index + 1 == arguments.size()) {
                if (is_matrix_argument(arguments[optional_index])) {
                    parameter_value = evaluate_expression_value(this,
                                                               impl_.get(),
                                                               arguments[optional_index],
                                                               false);
                    has_parameter = true;
                } else {
                    event_expression = arguments[optional_index];
                    has_event = true;
                }
                ++optional_index;
            } else {
                event_expression = arguments[optional_index];
                has_event = true;
                ++optional_index;
                parameter_value =
                    evaluate_expression_value(this, impl_.get(), arguments[optional_index], false);
                has_parameter = true;
                ++optional_index;
            }
        }

        if (optional_index != arguments.size()) {
            throw std::runtime_error(
                command +
                " received too many optional arguments");
        }

        const auto evaluate_rhs = build_scoped_scalar_evaluator(arguments[0]);
        std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)> evaluate_event;
        if (has_event) {
            evaluate_event = build_scoped_scalar_evaluator(event_expression);
        }
        const ODESolver solver(
            [evaluate_rhs,
             has_parameter,
             parameter_value,
             &make_scalar_stored,
             &append_parameter_assignments](double x_value, double y_value) {
                std::vector<std::pair<std::string, StoredValue>> assignments;
                assignments.reserve(has_parameter ? 4 : 2);
                assignments.push_back({"x", make_scalar_stored(x_value)});
                assignments.push_back({"y", make_scalar_stored(y_value)});
                if (has_parameter) {
                    append_parameter_assignments(parameter_value, &assignments);
                }
                return evaluate_rhs(assignments);
            },
            has_event
                ? ODESolver::EventFunction(
                      [evaluate_event,
                       has_parameter,
                       parameter_value,
                       &make_scalar_stored,
                       &append_parameter_assignments](double x_value, double y_value) {
                          std::vector<std::pair<std::string, StoredValue>> assignments;
                          assignments.reserve(has_parameter ? 4 : 2);
                          assignments.push_back({"x", make_scalar_stored(x_value)});
                          assignments.push_back({"y", make_scalar_stored(y_value)});
                          if (has_parameter) {
                              append_parameter_assignments(parameter_value, &assignments);
                          }
                          return evaluate_event(assignments);
                      })
                : ODESolver::EventFunction());
        const double x0 = x0_parser.parse();
        const double y0 = y0_parser.parse();
        const double x1 = x1_parser.parse();

        if (command == "ode") {
            *output = format_decimal(normalize_result(
                solver.solve(x0, y0, x1, steps)));
            return true;
        }

        const std::vector<ODEPoint> points =
            solver.solve_trajectory(x0, y0, x1, steps);
        matrix::Matrix table(points.size(), 2, 0.0);
        for (std::size_t i = 0; i < points.size(); ++i) {
            table.at(i, 0) = normalize_result(points[i].x);
            table.at(i, 1) = normalize_result(points[i].y);
        }
        *output = matrix_literal_expression(table);
        return true;
    };

    auto handle_ode_system_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        const std::vector<std::string> arguments =
            split_top_level_arguments(inside);
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                command +
                " expects rhs_vector, x0, y0_vector, x1, optional steps, optional event, and optional params");
        }

        int steps = command == "ode_system" ? 100 : 10;
        std::size_t optional_index = 4;
        int parsed_steps = steps;
        if (optional_index < arguments.size() &&
            try_parse_positive_step_argument(arguments[optional_index], &parsed_steps)) {
            steps = parsed_steps;
            ++optional_index;
        }

        std::string event_expression;
        bool has_event = false;
        StoredValue parameter_value;
        bool has_parameter = false;
        if (optional_index < arguments.size()) {
            if (optional_index + 1 == arguments.size()) {
                if (is_matrix_argument(arguments[optional_index])) {
                    parameter_value = evaluate_expression_value(this,
                                                               impl_.get(),
                                                               arguments[optional_index],
                                                               false);
                    has_parameter = true;
                } else {
                    event_expression = arguments[optional_index];
                    has_event = true;
                }
                ++optional_index;
            } else {
                event_expression = arguments[optional_index];
                has_event = true;
                ++optional_index;
                parameter_value =
                    evaluate_expression_value(this, impl_.get(), arguments[optional_index], false);
                has_parameter = true;
                ++optional_index;
            }
        }

        if (optional_index != arguments.size()) {
            throw std::runtime_error(
                command +
                " received too many optional arguments");
        }

        const double x0 = parse_decimal_argument(arguments[1]);
        const double x1 = parse_decimal_argument(arguments[3]);
        const std::vector<double> initial_state =
            matrix_to_vector_values(parse_matrix_argument(arguments[2], command),
                                    command);
        const auto evaluate_rhs_matrix =
            build_scoped_matrix_evaluator(arguments[0]);
        std::function<double(const std::vector<std::pair<std::string, StoredValue>>&)> evaluate_event;
        if (has_event) {
            evaluate_event = build_scoped_scalar_evaluator(event_expression);
        }

        const ODESystemSolver solver(
            [evaluate_rhs_matrix,
             has_parameter,
             parameter_value,
             &make_scalar_stored,
             &vector_to_column_matrix,
             &append_parameter_assignments](double x_value,
                                            const std::vector<double>& y_value) {
                std::vector<std::pair<std::string, StoredValue>> assignments;
                assignments.reserve(y_value.size() + (has_parameter ? 4 : 2));

                assignments.push_back({"x", make_scalar_stored(x_value)});

                StoredValue y_matrix_stored;
                y_matrix_stored.is_matrix = true;
                y_matrix_stored.matrix = vector_to_column_matrix(y_value);
                assignments.push_back({"y", y_matrix_stored});

                for (std::size_t i = 0; i < y_value.size(); ++i) {
                    assignments.push_back({"y" + std::to_string(i + 1),
                                           make_scalar_stored(y_value[i])});
                }
                if (has_parameter) {
                    append_parameter_assignments(parameter_value, &assignments);
                }

                const matrix::Matrix rhs_matrix = evaluate_rhs_matrix(assignments);
                if (!rhs_matrix.is_vector()) {
                    throw std::runtime_error("ODE system right-hand side must evaluate to a vector");
                }
                const std::size_t result_size =
                    rhs_matrix.rows == 1 ? rhs_matrix.cols : rhs_matrix.rows;
                if (result_size != y_value.size()) {
                    throw std::runtime_error("ODE system right-hand side dimension mismatch");
                }

                std::vector<double> result(result_size, 0.0);
                for (std::size_t i = 0; i < result_size; ++i) {
                    result[i] = rhs_matrix.rows == 1 ? rhs_matrix.at(0, i) : rhs_matrix.at(i, 0);
                }
                return result;
            },
            has_event
                ? ODESystemSolver::EventFunction(
                      [evaluate_event,
                       has_parameter,
                       parameter_value,
                       &make_scalar_stored,
                       &vector_to_column_matrix,
                       &append_parameter_assignments](double x_value,
                                                      const std::vector<double>& y_value) {
                          std::vector<std::pair<std::string, StoredValue>> assignments;
                          assignments.reserve(y_value.size() + (has_parameter ? 4 : 2));
                          assignments.push_back({"x", make_scalar_stored(x_value)});

                          StoredValue y_matrix_stored;
                          y_matrix_stored.is_matrix = true;
                          y_matrix_stored.matrix = vector_to_column_matrix(y_value);
                          assignments.push_back({"y", y_matrix_stored});

                          for (std::size_t i = 0; i < y_value.size(); ++i) {
                              assignments.push_back({"y" + std::to_string(i + 1),
                                                     make_scalar_stored(y_value[i])});
                          }
                          if (has_parameter) {
                              append_parameter_assignments(parameter_value, &assignments);
                          }
                          return evaluate_event(assignments);
                      })
                : ODESystemSolver::EventFunction());

        if (command == "ode_system") {
            const std::vector<double> final_state =
                solver.solve(x0, initial_state, x1, steps);
            *output = matrix::Matrix::vector(final_state).to_string();
            return true;
        }

        const std::vector<ODESystemPoint> points =
            solver.solve_trajectory(x0, initial_state, x1, steps);
        matrix::Matrix table(points.size(), initial_state.size() + 1, 0.0);
        for (std::size_t row = 0; row < points.size(); ++row) {
            table.at(row, 0) = normalize_result(points[row].x);
            for (std::size_t col = 0; col < points[row].y.size(); ++col) {
                table.at(row, col + 1) = normalize_result(points[row].y[col]);
            }
        }
        *output = matrix_literal_expression(table);
        return true;
    };

    auto handle_decomposition_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 1 || !is_matrix_argument(arguments[0])) {
            throw std::runtime_error(command + " expects exactly one matrix argument");
        }

        const matrix::Matrix matrix_value = evaluate_expression_value(this,
                                                                      impl_.get(),
                                                                      arguments[0],
                                                                      false).matrix;
        if (command == "svd") {
            *output = "U: " + matrix::svd_u(matrix_value).to_string() +
                      "\nS: " + matrix::svd_s(matrix_value).to_string() +
                      "\nVt: " + matrix::svd_vt(matrix_value).to_string();
            return true;
        }

        try {
            *output = "values: " + matrix::eigenvalues(matrix_value).to_string() +
                      "\nvectors: " + matrix::eigenvectors(matrix_value).to_string();
            return true;
        } catch (const std::exception&) {
            if (matrix_value.rows == 2 && matrix_value.cols == 2) {
                const double trace = matrix_value.at(0, 0) + matrix_value.at(1, 1);
                const double det = matrix::determinant(matrix_value);
                const double discriminant = trace * trace - 4.0 * det;
                if (discriminant < 0.0) {
                    const double real = trace * 0.5;
                    const double imag = mymath::sqrt(-discriminant) * 0.5;
                    std::ostringstream out;
                    out << "values: [complex(" << format_decimal(real) << ", "
                        << format_decimal(imag) << "), complex("
                        << format_decimal(real) << ", " << format_decimal(-imag)
                        << ")]\nvectors: unavailable for complex eigenvalues";
                    *output = out.str();
                    return true;
                }
            }
            throw;
        }
    };

    auto handle_extrema_command = [&](const std::string& command,
                         const std::string& inside,
                         std::string* output) {
        (void)command;
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3) {
            throw std::runtime_error("extrema expects exactly three arguments");
        }

        const FunctionAnalysis analysis = build_analysis(arguments[0]);
        DecimalParser left_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser right_parser(arguments[2], &impl_->variables, &impl_->functions);
        const std::vector<ExtremumPoint> points =
            analysis.solve_extrema(left_parser.parse(), right_parser.parse());

        if (points.empty()) {
            *output = "No extrema found in the given interval.";
            return true;
        }

        std::ostringstream out;
        auto display_value = [](double value) {
            if (is_integer_double(value, 1e-6)) {
                return format_decimal(static_cast<double>(round_to_long_long(value)));
            }
            return format_decimal(normalize_result(value));
        };
        for (std::size_t i = 0; i < points.size(); ++i) {
            if (i != 0) {
                out << '\n';
            }
            out << (points[i].is_maximum ? "max" : "min")
                << ": x = " << display_value(points[i].x)
                << ", f(x) = " << display_value(points[i].value);
        }
        *output = out.str();
        return true;
    };

    const std::vector<FunctionCommandRegistration> command_registry = {
        { {"poly_add", "poly_sub", "poly_mul", "poly_div", "roots"},
          handle_polynomial_command },
        { {"taylor", "pade", "puiseux"},
          handle_series_module_command },
        { {"series_sum", "summation"},
          handle_series_sum_command },
        { {"fourier", "ifourier", "inverse_fourier", "laplace", "ilaplace",
           "inverse_laplace", "ztrans", "z_transform", "iztrans", "inverse_z"},
          handle_transform_command },
        { {"double_integral", "double_integral_cyl", "double_integral_polar",
           "triple_integral", "triple_integral_cyl", "triple_integral_sph"},
          handle_integration_command },
        { {"solve", "bisect", "secant", "fixed_point"},
          handle_rootfinding_command },
        { {"lp_max", "lp_min", "ilp_max", "ilp_min", "milp_max", "milp_min",
           "bip_max", "binary_max", "bip_min", "binary_min"},
          handle_optimization_command },
        { {"simplify"}, handle_simplify_command },
        { {"gradient"}, handle_gradient_command },
        { {"hessian"}, handle_hessian_command },
        { {"critical"}, handle_critical_command },
        { {"jacobian"}, handle_jacobian_command },
        { {"diff"}, handle_diff_command },
        { {"limit"}, handle_limit_command },
        { {"integral"}, handle_integral_command },
        { {"ode", "ode_table"}, handle_ode_command },
        { {"ode_system", "ode_system_table"}, handle_ode_system_command },
        { {"eig", "svd"}, handle_decomposition_command },
        { {"extrema"}, handle_extrema_command },
    };

    for (const FunctionCommandRegistration& registration : command_registry) {
        for (const std::string& command_name : registration.names) {
            std::string inside;
            if (!split_named_call(trimmed, command_name, &inside)) {
                continue;
            }
            return registration.handler(command_name, inside, output);
        }
    }

    return false;
}

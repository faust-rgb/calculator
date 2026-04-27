#include "calculator_internal_types.h"

#include "function_analysis.h"
#include "matrix.h"
#include "multivariable_integrator.h"
#include "mymath.h"
#include "ode_solver.h"
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
        if (lhs.size() != rhs.size()) {
            throw std::runtime_error("objective and solution dimension mismatch");
        }
        long double total = 0.0L;
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            total += static_cast<long double>(lhs[i]) *
                     static_cast<long double>(rhs[i]);
        }
        return static_cast<double>(total);
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

    std::string inside;
    if (split_named_call(trimmed, "poly_add", &inside) ||
        split_named_call(trimmed, "poly_sub", &inside) ||
        split_named_call(trimmed, "poly_mul", &inside) ||
        split_named_call(trimmed, "poly_div", &inside)) {
        // 这些命令只处理“已经定义好的多项式函数名”，
        // 例如 poly_mul(p, q)，而不是直接吃任意表达式。
        // 这样接口更清晰，也便于后续缓存或扩展更多多项式算法。
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 2) {
            throw std::runtime_error("polynomial operations expect exactly two arguments");
        }

        std::string lhs_variable;
        std::string rhs_variable;
        std::vector<double> lhs_coefficients;
        std::vector<double> rhs_coefficients;
        build_polynomial(arguments[0], &lhs_variable, &lhs_coefficients);
        build_polynomial(arguments[1], &rhs_variable, &rhs_coefficients);

        const std::string output_variable = lhs_variable;
        if (trimmed.rfind("poly_add", 0) == 0) {
            *output = polynomial_to_string(
                polynomial_add(lhs_coefficients, rhs_coefficients),
                output_variable);
            return true;
        }
        if (trimmed.rfind("poly_sub", 0) == 0) {
            *output = polynomial_to_string(
                polynomial_subtract(lhs_coefficients, rhs_coefficients),
                output_variable);
            return true;
        }
        if (trimmed.rfind("poly_mul", 0) == 0) {
            *output = polynomial_to_string(
                polynomial_multiply(lhs_coefficients, rhs_coefficients),
                output_variable);
            return true;
        }

        const PolynomialDivisionResult division =
            polynomial_divide(lhs_coefficients, rhs_coefficients);
        *output = "quotient: " +
                  polynomial_to_string(division.quotient, output_variable) +
                  ", remainder: " +
                  polynomial_to_string(division.remainder, output_variable);
        return true;
    }

    if (split_named_call(trimmed, "roots", &inside)) {
        // 求根目前限定为“多项式实根”。
        // 如果用户传进来的自定义函数不能化成多项式，会在 build_polynomial
        // 阶段直接报错，避免对非多项式做错误的根求解。
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 1) {
            throw std::runtime_error("roots expects exactly one argument");
        }

        std::string variable_name;
        std::vector<double> coefficients;
        build_polynomial(arguments[0], &variable_name, &coefficients);
        const std::vector<double> roots = polynomial_real_roots(coefficients);
        if (roots.empty()) {
            *output = "No real roots.";
            return true;
        }

        std::ostringstream out;
        for (std::size_t i = 0; i < roots.size(); ++i) {
            if (i != 0) {
                out << ", ";
            }
            const double root =
                is_integer_double(roots[i], 1e-6)
                    ? static_cast<double>(round_to_long_long(roots[i]))
                    : roots[i];
            out << format_symbolic_scalar(root);
        }
        *output = out.str();
        return true;
    }

    if (split_named_call(trimmed, "taylor", &inside)) {
        // Taylor 多项式按定义生成：
        //   f(a) + f'(a)(x-a) + f''(a)/2!(x-a)^2 + ...
        //
        // 这里采取一个务实实现：
        // - 导数本身走符号求导
        // - 每一阶导数在 a 处的值走现有数值求值
        //
        // 这样既能得到可读的展开式，也不需要再单独实现一套
        // “高阶导数闭式系数提取器”。
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3) {
            throw std::runtime_error("taylor expects exactly three arguments");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], true, &variable_name, &expression);
        DecimalParser center_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser degree_parser(arguments[2], &impl_->variables, &impl_->functions);
        const double center = center_parser.parse();
        const double degree_value = degree_parser.parse();
        if (!is_integer_double(degree_value) || degree_value < 0.0) {
            throw std::runtime_error("taylor degree must be a non-negative integer");
        }

        const int degree = static_cast<int>(round_to_long_long(degree_value));
        const std::vector<double> coefficients =
            build_taylor_coefficients(expression, variable_name, center, degree);
        *output = taylor_series_to_string(coefficients, variable_name, center);
        return true;
    }

    if (split_named_call(trimmed, "pade", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "pade expects expr, m, n or expr, center, m, n");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], true, &variable_name, &expression);

        const bool explicit_center = arguments.size() == 4;
        const double center = explicit_center
                                  ? parse_decimal_argument(arguments[1])
                                  : 0.0;
        const double numerator_degree_value = parse_decimal_argument(
            arguments[explicit_center ? 2 : 1]);
        const double denominator_degree_value = parse_decimal_argument(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(numerator_degree_value) ||
            numerator_degree_value < 0.0 ||
            !is_integer_double(denominator_degree_value) ||
            denominator_degree_value < 0.0) {
            throw std::runtime_error(
                "pade degrees must be non-negative integers");
        }

        const int numerator_degree =
            static_cast<int>(round_to_long_long(numerator_degree_value));
        const int denominator_degree =
            static_cast<int>(round_to_long_long(denominator_degree_value));
        if (numerator_degree == 0 && denominator_degree == 0) {
            throw std::runtime_error("pade requires at least one non-zero degree");
        }

        const std::vector<double> coefficients = build_taylor_coefficients(
            expression,
            variable_name,
            center,
            numerator_degree + denominator_degree);
        auto coefficient_at = [&](int index) {
            if (index < 0 ||
                index >= static_cast<int>(coefficients.size())) {
                return 0.0;
            }
            return coefficients[static_cast<std::size_t>(index)];
        };

        std::vector<double> denominator(denominator_degree + 1, 0.0);
        denominator[0] = 1.0;
        if (denominator_degree > 0) {
            std::vector<std::vector<double>> matrix(
                static_cast<std::size_t>(denominator_degree),
                std::vector<double>(static_cast<std::size_t>(denominator_degree), 0.0));
            std::vector<double> rhs(static_cast<std::size_t>(denominator_degree), 0.0);
            for (int row = 0; row < denominator_degree; ++row) {
                for (int col = 0; col < denominator_degree; ++col) {
                    matrix[static_cast<std::size_t>(row)]
                          [static_cast<std::size_t>(col)] =
                        coefficient_at(numerator_degree + row - col);
                }
                rhs[static_cast<std::size_t>(row)] =
                    -coefficient_at(numerator_degree + row + 1);
            }
            const std::vector<double> solved = solve_dense_linear_system(
                matrix, rhs, "pade");
            for (int i = 0; i < denominator_degree; ++i) {
                denominator[static_cast<std::size_t>(i + 1)] =
                    solved[static_cast<std::size_t>(i)];
            }
        }

        std::vector<double> numerator(numerator_degree + 1, 0.0);
        for (int i = 0; i <= numerator_degree; ++i) {
            double value = 0.0;
            for (int j = 0; j <= denominator_degree && j <= i; ++j) {
                value += denominator[static_cast<std::size_t>(j)] *
                         coefficient_at(i - j);
            }
            numerator[static_cast<std::size_t>(i)] = value;
        }

        const std::string base = shifted_series_base(variable_name, center);
        const std::string numerator_text =
            polynomial_to_string(numerator, base);
        const std::string denominator_text =
            polynomial_to_string(denominator, base);
        if (denominator_text == "1") {
            *output = simplify_symbolic_text(numerator_text);
        } else {
            *output = simplify_symbolic_text(
                "(" + numerator_text + ") / (" + denominator_text + ")");
        }
        return true;
    }

    if (split_named_call(trimmed, "puiseux", &inside)) {
        const std::vector<std::string> arguments = split_top_level_arguments(inside);
        if (arguments.size() != 3 && arguments.size() != 4) {
            throw std::runtime_error(
                "puiseux expects expr, degree, denominator or expr, center, degree, denominator");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], true, &variable_name, &expression);

        const bool explicit_center = arguments.size() == 4;
        const double center = explicit_center
                                  ? parse_decimal_argument(arguments[1])
                                  : 0.0;
        const double degree_value = parse_decimal_argument(
            arguments[explicit_center ? 2 : 1]);
        const double denominator_value = parse_decimal_argument(
            arguments[explicit_center ? 3 : 2]);
        if (!is_integer_double(degree_value) || degree_value < 0.0) {
            throw std::runtime_error(
                "puiseux degree must be a non-negative integer");
        }
        if (!is_integer_double(denominator_value) || denominator_value <= 0.0) {
            throw std::runtime_error(
                "puiseux denominator must be a positive integer");
        }

        const int degree = static_cast<int>(round_to_long_long(degree_value));
        const int denominator =
            static_cast<int>(round_to_long_long(denominator_value));
        const std::string auxiliary_variable = "puiseux_t";
        const std::string replacement_text =
            mymath::is_near_zero(center, 1e-10)
                ? auxiliary_variable + " ^ " + std::to_string(denominator)
                : format_symbolic_scalar(center) + " + " +
                      auxiliary_variable + " ^ " +
                      std::to_string(denominator);
        const SymbolicExpression substituted = expression.substitute(
            variable_name,
            SymbolicExpression::parse(replacement_text));
        const std::vector<double> coefficients = build_taylor_coefficients(
            substituted, auxiliary_variable, 0.0, degree);
        *output = generalized_series_to_string(
            coefficients, variable_name, center, denominator);
        return true;
    }

    if (split_named_call(trimmed, "series_sum", &inside) ||
        split_named_call(trimmed, "summation", &inside)) {
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
    }

    if (split_named_call(trimmed, "simplify", &inside)) {
        const std::string argument = trim_copy(inside);
        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(argument, false, &variable_name, &expression);
        *output = expression.simplify().to_string();
        return true;
    }

    std::string transform_inside;
    std::string transform_command;
    if (split_named_call(trimmed, "fourier", &transform_inside)) {
        transform_command = "fourier";
    } else if (split_named_call(trimmed, "ifourier", &transform_inside) ||
               split_named_call(trimmed, "inverse_fourier", &transform_inside)) {
        transform_command = "ifourier";
    } else if (split_named_call(trimmed, "laplace", &transform_inside)) {
        transform_command = "laplace";
    } else if (split_named_call(trimmed, "ilaplace", &transform_inside) ||
               split_named_call(trimmed, "inverse_laplace", &transform_inside)) {
        transform_command = "ilaplace";
    } else if (split_named_call(trimmed, "ztrans", &transform_inside) ||
               split_named_call(trimmed, "z_transform", &transform_inside)) {
        transform_command = "ztrans";
    } else if (split_named_call(trimmed, "iztrans", &transform_inside) ||
               split_named_call(trimmed, "inverse_z", &transform_inside)) {
        transform_command = "iztrans";
    }
    if (!transform_command.empty()) {
        const std::vector<std::string> arguments =
            split_top_level_arguments(transform_inside);
        if (arguments.size() != 1 && arguments.size() != 3) {
            throw std::runtime_error(
                transform_command +
                " expects either one symbolic expression or expression plus input/output variable names");
        }

        std::string variable_name;
        SymbolicExpression expression;
        resolve_symbolic_expression(arguments[0], false, &variable_name, &expression);

        std::string input_variable;
        std::string output_variable;
        if (arguments.size() == 3) {
            input_variable = trim_copy(arguments[1]);
            output_variable = trim_copy(arguments[2]);
            if (!is_identifier_text(input_variable) || !is_identifier_text(output_variable)) {
                throw std::runtime_error(transform_command + " variable names must be identifiers");
            }
        } else if (transform_command == "fourier") {
            input_variable = variable_name.empty() ? "t" : variable_name;
            output_variable = "w";
        } else if (transform_command == "ifourier") {
            input_variable = variable_name.empty() ? "w" : variable_name;
            output_variable = "t";
        } else if (transform_command == "laplace") {
            input_variable = variable_name.empty() ? "t" : variable_name;
            output_variable = "s";
        } else if (transform_command == "ilaplace") {
            input_variable = variable_name.empty() ? "s" : variable_name;
            output_variable = "t";
        } else if (transform_command == "ztrans") {
            input_variable = variable_name.empty() ? "n" : variable_name;
            output_variable = "z";
        } else {
            input_variable = variable_name.empty() ? "z" : variable_name;
            output_variable = "n";
        }

        if (transform_command == "fourier") {
            *output = expression.fourier_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else if (transform_command == "ifourier") {
            *output = expression.inverse_fourier_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else if (transform_command == "laplace") {
            *output = expression.laplace_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else if (transform_command == "ilaplace") {
            *output = expression.inverse_laplace_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else if (transform_command == "ztrans") {
            *output = expression.z_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        } else {
            *output = expression.inverse_z_transform(input_variable, output_variable)
                          .simplify()
                          .to_string();
        }
        return true;
    }

    if (split_named_call(trimmed, "gradient", &inside)) {
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
    }

    if (split_named_call(trimmed, "hessian", &inside)) {
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
    }

    if (split_named_call(trimmed, "critical", &inside)) {
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
    }

    if (split_named_call(trimmed, "jacobian", &inside)) {
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
    }

    if (split_named_call(trimmed, "diff", &inside)) {
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
    }

    if (split_named_call(trimmed, "limit", &inside)) {
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
    }

    if (split_named_call(trimmed, "integral", &inside)) {
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
    }

    std::string multivariable_inside;
    std::string multivariable_command;
    if (split_named_call(trimmed, "double_integral", &multivariable_inside)) {
        multivariable_command = "double_integral";
    } else if (split_named_call(trimmed, "double_integral_cyl", &multivariable_inside)) {
        multivariable_command = "double_integral_cyl";
    } else if (split_named_call(trimmed, "double_integral_polar", &multivariable_inside)) {
        multivariable_command = "double_integral_polar";
    } else if (split_named_call(trimmed, "triple_integral", &multivariable_inside)) {
        multivariable_command = "triple_integral";
    } else if (split_named_call(trimmed, "triple_integral_cyl", &multivariable_inside)) {
        multivariable_command = "triple_integral_cyl";
    } else if (split_named_call(trimmed, "triple_integral_sph", &multivariable_inside)) {
        multivariable_command = "triple_integral_sph";
    }
    if (!multivariable_command.empty()) {
        const std::vector<std::string> arguments =
            split_top_level_arguments(multivariable_inside);

        if (multivariable_command == "double_integral") {
            if (arguments.size() != 5 && arguments.size() != 7) {
                throw std::runtime_error(
                    "double_integral expects expr, x0, x1, y0, y1, and optional nx, ny");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            const std::vector<int> subdivisions =
                parse_subdivisions(arguments, 5, {32, 32});
            const MultivariableIntegrator integrator(
                [evaluate_expression](const std::vector<double>& point) {
                    return evaluate_expression(
                        {{"x", point[0]}, {"y", point[1]}});
                });
            *output = format_decimal(normalize_result(integrator.integrate(
                {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
                 {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])}},
                subdivisions)));
            return true;
        }

        if (multivariable_command == "double_integral_cyl" ||
            multivariable_command == "double_integral_polar") {
            if (arguments.size() != 5 && arguments.size() != 7) {
                throw std::runtime_error(
                    multivariable_command +
                    " expects expr, r0, r1, theta0, theta1, and optional nr, ntheta");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            const std::vector<int> subdivisions =
                parse_subdivisions(arguments, 5, {32, 32});
            const MultivariableIntegrator integrator(
                [evaluate_expression](const std::vector<double>& point) {
                    const double r = point[0];
                    const double theta = point[1];
                    const double x = r * mymath::cos(theta);
                    const double y = r * mymath::sin(theta);
                    return evaluate_expression(
                               {{"r", r},
                                {"theta", theta},
                                {"x", x},
                                {"y", y}}) *
                           r;
                });
            *output = format_decimal(normalize_result(integrator.integrate(
                {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
                 {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])}},
                subdivisions)));
            return true;
        }

        if (multivariable_command == "triple_integral") {
            if (arguments.size() != 7 && arguments.size() != 10) {
                throw std::runtime_error(
                    "triple_integral expects expr, x0, x1, y0, y1, z0, z1, and optional nx, ny, nz");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            const std::vector<int> subdivisions =
                parse_subdivisions(arguments, 7, {16, 16, 16});
            const MultivariableIntegrator integrator(
                [evaluate_expression](const std::vector<double>& point) {
                    return evaluate_expression(
                        {{"x", point[0]}, {"y", point[1]}, {"z", point[2]}});
                });
            *output = format_decimal(normalize_result(integrator.integrate(
                {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
                 {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])},
                 {parse_decimal_argument(arguments[5]), parse_decimal_argument(arguments[6])}},
                subdivisions)));
            return true;
        }

        if (multivariable_command == "triple_integral_cyl") {
            if (arguments.size() != 7 && arguments.size() != 10) {
                throw std::runtime_error(
                    "triple_integral_cyl expects expr, r0, r1, theta0, theta1, z0, z1, and optional nr, ntheta, nz");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            const std::vector<int> subdivisions =
                parse_subdivisions(arguments, 7, {16, 16, 16});
            const MultivariableIntegrator integrator(
                [evaluate_expression](const std::vector<double>& point) {
                    const double r = point[0];
                    const double theta = point[1];
                    const double z = point[2];
                    const double x = r * mymath::cos(theta);
                    const double y = r * mymath::sin(theta);
                    return evaluate_expression(
                               {{"r", r},
                                {"theta", theta},
                                {"z", z},
                                {"x", x},
                                {"y", y}}) *
                           r;
                });
            *output = format_decimal(normalize_result(integrator.integrate(
                {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
                 {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])},
                 {parse_decimal_argument(arguments[5]), parse_decimal_argument(arguments[6])}},
                subdivisions)));
            return true;
        }

        if (arguments.size() != 7 && arguments.size() != 10) {
            throw std::runtime_error(
                "triple_integral_sph expects expr, rho0, rho1, theta0, theta1, phi0, phi1, and optional nrho, ntheta, nphi");
        }
        const auto evaluate_expression =
            build_scoped_decimal_evaluator(arguments[0]);
        const std::vector<int> subdivisions =
            parse_subdivisions(arguments, 7, {16, 16, 16});
        const MultivariableIntegrator integrator(
            [evaluate_expression](const std::vector<double>& point) {
                const double rho = point[0];
                const double theta = point[1];
                const double phi = point[2];
                const double sin_phi = mymath::sin(phi);
                const double cos_phi = mymath::cos(phi);
                const double cos_theta = mymath::cos(theta);
                const double sin_theta = mymath::sin(theta);
                const double x = rho * sin_phi * cos_theta;
                const double y = rho * sin_phi * sin_theta;
                const double z = rho * cos_phi;
                return evaluate_expression(
                           {{"rho", rho},
                            {"theta", theta},
                            {"phi", phi},
                            {"x", x},
                            {"y", y},
                            {"z", z}}) *
                       rho * rho * sin_phi;
            });
        *output = format_decimal(normalize_result(integrator.integrate(
            {{parse_decimal_argument(arguments[1]), parse_decimal_argument(arguments[2])},
             {parse_decimal_argument(arguments[3]), parse_decimal_argument(arguments[4])},
             {parse_decimal_argument(arguments[5]), parse_decimal_argument(arguments[6])}},
            subdivisions)));
        return true;
    }

    std::string ode_inside;
    std::string ode_command_name;
    if (split_named_call(trimmed, "ode", &ode_inside)) {
        ode_command_name = "ode";
    } else if (split_named_call(trimmed, "ode_table", &ode_inside)) {
        ode_command_name = "ode_table";
    }
    if (!ode_command_name.empty()) {
        const std::vector<std::string> arguments = split_top_level_arguments(ode_inside);
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                ode_command_name +
                " expects rhs, x0, y0, x1, optional steps, optional event, and optional params");
        }

        DecimalParser x0_parser(arguments[1], &impl_->variables, &impl_->functions);
        DecimalParser y0_parser(arguments[2], &impl_->variables, &impl_->functions);
        DecimalParser x1_parser(arguments[3], &impl_->variables, &impl_->functions);
        int steps = ode_command_name == "ode" ? 100 : 10;

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
                ode_command_name +
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

        if (ode_command_name == "ode") {
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
    }

    std::string ode_system_inside;
    std::string ode_system_command_name;
    if (split_named_call(trimmed, "ode_system", &ode_system_inside)) {
        ode_system_command_name = "ode_system";
    } else if (split_named_call(trimmed, "ode_system_table", &ode_system_inside)) {
        ode_system_command_name = "ode_system_table";
    }
    if (!ode_system_command_name.empty()) {
        const std::vector<std::string> arguments =
            split_top_level_arguments(ode_system_inside);
        if (arguments.size() < 4 || arguments.size() > 7) {
            throw std::runtime_error(
                ode_system_command_name +
                " expects rhs_vector, x0, y0_vector, x1, optional steps, optional event, and optional params");
        }

        int steps = ode_system_command_name == "ode_system" ? 100 : 10;
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
                ode_system_command_name +
                " received too many optional arguments");
        }

        const double x0 = parse_decimal_argument(arguments[1]);
        const double x1 = parse_decimal_argument(arguments[3]);
        const std::vector<double> initial_state =
            matrix_to_vector_values(parse_matrix_argument(arguments[2], ode_system_command_name),
                                    ode_system_command_name);
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

        if (ode_system_command_name == "ode_system") {
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
    }

    std::string planning_inside;
    std::string planning_command;
    if (split_named_call(trimmed, "lp_max", &planning_inside)) {
        planning_command = "lp_max";
    } else if (split_named_call(trimmed, "lp_min", &planning_inside)) {
        planning_command = "lp_min";
    } else if (split_named_call(trimmed, "ilp_max", &planning_inside)) {
        planning_command = "ilp_max";
    } else if (split_named_call(trimmed, "ilp_min", &planning_inside)) {
        planning_command = "ilp_min";
    } else if (split_named_call(trimmed, "milp_max", &planning_inside)) {
        planning_command = "milp_max";
    } else if (split_named_call(trimmed, "milp_min", &planning_inside)) {
        planning_command = "milp_min";
    } else if (split_named_call(trimmed, "bip_max", &planning_inside) ||
               split_named_call(trimmed, "binary_max", &planning_inside)) {
        planning_command = "bip_max";
    } else if (split_named_call(trimmed, "bip_min", &planning_inside) ||
               split_named_call(trimmed, "binary_min", &planning_inside)) {
        planning_command = "bip_min";
    }
    if (!planning_command.empty()) {
        const std::vector<std::string> arguments =
            split_top_level_arguments(planning_inside);
        const std::vector<double> objective = matrix_to_vector_values(
            parse_matrix_argument(arguments[0], planning_command), planning_command);
        const std::size_t variable_count = objective.size();
        const double planning_tolerance = 1e-8;
        const bool is_binary_program =
            planning_command == "bip_max" || planning_command == "bip_min";
        const bool is_mixed_integer =
            planning_command == "milp_max" || planning_command == "milp_min";
        const bool is_integer_program =
            planning_command == "ilp_max" || planning_command == "ilp_min";
        const bool maximize =
            planning_command == "lp_max" || planning_command == "ilp_max" ||
            planning_command == "milp_max" || planning_command == "bip_max";

        std::size_t argument_index = 1;
        const matrix::Matrix inequality_matrix =
            parse_matrix_argument(arguments[argument_index++], planning_command);
        const std::vector<double> inequality_rhs = matrix_to_vector_values(
            parse_matrix_argument(arguments[argument_index++], planning_command), planning_command);

        matrix::Matrix equality_matrix(0, variable_count, 0.0);
        std::vector<double> equality_rhs;
        std::vector<double> lower_bounds(variable_count, 0.0);
        std::vector<double> upper_bounds(variable_count, 0.0);
        std::vector<double> integrality(variable_count, 0.0);

        if (is_binary_program) {
            if (arguments.size() != 3 && arguments.size() != 5) {
                throw std::runtime_error(
                    planning_command +
                    " expects objective_vector, A, b, and optional Aeq, beq");
            }
            if (arguments.size() == 5) {
                equality_matrix =
                    parse_matrix_argument(arguments[argument_index++], planning_command);
                equality_rhs = matrix_to_vector_values(
                    parse_matrix_argument(arguments[argument_index++], planning_command),
                    planning_command);
            }
            std::fill(lower_bounds.begin(), lower_bounds.end(), 0.0);
            std::fill(upper_bounds.begin(), upper_bounds.end(), 1.0);
            std::fill(integrality.begin(), integrality.end(), 1.0);
        } else if (is_mixed_integer) {
            if (arguments.size() != 6 && arguments.size() != 8) {
                throw std::runtime_error(
                    planning_command +
                    " expects objective_vector, A, b, lower_bounds, upper_bounds, integrality, and optional Aeq, beq");
            }
            if (arguments.size() == 8) {
                equality_matrix =
                    parse_matrix_argument(arguments[argument_index++], planning_command);
                equality_rhs = matrix_to_vector_values(
                    parse_matrix_argument(arguments[argument_index++], planning_command),
                    planning_command);
            }
            lower_bounds = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
            upper_bounds = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
            integrality = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
        } else {
            if (arguments.size() != 5 && arguments.size() != 7) {
                throw std::runtime_error(
                    planning_command +
                    " expects objective_vector, A, b, lower_bounds, upper_bounds, and optional Aeq, beq");
            }
            if (arguments.size() == 7) {
                equality_matrix =
                    parse_matrix_argument(arguments[argument_index++], planning_command);
                equality_rhs = matrix_to_vector_values(
                    parse_matrix_argument(arguments[argument_index++], planning_command),
                    planning_command);
            }
            lower_bounds = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
            upper_bounds = matrix_to_vector_values(
                parse_matrix_argument(arguments[argument_index++], planning_command),
                planning_command);
            if (is_integer_program) {
                std::fill(integrality.begin(), integrality.end(), 1.0);
            }
        }

        if (argument_index != arguments.size()) {
            throw std::runtime_error(planning_command + " received an invalid argument count");
        }

        if (inequality_matrix.cols != variable_count ||
            inequality_rhs.size() != inequality_matrix.rows ||
            equality_matrix.cols != variable_count ||
            equality_rhs.size() != equality_matrix.rows ||
            lower_bounds.size() != variable_count ||
            upper_bounds.size() != variable_count ||
            integrality.size() != variable_count) {
            throw std::runtime_error(planning_command + " dimension mismatch");
        }

        for (std::size_t i = 0; i < variable_count; ++i) {
            if (lower_bounds[i] > upper_bounds[i]) {
                throw std::runtime_error(planning_command + " requires lower_bounds <= upper_bounds");
            }
            if (integrality[i] != 0.0 && !is_integer_double(lower_bounds[i])) {
                throw std::runtime_error(planning_command + " requires integer lower bounds for integer variables");
            }
            if (integrality[i] != 0.0 && !is_integer_double(upper_bounds[i])) {
                throw std::runtime_error(planning_command + " requires integer upper bounds for integer variables");
            }
        }

        std::vector<double> transformed_objective = objective;
        if (!maximize) {
            for (double& value : transformed_objective) {
                value = -value;
            }
        }

        std::vector<std::size_t> integer_indices;
        std::vector<std::size_t> continuous_indices;
        for (std::size_t i = 0; i < variable_count; ++i) {
            if (mymath::is_near_zero(integrality[i], planning_tolerance)) {
                continuous_indices.push_back(i);
            } else {
                integer_indices.push_back(i);
            }
        }

        if (integer_indices.empty()) {
            std::vector<double> best_solution;
            double best_value = 0.0;
            std::string planning_diagnostic;
            if (!solve_linear_box_problem(transformed_objective,
                                          inequality_matrix,
                                          inequality_rhs,
                                          equality_matrix,
                                          equality_rhs,
                                          lower_bounds,
                                          upper_bounds,
                                          planning_tolerance,
                                          &best_solution,
                                          &best_value,
                                          &planning_diagnostic)) {
                std::string message =
                    planning_command + " found no feasible bounded solution";
                if (!planning_diagnostic.empty()) {
                    message += " (" + planning_diagnostic + ")";
                }
                throw std::runtime_error(message);
            }
            *output = format_planning_result(best_solution,
                                             dot_product(objective, best_solution));
            return true;
        }

        static constexpr std::size_t kMaxIntegerAssignments = 1000000;
        static constexpr std::size_t kMaxIntegerSearchNodes = 2000000;
        std::size_t estimated_integer_assignments = 1;
        std::vector<long long> integer_lower(variable_count, 0);
        std::vector<long long> integer_upper(variable_count, 0);
        for (std::size_t index : integer_indices) {
            integer_lower[index] = round_to_long_long(lower_bounds[index]);
            integer_upper[index] = round_to_long_long(upper_bounds[index]);
            const long double width =
                static_cast<long double>(integer_upper[index]) -
                static_cast<long double>(integer_lower[index]) + 1.0L;
            if (width <= 0.0L ||
                width > static_cast<long double>(kMaxIntegerAssignments) ||
                estimated_integer_assignments >
                    kMaxIntegerAssignments / static_cast<std::size_t>(width)) {
                estimated_integer_assignments = kMaxIntegerAssignments + 1;
            } else {
                estimated_integer_assignments *= static_cast<std::size_t>(width);
            }
        }
        if (estimated_integer_assignments > kMaxIntegerAssignments) {
            throw std::runtime_error(
                planning_command + " integer search limit exceeded: more than " +
                std::to_string(kMaxIntegerAssignments) +
                " possible integer assignments");
        }

        bool found = false;
        double best_value = 0.0;
        std::vector<double> best_solution(variable_count, 0.0);
        std::vector<long long> current_integer_values(variable_count, 0);
        std::size_t visited_integer_nodes = 0;

        std::function<void(std::size_t, long double)> search_integer =
            [&](std::size_t depth, long double current_objective) {
                ++visited_integer_nodes;
                if (visited_integer_nodes > kMaxIntegerSearchNodes) {
                    throw std::runtime_error(
                        planning_command + " integer search node limit exceeded after " +
                        std::to_string(kMaxIntegerSearchNodes) + " nodes");
                }

                for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                    long double assigned_total = 0.0L;
                    for (std::size_t assigned_depth = 0; assigned_depth < depth; ++assigned_depth) {
                        const std::size_t col = integer_indices[assigned_depth];
                        assigned_total += static_cast<long double>(inequality_matrix.at(row, col)) *
                                          static_cast<long double>(current_integer_values[col]);
                    }

                    long double minimum_possible = assigned_total;
                    for (std::size_t remaining_depth = depth;
                         remaining_depth < integer_indices.size();
                         ++remaining_depth) {
                        const std::size_t col = integer_indices[remaining_depth];
                        const long double coefficient =
                            static_cast<long double>(inequality_matrix.at(row, col));
                        minimum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(integer_lower[col])
                                : coefficient * static_cast<long double>(integer_upper[col]);
                    }
                    for (std::size_t col : continuous_indices) {
                        const long double coefficient =
                            static_cast<long double>(inequality_matrix.at(row, col));
                        minimum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(lower_bounds[col])
                                : coefficient * static_cast<long double>(upper_bounds[col]);
                    }
                    if (minimum_possible >
                        static_cast<long double>(inequality_rhs[row]) + planning_tolerance) {
                        return;
                    }
                }

                for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                    long double assigned_total = 0.0L;
                    for (std::size_t assigned_depth = 0; assigned_depth < depth; ++assigned_depth) {
                        const std::size_t col = integer_indices[assigned_depth];
                        assigned_total += static_cast<long double>(equality_matrix.at(row, col)) *
                                          static_cast<long double>(current_integer_values[col]);
                    }

                    long double minimum_possible = assigned_total;
                    long double maximum_possible = assigned_total;
                    for (std::size_t remaining_depth = depth;
                         remaining_depth < integer_indices.size();
                         ++remaining_depth) {
                        const std::size_t col = integer_indices[remaining_depth];
                        const long double coefficient =
                            static_cast<long double>(equality_matrix.at(row, col));
                        minimum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(integer_lower[col])
                                : coefficient * static_cast<long double>(integer_upper[col]);
                        maximum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(integer_upper[col])
                                : coefficient * static_cast<long double>(integer_lower[col]);
                    }
                    for (std::size_t col : continuous_indices) {
                        const long double coefficient =
                            static_cast<long double>(equality_matrix.at(row, col));
                        minimum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(lower_bounds[col])
                                : coefficient * static_cast<long double>(upper_bounds[col]);
                        maximum_possible +=
                            coefficient >= 0.0L
                                ? coefficient * static_cast<long double>(upper_bounds[col])
                                : coefficient * static_cast<long double>(lower_bounds[col]);
                    }

                    const long double target = static_cast<long double>(equality_rhs[row]);
                    if (target < minimum_possible - planning_tolerance ||
                        target > maximum_possible + planning_tolerance) {
                        return;
                    }
                }

                long double optimistic_objective = current_objective;
                for (std::size_t remaining_depth = depth;
                     remaining_depth < integer_indices.size();
                     ++remaining_depth) {
                    const std::size_t col = integer_indices[remaining_depth];
                    const long double coefficient =
                        static_cast<long double>(transformed_objective[col]);
                    optimistic_objective +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(integer_upper[col])
                            : coefficient * static_cast<long double>(integer_lower[col]);
                }
                for (std::size_t col : continuous_indices) {
                    const long double coefficient =
                        static_cast<long double>(transformed_objective[col]);
                    optimistic_objective +=
                        coefficient >= 0.0L
                            ? coefficient * static_cast<long double>(upper_bounds[col])
                            : coefficient * static_cast<long double>(lower_bounds[col]);
                }
                if (found &&
                    optimistic_objective <=
                        static_cast<long double>(best_value) + planning_tolerance) {
                    return;
                }

                if (depth == integer_indices.size()) {
                    std::vector<double> candidate(variable_count, 0.0);
                    for (std::size_t col = 0; col < variable_count; ++col) {
                        candidate[col] = lower_bounds[col];
                    }
                    for (std::size_t col : integer_indices) {
                        candidate[col] = static_cast<double>(current_integer_values[col]);
                    }

                    if (continuous_indices.empty()) {
                        bool feasible = true;
                        for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                            long double total = 0.0L;
                            for (std::size_t col = 0; col < variable_count; ++col) {
                                total += static_cast<long double>(inequality_matrix.at(row, col)) *
                                         static_cast<long double>(candidate[col]);
                            }
                            if (total >
                                static_cast<long double>(inequality_rhs[row]) + planning_tolerance) {
                                feasible = false;
                                break;
                            }
                        }
                        if (!feasible) {
                            return;
                        }
                        for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                            long double total = 0.0L;
                            for (std::size_t col = 0; col < variable_count; ++col) {
                                total += static_cast<long double>(equality_matrix.at(row, col)) *
                                         static_cast<long double>(candidate[col]);
                            }
                            if (mymath::abs(static_cast<double>(total - equality_rhs[row])) >
                                planning_tolerance) {
                                feasible = false;
                                break;
                            }
                        }
                        if (!feasible) {
                            return;
                        }

                        const double objective_value = dot_product(transformed_objective, candidate);
                        if (!found || objective_value > best_value + planning_tolerance) {
                            found = true;
                            best_value = objective_value;
                            best_solution = candidate;
                        }
                        return;
                    }

                    matrix::Matrix reduced_inequality(inequality_matrix.rows,
                                                      continuous_indices.size(),
                                                      0.0);
                    std::vector<double> reduced_inequality_rhs(inequality_rhs.size(), 0.0);
                    for (std::size_t row = 0; row < inequality_matrix.rows; ++row) {
                        long double rhs_adjustment = static_cast<long double>(inequality_rhs[row]);
                        for (std::size_t col : integer_indices) {
                            rhs_adjustment -=
                                static_cast<long double>(inequality_matrix.at(row, col)) *
                                static_cast<long double>(candidate[col]);
                        }
                        reduced_inequality_rhs[row] = static_cast<double>(rhs_adjustment);
                        for (std::size_t reduced_col = 0;
                             reduced_col < continuous_indices.size();
                             ++reduced_col) {
                            reduced_inequality.at(row, reduced_col) =
                                inequality_matrix.at(row, continuous_indices[reduced_col]);
                        }
                    }

                    matrix::Matrix reduced_equality(equality_matrix.rows,
                                                    continuous_indices.size(),
                                                    0.0);
                    std::vector<double> reduced_equality_rhs(equality_rhs.size(), 0.0);
                    for (std::size_t row = 0; row < equality_matrix.rows; ++row) {
                        long double rhs_adjustment = static_cast<long double>(equality_rhs[row]);
                        for (std::size_t col : integer_indices) {
                            rhs_adjustment -=
                                static_cast<long double>(equality_matrix.at(row, col)) *
                                static_cast<long double>(candidate[col]);
                        }
                        reduced_equality_rhs[row] = static_cast<double>(rhs_adjustment);
                        for (std::size_t reduced_col = 0;
                             reduced_col < continuous_indices.size();
                             ++reduced_col) {
                            reduced_equality.at(row, reduced_col) =
                                equality_matrix.at(row, continuous_indices[reduced_col]);
                        }
                    }

                    std::vector<double> reduced_objective(continuous_indices.size(), 0.0);
                    std::vector<double> reduced_lower(continuous_indices.size(), 0.0);
                    std::vector<double> reduced_upper(continuous_indices.size(), 0.0);
                    for (std::size_t reduced_col = 0;
                         reduced_col < continuous_indices.size();
                         ++reduced_col) {
                        const std::size_t original_col = continuous_indices[reduced_col];
                        reduced_objective[reduced_col] = transformed_objective[original_col];
                        reduced_lower[reduced_col] = lower_bounds[original_col];
                        reduced_upper[reduced_col] = upper_bounds[original_col];
                    }

                    std::vector<double> reduced_solution;
                    double reduced_objective_value = 0.0;
                    std::string reduced_diagnostic;
                    if (!solve_linear_box_problem(reduced_objective,
                                                  reduced_inequality,
                                                  reduced_inequality_rhs,
                                                  reduced_equality,
                                                  reduced_equality_rhs,
                                                  reduced_lower,
                                                  reduced_upper,
                                                  planning_tolerance,
                                                  &reduced_solution,
                                                  &reduced_objective_value,
                                                  &reduced_diagnostic)) {
                        return;
                    }

                    for (std::size_t reduced_col = 0;
                         reduced_col < continuous_indices.size();
                         ++reduced_col) {
                        candidate[continuous_indices[reduced_col]] = reduced_solution[reduced_col];
                    }

                    const double objective_value = dot_product(transformed_objective, candidate);
                    if (!found || objective_value > best_value + planning_tolerance) {
                        found = true;
                        best_value = objective_value;
                        best_solution = candidate;
                    }
                    return;
                }

                const std::size_t current_col = integer_indices[depth];
                const bool descending = transformed_objective[current_col] >= 0.0;
                if (descending) {
                    for (long long value = integer_upper[current_col];
                         value >= integer_lower[current_col];
                         --value) {
                        current_integer_values[current_col] = value;
                        search_integer(
                            depth + 1,
                            current_objective +
                                static_cast<long double>(transformed_objective[current_col]) *
                                    static_cast<long double>(value));
                    }
                } else {
                    for (long long value = integer_lower[current_col];
                         value <= integer_upper[current_col];
                         ++value) {
                        current_integer_values[current_col] = value;
                        search_integer(
                            depth + 1,
                            current_objective +
                                static_cast<long double>(transformed_objective[current_col]) *
                                    static_cast<long double>(value));
                    }
                }
            };

        search_integer(0, 0.0L);
        if (!found) {
            throw std::runtime_error(planning_command + " found no feasible mixed-integer solution");
        }

        *output = format_planning_result(best_solution,
                                         dot_product(objective, best_solution));
        return true;
    }

    std::string root_inside;
    std::string root_command_name;
    if (split_named_call(trimmed, "solve", &root_inside)) {
        root_command_name = "solve";
    } else if (split_named_call(trimmed, "bisect", &root_inside)) {
        root_command_name = "bisect";
    } else if (split_named_call(trimmed, "secant", &root_inside)) {
        root_command_name = "secant";
    } else if (split_named_call(trimmed, "fixed_point", &root_inside)) {
        root_command_name = "fixed_point";
    }
    if (!root_command_name.empty()) {
        const std::vector<std::string> arguments = split_top_level_arguments(root_inside);
        if (root_command_name == "solve") {
            if (arguments.size() == 2 &&
                !is_matrix_argument(arguments[0]) &&
                !is_matrix_argument(arguments[1])) {
                const auto evaluate_expression =
                    build_scoped_decimal_evaluator(arguments[0]);
                double x = parse_decimal_argument(arguments[1]);
                for (int iteration = 0; iteration < 64; ++iteration) {
                    const double fx = evaluate_expression({{"x", x}});
                    if (mymath::abs(fx) <= root_function_tolerance(fx)) {
                        *output = format_decimal(normalize_result(x));
                        return true;
                    }
                    const double h = root_derivative_step(x);
                    const long double derivative =
                        (static_cast<long double>(evaluate_expression({{"x", x + h}})) -
                         static_cast<long double>(evaluate_expression({{"x", x - h}}))) /
                        (2.0L * static_cast<long double>(h));
                    if (mymath::abs_long_double(derivative) <=
                        1e-12L * std::max(1.0L, mymath::abs_long_double(static_cast<long double>(fx)))) {
                        throw std::runtime_error("solve failed because the derivative vanished");
                    }
                    const double next = static_cast<double>(
                        static_cast<long double>(x) -
                        static_cast<long double>(fx) / derivative);
                    if (mymath::abs(next - x) <=
                        root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
                        *output = format_decimal(normalize_result(next));
                        return true;
                    }
                    x = next;
                }
                *output = format_decimal(normalize_result(x));
                return true;
            }
        } else if (root_command_name == "bisect") {
            if (arguments.size() != 3 || is_matrix_argument(arguments[0])) {
                throw std::runtime_error("bisect expects expression, a, b");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            double left = parse_decimal_argument(arguments[1]);
            double right = parse_decimal_argument(arguments[2]);
            if (left > right) {
                std::swap(left, right);
            }
            double left_value = evaluate_expression({{"x", left}});
            double right_value = evaluate_expression({{"x", right}});
            if (left_value * right_value > 0.0) {
                throw std::runtime_error("bisect requires f(a) and f(b) to have opposite signs");
            }
            for (int iteration = 0; iteration < 100; ++iteration) {
                const double mid = 0.5 * (left + right);
                const double mid_value = evaluate_expression({{"x", mid}});
                if (mymath::abs(mid_value) <= root_function_tolerance(mid_value) ||
                    mymath::abs(right - left) <=
                        root_position_tolerance(std::max(mymath::abs(left), mymath::abs(right)))) {
                    *output = format_decimal(normalize_result(mid));
                    return true;
                }
                if ((left_value < 0.0 && mid_value > 0.0) ||
                    (left_value > 0.0 && mid_value < 0.0)) {
                    right = mid;
                    right_value = mid_value;
                } else {
                    left = mid;
                    left_value = mid_value;
                }
            }
            *output = format_decimal(normalize_result(0.5 * (left + right)));
            return true;
        } else if (root_command_name == "secant") {
            if (arguments.size() != 3 || is_matrix_argument(arguments[0])) {
                throw std::runtime_error("secant expects expression, x0, x1");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            double x0 = parse_decimal_argument(arguments[1]);
            double x1 = parse_decimal_argument(arguments[2]);
            for (int iteration = 0; iteration < 64; ++iteration) {
                const double f0 = evaluate_expression({{"x", x0}});
                const double f1 = evaluate_expression({{"x", x1}});
                const long double denominator =
                    static_cast<long double>(f1) - static_cast<long double>(f0);
                if (mymath::abs_long_double(denominator) <=
                    1e-12L * std::max({1.0L,
                                       mymath::abs_long_double(static_cast<long double>(f0)),
                                       mymath::abs_long_double(static_cast<long double>(f1))})) {
                    throw std::runtime_error("secant failed because consecutive function values matched");
                }
                const double next = static_cast<double>(
                    static_cast<long double>(x1) -
                    static_cast<long double>(f1) *
                        (static_cast<long double>(x1) - static_cast<long double>(x0)) /
                        denominator);
                if (mymath::abs(next - x1) <=
                    root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x1)))) {
                    *output = format_decimal(normalize_result(next));
                    return true;
                }
                x0 = x1;
                x1 = next;
            }
            *output = format_decimal(normalize_result(x1));
            return true;
        } else if (root_command_name == "fixed_point") {
            if (arguments.size() != 2 || is_matrix_argument(arguments[0])) {
                throw std::runtime_error("fixed_point expects expression, x0");
            }
            const auto evaluate_expression =
                build_scoped_decimal_evaluator(arguments[0]);
            double x = parse_decimal_argument(arguments[1]);
            for (int iteration = 0; iteration < 128; ++iteration) {
                const double next = evaluate_expression({{"x", x}});
                if (mymath::abs(next - x) <=
                    root_position_tolerance(std::max(mymath::abs(next), mymath::abs(x)))) {
                    *output = format_decimal(normalize_result(next));
                    return true;
                }
                x = next;
            }
            *output = format_decimal(normalize_result(x));
            return true;
        }
    }

    std::string decomposition_inside;
    std::string decomposition_command;
    if (split_named_call(trimmed, "eig", &decomposition_inside)) {
        decomposition_command = "eig";
    } else if (split_named_call(trimmed, "svd", &decomposition_inside)) {
        decomposition_command = "svd";
    }
    if (!decomposition_command.empty()) {
        const std::vector<std::string> arguments = split_top_level_arguments(decomposition_inside);
        if (arguments.size() != 1 || !is_matrix_argument(arguments[0])) {
            throw std::runtime_error(decomposition_command + " expects exactly one matrix argument");
        }

        const matrix::Matrix matrix_value = evaluate_expression_value(this,
                                                                      impl_.get(),
                                                                      arguments[0],
                                                                      false).matrix;
        if (decomposition_command == "svd") {
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
    }

    if (split_named_call(trimmed, "extrema", &inside)) {
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
    }

    return false;
}

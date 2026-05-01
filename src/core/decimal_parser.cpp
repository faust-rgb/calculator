#include "calculator_internal_types.h"
#include "base_parser.h"
#include "matrix.h"
#include "mymath.h"
#include "statistics/calculator_statistics.h"
#include <algorithm>
#include <map>

class DecimalParserImpl : public BaseParser {
public:
    using ScalarFunction = std::function<double(const std::vector<double>&)>;

    DecimalParserImpl(std::string source,
                      const VariableResolver& variables,
                      const std::map<std::string, CustomFunction>* functions,
                      const std::map<std::string, ScalarFunction>* scalar_functions = nullptr,
                      HasScriptFunctionCallback has_script_function = {},
                      InvokeScriptFunctionDecimalCallback invoke_script_function = {})
        : BaseParser(std::move(source)),
          variables_(variables),
          functions_(functions),
          scalar_functions_(scalar_functions),
          has_script_function_(std::move(has_script_function)),
          invoke_script_function_(std::move(invoke_script_function)) {}

    double parse() {
        double value = parse_comparison();
        skip_spaces();
        if (!is_at_end()) {
            throw SyntaxError("unexpected token near: " + source_.substr(pos_, 1));
        }
        return value;
    }

private:
    double parse_comparison() {
        double value = parse_expression();
        while (true) {
            skip_spaces();
            if (match_string("==")) {
                value = mymath::is_near_zero(value - parse_expression(), 1e-10) ? 1.0 : 0.0;
            } else if (match_string("!=")) {
                value = mymath::is_near_zero(value - parse_expression(), 1e-10) ? 0.0 : 1.0;
            } else if (match_string("<=")) {
                value = value <= parse_expression() ? 1.0 : 0.0;
            } else if (match_string(">=")) {
                value = value >= parse_expression() ? 1.0 : 0.0;
            } else if (match('<')) {
                value = value < parse_expression() ? 1.0 : 0.0;
            } else if (match('>')) {
                value = value > parse_expression() ? 1.0 : 0.0;
            } else {
                break;
            }
        }
        return value;
    }

    double parse_expression() {
        double value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value += parse_term();
            } else if (match('-')) {
                value -= parse_term();
            } else {
                break;
            }
        }
        return value;
    }

    double parse_term() {
        double value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value *= parse_unary();
            } else if (match('/')) {
                const double divisor = parse_unary();
                if (divisor == 0.0) {
                    throw MathError("division by zero");
                }
                value /= divisor;
            } else {
                break;
            }
        }
        return value;
    }

    double parse_power() {
        double value = parse_primary();
        skip_spaces();
        if (match('^')) {
            const double exponent = parse_unary();
            return mymath::pow(value, exponent);
        }
        return value;
    }

    double parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            return -parse_unary();
        }
        return parse_power();
    }

    double parse_primary() {
        skip_spaces();
        if (match('(')) {
            const double value = parse_expression();
            skip_spaces();
            expect(')');
            return value;
        }

        if (peek_is_alpha()) {
            const std::string name = parse_identifier();

            skip_spaces();
            if (!peek('(')) {
                return lookup_variable(name);
            }

            skip_spaces();
            expect('(');
            const std::vector<double> arguments = parse_argument_list();
            expect(')');
            return apply_function(name, arguments);
        }

        return parse_number();
    }

    std::vector<double> parse_argument_list() {
        std::vector<double> arguments;
        skip_spaces();
        if (peek(')')) {
            return arguments;
        }

        while (true) {
            arguments.push_back(parse_expression());
            skip_spaces();
            if (!match(',')) {
                break;
            }
        }
        return arguments;
    }

    double parse_number() {
        skip_spaces();

        if (!is_at_end() &&
            source_[pos_] == '0' &&
            pos_ + 1 < source_.size()) {
            int base = 10;
            if (prefixed_base(source_[pos_ + 1], &base)) {
                const std::size_t start = pos_;
                pos_ += 2;
                while (!is_at_end()) {
                    const int digit = digit_value(source_[pos_]);
                    if (digit < 0 || digit >= base) {
                        break;
                    }
                    ++pos_;
                }
                return static_cast<double>(
                    parse_prefixed_integer_token(source_.substr(start, pos_ - start)));
            }
        }

        const std::size_t start = pos_;
        bool has_digit = false;
        bool seen_dot = false;

        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isdigit(static_cast<unsigned char>(ch))) {
                has_digit = true;
                ++pos_;
            } else if (ch == '.' && !seen_dot) {
                seen_dot = true;
                ++pos_;
            } else {
                break;
            }
        }

        if (!is_at_end() && (source_[pos_] == 'e' || source_[pos_] == 'E')) {
            const std::size_t exponent_pos = pos_;
            ++pos_;
            if (!is_at_end() && (source_[pos_] == '+' || source_[pos_] == '-')) {
                ++pos_;
            }
            const std::size_t exponent_digits = pos_;
            while (!is_at_end() &&
                   std::isdigit(static_cast<unsigned char>(source_[pos_]))) {
                ++pos_;
            }
            if (exponent_digits == pos_) {
                pos_ = exponent_pos;
            }
        }

        if (!has_digit) {
            throw SyntaxError("expected number");
        }

        return std::stod(source_.substr(start, pos_ - start));
    }

    double lookup_variable(const std::string& name) const {
        const StoredValue* found = variables_.lookup(name);
        if (found) {
            if (found->is_matrix || found->is_complex ||
                found->is_string || found->has_symbolic_text) {
                throw MathError("unsupported variable type in numeric expression: " + name);
            }
            return found->exact ? rational_to_double(found->rational)
                                   : found->decimal;
        }

        throw UndefinedError("unknown variable: " + name);
    }

    double apply_function(const std::string& name, const std::vector<double>& arguments) {
        const auto require_integer_argument =
            [&name](double value, const std::string& label) -> long long {
                if (!is_integer_double(value)) {
                    throw MathError(name + " requires integer " + label);
                }
                return round_to_long_long(value);
            };

        const auto it = functions_->find(name);
        if (it != functions_->end()) {
            if (arguments.size() != 1) {
                throw MathError("custom function " + name + " expects 1 argument");
            }
            std::map<std::string, StoredValue> snapshot = variables_.snapshot();
            StoredValue arg_value;
            arg_value.decimal = arguments[0];
            snapshot[it->second.parameter_name] = arg_value;
            DecimalParser parser(it->second.expression,
                                 VariableResolver(&snapshot, nullptr),
                                 functions_,
                                 scalar_functions_,
                                 has_script_function_,
                                 invoke_script_function_);
            return parser.parse();
        }

        if (has_script_function_ && has_script_function_(name)) {
            return invoke_script_function_(name, arguments);
        }

        if (scalar_functions_) {
            const auto it_ext = scalar_functions_->find(name);
            if (it_ext != scalar_functions_->end()) {
                return it_ext->second(arguments);
            }
        }

        throw UndefinedError("unknown function: " + name);
    }

    VariableResolver variables_;
    const std::map<std::string, CustomFunction>* functions_;
    const std::map<std::string, ScalarFunction>* scalar_functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionDecimalCallback invoke_script_function_;
};

double parse_decimal_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, DecimalParser::ScalarFunction>* scalar_functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionDecimalCallback invoke_script_function) {
    DecimalParserImpl parser(expression,
                             variables,
                             functions,
                             scalar_functions,
                             std::move(has_script_function),
                             std::move(invoke_script_function));
    return parser.parse();
}

DecimalParser::DecimalParser(
    std::string source,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    const std::map<std::string, ScalarFunction>* scalar_functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionDecimalCallback invoke_script_function)
    : source_(std::move(source)),
      variables_(variables),
      functions_(functions),
      scalar_functions_(scalar_functions),
      has_script_function_(std::move(has_script_function)),
      invoke_script_function_(std::move(invoke_script_function)) {}

double DecimalParser::parse() {
    DecimalParserImpl parser(source_,
                             variables_,
                             functions_,
                             scalar_functions_,
                             std::move(has_script_function_),
                             std::move(invoke_script_function_));
    return parser.parse();
}

bool try_evaluate_matrix_expression(const std::string& expression,
                                    const VariableResolver& variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const std::map<std::string, DecimalParser::ScalarFunction>* scalar_functions,
                                    const std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>>* matrix_functions,
                                    const std::map<std::string, matrix::ValueFunction>* value_functions,
                                    const HasScriptFunctionCallback& has_script_function,
                                    const InvokeScriptFunctionDecimalCallback& invoke_script_function,
                                    matrix::Value* value) {
    const matrix::ScalarEvaluator scalar_evaluator =
        [variables, functions, scalar_functions, has_script_function, invoke_script_function](const std::string& text) {
            DecimalParserImpl parser(text,
                                     variables,
                                     functions,
                                     scalar_functions,
                                     has_script_function,
                                     invoke_script_function);
            const double scalar_value = parser.parse();
            return mymath::is_near_zero(scalar_value, 1e-10) ? 0.0 : scalar_value;
        };
    const matrix::MatrixLookup matrix_lookup =
        [variables](const std::string& name, matrix::Matrix* matrix_value) {
            const StoredValue* found = variables.lookup(name);
            if (!found || !found->is_matrix) {
                return false;
            }
            *matrix_value = found->matrix;
            return true;
        };
    const matrix::ComplexLookup complex_lookup =
        [variables](const std::string& name, matrix::ComplexNumber* complex_value) {
            const StoredValue* found = variables.lookup(name);
            if (!found || !found->is_complex) {
                return false;
            }
            *complex_value = found->complex;
            return true;
        };
    return matrix::try_evaluate_expression(expression,
                                           scalar_evaluator,
                                           matrix_lookup,
                                           complex_lookup,
                                           matrix_functions,
                                           value_functions,
                                           value);
}

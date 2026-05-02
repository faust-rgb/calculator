#include "symbolic_render_parser.h"
#include "calculator_internal_types.h"
#include "parser/base_parser.h"
#include "math/helpers/base_conversions.h"
#include "mymath.h"
#include "symbolic_expression.h"
#include <algorithm>
#include <cctype>
#include <sstream>

class SymbolicRenderParserImpl : public BaseParser {
public:
    SymbolicRenderParserImpl(std::string_view source,
                             const VariableResolver& variables,
                             const std::map<std::string, CustomFunction>* functions,
                             int depth = 0)
        : BaseParser(source),
          variables_(variables),
          functions_(functions),
          depth_(depth) {}

    bool parse(std::string* output, bool* used_symbolic_constant) {
        if (depth_ > 12) {
            return false;
        }
        try {
            used_symbolic_constant_ = false;
            const std::string text = parse_expression();
            skip_spaces();
            if (pos_ != source_.size()) {
                return false;
            }
            SymbolicExpression expression = SymbolicExpression::parse(text);
            *output = expression.to_string();
            *used_symbolic_constant = used_symbolic_constant_;
            return true;
        } catch (const std::exception&) {
            return false;
        }
    }

private:
    std::string parse_expression() {
        std::string value = parse_term();
        while (true) {
            skip_spaces();
            if (match('+')) {
                value = "(" + value + " + " + parse_term() + ")";
            } else if (match('-')) {
                value = "(" + value + " - " + parse_term() + ")";
            } else {
                return value;
            }
        }
    }

    std::string parse_term() {
        std::string value = parse_unary();
        while (true) {
            skip_spaces();
            if (match('*')) {
                value = "(" + value + " * " + parse_unary() + ")";
            } else if (match('/')) {
                value = "(" + value + " / " + parse_unary() + ")";
            } else {
                return value;
            }
        }
    }

    std::string parse_power() {
        std::string value = parse_primary();
        skip_spaces();
        if (match('^')) {
            value = "(" + value + " ^ " + parse_unary() + ")";
        }
        return value;
    }

    std::string parse_unary() {
        skip_spaces();
        if (match('+')) {
            return parse_unary();
        }
        if (match('-')) {
            return "(-" + parse_unary() + ")";
        }
        return parse_power();
    }

    std::string parse_primary() {
        skip_spaces();
        if (match('(')) {
            const std::string value = parse_expression();
            skip_spaces();
            expect(')');
            return "(" + value + ")";
        }
        if (peek_is_alpha()) {
            const std::string name(parse_identifier());
            skip_spaces();
            if (peek() != '(') {
                return render_identifier(name);
            }

            expect('(');
            const std::vector<std::string> arguments = parse_argument_list();
            expect(')');
            return render_function(name, arguments);
        }
        return parse_number_token();
    }

    std::vector<std::string> parse_argument_list() {
        std::vector<std::string> arguments;
        skip_spaces();
        if (peek() == ')') {
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

    std::string render_identifier(const std::string& name) {
        if (name == "pi" || name == "e") {
            used_symbolic_constant_ = true;
            return name;
        }

        double builtin_constant = 0.0;
        if (lookup_builtin_constant(name, &builtin_constant)) {
            return format_symbolic_scalar(builtin_constant);
        }

        const StoredValue* found = variables_.lookup(name);
        if (!found) {
            return name;
        }
        if (found->is_matrix || found->is_complex || found->is_string) {
            throw std::runtime_error("unsupported symbolic variable");
        }
        if (found->has_symbolic_text) {
            used_symbolic_constant_ = true;
        }
        if (found->has_symbolic_text) return "(" + found->symbolic_text + ")";
        if (found->exact) return "(" + found->rational.to_string() + ")";
        if (found->has_precise_decimal_text) return "(" + found->precise_decimal_text + ")";
        return "(" + format_decimal(normalize_display_decimal(found->decimal)) + ")";
    }

    std::string render_function(const std::string& name,
                                const std::vector<std::string>& arguments) {
        if (name == "pow") {
            if (arguments.size() != 2) {
                throw std::runtime_error("pow expects two arguments");
            }
            return "((" + arguments[0] + ") ^ (" + arguments[1] + "))";
        }
        if (name == "root") {
            if (arguments.size() != 2) {
                throw std::runtime_error("root expects two arguments");
            }
            return "((" + arguments[0] + ") ^ (1 / (" + arguments[1] + ")))";
        }
        const auto function_it = functions_->find(name);
        if (function_it != functions_->end()) {
            if (arguments.size() != function_it->second.parameter_names.size()) {
                throw std::runtime_error("custom function " + name + " expects " +
                                         std::to_string(function_it->second.parameter_names.size()) +
                                         " arguments, but got " + std::to_string(arguments.size()));
            }
            std::map<std::string, StoredValue> scoped_variables = variables_.snapshot();
            for (std::size_t i = 0; i < arguments.size(); ++i) {
                StoredValue parameter_value;
                parameter_value.has_symbolic_text = true;
                parameter_value.symbolic_text = arguments[i];
                scoped_variables[function_it->second.parameter_names[i]] = parameter_value;
            }

            SymbolicRenderParserImpl nested(function_it->second.expression,
                                            VariableResolver(&scoped_variables, nullptr),
                                            functions_,
                                            depth_ + 1);
            std::string expanded;
            bool nested_symbolic = false;
            if (!nested.parse(&expanded, &nested_symbolic)) {
                throw std::runtime_error("unable to expand custom function symbolically");
            }
            used_symbolic_constant_ = used_symbolic_constant_ || nested_symbolic;
            return "(" + expanded + ")";
        }

        static const auto is_supported_symbolic_unary_function = [](const std::string& n) {
            return n == "sin" || n == "cos" || n == "tan" ||
                   n == "asin" || n == "acos" || n == "atan" ||
                   n == "exp" || n == "ln" || n == "log10" ||
                   n == "sqrt" || n == "abs" || n == "sign" ||
                   n == "floor" || n == "ceil" || n == "cbrt" ||
                   n == "step" || n == "delta";
        };

        if (!is_supported_symbolic_unary_function(name) || arguments.size() != 1) {
            throw std::runtime_error("unsupported symbolic function");
        }
        return name + "(" + arguments[0] + ")";
    }

    std::string parse_number_token() {
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
                return format_decimal(static_cast<double>(
                    parse_prefixed_integer_token(std::string(source_.substr(start, pos_ - start)))));
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
            throw std::runtime_error("expected number");
        }
        return format_decimal(std::stod(std::string(source_.substr(start, pos_ - start))));
    }

    VariableResolver variables_;
    const std::map<std::string, CustomFunction>* functions_;
    int depth_ = 0;
    bool used_symbolic_constant_ = false;
};

bool try_symbolic_constant_expression(const std::string& expression,
                                      const VariableResolver& variables,
                                      const std::map<std::string, CustomFunction>* functions,
                                      std::string* output) {
    SymbolicRenderParserImpl parser(expression, variables, functions);
    bool used_symbolic_constant = false;
    if (!parser.parse(output, &used_symbolic_constant)) {
        return false;
    }
    return used_symbolic_constant;
}

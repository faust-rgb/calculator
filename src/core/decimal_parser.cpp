#include "calculator_internal_types.h"
#include "base_parser.h"
#include "matrix.h"
#include "mymath.h"
#include "statistics/calculator_statistics.h"
#include <algorithm>
#include <map>

class DecimalParserImpl : public BaseParser {
public:
    DecimalParserImpl(std::string source,
                      const VariableResolver& variables,
                      const std::map<std::string, CustomFunction>* functions,
                      HasScriptFunctionCallback has_script_function = {},
                      InvokeScriptFunctionDecimalCallback invoke_script_function = {})
        : BaseParser(std::move(source)),
          variables_(variables),
          functions_(functions),
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
                                 has_script_function_,
                                 invoke_script_function_);
            return parser.parse();
        }

        if (has_script_function_ && has_script_function_(name)) {
            return invoke_script_function_(name, arguments);
        }

        if (name == "sin") {
            if (arguments.size() != 1) throw MathError("sin expects 1 argument");
            return mymath::sin(arguments[0]);
        }
        if (name == "cos") {
            if (arguments.size() != 1) throw MathError("cos expects 1 argument");
            return mymath::cos(arguments[0]);
        }
        if (name == "tan") {
            if (arguments.size() != 1) throw MathError("tan expects 1 argument");
            return mymath::tan(arguments[0]);
        }
        if (name == "asin") {
            if (arguments.size() != 1) throw MathError("asin expects 1 argument");
            return mymath::asin(arguments[0]);
        }
        if (name == "acos") {
            if (arguments.size() != 1) throw MathError("acos expects 1 argument");
            return mymath::acos(arguments[0]);
        }
        if (name == "atan") {
            if (arguments.size() == 1) return mymath::atan(arguments[0]);
            if (arguments.size() == 2) return mymath::atan2(arguments[0], arguments[1]);
            throw MathError("atan expects 1 or 2 arguments");
        }
        if (name == "sec") {
            if (arguments.size() != 1) throw MathError("sec expects 1 argument");
            return mymath::sec(arguments[0]);
        }
        if (name == "csc") {
            if (arguments.size() != 1) throw MathError("csc expects 1 argument");
            return mymath::csc(arguments[0]);
        }
        if (name == "cot") {
            if (arguments.size() != 1) throw MathError("cot expects 1 argument");
            return mymath::cot(arguments[0]);
        }
        if (name == "asec") {
            if (arguments.size() != 1) throw MathError("asec expects 1 argument");
            return mymath::asec(arguments[0]);
        }
        if (name == "acsc") {
            if (arguments.size() != 1) throw MathError("acsc expects 1 argument");
            return mymath::acsc(arguments[0]);
        }
        if (name == "acot") {
            if (arguments.size() != 1) throw MathError("acot expects 1 argument");
            return mymath::acot(arguments[0]);
        }
        if (name == "sinh") {
            if (arguments.size() != 1) throw MathError("sinh expects 1 argument");
            return mymath::sinh(arguments[0]);
        }
        if (name == "cosh") {
            if (arguments.size() != 1) throw MathError("cosh expects 1 argument");
            return mymath::cosh(arguments[0]);
        }
        if (name == "tanh") {
            if (arguments.size() != 1) throw MathError("tanh expects 1 argument");
            return mymath::tanh(arguments[0]);
        }
        if (name == "asinh") {
            if (arguments.size() != 1) throw MathError("asinh expects 1 argument");
            return mymath::asinh(arguments[0]);
        }
        if (name == "acosh") {
            if (arguments.size() != 1) throw MathError("acosh expects 1 argument");
            return mymath::acosh(arguments[0]);
        }
        if (name == "atanh") {
            if (arguments.size() != 1) throw MathError("atanh expects 1 argument");
            return mymath::atanh(arguments[0]);
        }
        if (name == "exp") {
            if (arguments.size() != 1) throw MathError("exp expects 1 argument");
            return mymath::exp(arguments[0]);
        }
        if (name == "exp2") {
            if (arguments.size() != 1) throw MathError("exp2 expects 1 argument");
            return mymath::pow(2.0, arguments[0]);
        }
        if (name == "ln") {
            if (arguments.size() != 1) throw MathError("ln expects 1 argument");
            return mymath::ln(arguments[0]);
        }
        if (name == "log") {
            if (arguments.size() == 1) return mymath::ln(arguments[0]);
            if (arguments.size() == 2) {
                if (arguments[1] <= 0.0 || mymath::is_near_zero(arguments[1] - 1.0)) {
                    throw MathError("log base must be positive and not equal to 1");
                }
                return mymath::ln(arguments[0]) / mymath::ln(arguments[1]);
            }
            throw MathError("log expects 1 or 2 arguments");
        }
        if (name == "log2") {
            if (arguments.size() != 1) throw MathError("log2 expects 1 argument");
            return mymath::ln(arguments[0]) / mymath::ln(2.0);
        }
        if (name == "log10") {
            if (arguments.size() != 1) throw MathError("log10 expects 1 argument");
            return mymath::ln(arguments[0]) / mymath::ln(10.0);
        }
        if (name == "sqrt") {
            if (arguments.size() != 1) throw MathError("sqrt expects 1 argument");
            return mymath::sqrt(arguments[0]);
        }
        if (name == "cbrt") {
            if (arguments.size() != 1) throw MathError("cbrt expects 1 argument");
            return mymath::cbrt(arguments[0]);
        }
        if (name == "root") {
            if (arguments.size() != 2) throw MathError("root expects 2 arguments");
            return mymath::root(arguments[0], arguments[1]);
        }
        if (name == "pow") {
            if (arguments.size() != 2) throw MathError("pow expects 2 arguments");
            return mymath::pow(arguments[0], arguments[1]);
        }
        if (name == "abs") {
            if (arguments.size() != 1) throw MathError("abs expects 1 argument");
            return mymath::abs(arguments[0]);
        }
        if (name == "sign") {
            if (arguments.size() != 1) throw MathError("sign expects 1 argument");
            if (mymath::is_near_zero(arguments[0], 1e-12)) {
                return 0.0;
            }
            return arguments[0] > 0.0 ? 1.0 : -1.0;
        }
        if (name == "floor") {
            if (arguments.size() != 1) throw MathError("floor expects 1 argument");
            return static_cast<double>(floor_to_long_long(arguments[0]));
        }
        if (name == "ceil") {
            if (arguments.size() != 1) throw MathError("ceil expects 1 argument");
            return static_cast<double>(ceil_to_long_long(arguments[0]));
        }
        if (name == "round") {
            if (arguments.size() != 1) throw MathError("round expects 1 argument");
            return static_cast<double>(round_to_long_long(arguments[0]));
        }
        if (name == "trunc") {
            if (arguments.size() != 1) throw MathError("trunc expects 1 argument");
            return static_cast<double>(trunc_to_long_long(arguments[0]));
        }
        if (name == "min") {
            if (arguments.empty()) throw MathError("min expects at least 1 argument");
            double res = arguments[0];
            for (std::size_t i = 1; i < arguments.size(); ++i) res = std::min(res, arguments[i]);
            return res;
        }
        if (name == "max") {
            if (arguments.empty()) throw MathError("max expects at least 1 argument");
            double res = arguments[0];
            for (std::size_t i = 1; i < arguments.size(); ++i) res = std::max(res, arguments[i]);
            return res;
        }
        if (name == "clamp") {
            if (arguments.size() != 3) throw MathError("clamp expects 3 arguments");
            const double low = std::min(arguments[1], arguments[2]);
            const double high = std::max(arguments[1], arguments[2]);
            return std::clamp(arguments[0], low, high);
        }
        if (name == "gamma") {
            if (arguments.size() != 1) throw MathError("gamma expects 1 argument");
            return mymath::gamma(arguments[0]);
        }
        if (name == "beta") {
            if (arguments.size() != 2) throw MathError("beta expects 2 arguments");
            return mymath::beta(arguments[0], arguments[1]);
        }
        if (name == "zeta") {
            if (arguments.size() != 1) throw MathError("zeta expects 1 argument");
            return mymath::zeta(arguments[0]);
        }
        if (name == "erf") {
            if (arguments.size() != 1) throw MathError("erf expects 1 argument");
            return mymath::erf(arguments[0]);
        }
        if (name == "erfc") {
            if (arguments.size() != 1) throw MathError("erfc expects 1 argument");
            return mymath::erfc(arguments[0]);
        }
        if (name == "bessel" || name == "bessel_j") {
            if (arguments.size() != 2) throw MathError("bessel expects 2 arguments");
            return mymath::bessel_j(
                static_cast<int>(require_integer_argument(arguments[0], "order")),
                arguments[1]);
        }
        if (name == "gcd") {
            if (arguments.size() != 2) throw MathError("gcd expects 2 arguments");
            return static_cast<double>(gcd_ll(round_to_long_long(arguments[0]), round_to_long_long(arguments[1])));
        }
        if (name == "lcm") {
            if (arguments.size() != 2) throw MathError("lcm expects 2 arguments");
            return static_cast<double>(lcm_ll(round_to_long_long(arguments[0]), round_to_long_long(arguments[1])));
        }
        if (name == "mod") {
            if (arguments.size() != 2) throw MathError("mod expects 2 arguments");
            const long long lhs = require_integer_argument(arguments[0], "lhs");
            const long long rhs = require_integer_argument(arguments[1], "rhs");
            if (rhs == 0) {
                throw MathError("mod by zero");
            }
            return static_cast<double>(lhs % rhs);
        }
        if (name == "factorial") {
            if (arguments.size() != 1) throw MathError("factorial expects 1 argument");
            return factorial_value(require_integer_argument(arguments[0], "argument"));
        }
        if (name == "nCr" || name == "binom") {
            if (arguments.size() != 2) throw MathError("combination expects 2 arguments");
            return combination_value(require_integer_argument(arguments[0], "n"),
                                     require_integer_argument(arguments[1], "r"));
        }
        if (name == "nPr") {
            if (arguments.size() != 2) throw MathError("permutation expects 2 arguments");
            return permutation_value(require_integer_argument(arguments[0], "n"),
                                     require_integer_argument(arguments[1], "r"));
        }
        if (name == "fib") {
            if (arguments.size() != 1) throw MathError("fib expects 1 argument");
            return fibonacci_value(round_to_long_long(arguments[0]));
        }
        if (name == "is_prime") {
            if (arguments.size() != 1) throw MathError("is_prime expects 1 argument");
            return is_prime_ll(require_integer_argument(arguments[0], "argument")) ? 1.0 : 0.0;
        }
        if (name == "next_prime") {
            if (arguments.size() != 1) throw MathError("next_prime expects 1 argument");
            return static_cast<double>(next_prime_ll(require_integer_argument(arguments[0], "argument")));
        }
        if (name == "prev_prime") {
            if (arguments.size() != 1) throw MathError("prev_prime expects 1 argument");
            return static_cast<double>(prev_prime_ll(require_integer_argument(arguments[0], "argument")));
        }
        if (name == "euler_phi") {
            if (arguments.size() != 1) throw MathError("euler_phi expects 1 argument");
            return static_cast<double>(euler_phi_ll(require_integer_argument(arguments[0], "argument")));
        }
        if (name == "phi") {
            if (arguments.size() != 1) throw MathError("phi expects 1 argument");
            return static_cast<double>(euler_phi_ll(require_integer_argument(arguments[0], "argument")));
        }
        if (name == "mobius") {
            if (arguments.size() != 1) throw MathError("mobius expects 1 argument");
            return static_cast<double>(mobius_ll(require_integer_argument(arguments[0], "argument")));
        }
        if (name == "prime_pi") {
            if (arguments.size() != 1) throw MathError("prime_pi expects 1 argument");
            return static_cast<double>(prime_pi_ll(require_integer_argument(arguments[0], "argument")));
        }
        if (name == "egcd") {
            if (arguments.size() != 2) throw MathError("egcd expects 2 arguments");
            long long x = 0;
            long long y = 0;
            return static_cast<double>(extended_gcd_ll(
                require_integer_argument(arguments[0], "a"),
                require_integer_argument(arguments[1], "b"),
                &x,
                &y));
        }
        if (name == "rand") {
            if (!arguments.empty()) throw MathError("rand expects no arguments");
            return stats_ops::apply_probability(name, arguments);
        }
        if (name == "randn") {
            if (!arguments.empty()) throw MathError("randn expects no arguments");
            return stats_ops::apply_probability(name, arguments);
        }
        if (name == "randint") {
            if (arguments.size() != 2) throw MathError("randint expects 2 arguments");
            return stats_ops::apply_probability(name, arguments);
        }
        if (name == "pdf_normal") {
            if (arguments.size() != 3) throw MathError("pdf_normal expects 3 arguments");
            return stats_ops::apply_probability("pdf_normal", arguments);
        }
        if (name == "cdf_normal") {
            if (arguments.size() != 3) throw MathError("cdf_normal expects 3 arguments");
            return stats_ops::apply_probability("cdf_normal", arguments);
        }
        if (name == "deg") {
            if (arguments.size() != 1) throw MathError("deg expects 1 argument");
            return radians_to_degrees(arguments[0]);
        }
        if (name == "rad") {
            if (arguments.size() != 1) throw MathError("rad expects 1 argument");
            return degrees_to_radians(arguments[0]);
        }
        if (name == "deg2rad") {
            if (arguments.size() != 1) throw MathError("deg2rad expects 1 argument");
            return degrees_to_radians(arguments[0]);
        }
        if (name == "rad2deg") {
            if (arguments.size() != 1) throw MathError("rad2deg expects 1 argument");
            return radians_to_degrees(arguments[0]);
        }
        if (name == "sin_deg") {
            if (arguments.size() != 1) throw MathError("sin_deg expects 1 argument");
            return mymath::sin(degrees_to_radians(arguments[0]));
        }
        if (name == "cos_deg") {
            if (arguments.size() != 1) throw MathError("cos_deg expects 1 argument");
            return mymath::cos(degrees_to_radians(arguments[0]));
        }
        if (name == "celsius") {
            if (arguments.size() != 1) throw MathError("celsius expects 1 argument");
            return fahrenheit_to_celsius(arguments[0]);
        }
        if (name == "fahrenheit") {
            if (arguments.size() != 1) throw MathError("fahrenheit expects 1 argument");
            return celsius_to_fahrenheit(arguments[0]);
        }
        if (name == "kelvin") {
            if (arguments.size() != 1) throw MathError("kelvin expects 1 argument");
            return arguments[0] + 273.15;
        }
        if (name == "c2f") {
            if (arguments.size() != 1) throw MathError("c2f expects 1 argument");
            return celsius_to_fahrenheit(arguments[0]);
        }
        if (name == "f2c") {
            if (arguments.size() != 1) throw MathError("f2c expects 1 argument");
            return fahrenheit_to_celsius(arguments[0]);
        }
        if (name == "sum") {
            if (arguments.empty()) throw MathError("sum expects at least 1 argument");
            long double sum = 0.0L;
            long double compensation = 0.0L;
            for (double arg : arguments) {
                const long double adjusted = static_cast<long double>(arg) - compensation;
                const long double next = sum + adjusted;
                compensation = (next - sum) - adjusted;
                sum = next;
            }
            return static_cast<double>(sum);
        }
        if (name == "mean" || name == "avg") {
            return stats::mean(arguments);
        }
        if (name == "median") {
            return stats::median(arguments);
        }
        if (name == "mode") {
            return stats::mode(arguments);
        }
        if (name == "percentile") {
            if (arguments.size() < 2) throw MathError("percentile expects percentage and data");
            std::vector<double> data(arguments.begin() + 1, arguments.end());
            return stats::percentile(data, arguments[0]);
        }
        if (name == "quartile") {
            if (arguments.size() < 2) throw MathError("quartile expects q and data");
            if (!is_integer_double(arguments[0])) {
                throw MathError("quartile q must be an integer");
            }
            std::vector<double> data(arguments.begin() + 1, arguments.end());
            return stats::quartile(data, static_cast<int>(arguments[0]));
        }
        if (name == "var") {
            return stats::variance(arguments);
        }
        if (name == "std") {
            return stats::stddev(arguments);
        }
        if (name == "sample_var") {
            return stats::sample_variance(arguments);
        }
        if (name == "sample_std") {
            return stats::sample_stddev(arguments);
        }
        if (name == "skewness") {
            return stats::skewness(arguments);
        }
        if (name == "kurtosis") {
            return stats::kurtosis(arguments);
        }
        if (name == "cov" || name == "corr" || name == "slope" || name == "intercept") {
            return stats_ops::apply_statistic(name, arguments);
        }
        if (name == "bernoulli") {
            return stats_ops::apply_probability(name, arguments);
        }
        if (name == "and") {
            if (arguments.size() != 2) throw MathError("and expects 2 arguments");
            return static_cast<double>(require_integer_argument(arguments[0], "lhs") &
                                       require_integer_argument(arguments[1], "rhs"));
        }
        if (name == "or") {
            if (arguments.size() != 2) throw MathError("or expects 2 arguments");
            return static_cast<double>(require_integer_argument(arguments[0], "lhs") |
                                       require_integer_argument(arguments[1], "rhs"));
        }
        if (name == "xor") {
            if (arguments.size() != 2) throw MathError("xor expects 2 arguments");
            return static_cast<double>(require_integer_argument(arguments[0], "lhs") ^
                                       require_integer_argument(arguments[1], "rhs"));
        }
        if (name == "not") {
            if (arguments.size() != 1) throw MathError("not expects 1 argument");
            return static_cast<double>(~require_integer_argument(arguments[0], "argument"));
        }
        if (name == "shl") {
            if (arguments.size() != 2) throw MathError("shl expects 2 arguments");
            const long long count = require_integer_argument(arguments[1], "shift");
            if (count < 0) {
                throw MathError("shift count cannot be negative");
            }
            return static_cast<double>(require_integer_argument(arguments[0], "value") << count);
        }
        if (name == "shr") {
            if (arguments.size() != 2) throw MathError("shr expects 2 arguments");
            const long long count = require_integer_argument(arguments[1], "shift");
            if (count < 0) {
                throw MathError("shift count cannot be negative");
            }
            return static_cast<double>(require_integer_argument(arguments[0], "value") >> count);
        }
        if (name == "rol") {
            if (arguments.size() != 2) throw MathError("rol expects 2 arguments");
            return static_cast<double>(from_unsigned_bits(rotate_left_bits(
                to_unsigned_bits(require_integer_argument(arguments[0], "value")),
                normalize_rotation_count(require_integer_argument(arguments[1], "shift")))));
        }
        if (name == "ror") {
            if (arguments.size() != 2) throw MathError("ror expects 2 arguments");
            return static_cast<double>(from_unsigned_bits(rotate_right_bits(
                to_unsigned_bits(require_integer_argument(arguments[0], "value")),
                normalize_rotation_count(require_integer_argument(arguments[1], "shift")))));
        }
        if (name == "popcount") {
            if (arguments.size() != 1) throw MathError("popcount expects 1 argument");
            return static_cast<double>(popcount_bits(
                to_unsigned_bits(require_integer_argument(arguments[0], "argument"))));
        }
        if (name == "bitlen") {
            if (arguments.size() != 1) throw MathError("bitlen expects 1 argument");
            return static_cast<double>(bit_length_bits(
                to_unsigned_bits(require_integer_argument(arguments[0], "argument"))));
        }
        if (name == "ctz") {
            if (arguments.size() != 1) throw MathError("ctz expects 1 argument");
            return static_cast<double>(trailing_zero_count_bits(
                to_unsigned_bits(require_integer_argument(arguments[0], "argument"))));
        }
        if (name == "clz") {
            if (arguments.size() != 1) throw MathError("clz expects 1 argument");
            return static_cast<double>(leading_zero_count_bits(
                to_unsigned_bits(require_integer_argument(arguments[0], "argument"))));
        }
        if (name == "parity") {
            if (arguments.size() != 1) throw MathError("parity expects 1 argument");
            return static_cast<double>(parity_bits(
                to_unsigned_bits(require_integer_argument(arguments[0], "argument"))));
        }
        if (name == "reverse_bits") {
            if (arguments.size() != 1) throw MathError("reverse_bits expects 1 argument");
            return static_cast<double>(from_unsigned_bits(reverse_bits(
                to_unsigned_bits(require_integer_argument(arguments[0], "argument")))));
        }
        if (name == "step" || name == "heaviside") {
            if (arguments.size() != 1) throw MathError(name + " expects 1 argument");
            return arguments[0] < 0.0 ? 0.0 : 1.0;
        }
        if (name == "delta" || name == "impulse") {
            if (arguments.size() != 1) throw MathError(name + " expects 1 argument");
            return mymath::is_near_zero(arguments[0], 1e-12) ? 1.0 : 0.0;
        }
        if (name == "poisson_pmf" || name == "poisson_cdf" ||
            name == "binom_pmf" || name == "binom_cdf") {
            return stats_ops::apply_probability(name, arguments);
        }
        if (name == "rat") {
            if (arguments.size() == 1) return arguments[0];
            if (arguments.size() == 2) return arguments[0];
            throw MathError("rat expects 1 or 2 arguments");
        }

        throw UndefinedError("unknown function: " + name);
    }

    VariableResolver variables_;
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionDecimalCallback invoke_script_function_;
};

double parse_decimal_expression(
    const std::string& expression,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionDecimalCallback invoke_script_function) {
    DecimalParserImpl parser(expression,
                             variables,
                             functions,
                             std::move(has_script_function),
                             std::move(invoke_script_function));
    return parser.parse();
}

DecimalParser::DecimalParser(
    std::string source,
    const VariableResolver& variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionDecimalCallback invoke_script_function)
    : source_(std::move(source)),
      variables_(variables),
      functions_(functions),
      has_script_function_(std::move(has_script_function)),
      invoke_script_function_(std::move(invoke_script_function)) {}

double DecimalParser::parse() {
    DecimalParserImpl parser(source_,
                             variables_,
                             functions_,
                             std::move(has_script_function_),
                             std::move(invoke_script_function_));
    return parser.parse();
}

bool try_evaluate_matrix_expression(const std::string& expression,
                                    const VariableResolver& variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const HasScriptFunctionCallback& has_script_function,
                                    const InvokeScriptFunctionDecimalCallback& invoke_script_function,
                                    matrix::Value* value) {
    const matrix::ScalarEvaluator scalar_evaluator =
        [variables, functions, has_script_function, invoke_script_function](const std::string& text) {
            DecimalParserImpl parser(text,
                                     variables,
                                     functions,
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
                                           value);
}

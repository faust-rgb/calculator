#include "calculator_internal_types.h"

#include "matrix.h"
#include "mymath.h"
#include "statistics/calculator_statistics.h"

#include <algorithm>
#include <cctype>
#include <map>
#include <random>
#include <stdexcept>
#include <utility>

class DecimalParserImpl {
public:
    DecimalParserImpl(std::string source,
                      const std::map<std::string, StoredValue>* variables,
                      const std::map<std::string, CustomFunction>* functions,
                      HasScriptFunctionCallback has_script_function = {},
                      InvokeScriptFunctionDecimalCallback invoke_script_function = {})
        : source_(std::move(source)),
          variables_(variables),
          functions_(functions),
          has_script_function_(std::move(has_script_function)),
          invoke_script_function_(std::move(invoke_script_function)) {}

    double parse() {
        double value = parse_comparison();
        skip_spaces();
        if (!is_at_end()) {
            throw std::runtime_error("unexpected token near: " + source_.substr(pos_, 1));
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
                if (mymath::is_near_zero(divisor)) {
                    throw std::runtime_error("division by zero");
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

        // 先尝试解析 0b/0o/0x 这类前缀整数。
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
            throw std::runtime_error("expected number");
        }

        return parse_decimal(source_.substr(start, pos_ - start));
    }

    static double parse_decimal(const std::string& token) {
        return std::stod(token);
    }

    double apply_function(const std::string& name, const std::vector<double>& arguments) {
        // 这里是“普通浮点路径”的函数分发中心。
        // 需要字符串结果的功能，例如 factor/bin/hex，不从这里返回。
        const BuiltinFunction* builtin = find_builtin_function(name);
        if (builtin != nullptr) {
            return (*builtin)(arguments);
        }

        const auto function_it = functions_->find(name);
        if (function_it != functions_->end()) {
            if (arguments.size() != 1) {
                throw std::runtime_error("custom function " + name +
                                         " expects exactly one argument");
            }

            std::map<std::string, StoredValue> scoped_variables = *variables_;
            StoredValue parameter_value;
            parameter_value.decimal = arguments[0];
            parameter_value.exact = false;
            scoped_variables[function_it->second.parameter_name] = parameter_value;

            DecimalParserImpl nested_parser(function_it->second.expression,
                                            &scoped_variables,
                                            functions_,
                                            has_script_function_,
                                            invoke_script_function_);
            return nested_parser.parse();
        }

        if (has_script_function_ && has_script_function_(name)) {
            return invoke_script_function_(name, arguments);
        }

        if (arguments.size() != 1) {
            throw std::runtime_error("function " + name + " expects exactly one argument");
        }

        const double argument = arguments[0];
        if (name == "abs") {
            return mymath::abs(argument);
        }
        if (name == "not") {
            // 位运算函数只接受整数，浮点路径里也保持这个约束。
            if (!is_integer_double(argument)) {
                throw std::runtime_error("not only accepts integers");
            }
            return static_cast<double>(~round_to_long_long(argument));
        }
        if (name == "sign") {
            if (mymath::is_near_zero(argument)) {
                return 0.0;
            }
            return argument > 0.0 ? 1.0 : -1.0;
        }
        if (name == "floor") {
            return static_cast<double>(floor_to_long_long(argument));
        }
        if (name == "ceil") {
            return static_cast<double>(ceil_to_long_long(argument));
        }
        if (name == "round") {
            return static_cast<double>(round_to_long_long(argument));
        }
        if (name == "trunc") {
            return static_cast<double>(trunc_to_long_long(argument));
        }
        if (name == "cbrt") {
            return mymath::cbrt(argument);
        }
        if (name == "asinh") {
            return mymath::asinh(argument);
        }
        if (name == "acosh") {
            return mymath::acosh(argument);
        }
        if (name == "atanh") {
            return mymath::atanh(argument);
        }
        if (name == "sinh") {
            return mymath::sinh(argument);
        }
        if (name == "cosh") {
            return mymath::cosh(argument);
        }
        if (name == "tanh") {
            return mymath::tanh(argument);
        }
        if (name == "sec") {
            return mymath::sec(argument);
        }
        if (name == "csc") {
            return mymath::csc(argument);
        }
        if (name == "cot") {
            return mymath::cot(argument);
        }
        if (name == "sin") {
            return mymath::sin(argument);
        }
        if (name == "cos") {
            return mymath::cos(argument);
        }
        if (name == "tan") {
            return mymath::tan(argument);
        }
        if (name == "atan") {
            return mymath::atan(argument);
        }
        if (name == "asin") {
            return mymath::asin(argument);
        }
        if (name == "acos") {
            return mymath::acos(argument);
        }
        if (name == "asec") {
            return mymath::asec(argument);
        }
        if (name == "acsc") {
            return mymath::acsc(argument);
        }
        if (name == "acot") {
            return mymath::acot(argument);
        }
        if (name == "ln") {
            return mymath::ln(argument);
        }
        if (name == "log2") {
            return mymath::ln(argument) / mymath::ln(2.0);
        }
        if (name == "log10") {
            return mymath::log10(argument);
        }
        if (name == "exp") {
            return mymath::exp(argument);
        }
        if (name == "exp2") {
            return mymath::exp(argument * mymath::ln(2.0));
        }
        if (name == "gamma") {
            return mymath::gamma(argument);
        }
        if (name == "erf") {
            return mymath::erf(argument);
        }
        if (name == "erfc") {
            return mymath::erfc(argument);
        }
        if (name == "sqrt") {
            return mymath::sqrt(argument);
        }
        if (name == "step" || name == "u" || name == "heaviside") {
            return argument >= 0.0 ? 1.0 : 0.0;
        }
        if (name == "delta" || name == "impulse") {
            return mymath::is_near_zero(argument, 1e-10) ? 1.0 : 0.0;
        }
        if (name == "deg") {
            return radians_to_degrees(argument);
        }
        if (name == "rad") {
            return degrees_to_radians(argument);
        }
        if (name == "deg2rad") {
            return degrees_to_radians(argument);
        }
        if (name == "rad2deg") {
            return radians_to_degrees(argument);
        }
        if (name == "sin_deg") {
            return mymath::sin(degrees_to_radians(argument));
        }
        if (name == "cos_deg") {
            return mymath::cos(degrees_to_radians(argument));
        }
        if (name == "celsius") {
            return fahrenheit_to_celsius(argument);
        }
        if (name == "fahrenheit") {
            return celsius_to_fahrenheit(argument);
        }
        if (name == "kelvin") {
            return argument + 273.15;
        }
        if (name == "c2f") {
            return celsius_to_fahrenheit(argument);
        }
        if (name == "f2c") {
            return fahrenheit_to_celsius(argument);
        }

        throw std::runtime_error("unknown function: " + name);
    }

    using BuiltinFunction = double (*)(const std::vector<double>&);

    struct BuiltinFunctionEntry {
        const char* name;
        BuiltinFunction function;
    };

    static const BuiltinFunction* find_builtin_function(const std::string& name) {
        static const BuiltinFunctionEntry functions[] = {
            {"and", apply_and},
            {"avg", apply_avg},
            {"bessel", apply_bessel},
            {"beta", apply_beta},
            {"binom", apply_ncr},
            {"binom_cdf", apply_binom_cdf},
            {"binom_pmf", apply_binom_pmf},
            {"bitlen", apply_bitlen},
            {"cdf_normal", apply_cdf_normal},
            {"clamp", apply_clamp},
            {"clz", apply_clz},
            {"ctz", apply_ctz},
            {"egcd", apply_egcd},
            {"euler_phi", apply_euler_phi},
            {"factorial", apply_factorial},
            {"fib", apply_fib},
            {"gcd", apply_gcd},
            {"is_prime", apply_is_prime},
            {"kurtosis", apply_kurtosis},
            {"lcm", apply_lcm},
            {"log", apply_log},
            {"max", apply_max},
            {"mean", apply_mean},
            {"median", apply_median},
            {"min", apply_min},
            {"mobius", apply_mobius},
            {"mod", apply_mod},
            {"mode", apply_mode},
            {"nCr", apply_ncr},
            {"next_prime", apply_next_prime},
            {"nPr", apply_npr},
            {"or", apply_or},
            {"parity", apply_parity},
            {"pdf_normal", apply_pdf_normal},
            {"percentile", apply_percentile},
            {"phi", apply_euler_phi},
            {"poisson_cdf", apply_poisson_cdf},
            {"poisson_pmf", apply_poisson_pmf},
            {"popcount", apply_popcount},
            {"pow", apply_pow},
            {"prev_prime", apply_prev_prime},
            {"prime_pi", apply_prime_pi},
            {"quartile", apply_quartile},
            {"rand", apply_rand},
            {"randint", apply_randint},
            {"randn", apply_randn},
            {"reverse_bits", apply_reverse_bits},
            {"rol", apply_rol},
            {"root", apply_root},
            {"ror", apply_ror},
            {"shl", apply_shl},
            {"shr", apply_shr},
            {"skew", apply_skewness},
            {"skewness", apply_skewness},
            {"std", apply_stddev},
            {"sum", apply_sum},
            {"var", apply_variance},
            {"xor", apply_xor},
            {"zeta", apply_zeta},
        };

        for (const BuiltinFunctionEntry& entry : functions) {
            if (name == entry.name) {
                return &entry.function;
            }
        }
        return nullptr;
    }

    static double apply_gcd(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("gcd expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("gcd only accepts integers");
        }
        return static_cast<double>(
            gcd_ll(round_to_long_long(arguments[0]), round_to_long_long(arguments[1])));
    }

    static double apply_pow(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("pow expects exactly two arguments");
        }
        return mymath::pow(arguments[0], arguments[1]);
    }

    static double apply_root(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("root expects exactly two arguments");
        }
        return mymath::root(arguments[0], arguments[1]);
    }

    static double apply_lcm(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("lcm expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("lcm only accepts integers");
        }
        return static_cast<double>(
            lcm_ll(round_to_long_long(arguments[0]), round_to_long_long(arguments[1])));
    }

    static double apply_mod(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("mod expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("mod only accepts integers");
        }

        const long long lhs = round_to_long_long(arguments[0]);
        const long long rhs = round_to_long_long(arguments[1]);
        if (rhs == 0) {
            throw std::runtime_error("mod divisor cannot be zero");
        }
        return static_cast<double>(lhs % rhs);
    }

    static double apply_min(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("min expects exactly two arguments");
        }
        return arguments[0] < arguments[1] ? arguments[0] : arguments[1];
    }

    static double apply_max(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("max expects exactly two arguments");
        }
        return arguments[0] > arguments[1] ? arguments[0] : arguments[1];
    }

    static double apply_clamp(const std::vector<double>& arguments) {
        if (arguments.size() != 3) {
            throw std::runtime_error("clamp expects exactly three arguments");
        }
        double lower = arguments[1];
        double upper = arguments[2];
        if (lower > upper) {
            std::swap(lower, upper);
        }
        if (arguments[0] < lower) {
            return lower;
        }
        if (arguments[0] > upper) {
            return upper;
        }
        return arguments[0];
    }

    static double apply_log(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("log expects exactly two arguments");
        }
        if (mymath::is_near_zero(arguments[1] - 1.0)) {
            throw std::runtime_error("log base cannot be 1");
        }
        return mymath::ln(arguments[0]) / mymath::ln(arguments[1]);
    }

    static double apply_sum(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("sum expects at least one argument");
        }
        long double total = 0.0L;
        long double compensation = 0.0L;
        for (double value : arguments) {
            const long double adjusted = static_cast<long double>(value) - compensation;
            const long double next = total + adjusted;
            compensation = (next - total) - adjusted;
            total = next;
        }
        return static_cast<double>(total);
    }

    static double apply_avg(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("avg", arguments);
    }

    static double apply_mean(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("mean", arguments);
    }

    static double apply_median(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("median", arguments);
    }

    static double apply_mode(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("mode", arguments);
    }

    static double apply_variance(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("var", arguments);
    }

    static double apply_stddev(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("std", arguments);
    }

    static double apply_skewness(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("skewness", arguments);
    }

    static double apply_kurtosis(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("kurtosis", arguments);
    }

    static double apply_percentile(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("percentile", arguments);
    }

    static double apply_quartile(const std::vector<double>& arguments) {
        return stats_ops::apply_statistic("quartile", arguments);
    }

    static double apply_factorial(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("factorial", arguments);
    }

    static double apply_ncr(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("nCr", arguments);
    }

    static double apply_npr(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("nPr", arguments);
    }

    static double apply_rand(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("rand", arguments);
    }

    static double apply_randn(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("randn", arguments);
    }

    static double apply_randint(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("randint", arguments);
    }

    static double apply_pdf_normal(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("pdf_normal", arguments);
    }

    static double apply_cdf_normal(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("cdf_normal", arguments);
    }

    static double apply_poisson_pmf(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("poisson_pmf", arguments);
    }

    static double apply_poisson_cdf(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("poisson_cdf", arguments);
    }

    static double apply_binom_pmf(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("binom_pmf", arguments);
    }

    static double apply_binom_cdf(const std::vector<double>& arguments) {
        return stats_ops::apply_probability("binom_cdf", arguments);
    }

    static double apply_fib(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("fib expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("fib only accepts integers");
        }
        return fibonacci_value(round_to_long_long(arguments[0]));
    }

    static double apply_is_prime(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("is_prime expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("is_prime only accepts integers");
        }
        return is_prime_ll(round_to_long_long(arguments[0])) ? 1.0 : 0.0;
    }

    static double apply_next_prime(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("next_prime expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("next_prime only accepts integers");
        }
        return static_cast<double>(next_prime_ll(round_to_long_long(arguments[0])));
    }

    static double apply_prev_prime(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("prev_prime expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("prev_prime only accepts integers");
        }
        return static_cast<double>(prev_prime_ll(round_to_long_long(arguments[0])));
    }

    static double apply_prime_pi(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("prime_pi expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("prime_pi only accepts integers");
        }
        return static_cast<double>(prime_pi_ll(round_to_long_long(arguments[0])));
    }

    static double apply_euler_phi(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("euler_phi expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("euler_phi only accepts integers");
        }
        return static_cast<double>(euler_phi_ll(round_to_long_long(arguments[0])));
    }

    static double apply_mobius(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("mobius expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("mobius only accepts integers");
        }
        return static_cast<double>(mobius_ll(round_to_long_long(arguments[0])));
    }

    static double apply_egcd(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("egcd expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("egcd only accepts integers");
        }
        long long x = 0;
        long long y = 0;
        return static_cast<double>(extended_gcd_ll(round_to_long_long(arguments[0]),
                                                  round_to_long_long(arguments[1]),
                                                  &x,
                                                  &y));
    }

    static double apply_beta(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("beta expects exactly two arguments");
        }
        return mymath::beta(arguments[0], arguments[1]);
    }

    static double apply_zeta(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("zeta expects exactly one argument");
        }
        return mymath::zeta(arguments[0]);
    }

    static double apply_bessel(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("bessel expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("bessel order must be an integer");
        }
        return mymath::bessel_j(static_cast<int>(round_to_long_long(arguments[0])),
                                arguments[1]);
    }

    static double apply_and(const std::vector<double>& arguments) {
        // 位运算统一通过整数检查后再执行，避免隐式截断带来困惑。
        if (arguments.size() != 2) {
            throw std::runtime_error("and expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("and only accepts integers");
        }
        return static_cast<double>(
            round_to_long_long(arguments[0]) & round_to_long_long(arguments[1]));
    }

    static double apply_or(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("or expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("or only accepts integers");
        }
        return static_cast<double>(
            round_to_long_long(arguments[0]) | round_to_long_long(arguments[1]));
    }

    static double apply_xor(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("xor expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("xor only accepts integers");
        }
        return static_cast<double>(
            round_to_long_long(arguments[0]) ^ round_to_long_long(arguments[1]));
    }

    static double apply_shl(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("shl expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("shl only accepts integers");
        }
        const long long shift = round_to_long_long(arguments[1]);
        if (shift < 0) {
            throw std::runtime_error("shift count cannot be negative");
        }
        return static_cast<double>(round_to_long_long(arguments[0]) << shift);
    }

    static double apply_shr(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("shr expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("shr only accepts integers");
        }
        const long long shift = round_to_long_long(arguments[1]);
        if (shift < 0) {
            throw std::runtime_error("shift count cannot be negative");
        }
        return static_cast<double>(round_to_long_long(arguments[0]) >> shift);
    }

    static double apply_rol(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("rol expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("rol only accepts integers");
        }
        const std::uint64_t value = to_unsigned_bits(round_to_long_long(arguments[0]));
        const unsigned count = normalize_rotation_count(round_to_long_long(arguments[1]));
        return static_cast<double>(from_unsigned_bits(rotate_left_bits(value, count)));
    }

    static double apply_ror(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("ror expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("ror only accepts integers");
        }
        const std::uint64_t value = to_unsigned_bits(round_to_long_long(arguments[0]));
        const unsigned count = normalize_rotation_count(round_to_long_long(arguments[1]));
        return static_cast<double>(from_unsigned_bits(rotate_right_bits(value, count)));
    }

    static double apply_popcount(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("popcount expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("popcount only accepts integers");
        }
        return static_cast<double>(
            popcount_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_bitlen(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("bitlen expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("bitlen only accepts integers");
        }
        return static_cast<double>(
            bit_length_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_ctz(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("ctz expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("ctz only accepts integers");
        }
        return static_cast<double>(
            trailing_zero_count_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_clz(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("clz expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("clz only accepts integers");
        }
        return static_cast<double>(
            leading_zero_count_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_parity(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("parity expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("parity only accepts integers");
        }
        return static_cast<double>(
            parity_bits(to_unsigned_bits(round_to_long_long(arguments[0]))));
    }

    static double apply_reverse_bits(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("reverse_bits expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("reverse_bits only accepts integers");
        }
        return static_cast<double>(
            from_unsigned_bits(reverse_bits(to_unsigned_bits(round_to_long_long(arguments[0])))));
    }

    double lookup_variable(const std::string& name) const {
        const auto it = variables_->find(name);
        if (it == variables_->end()) {
            double constant_value = 0.0;
            if (lookup_builtin_constant(name, &constant_value)) {
                return constant_value;
            }
            throw std::runtime_error("unknown variable: " + name);
        }
        if (it->second.is_matrix) {
            throw std::runtime_error("matrix variable " + name + " cannot be used as a scalar");
        }
        if (it->second.is_complex) {
            throw std::runtime_error("complex variable " + name + " cannot be used as a real scalar");
        }
        if (it->second.is_string) {
            throw std::runtime_error("string variable " + name + " cannot be used as a number");
        }
        return it->second.exact ? rational_to_double(it->second.rational)
                                : it->second.decimal;
    }

    std::string parse_identifier() {
        const std::size_t start = pos_;
        while (!is_at_end()) {
            const char ch = source_[pos_];
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
                ++pos_;
            } else {
                break;
            }
        }
        return source_.substr(start, pos_ - start);
    }

    bool peek_is_alpha() const {
        return !is_at_end() &&
               std::isalpha(static_cast<unsigned char>(source_[pos_]));
    }

    bool peek(char expected) const {
        return !is_at_end() && source_[pos_] == expected;
    }

    void skip_spaces() {
        while (!is_at_end() &&
               std::isspace(static_cast<unsigned char>(source_[pos_]))) {
            ++pos_;
        }
    }

    bool match(char expected) {
        if (is_at_end() || source_[pos_] != expected) {
            return false;
        }
        ++pos_;
        return true;
    }

    bool match_string(const std::string& text) {
        if (source_.compare(pos_, text.size(), text) != 0) {
            return false;
        }
        pos_ += text.size();
        return true;
    }

    void expect(char expected) {
        if (!match(expected)) {
            throw std::runtime_error(std::string("expected '") + expected + "'");
        }
    }

    bool is_at_end() const {
        return pos_ >= source_.size();
    }

    std::string source_;
    std::size_t pos_ = 0;
    const std::map<std::string, StoredValue>* variables_;
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionDecimalCallback invoke_script_function_;
};

double parse_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables,
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
    const std::map<std::string, StoredValue>* variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function,
    InvokeScriptFunctionDecimalCallback invoke_script_function)
    : source_(std::move(source)),
      variables_(variables),
      functions_(functions),
      has_script_function_(std::move(has_script_function)),
      invoke_script_function_(std::move(invoke_script_function)) {}

double DecimalParser::parse() {
    return parse_decimal_expression(source_,
                                    variables_,
                                    functions_,
                                    std::move(has_script_function_),
                                    std::move(invoke_script_function_));
}

bool try_evaluate_matrix_expression(const std::string& expression,
                                    const std::map<std::string, StoredValue>* variables,
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
            const auto it = variables->find(name);
            if (it == variables->end() || !it->second.is_matrix) {
                return false;
            }
            *matrix_value = it->second.matrix;
            return true;
        };
    const matrix::ComplexLookup complex_lookup =
        [variables](const std::string& name, matrix::ComplexNumber* complex_value) {
            const auto it = variables->find(name);
            if (it == variables->end() || !it->second.is_complex) {
                return false;
            }
            *complex_value = it->second.complex;
            return true;
        };
    return matrix::try_evaluate_expression(expression,
                                           scalar_evaluator,
                                           matrix_lookup,
                                           complex_lookup,
                                           value);
}

// ExactParser 只处理能够保持为有理数的表达式。
// 当遇到 sin、pi 或非整数指数这类无法精确表示为分数的情况时，
// 它会抛出 ExactModeUnsupported，调用方再回退到普通浮点模式。

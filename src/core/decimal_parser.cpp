#include "calculator_internal_types.h"

#include "matrix.h"
#include "mymath.h"

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
        if (name == "pow") {
            return apply_pow(arguments);
        }
        if (name == "root") {
            return apply_root(arguments);
        }
        if (name == "and") {
            return apply_and(arguments);
        }
        if (name == "or") {
            return apply_or(arguments);
        }
        if (name == "xor") {
            return apply_xor(arguments);
        }
        if (name == "shl") {
            return apply_shl(arguments);
        }
        if (name == "shr") {
            return apply_shr(arguments);
        }
        if (name == "rol") {
            return apply_rol(arguments);
        }
        if (name == "ror") {
            return apply_ror(arguments);
        }
        if (name == "gcd") {
            return apply_gcd(arguments);
        }
        if (name == "lcm") {
            return apply_lcm(arguments);
        }
        if (name == "mod") {
            return apply_mod(arguments);
        }
        if (name == "min") {
            return apply_min(arguments);
        }
        if (name == "max") {
            return apply_max(arguments);
        }
        if (name == "clamp") {
            return apply_clamp(arguments);
        }
        if (name == "log") {
            return apply_log(arguments);
        }
        if (name == "sum") {
            return apply_sum(arguments);
        }
        if (name == "mean") {
            return apply_mean(arguments);
        }
        if (name == "avg") {
            return apply_avg(arguments);
        }
        if (name == "median") {
            return apply_median(arguments);
        }
        if (name == "mode") {
            return apply_mode(arguments);
        }
        if (name == "var") {
            return apply_variance(arguments);
        }
        if (name == "std") {
            return apply_stddev(arguments);
        }
        if (name == "percentile") {
            return apply_percentile(arguments);
        }
        if (name == "quartile") {
            return apply_quartile(arguments);
        }
        if (name == "factorial") {
            return apply_factorial(arguments);
        }
        if (name == "nCr") {
            return apply_ncr(arguments);
        }
        if (name == "binom") {
            return apply_ncr(arguments);
        }
        if (name == "nPr") {
            return apply_npr(arguments);
        }
        if (name == "fib") {
            return apply_fib(arguments);
        }
        if (name == "is_prime") {
            return apply_is_prime(arguments);
        }
        if (name == "next_prime") {
            return apply_next_prime(arguments);
        }
        if (name == "rand") {
            return apply_rand(arguments);
        }
        if (name == "randn") {
            return apply_randn(arguments);
        }
        if (name == "randint") {
            return apply_randint(arguments);
        }
        if (name == "beta") {
            return apply_beta(arguments);
        }
        if (name == "zeta") {
            return apply_zeta(arguments);
        }
        if (name == "bessel") {
            return apply_bessel(arguments);
        }
        if (name == "pdf_normal") {
            return apply_pdf_normal(arguments);
        }
        if (name == "cdf_normal") {
            return apply_cdf_normal(arguments);
        }
        if (name == "popcount") {
            return apply_popcount(arguments);
        }
        if (name == "bitlen") {
            return apply_bitlen(arguments);
        }
        if (name == "ctz") {
            return apply_ctz(arguments);
        }
        if (name == "clz") {
            return apply_clz(arguments);
        }
        if (name == "parity") {
            return apply_parity(arguments);
        }
        if (name == "reverse_bits") {
            return apply_reverse_bits(arguments);
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
        double total = 0.0;
        for (double value : arguments) {
            total += value;
        }
        return total;
    }

    static double apply_avg(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("avg expects at least one argument");
        }
        return apply_sum(arguments) / static_cast<double>(arguments.size());
    }

    static double apply_mean(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("mean expects at least one argument");
        }
        return apply_sum(arguments) / static_cast<double>(arguments.size());
    }

    static double apply_median(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("median expects at least one argument");
        }

        std::vector<double> sorted = arguments;
        std::sort(sorted.begin(), sorted.end());
        const std::size_t middle = sorted.size() / 2;
        if (sorted.size() % 2 == 1) {
            return sorted[middle];
        }
        return (sorted[middle - 1] + sorted[middle]) / 2.0;
    }

    static double apply_mode(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("mode expects at least one argument");
        }
        std::vector<double> sorted = arguments;
        std::sort(sorted.begin(), sorted.end());
        double best_value = sorted.front();
        int best_count = 1;
        double current_value = sorted.front();
        int current_count = 1;
        for (std::size_t i = 1; i < sorted.size(); ++i) {
            if (mymath::is_near_zero(sorted[i] - current_value, 1e-10)) {
                ++current_count;
                continue;
            }
            if (current_count > best_count) {
                best_count = current_count;
                best_value = current_value;
            }
            current_value = sorted[i];
            current_count = 1;
        }
        if (current_count > best_count) {
            best_value = current_value;
        }
        return best_value;
    }

    static double apply_variance(const std::vector<double>& arguments) {
        if (arguments.empty()) {
            throw std::runtime_error("var expects at least one argument");
        }
        const double mean = apply_mean(arguments);
        double sum = 0.0;
        for (double value : arguments) {
            const double delta = value - mean;
            sum += delta * delta;
        }
        return sum / static_cast<double>(arguments.size());
    }

    static double apply_stddev(const std::vector<double>& arguments) {
        return mymath::sqrt(apply_variance(arguments));
    }

    static double apply_percentile(const std::vector<double>& arguments) {
        if (arguments.size() < 2) {
            throw std::runtime_error("percentile expects p followed by at least one value");
        }
        const double p = arguments[0];
        if (p < 0.0 || p > 100.0) {
            throw std::runtime_error("percentile p must be in [0, 100]");
        }
        std::vector<double> values(arguments.begin() + 1, arguments.end());
        std::sort(values.begin(), values.end());
        if (values.size() == 1) {
            return values.front();
        }
        const double position =
            p * static_cast<double>(values.size() - 1) / 100.0;
        const long long lower_index = floor_to_long_long(position);
        const long long upper_index = ceil_to_long_long(position);
        if (lower_index == upper_index) {
            return values[static_cast<std::size_t>(lower_index)];
        }
        const double fraction = position - static_cast<double>(lower_index);
        const double lower = values[static_cast<std::size_t>(lower_index)];
        const double upper = values[static_cast<std::size_t>(upper_index)];
        return lower + (upper - lower) * fraction;
    }

    static double apply_quartile(const std::vector<double>& arguments) {
        if (arguments.size() < 2) {
            throw std::runtime_error("quartile expects q followed by at least one value");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("quartile q must be an integer");
        }
        const long long q = round_to_long_long(arguments[0]);
        if (q < 0 || q > 4) {
            throw std::runtime_error("quartile q must be between 0 and 4");
        }
        std::vector<double> percentile_arguments;
        percentile_arguments.reserve(arguments.size());
        percentile_arguments.push_back(static_cast<double>(q * 25));
        percentile_arguments.insert(percentile_arguments.end(),
                                    arguments.begin() + 1,
                                    arguments.end());
        return apply_percentile(percentile_arguments);
    }

    static double apply_factorial(const std::vector<double>& arguments) {
        if (arguments.size() != 1) {
            throw std::runtime_error("factorial expects exactly one argument");
        }
        if (!is_integer_double(arguments[0])) {
            throw std::runtime_error("factorial only accepts integers");
        }
        return factorial_value(round_to_long_long(arguments[0]));
    }

    static double apply_ncr(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("nCr expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("nCr only accepts integers");
        }
        return combination_value(round_to_long_long(arguments[0]),
                                 round_to_long_long(arguments[1]));
    }

    static double apply_npr(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("nPr expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("nPr only accepts integers");
        }
        return permutation_value(round_to_long_long(arguments[0]),
                                 round_to_long_long(arguments[1]));
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

    static double apply_rand(const std::vector<double>& arguments) {
        if (!arguments.empty()) {
            throw std::runtime_error("rand expects no arguments");
        }
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return distribution(global_rng());
    }

    static double apply_randn(const std::vector<double>& arguments) {
        if (!arguments.empty()) {
            throw std::runtime_error("randn expects no arguments");
        }
        std::normal_distribution<double> distribution(0.0, 1.0);
        return distribution(global_rng());
    }

    static double apply_randint(const std::vector<double>& arguments) {
        if (arguments.size() != 2) {
            throw std::runtime_error("randint expects exactly two arguments");
        }
        if (!is_integer_double(arguments[0]) || !is_integer_double(arguments[1])) {
            throw std::runtime_error("randint only accepts integers");
        }
        long long left = round_to_long_long(arguments[0]);
        long long right = round_to_long_long(arguments[1]);
        if (left > right) {
            std::swap(left, right);
        }
        std::uniform_int_distribution<long long> distribution(left, right);
        return static_cast<double>(distribution(global_rng()));
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

    static double apply_pdf_normal(const std::vector<double>& arguments) {
        if (arguments.size() != 3) {
            throw std::runtime_error("pdf_normal expects exactly three arguments");
        }
        return normal_pdf(arguments[0], arguments[1], arguments[2]);
    }

    static double apply_cdf_normal(const std::vector<double>& arguments) {
        if (arguments.size() != 3) {
            throw std::runtime_error("cdf_normal expects exactly three arguments");
        }
        return normal_cdf(arguments[0], arguments[1], arguments[2]);
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
    return matrix::try_evaluate_expression(expression,
                                           scalar_evaluator,
                                           matrix_lookup,
                                           value);
}

// ExactParser 只处理能够保持为有理数的表达式。
// 当遇到 sin、pi 或非整数指数这类无法精确表示为分数的情况时，
// 它会抛出 ExactModeUnsupported，调用方再回退到普通浮点模式。

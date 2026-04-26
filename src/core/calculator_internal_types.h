#ifndef CALCULATOR_INTERNAL_TYPES_H
#define CALCULATOR_INTERNAL_TYPES_H

#include "calculator.h"

#include "matrix.h"
#include "mymath.h"
#include "script_ast.h"

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

constexpr double kDisplayZeroEps = mymath::kDoubleDenormMin;
constexpr double kDisplayIntegerEps = 1e-9;

class ExactModeUnsupported : public std::runtime_error {
public:
    explicit ExactModeUnsupported(const std::string& message)
        : std::runtime_error(message) {}
};

class PreciseDecimalUnsupported : public std::runtime_error {
public:
    explicit PreciseDecimalUnsupported(const std::string& message)
        : std::runtime_error(message) {}
};

struct Rational {
    long long numerator = 0;
    long long denominator = 1;

    Rational() = default;
    Rational(long long num, long long den);
    void normalize();
    bool is_integer() const;
    std::string to_string() const;
};

Rational operator+(const Rational& lhs, const Rational& rhs);
Rational operator-(const Rational& lhs, const Rational& rhs);
Rational operator*(const Rational& lhs, const Rational& rhs);
Rational operator/(const Rational& lhs, const Rational& rhs);

struct PreciseDecimal {
    std::string digits = "0";
    int scale = 0;
    bool negative = false;

    void normalize();
    bool is_zero() const;
    std::string to_string() const;
    double to_double() const;
    static PreciseDecimal from_digits(std::string raw_digits,
                                      int raw_scale,
                                      bool is_negative);
    static PreciseDecimal from_integer_string(const std::string& integer_text,
                                              bool is_negative);
    static PreciseDecimal from_decimal_literal(const std::string& token);
};

struct StoredValue {
    bool is_matrix = false;
    bool is_string = false;
    bool has_symbolic_text = false;
    bool has_precise_decimal_text = false;
    bool exact = false;
    Rational rational;
    double decimal = 0.0;
    std::string string_value;
    std::string symbolic_text;
    std::string precise_decimal_text;
    matrix::Matrix matrix;
};

struct CustomFunction {
    std::string parameter_name;
    std::string expression;
};

struct ScriptFunction {
    std::vector<std::string> parameter_names;
    std::shared_ptr<const script::BlockStatement> body;
};

using HasScriptFunctionCallback = std::function<bool(const std::string&)>;
using InvokeScriptFunctionDecimalCallback =
    std::function<double(const std::string&, const std::vector<double>&)>;

class DecimalParser {
public:
    DecimalParser(std::string source,
                  const std::map<std::string, StoredValue>* variables,
                  const std::map<std::string, CustomFunction>* functions,
                  HasScriptFunctionCallback has_script_function = {},
                  InvokeScriptFunctionDecimalCallback invoke_script_function = {});

    double parse();

private:
    std::string source_;
    const std::map<std::string, StoredValue>* variables_;
    const std::map<std::string, CustomFunction>* functions_;
    HasScriptFunctionCallback has_script_function_;
    InvokeScriptFunctionDecimalCallback invoke_script_function_;
};

struct Calculator::Impl {
    std::map<std::string, StoredValue> variables;
    std::map<std::string, CustomFunction> functions;
    std::map<std::string, ScriptFunction> script_functions;
    std::vector<std::map<std::string, StoredValue>> local_scopes;
    bool symbolic_constants_mode = false;
    bool hex_prefix_mode = false;
    bool hex_uppercase_mode = true;
};

struct ScriptSignal {
    enum class Kind {
        kNone,
        kReturn,
        kBreak,
        kContinue,
    };

    Kind kind = Kind::kNone;
    bool has_value = false;
    StoredValue value;

    static ScriptSignal make_return(const StoredValue& return_value);
    static ScriptSignal make_break();
    static ScriptSignal make_continue();
};

struct HexFormatOptions {
    bool prefix = false;
    bool uppercase = true;
};

long long gcd_ll(long long a, long long b);
long long lcm_ll(long long a, long long b);
bool is_integer_double(double x, double eps = 1e-10);
long long round_to_long_long(double x);
long long trunc_to_long_long(double x);
long long floor_to_long_long(double x);
long long ceil_to_long_long(double x);
double normalize_display_decimal(double value);
std::mt19937_64& global_rng();
bool lookup_builtin_constant(const std::string& name, double* value);
double degrees_to_radians(double value);
double radians_to_degrees(double value);
double celsius_to_fahrenheit(double value);
double fahrenheit_to_celsius(double value);
double normal_pdf(double x, double mean, double sigma);
double normal_cdf(double x, double mean, double sigma);
bool is_prime_ll(long long value);
long long next_prime_ll(long long value);
long long prev_prime_ll(long long value);
long long euler_phi_ll(long long value);
long long mobius_ll(long long value);
long long prime_pi_ll(long long value);
long long extended_gcd_ll(long long a, long long b, long long* x, long long* y);
double fibonacci_value(long long n);
std::uint64_t to_unsigned_bits(long long value);
long long from_unsigned_bits(std::uint64_t value);
unsigned normalize_rotation_count(long long count);
std::uint64_t rotate_left_bits(std::uint64_t value, unsigned count);
std::uint64_t rotate_right_bits(std::uint64_t value, unsigned count);
int popcount_bits(std::uint64_t value);
int bit_length_bits(std::uint64_t value);
int trailing_zero_count_bits(std::uint64_t value);
int leading_zero_count_bits(std::uint64_t value);
int parity_bits(std::uint64_t value);
std::uint64_t reverse_bits(std::uint64_t value);
std::string factor_integer(long long value);
int digit_value(char ch);
bool prefixed_base(char prefix, int* base);
long long parse_prefixed_integer_token(const std::string& token);
Rational pow_rational(Rational base, long long exponent);
Rational abs_rational(Rational value);
std::string format_decimal(double value);
std::string format_symbolic_number(double value);
double root_position_tolerance(double value);
double root_function_tolerance(double value);
double root_derivative_step(double value);
PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs,
                                   const PreciseDecimal& rhs);
PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs,
                                        const PreciseDecimal& rhs);
PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs,
                                        const PreciseDecimal& rhs);
PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs,
                                      const PreciseDecimal& rhs);
double rational_to_double(const Rational& value);
bool is_valid_variable_name(const std::string& name);
std::string trim_copy(const std::string& text);
bool is_identifier_text(const std::string& text);
bool is_string_literal(const std::string& text);
std::string parse_string_literal_value(const std::string& text);
std::string encode_state_field(const std::string& text);
std::string decode_state_field(const std::string& text);
bool split_assignment(const std::string& expression,
                      std::string* lhs,
                      std::string* rhs);
bool split_named_call(const std::string& expression,
                      const std::string& name,
                      std::string* inside);
bool split_named_call_with_arguments(const std::string& expression,
                                     const std::string& name,
                                     std::vector<std::string>* arguments);
std::vector<std::string> split_top_level_arguments(const std::string& text);
std::string expand_inline_function_commands(Calculator* calculator,
                                            const std::string& expression);
std::string stored_value_precise_decimal_text(const StoredValue& value);
PreciseDecimal parse_precise_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables);
std::string format_stored_value(const StoredValue& value, bool symbolic_constants_mode);
std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode);
std::string format_symbolic_scalar(double value);
double factorial_int(int n);
double factorial_value(long long n);
Rational factorial_rational(long long n);
double combination_value(long long n, long long r);
Rational combination_rational(long long n, long long r);
double permutation_value(long long n, long long r);
Rational permutation_rational(long long n, long long r);
std::string taylor_series_to_string(const std::vector<double>& coefficients,
                                    const std::string& variable_name,
                                    double center);
std::string shifted_series_base(const std::string& variable_name, double center);
std::string generalized_series_to_string(const std::vector<double>& coefficients,
                                         const std::string& variable_name,
                                         double center,
                                         int denominator);
std::vector<double> solve_dense_linear_system(std::vector<std::vector<double>> matrix,
                                              std::vector<double> rhs,
                                              const std::string& context);
bool is_reserved_function_name(const std::string& name);
bool split_function_definition(const std::string& expression,
                               std::string* function_name,
                               std::string* parameter_name,
                               std::string* body);
double parse_decimal_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function = {},
    InvokeScriptFunctionDecimalCallback invoke_script_function = {});
bool try_evaluate_matrix_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables,
    const std::map<std::string, CustomFunction>* functions,
    const HasScriptFunctionCallback& has_script_function,
    const InvokeScriptFunctionDecimalCallback& invoke_script_function,
    matrix::Value* value);
Rational parse_exact_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables,
    const std::map<std::string, CustomFunction>* functions,
    HasScriptFunctionCallback has_script_function = {});
bool try_base_conversion_expression(const std::string& expression,
                                    const std::map<std::string, StoredValue>* variables,
                                    const std::map<std::string, CustomFunction>* functions,
                                    const HexFormatOptions& hex_options,
                                    std::string* output);
std::string matrix_literal_expression(const matrix::Matrix& value);
bool try_symbolic_constant_expression(
    const std::string& expression,
    const std::map<std::string, StoredValue>* variables,
    const std::map<std::string, CustomFunction>* functions,
    std::string* output);
std::map<std::string, StoredValue> visible_variables(const Calculator::Impl* impl);
bool has_visible_script_function(const Calculator::Impl* impl, const std::string& name);
void assign_visible_variable(Calculator::Impl* impl,
                             const std::string& name,
                             const StoredValue& value);
StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode);
double invoke_script_function_decimal(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& name,
                                      const std::vector<double>& arguments);
ScriptSignal execute_script_statement(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const script::Statement& statement,
                                      bool exact_mode,
                                      std::string* last_output,
                                      bool create_scope);
std::string render_script_block(const script::BlockStatement& block, int indent);

#endif

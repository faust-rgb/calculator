#include "function_registry.h"

#include <algorithm>
#include <map>
#include <stdexcept>

namespace runtime {
namespace {

bool exact_sqrt_bigint(const numeric::BigInt& value, numeric::BigInt* root) {
    if (value.sign() < 0) {
        return false;
    }
    numeric::BigInt low(0);
    numeric::BigInt high = value + numeric::BigInt(1);
    while (low + numeric::BigInt(1) < high) {
        numeric::BigInt mid = (low + high);
        mid.div_uint32(2);
        const numeric::BigInt square = mid * mid;
        if (square <= value) {
            low = mid;
        } else {
            high = mid;
        }
    }
    if (low * low == value) {
        *root = low;
        return true;
    }
    return false;
}

int small_nonnegative_integer(const numeric::BigInt& value, const std::string& context) {
    if (value.sign() < 0) {
        throw std::runtime_error(context + " requires a non-negative integer");
    }
    const std::string text = value.to_string();
    if (text.size() > 6) {
        throw std::runtime_error(context + " argument is too large for the initial 2.0 evaluator");
    }
    int result = 0;
    for (char ch : text) {
        result = result * 10 + (ch - '0');
    }
    return result;
}

long long small_integer(const numeric::BigInt& value, const std::string& context) {
    const std::string text = value.to_string();
    if (text.size() > 18) {
        throw std::runtime_error(context + " argument is too large for the initial 2.0 evaluator");
    }
    long long result = 0;
    std::size_t pos = 0;
    bool negative = false;
    if (text[pos] == '-') {
        negative = true;
        ++pos;
    }
    for (; pos < text.size(); ++pos) {
        result = result * 10 + (text[pos] - '0');
    }
    return negative ? -result : result;
}

numeric::BigInt require_integer(const numeric::Number& value, const std::string& context) {
    if (!value.is_integer()) {
        throw std::runtime_error(context + " requires integer arguments");
    }
    return value.as_integer();
}

numeric::Number integer_div_value(const std::vector<numeric::Number>& args,
                                  const numeric::PrecisionContext&) {
    const numeric::BigInt lhs = require_integer(args[0], "idiv");
    const numeric::BigInt rhs = require_integer(args[1], "idiv");
    if (rhs.is_zero()) {
        throw std::runtime_error("division by zero");
    }
    return numeric::Number(lhs / rhs);
}

numeric::Number absolute_value(const std::vector<numeric::Number>& args,
                               const numeric::PrecisionContext&) {
    const numeric::Number& value = args[0];
    if (value.is_integer()) {
        return numeric::Number(value.as_integer().abs());
    }
    if (value.is_rational()) {
        return numeric::Number(numeric::Rational(value.as_rational().numerator().abs(),
                                                value.as_rational().denominator()));
    }
    if (value.is_decimal()) {
        return numeric::Number(numeric::BigDecimal(value.as_decimal().coefficient().abs(),
                                                  value.as_decimal().scale()));
    }
    throw std::runtime_error("abs for complex values is not migrated yet");
}

numeric::Number sign_value(const std::vector<numeric::Number>& args,
                           const numeric::PrecisionContext&) {
    const numeric::Number& value = args[0];
    if (value.is_integer()) {
        return numeric::Number(value.as_integer().sign());
    }
    if (value.is_rational()) {
        return numeric::Number(value.as_rational().numerator().sign());
    }
    if (value.is_decimal()) {
        return numeric::Number(value.as_decimal().coefficient().sign());
    }
    throw std::runtime_error("sign for complex values is not defined");
}

numeric::Number step_value(const std::vector<numeric::Number>& args,
                           const numeric::PrecisionContext& context) {
    return numeric::Number(args[0].compare(numeric::Number(0), context) >= 0 ? 1 : 0);
}

numeric::Number delta_value(const std::vector<numeric::Number>& args,
                            const numeric::PrecisionContext& context) {
    return numeric::Number(args[0].compare(numeric::Number(0), context) == 0 ? 1 : 0);
}

numeric::Number min_value(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext& context) {
    numeric::Number result = args[0];
    for (std::size_t i = 1; i < args.size(); ++i) {
        if (args[i].compare(result, context) < 0) {
            result = args[i];
        }
    }
    return result;
}

numeric::Number max_value(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext& context) {
    numeric::Number result = args[0];
    for (std::size_t i = 1; i < args.size(); ++i) {
        if (args[i].compare(result, context) > 0) {
            result = args[i];
        }
    }
    return result;
}

numeric::Number clamp_value(const std::vector<numeric::Number>& args,
                            const numeric::PrecisionContext& context) {
    numeric::Number lo = args[1];
    numeric::Number hi = args[2];
    if (lo.compare(hi, context) > 0) {
        const numeric::Number tmp = lo;
        lo = hi;
        hi = tmp;
    }
    if (args[0].compare(lo, context) < 0) {
        return lo;
    }
    if (args[0].compare(hi, context) > 0) {
        return hi;
    }
    return args[0];
}

numeric::Number sum_value(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext& context) {
    numeric::Number result(0);
    for (const numeric::Number& arg : args) {
        result = numeric::add(result, arg, context);
    }
    return result;
}

numeric::Number mean_value(const std::vector<numeric::Number>& args,
                           const numeric::PrecisionContext& context) {
    return numeric::divide(sum_value(args, context),
                           numeric::Number(static_cast<long long>(args.size())),
                           context);
}

numeric::BigInt floor_rational(const numeric::Rational& value) {
    numeric::BigInt quotient = value.numerator() / value.denominator();
    const numeric::BigInt remainder = value.numerator() % value.denominator();
    if (value.numerator().sign() < 0 && !remainder.is_zero()) {
        quotient = quotient - numeric::BigInt(1);
    }
    return quotient;
}

numeric::BigInt ceil_rational(const numeric::Rational& value) {
    numeric::BigInt quotient = value.numerator() / value.denominator();
    const numeric::BigInt remainder = value.numerator() % value.denominator();
    if (value.numerator().sign() > 0 && !remainder.is_zero()) {
        quotient = quotient + numeric::BigInt(1);
    }
    return quotient;
}

numeric::Rational rational_from_number(const numeric::Number& value) {
    if (value.is_integer()) {
        return numeric::Rational(value.as_integer());
    }
    if (value.is_rational()) {
        return value.as_rational();
    }
    if (value.is_decimal()) {
        return numeric::Rational(value.as_decimal().coefficient(),
                                 numeric::BigDecimal::pow10(value.as_decimal().scale()));
    }
    throw std::runtime_error("complex value is not a real scalar");
}

numeric::Number floor_value(const std::vector<numeric::Number>& args,
                            const numeric::PrecisionContext&) {
    return numeric::Number(floor_rational(rational_from_number(args[0])));
}

numeric::Number ceil_value(const std::vector<numeric::Number>& args,
                           const numeric::PrecisionContext&) {
    return numeric::Number(ceil_rational(rational_from_number(args[0])));
}

numeric::Number round_value(const std::vector<numeric::Number>& args,
                            const numeric::PrecisionContext&) {
    const numeric::Rational value = rational_from_number(args[0]);
    numeric::BigInt quotient = value.numerator() / value.denominator();
    const numeric::BigInt remainder = value.numerator() % value.denominator();
    if (remainder.is_zero()) {
        return numeric::Number(quotient);
    }
    if (remainder.abs() * numeric::BigInt(2) >= value.denominator()) {
        quotient = quotient + numeric::BigInt(value.numerator().sign() < 0 ? -1 : 1);
    }
    return numeric::Number(quotient);
}

numeric::Number trunc_value(const std::vector<numeric::Number>& args,
                            const numeric::PrecisionContext&) {
    const numeric::Rational value = rational_from_number(args[0]);
    return numeric::Number(value.numerator() / value.denominator());
}

numeric::Number median_value(std::vector<numeric::Number> args,
                             const numeric::PrecisionContext& context) {
    std::sort(args.begin(), args.end(), [&context](const numeric::Number& lhs,
                                                   const numeric::Number& rhs) {
        return lhs.compare(rhs, context) < 0;
    });
    const std::size_t mid = args.size() / 2;
    if (args.size() % 2 == 1) {
        return args[mid];
    }
    return numeric::divide(numeric::add(args[mid - 1], args[mid], context),
                           numeric::Number(2),
                           context);
}

numeric::Number variance_value(const std::vector<numeric::Number>& args,
                               const numeric::PrecisionContext& context) {
    const numeric::Number mean = mean_value(args, context);
    numeric::Number total(0);
    for (const numeric::Number& arg : args) {
        const numeric::Number delta = numeric::subtract(arg, mean, context);
        total = numeric::add(total, numeric::multiply(delta, delta, context), context);
    }
    return numeric::divide(total, numeric::Number(static_cast<long long>(args.size())), context);
}

numeric::Number mode_value(const std::vector<numeric::Number>& args,
                           const numeric::PrecisionContext& context) {
    numeric::Number best = args[0];
    std::size_t best_count = 0;
    for (const numeric::Number& candidate : args) {
        std::size_t count = 0;
        for (const numeric::Number& arg : args) {
            if (arg.compare(candidate, context) == 0) {
                ++count;
            }
        }
        if (count > best_count ||
            (count == best_count && candidate.compare(best, context) < 0)) {
            best = candidate;
            best_count = count;
        }
    }
    return best;
}

numeric::Number mod_value(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext&) {
    const numeric::BigInt lhs = require_integer(args[0], "mod");
    const numeric::BigInt rhs = require_integer(args[1], "mod");
    if (rhs.is_zero()) {
        throw std::runtime_error("division by zero");
    }
    return numeric::Number(lhs % rhs);
}

numeric::Number gcd_value(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext&) {
    numeric::BigInt result = require_integer(args[0], "gcd").abs();
    for (std::size_t i = 1; i < args.size(); ++i) {
        result = numeric::BigInt::gcd(result, require_integer(args[i], "gcd"));
    }
    return numeric::Number(result);
}

numeric::Number lcm_value(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext&) {
    numeric::BigInt result = require_integer(args[0], "lcm").abs();
    for (std::size_t i = 1; i < args.size(); ++i) {
        const numeric::BigInt next = require_integer(args[i], "lcm").abs();
        if (result.is_zero() || next.is_zero()) {
            result = numeric::BigInt(0);
        } else {
            result = (result / numeric::BigInt::gcd(result, next)) * next;
        }
    }
    return numeric::Number(result);
}

numeric::Number fibonacci_value(const std::vector<numeric::Number>& args,
                                const numeric::PrecisionContext&) {
    const int n = small_nonnegative_integer(require_integer(args[0], "fib"), "fib");
    numeric::BigInt a(0);
    numeric::BigInt b(1);
    for (int i = 0; i < n; ++i) {
        const numeric::BigInt next = a + b;
        a = b;
        b = next;
    }
    return numeric::Number(a);
}

bool is_prime_small(long long value) {
    if (value < 2) {
        return false;
    }
    if (value == 2) {
        return true;
    }
    if (value % 2 == 0) {
        return false;
    }
    for (long long divisor = 3; divisor <= value / divisor; divisor += 2) {
        if (value % divisor == 0) {
            return false;
        }
    }
    return true;
}

numeric::Number is_prime_value(const std::vector<numeric::Number>& args,
                               const numeric::PrecisionContext&) {
    const long long n = small_integer(require_integer(args[0], "is_prime"), "is_prime");
    return numeric::Number(is_prime_small(n) ? 1 : 0);
}

numeric::Number next_prime_value(const std::vector<numeric::Number>& args,
                                 const numeric::PrecisionContext&) {
    long long n = small_integer(require_integer(args[0], "next_prime"), "next_prime") + 1;
    while (!is_prime_small(n)) {
        ++n;
    }
    return numeric::Number(n);
}

numeric::Number prev_prime_value(const std::vector<numeric::Number>& args,
                                 const numeric::PrecisionContext&) {
    long long n = small_integer(require_integer(args[0], "prev_prime"), "prev_prime") - 1;
    while (n >= 2 && !is_prime_small(n)) {
        --n;
    }
    return numeric::Number(n >= 2 ? n : 0);
}

numeric::Number prime_pi_value(const std::vector<numeric::Number>& args,
                               const numeric::PrecisionContext&) {
    const int n = small_nonnegative_integer(require_integer(args[0], "prime_pi"), "prime_pi");
    int count = 0;
    for (int value = 2; value <= n; ++value) {
        if (is_prime_small(value)) {
            ++count;
        }
    }
    return numeric::Number(count);
}

numeric::Number euler_phi_value(const std::vector<numeric::Number>& args,
                                const numeric::PrecisionContext&) {
    long long n = small_integer(require_integer(args[0], "euler_phi"), "euler_phi");
    if (n <= 0) {
        throw std::runtime_error("euler_phi requires a positive integer");
    }
    long long result = n;
    for (long long p = 2; p <= n / p; ++p) {
        if (n % p == 0) {
            while (n % p == 0) {
                n /= p;
            }
            result -= result / p;
        }
    }
    if (n > 1) {
        result -= result / n;
    }
    return numeric::Number(result);
}

numeric::BigInt factorial_bigint(int n) {
    numeric::BigInt result(1);
    for (int i = 2; i <= n; ++i) {
        result = result * numeric::BigInt(i);
    }
    return result;
}

numeric::Number factorial_value(const std::vector<numeric::Number>& args,
                                const numeric::PrecisionContext&) {
    const int n = small_nonnegative_integer(require_integer(args[0], "factorial"),
                                            "factorial");
    return numeric::Number(factorial_bigint(n));
}

numeric::Number combination_value(const std::vector<numeric::Number>& args,
                                  const numeric::PrecisionContext&) {
    const int n = small_nonnegative_integer(require_integer(args[0], "nCr"), "nCr");
    const int r = small_nonnegative_integer(require_integer(args[1], "nCr"), "nCr");
    if (r > n) {
        return numeric::Number(0);
    }
    const int k = r < n - r ? r : n - r;
    numeric::BigInt result(1);
    for (int i = 1; i <= k; ++i) {
        result = (result * numeric::BigInt(n - k + i)) / numeric::BigInt(i);
    }
    return numeric::Number(result);
}

numeric::Number permutation_value(const std::vector<numeric::Number>& args,
                                  const numeric::PrecisionContext&) {
    const int n = small_nonnegative_integer(require_integer(args[0], "nPr"), "nPr");
    const int r = small_nonnegative_integer(require_integer(args[1], "nPr"), "nPr");
    if (r > n) {
        return numeric::Number(0);
    }
    numeric::BigInt result(1);
    for (int i = 0; i < r; ++i) {
        result = result * numeric::BigInt(n - i);
    }
    return numeric::Number(result);
}

numeric::Number temperature_affine(const numeric::Number& value,
                                   long long numerator,
                                   long long denominator,
                                   long long offset,
                                   const numeric::PrecisionContext& context) {
    return numeric::add(numeric::divide(numeric::multiply(value, numeric::Number(numerator), context),
                                        numeric::Number(denominator),
                                        context),
                        numeric::Number(offset),
                        context);
}

numeric::Number celsius_value(const std::vector<numeric::Number>& args,
                              const numeric::PrecisionContext& context) {
    return temperature_affine(numeric::subtract(args[0], numeric::Number(32), context),
                              5,
                              9,
                              0,
                              context);
}

numeric::Number fahrenheit_value(const std::vector<numeric::Number>& args,
                                 const numeric::PrecisionContext& context) {
    return temperature_affine(args[0], 9, 5, 32, context);
}

numeric::Number kelvin_value(const std::vector<numeric::Number>& args,
                             const numeric::PrecisionContext& context) {
    return numeric::add(args[0], numeric::Number(numeric::Rational(5463, 20)), context);
}

numeric::Number pow_integer(const numeric::Number& base,
                            const numeric::Number& exponent,
                            const numeric::PrecisionContext& context) {
    if (!exponent.is_integer()) {
        throw std::runtime_error("pow currently requires an integer exponent");
    }
    const std::string exponent_text = exponent.as_integer().to_string();
    if (exponent_text.size() > 8) {
        throw std::runtime_error("pow exponent is too large for the initial 2.0 evaluator");
    }
    int exp = 0;
    std::size_t pos = 0;
    bool negative = false;
    if (exponent_text[pos] == '-') {
        negative = true;
        ++pos;
    }
    for (; pos < exponent_text.size(); ++pos) {
        exp = exp * 10 + (exponent_text[pos] - '0');
    }
    if (negative) {
        exp = -exp;
    }

    numeric::Number result(1);
    const int count = exp < 0 ? -exp : exp;
    for (int i = 0; i < count; ++i) {
        result = numeric::multiply(result, base, context);
    }
    if (exp < 0) {
        result = numeric::divide(numeric::Number(1), result, context);
    }
    return result;
}

numeric::Number sqrt_exact(const std::vector<numeric::Number>& args,
                           const numeric::PrecisionContext&) {
    const numeric::Number& value = args[0];
    if (value.is_integer()) {
        numeric::BigInt root;
        if (exact_sqrt_bigint(value.as_integer(), &root)) {
            return numeric::Number(root);
        }
    }
    if (value.is_rational()) {
        numeric::BigInt numerator_root;
        numeric::BigInt denominator_root;
        if (exact_sqrt_bigint(value.as_rational().numerator(), &numerator_root) &&
            exact_sqrt_bigint(value.as_rational().denominator(), &denominator_root)) {
            return numeric::Number(numeric::Rational(numerator_root, denominator_root));
        }
    }
    throw std::runtime_error("sqrt has no exact numeric result for this value yet");
}

bool number_text_is(const numeric::Number& value, const std::string& text) {
    return value.to_string() == text;
}

numeric::Number sin_exact(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext&) {
    if (number_text_is(args[0], "0")) {
        return numeric::Number(0);
    }
    throw std::runtime_error("sin has no exact numeric result for this value yet");
}

numeric::Number cos_exact(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext&) {
    if (number_text_is(args[0], "0")) {
        return numeric::Number(1);
    }
    throw std::runtime_error("cos has no exact numeric result for this value yet");
}

numeric::Number tan_exact(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext&) {
    if (number_text_is(args[0], "0")) {
        return numeric::Number(0);
    }
    throw std::runtime_error("tan has no exact numeric result for this value yet");
}

numeric::Number exp_exact(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext&) {
    if (number_text_is(args[0], "0")) {
        return numeric::Number(1);
    }
    throw std::runtime_error("exp has no exact numeric result for this value yet");
}

numeric::Number ln_exact(const std::vector<numeric::Number>& args,
                         const numeric::PrecisionContext&) {
    if (number_text_is(args[0], "1")) {
        return numeric::Number(0);
    }
    throw std::runtime_error("ln has no exact numeric result for this value yet");
}

numeric::Number sinh_exact(const std::vector<numeric::Number>& args,
                           const numeric::PrecisionContext&) {
    if (number_text_is(args[0], "0")) {
        return numeric::Number(0);
    }
    throw std::runtime_error("sinh has no exact numeric result for this value yet");
}

numeric::Number cosh_exact(const std::vector<numeric::Number>& args,
                           const numeric::PrecisionContext&) {
    if (number_text_is(args[0], "0")) {
        return numeric::Number(1);
    }
    throw std::runtime_error("cosh has no exact numeric result for this value yet");
}

numeric::Number tanh_exact(const std::vector<numeric::Number>& args,
                           const numeric::PrecisionContext&) {
    if (number_text_is(args[0], "0")) {
        return numeric::Number(0);
    }
    throw std::runtime_error("tanh has no exact numeric result for this value yet");
}

numeric::Number inverse_zero_exact(const std::vector<numeric::Number>& args,
                                   const numeric::PrecisionContext&) {
    if (number_text_is(args[0], "0")) {
        return numeric::Number(0);
    }
    throw std::runtime_error("inverse function has no exact numeric result for this value yet");
}

numeric::Number real_part(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext& context) {
    if (!args[0].is_complex()) {
        return args[0];
    }
    return numeric::Number(args[0].to_complex(context).real());
}

numeric::Number imag_part(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext& context) {
    if (!args[0].is_complex()) {
        return numeric::Number(0);
    }
    return numeric::Number(args[0].to_complex(context).imag());
}

numeric::Number conjugate(const std::vector<numeric::Number>& args,
                          const numeric::PrecisionContext& context) {
    if (!args[0].is_complex()) {
        return args[0];
    }
    return numeric::Number(numeric::conj(args[0].to_complex(context)));
}

FunctionRegistry make_builtins() {
    FunctionRegistry registry;
    registry.register_function({"sqrt", {}, "elementary", "Exact square root when available", {"sqrt(x)"}, "x'/(2*sqrt(x))", "", 1, 1, sqrt_exact});
    registry.register_function({"sin", {}, "elementary", "Exact sine for known values", {"sin(x)"}, "cos(x)*x'", "-cos(x)", 1, 1, sin_exact});
    registry.register_function({"cos", {}, "elementary", "Exact cosine for known values", {"cos(x)"}, "-sin(x)*x'", "sin(x)", 1, 1, cos_exact});
    registry.register_function({"tan", {}, "elementary", "Exact tangent for known values", {"tan(x)"}, "sec(x)^2*x'", "-ln(cos(x))", 1, 1, tan_exact});
    registry.register_function({"exp", {}, "elementary", "Exact exponential for known values", {"exp(x)"}, "exp(x)*x'", "exp(x)", 1, 1, exp_exact});
    registry.register_function({"ln", {"log"}, "elementary", "Exact natural logarithm for known values", {"ln(x)"}, "x'/x", "x*ln(x)-x", 1, 1, ln_exact});
    registry.register_function({"sinh", {}, "elementary", "Exact hyperbolic sine for known values", {"sinh(x)"}, "cosh(x)*x'", "cosh(x)", 1, 1, sinh_exact});
    registry.register_function({"cosh", {}, "elementary", "Exact hyperbolic cosine for known values", {"cosh(x)"}, "sinh(x)*x'", "sinh(x)", 1, 1, cosh_exact});
    registry.register_function({"tanh", {}, "elementary", "Exact hyperbolic tangent for known values", {"tanh(x)"}, "sech(x)^2*x'", "", 1, 1, tanh_exact});
    registry.register_function({"asin", {}, "elementary", "Exact inverse sine for known values", {"asin(x)"}, "x'/sqrt(1-x^2)", "", 1, 1, inverse_zero_exact});
    registry.register_function({"atan", {}, "elementary", "Exact inverse tangent for known values", {"atan(x)"}, "x'/(1+x^2)", "", 1, 1, inverse_zero_exact});
    registry.register_function({"asinh", {}, "elementary", "Exact inverse hyperbolic sine for known values", {"asinh(x)"}, "x'/sqrt(x^2+1)", "", 1, 1, inverse_zero_exact});
    registry.register_function({"atanh", {}, "elementary", "Exact inverse hyperbolic tangent for known values", {"atanh(x)"}, "x'/(1-x^2)", "", 1, 1, inverse_zero_exact});
    registry.register_function({"abs", {}, "arithmetic", "Absolute value for real exact numbers", {"abs(x)"}, "sign(x)*x'", "", 1, 1, absolute_value});
    registry.register_function({"sign", {}, "arithmetic", "Sign of a real exact number", {"sign(x)"}, "0 away from discontinuities", "", 1, 1, sign_value});
    registry.register_function({"step", {"heaviside"}, "arithmetic", "Unit step", {"step(x)"}, "delta(x)*x'", "", 1, 1, step_value});
    registry.register_function({"delta", {"impulse"}, "arithmetic", "Discrete impulse shorthand", {"delta(x)"}, "", "", 1, 1, delta_value});
    registry.register_function({"floor", {}, "arithmetic", "Exact floor for rational/decimal values", {"floor(x)"}, "", "", 1, 1, floor_value});
    registry.register_function({"ceil", {}, "arithmetic", "Exact ceiling for rational/decimal values", {"ceil(x)"}, "", "", 1, 1, ceil_value});
    registry.register_function({"round", {}, "arithmetic", "Round half away from zero", {"round(x)"}, "", "", 1, 1, round_value});
    registry.register_function({"trunc", {}, "arithmetic", "Truncate toward zero", {"trunc(x)"}, "", "", 1, 1, trunc_value});
    registry.register_function({"min", {}, "arithmetic", "Minimum", {"min(a, b)"}, "", "", 2, 64, min_value});
    registry.register_function({"max", {}, "arithmetic", "Maximum", {"max(a, b)"}, "", "", 2, 64, max_value});
    registry.register_function({"clamp", {}, "arithmetic", "Clamp into a range", {"clamp(x, lo, hi)"}, "", "", 3, 3, clamp_value});
    registry.register_function({"sum", {}, "aggregate", "Exact aggregate sum", {"sum(a, b, ...)"}, "", "", 1, 256, sum_value});
    registry.register_function({"mean", {"avg"}, "aggregate", "Exact aggregate mean", {"mean(a, b, ...)"}, "", "", 1, 256, mean_value});
    registry.register_function({"median", {}, "statistics", "Exact median", {"median(a, b, ...)"}, "", "", 1, 256, median_value});
    registry.register_function({"mode", {}, "statistics", "Mode with lowest-value tie break", {"mode(a, b, ...)"}, "", "", 1, 256, mode_value});
    registry.register_function({"var", {"variance"}, "statistics", "Population variance", {"var(a, b, ...)"}, "", "", 1, 256, variance_value});
    registry.register_function({"idiv", {}, "integer", "Integer quotient", {"idiv(a, b)"}, "", "", 2, 2, integer_div_value});
    registry.register_function({"mod", {}, "integer", "Integer remainder", {"mod(a, b)"}, "", "", 2, 2, mod_value});
    registry.register_function({"gcd", {}, "integer", "Greatest common divisor", {"gcd(a, b)"}, "", "", 2, 64, gcd_value});
    registry.register_function({"lcm", {}, "integer", "Least common multiple", {"lcm(a, b)"}, "", "", 2, 64, lcm_value});
    registry.register_function({"factorial", {}, "integer", "Exact factorial", {"factorial(n)"}, "", "", 1, 1, factorial_value});
    registry.register_function({"nCr", {"binom"}, "integer", "Exact combinations", {"nCr(n, r)"}, "", "", 2, 2, combination_value});
    registry.register_function({"nPr", {}, "integer", "Exact permutations", {"nPr(n, r)"}, "", "", 2, 2, permutation_value});
    registry.register_function({"fib", {"fibonacci"}, "integer", "Exact Fibonacci number", {"fib(n)"}, "", "", 1, 1, fibonacci_value});
    registry.register_function({"is_prime", {}, "integer", "Primality test for migrated small integers", {"is_prime(n)"}, "", "", 1, 1, is_prime_value});
    registry.register_function({"next_prime", {}, "integer", "Next prime after n", {"next_prime(n)"}, "", "", 1, 1, next_prime_value});
    registry.register_function({"prev_prime", {}, "integer", "Previous prime before n", {"prev_prime(n)"}, "", "", 1, 1, prev_prime_value});
    registry.register_function({"prime_pi", {}, "integer", "Count primes <= n", {"prime_pi(n)"}, "", "", 1, 1, prime_pi_value});
    registry.register_function({"euler_phi", {"phi"}, "integer", "Euler totient", {"euler_phi(n)"}, "", "", 1, 1, euler_phi_value});
    registry.register_function({"celsius", {"f2c"}, "unit", "Fahrenheit to Celsius", {"celsius(f)"}, "", "", 1, 1, celsius_value});
    registry.register_function({"fahrenheit", {"c2f"}, "unit", "Celsius to Fahrenheit", {"fahrenheit(c)"}, "", "", 1, 1, fahrenheit_value});
    registry.register_function({"kelvin", {}, "unit", "Celsius to Kelvin", {"kelvin(c)"}, "", "", 1, 1, kelvin_value});
    registry.register_function({
        "pow",
        {},
        "elementary",
        "Exact integer power",
        {"pow(x, n)"},
        "x^n*(n'*ln(x)+n*x'/x)",
        "",
        2,
        2,
        [](const std::vector<numeric::Number>& args,
           const numeric::PrecisionContext& context) {
            return pow_integer(args[0], args[1], context);
        },
    });
    registry.register_function({"real", {}, "complex", "Real part", {"real(z)"}, "", "", 1, 1, real_part});
    registry.register_function({"imag", {}, "complex", "Imaginary part", {"imag(z)"}, "", "", 1, 1, imag_part});
    registry.register_function({"conj", {}, "complex", "Complex conjugate", {"conj(z)"}, "", "", 1, 1, conjugate});
    return registry;
}

}  // namespace

const FunctionRegistry& FunctionRegistry::builtins() {
    static const FunctionRegistry registry = make_builtins();
    return registry;
}

void FunctionRegistry::register_function(const FunctionSpec& spec) {
    functions_[spec.name] = spec;
    for (const std::string& alias : spec.aliases) {
        aliases_[alias] = spec.name;
    }
}

const FunctionSpec* FunctionRegistry::find(const std::string& name) const {
    auto it = functions_.find(name);
    if (it == functions_.end()) {
        const auto alias_it = aliases_.find(name);
        if (alias_it != aliases_.end()) {
            it = functions_.find(alias_it->second);
        }
    }
    if (it == functions_.end()) {
        return nullptr;
    }
    return &it->second;
}

std::vector<FunctionSpec> FunctionRegistry::functions() const {
    std::vector<FunctionSpec> result;
    result.reserve(functions_.size());
    for (const auto& [_, spec] : functions_) {
        result.push_back(spec);
    }
    return result;
}

numeric::Number FunctionRegistry::evaluate_numeric(
    const std::string& name,
    const std::vector<numeric::Number>& args,
    const numeric::PrecisionContext& context) const {
    const FunctionSpec* spec = find(name);
    if (spec == nullptr || !spec->numeric_eval) {
        throw std::runtime_error("unknown numeric function: " + name);
    }
    if (args.size() < spec->min_arity || args.size() > spec->max_arity) {
        throw std::runtime_error(name + " received the wrong number of arguments");
    }
    return spec->numeric_eval(args, context);
}

}  // namespace runtime

#include "functions.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace numeric {
namespace {

PrecisionContext work_context(const PrecisionContext& context, int guard = 12) {
    PrecisionContext result = context;
    result.digits = std::max(context.digits + guard, context.digits);
    result.max_iterations = std::max(context.max_iterations, 1000);
    return result;
}

Number integer(long long value) {
    return Number(value);
}

Number decimal_power_of_ten(int exponent) {
    if (exponent >= 0) {
        return Number(BigDecimal(BigDecimal::pow10(exponent), 0));
    }
    return Number(BigDecimal(BigInt(1), -exponent));
}

Number tolerance(const PrecisionContext& context) {
    return decimal_power_of_ten(-(context.digits + 4));
}

Number round_to_context(const Number& value) {
    const PrecisionContext& context = value.context();
    if (value.is_complex()) {
        return value;
    }
    if (abs(value).compare(tolerance(context)) <= 0) {
        return Number(0).with_context(context);
    }
    BigDecimal decimal = value.to_decimal();
    if (decimal.scale() <= context.digits) {
        return Number(decimal).with_context(context);
    }
    const int drop = decimal.scale() - context.digits;
    BigInt divisor = BigDecimal::pow10(drop);
    BigInt quotient = decimal.coefficient() / divisor;
    BigInt remainder = decimal.coefficient().abs() % divisor;
    if (remainder * BigInt(2) >= divisor) {
        quotient = quotient + BigInt(decimal.coefficient().sign() < 0 ? -1 : 1);
    }
    return Number(BigDecimal(quotient, context.digits)).with_context(context);
}

bool exact_nth_root_bigint(const BigInt& value, long long degree, BigInt* root) {
    if (degree <= 0 || value.sign() < 0) {
        return false;
    }
    BigInt low(0);
    BigInt high = value + BigInt(1);
    while (low + BigInt(1) < high) {
        BigInt mid = low + high;
        mid.div_uint32(2);
        const BigInt powered = BigInt::pow(mid, static_cast<unsigned int>(degree));
        if (powered <= value) {
            low = mid;
        } else {
            high = mid;
        }
    }
    if (BigInt::pow(low, static_cast<unsigned int>(degree)) == value) {
        *root = low;
        return true;
    }
    return false;
}

bool less_abs_than(const Number& value, const Number& threshold) {
    return abs(value).compare(threshold) <= 0;
}

long long small_integer_value(const Number& value, const std::string& name) {
    if (!value.is_integer()) {
        throw std::runtime_error(name + " requires an integer argument");
    }
    const std::string text = value.as_integer().to_string();
    if (text.size() > 9 || text == "-2147483648") {
        throw std::runtime_error(name + " integer argument is too large");
    }
    long long result = 0;
    std::size_t pos = 0;
    bool negative = false;
    if (!text.empty() && text[0] == '-') {
        negative = true;
        pos = 1;
    }
    for (; pos < text.size(); ++pos) {
        result = result * 10 + (text[pos] - '0');
    }
    return negative ? -result : result;
}

Number atan_series(Number x) {
    const PrecisionContext& context = x.context();
    PrecisionContext work = work_context(context);
    const Number x2 = multiply(x, x).with_context(work);
    Number term = x;
    term.set_context(work);
    Number sum = term;
    const Number eps = tolerance(work);
    for (int n = 1; n < work.max_iterations; ++n) {
        term = multiply(term, x2).with_context(work);
        const Number addend = divide(term, integer(2 * n + 1)).with_context(work);
        sum = (n % 2 == 0) ? add(sum, addend).with_context(work)
                           : subtract(sum, addend).with_context(work);
        if (less_abs_than(addend, eps)) {
            break;
        }
    }
    return round_to_context(sum.with_context(context));
}

Number reduce_radians(Number value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number two_pi = multiply(integer(2), pi()).with_context(work);
    const Number pi_value = pi().with_context(work);
    int guard = 0;
    while (value.compare(pi_value) > 0 && guard++ < 10000) {
        value = subtract(value, two_pi).with_context(work);
    }
    guard = 0;
    const Number neg_pi = negate(pi_value);
    while (value.compare(neg_pi) < 0 && guard++ < 10000) {
        value = add(value, two_pi).with_context(work);
    }
    return value.with_context(context);
}

Number exp_small(Number value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    Number term(1);
    term.set_context(work);
    Number sum(1);
    sum.set_context(work);
    const Number eps = tolerance(work);
    for (int n = 1; n < work.max_iterations; ++n) {
        term = divide(multiply(term, value).with_context(work), integer(n)).with_context(work);
        sum = add(sum, term).with_context(work);
        if (less_abs_than(term, eps)) {
            break;
        }
    }
    return round_to_context(sum.with_context(context));
}

Number ln2() {
    PrecisionContext work = work_context(default_precision());
    const Number z = divide(integer(1), integer(3)).with_context(work);
    const Number z2 = multiply(z, z).with_context(work);
    Number term = z;
    term.set_context(work);
    Number sum = term;
    sum.set_context(work);
    const Number eps = tolerance(work);
    for (int n = 1; n < work.max_iterations; ++n) {
        term = multiply(term, z2).with_context(work);
        const Number addend = divide(term, integer(2 * n + 1)).with_context(work);
        sum = add(sum, addend).with_context(work);
        if (less_abs_than(addend, eps)) {
            break;
        }
    }
    return round_to_context(multiply(integer(2), sum).with_context(work));
}

}  // namespace

bool is_zero(const Number& value) {
    if (value.is_integer()) {
        return value.as_integer().is_zero();
    }
    if (value.is_rational()) {
        return value.as_rational().numerator().is_zero();
    }
    if (value.is_decimal()) {
        return value.as_decimal().coefficient().is_zero();
    }
    return value.as_complex().real().coefficient().is_zero() &&
           value.as_complex().imag().coefficient().is_zero();
}

Number abs(const Number& value) {
    if (value.is_integer()) {
        return Number(value.as_integer().abs()).with_context(value.context());
    }
    if (value.is_rational()) {
        return Number(Rational(value.as_rational().numerator().abs(),
                               value.as_rational().denominator())).with_context(value.context());
    }
    if (value.is_decimal()) {
        return Number(BigDecimal(value.as_decimal().coefficient().abs(),
                                 value.as_decimal().scale())).with_context(value.context());
    }
    throw std::runtime_error("complex abs is not implemented in arbitrary precision scalar functions");
}

Number negate(const Number& value) {
    return subtract(integer(0), value);
}

Number sqrt(const Number& value) {
    const PrecisionContext& context = value.context();
    if (value.compare(integer(0)) < 0) {
        throw std::runtime_error("sqrt requires a non-negative real argument");
    }
    if (is_zero(value)) {
        return integer(0);
    }
    PrecisionContext work = work_context(context);
    Number x = value.compare(integer(1)) > 0 ? value : integer(1);
    x.set_context(work);
    const Number eps = tolerance(context);
    for (int i = 0; i < work.max_iterations; ++i) {
        const Number next = divide(add(x, divide(value, x).with_context(work)).with_context(work), integer(2)).with_context(work);
        if (less_abs_than(subtract(next, x).with_context(work), eps)) {
            return round_to_context(next.with_context(context));
        }
        x = next;
    }
    return round_to_context(x.with_context(context));
}

Number cbrt(const Number& value) {
    return root(value, integer(3));
}

Number root(const Number& value, const Number& degree) {
    const PrecisionContext& context = value.context();
    const long long n = small_integer_value(degree, "root");
    if (n == 0) {
        throw std::runtime_error("root degree cannot be zero");
    }
    if (n < 0) {
        return divide(integer(1), root(value, integer(-n)));
    }
    if (value.compare(integer(0)) < 0 && n % 2 == 0) {
        throw std::runtime_error("even root of negative value is not a real scalar");
    }
    if (is_zero(value)) {
        return integer(0);
    }
    if (value.is_integer()) {
        const bool negative_integer = value.as_integer().sign() < 0;
        BigInt root_value;
        if (exact_nth_root_bigint(value.as_integer().abs(), n, &root_value)) {
            if (negative_integer) {
                root_value = -root_value;
            }
            return Number(root_value).with_context(context);
        }
    }
    if (value.is_rational()) {
        const bool negative_rational = value.as_rational().numerator().sign() < 0;
        BigInt numerator_root;
        BigInt denominator_root;
        if (exact_nth_root_bigint(value.as_rational().numerator().abs(), n, &numerator_root) &&
            exact_nth_root_bigint(value.as_rational().denominator(), n, &denominator_root)) {
            if (negative_rational) {
                numerator_root = -numerator_root;
            }
            return Number(Rational(numerator_root, denominator_root)).with_context(context);
        }
    }
    PrecisionContext work = work_context(context);
    const bool negative = value.compare(integer(0)) < 0;
    const Number magnitude_exact = negative ? negate(value) : value;
    const Number magnitude = Number(magnitude_exact.to_decimal()).with_context(work);
    Number x = magnitude.compare(integer(1)) > 0 ? magnitude : integer(1);
    x.set_context(work);
    const Number n_number(n);
    const Number n_minus_one(n - 1);
    const Number eps = tolerance(context);
    for (int i = 0; i < work.max_iterations; ++i) {
        Number denominator(1);
        for (long long p = 0; p < n - 1; ++p) {
            denominator = multiply(denominator, x).with_context(work);
        }
        const Number next = divide(add(multiply(n_minus_one, x).with_context(work),
                                       divide(magnitude, denominator).with_context(work)).with_context(work),
                                   n_number).with_context(work);
        if (less_abs_than(subtract(next, x).with_context(work), eps)) {
            return negative ? negate(round_to_context(next.with_context(context)))
                            : round_to_context(next.with_context(context));
        }
        x = next;
    }
    return negative ? negate(round_to_context(x.with_context(context)))
                    : round_to_context(x.with_context(context));
}

Number pow(const Number& base, const Number& exponent) {
    const PrecisionContext& context = base.context();
    if (exponent.is_integer()) {
        long long n = small_integer_value(exponent, "pow");
        bool negative = n < 0;
        if (negative) {
            n = -n;
        }
        Number result(1);
        result.set_context(context);
        Number factor = base;
        factor.set_context(context);
        while (n > 0) {
            if ((n & 1LL) != 0) {
                result = multiply(result, factor);
            }
            factor = multiply(factor, factor);
            n >>= 1;
        }
        if (negative) {
            result = divide(integer(1), result);
        }
        return round_to_context(result);
    }
    return exp(multiply(exponent, ln(base)));
}

Number pi() {
    PrecisionContext work = work_context(default_precision(), 18);
    const Number atan_1_5 = atan_series(divide(integer(1), integer(5)).with_context(work));
    const Number atan_1_239 = atan_series(divide(integer(1), integer(239)).with_context(work));
    return round_to_context(multiply(integer(4),
                                      subtract(multiply(integer(4), atan_1_5).with_context(work),
                                               atan_1_239).with_context(work)).with_context(work));
}

Number exp(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    Number x = value;
    x.set_context(work);
    int halvings = 0;
    const Number one(1);
    while (abs(x).compare(one) > 0 && halvings < 64) {
        x = divide(x, integer(2)).with_context(work);
        ++halvings;
    }
    Number result = exp_small(x);
    result.set_context(work);
    for (int i = 0; i < halvings; ++i) {
        result = multiply(result, result);
    }
    return round_to_context(result.with_context(context));
}

Number ln(const Number& value) {
    const PrecisionContext& context = value.context();
    if (value.compare(integer(0)) <= 0) {
        throw std::runtime_error("ln requires a positive real argument");
    }
    PrecisionContext work = work_context(context);
    Number x = value;
    x.set_context(work);
    int twos = 0;
    const Number low = divide(integer(3), integer(4)).with_context(work);
    const Number high = divide(integer(3), integer(2)).with_context(work);
    while (x.compare(high) > 0 && twos < 10000) {
        x = divide(x, integer(2)).with_context(work);
        ++twos;
    }
    while (x.compare(low) < 0 && twos > -10000) {
        x = multiply(x, integer(2)).with_context(work);
        --twos;
    }
    const Number z = divide(subtract(x, integer(1)).with_context(work),
                            add(x, integer(1)).with_context(work)).with_context(work);
    const Number z2 = multiply(z, z).with_context(work);
    Number term = z;
    term.set_context(work);
    Number sum = term;
    sum.set_context(work);
    const Number eps = tolerance(context);
    for (int n = 1; n < work.max_iterations; ++n) {
        term = multiply(term, z2).with_context(work);
        const Number addend = divide(term, integer(2 * n + 1)).with_context(work);
        sum = add(sum, addend).with_context(work);
        if (less_abs_than(addend, eps)) {
            break;
        }
    }
    Number result = multiply(integer(2), sum).with_context(work);
    if (twos != 0) {
        result = add(result, multiply(integer(twos), ln2()).with_context(work)).with_context(work);
    }
    return round_to_context(result.with_context(context));
}

Number sin(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number x = reduce_radians(Number(value.to_decimal()).with_context(work));
    const Number x2 = multiply(x, x).with_context(work);
    Number term = x;
    term.set_context(work);
    Number sum = term;
    sum.set_context(work);
    const Number eps = tolerance(context);
    for (int n = 1; n < work.max_iterations; ++n) {
        term = divide(multiply(term, x2).with_context(work), integer((2 * n) * (2 * n + 1))).with_context(work);
        sum = (n % 2 == 0) ? add(sum, term).with_context(work) : subtract(sum, term).with_context(work);
        if (less_abs_than(term, eps)) {
            break;
        }
    }
    return round_to_context(sum.with_context(context));
}

Number cos(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number x = reduce_radians(Number(value.to_decimal()).with_context(work));
    const Number x2 = multiply(x, x).with_context(work);
    Number term(1);
    term.set_context(work);
    Number sum(1);
    sum.set_context(work);
    const Number eps = tolerance(context);
    for (int n = 1; n < work.max_iterations; ++n) {
        term = divide(multiply(term, x2).with_context(work), integer((2 * n - 1) * (2 * n))).with_context(work);
        sum = (n % 2 == 0) ? add(sum, term).with_context(work) : subtract(sum, term).with_context(work);
        if (less_abs_than(term, eps)) {
            break;
        }
    }
    return round_to_context(sum.with_context(context));
}

Number tan(const Number& value) {
    return divide(sin(value), cos(value));
}

Number atan(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    if (abs(value).compare(integer(1)) <= 0) {
        return atan_series(value);
    }
    const Number half_pi = divide(pi(), integer(2));
    const Number reduced = atan_series(divide(integer(1), abs(value)).with_context(work));
    Number result = subtract(half_pi, reduced).with_context(work);
    if (value.compare(integer(0)) < 0) {
        result = negate(result);
    }
    return round_to_context(result.with_context(context));
}

Number asin(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number one(1);
    const Number inside = subtract(one, multiply(value, value)).with_context(work);
    return atan(divide(value, sqrt(inside)).with_context(work));
}

Number acos(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    return round_to_context(subtract(divide(pi(), integer(2)).with_context(work),
                                     asin(value).with_context(work)).with_context(context));
}

Number sinh(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number e = exp(value).with_context(work);
    const Number inv = divide(integer(1), e).with_context(work);
    return round_to_context(divide(subtract(e, inv).with_context(work), integer(2)).with_context(context));
}

Number cosh(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number e = exp(value).with_context(work);
    const Number inv = divide(integer(1), e).with_context(work);
    return round_to_context(divide(add(e, inv).with_context(work), integer(2)).with_context(context));
}

Number tanh(const Number& value) {
    return divide(sinh(value), cosh(value));
}

Number asinh(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number x2 = multiply(value, value).with_context(work);
    const Number inner = add(integer(1), x2).with_context(work);
    return ln(add(value, sqrt(inner)).with_context(work));
}

Number acosh(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    if (value.compare(integer(1)) < 0) {
        throw std::runtime_error("acosh requires argument >= 1");
    }
    const Number x2 = multiply(value, value).with_context(work);
    const Number inner = subtract(x2, integer(1)).with_context(work);
    return ln(add(value, sqrt(inner)).with_context(work));
}

Number atanh(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    if (abs(value).compare(integer(1)) >= 0) {
        throw std::runtime_error("atanh requires |argument| < 1");
    }
    const Number one = integer(1);
    const Number numerator = add(one, value).with_context(work);
    const Number denominator = subtract(one, value).with_context(work);
    return divide(ln(numerator), integer(2));
}

Number sec(const Number& value) {
    return divide(integer(1), cos(value));
}

Number csc(const Number& value) {
    return divide(integer(1), sin(value));
}

Number cot(const Number& value) {
    return divide(cos(value), sin(value));
}

Number asec(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    if (abs(value).compare(integer(1)) < 0) {
        throw std::runtime_error("asec requires |argument| >= 1");
    }
    return acos(divide(integer(1), value));
}

Number acsc(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    if (abs(value).compare(integer(1)) < 0) {
        throw std::runtime_error("acsc requires |argument| >= 1");
    }
    return asin(divide(integer(1), value));
}

Number acot(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number half_pi = divide(pi(), integer(2));
    return subtract(half_pi, atan(value));
}

Number erf(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number x2 = multiply(value, value).with_context(work);
    Number term = value;
    term.set_context(work);
    Number sum = term;
    sum.set_context(work);
    const Number eps = tolerance(context);
    for (int n = 1; n < work.max_iterations; ++n) {
        term = divide(multiply(term, x2).with_context(work), integer(n)).with_context(work);
        const Number addend = divide(term, integer(2 * n + 1)).with_context(work);
        sum = (n % 2 == 0) ? add(sum, addend).with_context(work) : subtract(sum, addend).with_context(work);
        if (less_abs_than(addend, eps)) {
            break;
        }
    }
    const Number factor = divide(integer(2), sqrt(pi())).with_context(work);
    return round_to_context(multiply(factor, sum).with_context(context));
}

Number erfc(const Number& value) {
    return subtract(integer(1), erf(value));
}

Number log10(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);
    const Number ln10 = ln(integer(10)).with_context(work);
    return round_to_context(divide(ln(value), ln10).with_context(context));
}

namespace {

// Lanczos coefficients for gamma function
const Number kLanczosCoefficients[] = {
    Number(BigDecimal::from_string("0.99999999999980993")),
    Number(BigDecimal::from_string("676.5203681218851")),
    Number(BigDecimal::from_string("-1259.1392167224028")),
    Number(BigDecimal::from_string("771.32342877765313")),
    Number(BigDecimal::from_string("-176.61502916214059")),
    Number(BigDecimal::from_string("12.507343278686905")),
    Number(BigDecimal::from_string("-0.13857109526572012")),
    Number(BigDecimal::from_string("9.9843695780195716e-6")),
    Number(BigDecimal::from_string("1.5056327351493116e-7")),
};

Number log_gamma_positive(const Number& x) {
    const PrecisionContext& context = x.context();
    PrecisionContext work = work_context(context);
    const Number z = subtract(x, integer(1)).with_context(work);
    Number series = kLanczosCoefficients[0];
    series.set_context(work);
    for (int i = 1; i < 9; ++i) {
        series = add(series, divide(kLanczosCoefficients[i], add(z, integer(i)).with_context(work))).with_context(work);
    }
    const Number t = add(z, Number(BigDecimal::from_string("7.5"))).with_context(work);
    const Number half = Number(BigDecimal::from_string("0.5"));
    const Number two_pi = multiply(integer(2), pi()).with_context(work);
    Number result = add(
        multiply(half, ln(two_pi)).with_context(work),
        multiply(add(z, half).with_context(work), ln(t)).with_context(work)).with_context(work);
    result = subtract(result, t).with_context(work);
    result = add(result, ln(series)).with_context(work);
    return round_to_context(result.with_context(context));
}

}  // namespace

Number gamma(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);

    // Check for non-positive integers
    if (is_integer_value(value)) {
        BigInt int_val = value.is_integer() ? value.as_integer() :
            (value.is_rational() ? value.as_rational().numerator() :
             value.to_decimal().coefficient());
        if (int_val.sign() <= 0) {
            throw std::runtime_error("gamma is undefined for non-positive integers");
        }
    }

    Number x = value;

    // Reflection formula for x < 0.5
    if (x.compare(integer(0)) < 0) {
        const Number one = integer(1);
        const Number one_minus_x = subtract(one, x).with_context(work);
        const Number sin_pi_x = sin(multiply(pi(), x).with_context(work)).with_context(work);
        if (is_near_zero(sin_pi_x)) {
            throw std::runtime_error("gamma is undefined at this input");
        }
        return divide(pi(), multiply(sin_pi_x, gamma(one_minus_x)).with_context(work));
    }

    // For x >= 0.5, use Lanczos approximation
    return exp(log_gamma_positive(x));
}

Number beta(const Number& a, const Number& b) {
    const PrecisionContext& context = a.context();
    PrecisionContext work = work_context(context);
    const Number ga = gamma(a).with_context(work);
    const Number gb = gamma(b).with_context(work);
    const Number gab = gamma(add(a, b)).with_context(work);
    return round_to_context(divide(multiply(ga, gb), gab).with_context(context));
}

Number zeta(const Number& value) {
    const PrecisionContext& context = value.context();
    PrecisionContext work = work_context(context);

    // zeta(1) is undefined (pole)
    if (value.compare(integer(1)) == 0) {
        throw std::runtime_error("zeta is undefined at s = 1");
    }

    // For s < 0, use reflection formula
    if (value.compare(integer(0)) < 0) {
        const Number one = integer(1);
        const Number one_minus_s = subtract(one, value).with_context(work);
        const Number two_pow_s = pow(integer(2), value);
        const Number pi_pow = pow(pi(), subtract(value, one)).with_context(work);
        const Number sin_part = sin(multiply(pi(), divide(value, integer(2))).with_context(work)).with_context(work);
        const Number gamma_part = gamma(one_minus_s);
        const Number zeta_part = zeta(one_minus_s);
        return round_to_context(
            multiply(multiply(multiply(two_pow_s, pi_pow), sin_part),
                     multiply(gamma_part, zeta_part)).with_context(context));
    }

    // zeta(0) = -0.5
    if (is_zero(value)) {
        return Number(BigDecimal::from_string("-0.5")).with_context(context);
    }

    // Use Dirichlet eta series
    const Number one = integer(1);
    const Number denominator = subtract(one, pow(integer(2), subtract(one, value))).with_context(work);
    if (is_near_zero(denominator)) {
        throw std::runtime_error("zeta is numerically unstable near s = 1");
    }

    Number eta = integer(0);
    eta.set_context(work);
    const Number eps = tolerance(context);
    for (int n = 1; n < work.max_iterations; ++n) {
        Number term = divide(integer(n % 2 == 1 ? 1 : -1),
                            pow(integer(n), value)).with_context(work);
        eta = add(eta, term).with_context(work);
        if (abs(term).compare(eps) < 0) {
            break;
        }
    }
    return round_to_context(divide(eta, denominator).with_context(context));
}

Number bessel_j(int order, const Number& x) {
    const PrecisionContext& context = x.context();
    PrecisionContext work = work_context(context);

    if (is_zero(x)) {
        return integer(order == 0 ? 1 : 0);
    }

    const Number abs_x = abs(x);
    const Number half_x = divide(x, integer(2)).with_context(work);

    // For large x, use asymptotic expansion
    if (abs_x.compare(integer(50)) > 0) {
        const Number phase = subtract(subtract(abs_x, multiply(integer(order), divide(pi(), integer(2))).with_context(work)).with_context(work),
                                      divide(pi(), integer(4))).with_context(work);
        const Number amplitude = sqrt(divide(integer(2), multiply(pi(), abs_x)).with_context(work)).with_context(work);
        Number result = multiply(amplitude, cos(phase)).with_context(work);
        if (x.compare(integer(0)) < 0 && order % 2 != 0) {
            result = negate(result);
        }
        return round_to_context(result.with_context(context));
    }

    // Series expansion
    Number sum = integer(0);
    sum.set_context(work);
    Number term = divide(pow(half_x, integer(order)).with_context(work), gamma(add(integer(order), integer(1))).with_context(work)).with_context(work);
    const Number eps = tolerance(context);
    const Number half_x_sq = multiply(half_x, half_x).with_context(work);

    for (int k = 0; k < work.max_iterations; ++k) {
        sum = add(sum, term).with_context(work);
        if (abs(term).compare(eps) < 0) {
            break;
        }
        term = multiply(term, divide(half_x_sq,
                                     multiply(integer(k + 1), integer(k + order + 1))).with_context(work)).with_context(work);
        term = negate(term);
    }
    return round_to_context(sum.with_context(context));
}

bool is_near_zero(const Number& value) {
    return abs(value).compare(tolerance(value.context())) <= 0;
}

bool is_integer_value(const Number& value) {
    if (value.is_integer()) {
        return true;
    }
    if (value.is_rational()) {
        return value.as_rational().denominator() == BigInt(1);
    }
    if (value.is_decimal()) {
        const BigDecimal& dec = value.as_decimal();
        if (dec.scale() <= 0) {
            return true;
        }
        BigInt remainder = dec.coefficient().abs() % BigDecimal::pow10(dec.scale());
        return remainder.is_zero();
    }
    return false;
}

bool approximate_fraction(const Number& value,
                          BigInt* numerator,
                          BigInt* denominator,
                          long long max_denominator) {
    if (value.is_integer()) {
        *numerator = value.as_integer();
        *denominator = BigInt(1);
        return true;
    }
    if (value.is_rational()) {
        const Rational& r = value.as_rational();
        if (r.denominator() <= BigInt(max_denominator)) {
            *numerator = r.numerator();
            *denominator = r.denominator();
            return true;
        }
        return false;
    }
    if (value.is_decimal()) {
        const BigDecimal& dec = value.as_decimal();
        BigInt num = dec.coefficient();
        BigInt den = BigDecimal::pow10(dec.scale());
        BigInt g = BigInt::gcd(num.abs(), den);
        num = num / g;
        den = den / g;
        if (den <= BigInt(max_denominator)) {
            *numerator = num;
            *denominator = den;
            return true;
        }
    }
    return false;
}

bool is_near_zero(double value, double eps) {
    return std::abs(value) < eps;
}

bool is_integer_value(double value, double eps) {
    return std::abs(value - std::round(value)) < eps;
}

bool approximate_fraction(double value,
                          long long* numerator,
                          long long* denominator,
                          long long max_denominator,
                          double eps) {
    // Simple continued fraction algorithm
    double x = std::abs(value);
    long long a0 = static_cast<long long>(std::floor(x));
    long long h0 = 1, h1 = a0;
    long long k0 = 0, k1 = 1;
    double remainder = x - a0;

    while (std::abs(remainder) > eps && k1 <= max_denominator) {
        double reciprocal = 1.0 / remainder;
        long long a = static_cast<long long>(std::floor(reciprocal));
        remainder = reciprocal - a;

        long long h2 = a * h1 + h0;
        long long k2 = a * k1 + k0;

        if (k2 > max_denominator) {
            break;
        }

        h0 = h1; h1 = h2;
        k0 = k1; k1 = k2;
    }

    if (k1 <= max_denominator) {
        *numerator = (value < 0) ? -h1 : h1;
        *denominator = k1;
        return true;
    }
    return false;
}

}  // namespace numeric

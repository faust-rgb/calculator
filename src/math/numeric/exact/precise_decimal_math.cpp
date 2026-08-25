// ============================================================================
// PreciseDecimal 数学函数实现
// ============================================================================
//
// 本文件实现 PreciseDecimal 的数学函数：
// - 基本数学函数：abs, sqrt, pow, floor, ceil, round
// - 三角函数：sin, cos, tan, asin, acos, atan
// - 双曲函数：sinh, cosh, tanh, asinh, acosh, atanh
// - 指数对数函数：exp, ln, log10, exp2, log2
// - 常量：pi, e, ln2
// ============================================================================

#include "precise_decimal.h"
#include "math/numeric/rational/rational.h"
#include "core/value/stored_value.h"
#include "core/common/calculator_exceptions.h"
#include "math/mymath.h"
#include "math/runtime/precision/default_precision.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace precise {

// ============================================================================
// 辅助函数
// ============================================================================

PreciseDecimal one() { return PreciseDecimal(1LL); }
PreciseDecimal two() { return PreciseDecimal(2LL); }
PreciseDecimal half() { return PreciseDecimal("0.5"); }

PreciseDecimal decimal_from_uint(uint32_t value) {
    return PreciseDecimal(static_cast<long long>(value));
}

PreciseDecimal decimal_power_of_10(int exponent) {
    PreciseDecimal result(1LL);
    if (exponent >= 0) {
        result.data = multiply_bigint_by_power_of_10(result.data, exponent);
    } else {
        result.scale = -exponent;
    }
    result.normalize();
    return result;
}

PreciseDecimal scale_epsilon(int extra_digits = 0) {
    return PreciseDecimal("1e-" + std::to_string(PrecisionContext::get_default_scale() + extra_digits));
}

bool try_integral_exponent(const PreciseDecimal& value, long long* exponent) {
    if (!value.is_integer_value()) return false;
    *exponent = static_cast<long long>(value.to_double());
    return true;
}

bool try_thirds_exponent(const PreciseDecimal& value, long long* thirds) {
    PreciseDecimal scaled = value * PreciseDecimal(3LL);
    PreciseDecimal rounded = precise::round(scaled);
    if (precise::abs(scaled - rounded) > scale_epsilon(-2)) return false;
    *thirds = static_cast<long long>(rounded.to_double());
    return true;
}

int series_iteration_limit(int minimum) {
    int scale = PrecisionContext::get_default_scale();
    return std::max(minimum, scale + scale / 4 + 10);
}

void trim_fraction_scale(PreciseDecimal* value, int max_scale) {
    if (!value || value->is_inf || value->is_nan || value->scale <= max_scale) {
        return;
    }

    const int excess = value->scale - max_scale;
    value->data = divide_bigint_by_pow10(std::move(value->data), excess);
    value->scale = max_scale;
    ScopedNormalizationEnable enable;
    value->normalize();
}

// ============================================================================
// 常量
// ============================================================================

PreciseDecimal compute_pi_machin(int target_scale) {
    const int work_scale = std::max(target_scale, PrecisionContext::get_default_scale());
    ScopedPrecision guard(work_scale + 16);
    const PreciseDecimal epsilon = scale_epsilon(6);

    // Machin-like formula: pi/4 = 4*arctan(1/5) - arctan(1/239)
    auto compute_atan_inv_int = [&](long long denom) -> PreciseDecimal {
        PreciseDecimal d(denom);
        PreciseDecimal d2 = d * d;
        PreciseDecimal term = one() / d;
        PreciseDecimal sum = term;
        int limit = std::max(50, work_scale * 2 + 50);
        for (int i = 1; i < limit; ++i) {
            term = -term / d2;
            PreciseDecimal add = term / decimal_from_uint(static_cast<uint32_t>(2 * i + 1));
            sum += add;
            if (precise::abs(add) < epsilon) break;
        }
        return sum;
    };

    PreciseDecimal atan_1_5 = compute_atan_inv_int(5);
    PreciseDecimal atan_1_239 = compute_atan_inv_int(239);
    PreciseDecimal p = PreciseDecimal(4LL) * (PreciseDecimal(4LL) * atan_1_5 - atan_1_239);
    p.normalize();
    return p;
}

PreciseDecimal compute_e_series(int target_scale) {
    const int work_scale = std::max(target_scale, PrecisionContext::get_default_scale());
    ScopedPrecision guard(work_scale + 16);
    const PreciseDecimal epsilon = scale_epsilon(6);

    PreciseDecimal term = one();
    PreciseDecimal sum = two();
    int limit = std::max(50, work_scale + 50);
    for (int i = 2; i < limit; ++i) {
        term = term / decimal_from_uint(static_cast<uint32_t>(i));
        sum += term;
        if (precise::abs(term) < epsilon) break;
    }
    sum.normalize();
    return sum;
}

PreciseDecimal compute_ln2_series(int target_scale) {
    ScopedPrecision guard(16);
    // ln(2) = ln((1 + 1/3) / (1 - 1/3)) = 2 * (1/3 + 1/(3*3^3) + 1/(5*3^5) + ...)
    const int work_scale = std::max(target_scale, PrecisionContext::get_default_scale());
    const PreciseDecimal epsilon = scale_epsilon(6);

    PreciseDecimal one_third = one() / PreciseDecimal(3LL);
    PreciseDecimal term = one_third;
    PreciseDecimal sum = one_third;
    int limit = std::max(50, static_cast<int>(work_scale * 1.2) + 20);
    PreciseDecimal nine = decimal_from_uint(9);
    for (int i = 1; i < limit; ++i) {
        term = term / nine;
        PreciseDecimal add = term / decimal_from_uint(static_cast<uint32_t>(2 * i + 1));
        sum += add;
        if (precise::abs(add) < epsilon) break;
    }
    sum = two() * sum;
    sum.normalize();
    return sum;
}

PreciseDecimal pi() {
    int scale = PrecisionContext::get_default_scale();
    if (scale > 230) {
        return compute_pi_machin(scale);
    }
    static const PreciseDecimal p(
        "3.141592653589793238462643383279502884197169399375105820974944592307816406286"
        "208998628034825342117067982148086513282306647093844609550582231725359408128481"
        "117450284102701938521105559644622948954930381964428810975665933446128475648233");
    return p;
}

PreciseDecimal two_pi() {
    int scale = PrecisionContext::get_default_scale();
    if (scale > 230) {
        return two() * pi();
    }
    static const PreciseDecimal tp(
        "6.283185307179586476925286766559005768394338798750211641949889184615632812572"
        "417997256069650684234135964296173026564613294187689219101164434507187816256962"
        "234900568205403877042211119289245897909860763928857621951331866892256951296466");
    return tp;
}

PreciseDecimal half_pi() {
    int scale = PrecisionContext::get_default_scale();
    if (scale > 230) {
        return pi() / two();
    }
    static const PreciseDecimal hp(
        "1.5707963267948966192313216916397514420985846996875529104874722961539082031431"
        "0449931401741267105853399107404325664115332354692223077529111586285977040642405"
        "587251420513509692605527798223114744974651909822144054878329667323064237824116");
    return hp;
}

PreciseDecimal e() {
    int scale = PrecisionContext::get_default_scale();
    if (scale > 230) {
        return compute_e_series(scale);
    }
    static const PreciseDecimal val_e(
        "2.718281828459045235360287471352662497757247093699959574966967627724076630353"
        "547594571382178525166427427466391932003059921817413596629043572900334295260595"
        "630738132328627943490763233829880753195251019011573834187930702154089149934884");
    return val_e;
}

PreciseDecimal ln2() {
    int scale = PrecisionContext::get_default_scale();
    if (scale > 230) {
        return compute_ln2_series(scale);
    }
    static const PreciseDecimal val_ln2(
        "0.693147180559945309417232121458176568075500134360255254120680009493393621969"
        "694715605863326996418687542001481020570685733685520235758130557032670751635075"
        "96519478973694723683769222563258791851371881657850144027689067671439205263563");
    return val_ln2;
}

// ============================================================================
// 基本数学函数
// ============================================================================

PreciseDecimal abs(const PreciseDecimal& val) {
    PreciseDecimal res = val;
    res.negative = false;
    return res;
}

PreciseDecimal floor(const PreciseDecimal& val) {
    if (val.scale <= 0) return val;

    PreciseDecimal res;
    res.data = divide_bigint_by_pow10(val.data, val.scale);
    res.scale = 0;
    res.negative = val.negative;
    res.normalize();

    if (val.negative && !is_bigint_multiple_of_pow10(val.data, val.scale)) {
        res -= PreciseDecimal(1LL);
    }
    return res;
}

PreciseDecimal ceil(const PreciseDecimal& val) {
    if (val.scale <= 0) return val;
    if (is_bigint_multiple_of_pow10(val.data, val.scale)) {
        PreciseDecimal res = val;
        res.normalize();
        return res;
    }
    PreciseDecimal f = precise::floor(val);
    return f + PreciseDecimal(1LL);
}

PreciseDecimal round(const PreciseDecimal& val) {
    if (val.is_zero()) return val;
    PreciseDecimal half_val = half();
    if (val.negative) return precise::ceil(val - half_val);
    return precise::floor(val + half_val);
}

// ============================================================================
// 辅助函数实现
// ============================================================================

PreciseDecimal sqrt_initial_guess(const PreciseDecimal& val) {
    std::string digits = bigint_to_string(val.data);
    int decimal_exponent = static_cast<int>(digits.size()) - val.scale - 1;
    int mantissa_digits = std::min<int>(18, digits.size());
    std::string mantissa_text = digits.substr(0, static_cast<std::size_t>(mantissa_digits));

    PreciseDecimal mantissa = PreciseDecimal::from_digits(mantissa_text, mantissa_digits - 1, false);
    if (decimal_exponent % 2 != 0) {
        mantissa *= PreciseDecimal(10LL);
        --decimal_exponent;
    }

    PreciseDecimal x = mantissa > one() ? mantissa : one();
    const PreciseDecimal one_half = half();
    for (int i = 0; i < 8; ++i) {
        x = one_half * (x + mantissa / x);
    }

    return x * decimal_power_of_10(decimal_exponent / 2);
}

PreciseDecimal sqrt(const PreciseDecimal& val) {
    if (val.is_zero()) return PreciseDecimal(0LL);
    if (val.negative) throw PreciseDecimalUnsupported("sqrt of negative number");

    int target_scale = PrecisionContext::get_default_scale();

    PreciseDecimal x = sqrt_initial_guess(val);
    const PreciseDecimal one_half = half();
    const PreciseDecimal epsilon = scale_epsilon(4);

    NormalizationSuppressor suppressor;

    const int iterations = std::max(16, target_scale / 4 + 20);
    for (int i = 0; i < iterations; ++i) {
        PreciseDecimal next = one_half * (x + val / x);
        if (precise::abs(next - x) <= epsilon * std::max(one(), precise::abs(next))) {
            ScopedNormalizationEnable enable;
            next.normalize();
            return next;
        }
        x = next;
    }
    ScopedNormalizationEnable enable;
    x.normalize();
    return x;
}

PreciseDecimal cbrt_precise_decimal(const PreciseDecimal& val) {
    if (val.is_zero()) return PreciseDecimal(0LL);
    PreciseDecimal x = sqrt_initial_guess(precise::abs(val));
    if (x.is_zero()) x = one();
    const PreciseDecimal abs_val = precise::abs(val);
    for (int i = 0; i < std::max(14, PrecisionContext::get_default_scale() / 6 + 8); ++i) {
        x = (two() * x + abs_val / (x * x)) / PreciseDecimal(3LL);
    }
    return val.negative ? -x : x;
}

PreciseDecimal cbrt(const PreciseDecimal& val) {
    return cbrt_precise_decimal(val);
}

// ============================================================================
// 幂运算
// ============================================================================

PreciseDecimal pow_precise_decimal(const PreciseDecimal& base, long long exponent) {
    if (exponent < 0) throw PreciseDecimalUnsupported("negative exponent not supported for precise decimal");
    if (exponent == 0) return PreciseDecimal::from_integer_string("1", false);
    if (exponent == 1) return base;

    if (base.is_zero()) return PreciseDecimal(0LL);
    if (base == one()) return one();
    if (base == two()) {
        BigIntData result = {1};
        for (long long i = 0; i < exponent; ++i) {
            result = multiply_bigint_by_uint32(result, 2);
        }
        PreciseDecimal res;
        res.data = result;
        res.scale = 0;
        res.negative = false;
        return res;
    }

    PreciseDecimal res = PreciseDecimal::from_integer_string("1", false);
    PreciseDecimal b = base;
    long long exp = exponent;

    while (exp > 0) {
        if (exp & 1) {
            res = multiply_precise_decimal(res, b);
        }
        b = multiply_precise_decimal(b, b);
        exp >>= 1;
    }
    return res;
}

PreciseDecimal pow(const PreciseDecimal& base, long long exp) {
    return pow_precise_decimal(base, exp);
}

PreciseDecimal pow(const PreciseDecimal& base, const PreciseDecimal& exp) {
    if (exp.is_zero()) return one();
    if (base.is_zero()) {
        if (exp.negative) throw std::domain_error("0^negative is undefined");
        return PreciseDecimal(0LL);
    }
    if (base == one()) return one();

    long long integer_exp = 0;
    if (try_integral_exponent(exp, &integer_exp)) {
        if (integer_exp < 0) {
            return one() / pow_precise_decimal(base, -integer_exp);
        }
        return pow_precise_decimal(base, integer_exp);
    }

    if (base.negative) {
        PreciseDecimal abs_base = precise::abs(base);
        long long thirds = 0;
        if (try_thirds_exponent(exp, &thirds) && thirds % 3 != 0) {
            PreciseDecimal result = precise::exp(exp * precise::ln(abs_base));
            if (thirds % 2 != 0) {
                result.negative = true;
            }
            return result;
        }

        throw std::domain_error("negative base with non-integer exponent");
    }

    return precise::exp(precise::ln(base) * exp);
}

// ============================================================================
// 指数和对数函数
// ============================================================================

PreciseDecimal exp2(const PreciseDecimal& x) {
    if (x.is_zero()) return one();

    ScopedPrecision guard(8);

    PreciseDecimal int_part = floor(x);
    PreciseDecimal frac_part = x - int_part;

    if (int_part.negative) {
        return one() / exp2(-x);
    }

    long long n = static_cast<long long>(int_part.to_double());
    PreciseDecimal int_result = one();
    if (n > 0) {
        BigIntData power_of_2 = {1};
        long long exp = n;
        BigIntData base_pow = {2};
        while (exp > 0) {
            if (exp & 1) {
                power_of_2 = multiply_bigint(power_of_2, base_pow);
            }
            base_pow = multiply_bigint(base_pow, base_pow);
            exp >>= 1;
        }
        int_result.data = power_of_2;
        int_result.scale = 0;
        int_result.normalize();
    }

    PreciseDecimal y = frac_part;
    PreciseDecimal ln2_val = ln2();

    NormalizationSuppressor suppressor;

    PreciseDecimal term = one();
    PreciseDecimal sum = one();
    PreciseDecimal y_ln2 = y * ln2_val;
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(150);

    for (int i = 1; i < limit; ++i) {
        term = term * y_ln2 / decimal_from_uint(static_cast<uint32_t>(i));
        sum += term;
        if (precise::abs(term) < epsilon) break;
    }

    g_suppress_normalization = false;
    sum.normalize();

    return int_result * sum;
}

PreciseDecimal log2(const PreciseDecimal& x) {
    if (x <= PreciseDecimal(0LL)) throw std::domain_error("log2 of non-positive number");
    if (x == one()) return PreciseDecimal(0LL);
    if (x == two()) return one();

    return ln(x) / ln2();
}

PreciseDecimal ln_agm(const PreciseDecimal& x) {
    if (x <= PreciseDecimal(0LL)) throw std::domain_error("ln of non-positive number");
    if (x == one()) return PreciseDecimal(0LL);
    if (x < one()) {
        return -ln_agm(one() / x);
    }

    const int user_scale = PrecisionContext::get_default_scale();
    ScopedPrecision guard(user_scale / 4 + 16);
    const int work_scale = PrecisionContext::get_default_scale();
    const PreciseDecimal epsilon = scale_epsilon(8);

    int m = static_cast<int>(user_scale * 1.7) + 8;

    PreciseDecimal two_pow_m;
    {
        BigIntData power_of_2 = {1};
        long long exp_val = m;
        BigIntData base_pow = {2};
        while (exp_val > 0) {
            if (exp_val & 1) {
                power_of_2 = multiply_bigint(power_of_2, base_pow);
            }
            base_pow = multiply_bigint(base_pow, base_pow);
            exp_val >>= 1;
        }
        two_pow_m.data = power_of_2;
        two_pow_m.scale = 0;
        two_pow_m.normalize();
    }
    PreciseDecimal s = x * two_pow_m;

    PreciseDecimal a = one();
    PreciseDecimal b = PreciseDecimal(4LL) / s;

    const int max_iter = std::max(30, work_scale / 4 + 10);
    for (int i = 0; i < max_iter; ++i) {
        PreciseDecimal next_a = half() * (a + b);
        PreciseDecimal next_b = precise::sqrt(a * b);
        PreciseDecimal diff = precise::abs(next_a - next_b);
        a = next_a;
        b = next_b;
        if (diff < epsilon) break;
    }

    PreciseDecimal pi_val = pi();
    PreciseDecimal ln2_val = ln2();
    PreciseDecimal res = pi_val / (two() * a) - PreciseDecimal(static_cast<long long>(m)) * ln2_val;
    res.normalize();
    return res;
}

PreciseDecimal ln_near_one(const PreciseDecimal& x) {
    const PreciseDecimal z = (x - one()) / (x + one());
    const PreciseDecimal z2 = z * z;

    NormalizationSuppressor suppressor;

    PreciseDecimal term = z;
    PreciseDecimal sum = z;
    const PreciseDecimal epsilon = scale_epsilon(2);
    const int limit = series_iteration_limit(200);

    for (int n = 3; n < limit * 2; n += 2) {
        term *= z2;
        PreciseDecimal add = term / decimal_from_uint(static_cast<uint32_t>(n));
        sum += add;
        if (precise::abs(add) < epsilon) break;
    }

    g_suppress_normalization = false;
    sum.normalize();
    return two() * sum;
}

PreciseDecimal ln(const PreciseDecimal& x) {
    if (x <= PreciseDecimal(0LL)) throw std::domain_error("ln of non-positive number");
    if (x == one()) return PreciseDecimal(0LL);

    int scale = PrecisionContext::get_default_scale();
    if (scale > 100) {
        return ln_agm(x);
    }

    ScopedPrecision guard(10);

    const PreciseDecimal ln2_val = ln2();
    const PreciseDecimal two_val = two();
    const PreciseDecimal half_val = half();
    const PreciseDecimal lower("0.7071067811865475");
    const PreciseDecimal upper("1.414213562373095");

    PreciseDecimal reduced = x;
    int k = 0;

    while (reduced > upper) {
        reduced *= half_val;
        ++k;
    }
    while (reduced < lower) {
        reduced *= two_val;
        --k;
    }

    PreciseDecimal result = ln_near_one(reduced);
    if (k != 0) {
        result = result + PreciseDecimal(static_cast<long long>(k)) * ln2_val;
    }
    result.normalize();
    return result;
}

PreciseDecimal ln10() {
    int scale = PrecisionContext::get_default_scale();
    if (scale > 230) {
        return ln(PreciseDecimal(10LL));
    }
    static const PreciseDecimal val_ln10(
        "2.302585092994045684017991454684364207601101488628772976033327900967572609677"
        "352480235997205089598298341967784042286248633409525465082806756666287369098781"
        "68948290720832555468084287239488982714848914692334149091773093118831437690469");
    return val_ln10;
}

PreciseDecimal log10(const PreciseDecimal& x) {
    return ln(x) / ln10();
}

PreciseDecimal exp(const PreciseDecimal& x) {
    if (x.is_zero()) return one();

    const long double x_ld = static_cast<long double>(x);
    if (x_ld >= mymath::kLnDoubleMax) return PreciseDecimal::infinity();
    if (x_ld <= mymath::kLnDoubleDenormMin) return PreciseDecimal(0LL);

    const long double x_magnitude = mymath::abs(x_ld);
    const int extra_scale = std::max(8, static_cast<int>(x_magnitude / std::log(10.0L)) + 8);
    ScopedPrecision guard(extra_scale);
    const int work_scale = PrecisionContext::get_default_scale() + 8;

    if (x.negative) {
        PreciseDecimal abs_x = precise::abs(x);
        return one() / exp(abs_x);
    }

    int scale = PrecisionContext::get_default_scale();
    if (scale > 100) {
        const PreciseDecimal ln2_val = ln2();
        PreciseDecimal k_val = precise::floor(x / ln2_val);
        long long k = static_cast<long long>(k_val.to_double());
        PreciseDecimal r = x - k_val * ln2_val;

        NormalizationSuppressor suppressor;
        PreciseDecimal term = one();
        PreciseDecimal sum = one();
        const PreciseDecimal epsilon = scale_epsilon();
        const int limit = series_iteration_limit(100);
        for (int i = 1; i < limit; ++i) {
            term = term * r / decimal_from_uint(static_cast<uint32_t>(i));
            trim_fraction_scale(&term, work_scale);
            sum += term;
            trim_fraction_scale(&sum, work_scale);
            if (precise::abs(term) < epsilon) break;
        }

        ScopedNormalizationEnable enable;
        sum.normalize();

        if (k != 0) {
            BigIntData power_of_2 = {1};
            long long exp = k > 0 ? k : -k;
            BigIntData base_pow = {2};
            while (exp > 0) {
                if (exp & 1) {
                    power_of_2 = multiply_bigint(power_of_2, base_pow);
                }
                base_pow = multiply_bigint(base_pow, base_pow);
                exp >>= 1;
            }
            PreciseDecimal two_pow_k;
            two_pow_k.data = power_of_2;
            two_pow_k.scale = 0;
            two_pow_k.normalize();

            if (k > 0) {
                return sum * two_pow_k;
            } else {
                return sum / two_pow_k;
            }
        }
        return sum;
    }

    int k = 0;
    PreciseDecimal r = x;
    const PreciseDecimal threshold("0.01");
    while (r > threshold) {
        r /= two();
        k++;
    }

    NormalizationSuppressor suppressor;

    PreciseDecimal term = one();
    PreciseDecimal sum = one();
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(150);
    for (int i = 1; i < limit; ++i) {
        term = term * r / decimal_from_uint(static_cast<uint32_t>(i));
        trim_fraction_scale(&term, work_scale);
        sum += term;
        trim_fraction_scale(&sum, work_scale);
        if (precise::abs(term) < epsilon) break;
    }
    for (int i = 0; i < k; ++i) {
        sum = sum * sum;
        trim_fraction_scale(&sum, work_scale);
    }

    ScopedNormalizationEnable enable;
    sum.normalize();
    return sum;
}

// ============================================================================
// 三角函数
// ============================================================================

PreciseDecimal reduce_mod_positive(const PreciseDecimal& value, const PreciseDecimal& modulus) {
    PreciseDecimal q = precise::floor(value / modulus);
    PreciseDecimal r = value - q * modulus;
    if (r < PreciseDecimal(0LL)) r += modulus;
    return r;
}

PreciseDecimal reduce_trig_argument(const PreciseDecimal& x) {
    int extra_guard = 8;
    if (precise::abs(x) > PreciseDecimal(100LL)) {
        long double x_d = mymath::abs(static_cast<long double>(x));
        int log10_x = static_cast<int>(std::log10(std::max(1.0L, x_d)));
        extra_guard += log10_x;
    }
    ScopedPrecision guard(extra_guard);
    const PreciseDecimal p = precise::pi();
    const PreciseDecimal two_pi_val = precise::two_pi();
    PreciseDecimal r = reduce_mod_positive(x, two_pi_val);
    if (r > p) r -= two_pi_val;
    r.normalize();
    return r;
}

PreciseDecimal reduce_sin_argument(const PreciseDecimal& x, bool* negate) {
    PreciseDecimal r = reduce_trig_argument(x);
    const PreciseDecimal p = precise::pi();
    const PreciseDecimal half_p = precise::half_pi();
    *negate = false;
    if (r > half_p) {
        r = p - r;
    } else if (r < -half_p) {
        r = -p - r;
    }
    return r;
}

PreciseDecimal reduce_cos_argument(const PreciseDecimal& x, bool* negate) {
    PreciseDecimal r = reduce_trig_argument(x);
    const PreciseDecimal p = precise::pi();
    const PreciseDecimal half_p = precise::half_pi();
    *negate = false;
    if (r > half_p) {
        r = p - r;
        *negate = true;
    } else if (r < -half_p) {
        r = -p - r;
        *negate = true;
    }
    return r;
}

PreciseDecimal sin(const PreciseDecimal& x) {
    if (x.is_zero()) return PreciseDecimal(0LL);

    ScopedPrecision guard(8);
    bool negate = false;
    PreciseDecimal r = reduce_sin_argument(x, &negate);

    int scale = PrecisionContext::get_default_scale();
    const int work_scale = scale + 8;

    if (scale > 100) {
        int k = 0;
        const PreciseDecimal threshold("0.01");
        while (precise::abs(r) > threshold && k < 20) {
            r *= half();
            ++k;
        }

        NormalizationSuppressor suppressor;

        PreciseDecimal r2 = r * r;
        trim_fraction_scale(&r2, work_scale);
        PreciseDecimal sin_r = r;
        PreciseDecimal cos_r = one();
        PreciseDecimal term_s = r;
        PreciseDecimal term_c = one();
        const PreciseDecimal epsilon = scale_epsilon();
        const int limit = series_iteration_limit(80);

        for (int i = 1; i < limit; ++i) {
            term_c = -term_c * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i - 1) * (2 * i)));
            term_s = -term_s * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i) * (2 * i + 1)));
            trim_fraction_scale(&term_c, work_scale);
            trim_fraction_scale(&term_s, work_scale);
            cos_r += term_c;
            sin_r += term_s;
            trim_fraction_scale(&cos_r, work_scale);
            trim_fraction_scale(&sin_r, work_scale);
            if (precise::abs(term_s) < epsilon && precise::abs(term_c) < epsilon) break;
        }

        for (int i = 0; i < k; ++i) {
            PreciseDecimal new_sin = two() * sin_r * cos_r;
            PreciseDecimal new_cos = cos_r * cos_r - sin_r * sin_r;
            trim_fraction_scale(&new_sin, work_scale);
            trim_fraction_scale(&new_cos, work_scale);
            sin_r = new_sin;
            cos_r = new_cos;
        }

        ScopedNormalizationEnable enable;
        sin_r.normalize();
        return negate ? -sin_r : sin_r;
    }

    NormalizationSuppressor suppressor;

    PreciseDecimal term = r;
    PreciseDecimal sum = r;
    PreciseDecimal r2 = r * r;
    trim_fraction_scale(&r2, work_scale);
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(150);
    for (int i = 1; i < limit; ++i) {
        term = -term * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i) * (2 * i + 1)));
        trim_fraction_scale(&term, work_scale);
        sum += term;
        trim_fraction_scale(&sum, work_scale);
        if (precise::abs(term) < epsilon) break;
    }

    ScopedNormalizationEnable enable;
    sum.normalize();
    return negate ? -sum : sum;
}

PreciseDecimal cos(const PreciseDecimal& x) {
    if (x.is_zero()) return one();

    ScopedPrecision guard(8);
    bool negate = false;
    PreciseDecimal r = reduce_cos_argument(x, &negate);

    int scale = PrecisionContext::get_default_scale();
    const int work_scale = scale + 8;

    if (scale > 100) {
        int k = 0;
        const PreciseDecimal threshold("0.01");
        while (precise::abs(r) > threshold && k < 20) {
            r *= half();
            ++k;
        }

        NormalizationSuppressor suppressor;

        PreciseDecimal r2 = r * r;
        trim_fraction_scale(&r2, work_scale);
        PreciseDecimal sin_r = r;
        PreciseDecimal cos_r = one();
        PreciseDecimal term_s = r;
        PreciseDecimal term_c = one();
        const PreciseDecimal epsilon = scale_epsilon();
        const int limit = series_iteration_limit(80);

        for (int i = 1; i < limit; ++i) {
            term_c = -term_c * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i - 1) * (2 * i)));
            term_s = -term_s * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i) * (2 * i + 1)));
            trim_fraction_scale(&term_c, work_scale);
            trim_fraction_scale(&term_s, work_scale);
            cos_r += term_c;
            sin_r += term_s;
            trim_fraction_scale(&cos_r, work_scale);
            trim_fraction_scale(&sin_r, work_scale);
            if (precise::abs(term_s) < epsilon && precise::abs(term_c) < epsilon) break;
        }

        for (int i = 0; i < k; ++i) {
            PreciseDecimal new_sin = two() * sin_r * cos_r;
            PreciseDecimal new_cos = cos_r * cos_r - sin_r * sin_r;
            trim_fraction_scale(&new_sin, work_scale);
            trim_fraction_scale(&new_cos, work_scale);
            sin_r = new_sin;
            cos_r = new_cos;
        }

        ScopedNormalizationEnable enable;
        cos_r.normalize();
        return negate ? -cos_r : cos_r;
    }

    NormalizationSuppressor suppressor;

    PreciseDecimal term = one();
    PreciseDecimal sum = one();
    PreciseDecimal r2 = r * r;
    trim_fraction_scale(&r2, work_scale);
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(150);
    for (int i = 1; i < limit; ++i) {
        term = -term * r2 / decimal_from_uint(static_cast<uint32_t>((2 * i - 1) * (2 * i)));
        trim_fraction_scale(&term, work_scale);
        sum += term;
        trim_fraction_scale(&sum, work_scale);
        if (precise::abs(term) < epsilon) break;
    }

    ScopedNormalizationEnable enable;
    sum.normalize();
    return negate ? -sum : sum;
}

PreciseDecimal tan(const PreciseDecimal& x) {
    PreciseDecimal p = pi();
    PreciseDecimal half_p = p / two();
    PreciseDecimal r = x - precise::floor(x / p) * p;
    if (precise::abs(r - half_p) < scale_epsilon(4)) {
        throw std::domain_error("tan undefined near singularity");
    }
    PreciseDecimal s = sin(x);
    PreciseDecimal c = cos(x);
    if (c.is_zero()) throw std::domain_error("tan undefined (cos is zero)");
    return s / c;
}

PreciseDecimal asin(const PreciseDecimal& x) {
    if (precise::abs(x) > one()) throw std::domain_error("asin out of range");
    if (x == one()) return pi() / two();
    if (x == -one()) return -pi() / two();
    return atan(x / precise::sqrt(one() - x * x));
}

PreciseDecimal acos(const PreciseDecimal& x) {
    return pi() / two() - asin(x);
}

PreciseDecimal atan(const PreciseDecimal& x) {
    if (x.is_zero()) return PreciseDecimal(0LL);

    ScopedPrecision guard(8);

    if (x == one()) return pi() / PreciseDecimal(4LL);
    if (x == -one()) return -pi() / PreciseDecimal(4LL);
    if (precise::abs(x) > one()) {
        if (x.negative) return -pi() / two() - atan(one() / x);
        return pi() / two() - atan(one() / x);
    }

    NormalizationSuppressor suppressor;

    PreciseDecimal term = x;
    PreciseDecimal sum = x;
    const PreciseDecimal x2 = x * x;
    const PreciseDecimal epsilon = scale_epsilon();
    const int limit = series_iteration_limit(200);
    for (int i = 1; i < limit; ++i) {
        term = -term * x2;
        PreciseDecimal add = term / decimal_from_uint(static_cast<uint32_t>(2 * i + 1));
        sum += add;
        if (precise::abs(add) < epsilon) break;
    }

    g_suppress_normalization = false;
    sum.normalize();
    return sum;
}

// ============================================================================
// 双曲函数
// ============================================================================

PreciseDecimal sinh(const PreciseDecimal& x) {
    PreciseDecimal ex = exp(x);
    return (ex - one() / ex) / two();
}

PreciseDecimal cosh(const PreciseDecimal& x) {
    PreciseDecimal ex = exp(x);
    return (ex + one() / ex) / two();
}

PreciseDecimal tanh(const PreciseDecimal& x) {
    PreciseDecimal ex = exp(x);
    PreciseDecimal einv = one() / ex;
    return (ex - einv) / (ex + einv);
}

PreciseDecimal asinh(const PreciseDecimal& x) {
    return ln(x + precise::sqrt(x * x + one()));
}

PreciseDecimal acosh(const PreciseDecimal& x) {
    if (x < one()) throw std::domain_error("acosh out of range");
    return ln(x + precise::sqrt(x * x - one()));
}

PreciseDecimal atanh(const PreciseDecimal& x) {
    if (precise::abs(x) >= one()) throw std::domain_error("atanh out of range");
    return ln((one() + x) / (one() - x)) / two();
}

// ============================================================================
// 算术运算
// ============================================================================

void align_precise_scales(PreciseDecimal* lhs, PreciseDecimal* rhs) {
    if (lhs->scale == rhs->scale) return;
    if (lhs->scale < rhs->scale) {
        lhs->data = multiply_bigint_by_power_of_10(lhs->data, rhs->scale - lhs->scale);
        lhs->scale = rhs->scale;
    } else {
        rhs->data = multiply_bigint_by_power_of_10(rhs->data, lhs->scale - rhs->scale);
        rhs->scale = lhs->scale;
    }
}

PreciseDecimal add_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (lhs.is_inf && rhs.is_inf) {
        if (lhs.negative != rhs.negative) {
            return PreciseDecimal::nan();
        }
        return lhs;
    }
    if (lhs.is_inf) return lhs;
    if (rhs.is_inf) return rhs;

    if (lhs.is_nan || rhs.is_nan) return PreciseDecimal::nan();

    if (lhs.negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return subtract_precise_decimal(lhs, rhs_flipped);
    }

    PreciseDecimal l = lhs;
    PreciseDecimal r = rhs;
    align_precise_scales(&l, &r);

    PreciseDecimal res;
    res.scale = l.scale;
    res.negative = lhs.negative;

    size_t max_len = std::max(l.data.size(), r.data.size()) + 1;
    res.data.resize(max_len);
    size_t new_size = add_raw(l.data.data(), l.data.size(), r.data.data(), r.data.size(), res.data.data());
    res.data.resize(new_size);
    res.normalize();
    return res;
}

PreciseDecimal subtract_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (lhs.is_inf && rhs.is_inf) {
        if (lhs.negative == rhs.negative) {
            return PreciseDecimal::nan();
        }
        return lhs;
    }
    if (lhs.is_inf) return lhs;
    if (rhs.is_inf) {
        PreciseDecimal res = rhs;
        res.negative = !res.negative;
        return res;
    }

    if (lhs.is_nan || rhs.is_nan) return PreciseDecimal::nan();

    if (lhs.negative != rhs.negative) {
        PreciseDecimal rhs_flipped = rhs;
        rhs_flipped.negative = !rhs_flipped.negative;
        return add_precise_decimal(lhs, rhs_flipped);
    }

    PreciseDecimal l = lhs;
    PreciseDecimal r = rhs;
    align_precise_scales(&l, &r);

    int cmp = compare_raw(l.data.data(), l.data.size(), r.data.data(), r.data.size());
    if (cmp == 0) return {};

    PreciseDecimal res;
    res.scale = l.scale;
    res.data.resize(std::max(l.data.size(), r.data.size()));

    if (cmp > 0) {
        size_t new_size = sub_raw(l.data.data(), l.data.size(), r.data.data(), r.data.size(), res.data.data());
        res.data.resize(new_size);
        res.negative = lhs.negative;
    } else {
        size_t new_size = sub_raw(r.data.data(), r.data.size(), l.data.data(), l.data.size(), res.data.data());
        res.data.resize(new_size);
        res.negative = !lhs.negative;
    }
    res.normalize();
    return res;
}

PreciseDecimal multiply_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (lhs.is_nan || rhs.is_nan) return PreciseDecimal::nan();
    if (lhs.is_inf || rhs.is_inf) {
        if (lhs.is_zero() || rhs.is_zero()) {
            return PreciseDecimal::nan();
        }
        PreciseDecimal res = lhs.is_inf ? lhs : rhs;
        res.negative = lhs.negative != rhs.negative;
        return res;
    }

    PreciseDecimal res;
    res.data = multiply_bigint(lhs.data, rhs.data);
    res.scale = lhs.scale + rhs.scale;
    res.negative = lhs.negative != rhs.negative;
    res.normalize();
    return res;
}

PreciseDecimal divide_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (lhs.is_nan || rhs.is_nan) return PreciseDecimal::nan();
    if (lhs.is_inf && rhs.is_inf) {
        return PreciseDecimal::nan();
    }
    if (lhs.is_inf) {
        PreciseDecimal res = lhs;
        res.negative = lhs.negative != rhs.negative;
        return res;
    }
    if (rhs.is_inf) {
        return PreciseDecimal(0LL);
    }

    if (rhs.is_zero()) throw std::runtime_error("division by zero");
    if (lhs.is_zero()) return {};

    int target_scale = PrecisionContext::get_default_scale();
    int numerator_shift = target_scale + rhs.scale - lhs.scale;
    BigIntData numerator = lhs.data;
    if (numerator_shift >= 0) {
        numerator = multiply_bigint_by_power_of_10(numerator, numerator_shift);
    } else {
        BigIntData divisor = multiply_bigint_by_power_of_10({1}, -numerator_shift);
        BigIntData truncated, ignored_remainder;
        div_bigint(numerator, divisor, &truncated, &ignored_remainder);
        if (!ignored_remainder.empty() &&
            compare_bigint(multiply_bigint_by_uint32(ignored_remainder, 2), divisor) >= 0) {
            truncated = add_bigint(truncated, BigIntData{1});
        }
        numerator = truncated;
    }
    BigIntData denominator = rhs.data;

    BigIntData q, r;
    div_bigint(numerator, denominator, &q, &r);

    // Round to nearest, ties away from zero instead of silently truncating.
    if (!r.empty() && compare_bigint(multiply_bigint_by_uint32(r, 2), denominator) >= 0) {
        q = add_bigint(q, BigIntData{1});
    }

    PreciseDecimal res;
    res.data = q;
    res.scale = target_scale;
    res.negative = lhs.negative != rhs.negative;
    res.normalize();
    return res;
}

int compare_precise_decimal(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (lhs.is_zero() && rhs.is_zero()) return 0;
    if (lhs.negative != rhs.negative) return lhs.negative ? -1 : 1;

    int res = 0;
    if (lhs.scale == rhs.scale) {
        res = compare_bigint(lhs.data, rhs.data);
    } else if (lhs.scale < rhs.scale) {
        BigIntData shifted_lhs = multiply_bigint_by_power_of_10(lhs.data, rhs.scale - lhs.scale);
        res = compare_bigint(shifted_lhs, rhs.data);
    } else {
        BigIntData shifted_rhs = multiply_bigint_by_power_of_10(rhs.data, lhs.scale - rhs.scale);
        res = compare_bigint(lhs.data, shifted_rhs);
    }

    return lhs.negative ? -res : res;
}

// ============================================================================
// 辅助输出函数
// ============================================================================

std::string rational_to_precise_decimal_text(const Rational& value) {
    PreciseDecimal numerator = PreciseDecimal::from_integer_string(
        std::to_string(value.numerator < 0 ? -value.numerator : value.numerator),
        value.numerator < 0);
    const PreciseDecimal denominator =
        PreciseDecimal::from_integer_string(std::to_string(value.denominator), false);
    return divide_precise_decimal(numerator, denominator).to_string();
}

std::string stored_value_precise_decimal_text(const StoredValue& value) {
    if (value.exact) return rational_to_precise_decimal_text(value.rational);
    if (value.has_precise_decimal_text) return value.precise_decimal_text;
    return format_decimal(normalize_display_decimal(value.decimal.to_long_double()), 15);
}

} // namespace precise

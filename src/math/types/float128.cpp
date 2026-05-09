/**
 * @file float128.cpp
 * @brief Implementation of 128-bit high-precision floating-point arithmetic
 *
 * This file implements mathematical operations for float128_t type using
 * Double-Long-Double arithmetic, providing approximately 36-38 significant digits.
 *
 * Key techniques:
 * - Error-Free Transformation (EFT)
 * - Kahan compensated summation
 * - Taylor series with range reduction
 */

#include "float128.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <climits>

namespace mymath {

namespace {

// Error-Free Transformation: a + b = s + e
inline float128_t two_sum(long double a, long double b) {
    long double s = a + b;
    long double v = s - a;
    long double e = (a - (s - v)) + (b - v);
    return {s, e};
}

// Kahan summation helper for improved precision in series summation
class KahanAccumulator {
public:
    KahanAccumulator() : sum_{0.0L, 0.0L}, compensation_{0.0L, 0.0L} {}

    void add(float128_t term) {
        // Kahan summation: y = term - compensation
        float128_t y = term - compensation_;
        // t = sum + y
        float128_t t = sum_ + y;
        // compensation = (t - sum) - y (new compensation value)
        compensation_ = (t - sum_) - y;
        sum_ = t;
    }

    float128_t get() const { return sum_; }

private:
    float128_t sum_;
    float128_t compensation_;
};

// Error-Free Transformation: a * b = p + e
// Uses Dekker's algorithm for multiplication decomposition
inline float128_t two_prod(long double a, long double b) {
    long double p = a * b;
#ifdef __FMA__
    long double e = __builtin_fmal(a, b, -p);
#else
    // Dekker's algorithm
    constexpr long double kSplitLD = 4294967297.0L; // 2^32 + 1 for 64-bit mantissa

    auto split = [](long double x) {
        long double t = x * kSplitLD;
        long double x_hi = t - (t - x);
        long double x_lo = x - x_hi;
        return std::make_pair(x_hi, x_lo);
    };

    auto a_split = split(a);
    auto b_split = split(b);

    long double e = ((a_split.first * b_split.first - p) +
                     a_split.first * b_split.second +
                     a_split.second * b_split.first) +
                    a_split.second * b_split.second;
#endif
    return {p, e};
}

// Forward declarations for functions needed from core layer
// These will be resolved at link time
long double infinity_impl();
long double quiet_nan_impl();

} // namespace

float128_t operator+(float128_t a, float128_t b) {
    float128_t s = two_sum(a.hi, b.hi);
    s.lo += a.lo + b.lo;
    return two_sum(s.hi, s.lo);
}

float128_t operator-(float128_t a, float128_t b) {
    return a + (-b);
}

float128_t operator-(float128_t a) {
    return {-a.hi, -a.lo};
}

float128_t operator*(float128_t a, float128_t b) {
    float128_t p = two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return two_sum(p.hi, p.lo);
}

float128_t operator/(float128_t a, float128_t b) {
    long double q1 = a.hi / b.hi;
    float128_t r = a - (b * q1);
    long double q2 = r.hi / b.hi;
    return two_sum(q1, q2);
}

bool operator==(float128_t a, float128_t b) {
    return a.hi == b.hi && a.lo == b.lo;
}

bool operator!=(float128_t a, float128_t b) {
    return !(a == b);
}

bool operator<(float128_t a, float128_t b) {
    return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
}

namespace precise128 {

// ============================================================================
// Basic utility functions
// ============================================================================

bool isnan(float128_t a) {
    return std::isnan(a.hi) || std::isnan(a.lo);
}

bool isinf(float128_t a) {
    return std::isinf(a.hi);
}

bool isfinite(float128_t a) {
    return std::isfinite(a.hi) && std::isfinite(a.lo);
}

float128_t floor(float128_t a) {
    long double f = std::floor(a.hi);
    if (a.hi == f && a.lo != 0) {
        if (a.lo < 0) {
            return {f - 1.0L, 0.0L};
        }
    }
    return {f, 0.0L};
}

float128_t ceil(float128_t a) {
    long double c = std::ceil(a.hi);
    if (a.hi == c && a.lo != 0) {
        if (a.lo > 0) {
            return {c + 1.0L, 0.0L};
        }
    }
    return {c, 0.0L};
}

float128_t round(float128_t a) {
    long double r = std::round(a.hi);
    if (a.hi == r && a.lo != 0) {
        if (a.lo >= 0.5L) {
            return {r + 1.0L, 0.0L};
        } else if (a.lo <= -0.5L) {
            return {r - 1.0L, 0.0L};
        }
    }
    return {r, 0.0L};
}

float128_t trunc(float128_t a) {
    return {std::trunc(a.hi), 0.0L};
}

float128_t fmod(float128_t a, float128_t b) {
    float128_t q = a / b;
    float128_t q_trunc = precise128::trunc(q);
    return a - b * q_trunc;
}

float128_t remainder(float128_t a, float128_t b) {
    if (precise128::abs(a).hi > 1e15L) {
        float128_t b_scaled = b * float128_t(1e15L);
        float128_t a_reduced = a;
        while (precise128::abs(a_reduced).hi > precise128::abs(b_scaled).hi) {
            float128_t q = a_reduced / b_scaled;
            float128_t q_trunc = precise128::trunc(q);
            a_reduced = a_reduced - b_scaled * q_trunc;
        }
        float128_t q = a_reduced / b;
        float128_t q_round = precise128::round(q);
        return a_reduced - b * q_round;
    }
    float128_t q = a / b;
    float128_t q_round = precise128::round(q);
    return a - b * q_round;
}

float128_t hypot(float128_t a, float128_t b) {
    float128_t abs_a = precise128::abs(a);
    float128_t abs_b = precise128::abs(b);
    if (abs_a.hi == 0) return abs_b;
    if (abs_b.hi == 0) return abs_a;
    if (abs_a > abs_b) {
        float128_t ratio = abs_b / abs_a;
        return abs_a * precise128::sqrt(float128_t(1.0L) + ratio * ratio);
    }
    float128_t ratio = abs_a / abs_b;
    return abs_b * precise128::sqrt(float128_t(1.0L) + ratio * ratio);
}

float128_t normalize_angle(float128_t x) {
    float128_t two_pi_val = precise128::two_pi();
    float128_t result = precise128::fmod(x, two_pi_val);
    float128_t p = precise128::pi();
    if (result > p) {
        result = result - two_pi_val;
    } else if (result < -p) {
        result = result + two_pi_val;
    }
    return result;
}

// ============================================================================
// Basic math functions
// ============================================================================

float128_t abs(float128_t a) {
    return (a.hi < 0 || (a.hi == 0 && a.lo < 0)) ? -a : a;
}

float128_t sqrt(float128_t a) {
    if (precise128::isnan(a)) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    if (a.hi == 0) return {0, 0};
    if (a.hi < 0) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};

    long double x = std::sqrt(a.hi);
    float128_t res(x);
    for (int i = 0; i < 5; ++i) {
        float128_t next_res = (res + a / res) * 0.5L;
        float128_t diff = precise128::abs(next_res - res);
        if (diff.hi < 1e-35L * precise128::abs(next_res).hi) {
            break;
        }
        res = next_res;
    }
    return res;
}

float128_t cbrt(float128_t a) {
    if (a.hi == 0) return {0, 0};

    long double x = std::cbrt(a.hi);
    float128_t res(x);
    for (int i = 0; i < 5; ++i) {
        float128_t x2 = res * res;
        float128_t x3 = x2 * res;
        float128_t next_res = res - (x3 - a) / (float128_t(3.0L) * x2);
        float128_t diff = precise128::abs(next_res - res);
        if (diff.hi < 1e-35L * precise128::abs(next_res).hi) {
            break;
        }
        res = next_res;
    }
    return res;
}

// High-precision exp using Taylor series and range reduction
float128_t exp_high_prec(float128_t a) {
    if (a.hi > 11356.0L) return {std::numeric_limits<long double>::infinity(), 0.0L};
    if (a.hi < -11399.0L) return {0.0L, 0.0L};

    static const float128_t ln2 = from_string("0.6931471805599453094172321214581765680755");
    static const float128_t inv_ln2 = from_string("1.4426950408889634073599246810018921374266");

    float128_t k_f = precise128::round(a * inv_ln2);
    int k = static_cast<int>(k_f.hi);
    float128_t r = a - k_f * ln2;

    KahanAccumulator accum;
    float128_t term = 1.0L;
    accum.add(term);
    for (int i = 1; i <= 30; ++i) {
        term = term * r / float128_t(static_cast<long double>(i));
        accum.add(term);
        if (precise128::abs(term).hi < 1e-40L) break;
    }

    return precise128::ldexp(accum.get(), k);
}

float128_t exp(float128_t a) {
    return exp_high_prec(a);
}

float128_t ln(float128_t a) {
    if (a.hi <= 0) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};

    long double guess_hi = std::log(a.hi);
    float128_t res(guess_hi);

    for (int i = 0; i < 6; ++i) {
        float128_t ex = exp_high_prec(res);
        float128_t next_res = res + float128_t(2.0L) * (a - ex) / (a + ex);
        float128_t diff = precise128::abs(next_res - res);
        if (diff.hi < 1e-35L * precise128::abs(next_res).hi) {
            break;
        }
        res = next_res;
    }
    return res;
}

float128_t sin(float128_t a) {
    if (precise128::isinf(a) || precise128::isnan(a)) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};

    float128_t x = precise128::remainder(a, precise128::two_pi());

    float128_t half_pi = precise128::half_pi();
    float128_t quarter_pi = half_pi * 0.5L;

    bool negate_result = false;
    if (x.hi > half_pi.hi) {
        x = precise128::pi() - x;
    } else if (x.hi < -half_pi.hi) {
        x = x + precise128::pi();
        negate_result = true;
    }

    bool use_cos = false;
    if (x.hi > quarter_pi.hi) {
        x = half_pi - x;
        use_cos = true;
    } else if (x.hi < -quarter_pi.hi) {
        x = half_pi + x;
        use_cos = true;
        negate_result = !negate_result;
    }

    float128_t result;
    if (use_cos) {
        KahanAccumulator accum;
        float128_t x2 = x * x;
        float128_t term = 1.0L;
        accum.add(term);
        for (int i = 2; i <= 50; i += 2) {
            term = -term * x2 / float128_t(static_cast<long double>(i * (i - 1)));
            accum.add(term);
            if (precise128::abs(term).hi < 1e-40L) break;
        }
        result = accum.get();
    } else {
        KahanAccumulator accum;
        float128_t x2 = x * x;
        float128_t term = x;
        accum.add(term);
        for (int i = 3; i <= 51; i += 2) {
            term = -term * x2 / float128_t(static_cast<long double>(i * (i - 1)));
            accum.add(term);
            if (precise128::abs(term).hi < 1e-40L) break;
        }
        result = accum.get();
    }

    return negate_result ? -result : result;
}

float128_t cos(float128_t a) {
    if (precise128::isinf(a) || precise128::isnan(a)) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};

    float128_t x = precise128::remainder(a, precise128::two_pi());

    float128_t half_pi = precise128::half_pi();
    float128_t quarter_pi = half_pi * 0.5L;

    bool negate_result = false;
    if (x.hi > half_pi.hi) {
        x = precise128::pi() - x;
        negate_result = true;
    } else if (x.hi < -half_pi.hi) {
        x = x + precise128::pi();
        negate_result = true;
    }

    bool use_sin = false;
    if (x.hi > quarter_pi.hi) {
        x = half_pi - x;
        use_sin = true;
    } else if (x.hi < -quarter_pi.hi) {
        x = half_pi + x;
        use_sin = true;
    }

    float128_t result;
    if (use_sin) {
        KahanAccumulator accum;
        float128_t x2 = x * x;
        float128_t term = x;
        accum.add(term);
        for (int i = 3; i <= 51; i += 2) {
            term = -term * x2 / float128_t(static_cast<long double>(i * (i - 1)));
            accum.add(term);
            if (precise128::abs(term).hi < 1e-40L) break;
        }
        result = accum.get();
    } else {
        KahanAccumulator accum;
        float128_t x2 = x * x;
        float128_t term = 1.0L;
        accum.add(term);
        for (int i = 2; i <= 50; i += 2) {
            term = -term * x2 / float128_t(static_cast<long double>(i * (i - 1)));
            accum.add(term);
            if (precise128::abs(term).hi < 1e-40L) break;
        }
        result = accum.get();
    }

    return negate_result ? -result : result;
}

float128_t tan(float128_t a) {
    return precise128::sin(a) / precise128::cos(a);
}

float128_t asin(float128_t a) {
    if (precise128::abs(a).hi > 1.0L) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    if (a.hi == 1.0L && a.lo >= 0) return precise128::pi() * 0.5L;
    if (a.hi == -1.0L && a.lo <= 0) return precise128::pi() * -0.5L;

    float128_t one_minus_x2 = float128_t(1.0L) - a * a;
    return precise128::atan(a / precise128::sqrt(one_minus_x2));
}

float128_t acos(float128_t a) {
    if (precise128::abs(a).hi > 1.0L) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    return precise128::pi() * 0.5L - precise128::asin(a);
}

float128_t sec(float128_t a) {
    return float128_t(1.0L) / precise128::cos(a);
}

float128_t csc(float128_t a) {
    return float128_t(1.0L) / precise128::sin(a);
}

float128_t cot(float128_t a) {
    return precise128::cos(a) / precise128::sin(a);
}

float128_t asec(float128_t a) {
    if (precise128::abs(a).hi < 1.0L) {
        return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    }
    return precise128::acos(float128_t(1.0L) / a);
}

float128_t acsc(float128_t a) {
    if (precise128::abs(a).hi < 1.0L) {
        return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    }
    return precise128::asin(float128_t(1.0L) / a);
}

float128_t acot(float128_t a) {
    return precise128::pi() * 0.5L - precise128::atan(a);
}

// ============================================================================
// Hyperbolic functions
// ============================================================================

float128_t sinh(float128_t a) {
    float128_t e_x = precise128::exp(a);
    float128_t e_neg_x = float128_t(1.0L) / e_x;
    return (e_x - e_neg_x) * 0.5L;
}

float128_t cosh(float128_t a) {
    float128_t e_x = precise128::exp(a);
    float128_t e_neg_x = float128_t(1.0L) / e_x;
    return (e_x + e_neg_x) * 0.5L;
}

float128_t tanh(float128_t a) {
    if (a.hi > 20.0L) {
        return {1.0L, 0.0L};
    }
    if (a.hi < -20.0L) {
        return {-1.0L, 0.0L};
    }
    float128_t e2x = precise128::exp(a * 2.0L);
    return (e2x - float128_t(1.0L)) / (e2x + float128_t(1.0L));
}

float128_t asinh(float128_t a) {
    return precise128::ln(a + precise128::sqrt(a * a + float128_t(1.0L)));
}

float128_t acosh(float128_t a) {
    if (a.hi < 1.0L) {
        return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    }
    return precise128::ln(a + precise128::sqrt(a - float128_t(1.0L)) * precise128::sqrt(a + float128_t(1.0L)));
}

float128_t atanh(float128_t a) {
    if (precise128::abs(a).hi >= 1.0L) {
        return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    }
    return precise128::ln((float128_t(1.0L) + a) / (float128_t(1.0L) - a)) * 0.5L;
}

float128_t log10(float128_t a) {
    static const float128_t ln10 = from_string("2.302585092994045684017991454684364207601");
    return precise128::ln(a) / ln10;
}

float128_t log2(float128_t a) {
    static const float128_t ln2 = from_string("0.6931471805599453094172321214581765680755");
    return precise128::ln(a) / ln2;
}

float128_t log1p(float128_t a) {
    if (precise128::abs(a).hi < 0.125L) {
        KahanAccumulator accum;
        float128_t term = a;
        accum.add(term);
        for (int i = 2; i <= 30; ++i) {
            term = -term * a / float128_t(static_cast<long double>(i));
            accum.add(term);
            if (precise128::abs(term).hi < 1e-40L) break;
        }
        return accum.get();
    }
    return precise128::ln(float128_t(1.0L) + a);
}

float128_t pow(float128_t base, float128_t exponent) {
    if (base.hi == 0) {
        if (exponent.hi == 0 && exponent.lo == 0) {
            return {1.0L, 0.0L};
        }
        if (exponent.hi < 0 || (exponent.hi == 0 && exponent.lo < 0)) {
            return {std::numeric_limits<long double>::infinity(), 0.0L};
        }
        return {0, 0};
    }
    if (base.hi < 0) {
        float128_t exp_rounded = precise128::round(exponent);
        if (precise128::abs(exponent - exp_rounded).hi < 1e-35L) {
            long double exp_int = exp_rounded.hi;
            float128_t abs_base = precise128::abs(base);
            float128_t result = precise128::exp(exp_rounded * precise128::ln(abs_base));
            long long exp_ll = static_cast<long long>(exp_int);
            if (exp_ll % 2 != 0) {
                result = -result;
            }
            return result;
        }
        return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    }
    return precise128::exp(exponent * precise128::ln(base));
}

float128_t atan(float128_t a) {
    long double x = std::atan(a.hi);
    float128_t res(x);
    for (int i = 0; i < 6; ++i) {
        float128_t ty = precise128::tan(res);
        float128_t next_res = res + (a - ty) / (float128_t(1.0L) + ty * ty);
        float128_t diff = precise128::abs(next_res - res);
        if (diff.hi < 1e-35L * (precise128::abs(next_res).hi + 1.0L)) {
            break;
        }
        res = next_res;
    }
    return res;
}

float128_t atan2(float128_t y, float128_t x) {
    if (x.hi > 0) return precise128::atan(y / x);
    if (x.hi < 0) {
        if (y.hi >= 0) return precise128::atan(y / x) + precise128::pi();
        return precise128::atan(y / x) - precise128::pi();
    }
    if (y.hi > 0) return precise128::pi() * 0.5L;
    if (y.hi < 0) return precise128::pi() * -0.5L;
    return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
}

float128_t pi() {
    static const float128_t value = from_string("3.141592653589793238462643383279502884197");
    return value;
}

float128_t two_pi() {
    static const float128_t value = from_string("6.283185307179586476925286766559005768394");
    return value;
}

float128_t half_pi() {
    static const float128_t value = from_string("1.570796326794896619231321691639751442099");
    return value;
}

float128_t sqrt_pi() {
    static const float128_t value = from_string("1.772453850905516027298167483341145182798");
    return value;
}

float128_t e() {
    static const float128_t value = from_string("2.718281828459045235360287471352662497757");
    return value;
}

float128_t infinity() {
    return {std::numeric_limits<long double>::infinity(), 0.0L};
}

// ============================================================================
// Additional math functions
// ============================================================================

float128_t exp2(float128_t a) {
    static const float128_t ln2 = from_string("0.6931471805599453094172321214581765680755");
    return precise128::exp(a * ln2);
}

float128_t expm1(float128_t a) {
    if (precise128::abs(a).hi < 0.125L) {
        KahanAccumulator accum;
        float128_t term = a;
        accum.add(term);
        for (int i = 2; i <= 35; ++i) {
            term = term * a / float128_t(static_cast<long double>(i));
            accum.add(term);
            if (precise128::abs(term).hi < 1e-40L) break;
        }
        return accum.get();
    }
    return precise128::exp(a) - float128_t(1.0L);
}

float128_t ldexp(float128_t a, int exp) {
    if (exp == 0) return a;
    long double scale = std::ldexp(1.0L, exp);
    return a * float128_t(scale);
}

int ilogb(float128_t a) {
    if (a.hi == 0) return FP_ILOGB0;
    if (precise128::isnan(a)) return FP_ILOGBNAN;
    if (precise128::isinf(a)) return INT_MAX;
    return std::ilogb(a.hi);
}

float128_t scalbn(float128_t a, int n) {
    return ldexp(a, n);
}

float128_t logb(float128_t a) {
    if (a.hi == 0) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    if (precise128::isnan(a)) return a;
    if (precise128::isinf(a)) return {std::numeric_limits<long double>::infinity(), 0.0L};
    return float128_t(std::logb(a.hi));
}

float128_t nextafter(float128_t a, float128_t b) {
    if (precise128::isnan(a) || precise128::isnan(b)) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    if (a == b) return b;

    long double next_hi = std::nextafter(a.hi, b.hi);
    return {next_hi, 0.0L};
}

float128_t copysign(float128_t a, float128_t b) {
    if (a.hi == 0) return a;
    bool a_neg = a.hi < 0;
    bool b_neg = b.hi < 0;
    if (a_neg == b_neg) return a;
    return -a;
}

float128_t fma(float128_t a, float128_t b, float128_t c) {
    float128_t p = a * b;
    return p + c;
}

float128_t fma_fast(float128_t a, float128_t b, float128_t c) {
#ifdef __FMA__
    long double result = __builtin_fmal(a.hi, b.hi, c.hi);
    return {result, 0.0L};
#else
    return fma(a, b, c);
#endif
}

float128_t fdim(float128_t a, float128_t b) {
    if (precise128::isnan(a) || precise128::isnan(b)) return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    float128_t diff = a - b;
    if (diff.hi < 0) return {0.0L, 0.0L};
    return diff;
}

float128_t fmax(float128_t a, float128_t b) {
    if (precise128::isnan(a)) return b;
    if (precise128::isnan(b)) return a;
    return (a > b) ? a : b;
}

float128_t fmin(float128_t a, float128_t b) {
    if (precise128::isnan(a)) return b;
    if (precise128::isnan(b)) return a;
    return (a < b) ? a : b;
}

bool is_near_zero(float128_t x, float128_t eps) {
    return precise128::abs(x) <= eps;
}

} // namespace precise128

// ============================================================================
// Stream operators and string conversion
// ============================================================================

std::ostream& operator<<(std::ostream& os, float128_t a) {
    os << to_string(a, os.precision() > 0 ? static_cast<int>(os.precision()) : 32);
    return os;
}

std::istream& operator>>(std::istream& is, float128_t& a) {
    std::string s;
    is >> s;
    a = from_string(s);
    return is;
}

std::string to_string(float128_t a, int precision) {
    if (precise128::isnan(a)) return "nan";
    if (precise128::isinf(a)) return a.hi > 0 ? "inf" : "-inf";
    if (a.hi == 0 && a.lo == 0) return "0";

    std::string s;
    if (a.hi < 0 || (a.hi == 0 && a.lo < 0)) {
        s += "-";
        a = -a;
    }

    float128_t integral = precise128::floor(a);
    float128_t fractional = a - integral;

    if (integral.hi < 1e18L) {
        unsigned long long int_val = static_cast<unsigned long long>(integral.hi);
        s += std::to_string(int_val);
    } else {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(0) << integral.hi;
        s += oss.str();
    }

    if (precision > 0) {
        s += ".";
        for (int i = 0; i < precision; ++i) {
            fractional = fractional * 10.0L;
            float128_t digit_f = precise128::floor(fractional);
            int digit = static_cast<int>(digit_f.hi);
            if (digit > 9) digit = 9;
            if (digit < 0) digit = 0;
            s += std::to_string(digit);
            fractional = fractional - digit_f;
        }
    }

    if (precision > 0) {
        while (s.back() == '0') s.pop_back();
        if (s.back() == '.') s.pop_back();
    }

    return s;
}

float128_t from_string(const std::string& s) {
    if (s == "nan") return {std::numeric_limits<long double>::quiet_NaN(), std::numeric_limits<long double>::quiet_NaN()};
    if (s == "inf" || s == "+inf") return {std::numeric_limits<long double>::infinity(), 0.0L};
    if (s == "-inf") return {-std::numeric_limits<long double>::infinity(), 0.0L};

    float128_t result(0.0L);
    float128_t sign(1.0L);
    size_t i = 0;

    if (i < s.size() && s[i] == '-') {
        sign = -1.0L;
        i++;
    } else if (i < s.size() && s[i] == '+') {
        i++;
    }

    while (i < s.size() && std::isdigit(s[i])) {
        result = result * 10.0L + float128_t(static_cast<long double>(s[i] - '0'));
        i++;
    }

    if (i < s.size() && s[i] == '.') {
        i++;
        float128_t divisor(10.0L);
        while (i < s.size() && std::isdigit(s[i])) {
            result += float128_t(static_cast<long double>(s[i] - '0')) / divisor;
            divisor = divisor * 10.0L;
            i++;
        }
    }

    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
        i++;
        int exp_sign = 1;
        if (i < s.size() && s[i] == '-') {
            exp_sign = -1;
            i++;
        } else if (i < s.size() && s[i] == '+') {
            i++;
        }
        int exp_val = 0;
        while (i < s.size() && std::isdigit(s[i])) {
            exp_val = exp_val * 10 + (s[i] - '0');
            i++;
        }
        if (exp_val != 0) {
            result = result * precise128::pow(float128_t(10.0L), float128_t(exp_val * exp_sign));
        }
    }

    return result * sign;
}

} // namespace mymath

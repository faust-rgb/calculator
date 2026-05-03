#include "integer_helpers.h"
#include "mymath.h"
#include <stdexcept>
#include <sstream>

long long gcd_ll(long long a, long long b) {
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;
    if (a == 0) return b;
    if (b == 0) return a;
    while (b != 0) {
        const long long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

long long lcm_ll(long long a, long long b) {
    if (a == 0 || b == 0) return 0;
    const long long result = (a / gcd_ll(a, b)) * b;
    return result < 0 ? -result : result;
}

bool is_integer_double(double x, double eps) {
    return mymath::is_integer(x, eps);
}

long long round_to_long_long(double x) {
    return static_cast<long long>(x >= 0.0 ? x + 0.5 : x - 0.5);
}

long long trunc_to_long_long(double x) {
    return static_cast<long long>(x);
}

long long floor_to_long_long(double x) {
    long long truncated = static_cast<long long>(x);
    if (x < 0.0 && static_cast<double>(truncated) != x) {
        --truncated;
    }
    return truncated;
}

long long ceil_to_long_long(double x) {
    long long truncated = static_cast<long long>(x);
    if (x > 0.0 && static_cast<double>(truncated) != x) {
        ++truncated;
    }
    return truncated;
}

bool is_prime_ll(long long value) {
    if (value <= 1) return false;
    if (value <= 3) return true;
    if (value % 2 == 0 || value % 3 == 0) return false;
    for (long long i = 5; i * i <= value; i += 6) {
        if (value % i == 0 || value % (i + 2) == 0) return false;
    }
    return true;
}

long long next_prime_ll(long long value) {
    long long candidate = value <= 2 ? 2 : value + 1;
    if (candidate % 2 == 0 && candidate != 2) ++candidate;
    while (!is_prime_ll(candidate)) {
        candidate += (candidate == 2 ? 1 : 2);
    }
    return candidate;
}

long long prev_prime_ll(long long value) {
    if (value <= 2) throw std::runtime_error("prev_prime requires n > 2");
    long long candidate = value - 1;
    if (candidate % 2 == 0) --candidate;
    while (candidate >= 2 && !is_prime_ll(candidate)) {
        candidate -= 2;
    }
    if (candidate < 2) throw std::runtime_error("prev_prime requires n > 2");
    return candidate;
}

long long euler_phi_ll(long long value) {
    if (value <= 0) throw std::runtime_error("euler_phi only accepts positive integers");
    long long n = value;
    long long result = value;
    for (long long p = 2; p * p <= n; ++p) {
        if (n % p != 0) continue;
        while (n % p == 0) n /= p;
        result -= result / p;
    }
    if (n > 1) result -= result / n;
    return result;
}

long long mobius_ll(long long value) {
    if (value <= 0) throw std::runtime_error("mobius only accepts positive integers");
    long long n = value;
    int prime_factor_count = 0;
    for (long long p = 2; p * p <= n; ++p) {
        if (n % p != 0) continue;
        n /= p;
        ++prime_factor_count;
        if (n % p == 0) return 0;
        while (n % p == 0) n /= p;
    }
    if (n > 1) ++prime_factor_count;
    return prime_factor_count % 2 == 0 ? 1 : -1;
}

long long prime_pi_ll(long long value) {
    if (value < 2) return 0;
    long long count = 0;
    for (long long n = 2; n <= value; ++n) {
        if (is_prime_ll(n)) ++count;
    }
    return count;
}

long long extended_gcd_ll(long long a, long long b, long long* x, long long* y) {
    long long old_r = a, r = b;
    long long old_s = 1, s = 0;
    long long old_t = 0, t = 1;

    while (r != 0) {
        const long long quotient = old_r / r;
        const long long next_r = old_r - quotient * r;
        old_r = r; r = next_r;

        const long long next_s = old_s - quotient * s;
        old_s = s; s = next_s;

        const long long next_t = old_t - quotient * t;
        old_t = t; t = next_t;
    }

    if (old_r < 0) {
        old_r = -old_r; old_s = -old_s; old_t = -old_t;
    }
    *x = old_s; *y = old_t;
    return old_r;
}

std::string factor_integer(long long value) {
    if (value == 0) return "0";
    if (value == 1) return "1";
    if (value == -1) return "-1";

    std::ostringstream out;
    bool first = true;
    if (value < 0) {
        out << "-1";
        first = false;
        value = -value;
    }

    for (long long p = 2; p * p <= value; ++p) {
        int exponent = 0;
        while (value % p == 0) {
            value /= p;
            ++exponent;
        }
        if (exponent == 0) continue;
        if (!first) out << " * ";
        first = false;
        out << p;
        if (exponent > 1) out << "^" << exponent;
    }

    if (value > 1) {
        if (!first) out << " * ";
        out << value;
    }

    return out.str();
}

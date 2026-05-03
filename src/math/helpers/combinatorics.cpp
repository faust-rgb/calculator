#include "combinatorics.h"
#include "mymath.h"
#include <stdexcept>
#include <algorithm>

double fibonacci_value(long long n) {
    if (n < 0) throw std::runtime_error("fib only accepts non-negative integers");
    if (n > 186) throw std::runtime_error("fib is limited to n <= 186 to avoid overflow");
    if (n == 0) return 0.0;
    long long a = 0, b = 1;
    for (long long i = 1; i < n; ++i) {
        const long long next = a + b;
        a = b; b = next;
    }
    return static_cast<double>(b);
}

double factorial_value(long long n) {
    if (n < 0) throw std::runtime_error("factorial only accepts non-negative integers");
    if (n > 170) throw std::runtime_error("factorial is limited to n <= 170 to avoid overflow");
    long double result = 1.0L;
    for (long long i = 2; i <= n; ++i) {
        result *= static_cast<long double>(i);
    }
    return static_cast<double>(result);
}

Rational factorial_rational(long long n) {
    return Rational(static_cast<long long>(factorial_value(n)), 1);
}

double combination_value(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("combination requires 0 <= r <= n");
    if (n > 170) throw std::runtime_error("nCr is limited to n <= 170 to avoid overflow");
    r = std::min(r, n - r);
    long double result = 1.0L;
    for (long long i = 1; i <= r; ++i) {
        result *= static_cast<long double>(n - r + i);
        result /= static_cast<long double>(i);
    }
    return static_cast<double>(result);
}

Rational combination_rational(long long n, long long r) {
    return Rational(static_cast<long long>(combination_value(n, r)), 1);
}

double permutation_value(long long n, long long r) {
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("permutation requires 0 <= r <= n");
    if (n > 170) throw std::runtime_error("nPr is limited to n <= 170 to avoid overflow");
    long double result = 1.0L;
    for (long long i = 0; i < r; ++i) {
        result *= static_cast<long double>(n - i);
    }
    return static_cast<double>(result);
}

Rational permutation_rational(long long n, long long r) {
    return Rational(static_cast<long long>(permutation_value(n, r)), 1);
}

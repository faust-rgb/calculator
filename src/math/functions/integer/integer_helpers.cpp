/**
 * @file integer_helpers.cpp
 * @brief 整数辅助函数实现
 *
 * 本文件提供整数相关的辅助函数，包括：
 * - 最大公约数和最小公倍数
 * - 浮点数到整数的各种转换方式
 * - 素数判断和素数查找
 * - 数论函数（欧拉函数、莫比乌斯函数、素数计数函数）
 * - 扩展欧几里得算法
 * - 整数因式分解
 */

#include "integer_helpers.h"
#include "mymath.h"
#include <stdexcept>
#include <sstream>

namespace mymath {

int abs(int x) { return x < 0 ? -x : x; }
long abs(long x) { return x < 0 ? -x : x; }
long long abs(long long x) { return x < 0 ? -x : x; }

long long gcd(long long a, long long b) {
    while (b != 0) {
        const long long remainder = a % b;
        a = b;
        b = remainder;
    }
    return a < 0 ? -a : a;
}

}  // namespace mymath

/**
 * @brief 计算两个整数的最大公约数（GCD）
 * @param a 第一个整数
 * @param b 第二个整数
 * @return a 和 b 的最大公约数
 *
 * 使用欧几里得算法实现，自动处理负数参数。
 */
long long gcd_ll(long long a, long long b) {
    constexpr long long kLongLongMax = 9223372036854775807LL;
    unsigned long long ua = a < 0 ? 0ULL - static_cast<unsigned long long>(a)
                                  : static_cast<unsigned long long>(a);
    unsigned long long ub = b < 0 ? 0ULL - static_cast<unsigned long long>(b)
                                  : static_cast<unsigned long long>(b);
    while (ub != 0) {
        const unsigned long long t = ua % ub;
        ua = ub;
        ub = t;
    }
    if (ua > static_cast<unsigned long long>(kLongLongMax)) {
        throw std::overflow_error("gcd exceeds signed range");
    }
    return static_cast<long long>(ua);
}

/**
 * @brief 计算两个整数的最小公倍数（LCM）
 * @param a 第一个整数
 * @param b 第二个整数
 * @return a 和 b 的最小公倍数
 *
 * 使用公式 LCM(a, b) = |a * b| / GCD(a, b) 计算。
 */
long long lcm_ll(long long a, long long b) {
    if (a == 0 || b == 0) return 0;
    const long long result = (a / gcd_ll(a, b)) * b;
    return result < 0 ? -result : result;
}

/**
 * @brief 检查 long double 值是否为整数
 * @param x 待检查的值
 * @param eps 允许的误差阈值
 * @return true 如果 x 接近整数
 */
bool is_integer_double(long double x, long double eps) {
    return mymath::is_integer(x, eps);
}

/**
 * @brief 将 long double 四舍五入为 long long
 * @param x 待转换的值
 * @return 四舍五入后的整数值
 */
long long round_to_long_long(long double x) {
    return static_cast<long long>(x >= 0.0L ? x + 0.5 : x - 0.5);
}

/**
 * @brief 将 long double 向零截断为 long long
 * @param x 待转换的值
 * @return 向零截断后的整数值
 */
long long trunc_to_long_long(long double x) {
    return static_cast<long long>(x);
}

/**
 * @brief 将 long double 向下取整为 long long
 * @param x 待转换的值
 * @return 不大于 x 的最大整数
 */
long long floor_to_long_long(long double x) {
    long long truncated = static_cast<long long>(x);
    if (x < 0.0L && static_cast<long double>(truncated) != x) {
        --truncated;
    }
    return truncated;
}

/**
 * @brief 将 long double 向上取整为 long long
 * @param x 待转换的值
 * @return 不小于 x 的最小整数
 */
long long ceil_to_long_long(long double x) {
    long long truncated = static_cast<long long>(x);
    if (x > 0.0L && static_cast<long double>(truncated) != x) {
        ++truncated;
    }
    return truncated;
}

/**
 * @brief 判断一个整数是否为素数
 * @param value 待判断的值
 * @return true 如果是素数，false 否则
 *
 * 使用试除法，检查到 sqrt(value) 为止。
 */
bool is_prime_ll(long long value) {
    if (value <= 1) return false;
    if (value <= 3) return true;
    if (value % 2 == 0 || value % 3 == 0) return false;
    for (long long i = 5; i * i <= value; i += 6) {
        if (value % i == 0 || value % (i + 2) == 0) return false;
    }
    return true;
}

/**
 * @brief 找到大于给定值的最小素数
 * @param value 起始值
 * @return 大于 value 的最小素数
 */
long long next_prime_ll(long long value) {
    long long candidate = value <= 2 ? 2 : value + 1;
    if (candidate % 2 == 0 && candidate != 2) ++candidate;
    while (!is_prime_ll(candidate)) {
        candidate += (candidate == 2 ? 1 : 2);
    }
    return candidate;
}

/**
 * @brief 找到小于给定值的最大素数
 * @param value 起始值
 * @return 小于 value 的最大素数
 * @throws std::runtime_error 如果 value <= 2
 */
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

/**
 * @brief 计算欧拉函数 φ(n)
 * @param value 正整数
 * @return 小于等于 n 且与 n 互质的正整数个数
 * @throws std::runtime_error 如果 value 不是正整数
 *
 * 欧拉函数 φ(n) 表示 1 到 n 中与 n 互质的整数个数。
 */
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

/**
 * @brief 计算莫比乌斯函数 μ(n)
 * @param value 正整数
 * @return μ(n) 的值（-1, 0, 或 1）
 * @throws std::runtime_error 如果 value 不是正整数
 *
 * 莫比乌斯函数定义：
 * - μ(n) = 1，如果 n 是偶数个不同素数的乘积
 * - μ(n) = -1，如果 n 是奇数个不同素数的乘积
 * - μ(n) = 0，如果 n 有平方因子
 */
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

/**
 * @brief 计算素数计数函数 π(n)
 * @param value 非负整数
 * @return 小于等于 n 的素数个数
 */
long long prime_pi_ll(long long value) {
    if (value < 2) return 0;
    long long count = 0;
    for (long long n = 2; n <= value; ++n) {
        if (is_prime_ll(n)) ++count;
    }
    return count;
}

/**
 * @brief 扩展欧几里得算法
 * @param a 第一个整数
 * @param b 第二个整数
 * @param x 输出参数，存储贝祖系数 x
 * @param y 输出参数，存储贝祖系数 y
 * @return a 和 b 的最大公约数
 *
 * 计算 gcd(a, b) 并找到整数 x, y 满足 ax + by = gcd(a, b)。
 */
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

/**
 * @brief 整数因式分解
 * @param value 待分解的整数
 * @return 因式分解的字符串表示（如 "2^2 * 3 * 5"）
 *
 * 将整数分解为素因子的乘积形式，返回易读的字符串表示。
 */
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

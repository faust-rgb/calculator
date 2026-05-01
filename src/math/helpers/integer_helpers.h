#ifndef INTEGER_HELPERS_H
#define INTEGER_HELPERS_H

#include <string>

/** @brief 最大公约数 */
long long gcd_ll(long long a, long long b);

/** @brief 最小公倍数 */
long long lcm_ll(long long a, long long b);

/** @brief 检查 double 是否为整数 */
bool is_integer_double(double x, double eps = 1e-10);

/** @brief 四舍五入 */
long long round_to_long_long(double x);

/** @brief 向零截断 */
long long trunc_to_long_long(double x);

/** @brief 向下取整 */
long long floor_to_long_long(double x);

/** @brief 向上取整 */
long long ceil_to_long_long(double x);

/** @brief 素数判断 */
bool is_prime_ll(long long value);

/** @brief 下一个素数 */
long long next_prime_ll(long long value);

/** @brief 上一个素数 */
long long prev_prime_ll(long long value);

/** @brief 欧拉函数 */
long long euler_phi_ll(long long value);

/** @brief 莫比乌斯函数 */
long long mobius_ll(long long value);

/** @brief 素数计数函数 */
long long prime_pi_ll(long long value);

/** @brief 扩展欧几里得 */
long long extended_gcd_ll(long long a, long long b, long long* x, long long* y);

/** @brief 整数因式分解 */
std::string factor_integer(long long value);

#endif

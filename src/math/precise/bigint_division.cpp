// ============================================================================
// 大整数除法算法实现
// ============================================================================
//
// 本文件实现高效的大整数除法算法：
// - Knuth 算法 D (经典除法)
// - Newton-Raphson 除法 (大规模优化)
// - Barrett 约减 (重复除法优化)
// ============================================================================

#include "precise_decimal.h"

#include <algorithm>
#include <stdexcept>

namespace precise {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;
constexpr int kBaseDigits = 9;

// 前向声明
void div_bigint(const BigIntData& num, const BigIntData& den, BigIntData* quotient, BigIntData* remainder);

// ============================================================================
// 辅助函数
// ============================================================================

int count_trailing_zeros(const BigIntData& v) {
    if (v.empty() || (v.size() == 1 && v[0] == 0)) return 0;
    int count = 0;
    std::size_t i = 0;
    while (i < v.size() && v[i] == 0) {
        count += kBaseDigits;
        i++;
    }
    if (i < v.size()) {
        uint32_t val = v[i];
        while (val > 0 && val % 10 == 0) {
            count++;
            val /= 10;
        }
    }
    return count;
}

BigIntData divide_bigint_by_pow10(BigIntData v, int n) {
    if (n <= 0) return v;
    int chunk_shift = n / kBaseDigits;
    int digit_shift = n % kBaseDigits;

    if (chunk_shift >= static_cast<int>(v.size())) return {0};

    if (chunk_shift > 0) {
        v.erase(v.begin(), v.begin() + chunk_shift);
    }

    if (digit_shift > 0) {
        uint32_t divisor = 1;
        for (int i = 0; i < digit_shift; ++i) divisor *= 10;

        uint32_t rem = 0;
        for (int i = static_cast<int>(v.size()) - 1; i >= 0; --i) {
            uint64_t current = v[i] + static_cast<uint64_t>(rem) * kBase;
            v[i] = static_cast<uint32_t>(current / divisor);
            rem = static_cast<uint32_t>(current % divisor);
        }
    }
    while (v.size() > 1 && v.back() == 0) v.pop_back();
    if (v.empty()) v = {0};
    return v;
}

bool is_bigint_multiple_of_pow10(const BigIntData& v, int n) {
    if (n <= 0) return true;
    return count_trailing_zeros(v) >= n;
}

BigIntData multiply_bigint_by_power_of_10(BigIntData v, int n) {
    if (n <= 0) return v;
    if (v.size() == 1 && v[0] == 0) return v;

    int chunk_shift = n / kBaseDigits;
    int digit_shift = n % kBaseDigits;

    if (digit_shift > 0) {
        uint32_t multiplier = 1;
        for (int i = 0; i < digit_shift; ++i) multiplier *= 10;
        v = multiply_bigint_by_uint32(v, multiplier);
    }

    if (chunk_shift > 0) {
        v.insert(v.begin(), static_cast<std::size_t>(chunk_shift), 0u);
    }
    return v;
}

// ============================================================================
// Newton-Raphson 除法
// ============================================================================

BigIntData reciprocal_newton(const BigIntData& d, std::size_t precision_chunks) {
    if (d.empty() || (d.size() == 1 && d[0] == 0)) {
        throw std::runtime_error("division by zero in reciprocal_newton");
    }

    uint64_t d_high = 0;
    if (d.size() >= 2) {
        d_high = (static_cast<uint64_t>(d[d.size() - 1]) << 32) | d[d.size() - 2];
    } else {
        d_high = d[0];
    }

    uint64_t x0_val = (static_cast<uint64_t>(1) << 63) / (d_high >> 1);
    if (x0_val == 0) x0_val = 1;

    BigIntData x = {static_cast<uint32_t>(x0_val), static_cast<uint32_t>(x0_val >> 32)};

    std::size_t current_precision = 2;

    while (current_precision < precision_chunks + 2) {
        std::size_t needed_precision = std::min(current_precision * 2, precision_chunks + 2);

        std::size_t d_truncate = std::min(d.size(), needed_precision + 2);
        BigIntData d_trunc(d.begin(), d.begin() + static_cast<std::ptrdiff_t>(d_truncate));

        BigIntData d_times_x = multiply_bigint(d_trunc, x);
        if (d_times_x.size() > needed_precision + current_precision) {
            d_times_x.resize(needed_precision + current_precision);
        }

        BigIntData two_pow_2k(2 * current_precision + 1, 0);
        two_pow_2k[2 * current_precision] = 1;

        BigIntData correction;
        if (compare_bigint(two_pow_2k, d_times_x) > 0) {
            correction = subtract_bigint(two_pow_2k, d_times_x);
        } else {
            correction = subtract_bigint(d_times_x, two_pow_2k);
            current_precision *= 2;
            continue;
        }

        BigIntData x_new = multiply_bigint(x, correction);

        if (x_new.size() > current_precision) {
            x.assign(x_new.begin() + static_cast<std::ptrdiff_t>(current_precision), x_new.end());
        } else {
            x = {0};
        }

        while (x.size() < needed_precision) {
            x.push_back(0);
        }

        current_precision = needed_precision;
    }

    if (x.size() > precision_chunks) {
        x.resize(precision_chunks);
    }

    return x;
}

// ============================================================================
// Barrett 约减器
// ============================================================================

struct BarrettReducer {
    BigIntData divisor;
    BigIntData reciprocal;
    std::size_t shift_chunks;

    explicit BarrettReducer(const BigIntData& d, bool use_newton = false) : divisor(d) {
        if (d.empty() || (d.size() == 1 && d[0] == 0)) {
            throw std::runtime_error("BarrettReducer: divisor cannot be zero");
        }

        std::size_t n = divisor.size();
        shift_chunks = 2 * n;

        if (use_newton || n > 256) {
            reciprocal = reciprocal_newton(divisor, shift_chunks);
        } else {
            BigIntData two_pow_2n(shift_chunks + 1, 0);
            two_pow_2n[shift_chunks] = 1;

            BigIntData rem;
            div_bigint(two_pow_2n, divisor, &reciprocal, &rem);
        }
    }

    BigIntData reduce(const BigIntData& x) const {
        if (compare_bigint(x, divisor) < 0) {
            return x;
        }

        BigIntData x_times_r = multiply_bigint(x, reciprocal);

        if (shift_chunks >= x_times_r.size()) {
            return x;
        }

        BigIntData q(x_times_r.begin() + static_cast<std::ptrdiff_t>(shift_chunks), x_times_r.end());

        BigIntData q_times_d = multiply_bigint(q, divisor);
        BigIntData r;
        if (compare_bigint(x, q_times_d) >= 0) {
            r = subtract_bigint(x, q_times_d);
        } else {
            q = subtract_bigint(q, BigIntData{1});
            q_times_d = multiply_bigint(q, divisor);
            r = subtract_bigint(x, q_times_d);
        }

        while (compare_bigint(r, divisor) >= 0) {
            r = subtract_bigint(r, divisor);
        }

        return r;
    }

    BigIntData divide(const BigIntData& x) const {
        if (compare_bigint(x, divisor) < 0) {
            return {0};
        }

        BigIntData x_times_r = multiply_bigint(x, reciprocal);

        if (shift_chunks >= x_times_r.size()) {
            return {0};
        }

        BigIntData q(x_times_r.begin() + static_cast<std::ptrdiff_t>(shift_chunks), x_times_r.end());

        BigIntData q_times_d = multiply_bigint(q, divisor);
        BigIntData r;
        if (compare_bigint(x, q_times_d) >= 0) {
            r = subtract_bigint(x, q_times_d);
        } else {
            q = subtract_bigint(q, BigIntData{1});
            q_times_d = multiply_bigint(q, divisor);
        }

        while (compare_bigint(multiply_bigint(add_bigint(q, BigIntData{1}), divisor), x) <= 0) {
            q = add_bigint(q, BigIntData{1});
        }

        return q;
    }
};

// ============================================================================
// Newton-Raphson 除法入口
// ============================================================================

void div_bigint_newton(const BigIntData& num,
                       const BigIntData& den,
                       BigIntData* quotient,
                       BigIntData* remainder) {
    if (den.empty() || (den.size() == 1 && den[0] == 0)) {
        throw std::runtime_error("division by zero in div_bigint_newton");
    }

    constexpr std::size_t NEWTON_THRESHOLD = 256;

    if (den.size() < NEWTON_THRESHOLD || num.size() < NEWTON_THRESHOLD) {
        div_bigint(num, den, quotient, remainder);
        return;
    }

    BarrettReducer reducer(den, true);
    *quotient = reducer.divide(num);
    *remainder = reducer.reduce(num);
}

// ============================================================================
// Knuth 算法 D (经典除法)
// ============================================================================

void div_bigint(const BigIntData& num,
                const BigIntData& den,
                BigIntData* quotient,
                BigIntData* remainder) {
    std::size_t m = den.size();
    while (m > 0 && den[m - 1] == 0) --m;
    if (m == 0) {
        throw std::runtime_error("division by zero in div_bigint");
    }

    if (num.empty() || (num.size() == 1 && num[0] == 0)) {
        *quotient = {0};
        *remainder = {0};
        return;
    }

    std::size_t n = num.size();
    while (n > 1 && num[n - 1] == 0) --n;

    if (m > n || (m == n && compare_bigint(
            BigIntData(num.begin(), num.begin() + static_cast<std::ptrdiff_t>(n)),
            BigIntData(den.begin(), den.begin() + static_cast<std::ptrdiff_t>(m))) < 0)) {
        *quotient = {0};
        *remainder = BigIntData(num.begin(), num.begin() + static_cast<std::ptrdiff_t>(n));
        while (remainder->size() > 1 && remainder->back() == 0) remainder->pop_back();
        return;
    }

    if (m == 1) {
        BigIntData q(n, 0);
        uint64_t rem = 0;
        uint32_t d = den[0];
        for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
            uint64_t cur = rem * kBase + num[i];
            q[i] = static_cast<uint32_t>(cur / d);
            rem = cur % d;
        }
        while (q.size() > 1 && q.back() == 0) q.pop_back();
        *quotient = std::move(q);
        *remainder = {static_cast<uint32_t>(rem)};
        return;
    }

    uint32_t d = static_cast<uint32_t>(kBase / (static_cast<uint64_t>(den[m - 1]) + 1));
    if (d == 0) d = 1;

    BigIntData u(n + 1, 0);
    BigIntData v(m, 0);

    uint64_t carry = 0;
    for (std::size_t i = 0; i < n; ++i) {
        uint64_t cur = static_cast<uint64_t>(num[i]) * d + carry;
        u[i] = static_cast<uint32_t>(cur % kBase);
        carry = cur / kBase;
    }
    u[n] = static_cast<uint32_t>(carry);

    carry = 0;
    for (std::size_t i = 0; i < m; ++i) {
        uint64_t cur = static_cast<uint64_t>(den[i]) * d + carry;
        v[i] = static_cast<uint32_t>(cur % kBase);
        carry = cur / kBase;
    }

    BigIntData q(n - m + 1, 0);

    for (int j = static_cast<int>(n - m); j >= 0; --j) {
        uint64_t numerator_val = static_cast<uint64_t>(u[j + m]) * kBase + u[j + m - 1];
        uint64_t q_hat = numerator_val / v[m - 1];
        uint64_t r_hat = numerator_val % v[m - 1];

        while (q_hat >= kBase ||
               (q_hat * v[m - 2] > r_hat * kBase + (j + static_cast<int>(m) - 2 >= 0 ? u[j + m - 2] : 0))) {
            --q_hat;
            r_hat += v[m - 1];
            if (r_hat >= kBase) break;
        }

        int64_t borrow = 0;
        for (std::size_t i = 0; i < m; ++i) {
            uint64_t product = q_hat * v[i];
            int64_t diff = static_cast<int64_t>(u[j + i]) - static_cast<int64_t>(product % kBase) - borrow;
            borrow = static_cast<int64_t>(product / kBase);
            if (diff < 0) {
                diff += kBase;
                ++borrow;
            }
            u[j + i] = static_cast<uint32_t>(diff);
        }
        int64_t diff = static_cast<int64_t>(u[j + m]) - borrow;
        u[j + m] = static_cast<uint32_t>(diff);

        q[j] = static_cast<uint32_t>(q_hat);
        if (diff < 0) {
            --q[j];
            uint64_t c = 0;
            for (std::size_t i = 0; i < m; ++i) {
                uint64_t sum = static_cast<uint64_t>(u[j + i]) + v[i] + c;
                u[j + i] = static_cast<uint32_t>(sum % kBase);
                c = sum / kBase;
            }
            u[j + m] += static_cast<uint32_t>(c);
        }
    }

    while (q.size() > 1 && q.back() == 0) q.pop_back();
    *quotient = std::move(q);

    BigIntData rem(m, 0);
    uint64_t r = 0;
    for (int i = static_cast<int>(m) - 1; i >= 0; --i) {
        uint64_t cur = r * kBase + u[i];
        rem[i] = static_cast<uint32_t>(cur / d);
        r = cur % d;
    }
    while (rem.size() > 1 && rem.back() == 0) rem.pop_back();
    *remainder = std::move(rem);
}

} // namespace precise
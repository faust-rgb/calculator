// ============================================================================
// 快速除法与开方优化 (Newton-Raphson + Barrett)
// ============================================================================
//
// 本文件实现高效的除法和开方算法：
// - 使用硬件浮点数计算初值（提供 15-18 位有效数字）
// - Barrett 约减用于重复除法
// - Newton-Raphson 迭代优化
//
// 优化目标：
// - 除法：O(M(n)) 而非 O(n²)
// - 开方：减少 2-3 次迭代
// ============================================================================

#include "precise_decimal.h"
#include "core/common/calculator_exceptions.h"
#include "math/mymath.h"

#include <algorithm>
#include <cstdint>
#include <memory>

namespace precise {

// ============================================================================
// 常量定义
// ============================================================================

constexpr uint32_t kBase = 1000000000;
constexpr uint64_t kBase64 = 1000000000ULL;
constexpr int kBaseDigits = 9;

static inline PreciseDecimal make_power_of_10(int exponent) {
    PreciseDecimal result(1LL);
    if (exponent >= 0) {
        result.data = multiply_bigint_by_power_of_10(result.data, exponent);
    } else {
        result.scale = -exponent;
    }
    result.normalize();
    return result;
}

// 提取大数的高位尾数及十进制指数（避免极值尺度下的浮点溢出）
static void extract_mantissa_and_exponent(const BigIntData& data, int scale, long double* out_mantissa, int* out_exponent) {
    if (data.empty() || (data.size() == 1 && data[0] == 0)) {
        *out_mantissa = 0.0L;
        *out_exponent = 0;
        return;
    }
    std::string digits = bigint_to_string(data);
    int total_digits = static_cast<int>(digits.size());
    int exponent = total_digits - 1 - scale;
    int mantissa_len = std::min<int>(18, total_digits);
    std::string mantissa_str = digits.substr(0, static_cast<std::size_t>(mantissa_len));
    long double mantissa = std::stold(mantissa_str);
    for (int i = 1; i < mantissa_len; ++i) {
        mantissa /= 10.0L;
    }
    *out_mantissa = mantissa;
    *out_exponent = exponent;
}

// ============================================================================
// 硬件浮点初值计算
// ============================================================================
//
// 使用 long double 提供约 18 位有效数字的初值
// 这可以将 Newton-Raphson 迭代次数减少 2-3 次
// ============================================================================

// 计算 sqrt(x) 的硬件初值
long double compute_sqrt_initial_guess(const BigIntData& data, int scale) {
    long double mantissa = 0.0L;
    int exp = 0;
    extract_mantissa_and_exponent(data, scale, &mantissa, &exp);
    if (mantissa <= 0.0L) {
        return 0.0L;
    }

    if (exp % 2 != 0) {
        mantissa *= 10.0L;
        --exp;
    }

    if (exp > 300 || exp < -300) {
        return mymath::sqrt(mantissa);
    }

    long double approx = mymath::sqrt(mantissa) * mymath::pow(10.0L, static_cast<long double>(exp / 2));
    if (!mymath::isfinite(approx) || approx <= 0.0L) {
        return 1.0L;
    }
    return approx;
}

// 计算 1/x 的硬件初值
long double compute_reciprocal_initial_guess(const BigIntData& data, int scale) {
    long double mantissa = 0.0L;
    int exp = 0;
    extract_mantissa_and_exponent(data, scale, &mantissa, &exp);
    if (mantissa <= 0.0L) {
        return 0.0L;
    }

    if (exp > 300 || exp < -300) {
        return 1.0L / mantissa;
    }

    long double approx = (1.0L / mantissa) * mymath::pow(10.0L, static_cast<long double>(-exp));
    if (!mymath::isfinite(approx) || approx <= 0.0L) {
        return 1.0L;
    }
    return approx;
}

// ============================================================================
// Barrett 约减器（优化版）
// ============================================================================
//
// 用于对同一除数进行多次除法的优化
// 预计算准倒数，将除法转换为两次乘法和一次位移
// ============================================================================

class BarrettReducerOptimized {
public:
    // 构造函数：预计算准倒数
    explicit BarrettReducerOptimized(const BigIntData& divisor, int divisor_scale = 0)
        : divisor_(divisor), divisor_scale_(divisor_scale) {

        if (divisor.empty() || (divisor.size() == 1 && divisor[0] == 0)) {
            valid_ = false;
            return;
        }

        valid_ = true;
        n_ = divisor.size();

        // 计算准倒数：floor(2^(2n*kBaseDigits) / divisor)
        // 使用 Newton-Raphson 方法计算倒数

        // 首先获取硬件初值
        long double approx = compute_reciprocal_initial_guess(divisor, divisor_scale);

        // 将初值转换为大整数形式
        // 准倒数的精度需要 2n 个 limb
        size_t target_precision = 2 * n_ + 2;

        // 从硬件初值构建初始近似
        reciprocal_ = long_double_to_bigint(approx, target_precision);

        // 使用 Newton-Raphson 精化倒数
        // x_{k+1} = x_k * (2 - d * x_k)
        refine_reciprocal(target_precision);
    }

    // 快速除法
    BigIntData divide(const BigIntData& num, int num_scale, int target_scale) const {
        if (!valid_) {
            throw std::runtime_error("Barrett division by zero");
        }

        // 对齐小数点
        int scale_diff = num_scale - divisor_scale_;
        BigIntData adjusted_num = num;

        if (scale_diff > 0) {
            // num 有更多小数位，需要放大
            adjusted_num = multiply_bigint_by_power_of_10(adjusted_num, scale_diff);
        } else if (scale_diff < 0) {
            // divisor 有更多小数位，需要调整目标精度
            target_scale -= scale_diff;
        }

        // 使用 Barrett 约减计算商
        return barrett_divide(adjusted_num, target_scale);
    }

    // 快速取模
    BigIntData mod(const BigIntData& num) const {
        if (!valid_) {
            throw std::runtime_error("Barrett mod by zero");
        }

        BigIntData q = barrett_divide(num, 0);
        BigIntData q_times_d = multiply_bigint(q, divisor_);
        BigIntData r = subtract_bigint(num, q_times_d);

        // 规范化
        while (r.size() > 1 && r.back() == 0) {
            r.pop_back();
        }

        return r;
    }

    bool is_valid() const { return valid_; }
    const BigIntData& get_divisor() const { return divisor_; }

private:
    BigIntData divisor_;
    int divisor_scale_;
    BigIntData reciprocal_;
    size_t n_;
    bool valid_;

    // 将 long double 转换为大整数
    BigIntData long_double_to_bigint(long double val, size_t precision) {
        BigIntData result;
        result.reserve(precision);

        // 提取整数部分
        if (val < 1.0L) {
            // 对于小于 1 的数，需要计算 1/val 的精度
            val = 1.0L / val;
            // 存储为倒数形式（简化处理）
        }

        // 将 val 转换为整数形式
        // 乘以足够的 10 的幂次
        int exponent = static_cast<int>(precision) * kBaseDigits;
        long double scaled = val * mymath::pow(10.0L, static_cast<long double>(exponent));

        // 转换为 limbs
        while (scaled > 0.5L && result.size() < precision) {
            long double digit = mymath::fmod(scaled, static_cast<long double>(kBase));
            result.push_back(static_cast<uint32_t>(digit + 0.5L));
            scaled /= kBase;
            scaled = mymath::floor(scaled);
        }

        if (result.empty()) {
            result.push_back(0);
        }

        return result;
    }

    // 使用 Newton-Raphson 精化倒数
    void refine_reciprocal(size_t target_precision) {
        // 迭代公式：x_{k+1} = x_k * (2 - d * x_k)
        // 收敛速度：每次迭代精度翻倍

        size_t current_precision = 2;  // 初始精度（来自硬件初值）
        int max_iterations = 10;  // 最大迭代次数

        for (int iter = 0; iter < max_iterations && current_precision < target_precision; ++iter) {
            // 计算 d * x_k
            BigIntData d_times_x = multiply_bigint(divisor_, reciprocal_);

            // 计算 2 - d * x_k（在适当的精度下）
            // 注意：这里需要处理符号
            BigIntData two_pow_p(current_precision * 2 + 1, 0);
            two_pow_p[current_precision * 2] = 1;  // 2^p 的近似

            BigIntData correction;
            if (compare_bigint(two_pow_p, d_times_x) > 0) {
                correction = subtract_bigint(two_pow_p, d_times_x);
            } else {
                // 精度已足够，停止迭代
                break;
            }

            // 计算 x_{k+1} = x_k * correction
            BigIntData new_reciprocal = multiply_bigint(reciprocal_, correction);

            // 截断到目标精度
            if (new_reciprocal.size() > target_precision) {
                reciprocal_.assign(
                    new_reciprocal.begin() + static_cast<std::ptrdiff_t>(new_reciprocal.size() - target_precision),
                    new_reciprocal.end()
                );
            } else {
                reciprocal_ = std::move(new_reciprocal);
            }

            current_precision *= 2;
        }

        // 规范化
        while (reciprocal_.size() > 1 && reciprocal_.back() == 0) {
            reciprocal_.pop_back();
        }
    }

    // Barrett 除法核心
    BigIntData barrett_divide(const BigIntData& num, int target_scale) const {
        // q ≈ floor(num * reciprocal / 2^(2n*kBaseDigits))
        (void)target_scale;  // 目前未使用，未来可用于调整结果精度
        BigIntData num_times_r = multiply_bigint(num, reciprocal_);

        // 右移 2n 个 limb
        size_t shift = 2 * n_;
        if (num_times_r.size() > shift) {
            BigIntData q(num_times_r.begin() + static_cast<std::ptrdiff_t>(shift),
                         num_times_r.end());

            // 修正商
            BigIntData q_times_d = multiply_bigint(q, divisor_);
            BigIntData r;

            if (compare_bigint(num, q_times_d) >= 0) {
                r = subtract_bigint(num, q_times_d);
            } else {
                // 商太大，减 1
                q = subtract_bigint(q, BigIntData{1});
                q_times_d = multiply_bigint(q, divisor_);
                r = subtract_bigint(num, q_times_d);
            }

            // 最终修正
            while (compare_bigint(r, divisor_) >= 0) {
                q = add_bigint(q, BigIntData{1});
                r = subtract_bigint(r, divisor_);
            }

            return q;
        }

        return {0};
    }
};

// ============================================================================
// 优化的 sqrt 实现
// ============================================================================
//
// 使用硬件初值 + Newton-Raphson 迭代
// 初值精度约 18 位，可减少 2-3 次迭代
// ============================================================================

PreciseDecimal sqrt_optimized(const PreciseDecimal& val) {
    if (val.is_zero()) {
        return PreciseDecimal(0LL);
    }
    if (val.negative) {
        throw MathError("sqrt of negative number");
    }

    int target_scale = PrecisionContext::get_default_scale();

    // 基于尾数和指数构建安全初值（防止超大/超小尺度溢出）
    long double mantissa = 0.0L;
    int exp = 0;
    extract_mantissa_and_exponent(val.data, val.scale, &mantissa, &exp);
    if (exp % 2 != 0) {
        mantissa *= 10.0L;
        --exp;
    }
    long double sqrt_m = (mantissa > 0.0L) ? mymath::sqrt(mantissa) : 1.0L;
    PreciseDecimal x = PreciseDecimal(sqrt_m) * make_power_of_10(exp / 2);
    if (x.is_zero()) {
        x = make_power_of_10(exp / 2);
    }

    const PreciseDecimal one_half("0.5");
    const PreciseDecimal epsilon = PreciseDecimal(
        "1e-" + std::to_string(target_scale + 8));

    NormalizationSuppressor suppressor;

    // Newton-Raphson 迭代：x_{k+1} = 0.5 * (x_k + val / x_k)
    // 由于初值有约 18 位精度，迭代次数大大减少
    int max_iterations = std::max(5, target_scale / 20 + 2);

    for (int i = 0; i < max_iterations; ++i) {
        PreciseDecimal next = one_half * (x + val / x);

        PreciseDecimal diff = precise::abs(next - x);
        if (diff.is_zero() || diff < epsilon * std::max(PreciseDecimal(1LL), precise::abs(next))) {
            g_suppress_normalization = false;
            next.normalize();
            return next;
        }

        x = next;
    }

    g_suppress_normalization = false;
    x.normalize();
    return x;
}

// ============================================================================
// 优化的倒数实现
// ============================================================================

PreciseDecimal reciprocal_optimized(const PreciseDecimal& val) {
    if (val.is_zero()) {
        throw std::runtime_error("reciprocal of zero");
    }

    int target_scale = PrecisionContext::get_default_scale();

    // 基于尾数和指数构建安全初值（防止超大/超小尺度溢出）
    long double mantissa = 0.0L;
    int exp = 0;
    extract_mantissa_and_exponent(val.data, val.scale, &mantissa, &exp);
    long double recip_m = (mantissa > 0.0L) ? (1.0L / mantissa) : 1.0L;
    PreciseDecimal x = PreciseDecimal(recip_m) * make_power_of_10(-exp);
    if (x.is_zero()) {
        x = make_power_of_10(-exp);
    }

    const PreciseDecimal two(2LL);
    const PreciseDecimal epsilon = PreciseDecimal(
        "1e-" + std::to_string(target_scale + 8));

    NormalizationSuppressor suppressor;

    // Newton-Raphson 迭代：x_{k+1} = x_k * (2 - val * x_k)
    int max_iterations = std::max(5, target_scale / 20 + 2);

    for (int i = 0; i < max_iterations; ++i) {
        PreciseDecimal next = x * (two - val * x);

        PreciseDecimal diff = precise::abs(next - x);
        if (diff.is_zero() || diff < epsilon * std::max(PreciseDecimal(1LL), precise::abs(next))) {
            g_suppress_normalization = false;
            next.normalize();
            return next;
        }

        x = next;
    }

    g_suppress_normalization = false;
    x.normalize();
    return x;
}

// ============================================================================
// 优化的除法实现（使用 Barrett 约减）
// ============================================================================

PreciseDecimal divide_optimized(const PreciseDecimal& lhs, const PreciseDecimal& rhs) {
    if (rhs.is_zero()) {
        throw std::runtime_error("division by zero");
    }
    if (lhs.is_zero()) {
        return PreciseDecimal(0LL);
    }

    // The generic implementation applies the same rounding contract as the
    // public division operator. Keep the Barrett path opt-in until it does so
    // as well; truncation here is worse than the performance difference.
    return divide_precise_decimal(lhs, rhs);

    // 对于小规模除法，使用标准算法
    if (rhs.data.size() <= 32) {
        return divide_precise_decimal(lhs, rhs);
    }

    // 使用 Barrett 约减
    BarrettReducerOptimized reducer(rhs.data, rhs.scale);

    if (!reducer.is_valid()) {
        return divide_precise_decimal(lhs, rhs);
    }

    int target_scale = PrecisionContext::get_default_scale();

    // 调整被除数的精度
    int numerator_shift = target_scale + rhs.scale - lhs.scale;
    BigIntData numerator = lhs.data;

    if (numerator_shift >= 0) {
        numerator = multiply_bigint_by_power_of_10(numerator, numerator_shift);
    } else {
        BigIntData divisor = multiply_bigint_by_power_of_10({1}, -numerator_shift);
        BigIntData truncated, ignored;
        div_bigint(numerator, divisor, &truncated, &ignored);
        numerator = truncated;
    }

    // 使用 Barrett 除法
    BigIntData q = reducer.divide(numerator, lhs.scale, target_scale);

    PreciseDecimal res;
    res.data = q;
    res.scale = target_scale;
    res.negative = lhs.negative != rhs.negative;
    res.normalize();

    return res;
}

// ============================================================================
// 缓存的 Barrett 约减器
// ============================================================================
//
// 用于需要多次使用同一除数的场景
// ============================================================================

class CachedBarrettReducer {
public:
    static CachedBarrettReducer& instance() {
        static thread_local CachedBarrettReducer inst;
        return inst;
    }

    // 设置当前除数
    void set_divisor(const BigIntData& divisor, int scale = 0) {
        if (!cached_divisor_valid_ ||
            compare_bigint(cached_divisor_, divisor) != 0 ||
            cached_scale_ != scale) {
            cached_divisor_ = divisor;
            cached_divisor_valid_ = true;
            cached_scale_ = scale;
            reducer_ = std::make_unique<BarrettReducerOptimized>(divisor, scale);
        }
    }

    // 执行除法
    BigIntData divide(const BigIntData& num, int num_scale, int target_scale) {
        if (reducer_ && reducer_->is_valid()) {
            return reducer_->divide(num, num_scale, target_scale);
        }
        throw std::runtime_error("No cached divisor set");
    }

private:
    CachedBarrettReducer() : cached_divisor_valid_(false) {}

    BigIntData cached_divisor_;
    bool cached_divisor_valid_;
    int cached_scale_ = 0;
    std::unique_ptr<BarrettReducerOptimized> reducer_;
};

// ============================================================================
// 导出函数
// ============================================================================

// 使用缓存的 Barrett 约减进行批量除法
void set_cached_divisor(const BigIntData& divisor, int scale) {
    CachedBarrettReducer::instance().set_divisor(divisor, scale);
}

BigIntData divide_with_cached_divisor(const BigIntData& num, int num_scale, int target_scale) {
    return CachedBarrettReducer::instance().divide(num, num_scale, target_scale);
}

} // namespace precise

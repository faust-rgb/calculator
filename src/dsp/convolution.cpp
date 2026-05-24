/**
 * @file convolution.cpp
 * @brief 卷积与相关函数实现
 *
 * 本文件实现了信号处理中的卷积和相关函数，包括：
 * - 线性卷积：直接法和 FFT 快速法
 * - 循环卷积：使用 FFT 实现的圆周卷积
 * - 互相关函数：测量两个信号的相似性
 * - 自相关函数：分析信号自身的周期性
 *
 * 对于小信号使用直接卷积，大信号自动切换到 FFT 快速卷积以提高效率。
 *
 * @author Calculator Development Team
 * @date 2024
 */

#include "signal_processing.h"
#include "math/mymath.h"
#include "app/scalar_type.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace signal {

// Use Scalar for internal high-precision real computations
using Scalar = mymath::Scalar;

// 使用 mymath 命名空间中的常量
using mymath::kPi;

// ============================================================================
// 线性卷积（直接法）
// ============================================================================

/**
 * @brief 计算两个信号的线性卷积
 *
 * 实现线性卷积运算，输出长度为 M + N - 1。
 * 对于小信号（N * M <= 1024）使用直接卷积，大信号自动使用 FFT 快速卷积。
 *
 * @param signal1 第一个输入信号
 * @param signal2 第二个输入信号
 * @return 卷积结果向量
 */
std::vector<Scalar> convolve(const std::vector<Scalar>& signal1,
                               const std::vector<Scalar>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();

    if (n == 0 || m == 0) {
        return {};
    }

    const std::size_t out_len = n + m - 1;
    std::vector<Scalar> result(out_len, Scalar(0.0L));

    // 对于小信号使用直接卷积（用 Scalar 进行高精度计算）
    if (n * m <= 1024) {
        std::vector<Scalar> result_scalar(out_len, Scalar(0.0L));
        std::vector<Scalar> s1_scalar(n), s2_scalar(m);
        for (std::size_t i = 0; i < n; ++i) s1_scalar[i] = Scalar(signal1[i]);
        for (std::size_t i = 0; i < m; ++i) s2_scalar[i] = Scalar(signal2[i]);

        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < m; ++j) {
                result_scalar[i + j] += s1_scalar[i] * s2_scalar[j];
            }
        }

        for (std::size_t i = 0; i < out_len; ++i) {
            result[i] = (result_scalar[i]);
        }
        return result;
    }

    // 大信号使用 FFT 卷积
    return fft_convolve(signal1, signal2);
}

// ============================================================================
// FFT 快速卷积
// ============================================================================

std::vector<Scalar> fast_convolve(const std::vector<Scalar>& signal1,
                                   const std::vector<Scalar>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();
    if (n == 0 || m == 0) return {};

    const std::size_t out_len = n + m - 1;
    const std::size_t fft_len = next_power_of_two(out_len);

    std::vector<Scalar> x1(fft_len, 0.0L);
    std::vector<Scalar> x2(fft_len, 0.0L);
    std::copy(signal1.begin(), signal1.end(), x1.begin());
    std::copy(signal2.begin(), signal2.end(), x2.begin());

    // 使用实数 FFT，输出长度为 fft_len/2 + 1
    std::vector<Complex> X1 = rfft(x1);
    std::vector<Complex> X2 = rfft(x2);

    // Use Scalar for multiplication precision
    for (std::size_t i = 0; i < X1.size(); ++i) {
        const Scalar re1(X1[i].real()), im1(X1[i].imag());
        const Scalar re2(X2[i].real()), im2(X2[i].imag());
        // Complex multiplication: (re1 + im1*i) * (re2 + im2*i)
        const Scalar new_re = re1 * re2 - im1 * im2;
        const Scalar new_im = re1 * im2 + im1 * re2;
        X1[i] = Complex((new_re), (new_im));
    }

    std::vector<Scalar> full_result = irfft(X1, fft_len);
    full_result.resize(out_len);
    return full_result;
}

std::vector<Scalar> fft_convolve(const std::vector<Scalar>& signal1,
                                  const std::vector<Scalar>& signal2) {
    return fast_convolve(signal1, signal2);
}

std::vector<Scalar> circular_convolve(const std::vector<Scalar>& signal1,
                                       const std::vector<Scalar>& signal2,
                                       std::size_t n) {
    const std::size_t len1 = signal1.size();
    const std::size_t len2 = signal2.size();
    if (len1 == 0 || len2 == 0) return {};

    // 改进：默认使用 max(len1, len2) 而非 lcm，避免内存爆炸
    if (n == 0) n = std::max(len1, len2);

    std::vector<Scalar> x1(n, 0.0L);
    std::vector<Scalar> x2(n, 0.0L);
    // 处理混叠 - use Scalar for precision
    std::vector<Scalar> x1_scalar(n, Scalar(0.0L));
    std::vector<Scalar> x2_scalar(n, Scalar(0.0L));
    for (std::size_t i = 0; i < len1; ++i) x1_scalar[i % n] += Scalar(signal1[i]);
    for (std::size_t i = 0; i < len2; ++i) x2_scalar[i % n] += Scalar(signal2[i]);
    for (std::size_t i = 0; i < n; ++i) {
        x1[i] = (x1_scalar[i]);
        x2[i] = (x2_scalar[i]);
    }

    std::vector<Complex> X1 = rfft(x1);
    std::vector<Complex> X2 = rfft(x2);

    // Use Scalar for multiplication
    for (std::size_t i = 0; i < X1.size(); ++i) {
        const Scalar re1(X1[i].real()), im1(X1[i].imag());
        const Scalar re2(X2[i].real()), im2(X2[i].imag());
        // Complex multiplication
        const Scalar new_re = re1 * re2 - im1 * im2;
        const Scalar new_im = re1 * im2 + im1 * re2;
        X1[i] = Complex((new_re), (new_im));
    }

    return irfft(X1, n);
}

std::vector<Scalar> xcorr(const std::vector<Scalar>& signal1,
                           const std::vector<Scalar>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();
    if (n == 0 || m == 0) return {};

    const std::size_t out_len = n + m - 1;
    const std::size_t fft_len = next_power_of_two(out_len);

    std::vector<Scalar> x1(fft_len, 0.0L);
    std::vector<Scalar> x2(fft_len, 0.0L);
    std::copy(signal1.begin(), signal1.end(), x1.begin());
    std::copy(signal2.begin(), signal2.end(), x2.begin());

    std::vector<Complex> X1 = rfft(x1);
    std::vector<Complex> X2 = rfft(x2);

    // 互相关：X1 * conj(X2) - use Scalar for precision
    for (std::size_t i = 0; i < X1.size(); ++i) {
        const Scalar re1(X1[i].real()), im1(X1[i].imag());
        const Scalar re2(X2[i].real()), im2(X2[i].imag());
        // conj(X2) = re2 - im2*i
        // (re1 + im1*i) * (re2 - im2*i) = (re1*re2 + im1*im2) + (im1*re2 - re1*im2)*i
        const Scalar new_re = re1 * re2 + im1 * im2;
        const Scalar new_im = im1 * re2 - re1 * im2;
        X1[i] = Complex((new_re), (new_im));
    }

    std::vector<Scalar> result = irfft(X1, fft_len);
    result.resize(out_len);
    return result;
}

// ============================================================================
// 自相关函数
// ============================================================================

std::vector<Scalar> autocorr(const std::vector<Scalar>& signal) {
    return xcorr(signal, signal);
}

// ============================================================================
// 归一化互相关
// ============================================================================

std::vector<Scalar> normalized_xcorr(const std::vector<Scalar>& signal1,
                                      const std::vector<Scalar>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();

    if (n == 0 || m == 0) {
        return {};
    }

    // 计算能量 - use Scalar for precision
    Scalar energy1 = Scalar(0.0L);
    Scalar energy2 = Scalar(0.0L);
    for (std::size_t i = 0; i < n; ++i) {
        const Scalar val(signal1[i]);
        energy1 += val * val;
    }
    for (std::size_t i = 0; i < m; ++i) {
        const Scalar val(signal2[i]);
        energy2 += val * val;
    }

    if (energy1 == Scalar(0.0L) || energy2 == Scalar(0.0L)) {
        return std::vector<Scalar>(n + m - 1, 0.0L);
    }

    // 计算互相关
    std::vector<Scalar> corr = xcorr(signal1, signal2);

    // 归一化
    const Scalar norm_factor = mymath::sqrt(energy1 * energy2);
    for (Scalar& value : corr) {
        value = (Scalar(value) / norm_factor);
    }

    return corr;
}

}  // namespace signal

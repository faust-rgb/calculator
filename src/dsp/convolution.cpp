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

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace signal {

// 数学常量
constexpr long double kPi = 3.14159265358979323846;

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
std::vector<long double> convolve(const std::vector<long double>& signal1,
                               const std::vector<long double>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();

    if (n == 0 || m == 0) {
        return {};
    }

    const std::size_t out_len = n + m - 1;
    std::vector<long double> result(out_len, 0.0L);

    // 对于小信号使用直接卷积
    if (n * m <= 1024) {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < m; ++j) {
                result[i + j] += signal1[i] * signal2[j];
            }
        }
        return result;
    }

    // 大信号使用 FFT 卷积
    return fft_convolve(signal1, signal2);
}

// ============================================================================
// FFT 快速卷积
// ============================================================================

std::vector<long double> fast_convolve(const std::vector<long double>& signal1,
                                   const std::vector<long double>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();
    if (n == 0 || m == 0) return {};

    const std::size_t out_len = n + m - 1;
    const std::size_t fft_len = next_power_of_two(out_len);

    std::vector<long double> x1(fft_len, 0.0L);
    std::vector<long double> x2(fft_len, 0.0L);
    std::copy(signal1.begin(), signal1.end(), x1.begin());
    std::copy(signal2.begin(), signal2.end(), x2.begin());

    // 使用实数 FFT，输出长度为 fft_len/2 + 1
    std::vector<Complex> X1 = rfft(x1);
    std::vector<Complex> X2 = rfft(x2);

    for (std::size_t i = 0; i < X1.size(); ++i) {
        X1[i] *= X2[i];
    }

    std::vector<long double> full_result = irfft(X1, fft_len);
    full_result.resize(out_len);
    return full_result;
}

std::vector<long double> fft_convolve(const std::vector<long double>& signal1,
                                  const std::vector<long double>& signal2) {
    return fast_convolve(signal1, signal2);
}

std::vector<long double> circular_convolve(const std::vector<long double>& signal1,
                                       const std::vector<long double>& signal2,
                                       std::size_t n) {
    const std::size_t len1 = signal1.size();
    const std::size_t len2 = signal2.size();
    if (len1 == 0 || len2 == 0) return {};

    // 改进：默认使用 max(len1, len2) 而非 lcm，避免内存爆炸
    if (n == 0) n = std::max(len1, len2);

    std::vector<long double> x1(n, 0.0L);
    std::vector<long double> x2(n, 0.0L);
    // 处理混叠
    for (std::size_t i = 0; i < len1; ++i) x1[i % n] += signal1[i];
    for (std::size_t i = 0; i < len2; ++i) x2[i % n] += signal2[i];

    std::vector<Complex> X1 = rfft(x1);
    std::vector<Complex> X2 = rfft(x2);
    for (std::size_t i = 0; i < X1.size(); ++i) X1[i] *= X2[i];

    return irfft(X1, n);
}

std::vector<long double> xcorr(const std::vector<long double>& signal1,
                           const std::vector<long double>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();
    if (n == 0 || m == 0) return {};

    const std::size_t out_len = n + m - 1;
    const std::size_t fft_len = next_power_of_two(out_len);

    std::vector<long double> x1(fft_len, 0.0L);
    std::vector<long double> x2(fft_len, 0.0L);
    std::copy(signal1.begin(), signal1.end(), x1.begin());
    std::copy(signal2.begin(), signal2.end(), x2.begin());

    std::vector<Complex> X1 = rfft(x1);
    std::vector<Complex> X2 = rfft(x2);

    // 互相关：X1 * conj(X2)
    for (std::size_t i = 0; i < X1.size(); ++i) {
        X1[i] *= mymath::conj(X2[i]);
    }

    std::vector<long double> result = irfft(X1, fft_len);
    result.resize(out_len);
    return result;
}

// ============================================================================
// 自相关函数
// ============================================================================

std::vector<long double> autocorr(const std::vector<long double>& signal) {
    return xcorr(signal, signal);
}

// ============================================================================
// 归一化互相关
// ============================================================================

std::vector<long double> normalized_xcorr(const std::vector<long double>& signal1,
                                      const std::vector<long double>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();

    if (n == 0 || m == 0) {
        return {};
    }

    // 计算能量
    long double energy1 = 0.0L;
    long double energy2 = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
        energy1 += signal1[i] * signal1[i];
    }
    for (std::size_t i = 0; i < m; ++i) {
        energy2 += signal2[i] * signal2[i];
    }

    if (energy1 == 0.0L || energy2 == 0.0L) {
        return std::vector<long double>(n + m - 1, 0.0L);
    }

    // 计算互相关
    std::vector<long double> corr = xcorr(signal1, signal2);

    // 归一化
    const long double norm_factor = mymath::sqrt(energy1 * energy2);
    for (long double& value : corr) {
        value /= norm_factor;
    }

    return corr;
}

}  // namespace signal

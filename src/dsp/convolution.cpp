/**
 * @file convolution.cpp
 * @brief 卷积与相关函数实现
 *
 * 实现：
 * - 线性卷积
 * - FFT 快速卷积
 * - 循环卷积
 * - 互相关与自相关
 */

#include "signal_processing.h"
#include "math/mymath.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace signal {

// 数学常量
constexpr double kPi = 3.14159265358979323846;

// ============================================================================
// 线性卷积（直接法）
// ============================================================================

std::vector<double> convolve(const std::vector<double>& signal1,
                               const std::vector<double>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();

    if (n == 0 || m == 0) {
        return {};
    }

    const std::size_t out_len = n + m - 1;
    std::vector<double> result(out_len, 0.0);

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

std::vector<double> fast_convolve(const std::vector<double>& signal1,
                                   const std::vector<double>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();
    if (n == 0 || m == 0) return {};

    const std::size_t out_len = n + m - 1;
    const std::size_t fft_len = next_power_of_two(out_len);

    std::vector<double> x1(fft_len, 0.0);
    std::vector<double> x2(fft_len, 0.0);
    std::copy(signal1.begin(), signal1.end(), x1.begin());
    std::copy(signal2.begin(), signal2.end(), x2.begin());

    // 使用实数 FFT，输出长度为 fft_len/2 + 1
    std::vector<Complex> X1 = rfft(x1);
    std::vector<Complex> X2 = rfft(x2);

    for (std::size_t i = 0; i < X1.size(); ++i) {
        X1[i] *= X2[i];
    }

    std::vector<double> full_result = irfft(X1, fft_len);
    full_result.resize(out_len);
    return full_result;
}

std::vector<double> fft_convolve(const std::vector<double>& signal1,
                                  const std::vector<double>& signal2) {
    return fast_convolve(signal1, signal2);
}

std::vector<double> circular_convolve(const std::vector<double>& signal1,
                                       const std::vector<double>& signal2,
                                       std::size_t n) {
    const std::size_t len1 = signal1.size();
    const std::size_t len2 = signal2.size();
    if (len1 == 0 || len2 == 0) return {};

    // 改进：默认使用 max(len1, len2) 而非 lcm，避免内存爆炸
    if (n == 0) n = std::max(len1, len2);

    std::vector<double> x1(n, 0.0);
    std::vector<double> x2(n, 0.0);
    for (std::size_t i = 0; i < std::min(len1, n); ++i) x1[i] = signal1[i];
    for (std::size_t i = 0; i < std::min(len2, n); ++i) x2[i] = signal2[i];

    std::vector<Complex> X1 = rfft(x1);
    std::vector<Complex> X2 = rfft(x2);
    for (std::size_t i = 0; i < X1.size(); ++i) X1[i] *= X2[i];

    return irfft(X1, n);
}

std::vector<double> xcorr(const std::vector<double>& signal1,
                           const std::vector<double>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();
    if (n == 0 || m == 0) return {};

    const std::size_t out_len = n + m - 1;
    const std::size_t fft_len = next_power_of_two(out_len);

    std::vector<double> x1(fft_len, 0.0);
    std::vector<double> x2(fft_len, 0.0);
    std::copy(signal1.begin(), signal1.end(), x1.begin());
    std::copy(signal2.begin(), signal2.end(), x2.begin());

    std::vector<Complex> X1 = rfft(x1);
    std::vector<Complex> X2 = rfft(x2);

    // 互相关：X1 * conj(X2)
    for (std::size_t i = 0; i < X1.size(); ++i) {
        X1[i] *= mymath::conj(X2[i]);
    }

    std::vector<double> result = irfft(X1, fft_len);
    result.resize(out_len);
    return result;
}

// ============================================================================
// 自相关函数
// ============================================================================

std::vector<double> autocorr(const std::vector<double>& signal) {
    return xcorr(signal, signal);
}

// ============================================================================
// 归一化互相关
// ============================================================================

std::vector<double> normalized_xcorr(const std::vector<double>& signal1,
                                      const std::vector<double>& signal2) {
    const std::size_t n = signal1.size();
    const std::size_t m = signal2.size();

    if (n == 0 || m == 0) {
        return {};
    }

    // 计算能量
    double energy1 = 0.0;
    double energy2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        energy1 += signal1[i] * signal1[i];
    }
    for (std::size_t i = 0; i < m; ++i) {
        energy2 += signal2[i] * signal2[i];
    }

    if (energy1 == 0.0 || energy2 == 0.0) {
        return std::vector<double>(n + m - 1, 0.0);
    }

    // 计算互相关
    std::vector<double> corr = xcorr(signal1, signal2);

    // 归一化
    const double norm_factor = mymath::sqrt(energy1 * energy2);
    for (double& value : corr) {
        value /= norm_factor;
    }

    return corr;
}

}  // namespace signal

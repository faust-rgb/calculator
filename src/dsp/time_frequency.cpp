/**
 * @file time_frequency.cpp
 * @brief 时频分析实现
 *
 * 实现：
 * - 功率谱密度（周期图法、Welch 方法）
 * - 短时傅里叶变换（STFT）
 * - 语谱图
 */

#include "signal_processing.h"
#include "math/mymath.h"
#include "core/common/scalar_type.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace signal {

// Use Scalar for internal high-precision real computations
using Scalar = mymath::Scalar;

// 使用 mymath 命名空间中的常量
using mymath::kPi;

// ============================================================================
// 功率谱密度（周期图法）
// ============================================================================

std::vector<Scalar> periodogram(const std::vector<Scalar>& signal, Scalar sample_rate) {
    if (signal.empty()) {
        return {};
    }

    const std::size_t n = signal.size();

    // 计算 FFT
    std::vector<Complex> spectrum = rfft(signal);

    // 计算功率谱 - use Scalar for precision
    std::vector<Scalar> psd(spectrum.size());
    const Scalar scale = Scalar(2.0L) / (Scalar((n)) * Scalar(sample_rate));

    for (std::size_t i = 0; i < spectrum.size(); ++i) {
        const Scalar re(spectrum[i].real());
        const Scalar im(spectrum[i].imag());
        const Scalar norm_val = re * re + im * im;
        psd[i] = (scale * norm_val / Scalar((n)));
    }

    // 处理 DC 和 Nyquist (不需要乘2)
    if (!psd.empty()) psd[0] /= 2.0;
    if (n % 2 == 0 && psd.size() > 1) psd.back() /= 2.0;

    return psd;
}

// ============================================================================
// 功率谱密度（Welch 方法）
// ============================================================================

std::vector<Scalar> pwelch(const std::vector<Scalar>& signal,
                            std::size_t nfft,
                            WindowType window_type,
                            std::size_t noverlap,
                            Scalar sample_rate) {
    if (signal.empty()) {
        return {};
    }

    const std::size_t n = signal.size();

    // 默认参数
    if (nfft == 0) {
        nfft = std::min(n, static_cast<std::size_t>(256));
    }
    if (nfft > n) {
        nfft = n;
    }

    // 默认重叠 50%
    if (noverlap == 0) {
        noverlap = nfft / 2;
    }
    if (noverlap >= nfft) {
        throw std::runtime_error("Overlap must be less than segment length");
    }

    // 生成窗函数
    std::vector<Scalar> win = window(window_type, nfft);

    // Use Scalar for window power calculation
    Scalar win_power = Scalar(0.0L);
    for (std::size_t i = 0; i < nfft; ++i) {
        const Scalar w(win[i]);
        win_power += w * w;
    }

    // 计算分段数
    const std::size_t step = nfft - noverlap;
    const std::size_t n_segments = (n - nfft) / step + 1;

    if (n_segments == 0) {
        return periodogram(signal, sample_rate);
    }

    // 初始化 PSD 累加器 - use Scalar for precision
    std::vector<Scalar> psd_sum(nfft / 2 + 1, Scalar(0.0L));

    // 对每个分段计算周期图
    for (std::size_t seg = 0; seg < n_segments; ++seg) {
        const std::size_t start = seg * step;

        // 提取并加窗
        std::vector<Scalar> segment(nfft);
        for (std::size_t i = 0; i < nfft; ++i) {
            segment[i] = signal[start + i] * win[i];
        }

        // FFT
        std::vector<Complex> spectrum = rfft(segment);

        // 累加功率
        for (std::size_t k = 0; k < spectrum.size(); ++k) {
            const Scalar re(spectrum[k].real());
            const Scalar im(spectrum[k].imag());
            psd_sum[k] += re * re + im * im;
        }
    }

    // 平均并归一化
    std::vector<Scalar> psd(psd_sum.size());
    const Scalar scale = Scalar(2.0L) / (Scalar((n_segments)) *
                                          Scalar(sample_rate) * win_power);

    for (std::size_t k = 0; k < psd.size(); ++k) {
        psd[k] = (psd_sum[k] * scale);
    }

    // 处理 DC 和 Nyquist
    if (!psd.empty()) psd[0] /= 2.0;
    if (nfft % 2 == 0 && psd.size() > 1) psd.back() /= 2.0;

    return psd;
}

// ============================================================================
// 短时傅里叶变换（STFT）
// ============================================================================

STFTResult stft(const std::vector<Scalar>& signal,
                std::size_t nfft,
                WindowType window_type,
                std::size_t noverlap) {
    STFTResult result;

    if (signal.empty()) {
        return result;
    }

    const std::size_t n = signal.size();

    // 默认参数
    if (nfft == 0) {
        nfft = std::min(n, static_cast<std::size_t>(256));
    }

    if (noverlap == 0) {
        noverlap = nfft / 2;
    }
    if (noverlap >= nfft) {
        throw std::runtime_error("Overlap must be less than FFT length");
    }

    // 生成窗函数
    std::vector<Scalar> win = window(window_type, nfft);

    // 计算分段数
    const std::size_t step = nfft - noverlap;
    const std::size_t n_frames = (n - noverlap) / step;

    if (n_frames == 0) {
        return result;
    }

    // 频率轴
    const std::size_t n_freqs = nfft / 2 + 1;
    result.freq_vector.resize(n_freqs);
    for (std::size_t k = 0; k < n_freqs; ++k) {
        result.freq_vector[k] = (k) / (nfft);
    }

    // 时间轴
    result.time_vector.resize(n_frames);
    for (std::size_t frame = 0; frame < n_frames; ++frame) {
        result.time_vector[frame] = (frame * step + nfft / 2);
    }

    // STFT 矩阵
    result.stft_matrix.resize(n_frames);
    result.window_type = window_type;

    for (std::size_t frame = 0; frame < n_frames; ++frame) {
        const std::size_t start = frame * step;

        // 提取并加窗
        std::vector<Scalar> frame_data(nfft);
        for (std::size_t i = 0; i < nfft; ++i) {
            frame_data[i] = signal[start + i] * win[i];
        }

        // FFT
        result.stft_matrix[frame] = rfft(frame_data);
    }

    return result;
}

// ============================================================================
// 逆短时傅里叶变换
// ============================================================================

std::vector<Scalar> istft(const STFTResult& stft_result) {
    if (stft_result.stft_matrix.empty()) {
        return {};
    }

    const std::size_t n_frames = stft_result.stft_matrix.size();
    const std::size_t n_freqs = stft_result.stft_matrix[0].size();
    const std::size_t nfft = (n_freqs - 1) * 2;

    // 假设 50% 重叠
    const std::size_t step = nfft / 2;
    const std::size_t n_samples = (n_frames - 1) * step + nfft;

    // Use Scalar for precision in accumulation
    std::vector<Scalar> signal_scalar(n_samples, Scalar(0.0L));
    std::vector<Scalar> window_sum(n_samples, Scalar(0.0L));

    // 使用 STFT 使用的窗函数进行重叠相加
    std::vector<Scalar> win = window(stft_result.window_type, nfft);

    for (std::size_t frame = 0; frame < n_frames; ++frame) {
        // IFFT
        std::vector<Scalar> frame_time = irfft(stft_result.stft_matrix[frame], nfft);

        // 重叠相加
        const std::size_t start = frame * step;
        for (std::size_t i = 0; i < nfft; ++i) {
            const Scalar frame_val(frame_time[i]);
            const Scalar win_val(win[i]);
            signal_scalar[start + i] += frame_val * win_val;
            window_sum[start + i] += win_val * win_val;
        }
    }

    // 归一化
    std::vector<Scalar> signal(n_samples);
    for (std::size_t i = 0; i < n_samples; ++i) {
        if (window_sum[i] > Scalar(1e-10L)) {
            signal[i] = (signal_scalar[i] / window_sum[i]);
        } else {
            signal[i] = 0.0L;
        }
    }

    return signal;
}

// ============================================================================
// 语谱图
// ============================================================================

std::vector<std::vector<Scalar>> spectrogram(const std::vector<Scalar>& signal,
                                               std::size_t nfft,
                                               WindowType window_type,
                                               std::size_t noverlap) {
    STFTResult stft_data = stft(signal, nfft, window_type, noverlap);

    if (stft_data.stft_matrix.empty()) {
        return {};
    }

    const std::size_t n_frames = stft_data.stft_matrix.size();
    const std::size_t n_freqs = stft_data.stft_matrix[0].size();

    std::vector<std::vector<Scalar>> spec(n_freqs, std::vector<Scalar>(n_frames));

    // 转置并计算功率 - use Scalar for precision
    for (std::size_t frame = 0; frame < n_frames; ++frame) {
        for (std::size_t freq = 0; freq < n_freqs; ++freq) {
            const Scalar re(stft_data.stft_matrix[frame][freq].real());
            const Scalar im(stft_data.stft_matrix[frame][freq].imag());
            spec[freq][frame] = (re * re + im * im);
        }
    }

    return spec;
}

}  // namespace signal

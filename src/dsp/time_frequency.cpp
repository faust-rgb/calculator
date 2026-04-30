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
#include "../math/mymath.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace signal {

// 数学常量
constexpr double kPi = 3.14159265358979323846;

// ============================================================================
// 功率谱密度（周期图法）
// ============================================================================

std::vector<double> periodogram(const std::vector<double>& signal) {
    if (signal.empty()) {
        return {};
    }

    const std::size_t n = signal.size();

    // 计算 FFT
    std::vector<Complex> spectrum = rfft(signal);

    // 计算功率谱
    std::vector<double> psd(spectrum.size());
    const double scale = 1.0 / static_cast<double>(n * n);

    for (std::size_t i = 0; i < spectrum.size(); ++i) {
        psd[i] = scale * mymath::norm(spectrum[i]);
    }

    return psd;
}

// ============================================================================
// 功率谱密度（Welch 方法）
// ============================================================================

std::vector<double> pwelch(const std::vector<double>& signal,
                            std::size_t nfft,
                            WindowType window_type,
                            std::size_t noverlap) {
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
    std::vector<double> win = window(window_type, nfft);
    const double win_energy = mymath::sqrt(std::inner_product(win.begin(), win.end(), win.begin(), 0.0));

    // 计算分段数
    const std::size_t step = nfft - noverlap;
    const std::size_t n_segments = (n - nfft) / step + 1;

    if (n_segments == 0) {
        return periodogram(signal);
    }

    // 初始化 PSD 累加器
    std::vector<double> psd_sum(nfft / 2 + 1, 0.0);

    // 对每个分段计算周期图
    for (std::size_t seg = 0; seg < n_segments; ++seg) {
        const std::size_t start = seg * step;

        // 提取并加窗
        std::vector<double> segment(nfft);
        for (std::size_t i = 0; i < nfft; ++i) {
            segment[i] = signal[start + i] * win[i];
        }

        // FFT
        std::vector<Complex> spectrum = rfft(segment);

        // 累加功率
        for (std::size_t k = 0; k < spectrum.size(); ++k) {
            psd_sum[k] += mymath::norm(spectrum[k]);
        }
    }

    // 平均并归一化
    std::vector<double> psd(psd_sum.size());
    const double scale = 1.0 / (static_cast<double>(n_segments) * win_energy * win_energy);

    for (std::size_t k = 0; k < psd.size(); ++k) {
        psd[k] = psd_sum[k] * scale;
    }

    return psd;
}

// ============================================================================
// 短时傅里叶变换（STFT）
// ============================================================================

STFTResult stft(const std::vector<double>& signal,
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
    std::vector<double> win = window(window_type, nfft);

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
        result.freq_vector[k] = static_cast<double>(k) / static_cast<double>(nfft);
    }

    // 时间轴
    result.time_vector.resize(n_frames);
    for (std::size_t frame = 0; frame < n_frames; ++frame) {
        result.time_vector[frame] = static_cast<double>(frame * step + nfft / 2);
    }

    // STFT 矩阵
    result.stft_matrix.resize(n_frames);

    for (std::size_t frame = 0; frame < n_frames; ++frame) {
        const std::size_t start = frame * step;

        // 提取并加窗
        std::vector<double> frame_data(nfft);
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

std::vector<double> istft(const STFTResult& stft_result) {
    if (stft_result.stft_matrix.empty()) {
        return {};
    }

    const std::size_t n_frames = stft_result.stft_matrix.size();
    const std::size_t n_freqs = stft_result.stft_matrix[0].size();
    const std::size_t nfft = (n_freqs - 1) * 2;

    // 假设 50% 重叠
    const std::size_t step = nfft / 2;
    const std::size_t n_samples = (n_frames - 1) * step + nfft;

    std::vector<double> signal(n_samples, 0.0);
    std::vector<double> window_sum(n_samples, 0.0);

    // 使用汉宁窗进行重叠相加
    std::vector<double> win = hanning_window(nfft);

    for (std::size_t frame = 0; frame < n_frames; ++frame) {
        // IFFT
        std::vector<double> frame_time = irfft(stft_result.stft_matrix[frame], nfft);

        // 重叠相加
        const std::size_t start = frame * step;
        for (std::size_t i = 0; i < nfft; ++i) {
            signal[start + i] += frame_time[i] * win[i];
            window_sum[start + i] += win[i] * win[i];
        }
    }

    // 归一化
    for (std::size_t i = 0; i < n_samples; ++i) {
        if (window_sum[i] > 1e-10) {
            signal[i] /= window_sum[i];
        }
    }

    return signal;
}

// ============================================================================
// 语谱图
// ============================================================================

std::vector<std::vector<double>> spectrogram(const std::vector<double>& signal,
                                               std::size_t nfft,
                                               WindowType window_type,
                                               std::size_t noverlap) {
    STFTResult stft_data = stft(signal, nfft, window_type, noverlap);

    if (stft_data.stft_matrix.empty()) {
        return {};
    }

    const std::size_t n_frames = stft_data.stft_matrix.size();
    const std::size_t n_freqs = stft_data.stft_matrix[0].size();

    std::vector<std::vector<double>> spec(n_freqs, std::vector<double>(n_frames));

    // 转置并计算功率
    for (std::size_t frame = 0; frame < n_frames; ++frame) {
        for (std::size_t freq = 0; freq < n_freqs; ++freq) {
            spec[freq][frame] = mymath::norm(stft_data.stft_matrix[frame][freq]);
        }
    }

    return spec;
}

}  // namespace signal

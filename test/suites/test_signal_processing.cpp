/**
 * @file test_signal_processing.cpp
 * @brief 信号处理模块测试实现
 *
 * 测试信号处理模块的各项功能：
 * - FFT 算法（基-2、混合基、Bluestein、实数 FFT）
 * - 卷积与相关分析
 * - 窗函数
 * - 滤波器设计
 * - 时频分析
 */

#include "dsp/signal_processing.h"
#include "math/mymath.h"
#include "test_helpers.h"
#include <iostream>
#include <vector>

namespace test_suites {

// 容差
constexpr double kEps = 1e-10;

// ============================================================================
// FFT 测试
// ============================================================================

static void test_fft_radix2(int& passed, int& failed) {
    std::cout << "  Testing FFT radix-2..." << std::endl;

    // 测试简单信号
    std::vector<signal::Complex> sig = {1.0, 0.0, 0.0, 0.0};
    std::vector<signal::Complex> spectrum = signal::fft_radix2(sig);

    bool ok = (spectrum.size() == 4);
    if (ok) {
        for (const auto& c : spectrum) {
            if (!test_helpers::nearly_equal(c.real(), 1.0, kEps) ||
                !test_helpers::nearly_equal(c.imag(), 0.0, kEps)) {
                ok = false;
                break;
            }
        }
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: FFT radix-2 basic test" << std::endl;
    }
}

static void test_fft_ifft_roundtrip(int& passed, int& failed) {
    std::cout << "  Testing FFT/IFFT roundtrip..." << std::endl;

    std::vector<signal::Complex> original = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

    std::vector<signal::Complex> spectrum = signal::fft(original);
    std::vector<signal::Complex> reconstructed = signal::ifft(spectrum);

    bool ok = (reconstructed.size() == original.size());
    if (ok) {
        for (std::size_t i = 0; i < original.size(); ++i) {
            if (!test_helpers::nearly_equal(reconstructed[i].real(), original[i].real(), kEps) ||
                !test_helpers::nearly_equal(reconstructed[i].imag(), 0.0, kEps)) {
                ok = false;
                break;
            }
        }
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: FFT/IFFT roundtrip test" << std::endl;
    }
}

static void test_fft_non_power_of_two(int& passed, int& failed) {
    std::cout << "  Testing FFT non-power-of-two..." << std::endl;

    std::vector<signal::Complex> sig = {1.0, 2.0, 3.0, 4.0, 5.0};

    std::vector<signal::Complex> spectrum = signal::fft(sig);
    std::vector<signal::Complex> reconstructed = signal::ifft(spectrum);

    bool ok = (reconstructed.size() == sig.size());
    if (ok) {
        for (std::size_t i = 0; i < sig.size(); ++i) {
            if (!test_helpers::nearly_equal(reconstructed[i].real(), sig[i].real(), kEps)) {
                ok = false;
                break;
            }
        }
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: FFT non-power-of-two test" << std::endl;
    }
}

static void test_rfft(int& passed, int& failed) {
    std::cout << "  Testing RFFT..." << std::endl;

    std::vector<double> sig = {1.0, 1.0, 0.0, 0.0};

    std::vector<signal::Complex> spectrum = signal::rfft(sig);

    bool ok = (spectrum.size() == 3u);
    if (ok) {
        ok = test_helpers::nearly_equal(spectrum[0].real(), 2.0, kEps);
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: RFFT test" << std::endl;
    }
}

static void test_fftshift(int& passed, int& failed) {
    std::cout << "  Testing FFT shift..." << std::endl;

    std::vector<signal::Complex> spectrum = {1.0, 2.0, 3.0, 4.0};

    std::vector<signal::Complex> shifted = signal::fftshift(spectrum);
    std::vector<signal::Complex> unshifted = signal::ifftshift(shifted);

    bool ok = (shifted.size() == 4u);
    if (ok) {
        ok = test_helpers::nearly_equal(shifted[0].real(), 3.0, kEps) &&
             test_helpers::nearly_equal(shifted[1].real(), 4.0, kEps) &&
             test_helpers::nearly_equal(shifted[2].real(), 1.0, kEps) &&
             test_helpers::nearly_equal(shifted[3].real(), 2.0, kEps);
    }

    if (ok) {
        for (std::size_t i = 0; i < spectrum.size(); ++i) {
            if (!test_helpers::nearly_equal(unshifted[i].real(), spectrum[i].real(), kEps)) {
                ok = false;
                break;
            }
        }
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: FFT shift test" << std::endl;
    }
}

// ============================================================================
// 卷积测试
// ============================================================================

static void test_convolve(int& passed, int& failed) {
    std::cout << "  Testing convolution..." << std::endl;

    std::vector<double> signal1 = {1.0, 2.0, 3.0};
    std::vector<double> signal2 = {1.0, 1.0};

    std::vector<double> result = signal::convolve(signal1, signal2);

    bool ok = (result.size() == 4u);
    if (ok) {
        ok = test_helpers::nearly_equal(result[0], 1.0, kEps) &&
             test_helpers::nearly_equal(result[1], 3.0, kEps) &&
             test_helpers::nearly_equal(result[2], 5.0, kEps) &&
             test_helpers::nearly_equal(result[3], 3.0, kEps);
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Convolution test" << std::endl;
    }
}

static void test_circular_convolve(int& passed, int& failed) {
    std::cout << "  Testing circular convolution..." << std::endl;

    std::vector<double> signal1 = {1.0, 2.0, 3.0};
    std::vector<double> signal2 = {1.0, 1.0, 0.0};

    std::vector<double> result = signal::circular_convolve(signal1, signal2, 3);

    bool ok = (result.size() == 3u);

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Circular convolution test" << std::endl;
    }
}

static void test_xcorr(int& passed, int& failed) {
    std::cout << "  Testing cross-correlation..." << std::endl;

    std::vector<double> signal1 = {1.0, 2.0, 3.0};
    std::vector<double> signal2 = {1.0, 2.0};

    std::vector<double> result = signal::xcorr(signal1, signal2);

    bool ok = (result.size() == 4u);

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Cross-correlation test" << std::endl;
    }
}

static void test_autocorr(int& passed, int& failed) {
    std::cout << "  Testing auto-correlation..." << std::endl;

    std::vector<double> sig = {1.0, 2.0, 3.0};

    std::vector<double> result = signal::autocorr(sig);

    bool ok = (result.size() == 5u);
    if (ok) {
        // 自相关在零延迟处应该最大
        // 对于信号 {1, 2, 3}，能量是 1+4+9=14
        // xcorr 返回的结果中，零延迟在索引 0
        ok = test_helpers::nearly_equal(result[0], 14.0, 1e-6);
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Auto-correlation test (got " << result[0] << ", expected 14)" << std::endl;
    }
}

// ============================================================================
// 窗函数测试
// ============================================================================

static void test_hanning_window(int& passed, int& failed) {
    std::cout << "  Testing Hanning window..." << std::endl;

    std::vector<double> win = signal::hanning_window(8);

    bool ok = (win.size() == 8u);
    if (ok) {
        // 检查对称性
        for (std::size_t i = 0; i < 4; ++i) {
            if (!test_helpers::nearly_equal(win[i], win[7 - i], kEps)) {
                ok = false;
                break;
            }
        }
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Hanning window test" << std::endl;
    }
}

static void test_hamming_window(int& passed, int& failed) {
    std::cout << "  Testing Hamming window..." << std::endl;

    std::vector<double> win = signal::hamming_window(8);

    bool ok = (win.size() == 8u);
    if (ok) {
        // 检查对称性
        for (std::size_t i = 0; i < 4; ++i) {
            if (!test_helpers::nearly_equal(win[i], win[7 - i], kEps)) {
                ok = false;
                break;
            }
        }
    }
    if (ok) {
        // 端点不应该为零
        ok = (win[0] > 0.0 && win[7] > 0.0);
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Hamming window test" << std::endl;
    }
}

static void test_kaiser_window(int& passed, int& failed) {
    std::cout << "  Testing Kaiser window..." << std::endl;

    std::vector<double> win = signal::kaiser_window(16, 5.0);

    bool ok = (win.size() == 16u);
    if (ok) {
        // 检查对称性
        for (std::size_t i = 0; i < 8; ++i) {
            if (!test_helpers::nearly_equal(win[i], win[15 - i], kEps)) {
                ok = false;
                break;
            }
        }
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Kaiser window test" << std::endl;
    }
}

static void test_window_function(int& passed, int& failed) {
    std::cout << "  Testing generic window function..." << std::endl;

    std::vector<double> hann = signal::window(signal::WindowType::Hanning, 8);
    std::vector<double> rect = signal::window(signal::WindowType::Rectangular, 8);

    bool ok = (hann.size() == 8u && rect.size() == 8u);
    if (ok) {
        // 矩形窗应该全为 1
        for (double v : rect) {
            if (!test_helpers::nearly_equal(v, 1.0, kEps)) {
                ok = false;
                break;
            }
        }
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Generic window function test" << std::endl;
    }
}

// ============================================================================
// 滤波器测试
// ============================================================================

static void test_fir_lowpass(int& passed, int& failed) {
    std::cout << "  Testing FIR lowpass filter design..." << std::endl;

    signal::FilterCoefficients fir = signal::design_fir(16, 0.2, signal::FilterType::LowPass);

    bool ok = !fir.b.empty() && fir.a.size() == 1u;
    if (ok) {
        ok = test_helpers::nearly_equal(fir.a[0], 1.0, kEps);
    }
    if (ok) {
        // 检查对称性
        std::size_t n = fir.b.size();
        for (std::size_t i = 0; i < n / 2; ++i) {
            if (!test_helpers::nearly_equal(fir.b[i], fir.b[n - 1 - i], kEps)) {
                ok = false;
                break;
            }
        }
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: FIR lowpass test" << std::endl;
    }
}

static void test_filter_application(int& passed, int& failed) {
    std::cout << "  Testing filter application..." << std::endl;

    std::vector<double> b = {0.5, 0.5};
    std::vector<double> a = {1.0};
    std::vector<double> sig = {1.0, 2.0, 3.0, 4.0, 5.0};

    std::vector<double> filtered = signal::filter(b, a, sig);

    bool ok = (filtered.size() == sig.size());
    if (ok) {
        ok = test_helpers::nearly_equal(filtered[0], 0.5, kEps) &&
             test_helpers::nearly_equal(filtered[1], 1.5, kEps) &&
             test_helpers::nearly_equal(filtered[2], 2.5, kEps);
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Filter application test" << std::endl;
    }
}

static void test_freqz(int& passed, int& failed) {
    std::cout << "  Testing frequency response..." << std::endl;

    std::vector<double> b = {1.0, 0.0};
    std::vector<double> a = {1.0, -0.5};

    std::vector<signal::Complex> response = signal::freqz(b, a, 64);

    bool ok = (response.size() == 64u);
    if (ok) {
        // 直流响应应该是 1/(1-0.5) = 2
        ok = mymath::abs(response[0]) > 1.5 && mymath::abs(response[0]) < 3.0;
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Frequency response test" << std::endl;
    }
}

// ============================================================================
// 时频分析测试
// ============================================================================

static void test_periodogram(int& passed, int& failed) {
    std::cout << "  Testing periodogram..." << std::endl;

    std::vector<double> sig(64);
    for (std::size_t i = 0; i < sig.size(); ++i) {
        sig[i] = mymath::sin(2.0 * mymath::kPi * 8.0 * static_cast<double>(i) / 64.0);
    }

    std::vector<double> psd = signal::periodogram(sig);

    bool ok = !psd.empty();
    if (ok) {
        // 找峰值位置
        std::size_t peak_idx = 0;
        double max_val = 0.0;
        for (std::size_t i = 0; i < psd.size(); ++i) {
            if (psd[i] > max_val) {
                max_val = psd[i];
                peak_idx = i;
            }
        }
        // 峰值应该在频率索引 8 附近
        ok = (peak_idx >= 6 && peak_idx <= 10);
    }

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Periodogram test" << std::endl;
    }
}

static void test_pwelch(int& passed, int& failed) {
    std::cout << "  Testing Welch PSD..." << std::endl;

    std::vector<double> sig(256);
    for (std::size_t i = 0; i < sig.size(); ++i) {
        sig[i] = mymath::sin(2.0 * mymath::kPi * 32.0 * static_cast<double>(i) / 256.0);
    }

    std::vector<double> psd = signal::pwelch(sig, 64);

    bool ok = !psd.empty();

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Welch PSD test" << std::endl;
    }
}

static void test_stft(int& passed, int& failed) {
    std::cout << "  Testing STFT..." << std::endl;

    std::vector<double> sig(128);
    for (std::size_t i = 0; i < sig.size(); ++i) {
        sig[i] = mymath::sin(2.0 * mymath::kPi * 16.0 * static_cast<double>(i) / 128.0);
    }

    signal::STFTResult result = signal::stft(sig, 32);

    bool ok = !result.stft_matrix.empty() &&
              !result.time_vector.empty() &&
              !result.freq_vector.empty();

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: STFT test" << std::endl;
    }
}

static void test_spectrogram(int& passed, int& failed) {
    std::cout << "  Testing spectrogram..." << std::endl;

    std::vector<double> sig(128);
    for (std::size_t i = 0; i < sig.size(); ++i) {
        sig[i] = mymath::sin(2.0 * mymath::kPi * 16.0 * static_cast<double>(i) / 128.0);
    }

    std::vector<std::vector<double>> spec = signal::spectrogram(sig, 32);

    bool ok = !spec.empty() && !spec[0].empty();

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Spectrogram test" << std::endl;
    }
}

// ============================================================================
// 辅助函数测试
// ============================================================================

static void test_is_power_of_two(int& passed, int& failed) {
    std::cout << "  Testing is_power_of_two..." << std::endl;

    bool ok = signal::is_power_of_two(1) &&
              signal::is_power_of_two(2) &&
              signal::is_power_of_two(4) &&
              signal::is_power_of_two(8) &&
              signal::is_power_of_two(1024) &&
              !signal::is_power_of_two(0) &&
              !signal::is_power_of_two(3) &&
              !signal::is_power_of_two(5) &&
              !signal::is_power_of_two(6) &&
              !signal::is_power_of_two(7);

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: is_power_of_two test" << std::endl;
    }
}

static void test_next_power_of_two(int& passed, int& failed) {
    std::cout << "  Testing next_power_of_two..." << std::endl;

    bool ok = signal::next_power_of_two(0) == 1u &&
              signal::next_power_of_two(1) == 1u &&
              signal::next_power_of_two(2) == 2u &&
              signal::next_power_of_two(3) == 4u &&
              signal::next_power_of_two(5) == 8u &&
              signal::next_power_of_two(9) == 16u;

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: next_power_of_two test" << std::endl;
    }
}

static void test_gcd_lcm(int& passed, int& failed) {
    std::cout << "  Testing GCD and LCM..." << std::endl;

    bool ok = signal::gcd(12, 8) == 4u &&
              signal::gcd(17, 13) == 1u &&
              signal::lcm(12, 8) == 24u &&
              signal::lcm(3, 5) == 15u;

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: GCD/LCM test" << std::endl;
    }
}

static void test_window_type_conversion(int& passed, int& failed) {
    std::cout << "  Testing window type conversion..." << std::endl;

    bool ok = signal::window_type_to_string(signal::WindowType::Hanning) == "hanning" &&
              signal::window_type_to_string(signal::WindowType::Hamming) == "hamming" &&
              signal::window_type_to_string(signal::WindowType::Kaiser) == "kaiser" &&
              signal::string_to_window_type("hanning") == signal::WindowType::Hanning &&
              signal::string_to_window_type("hamming") == signal::WindowType::Hamming &&
              signal::string_to_window_type("kaiser") == signal::WindowType::Kaiser;

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Window type conversion test" << std::endl;
    }
}

static void test_filter_type_conversion(int& passed, int& failed) {
    std::cout << "  Testing filter type conversion..." << std::endl;

    bool ok = signal::filter_type_to_string(signal::FilterType::LowPass) == "lowpass" &&
              signal::filter_type_to_string(signal::FilterType::HighPass) == "highpass" &&
              signal::filter_type_to_string(signal::FilterType::BandPass) == "bandpass" &&
              signal::string_to_filter_type("lowpass") == signal::FilterType::LowPass &&
              signal::string_to_filter_type("highpass") == signal::FilterType::HighPass &&
              signal::string_to_filter_type("bandpass") == signal::FilterType::BandPass;

    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "    FAILED: Filter type conversion test" << std::endl;
    }
}

// ============================================================================
// 测试入口
// ============================================================================

int run_signal_processing_tests(int& passed, int& failed) {
    std::cout << "Running Signal Processing Tests..." << std::endl;

    // FFT 测试
    test_fft_radix2(passed, failed);
    test_fft_ifft_roundtrip(passed, failed);
    test_fft_non_power_of_two(passed, failed);
    test_rfft(passed, failed);
    test_fftshift(passed, failed);

    // 卷积测试
    test_convolve(passed, failed);
    test_circular_convolve(passed, failed);
    test_xcorr(passed, failed);
    test_autocorr(passed, failed);

    // 窗函数测试
    test_hanning_window(passed, failed);
    test_hamming_window(passed, failed);
    test_kaiser_window(passed, failed);
    test_window_function(passed, failed);

    // 滤波器测试
    test_fir_lowpass(passed, failed);
    test_filter_application(passed, failed);
    test_freqz(passed, failed);

    // 时频分析测试
    test_periodogram(passed, failed);
    test_pwelch(passed, failed);
    test_stft(passed, failed);
    test_spectrogram(passed, failed);

    // 辅助函数测试
    test_is_power_of_two(passed, failed);
    test_next_power_of_two(passed, failed);
    test_gcd_lcm(passed, failed);
    test_window_type_conversion(passed, failed);
    test_filter_type_conversion(passed, failed);

    return 0;
}

}  // namespace test_suites

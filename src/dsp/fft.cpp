/**
 * @file fft.cpp
 * @brief FFT 高级算法实现
 *
 * 实现：
 * - 基-2 FFT（Cooley-Tukey）
 * - 混合基 FFT
 * - Bluestein FFT
 * - 实数 FFT
 */

#include "signal_processing.h"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <vector>

namespace signal {

// 数学常量
constexpr double kPi = 3.14159265358979323846;

// ============================================================================
// 辅助函数
// ============================================================================

bool is_power_of_two(std::size_t n) {
    return n > 0 && (n & (n - 1)) == 0;
}

std::size_t next_power_of_two(std::size_t n) {
    if (n == 0) return 1;
    std::size_t result = 1;
    while (result < n) {
        result <<= 1;
    }
    return result;
}

std::size_t gcd(std::size_t a, std::size_t b) {
    while (b != 0) {
        std::size_t temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

std::size_t lcm(std::size_t a, std::size_t b) {
    return a / gcd(a, b) * b;
}

// ============================================================================
// 基-2 FFT（Cooley-Tukey 算法）
// ============================================================================

std::vector<Complex> fft_radix2(const std::vector<Complex>& input) {
    const std::size_t n = input.size();
    if (n == 0) {
        return {};
    }
    if (!is_power_of_two(n)) {
        throw std::runtime_error("fft_radix2 requires input length to be power of 2");
    }

    std::vector<Complex> output = input;

    // 位反转排列
    std::size_t reversed = 0;
    for (std::size_t i = 1; i < n; ++i) {
        std::size_t bit = n >> 1;
        while ((reversed & bit) != 0) {
            reversed ^= bit;
            bit >>= 1;
        }
        reversed ^= bit;
        if (i < reversed) {
            std::swap(output[i], output[reversed]);
        }
    }

    // Cooley-Tukey 迭代
    for (std::size_t length = 2; length <= n; length <<= 1) {
        const double angle = -2.0 * kPi / static_cast<double>(length);
        const Complex step = mymath::polar(1.0, angle);
        for (std::size_t start = 0; start < n; start += length) {
            Complex twiddle = 1.0;
            const std::size_t half = length >> 1;
            for (std::size_t offset = 0; offset < half; ++offset) {
                const Complex even = output[start + offset];
                const Complex odd = output[start + offset + half] * twiddle;
                output[start + offset] = even + odd;
                output[start + offset + half] = even - odd;
                twiddle *= step;
            }
        }
    }

    return output;
}

// ============================================================================
// 混合基 FFT
// ============================================================================

// 找到最优分解因子
static std::vector<std::size_t> find_factors(std::size_t n) {
    std::vector<std::size_t> factors;

    // 优先提取 2 的因子
    while ((n % 2) == 0 && n > 1) {
        factors.push_back(2);
        n /= 2;
    }

    // 提取其他小因子（3, 5）
    for (std::size_t p = 3; p <= 5 && n > 1; p += 2) {
        while ((n % p) == 0) {
            factors.push_back(p);
            n /= p;
        }
    }

    // 剩余部分
    if (n > 1) {
        factors.push_back(n);
    }

    return factors;
}

// 小因子 FFT（用于混合基）
static void small_factor_fft(std::vector<Complex>& data,
                              std::size_t start,
                              std::size_t factor,
                              std::size_t stride,
                              std::size_t total_length) {
    if (factor == 2) {
        // 基-2
        const double angle = -2.0 * kPi / static_cast<double>(total_length);
        Complex twiddle = mymath::polar(1.0, angle * static_cast<double>(start));
        const Complex step = mymath::polar(1.0, angle * static_cast<double>(stride));

        for (std::size_t i = 0; i < stride; ++i) {
            const Complex even = data[start + i];
            const Complex odd = data[start + i + stride] * twiddle;
            data[start + i] = even + odd;
            data[start + i + stride] = even - odd;
            twiddle *= step;
        }
    } else if (factor == 3) {
        // 基-3
        const double angle = -2.0 * kPi / static_cast<double>(total_length * factor);
        const Complex w1 = mymath::polar(1.0, angle);
        const Complex w2 = mymath::polar(1.0, 2.0 * angle);

        const double base_angle = -2.0 * kPi / static_cast<double>(total_length);
        Complex twiddle = mymath::polar(1.0, base_angle * static_cast<double>(start));
        const Complex step = mymath::polar(1.0, base_angle * static_cast<double>(stride));

        for (std::size_t i = 0; i < stride; ++i) {
            Complex sum = data[start + i];
            Complex sum1 = data[start + i + stride] * w1;
            Complex sum2 = data[start + i + 2 * stride] * w2;

            data[start + i] = twiddle * (sum + sum1 + sum2);
            data[start + i + stride] = twiddle * (sum + w2 * sum1 + w1 * sum2);
            data[start + i + 2 * stride] = twiddle * (sum + w1 * sum1 + w2 * sum2);
            twiddle *= step;
        }
    } else if (factor == 4) {
        // 基-4
        const double angle = -2.0 * kPi / static_cast<double>(total_length);
        Complex twiddle = mymath::polar(1.0, angle * static_cast<double>(start));
        const Complex step = mymath::polar(1.0, angle * static_cast<double>(stride));
        const Complex j = mymath::polar(1.0, kPi / 2.0);  // j = sqrt(-1)

        for (std::size_t i = 0; i < stride; ++i) {
            const Complex a = data[start + i];
            const Complex b = data[start + i + stride] * twiddle;
            const Complex c = data[start + i + 2 * stride] * twiddle * twiddle;
            const Complex d = data[start + i + 3 * stride] * twiddle * twiddle * twiddle;

            const Complex apc = a + c;
            const Complex amc = a - c;
            const Complex bpd = b + d;
            const Complex bmd_j = (b - d) * j;

            data[start + i] = twiddle * (apc + bpd);
            data[start + i + stride] = twiddle * (amc - bmd_j);
            data[start + i + 2 * stride] = twiddle * (apc - bpd);
            data[start + i + 3 * stride] = twiddle * (amc + bmd_j);
            twiddle *= step;
        }
    } else if (factor == 5) {
        // 基-5
        const double angle = -2.0 * kPi / 5.0;
        const Complex w1 = mymath::polar(1.0, angle);
        const Complex w2 = mymath::polar(1.0, 2.0 * angle);
        const Complex w3 = mymath::polar(1.0, 3.0 * angle);
        const Complex w4 = mymath::polar(1.0, 4.0 * angle);

        const double base_angle = -2.0 * kPi / static_cast<double>(total_length);
        Complex twiddle = mymath::polar(1.0, base_angle * static_cast<double>(start));
        const Complex step = mymath::polar(1.0, base_angle * static_cast<double>(stride));

        for (std::size_t i = 0; i < stride; ++i) {
            Complex a = data[start + i];
            Complex b = data[start + i + stride];
            Complex c = data[start + i + 2 * stride];
            Complex d = data[start + i + 3 * stride];
            Complex e = data[start + i + 4 * stride];

            Complex sum = a + b + c + d + e;
            Complex s1 = a + w1 * b + w2 * c + w3 * d + w4 * e;
            Complex s2 = a + w2 * b + w4 * c + w1 * d + w3 * e;
            Complex s3 = a + w3 * b + w1 * c + w4 * d + w2 * e;
            Complex s4 = a + w4 * b + w3 * c + w2 * d + w1 * e;

            data[start + i] = twiddle * sum;
            data[start + i + stride] = twiddle * s1;
            data[start + i + 2 * stride] = twiddle * s2;
            data[start + i + 3 * stride] = twiddle * s3;
            data[start + i + 4 * stride] = twiddle * s4;
            twiddle *= step;
        }
    }
}

// 混合基 FFT 递归实现
static void mixed_radix_fft_recursive(std::vector<Complex>& data,
                                       std::size_t start,
                                       std::size_t stride,
                                       std::size_t total_length,
                                       const std::vector<std::size_t>& factors,
                                       std::size_t factor_index) {
    if (factor_index >= factors.size()) {
        return;
    }

    const std::size_t factor = factors[factor_index];
    const std::size_t sub_length = total_length / factor;

    // 先对每个子序列递归处理
    for (std::size_t k = 0; k < factor; ++k) {
        mixed_radix_fft_recursive(data, start + k * stride, stride * factor,
                                   sub_length, factors, factor_index + 1);
    }

    // 然后应用当前层的蝶形运算
    small_factor_fft(data, start, factor, stride, total_length);
}

std::vector<Complex> fft_mixed_radix(const std::vector<Complex>& input) {
    const std::size_t n = input.size();
    if (n == 0) {
        return {};
    }
    if (n == 1) {
        return input;
    }

    // 如果是 2 的幂次，直接用基-2
    if (is_power_of_two(n)) {
        return fft_radix2(input);
    }

    // 找到因子分解
    std::vector<std::size_t> factors = find_factors(n);

    // 检查是否有大因子（不适合混合基）
    for (std::size_t f : factors) {
        if (f > 7) {
            // 大因子，使用 Bluestein 更高效
            return fft_bluestein(input);
        }
    }

    std::vector<Complex> output = input;
    mixed_radix_fft_recursive(output, 0, 1, n, factors, 0);

    return output;
}

// ============================================================================
// Bluestein FFT（Chirp-Z 变换）
// ============================================================================

std::vector<Complex> fft_bluestein(const std::vector<Complex>& input) {
    const std::size_t n = input.size();
    if (n == 0) {
        return {};
    }
    if (n == 1) {
        return input;
    }

    // 找到合适的 FFT 长度（至少 2n-1，且为 2 的幂）
    const std::size_t m = next_power_of_two(2 * n - 1);

    // 线程本地缓存以优化重复调用
    thread_local std::map<std::size_t, std::vector<Complex>> chirp_fft_cache;
    
    std::vector<Complex> chirp_fft;
    auto cache_it = chirp_fft_cache.find(n);
    
    if (cache_it != chirp_fft_cache.end()) {
        chirp_fft = cache_it->second;
    } else {
        // 生成 chirp 信号
        std::vector<Complex> chirp(m, 0.0);
        for (std::size_t k = 0; k < n; ++k) {
            const double angle = -kPi * static_cast<double>(k * k) / static_cast<double>(n);
            chirp[k] = mymath::polar(1.0, angle);
        }
        for (std::size_t k = m - n + 1; k < m; ++k) {
            const std::size_t kk = k - m;  // 负索引
            const double angle = -kPi * static_cast<double>(kk * kk) / static_cast<double>(n);
            chirp[k] = mymath::polar(1.0, angle);
        }
        chirp_fft = fft_radix2(chirp);
        chirp_fft_cache[n] = chirp_fft;
    }

    // 输入信号乘以 chirp
    std::vector<Complex> y(m, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
        const double angle = kPi * static_cast<double>(k * k) / static_cast<double>(n);
        y[k] = input[k] * mymath::polar(1.0, angle);
    }

    // FFT(y)
    std::vector<Complex> y_fft = fft_radix2(y);

    // 逐点相乘
    for (std::size_t k = 0; k < m; ++k) {
        y_fft[k] *= chirp_fft[k];
    }

    // IFFT
    std::vector<Complex> result = ifft(y_fft);

    // 取前 n 个点并乘以 chirp
    std::vector<Complex> output(n);
    for (std::size_t k = 0; k < n; ++k) {
        const double angle = kPi * static_cast<double>(k * k) / static_cast<double>(n);
        output[k] = result[k] * mymath::polar(1.0, angle);
    }

    return output;
}

// ============================================================================
// 自动选择最优 FFT
// ============================================================================

std::vector<Complex> fft(const std::vector<Complex>& input) {
    const std::size_t n = input.size();
    if (n == 0) {
        return {};
    }
    if (n == 1) {
        return input;
    }

    if (is_power_of_two(n)) {
        return fft_radix2(input);
    }

    // 检查因子分解是否适合混合基
    std::vector<std::size_t> factors = find_factors(n);
    bool suitable_for_mixed = true;
    for (std::size_t f : factors) {
        if (f > 7) {
            suitable_for_mixed = false;
            break;
        }
    }

    if (suitable_for_mixed) {
        return fft_mixed_radix(input);
    }

    return fft_bluestein(input);
}

// ============================================================================
// 逆 FFT
// ============================================================================

std::vector<Complex> ifft(const std::vector<Complex>& spectrum) {
    const std::size_t n = spectrum.size();
    if (n == 0) {
        return {};
    }

    // 共轭输入
    std::vector<Complex> conjugated(n);
    for (std::size_t i = 0; i < n; ++i) {
        conjugated[i] = mymath::conj(spectrum[i]);
    }

    // 正向 FFT
    std::vector<Complex> result = fft(conjugated);

    // 共轭并缩放
    const double scale = 1.0 / static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = mymath::conj(result[i]) * scale;
    }

    return result;
}

// ============================================================================
// 实数 FFT
// ============================================================================

std::vector<Complex> rfft(const std::vector<double>& input) {
    const std::size_t n = input.size();
    if (n == 0) return {};
    if (n == 1) return {Complex(input[0], 0.0)};

    if (n % 2 == 0) {
        // 针对偶数长度的高效打包算法：将实数打包为复数
        std::vector<Complex> packed(n / 2);
        for (std::size_t i = 0; i < n / 2; ++i) {
            packed[i] = Complex(input[2 * i], input[2 * i + 1]);
        }

        std::vector<Complex> packed_fft = fft(packed);

        // 解包利用 Hermitian 对称性得到正频率部分 (0 到 n/2)
        std::vector<Complex> result(n / 2 + 1);
        
        // 处理 DC 和 Nyquist (如果是偶数)
        const Complex f0 = packed_fft[0];
        result[0] = Complex(f0.real() + f0.imag(), 0.0);
        result[n / 2] = Complex(f0.real() - f0.imag(), 0.0);

        for (std::size_t k = 1; k < n / 2; ++k) {
            const Complex fk = packed_fft[k];
            const Complex fnk = mymath::conj(packed_fft[n / 2 - k]);
            
            const Complex fe = 0.5 * (fk + fnk);
            const Complex fo = Complex(0.0, -0.5) * (fk - fnk);
            
            const double angle = -2.0 * kPi * static_cast<double>(k) / static_cast<double>(n);
            const Complex twiddle = mymath::polar(1.0, angle);
            
            result[k] = fe + twiddle * fo;
        }
        return result;
    }

    // 奇数长度：回退到标准 FFT
    std::vector<Complex> complex_input(n);
    for (std::size_t i = 0; i < n; ++i) {
        complex_input[i] = Complex(input[i], 0.0);
    }
    std::vector<Complex> full_fft = fft(complex_input);
    std::vector<Complex> result(n / 2 + 1);
    std::copy(full_fft.begin(), full_fft.begin() + (n / 2 + 1), result.begin());
    return result;
}

std::vector<double> irfft(const std::vector<Complex>& spectrum, std::size_t n) {
    if (spectrum.empty()) return {};
    if (n == 0) return {};
    
    if (n % 2 == 0 && spectrum.size() == n / 2 + 1) {
        // 逆过程：解包并使用复数 IFFT
        std::vector<Complex> packed(n / 2);
        
        for (std::size_t k = 0; k < n / 2; ++k) {
            const Complex rk = spectrum[k];
            const Complex rnk = mymath::conj(spectrum[n / 2 - k]);
            
            const Complex fe = 0.5 * (rk + rnk);
            const Complex fo = 0.5 * (rk - rnk);
            
            const double angle = 2.0 * kPi * static_cast<double>(k) / static_cast<double>(n);
            const Complex twiddle_inv = mymath::polar(1.0, angle);
            
            packed[k] = fe + Complex(0.0, 1.0) * twiddle_inv * fo;
        }
        
        std::vector<Complex> time_complex = ifft(packed);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n / 2; ++i) {
            result[2 * i] = time_complex[i].real();
            result[2 * i + 1] = time_complex[i].imag();
        }
        return result;
    }

    // 回退到标准 IFFT
    std::vector<Complex> full_spectrum(n);
    for (std::size_t i = 0; i < spectrum.size(); ++i) full_spectrum[i] = spectrum[i];
    for (std::size_t i = spectrum.size(); i < n; ++i) full_spectrum[i] = mymath::conj(spectrum[n - i]);
    
    std::vector<Complex> time_signal = ifft(full_spectrum);
    std::vector<double> result(n);
    for (std::size_t i = 0; i < n; ++i) result[i] = time_signal[i].real();
    return result;
}

// ============================================================================
// FFT 频率轴和频移
// ============================================================================

std::vector<double> fft_frequencies(std::size_t n, double sample_rate) {
    std::vector<double> freqs(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i <= n / 2) {
            freqs[i] = static_cast<double>(i) * sample_rate / static_cast<double>(n);
        } else {
            freqs[i] = static_cast<double>(i - n) * sample_rate / static_cast<double>(n);
        }
    }
    return freqs;
}

std::vector<Complex> fftshift(const std::vector<Complex>& spectrum) {
    const std::size_t n = spectrum.size();
    if (n == 0) {
        return {};
    }

    std::vector<Complex> shifted(n);
    const std::size_t half = n / 2;

    for (std::size_t i = 0; i < half; ++i) {
        shifted[i] = spectrum[i + half + (n % 2)];
    }
    for (std::size_t i = half; i < n; ++i) {
        shifted[i] = spectrum[i - half];
    }

    return shifted;
}

std::vector<Complex> ifftshift(const std::vector<Complex>& spectrum) {
    const std::size_t n = spectrum.size();
    if (n == 0) {
        return {};
    }

    std::vector<Complex> shifted(n);
    const std::size_t half = n / 2;

    for (std::size_t i = 0; i < half + (n % 2); ++i) {
        shifted[i] = spectrum[i + half];
    }
    for (std::size_t i = half + (n % 2); i < n; ++i) {
        shifted[i] = spectrum[i - half - (n % 2)];
    }

    return shifted;
}

}  // namespace signal

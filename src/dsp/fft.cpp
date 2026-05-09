/**
 * @file fft.cpp
 * @brief FFT 高级算法实现
 *
 * 实现：
 * - 基-2 FFT（Cooley-Tukey）
 * - 混合基 FFT（迭代实现）
 * - Bluestein FFT（带 LRU 缓存）
 * - 实数 FFT
 */

#include "signal_processing.h"
#include "math/mymath.h"
#include "core/scalar_type.h"

#include <algorithm>
#include <list>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace signal {

// Use Scalar for internal high-precision real computations
using Scalar = mymath::Scalar;

// 使用 mymath 命名空间中的常量
using mymath::kPi;
static const Scalar kPiScalar = mymath::precise128::pi();

// ============================================================================
// LRU 缓存（用于 Bluestein FFT）
// ============================================================================

/**
 * @brief Chirp FFT 缓存的 LRU 缓存管理器
 *
 * 解决原 thread_local std::map 无上限导致的内存膨胀问题。
 * 最大缓存 16 个不同长度的 Chirp FFT 结果。
 */
class ChirpFftLruCache {
public:
    static constexpr std::size_t kMaxCapacity = 16;

    std::vector<Complex> get(std::size_t n) {
        auto it = index_.find(n);
        if (it != index_.end()) {
            // 移动到前端（最近使用）
            entries_.splice(entries_.begin(), entries_, it->second);
            return it->second->second;
        }
        return {};  // 未命中
    }

    void put(std::size_t n, const std::vector<Complex>& chirp_fft) {
        auto it = index_.find(n);
        if (it != index_.end()) {
            // 已存在，更新并移到前端
            it->second->second = chirp_fft;
            entries_.splice(entries_.begin(), entries_, it->second);
            return;
        }

        // 淘汰最久未使用
        if (entries_.size() >= kMaxCapacity) {
            index_.erase(entries_.back().first);
            entries_.pop_back();
        }

        entries_.emplace_front(n, chirp_fft);
        index_[n] = entries_.begin();
    }

private:
    std::list<std::pair<std::size_t, std::vector<Complex>>> entries_;
    std::unordered_map<std::size_t,
                       std::list<std::pair<std::size_t, std::vector<Complex>>>::iterator> index_;
};

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

// Helper to compute twiddle factor with high precision
static inline Complex twiddle_factor(Scalar angle) {
    return Complex(mymath::cos(angle), mymath::sin(angle));
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

    // Cooley-Tukey 迭代 - 使用预计算旋转因子表提高数值稳定性
    for (std::size_t length = 2; length <= n; length <<= 1) {
        const Scalar base_angle = -Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(length));
        const std::size_t half = length >> 1;

        // 预计算旋转因子表，避免累乘误差
        std::vector<Complex> twiddles(half);
        for (std::size_t k = 0; k < half; ++k) {
            twiddles[k] = twiddle_factor(base_angle * Scalar(static_cast<long long>(k)));
        }

        for (std::size_t start = 0; start < n; start += length) {
            for (std::size_t offset = 0; offset < half; ++offset) {
                const Complex even = output[start + offset];
                const Complex odd = output[start + offset + half] * twiddles[offset];
                output[start + offset] = even + odd;
                output[start + offset + half] = even - odd;
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

// 小因子 FFT（用于混合基）- 使用预计算旋转因子表
static void small_factor_fft(std::vector<Complex>& data,
                              std::size_t start,
                              std::size_t factor,
                              std::size_t stride,
                              std::size_t total_length,
                              const std::vector<Complex>& twiddles_base) {
    if (factor == 2) {
        // 基-2
        const Scalar base_angle = -Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(total_length));
        Complex twiddle = twiddle_factor(base_angle * Scalar(static_cast<long long>(start)));
        const Complex step = twiddle_factor(base_angle * Scalar(static_cast<long long>(stride)));

        for (std::size_t i = 0; i < stride; ++i) {
            const Complex even = data[start + i];
            const Complex odd = data[start + i + stride] * twiddle;
            data[start + i] = even + odd;
            data[start + i + stride] = even - odd;
            twiddle *= step;
        }
    } else if (factor == 3) {
        // 基-3
        const Scalar angle = -Scalar(2.0L) * kPiScalar / Scalar(3.0L);
        const Complex w1 = twiddle_factor(angle);
        const Complex w2 = twiddle_factor(2.0L * angle);

        const Scalar base_angle = -Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(total_length));
        Complex twiddle = twiddle_factor(base_angle * Scalar(static_cast<long long>(start)));
        const Complex step = twiddle_factor(base_angle * Scalar(static_cast<long long>(stride)));

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
        const Scalar base_angle = -Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(total_length));
        Complex twiddle = twiddle_factor(base_angle * Scalar(static_cast<long long>(start)));
        const Complex step = twiddle_factor(base_angle * Scalar(static_cast<long long>(stride)));
        const Complex j = twiddle_factor(kPiScalar / Scalar(2.0L));  // j = sqrt(-1)

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
        const Scalar angle = -Scalar(2.0L) * kPiScalar / Scalar(5.0L);
        const Complex w1 = twiddle_factor(angle);
        const Complex w2 = twiddle_factor(2.0L * angle);
        const Complex w3 = twiddle_factor(3.0L * angle);
        const Complex w4 = twiddle_factor(4.0L * angle);

        const Scalar base_angle = -Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(total_length));
        Complex twiddle = twiddle_factor(base_angle * Scalar(static_cast<long long>(start)));
        const Complex step = twiddle_factor(base_angle * Scalar(static_cast<long long>(stride)));

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

// 混合基 FFT 迭代实现（避免递归栈溢出）
static void mixed_radix_fft_iterative(std::vector<Complex>& data,
                                       std::size_t n,
                                       const std::vector<std::size_t>& factors) {
    if (factors.empty() || n <= 1) return;

    // 计算各层的 stride 和 total_length
    // 使用迭代方式从最内层向外处理
    const std::size_t num_stages = factors.size();

    // 预计算每层的参数
    std::vector<std::size_t> stage_lengths(num_stages + 1);
    stage_lengths[0] = n;
    for (std::size_t i = 0; i < num_stages; ++i) {
        stage_lengths[i + 1] = stage_lengths[i] / factors[i];
    }

    // 预计算旋转因子表（空表，small_factor_fft 会自行计算）
    std::vector<Complex> twiddles_base;

    // 从最内层向外迭代处理
    // 使用栈模拟递归，但按层处理避免深度递归
    for (std::size_t stage = 0; stage < num_stages; ++stage) {
        const std::size_t factor = factors[stage];
        const std::size_t total_length = stage_lengths[stage];
        const std::size_t sub_length = stage_lengths[stage + 1];

        // 计算当前层的 stride
        std::size_t stride = 1;
        for (std::size_t j = stage + 1; j < num_stages; ++j) {
            stride *= factors[j];
        }

        // 对当前层的所有块应用蝶形运算
        const std::size_t num_blocks = n / total_length;
        for (std::size_t block = 0; block < num_blocks; ++block) {
            const std::size_t start = block * total_length;
            small_factor_fft(data, start, factor, stride, total_length, twiddles_base);
        }
    }
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
    mixed_radix_fft_iterative(output, n, factors);

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
    const std::size_t m = next_power_of_two(2 * static_cast<long long>(n) - 1);

    // 使用 LRU 缓存替代无上限的 thread_local std::map
    static thread_local ChirpFftLruCache chirp_fft_cache;

    std::vector<Complex> chirp_fft;
    chirp_fft = chirp_fft_cache.get(n);

    if (chirp_fft.empty()) {
        // 生成 chirp 信号
        std::vector<Complex> chirp(m, Complex(0.0L, 0.0L));
        for (std::size_t k = 0; k < n; ++k) {
            const Scalar angle = -kPiScalar * Scalar((k * static_cast<long long>(k))) / Scalar(static_cast<long long>(n));
            chirp[k] = twiddle_factor(angle);
        }
        for (std::size_t k = m - n + 1; k < m; ++k) {
            const std::size_t kk = k - m;  // 负索引
            const Scalar angle = -kPiScalar * Scalar((kk * kk)) / Scalar(static_cast<long long>(n));
            chirp[k] = twiddle_factor(angle);
        }
        chirp_fft = fft_radix2(chirp);
        chirp_fft_cache.put(n, chirp_fft);
    }

    // 输入信号乘以 chirp
    std::vector<Complex> y(m, Complex(0.0L, 0.0L));
    for (std::size_t k = 0; k < n; ++k) {
        const Scalar angle = kPiScalar * Scalar((k * static_cast<long long>(k))) / Scalar(static_cast<long long>(n));
        y[k] = input[k] * twiddle_factor(angle);
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
        const Scalar angle = kPiScalar * Scalar((k * static_cast<long long>(k))) / Scalar(static_cast<long long>(n));
        output[k] = result[k] * twiddle_factor(angle);
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
    const Scalar scale = Scalar(1.0L) / Scalar(static_cast<long long>(n));
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = mymath::conj(result[i]) * scale;
    }

    return result;
}

// ============================================================================
// 实数 FFT
// ============================================================================

std::vector<Complex> rfft(const std::vector<Scalar>& input) {
    const std::size_t n = input.size();
    if (n == 0) return {};
    if (n == 1) return {Complex(input[0], 0.0L)};

    if (n % 2 == 0) {
        // 针对偶数长度的高效打包算法：将实数打包为复数
        std::vector<Complex> packed(n / 2);
        for (std::size_t i = 0; i < n / 2; ++i) {
            packed[i] = Complex(input[2 * static_cast<long long>(i)], input[2 * static_cast<long long>(i) + 1]);
        }

        std::vector<Complex> packed_fft = fft(packed);

        // 解包利用 Hermitian 对称性得到正频率部分 (0 到 n/2)
        std::vector<Complex> result(n / 2 + 1);

        // 处理 DC 和 Nyquist (如果是偶数)
        const Complex f0 = packed_fft[0];
        result[0] = Complex(f0.real() + f0.imag(), 0.0L);
        result[n / 2] = Complex(f0.real() - f0.imag(), 0.0L);

        for (std::size_t k = 1; k < n / 2; ++k) {
            const Complex fk = packed_fft[k];
            const Complex fnk = mymath::conj(packed_fft[n / 2 - k]);

            const Complex fe = 0.5L * (fk + fnk);
            const Complex fo = Complex(0.0L, -0.5L) * (fk - fnk);

            const Scalar angle_scalar = -Scalar(2.0L) * kPiScalar * Scalar(static_cast<long long>(k)) / Scalar(static_cast<long long>(n));
            const Complex twiddle = twiddle_factor((angle_scalar));

            result[k] = fe + twiddle * fo;
        }
        return result;
    }

    // 奇数长度：回退到标准 FFT
    std::vector<Complex> complex_input(n);
    for (std::size_t i = 0; i < n; ++i) {
        complex_input[i] = Complex(input[i], 0.0L);
    }
    std::vector<Complex> full_fft = fft(complex_input);
    std::vector<Complex> result(n / 2 + 1);
    std::copy(full_fft.begin(), full_fft.begin() + (n / 2 + 1), result.begin());
    return result;
}

std::vector<Scalar> irfft(const std::vector<Complex>& spectrum, std::size_t n) {
    if (spectrum.empty()) return {};
    if (n == 0) return {};

    if (n % 2 == 0 && spectrum.size() == n / 2 + 1) {
        // 逆过程：解包并使用复数 IFFT
        std::vector<Complex> packed(n / 2);

        for (std::size_t k = 0; k < n / 2; ++k) {
            const Complex rk = spectrum[k];
            const Complex rnk = mymath::conj(spectrum[n / 2 - k]);

            const Complex fe = 0.5L * (rk + rnk);
            const Complex fo = 0.5L * (rk - rnk);

            const Scalar angle_scalar = Scalar(2.0L) * kPiScalar * Scalar(static_cast<long long>(k)) / Scalar(static_cast<long long>(n));
            const Complex twiddle_inv = twiddle_factor((angle_scalar));

            packed[k] = fe + Complex(0.0L, 1.0L) * twiddle_inv * fo;
        }

        std::vector<Complex> time_complex = ifft(packed);
        std::vector<Scalar> result(n);
        for (std::size_t i = 0; i < n / 2; ++i) {
            result[2 * static_cast<long long>(i)] = time_complex[i].real();
            result[2 * static_cast<long long>(i) + 1] = time_complex[i].imag();
        }
        return result;
    }

    // 回退到标准 IFFT
    std::vector<Complex> full_spectrum(n);
    for (std::size_t i = 0; i < spectrum.size(); ++i) full_spectrum[i] = spectrum[i];
    for (std::size_t i = spectrum.size(); i < n; ++i) full_spectrum[i] = mymath::conj(spectrum[n - i]);

    std::vector<Complex> time_signal = ifft(full_spectrum);
    std::vector<Scalar> result(n);
    for (std::size_t i = 0; i < n; ++i) result[i] = time_signal[i].real();
    return result;
}

// ============================================================================
// FFT 频率轴和频移
// ============================================================================

std::vector<Scalar> fft_frequencies(std::size_t n, Scalar sample_rate) {
    std::vector<Scalar> freqs(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i <= n / 2) {
            freqs[i] = (i) * sample_rate / (n);
        } else {
            freqs[i] = (i - n) * sample_rate / (n);
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

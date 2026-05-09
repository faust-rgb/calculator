/**
 * @file window_functions.cpp
 * @brief 窗函数实现
 *
 * 实现各种窗函数：
 * - 矩形窗、汉宁窗、汉明窗
 * - 布莱克曼窗、布莱克曼-哈里斯窗
 * - 巴特利特窗、凯撒窗、高斯窗
 * - 平顶窗、Tukey 窗
 */

#include "signal_processing.h"
#include "math/mymath.h"
#include "core/common/scalar_type.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace signal {

// Use Scalar for internal high-precision real computations
using Scalar = mymath::Scalar;

// 使用 mymath 命名空间中的常量
using mymath::kPi;
static const Scalar kPiScalar = mymath::precise128::pi();

// ============================================================================
// 窗函数类型转换
// ============================================================================

std::string window_type_to_string(WindowType type) {
    switch (type) {
        case WindowType::Rectangular:    return "rectangular";
        case WindowType::Hanning:        return "hanning";
        case WindowType::Hamming:        return "hamming";
        case WindowType::Blackman:       return "blackman";
        case WindowType::BlackmanHarris: return "blackman-harris";
        case WindowType::Bartlett:       return "bartlett";
        case WindowType::Kaiser:         return "kaiser";
        case WindowType::Gaussian:       return "gaussian";
        case WindowType::Flattop:        return "flattop";
        case WindowType::Tukey:          return "tukey";
        default:                         return "unknown";
    }
}

WindowType string_to_window_type(const std::string& str) {
    if (str == "rectangular" || str == "rect") {
        return WindowType::Rectangular;
    } else if (str == "hanning" || str == "hann") {
        return WindowType::Hanning;
    } else if (str == "hamming") {
        return WindowType::Hamming;
    } else if (str == "blackman") {
        return WindowType::Blackman;
    } else if (str == "blackman-harris" || str == "blackmanharris") {
        return WindowType::BlackmanHarris;
    } else if (str == "bartlett" || str == "triangular") {
        return WindowType::Bartlett;
    } else if (str == "kaiser") {
        return WindowType::Kaiser;
    } else if (str == "gaussian" || str == "gauss") {
        return WindowType::Gaussian;
    } else if (str == "flattop" || str == "flat") {
        return WindowType::Flattop;
    } else if (str == "tukey") {
        return WindowType::Tukey;
    }
    throw std::runtime_error("Unknown window type: " + str);
}

// ============================================================================
// 矩形窗
// ============================================================================

std::vector<Scalar> rectangular_window(std::size_t length) {
    return std::vector<Scalar>(length, 1.0L);
}

// ============================================================================
// 汉宁窗（Hanning / Hann）
// ============================================================================

std::vector<Scalar> hanning_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<Scalar> w(length);
    const Scalar factor = Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(length - 1));

    for (std::size_t i = 0; i < length; ++i) {
        const Scalar x = factor * Scalar(static_cast<long long>(i));
        w[i] = (Scalar(0.5L) * (Scalar(1.0L) - mymath::precise128::cos(x)));
    }

    return w;
}

// ============================================================================
// 汉明窗（Hamming）
// ============================================================================

std::vector<Scalar> hamming_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<Scalar> w(length);
    const Scalar factor = Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(length - 1));

    // Hamming 窗系数：a0 = 0.54, a1 = 0.46
    constexpr Scalar a0 = Scalar(0.54L);
    constexpr Scalar a1 = Scalar(0.46L);

    for (std::size_t i = 0; i < length; ++i) {
        const Scalar x = factor * Scalar(static_cast<long long>(i));
        w[i] = (a0 - a1 * mymath::precise128::cos(x));
    }

    return w;
}

// ============================================================================
// 布莱克曼窗（Blackman）
// ============================================================================

std::vector<Scalar> blackman_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<Scalar> w(length);
    const Scalar factor = Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(length - 1));

    // Blackman 窗系数
    constexpr Scalar a0 = Scalar(0.42L);
    constexpr Scalar a1 = Scalar(0.5L);
    constexpr Scalar a2 = Scalar(0.08L);

    for (std::size_t i = 0; i < length; ++i) {
        const Scalar x = factor * Scalar(static_cast<long long>(i));
        w[i] = (a0 - a1 * mymath::precise128::cos(x) +
                                        a2 * mymath::precise128::cos(Scalar(2.0L) * x));
    }

    return w;
}

// ============================================================================
// 布莱克曼-哈里斯窗（Blackman-Harris）
// ============================================================================

std::vector<Scalar> blackman_harris_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<Scalar> w(length);
    const Scalar factor = Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(length - 1));

    // 4 项 Blackman-Harris 窗系数
    constexpr Scalar a0 = Scalar(0.35875L);
    constexpr Scalar a1 = Scalar(0.48829L);
    constexpr Scalar a2 = Scalar(0.14128L);
    constexpr Scalar a3 = Scalar(0.01168L);

    for (std::size_t i = 0; i < length; ++i) {
        const Scalar x = factor * Scalar(static_cast<long long>(i));
        w[i] = (a0 - a1 * mymath::precise128::cos(x) +
                                        a2 * mymath::precise128::cos(Scalar(2.0L) * x) -
                                        a3 * mymath::precise128::cos(Scalar(3.0L) * x));
    }

    return w;
}

// ============================================================================
// 巴特利特窗（Bartlett / 三角窗）
// ============================================================================

std::vector<Scalar> bartlett_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<Scalar> w(length);
    const Scalar half = Scalar(static_cast<long long>(length - 1)) / Scalar(2.0L);

    for (std::size_t i = 0; i < length; ++i) {
        const Scalar idx = Scalar(static_cast<long long>(i));
        w[i] = (Scalar(1.0L) - mymath::precise128::abs(idx - half) / half);
    }

    return w;
}

// ============================================================================
// 凯撒窗（Kaiser）
// ============================================================================

// 第一类修正贝塞尔函数 I0
static Scalar bessel_i0_scalar(Scalar x) {
    // 使用级数展开
    Scalar result = Scalar(1.0L);
    Scalar term = Scalar(1.0L);
    Scalar x_squared_over_4 = x * x / Scalar(4.0L);

    for (int k = 1; k <= 25; ++k) {
        term *= x_squared_over_4 / Scalar((k * k));
        result += term;
        if (term < Scalar(1e-12L) * result) {
            break;
        }
    }

    return result;
}

std::vector<Scalar> kaiser_window(std::size_t length, Scalar beta) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<Scalar> w(length);
    const Scalar beta_scalar = Scalar(beta);
    const Scalar i0_beta = bessel_i0_scalar(beta_scalar);
    const Scalar half = Scalar(static_cast<long long>(length - 1)) / Scalar(2.0L);

    for (std::size_t i = 0; i < length; ++i) {
        const Scalar x = (Scalar(static_cast<long long>(i)) - half) / half;
        w[i] = (bessel_i0_scalar(beta_scalar *
                                    mymath::precise128::sqrt(Scalar(1.0L) - x * x)) / i0_beta);
    }

    return w;
}

// ============================================================================
// 高斯窗（Gaussian）
// ============================================================================

std::vector<Scalar> gaussian_window(std::size_t length, Scalar sigma) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<Scalar> w(length);
    const Scalar half = Scalar(static_cast<long long>(length - 1)) / Scalar(2.0L);
    const Scalar sigma_scalar = Scalar(sigma);
    const Scalar sigma_scaled = sigma_scalar * half;

    for (std::size_t i = 0; i < length; ++i) {
        const Scalar x = Scalar(static_cast<long long>(i)) - half;
        const Scalar ratio = x / sigma_scaled;
        w[i] = (mymath::precise128::exp(Scalar(-0.5L) * ratio * ratio));
    }

    return w;
}

// ============================================================================
// 平顶窗（Flattop）
// ============================================================================

std::vector<Scalar> flattop_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<Scalar> w(length);
    const Scalar factor = Scalar(2.0L) * kPiScalar / Scalar(static_cast<long long>(length - 1));

    // 平顶窗系数
    constexpr Scalar a0 = Scalar(0.21557895L);
    constexpr Scalar a1 = Scalar(0.41663158L);
    constexpr Scalar a2 = Scalar(0.277263158L);
    constexpr Scalar a3 = Scalar(0.083578947L);
    constexpr Scalar a4 = Scalar(0.006947368L);

    for (std::size_t i = 0; i < length; ++i) {
        const Scalar x = factor * Scalar(static_cast<long long>(i));
        w[i] = (a0 - a1 * mymath::precise128::cos(x) +
                                        a2 * mymath::precise128::cos(Scalar(2.0L) * x) -
                                        a3 * mymath::precise128::cos(Scalar(3.0L) * x) +
                                        a4 * mymath::precise128::cos(Scalar(4.0L) * x));
    }

    return w;
}

// ============================================================================
// Tukey 窗
// ============================================================================

std::vector<Scalar> tukey_window(std::size_t length, Scalar alpha) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    // 限制 alpha 在 [0, 1] 范围内
    Scalar alpha_scalar = Scalar(alpha);
    alpha_scalar = mymath::precise128::fmax(Scalar(0.0L), mymath::precise128::fmin(Scalar(1.0L), alpha_scalar));

    std::vector<Scalar> w(length);

    if (alpha_scalar == Scalar(0.0L)) {
        // 退化为矩形窗
        return rectangular_window(length);
    } else if (alpha_scalar == Scalar(1.0L)) {
        // 退化为汉宁窗
        return hanning_window(length);
    }

    const Scalar width = alpha_scalar * Scalar(static_cast<long long>(length - 1)) / Scalar(2.0L);

    for (std::size_t i = 0; i < length; ++i) {
        const Scalar x = Scalar(static_cast<long long>(i));

        if (x < width) {
            // 上升沿
            w[i] = (Scalar(0.5L) * (Scalar(1.0L) +
                    mymath::precise128::cos(kPiScalar * (x / width - Scalar(1.0L)))));
        } else if (x > Scalar(static_cast<long long>(length - 1)) - width) {
            // 下降沿
            w[i] = (Scalar(0.5L) * (Scalar(1.0L) +
                    mymath::precise128::cos(kPiScalar * ((x - Scalar(static_cast<long long>(length - 1))) / width +
                    Scalar(1.0L)))));
        } else {
            // 平坦部分
            w[i] = 1.0L;
        }
    }

    return w;
}

// ============================================================================
// 通用窗函数生成
// ============================================================================

std::vector<Scalar> window(WindowType type, std::size_t length, Scalar param) {
    switch (type) {
        case WindowType::Rectangular:
            return rectangular_window(length);
        case WindowType::Hanning:
            return hanning_window(length);
        case WindowType::Hamming:
            return hamming_window(length);
        case WindowType::Blackman:
            return blackman_window(length);
        case WindowType::BlackmanHarris:
            return blackman_harris_window(length);
        case WindowType::Bartlett:
            return bartlett_window(length);
        case WindowType::Kaiser:
            return kaiser_window(length, param > 0.0L ? param : 5.0);
        case WindowType::Gaussian:
            return gaussian_window(length, param > 0.0L ? param : 0.4);
        case WindowType::Flattop:
            return flattop_window(length);
        case WindowType::Tukey:
            return tukey_window(length, param > 0.0L ? param : 0.5);
        default:
            throw std::runtime_error("Unknown window type");
    }
}

// ============================================================================
// 应用窗函数
// ============================================================================

std::vector<Scalar> apply_window(const std::vector<Scalar>& signal,
                                  WindowType window_type,
                                  Scalar param) {
    if (signal.empty()) {
        return {};
    }

    std::vector<Scalar> win = window(window_type, signal.size(), param);
    std::vector<Scalar> result(signal.size());

    // Use Scalar for multiplication to get higher precision
    for (std::size_t i = 0; i < signal.size(); ++i) {
        result[i] = (Scalar(signal[i]) * Scalar(win[i]));
    }

    return result;
}

}  // namespace signal
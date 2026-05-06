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

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace signal {

// 数学常量
constexpr long double kPi = 3.14159265358979323846;

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

std::vector<long double> rectangular_window(std::size_t length) {
    return std::vector<long double>(length, 1.0L);
}

// ============================================================================
// 汉宁窗（Hanning / Hann）
// ============================================================================

std::vector<long double> hanning_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<long double> w(length);
    const long double factor = 2.0 * kPi / static_cast<long double>(length - 1);

    for (std::size_t i = 0; i < length; ++i) {
        w[i] = 0.5 * (1.0L - mymath::cos(factor * static_cast<long double>(i)));
    }

    return w;
}

// ============================================================================
// 汉明窗（Hamming）
// ============================================================================

std::vector<long double> hamming_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<long double> w(length);
    const long double factor = 2.0 * kPi / static_cast<long double>(length - 1);

    // Hamming 窗系数：a0 = 0.54, a1 = 0.46
    for (std::size_t i = 0; i < length; ++i) {
        w[i] = 0.54 - 0.46 * mymath::cos(factor * static_cast<long double>(i));
    }

    return w;
}

// ============================================================================
// 布莱克曼窗（Blackman）
// ============================================================================

std::vector<long double> blackman_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<long double> w(length);
    const long double factor = 2.0 * kPi / static_cast<long double>(length - 1);

    // Blackman 窗系数
    constexpr long double a0 = 0.42;
    constexpr long double a1 = 0.5;
    constexpr long double a2 = 0.08;

    for (std::size_t i = 0; i < length; ++i) {
        const long double x = factor * static_cast<long double>(i);
        w[i] = a0 - a1 * mymath::cos(x) + a2 * mymath::cos(2.0 * x);
    }

    return w;
}

// ============================================================================
// 布莱克曼-哈里斯窗（Blackman-Harris）
// ============================================================================

std::vector<long double> blackman_harris_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<long double> w(length);
    const long double factor = 2.0 * kPi / static_cast<long double>(length - 1);

    // 4 项 Blackman-Harris 窗系数
    constexpr long double a0 = 0.35875;
    constexpr long double a1 = 0.48829;
    constexpr long double a2 = 0.14128;
    constexpr long double a3 = 0.01168;

    for (std::size_t i = 0; i < length; ++i) {
        const long double x = factor * static_cast<long double>(i);
        w[i] = a0 - a1 * mymath::cos(x) + a2 * mymath::cos(2.0 * x) - a3 * mymath::cos(3.0 * x);
    }

    return w;
}

// ============================================================================
// 巴特利特窗（Bartlett / 三角窗）
// ============================================================================

std::vector<long double> bartlett_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<long double> w(length);
    const long double half = static_cast<long double>(length - 1) / 2.0;

    for (std::size_t i = 0; i < length; ++i) {
        w[i] = 1.0L - mymath::abs(static_cast<long double>(i) - half) / half;
    }

    return w;
}

// ============================================================================
// 凯撒窗（Kaiser）
// ============================================================================

// 第一类修正贝塞尔函数 I0
static long double bessel_i0(long double x) {
    // 使用级数展开
    long double result = 1.0L;
    long double term = 1.0L;
    long double x_squared_over_4 = x * x / 4.0;

    for (int k = 1; k <= 25; ++k) {
        term *= x_squared_over_4 / static_cast<long double>(k * k);
        result += term;
        if (term < 1e-12 * result) {
            break;
        }
    }

    return result;
}

std::vector<long double> kaiser_window(std::size_t length, long double beta) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<long double> w(length);
    const long double i0_beta = bessel_i0(beta);
    const long double half = static_cast<long double>(length - 1) / 2.0;

    for (std::size_t i = 0; i < length; ++i) {
        const long double x = (static_cast<long double>(i) - half) / half;
        w[i] = bessel_i0(beta * mymath::sqrt(1.0L - x * x)) / i0_beta;
    }

    return w;
}

// ============================================================================
// 高斯窗（Gaussian）
// ============================================================================

std::vector<long double> gaussian_window(std::size_t length, long double sigma) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<long double> w(length);
    const long double half = static_cast<long double>(length - 1) / 2.0;
    const long double sigma_scaled = sigma * half;

    for (std::size_t i = 0; i < length; ++i) {
        const long double x = static_cast<long double>(i) - half;
        w[i] = mymath::exp(-0.5 * (x / sigma_scaled) * (x / sigma_scaled));
    }

    return w;
}

// ============================================================================
// 平顶窗（Flattop）
// ============================================================================

std::vector<long double> flattop_window(std::size_t length) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    std::vector<long double> w(length);
    const long double factor = 2.0 * kPi / static_cast<long double>(length - 1);

    // 平顶窗系数
    constexpr long double a0 = 0.21557895;
    constexpr long double a1 = 0.41663158;
    constexpr long double a2 = 0.277263158;
    constexpr long double a3 = 0.083578947;
    constexpr long double a4 = 0.006947368;

    for (std::size_t i = 0; i < length; ++i) {
        const long double x = factor * static_cast<long double>(i);
        w[i] = a0 - a1 * mymath::cos(x) + a2 * mymath::cos(2.0 * x)
               - a3 * mymath::cos(3.0 * x) + a4 * mymath::cos(4.0 * x);
    }

    return w;
}

// ============================================================================
// Tukey 窗
// ============================================================================

std::vector<long double> tukey_window(std::size_t length, long double alpha) {
    if (length == 0) {
        return {};
    }
    if (length == 1) {
        return {1.0L};
    }

    // 限制 alpha 在 [0, 1] 范围内
    alpha = std::max(0.0L, std::min(1.0L, alpha));

    std::vector<long double> w(length);

    if (alpha == 0.0L) {
        // 退化为矩形窗
        return rectangular_window(length);
    } else if (alpha == 1.0L) {
        // 退化为汉宁窗
        return hanning_window(length);
    }

    const long double width = alpha * static_cast<long double>(length - 1) / 2.0;

    for (std::size_t i = 0; i < length; ++i) {
        const long double x = static_cast<long double>(i);

        if (x < width) {
            // 上升沿
            w[i] = 0.5 * (1.0L + mymath::cos(kPi * (x / width - 1.0L)));
        } else if (x > static_cast<long double>(length - 1) - width) {
            // 下降沿
            w[i] = 0.5 * (1.0L + mymath::cos(kPi * ((x - static_cast<long double>(length - 1)) / width + 1.0L)));
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

std::vector<long double> window(WindowType type, std::size_t length, long double param) {
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

std::vector<long double> apply_window(const std::vector<long double>& signal,
                                  WindowType window_type,
                                  long double param) {
    if (signal.empty()) {
        return {};
    }

    std::vector<long double> win = window(window_type, signal.size(), param);
    std::vector<long double> result(signal.size());

    for (std::size_t i = 0; i < signal.size(); ++i) {
        result[i] = signal[i] * win[i];
    }

    return result;
}

}  // namespace signal
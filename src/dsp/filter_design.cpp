/**
 * @file filter_design.cpp
 * @brief 滤波器设计实现
 *
 * 实现：
 * - FIR 滤波器设计（窗函数法）
 * - IIR 滤波器设计（巴特沃斯、切比雪夫、椭圆）
 * - 滤波器应用
 * - 频率响应计算
 */

#include "signal_processing.h"
#include "polynomial/polynomial.h"
#include "math/mymath.h"
#include "types/scalar_type.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace signal {

// Use Scalar for internal high-precision real computations
using Scalar = mymath::Scalar;

// 使用 mymath 命名空间中的常量
using mymath::kPi;
static const Scalar kPiScalar = mymath::constants::pi<Scalar>();

// ============================================================================
// 滤波器类型转换
// ============================================================================

std::string filter_type_to_string(FilterType type) {
    switch (type) {
        case FilterType::LowPass:  return "lowpass";
        case FilterType::HighPass: return "highpass";
        case FilterType::BandPass: return "bandpass";
        case FilterType::BandStop: return "bandstop";
        default:                   return "unknown";
    }
}

FilterType string_to_filter_type(const std::string& str) {
    if (str == "lowpass" || str == "lp") {
        return FilterType::LowPass;
    } else if (str == "highpass" || str == "hp") {
        return FilterType::HighPass;
    } else if (str == "bandpass" || str == "bp") {
        return FilterType::BandPass;
    } else if (str == "bandstop" || str == "bs" || str == "notch") {
        return FilterType::BandStop;
    }
    throw std::runtime_error("Unknown filter type: " + str);
}

// ============================================================================
// 理想低通滤波器脉冲响应（sinc 函数）
// ============================================================================

static Scalar sinc_scalar(Scalar x) {
    const Scalar x_s(x);
    if (mymath::abs(x_s) < Scalar(1e-10L)) {
        return 1.0L;
    }
    return (mymath::sin(kPiScalar * x_s) / (kPiScalar * x_s));
}

// ============================================================================
// FIR 滤波器设计（窗函数法）
// ============================================================================

FilterCoefficients design_fir(int order,
                               Scalar cutoff,
                               FilterType type,
                               WindowType window_type) {
    if (order < 1) {
        throw std::runtime_error("Filter order must be at least 1");
    }
    if (cutoff <= 0.0L || cutoff >= 1.0L) {
        throw std::runtime_error("Cutoff frequency must be in (0, 1)");
    }

    // 滤波器长度 = 阶数 + 1
    const std::size_t length = static_cast<std::size_t>(order + 1);
    const int half = order / 2;

    // 归一化截止频率（0 到 0.5 对应 0 到 Nyquist）
    const Scalar fc = cutoff / 2.0;

    // 生成窗函数
    std::vector<Scalar> win = window(window_type, length);

    // 计算理想脉冲响应
    std::vector<Scalar> h(length);

    if (type == FilterType::LowPass) {
        // 低通：h[n] = 2*fc*sinc(2*fc*(n-M/2))
        for (int n = 0; n <= order; ++n) {
            const Scalar x = (n - half);
            h[n] = 2.0 * fc * sinc_scalar(2.0 * fc * x);
        }
    } else if (type == FilterType::HighPass) {
        // 高通：h[n] = delta(n-M/2) - 2*fc*sinc(2*fc*(n-M/2))
        for (int n = 0; n <= order; ++n) {
            const Scalar x = (n - half);
            h[n] = (n == half ? 1.0L : 0.0L) - 2.0 * fc * sinc_scalar(2.0 * fc * x);
        }
    } else {
        throw std::runtime_error("Use design_fir_band for bandpass/bandstop filters");
    }

    // 应用窗函数
    for (std::size_t n = 0; n < length; ++n) {
        h[n] *= win[n];
    }

    return FilterCoefficients{h, {1.0L}};
}

FilterCoefficients design_fir_band(int order,
                                    Scalar cutoff_low,
                                    Scalar cutoff_high,
                                    FilterType type,
                                    WindowType window_type) {
    if (order < 1) {
        throw std::runtime_error("Filter order must be at least 1");
    }
    if (cutoff_low <= 0.0L || cutoff_low >= cutoff_high || cutoff_high >= 1.0L) {
        throw std::runtime_error("Invalid cutoff frequencies");
    }

    const std::size_t length = static_cast<std::size_t>(order + 1);
    const int half = order / 2;

    const Scalar fc1 = cutoff_low / 2.0;
    const Scalar fc2 = cutoff_high / 2.0;

    std::vector<Scalar> win = window(window_type, length);
    std::vector<Scalar> h(length);

    if (type == FilterType::BandPass) {
        // 带通 = 高通(fc1) - 高通(fc2) = 低通(fc2) - 低通(fc1)
        for (int n = 0; n <= order; ++n) {
            const Scalar x = (n - half);
            h[n] = 2.0 * (fc2 - fc1) * sinc_scalar((fc2 - fc1) * x) *
                   (mymath::cos((kPiScalar) * (fc2 + fc1) * x));
        }
    } else if (type == FilterType::BandStop) {
        // 带阻 = 低通(fc1) + 高通(fc2)
        for (int n = 0; n <= order; ++n) {
            const Scalar x = (n - half);
            h[n] = (n == half ? 1.0L : 0.0L) -
                   2.0 * (fc2 - fc1) * sinc_scalar((fc2 - fc1) * x) *
                   (mymath::cos((kPiScalar) * (fc2 + fc1) * x));
        }
    } else {
        throw std::runtime_error("Use design_fir for lowpass/highpass filters");
    }

    for (std::size_t n = 0; n < length; ++n) {
        h[n] *= win[n];
    }

    return FilterCoefficients{h, {1.0L}};
}

// ============================================================================
// IIR 滤波器设计辅助函数
// ============================================================================

// 双线性变换：将模拟极点转换为数字极点
static Complex bilinear_transform(const Complex& s_pole, Scalar fs) {
    return (1.0L + s_pole / (2.0 * fs)) / (1.0L - s_pole / (2.0 * fs));
}

// 将模拟传递函数转换为数字传递函数（双线性变换）
static void bilinear(const std::vector<Complex>& s_zeros,
                     const std::vector<Complex>& s_poles,
                     Scalar fs,
                     std::vector<Scalar>& b,
                     std::vector<Scalar>& a) {
    const std::size_t n_zeros = s_zeros.size();
    const std::size_t n_poles = s_poles.size();

    std::vector<Complex> z_zeros(n_zeros);
    std::vector<Complex> z_poles(n_poles);

    for (std::size_t i = 0; i < n_zeros; ++i) {
        z_zeros[i] = bilinear_transform(s_zeros[i], fs);
    }
    for (std::size_t i = 0; i < n_poles; ++i) {
        z_poles[i] = bilinear_transform(s_poles[i], fs);
    }

    // 从极点/零点展开多项式系数
    auto poly_from_roots = [](const std::vector<Complex>& roots) -> std::vector<Complex> {
        std::vector<Complex> coeffs = {Scalar(1.0L)};
        for (const Complex& r : roots) {
            std::vector<Complex> new_coeffs(coeffs.size() + 1, Scalar(0.0L));
            for (std::size_t i = 0; i < coeffs.size(); ++i) {
                new_coeffs[i] += coeffs[i];
                new_coeffs[i + 1] -= r * coeffs[i];
            }
            coeffs = new_coeffs;
        }
        return coeffs;
    };

    std::vector<Complex> b_complex = poly_from_roots(z_zeros);
    std::vector<Complex> a_complex = poly_from_roots(z_poles);

    // 归一化
    Complex gain = a_complex[0];
    for (Complex& c : b_complex) {
        c /= gain;
    }
    for (Complex& c : a_complex) {
        c /= gain;
    }

    // 提取实部（应该都是实数）
    b.resize(b_complex.size());
    a.resize(a_complex.size());
    for (std::size_t i = 0; i < b.size(); ++i) {
        b[i] = b_complex[i].real();
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        a[i] = a_complex[i].real();
    }
}

// ============================================================================
// 巴特沃斯滤波器设计
// ============================================================================

FilterCoefficients design_butterworth(int order,
                                       Scalar cutoff,
                                       FilterType type) {
    if (order < 1) {
        throw std::runtime_error("Filter order must be at least 1");
    }
    if (cutoff <= 0.0L || cutoff >= 1.0L) {
        throw std::runtime_error("Cutoff frequency must be in (0, 1)");
    }

    // 预扭曲频率（双线性变换）
    const Scalar fs = 1.0L;  // 归一化采样率
    // Use Scalar for tan computation
    const Scalar cutoff_s(cutoff);
    const Scalar tan_arg = (kPiScalar) * cutoff_s / Scalar(2.0L);
    const Scalar wc_scalar = Scalar(2.0L) * fs * mymath::tan(tan_arg);
    const Scalar wc = (wc_scalar);

    // 计算模拟巴特沃斯极点
    std::vector<Complex> s_poles;
    for (int k = 0; k < order; ++k) {
        const Scalar angle_scalar = kPiScalar * (Scalar(2.0L) * Scalar((k)) + Scalar(1.0L)) /
                                    (Scalar(2.0L) * Scalar((order)));
        const Scalar angle = (angle_scalar);
        s_poles.push_back(wc * mymath::polar(Scalar(1.0L), angle + (kPiScalar) / 2.0L));
    }

    std::vector<Complex> s_zeros;  // 巴特沃斯没有有限零点

    std::vector<Scalar> b, a;
    bilinear(s_zeros, s_poles, fs, b, a);

    // 根据滤波器类型调整
    if (type == FilterType::HighPass) {
        // 高通：z -> -z 变换
        for (std::size_t i = 0; i < b.size(); ++i) {
            if (i % 2 == 1) {
                b[i] = -b[i];
            }
        }
        for (std::size_t i = 0; i < a.size(); ++i) {
            if (i % 2 == 1) {
                a[i] = -a[i];
            }
        }
    }

    return FilterCoefficients{b, a};
}

// ============================================================================
// 切比雪夫 I 型滤波器设计
// ============================================================================

FilterCoefficients design_chebyshev1(int order,
                                      Scalar cutoff,
                                      Scalar ripple,
                                      FilterType type) {
    if (order < 1) {
        throw std::runtime_error("Filter order must be at least 1");
    }
    if (cutoff <= 0.0L || cutoff >= 1.0L) {
        throw std::runtime_error("Cutoff frequency must be in (0, 1)");
    }
    if (ripple <= 0.0L) {
        throw std::runtime_error("Ripple must be positive");
    }

    const Scalar fs = 1.0L;
    // Use Scalar for tan computation
    const Scalar cutoff_s(cutoff);
    const Scalar tan_arg = (kPiScalar) * cutoff_s / Scalar(2.0L);
    const Scalar wc_scalar = Scalar(2.0L) * fs * mymath::tan(tan_arg);
    const Scalar wc = (wc_scalar);

    // 从波纹计算 epsilon - use Scalar for precision
    const Scalar ripple_s(ripple);
    const Scalar epsilon = mymath::sqrt(mymath::pow(Scalar(10.0L), ripple_s / Scalar(10.0L)) - Scalar(1.0L));

    // 计算切比雪夫极点
    const Scalar mu = mymath::asinh(Scalar(1.0L) / epsilon) / Scalar((order));

    std::vector<Complex> s_poles;
    for (int k = 0; k < order; ++k) {
        const Scalar theta = kPiScalar * (Scalar(2.0L) * Scalar((k)) + Scalar(1.0L)) /
                              (Scalar(2.0L) * Scalar((order)));
        const Scalar sigma_val = -mymath::sinh(mu) * mymath::sin(theta);
        const Scalar omega_val = mymath::cosh(mu) * mymath::cos(theta);
        s_poles.push_back(wc * Complex((sigma_val), (omega_val)));
    }

    std::vector<Complex> s_zeros;
    std::vector<Scalar> b, a;
    bilinear(s_zeros, s_poles, fs, b, a);

    // 直流增益调整
    Scalar dc_gain = 0.0L;
    for (Scalar coeff : b) {
        dc_gain += coeff;
    }
    Scalar dc_denom = 0.0L;
    for (Scalar coeff : a) {
        dc_denom += coeff;
    }
    dc_gain /= dc_denom;

    if (type == FilterType::LowPass) {
        // 归一化直流增益为 1
        for (Scalar& coeff : b) {
            coeff /= dc_gain;
        }
    } else if (type == FilterType::HighPass) {
        // 高通变换
        for (std::size_t i = 0; i < b.size(); ++i) {
            if (i % 2 == 1) {
                b[i] = -b[i];
            }
        }
        for (std::size_t i = 0; i < a.size(); ++i) {
            if (i % 2 == 1) {
                a[i] = -a[i];
            }
        }
    }

    return FilterCoefficients{b, a};
}

// ============================================================================
// 椭圆滤波器设计（简化实现）
// ============================================================================

FilterCoefficients design_elliptic(int order,
                                    Scalar cutoff,
                                    Scalar ripple,
                                    Scalar stopband_atten,
                                    FilterType type) {
    (void)order;
    (void)cutoff;
    (void)ripple;
    (void)stopband_atten;
    (void)type;
    throw std::runtime_error("Elliptic filter design is not fully implemented yet (requires Jacobi elliptic functions).");
}

// ============================================================================
// 二阶节 (SOS) 滤波器实现
// ============================================================================

std::vector<Scalar> sosfilter(const std::vector<SOS>& sos, const std::vector<Scalar>& signal) {
    if (signal.empty()) return {};
    if (sos.empty()) return signal;

    std::vector<Scalar> output = signal;

    // 用于存储每个节的状态 (Direct Form II) - use Scalar for precision
    struct State { Scalar w1 = Scalar(0.0L), w2 = Scalar(0.0L); };
    std::vector<State> states(sos.size());

    for (Scalar& x : output) {
        Scalar val = Scalar(x);
        for (std::size_t i = 0; i < sos.size(); ++i) {
            const SOS& s = sos[i];
            State& st = states[i];

            // w[n] = x[n] - a1*w[n-1] - a2*w[n-2]
            Scalar w = val - Scalar(s.a1) * st.w1 - Scalar(s.a2) * st.w2;
            // y[n] = b0*w[n] + b1*w[n-1] + b2*w[n-2]
            val = Scalar(s.b0) * w + Scalar(s.b1) * st.w1 + Scalar(s.b2) * st.w2;

            st.w2 = st.w1;
            st.w1 = w;
        }
        x = (val);
    }

    return output;
}

std::vector<SOS> tf2sos(const std::vector<Scalar>& b, const std::vector<Scalar>& a) {
    if (a.empty() || b.empty()) return {};

    const Scalar a0 = a[0];
    const Scalar b0 = b[0];
    if (mymath::is_near_zero(a0, 1e-12L)) {
        throw std::runtime_error("tf2sos: leading denominator coefficient a[0] must not be zero");
    }

    if (a.size() <= 3 && b.size() <= 3) {
        SOS s;
        s.b0 = b[0] / a0;
        s.b1 = (b.size() > 1 ? b[1] / a0 : 0.0L);
        s.b2 = (b.size() > 2 ? b[2] / a0 : 0.0L);
        s.a1 = (a.size() > 1 ? a[1] / a0 : 0.0L);
        s.a2 = (a.size() > 2 ? a[2] / a0 : 0.0L);
        return {s};
    }

    std::vector<Scalar> a_poly(a.rbegin(), a.rend());
    std::vector<Scalar> b_poly(b.rbegin(), b.rend());

    std::vector<mymath::complex<Scalar>> poles;
    if (a.size() > 1) {
        poles = polynomial_complex_roots(a_poly);
    }
    std::vector<mymath::complex<Scalar>> zeros;
    if (b.size() > 1) {
        zeros = polynomial_complex_roots(b_poly);
    }

    struct SectionPair {
        Scalar c1 = 0;
        Scalar c2 = 0;
    };

    auto pair_roots = [](std::vector<mymath::complex<Scalar>>& roots) -> std::vector<SectionPair> {
        std::vector<SectionPair> sections;
        std::vector<bool> used(roots.size(), false);

        // 先配对共轭复根
        for (std::size_t i = 0; i < roots.size(); ++i) {
            if (used[i]) continue;
            if (mymath::abs(roots[i].imag()) > Scalar(1e-7L)) {
                int best_j = -1;
                Scalar best_diff = Scalar(1e-4L);
                for (std::size_t j = i + 1; j < roots.size(); ++j) {
                    if (used[j]) continue;
                    Scalar diff = mymath::abs(roots[i].real() - roots[j].real()) +
                                  mymath::abs(roots[i].imag() + roots[j].imag());
                    if (diff < best_diff) {
                        best_diff = diff;
                        best_j = static_cast<int>(j);
                    }
                }
                if (best_j >= 0) {
                    used[i] = true;
                    used[best_j] = true;
                    SectionPair sp;
                    sp.c1 = -Scalar(2.0L) * roots[i].real();
                    sp.c2 = roots[i].real() * roots[i].real() + roots[i].imag() * roots[i].imag();
                    sections.push_back(sp);
                }
            }
        }

        // 再配对实根
        std::vector<Scalar> real_roots;
        for (std::size_t i = 0; i < roots.size(); ++i) {
            if (!used[i]) {
                real_roots.push_back(roots[i].real());
            }
        }

        for (std::size_t i = 0; i + 1 < real_roots.size(); i += 2) {
            SectionPair sp;
            sp.c1 = -(real_roots[i] + real_roots[i + 1]);
            sp.c2 = real_roots[i] * real_roots[i + 1];
            sections.push_back(sp);
        }
        if (real_roots.size() % 2 == 1) {
            SectionPair sp;
            sp.c1 = -real_roots.back();
            sp.c2 = Scalar(0.0L);
            sections.push_back(sp);
        }
        return sections;
    };

    std::vector<SectionPair> a_sections = pair_roots(poles);
    std::vector<SectionPair> b_sections = pair_roots(zeros);

    std::size_t num_sections = std::max(a_sections.size(), b_sections.size());
    if (num_sections == 0) num_sections = 1;

    std::vector<SOS> result(num_sections);
    const Scalar overall_gain = b0 / a0;

    for (std::size_t i = 0; i < num_sections; ++i) {
        result[i].b0 = (i == 0) ? overall_gain : Scalar(1.0L);
        if (i < b_sections.size()) {
            result[i].b1 = (i == 0) ? overall_gain * b_sections[i].c1 : b_sections[i].c1;
            result[i].b2 = (i == 0) ? overall_gain * b_sections[i].c2 : b_sections[i].c2;
        } else {
            result[i].b1 = Scalar(0.0L);
            result[i].b2 = Scalar(0.0L);
        }

        if (i < a_sections.size()) {
            result[i].a1 = a_sections[i].c1;
            result[i].a2 = a_sections[i].c2;
        } else {
            result[i].a1 = Scalar(0.0L);
            result[i].a2 = Scalar(0.0L);
        }
    }
    return result;
}

// ============================================================================
// 滤波器应用
// ============================================================================

std::vector<Scalar> filter(const std::vector<Scalar>& b,
                            const std::vector<Scalar>& a,
                            const std::vector<Scalar>& signal) {
    if (signal.empty()) return {};

    // 优化：对于高阶 IIR 滤波器（a.size > 3），尝试使用 SOS
    // 但由于目前 tf2sos 还没实现全自动配对，先保持直接型，
    // 但增加了数值保护和归一化。

    const std::size_t n = signal.size();
    const std::size_t nb = b.size();
    const std::size_t na = a.size();

    std::vector<Scalar> output(n, 0.0L);
    const Scalar a0 = Scalar(a[0]);

    // 使用 Scalar 进行高精度计算
    std::vector<Scalar> b_scalar(nb);
    std::vector<Scalar> a_scalar(na);
    for (std::size_t i = 0; i < nb; ++i) b_scalar[i] = Scalar(b[i]);
    for (std::size_t i = 0; i < na; ++i) a_scalar[i] = Scalar(a[i]);

    // 使用状态空间或缓冲区以提高性能并减少内存访问
    std::vector<Scalar> w(std::max(nb, na), Scalar(0.0L));

    for (std::size_t i = 0; i < n; ++i) {
        Scalar val = Scalar(signal[i]);

        // 采用 Direct Form II
        Scalar wn = val;
        for (std::size_t j = 1; j < na; ++j) {
            wn -= (a_scalar[j] / a0) * w[j - 1];
        }

        Scalar yn = (b_scalar[0] / a0) * wn;
        for (std::size_t j = 1; j < nb; ++j) {
            yn += (b_scalar[j] / a0) * w[j - 1];
        }

        // 更新状态
        for (std::size_t j = std::max(nb, na) - 1; j > 0; --j) {
            w[j] = w[j - 1];
        }
        w[0] = wn;

        output[i] = (yn);
    }

    return output;
}

// ============================================================================
// 零相位滤波
// ============================================================================

std::vector<Scalar> filtfilt(const std::vector<Scalar>& b,
                              const std::vector<Scalar>& a,
                              const std::vector<Scalar>& signal) {
    if (signal.empty()) {
        return {};
    }

    // 计算填充长度
    std::size_t padlen = 3 * std::max(b.size(), a.size());
    if (signal.size() <= padlen) {
        padlen = signal.size() - 1;
    }

    std::vector<Scalar> padded_signal;
    padded_signal.reserve(signal.size() + 2 * padlen);

    // 奇对称延拓 (Odd Extension) 左侧
    for (std::size_t i = padlen; i > 0; --i) {
        padded_signal.push_back(2.0 * signal[0] - signal[i]);
    }

    // 原始信号
    padded_signal.insert(padded_signal.end(), signal.begin(), signal.end());

    // 奇对称延拓 (Odd Extension) 右侧
    for (std::size_t i = 1; i <= padlen; ++i) {
        padded_signal.push_back(2.0 * signal.back() - signal[signal.size() - 1 - i]);
    }

    // 前向滤波
    std::vector<Scalar> forward = filter(b, a, padded_signal);

    // 反向
    std::vector<Scalar> reversed(forward.rbegin(), forward.rend());

    // 反向滤波
    std::vector<Scalar> backward = filter(b, a, reversed);

    // 再反向恢复原顺序，并截取有效部分
    std::vector<Scalar> result(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        result[i] = backward[backward.size() - 1 - (padlen + i)];
    }

    return result;
}

// ============================================================================
// 频率响应
// ============================================================================

std::vector<Complex> freqz(const std::vector<Scalar>& b,
                            const std::vector<Scalar>& a,
                            std::size_t n) {
    if (b.empty()) {
        throw std::runtime_error("Numerator coefficients cannot be empty");
    }
    if (a.empty()) {
        throw std::runtime_error("Denominator coefficients cannot be empty");
    }

    std::vector<Complex> response(n);

    // Convert coefficients to Scalar for high precision computation
    std::vector<Scalar> b_scalar(b.size()), a_scalar(a.size());
    for (std::size_t i = 0; i < b.size(); ++i) b_scalar[i] = Scalar(b[i]);
    for (std::size_t i = 0; i < a.size(); ++i) a_scalar[i] = Scalar(a[i]);

    for (std::size_t k = 0; k < n; ++k) {
        const Scalar angle = kPiScalar * Scalar(static_cast<long long>(k)) / Scalar(static_cast<long long>(n));
        const Scalar omega = (angle);
        const Complex z = Complex(mymath::cos(omega), mymath::sin(omega));

        // 计算分子
        Complex num = Scalar(0.0L);
        Complex z_pow = Scalar(1.0L);
        for (Scalar coeff : b) {
            num += coeff * z_pow;
            z_pow *= z;
        }

        // 计算分母
        Complex den = Scalar(0.0L);
        z_pow = Scalar(1.0L);
        for (Scalar coeff : a) {
            den += coeff * z_pow;
            z_pow *= z;
        }

        response[k] = num / den;
    }

    return response;
}

// ============================================================================
// 群延迟
// ============================================================================

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
#endif
std::vector<Scalar> grpdelay(const std::vector<Scalar>& b,
                              const std::vector<Scalar>& a,
                              std::size_t n) {
    std::vector<Complex> h = freqz(b, a, n);
    std::vector<Scalar> gd(n);

    if (n < 2) return std::vector<Scalar>(n, 0.0L);

    // 提取相位并展开 (unwrap)
    std::vector<Scalar> phase(n);
    for (std::size_t k = 0; k < n; ++k) {
        phase[k] = mymath::arg(h[k]);
    }

    for (std::size_t k = 1; k < n; ++k) {
        Scalar dp = phase[k] - phase[k - 1];
        while (dp > kPi) { phase[k] -= 2.0 * kPi; dp -= 2.0 * kPi; }
        while (dp < -kPi) { phase[k] += 2.0 * kPi; dp += 2.0 * kPi; }
    }

    // 数值微分计算群延迟
    const Scalar dw = kPi / Scalar(static_cast<long long>(n));

    for (std::size_t k = 0; k < n; ++k) {
        if (k == 0) {
            // 前向差分
            const Scalar dphi = phase[1] - phase[0];
            gd[k] = -dphi / dw;
        } else if (k == n - 1) {
            // 后向差分
            const Scalar dphi = phase[k] - phase[k - 1];
            gd[k] = -dphi / dw;
        } else {
            // 中心差分
            const Scalar dphi = phase[k + 1] - phase[k - 1];
            gd[k] = -dphi / (2.0 * dw);
        }
    }

    return gd;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

}  // namespace signal

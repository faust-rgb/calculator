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
#include "math/mymath.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace signal {

// 数学常量
constexpr double kPi = 3.14159265358979323846;

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

static double sinc(double x) {
    if (mymath::abs(x) < 1e-10) {
        return 1.0;
    }
    return mymath::sin(kPi * x) / (kPi * x);
}

// ============================================================================
// FIR 滤波器设计（窗函数法）
// ============================================================================

FilterCoefficients design_fir(int order,
                               double cutoff,
                               FilterType type,
                               WindowType window_type) {
    if (order < 1) {
        throw std::runtime_error("Filter order must be at least 1");
    }
    if (cutoff <= 0.0 || cutoff >= 1.0) {
        throw std::runtime_error("Cutoff frequency must be in (0, 1)");
    }

    // 滤波器长度 = 阶数 + 1
    const std::size_t length = static_cast<std::size_t>(order + 1);
    const int half = order / 2;

    // 归一化截止频率（0 到 0.5 对应 0 到 Nyquist）
    const double fc = cutoff / 2.0;

    // 生成窗函数
    std::vector<double> win = window(window_type, length);

    // 计算理想脉冲响应
    std::vector<double> h(length);

    if (type == FilterType::LowPass) {
        // 低通：h[n] = 2*fc*sinc(2*fc*(n-M/2))
        for (int n = 0; n <= order; ++n) {
            const double x = static_cast<double>(n - half);
            h[n] = 2.0 * fc * sinc(2.0 * fc * x);
        }
    } else if (type == FilterType::HighPass) {
        // 高通：h[n] = delta(n-M/2) - 2*fc*sinc(2*fc*(n-M/2))
        for (int n = 0; n <= order; ++n) {
            const double x = static_cast<double>(n - half);
            h[n] = (n == half ? 1.0 : 0.0) - 2.0 * fc * sinc(2.0 * fc * x);
        }
    } else {
        throw std::runtime_error("Use design_fir_band for bandpass/bandstop filters");
    }

    // 应用窗函数
    for (std::size_t n = 0; n < length; ++n) {
        h[n] *= win[n];
    }

    return FilterCoefficients{h, {1.0}};
}

FilterCoefficients design_fir_band(int order,
                                    double cutoff_low,
                                    double cutoff_high,
                                    FilterType type,
                                    WindowType window_type) {
    if (order < 1) {
        throw std::runtime_error("Filter order must be at least 1");
    }
    if (cutoff_low <= 0.0 || cutoff_low >= cutoff_high || cutoff_high >= 1.0) {
        throw std::runtime_error("Invalid cutoff frequencies");
    }

    const std::size_t length = static_cast<std::size_t>(order + 1);
    const int half = order / 2;

    const double fc1 = cutoff_low / 2.0;
    const double fc2 = cutoff_high / 2.0;

    std::vector<double> win = window(window_type, length);
    std::vector<double> h(length);

    if (type == FilterType::BandPass) {
        // 带通 = 高通(fc1) - 高通(fc2) = 低通(fc2) - 低通(fc1)
        for (int n = 0; n <= order; ++n) {
            const double x = static_cast<double>(n - half);
            h[n] = 2.0 * (fc2 - fc1) * sinc((fc2 - fc1) * x) *
                   mymath::cos(kPi * (fc2 + fc1) * x);
        }
    } else if (type == FilterType::BandStop) {
        // 带阻 = 低通(fc1) + 高通(fc2)
        for (int n = 0; n <= order; ++n) {
            const double x = static_cast<double>(n - half);
            h[n] = (n == half ? 1.0 : 0.0) -
                   2.0 * (fc2 - fc1) * sinc((fc2 - fc1) * x) *
                   mymath::cos(kPi * (fc2 + fc1) * x);
        }
    } else {
        throw std::runtime_error("Use design_fir for lowpass/highpass filters");
    }

    for (std::size_t n = 0; n < length; ++n) {
        h[n] *= win[n];
    }

    return FilterCoefficients{h, {1.0}};
}

// ============================================================================
// IIR 滤波器设计辅助函数
// ============================================================================

// 双线性变换：将模拟极点转换为数字极点
static Complex bilinear_transform(const Complex& s_pole, double fs) {
    return (1.0 + s_pole / (2.0 * fs)) / (1.0 - s_pole / (2.0 * fs));
}

// 将模拟传递函数转换为数字传递函数（双线性变换）
static void bilinear(const std::vector<Complex>& s_zeros,
                     const std::vector<Complex>& s_poles,
                     double fs,
                     std::vector<double>& b,
                     std::vector<double>& a) {
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
        std::vector<Complex> coeffs = {1.0};
        for (const Complex& r : roots) {
            std::vector<Complex> new_coeffs(coeffs.size() + 1, 0.0);
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
                                       double cutoff,
                                       FilterType type) {
    if (order < 1) {
        throw std::runtime_error("Filter order must be at least 1");
    }
    if (cutoff <= 0.0 || cutoff >= 1.0) {
        throw std::runtime_error("Cutoff frequency must be in (0, 1)");
    }

    // 预扭曲频率（双线性变换）
    const double fs = 1.0;  // 归一化采样率
    const double wc = 2.0 * fs * mymath::tan(kPi * cutoff / 2.0);

    // 计算模拟巴特沃斯极点
    std::vector<Complex> s_poles;
    for (int k = 0; k < order; ++k) {
        const double angle = kPi * (2.0 * static_cast<double>(k) + 1.0) / (2.0 * static_cast<double>(order));
        s_poles.push_back(wc * mymath::polar(1.0, angle + kPi / 2.0));
    }

    std::vector<Complex> s_zeros;  // 巴特沃斯没有有限零点

    std::vector<double> b, a;
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
                                      double cutoff,
                                      double ripple,
                                      FilterType type) {
    if (order < 1) {
        throw std::runtime_error("Filter order must be at least 1");
    }
    if (cutoff <= 0.0 || cutoff >= 1.0) {
        throw std::runtime_error("Cutoff frequency must be in (0, 1)");
    }
    if (ripple <= 0.0) {
        throw std::runtime_error("Ripple must be positive");
    }

    const double fs = 1.0;
    const double wc = 2.0 * fs * mymath::tan(kPi * cutoff / 2.0);

    // 从波纹计算 epsilon
    const double epsilon = mymath::sqrt(mymath::pow(10.0, ripple / 10.0) - 1.0);

    // 计算切比雪夫极点
    const double mu = mymath::asinh(1.0 / epsilon) / static_cast<double>(order);

    std::vector<Complex> s_poles;
    for (int k = 0; k < order; ++k) {
        const double theta = kPi * (2.0 * static_cast<double>(k) + 1.0) / (2.0 * static_cast<double>(order));
        const double sigma = -mymath::sinh(mu) * mymath::sin(theta);
        const double omega = mymath::cosh(mu) * mymath::cos(theta);
        s_poles.push_back(wc * Complex(sigma, omega));
    }

    std::vector<Complex> s_zeros;
    std::vector<double> b, a;
    bilinear(s_zeros, s_poles, fs, b, a);

    // 直流增益调整
    double dc_gain = 0.0;
    for (double coeff : b) {
        dc_gain += coeff;
    }
    double dc_denom = 0.0;
    for (double coeff : a) {
        dc_denom += coeff;
    }
    dc_gain /= dc_denom;

    if (type == FilterType::LowPass) {
        // 归一化直流增益为 1
        for (double& coeff : b) {
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
                                    double cutoff,
                                    double ripple,
                                    double,
                                    FilterType type) {
    // 椭圆滤波器设计较复杂，这里使用简化的近似方法
    // 实际实现需要计算椭圆积分和雅可比椭圆函数

    if (order < 1) {
        throw std::runtime_error("Filter order must be at least 1");
    }
    if (cutoff <= 0.0 || cutoff >= 1.0) {
        throw std::runtime_error("Cutoff frequency must be in (0, 1)");
    }

    // 简化：使用切比雪夫作为近似
    // 实际椭圆滤波器需要更复杂的计算
    return design_chebyshev1(order, cutoff, ripple, type);
}

// ============================================================================
// 二阶节 (SOS) 滤波器实现
// ============================================================================

std::vector<double> sosfilter(const std::vector<SOS>& sos, const std::vector<double>& signal) {
    if (signal.empty()) return {};
    if (sos.empty()) return signal;

    std::vector<double> output = signal;
    
    // 用于存储每个节的状态 (Direct Form II)
    struct State { double w1 = 0.0, w2 = 0.0; };
    std::vector<State> states(sos.size());

    for (double& x : output) {
        double val = x;
        for (std::size_t i = 0; i < sos.size(); ++i) {
            const SOS& s = sos[i];
            State& st = states[i];
            
            // w[n] = x[n] - a1*w[n-1] - a2*w[n-2]
            double w = val - s.a1 * st.w1 - s.a2 * st.w2;
            // y[n] = b0*w[n] + b1*w[n-1] + b2*w[n-2]
            val = s.b0 * w + s.b1 * st.w1 + s.b2 * st.w2;
            
            st.w2 = st.w1;
            st.w1 = w;
        }
        x = val;
    }

    return output;
}

std::vector<SOS> tf2sos(const std::vector<double>& b, const std::vector<double>& a) {
    // 这是一个简化的实现。
    // 对于真正的 tf2sos，我们需要通过多项式求根找到极点和零点，然后配对。
    // 在本仓库中，由于 IIR 设计函数已知极点，建议在设计函数中直接生成 SOS。
    // 这里暂时返回一个单节（如果是低阶）或抛出异常。
    
    if (a.size() <= 3 && b.size() <= 3) {
        SOS s;
        const double a0 = a[0];
        s.b0 = b[0] / a0;
        s.b1 = (b.size() > 1 ? b[1] / a0 : 0.0);
        s.b2 = (b.size() > 2 ? b[2] / a0 : 0.0);
        s.a1 = (a.size() > 1 ? a[1] / a0 : 0.0);
        s.a2 = (a.size() > 2 ? a[2] / a0 : 0.0);
        return {s};
    }
    
    throw std::runtime_error("tf2sos for high-order filters requires polynomial root finding (not yet fully implemented)");
}

// ============================================================================
// 滤波器应用
// ============================================================================

std::vector<double> filter(const std::vector<double>& b,
                            const std::vector<double>& a,
                            const std::vector<double>& signal) {
    if (signal.empty()) return {};
    
    // 优化：对于高阶 IIR 滤波器（a.size > 3），尝试使用 SOS
    // 但由于目前 tf2sos 还没实现全自动配对，先保持直接型，
    // 但增加了数值保护和归一化。

    const std::size_t n = signal.size();
    const std::size_t nb = b.size();
    const std::size_t na = a.size();
    
    std::vector<double> output(n, 0.0);
    const double a0 = a[0];
    
    // 使用状态空间或缓冲区以提高性能并减少内存访问
    std::vector<double> w(std::max(nb, na), 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        double val = signal[i];
        
        // 采用 Direct Form II
        double wn = val;
        for (std::size_t j = 1; j < na; ++j) {
            wn -= (a[j] / a0) * w[j - 1];
        }
        
        double yn = (b[0] / a0) * wn;
        for (std::size_t j = 1; j < nb; ++j) {
            yn += (b[j] / a0) * w[j - 1];
        }
        
        // 更新状态
        for (std::size_t j = std::max(nb, na) - 1; j > 0; --j) {
            w[j] = w[j - 1];
        }
        w[0] = wn;
        
        output[i] = yn;
    }

    return output;
}

// ============================================================================
// 零相位滤波
// ============================================================================

std::vector<double> filtfilt(const std::vector<double>& b,
                              const std::vector<double>& a,
                              const std::vector<double>& signal) {
    if (signal.empty()) {
        return {};
    }

    // 前向滤波
    std::vector<double> forward = filter(b, a, signal);

    // 反向
    std::vector<double> reversed(forward.rbegin(), forward.rend());

    // 反向滤波
    std::vector<double> backward = filter(b, a, reversed);

    // 再反向恢复原顺序
    return std::vector<double>(backward.rbegin(), backward.rend());
}

// ============================================================================
// 频率响应
// ============================================================================

std::vector<Complex> freqz(const std::vector<double>& b,
                            const std::vector<double>& a,
                            std::size_t n) {
    if (b.empty()) {
        throw std::runtime_error("Numerator coefficients cannot be empty");
    }
    if (a.empty()) {
        throw std::runtime_error("Denominator coefficients cannot be empty");
    }

    std::vector<Complex> response(n);

    for (std::size_t k = 0; k < n; ++k) {
        const double omega = kPi * static_cast<double>(k) / static_cast<double>(n);
        const Complex z = mymath::polar(1.0, omega);

        // 计算分子
        Complex num = 0.0;
        Complex z_pow = 1.0;
        for (double coeff : b) {
            num += coeff * z_pow;
            z_pow *= z;
        }

        // 计算分母
        Complex den = 0.0;
        z_pow = 1.0;
        for (double coeff : a) {
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

std::vector<double> grpdelay(const std::vector<double>& b,
                              const std::vector<double>& a,
                              std::size_t n) {
    std::vector<Complex> h = freqz(b, a, n);
    std::vector<double> gd(n);

    // 数值微分计算群延迟
    const double dw = kPi / static_cast<double>(n);

    for (std::size_t k = 0; k < n; ++k) {
        if (k == 0) {
            // 前向差分
            const double dphi = mymath::arg(h[1]) - mymath::arg(h[0]);
            gd[k] = -dphi / dw;
        } else if (k == n - 1) {
            // 后向差分
            const double dphi = mymath::arg(h[k]) - mymath::arg(h[k - 1]);
            gd[k] = -dphi / dw;
        } else {
            // 中心差分
            const double dphi = mymath::arg(h[k + 1]) - mymath::arg(h[k - 1]);
            gd[k] = -dphi / (2.0 * dw);
        }
    }

    return gd;
}

}  // namespace signal

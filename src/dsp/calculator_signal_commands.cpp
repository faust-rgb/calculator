/**
 * @file calculator_signal_commands.cpp
 * @brief 信号处理命令实现
 */

#include "calculator_signal_commands.h"
#include "parser/unified_expression_parser.h"
#include "execution/variable_resolver.h"
#include "parser/parser_utils.h"
#include "core/string_utils.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace signal_cmds {

// 数学常量
constexpr double kPi = 3.14159265358979323846;

// ============================================================================
// 辅助函数
// ============================================================================

// 使用统一解析器解析向量
static std::vector<double> parse_vector(const SignalContext& ctx, const std::string& str) {
    if (utils::trim_view(str).empty()) return {};

    UnifiedExpressionParser parser(VariableResolver(ctx.variables, nullptr), ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);
    StoredValue val = parser.evaluate_stored(str);

    std::vector<double> result;
    if (val.is_matrix) {
        for (double d : val.matrix.data) result.push_back(d);
    } else if (val.is_list && val.list_value) {
        for (const auto& item : *val.list_value) {
            result.push_back(item.exact ? rational_to_double(item.rational) : item.decimal);
        }
    } else {
        result.push_back(val.exact ? rational_to_double(val.rational) : val.decimal);
    }
    return result;
}

// 格式化向量输出
static std::string format_vector(const std::vector<double>& vec, int precision = 6) {
    std::ostringstream oss;
    oss << std::setprecision(precision) << std::fixed;
    oss << "[";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << vec[i];
    }
    oss << "]";
    return oss.str();
}

// 格式化复数向量输出
static std::string format_complex_vector(const std::vector<signal::Complex>& vec, int precision = 6) {
    std::ostringstream oss;
    oss << std::setprecision(precision) << std::fixed;
    oss << "[";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) oss << ", ";
        double re = vec[i].real();
        double im = vec[i].imag();
        oss << re;
        if (im >= 0) oss << "+";
        oss << im << "j";
    }
    oss << "]";
    return oss.str();
}

// 解析参数列表（使用 parser_utils）
static std::vector<std::string> parse_args(const std::string& args) {
    auto result = parser_utils::split_top_level(args, ',');
    for (auto& arg : result) {
        arg = utils::trim_copy(arg);
    }
    return result;
}

// ============================================================================
// FFT 命令
// ============================================================================

bool handle_fft_command(const SignalContext& ctx,
                        const std::string& args,
                        std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.empty()) {
        *output = "Error: fft requires a signal argument";
        return false;
    }

    std::vector<double> signal = parse_vector(ctx, arg_list[0]);
    if (signal.empty()) {
        *output = "Error: empty signal";
        return false;
    }

    // 转换为复数
    std::vector<signal::Complex> complex_signal(signal.size());
    for (std::size_t i = 0; i < signal.size(); ++i) {
        complex_signal[i] = signal::Complex(signal[i], 0.0);
    }

    // 计算 FFT
    std::vector<signal::Complex> spectrum = signal::fft(complex_signal);

    *output = format_complex_vector(spectrum);
    return true;
}

bool handle_ifft_command(const SignalContext& ctx,
                         const std::string& args,
                         std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.empty()) {
        *output = "Error: ifft requires a spectrum argument";
        return false;
    }

    // 解析复数频谱（简化处理：目前暂只支持实数输入向量）
    std::vector<double> spectrum_real = parse_vector(ctx, arg_list[0]);
    if (spectrum_real.empty()) {
        *output = "Error: empty spectrum";
        return false;
    }

    std::vector<signal::Complex> spectrum(spectrum_real.size());
    for (std::size_t i = 0; i < spectrum_real.size(); ++i) {
        spectrum[i] = signal::Complex(spectrum_real[i], 0.0);
    }

    std::vector<signal::Complex> time_signal = signal::ifft(spectrum);

    // 提取实部
    std::vector<double> result(time_signal.size());
    for (std::size_t i = 0; i < time_signal.size(); ++i) {
        result[i] = time_signal[i].real();
    }

    *output = format_vector(result);
    return true;
}

bool handle_rfft_command(const SignalContext& ctx,
                         const std::string& args,
                         std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.empty()) {
        *output = "Error: rfft requires a signal argument";
        return false;
    }

    std::vector<double> signal = parse_vector(ctx, arg_list[0]);
    if (signal.empty()) {
        *output = "Error: empty signal";
        return false;
    }

    std::vector<signal::Complex> spectrum = signal::rfft(signal);

    *output = format_complex_vector(spectrum);
    return true;
}

// ============================================================================
// 卷积与相关命令
// ============================================================================

bool handle_conv_command(const SignalContext& ctx,
                         const std::string& args,
                         std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.size() < 2) {
        *output = "Error: conv requires two signal arguments";
        return false;
    }

    std::vector<double> signal1 = parse_vector(ctx, arg_list[0]);
    std::vector<double> signal2 = parse_vector(ctx, arg_list[1]);

    if (signal1.empty() || signal2.empty()) {
        *output = "Error: empty signal";
        return false;
    }

    std::vector<double> result = signal::convolve(signal1, signal2);

    *output = format_vector(result);
    return true;
}

bool handle_cconv_command(const SignalContext& ctx,
                          const std::string& args,
                          std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.size() < 2) {
        *output = "Error: cconv requires two signal arguments";
        return false;
    }

    std::vector<double> signal1 = parse_vector(ctx, arg_list[0]);
    std::vector<double> signal2 = parse_vector(ctx, arg_list[1]);

    if (signal1.empty() || signal2.empty()) {
        *output = "Error: empty signal";
        return false;
    }

    std::size_t n = 0;
    if (arg_list.size() >= 3) {
        try {
            UnifiedExpressionParser parser(VariableResolver(ctx.variables, nullptr), ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);
            n = static_cast<std::size_t>(parser.evaluate(arg_list[2]));
        } catch (...) {
            // 使用默认值
        }
    }

    std::vector<double> result = signal::circular_convolve(signal1, signal2, n);

    *output = format_vector(result);
    return true;
}

bool handle_xcorr_command(const SignalContext& ctx,
                          const std::string& args,
                          std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.size() < 2) {
        *output = "Error: xcorr requires two signal arguments";
        return false;
    }

    std::vector<double> signal1 = parse_vector(ctx, arg_list[0]);
    std::vector<double> signal2 = parse_vector(ctx, arg_list[1]);

    if (signal1.empty() || signal2.empty()) {
        *output = "Error: empty signal";
        return false;
    }

    std::vector<double> result = signal::xcorr(signal1, signal2);

    *output = format_vector(result);
    return true;
}

bool handle_autocorr_command(const SignalContext& ctx,
                             const std::string& args,
                             std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.empty()) {
        *output = "Error: autocorr requires a signal argument";
        return false;
    }

    std::vector<double> signal = parse_vector(ctx, arg_list[0]);
    if (signal.empty()) {
        *output = "Error: empty signal";
        return false;
    }

    std::vector<double> result = signal::autocorr(signal);

    *output = format_vector(result);
    return true;
}

// ============================================================================
// 窗函数命令
// ============================================================================

bool handle_window_command(const SignalContext& ctx,
                           const std::string& args,
                           std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.size() < 2) {
        *output = "Error: window requires type and length arguments";
        return false;
    }

    UnifiedExpressionParser parser(VariableResolver(ctx.variables, nullptr), ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);
    
    std::string type_str = arg_list[0];
    if (type_str.front() == '"' || type_str.front() == '\'') type_str = type_str.substr(1, type_str.size() - 2);

    std::size_t length = 0;
    try {
        length = static_cast<std::size_t>(parser.evaluate(arg_list[1]));
    } catch (...) {
        *output = "Error: invalid length";
        return false;
    }

    double param = 0.0;
    if (arg_list.size() >= 3) {
        try {
            param = parser.evaluate(arg_list[2]);
        } catch (...) {
            // 使用默认值
        }
    }

    try {
        signal::WindowType type = signal::string_to_window_type(type_str);
        std::vector<double> result = signal::window(type, length, param);
        *output = format_vector(result);
        return true;
    } catch (const std::exception& e) {
        *output = std::string("Error: ") + e.what();
        return false;
    }
}

// ============================================================================
// 滤波器命令
// ============================================================================

bool handle_filter_command(const SignalContext& ctx,
                           const std::string& args,
                           std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.size() < 3) {
        *output = "Error: filter requires b, a, and signal arguments";
        return false;
    }

    std::vector<double> b = parse_vector(ctx, arg_list[0]);
    std::vector<double> a = parse_vector(ctx, arg_list[1]);
    std::vector<double> signal = parse_vector(ctx, arg_list[2]);

    if (b.empty() || a.empty() || signal.empty()) {
        *output = "Error: empty coefficients or signal";
        return false;
    }

    try {
        std::vector<double> result = signal::filter(b, a, signal);
        *output = format_vector(result);
        return true;
    } catch (const std::exception& e) {
        *output = std::string("Error: ") + e.what();
        return false;
    }
}

bool handle_fir_design_command(const SignalContext& ctx,
                               const std::string& args,
                               std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.size() < 3) {
        *output = "Error: fir_design requires order, cutoff, and type arguments";
        return false;
    }

    UnifiedExpressionParser parser(VariableResolver(ctx.variables, nullptr), ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);

    int order = 0;
    double cutoff = 0.0;

    try {
        order = static_cast<int>(parser.evaluate(arg_list[0]));
        cutoff = parser.evaluate(arg_list[1]);
    } catch (...) {
        *output = "Error: invalid order or cutoff";
        return false;
    }

    std::string type_str = arg_list[2];
    if (type_str.front() == '"' || type_str.front() == '\'') type_str = type_str.substr(1, type_str.size() - 2);

    signal::FilterType type;
    try {
        type = signal::string_to_filter_type(type_str);
    } catch (...) {
        *output = "Error: unknown filter type";
        return false;
    }

    signal::WindowType win_type = signal::WindowType::Hamming;
    if (arg_list.size() >= 4 && !std::isdigit(static_cast<unsigned char>(arg_list[3][0]))) {
        try {
            std::string win_str = arg_list[3];
            if (win_str.front() == '"' || win_str.front() == '\'') win_str = win_str.substr(1, win_str.size() - 2);
            win_type = signal::string_to_window_type(win_str);
        } catch (...) {
            // 使用默认值
        }
    }

    try {
        signal::FilterCoefficients coeffs;
        if (type == signal::FilterType::BandPass || type == signal::FilterType::BandStop) {
            if (arg_list.size() < 4) {
                *output = "Error: band filters require two cutoff frequencies";
                return false;
            }
            double cutoff2 = parser.evaluate(arg_list[3]);
            coeffs = signal::design_fir_band(order, cutoff, cutoff2, type, win_type);
        } else {
            coeffs = signal::design_fir(order, cutoff, type, win_type);
        }

        std::ostringstream oss;
        oss << "b = " << format_vector(coeffs.b) << "\n";
        oss << "a = " << format_vector(coeffs.a);
        *output = oss.str();
        return true;
    } catch (const std::exception& e) {
        *output = std::string("Error: ") + e.what();
        return false;
    }
}

bool handle_iir_design_command(const SignalContext& ctx,
                               const std::string& args,
                               std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.size() < 3) {
        *output = "Error: iir_design requires order, cutoff, and type arguments";
        return false;
    }

    UnifiedExpressionParser parser(VariableResolver(ctx.variables, nullptr), ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);

    int order = 0;
    double cutoff = 0.0;

    try {
        order = static_cast<int>(parser.evaluate(arg_list[0]));
        cutoff = parser.evaluate(arg_list[1]);
    } catch (...) {
        *output = "Error: invalid order or cutoff";
        return false;
    }

    std::string type_str = arg_list[2];
    if (type_str.front() == '"' || type_str.front() == '\'') type_str = type_str.substr(1, type_str.size() - 2);

    signal::FilterType type;
    try {
        type = signal::string_to_filter_type(type_str);
    } catch (...) {
        *output = "Error: unknown filter type";
        return false;
    }

    try {
        signal::FilterCoefficients coeffs = signal::design_butterworth(order, cutoff, type);

        std::ostringstream oss;
        oss << "b = " << format_vector(coeffs.b) << "\n";
        oss << "a = " << format_vector(coeffs.a);
        *output = oss.str();
        return true;
    } catch (const std::exception& e) {
        *output = std::string("Error: ") + e.what();
        return false;
    }
}

bool handle_freqz_command(const SignalContext& ctx,
                          const std::string& args,
                          std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.size() < 2) {
        *output = "Error: freqz requires b and a arguments";
        return false;
    }

    std::vector<double> b = parse_vector(ctx, arg_list[0]);
    std::vector<double> a = parse_vector(ctx, arg_list[1]);

    if (b.empty() || a.empty()) {
        *output = "Error: empty coefficients";
        return false;
    }

    std::size_t n = 512;
    if (arg_list.size() >= 3) {
        try {
            UnifiedExpressionParser parser(VariableResolver(ctx.variables, nullptr), ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);
            n = static_cast<std::size_t>(parser.evaluate(arg_list[2]));
        } catch (...) {
            // 使用默认值
        }
    }

    try {
        std::vector<signal::Complex> response = signal::freqz(b, a, n);
        *output = format_complex_vector(response);
        return true;
    } catch (const std::exception& e) {
        *output = std::string("Error: ") + e.what();
        return false;
    }
}

// ============================================================================
// 时频分析命令
// ============================================================================

bool handle_psd_command(const SignalContext& ctx,
                        const std::string& args,
                        std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.empty()) {
        *output = "Error: psd requires a signal argument";
        return false;
    }

    std::vector<double> signal = parse_vector(ctx, arg_list[0]);
    if (signal.empty()) {
        *output = "Error: empty signal";
        return false;
    }

    std::size_t nfft = 256;
    if (arg_list.size() >= 2) {
        try {
            UnifiedExpressionParser parser(VariableResolver(ctx.variables, nullptr), ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);
            nfft = static_cast<std::size_t>(parser.evaluate(arg_list[1]));
        } catch (...) {
            // 使用默认值
        }
    }

    std::vector<double> result = signal::pwelch(signal, nfft);
    *output = format_vector(result);
    return true;
}

bool handle_stft_command(const SignalContext& ctx,
                         const std::string& args,
                         std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.empty()) {
        *output = "Error: stft requires a signal argument";
        return false;
    }

    std::vector<double> signal = parse_vector(ctx, arg_list[0]);
    if (signal.empty()) {
        *output = "Error: empty signal";
        return false;
    }

    std::size_t nfft = 256;
    if (arg_list.size() >= 2) {
        try {
            UnifiedExpressionParser parser(VariableResolver(ctx.variables, nullptr), ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);
            nfft = static_cast<std::size_t>(parser.evaluate(arg_list[1]));
        } catch (...) {
            // 使用默认值
        }
    }

    signal::STFTResult result = signal::stft(signal, nfft);

    std::ostringstream oss;
    oss << "STFT result:\n";
    oss << "Frames: " << result.stft_matrix.size() << "\n";
    oss << "Frequencies: " << result.freq_vector.size() << "\n";
    oss << "Time vector: " << format_vector(result.time_vector) << "\n";
    oss << "Freq vector: " << format_vector(result.freq_vector);

    *output = oss.str();
    return true;
}

bool handle_spectrogram_command(const SignalContext& ctx,
                                const std::string& args,
                                std::string* output) {
    auto arg_list = parse_args(args);
    if (arg_list.empty()) {
        *output = "Error: spectrogram requires a signal argument";
        return false;
    }

    std::vector<double> signal = parse_vector(ctx, arg_list[0]);
    if (signal.empty()) {
        *output = "Error: empty signal";
        return false;
    }

    std::size_t nfft = 256;
    if (arg_list.size() >= 2) {
        try {
            UnifiedExpressionParser parser(VariableResolver(ctx.variables, nullptr), ctx.functions, ctx.scalar_functions, nullptr, nullptr, ctx.has_script_function, ctx.invoke_script_function);
            nfft = static_cast<std::size_t>(parser.evaluate(arg_list[1]));
        } catch (...) {
            // 使用默认值
        }
    }

    std::vector<std::vector<double>> spec = signal::spectrogram(signal, nfft);

    std::ostringstream oss;
    oss << "Spectrogram: " << spec.size() << " frequency bins x ";
    if (!spec.empty()) {
        oss << spec[0].size() << " time frames";
    }

    *output = oss.str();
    return true;
}

// ============================================================================
// 命令注册
// ============================================================================

bool is_signal_command(const std::string& command) {
    static const std::vector<std::string> commands = {
        "fft", "ifft", "rfft", "irfft",
        "conv", "cconv", "xcorr", "autocorr",
        "window",
        "filter", "fir_design", "iir_design", "freqz",
        "psd", "stft", "spectrogram", "residue"
    };

    return std::find(commands.begin(), commands.end(), command) != commands.end();
}

bool handle_signal_command(const SignalContext& ctx,
                           const std::string& command,
                           const std::string& args,
                           std::string* output) {
    if (command == "fft") return handle_fft_command(ctx, args, output);
    if (command == "ifft") return handle_ifft_command(ctx, args, output);
    if (command == "rfft") return handle_rfft_command(ctx, args, output);

    if (command == "conv") return handle_conv_command(ctx, args, output);
    if (command == "cconv") return handle_cconv_command(ctx, args, output);
    if (command == "xcorr") return handle_xcorr_command(ctx, args, output);
    if (command == "autocorr") return handle_autocorr_command(ctx, args, output);

    if (command == "window") return handle_window_command(ctx, args, output);

    if (command == "filter") return handle_filter_command(ctx, args, output);
    if (command == "fir_design") return handle_fir_design_command(ctx, args, output);
    if (command == "iir_design") return handle_iir_design_command(ctx, args, output);
    if (command == "freqz") return handle_freqz_command(ctx, args, output);

    if (command == "psd") return handle_psd_command(ctx, args, output);
    if (command == "stft") return handle_stft_command(ctx, args, output);
    if (command == "spectrogram") return handle_spectrogram_command(ctx, args, output);

    *output = "Unknown signal command: " + command;
    return false;
}

std::string signal_help_text() {
    return R"(
Signal Processing Commands:
===========================

FFT Commands:
  fft(signal)           - Fast Fourier Transform
  ifft(spectrum)        - Inverse FFT
  rfft(signal)          - Real FFT (optimized for real signals)

Convolution & Correlation:
  conv(s1, s2)          - Linear convolution
  cconv(s1, s2 [, n])   - Circular convolution
  xcorr(s1, s2)         - Cross-correlation
  autocorr(signal)      - Auto-correlation

Window Functions:
  window(type, length [, param])
    Types: rectangular, hanning, hamming, blackman, blackman-harris,
           bartlett, kaiser, gaussian, flattop, tukey

Filter Design:
  fir_design(order, cutoff, type [, window])
    Types: lowpass, highpass, bandpass, bandstop
  iir_design(order, cutoff, type)
    Butterworth filter design
  filter(b, a, signal)  - Apply digital filter
  freqz(b, a [, n])     - Frequency response

Time-Frequency Analysis:
  psd(signal [, nfft])  - Power spectral density (Welch method)
  stft(signal [, nfft]) - Short-time Fourier transform
  spectrogram(signal [, nfft]) - Spectrogram

Examples:
  fft([1, 2, 3, 4])
  conv([1, 2], [1, 1, 1])
  window(hamming, 64)
  fir_design(32, 0.2, lowpass)
  psd([1, 2, 3, 4, 5, 6, 7, 8], 8)
)";
}

}  // namespace signal_cmds

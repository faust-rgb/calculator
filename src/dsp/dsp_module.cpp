/**
 * @file dsp_module.cpp
 * @brief 数字信号处理模块实现
 *
 * 本文件实现了 DSP 模块的命令处理和函数注册：
 * - residue: 计算有理函数在指定点的留数（命令形式）
 * - 窗函数: hann, hamming, blackman, bartlett, kaiser, gaussian, flattop, tukey
 * - 变换与相关: fft, ifft, dft, idft, conv, convolve, xcorr, autocorr
 * - 时频与滤波: filter, freqz, fir1, periodogram, pwelch, residue_pf
 */

#include "dsp_module.h"
#include "core/services/core_manager_interfaces.h"
#include "residue.h"
#include "signal_processing.h"
#include "matrix_dsp.h"
#include "matrix_internal.h"
#include "math/mymath.h"
#include "core/services/string_utils.h"
#include <stdexcept>
#include <memory>

namespace {

using matrix::Matrix;
using matrix::filter;
using matrix::freqz;
using matrix::residue;
using Scalar = mymath::Scalar;

Matrix require_matrix(const StoredValue& val, const std::string& func_name) {
    if (!val.is_matrix_type() || !val.get_matrix_ptr()) {
        throw std::runtime_error(func_name + " expects a matrix argument");
    }
    return *val.get_matrix_ptr();
}

std::vector<Scalar> matrix_to_vector(const Matrix& m) {
    std::vector<Scalar> vec;
    vec.reserve(m.rows * m.cols);
    for (std::size_t r = 0; r < m.rows; ++r) {
        for (std::size_t c = 0; c < m.cols; ++c) {
            vec.push_back(m.at(r, c));
        }
    }
    return vec;
}

Matrix vector_to_matrix(const std::vector<Scalar>& vec) {
    Matrix m(1, vec.size(), 0.0L);
    for (std::size_t i = 0; i < vec.size(); ++i) {
        m.at(0, i) = vec[i];
    }
    return m;
}

std::size_t require_positive_length(const StoredValue& value, const std::string& name) {
    if (value.is_matrix_type() || value.is_complex_type() ||
        !mymath::is_integer(value.as_scalar()) || value.as_scalar() <= Scalar(0)) {
        throw std::runtime_error(name + " expects a positive integer length");
    }
    return static_cast<std::size_t>(static_cast<long double>(value.as_scalar()) + 0.5L);
}

StoredValue make_dsp_window(const std::vector<StoredValue>& args, signal::WindowType type, const std::string& name) {
    if (args.empty() || args.size() > 2) throw std::runtime_error(name + " expects length and optional parameter");
    const std::size_t n = require_positive_length(args[0], name);
    std::vector<Scalar> win;
    if (type == signal::WindowType::Kaiser) {
        Scalar beta = (args.size() == 2) ? args[1].as_scalar() : Scalar(5.0L);
        win = signal::kaiser_window(n, beta);
    } else if (type == signal::WindowType::Gaussian) {
        Scalar sigma = (args.size() == 2) ? args[1].as_scalar() : Scalar(0.5L);
        win = signal::gaussian_window(n, sigma);
    } else if (type == signal::WindowType::Tukey) {
        Scalar alpha = (args.size() == 2) ? args[1].as_scalar() : Scalar(0.5L);
        win = signal::tukey_window(n, alpha);
    } else {
        win = signal::window(type, n);
    }
    return StoredValue(vector_to_matrix(win));
}

StoredValue make_dsp_transform_result(const std::vector<StoredValue>& args,
                                      const std::string& name, bool inverse) {
    if (args.size() != 1) throw std::runtime_error(name + " expects one sequence argument");
    const Matrix input = require_matrix(args[0], name);
    return StoredValue(matrix::internal::complex_sequence_to_matrix(
        matrix::internal::discrete_fourier_transform(
            matrix::internal::as_complex_sequence(input, name), inverse),
        inverse));
}

StoredValue make_dsp_convolution_result(const std::vector<StoredValue>& args,
                                        const std::string& name) {
    if (args.size() != 2) throw std::runtime_error(name + " expects two sequence arguments");
    return StoredValue(matrix::internal::complex_sequence_to_matrix(
        matrix::internal::convolve_sequences(
            matrix::internal::as_complex_sequence(require_matrix(args[0], name), name),
            matrix::internal::as_complex_sequence(require_matrix(args[1], name), name)),
        true));
}

} // namespace

// ============================================================================
// 命令接口实现
// ============================================================================

std::vector<std::string> DspModule::get_commands() const {
    return {"residue"};
}

std::string DspModule::execute_args_view(std::string_view command,
                                    const std::vector<std::string_view>& args,
                                    ModuleServices& module_services) {
    std::vector<std::string> string_args;
    for (const auto& arg : args) string_args.emplace_back(arg);

    return dsp_ops::handle_residue_command(std::string(command), string_args, module_services);
}

// ============================================================================
// 函数接口实现 — DSP 函数注册
// ============================================================================

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
DspModule::get_functions_map() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    // 窗函数
    funcs["hann"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, signal::WindowType::Hanning, "hann"); };
    funcs["hanning"] = funcs["hann"];
    funcs["hamming"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, signal::WindowType::Hamming, "hamming"); };
    funcs["blackman"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, signal::WindowType::Blackman, "blackman"); };
    funcs["bartlett"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, signal::WindowType::Bartlett, "bartlett"); };
    funcs["kaiser"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, signal::WindowType::Kaiser, "kaiser"); };
    funcs["gaussian"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, signal::WindowType::Gaussian, "gaussian"); };
    funcs["flattop"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, signal::WindowType::Flattop, "flattop"); };
    funcs["tukey"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, signal::WindowType::Tukey, "tukey"); };

    // 变换与卷积
    funcs["dft"] = [](const std::vector<StoredValue>& args) { return make_dsp_transform_result(args, "dft", false); };
    funcs["fft"] = funcs["dft"];
    funcs["idft"] = [](const std::vector<StoredValue>& args) { return make_dsp_transform_result(args, "idft", true); };
    funcs["ifft"] = funcs["idft"];
    funcs["convolve"] = [](const std::vector<StoredValue>& args) { return make_dsp_convolution_result(args, "convolve"); };
    funcs["conv"] = funcs["convolve"];

    // 互相关与自相关
    funcs["xcorr"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("xcorr expects 2 sequence arguments");
        Matrix m1 = require_matrix(args[0], "xcorr");
        Matrix m2 = require_matrix(args[1], "xcorr");
        return StoredValue(vector_to_matrix(signal::xcorr(matrix_to_vector(m1), matrix_to_vector(m2))));
    };
    funcs["autocorr"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("autocorr expects 1 sequence argument");
        Matrix m = require_matrix(args[0], "autocorr");
        return StoredValue(vector_to_matrix(signal::autocorr(matrix_to_vector(m))));
    };

    // 功率谱密度
    funcs["periodogram"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty() || args.size() > 2) throw std::runtime_error("periodogram expects 1 or 2 arguments");
        Matrix m = require_matrix(args[0], "periodogram");
        Scalar fs = (args.size() == 2) ? args[1].as_scalar() : Scalar(1.0L);
        return StoredValue(vector_to_matrix(signal::periodogram(matrix_to_vector(m), fs)));
    };
    funcs["pwelch"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty() || args.size() > 4) throw std::runtime_error("pwelch expects 1 to 4 arguments");
        Matrix m = require_matrix(args[0], "pwelch");
        std::size_t nfft = (args.size() >= 2) ? require_positive_length(args[1], "pwelch") : 0;
        std::size_t noverlap = (args.size() >= 3) ? static_cast<std::size_t>(static_cast<long double>(args[2].as_scalar()) + 0.5L) : 0;
        Scalar fs = (args.size() >= 4) ? args[3].as_scalar() : Scalar(1.0L);
        return StoredValue(vector_to_matrix(signal::pwelch(matrix_to_vector(m), nfft, signal::WindowType::Hanning, noverlap, fs)));
    };

    // 数字滤波器
    funcs["filter"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 3) throw std::runtime_error("filter expects 3 arguments");
        Matrix b = require_matrix(args[0], "filter");
        Matrix a = require_matrix(args[1], "filter");
        Matrix x = require_matrix(args[2], "filter");
        return StoredValue(filter(b, a, x));
    };

    // FIR 滤波器设计
    funcs["fir1"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2 || args.size() > 3) throw std::runtime_error("fir1(order, cutoff, [type]) expects 2 or 3 arguments");
        int order = static_cast<int>(static_cast<long double>(args[0].as_scalar()) + 0.5L);
        Scalar cutoff = args[1].as_scalar();
        signal::FilterType ftype = signal::FilterType::LowPass;
        if (args.size() == 3 && args[2].is_string_type()) {
            ftype = signal::string_to_filter_type(args[2].as_string());
        }
        auto fir = signal::design_fir(order, cutoff, ftype);
        return StoredValue(vector_to_matrix(fir.b));
    };

    // 频率响应
    funcs["freqz"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2 || args.size() > 3) throw std::runtime_error("freqz expects 2 or 3 arguments");
        Matrix b = require_matrix(args[0], "freqz");
        Matrix a = require_matrix(args[1], "freqz");
        std::size_t n = 512;
        if (args.size() == 3) {
            Matrix nm = require_matrix(args[2], "freqz");
            if (nm.rows > 0 && nm.cols > 0) {
                n = static_cast<std::size_t>(static_cast<long double>(nm.at(0, 0)) + 0.5);
            }
        }
        return StoredValue(freqz(b, a, n));
    };

    // 部分分式展开（矩阵版）
    funcs["residue_pf"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("residue_pf expects 2 arguments");
        Matrix b = require_matrix(args[0], "residue_pf");
        Matrix a = require_matrix(args[1], "residue_pf");
        return StoredValue(residue(b, a));
    };

    return funcs;
}

std::vector<std::string> DspModule::get_function_names() const {
    static const std::vector<std::string> names = {
        "autocorr", "bartlett", "blackman", "conv", "convolve", "dft", "fft",
        "filter", "fir1", "flattop", "freqz", "gaussian", "hamming", "hann",
        "hanning", "idft", "ifft", "kaiser", "periodogram", "pwelch", "residue_pf",
        "tukey", "xcorr"
    };
    return names;
}

std::string DspModule::get_help_snippet(const std::string& topic) const {
    if (topic == "dsp") {
        return "DSP commands:\n"
               "  residue(expr, var, point) - Compute residue of rational function at point\n"
               "DSP functions:\n"
               "  filter(b, a, x) - Digital filter\n"
               "  fir1(order, cutoff, [type]) - Window-based FIR filter design\n"
               "  freqz(b, a, n) - Frequency response\n"
               "  fft(x), ifft(x), dft(x), idft(x) - Fourier transforms\n"
               "  conv(x, y), convolve(x, y) - Convolution\n"
               "  xcorr(x, y), autocorr(x) - Cross- and auto-correlation\n"
               "  periodogram(x, [fs]), pwelch(x, [nfft], [overlap], [fs]) - Power spectral density\n"
               "  hann(n), hamming(n), blackman(n), bartlett(n), kaiser(n, [beta]), gaussian(n, [sigma]), flattop(n), tukey(n, [alpha]) - Window functions\n"
               "  residue_pf(b, a) - Partial fraction expansion";
    }
    return "";
}

/**
 * @file calculator_signal_commands.h
 * @brief 信号处理命令接口
 *
 * 提供信号处理命令的计算逻辑，包括：
 * - FFT 相关命令 (fft, ifft, rfft, irfft)
 * - 卷积相关命令 (conv, xcorr, autocorr)
 * - 窗函数命令 (window)
 * - 滤波器命令 (filter, fir_design, iir_design)
 * - 时频分析命令 (psd, stft, spectrogram)
 */

#ifndef CALCULATOR_SIGNAL_COMMANDS_H
#define CALCULATOR_SIGNAL_COMMANDS_H

#include "calculator_internal_types.h"
#include "signal_processing.h"

#include <functional>
#include <map>
#include <string>
#include <vector>

namespace signal_cmds {

// ============================================================================
// 信号数据结构
// ============================================================================

/**
 * @brief 信号数据（用于命令间传递）
 */
struct SignalData {
    std::vector<double> samples;
    double sample_rate = 1.0;
};

/**
 * @brief 复数信号数据
 */
struct ComplexSignalData {
    std::vector<signal::Complex> samples;
    double sample_rate = 1.0;
};

/**
 * @brief 信号处理上下文
 */
struct SignalContext {
    const std::map<std::string, StoredValue>* variables;
    const std::map<std::string, CustomFunction>* functions;
    const std::map<std::string, std::function<double(const std::vector<double>&)>>* scalar_functions;
    std::function<bool(const std::string&)> has_script_function;
    std::function<double(const std::string&, const std::vector<double>&)> invoke_script_function;
    std::function<SignalData(const std::string&, std::string*)> resolve_signal;
};

// ============================================================================
// FFT 命令
// ============================================================================

/**
 * @brief FFT 命令处理
 * @param ctx 信号处理上下文
 * @param args 参数字符串
 * @param output 输出字符串
 * @return 是否成功处理
 */
bool handle_fft_command(const SignalContext& ctx,
                        const std::string& args,
                        std::string* output);

/**
 * @brief IFFT 命令处理
 */
bool handle_ifft_command(const SignalContext& ctx,
                         const std::string& args,
                         std::string* output);

/**
 * @brief 实数 FFT 命令处理
 */
bool handle_rfft_command(const SignalContext& ctx,
                         const std::string& args,
                         std::string* output);

// ============================================================================
// 卷积与相关命令
// ============================================================================

/**
 * @brief 卷积命令处理
 */
bool handle_conv_command(const SignalContext& ctx,
                         const std::string& args,
                         std::string* output);

/**
 * @brief 循环卷积命令处理
 */
bool handle_cconv_command(const SignalContext& ctx,
                          const std::string& args,
                          std::string* output);

/**
 * @brief 互相关命令处理
 */
bool handle_xcorr_command(const SignalContext& ctx,
                          const std::string& args,
                          std::string* output);

/**
 * @brief 自相关命令处理
 */
bool handle_autocorr_command(const SignalContext& ctx,
                             const std::string& args,
                             std::string* output);

// ============================================================================
// 窗函数命令
// ============================================================================

/**
 * @brief 窗函数命令处理
 */
bool handle_window_command(const SignalContext& ctx,
                           const std::string& args,
                           std::string* output);

// ============================================================================
// 滤波器命令
// ============================================================================

/**
 * @brief 滤波命令处理
 */
bool handle_filter_command(const SignalContext& ctx,
                           const std::string& args,
                           std::string* output);

/**
 * @brief FIR 滤波器设计命令处理
 */
bool handle_fir_design_command(const SignalContext& ctx,
                               const std::string& args,
                               std::string* output);

/**
 * @brief IIR 滤波器设计命令处理
 */
bool handle_iir_design_command(const SignalContext& ctx,
                               const std::string& args,
                               std::string* output);

/**
 * @brief 频率响应命令处理
 */
bool handle_freqz_command(const SignalContext& ctx,
                          const std::string& args,
                          std::string* output);

// ============================================================================
// 时频分析命令
// ============================================================================

/**
 * @brief 功率谱密度命令处理
 */
bool handle_psd_command(const SignalContext& ctx,
                        const std::string& args,
                        std::string* output);

/**
 * @brief STFT 命令处理
 */
bool handle_stft_command(const SignalContext& ctx,
                         const std::string& args,
                         std::string* output);

/**
 * @brief 语谱图命令处理
 */
bool handle_spectrogram_command(const SignalContext& ctx,
                                const std::string& args,
                                std::string* output);

// ============================================================================
// 命令注册
// ============================================================================

/**
 * @brief 检查是否为信号处理命令
 */
bool is_signal_command(const std::string& command);

/**
 * @brief 处理信号处理命令
 */
bool handle_signal_command(const SignalContext& ctx,
                           const std::string& command,
                           const std::string& args,
                           std::string* output);

/**
 * @brief 获取信号处理帮助文本
 */
std::string signal_help_text();

}  // namespace signal_cmds

#endif  // CALCULATOR_SIGNAL_COMMANDS_H

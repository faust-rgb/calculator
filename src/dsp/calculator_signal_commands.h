/**
 * @file calculator_signal_commands.h
 * @brief 信号处理命令接口
 *
 * 本头文件定义了信号处理命令的接口，包括：
 * - FFT 相关命令 (fft, ifft, rfft, irfft)：快速傅里叶变换及其逆变换
 * - 卷积相关命令 (conv, xcorr, autocorr)：线性卷积、循环卷积和相关函数
 * - 窗函数命令 (window)：生成各类窗函数
 * - 滤波器命令 (filter, fir_design, iir_design)：数字滤波器设计和应用
 * - 时频分析命令 (psd, stft, spectrogram)：功率谱和时频分析
 *
 * 这些命令为计算器提供完整的数字信号处理功能支持。
 *
 * @author Calculator Development Team
 * @date 2024
 */

#ifndef CALCULATOR_SIGNAL_COMMANDS_H
#define CALCULATOR_SIGNAL_COMMANDS_H

#include "core/types/module_types.h"
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
 * @struct SignalData
 * @brief 信号数据结构
 *
 * 用于在命令间传递实数信号数据，包含采样值和采样率信息。
 */
struct SignalData {
    std::vector<Scalar> samples;
    Scalar sample_rate = 1.0L;
};

/**
 * @struct ComplexSignalData
 * @brief 复数信号数据结构
 *
 * 用于在命令间传递复数信号数据，包含复数采样值和采样率信息。
 */
struct ComplexSignalData {
    std::vector<signal::Complex> samples;
    Scalar sample_rate = 1.0L;
};

/**
 * @struct SignalContext
 * @brief 信号处理上下文结构
 *
 * 包含信号处理命令执行所需的上下文信息，
 * 如变量表、函数表、标量函数表以及脚本函数回调等。
 */
struct SignalContext {
    const std::map<std::string, StoredValue>* variables;
    const std::map<std::string, CustomFunction>* functions;
    const std::map<std::string, std::function<Scalar(const std::vector<Scalar>&)>>* scalar_functions;
    std::function<bool(const std::string&)> has_script_function;
    std::function<Scalar(const std::string&, const std::vector<Scalar>&)> invoke_script_function;
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

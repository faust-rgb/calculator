#ifndef SIGNAL_PROCESSING_H
#define SIGNAL_PROCESSING_H

/**
 * @file signal_processing.h
 * @brief 信号处理库
 *
 * 提供完整的数字信号处理功能，包括：
 * - 高级 FFT 算法（混合基、Bluestein、实数 FFT）
 * - 卷积与相关分析
 * - 窗函数
 * - 滤波器设计
 * - 时频分析（PSD、STFT）
 */

#include "mymath_complex.h"

#include <cstddef>
#include <string>
#include <vector>

namespace signal {

// ============================================================================
// 复数类型定义
// ============================================================================

using Complex = mymath::complex<double>;

// ============================================================================
// 窗函数类型定义
// ============================================================================

/**
 * @enum WindowType
 * @brief 窗函数类型
 */
enum class WindowType {
    Rectangular,    ///< 矩形窗
    Hanning,        ///< 汉宁窗
    Hamming,        ///< 汉明窗
    Blackman,       ///< 布莱克曼窗
    BlackmanHarris, ///< 布莱克曼-哈里斯窗
    Bartlett,       ///< 巴特利特窗（三角窗）
    Kaiser,         ///< 凯撒窗
    Gaussian,       ///< 高斯窗
    Flattop,        ///< 平顶窗
    Tukey           ///< Tukey 窗
};

// ============================================================================
// FFT 相关结构体
// ============================================================================

/**
 * @struct FFTResult
 * @brief FFT 结果结构
 */
struct FFTResult {
    std::vector<Complex> spectrum;  ///< 频谱数据
    std::vector<double> frequencies; ///< 频率轴（可选）
    double sample_rate = 1.0;       ///< 采样率
};

/**
 * @struct STFTResult
 * @brief 短时傅里叶变换结果
 */
struct STFTResult {
    std::vector<std::vector<Complex>> stft_matrix; ///< STFT 矩阵（时间 × 频率）
    std::vector<double> time_vector;               ///< 时间轴
    std::vector<double> freq_vector;               ///< 频率轴
    WindowType window_type = WindowType::Hanning;  ///< 使用的窗函数
};

// ============================================================================
// FFT 高级算法
// ============================================================================

/**
 * @brief 基-2 FFT（Cooley-Tukey 算法）
 * @param input 输入信号
 * @return 频谱
 * @throw std::runtime_error 当输入长度不是 2 的幂次时抛出
 */
std::vector<Complex> fft_radix2(const std::vector<Complex>& input);

/**
 * @brief 混合基 FFT（支持任意长度）
 * @param input 输入信号
 * @return 频谱
 *
 * 使用 Cooley-Tukey 混合基算法，自动分解为最优因子组合。
 */
std::vector<Complex> fft_mixed_radix(const std::vector<Complex>& input);

/**
 * @brief Bluestein FFT（支持任意长度的 Chirp-Z 变换）
 * @param input 输入信号
 * @return 频谱
 *
 * 使用 Bluestein 算法，支持任意长度 N 的 FFT，
 * 时间复杂度 O(N log N)。
 */
std::vector<Complex> fft_bluestein(const std::vector<Complex>& input);

/**
 * @brief 自动选择最优 FFT 算法
 * @param input 输入信号
 * @return 频谱
 *
 * 根据输入长度自动选择最优算法：
 * - 2 的幂次：基-2 FFT
 * - 高合成数：混合基 FFT
 * - 其他：Bluestein FFT
 */
std::vector<Complex> fft(const std::vector<Complex>& input);

/**
 * @brief 实数 FFT（针对实数信号优化的 FFT）
 * @param input 实数输入信号
 * @return 频谱（仅包含正频率部分，利用 Hermitian 对称性）
 */
std::vector<Complex> rfft(const std::vector<double>& input);

/**
 * @brief 实数 IFFT
 * @param spectrum 频谱（正频率部分）
 * @param n 原始信号长度
 * @return 实数信号
 */
std::vector<double> irfft(const std::vector<Complex>& spectrum, std::size_t n);

/**
 * @brief 逆 FFT
 * @param spectrum 频谱
 * @return 时域信号
 */
std::vector<Complex> ifft(const std::vector<Complex>& spectrum);

/**
 * @brief FFT 频率轴生成
 * @param n FFT 点数
 * @param sample_rate 采样率
 * @return 频率向量
 */
std::vector<double> fft_frequencies(std::size_t n, double sample_rate);

/**
 * @brief FFT 频移（将零频率移到中心）
 * @param spectrum 频谱
 * @return 频移后的频谱
 */
std::vector<Complex> fftshift(const std::vector<Complex>& spectrum);

/**
 * @brief 逆 FFT 频移
 * @param spectrum 频移后的频谱
 * @return 原始频谱
 */
std::vector<Complex> ifftshift(const std::vector<Complex>& spectrum);

// ============================================================================
// 卷积与相关
// ============================================================================

/**
 * @brief 线性卷积
 * @param signal1 第一个信号
 * @param signal2 第二个信号
 * @return 卷积结果（长度 = n + m - 1）
 */
std::vector<double> convolve(const std::vector<double>& signal1,
                              const std::vector<double>& signal2);

/**
 * @brief FFT 快速卷积
 * @param signal1 第一个信号
 * @param signal2 第二个信号
 * @return 卷积结果
 */
std::vector<double> fft_convolve(const std::vector<double>& signal1,
                                  const std::vector<double>& signal2);

/**
 * @brief 循环卷积（圆周卷积）
 * @param signal1 第一个信号
 * @param signal2 第二个信号
 * @param n 卷积长度（默认为两信号长度的最小公倍数）
 * @return 循环卷积结果
 */
std::vector<double> circular_convolve(const std::vector<double>& signal1,
                                       const std::vector<double>& signal2,
                                       std::size_t n = 0);

/**
 * @brief 互相关函数
 * @param signal1 第一个信号
 * @param signal2 第二个信号
 * @return 互相关结果
 */
std::vector<double> xcorr(const std::vector<double>& signal1,
                          const std::vector<double>& signal2);

/**
 * @brief 自相关函数
 * @param signal 输入信号
 * @return 自相关结果
 */
std::vector<double> autocorr(const std::vector<double>& signal);

/**
 * @brief 归一化互相关
 * @param signal1 第一个信号
 * @param signal2 第二个信号
 * @return 归一化互相关（范围 [-1, 1]）
 */
std::vector<double> normalized_xcorr(const std::vector<double>& signal1,
                                      const std::vector<double>& signal2);

// ============================================================================
// 窗函数
// ============================================================================

/**
 * @brief 生成窗函数
 * @param type 窗函数类型
 * @param length 窗长度
 * @param param 可选参数（如 Kaiser 窗的 beta 值）
 * @return 窗函数系数
 */
std::vector<double> window(WindowType type, std::size_t length, double param = 0.0);

/**
 * @brief 矩形窗
 */
std::vector<double> rectangular_window(std::size_t length);

/**
 * @brief 汉宁窗
 */
std::vector<double> hanning_window(std::size_t length);

/**
 * @brief 汉明窗
 */
std::vector<double> hamming_window(std::size_t length);

/**
 * @brief 布莱克曼窗
 */
std::vector<double> blackman_window(std::size_t length);

/**
 * @brief 布莱克曼-哈里斯窗
 */
std::vector<double> blackman_harris_window(std::size_t length);

/**
 * @brief 巴特利特窗（三角窗）
 */
std::vector<double> bartlett_window(std::size_t length);

/**
 * @brief 凯撒窗
 * @param length 窗长度
 * @param beta 形状参数（默认 5.0，控制旁瓣衰减）
 */
std::vector<double> kaiser_window(std::size_t length, double beta = 5.0);

/**
 * @brief 高斯窗
 * @param length 窗长度
 * @param sigma 标准差（相对于窗长度的比例，默认 0.4）
 */
std::vector<double> gaussian_window(std::size_t length, double sigma = 0.4);

/**
 * @brief 平顶窗
 */
std::vector<double> flattop_window(std::size_t length);

/**
 * @brief Tukey 窗
 * @param length 窗长度
 * @param alpha 锥度比例（0 = 矩形窗，1 = 汉宁窗）
 */
std::vector<double> tukey_window(std::size_t length, double alpha = 0.5);

/**
 * @brief 应用窗函数到信号
 * @param signal 输入信号
 * @param window_type 窗函数类型
 * @param param 可选参数
 * @return 加窗后的信号
 */
std::vector<double> apply_window(const std::vector<double>& signal,
                                  WindowType window_type,
                                  double param = 0.0);

// ============================================================================
// 滤波器设计
// ============================================================================

/**
 * @enum FilterType
 * @brief 滤波器类型
 */
enum class FilterType {
    LowPass,    ///< 低通滤波器
    HighPass,   ///< 高通滤波器
    BandPass,   ///< 带通滤波器
    BandStop    ///< 带阻滤波器
};
/**
 * @struct FilterCoefficients
 * @brief 滤波器系数
 */
struct FilterCoefficients {
    std::vector<double> b;  ///< 分子系数
    std::vector<double> a;  ///< 分母系数
};

/**
 * @struct SOS
 * @brief 二阶节 (Second-Order Section) 系数
 */
struct SOS {
    double b0, b1, b2;  ///< 分子系数
    double a1, a2;      ///< 分母系数 (a0 默认为 1)
};

/**
 * @brief 将 [b, a] 系数转换为 SOS 形式
 */
std::vector<SOS> tf2sos(const std::vector<double>& b, const std::vector<double>& a);

/**
 * @brief 使用 SOS 形式应用滤波器
 */
std::vector<double> sosfilter(const std::vector<SOS>& sos, const std::vector<double>& signal);

/**
 * @brief 线性卷积（针对实数信号优化的 FFT 卷积）
 */
std::vector<double> fast_convolve(const std::vector<double>& signal1,
                                   const std::vector<double>& signal2);

/**
 * @brief 设计 FIR 滤波器（窗函数法）
 * @param order 滤波器阶数
 * @param cutoff 截止频率（归一化，0-1）
 * @param type 滤波器类型
 * @param window_type 窗函数类型
 * @return 滤波器系数
 */
FilterCoefficients design_fir(int order,
                               double cutoff,
                               FilterType type,
                               WindowType window_type = WindowType::Hamming);

/**
 * @brief 设计带通/带阻 FIR 滤波器
 * @param order 滤波器阶数
 * @param cutoff_low 低截止频率
 * @param cutoff_high 高截止频率
 * @param type 滤波器类型（BandPass 或 BandStop）
 * @param window_type 窗函数类型
 * @return 滤波器系数
 */
FilterCoefficients design_fir_band(int order,
                                    double cutoff_low,
                                    double cutoff_high,
                                    FilterType type,
                                    WindowType window_type = WindowType::Hamming);

/**
 * @brief 设计巴特沃斯 IIR 滤波器
 * @param order 滤波器阶数
 * @param cutoff 截止频率（归一化）
 * @param type 滤波器类型
 * @return 滤波器系数
 */
FilterCoefficients design_butterworth(int order,
                                       double cutoff,
                                       FilterType type);

/**
 * @brief 设计切比雪夫 I 型 IIR 滤波器
 * @param order 滤波器阶数
 * @param cutoff 截止频率
 * @param ripple 通带波纹（dB）
 * @param type 滤波器类型
 * @return 滤波器系数
 */
FilterCoefficients design_chebyshev1(int order,
                                      double cutoff,
                                      double ripple,
                                      FilterType type);

/**
 * @brief 设计椭圆 IIR 滤波器
 * @param order 滤波器阶数
 * @param cutoff 截止频率
 * @param ripple 通带波纹（dB）
 * @param stopband_atten 阻带衰减（dB）
 * @param type 滤波器类型
 * @return 滤波器系数
 */
FilterCoefficients design_elliptic(int order,
                                    double cutoff,
                                    double ripple,
                                    double stopband_atten,
                                    FilterType type);

/**
 * @brief 应用数字滤波器
 * @param b 分子系数
 * @param a 分母系数
 * @param signal 输入信号
 * @return 滤波后的信号
 */
std::vector<double> filter(const std::vector<double>& b,
                            const std::vector<double>& a,
                            const std::vector<double>& signal);

/**
 * @brief 零相位滤波（前向-后向滤波）
 * @param b 分子系数
 * @param a 分母系数
 * @param signal 输入信号
 * @return 滤波后的信号（零相位失真）
 */
std::vector<double> filtfilt(const std::vector<double>& b,
                              const std::vector<double>& a,
                              const std::vector<double>& signal);

/**
 * @brief 计算滤波器频率响应
 * @param b 分子系数
 * @param a 分母系数
 * @param n 频率点数
 * @return 频率响应（复数）
 */
std::vector<Complex> freqz(const std::vector<double>& b,
                            const std::vector<double>& a,
                            std::size_t n = 512);

/**
 * @brief 计算滤波器群延迟
 * @param b 分子系数
 * @param a 分母系数
 * @param n 频率点数
 * @return 群延迟（采样点数）
 */
std::vector<double> grpdelay(const std::vector<double>& b,
                              const std::vector<double>& a,
                              std::size_t n = 512);

// ============================================================================
// 时频分析
// ============================================================================

/**
 * @brief 计算功率谱密度（Welch 方法）
 * @param signal 输入信号
 * @param nfft FFT 点数
 * @param window_type 窗函数类型
 * @param noverlap 重叠点数
 * @param sample_rate 采样率
 * @return 功率谱密度
 */
std::vector<double> pwelch(const std::vector<double>& signal,
                            std::size_t nfft,
                            WindowType window_type = WindowType::Hanning,
                            std::size_t noverlap = 0,
                            double sample_rate = 1.0);

/**
 * @brief 计算功率谱密度（周期图法）
 * @param signal 输入信号
 * @param sample_rate 采样率
 * @return 功率谱密度
 */
std::vector<double> periodogram(const std::vector<double>& signal, double sample_rate = 1.0);

/**
 * @brief 短时傅里叶变换
 * @param signal 输入信号
 * @param nfft FFT 点数
 * @param window_type 窗函数类型
 * @param noverlap 重叠点数
 * @return STFT 结果
 */
STFTResult stft(const std::vector<double>& signal,
                std::size_t nfft,
                WindowType window_type = WindowType::Hanning,
                std::size_t noverlap = 0);

/**
 * @brief 逆短时傅里叶变换
 * @param stft_result STFT 结果
 * @return 重构的时域信号
 */
std::vector<double> istft(const STFTResult& stft_result);

/**
 * @brief 计算语谱图
 * @param signal 输入信号
 * @param nfft FFT 点数
 * @param window_type 窗函数类型
 * @param noverlap 重叠点数
 * @return 语谱图（功率谱矩阵）
 */
std::vector<std::vector<double>> spectrogram(const std::vector<double>& signal,
                                              std::size_t nfft,
                                              WindowType window_type = WindowType::Hanning,
                                              std::size_t noverlap = 0);

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 检查是否为 2 的幂
 */
bool is_power_of_two(std::size_t n);

/**
 * @brief 计算下一个 2 的幂
 */
std::size_t next_power_of_two(std::size_t n);

/**
 * @brief 计算两个数的最大公约数
 */
std::size_t gcd(std::size_t a, std::size_t b);

/**
 * @brief 计算两个数的最小公倍数
 */
std::size_t lcm(std::size_t a, std::size_t b);

/**
 * @brief 窗函数类型转字符串
 */
std::string window_type_to_string(WindowType type);

/**
 * @brief 字符串转窗函数类型
 */
WindowType string_to_window_type(const std::string& str);

/**
 * @brief 滤波器类型转字符串
 */
std::string filter_type_to_string(FilterType type);

/**
 * @brief 字符串转滤波器类型
 */
FilterType string_to_filter_type(const std::string& str);

}  // namespace signal

#endif  // SIGNAL_PROCESSING_H

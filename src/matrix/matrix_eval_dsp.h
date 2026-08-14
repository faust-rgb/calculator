// ============================================================================
// matrix_eval_dsp.h - 矩阵/向量 DSP 与变换函数求值器
// ============================================================================

#ifndef MATRIX_EVAL_DSP_H
#define MATRIX_EVAL_DSP_H

#include "matrix.h"
#include <string>
#include <vector>
#include <functional>

namespace matrix {
namespace internal {

/**
 * @brief 尝试求值 DSP 与信号变换相关函数
 *
 * 支持：dft, fft, idft, ifft, convolve, conv, freqz, hann, hanning, hamming, blackman
 *
 * @param name 函数名
 * @param arguments 参数表达式列表
 * @param scalar_eval 标量求值回调
 * @param require_matrix 矩阵提取回调
 * @param result 输出结果
 * @return 如果成功分发求值则返回 true
 */
bool try_evaluate_dsp_function(
    const std::string& name,
    const std::vector<std::string>& arguments,
    const ScalarEvaluator& scalar_eval,
    const std::function<Matrix(const std::string&, const std::string&)>& require_matrix,
    Value* result);

} // namespace internal
} // namespace matrix

#endif // MATRIX_EVAL_DSP_H

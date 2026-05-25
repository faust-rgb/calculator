/**
 * @file matrix_dsp.h
 * @brief 矩阵信号处理运算函数声明
 *
 * 从 matrix.h 中分离出的 DSP 相关函数接口，包括：
 * - filter: 数字滤波器
 * - freqz: 频率响应计算
 * - residue: 部分分式展开（留数分解）
 *
 * 这些函数仅由需要信号处理功能的模块使用，
 * 不应被仅需基础矩阵数据结构的模块包含。
 */

#ifndef MATRIX_DSP_H
#define MATRIX_DSP_H

#include "matrix.h"

namespace matrix {

// ============================================================================
// 信号处理运算 (泛型)
// ============================================================================

/** @brief 数字滤波器：计算差分方程 y[n] = (b0*x[n] + b1*x[n-1] + ...) - (a1*y[n-1] + ...) */
template <typename T> TMatrix<T> filter(const TMatrix<T>& b, const TMatrix<T>& a, const TMatrix<T>& x);

/** @brief 频率响应计算：在 n 个等间隔频率点上计算 H(w) */
template <typename T> TMatrix<T> freqz(const TMatrix<T>& b, const TMatrix<T>& a, std::size_t n = 512);

/** @brief 部分分式展开（留数分解）：计算 B(s)/A(s) 的留数、极点和直接项 */
template <typename T> TMatrix<T> residue(const TMatrix<T>& b, const TMatrix<T>& a);

}  // namespace matrix

#endif  // MATRIX_DSP_H
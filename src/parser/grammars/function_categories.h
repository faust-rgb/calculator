// ============================================================================
// 函数分类定义
// ============================================================================
//
// 统一定义矩阵函数和复数函数列表，避免多处硬编码。
// 供表达式分析、编译等模块使用。
//
// ============================================================================

#ifndef PARSER_FUNCTION_CATEGORIES_H
#define PARSER_FUNCTION_CATEGORIES_H

#include <string_view>

// ============================================================================
// 矩阵函数列表
// ============================================================================

// 返回矩阵或操作矩阵的函数
inline constexpr std::string_view kMatrixFunctions[] = {
    "mat", "vec", "zeros", "identity", "diag", "eye", "linspace", "logspace",
    "get", "set", "resize", "append_row", "append_col", "transpose", "inverse",
    "pinv", "dot", "outer", "kron", "hadamard", "null", "least_squares",
    "qr_q", "qr_r", "lu_l", "lu_u", "lu_p", "svd_u", "svd_s", "svd_vt",
    "solve", "norm", "cond", "trace", "det", "rank", "rref", "eigvals",
    "eigvecs", "reshape", "cholesky", "hessenberg", "schur", "poly_eval",
    "poly_deriv", "poly_integ", "poly_compose", "poly_gcd", "poly_fit",
    "polynomial_fit", "lagrange", "linear_regression", "dft", "fft",
    "idft", "ifft", "convolve", "hann", "hanning", "hamming", "blackman",
    "divisors", "extended_gcd", "xgcd", "randmat", "random_matrix"
};

// ============================================================================
// 复数函数列表
// ============================================================================

// 返回复数或操作复数的函数
inline constexpr std::string_view kComplexFunctions[] = {
    "complex", "polar", "real", "imag", "arg", "conj"
};

// ============================================================================
// 检测函数
// ============================================================================

inline bool is_matrix_function(std::string_view name) {
    for (const auto& func : kMatrixFunctions) {
        if (func == name) return true;
    }
    return false;
}

inline bool is_complex_function(std::string_view name) {
    for (const auto& func : kComplexFunctions) {
        if (func == name) return true;
    }
    return false;
}

inline bool is_matrix_or_complex_function(std::string_view name) {
    return is_matrix_function(name) || is_complex_function(name);
}

#endif // PARSER_FUNCTION_CATEGORIES_H

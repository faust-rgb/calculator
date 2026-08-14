// ============================================================================
// matrix_eval_dsp.cpp - 矩阵/向量 DSP 与变换函数求值实现
// ============================================================================

#include "matrix_eval_dsp.h"
#include "matrix_dsp.h"
#include "matrix_internal.h"
#include "mymath.h"
#include <stdexcept>

namespace matrix {
namespace internal {

namespace {

Matrix make_dsp_window(std::size_t n, const std::string& name) {
    if (n == 0) {
        throw std::runtime_error(name + " requires a positive length");
    }
    Matrix result(1, n, 0.0L);
    if (n == 1) {
        result.at(0, 0) = 1.0L;
        return result;
    }

    for (std::size_t i = 0; i < n; ++i) {
        const Scalar phase =
            2.0 * mymath::kPi * static_cast<long double>(i) / static_cast<long double>(n - 1);
        if (name == "hann" || name == "hanning") {
            result.at(0, i) = 0.5 - 0.5 * mymath::cos(phase);
        } else if (name == "hamming") {
            result.at(0, i) = 0.54 - 0.46 * mymath::cos(phase);
        } else if (name == "blackman") {
            result.at(0, i) = 0.42 - 0.5 * mymath::cos(phase) +
                              0.08 * mymath::cos(2.0 * phase);
        }
    }
    return result;
}

std::size_t parse_size_arg(const std::string& expression,
                           const ScalarEvaluator& scalar_eval) {
    const Scalar value = scalar_eval(expression);
    if (!mymath::is_integer(value) || value < Scalar(0.0L)) {
        throw std::runtime_error("dimensions must be non-negative integers");
    }
    return static_cast<std::size_t>(static_cast<long double>(value) + 0.5);
}

} // namespace

bool try_evaluate_dsp_function(
    const std::string& name,
    const std::vector<std::string>& arguments,
    const ScalarEvaluator& scalar_eval,
    const std::function<Matrix(const std::string&, const std::string&)>& require_matrix,
    Value* result) {
    if (!result) return false;

    if (name == "hann" || name == "hanning" ||
        name == "hamming" || name == "blackman") {
        if (arguments.size() != 1) {
            throw std::runtime_error(name + " expects length");
        }
        *result = Value::from_matrix(make_dsp_window(
            parse_size_arg(arguments[0], scalar_eval),
            name));
        return true;
    }

    if (name == "freqz") {
        if (arguments.size() != 2 && arguments.size() != 3) {
            throw std::runtime_error("freqz expects b, a, and optional n");
        }
        std::size_t n = 512;
        if (arguments.size() == 3) {
            n = parse_size_arg(arguments[2], scalar_eval);
        }
        *result = Value::from_matrix(freqz(
            require_matrix(arguments[0], "freqz"),
            require_matrix(arguments[1], "freqz"),
            n));
        return true;
    }

    if (name == "dft" || name == "fft") {
        if (arguments.size() != 1) {
            throw std::runtime_error(name + " expects exactly one sequence argument");
        }
        *result = Value::from_matrix(complex_sequence_to_matrix(
            discrete_fourier_transform(
                as_complex_sequence(require_matrix(arguments[0], name), name),
                false),
            false));
        return true;
    }

    if (name == "idft" || name == "ifft") {
        if (arguments.size() != 1) {
            throw std::runtime_error(name + " expects exactly one sequence argument");
        }
        *result = Value::from_matrix(complex_sequence_to_matrix(
            discrete_fourier_transform(
                as_complex_sequence(require_matrix(arguments[0], name), name),
                true),
            true));
        return true;
    }

    if (name == "convolve" || name == "conv") {
        if (arguments.size() != 2) {
            throw std::runtime_error(name + " expects exactly two sequence arguments");
        }
        *result = Value::from_matrix(complex_sequence_to_matrix(
            convolve_sequences(
                as_complex_sequence(require_matrix(arguments[0], name), name),
                as_complex_sequence(require_matrix(arguments[1], name), name)),
            true));
        return true;
    }

    return false;
}

} // namespace internal
} // namespace matrix

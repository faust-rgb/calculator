#include "matrix.h"
#include "matrix_internal.h"
#include "dsp/signal_processing.h"
#include "mymath.h"
#include "mymath_complex.h"
#include "string_utils.h"
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace matrix {
namespace internal {

int& mutable_display_precision() {
    static int precision = 12;
    return precision;
}

int clamp_display_precision(int precision) {
    return std::clamp(precision, 1, 17);
}

std::string format_number(double value) {
    if (mymath::is_near_zero(value, 1e-10)) {
        value = 0.0;
    }

    std::ostringstream out;
    out << std::setprecision(mutable_display_precision()) << value;
    return out.str();
}

std::string trim_copy(const std::string& text) {
    return utils::trim_copy(text);
}

std::string format_complex(ComplexNumber value) {
    value = normalize_complex(value);
    return "complex(" + format_number(value.real) + ", " +
           format_number(value.imag) + ")";
}

long double abs_ld(long double value) {
    return mymath::abs_long_double(value);
}

long double sqrt_ld(long double value) {
    return static_cast<long double>(mymath::sqrt(static_cast<double>(value)));
}

std::size_t floor_to_size_t(double value) {
    const long long truncated = static_cast<long long>(value);
    if (value < 0.0 && static_cast<double>(truncated) != value) {
        return static_cast<std::size_t>(truncated - 1);
    }
    return static_cast<std::size_t>(truncated);
}

std::size_t ceil_to_size_t(double value) {
    const long long truncated = static_cast<long long>(value);
    if (value > 0.0 && static_cast<double>(truncated) != value) {
        return static_cast<std::size_t>(truncated + 1);
    }
    return static_cast<std::size_t>(truncated);
}

void require_same_shape(const Matrix& lhs, const Matrix& rhs, const std::string& op_name) {
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) {
        throw std::runtime_error(op_name + " requires matrices of the same shape");
    }
}

void swap_rows(Matrix* matrix, std::size_t lhs, std::size_t rhs) {
    if (lhs == rhs) {
        return;
    }
    for (std::size_t col = 0; col < matrix->cols; ++col) {
        const double temp = matrix->at(lhs, col);
        matrix->at(lhs, col) = matrix->at(rhs, col);
        matrix->at(rhs, col) = temp;
    }
}

double vector_norm_squared(const std::vector<double>& values) {
    long double sum = 0.0L;
    long double compensation = 0.0L;
    for (double value : values) {
        const long double value_ld = static_cast<long double>(value);
        const long double term = value_ld * value_ld;
        const long double adjusted = term - compensation;
        const long double next = sum + adjusted;
        compensation = (next - sum) - adjusted;
        sum = next;
    }
    return static_cast<double>(sum);
}

double max_abs_entry(const Matrix& matrix) {
    double max_value = 0.0;
    for (double value : matrix.data) {
        const double magnitude = mymath::abs(value);
        if (magnitude > max_value) {
            max_value = magnitude;
        }
    }
    return max_value;
}

double matrix_tolerance(double scale) {
    return std::max(kMatrixPivotAbsoluteEps, scale * kMatrixPivotRelativeEps);
}

double matrix_tolerance(const Matrix& matrix) {
    return matrix_tolerance(max_abs_entry(matrix));
}

std::size_t vector_length(const Matrix& matrix, const std::string& func_name) {
    if (!matrix.is_vector()) {
        throw std::runtime_error(func_name + " only accepts vectors");
    }
    return matrix.rows == 1 ? matrix.cols : matrix.rows;
}

double vector_entry(const Matrix& matrix, std::size_t index) {
    return matrix.rows == 1 ? matrix.at(0, index) : matrix.at(index, 0);
}

double off_diagonal_magnitude(const Matrix& matrix) {
    long double sum = 0.0L;
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t col = 0; col < matrix.cols; ++col) {
            if (row != col) {
                sum += static_cast<long double>(mymath::abs(matrix.at(row, col)));
            }
        }
    }
    return static_cast<double>(sum);
}

std::vector<double> matrix_column(const Matrix& matrix, std::size_t col) {
    std::vector<double> values(matrix.rows, 0.0);
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        values[row] = matrix.at(row, col);
    }
    return values;
}

void set_matrix_column(Matrix* matrix, std::size_t col, const std::vector<double>& values) {
    for (std::size_t row = 0; row < matrix->rows; ++row) {
        matrix->at(row, col) = values[row];
    }
}

double dot_vectors(const std::vector<double>& lhs, const std::vector<double>& rhs) {
    long double sum = 0.0L;
    long double compensation = 0.0L;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const long double term =
            static_cast<long double>(lhs[i]) * static_cast<long double>(rhs[i]);
        const long double adjusted = term - compensation;
        const long double next = sum + adjusted;
        compensation = (next - sum) - adjusted;
        sum = next;
    }
    return static_cast<double>(sum);
}

bool orthonormalize(std::vector<double>* values,
                    const std::vector<std::vector<double>>& basis) {
    for (const std::vector<double>& existing : basis) {
        const double projection = dot_vectors(*values, existing);
        for (std::size_t i = 0; i < values->size(); ++i) {
            (*values)[i] -= projection * existing[i];
        }
    }

    const double magnitude = mymath::sqrt(vector_norm_squared(*values));
    if (magnitude <= matrix_tolerance(magnitude)) {
        return false;
    }
    for (double& value : *values) {
        value /= magnitude;
    }
    return true;
}

std::vector<double> standard_basis_vector(std::size_t size, std::size_t index) {
    std::vector<double> values(size, 0.0);
    values[index] = 1.0;
    return values;
}

std::vector<double> as_vector_values(const Matrix& matrix, const std::string& func_name) {
    const std::size_t length = vector_length(matrix, func_name);
    std::vector<double> values(length, 0.0);
    for (std::size_t i = 0; i < length; ++i) {
        values[i] = vector_entry(matrix, i);
    }
    return values;
}

void require_nonempty_values(const std::vector<double>& values,
                             const std::string& func_name) {
    if (values.empty()) {
        throw std::runtime_error(func_name + " requires at least one value");
    }
}

std::vector<double> sort_values(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    return values;
}

std::pair<double, double> linear_regression_fit(const std::vector<double>& x,
                                                const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("linreg requires sample vectors of the same non-zero length");
    }

    long double sum_x = 0.0L;
    long double sum_y = 0.0L;
    long double sum_xx = 0.0L;
    long double sum_xy = 0.0L;
    for (std::size_t i = 0; i < x.size(); ++i) {
        sum_x += static_cast<long double>(x[i]);
        sum_y += static_cast<long double>(y[i]);
        sum_xx += static_cast<long double>(x[i]) * static_cast<long double>(x[i]);
        sum_xy += static_cast<long double>(x[i]) * static_cast<long double>(y[i]);
    }

    const long double n = static_cast<long double>(x.size());
    const long double denominator = n * sum_xx - sum_x * sum_x;
    if (mymath::is_near_zero(static_cast<double>(denominator), 1e-12)) {
        throw std::runtime_error("linreg requires non-constant x values");
    }

    const long double slope = (n * sum_xy - sum_x * sum_y) / denominator;
    const long double intercept = (sum_y - slope * sum_x) / n;
    return {static_cast<double>(slope), static_cast<double>(intercept)};
}

ComplexNumber normalize_complex(ComplexNumber value) {
    if (mymath::is_near_zero(value.real, 1e-12)) {
        value.real = 0.0;
    }
    if (mymath::is_near_zero(value.imag, 1e-12)) {
        value.imag = 0.0;
    }
    return value;
}

bool is_complex_vector(const Matrix& matrix) {
    return (matrix.rows == 1 && matrix.cols == 2) ||
           (matrix.rows == 2 && matrix.cols == 1);
}

double complex_real(const Matrix& matrix) {
    return matrix.at(0, 0);
}

double complex_imag(const Matrix& matrix) {
    return matrix.rows == 1 ? matrix.at(0, 1) : matrix.at(1, 0);
}

Matrix complex_value(double real, double imag) {
    Matrix result(1, 2, 0.0);
    result.at(0, 0) = real;
    result.at(0, 1) = imag;
    return result;
}

ComplexNumber complex_from_matrix(const Matrix& matrix) {
    if (!is_complex_vector(matrix)) {
        throw std::runtime_error("complex value requires a 2-element vector");
    }
    return normalize_complex({complex_real(matrix), complex_imag(matrix)});
}

bool try_complex_from_value(const Value& value, ComplexNumber* complex) {
    if (complex == nullptr) {
        return false;
    }
    if (value.is_complex) {
        *complex = normalize_complex(value.complex);
        return true;
    }
    if (value.is_matrix) {
        if (is_complex_vector(value.matrix)) {
            *complex = complex_from_matrix(value.matrix);
            return true;
        }
        if (value.matrix.rows == 1 && value.matrix.cols == 1) {
            *complex = {value.matrix.at(0, 0), 0.0};
            return true;
        }
        return false;
    }
    *complex = {value.scalar, 0.0};
    return true;
}

std::vector<ComplexSample> as_complex_sequence(const Matrix& matrix,
                                               const std::string& func_name) {
    std::vector<ComplexSample> result;
    if (matrix.is_vector()) {
        const std::size_t length = vector_length(matrix, func_name);
        result.reserve(length);
        for (std::size_t i = 0; i < length; ++i) {
            result.push_back({vector_entry(matrix, i), 0.0});
        }
        return result;
    }

    if (matrix.cols == 2) {
        result.reserve(matrix.rows);
        for (std::size_t row = 0; row < matrix.rows; ++row) {
            const ComplexNumber value =
                normalize_complex({matrix.at(row, 0), matrix.at(row, 1)});
            result.push_back({value.real, value.imag});
        }
        return result;
    }

    throw std::runtime_error(func_name + " requires a real vector or Nx2 complex sample matrix");
}

Matrix complex_sequence_to_matrix(const std::vector<ComplexSample>& values,
                                  bool prefer_real_vector) {
    bool all_real = prefer_real_vector;
    for (const ComplexSample& value : values) {
        if (!mymath::is_near_zero(value.imag, 1e-12)) {
            all_real = false;
            break;
        }
    }

    if (all_real) {
        Matrix result(1, values.size(), 0.0);
        for (std::size_t i = 0; i < values.size(); ++i) {
            result.at(0, i) = values[i].real;
        }
        return result;
    }

    Matrix result(values.size(), 2, 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const ComplexNumber value = normalize_complex({values[i].real, values[i].imag});
        result.at(i, 0) = value.real;
        result.at(i, 1) = value.imag;
    }
    return result;
}

std::vector<ComplexSample> discrete_fourier_transform(
    const std::vector<ComplexSample>& input,
    bool inverse) {
    std::vector<signal::Complex> values;
    values.reserve(input.size());
    for (const ComplexSample& sample : input) {
        values.emplace_back(sample.real, sample.imag);
    }

    const std::vector<signal::Complex> transformed =
        inverse ? signal::ifft(values) : signal::fft(values);

    std::vector<ComplexSample> result;
    result.reserve(transformed.size());
    for (const signal::Complex& value : transformed) {
        const ComplexNumber normalized = normalize_complex({value.real(), value.imag()});
        result.push_back({normalized.real, normalized.imag});
    }
    return result;
}

std::vector<ComplexSample> convolve_sequences(
    const std::vector<ComplexSample>& lhs,
    const std::vector<ComplexSample>& rhs) {
    if (lhs.empty() || rhs.empty()) {
        return {};
    }

    std::vector<ComplexSample> result(lhs.size() + rhs.size() - 1, {0.0, 0.0});
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        const mymath::complex<double> a(lhs[i].real, lhs[i].imag);
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            const mymath::complex<double> b(rhs[j].real, rhs[j].imag);
            const mymath::complex<double> product = a * b;
            result[i + j].real += product.real();
            result[i + j].imag += product.imag();
        }
    }

    for (ComplexSample& sample : result) {
        const ComplexNumber normalized = normalize_complex({sample.real, sample.imag});
        sample.real = normalized.real;
        sample.imag = normalized.imag;
    }
    return result;
}

} // namespace internal

void set_display_precision(int precision) {
    internal::mutable_display_precision() = internal::clamp_display_precision(precision);
}

} // namespace matrix

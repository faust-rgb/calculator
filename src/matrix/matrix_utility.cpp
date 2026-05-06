#include "matrix.h"
#include "matrix_internal.h"
#include "dsp/signal_processing.h"
#include "mymath.h"
#include "mymath_complex.h"
#include "string_utils.h"
#include "precise/precise_decimal.h"
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

template <typename T>
T t_abs(T v) {
    if constexpr (std::is_same_v<T, PreciseDecimal>) return precise::abs(v);
    else return static_cast<T>(mymath::abs(static_cast<double>(v)));
}

template <typename T>
T t_sqrt(T v) {
    if constexpr (std::is_same_v<T, PreciseDecimal>) return precise::sqrt(v);
    else return static_cast<T>(mymath::sqrt(static_cast<double>(v)));
}

template <typename T>
std::string format_number(T value) {
    if constexpr (std::is_same_v<T, PreciseDecimal>) return value.to_string();
    else {
        if (mymath::is_near_zero(static_cast<double>(value), 1e-10)) value = T(0);
        std::ostringstream out;
        out << std::setprecision(mutable_display_precision()) << static_cast<double>(value);
        return out.str();
    }
}

template <typename T>
std::string format_complex(TComplex<T> value) {
    value = normalize_complex(value);
    return "complex(" + format_number(value.real) + ", " + format_number(value.imag) + ")";
}

std::string trim_copy(const std::string& text) { return utils::trim_copy(text); }

template <typename T>
void require_same_shape(const TMatrix<T>& lhs, const TMatrix<T>& rhs, const std::string& op_name) {
    if (lhs.rows != rhs.rows || lhs.cols != rhs.cols) throw std::runtime_error(op_name + " requires same shape");
}

template <typename T>
void swap_rows(TMatrix<T>* m, std::size_t lhs, std::size_t rhs) {
    if (lhs == rhs) return;
    for (std::size_t j = 0; j < m->cols; ++j) std::swap(m->at(lhs, j), m->at(rhs, j));
}

template <typename T>
T vector_norm_squared(const std::vector<T>& values) {
    T sum = T(0);
    for (const T& v : values) sum += v * v;
    return sum;
}

template <>
double vector_norm_squared<double>(const std::vector<double>& values) {
    long double sum = 0.0L, compensation = 0.0L;
    for (double v : values) {
        long double val = static_cast<long double>(v);
        long double term = val * val - compensation;
        long double next = sum + term;
        compensation = (next - sum) - term;
        sum = next;
    }
    return static_cast<double>(sum);
}

template <typename T>
T max_abs_entry(const TMatrix<T>& m) {
    T max_v = T(0);
    for (const T& v : m.data) {
        T mag = t_abs(v);
        if (mag > max_v) max_v = mag;
    }
    return max_v;
}

template <typename T>
T matrix_tolerance(T scale) {
    if constexpr (std::is_same_v<T, PreciseDecimal>) return scale * PreciseDecimal("1e-35");
    else return std::max(kMatrixPivotAbsoluteEps, scale * kMatrixPivotRelativeEps);
}

template <typename T>
T matrix_tolerance(const TMatrix<T>& m) { return matrix_tolerance(max_abs_entry(m)); }

template <typename T>
std::size_t vector_length(const TMatrix<T>& m, const std::string& name) {
    if (!m.is_vector()) throw std::runtime_error(name + " requires a vector");
    return m.rows == 1 ? m.cols : m.rows;
}

template <typename T>
T vector_entry(const TMatrix<T>& m, std::size_t i) {
    return m.rows == 1 ? m.at(0, i) : m.at(i, 0);
}

template <typename T>
T off_diagonal_magnitude(const TMatrix<T>& m) {
    T sum = T(0);
    for (std::size_t i = 0; i < m.rows; ++i) {
        for (std::size_t j = 0; j < m.cols; ++j) {
            if (i != j) sum += t_abs(m.at(i, j));
        }
    }
    return sum;
}

template <typename T>
std::vector<T> matrix_column(const TMatrix<T>& m, std::size_t j) {
    std::vector<T> col(m.rows);
    for (std::size_t i = 0; i < m.rows; ++i) col[i] = m.at(i, j);
    return col;
}

template <typename T>
void set_matrix_column(TMatrix<T>* m, std::size_t j, const std::vector<T>& values) {
    for (std::size_t i = 0; i < m->rows; ++i) m->at(i, j) = values[i];
}

template <typename T>
bool orthonormalize(std::vector<T>* v, const std::vector<std::vector<T>>& basis) {
    for (const auto& b : basis) {
        T proj = T(0);
        for (std::size_t i = 0; i < v->size(); ++i) proj += (*v)[i] * b[i];
        for (std::size_t i = 0; i < v->size(); ++i) (*v)[i] -= proj * b[i];
    }
    T mag = t_sqrt(vector_norm_squared(*v));
    if (mag <= matrix_tolerance(mag)) return false;
    for (T& x : *v) x /= mag;
    return true;
}

template <typename T>
std::vector<T> standard_basis_vector(std::size_t n, std::size_t i) {
    std::vector<T> v(n, T(0)); v[i] = T(1); return v;
}

template <typename T>
std::vector<T> as_vector_values(const TMatrix<T>& m, const std::string& name) {
    std::size_t n = vector_length(m, name);
    std::vector<T> v(n);
    for (std::size_t i = 0; i < n; ++i) v[i] = vector_entry(m, i);
    return v;
}

template <typename T>
std::pair<T, T> linear_regression_fit(const std::vector<T>& x, const std::vector<T>& y) {
    if (x.size() != y.size() || x.empty()) throw std::runtime_error("invalid data for linreg");
    T sum_x = T(0), sum_y = T(0), sum_xx = T(0), sum_xy = T(0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        sum_x += x[i]; sum_y += y[i];
        sum_xx += x[i] * x[i]; sum_xy += x[i] * y[i];
    }
    T n = T(static_cast<long long>(x.size()));
    T den = n * sum_xx - sum_x * sum_x;
    if (t_abs(den) < matrix_tolerance(T(1))) throw std::runtime_error("linreg requires non-constant x");
    T slope = (n * sum_xy - sum_x * sum_y) / den;
    T intercept = (sum_y - slope * sum_x) / n;
    return {slope, intercept};
}

template <typename T>
TComplex<T> normalize_complex(TComplex<T> v) {
    T eps = matrix_tolerance(T(1)) * T(100);
    if (t_abs(v.real) < eps) v.real = T(0);
    if (t_abs(v.imag) < eps) v.imag = T(0);
    return v;
}

template <typename T> bool is_complex_vector(const TMatrix<T>& m) { return (m.rows==1 && m.cols==2) || (m.rows==2 && m.cols==1); }
template <typename T> T complex_real(const TMatrix<T>& m) { return m.at(0,0); }
template <typename T> T complex_imag(const TMatrix<T>& m) { return m.rows==1 ? m.at(0,1) : m.at(1,0); }
template <typename T> TMatrix<T> complex_value(T re, T im) { TMatrix<T> r(1,2); r.at(0,0)=re; r.at(0,1)=im; return r; }

template <typename T>
TComplex<T> complex_from_matrix(const TMatrix<T>& m) {
    if (!is_complex_vector(m)) throw std::runtime_error("2-element vector required");
    return normalize_complex<T>({complex_real(m), complex_imag(m)});
}

template <typename T>
bool try_complex_from_value(const TValue<T>& val, TComplex<T>* c) {
    if (!c) return false;
    if (val.is_complex) { *c = normalize_complex(val.complex); return true; }
    if (val.is_matrix) {
        if (is_complex_vector(val.matrix)) { *c = complex_from_matrix(val.matrix); return true; }
        if (val.matrix.rows==1 && val.matrix.cols==1) { *c = {val.matrix.at(0,0), T(0)}; return true; }
        return false;
    }
    *c = {val.scalar, T(0)}; return true;
}

template <typename T>
std::vector<TComplexSample<T>> as_complex_sequence(const TMatrix<T>& m, const std::string& name) {
    std::vector<TComplexSample<T>> res;
    if (m.is_vector()) {
        std::size_t n = vector_length(m, name);
        for (std::size_t i = 0; i < n; ++i) res.push_back({vector_entry(m, i), T(0)});
        return res;
    }
    if (m.cols == 2) {
        for (std::size_t i = 0; i < m.rows; ++i) res.push_back({m.at(i,0), m.at(i,1)});
        return res;
    }
    throw std::runtime_error(name + " requires vector or Nx2 matrix");
}

template <typename T>
TMatrix<T> complex_sequence_to_matrix(const std::vector<TComplexSample<T>>& values, bool prefer_real) {
    bool all_real = prefer_real;
    T eps = matrix_tolerance(T(1)) * T(100);
    for (const auto& v : values) if (t_abs(v.imag) > eps) { all_real = false; break; }
    if (all_real) {
        TMatrix<T> r(1, values.size());
        for (std::size_t i = 0; i < values.size(); ++i) r.at(0,i) = values[i].real;
        return r;
    }
    TMatrix<T> r(values.size(), 2);
    for (std::size_t i = 0; i < values.size(); ++i) { r.at(i,0) = values[i].real; r.at(i,1) = values[i].imag; }
    return r;
}

template <typename T>
std::vector<TComplexSample<T>> discrete_fourier_transform(const std::vector<TComplexSample<T>>& input, bool inverse) {
    if constexpr (std::is_same_v<T, double>) {
        std::vector<signal::Complex> v;
        for (const auto& s : input) v.emplace_back(s.real, s.imag);
        auto trans = inverse ? signal::ifft(v) : signal::fft(v);
        std::vector<TComplexSample<T>> res;
        for (const auto& val : trans) res.push_back({val.real(), val.imag()});
        return res;
    } else {
        throw std::runtime_error("DFT not supported for this type yet");
    }
}

template <typename T>
std::vector<TComplexSample<T>> convolve_sequences(const std::vector<TComplexSample<T>>& lhs, const std::vector<TComplexSample<T>>& rhs) {
    if (lhs.empty() || rhs.empty()) return {};
    std::vector<TComplexSample<T>> res(lhs.size() + rhs.size() - 1, {T(0), T(0)});
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
            T ac = lhs[i].real * rhs[j].real;
            T bd = lhs[i].imag * rhs[j].imag;
            T ad = lhs[i].real * rhs[j].imag;
            T bc = lhs[i].imag * rhs[j].real;
            res[i+j].real += (ac - bd);
            res[i+j].imag += (ad + bc);
        }
    }
    return res;
}

} // namespace internal

void set_display_precision(int precision) {
    internal::mutable_display_precision() = internal::clamp_display_precision(precision);
}

} // namespace matrix

#define INSTANTIATE_UTIL(TYPE) \
    template TYPE matrix::internal::t_abs(TYPE); \
    template TYPE matrix::internal::t_sqrt(TYPE); \
    template std::string matrix::internal::format_number(TYPE); \
    template std::string matrix::internal::format_complex(matrix::TComplex<TYPE>); \
    template void matrix::internal::require_same_shape(const matrix::TMatrix<TYPE>&, const matrix::TMatrix<TYPE>&, const std::string&); \
    template void matrix::internal::swap_rows(matrix::TMatrix<TYPE>*, std::size_t, std::size_t); \
    template TYPE matrix::internal::vector_norm_squared(const std::vector<TYPE>&); \
    template TYPE matrix::internal::max_abs_entry(const matrix::TMatrix<TYPE>&); \
    template TYPE matrix::internal::matrix_tolerance(TYPE); \
    template TYPE matrix::internal::matrix_tolerance(const matrix::TMatrix<TYPE>&); \
    template std::size_t matrix::internal::vector_length(const matrix::TMatrix<TYPE>&, const std::string&); \
    template TYPE matrix::internal::vector_entry(const matrix::TMatrix<TYPE>&, std::size_t); \
    template TYPE matrix::internal::off_diagonal_magnitude(const matrix::TMatrix<TYPE>&); \
    template std::vector<TYPE> matrix::internal::matrix_column(const matrix::TMatrix<TYPE>&, std::size_t); \
    template void matrix::internal::set_matrix_column(matrix::TMatrix<TYPE>*, std::size_t, const std::vector<TYPE>&); \
    template bool matrix::internal::orthonormalize(std::vector<TYPE>*, const std::vector<std::vector<TYPE>>&); \
    template std::vector<TYPE> matrix::internal::standard_basis_vector(std::size_t, std::size_t); \
    template std::vector<TYPE> matrix::internal::as_vector_values(const matrix::TMatrix<TYPE>&, const std::string&); \
    template std::pair<TYPE, TYPE> matrix::internal::linear_regression_fit(const std::vector<TYPE>&, const std::vector<TYPE>&); \
    template matrix::TComplex<TYPE> matrix::internal::normalize_complex(matrix::TComplex<TYPE>); \
    template bool matrix::internal::is_complex_vector(const matrix::TMatrix<TYPE>&); \
    template TYPE matrix::internal::complex_real(const matrix::TMatrix<TYPE>&); \
    template TYPE matrix::internal::complex_imag(const matrix::TMatrix<TYPE>&); \
    template matrix::TMatrix<TYPE> matrix::internal::complex_value(TYPE, TYPE); \
    template matrix::TComplex<TYPE> matrix::internal::complex_from_matrix(const matrix::TMatrix<TYPE>&); \
    template bool matrix::internal::try_complex_from_value(const matrix::TValue<TYPE>&, matrix::TComplex<TYPE>*); \
    template std::vector<matrix::internal::TComplexSample<TYPE>> matrix::internal::as_complex_sequence(const matrix::TMatrix<TYPE>&, const std::string&); \
    template matrix::TMatrix<TYPE> matrix::internal::complex_sequence_to_matrix(const std::vector<matrix::internal::TComplexSample<TYPE>>&, bool); \
    template std::vector<matrix::internal::TComplexSample<TYPE>> matrix::internal::discrete_fourier_transform(const std::vector<matrix::internal::TComplexSample<TYPE>>&, bool); \
    template std::vector<matrix::internal::TComplexSample<TYPE>> matrix::internal::convolve_sequences(const std::vector<matrix::internal::TComplexSample<TYPE>>&, const std::vector<matrix::internal::TComplexSample<TYPE>>&);

INSTANTIATE_UTIL(double)
INSTANTIATE_UTIL(PreciseDecimal)

/**
 * @file matrix_module.cpp
 * @brief 矩阵模块实现
 *
 * 所有函数通过统一的 StoredValue 接口注册，
 * 不再依赖 ScalarEvaluator/ValueFunction 遗留类型。
 */

#include "matrix_module.h"
#include "core/services/service_locator.h"
#include "core/services/core_manager_interfaces.h"
#include "matrix.h"
#include "matrix_internal.h"
#include "polynomial/polynomial.h"
#include "matrix_dsp.h"
#include "mymath.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "core/common/calculator_exceptions.h"
#include "math/functions/integer/integer_helpers.h"
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <random>
#include <cmath>
#include <tuple>

namespace {

using namespace matrix;
using Scalar = mymath::Scalar;
using ComplexNumber = mymath::complex<Scalar>;

Matrix require_matrix(const StoredValue& val, const std::string& func_name) {
    if (!val.is_matrix || !val.matrix_ptr) {
        throw std::runtime_error(func_name + " expects a matrix argument");
    }
    return *val.matrix_ptr;
}

bool try_complex_from_stored(const StoredValue& v, ComplexNumber* z) {
    if (v.is_complex) {
        *z = v.complex;
        return true;
    }
    if (!v.is_matrix) {
        z->real(v.decimal);
        z->imag(0.0L);
        return true;
    }
    return false;
}

ComplexNumber require_complex_argument(const StoredValue& val, const std::string& func_name) {
    ComplexNumber z;
    if (try_complex_from_stored(val, &z)) return z;
    throw std::runtime_error(func_name + " expects a scalar or complex argument");
}

std::size_t parse_index_argument(const StoredValue& val, const std::string& func_name) {
    if (val.is_matrix || val.is_complex) {
        throw std::runtime_error(func_name + " requires non-negative integer index");
    }
    const Scalar dec = val.get_decimal();
    if (!mymath::is_integer(dec) || dec < Scalar(0.0L)) {
        throw std::runtime_error(func_name + " requires non-negative integer index");
    }
    return static_cast<std::size_t>(static_cast<long double>(dec) + 0.5);
}

StoredValue make_matrix_result(Matrix m) {
    StoredValue res;
    res.is_matrix = true;
    res.matrix_ptr = std::make_shared<Matrix>(std::move(m));
    return res;
}

StoredValue make_scalar_result(Scalar s) {
    return StoredValue(s);
}

StoredValue make_complex_result(ComplexNumber z) {
    return StoredValue(z);
}

std::vector<Scalar> require_vector(const StoredValue& value, const std::string& name) {
    const Matrix matrix = require_matrix(value, name);
    std::vector<Scalar> result;
    result.reserve(matrix.data.size());
    for (const Scalar entry : matrix.data) result.push_back(entry);
    return result;
}

std::size_t require_positive_length(const StoredValue& value, const std::string& name) {
    if (value.is_matrix || value.is_complex || !mymath::is_integer(value.get_decimal()) ||
        value.get_decimal() <= Scalar(0)) {
        throw std::runtime_error(name + " expects a positive integer length");
    }
    return static_cast<std::size_t>(static_cast<long double>(value.get_decimal()) + 0.5L);
}

StoredValue make_dsp_window(const std::vector<StoredValue>& args, const std::string& name) {
    if (args.size() != 1) throw std::runtime_error(name + " expects length");
    const std::size_t n = require_positive_length(args[0], name);
    Matrix result(1, n, 0.0L);
    if (n == 1) {
        result.at(0, 0) = Scalar(1);
        return make_matrix_result(std::move(result));
    }
    for (std::size_t i = 0; i < n; ++i) {
        const Scalar phase = Scalar(2) * mymath::kPi * Scalar(i) / Scalar(n - 1);
        if (name == "hann" || name == "hanning") {
            result.at(0, i) = Scalar(0.5L) - Scalar(0.5L) * mymath::cos(phase);
        } else if (name == "hamming") {
            result.at(0, i) = Scalar(0.54L) - Scalar(0.46L) * mymath::cos(phase);
        } else {
            result.at(0, i) = Scalar(0.42L) - Scalar(0.5L) * mymath::cos(phase) +
                               Scalar(0.08L) * mymath::cos(Scalar(2) * phase);
        }
    }
    return make_matrix_result(std::move(result));
}

StoredValue make_transform_result(const std::vector<StoredValue>& args,
                                  const std::string& name, bool inverse) {
    if (args.size() != 1) throw std::runtime_error(name + " expects one sequence argument");
    const Matrix input = require_matrix(args[0], name);
    return make_matrix_result(matrix::internal::complex_sequence_to_matrix(
        matrix::internal::discrete_fourier_transform(
            matrix::internal::as_complex_sequence(input, name), inverse),
        inverse));
}

StoredValue make_convolution_result(const std::vector<StoredValue>& args,
                                    const std::string& name) {
    if (args.size() != 2) throw std::runtime_error(name + " expects two sequence arguments");
    return make_matrix_result(matrix::internal::complex_sequence_to_matrix(
        matrix::internal::convolve_sequences(
            matrix::internal::as_complex_sequence(require_matrix(args[0], name), name),
            matrix::internal::as_complex_sequence(require_matrix(args[1], name), name)),
        true));
}

std::string format_eigenvalue_matrix(const Matrix& values) {
    if (values.rows <= 1 || values.cols != 2) {
        return values.to_string();
    }
    std::ostringstream out;
    out << "[";
    for (std::size_t row = 0; row < values.rows; ++row) {
        if (row != 0) out << ", ";
        out << matrix::internal::format_complex<Scalar>({values.at(row, 0), values.at(row, 1)});
    }
    out << "]";
    return out.str();
}

} // namespace

// ============================================================================
// 命令接口实现
// ============================================================================

std::vector<std::string> MatrixModule::get_commands() const {
    return {"eig", "svd", "lu_p"};
}

std::string MatrixModule::execute_args_view(std::string_view command,
                                            const std::vector<std::string_view>& args,
                                            ServiceLocator& locator) {
    using namespace module_helpers;
    require_args_count(args, 1, 1, command);

    auto& services = *locator.resolve<CoreServices>();
    const StoredValue val = services.evaluation.evaluate_value(std::string(args[0]), false);
    const Matrix matrix_value = require_matrix(val, std::string(command));

    if (command == "svd") {
        return "U: " + svd_u<Scalar>(matrix_value).to_string() +
               "\nS: " + svd_s<Scalar>(matrix_value).to_string() +
               "\nVt: " + svd_vt<Scalar>(matrix_value).to_string();
    }

    if (command == "lu_p") {
        return lu_p<Scalar>(matrix_value).to_string();
    }

    // eig
    try {
        const Matrix values = eigenvalues<Scalar>(matrix_value);
        if (values.rows > 1 && values.cols == 2) {
            return "values: " + format_eigenvalue_matrix(values) +
                   "\nvectors: unavailable for complex eigenvalues";
        }
        return "values: " + values.to_string() +
               "\nvectors: " + eigenvectors<Scalar>(matrix_value).to_string();
    } catch (const std::exception&) {
        if (matrix_value.rows == 2 && matrix_value.cols == 2) {
            const Scalar tr = matrix_value.at(0, 0) + matrix_value.at(1, 1);
            const Scalar det = determinant<Scalar>(matrix_value);
            const Scalar discriminant = tr * tr - Scalar(4.0L) * det;
            if (discriminant < Scalar(0.0L)) {
                const Scalar real = tr * Scalar(0.5L);
                const Scalar imag = mymath::sqrt(-discriminant) * Scalar(0.5L);
                std::ostringstream out;
                out << "values: [complex(" << format_decimal(real) << ", "
                    << format_decimal(imag) << "), complex("
                    << format_decimal(real) << ", " << format_decimal(-imag)
                    << ")]\nvectors: unavailable for complex eigenvalues";
                return out.str();
            }
        }
        throw;
    }
}

// ============================================================================
// 函数注册 — 全部使用 StoredValue 统一接口
// ============================================================================

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
MatrixModule::get_functions_map() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    // Matrix construction functions are registered here so the unified AST
    // does not need the legacy matrix parser for ordinary matrix expressions.
    funcs["vec"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty()) throw std::runtime_error("vec expects at least one element");
        if (args.size() == 1 && args[0].is_matrix) {
            return make_matrix_result(vectorize(require_matrix(args[0], "vec")));
        }
        std::vector<Scalar> values;
        values.reserve(args.size());
        for (const auto& arg : args) {
            if (!arg.is_scalar()) throw std::runtime_error("vec elements must be scalar");
            values.push_back(arg.get_decimal());
        }
        return make_matrix_result(Matrix::vector(values));
    };
    funcs["mat"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2) throw std::runtime_error("mat expects rows, cols, and optional elements");
        const auto size_arg = [](const StoredValue& value) -> std::size_t {
            const Scalar n = value.get_decimal();
            if (!mymath::is_integer(n) || n < 0) throw std::runtime_error("matrix dimensions must be non-negative integers");
            return static_cast<std::size_t>(round_to_long_long(n));
        };
        const std::size_t rows = size_arg(args[0]);
        const std::size_t cols = size_arg(args[1]);
        if (args.size() != rows * cols + 2) throw std::runtime_error("mat element count does not match the requested shape");
        Matrix result(rows, cols, 0.0L);
        for (std::size_t i = 0; i < rows * cols; ++i) {
            if (!args[i + 2].is_scalar()) throw std::runtime_error("mat elements must be scalar");
            result.data[i] = args[i + 2].get_decimal();
        }
        return make_matrix_result(std::move(result));
    };
    funcs["zeros"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("zeros expects exactly two arguments");
        return make_matrix_result(Matrix::zero(parse_index_argument(args[0], "zeros"),
                                               parse_index_argument(args[1], "zeros")));
    };
    funcs["eye"] = funcs["identity"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("eye expects exactly one argument");
        return make_matrix_result(Matrix::identity(static_cast<std::size_t>(round_to_long_long(args[0].get_decimal()))));
    };

    funcs["hann"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, "hann"); };
    funcs["hanning"] = funcs["hann"];
    funcs["hamming"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, "hamming"); };
    funcs["blackman"] = [](const std::vector<StoredValue>& args) { return make_dsp_window(args, "blackman"); };
    funcs["dft"] = [](const std::vector<StoredValue>& args) { return make_transform_result(args, "dft", false); };
    funcs["fft"] = funcs["dft"];
    funcs["idft"] = [](const std::vector<StoredValue>& args) { return make_transform_result(args, "idft", true); };
    funcs["ifft"] = funcs["idft"];
    funcs["convolve"] = [](const std::vector<StoredValue>& args) { return make_convolution_result(args, "convolve"); };
    funcs["conv"] = funcs["convolve"];

    funcs["divisors"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1 || args[0].is_matrix || args[0].is_complex) {
            throw std::runtime_error("divisors expects one integer");
        }
        const Scalar raw = args[0].get_decimal();
        if (!mymath::is_integer(raw) || raw <= Scalar(0)) throw std::runtime_error("divisors expects a positive integer");
        const long long value = static_cast<long long>(raw);
        std::vector<Scalar> values;
        for (long long d = 1; d <= value / d; ++d) {
            if (value % d != 0) continue;
            values.push_back(Scalar(d));
            if (d != value / d) values.push_back(Scalar(value / d));
        }
        std::sort(values.begin(), values.end());
        return make_matrix_result(Matrix::vector(std::move(values)));
    };
    funcs["extended_gcd"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2 || args[0].is_matrix || args[1].is_matrix ||
            args[0].is_complex || args[1].is_complex) {
            throw std::runtime_error("extended_gcd expects two integers");
        }
        if (!mymath::is_integer(args[0].get_decimal()) || !mymath::is_integer(args[1].get_decimal())) {
            throw std::runtime_error("extended_gcd expects two integers");
        }
        long long a = static_cast<long long>(args[0].get_decimal());
        long long b = static_cast<long long>(args[1].get_decimal());
        long long old_r = a, r = b, old_s = 1, s = 0, old_t = 0, t = 1;
        while (r != 0) {
            const long long q = old_r / r;
            std::tie(old_r, r) = std::make_pair(r, old_r - q * r);
            std::tie(old_s, s) = std::make_pair(s, old_s - q * s);
            std::tie(old_t, t) = std::make_pair(t, old_t - q * t);
        }
        return make_matrix_result(Matrix::vector({Scalar(old_r), Scalar(old_s), Scalar(old_t)}));
    };
    funcs["xgcd"] = funcs["extended_gcd"];

    // --- Matrix-only functions (1 matrix arg → matrix result) ---
    auto mat_func_1 = [](auto f, const std::string& name) {
        return [f, name](const std::vector<StoredValue>& args) -> StoredValue {
            if (args.size() != 1) throw std::runtime_error(name + " expects 1 argument");
            return make_matrix_result(f(require_matrix(args[0], name)));
        };
    };

    funcs["transpose"] = mat_func_1(transpose<Scalar>, "transpose");
    funcs["inverse"] = mat_func_1(inverse<Scalar>, "inverse");
    funcs["inv"] = funcs["inverse"];
    funcs["pinv"] = mat_func_1(pseudo_inverse<Scalar>, "pinv");
    funcs["null"] = mat_func_1(nullspace<Scalar>, "null");
    funcs["qr_q"] = mat_func_1(qr_q<Scalar>, "qr_q");
    funcs["qr_r"] = mat_func_1(qr_r<Scalar>, "qr_r");
    funcs["lu_l"] = mat_func_1(lu_l<Scalar>, "lu_l");
    funcs["lu_u"] = mat_func_1(lu_u<Scalar>, "lu_u");
    funcs["lu_p"] = mat_func_1(lu_p<Scalar>, "lu_p");
    funcs["svd_u"] = mat_func_1(svd_u<Scalar>, "svd_u");
    funcs["svd_s"] = mat_func_1(svd_s<Scalar>, "svd_s");
    funcs["svd_vt"] = mat_func_1(svd_vt<Scalar>, "svd_vt");
    funcs["cholesky"] = mat_func_1(cholesky<Scalar>, "cholesky");
    funcs["hessenberg"] = mat_func_1(hessenberg<Scalar>, "hessenberg");
    funcs["schur"] = mat_func_1(schur<Scalar>, "schur");
    funcs["diag"] = mat_func_1(diag<Scalar>, "diag");
    funcs["rref"] = mat_func_1(rref<Scalar>, "rref");
    funcs["eigvals"] = mat_func_1(eigenvalues<Scalar>, "eigvals");
    funcs["expm"] = mat_func_1(matrix_exponential<Scalar>, "expm");
    funcs["eigvecs"] = mat_func_1(eigenvectors<Scalar>, "eigvecs");

    funcs["resize"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 3) throw std::runtime_error("resize expects matrix, rows, cols");
        Matrix result = require_matrix(args[0], "resize");
        const auto rows = parse_index_argument(args[1], "resize");
        const auto cols = parse_index_argument(args[2], "resize");
        result.resize(rows, cols, Scalar(0.0L));
        return make_matrix_result(std::move(result));
    };
    funcs["append_row"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2) throw std::runtime_error("append_row expects a matrix and values");
        Matrix result = require_matrix(args[0], "append_row");
        std::vector<Scalar> values;
        for (std::size_t i = 1; i < args.size(); ++i) {
            if (!args[i].is_scalar()) throw std::runtime_error("append_row values must be scalar");
            values.push_back(args[i].get_decimal());
        }
        result.append_row(values);
        return make_matrix_result(std::move(result));
    };
    funcs["append_col"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < 2) throw std::runtime_error("append_col expects a matrix and values");
        Matrix result = require_matrix(args[0], "append_col");
        std::vector<Scalar> values;
        for (std::size_t i = 1; i < args.size(); ++i) {
            if (!args[i].is_scalar()) throw std::runtime_error("append_col values must be scalar");
            values.push_back(args[i].get_decimal());
        }
        result.append_col(values);
        return make_matrix_result(std::move(result));
    };
    funcs["randmat"] = funcs["random_matrix"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2 && args.size() != 4) throw std::runtime_error("randmat expects rows, cols, optional min, max");
        const auto rows = parse_index_argument(args[0], "randmat");
        const auto cols = parse_index_argument(args[1], "randmat");
        const Scalar lo = args.size() == 4 ? args[2].get_decimal() : Scalar(0.0L);
        const Scalar hi = args.size() == 4 ? args[3].get_decimal() : Scalar(1.0L);
        if (lo > hi) throw std::runtime_error("randmat requires min <= max");
        std::mt19937_64 engine(std::random_device{}());
        std::uniform_real_distribution<long double> distribution(static_cast<long double>(lo), static_cast<long double>(hi));
        Matrix result(rows, cols);
        for (auto& cell : result.data) cell = Scalar(distribution(engine));
        return make_matrix_result(std::move(result));
    };

    funcs["charpoly"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.empty() || args.size() > 2) throw std::runtime_error("charpoly expects 1 or 2 arguments");
        Matrix m = require_matrix(args[0], "charpoly");
        Matrix poly = characteristic_polynomial<Scalar>(m);
        if (args.size() == 2 && args[1].is_string) {
            std::string var = args[1].string_value;
            std::string expr_str;
            for (std::size_t i = poly.cols; i > 0; --i) {
                std::size_t deg = i - 1;
                Scalar c = poly.at(0, deg);
                if (c == Scalar(0)) continue;
                std::string c_str = c.to_string();
                bool neg = (!c_str.empty() && c_str[0] == '-');
                if (neg) {
                    c_str = c_str.substr(1);
                    expr_str += expr_str.empty() ? "-" : " - ";
                } else if (!expr_str.empty()) {
                    expr_str += " + ";
                }
                if (deg == 0) {
                    expr_str += c_str;
                } else if (deg == 1) {
                    if (c_str == "1") expr_str += var;
                    else expr_str += c_str + "*" + var;
                } else {
                    if (c_str == "1") expr_str += var + "^" + std::to_string(deg);
                    else expr_str += c_str + "*" + var + "^" + std::to_string(deg);
                }
            }
            StoredValue res;
            res.is_string = true;
            res.string_value = expr_str.empty() ? "0" : expr_str;
            return res;
        }
        return make_matrix_result(std::move(poly));
    };

    // --- Matrix-only functions (1 matrix arg → scalar result) ---
    auto mat_scalar_1 = [](auto f, const std::string& name) {
        return [f, name](const std::vector<StoredValue>& args) -> StoredValue {
            if (args.size() != 1) throw std::runtime_error(name + " expects 1 argument");
            return make_scalar_result(f(require_matrix(args[0], name)));
        };
    };

    funcs["norm"] = mat_scalar_1(norm<Scalar>, "norm");
    funcs["trace"] = mat_scalar_1(trace<Scalar>, "trace");
    funcs["det"] = mat_scalar_1(determinant<Scalar>, "det");
    funcs["rank"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("rank expects 1 argument");
        return make_scalar_result(static_cast<long double>(rank<Scalar>(require_matrix(args[0], "rank"))));
    };
    funcs["cond"] = mat_scalar_1(condition_number<Scalar>, "cond");

    // --- Matrix-only functions (2 matrix args) ---
    funcs["outer"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("outer expects 2 arguments");
        return make_matrix_result(outer<Scalar>(require_matrix(args[0], "outer"), require_matrix(args[1], "outer")));
    };
    funcs["kron"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("kron expects 2 arguments");
        return make_matrix_result(kronecker<Scalar>(require_matrix(args[0], "kron"), require_matrix(args[1], "kron")));
    };
    funcs["hadamard"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("hadamard expects 2 arguments");
        return make_matrix_result(hadamard<Scalar>(require_matrix(args[0], "hadamard"), require_matrix(args[1], "hadamard")));
    };
    funcs["least_squares"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("least_squares expects 2 arguments");
        return make_matrix_result(least_squares<Scalar>(require_matrix(args[0], "least_squares"), require_matrix(args[1], "least_squares")));
    };
    funcs["solve"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("solve expects 2 arguments");
        return make_matrix_result(solve<Scalar>(require_matrix(args[0], "solve"), require_matrix(args[1], "solve")));
    };
    funcs["dot"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("dot expects 2 arguments");
        return make_scalar_result(dot<Scalar>(require_matrix(args[0], "dot"), require_matrix(args[1], "dot")));
    };

    // --- reshape (1 matrix + 2 scalar args) ---
    funcs["reshape"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 3) throw std::runtime_error("reshape expects 3 arguments");
        auto rows = static_cast<std::size_t>(static_cast<long double>(args[1].decimal) + 0.5);
        auto cols = static_cast<std::size_t>(static_cast<long double>(args[2].decimal) + 0.5);
        return make_matrix_result(reshape<Scalar>(require_matrix(args[0], "reshape"), rows, cols));
    };

    // --- get (1 matrix + 1 or 2 index args → scalar) ---
    funcs["get"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2 && args.size() != 3)
            throw std::runtime_error("get expects 2 or 3 arguments");
        Matrix m = require_matrix(args[0], "get");
        if (args.size() == 2) {
            return make_scalar_result(get<Scalar>(m, parse_index_argument(args[1], "get")));
        }
        return make_scalar_result(get<Scalar>(m,
            parse_index_argument(args[1], "get"),
            parse_index_argument(args[2], "get")));
    };

    // --- set (1 matrix + 2 or 3 args → matrix) ---
    funcs["set"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 3 && args.size() != 4)
            throw std::runtime_error("set expects 3 or 4 arguments");
        Matrix m = require_matrix(args[0], "set");
        if (args.size() == 3) {
            return make_matrix_result(set<Scalar>(m,
                parse_index_argument(args[1], "set"),
                args[2].decimal));
        }
        return make_matrix_result(set<Scalar>(m,
            parse_index_argument(args[1], "set"),
            parse_index_argument(args[2], "set"),
            args[3].decimal));
    };

    funcs["poly_eval"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("poly_eval expects coefficient vector and x");
        return make_scalar_result(polynomial_evaluate(require_vector(args[0], "poly_eval"), args[1].get_decimal()));
    };
    funcs["poly_deriv"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("poly_deriv expects one coefficient vector");
        return make_matrix_result(Matrix::vector(polynomial_derivative(require_vector(args[0], "poly_deriv"))));
    };
    funcs["poly_integ"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("poly_integ expects one coefficient vector");
        return make_matrix_result(Matrix::vector(polynomial_integral(require_vector(args[0], "poly_integ"))));
    };
    funcs["poly_compose"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("poly_compose expects two coefficient vectors");
        return make_matrix_result(Matrix::vector(polynomial_compose(require_vector(args[0], "poly_compose"),
                                                                     require_vector(args[1], "poly_compose"))));
    };
    funcs["poly_gcd"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("poly_gcd expects two coefficient vectors");
        return make_matrix_result(Matrix::vector(polynomial_gcd(require_vector(args[0], "poly_gcd"),
                                                                require_vector(args[1], "poly_gcd"))));
    };
    funcs["poly_fit"] = funcs["polynomial_fit"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 3) throw std::runtime_error("poly_fit expects x, y, degree");
        const Scalar degree = args[2].get_decimal();
        if (!mymath::is_integer(degree) || degree < 0) throw std::runtime_error("poly_fit degree must be a non-negative integer");
        return make_matrix_result(Matrix::vector(polynomial_fit(require_vector(args[0], "poly_fit"),
                                                                require_vector(args[1], "poly_fit"),
                                                                static_cast<int>(degree))));
    };
    funcs["lagrange"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 3) throw std::runtime_error("lagrange expects x, y, xi");
        return make_scalar_result(internal::lagrange_interpolate(require_vector(args[0], "lagrange"),
                                                                require_vector(args[1], "lagrange"), args[2].get_decimal()));
    };
    funcs["linear_regression"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("linear_regression expects two vectors");
        const auto fit = internal::linear_regression_fit(require_vector(args[0], "linear_regression"),
                                                        require_vector(args[1], "linear_regression"));
        return make_matrix_result(Matrix::vector(std::vector<Scalar>{fit.first, fit.second}));
    };

    // --- Complex construction ---
    funcs["complex"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("complex expects 2 arguments");
        return make_complex_result(ComplexNumber(args[0].decimal, args[1].decimal));
    };

    funcs["polar"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 2) throw std::runtime_error("polar expects 2 arguments");
        Scalar r = args[0].decimal, theta = args[1].decimal;
        return make_complex_result(ComplexNumber(r * mymath::cos(theta), r * mymath::sin(theta)));
    };

    // --- Complex/scalar extraction ---
    funcs["real"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("real expects 1 argument");
        return make_scalar_result(require_complex_argument(args[0], "real").real());
    };

    funcs["imag"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("imag expects 1 argument");
        return make_scalar_result(require_complex_argument(args[0], "imag").imag());
    };

    funcs["arg"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("arg expects 1 argument");
        const ComplexNumber z = require_complex_argument(args[0], "arg");
        const Scalar real = z.real(), imag = z.imag();
        const Scalar eps = matrix::internal::matrix_epsilon<Scalar>();
        if (mymath::is_near_zero(real, eps)) {
            if (mymath::is_near_zero(imag, eps)) return make_scalar_result(Scalar(0.0L));
            return make_scalar_result(imag > Scalar(0.0L) ? Scalar(mymath::kPi / 2.0) : Scalar(-mymath::kPi / 2.0));
        }
        Scalar angle = mymath::atan(imag / real);
        if (real < Scalar(0.0L)) {
            angle += imag >= Scalar(0.0L) ? Scalar(mymath::kPi) : Scalar(-mymath::kPi);
        }
        return make_scalar_result(angle);
    };

    funcs["conj"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("conj expects 1 argument");
        const ComplexNumber z = require_complex_argument(args[0], "conj");
        return make_complex_result(ComplexNumber(z.real(), -z.imag()));
    };

    // --- Polymorphic functions (scalar/complex/matrix) ---
    funcs["abs"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("abs expects 1 argument");
        if (args[0].is_matrix) {
            return make_scalar_result(norm<Scalar>(require_matrix(args[0], "abs")));
        }
        ComplexNumber z;
        if (try_complex_from_stored(args[0], &z) && args[0].is_complex) {
            return make_scalar_result(mymath::sqrt(z.real() * z.real() + z.imag() * z.imag()));
        }
        return make_scalar_result(mymath::abs(args[0].decimal));
    };

    funcs["exp"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("exp expects 1 argument");
        ComplexNumber z;
        if (try_complex_from_stored(args[0], &z) && args[0].is_complex) {
            Scalar m = mymath::exp(z.real());
            return make_complex_result(ComplexNumber(m * mymath::cos(z.imag()), m * mymath::sin(z.imag())));
        }
        return make_scalar_result(mymath::exp(args[0].decimal));
    };

    funcs["ln"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("ln expects 1 argument");
        ComplexNumber z;
        if (try_complex_from_stored(args[0], &z) && args[0].is_complex) {
            Scalar r = z.real(), i = z.imag();
            return make_complex_result(ComplexNumber(Scalar(0.5L) * mymath::ln(r * r + i * i), mymath::atan2(i, r)));
        }
        return make_scalar_result(mymath::ln(args[0].decimal));
    };

    funcs["sin"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("sin expects 1 argument");
        ComplexNumber z;
        if (try_complex_from_stored(args[0], &z) && args[0].is_complex) {
            Scalar r = z.real(), i = z.imag();
            return make_complex_result(ComplexNumber(mymath::sin(r) * mymath::cosh(i), mymath::cos(r) * mymath::sinh(i)));
        }
        return make_scalar_result(mymath::sin(args[0].decimal));
    };

    funcs["cos"] = [](const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() != 1) throw std::runtime_error("cos expects 1 argument");
        ComplexNumber z;
        if (try_complex_from_stored(args[0], &z) && args[0].is_complex) {
            Scalar r = z.real(), i = z.imag();
            return make_complex_result(ComplexNumber(mymath::cos(r) * mymath::cosh(i), -mymath::sin(r) * mymath::sinh(i)));
        }
        return make_scalar_result(mymath::cos(args[0].decimal));
    };

    return funcs;
}

std::vector<std::string> MatrixModule::get_function_names() const {
    std::vector<std::string> names;
    auto funcs = get_functions_map();
    for (const auto& [name, _] : funcs) names.push_back(name);

    // Matrix creation commands (handled separately by the expression parser)
    std::vector<std::string> creation = {"vec", "mat", "zeros", "eye", "identity", "randmat"};
    names.insert(names.end(), creation.begin(), creation.end());

    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
}

std::string MatrixModule::get_help_snippet(const std::string& topic) const {
    if (topic == "matrix") {
        return "Matrix guide:\n"
               "  Create:  [a,b;c,d] vec mat zeros eye identity randmat\n"
               "  Shape:   resize append_row append_col transpose\n"
               "  Anal.:   norm trace det rank rref eigvals solve cond diag\n"
               "  Complex: real imag arg conj polar complex abs";
    }
    return "";
}

REGISTER_CALCULATOR_MODULE(MatrixModule)

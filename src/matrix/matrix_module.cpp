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
#include "matrix_dsp.h"
#include "mymath.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "core/common/calculator_exceptions.h"
#include "math/helpers/integer_helpers.h"
#include <stdexcept>
#include <algorithm>
#include <sstream>

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
    if (!mymath::is_integer(val.decimal) || val.decimal < Scalar(0.0L)) {
        throw std::runtime_error(func_name + " requires non-negative integer index");
    }
    return static_cast<std::size_t>(static_cast<long double>(val.decimal) + 0.5);
}

StoredValue make_matrix_result(Matrix m) {
    StoredValue res;
    res.is_matrix = true;
    res.matrix_ptr = std::make_shared<Matrix>(std::move(m));
    return res;
}

StoredValue make_scalar_result(Scalar s) {
    StoredValue res;
    res.decimal = s;
    return res;
}

StoredValue make_complex_result(ComplexNumber z) {
    StoredValue res;
    res.is_complex = true;
    res.complex = z;
    return res;
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

    // --- Matrix-only functions (1 matrix arg → matrix result) ---
    auto mat_func_1 = [](auto f, const std::string& name) {
        return [f, name](const std::vector<StoredValue>& args) -> StoredValue {
            if (args.size() != 1) throw std::runtime_error(name + " expects 1 argument");
            return make_matrix_result(f(require_matrix(args[0], name)));
        };
    };

    funcs["transpose"] = mat_func_1(transpose<Scalar>, "transpose");
    funcs["inverse"] = mat_func_1(inverse<Scalar>, "inverse");
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
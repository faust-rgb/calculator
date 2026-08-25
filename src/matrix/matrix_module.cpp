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
#include "mymath.h"
#include "core/services/string_utils.h"
#include "core/services/format_utils.h"
#include "core/common/calculator_exceptions.h"
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <random>
#include <cmath>

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

    return funcs;
}

std::vector<std::string> MatrixModule::get_function_names() const {
    static const std::vector<std::string> names = {
        "append_col", "append_row", "charpoly", "cholesky", "cond", "det",
        "diag", "dot", "eigvals", "eigvecs", "expm", "eye", "get", "hadamard",
        "hessenberg", "identity", "inv", "inverse", "kron", "least_squares",
        "lu_l", "lu_p", "lu_u", "mat", "norm", "null", "outer", "pinv",
        "qr_q", "qr_r", "randmat", "random_matrix", "rank", "reshape",
        "resize", "rref", "schur", "set", "solve", "svd_s", "svd_u",
        "svd_vt", "trace", "transpose", "vec", "zeros"
    };
    return names;
}

std::string MatrixModule::get_help_snippet(const std::string& topic) const {
    if (topic == "matrix") {
        return "Matrix guide:\n"
               "  Create:  [a,b;c,d] vec mat zeros eye identity randmat\n"
               "  Shape:   resize append_row append_col transpose reshape\n"
               "  Anal.:   norm trace det rank rref eigvals eigvecs solve cond diag\n"
               "  Decomp:  lu_l lu_u lu_p qr_q qr_r svd_u svd_s svd_vt cholesky schur hessenberg\n"
               "  Access:  get set";
    }
    return "";
}

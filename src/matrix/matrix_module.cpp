#include "matrix_module.h"
#include "matrix.h"
#include "matrix_internal.h"
#include "mymath.h"
#include "../core/string_utils.h"
#include "../core/calculator_exceptions.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace {

using namespace matrix;

// 辅助函数：解析索引参数
std::size_t parse_index_argument(const std::string& expression,
                                 const ScalarEvaluator& scalar_evaluator,
                                 const std::string& func_name) {
    const double value = scalar_evaluator(expression);
    if (!mymath::is_integer(value) || value < 0.0) {
        throw std::runtime_error(func_name + " requires non-negative integer index");
    }
    return static_cast<std::size_t>(value >= 0.0 ? value + 0.5 : value - 0.5);
}

// 辅助函数：要求参数为矩阵
Matrix require_matrix(const std::string& expression,
                      const std::string& func_name,
                      const ScalarEvaluator& scalar_evaluator,
                      const MatrixLookup& matrix_lookup,
                      const ComplexLookup& complex_lookup,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* matrix_functions) {
    Value value;
    if (!try_evaluate_expression(expression, scalar_evaluator, matrix_lookup, complex_lookup, matrix_functions, nullptr, &value) ||
        !value.is_matrix) {
        throw std::runtime_error(func_name + " expects a matrix argument");
    }
    return value.matrix;
}

// 辅助函数：尝试从 Value 获取复数
bool try_complex_from_value(const Value& v, ComplexNumber* z) {
    if (v.is_complex) {
        *z = v.complex;
        return true;
    }
    if (!v.is_matrix) {
        z->real = v.scalar;
        z->imag = 0.0;
        return true;
    }
    return false;
}

// 辅助函数：要求参数为复数
ComplexNumber require_complex_argument(const std::string& expression,
                                       const std::string& func_name,
                                       const ScalarEvaluator& scalar_evaluator,
                                       const MatrixLookup& matrix_lookup,
                                       const ComplexLookup& complex_lookup,
                                       const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* matrix_functions) {
    Value value;
    if (try_evaluate_expression(expression, scalar_evaluator, matrix_lookup, complex_lookup, matrix_functions, nullptr, &value)) {
        ComplexNumber z;
        if (try_complex_from_value(value, &z)) {
            return z;
        }
    }
    throw std::runtime_error(func_name + " expects a scalar or complex argument");
}

// 格式化特征值矩阵（复数形式）
std::string format_eigenvalue_matrix(const Matrix& values) {
    if (values.rows <= 1 || values.cols != 2) {
        return values.to_string();
    }
    std::ostringstream out;
    out << "[";
    for (std::size_t row = 0; row < values.rows; ++row) {
        if (row != 0) {
            out << ", ";
        }
        out << internal::format_complex({values.at(row, 0), values.at(row, 1)});
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


std::string MatrixModule::execute_args(const std::string& command,
                                       const std::vector<std::string>& args,
                                       const CoreServices& services) {
    // 命令已由路由层验证，无需再检查

    if (args.size() != 1) {
        throw std::runtime_error(command + " expects exactly one matrix argument");
    }

    const std::string& arg = args[0];
    if (!services.is_matrix_argument(arg)) {
        throw std::runtime_error(command + " expects a matrix argument");
    }

    const Matrix matrix_value = services.parse_matrix_argument(arg, command);

    if (command == "svd") {
        return "U: " + svd_u(matrix_value).to_string() +
               "\nS: " + svd_s(matrix_value).to_string() +
               "\nVt: " + svd_vt(matrix_value).to_string();
    }

    if (command == "lu_p") {
        return lu_p(matrix_value).to_string();
    }

    // command == "eig"
    try {
        const Matrix values = eigenvalues(matrix_value);
        if (values.rows > 1 && values.cols == 2) {
            return "values: " + format_eigenvalue_matrix(values) +
                   "\nvectors: unavailable for complex eigenvalues";
        }
        return "values: " + values.to_string() +
               "\nvectors: " + eigenvectors(matrix_value).to_string();
    } catch (const std::exception&) {
        if (matrix_value.rows == 2 && matrix_value.cols == 2) {
            const double trace = matrix_value.at(0, 0) + matrix_value.at(1, 1);
            const double det = determinant(matrix_value);
            const double discriminant = trace * trace - 4.0 * det;
            if (discriminant < 0.0) {
                const double real = trace * 0.5;
                const double imag = mymath::sqrt(-discriminant) * 0.5;
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

std::map<std::string, std::function<matrix::Matrix(const std::vector<matrix::Matrix>&)>>
MatrixModule::get_matrix_functions() const {
    std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>> funcs;

    // 矩阵转置
    funcs["transpose"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("transpose expects 1 argument");
        return transpose(args[0]);
    };

    // 矩阵求逆
    funcs["inverse"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("inverse expects 1 argument");
        return inverse(args[0]);
    };

    // 伪逆
    funcs["pinv"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("pinv expects 1 argument");
        return pseudo_inverse(args[0]);
    };

    // 外积
    funcs["outer"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 2) throw std::runtime_error("outer expects 2 arguments");
        return outer(args[0], args[1]);
    };

    // Kronecker 积
    funcs["kron"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 2) throw std::runtime_error("kron expects 2 arguments");
        return kronecker(args[0], args[1]);
    };

    // Hadamard 积
    funcs["hadamard"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 2) throw std::runtime_error("hadamard expects 2 arguments");
        return hadamard(args[0], args[1]);
    };

    // 零空间
    funcs["null"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("null expects 1 argument");
        return nullspace(args[0]);
    };

    // 最小二乘
    funcs["least_squares"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 2) throw std::runtime_error("least_squares expects 2 arguments");
        return least_squares(args[0], args[1]);
    };

    // QR 分解
    funcs["qr_q"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("qr_q expects 1 argument");
        return qr_q(args[0]);
    };
    funcs["qr_r"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("qr_r expects 1 argument");
        return qr_r(args[0]);
    };

    // LU 分解
    funcs["lu_l"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("lu_l expects 1 argument");
        return lu_l(args[0]);
    };
    funcs["lu_u"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("lu_u expects 1 argument");
        return lu_u(args[0]);
    };
    funcs["lu_p"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("lu_p expects 1 argument");
        return lu_p(args[0]);
    };

    // SVD 分解
    funcs["svd_u"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("svd_u expects 1 argument");
        return svd_u(args[0]);
    };
    funcs["svd_s"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("svd_s expects 1 argument");
        return svd_s(args[0]);
    };
    funcs["svd_vt"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("svd_vt expects 1 argument");
        return svd_vt(args[0]);
    };

    // Cholesky 分解
    funcs["cholesky"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("cholesky expects 1 argument");
        return cholesky(args[0]);
    };

    // Hessenberg 形式
    funcs["hessenberg"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("hessenberg expects 1 argument");
        return hessenberg(args[0]);
    };

    // Schur 分解
    funcs["schur"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("schur expects 1 argument");
        return schur(args[0]);
    };

    // 滤波器
    funcs["filter"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 3) throw std::runtime_error("filter expects 3 arguments");
        return filter(args[0], args[1], args[2]);
    };

    // 频率响应
    funcs["freqz"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() < 2 || args.size() > 3) throw std::runtime_error("freqz expects 2 or 3 arguments");
        std::size_t n = 512;
        if (args.size() == 3 && args[2].rows > 0 && args[2].cols > 0) {
            n = static_cast<std::size_t>(args[2].at(0, 0) + 0.5);
        }
        return freqz(args[0], args[1], n);
    };

    // 部分分式展开
    funcs["residue"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 2) throw std::runtime_error("residue expects 2 arguments");
        return residue(args[0], args[1]);
    };

    // 对角矩阵
    funcs["diag"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 1) throw std::runtime_error("diag expects 1 argument");
        return diag(args[0]);
    };

    // 矩阵重塑
    funcs["reshape"] = [](const std::vector<Matrix>& args) -> Matrix {
        if (args.size() != 3) throw std::runtime_error("reshape expects 3 arguments");
        const std::size_t rows = static_cast<std::size_t>(args[1].at(0, 0) + 0.5);
        const std::size_t cols = static_cast<std::size_t>(args[2].at(0, 0) + 0.5);
        return reshape(args[0], rows, cols);
    };

    return funcs;
}

std::map<std::string, CalculatorModule::ValueFunction>
MatrixModule::get_value_functions() const {
    std::map<std::string, ValueFunction> funcs;

    // complex - 构造复数
    funcs["complex"] = [](const std::vector<std::string>& args,
                          const ScalarEvaluator& se,
                          const MatrixLookup&,
                          const ComplexLookup&,
                          const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>*) -> Value {
        if (args.size() != 2) throw std::runtime_error("complex expects 2 arguments");
        return Value::from_complex(se(args[0]), se(args[1]));
    };

    // polar - 极坐标构造复数
    funcs["polar"] = [](const std::vector<std::string>& args,
                        const ScalarEvaluator& se,
                        const MatrixLookup&,
                        const ComplexLookup&,
                        const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>*) -> Value {
        if (args.size() != 2) throw std::runtime_error("polar expects 2 arguments");
        const double r = se(args[0]);
        const double theta = se(args[1]);
        return Value::from_complex(r * mymath::cos(theta), r * mymath::sin(theta));
    };

    // real - 取实部
    funcs["real"] = [](const std::vector<std::string>& args,
                       const ScalarEvaluator& se,
                       const MatrixLookup& ml,
                       const ComplexLookup& cl,
                       const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("real expects 1 argument");
        const ComplexNumber z = require_complex_argument(args[0], "real", se, ml, cl, mf);
        return Value::from_scalar(z.real);
    };

    // imag - 取虚部
    funcs["imag"] = [](const std::vector<std::string>& args,
                       const ScalarEvaluator& se,
                       const MatrixLookup& ml,
                       const ComplexLookup& cl,
                       const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("imag expects 1 argument");
        const ComplexNumber z = require_complex_argument(args[0], "imag", se, ml, cl, mf);
        return Value::from_scalar(z.imag);
    };

    // arg - 复数辐角
    funcs["arg"] = [](const std::vector<std::string>& args,
                      const ScalarEvaluator& se,
                      const MatrixLookup& ml,
                      const ComplexLookup& cl,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("arg expects 1 argument");
        const ComplexNumber z = require_complex_argument(args[0], "arg", se, ml, cl, mf);
        const double real = z.real;
        const double imag = z.imag;
        if (mymath::is_near_zero(real, 1e-12)) {
            if (mymath::is_near_zero(imag, 1e-12)) {
                return Value::from_scalar(0.0);
            }
            return Value::from_scalar(imag > 0.0 ? mymath::kPi / 2.0 : -mymath::kPi / 2.0);
        }
        double angle = mymath::atan(imag / real);
        if (real < 0.0) {
            angle += imag >= 0.0 ? mymath::kPi : -mymath::kPi;
        }
        return Value::from_scalar(angle);
    };

    // conj - 共轭复数
    funcs["conj"] = [](const std::vector<std::string>& args,
                       const ScalarEvaluator& se,
                       const MatrixLookup& ml,
                       const ComplexLookup& cl,
                       const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("conj expects 1 argument");
        const ComplexNumber z = require_complex_argument(args[0], "conj", se, ml, cl, mf);
        return Value::from_complex(z.real, -z.imag);
    };

    // abs - 多态绝对值
    funcs["abs"] = [](const std::vector<std::string>& args,
                      const ScalarEvaluator& se,
                      const MatrixLookup& ml,
                      const ComplexLookup& cl,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("abs expects 1 argument");
        Value v;
        if (try_evaluate_expression(args[0], se, ml, cl, mf, nullptr, &v)) {
            ComplexNumber z;
            if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                const double r = z.real, i = z.imag;
                return Value::from_scalar(mymath::sqrt(r * r + i * i));
            } else if (v.is_matrix) {
                return Value::from_scalar(norm(v.matrix));
            } else {
                return Value::from_scalar(mymath::abs(v.scalar));
            }
        }
        return Value::from_scalar(se("abs(" + args[0] + ")"));
    };

    // exp - 多态指数
    funcs["exp"] = [](const std::vector<std::string>& args,
                      const ScalarEvaluator& se,
                      const MatrixLookup& ml,
                      const ComplexLookup& cl,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("exp expects 1 argument");
        Value v;
        if (try_evaluate_expression(args[0], se, ml, cl, mf, nullptr, &v)) {
            ComplexNumber z;
            if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                const double r = z.real, i = z.imag, m = mymath::exp(r);
                return Value::from_complex(m * mymath::cos(i), m * mymath::sin(i));
            } else if (!v.is_matrix) {
                return Value::from_scalar(mymath::exp(v.scalar));
            }
        }
        return Value::from_scalar(se("exp(" + args[0] + ")"));
    };

    // ln - 多态对数
    funcs["ln"] = [](const std::vector<std::string>& args,
                     const ScalarEvaluator& se,
                     const MatrixLookup& ml,
                     const ComplexLookup& cl,
                     const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("ln expects 1 argument");
        Value v;
        if (try_evaluate_expression(args[0], se, ml, cl, mf, nullptr, &v)) {
            ComplexNumber z;
            if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                const double r = z.real, i = z.imag;
                return Value::from_complex(0.5 * mymath::ln(r * r + i * i), mymath::atan2(i, r));
            } else if (!v.is_matrix) {
                return Value::from_scalar(mymath::ln(v.scalar));
            }
        }
        return Value::from_scalar(se("ln(" + args[0] + ")"));
    };

    // sin - 多态正弦
    funcs["sin"] = [](const std::vector<std::string>& args,
                      const ScalarEvaluator& se,
                      const MatrixLookup& ml,
                      const ComplexLookup& cl,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("sin expects 1 argument");
        Value v;
        if (try_evaluate_expression(args[0], se, ml, cl, mf, nullptr, &v)) {
            ComplexNumber z;
            if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                const double r = z.real, i = z.imag;
                return Value::from_complex(mymath::sin(r) * mymath::cosh(i), mymath::cos(r) * mymath::sinh(i));
            } else if (!v.is_matrix) {
                return Value::from_scalar(mymath::sin(v.scalar));
            }
        }
        return Value::from_scalar(se("sin(" + args[0] + ")"));
    };

    // cos - 多态余弦
    funcs["cos"] = [](const std::vector<std::string>& args,
                      const ScalarEvaluator& se,
                      const MatrixLookup& ml,
                      const ComplexLookup& cl,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("cos expects 1 argument");
        Value v;
        if (try_evaluate_expression(args[0], se, ml, cl, mf, nullptr, &v)) {
            ComplexNumber z;
            if (try_complex_from_value(v, &z) && (v.is_complex || v.is_matrix)) {
                const double r = z.real, i = z.imag;
                return Value::from_complex(mymath::cos(r) * mymath::cosh(i), -mymath::sin(r) * mymath::sinh(i));
            } else if (!v.is_matrix) {
                return Value::from_scalar(mymath::cos(v.scalar));
            }
        }
        return Value::from_scalar(se("cos(" + args[0] + ")"));
    };

    // norm - 矩阵范数
    funcs["norm"] = [](const std::vector<std::string>& args,
                       const ScalarEvaluator& se,
                       const MatrixLookup& ml,
                       const ComplexLookup& cl,
                       const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("norm expects 1 argument");
        const Matrix m = require_matrix(args[0], "norm", se, ml, cl, mf);
        return Value::from_scalar(norm(m));
    };

    // dot - 点积
    funcs["dot"] = [](const std::vector<std::string>& args,
                      const ScalarEvaluator& se,
                      const MatrixLookup& ml,
                      const ComplexLookup& cl,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 2) throw std::runtime_error("dot expects 2 arguments");
        const Matrix a = require_matrix(args[0], "dot", se, ml, cl, mf);
        const Matrix b = require_matrix(args[1], "dot", se, ml, cl, mf);
        return Value::from_scalar(dot(a, b));
    };

    // cond - 条件数
    funcs["cond"] = [](const std::vector<std::string>& args,
                       const ScalarEvaluator& se,
                       const MatrixLookup& ml,
                       const ComplexLookup& cl,
                       const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("cond expects 1 argument");
        const Matrix m = require_matrix(args[0], "cond", se, ml, cl, mf);
        return Value::from_scalar(condition_number(m));
    };

    // trace - 迹
    funcs["trace"] = [](const std::vector<std::string>& args,
                        const ScalarEvaluator& se,
                        const MatrixLookup& ml,
                        const ComplexLookup& cl,
                        const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("trace expects 1 argument");
        const Matrix m = require_matrix(args[0], "trace", se, ml, cl, mf);
        return Value::from_scalar(trace(m));
    };

    // det - 行列式
    funcs["det"] = [](const std::vector<std::string>& args,
                      const ScalarEvaluator& se,
                      const MatrixLookup& ml,
                      const ComplexLookup& cl,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("det expects 1 argument");
        const Matrix m = require_matrix(args[0], "det", se, ml, cl, mf);
        return Value::from_scalar(determinant(m));
    };

    // rank - 秩
    funcs["rank"] = [](const std::vector<std::string>& args,
                       const ScalarEvaluator& se,
                       const MatrixLookup& ml,
                       const ComplexLookup& cl,
                       const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("rank expects 1 argument");
        const Matrix m = require_matrix(args[0], "rank", se, ml, cl, mf);
        return Value::from_scalar(static_cast<double>(rank(m)));
    };

    // rref - 行最简形
    funcs["rref"] = [](const std::vector<std::string>& args,
                       const ScalarEvaluator& se,
                       const MatrixLookup& ml,
                       const ComplexLookup& cl,
                       const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("rref expects 1 argument");
        const Matrix m = require_matrix(args[0], "rref", se, ml, cl, mf);
        return Value::from_matrix(rref(m));
    };

    // eigvals - 特征值
    funcs["eigvals"] = [](const std::vector<std::string>& args,
                          const ScalarEvaluator& se,
                          const MatrixLookup& ml,
                          const ComplexLookup& cl,
                          const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 1) throw std::runtime_error("eigvals expects 1 argument");
        const Matrix m = require_matrix(args[0], "eigvals", se, ml, cl, mf);
        return Value::from_matrix(eigenvalues(m));
    };

    // solve - 线性方程组求解
    funcs["solve"] = [](const std::vector<std::string>& args,
                        const ScalarEvaluator& se,
                        const MatrixLookup& ml,
                        const ComplexLookup& cl,
                        const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 2) throw std::runtime_error("solve expects 2 arguments");
        const Matrix A = require_matrix(args[0], "solve", se, ml, cl, mf);
        const Matrix b = require_matrix(args[1], "solve", se, ml, cl, mf);
        return Value::from_matrix(solve(A, b));
    };

    // get - 获取矩阵元素
    funcs["get"] = [](const std::vector<std::string>& args,
                      const ScalarEvaluator& se,
                      const MatrixLookup& ml,
                      const ComplexLookup& cl,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 2 && args.size() != 3) {
            throw std::runtime_error("get expects 2 or 3 arguments");
        }
        Matrix m = require_matrix(args[0], "get", se, ml, cl, mf);
        if (args.size() == 2) {
            return Value::from_scalar(get(m, parse_index_argument(args[1], se, "get")));
        }
        return Value::from_scalar(get(m,
            parse_index_argument(args[1], se, "get"),
            parse_index_argument(args[2], se, "get")));
    };

    // set - 设置矩阵元素
    funcs["set"] = [](const std::vector<std::string>& args,
                      const ScalarEvaluator& se,
                      const MatrixLookup& ml,
                      const ComplexLookup& cl,
                      const std::map<std::string, std::function<Matrix(const std::vector<Matrix>&)>>* mf) -> Value {
        if (args.size() != 3 && args.size() != 4) {
            throw std::runtime_error("set expects 3 or 4 arguments");
        }
        Matrix m = require_matrix(args[0], "set", se, ml, cl, mf);
        if (args.size() == 3) {
            return Value::from_matrix(set(m, parse_index_argument(args[1], se, "set"), se(args[2])));
        }
        return Value::from_matrix(set(m,
            parse_index_argument(args[1], se, "set"),
            parse_index_argument(args[2], se, "set"),
            se(args[3])));
    };

    return funcs;
}

std::vector<std::string> MatrixModule::get_functions() const {
    std::vector<std::string> names;
    auto mat_funcs = get_matrix_functions();
    for (const auto& [name, _] : mat_funcs) names.push_back(name);
    auto val_funcs = get_value_functions();
    for (const auto& [name, _] : val_funcs) names.push_back(name);
    
    // Add matrix creation functions
    std::vector<std::string> creation = { "vec", "mat", "zeros", "eye", "identity", "randmat" };
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

// ============================================================================
// critical_point_solver.cpp - 临界点求解算法
// ============================================================================
//
// 本文件实现了 critical 命令的核心算法逻辑，包括：
// - 单变量函数的临界点求解（二分法 + Newton-Raphson）
// - 多变量函数的临界点求解（梯度下降 + Hessian矩阵分析）
// - Groebner 基分析用于多项式系统
//
// 这些算法从 analysis_module.cpp 中拆分出来，使模块类更专注于路由。
// ============================================================================

#include "critical_point_solver.h"
#include "math/base/precision_constants.h"
#include "symbolic/algebra/groebner/groebner_basis.h"
#include "symbolic/solver/symbolic_solver.h"
#include "math/mymath.h"
#include <algorithm>
#include <cmath>

namespace analysis {

using Scalar = mymath::Scalar;

// ============================================================================
// AST 数值求值
// ============================================================================

Scalar evaluate_ast_node(const std::shared_ptr<SymbolicExpression::Node>& node,
                         const std::map<std::string, Scalar>& point) {
    if (!node) return Scalar(0);

    switch (node->type) {
        case NodeType::kNumber:
            return Scalar(node->number_value);
        case NodeType::kVariable: {
            auto it = point.find(node->text);
            if (it != point.end()) return it->second;
            throw std::runtime_error("unknown variable: " + node->text);
        }
        case NodeType::kPi:
            return Scalar(mymath::pi());
        case NodeType::kE:
            return Scalar(mymath::e());
        case NodeType::kAdd:
            return evaluate_ast_node(node->left, point) + evaluate_ast_node(node->right, point);
        case NodeType::kSubtract:
            return evaluate_ast_node(node->left, point) - evaluate_ast_node(node->right, point);
        case NodeType::kMultiply:
            return evaluate_ast_node(node->left, point) * evaluate_ast_node(node->right, point);
        case NodeType::kDivide:
            return evaluate_ast_node(node->left, point) / evaluate_ast_node(node->right, point);
        case NodeType::kPower:
            return mymath::pow(evaluate_ast_node(node->left, point),
                              evaluate_ast_node(node->right, point));
        case NodeType::kNegate:
            return -evaluate_ast_node(node->left, point);
        case NodeType::kFunction: {
            Scalar arg = evaluate_ast_node(node->left, point);
            const std::string& fname = node->text;
            if (fname == "sin") return mymath::sin(arg);
            if (fname == "cos") return mymath::cos(arg);
            if (fname == "tan") return mymath::tan(arg);
            if (fname == "exp") return mymath::exp(arg);
            if (fname == "ln" || fname == "log") return mymath::log(arg);
            if (fname == "sqrt") return mymath::sqrt(arg);
            throw std::runtime_error("unsupported function in numeric eval: " + fname);
        }
        default:
            throw std::runtime_error("unsupported node type in numeric eval");
    }
}

// ============================================================================
// 单变量临界点求解
// ============================================================================

std::vector<Scalar> find_univariate_critical_points(
    const SymbolicExpression& derivative,
    const std::string& variable,
    const std::function<Scalar(Scalar)>& normalize_result) {

    std::vector<Scalar> critical_points;

    auto eval_derivative = [&](Scalar x) -> Scalar {
        try {
            return evaluate_ast_node(derivative.node_, {{variable, x}});
        } catch (...) {
            throw std::runtime_error("derivative is not numeric at this point");
        }
    };

    auto add_point = [&](Scalar x) {
        for (const auto& existing : critical_points) {
            if (mymath::abs(existing - x) < Scalar(1e-5L)) return;
        }
        critical_points.push_back(normalize_result(x));
    };

    const Scalar scan_min = Scalar(-100);
    const Scalar scan_max = Scalar(100);
    const int coarse_segments = 512;

    // 计算二阶导数
    SymbolicExpression second_deriv = derivative.derivative(variable).simplify();
    auto eval_second = [&](Scalar x) -> Scalar {
        try {
            return evaluate_ast_node(second_deriv.node_, {{variable, x}});
        } catch (...) {
            return Scalar(0);
        }
    };

    // 扫描区间寻找符号变化
    Scalar previous_x = scan_min;
    Scalar previous_value = eval_derivative(previous_x);

    for (int i = 1; i <= coarse_segments; ++i) {
        const Scalar current_x = scan_min + (scan_max - scan_min) * Scalar(i) / Scalar(coarse_segments);
        const Scalar current_value = eval_derivative(current_x);

        // 二分法寻找根
        if ((previous_value < Scalar(0) && current_value > Scalar(0)) ||
            (previous_value > Scalar(0) && current_value < Scalar(0))) {
            Scalar left = previous_x;
            Scalar right = current_x;
            Scalar left_value = previous_value;

            for (int iter = 0; iter < 80; ++iter) {
                const Scalar mid = (left + right) * Scalar(0.5L);
                const Scalar mid_value = eval_derivative(mid);

                if (mymath::isfinite(mid_value) &&
                    mymath::abs(mid_value) < precision::newton_tolerance<Scalar>()) {
                    left = right = mid;
                    break;
                }

                if ((left_value < Scalar(0) && mid_value > Scalar(0)) ||
                    (left_value > Scalar(0) && mid_value < Scalar(0))) {
                    right = mid;
                } else {
                    left = mid;
                    left_value = mid_value;
                }
            }
            add_point((left + right) * Scalar(0.5L));
        }
        previous_x = current_x;
        previous_value = current_value;
    }

    // Newton-Raphson 精化
    for (int i = 0; i <= coarse_segments; ++i) {
        const Scalar x = scan_min + (scan_max - scan_min) * Scalar(i) / Scalar(coarse_segments);
        const Scalar deriv_val = eval_derivative(x);
        const Scalar second_val = eval_second(x);

        if (mymath::abs(deriv_val) < precision::gradient_convergence_threshold<Scalar>() &&
            mymath::abs(second_val) > precision::sqrt_epsilon<Scalar>()) {
            Scalar refined_x = x;

            for (int iter = 0; iter < 20; ++iter) {
                const Scalar f_prime = eval_derivative(refined_x);
                const Scalar f_double_prime = eval_second(refined_x);

                if (mymath::isfinite(f_prime) &&
                    mymath::abs(f_prime) < precision::newton_tolerance<Scalar>()) break;
                if (mymath::abs(f_double_prime) < precision::gradient_convergence_threshold<Scalar>()) break;

                refined_x = refined_x - f_prime / f_double_prime;
            }

            if (mymath::isfinite(eval_derivative(refined_x)) &&
                mymath::abs(eval_derivative(refined_x)) < precision::newton_tolerance<Scalar>() &&
                mymath::abs(eval_second(refined_x)) > precision::sqrt_epsilon<Scalar>()) {
                add_point(refined_x);
            }
        }
    }

    return critical_points;
}

// ============================================================================
// 多变量临界点求解
// ============================================================================

std::vector<std::map<std::string, Scalar>> find_multivariate_critical_points(
    const std::vector<SymbolicExpression>& gradient,
    const std::vector<std::string>& variables,
    const std::function<Scalar(Scalar)>& normalize_result) {

    std::vector<std::map<std::string, Scalar>> critical_points;

    auto eval_gradient_at = [&](const SymbolicExpression& g,
                                const std::map<std::string, Scalar>& point) -> Scalar {
        return evaluate_ast_node(g.node_, point);
    };

    auto gradient_norm_at = [&](const std::map<std::string, Scalar>& point) -> Scalar {
        Scalar norm = Scalar(0);
        for (const auto& g : gradient) {
            const Scalar val = eval_gradient_at(g, point);
            norm = norm + val * val;
        }
        return mymath::sqrt(norm);
    };

    auto add_critical_point = [&](const std::map<std::string, Scalar>& point) {
        bool duplicate = false;
        for (const auto& existing : critical_points) {
            bool same = true;
            for (const auto& v : variables) {
                const auto it_existing = existing.find(v);
                const auto it_current = point.find(v);
                if (it_existing == existing.end() || it_current == point.end() ||
                    mymath::abs(it_existing->second - it_current->second) > Scalar(1e-2L)) {
                    same = false;
                    break;
                }
            }
            if (same) { duplicate = true; break; }
        }
        if (!duplicate) {
            std::map<std::string, Scalar> normalized;
            for (const auto& [k, v] : point) {
                normalized[k] = normalize_result(v);
            }
            critical_points.push_back(normalized);
        }
    };

    // 收集起始点
    std::vector<std::map<std::string, Scalar>> starting_points;

    // 检查是否为多项式系统，尝试 Groebner 基分析
    auto is_polynomial_node = [](const std::shared_ptr<SymbolicExpression::Node>& node) -> bool {
        std::function<bool(const std::shared_ptr<SymbolicExpression::Node>&)> check =
            [&](const std::shared_ptr<SymbolicExpression::Node>& n) -> bool {
            if (!n) return true;
            switch (n->type) {
                case NodeType::kNumber:
                case NodeType::kPi:
                case NodeType::kE:
                case NodeType::kVariable:
                    return true;
                case NodeType::kAdd:
                case NodeType::kSubtract:
                case NodeType::kMultiply:
                    return check(n->left) && check(n->right);
                case NodeType::kNegate:
                    return check(n->left);
                case NodeType::kPower: {
                    Scalar val;
                    if (SymbolicExpression(n->right).is_number(&val)) {
                        if (val >= Scalar(0) &&
                            mymath::abs(val - mymath::floor(val)) < precision::epsilon<Scalar>() * Scalar(100)) {
                            return check(n->left);
                        }
                    }
                    return false;
                }
                default:
                    return false;
            }
        };
        return check(node);
    };

    bool is_poly_system = true;
    for (const auto& g : gradient) {
        if (!is_polynomial_node(g.node_)) {
            is_poly_system = false;
            break;
        }
    }

    if (is_poly_system && variables.size() > 1) {
        try {
            auto basis = symbolic_groebner::compute_groebner_basis(gradient, variables);
            if (!basis.empty()) {
                const SymbolicExpression& last_poly = basis.back();
                symbolic_solver::SymbolicSolver solver;
                auto sol = solver.solve(last_poly, variables.back());
                for (const auto& val_expr : sol.values) {
                    Scalar val;
                    if (val_expr.is_number(&val)) {
                        std::map<std::string, Scalar> pt;
                        pt[variables.back()] = val;
                        for (size_t i = 0; i < variables.size() - 1; ++i) pt[variables[i]] = Scalar(0);
                        starting_points.push_back(pt);

                        std::map<std::string, Scalar> pt1;
                        pt1[variables.back()] = val;
                        for (size_t i = 0; i < variables.size() - 1; ++i) pt1[variables[i]] = Scalar(1);
                        starting_points.push_back(pt1);
                    }
                }
            }
        } catch (...) {
            // Groebner 基计算失败，继续使用其他方法
        }
    }

    // 添加原点和网格点作为起始点
    std::map<std::string, Scalar> origin;
    for (const auto& v : variables) origin[v] = Scalar(0);
    starting_points.push_back(origin);

    const std::vector<Scalar> grid_values = {Scalar(-10), Scalar(-1), Scalar(1), Scalar(10)};
    if (variables.size() <= 3) {
        for (Scalar v0 : grid_values) {
            for (Scalar v1 : grid_values) {
                std::map<std::string, Scalar> pt;
                pt[variables[0]] = v0;
                pt[variables[1]] = v1;
                if (variables.size() == 3) {
                    for (Scalar v2 : grid_values) {
                        pt[variables[2]] = v2;
                        starting_points.push_back(pt);
                    }
                } else {
                    starting_points.push_back(pt);
                }
            }
        }
    }

    // 添加随机起始点
    for (int i = 0; i < 5; ++i) {
        std::map<std::string, Scalar> pt;
        for (const auto& v : variables) {
            pt[v] = Scalar(static_cast<long double>(rand()) / RAND_MAX * 20.0L - 10.0L);
        }
        starting_points.push_back(pt);
    }

    // 计算 Hessian 矩阵
    std::vector<std::vector<SymbolicExpression>> symbolic_hessian(
        variables.size(), std::vector<SymbolicExpression>(variables.size()));
    for (size_t i = 0; i < variables.size(); ++i) {
        for (size_t j = 0; j < variables.size(); ++j) {
            symbolic_hessian[i][j] = gradient[i].derivative(variables[j]).simplify();
        }
    }

    // 从每个起始点进行优化
    for (const auto& start : starting_points) {
        try {
            std::map<std::string, Scalar> current = start;

            if (gradient_norm_at(current) < precision::gradient_convergence_threshold<Scalar>()) {
                add_critical_point(current);
                continue;
            }

            for (int iter = 0; iter < 100; ++iter) {
                std::vector<Scalar> rhs(variables.size(), Scalar(0));
                std::vector<std::vector<Scalar>> jac(
                    variables.size(), std::vector<Scalar>(variables.size(), Scalar(0)));

                Scalar current_norm = gradient_norm_at(current);
                if (current_norm < precision::newton_tolerance<Scalar>()) break;

                bool eval_ok = true;
                for (size_t row = 0; row < variables.size(); ++row) {
                    rhs[row] = Scalar(0) - eval_gradient_at(gradient[row], current);
                    for (size_t col = 0; col < variables.size(); ++col) {
                        try {
                            jac[row][col] = evaluate_ast_node(symbolic_hessian[row][col].node_, current);
                        } catch (...) {
                            eval_ok = false;
                            break;
                        }
                    }
                    if (!eval_ok) break;
                }

                if (!eval_ok) break;

                // 求解线性系统
                const std::vector<Scalar> delta = solve_linear_system(jac, rhs);

                // 线搜索
                Scalar alpha = Scalar(1.0L);
                bool line_search_success = false;

                for (int ls_iter = 0; ls_iter < 10; ++ls_iter) {
                    std::map<std::string, Scalar> next = current;
                    for (size_t i = 0; i < variables.size(); ++i) {
                        next[variables[i]] = next[variables[i]] + alpha * delta[i];
                    }

                    if (gradient_norm_at(next) < current_norm) {
                        current = next;
                        line_search_success = true;
                        break;
                    }
                    alpha *= Scalar(0.5L);
                }

                if (!line_search_success) {
                    for (size_t i = 0; i < variables.size(); ++i) {
                        current[variables[i]] = current[variables[i]] + alpha * delta[i];
                    }
                    if (alpha < precision::line_search_min_step<Scalar>()) break;
                }

                Scalar max_delta = Scalar(0);
                for (auto d : delta) max_delta = mymath::fmax(max_delta, mymath::abs(d * alpha));
                if (max_delta < precision::line_search_min_step<Scalar>()) break;
            }

            const Scalar grad_norm = gradient_norm_at(current);
            if (grad_norm < precision::gradient_convergence_threshold<Scalar>()) {
                add_critical_point(current);
            }
        } catch (const std::exception&) {
            // 忽略单个起始点的失败
        }
    }

    return critical_points;
}

// ============================================================================
// 线性系统求解
// ============================================================================

std::vector<Scalar> solve_linear_system(
    const std::vector<std::vector<Scalar>>& A,
    const std::vector<Scalar>& b) {

    const size_t n = A.size();
    if (n == 0 || b.size() != n) {
        throw std::runtime_error("invalid linear system dimensions");
    }

    // 创建增广矩阵
    std::vector<std::vector<Scalar>> aug(n, std::vector<Scalar>(n + 1));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) aug[i][j] = A[i][j];
        aug[i][n] = b[i];
    }

    // 高斯消元
    for (size_t col = 0; col < n; ++col) {
        // 寻找主元
        size_t max_row = col;
        for (size_t row = col + 1; row < n; ++row) {
            if (mymath::abs(aug[row][col]) > mymath::abs(aug[max_row][col])) {
                max_row = row;
            }
        }

        // 交换行
        std::swap(aug[col], aug[max_row]);

        // 检查奇异性
        if (mymath::abs(aug[col][col]) < precision::epsilon<Scalar>() * Scalar(1000)) {
            throw std::runtime_error("singular matrix in linear system");
        }

        // 消元
        for (size_t row = col + 1; row < n; ++row) {
            const Scalar factor = aug[row][col] / aug[col][col];
            for (size_t j = col; j <= n; ++j) {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    // 回代
    std::vector<Scalar> x(n);
    for (size_t i = n; i-- > 0;) {
        x[i] = aug[i][n];
        for (size_t j = i + 1; j < n; ++j) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }

    return x;
}

} // namespace analysis

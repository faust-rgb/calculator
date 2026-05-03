#include "symbolic/risch_algorithm.h"
#include "symbolic/symbolic_expression_internal.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

using namespace symbolic_expression_internal;

// 根查找
// ============================================================================

std::vector<RischAlgorithm::ComplexRoot> RischAlgorithm::find_all_roots(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& var_name) {

    std::vector<ComplexRoot> roots;

    std::size_t poly_degree = coeffs.size() - 1;
    if (poly_degree == 0) return roots;

    if (poly_degree == 1) {
        SymbolicExpression root = (make_negate(coeffs[0]) / coeffs[1]).simplify();
        roots.push_back(ComplexRoot::real(root));
        return roots;
    }

    // 尝试数值求解
    bool all_numeric = true;
    std::vector<double> num_coeffs;

    for (const auto& c : coeffs) {
        double val = 0.0;
        if (c.is_number(&val)) {
            num_coeffs.push_back(val);
        } else {
            all_numeric = false;
            break;
        }
    }

    if (all_numeric) {
        // 二次
        if (poly_degree == 2) {
            double a = num_coeffs[2], b = num_coeffs[1], c = num_coeffs[0];
            double delta = b * b - 4.0 * a * c;

            if (delta >= 0) {
                double sqrt_delta = std::sqrt(delta);
                roots.push_back(ComplexRoot::real(SymbolicExpression::number((-b + sqrt_delta) / (2.0 * a))));
                roots.push_back(ComplexRoot::real(SymbolicExpression::number((-b - sqrt_delta) / (2.0 * a))));
            } else {
                // 复数根：返回共轭对 (a + bi, a - bi) 作为单个条目
                double real_part = -b / (2.0 * a);
                double imag_part = std::sqrt(-delta) / (2.0 * a);
                // 存储为一个共轭对，表示 a+bi 和 a-bi
                roots.push_back(ComplexRoot::complex(
                    SymbolicExpression::number(real_part),
                    SymbolicExpression::number(imag_part),
                    true));
            }
            return roots;
        }

        // 三次（Cardano 公式）
        if (poly_degree == 3) {
            double a = num_coeffs[3], b = num_coeffs[2], c = num_coeffs[1], d = num_coeffs[0];

            double p = (3.0 * a * c - b * b) / (3.0 * a * a);
            double q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) / (27.0 * a * a * a);
            double disc = q * q / 4.0 + p * p * p / 27.0;

            if (disc > 0) {
                double u = std::sqrt(disc);
                double root1 = std::pow(-q / 2.0 + u, 1.0/3.0) + std::pow(-q / 2.0 - u, 1.0/3.0) - b / (3.0 * a);
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root1)));
            } else if (disc == 0) {
                double root1 = 3.0 * q / p - b / (3.0 * a);
                double root2 = -3.0 * q / (2.0 * p) - b / (3.0 * a);
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root1)));
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root2)));
            } else {
                double r = std::sqrt(-p * p * p / 27.0);
                double theta = std::acos(-q / (2.0 * r)) / 3.0;
                double root1 = 2.0 * std::pow(r, 1.0/3.0) * std::cos(theta) - b / (3.0 * a);
                double root2 = 2.0 * std::pow(r, 1.0/3.0) * std::cos(theta + 2.0 * M_PI / 3.0) - b / (3.0 * a);
                double root3 = 2.0 * std::pow(r, 1.0/3.0) * std::cos(theta + 4.0 * M_PI / 3.0) - b / (3.0 * a);
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root1)));
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root2)));
                roots.push_back(ComplexRoot::real(SymbolicExpression::number(root3)));
            }
            return roots;
        }

        // 更高次：使用牛顿法
        auto numeric_roots = find_numeric_roots_newton(num_coeffs);
        for (const auto& r : numeric_roots) {
            roots.push_back(ComplexRoot::real(r));
        }
        return roots;
    }

    // 符号系数：尝试整数和有理根
    auto int_roots = find_integer_roots(coeffs, var_name);
    for (const auto& r : int_roots) {
        roots.push_back(ComplexRoot::real(r));
    }

    if (roots.empty()) {
        auto rat_roots = find_rational_roots(coeffs, var_name);
        for (const auto& r : rat_roots) {
            roots.push_back(ComplexRoot::real(r));
        }
    }

    return roots;
}

std::vector<SymbolicExpression> RischAlgorithm::find_integer_roots(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& var_name) {

    std::vector<SymbolicExpression> roots;

    SymbolicExpression lc = coeffs.back();
    double lc_val = 1.0;
    if (!lc.is_number(&lc_val)) lc_val = 1.0;

    SymbolicExpression ct = coeffs.front();
    double ct_val = 0.0;
    if (!ct.is_number(&ct_val)) ct_val = 0.0;

    int max_search = 100;
    if (ct_val != 0.0 && std::abs(ct_val) < 1000) {
        max_search = static_cast<int>(std::abs(ct_val)) + 1;
    }

    for (int i = -max_search; i <= max_search; ++i) {
        if (i == 0 && coeffs.size() > 1 && !SymbolicPolynomial::coeff_is_zero(coeffs[0])) {
            continue;
        }

        SymbolicExpression val = SymbolicExpression::number(0.0);
        SymbolicExpression x_val = SymbolicExpression::number(i);
        SymbolicExpression power = SymbolicExpression::number(1.0);

        for (std::size_t j = 0; j < coeffs.size(); ++j) {
            val = (val + coeffs[j] * power).simplify();
            power = (power * x_val).simplify();
        }

        if (SymbolicPolynomial::coeff_is_zero(val)) {
            roots.push_back(SymbolicExpression::number(i));
        }
    }

    return roots;
}

std::vector<SymbolicExpression> RischAlgorithm::find_rational_roots(
    const std::vector<SymbolicExpression>& coeffs,
    const std::string& var_name) {

    std::vector<SymbolicExpression> roots;

    SymbolicExpression lc = coeffs.back();
    SymbolicExpression ct = coeffs.front();

    double lc_val = 1.0, ct_val = 0.0;
    if (!lc.is_number(&lc_val) || !ct.is_number(&ct_val)) {
        return roots;
    }

    auto get_divisors = [](double n) -> std::vector<int> {
        std::vector<int> divisors;
        int abs_n = static_cast<int>(std::abs(n) + 0.5);
        if (abs_n == 0) return divisors;
        for (int i = 1; i <= abs_n; ++i) {
            if (abs_n % i == 0) {
                divisors.push_back(i);
                divisors.push_back(-i);
            }
        }
        return divisors;
    };

    std::vector<int> p_divisors = get_divisors(ct_val);
    std::vector<int> q_divisors = get_divisors(lc_val);

    for (int p : p_divisors) {
        for (int q : q_divisors) {
            if (q == 0) continue;
            double rational = static_cast<double>(p) / static_cast<double>(q);

            SymbolicExpression val = SymbolicExpression::number(0.0);
            SymbolicExpression x_val = SymbolicExpression::number(rational);
            SymbolicExpression power = SymbolicExpression::number(1.0);

            for (std::size_t j = 0; j < coeffs.size(); ++j) {
                val = (val + coeffs[j] * power).simplify();
                power = (power * x_val).simplify();
            }

            if (SymbolicPolynomial::coeff_is_zero(val)) {
                roots.push_back(SymbolicExpression::number(rational));
            }
        }
    }

    return roots;
}

std::vector<SymbolicExpression> RischAlgorithm::find_numeric_roots_newton(
    const std::vector<double>& coeffs) {

    std::vector<SymbolicExpression> roots;
    if (coeffs.empty()) return roots;

    int deg = static_cast<int>(coeffs.size()) - 1;
    if (deg <= 0) return roots;

    auto eval_poly = [&coeffs](double x) -> double {
        double result = 0.0;
        double power = 1.0;
        for (double c : coeffs) {
            result += c * power;
            power *= x;
        }
        return result;
    };

    auto eval_deriv = [&coeffs, deg](double x) -> double {
        double result = 0.0;
        double power = 1.0;
        for (int i = 1; i <= deg; ++i) {
            result += i * coeffs[i] * power;
            power *= x;
        }
        return result;
    };

    std::vector<double> start_points = {-10.0, -5.0, -2.0, -1.0, -0.5, 0.5, 1.0, 2.0, 5.0, 10.0};

    for (double start : start_points) {
        double x = start;
        for (int iter = 0; iter < 50; ++iter) {
            double fx = eval_poly(x);
            double fpx = eval_deriv(x);

            if (std::abs(fpx) < 1e-12) break;

            double next_x = x - fx / fpx;

            if (std::abs(next_x - x) < 1e-10) {
                if (std::abs(eval_poly(next_x)) < 1e-9) {
                    bool already_found = false;
                    for (const auto& r : roots) {
                        double r_val = 0.0;
                        if (r.is_number(&r_val) && std::abs(r_val - next_x) < 1e-6) {
                            already_found = true;
                            break;
                        }
                    }
                    if (!already_found) {
                        roots.push_back(SymbolicExpression::number(next_x));
                    }
                }
                break;
            }
            x = next_x;
        }
    }

    return roots;
}

// ============================================================================
// Aberth-Ehrlich 方法找复数根 (Durand-Kerner 变体)
// ============================================================================

std::vector<std::pair<double, double>> RischAlgorithm::find_complex_roots_aberth(
    const std::vector<double>& coeffs) {

    std::vector<std::pair<double, double>> roots;
    if (coeffs.empty()) return roots;

    int n = static_cast<int>(coeffs.size()) - 1;
    if (n <= 0) return roots;

    // 使用 Durand-Kerner 方法
    auto eval_poly_complex = [&coeffs](const std::complex<double>& z) -> std::complex<double> {
        std::complex<double> result = 0.0;
        std::complex<double> power = 1.0;
        for (double c : coeffs) {
            result += c * power;
            power *= z;
        }
        return result;
    };

    // 初始猜测：在单位圆上均匀分布
    std::vector<std::complex<double>> z(n);
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * M_PI * i / n + 0.1;
        z[i] = std::complex<double>(std::cos(angle), std::sin(angle));
    }

    // Durand-Kerner 迭代
    for (int iter = 0; iter < 100; ++iter) {
        bool converged = true;
        for (int i = 0; i < n; ++i) {
            std::complex<double> p_z = eval_poly_complex(z[i]);
            std::complex<double> denom = 1.0;

            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    denom *= (z[i] - z[j]);
                }
            }

            if (std::abs(denom) < 1e-15) {
                denom = 1e-15;
            }

            std::complex<double> delta = p_z / denom;
            z[i] -= delta;

            if (std::abs(delta) > 1e-10) {
                converged = false;
            }
        }

        if (converged) break;
    }

    // 收集根并合并共轭对
    std::set<int> used;
    for (int i = 0; i < n; ++i) {
        if (used.count(i)) continue;

        double re = z[i].real();
        double im = z[i].imag();

        // 检查是否为实根
        if (std::abs(im) < 1e-9) {
            roots.push_back({re, 0.0});
            used.insert(i);
        } else {
            // 寻找共轭根
            bool found_conjugate = false;
            for (int j = i + 1; j < n; ++j) {
                if (!used.count(j) &&
                    std::abs(z[j].real() - re) < 1e-9 &&
                    std::abs(z[j].imag() + im) < 1e-9) {
                    // 共轭对
                    roots.push_back({re, std::abs(im)});
                    used.insert(i);
                    used.insert(j);
                    found_conjugate = true;
                    break;
                }
            }
            if (!found_conjugate) {
                roots.push_back({re, im});
                used.insert(i);
            }
        }
    }

    return roots;
}

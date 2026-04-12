/**
 * @file polynomial.cpp
 * @brief 多项式运算实现
 *
 * 实现多项式的基本运算（加减乘除）和实根查找。
 * 使用系数向量表示多项式，索引 i 对应 x^i 的系数。
 */

#include "polynomial.h"

#include "mymath.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace {

/** @brief 多项式运算的数值精度阈值 */
constexpr double kPolynomialEps = 1e-10;

/** @brief 根隔离过程的精度阈值 */
constexpr double kRootIsolationEps = 1e-8;

void trim_trailing_zeros(std::vector<double>* coefficients) {
    // 多项式内部采用“低次到高次”的系数表示：
    // coefficients[i] 对应 x^i 的系数。
    //
    // 例如：
    //   x^2 + 2x + 1  ->  [1, 2, 1]
    //
    // 算法过程中经常会出现最高项被消成接近 0 的情况，
    // 如果不在这里统一裁掉尾部 0，后面的次数判断、除法和求根都会出错。
    while (coefficients->size() > 1 &&
           mymath::is_near_zero(coefficients->back(), kPolynomialEps)) {
        coefficients->pop_back();
    }
    if (coefficients->empty()) {
        coefficients->push_back(0.0);
    }
}

std::string format_coefficient(double value) {
    if (mymath::is_integer(value, 1e-10)) {
        long long rounded =
            static_cast<long long>(value >= 0.0 ? value + 0.5 : value - 0.5);
        return std::to_string(rounded);
    }

    std::string text = std::to_string(value);
    while (!text.empty() && text.back() == '0') {
        text.pop_back();
    }
    if (!text.empty() && text.back() == '.') {
        text.pop_back();
    }
    return text;
}

double polynomial_evaluate(const std::vector<double>& coefficients, double x) {
    // 这里使用 Horner 形式求值：
    //   a0 + a1*x + a2*x^2 + ...
    // 可以改写成
    //   (...((an*x + a(n-1))*x + ...) * x + a0)
    //
    // 这样既减少乘法次数，也比直接逐项求幂更稳定。
    double result = 0.0;
    for (std::size_t i = coefficients.size(); i > 0; --i) {
        result = result * x + coefficients[i - 1];
    }
    return result;
}

std::vector<double> polynomial_derivative_coefficients(
    const std::vector<double>& coefficients) {
    if (coefficients.size() <= 1) {
        return {0.0};
    }

    std::vector<double> derivative(coefficients.size() - 1, 0.0);
    for (std::size_t i = 1; i < coefficients.size(); ++i) {
        derivative[i - 1] = coefficients[i] * static_cast<double>(i);
    }
    trim_trailing_zeros(&derivative);
    return derivative;
}

double polynomial_root_bound(const std::vector<double>& coefficients) {
    // 使用一个简单的 Cauchy 型根界：
    // 所有实根都会落在 [-B, B]，其中 B = 1 + max |ai/an|。
    //
    // 这样就可以把“无界的实轴搜索”收缩到一个有限区间，
    // 后续只需要在这个区间里结合导函数临界点做分段即可。
    const double leading = coefficients.back();
    double bound = 0.0;
    for (std::size_t i = 0; i + 1 < coefficients.size(); ++i) {
        const double ratio = mymath::abs(coefficients[i] / leading);
        if (ratio > bound) {
            bound = ratio;
        }
    }
    return 1.0 + bound;
}

double bisect_root(const std::vector<double>& coefficients, double left, double right) {
    // 进入这个函数前，调用方已经保证 [left, right] 两端函数值异号。
    // 因此直接用二分法做“可靠收敛”的根逼近。
    double left_value = polynomial_evaluate(coefficients, left);

    for (int i = 0; i < 100; ++i) {
        const double mid = (left + right) * 0.5;
        const double mid_value = polynomial_evaluate(coefficients, mid);
        if (mymath::is_near_zero(mid_value, kRootIsolationEps) ||
            mymath::abs(right - left) <= kRootIsolationEps) {
            return mid;
        }

        if ((left_value < 0.0 && mid_value > 0.0) ||
            (left_value > 0.0 && mid_value < 0.0)) {
            right = mid;
        } else {
            left = mid;
            left_value = mid_value;
        }
    }

    return (left + right) * 0.5;
}

void add_unique_root(std::vector<double>* roots, double candidate) {
    // 求根过程会同时从“导函数临界点命中”以及“分段变号”两条路径得到根。
    // 同一个根可能被重复发现，因此统一做一次近似去重。
    for (double existing : *roots) {
        if (mymath::abs(existing - candidate) <= 1e-6) {
            return;
        }
    }
    roots->push_back(candidate);
}

}  // namespace

std::vector<double> polynomial_add(const std::vector<double>& lhs,
                                   const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] += lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] += rhs[i];
    }
    trim_trailing_zeros(&result);
    return result;
}

std::vector<double> polynomial_subtract(const std::vector<double>& lhs,
                                        const std::vector<double>& rhs) {
    const std::size_t size = lhs.size() > rhs.size() ? lhs.size() : rhs.size();
    std::vector<double> result(size, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i] += lhs[i];
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result[i] -= rhs[i];
    }
    trim_trailing_zeros(&result);
    return result;
}

std::vector<double> polynomial_multiply(const std::vector<double>& lhs,
                                        const std::vector<double>& rhs) {
    std::vector<double> result(lhs.size() + rhs.size() - 1, 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            result[i + j] += lhs[i] * rhs[j];
        }
    }
    trim_trailing_zeros(&result);
    return result;
}

PolynomialDivisionResult polynomial_divide(const std::vector<double>& dividend,
                                           const std::vector<double>& divisor) {
    std::vector<double> normalized_dividend = dividend;
    std::vector<double> normalized_divisor = divisor;
    trim_trailing_zeros(&normalized_dividend);
    trim_trailing_zeros(&normalized_divisor);

    if (normalized_divisor.size() == 1 &&
        mymath::is_near_zero(normalized_divisor[0], kPolynomialEps)) {
        throw std::runtime_error("polynomial divisor cannot be zero");
    }

    if (normalized_dividend.size() < normalized_divisor.size()) {
        return {{0.0}, normalized_dividend};
    }

    std::vector<double> quotient(
        normalized_dividend.size() - normalized_divisor.size() + 1, 0.0);
    std::vector<double> remainder = normalized_dividend;

    // 这里实现的是标准多项式长除法：
    // 每一轮都拿“当前余式最高项 / 除式最高项”得到一个商项，
    // 再把对应倍数从余式里消掉，直到余式次数低于除式次数。
    while (remainder.size() >= normalized_divisor.size() &&
           !(remainder.size() == 1 && mymath::is_near_zero(remainder[0], kPolynomialEps))) {
        const std::size_t degree_diff =
            remainder.size() - normalized_divisor.size();
        const double factor =
            remainder.back() / normalized_divisor.back();
        quotient[degree_diff] = factor;

        for (std::size_t i = 0; i < normalized_divisor.size(); ++i) {
            remainder[degree_diff + i] -= factor * normalized_divisor[i];
        }
        trim_trailing_zeros(&remainder);
        if (remainder.size() < normalized_divisor.size()) {
            break;
        }
    }

    trim_trailing_zeros(&quotient);
    trim_trailing_zeros(&remainder);
    return {quotient, remainder};
}

std::vector<double> polynomial_real_roots(const std::vector<double>& coefficients) {
    std::vector<double> normalized = coefficients;
    trim_trailing_zeros(&normalized);

    if (normalized.size() <= 1) {
        throw std::runtime_error("constant polynomial does not have isolated roots");
    }

    if (normalized.size() == 2) {
        return {-normalized[0] / normalized[1]};
    }

    // 核心思路：
    // 1. 先求导函数的全部实根，这些点就是原多项式在实轴上的临界点
    // 2. 把根界区间 [-B, B] 用这些临界点切成若干段
    // 3. 在每一段里，多项式单调，因此最多只有一个实根
    // 4. 只要段两端变号，就可以安全地二分定位这个根
    //
    // 这种“递归导数 + 分段二分”的方式比对整个区间做稠密扫描更稳定，
    // 对当前只返回实根的需求也足够直接。
    const std::vector<double> derivative =
        polynomial_derivative_coefficients(normalized);
    std::vector<double> critical_points = polynomial_real_roots(derivative);
    std::sort(critical_points.begin(), critical_points.end());

    const double bound = polynomial_root_bound(normalized);
    std::vector<double> points;
    points.push_back(-bound);
    for (double point : critical_points) {
        points.push_back(point);
    }
    points.push_back(bound);

    std::vector<double> roots;
    for (double point : critical_points) {
        // 临界点本身如果恰好落在 x 轴上，往往对应重根。
        // 这类根不一定伴随普通的“左右变号”，所以需要单独补录。
        if (mymath::is_near_zero(polynomial_evaluate(normalized, point), 1e-6)) {
            add_unique_root(&roots, point);
        }
    }

    for (std::size_t i = 1; i < points.size(); ++i) {
        const double left = points[i - 1];
        const double right = points[i];
        const double left_value = polynomial_evaluate(normalized, left);
        const double right_value = polynomial_evaluate(normalized, right);

        if (mymath::is_near_zero(left_value, 1e-6)) {
            add_unique_root(&roots, left);
            continue;
        }
        if (mymath::is_near_zero(right_value, 1e-6)) {
            add_unique_root(&roots, right);
            continue;
        }

        if ((left_value < 0.0 && right_value > 0.0) ||
            (left_value > 0.0 && right_value < 0.0)) {
            add_unique_root(&roots, bisect_root(normalized, left, right));
        }
    }

    std::sort(roots.begin(), roots.end());
    return roots;
}

std::string polynomial_to_string(const std::vector<double>& coefficients,
                                 const std::string& variable_name) {
    // 这里负责把内部系数数组重新格式化成人可读的多项式字符串，
    // 例如 [1, 0, -3] -> "-3 * x ^ 2 + 1"。
    //
    // 输出格式保持和项目里现有符号表达式风格一致，
    // 方便 poly_*、roots 等命令直接复用。
    std::vector<double> normalized = coefficients;
    trim_trailing_zeros(&normalized);
    if (normalized.size() == 1 && mymath::is_near_zero(normalized[0], kPolynomialEps)) {
        return "0";
    }

    std::string result;
    bool first = true;
    for (std::size_t index = normalized.size(); index > 0; --index) {
        const std::size_t degree = index - 1;
        const double coefficient = normalized[degree];
        if (mymath::is_near_zero(coefficient, kPolynomialEps)) {
            continue;
        }

        const bool negative = coefficient < 0.0;
        const double abs_value = negative ? -coefficient : coefficient;

        std::string term;
        if (degree == 0) {
            term = format_coefficient(abs_value);
        } else {
            if (!mymath::is_near_zero(abs_value - 1.0, kPolynomialEps)) {
                term += format_coefficient(abs_value) + " * ";
            }
            term += variable_name;
            if (degree > 1) {
                term += " ^ " + std::to_string(degree);
            }
        }

        if (first) {
            result += negative ? "-" + term : term;
            first = false;
        } else {
            result += negative ? " - " + term : " + " + term;
        }
    }

    return result.empty() ? "0" : result;
}

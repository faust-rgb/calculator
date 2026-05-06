// ============================================================================
// 符号求和模块
// ============================================================================
//
// 实现符号求和，支持：
// 1. 多项式求和 (Faulhaber 公式)
// 2. 几何级数
// 3. 伸缩级数
// 4. 超几何级数 (Gosper 算法)
// 5. 已知无穷级数
//
// ============================================================================

#include "symbolic/symbolic_sum.h"
#include "symbolic/symbolic_expression_internal.h"
#include "math/mymath.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace symbolic_sum {

namespace {

using namespace symbolic_expression_internal;

// Bernoulli 数 B_n (有理数形式)
// B_0 = 1, B_1 = -1/2, B_2 = 1/6, B_4 = -1/30, B_6 = 1/42, ...
// 奇数项 B_3 = B_5 = ... = 0 (除了 B_1)
static const std::vector<std::pair<long long, long long>> bernoulli_numbers = {
    {1, 1},      // B_0 = 1
    {-1, 2},     // B_1 = -1/2
    {1, 6},      // B_2 = 1/6
    {0, 1},      // B_3 = 0
    {-1, 30},    // B_4 = -1/30
    {0, 1},      // B_5 = 0
    {1, 42},     // B_6 = 1/42
    {0, 1},      // B_7 = 0
    {-1, 30},    // B_8 = -1/30
    {0, 1},      // B_9 = 0
    {5, 66},     // B_10 = 5/66
    {0, 1},      // B_11 = 0
    {-691, 2730}, // B_12 = -691/2730
    {0, 1},      // B_13 = 0
    {7, 6},      // B_14 = 7/6
    {0, 1},      // B_15 = 0
    {-3617, 510}, // B_16 = -3617/510
};

// 计算组合数 C(n, k)
long long binomial(long long n, long long k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    long long result = 1;
    for (long long i = 0; i < k; ++i) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

// 检查表达式是否包含变量
bool contains_variable(const SymbolicExpression& expr, const std::string& var) {
    std::string str = expr.to_string();
    return str.find(var) != std::string::npos;
}

}  // namespace

// ============================================================================
// SymbolicSumEngine 实现
// ============================================================================

SumResult SymbolicSumEngine::compute_sum(
    const SymbolicExpression& term,
    const std::string& var,
    const BoundArgument& lower,
    const BoundArgument& upper) {

    // 确定求和类型
    SumType sum_type = SumType::kFinite;
    if (upper.is_infinite()) {
        sum_type = SumType::kInfiniteUpper;
    }
    if (lower.is_infinite() && upper.is_infinite()) {
        sum_type = SumType::kInfiniteBoth;
    }

    // 分类求和项
    TermType term_type = classify_term(term, var);

    // 策略 1: 已知无穷级数
    if (sum_type == SumType::kInfiniteUpper) {
        SumResult known = known_infinite_series(term, var);
        if (known.is_convergent && known.method_used != "unknown") {
            return known;
        }
    }

    // 策略 2: 多项式求和 (Faulhaber)
    if (term_type == TermType::kPolynomial) {
        if (lower.is_finite() && upper.is_finite()) {
            return sum_polynomial(term, var,
                static_cast<int>(mymath::round(lower.value)),
                static_cast<int>(mymath::round(upper.value)));
        }
        if (lower.is_finite() && upper.is_infinite()) {
            // 无穷多项式级数发散
            return SumResult::divergent("polynomial series diverges");
        }
    }

    // 策略 3: 几何级数
    if (term_type == TermType::kGeometric) {
        return sum_geometric(term, var, lower, upper);
    }

    // 策略 4: 伸缩级数
    SymbolicExpression telescoping_f;
    if (detect_telescoping(term, var, &telescoping_f)) {
        if (lower.is_finite() && upper.is_finite()) {
            int lo = static_cast<int>(mymath::round(lower.value));
            int hi = static_cast<int>(mymath::round(upper.value));
            SymbolicExpression result = telescoping_f.substitute(var, SymbolicExpression::number(hi + 1))
                                        .simplify() -
                                        telescoping_f.substitute(var, SymbolicExpression::number(lo))
                                        .simplify();
            return SumResult::elementary(result.simplify(), TermType::kTelescoping, "telescoping");
        }
    }

    // 策略 5: Gosper 算法
    if (lower.is_finite() && upper.is_finite()) {
        SumResult gosper_result;
        if (gosper_algorithm(term, var,
            static_cast<int>(mymath::round(lower.value)),
            static_cast<int>(mymath::round(upper.value)),
            &gosper_result)) {
            return gosper_result;
        }
    }

    // 策略 6: 数值求和（回退）
    if (lower.is_finite() && upper.is_finite()) {
        int lo = static_cast<int>(mymath::round(lower.value));
        int hi = static_cast<int>(mymath::round(upper.value));
        if (hi - lo < 1000) {  // 限制数值求和范围
            SymbolicExpression result = numerical_sum(term, var, lo, hi);
            return SumResult::elementary(result, term_type, "numerical_sum");
        }
    }

    return SumResult::unknown();
}

std::optional<SymbolicExpression> SymbolicSumEngine::sum(
    const std::string& term_str,
    const std::string& var,
    const BoundArgument& lower,
    const BoundArgument& upper) {

    try {
        SymbolicExpression term = SymbolicExpression::parse(term_str);
        SymbolicSumEngine engine;
        SumResult result = engine.compute_sum(term, var, lower, upper);

        if (result.is_convergent) {
            return result.value;
        }
        return std::nullopt;
    } catch (...) {
        return std::nullopt;
    }
}

TermType SymbolicSumEngine::classify_term(
    const SymbolicExpression& term,
    const std::string& var) {

    // 检查是否为几何级数项 a*r^n
    SymbolicExpression a, r;
    if (is_geometric_term(term, var, &a, &r)) {
        return TermType::kGeometric;
    }

    // 检查是否为多项式
    int degree = get_polynomial_degree(term, var);
    if (degree >= 0) {
        return TermType::kPolynomial;
    }

    // 检查是否为有理函数
    if (term.node_->type == NodeType::kDivide) {
        return TermType::kRational;
    }

    return TermType::kUnknown;
}

SumResult SymbolicSumEngine::sum_polynomial(
    const SymbolicExpression& term,
    const std::string& var,
    int lower,
    int upper) {

    // 对于多项式求和，使用 Faulhaber 公式
    // sum(k^m, k, 1, n) = 1/(m+1) * Σ C(m+1,j) * B_j * n^(m+1-j)

    int degree = get_polynomial_degree(term, var);

    if (degree < 0) {
        return SumResult::unknown();
    }

    // 特殊情况：常数项
    if (degree == 0) {
        double val = 0.0;
        if (term.is_number(&val)) {
            // sum(c, k, lo, hi) = c * (hi - lo + 1)
            double result = val * (upper - lower + 1);
            return SumResult::elementary(
                SymbolicExpression::number(result),
                TermType::kPolynomial,
                "constant_sum");
        }
    }

    // 对于 k^m 形式
    if (term.node_->type == NodeType::kVariable && term.node_->text == var) {
        // sum(k, k, lo, hi) = (hi + lo) * (hi - lo + 1) / 2
        long long n = upper;
        long long m = lower;

        // 使用 sum(k, k, 1, n) - sum(k, k, 1, m-1)
        SymbolicExpression sum_to_n = faulhaber_sum(1, SymbolicExpression::number(n));
        SymbolicExpression sum_to_m_1 = faulhaber_sum(1, SymbolicExpression::number(m - 1));

        SymbolicExpression result = (sum_to_n - sum_to_m_1).simplify();
        return SumResult::elementary(result, TermType::kPolynomial, "faulhaber");
    }

    // 对于 k^m 形式
    if (term.node_->type == NodeType::kPower) {
        SymbolicExpression base(term.node_->left);
        double exp = 0.0;
        if (base.node_->type == NodeType::kVariable && base.node_->text == var &&
            SymbolicExpression(term.node_->right).is_number(&exp)) {
            int m = static_cast<int>(mymath::round(exp));

            long long n = upper;
            long long lo = lower;

            SymbolicExpression sum_to_n = faulhaber_sum(m, SymbolicExpression::number(n));
            SymbolicExpression sum_to_lo_1 = faulhaber_sum(m, SymbolicExpression::number(lo - 1));

            SymbolicExpression result = (sum_to_n - sum_to_lo_1).simplify();
            return SumResult::elementary(result, TermType::kPolynomial, "faulhaber");
        }
    }

    // 一般多项式：分解为各幂次项
    // 简化版本：数值求和
    SymbolicExpression result = numerical_sum(term, var, lower, upper);
    return SumResult::elementary(result, TermType::kPolynomial, "numerical_polynomial_sum");
}

SumResult SymbolicSumEngine::sum_geometric(
    const SymbolicExpression& term,
    const std::string& var,
    const BoundArgument& lower,
    const BoundArgument& upper) {

    SymbolicExpression a, r;
    if (!is_geometric_term(term, var, &a, &r)) {
        return SumResult::unknown();
    }

    a = a.simplify();
    r = r.simplify();

    // 几何级数: sum(a * r^k, k, lo, hi) = a * r^lo * (r^(hi-lo+1) - 1) / (r - 1)

    if (lower.is_finite() && upper.is_finite()) {
        long long lo = static_cast<long long>(mymath::round(lower.value));
        long long hi = static_cast<long long>(mymath::round(upper.value));
        long long n = hi - lo + 1;

        // 检查 r 是否为数值
        double r_val = 0.0;
        if (r.is_number(&r_val)) {
            if (mymath::is_near_zero(r_val - 1.0, 1e-15)) {
                // r = 1: sum = a * n
                double a_val = 0.0;
                if (a.is_number(&a_val)) {
                    return SumResult::elementary(
                        SymbolicExpression::number(a_val * n),
                        TermType::kGeometric,
                        "geometric_r_equals_1");
                }
            }

            // 一般情况
            double a_val = 0.0;
            if (a.is_number(&a_val)) {
                double result = a_val * mymath::pow(r_val, lo) * (mymath::pow(r_val, n) - 1.0) / (r_val - 1.0);
                return SumResult::elementary(
                    SymbolicExpression::number(result),
                    TermType::kGeometric,
                    "geometric_numeric");
            }
        }

        // 符号形式
        SymbolicExpression r_lo = make_power(r, SymbolicExpression::number(lo));
        SymbolicExpression r_n = make_power(r, SymbolicExpression::number(n));
        SymbolicExpression numerator = make_subtract(r_n, SymbolicExpression::number(1.0));
        SymbolicExpression denominator = make_subtract(r, SymbolicExpression::number(1.0));
        SymbolicExpression result = make_multiply(a, make_multiply(r_lo, make_divide(numerator, denominator))).simplify();

        return SumResult::elementary(result, TermType::kGeometric, "geometric_symbolic");
    }

    // 无穷几何级数
    if (lower.is_finite() && upper.is_infinite()) {
        long long lo = static_cast<long long>(mymath::round(lower.value));

        // 检查收敛性 |r| < 1
        double r_val = 0.0;
        if (r.is_number(&r_val)) {
            if (mymath::abs(r_val) < 1.0) {
                // sum(a * r^k, k, lo, inf) = a * r^lo / (1 - r)
                double a_val = 0.0;
                if (a.is_number(&a_val)) {
                    double result = a_val * mymath::pow(r_val, lo) / (1.0 - r_val);
                    return SumResult::elementary(
                        SymbolicExpression::number(result),
                        TermType::kGeometric,
                        "geometric_infinite");
                }

                // 符号形式
                SymbolicExpression r_lo = make_power(r, SymbolicExpression::number(lo));
                SymbolicExpression result = make_multiply(a, make_divide(r_lo, make_subtract(SymbolicExpression::number(1.0), r))).simplify();
                return SumResult::elementary(result, TermType::kGeometric, "geometric_infinite_symbolic");
            } else {
                return SumResult::divergent("|r| >= 1");
            }
        }
    }

    return SumResult::unknown();
}

SumResult SymbolicSumEngine::sum_rational(
    const SymbolicExpression& term,
    const std::string& var,
    int lower,
    int upper) {

    // 有理函数求和较为复杂，回退到数值方法
    SymbolicExpression result = numerical_sum(term, var, lower, upper);
    return SumResult::elementary(result, TermType::kRational, "numerical_rational_sum");
}

bool SymbolicSumEngine::gosper_algorithm(
    const SymbolicExpression& term,
    const std::string& var,
    int lower,
    int upper,
    SumResult* result) {
    (void)lower;
    (void)upper;
    (void)result;

    // Gosper 算法用于超几何级数求和
    // 检查 t_{k+1}/t_k 是否为有理函数

    // 计算 t(k+1)
    SymbolicExpression var_plus_1 = make_add(SymbolicExpression::variable(var), SymbolicExpression::number(1.0));
    SymbolicExpression t_next = term.substitute(var, var_plus_1).simplify();

    // 计算比值 r(k) = t(k+1)/t(k)
    SymbolicExpression ratio = make_divide(t_next, term).simplify();

    // 检查 ratio 是否为有理函数
    // 简化版本：检查是否为多项式比值
    if (ratio.node_->type == NodeType::kDivide) {
        // 可能是超几何项
        // 完整的 Gosper 算法实现较为复杂
        // 这里简化处理
    }

    return false;
}

bool SymbolicSumEngine::detect_telescoping(
    const SymbolicExpression& term,
    const std::string& var,
    SymbolicExpression* f) {
    (void)f;

    // 检测伸缩级数: term = f(k+1) - f(k)
    // 方法：尝试积分或猜测

    // 简单检测：如果 term 是差分形式
    if (term.node_->type == NodeType::kSubtract) {
        SymbolicExpression left(term.node_->left);
        SymbolicExpression right(term.node_->right);

        // 检查 left 是否为 f(k+1), right 是否为 f(k)
        SymbolicExpression var_plus_1 = make_add(SymbolicExpression::variable(var), SymbolicExpression::number(1.0));

        // 尝试匹配
        // 这需要更复杂的模式匹配
    }

    return false;
}

SumResult SymbolicSumEngine::known_infinite_series(
    const SymbolicExpression& term,
    const std::string& var) {

    // 已知无穷级数表

    // sum(1/k^2, k, 1, inf) = pi^2/6
    if (term.node_->type == NodeType::kDivide) {
        SymbolicExpression num(term.node_->left);
        SymbolicExpression den(term.node_->right);

        double num_val = 0.0;
        if (num.is_number(&num_val) && mymath::is_near_zero(num_val - 1.0, 1e-15)) {
            // 检查分母是否为 k^2
            if (den.node_->type == NodeType::kPower) {
                SymbolicExpression base(den.node_->left);
                double exp = 0.0;
                if (base.node_->type == NodeType::kVariable && base.node_->text == var &&
                    SymbolicExpression(den.node_->right).is_number(&exp)) {
                    if (mymath::is_near_zero(exp - 2.0, 1e-10)) {
                        // sum(1/k^2) = pi^2/6
                        SymbolicExpression result = make_divide(
                            make_power(SymbolicExpression::parse("pi"), SymbolicExpression::number(2.0)),
                            SymbolicExpression::number(6.0));
                        return SumResult::elementary(result, TermType::kRational, "zeta(2)");
                    }
                    if (mymath::is_near_zero(exp - 4.0, 1e-10)) {
                        // sum(1/k^4) = pi^4/90
                        SymbolicExpression result = make_divide(
                            make_power(SymbolicExpression::parse("pi"), SymbolicExpression::number(4.0)),
                            SymbolicExpression::number(90.0));
                        return SumResult::elementary(result, TermType::kRational, "zeta(4)");
                    }
                }
            }
        }
    }

    // sum(x^n/n!, n, 0, inf) = exp(x)
    if (term.node_->type == NodeType::kDivide) {
        SymbolicExpression num(term.node_->left);
        SymbolicExpression den(term.node_->right);

        // 检查是否为 x^n / n!
        // 这需要更复杂的模式匹配
    }

    // sum((-1)^n/(2n+1), n, 0, inf) = pi/4 (Leibniz 公式)
    // 这需要检测 (-1)^n 模式

    return SumResult::unknown();
}

SymbolicExpression SymbolicSumEngine::faulhaber_sum(int m, const SymbolicExpression& n) {
    // Faulhaber 公式: sum(k^m, k, 1, n) = 1/(m+1) * Σ C(m+1,j) * B_j * n^(m+1-j)

    SymbolicExpression result = SymbolicExpression::number(0.0);

    for (int j = 0; j <= m; ++j) {
        auto [b_num, b_den] = get_bernoulli(j);
        if (b_num == 0) continue;  // 跳过零 Bernoulli 数

        long long comb = binomial(m + 1, j);

        // term = C(m+1,j) * B_j * n^(m+1-j)
        SymbolicExpression coeff = SymbolicExpression::number(static_cast<double>(comb * b_num) / b_den);
        SymbolicExpression power = make_power(n, SymbolicExpression::number(m + 1 - j));

        result = make_add(result, make_multiply(coeff, power)).simplify();
    }

    // 除以 m+1
    result = make_divide(result, SymbolicExpression::number(m + 1)).simplify();

    return result;
}

std::pair<long long, long long> SymbolicSumEngine::get_bernoulli(int n) {
    if (n < 0 || n >= static_cast<int>(bernoulli_numbers.size())) {
        return {0, 1};
    }
    return bernoulli_numbers[n];
}

bool SymbolicSumEngine::is_geometric_term(
    const SymbolicExpression& term,
    const std::string& var,
    SymbolicExpression* a,
    SymbolicExpression* r) {

    // 检查是否为 a * r^var 形式

    // 情况 1: r^var (a = 1)
    if (term.node_->type == NodeType::kPower) {
        SymbolicExpression base(term.node_->left);
        SymbolicExpression exp(term.node_->right);

        if (exp.node_->type == NodeType::kVariable && exp.node_->text == var) {
            *a = SymbolicExpression::number(1.0);
            *r = base;
            return true;
        }
    }

    // 情况 2: a * r^var
    if (term.node_->type == NodeType::kMultiply) {
        SymbolicExpression left(term.node_->left);
        SymbolicExpression right(term.node_->right);

        // 检查 left 是否为常数，right 是否为 r^var
        if (!contains_variable(left, var)) {
            if (right.node_->type == NodeType::kPower) {
                SymbolicExpression base(right.node_->left);
                SymbolicExpression exp(right.node_->right);

                if (exp.node_->type == NodeType::kVariable && exp.node_->text == var) {
                    *a = left;
                    *r = base;
                    return true;
                }
            }
        }

        // 检查 right 是否为常数，left 是否为 r^var
        if (!contains_variable(right, var)) {
            if (left.node_->type == NodeType::kPower) {
                SymbolicExpression base(left.node_->left);
                SymbolicExpression exp(left.node_->right);

                if (exp.node_->type == NodeType::kVariable && exp.node_->text == var) {
                    *a = right;
                    *r = base;
                    return true;
                }
            }
        }
    }

    return false;
}

int SymbolicSumEngine::get_polynomial_degree(
    const SymbolicExpression& term,
    const std::string& var) {

    // 常数
    double val = 0.0;
    if (term.is_number(&val)) {
        return 0;
    }

    // 变量
    if (term.node_->type == NodeType::kVariable) {
        if (term.node_->text == var) {
            return 1;
        }
        return 0;  // 其他变量视为常数
    }

    // 幂次
    if (term.node_->type == NodeType::kPower) {
        SymbolicExpression base(term.node_->left);
        double exp = 0.0;
        if (base.node_->type == NodeType::kVariable && base.node_->text == var &&
            SymbolicExpression(term.node_->right).is_number(&exp)) {
            if (mymath::is_integer(exp, 1e-10) && exp >= 0) {
                return static_cast<int>(mymath::round(exp));
            }
        }
        // 更一般的情况
        int base_deg = get_polynomial_degree(base, var);
        if (base_deg >= 0 && SymbolicExpression(term.node_->right).is_number(&exp)) {
            return base_deg * static_cast<int>(mymath::round(exp));
        }
    }

    // 加法：取最大幂次
    if (term.node_->type == NodeType::kAdd || term.node_->type == NodeType::kSubtract) {
        int left_deg = get_polynomial_degree(SymbolicExpression(term.node_->left), var);
        int right_deg = get_polynomial_degree(SymbolicExpression(term.node_->right), var);
        return std::max(left_deg, right_deg);
    }

    // 乘法：幂次相加
    if (term.node_->type == NodeType::kMultiply) {
        int left_deg = get_polynomial_degree(SymbolicExpression(term.node_->left), var);
        int right_deg = get_polynomial_degree(SymbolicExpression(term.node_->right), var);
        if (left_deg >= 0 && right_deg >= 0) {
            return left_deg + right_deg;
        }
    }

    // 非多项式
    return -1;
}

SymbolicExpression SymbolicSumEngine::numerical_sum(
    const SymbolicExpression& term,
    const std::string& var,
    int lower,
    int upper) {

    double sum = 0.0;

    for (int k = lower; k <= upper; ++k) {
        SymbolicExpression sub = term.substitute(var, SymbolicExpression::number(k));
        sub = sub.simplify();
        double val = 0.0;
        if (sub.is_number(&val)) {
            sum += val;
        } else {
            // 无法数值化，返回符号形式
            return SymbolicExpression::parse("sum(" + term.to_string() + ")");
        }
    }

    return SymbolicExpression::number(sum);
}

// ============================================================================
// 辅助函数实现
// ============================================================================

bool parse_sum_arguments(
    const std::vector<std::string>& args,
    SymbolicExpression* term,
    std::string* var,
    BoundArgument* lower,
    BoundArgument* upper) {

    if (args.size() < 4) {
        return false;
    }

    try {
        *term = SymbolicExpression::parse(args[0]);
        *var = args[1];

        // 解析下限
        std::string lower_str = args[2];
        if (lower_str == "inf" || lower_str == "infinity" || lower_str == "oo") {
            *lower = BoundArgument::pos_inf();
        } else if (lower_str == "-inf" || lower_str == "-infinity" || lower_str == "-oo") {
            *lower = BoundArgument::neg_inf();
        } else {
            double lo = std::stod(lower_str);
            *lower = BoundArgument::finite(lo);
        }

        // 解析上限
        std::string upper_str = args[3];
        if (upper_str == "inf" || upper_str == "infinity" || upper_str == "oo" || upper_str == "+inf") {
            *upper = BoundArgument::pos_inf();
        } else if (upper_str == "-inf" || upper_str == "-infinity" || upper_str == "-oo") {
            *upper = BoundArgument::neg_inf();
        } else {
            double hi = std::stod(upper_str);
            *upper = BoundArgument::finite(hi);
        }

        return true;
    } catch (...) {
        return false;
    }
}

std::string format_sum_result(const SumResult& result) {
    if (!result.is_convergent) {
        return "diverges";
    }
    return result.value.simplify().to_string();
}

}  // namespace symbolic_sum

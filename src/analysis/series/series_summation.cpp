#include "analysis/series/series_summation.h"
#include "analysis/modules/series_module.h"
#include "core/scalar_type.h"
#include "core/format_utils.h"
#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_expression_internal.h"
#include "math/mymath.h"
#include "math/helpers/integer_helpers.h"
#include "statistics/probability.h"
#include "expression_utils.h"
#include "string_utils.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace series_ops {
namespace summation {

bool detect_geometric_ratio_symbolic(const SymbolicExpression& summand,
                                      const std::string& index_name,
                                      Scalar* coefficient,
                                      Scalar* ratio) {
    SymbolicExpression n_plus_1 = SymbolicExpression::parse("(" + index_name + ") + 1");
    SymbolicExpression next_term = summand.substitute(index_name, n_plus_1);

    SymbolicExpression ratio_expr = (next_term / summand).simplify();

    if (ratio_expr.is_number(ratio)) {
        SymbolicExpression n_zero = SymbolicExpression::number(0.0L);
        SymbolicExpression term_at_zero = summand.substitute(index_name, n_zero).simplify();
        if (term_at_zero.is_number(coefficient)) {
            return true;
        }
        SymbolicExpression n_one = SymbolicExpression::number(1.0L);
        SymbolicExpression term_at_one = summand.substitute(index_name, n_one).simplify();
        Scalar val_one = 0.0L;
        if (term_at_one.is_number(&val_one) && !mymath::is_near_zero(*ratio, 1e-12)) {
            *coefficient = val_one / *ratio;
            return true;
        }
    }

    return false;
}

bool detect_arith_geo_series(const SymbolicExpression& summand,
                              const std::string& index_name,
                              std::vector<Scalar>* poly_coeffs,
                              Scalar* ratio) {
    SymbolicExpression n_plus_1 = SymbolicExpression::parse("(" + index_name + ") + 1");
    SymbolicExpression next_term = summand.substitute(index_name, n_plus_1);
    SymbolicExpression ratio_expr = (next_term / summand).simplify();

    Scalar r = 0.0L;
    if (ratio_expr.is_number(&r)) {
        return false;
    }

    constexpr int kNumPoints = 6;
    Scalar ratios_arr[kNumPoints];
    Scalar n_arr[kNumPoints];
    for (int i = 0; i < kNumPoints; ++i) {
        int n = 10 + i;
        SymbolicExpression n_val = SymbolicExpression::number(Scalar(n));
        SymbolicExpression term_n = summand.substitute(index_name, n_val).simplify();
        SymbolicExpression n_plus_val = SymbolicExpression::number(Scalar(n + 1));
        SymbolicExpression term_n1 = summand.substitute(index_name, n_plus_val).simplify();

        Scalar val_n = 0.0L, val_n1 = 0.0L;
        if (!term_n.is_number(&val_n) || !term_n1.is_number(&val_n1)) {
            return false;
        }
        if (mymath::is_near_zero(val_n, 1e-12)) {
            return false;
        }
        ratios_arr[i] = val_n1 / val_n;
        n_arr[i] = Scalar(n);
    }

    bool is_decreasing = true;
    for (int i = 1; i < kNumPoints; ++i) {
        if (ratios_arr[i] > ratios_arr[i-1] + 1e-6) {
            is_decreasing = false;
            break;
        }
    }

    Scalar ratio_diff1 = mymath::abs(ratios_arr[1] - ratios_arr[0]);
    Scalar ratio_diff2 = mymath::abs(ratios_arr[2] - ratios_arr[1]);
    bool is_converging = (ratio_diff2 < ratio_diff1 * 1.5) ||
                         (ratio_diff1 < 1e-3 && ratio_diff2 < 1e-3);

    if (!is_decreasing && !is_converging) {
        return false;
    }

    std::vector<std::vector<Scalar>> T(kNumPoints, std::vector<Scalar>(kNumPoints, 0.0L));

    for (int i = 0; i < kNumPoints; ++i) {
        T[i][0] = ratios_arr[i];
    }

    for (int j = 1; j < kNumPoints; ++j) {
        for (int i = 0; i < kNumPoints - j; ++i) {
            Scalar n_ipj = n_arr[i + j];
            Scalar n_i = n_arr[i];
            T[i][j] = T[i + 1][j - 1] + (T[i + 1][j - 1] - T[i][j - 1]) * n_ipj / (n_ipj - n_i);
        }
    }

    *ratio = T[0][kNumPoints - 1];

    poly_coeffs->clear();

    constexpr int kMaxPolyDegree = 10;
    for (int n = 0; n <= kMaxPolyDegree + 2; ++n) {
        SymbolicExpression n_val = SymbolicExpression::number(Scalar(n));
        SymbolicExpression term_n = summand.substitute(index_name, n_val).simplify();
        Scalar val = 0.0L;
        if (!term_n.is_number(&val)) {
            return false;
        }
        Scalar p_n = val / mymath::pow(*ratio, n);
        poly_coeffs->push_back(p_n);
    }

    std::vector<Scalar> diff = *poly_coeffs;
    for (int order = 1; order <= kMaxPolyDegree + 1; ++order) {
        std::vector<Scalar> new_diff;
        for (std::size_t i = 1; i < diff.size(); ++i) {
            new_diff.push_back(diff[i] - diff[i - 1]);
        }
        diff = new_diff;

        bool all_zero = true;
        for (Scalar d : diff) {
            if (!mymath::is_near_zero(d, 1e-8)) {
                all_zero = false;
                break;
            }
        }
        if (all_zero) {
            int degree = order - 1;
            if (static_cast<int>(poly_coeffs->size()) < degree + 1) {
                return false;
            }

            std::vector<Scalar> coeffs(degree + 1, 0.0L);

            if (degree == 0) {
                coeffs[0] = (*poly_coeffs)[0];
            } else if (degree == 1) {
                coeffs[0] = (*poly_coeffs)[0];
                coeffs[1] = (*poly_coeffs)[1] - (*poly_coeffs)[0];
            } else if (degree == 2) {
                coeffs[0] = (*poly_coeffs)[0];
                Scalar p1 = (*poly_coeffs)[1] - coeffs[0];
                Scalar p2 = (*poly_coeffs)[2] - coeffs[0];
                coeffs[2] = (p2 - 2 * p1) / 2.0;
                coeffs[1] = p1 - coeffs[2];
            } else {
                int n = degree + 1;
                std::vector<std::vector<Scalar>> A(n, std::vector<Scalar>(n, 0.0L));
                std::vector<Scalar> b(n, 0.0L);

                for (int i = 0; i < n; ++i) {
                    b[i] = (*poly_coeffs)[i];
                    Scalar power = 1.0L;
                    for (int j = 0; j < n; ++j) {
                        A[i][j] = power;
                        power *= Scalar(i);
                    }
                }

                for (int i = 0; i < n; ++i) {
                    int max_row = i;
                    for (int k = i + 1; k < n; ++k) {
                        if (mymath::abs(A[k][i]) > mymath::abs(A[max_row][i])) {
                            max_row = k;
                        }
                    }
                    std::swap(A[i], A[max_row]);
                    std::swap(b[i], b[max_row]);

                    if (mymath::abs(A[i][i]) < 1e-12) {
                        return false;
                    }

                    for (int k = i + 1; k < n; ++k) {
                        Scalar factor = A[k][i] / A[i][i];
                        for (int j = i; j < n; ++j) {
                            A[k][j] -= factor * A[i][j];
                        }
                        b[k] -= factor * b[i];
                    }
                }

                for (int i = n - 1; i >= 0; --i) {
                    coeffs[i] = b[i];
                    for (int j = i + 1; j < n; ++j) {
                        coeffs[i] -= A[i][j] * coeffs[j];
                    }
                    coeffs[i] /= A[i][i];
                }
            }

            *poly_coeffs = coeffs;
            return true;
        }

        if (diff.empty()) break;
    }

    return false;
}

std::string compute_power_times_geo_sum(int m, const std::string& r,
                                         const std::string& index_name) {
    if (m == 0) {
        std::ostringstream s0;
        s0 << "(1 - (" << r << ") ^ (" << index_name << " + 1)) / (1 - (" << r << "))";
        return s0.str();
    }

    if (m == 1) {
        std::ostringstream s1;
        s1 << "(" << r << ") * (1 - (" << index_name << " + 1) * ("
           << r << ") ^ " << index_name << " + " << index_name << " * ("
           << r << ") ^ (" << index_name << " + 1)) / (1 - (" << r << ")) ^ 2";
        return s1.str();
    }

    if (m == 2) {
        std::ostringstream s2;
        s2 << "(" << r << ") * (1 + (" << r << ") - ("
           << index_name << " + 1) ^ 2 * (" << r << ") ^ " << index_name
           << " + (2 * (" << index_name << ") ^ 2 + 2 * (" << index_name
           << ") - 1) * (" << r << ") ^ (" << index_name << " + 1) - ("
           << index_name << ") ^ 2 * (" << r << ") ^ (" << index_name
           << " + 2)) / (1 - (" << r << ")) ^ 3";
        return s2.str();
    }

    std::ostringstream sm;
    sm << "((" << index_name << ") ^ " << m << " * (" << r << ") ^ ("
       << index_name << " + 1)";

    for (int j = 0; j < m; ++j) {
        Scalar c = prob::nCr(m, j);
        if (mymath::is_near_zero(c, 1e-10)) continue;

        std::string s_j = compute_power_times_geo_sum(j, r, index_name);
        sm << " + (" << r << ") * (" << format_symbolic_scalar(c) << ") * (" << s_j << ")";
    }

    sm << ") / (1 - (" << r << "))";
    return sm.str();
}

std::string arith_geo_sum(const std::vector<Scalar>& poly_coeffs,
                           Scalar ratio,
                           const std::string& index_name,
                           const SymbolicExpression& upper_expr,
                           const std::string& lower) {
    std::string r = format_symbolic_scalar(ratio);

    if (mymath::is_near_zero(ratio - 1.0L, 1e-10)) {
        std::ostringstream result;
        for (std::size_t m = 0; m < poly_coeffs.size(); ++m) {
            if (mymath::is_near_zero(poly_coeffs[m], 1e-10)) continue;

            std::ostringstream faulhaber;
            faulhaber << "(1 / " << (m + 1) << ") * (";
            bool first = true;
            for (std::size_t j = 0; j <= m; ++j) {
                Scalar bj = prob::bernoulli(static_cast<int>(j));
                if (mymath::is_near_zero(bj, 1e-10)) continue;

                Scalar term_coeff = prob::nCr(m + 1, j) * bj;
                if (!first) faulhaber << " + ";
                faulhaber << format_symbolic_scalar(term_coeff) << " * ("
                          << index_name << " ^ " << (m + 1 - j) << ")";
                first = false;
            }
            faulhaber << ")";

            if (!mymath::is_near_zero(poly_coeffs[m], 1e-10)) {
                if (result.tellp() > 0) result << " + ";
                result << "(" << format_symbolic_scalar(poly_coeffs[m]) << ") * " << faulhaber.str();
            }
        }

        SymbolicExpression primitive = SymbolicExpression::parse(result.str()).simplify();
        SymbolicExpression lower_minus_one = SymbolicExpression::parse("(" + lower + ") - 1").simplify();
        return SymbolicExpression::parse(
                   "(" + primitive.substitute(index_name, upper_expr).to_string() +
                   ") - (" + primitive.substitute(index_name, lower_minus_one).to_string() + ")")
            .simplify()
            .to_string();
    }

    std::ostringstream result;
    for (std::size_t m = 0; m < poly_coeffs.size(); ++m) {
        if (mymath::is_near_zero(poly_coeffs[m], 1e-10)) continue;

        std::string s_m = compute_power_times_geo_sum(static_cast<int>(m), r, index_name);
        if (result.tellp() > 0) result << " + ";
        result << "(" << format_symbolic_scalar(poly_coeffs[m]) << ") * (" << s_m << ")";
    }

    if (result.tellp() == 0) {
        return "0";
    }

    SymbolicExpression primitive = SymbolicExpression::parse(result.str()).simplify();
    SymbolicExpression lower_minus_one = SymbolicExpression::parse("(" + lower + ") - 1").simplify();
    return SymbolicExpression::parse(
               "(" + primitive.substitute(index_name, upper_expr).to_string() +
               ") - (" + primitive.substitute(index_name, lower_minus_one).to_string() + ")")
        .simplify()
        .to_string();
}

bool detect_telescoping(const SymbolicExpression& summand,
                        const std::string& index_name,
                        SymbolicExpression* g_function) {
    if (summand.node_->type == NodeType::kSubtract) {
        SymbolicExpression left(summand.node_->left);
        SymbolicExpression right(summand.node_->right);

        SymbolicExpression n_plus_1 = SymbolicExpression::parse("(" + index_name + ") + 1");
        SymbolicExpression left_shifted = left.substitute(index_name, n_plus_1).simplify();

        if (left_shifted.to_string() == right.to_string() ||
            left_shifted.simplify().to_string() == right.simplify().to_string()) {
            *g_function = left;
            return true;
        }
    }

    if (summand.node_->type == NodeType::kDivide) {
        SymbolicExpression numerator(summand.node_->left);
        SymbolicExpression denominator(summand.node_->right);

        Scalar num_val = 0.0L;
        if (numerator.is_number(&num_val) && mymath::is_near_zero(num_val - 1.0L, 1e-10)) {
            if (denominator.node_->type == NodeType::kMultiply) {
                SymbolicExpression left_factor(denominator.node_->left);
                SymbolicExpression right_factor(denominator.node_->right);

                SymbolicExpression n_plus_1 = SymbolicExpression::parse("(" + index_name + ") + 1");
                SymbolicExpression left_shifted = left_factor.substitute(index_name, n_plus_1).simplify();

                if (left_shifted.to_string() == right_factor.to_string() ||
                    left_shifted.simplify().to_string() == right_factor.simplify().to_string()) {
                    *g_function = SymbolicExpression::parse("1 / (" + left_factor.to_string() + ")");
                    return true;
                }

                SymbolicExpression right_shifted = right_factor.substitute(index_name, n_plus_1).simplify();
                if (right_shifted.to_string() == left_factor.to_string() ||
                    right_shifted.simplify().to_string() == left_factor.simplify().to_string()) {
                    *g_function = SymbolicExpression::parse("1 / (" + right_factor.to_string() + ")");
                    return true;
                }
            }
        }
    }

    return false;
}

bool wynn_epsilon_acceleration(const std::vector<Scalar>& partial_sums,
                                int max_iterations,
                                Scalar tolerance,
                                Scalar* result) {
    if (partial_sums.size() < 3) {
        if (partial_sums.empty()) {
            *result = 0.0L;
            return true;
        }
        *result = partial_sums.back();
        return true;
    }

    const int n = static_cast<int>(partial_sums.size());

    int actual_max_iter = (max_iterations < 0) ? std::min(n / 2, 50) : max_iterations;
    Scalar actual_tolerance = (tolerance < 0) ? 1e-10 : tolerance;

    std::vector<Scalar> e_prev(n + 1, 0.0L);
    std::vector<Scalar> e_curr(n + 1, 0.0L);

    for (int i = 0; i < n; ++i) {
        e_prev[i + 1] = partial_sums[i];
    }

    Scalar best_estimate = partial_sums.back();
    Scalar best_change = mymath::kDoubleMax;
    int consecutive_no_improvement = 0;

    for (int k = 1; k <= actual_max_iter && k < n; ++k) {
        bool all_valid = true;
        for (int j = k; j <= n - k; ++j) {
            Scalar diff = e_prev[j + 1] - e_prev[j - 1];
            if (mymath::abs(diff) < 1e-30) {
                e_curr[j] = e_prev[j];
            } else {
                Scalar eps_diff = e_curr[j - 1] - e_prev[j];
                if (mymath::abs(eps_diff) < 1e-30) {
                    e_curr[j] = e_prev[j + 1];
                } else {
                    e_curr[j] = e_prev[j + 1] + 1.0L / eps_diff;
                }
            }
            if (!mymath::isfinite(e_curr[j])) {
                all_valid = false;
            }
        }

        if (k % 2 == 0 && all_valid) {
            int mid = (k + 1) / 2 + 1;
            if (mid <= n - k) {
                Scalar estimate = e_curr[mid];
                if (mymath::isfinite(estimate)) {
                    Scalar change = mymath::abs(estimate - best_estimate);
                    if (change < best_change) {
                        best_change = change;
                        best_estimate = estimate;
                        consecutive_no_improvement = 0;
                    } else {
                        consecutive_no_improvement++;
                    }

                    if (change < actual_tolerance) {
                        *result = estimate;
                        return true;
                    }

                    if (consecutive_no_improvement >= 3) {
                        break;
                    }
                }
            }
        }

        std::swap(e_prev, e_curr);
    }

    *result = best_estimate;
    return best_change < actual_tolerance * 10;
}

bool euler_transform(const std::vector<Scalar>& terms,
                     Scalar tolerance,
                     Scalar* result) {
    if (terms.empty()) {
        *result = 0.0L;
        return true;
    }

    const int n = static_cast<int>(terms.size());
    std::vector<Scalar> delta = terms;
    std::vector<Scalar> euler_sums;
    euler_sums.reserve(n);

    Scalar sum = 0.0L;
    Scalar pow2 = 1.0L;

    for (int k = 0; k < n && k < 30; ++k) {
        sum += delta[0] / pow2;
        euler_sums.push_back(sum);

        for (int i = 0; i < n - k - 1; ++i) {
            delta[i] = delta[i + 1] - delta[i];
        }

        pow2 *= 2.0;

        if (k > 2 && euler_sums.size() >= 2) {
            Scalar change = mymath::abs(euler_sums[k] - euler_sums[k - 1]);
            if (change < tolerance) {
                *result = sum;
                return true;
            }
        }
    }

    *result = sum;
    return euler_sums.size() >= 3;
}

bool levin_transform(const std::vector<Scalar>& partial_sums,
                     const std::vector<Scalar>& terms,
                     Scalar tolerance,
                     Scalar* result) {
    const int n = static_cast<int>(partial_sums.size());
    if (n < 4) {
        *result = partial_sums.empty() ? 0.0L : partial_sums.back();
        return true;
    }

    Scalar best_estimate = partial_sums.back();
    Scalar best_change = mymath::kDoubleMax;

    for (int k = 2; k < n - 1; ++k) {
        Scalar beta = 1.0L;

        Scalar num = partial_sums[k];
        Scalar den = 1.0L;

        for (int j = 1; j <= k && j < n; ++j) {
            Scalar coeff = Scalar(j) / Scalar(k + 1);
            Scalar term = mymath::pow(beta * coeff, j);
            if (j < static_cast<int>(terms.size())) {
                num += term * terms[k - j];
            }
            den += term;
        }

        if (mymath::abs(den) > 1e-15) {
            Scalar estimate = num / den;
            if (mymath::isfinite(estimate)) {
                Scalar change = mymath::abs(estimate - best_estimate);
                if (change < best_change) {
                    best_change = change;
                    best_estimate = estimate;
                }
                if (change < tolerance) {
                    *result = estimate;
                    return true;
                }
            }
        }
    }

    *result = best_estimate;
    return best_change < tolerance * 100;
}

bool is_alternating_series(const std::vector<Scalar>& terms) {
    if (terms.size() < 3) return false;

    int sign_changes = 0;
    for (std::size_t i = 1; i < terms.size(); ++i) {
        if ((terms[i] > 0) != (terms[i - 1] > 0)) {
            sign_changes++;
        }
    }

    return sign_changes > static_cast<int>(terms.size()) / 2;
}

bool numerical_infinite_sum(const SeriesContext& ctx,
                            const SymbolicExpression& summand,
                            const std::string& index_name,
                            const std::string& lower,
                            int max_terms,
                            Scalar* result) {
    Scalar lower_val = 0.0L;
    try {
        lower_val = ctx.parse_decimal(lower);
    } catch (...) {
        lower_val = 1.0L;
    }

    int actual_max_terms = (max_terms < 0) ? 2000 : max_terms;
    const Scalar tolerance = 1e-12;

    std::vector<Scalar> partial_sums;
    std::vector<Scalar> terms;
    partial_sums.reserve(actual_max_terms);
    terms.reserve(actual_max_terms);

    Scalar running_sum = 0.0L;
    Scalar prev_term = mymath::kDoubleMax;
    int terms_since_improvement = 0;
    Scalar best_term = mymath::kDoubleMax;

    for (int n = 0; n < actual_max_terms; ++n) {
        Scalar idx = lower_val + n;
        Scalar term = 0.0L;
        try {
            term = ctx.evaluate_at(summand, index_name, idx);
        } catch (...) {
            break;
        }

        if (!mymath::isfinite(term)) {
            break;
        }

        running_sum += term;
        partial_sums.push_back(running_sum);
        terms.push_back(term);

        if (n > 10 && mymath::abs(term) > mymath::abs(prev_term) * 1.1) {
            if (!is_alternating_series(terms)) {
                return false;
            }
        }

        if (mymath::abs(term) < tolerance * mymath::abs(running_sum)) {
            break;
        }

        if (mymath::abs(term) < best_term * 0.99) {
            best_term = mymath::abs(term);
            terms_since_improvement = 0;
        } else {
            terms_since_improvement++;
        }

        if (terms_since_improvement > 100 && n > 200) {
            break;
        }

        prev_term = term;
    }

    if (partial_sums.size() < 5) {
        return false;
    }

    bool is_alternating = is_alternating_series(terms);

    std::vector<std::pair<Scalar, bool>> estimates;

    Scalar wynn_result = 0.0L;
    bool wynn_ok = wynn_epsilon_acceleration(partial_sums, -1, tolerance, &wynn_result);
    if (wynn_ok) {
        estimates.push_back({wynn_result, true});
    }

    if (is_alternating) {
        Scalar euler_result = 0.0L;
        bool euler_ok = euler_transform(terms, tolerance, &euler_result);
        if (euler_ok) {
            estimates.push_back({euler_result, true});
        }
    }

    Scalar levin_result = 0.0L;
    bool levin_ok = levin_transform(partial_sums, terms, tolerance, &levin_result);
    if (levin_ok) {
        estimates.push_back({levin_result, true});
    }

    if (estimates.empty()) {
        *result = partial_sums.back();
        return false;
    }

    if (estimates.size() == 1) {
        *result = estimates[0].first;
        return true;
    }

    Scalar best_result = estimates[0].first;
    Scalar min_variance = mymath::kDoubleMax;

    for (const auto& est1 : estimates) {
        Scalar variance = 0.0L;
        for (const auto& est2 : estimates) {
            variance += mymath::abs(est1.first - est2.first);
        }
        if (variance < min_variance) {
            min_variance = variance;
            best_result = est1.first;
        }
    }

    *result = best_result;
    return true;
}

bool try_infinite_series_value(const SymbolicExpression& summand,
                                const std::string& index_name,
                                const std::string& lower,
                                std::string* result) {
    if (lower != "1") {
        return false;
    }

    for (int k = 1; k <= 5; ++k) {
        int two_k = 2 * k;
        std::string inv_pow_str = "1 / (" + index_name + " ^ " + std::to_string(two_k) + ")";
        SymbolicExpression pattern = SymbolicExpression::parse(inv_pow_str);
        SymbolicExpression diff = (summand - pattern).simplify();
        Scalar val = 0.0L;
        if (diff.is_number(&val) && mymath::is_near_zero(val, 1e-10)) {
            Scalar B = prob::bernoulli(two_k);
            Scalar coeff = mymath::pow(2.0, two_k - 1) * mymath::abs(B) / prob::factorial(two_k);
            std::ostringstream ss;
            if (!mymath::is_near_zero(coeff - 1.0L, 1e-10)) {
                ss << format_symbolic_scalar(coeff) << " * ";
            }
            ss << "pi ^ " << two_k;
            *result = ss.str();
            return true;
        }
    }

    std::string alt_harm = "(-1) ^ (" + index_name + " + 1) / " + index_name;
    SymbolicExpression alt_pattern = SymbolicExpression::parse(alt_harm);
    SymbolicExpression alt_diff = (summand - alt_pattern).simplify();
    Scalar alt_val = 0.0L;
    if (alt_diff.is_number(&alt_val) && mymath::is_near_zero(alt_val, 1e-10)) {
        *result = "ln(2)";
        return true;
    }

    std::string telescoping1 = "1 / (" + index_name + " * (" + index_name + " + 1))";
    SymbolicExpression teles_pattern1 = SymbolicExpression::parse(telescoping1);
    SymbolicExpression teles_diff1 = (summand - teles_pattern1).simplify();
    Scalar teles_val1 = 0.0L;
    if (teles_diff1.is_number(&teles_val1) && mymath::is_near_zero(teles_val1, 1e-10)) {
        *result = "1";
        return true;
    }

    std::string telescoping2 = "1 / (" + index_name + " * (" + index_name + " + 2))";
    SymbolicExpression teles_pattern2 = SymbolicExpression::parse(telescoping2);
    SymbolicExpression teles_diff2 = (summand - teles_pattern2).simplify();
    Scalar teles_val2 = 0.0L;
    if (teles_diff2.is_number(&teles_val2) && mymath::is_near_zero(teles_val2, 1e-10)) {
        *result = "3 / 4";
        return true;
    }

    return false;
}

std::string series_sum(const SeriesContext& ctx,
                       const std::string& expr,
                       const std::string& index_name,
                       const std::string& lower,
                       const std::string& upper) {
    SymbolicExpression summand = SymbolicExpression::parse(ctx.expand_inline(expr));
    SymbolicExpression upper_expression;
    const bool upper_is_infinite =
        upper == "inf" || upper == "oo" || upper == "infinity";
    if (!upper_is_infinite) {
        upper_expression = SymbolicExpression::parse(ctx.expand_inline(upper));
    }

    auto make_polynomial_sum_primitive =
        [&](const std::vector<Scalar>& coefficients) {
            std::vector<std::string> pieces;
            for (std::size_t p = 0; p < coefficients.size(); ++p) {
                if (mymath::is_near_zero(coefficients[p], 1e-10)) {
                    continue;
                }

                std::ostringstream poly_part;
                poly_part << "(1 / " << (p + 1) << ") * (";
                bool first = true;
                for (std::size_t j = 0; j <= p; ++j) {
                    Scalar bj = prob::bernoulli(static_cast<int>(j));
                    if (mymath::is_near_zero(bj, 1e-10)) continue;

                    Scalar term_coeff = prob::nCr(p + 1, j) * bj;
                    if (!first) poly_part << " + ";
                    poly_part << format_symbolic_scalar(term_coeff) << " * (" << index_name << " ^ " << (p + 1 - j) << ")";
                    first = false;
                }
                poly_part << ")";

                pieces.push_back("(" + format_symbolic_scalar(coefficients[p]) + ") * " + poly_part.str());
            }

            if (pieces.empty()) {
                return SymbolicExpression::number(0.0L);
            }
            std::ostringstream out;
            for (std::size_t i = 0; i < pieces.size(); ++i) {
                if (i != 0) {
                    out << " + ";
                }
                out << pieces[i];
            }
            return SymbolicExpression::parse(out.str()).simplify();
        };

    auto finite_sum_from_primitive =
        [&](const SymbolicExpression& primitive) {
            const SymbolicExpression lower_minus_one =
                SymbolicExpression::parse("(" + lower + ") - 1").simplify();
            return SymbolicExpression::parse(
                       "(" +
                       primitive.substitute(index_name, upper_expression).to_string() +
                       ") - (" +
                       primitive.substitute(index_name, lower_minus_one).to_string() +
                       ")")
                .simplify()
                .to_string();
        };

    std::vector<Scalar> polynomial_coefficients;
    if (summand.polynomial_coefficients(index_name, &polynomial_coefficients)) {
        if (upper_is_infinite) {
            bool all_zero = true;
            for (Scalar coefficient : polynomial_coefficients) {
                if (!mymath::is_near_zero(coefficient, 1e-10)) {
                    all_zero = false;
                    break;
                }
            }
            if (!all_zero) {
                throw std::runtime_error(
                    "series_sum does not support infinite polynomial sums (diverges)");
            }
            return "0";
        }

        const SymbolicExpression primitive =
            make_polynomial_sum_primitive(polynomial_coefficients);
        return finite_sum_from_primitive(primitive);
    }

    if (upper_is_infinite) {
        std::string special_value;
        if (try_infinite_series_value(summand, index_name, lower, &special_value)) {
            return ctx.simplify_symbolic(special_value);
        }
    }

    Scalar geometric_coefficient = 0.0L;
    Scalar geometric_ratio_value = 0.0L;
    bool is_geometric = detect_geometric_ratio_symbolic(
        summand, index_name, &geometric_coefficient, &geometric_ratio_value);

    if (!is_geometric) {
        auto geometric_ratio_numeric = [&](Scalar* coefficient, Scalar* ratio) {
            Scalar s[4];
            int offset = 0;
            while (offset < 10) {
                for (int i = 0; i < 4; ++i) {
                    s[i] = ctx.evaluate_at(summand, index_name, offset + i);
                }
                if (!mymath::is_near_zero(s[0], 1e-10)) {
                    break;
                }
                offset++;
            }
            if (offset == 10 || mymath::is_near_zero(s[0], 1e-10)) {
                return false;
            }

            const Scalar candidate = s[1] / s[0];
            if (!mymath::is_near_zero(s[2] - s[1] * candidate, 1e-8) ||
                !mymath::is_near_zero(s[3] - s[2] * candidate, 1e-8)) {
                return false;
            }

            if (mymath::is_near_zero(candidate, 1e-10)) {
                *ratio = 0.0L;
                *coefficient = s[0];
            } else {
                *ratio = candidate;
                *coefficient = s[0] / mymath::pow(candidate, offset);
            }
            return true;
        };

        is_geometric = geometric_ratio_numeric(&geometric_coefficient, &geometric_ratio_value);
    }

    std::vector<Scalar> arith_geo_poly;
    Scalar arith_geo_ratio = 0.0L;
    bool is_arith_geo = false;
    if (!is_geometric) {
        is_arith_geo = detect_arith_geo_series(summand, index_name, &arith_geo_poly, &arith_geo_ratio);
    }

    if (is_arith_geo && !upper_is_infinite) {
        return arith_geo_sum(arith_geo_poly, arith_geo_ratio, index_name, upper_expression, lower);
    }

    if (is_geometric) {
        const std::string coefficient_text =
            format_symbolic_scalar(geometric_coefficient);
        const std::string ratio_text = format_symbolic_scalar(geometric_ratio_value);

        if (upper_is_infinite) {
            if (mymath::abs(geometric_ratio_value) >= 1.0L - 1e-10) {
                throw std::runtime_error(
                    "series_sum infinite geometric series requires |r| < 1");
            }
            if (mymath::is_near_zero(geometric_ratio_value - 1.0L, 1e-10)) {
                throw std::runtime_error(
                    "series_sum infinite geometric series diverges for r = 1");
            }
            return ctx.simplify_symbolic(
                "(" + coefficient_text + ") * (" + ratio_text + ") ^ (" +
                lower + ") / (1 - (" + ratio_text + "))");
        }

        const std::string geometric_primitive_text =
            mymath::is_near_zero(geometric_ratio_value - 1.0L, 1e-10)
                ? "(" + coefficient_text + ") * (" + index_name + " + 1)"
                : "(" + coefficient_text + ") * (1 - (" + ratio_text +
                      ") ^ (" + index_name + " + 1)) / (1 - (" +
                      ratio_text + "))";
        const SymbolicExpression primitive =
            SymbolicExpression::parse(geometric_primitive_text).simplify();
        return finite_sum_from_primitive(primitive);
    }

    SymbolicExpression g_function;
    if (detect_telescoping(summand, index_name, &g_function)) {
        if (upper_is_infinite) {
            SymbolicExpression lower_expr = SymbolicExpression::parse(lower);
            return g_function.substitute(index_name, lower_expr).simplify().to_string();
        } else {
            SymbolicExpression lower_expr = SymbolicExpression::parse(lower);
            SymbolicExpression upper_plus_one = SymbolicExpression::parse("(" + upper + ") + 1");
            return SymbolicExpression::parse(
                "(" + g_function.substitute(index_name, lower_expr).to_string() +
                ") - (" + g_function.substitute(index_name, upper_plus_one).to_string() + ")")
                .simplify().to_string();
        }
    }

    if (upper_is_infinite) {
        Scalar numerical_result = 0.0L;
        if (numerical_infinite_sum(ctx, summand, index_name, lower, 1000, &numerical_result)) {
            std::ostringstream ss;
            ss << std::setprecision(15) << numerical_result;
            return ss.str() + " (numerical approximation)";
        }
    }

    throw std::runtime_error(
        "series_sum: unsupported series type. Supported: polynomials, geometric series, "
        "arithmetic-geometric series (n*r^n), zeta(2k) values, telescoping series, "
        "and numerical approximation for convergent infinite series");
}

}  // namespace summation
}  // namespace series_ops
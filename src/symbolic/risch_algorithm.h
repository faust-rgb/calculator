#ifndef RISCH_ALGORITHM_H
#define RISCH_ALGORITHM_H

#include "symbolic/symbolic_expression.h"
#include "symbolic/symbolic_polynomial.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <functional>

/**
 * @class RischAlgorithm
 * @brief 完整的 Risch 积分算法决策过程
 */
class RischAlgorithm {
public:
    // 复数根表示，用于 Rothstein-Trager 算法
    struct ComplexRoot {
        SymbolicExpression real_part;      // 实部
        SymbolicExpression imag_part;      // 虚部
        bool is_complex;                   // true 如果有非零虚部
        bool is_conjugate_pair;            // true 如果表示共轭对 a+bi 和 a-bi

        static ComplexRoot real(const SymbolicExpression& r) {
            return {r, SymbolicExpression::number(0.0), false, false};
        }
        static ComplexRoot complex(const SymbolicExpression& re, const SymbolicExpression& im, bool conjugate = true) {
            return {re, im, true, conjugate};
        }
    };

    struct DifferentialExtension {
        enum class Kind { kLogarithmic, kExponential, kAlgebraic, kNone };
        Kind kind;
        SymbolicExpression argument;
        SymbolicExpression derivation;
        std::string t_name;
        int dependency_depth;
        std::set<std::string> dependencies;
    };

    enum class IntegralType {
        kElementary,
        kNonElementary,
        kUnknown
    };

    struct IntegrationResult {
        bool success;
        IntegralType type;
        SymbolicExpression value;
        std::string message;

        static IntegrationResult elementary(const SymbolicExpression& value) {
            return {true, IntegralType::kElementary, value, ""};
        }
        static IntegrationResult non_elementary(const std::string& msg = "") {
            return {false, IntegralType::kNonElementary, SymbolicExpression::number(0.0), msg};
        }
        static IntegrationResult unknown(const std::string& msg = "") {
            return {false, IntegralType::kUnknown, SymbolicExpression::number(0.0), msg};
        }
    };

    static bool integrate(const SymbolicExpression& expression,
                         const std::string& variable_name,
                         SymbolicExpression* result);

    static IntegrationResult integrate_full(const SymbolicExpression& expression,
                                            const std::string& variable_name);

    static bool integrate_rational(const SymbolicPolynomial& numerator,
                                  const SymbolicPolynomial& denominator,
                                  const std::string& variable_name,
                                  SymbolicExpression* result,
                                  const std::vector<DifferentialExtension>& tower = {},
                                  int tower_index = -1,
                                  const std::string& main_var = "",
                                  const SymbolicExpression* t_prime = nullptr,
                                  DifferentialExtension::Kind kind = DifferentialExtension::Kind::kNone);

    static bool hermite_reduction(const SymbolicPolynomial& numerator,
                                 const SymbolicPolynomial& denominator,
                                 SymbolicExpression* rational_part,
                                 SymbolicPolynomial* reduced_numerator,
                                 SymbolicPolynomial* reduced_denominator,
                                 const std::vector<DifferentialExtension>& tower = {},
                                 int tower_index = -1,
                                 const std::string& main_var = "",
                                 const SymbolicExpression* t_prime = nullptr,
                                 DifferentialExtension::Kind kind = DifferentialExtension::Kind::kNone);

    static bool rothstein_trager(const SymbolicPolynomial& numerator,
                                const SymbolicPolynomial& denominator,
                                const std::string& variable_name,
                                SymbolicExpression* log_part,
                                const std::vector<DifferentialExtension>& tower = {},
                                int tower_index = -1,
                                const std::string& main_var = "",
                                const SymbolicExpression* t_prime = nullptr,
                                DifferentialExtension::Kind kind = DifferentialExtension::Kind::kNone);

    static IntegrationResult integrate_in_extension(const SymbolicExpression& expression,
                                      const std::vector<DifferentialExtension>& tower,
                                      int tower_index,
                                      const std::string& x_var);

private:
    static std::vector<ComplexRoot> find_all_roots(
        const std::vector<SymbolicExpression>& coeffs,
        const std::string& var_name);

    static std::vector<SymbolicExpression> find_integer_roots(
        const std::vector<SymbolicExpression>& coeffs,
        const std::string& var_name);

    static std::vector<SymbolicExpression> find_rational_roots(
        const std::vector<SymbolicExpression>& coeffs,
        const std::string& var_name);

    static std::vector<SymbolicExpression> find_numeric_roots_newton(
        const std::vector<double>& coeffs);

    static bool check_algebraic_independence(const SymbolicExpression& arg,
                                            DifferentialExtension::Kind kind,
                                            const std::vector<DifferentialExtension>& current_tower,
                                            const std::string& x_var,
                                            SymbolicExpression* substitution);

    static SymbolicExpression normalize_logarithm(const SymbolicExpression& expr,
                                                  const std::string& x_var);

    static SymbolicExpression normalize_exponential(const SymbolicExpression& expr,
                                                    const std::string& x_var);

    static std::vector<DifferentialExtension> build_differential_tower(
        const SymbolicExpression& expression,
        const std::string& variable_name);

    static void topological_sort_tower(std::vector<DifferentialExtension>& tower);

    static void collect_transcendental_extensions(
        const SymbolicExpression& expr,
        const std::string& x_var,
        std::vector<std::pair<SymbolicExpression, DifferentialExtension::Kind>>& extensions);

    static IntegrationResult solve_rde(const SymbolicExpression& f,
                         const SymbolicExpression& g,
                         const std::string& x_var,
                         const std::vector<DifferentialExtension>& tower,
                         int tower_index);

    static IntegrationResult solve_polynomial_rde(const SymbolicPolynomial& f_poly,
                                     const SymbolicPolynomial& g_poly,
                                     const std::string& x_var,
                                     const std::vector<DifferentialExtension>& tower,
                                     int tower_index);

    static int compute_rde_degree_bound(const SymbolicPolynomial& f,
                                        const SymbolicPolynomial& g,
                                        const std::vector<DifferentialExtension>& tower,
                                        int tower_index);

    static bool solve_parametric_rde(const SymbolicExpression& f,
                                     const std::vector<SymbolicExpression>& g_list,
                                     const std::string& x_var,
                                     const std::vector<DifferentialExtension>& tower,
                                     int tower_index,
                                     SymbolicExpression* y_out,
                                     std::vector<SymbolicExpression>* c_out);

    static bool solve_polynomial_parametric_rde(
        const SymbolicPolynomial& f_poly,
        const std::vector<SymbolicPolynomial>& g_polys,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index,
        SymbolicPolynomial* y_out,
        std::vector<SymbolicExpression>* c_out);

    static bool solve_linear_system(
        std::vector<std::vector<SymbolicExpression>>& matrix,
        std::vector<SymbolicExpression>& rhs,
        std::vector<SymbolicExpression>* solution);

    static bool solve_coefficient_identity_for_rde(
        const SymbolicPolynomial& D,
        const SymbolicPolynomial& F,
        const SymbolicPolynomial& G,
        const std::string& x_var,
        int max_deg,
        std::vector<SymbolicExpression>* unknowns);

    static SymbolicExpression convert_trig_to_exponential(const SymbolicExpression& expr);

    static IntegralType detect_non_elementary_pattern(const SymbolicExpression& expr,
                                                       const std::string& x_var);

    static bool is_non_elementary_integral(const SymbolicExpression& expr,
                                           const std::string& variable_name);

    static SymbolicExpression complex_log_to_real(const SymbolicExpression& expr,
                                                   const std::string& x_var);

    static bool handle_exponential_special_case(
        const SymbolicExpression& f,
        const SymbolicExpression& g,
        const std::string& x_var,
        const std::vector<DifferentialExtension>& tower,
        int tower_index,
        SymbolicExpression* result);

    static bool is_in_base_field(const SymbolicExpression& expr,
                                 const std::vector<DifferentialExtension>& tower,
                                 int tower_index);
};

#endif // RISCH_ALGORITHM_H

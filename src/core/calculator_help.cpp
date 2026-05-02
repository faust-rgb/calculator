#include "core/calculator.h"
#include "core/calculator_internal_types.h"
#include "core/calculator_module.h"
#include "command/command_registry.h"

#include <stdexcept>
#include <string>
#include <sstream>

namespace {

std::string build_fallback_help(const std::string& topic) {
    // 保留一些最基础的全局帮助，作为兜底
    if (topic == "commands") {
        return
        "Core Commands:\n"
        "  help, :help         Show this help message\n"
        "  :exact on|off       Toggle exact fraction mode\n"
        "  :precision n        Set decimal display precision (1..17)\n"
        "  :vars, :funcs       List variables and functions\n"
        "  :clear, :clearfunc  Management of variables and functions\n"
        "  :save, :load, :run  Persistence and script execution\n"
        "  exit, quit          Exit the calculator";
    }
    return "";
}

std::string supplemental_help(const std::string& topic) {
    if (topic == "commands") {
        return "Command aliases:\n"
               "  help, :help         Show help text\n"
               "  :run file.calc      Execute a script file\n"
               "  :load file          Load state from file\n"
               "  :save file          Save state to file";
    }
    if (topic == "functions") {
        return "Function index:\n"
               "  sqrt cbrt root sinh cosh tanh gamma\n"
               "  inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt\n"
               "  gcd lcm mod factor factorial nCr nPr\n"
               "  deg2rad rad2deg celsius fahrenheit kelvin bin oct hex base\n"
               "  sum avg median percentile quartile\n"
               "  def if elif else while for range return break continue print\n"
               "  poly_add poly_sub poly_mul poly_div roots\n"
               "  double_integral triple_integral triple_integral_sph\n"
               "  simplify diff integral gradient jacobian hessian critical\n"
               "  limit extrema ode ode_table ode_system ode_system_table\n"
               "  lp_max lp_min ilp_max ilp_min milp_max milp_min bip_max bip_min\n"
               "  step delta heaviside impulse\n"
               "  fourier ifourier laplace ilaplace ztrans iztrans\n"
               "  dft fft idft ifft conv convolve\n"
               "  taylor pade puiseux series_sum summation\n"
               "  and or xor not shl shr rol ror popcount bitlen ctz clz parity reverse_bits";
    }
    if (topic == "matrix") {
        return "Matrix guide:\n"
               "  vec mat zeros randmat random_matrix eye identity\n"
               "  resize append_row append_col transpose get set\n"
               "  inverse dot outer null least_squares qr_q qr_r lu_l lu_u svd_u svd_s svd_vt pinv kron hadamard\n"
               "  norm trace det rank rref eigvals eigvecs solve cond diag reshape cholesky schur hessenberg";
    }
    if (topic == "symbolic") {
        return "Symbolic guide:\n"
               "  simplify(expr) diff(expr, x) integral(expr, x)\n"
               "  gradient(expr, x, y) jacobian([f; g], x, y) hessian(expr, x, y)\n"
               "  taylor(expr, a, n) pade(expr, m, n) puiseux series_sum summation\n"
               "  step(t), delta(t) laplace(expr, t, s) fourier(expr, t, w) ztrans";
    }
    if (topic == "analysis") {
        return "Analysis guide:\n"
               "  solve(expr, guess) bisect(expr, left, right) secant(expr, x0, x1) fixed_point(expr, x0)\n"
               "  limit(expr, x0) extrema(f, left, right) integral(f, x0)\n"
               "  double_integral triple_integral ode(rhs, x0, y0, x1) ode_system";
    }
    if (topic == "planning") {
        return "Planning guide:\n"
               "  lp_max(c, A, b, lo, hi) lp_min(c, A, b, lo, hi)\n"
               "  ilp_max ilp_min milp_max milp_min bip_max bip_min";
    }
    if (topic == "variables") {
        return "Variable management:\n"
               "  :vars :clear Clear all variables :funcs :clearfunc";
    }
    if (topic == "programmer") {
        return "Programmer functions:\n"
               "  and or xor not shl shr rol ror popcount bitlen ctz clz bin oct hex base";
    }
    return "";
}

}  // namespace

std::string Calculator::help_text() const {
    std::ostringstream out;
    out << "Help topics:\n"
        << "  :help commands      Show command reference\n"
        << "  :help functions     Show supported functions\n"
        << "  :help matrix        Show matrix usage guide\n"
        << "  :help symbolic      Show symbolic algebra and transform help\n"
        << "  :help analysis      Show calculus, root solving, and ODE help\n"
        << "  :help planning      Show linear/integer planning help\n"
        << "  :help examples      Show example inputs\n"
        << "  :help exact         Show exact fraction mode help\n"
        << "  :help variables     Show variable and function usage help\n"
        << "  :help persistence   Show save/load help\n"
        << "  :help programmer    Show bitwise/base-conversion help\n"
        << "\n" << help_topic("commands")
        << "\n\n" << help_topic("matrix");
    return out.str();
}

std::string Calculator::help_topic(const std::string& topic) const {
    std::ostringstream out;
    bool found = false;

    // 从模块中收集对应主题的内容
    const auto it = impl_->help_topic_to_modules.find(topic);
    if (it != impl_->help_topic_to_modules.end()) {
        for (const auto& module : it->second) {
            std::string snippet = module->get_help_snippet(topic);
            if (!snippet.empty()) {
                if (found) out << "\n\n";
                out << "[" << module->name() << "]\n" << snippet;
                found = true;
            }
        }
    }

    if (!found) {
        std::string fallback = build_fallback_help(topic);
        if (!fallback.empty()) return fallback;
        throw std::runtime_error("unknown help topic: " + topic);
    }

    const std::string extra = supplemental_help(topic);
    if (!extra.empty()) {
        out << "\n\n" << extra;
    }

    return out.str();
}

/**
 * @file standard_math_module.cpp
 * @brief 标准数学函数模块实现
 *
 * 所有函数通过统一的 StoredValue 接口注册，
 * 使用 wrap_scalar 辅助函数将标量计算包装为 StoredValue 签名。
 */

#include "standard_math_module.h"
#include "math/mymath.h"
#include "math/functions/integer/integer_helpers.h"
#include "math/functions/conversion/unit_conversions.h"
#include "core/common/calculator_exceptions.h"
#include "core/types/module_types.h"
#include "core/services/string_utils.h"
#include <algorithm>
#include <map>

namespace {

using Scalar = mymath::Scalar;

auto wrap_scalar(std::function<Scalar(const std::vector<Scalar>&)> f,
                 const std::string& name, std::size_t min_args, std::size_t max_args) {
    return [f = std::move(f), name, min_args, max_args]
           (const std::vector<StoredValue>& args) -> StoredValue {
        if (args.size() < min_args || args.size() > max_args)
            throw MathError(name + " expects " + std::to_string(min_args) +
                            (min_args == max_args ? "" : " to " + std::to_string(max_args)) +
                            " argument(s)");
        std::vector<Scalar> sa;
        sa.reserve(args.size());
        for (const auto& a : args) sa.push_back(a.decimal);
        StoredValue sv;
        sv.decimal = f(sa);
        return sv;
    };
}

} // namespace

std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>>
StandardMathModule::get_functions_map() const {
    std::map<std::string, std::function<StoredValue(const std::vector<StoredValue>&)>> funcs;

    // Trigonometric
    funcs["sin"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::sin(a[0]); }, "sin", 1, 1);
    funcs["cos"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::cos(a[0]); }, "cos", 1, 1);
    funcs["tan"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::tan(a[0]); }, "tan", 1, 1);
    funcs["asin"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::asin(a[0]); }, "asin", 1, 1);
    funcs["acos"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::acos(a[0]); }, "acos", 1, 1);
    funcs["atan"] = wrap_scalar([](const std::vector<Scalar>& a) {
        if (a.size() == 1) return mymath::atan(a[0]);
        return mymath::atan2(a[0], a[1]);
    }, "atan", 1, 2);
    funcs["sec"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::sec(a[0]); }, "sec", 1, 1);
    funcs["csc"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::csc(a[0]); }, "csc", 1, 1);
    funcs["cot"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::cot(a[0]); }, "cot", 1, 1);
    funcs["asec"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::asec(a[0]); }, "asec", 1, 1);
    funcs["acsc"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::acsc(a[0]); }, "acsc", 1, 1);
    funcs["acot"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::acot(a[0]); }, "acot", 1, 1);

    // Hyperbolic
    funcs["sinh"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::sinh(a[0]); }, "sinh", 1, 1);
    funcs["cosh"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::cosh(a[0]); }, "cosh", 1, 1);
    funcs["tanh"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::tanh(a[0]); }, "tanh", 1, 1);
    funcs["asinh"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::asinh(a[0]); }, "asinh", 1, 1);
    funcs["acosh"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::acosh(a[0]); }, "acosh", 1, 1);
    funcs["atanh"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::atanh(a[0]); }, "atanh", 1, 1);

    // Log/Exp
    funcs["exp"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::exp(a[0]); }, "exp", 1, 1);
    funcs["exp2"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::pow(Scalar(2.0L), a[0]); }, "exp2", 1, 1);
    funcs["ln"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::ln(a[0]); }, "ln", 1, 1);
    funcs["log"] = wrap_scalar([](const std::vector<Scalar>& a) {
        if (a.size() == 1) return mymath::ln(a[0]);
        if (a[1] <= Scalar(0.0L) || mymath::is_near_zero(a[1] - Scalar(1.0L)))
            throw MathError("log base must be positive and not equal to 1");
        return mymath::ln(a[0]) / mymath::ln(a[1]);
    }, "log", 1, 2);
    funcs["log2"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::ln(a[0]) / mymath::ln(Scalar(2.0L)); }, "log2", 1, 1);
    funcs["log10"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::ln(a[0]) / mymath::ln(Scalar(10.0L)); }, "log10", 1, 1);

    // Roots/Power
    funcs["sqrt"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::sqrt(a[0]); }, "sqrt", 1, 1);
    funcs["cbrt"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::cbrt(a[0]); }, "cbrt", 1, 1);
    funcs["root"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::root(a[0], a[1]); }, "root", 2, 2);
    funcs["pow"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::pow(a[0], a[1]); }, "pow", 2, 2);

    // Basic
    funcs["abs"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::abs(a[0]); }, "abs", 1, 1);
    funcs["sign"] = wrap_scalar([](const std::vector<Scalar>& a) {
        if (a[0] == Scalar(0.0L)) return Scalar(0.0L);
        return a[0] > Scalar(0.0L) ? Scalar(1.0L) : Scalar(-1.0L);
    }, "sign", 1, 1);
    funcs["floor"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::floor(a[0]); }, "floor", 1, 1);
    funcs["ceil"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::ceil(a[0]); }, "ceil", 1, 1);
    funcs["round"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::round(a[0]); }, "round", 1, 1);
    funcs["trunc"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::trunc(a[0]); }, "trunc", 1, 1);
    funcs["min"] = wrap_scalar([](const std::vector<Scalar>& a) {
        Scalar res = a[0];
        for (std::size_t i = 1; i < a.size(); ++i) res = std::min(res, a[i]);
        return res;
    }, "min", 1, 255);
    funcs["max"] = wrap_scalar([](const std::vector<Scalar>& a) {
        Scalar res = a[0];
        for (std::size_t i = 1; i < a.size(); ++i) res = std::max(res, a[i]);
        return res;
    }, "max", 1, 255);
    funcs["clamp"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return std::clamp(a[0], std::min(a[1], a[2]), std::max(a[1], a[2]));
    }, "clamp", 3, 3);

    // Special Functions
    funcs["gamma"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::gamma(a[0]); }, "gamma", 1, 1);
    funcs["beta"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::beta(a[0], a[1]); }, "beta", 2, 2);
    funcs["zeta"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::zeta(a[0]); }, "zeta", 1, 1);
    funcs["erf"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::erf(a[0]); }, "erf", 1, 1);
    funcs["erfc"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::erfc(a[0]); }, "erfc", 1, 1);
    funcs["bessel"] = wrap_scalar([](const std::vector<Scalar>& a) {
        if (!is_integer_double(a[0])) throw MathError("bessel order must be an integer");
        return mymath::bessel_j(static_cast<int>(round_to_long_long(a[0])), a[1]);
    }, "bessel", 2, 2);
    funcs["bessel_j"] = funcs["bessel"];

    // Conversions
    funcs["deg"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::radians_to_degrees(a[0]); }, "deg", 1, 1);
    funcs["rad"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::degrees_to_radians(a[0]); }, "rad", 1, 1);
    funcs["deg2rad"] = funcs["rad"];
    funcs["rad2deg"] = funcs["deg"];
    funcs["sin_deg"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::sin(mymath::degrees_to_radians(a[0])); }, "sin_deg", 1, 1);
    funcs["cos_deg"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::cos(mymath::degrees_to_radians(a[0])); }, "cos_deg", 1, 1);
    funcs["celsius"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::fahrenheit_to_celsius(a[0]); }, "celsius", 1, 1);
    funcs["fahrenheit"] = wrap_scalar([](const std::vector<Scalar>& a) { return mymath::celsius_to_fahrenheit(a[0]); }, "fahrenheit", 1, 1);
    funcs["kelvin"] = wrap_scalar([](const std::vector<Scalar>& a) { return a[0] + Scalar(273.15L); }, "kelvin", 1, 1);
    funcs["c2f"] = funcs["fahrenheit"];
    funcs["f2c"] = funcs["celsius"];

    // Signals
    funcs["step"] = wrap_scalar([](const std::vector<Scalar>& a) { return a[0] < 0.0L ? 0.0L : 1.0L; }, "step", 1, 1);
    funcs["heaviside"] = funcs["step"];
    funcs["delta"] = wrap_scalar([](const std::vector<Scalar>& a) {
        return mymath::is_near_zero(a[0], Scalar(1e-12L)) ? Scalar(1.0L) : Scalar(0.0L);
    }, "delta", 1, 1);
    funcs["impulse"] = funcs["delta"];

    // Other
    funcs["rat"] = wrap_scalar([](const std::vector<Scalar>& a) { return a[0]; }, "rat", 1, 2);

    return funcs;
}

std::vector<std::string> StandardMathModule::get_function_names() const {
    std::vector<std::string> names;
    auto funcs = get_functions_map();
    for (const auto& [name, _] : funcs) names.push_back(name);
    return names;
}

std::string StandardMathModule::get_help_snippet(const std::string& topic) const {
    if (topic == "functions") {
        return "Common functions:\n"
               "  Trigonometric: sin cos tan sec csc cot asin acos atan ...\n"
               "  Exponential:   exp exp2 ln log log2 log10 pow gamma beta zeta erf ...\n"
               "  Roots:         sqrt cbrt root\n"
               "  Numeric:       abs sign floor ceil round trunc min max clamp sum ...\n"
               "  Aggregates:    mean avg median mode percentile quartile var std ...";
    }
    return "";
}

REGISTER_CALCULATOR_MODULE(StandardMathModule)

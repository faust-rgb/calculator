// ============================================================================
// 格式化工具函数实现
// ============================================================================
//
// 提供数值和存储值的格式化函数。
// ============================================================================

#include "format_utils.h"
#include "calculator_internal_types.h"
#include "math/helpers/integer_helpers.h"
#include "math/mymath.h"
#include "precise/rational.h"
#include "types/stored_value.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

// ============================================================================
// 常数和符号识别表格
// ============================================================================

namespace {

struct NamedConstant {
    double value;
    const char* name;
};

static const NamedConstant kNamedConstants[] = {
    {mymath::kPi, "pi"},
    {mymath::kE, "e"},
    {0.6931471805599453, "ln(2)"},
    {2.302585092994046, "ln(10)"},
    {1.618033988749895, "phi"}, // 黄金分割比
};

std::string format_rational_with_constant(const Rational& r, const std::string& name, bool reciprocal) {
    if (reciprocal) {
        // value = r / name => (r.num / r.den) / name
        if (r.numerator == 1 && r.denominator == 1) return "1 / " + name;
        if (r.denominator == 1) return std::to_string(r.numerator) + " / " + name;
        return std::to_string(r.numerator) + " / (" + std::to_string(r.denominator) + " * " + name + ")";
    } else {
        // value = r * name => (r.num / r.den) * name
        if (r.numerator == 1 && r.denominator == 1) return name;
        if (r.denominator == 1) return std::to_string(r.numerator) + " * " + name;
        if (r.numerator == 1) return name + " / " + std::to_string(r.denominator);
        return std::to_string(r.numerator) + " * " + name + " / " + std::to_string(r.denominator);
    }
}

std::string try_format_with_named_constants(double value, [[maybe_unused]] double eps) {
    const double abs_value = mymath::abs(value);
    Rational r;

    for (const auto& const_entry : kNamedConstants) {
        // 尝试 value = r * C
        if (try_make_simple_rational(abs_value / const_entry.value, 20, &r)) {
            return format_rational_with_constant(r, const_entry.name, false);
        }
        // 尝试 value = r / C
        if (try_make_simple_rational(abs_value * const_entry.value, 20, &r)) {
            return format_rational_with_constant(r, const_entry.name, true);
        }
    }
    return "";
}

std::string try_format_as_sqrt(double value, [[maybe_unused]] double eps) {
    const double abs_value = mymath::abs(value);
    const double squared = abs_value * abs_value;
    
    Rational r;
    // 尝试识别平方后是有理数的情况
    if (try_make_simple_rational(squared, 100, &r)) {
        long long n = r.numerator;
        long long d = r.denominator;

        // 完美平方检测
        auto is_perfect_square = [](long long x) {
            if (x < 0) return false;
            long long s = static_cast<long long>(mymath::sqrt(static_cast<double>(x)) + 0.5);
            return s * s == x;
        };

        if (is_perfect_square(n) && is_perfect_square(d)) {
            // 这其实应该在普通有理数识别中被捕获，但以防万一
            return ""; 
        }

        // 格式化为 sqrt(n) / sqrt(d) -> sqrt(n*d) / d
        long long inner = n * d;
        // 提取 inner 中的平方因子
        long long factor = 1;
        for (long long i = 2; i * i <= inner; ++i) {
            while (inner % (i * i) == 0) {
                factor *= i;
                inner /= (i * i);
            }
        }

        std::string res;
        if (inner == 1) {
            res = std::to_string(factor);
        } else {
            res = (factor == 1 ? "" : std::to_string(factor) + " * ") + "sqrt(" + std::to_string(inner) + ")";
        }

        if (d == 1) return res;
        
        // 如果 n=1, d=2, value = sqrt(1/2) = sqrt(2)/2
        // res 此时已经是 "sqrt(2)"，需要除以 d (d=2? 不，应该是 sqrt(d^2) = d)
        // 这里的逻辑：value = sqrt(n/d) = sqrt(n*d)/d
        return res + " / " + std::to_string(d);
    }
    return "";
}

} // namespace

// ============================================================================
// 符号识别入口
// ============================================================================

std::string try_format_symbolic_extended(double value, double eps) {
    if (!mymath::isfinite(value) || mymath::is_near_zero(value, eps)) {
        return "";
    }

    const bool negative = value < 0.0;
    const double abs_value = mymath::abs(value);

    // 1. 尝试命名常数比例 (pi, e, 等)
    std::string named_form = try_format_with_named_constants(abs_value, eps);
    if (!named_form.empty()) {
        return negative ? "-" + named_form : named_form;
    }

    // 2. 尝试平方根形式
    std::string sqrt_form = try_format_as_sqrt(abs_value, eps);
    if (!sqrt_form.empty()) {
        return negative ? "-" + sqrt_form : sqrt_form;
    }

    return "";
}

std::string try_format_as_pi_fraction(double value, double eps) {
    // 保持向前兼容，调用新的扩展识别
    return try_format_symbolic_extended(value, eps);
}


// ============================================================================
// 显示精度
// ============================================================================

namespace {

int& mutable_process_display_precision() {
    static int precision = kDefaultDisplayPrecision;
    return precision;
}

std::string format_dict_key(const std::string& key) {
    return "\"" + key + "\"";
}

} // namespace

int process_display_precision() {
    return mutable_process_display_precision();
}

void set_process_display_precision(int precision) {
    mutable_process_display_precision() =
        std::clamp(precision, kMinDisplayPrecision, kMaxDisplayPrecision);
}

// ============================================================================
// 数值规范化
// ============================================================================

double normalize_display_decimal(double value) {
    if (mymath::is_near_zero(value, kDisplayZeroEps)) {
        return 0.0;
    }
    if (mymath::abs(value) > kDisplayIntegerEps &&
        is_integer_double(value, kDisplayIntegerEps)) {
        return static_cast<double>(round_to_long_long(value));
    }
    return value;
}

// ============================================================================
// 数值格式化
// ============================================================================

std::string format_decimal(double value) {
    return format_decimal(value, process_display_precision());
}

std::string format_decimal(double value, int precision) {
    value = normalize_display_decimal(value);
    precision = std::clamp(precision, kMinDisplayPrecision, kMaxDisplayPrecision);
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    return out.str();
}

bool try_make_simple_rational(double value,
                              int max_denominator,
                              Rational* rational) {
    if (rational == nullptr || !mymath::isfinite(value)) {
        return false;
    }

    long long numerator = 0;
    long long denominator = 1;
    if (!mymath::approximate_fraction(value,
                                      &numerator,
                                      &denominator,
                                      max_denominator,
                                      1e-10)) {
        return false;
    }

    *rational = Rational(numerator, denominator);
    return true;
}

std::string format_symbolic_number(double value) {
    value = mymath::is_near_zero(value, kDisplayZeroEps) ? 0.0 : value;
    if (mymath::abs(value) > kDisplayIntegerEps &&
        is_integer_double(value, kDisplayIntegerEps)) {
        return std::to_string(round_to_long_long(value));
    }

    // 1. 尝试扩展符号识别（常数、根式）
    std::string extended_form = try_format_symbolic_extended(value, 1e-9);
    if (!extended_form.empty()) {
        return extended_form;
    }

    // 2. 尝试普通有理数识别
    Rational rational;
    if (try_make_simple_rational(value, 999, &rational)) {
        return rational.to_string();
    }

    return format_decimal(value);
}

std::string format_symbolic_scalar(double value) {
    return format_symbolic_number(value);
}

// ============================================================================
// 级数格式化辅助
// ============================================================================

std::string signed_center_text(double center) {
    if (mymath::is_near_zero(center, 1e-12)) {
        return "";
    }
    return center > 0.0
               ? " - " + format_symbolic_number(center)
               : " + " + format_symbolic_number(-center);
}

std::string power_term(const std::string& base, int numerator, int denominator) {
    if (numerator == 0) {
        return "";
    }
    if (denominator != 0 && numerator % denominator == 0) {
        numerator /= denominator;
        denominator = 1;
    }
    if (numerator == denominator) {
        return base;
    }
    if (denominator == 1) {
        return base + " ^ " + std::to_string(numerator);
    }
    if (numerator == 1) {
        return base + " ^ (1 / " + std::to_string(denominator) + ")";
    }
    return base + " ^ (" + std::to_string(numerator) + " / " +
           std::to_string(denominator) + ")";
}

std::string format_term(double coefficient, const std::string& factor) {
    const bool has_factor = !factor.empty();
    const double abs_coefficient = mymath::abs(coefficient);
    const bool omit_unit =
        has_factor && mymath::is_near_zero(abs_coefficient - 1.0, 1e-9);

    if (!has_factor) {
        return format_symbolic_number(coefficient);
    }
    const std::string coeff_text = format_symbolic_number(abs_coefficient);
    if (coeff_text == "1") {
        return coefficient < 0.0 ? "-" + factor : factor;
    }
    if (omit_unit) {
        return coefficient < 0.0 ? "-" + factor : factor;
    }
    return coefficient < 0.0 ? "-" + coeff_text + " * " + factor
                             : coeff_text + " * " + factor;
}

// ============================================================================
// 存储值格式化
// ============================================================================

std::string format_stored_value(const StoredValue& value, bool symbolic_constants_mode) {
    if (value.is_string) {
        return "\"" + value.string_value + "\"";
    }
    if (value.is_list) {
        std::ostringstream out;
        out << "[";
        if (value.list_value) {
            for (std::size_t i = 0; i < value.list_value->size(); ++i) {
                if (i != 0) out << ", ";
                out << format_stored_value((*value.list_value)[i], symbolic_constants_mode);
            }
        }
        out << "]";
        return out.str();
    }
    if (value.is_dict) {
        std::ostringstream out;
        out << "{";
        if (value.dict_value) {
            bool first = true;
            for (const auto& [key, dict_value] : *value.dict_value) {
                if (!first) out << ", ";
                first = false;
                out << format_dict_key(key) << ": "
                    << format_stored_value(dict_value, symbolic_constants_mode);
            }
        }
        out << "}";
        return out.str();
    }
    const std::string& symbolic_text = value.get_symbolic_text(symbolic_constants_mode);
    if (!symbolic_text.empty()) {
        return symbolic_text;
    }
    if (value.has_precise_decimal_text && !value.precise_decimal_text.empty()) {
        if (is_integer_double(value.decimal, kDisplayIntegerEps)) {
            return format_decimal(normalize_display_decimal(value.decimal));
        }
        return value.precise_decimal_text;
    }
    if (value.is_matrix) {
        return value.matrix.to_string();
    }
    if (value.is_complex) {
        return "complex(" + format_symbolic_number(value.complex.real) + ", " +
               format_symbolic_number(value.complex.imag) + ")";
    }
    if (value.exact) {
        return value.rational.to_string();
    }
    return symbolic_constants_mode ? format_symbolic_number(value.decimal)
                                   : format_decimal(normalize_display_decimal(value.decimal));
}

std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode) {
    if (value.is_string) {
        return value.string_value;
    }
    return format_stored_value(value, symbolic_constants_mode);
}

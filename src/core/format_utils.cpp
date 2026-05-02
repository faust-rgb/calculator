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
// 显示精度
// ============================================================================

namespace {

int& mutable_process_display_precision() {
    static int precision = kDefaultDisplayPrecision;
    return precision;
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
    (void)symbolic_constants_mode;
    if (value.is_string) {
        return "\"" + value.string_value + "\"";
    }
    if (value.has_symbolic_text && !value.symbolic_text.empty()) {
        return value.symbolic_text;
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
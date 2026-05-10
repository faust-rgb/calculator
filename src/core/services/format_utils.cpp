// ============================================================================
// 格式化工具函数实现
// ============================================================================
//
// 提供数值和存储值的格式化函数。
// ============================================================================

#include "format_utils.h"
#include "app/scalar_type.h"
#include "calculator_internal_types.h"
#include "math/helpers/integer_helpers.h"
#include "math/mymath.h"
#include "precise/rational.h"
#include "types/stored_value.h"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>

// ============================================================================
// 常数和符号识别表格
// ============================================================================

namespace {

using Scalar = mymath::Scalar;

/**
 * @struct NamedConstant
 * @brief 命名常量定义
 *
 * 存储常见数学常量的值和名称，用于符号格式化时识别。
 */
struct NamedConstant {
    long double value;    ///< 常量值
    const char* name;     ///< 常量名称
};

/// 预定义的命名常量表
static const NamedConstant kNamedConstants[] = {
    {mymath::kPi, "pi"},
    {mymath::kE, "e"},
    {0.6931471805599453, "ln(2)"},
    {2.302585092994046, "ln(10)"},
    {1.618033988749895, "phi"}, // 黄金分割比
};

/**
 * @brief 将有理数与命名常量组合格式化
 * @param r 有理数部分
 * @param name 常量名称
 * @param reciprocal 是否为倒数形式（r / name）
 * @return 格式化后的字符串
 */
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

/**
 * @brief 尝试用命名常量格式化数值
 * @param value 输入数值
 * @param eps 匹配误差阈值
 * @return 格式化后的字符串，如 "pi / 4"、"2 * e"；无法匹配时返回空字符串
 */
std::string try_format_with_named_constants(Scalar value, [[maybe_unused]] long double eps) {
    const Scalar abs_value = mymath::abs(value);
    Rational r;

    for (const auto& const_entry : kNamedConstants) {
        // 尝试 value = r * C
        if (try_make_simple_rational(abs_value / Scalar(const_entry.value), 20, &r)) {
            return format_rational_with_constant(r, const_entry.name, false);
        }
        // 尝试 value = r / C
        if (try_make_simple_rational(abs_value * Scalar(const_entry.value), 20, &r)) {
            return format_rational_with_constant(r, const_entry.name, true);
        }
    }
    return "";
}

/**
 * @brief 尝试将数值格式化为平方根形式
 * @param value 输入数值
 * @param eps 匹配误差阈值
 * @return 格式化后的字符串，如 "sqrt(2)"、"sqrt(3) / 2"；无法匹配时返回空字符串
 *
 * 检测形如 sqrt(n/d) 的值，其中 n/d 为简单有理数。
 * 提取平方因子以简化表达式，如 sqrt(8) = 2*sqrt(2)。
 */
std::string try_format_as_sqrt(Scalar value, [[maybe_unused]] long double eps) {
    const Scalar abs_value = mymath::abs(value);
    const Scalar squared = abs_value * abs_value;

    Rational r;
    // 尝试识别平方后是有理数的情况
    if (try_make_simple_rational(squared, 100, &r)) {
        long long n = r.numerator;
        long long d = r.denominator;

        // 完美平方检测
        auto is_perfect_square = [](long long x) {
            if (x < 0) return false;
            long long s = static_cast<long long>(mymath::sqrt(static_cast<long double>(x)) + 0.5);
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
        const long long common = std::gcd(factor, d);
        if (common > 1) {
            factor /= common;
            d /= common;
            if (inner == 1) {
                res = std::to_string(factor);
            } else {
                res = (factor == 1 ? "" : std::to_string(factor) + " * ") + "sqrt(" + std::to_string(inner) + ")";
            }
            if (d == 1) return res;
        }
        return res + " / " + std::to_string(d);
    }
    return "";
}

} // namespace

// ============================================================================
// 符号识别入口
// ============================================================================

/**
 * @brief 扩展符号格式化尝试
 * @param value 输入数值
 * @param eps 匹配误差阈值
 * @return 格式化后的字符串；无法匹配时返回空字符串
 *
 * 依次尝试：
 * 1. 命名常数比例形式（如 pi/4, 2*e）
 * 2. 平方根形式（如 sqrt(2), sqrt(3)/2）
 */
std::string try_format_symbolic_extended(Scalar value, long double eps) {
    if (!mymath::isfinite(value) || mymath::abs(value) < Scalar(eps)) {
        return "";
    }

    const bool negative = value < Scalar(0);
    const Scalar abs_value = mymath::abs(value);

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

/**
 * @brief 尝试将数值格式化为含 pi 的分数形式（向后兼容接口）
 */
std::string try_format_as_pi_fraction(Scalar value, long double eps) {
    return try_format_symbolic_extended(value, eps);
}


// ============================================================================
// 显示精度
// ============================================================================

namespace {

/**
 * @brief 获取可变的显示精度引用（内部使用）
 * @return 静态精度值的引用
 */
int& mutable_process_display_precision() {
    static int precision = kDefaultDisplayPrecision;
    return precision;
}

/**
 * @brief 格式化字典键
 * @param key 键名
 * @return 带引号的键字符串
 */
std::string format_dict_key(const std::string& key) {
    return "\"" + key + "\"";
}

} // namespace

/**
 * @brief 获取当前显示精度
 * @return 当前显示精度值
 */
int process_display_precision() {
    return mutable_process_display_precision();
}

/**
 * @brief 设置显示精度
 * @param precision 新的显示精度值，将被限制在 [kMinDisplayPrecision, kMaxDisplayPrecision] 范围内
 */
void set_process_display_precision(int precision) {
    mutable_process_display_precision() =
        std::clamp(precision, kMinDisplayPrecision, kMaxDisplayPrecision);
}

// ============================================================================
// 数值规范化
// ============================================================================

/**
 * @brief 规范化显示的小数值
 * @param value 原始数值
 * @return 规范化后的数值
 *
 * 处理以下边界情况：
 * - 接近零的值返回精确的 0
 * - 接近整数的值返回精确的整数
 */
Scalar normalize_display_decimal(Scalar value) {
    if (mymath::abs(value) < kDisplayZeroEps()) {
        return Scalar(0.0L);
    }
    
    // 只有在合理范围内才尝试整数转换
    static const Scalar max_ll(static_cast<long long>(9223372036854775807LL));
    static const Scalar min_ll(static_cast<long long>(-9223372036854775807LL - 1LL));
    
    if (mymath::abs(value) > kDisplayIntegerEps() &&
        value < max_ll && value > min_ll &&
        mymath::is_integer(value.to_long_double())) { // 这里使用 to_long_double 是安全的，因为已经在 long long 范围内
        return Scalar(static_cast<long long>(mymath::round(value.to_long_double())));
    }
    return value;
}

// ============================================================================
// 数值格式化
// ============================================================================

/**
 * @brief 格式化小数值（使用全局显示精度）
 * @param value 要格式化的数值
 * @return 格式化后的字符串
 */
std::string format_decimal(Scalar value) {
    return format_decimal(value, process_display_precision());
}

/**
 * @brief 格式化小数值
 * @param value 要格式化的数值
 * @param precision 显示精度
 * @return 格式化后的字符串
 */
std::string format_decimal(Scalar value, int precision) {
    Scalar normalized = normalize_display_decimal(value);
    
    if (normalized.is_zero()) {
        return "0";
    }

    // 对于 PreciseDecimal，如果它是高精度的，直接使用其 to_string
    // 除非它是一个简单的整数
    if (mymath::abs(normalized) < Scalar(1e-12L) || mymath::abs(normalized) > Scalar(1e15L)) {
        // 对于极大或极小的数，使用科学计数法或 PreciseDecimal 的原生格式
        return normalized.to_string();
    }

    // 默认回退到 long double 格式化用于短小数字（保持美观）
    precision = std::clamp(precision, kMinDisplayPrecision, kMaxDisplayPrecision);
    std::ostringstream out;
    out << std::setprecision(precision) << normalized.to_long_double();
    return out.str();
}

/**
 * @brief 尝试将 Scalar 转换为简单的有理数
 */
bool try_make_simple_rational(Scalar value,
                              int max_denominator,
                              Rational* rational) {
    if (rational == nullptr || !mymath::isfinite(value)) {
        return false;
    }

    long long numerator = 0;
    long long denominator = 1;
    if (!mymath::approximate_fraction(value.to_long_double(),
                                      &numerator,
                                      &denominator,
                                      max_denominator,
                                      1e-8)) {
        return false;
    }

    *rational = Rational(numerator, denominator);
    return true;
}

/**
 * @brief 格式化符号数值
 * @param value 要格式化的数值
 * @return 格式化后的字符串，优先使用符号形式（如 pi/4, sqrt(2)）
 */
std::string format_symbolic_number(Scalar value) {
    const Scalar zero_eps = kDisplayZeroEps();
    const Scalar int_eps = kDisplayIntegerEps();

    if (mymath::abs(value) < zero_eps) {
        return "0";
    }
    
    Scalar normalized = normalize_display_decimal(value);
    if (mymath::abs(normalized) > int_eps &&
        mymath::is_integer(normalized.to_long_double())) {
        return std::to_string(static_cast<long long>(mymath::round(normalized.to_long_double())));
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

/**
 * @brief 格式化符号标量
 */
std::string format_symbolic_scalar(Scalar value) {
    return format_symbolic_number(value);
}

// ============================================================================
// 级数格式化辅助
// ============================================================================

/**
 * @brief 生成带符号的中心文本
 */
std::string signed_center_text(Scalar center) {
    if (mymath::abs(center) < Scalar(1e-12)) {
        return "";
    }
    return center > Scalar(0)
               ? " - " + format_symbolic_number(center)
               : " + " + format_symbolic_number(-center);
}

/**
 * @brief 格式化幂次项
 */
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

/**
 * @brief 格式化级数项
 */
std::string format_term(Scalar coefficient, const std::string& factor) {
    const bool has_factor = !factor.empty();
    const Scalar abs_coefficient = mymath::abs(coefficient);
    const bool omit_unit =
        has_factor && mymath::abs(abs_coefficient - Scalar(1)) < Scalar(1e-9);

    if (!has_factor) {
        return format_symbolic_number(coefficient);
    }
    const std::string coeff_text = format_symbolic_number(abs_coefficient);
    if (coeff_text == "1") {
        return coefficient < Scalar(0) ? "-" + factor : factor;
    }
    if (omit_unit) {
        return coefficient < Scalar(0) ? "-" + factor : factor;
    }
    return coefficient < Scalar(0) ? "-" + coeff_text + " * " + factor
                             : coeff_text + " * " + factor;
}

// ============================================================================
// 存储值格式化
// ============================================================================

/**
 * @brief 格式化 StoredValue 用于显示
 * @param value 要格式化的存储值
 * @param symbolic_constants_mode 是否启用符号常量模式
 * @return 格式化后的字符串
 *
 * 根据值的类型选择合适的格式化方式：
 * - 字符串：加引号输出
 * - 列表：[elem1, elem2, ...]
 * - 字典：{key: value, ...}
 * - 矩阵：使用 Matrix::to_string()
 * - 复数：complex(real, imag)
 * - 精确值：使用有理数或符号形式
 * - 其他：使用小数格式
 */
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
        if (is_integer_double(value.decimal, 1e-9)) {
            return format_decimal(normalize_display_decimal(value.decimal));
        }
        // Preserve the original high-precision text
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

/**
 * @brief 格式化 StoredValue 用于 print 命令
 * @param value 要格式化的存储值
 * @param symbolic_constants_mode 是否启用符号常量模式
 * @return 格式化后的字符串
 *
 * 与 format_stored_value 的区别：
 * - 字符串值不加引号，直接输出内容
 * - 其他类型使用 format_stored_value
 */
std::string format_print_value(const StoredValue& value, bool symbolic_constants_mode) {
    if (value.is_string) {
        return value.string_value;
    }
    return format_stored_value(value, symbolic_constants_mode);
}

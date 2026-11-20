// ============================================================================
// 存储值类型 (Modern C++ std::variant 架构与轻量无状态代理)
// ============================================================================

#ifndef CORE_VALUE_STORED_VALUE_H
#define CORE_VALUE_STORED_VALUE_H

#include "math/numeric/rational/rational.h"
#include "types/scalar_type.h"
#include "math/numeric/types/complex.h"
#include "math/numeric/exact/precise_decimal.h"

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <variant>
#include <type_traits>
#include <cstddef>
#include <stdexcept>

namespace matrix {
template<typename T> struct TMatrix;
using Matrix = TMatrix<Scalar>;
}

/**
 * @struct StoredValue
 * @brief 计算器中存储的值结构体（基于 std::variant）
 */
struct StoredValue {
    using ListType = std::vector<StoredValue>;
    using DictType = std::map<std::string, StoredValue>;

    using VariantType = std::variant<
        std::monostate,                         // 0: 空值 (Nil)
        Scalar,                                 // 1: 标量
        Rational,                               // 2: 精确有理数
        mymath::complex<Scalar>,                // 3: 复数
        std::string,                            // 4: 字符串
        std::shared_ptr<matrix::Matrix>,        // 5: 矩阵
        std::shared_ptr<ListType>,              // 6: 列表
        std::shared_ptr<DictType>               // 7: 字典
    >;

    VariantType data = std::monostate{};

    enum class Kind { Nil, Scalar, Rational, Complex, String, Matrix, List, Dict };

    Kind kind() const {
        return static_cast<Kind>(data.index());
    }

    // 模式与延时解析元数据
    bool exact = false;
    Rational rational;
    mutable bool has_symbolic_text = false;
    bool has_precise_decimal_text = false;
    mutable bool symbolic_computed = false;
    mutable std::string source_expression;
    mutable std::string symbolic_text;
    std::string precise_decimal_text;
    std::shared_ptr<PreciseDecimal> precise_decimal_value;

    // 构造函数
    StoredValue() = default;
    StoredValue(Scalar val) : data(val) {}
    StoredValue(Rational val) : data(val), exact(true), rational(val) {}
    StoredValue(mymath::complex<Scalar> val) : data(val) {}
    StoredValue(std::string val) : data(std::move(val)) {}
    StoredValue(matrix::Matrix val);
    StoredValue(std::shared_ptr<matrix::Matrix> val) : data(std::move(val)) {}
    StoredValue(std::shared_ptr<ListType> val) : data(std::move(val)) {}
    StoredValue(std::shared_ptr<DictType> val) : data(std::move(val)) {}

    // 拷贝与移动操作全部默认（无需自引用指针维护）
    StoredValue(const StoredValue&) = default;
    StoredValue(StoredValue&&) noexcept = default;
    StoredValue& operator=(const StoredValue&) = default;
    StoredValue& operator=(StoredValue&&) noexcept = default;

    // 类型查询接口
    bool is_scalar_type() const {
        return std::holds_alternative<Scalar>(data) || std::holds_alternative<Rational>(data) ||
               (exact && std::holds_alternative<std::monostate>(data));
    }
    bool is_matrix_type() const {
        return std::holds_alternative<std::shared_ptr<matrix::Matrix>>(data) &&
               std::get<std::shared_ptr<matrix::Matrix>>(data) != nullptr;
    }
    bool is_complex_type() const { return std::holds_alternative<mymath::complex<Scalar>>(data); }
    bool is_string_type() const { return std::holds_alternative<std::string>(data); }
    bool is_list_type() const {
        return std::holds_alternative<std::shared_ptr<ListType>>(data) &&
               std::get<std::shared_ptr<ListType>>(data) != nullptr;
    }
    bool is_dict_type() const {
        return std::holds_alternative<std::shared_ptr<DictType>>(data) &&
               std::get<std::shared_ptr<DictType>>(data) != nullptr;
    }

    // 现代 API (Phase 1 方案规范)
    bool is_scalar() const { return is_scalar_type(); }
    bool is_rational() const { return std::holds_alternative<Rational>(data); }
    bool is_nil() const { return std::holds_alternative<std::monostate>(data); }
    bool is_matrix() const { return is_matrix_type(); }
    bool is_complex() const { return is_complex_type(); }
    bool is_string() const { return is_string_type(); }
    bool is_list() const { return is_list_type(); }
    bool is_dict() const { return is_dict_type(); }

    Scalar as_scalar() const {
        if (std::holds_alternative<Scalar>(data)) return std::get<Scalar>(data);
        if (std::holds_alternative<Rational>(data)) {
            const Rational& r = std::get<Rational>(data);
            return Scalar(r.numerator) / Scalar(r.denominator);
        }
        if (exact && std::holds_alternative<std::monostate>(data)) {
            return Scalar(rational.numerator) / Scalar(rational.denominator);
        }
        throw std::runtime_error("StoredValue does not contain a scalar");
    }
    const Rational& as_rational() const {
        if (is_rational()) return std::get<Rational>(data);
        if (exact && std::holds_alternative<std::monostate>(data)) return rational;
        throw std::runtime_error("StoredValue does not contain a rational");
    }
    const matrix::Matrix& as_matrix() const {
        const auto& ptr = std::get<std::shared_ptr<matrix::Matrix>>(data);
        if (!ptr) throw std::runtime_error("null matrix access");
        return *ptr;
    }
    const mymath::complex<Scalar>& as_complex() const { return std::get<mymath::complex<Scalar>>(data); }
    const std::string& as_string() const { return get_string_value(); }

    std::string to_latex() const;

    bool is_exact() const { return exact; }
    void set_exact(bool ex) { exact = ex; }

    // 访问器
    Scalar get_decimal() const {
        if (std::holds_alternative<Scalar>(data)) return std::get<Scalar>(data);
        if (std::holds_alternative<Rational>(data)) {
            const Rational& r = std::get<Rational>(data);
            return Scalar(r.numerator) / Scalar(r.denominator);
        }
        if (exact && std::holds_alternative<std::monostate>(data)) {
            return Scalar(rational.numerator) / Scalar(rational.denominator);
        }
        throw std::runtime_error("StoredValue does not contain a decimal value");
    }
    void set_decimal(Scalar val) {
        data = val;
        exact = false;
        rational = Rational();
    }

    void set_rational(Rational val) {
        data = val;
        exact = true;
        rational = val;
    }

    mymath::complex<Scalar> get_complex() const {
        if (std::holds_alternative<mymath::complex<Scalar>>(data)) return std::get<mymath::complex<Scalar>>(data);
        return mymath::complex<Scalar>(get_decimal(), Scalar(0.0L));
    }
    void set_complex(mymath::complex<Scalar> val) { data = val; }

    std::shared_ptr<matrix::Matrix> get_matrix_ptr() const {
        if (std::holds_alternative<std::shared_ptr<matrix::Matrix>>(data)) return std::get<std::shared_ptr<matrix::Matrix>>(data);
        return nullptr;
    }
    void set_matrix_ptr(std::shared_ptr<matrix::Matrix> val) { data = std::move(val); }

    const std::string& get_string_value() const {
        if (std::holds_alternative<std::string>(data)) return std::get<std::string>(data);
        static const std::string empty;
        return empty;
    }
    void set_string_value(std::string val) { data = std::move(val); }

    std::shared_ptr<ListType> get_list_value() const {
        if (std::holds_alternative<std::shared_ptr<ListType>>(data)) return std::get<std::shared_ptr<ListType>>(data);
        return nullptr;
    }
    void set_list_value(std::shared_ptr<ListType> val) { data = std::move(val); }

    std::shared_ptr<DictType> get_dict_value() const {
        if (std::holds_alternative<std::shared_ptr<DictType>>(data)) return std::get<std::shared_ptr<DictType>>(data);
        return nullptr;
    }
    void set_dict_value(std::shared_ptr<DictType> val) { data = std::move(val); }

    // 现代 Visitor 支持
    template <typename Visitor>
    decltype(auto) visit(Visitor&& visitor) const {
        return std::visit(std::forward<Visitor>(visitor), data);
    }

    const std::string& get_symbolic_text(bool need_symbolic) const {
        if (!need_symbolic || is_matrix_type() || is_complex_type() || is_string_type() || is_list_type() || is_dict_type()) {
            static const std::string empty;
            return empty;
        }
        if (!symbolic_computed && !source_expression.empty()) {
            symbolic_computed = true;
            if (symbolic_text.empty()) {
                symbolic_text = source_expression;
            }
            has_symbolic_text = true;
        }
        return symbolic_text;
    }

    void set_source_expression(const std::string& expr) {
        source_expression = expr;
        symbolic_computed = false;
    }
};

namespace precise {
std::string stored_value_precise_decimal_text(const StoredValue& value);
} // namespace precise

using precise::stored_value_precise_decimal_text;

#endif // CORE_VALUE_STORED_VALUE_H

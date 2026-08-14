// ============================================================================
// 存储值类型 (Modern C++ std::variant 架构与轻量无状态代理)
// ============================================================================

#ifndef TYPES_STORED_VALUE_H
#define TYPES_STORED_VALUE_H

#include "rational.h"
#include "types/scalar_type.h"
#include "math/types/complex.h"
#include "math/precise/precise_decimal.h"

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <variant>
#include <type_traits>
#include <cstddef>

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

    // 拷贝与移动操作全部默认（无需自引用指针维护）
    StoredValue(const StoredValue&) = default;
    StoredValue(StoredValue&&) noexcept = default;
    StoredValue& operator=(const StoredValue&) = default;
    StoredValue& operator=(StoredValue&&) noexcept = default;

    // 类型查询接口
    bool is_scalar_type() const { return std::holds_alternative<Scalar>(data) || std::holds_alternative<Rational>(data); }
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

    const Scalar& as_scalar() const {
        if (std::holds_alternative<Scalar>(data)) return std::get<Scalar>(data);
        static const Scalar zero(0.0L);
        return zero;
    }
    const Rational& as_rational() const { return std::get<Rational>(data); }
    const matrix::Matrix& as_matrix() const { return *std::get<std::shared_ptr<matrix::Matrix>>(data); }
    const mymath::complex<Scalar>& as_complex() const { return std::get<mymath::complex<Scalar>>(data); }
    const std::string& as_string() const { return get_string_value(); }

    bool is_exact() const { return exact; }
    void set_exact(bool ex) { exact = ex; }

    // 访问器
    Scalar get_decimal() const {
        if (std::holds_alternative<Scalar>(data)) return std::get<Scalar>(data);
        if (std::holds_alternative<Rational>(data)) return static_cast<Scalar>(rational_to_double(std::get<Rational>(data)));
        return Scalar(0.0L);
    }
    void set_decimal(Scalar val) { data = val; }

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
    std::shared_ptr<DictType> get_dict_value() const {
        if (std::holds_alternative<std::shared_ptr<DictType>>(data)) return std::get<std::shared_ptr<DictType>>(data);
        return nullptr;
    }

    // ========================================================================
    // 无状态轻量代理包装器（消除自引用 owner 指针）
    // ========================================================================

    struct DecimalProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator Scalar() const;
        explicit operator long double() const { return static_cast<long double>(get()); }
        Scalar get() const;
        long double to_long_double() const { return static_cast<long double>(get()); }
        DecimalProxy& operator=(Scalar v);
        DecimalProxy& operator+=(Scalar v);
        DecimalProxy& operator-=(Scalar v);
        DecimalProxy& operator*=(Scalar v);
        DecimalProxy& operator/=(Scalar v);

        template <typename Stream>
        friend Stream& operator<<(Stream& os, const DecimalProxy& dp) {
            os << dp.get();
            return os;
        }
    };

    struct ComplexProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator mymath::complex<Scalar>() const;
        ComplexProxy& operator=(mymath::complex<Scalar> v);
        Scalar real() const;
        Scalar imag() const;
        void real(Scalar r);
        void imag(Scalar i);
        mymath::complex<Scalar> get() const;
    };

    struct StringProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator const std::string&() const;
        StringProxy& operator=(std::string v);
        const std::string& get() const;

        auto begin() const { return get().begin(); }
        auto end() const { return get().end(); }
        bool empty() const { return get().empty(); }
        std::size_t size() const { return get().size(); }

        friend std::string operator+(const std::string& lhs, const StringProxy& rhs) { return lhs + rhs.get(); }
        friend std::string operator+(const char* lhs, const StringProxy& rhs) { return std::string(lhs) + rhs.get(); }
        friend std::string operator+(const StringProxy& lhs, const std::string& rhs) { return lhs.get() + rhs; }
        friend std::string operator+(const StringProxy& lhs, const char* rhs) { return lhs.get() + rhs; }
        friend bool operator==(const StringProxy& lhs, const StringProxy& rhs) { return lhs.get() == rhs.get(); }
        friend bool operator==(const StringProxy& lhs, const std::string& rhs) { return lhs.get() == rhs; }
        friend bool operator==(const StringProxy& lhs, const char* rhs) { return lhs.get() == rhs; }
    };

    struct MatrixPtrProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator std::shared_ptr<matrix::Matrix>() const;
        MatrixPtrProxy& operator=(std::shared_ptr<matrix::Matrix> v);
        matrix::Matrix* operator->() const;
        matrix::Matrix& operator*() const;
        explicit operator bool() const;
    };

    struct IsMatrixProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator bool() const;
        bool operator()() const { return operator bool(); }
        IsMatrixProxy& operator=(bool set_mat);
    };

    struct IsComplexProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator bool() const;
        bool operator()() const { return operator bool(); }
        IsComplexProxy& operator=(bool set_cplx);
    };

    struct IsStringProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator bool() const;
        bool operator()() const { return operator bool(); }
        IsStringProxy& operator=(bool set_str);
    };

    struct IsListProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator bool() const;
        bool operator()() const { return operator bool(); }
        IsListProxy& operator=(bool set_l);
    };

    struct IsDictProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator bool() const;
        bool operator()() const { return operator bool(); }
        IsDictProxy& operator=(bool set_d);
    };

    struct ListValueProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator std::shared_ptr<ListType>() const;
        ListValueProxy& operator=(std::shared_ptr<ListType> v);
        ListType* operator->() const;
        ListType& operator*() const;
        explicit operator bool() const;
    };

    struct DictValueProxy {
        StoredValue* owner();
        const StoredValue* owner() const;

        operator std::shared_ptr<DictType>() const;
        DictValueProxy& operator=(std::shared_ptr<DictType> v);
        DictType* operator->() const;
        DictType& operator*() const;
        explicit operator bool() const;
    };

    // 字段兼容代理实例（零指针开销）
    DecimalProxy decimal;
    ComplexProxy complex;
    StringProxy string_value;
    MatrixPtrProxy matrix_ptr;

    IsMatrixProxy is_matrix;
    IsComplexProxy is_complex;
    IsStringProxy is_string;
    IsListProxy is_list;
    IsDictProxy is_dict;

    ListValueProxy list_value;
    DictValueProxy dict_value;

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

// ============================================================================
// 代理内联实现（基于 offsetof 动态寻址父对象，零指针存储）
// ============================================================================

inline StoredValue* StoredValue::DecimalProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, decimal));
}
inline const StoredValue* StoredValue::DecimalProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, decimal));
}
inline StoredValue::DecimalProxy::operator Scalar() const { return owner()->get_decimal(); }
inline Scalar StoredValue::DecimalProxy::get() const { return owner()->get_decimal(); }
inline StoredValue::DecimalProxy& StoredValue::DecimalProxy::operator=(Scalar v) { owner()->set_decimal(v); return *this; }
inline StoredValue::DecimalProxy& StoredValue::DecimalProxy::operator+=(Scalar v) { owner()->set_decimal(owner()->get_decimal() + v); return *this; }
inline StoredValue::DecimalProxy& StoredValue::DecimalProxy::operator-=(Scalar v) { owner()->set_decimal(owner()->get_decimal() - v); return *this; }
inline StoredValue::DecimalProxy& StoredValue::DecimalProxy::operator*=(Scalar v) { owner()->set_decimal(owner()->get_decimal() * v); return *this; }
inline StoredValue::DecimalProxy& StoredValue::DecimalProxy::operator/=(Scalar v) { owner()->set_decimal(owner()->get_decimal() / v); return *this; }

inline StoredValue* StoredValue::ComplexProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, complex));
}
inline const StoredValue* StoredValue::ComplexProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, complex));
}
inline StoredValue::ComplexProxy::operator mymath::complex<Scalar>() const { return owner()->get_complex(); }
inline StoredValue::ComplexProxy& StoredValue::ComplexProxy::operator=(mymath::complex<Scalar> v) { owner()->set_complex(v); return *this; }
inline Scalar StoredValue::ComplexProxy::real() const { return owner()->get_complex().real(); }
inline Scalar StoredValue::ComplexProxy::imag() const { return owner()->get_complex().imag(); }
inline void StoredValue::ComplexProxy::real(Scalar r) {
    auto c = owner()->get_complex();
    c.real(r);
    owner()->set_complex(c);
}
inline void StoredValue::ComplexProxy::imag(Scalar i) {
    auto c = owner()->get_complex();
    c.imag(i);
    owner()->set_complex(c);
}
inline mymath::complex<Scalar> StoredValue::ComplexProxy::get() const { return owner()->get_complex(); }

inline StoredValue* StoredValue::StringProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, string_value));
}
inline const StoredValue* StoredValue::StringProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, string_value));
}
inline StoredValue::StringProxy::operator const std::string&() const { return owner()->get_string_value(); }
inline StoredValue::StringProxy& StoredValue::StringProxy::operator=(std::string v) { owner()->set_string_value(std::move(v)); return *this; }
inline const std::string& StoredValue::StringProxy::get() const { return owner()->get_string_value(); }

inline StoredValue* StoredValue::MatrixPtrProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, matrix_ptr));
}
inline const StoredValue* StoredValue::MatrixPtrProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, matrix_ptr));
}
inline StoredValue::MatrixPtrProxy::operator std::shared_ptr<matrix::Matrix>() const { return owner()->get_matrix_ptr(); }
inline StoredValue::MatrixPtrProxy& StoredValue::MatrixPtrProxy::operator=(std::shared_ptr<matrix::Matrix> v) { owner()->set_matrix_ptr(std::move(v)); return *this; }
inline matrix::Matrix* StoredValue::MatrixPtrProxy::operator->() const { return owner()->get_matrix_ptr().get(); }
inline matrix::Matrix& StoredValue::MatrixPtrProxy::operator*() const { return *owner()->get_matrix_ptr(); }
inline StoredValue::MatrixPtrProxy::operator bool() const { return owner()->get_matrix_ptr() != nullptr; }

inline StoredValue* StoredValue::IsMatrixProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, is_matrix));
}
inline const StoredValue* StoredValue::IsMatrixProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, is_matrix));
}

inline StoredValue* StoredValue::IsComplexProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, is_complex));
}
inline const StoredValue* StoredValue::IsComplexProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, is_complex));
}
inline StoredValue::IsComplexProxy::operator bool() const { return owner()->is_complex_type(); }
inline StoredValue::IsComplexProxy& StoredValue::IsComplexProxy::operator=(bool set_cplx) {
    if (set_cplx && !owner()->is_complex_type()) {
        owner()->data = mymath::complex<Scalar>(Scalar(0.0L), Scalar(0.0L));
    } else if (!set_cplx && owner()->is_complex_type()) {
        owner()->data = std::monostate{};
    }
    return *this;
}

inline StoredValue* StoredValue::IsStringProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, is_string));
}
inline const StoredValue* StoredValue::IsStringProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, is_string));
}
inline StoredValue::IsStringProxy::operator bool() const { return owner()->is_string_type(); }
inline StoredValue::IsStringProxy& StoredValue::IsStringProxy::operator=(bool set_str) {
    if (set_str && !owner()->is_string_type()) {
        owner()->data = std::string("");
    } else if (!set_str && owner()->is_string_type()) {
        owner()->data = std::monostate{};
    }
    return *this;
}

inline StoredValue* StoredValue::IsListProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, is_list));
}
inline const StoredValue* StoredValue::IsListProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, is_list));
}
inline StoredValue::IsListProxy::operator bool() const { return owner()->is_list_type(); }
inline StoredValue::IsListProxy& StoredValue::IsListProxy::operator=(bool set_l) {
    if (set_l && !owner()->is_list_type()) {
        owner()->data = std::make_shared<ListType>();
    } else if (!set_l && owner()->is_list_type()) {
        owner()->data = std::monostate{};
    }
    return *this;
}

inline StoredValue* StoredValue::IsDictProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, is_dict));
}
inline const StoredValue* StoredValue::IsDictProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, is_dict));
}
inline StoredValue::IsDictProxy::operator bool() const { return owner()->is_dict_type(); }
inline StoredValue::IsDictProxy& StoredValue::IsDictProxy::operator=(bool set_d) {
    if (set_d && !owner()->is_dict_type()) {
        owner()->data = std::make_shared<DictType>();
    } else if (!set_d && owner()->is_dict_type()) {
        owner()->data = std::monostate{};
    }
    return *this;
}

inline StoredValue* StoredValue::ListValueProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, list_value));
}
inline const StoredValue* StoredValue::ListValueProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, list_value));
}
inline StoredValue::ListValueProxy::operator std::shared_ptr<StoredValue::ListType>() const { return owner()->get_list_value(); }
inline StoredValue::ListValueProxy& StoredValue::ListValueProxy::operator=(std::shared_ptr<ListType> v) { owner()->data = std::move(v); return *this; }
inline StoredValue::ListType* StoredValue::ListValueProxy::operator->() const { return owner()->get_list_value().get(); }
inline StoredValue::ListType& StoredValue::ListValueProxy::operator*() const { return *owner()->get_list_value(); }
inline StoredValue::ListValueProxy::operator bool() const { return owner()->get_list_value() != nullptr; }

inline StoredValue* StoredValue::DictValueProxy::owner() {
    return reinterpret_cast<StoredValue*>(reinterpret_cast<char*>(this) - offsetof(StoredValue, dict_value));
}
inline const StoredValue* StoredValue::DictValueProxy::owner() const {
    return reinterpret_cast<const StoredValue*>(reinterpret_cast<const char*>(this) - offsetof(StoredValue, dict_value));
}
inline StoredValue::DictValueProxy::operator std::shared_ptr<StoredValue::DictType>() const { return owner()->get_dict_value(); }
inline StoredValue::DictValueProxy& StoredValue::DictValueProxy::operator=(std::shared_ptr<DictType> v) { owner()->data = std::move(v); return *this; }
inline StoredValue::DictType* StoredValue::DictValueProxy::operator->() const { return owner()->get_dict_value().get(); }
inline StoredValue::DictType& StoredValue::DictValueProxy::operator*() const { return *owner()->get_dict_value(); }
inline StoredValue::DictValueProxy::operator bool() const { return owner()->get_dict_value() != nullptr; }

namespace precise {
std::string stored_value_precise_decimal_text(const StoredValue& value);
} // namespace precise

using precise::stored_value_precise_decimal_text;

#endif // TYPES_STORED_VALUE_H

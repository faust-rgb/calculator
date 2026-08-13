// ============================================================================
// 存储值类型 (Modern C++ std::variant 架构)
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

    void init_proxies() {
        decimal.owner = this;
        complex.owner = this;
        string_value.owner = this;
        matrix_ptr.owner = this;
        is_matrix.owner = this;
        is_complex.owner = this;
        is_string.owner = this;
        is_list.owner = this;
        is_dict.owner = this;
        list_value.owner = this;
        dict_value.owner = this;
    }

    // 构造函数
    StoredValue() { init_proxies(); }
    StoredValue(Scalar val) : data(val) { init_proxies(); }
    StoredValue(Rational val) : data(val), exact(true), rational(val) { init_proxies(); }
    StoredValue(mymath::complex<Scalar> val) : data(val) { init_proxies(); }
    StoredValue(std::string val) : data(std::move(val)) { init_proxies(); }
    StoredValue(matrix::Matrix val);
    StoredValue(std::shared_ptr<matrix::Matrix> val) : data(std::move(val)) { init_proxies(); }

    // 拷贝与移动构造函数
    StoredValue(const StoredValue& other)
        : data(other.data), exact(other.exact), rational(other.rational),
          has_symbolic_text(other.has_symbolic_text), has_precise_decimal_text(other.has_precise_decimal_text),
          symbolic_computed(other.symbolic_computed), source_expression(other.source_expression),
          symbolic_text(other.symbolic_text), precise_decimal_text(other.precise_decimal_text),
          precise_decimal_value(other.precise_decimal_value) {
        init_proxies();
    }

    StoredValue(StoredValue&& other) noexcept
        : data(std::move(other.data)), exact(other.exact), rational(std::move(other.rational)),
          has_symbolic_text(other.has_symbolic_text), has_precise_decimal_text(other.has_precise_decimal_text),
          symbolic_computed(other.symbolic_computed), source_expression(std::move(other.source_expression)),
          symbolic_text(std::move(other.symbolic_text)), precise_decimal_text(std::move(other.precise_decimal_text)),
          precise_decimal_value(std::move(other.precise_decimal_value)) {
        init_proxies();
    }

    StoredValue& operator=(const StoredValue& other) {
        if (this != &other) {
            data = other.data;
            exact = other.exact;
            rational = other.rational;
            has_symbolic_text = other.has_symbolic_text;
            has_precise_decimal_text = other.has_precise_decimal_text;
            symbolic_computed = other.symbolic_computed;
            source_expression = other.source_expression;
            symbolic_text = other.symbolic_text;
            precise_decimal_text = other.precise_decimal_text;
            precise_decimal_value = other.precise_decimal_value;
            init_proxies();
        }
        return *this;
    }

    StoredValue& operator=(StoredValue&& other) noexcept {
        if (this != &other) {
            data = std::move(other.data);
            exact = other.exact;
            rational = std::move(other.rational);
            has_symbolic_text = other.has_symbolic_text;
            has_precise_decimal_text = other.has_precise_decimal_text;
            symbolic_computed = other.symbolic_computed;
            source_expression = std::move(other.source_expression);
            symbolic_text = std::move(other.symbolic_text);
            precise_decimal_text = std::move(other.precise_decimal_text);
            precise_decimal_value = std::move(other.precise_decimal_value);
            init_proxies();
        }
        return *this;
    }

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

    // 代理包装器实现旧接口向后兼容
    struct DecimalProxy {
        StoredValue* owner;
        operator Scalar() const { return owner->get_decimal(); }
        explicit operator long double() const { return static_cast<long double>(owner->get_decimal()); }
        Scalar get() const { return owner->get_decimal(); }
        long double to_long_double() const { return static_cast<long double>(owner->get_decimal()); }
        DecimalProxy& operator=(Scalar v) { owner->set_decimal(v); return *this; }
        DecimalProxy& operator+=(Scalar v) { owner->set_decimal(owner->get_decimal() + v); return *this; }
        DecimalProxy& operator-=(Scalar v) { owner->set_decimal(owner->get_decimal() - v); return *this; }
        DecimalProxy& operator*=(Scalar v) { owner->set_decimal(owner->get_decimal() * v); return *this; }
        DecimalProxy& operator/=(Scalar v) { owner->set_decimal(owner->get_decimal() / v); return *this; }

        template <typename Stream>
        friend Stream& operator<<(Stream& os, const DecimalProxy& dp) {
            os << dp.get();
            return os;
        }
    };

    struct ComplexProxy {
        StoredValue* owner;
        operator mymath::complex<Scalar>() const { return owner->get_complex(); }
        ComplexProxy& operator=(mymath::complex<Scalar> v) { owner->set_complex(v); return *this; }
        Scalar real() const { return owner->get_complex().real(); }
        Scalar imag() const { return owner->get_complex().imag(); }
        void real(Scalar r) {
            auto c = owner->get_complex();
            c.real(r);
            owner->set_complex(c);
        }
        void imag(Scalar i) {
            auto c = owner->get_complex();
            c.imag(i);
            owner->set_complex(c);
        }
        mymath::complex<Scalar> get() const { return owner->get_complex(); }
    };

    struct StringProxy {
        StoredValue* owner;
        operator const std::string&() const { return owner->get_string_value(); }
        StringProxy& operator=(std::string v) { owner->set_string_value(std::move(v)); return *this; }
        const std::string& get() const { return owner->get_string_value(); }

        auto begin() const { return owner->get_string_value().begin(); }
        auto end() const { return owner->get_string_value().end(); }
        bool empty() const { return owner->get_string_value().empty(); }
        std::size_t size() const { return owner->get_string_value().size(); }

        friend std::string operator+(const std::string& lhs, const StringProxy& rhs) { return lhs + rhs.get(); }
        friend std::string operator+(const char* lhs, const StringProxy& rhs) { return std::string(lhs) + rhs.get(); }
        friend std::string operator+(const StringProxy& lhs, const std::string& rhs) { return lhs.get() + rhs; }
        friend std::string operator+(const StringProxy& lhs, const char* rhs) { return lhs.get() + rhs; }
        friend bool operator==(const StringProxy& lhs, const StringProxy& rhs) { return lhs.get() == rhs.get(); }
        friend bool operator==(const StringProxy& lhs, const std::string& rhs) { return lhs.get() == rhs; }
        friend bool operator==(const StringProxy& lhs, const char* rhs) { return lhs.get() == rhs; }
    };

    struct MatrixPtrProxy {
        StoredValue* owner;
        operator std::shared_ptr<matrix::Matrix>() const { return owner->get_matrix_ptr(); }
        MatrixPtrProxy& operator=(std::shared_ptr<matrix::Matrix> v) { owner->set_matrix_ptr(std::move(v)); return *this; }
        matrix::Matrix* operator->() const { return owner->get_matrix_ptr().get(); }
        matrix::Matrix& operator*() const { return *owner->get_matrix_ptr(); }
        explicit operator bool() const { return owner->get_matrix_ptr() != nullptr; }
    };

    struct IsMatrixProxy {
        StoredValue* owner;
        operator bool() const;
        bool operator()() const { return operator bool(); }
        IsMatrixProxy& operator=(bool set_mat);
    };

    struct IsComplexProxy {
        StoredValue* owner;
        operator bool() const { return owner->is_complex_type(); }
        bool operator()() const { return operator bool(); }
        IsComplexProxy& operator=(bool set_cplx) {
            if (set_cplx && !owner->is_complex_type()) {
                owner->data = mymath::complex<Scalar>(Scalar(0.0L), Scalar(0.0L));
            } else if (!set_cplx && owner->is_complex_type()) {
                owner->data = std::monostate{};
            }
            return *this;
        }
    };

    struct IsStringProxy {
        StoredValue* owner;
        operator bool() const { return owner->is_string_type(); }
        bool operator()() const { return operator bool(); }
        IsStringProxy& operator=(bool set_str) {
            if (set_str && !owner->is_string_type()) {
                owner->data = std::string("");
            } else if (!set_str && owner->is_string_type()) {
                owner->data = std::monostate{};
            }
            return *this;
        }
    };

    struct IsListProxy {
        StoredValue* owner;
        operator bool() const { return owner->is_list_type(); }
        bool operator()() const { return operator bool(); }
        IsListProxy& operator=(bool set_l) {
            if (set_l && !owner->is_list_type()) {
                owner->data = std::make_shared<ListType>();
            } else if (!set_l && owner->is_list_type()) {
                owner->data = std::monostate{};
            }
            return *this;
        }
    };

    struct IsDictProxy {
        StoredValue* owner;
        operator bool() const { return owner->is_dict_type(); }
        bool operator()() const { return operator bool(); }
        IsDictProxy& operator=(bool set_d) {
            if (set_d && !owner->is_dict_type()) {
                owner->data = std::make_shared<DictType>();
            } else if (!set_d && owner->is_dict_type()) {
                owner->data = std::monostate{};
            }
            return *this;
        }
    };

    struct ListValueProxy {
        StoredValue* owner;
        operator std::shared_ptr<ListType>() const { return owner->get_list_value(); }
        ListValueProxy& operator=(std::shared_ptr<ListType> v) { owner->data = std::move(v); return *this; }
        ListType* operator->() const { return owner->get_list_value().get(); }
        ListType& operator*() const { return *owner->get_list_value(); }
        explicit operator bool() const { return owner->get_list_value() != nullptr; }
    };

    struct DictValueProxy {
        StoredValue* owner;
        operator std::shared_ptr<DictType>() const { return owner->get_dict_value(); }
        DictValueProxy& operator=(std::shared_ptr<DictType> v) { owner->data = std::move(v); return *this; }
        DictType* operator->() const { return owner->get_dict_value().get(); }
        DictType& operator*() const { return *owner->get_dict_value(); }
        explicit operator bool() const { return owner->get_dict_value() != nullptr; }
    };

    // 字段兼容代理
    DecimalProxy decimal{this};
    ComplexProxy complex{this};
    StringProxy string_value{this};
    MatrixPtrProxy matrix_ptr{this};

    IsMatrixProxy is_matrix{this};
    IsComplexProxy is_complex{this};
    IsStringProxy is_string{this};
    IsListProxy is_list{this};
    IsDictProxy is_dict{this};

    ListValueProxy list_value{this};
    DictValueProxy dict_value{this};

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

#endif // TYPES_STORED_VALUE_H

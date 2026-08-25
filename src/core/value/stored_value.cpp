// ============================================================================
// StoredValue 代理实现 (src/core/value/stored_value.cpp)
// ============================================================================

#include "core/value/stored_value.h"
#include "matrix/matrix.h"

StoredValue::StoredValue(matrix::Matrix val)
    : data(std::make_shared<matrix::Matrix>(std::move(val))) {
}

StoredValue::IsMatrixProxy::operator bool() const {
    return owner()->is_matrix_type();
}

StoredValue::IsMatrixProxy& StoredValue::IsMatrixProxy::operator=(bool set_mat) {
    if (set_mat && !owner()->is_matrix_type()) {
        owner()->data = std::make_shared<matrix::Matrix>();
    } else if (!set_mat && owner()->is_matrix_type()) {
        owner()->data = std::monostate{};
    }
    return *this;
}

std::string StoredValue::to_latex() const {
    if (is_matrix_type()) {
        auto m = get_matrix_ptr();
        return m ? m->to_latex() : "\\begin{pmatrix}\\end{pmatrix}";
    }
    if (is_complex_type()) {
        auto c = get_complex();
        Scalar re = c.real();
        Scalar im = c.imag();
        if (im == Scalar(0)) return re.to_string();
        if (re == Scalar(0)) {
            if (im == Scalar(1)) return "i";
            if (im == Scalar(-1)) return "-i";
            return im.to_string() + " i";
        }
        if (im > Scalar(0)) {
            if (im == Scalar(1)) return re.to_string() + " + i";
            return re.to_string() + " + " + im.to_string() + " i";
        } else {
            if (im == Scalar(-1)) return re.to_string() + " - i";
            return re.to_string() + " - " + (-im).to_string() + " i";
        }
    }
    if (is_rational()) {
        const auto& r = as_rational();
        if (r.denominator == 1) return std::to_string(r.numerator);
        if (r.numerator < 0) {
            return "-\\frac{" + std::to_string(-r.numerator) + "}{" + std::to_string(r.denominator) + "}";
        }
        return "\\frac{" + std::to_string(r.numerator) + "}{" + std::to_string(r.denominator) + "}";
    }
    if (exact && rational.denominator != 0) {
        if (rational.denominator == 1) return std::to_string(rational.numerator);
        if (rational.numerator < 0) {
            return "-\\frac{" + std::to_string(-rational.numerator) + "}{" + std::to_string(rational.denominator) + "}";
        }
        return "\\frac{" + std::to_string(rational.numerator) + "}{" + std::to_string(rational.denominator) + "}";
    }
    if (is_string_type()) {
        return "\\text{" + get_string_value() + "}";
    }
    if (is_list_type()) {
        auto l = get_list_value();
        if (!l) return "\\left[\\right]";
        std::string res = "\\left[ ";
        for (size_t i = 0; i < l->size(); ++i) {
            if (i > 0) res += ", ";
            res += (*l)[i].to_latex();
        }
        res += " \\right]";
        return res;
    }
    if (is_dict_type()) {
        auto d = get_dict_value();
        if (!d) return "\\left\\{\\right\\}";
        std::string res = "\\left\\{ ";
        bool first = true;
        for (const auto& [k, v] : *d) {
            if (!first) res += ", ";
            first = false;
            res += "\\text{" + k + "}: " + v.to_latex();
        }
        res += " \\right\\}";
        return res;
    }
    if (is_scalar()) {
        return get_decimal().to_string();
    }
    return "";
}

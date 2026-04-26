#include "calculator_internal_types.h"

#include "line_executor.h"

#include <cctype>
#include <stdexcept>

namespace {

bool has_identifier_char(char ch) {
    return std::isalnum(static_cast<unsigned char>(ch)) || ch == '_';
}

bool contains_v2_imaginary_unit(const std::string& text) {
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] != 'i') {
            continue;
        }
        const bool prev_identifier =
            i > 0 && (std::isalpha(static_cast<unsigned char>(text[i - 1])) ||
                      text[i - 1] == '_');
        const bool next_identifier =
            i + 1 < text.size() && has_identifier_char(text[i + 1]);
        if (!prev_identifier && !next_identifier) {
            return true;
        }
    }
    return false;
}

bool parse_positive_int(const std::string& text, int* value) {
    const std::string trimmed = trim_copy(text);
    if (trimmed.empty()) {
        return false;
    }
    int parsed = 0;
    for (char ch : trimmed) {
        if (!std::isdigit(static_cast<unsigned char>(ch))) {
            return false;
        }
        parsed = parsed * 10 + (ch - '0');
        if (parsed > 10000) {
            return false;
        }
    }
    if (parsed <= 0) {
        return false;
    }
    *value = parsed;
    return true;
}

}  // namespace

bool set_v2_precision(Calculator::Impl* impl,
                      const std::string& digits_text,
                      std::string* output) {
    int digits = 0;
    if (!parse_positive_int(digits_text, &digits)) {
        throw std::runtime_error(":precision expects a positive integer");
    }
    impl->v2_environment.precision().digits = digits;
    *output = "2.0 precision: " + std::to_string(digits);
    return true;
}

bool try_process_v2_line(Calculator::Impl* impl,
                         const std::string& expression,
                         std::string* output) {
    const std::string trimmed = trim_copy(expression);
    if (trimmed.empty()) {
        return false;
    }

    if (trimmed == ":v2vars") {
        *output = runtime::list_variables(impl->v2_environment);
        return true;
    }
    if (trimmed == ":v2clear") {
        impl->v2_environment.clear();
        *output = "Cleared all 2.0 variables.";
        return true;
    }
    if (trimmed.rfind(":v2precision ", 0) == 0) {
        return set_v2_precision(impl, trimmed.substr(13), output);
    }
    if (trimmed.rfind(":precision ", 0) == 0) {
        return set_v2_precision(impl, trimmed.substr(11), output);
    }
    if (trimmed.rfind(":v2 ", 0) == 0) {
        *output = runtime::evaluate_line(trimmed.substr(4), impl->v2_environment).to_string();
        return true;
    }

    if (!contains_v2_imaginary_unit(trimmed)) {
        return false;
    }

    *output = runtime::evaluate_line(trimmed, impl->v2_environment).to_string();
    return true;
}

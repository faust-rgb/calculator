#include "calculator_internal_types.h"

#include "symbolic_expression.h"

#include "mymath.h"

#include <algorithm>
#include <stdexcept>

void apply_calculator_display_precision(const Calculator::Impl* impl) {
    const int precision = impl == nullptr ? kDefaultDisplayPrecision : impl->display_precision;
    set_process_display_precision(precision);
    matrix::set_display_precision(precision);
    SymbolicExpression::set_display_precision(precision);
}

Calculator::Calculator() : impl_(new Impl()) {
    apply_calculator_display_precision(impl_.get());
}

Calculator::~Calculator() = default;

std::string Calculator::clear_variable(const std::string& name) {
    const auto it = impl_->variables.find(name);
    if (it == impl_->variables.end()) {
        throw std::runtime_error("unknown variable: " + name);
    }
    impl_->variables.erase(it);
    return "Cleared variable: " + name;
}

std::string Calculator::clear_all_variables() {
    impl_->variables.clear();
    return "Cleared all variables.";
}

std::string Calculator::set_hex_prefix_mode(bool enabled) {
    impl_->hex_prefix_mode = enabled;
    return std::string("Hex prefix mode: ") + (enabled ? "ON" : "OFF");
}

bool Calculator::hex_prefix_mode() const {
    return impl_->hex_prefix_mode;
}

std::string Calculator::set_hex_uppercase_mode(bool enabled) {
    impl_->hex_uppercase_mode = enabled;
    return std::string("Hex letter case: ") + (enabled ? "UPPER" : "LOWER");
}

bool Calculator::hex_uppercase_mode() const {
    return impl_->hex_uppercase_mode;
}

std::string Calculator::set_symbolic_constants_mode(bool enabled) {
    impl_->symbolic_constants_mode = enabled;
    return std::string("Symbolic constants mode: ") + (enabled ? "ON" : "OFF");
}

bool Calculator::symbolic_constants_mode() const {
    return impl_->symbolic_constants_mode;
}

std::string Calculator::set_display_precision(int precision) {
    if (precision < kMinDisplayPrecision || precision > kMaxDisplayPrecision) {
        throw std::runtime_error("display precision must be in the range 1..17");
    }
    impl_->display_precision = precision;
    apply_calculator_display_precision(impl_.get());
    return "Display precision: " + std::to_string(precision);
}

int Calculator::display_precision() const {
    return impl_->display_precision;
}

std::vector<std::string> Calculator::variable_names() const {
    std::vector<std::string> names;
    names.reserve(impl_->variables.size());
    for (const auto& [name, _] : impl_->variables) {
        names.push_back(name);
    }
    return names;
}

std::vector<std::string> Calculator::custom_function_names() const {
    std::vector<std::string> names;
    names.reserve(impl_->functions.size());
    for (const auto& [name, _] : impl_->functions) {
        names.push_back(name);
    }
    return names;
}

double Calculator::normalize_result(double value) {
    if (!mymath::isfinite(value)) {
        return value;
    }
    if (mymath::abs(value) < kDisplayZeroEps) {
        return 0.0;
    }
    if (mymath::abs(value) > kDisplayIntegerEps &&
        is_integer_double(value, kDisplayIntegerEps)) {
        return static_cast<double>(round_to_long_long(value));
    }
    return value;
}

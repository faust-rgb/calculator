#include "calculator_internal_types.h"

#include "mymath.h"

#include <algorithm>
#include <stdexcept>

namespace {

constexpr double kDisplayZeroEps = mymath::kDoubleDenormMin;
constexpr double kDisplayIntegerEps = 1e-9;

bool is_integer_double(double x, double eps = 1e-10) {
    return mymath::is_integer(x, eps);
}

long long round_to_long_long(double x) {
    return static_cast<long long>(x >= 0.0 ? x + 0.5 : x - 0.5);
}

}

Calculator::Calculator() : impl_(new Impl()) {}

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

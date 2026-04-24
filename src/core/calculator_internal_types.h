#ifndef CALCULATOR_INTERNAL_TYPES_H
#define CALCULATOR_INTERNAL_TYPES_H

#include "calculator.h"

#include "matrix.h"
#include "script_ast.h"

#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

struct Rational {
    long long numerator = 0;
    long long denominator = 1;

    Rational() = default;

    Rational(long long num, long long den) : numerator(num), denominator(den) {
        normalize();
    }

    void normalize() {
        if (denominator == 0) {
            throw std::runtime_error("division by zero");
        }
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }

        const long long divisor = std::gcd(numerator, denominator);
        numerator /= divisor;
        denominator /= divisor;
    }
};

struct StoredValue {
    bool is_matrix = false;
    bool is_string = false;
    bool has_symbolic_text = false;
    bool has_precise_decimal_text = false;
    bool exact = false;
    Rational rational;
    double decimal = 0.0;
    std::string string_value;
    std::string symbolic_text;
    std::string precise_decimal_text;
    matrix::Matrix matrix;
};

struct CustomFunction {
    std::string parameter_name;
    std::string expression;
};

struct ScriptFunction {
    std::vector<std::string> parameter_names;
    std::shared_ptr<const script::BlockStatement> body;
};

struct Calculator::Impl {
    std::map<std::string, StoredValue> variables;
    std::map<std::string, CustomFunction> functions;
    std::map<std::string, ScriptFunction> script_functions;
    std::vector<std::map<std::string, StoredValue>> local_scopes;
    bool symbolic_constants_mode = false;
    bool hex_prefix_mode = false;
    bool hex_uppercase_mode = true;
};

std::string render_script_block(const script::BlockStatement& block, int indent);

#endif

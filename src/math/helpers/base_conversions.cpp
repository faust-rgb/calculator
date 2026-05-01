#include "base_conversions.h"
#include <stdexcept>

int digit_value(char ch) {
    if (ch >= '0' && ch <= '9') return ch - '0';
    if (ch >= 'a' && ch <= 'f') return 10 + (ch - 'a');
    if (ch >= 'A' && ch <= 'F') return 10 + (ch - 'A');
    return -1;
}

bool prefixed_base(char prefix, int* base) {
    if (prefix == 'b' || prefix == 'B') {
        *base = 2; return true;
    }
    if (prefix == 'o' || prefix == 'O') {
        *base = 8; return true;
    }
    if (prefix == 'x' || prefix == 'X') {
        *base = 16; return true;
    }
    return false;
}

long long parse_prefixed_integer_token(const std::string& token) {
    if (token.size() < 3 || token[0] != '0') {
        throw std::runtime_error("invalid prefixed integer literal");
    }

    int base = 10;
    if (!prefixed_base(token[1], &base)) {
        throw std::runtime_error("invalid prefixed integer literal");
    }

    long long value = 0;
    bool has_digit = false;
    for (std::size_t i = 2; i < token.size(); ++i) {
        const int digit = digit_value(token[i]);
        if (digit < 0 || digit >= base) {
            throw std::runtime_error("invalid digit in prefixed integer literal");
        }
        has_digit = true;
        value = value * static_cast<long long>(base) + static_cast<long long>(digit);
    }

    if (!has_digit) {
        throw std::runtime_error("prefixed integer literal requires digits");
    }

    return value;
}

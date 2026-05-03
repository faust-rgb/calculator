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

std::string convert_to_base(long long value, int base, bool uppercase, bool prefix) {
    if (base < 2 || base > 16) {
        throw std::runtime_error("base must be in the range [2, 16]");
    }

    static const char upper_digits[] = "0123456789ABCDEF";
    static const char lower_digits[] = "0123456789abcdef";
    const char* digits = uppercase ? upper_digits : lower_digits;

    if (value == 0) {
        if (prefix) {
            if (base == 2) return "0b0";
            if (base == 8) return "0o0";
            if (base == 16) return "0x0";
        }
        return "0";
    }

    bool negative = value < 0;
    unsigned long long current = negative
                                     ? static_cast<unsigned long long>(-(value + 1)) + 1ULL
                                     : static_cast<unsigned long long>(value);

    std::string reversed;
    while (current > 0) {
        reversed.push_back(digits[current % static_cast<unsigned long long>(base)]);
        current /= static_cast<unsigned long long>(base);
    }

    std::string output;
    if (negative) {
        output.push_back('-');
    }
    if (prefix) {
        if (base == 2) output += "0b";
        else if (base == 8) output += "0o";
        else if (base == 16) output += "0x";
    }
    for (std::size_t i = reversed.size(); i > 0; --i) {
        output.push_back(reversed[i - 1]);
    }
    return output;
}

#include "bitwise_helpers.h"
#include <stdexcept>

constexpr unsigned kProgrammerBitWidth = 64;

std::uint64_t to_unsigned_bits(long long value) {
    return static_cast<std::uint64_t>(value);
}

long long from_unsigned_bits(std::uint64_t value) {
    return static_cast<long long>(value);
}

unsigned normalize_rotation_count(long long count) {
    if (count < 0) throw std::runtime_error("rotate count cannot be negative");
    return static_cast<unsigned>(count % static_cast<long long>(kProgrammerBitWidth));
}

std::uint64_t rotate_left_bits(std::uint64_t value, unsigned count) {
    if (count == 0) return value;
    return (value << count) | (value >> (kProgrammerBitWidth - count));
}

std::uint64_t rotate_right_bits(std::uint64_t value, unsigned count) {
    if (count == 0) return value;
    return (value >> count) | (value << (kProgrammerBitWidth - count));
}

int popcount_bits(std::uint64_t value) {
    int count = 0;
    while (value != 0) {
        value &= (value - 1);
        ++count;
    }
    return count;
}

int bit_length_bits(std::uint64_t value) {
    int length = 0;
    while (value != 0) {
        ++length;
        value >>= 1;
    }
    return length;
}

int trailing_zero_count_bits(std::uint64_t value) {
    if (value == 0) return static_cast<int>(kProgrammerBitWidth);
    int count = 0;
    while ((value & 1ULL) == 0ULL) {
        ++count;
        value >>= 1;
    }
    return count;
}

int leading_zero_count_bits(std::uint64_t value) {
    if (value == 0) return static_cast<int>(kProgrammerBitWidth);
    int count = 0;
    std::uint64_t mask = 1ULL << (kProgrammerBitWidth - 1);
    while ((value & mask) == 0ULL) {
        ++count;
        mask >>= 1;
    }
    return count;
}

int parity_bits(std::uint64_t value) {
    return popcount_bits(value) % 2;
}

std::uint64_t reverse_bits(std::uint64_t value) {
    std::uint64_t reversed = 0ULL;
    for (unsigned i = 0; i < kProgrammerBitWidth; ++i) {
        reversed = (reversed << 1) | (value & 1ULL);
        value >>= 1;
    }
    return reversed;
}

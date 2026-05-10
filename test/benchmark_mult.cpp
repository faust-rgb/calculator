#include <iostream>
#include <vector>
#include <chrono>
#include "precise/precise_decimal.h"

// Accessing internal functions by declaring them
namespace precise {
    std::vector<uint32_t> multiply_bigint_naive(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs);
    std::vector<uint32_t> multiply_bigint_karatsuba(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs);
    std::vector<uint32_t> multiply_bigint_ntt(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs);
}

// Since they are in anonymous namespace in precise_decimal.cpp, we can't easily access them
// unless we modify the file or use the public operator*.
// But wait, the grep showed they are NOT in anonymous namespace!
// Let's check src/precise/precise_decimal.cpp again for the namespace.

/*
File: src/precise/precise_decimal.cpp
L34- namespace {
...
L644: std::vector<uint32_t> multiply_bigint(const std::vector<uint32_t>& lhs, const std::vector<uint32_t>& rhs) {
*/
// Actually, they seem to be OUTSIDE the anonymous namespace based on the grep results.
// Let's verify.

int main() {
    std::vector<uint32_t> a(500, 123456789);
    std::vector<uint32_t> b(500, 987654321);

    auto start = std::chrono::high_resolution_clock::now();
    for(int i=0; i<10; ++i) {
        PreciseDecimal da, db;
        // We can't easily set data, but we can multiply large numbers
    }
    
    // Let's use strings to create large numbers
    std::string s1(2000, '9');
    std::string s2(2000, '9');
    PreciseDecimal p1(s1), p2(s2);
    
    start = std::chrono::high_resolution_clock::now();
    for(int i=0; i<100; ++i) {
        PreciseDecimal p3 = p1 * p2;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time for 100 mults (2000 digits): " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

    std::string s3(10000, '9');
    std::string s4(10000, '9');
    PreciseDecimal p4(s3), p5(s4);
    
    start = std::chrono::high_resolution_clock::now();
    for(int i=0; i<10; ++i) {
        PreciseDecimal p6 = p4 * p5;
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Time for 10 mults (10000 digits): " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

    return 0;
}
